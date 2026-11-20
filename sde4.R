library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(patchwork)
library(scales)
library(poisbinom)

setwd("~/Desktop/CME")

assoc <- read.table("cme_v2/associations_pre.txt", sep = "\t", header = T)
assoc <- assoc |> group_by(aid, host_bin, cme_id, host_segment) |> 
    slice_head(n = 1) |> ungroup(host_segment) |> 
    summarize(R = sum(reads), M = n(), mu = xcov[1], var = xvar[1]) |> 
    ungroup() |> mutate(var = ifelse(var <= mu, mu + 1e-10, var))
cmes <- unique(assoc$cme_id)

x <- read.table("cme_v2/host_abundance.txt", sep = "\t", header = T) |> rename(host_bin = bin)
x <- x |> select(aid, host_bin, x = early) |> mutate(x = x / 100) |> 
    pivot_wider(id_cols = aid, names_from = "host_bin", values_from = "x")
x[is.na(x)] <- 0
minval <- min(x[, 2:ncol(x)][x[, 2:ncol(x)] > 0])
x[, 2:ncol(x)] <- x[, 2:ncol(x)] + minval / 10
rows <- rowSums(x[, 2:ncol(x)])
x[, 2:ncol(x)] <- apply(x[, 2:ncol(x)], 2, \(x) x / rows)
x <- x |> pivot_longer(cols = -aid, names_to = "host_bin", values_to = "x")

abundance <- x |> pivot_wider(id_cols = aid, names_from = "host_bin", values_from = "x") |> 
    select(-aid)
kappa_dat <- list(n_sub = nrow(abundance),
                  n_bin = ncol(abundance),
                  x = as.matrix(abundance),
                  p = as.numeric(colMeans(abundance)))
kappa_mod <- cmdstan_model("kappa.stan")
kappa_fit <- kappa_mod$sample(
    data = kappa_dat,
    chains = 4,
    parallel_chains = 4)
kappa <- kappa_fit$draws(format = "df") |> pull(kappa) |> mean()

create_dat <- function(cme) {
    filt <- assoc |> filter(cme_id == cme)
    bins <- unique(filt$host_bin)
    M <- filt |> group_by(host_bin) |>
        summarize(M = max(M)) |> ungroup() 
    
    X <- x |> filter(host_bin %in% bins) |> group_by(aid) |> 
        summarize(X = sum(x))
    
    xcov <- assoc |> select(aid, host_bin, mu) |> 
        distinct() |> filter(host_bin %in% bins) |> 
        pivot_wider(id_cols = aid, names_from = host_bin, values_from = mu) |> 
        pivot_longer(-aid, names_to = "host_bin", values_to = "xcov")
    xvar <- assoc |> select(aid, host_bin, var) |> 
        distinct() |> filter(host_bin %in% bins) |> 
        pivot_wider(id_cols = aid, names_from = host_bin, values_from = var) |> 
        pivot_longer(-aid, names_to = "host_bin", values_to = "xvar")
    
    R <- assoc |> filter(host_bin %in% bins) |> 
        mutate(R = ifelse(cme_id == cme, R, 0)) |> 
        group_by(aid, host_bin) |> summarize(R = sum(R)) |> ungroup() |> 
        pivot_wider(id_cols = aid, names_from = host_bin, values_from = R) |> 
        pivot_longer(-aid, names_to = "host_bin", values_to = "R") |> 
        mutate(R = ifelse(is.na(R), 0, R))
    
    all <- full_join(xcov, xvar, by = c("aid", "host_bin")) |> 
        full_join(M, by = "host_bin") |> 
        full_join(X, by = "aid") |> 
        full_join(R, by = c("aid", "host_bin")) |> 
        filter(xcov > 0 & !is.na(xcov))
    
    list(n = nrow(all),
         M = all$M, 
         xcov = all$xcov,
         xvar = all$xvar,
         X = all$X,
         R = all$R,
         thresh = 0.0001,
         aid = all$aid,
         bin = all$host_bin)
}

create_df <- function(cme) {
    filt <- assoc |> filter(cme_id == cme)
    bins <- unique(filt$host_bin)
    M <- filt |> group_by(host_bin) |>
        summarize(M = max(M)) |> ungroup() 
    
    X <- x |> filter(host_bin %in% bins) |> group_by(aid) |> 
        summarize(X = sum(x))
    
    xcov <- assoc |> select(aid, host_bin, mu) |> 
        distinct() |> filter(host_bin %in% bins) |> 
        pivot_wider(id_cols = aid, names_from = host_bin, values_from = mu) |> 
        pivot_longer(-aid, names_to = "host_bin", values_to = "xcov")
    xvar <- assoc |> select(aid, host_bin, var) |> 
        distinct() |> filter(host_bin %in% bins) |> 
        pivot_wider(id_cols = aid, names_from = host_bin, values_from = var) |> 
        pivot_longer(-aid, names_to = "host_bin", values_to = "xvar")
    
    R <- assoc |> filter(host_bin %in% bins) |> 
        mutate(R = ifelse(cme_id == cme, R, 0)) |> 
        group_by(aid, host_bin) |> summarize(R = sum(R)) |> ungroup() |> 
        pivot_wider(id_cols = aid, names_from = host_bin, values_from = R) |> 
        pivot_longer(-aid, names_to = "host_bin", values_to = "R") |> 
        mutate(R = ifelse(is.na(R), 0, R))
    
    full_join(xcov, xvar, by = c("aid", "host_bin")) |> 
        full_join(M, by = "host_bin") |> 
        full_join(X, by = "aid") |> 
        full_join(R, by = c("aid", "host_bin")) |> 
        filter(xcov > 0 & !is.na(xcov))
}

fit_R0 <- function(cme) {
    # cme <- "c4"
    print(which(cme == cmes))
    dat <- create_dat(cme)
    
    P <- data.frame(aid = dat$aid, X = dat$X) |>
        distinct() |> pull(X) |> mean()
    
    dat$aid <- NULL
    dat$bin <- NULL
    
    fit <- mod$sample(
        data = dat,
        chains = 4,
        iter_warmup = 2000,
        iter_sampling = 2000,
        parallel_chains = 4,
        adapt_delta = 0.9)
    
    mcmc_trace(fit$draws(format = "df"))
    fit$draws(format = "df") |> pull(R0)
    data.frame(R0 = fit$draws(format = "df") |> pull(R0)) |>
        ggplot(aes(x = R0)) +
        geom_density() + theme_classic()
}

dosim <- function(simXi, testdat, estR0) {
    zstar <- 1 - 1 / (simXi * estR0)
    # smax <- log(exp(20 * zstar) + exp(20 * 0.01) ) / 20
    # smax <- log(exp(500 * zstar) + exp(500 * 1e-8) ) / 500
    # smax <- pmax(zstar, testdat$thresh)
    smax <- pmax(zstar, exp(-500))
    
    data.frame(sapply(1:100, \(x) rnbinom(length(simXi), mu = testdat$M * testdat$xcov * smax,
            size = testdat$M * testdat$xcov ^ 2 / (testdat$xvar - testdat$xcov)))) |> 
        mutate(aid = testdat$aid, bin = testdat$bin)
}

plotsim <- function(cme, estR0) {
    # cme <- "c1"
    # estR0 <- mean(ret[[cme]])
    # estR0 <- 0
    testdat <- create_dat(cme)
    
    P <- data.frame(aid = testdat$aid, X = testdat$X) |>
        distinct() |> pull(X) |> mean()
    
    # simX <- matrix(rbeta(42 * 10, kappa * P, kappa * (1 - P)), nrow = 42, ncol = 10)
    # simX <- simX |> data.frame() |> mutate(aid = unique(testdat$aid)) |> 
    #     right_join(data.frame(aid = testdat$aid), by = "aid")
    # 
    # sims <- bind_rows(apply(simX[, 1:10], 2, \(x) dosim(x, testdat, estR0))) |> mutate(id1 = rep(1:10, each = n() / 10))
    
    sims <- dosim(testdat$X, testdat, estR0)

    realtab <- data.frame(aid = testdat$aid, bin = testdat$bin, real = testdat$R) |> 
        group_by(aid, bin) |> summarize(detect = sum(real) > 0) |> 
        summarize(real = sum(detect)) |> 
        mutate(real = factor(real, levels = 0:max(real))) |> 
        pull(real) |> table()
    realtab <- realtab / sum(realtab)
    
    simdf <- sims |> pivot_longer(cols = starts_with("X"), names_to = "id2", values_to = "sim") |> 
        mutate(id2 = as.numeric(str_remove(id2, "X"))) |> 
        group_by(aid, id2, bin) |> summarize(detect = sum(sim) > 0) |> 
        summarize(sim = sum(detect)) |> 
        ungroup() |> mutate(sim = factor(sim, levels = 0:max(sim)))
    simtab <- rep(0, length(levels(simdf$sim)))
    names(simtab) <- levels(simdf$sim)
    # for (i in 1:10) {
    #     for (j in 1:10) {
    #         tab <- simdf |> filter(id1 == i, id2 == j) |> pull(sim) |> table()
    #         simtab <- simtab + tab
    #     }
    # }
    for (j in 1:100) {
        tab <- simdf |> filter(id2 == j) |> pull(sim) |> table()
        simtab <- simtab + tab
    }
    simtab <- simtab / sum(simtab)

    rl <- as.numeric(names(simtab[simtab > 0.01]))
    rl <- rl[length(rl)]
    rl <- max(rl, as.numeric(names(realtab)))
    
    p1 <- data.frame(count = as.numeric(names(realtab)), freq = as.numeric(realtab)) |> 
        ggplot(aes(x = count, y = freq)) +
        geom_col(width = 1, color = "black", fill = "white") +
        ylim(0, 1) + xlim(-0.5, rl + 0.5) +
        theme_classic() + xlab("data")
    p2 <- data.frame(count = as.numeric(names(simtab)), freq = as.numeric(simtab)) |> 
        ggplot(aes(x = count, y = freq)) +
        geom_col(width = 1, color = "black", fill = "white") +
        ylim(0, 1) + xlim(-0.5, rl + 0.5) +
        theme_classic() + xlab("simulation")
    
    # O <- 42 * realtab
    # E <- 42 * simtab
    # 
    # O <- c(O[1], "1" = sum(O[2:length(O)]))
    # E <- c(E[1], "1" = sum(E[2:length(E)]))
    # 
    # both <- rbind(data.frame(bin = 0:1, "var" = O) |> mutate(class = "O"), 
    #       data.frame(bin = 0:1, "var" = E) |> mutate(class = "E")) |> 
    #     pivot_wider(id_cols = bin, names_from = class, values_from = var)
    # both[is.na(both)] <- 0
    # both <- both |> filter(O > 0 | E > 0)
    # both$E <- both$E + min(E) / 10
    # 
    # pval <- pchisq(sum((both$O - both$E) ^ 2 / both$E),
    #        nrow(both) - 3 - 1,
    #        lower.tail = F)
    
    list(p1 / p2 + plot_annotation(title = cme) &
             theme(text = element_text(size = 12, family = "mono")),
         estR0, 1 / P)
}

mod <- cmdstan_model("R0.stan")

ret <- data.frame(sapply(cmes, fit_R0))
write_csv(ret, "ret2.csv")
ret <- read_csv("ret.csv")

# Ps <- sapply(cmes, \(cme) {
#     testdat <- create_dat(cme)
#     data.frame(aid = testdat$aid, X = testdat$X) |>
#         distinct() |> pull(X) |> mean()
#     })
# Ps <- data.frame(Ps)
# Ps$cme <- rownames(Ps)

allplot <- data.frame(X1 = apply(ret, 2, \(x) median(x))) |> mutate(P = Ps$Ps)
# allplot |> mutate(mean = X1, var = X2) |> 
#     arrange(-mean) |> mutate(idx = row_number()) |> 
#     mutate(ymin = mean - sqrt(var), ymax = mean + sqrt(var)) |> 
#     mutate(ymin = ifelse(ymin < 0, 0, ymin)) |> 
#     ggplot(aes(x = idx, y = mean)) +
#     geom_hline(yintercept = 1, linetype = 4, color = "red") +
#     geom_point(shape = 3, alpha = 0.5) +
#     geom_linerange(aes(ymin = ymin, ymax = ymax),
#                    alpha = 0.2) +
#     labs(x = "cme", y = "R0") +
#     theme_classic() +
#     theme(text = element_text(size = 12, family = "mono"))

tester <- "c1"
out <- plotsim(tester, 2)
out[[1]]
# # print(paste0("p-val = ", signif(out[[2]], 3)))
# print(paste0("R0 = ", signif(out[[2]], 3)))
# print(paste0("1 / pT = ", signif(out[[3]], 3)))

n <- 100
P <- seq(from = 0.01, to = 1, length.out = n)
R0 <- seq(from = 0.01, to = 20, length.out = n)
# P <- 10^(seq(from = log10(0.01), to = log10(1), length.out = n))
# R0 <- 10^(seq(from = -0.5, to = log10(25), length.out = n))
record <- data.frame(P = rep(P, each = n), 
                     R0 = rep(R0, n),
                     mean = 0, median = 0, p0 = 0)
idx <- 1
for (i in 1:n) {
    print(i)
    sims <- rbeta(1000, kappa * P[i], kappa * (1 - P[i]))
    for (j in 1:n) {
        zstar <- 1 - 1 / (sims * R0[j])
        smax <- log(exp(20 * zstar) + exp(20 * 0.01) ) / 20
        # smax <- sapply(zstar, \(x) max(x, 0))
        record$mean[idx] <- mean(smax)
        record$median[idx] <- median(smax)
        record$p0[idx] <- mean(smax == 0.001)
        idx <- idx + 1
    }
}
# record |>
#     ggplot(aes(x = (P), y = (R0))) +
#     geom_tile(aes(fill = mean)) +
#     geom_point(data = allplot |> select(R0 = X1, P),
#                color = "white", shape = 3, alpha = 0.7) +
#     scale_fill_viridis_c(option = 12, name = "infected\nfraction") +
#     scale_x_log10(
#         breaks = scales::trans_breaks("log10", \(x) 10^x),
#         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#     scale_y_log10(
#         breaks = scales::trans_breaks("log10", \(x) 10^x),
#         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#     labs(y = expression(R["0"]), x = expression(p["T"])) +
#     theme_classic() +
#     theme(text = element_text(size = 12, family = "mono")) +
#     annotation_logticks()
p1 <- record |>
    ggplot(aes(x = P, y = R0)) +
    geom_tile(aes(fill = mean)) +
    geom_point(data = allplot |> select(R0 = X1, P),
               color = "#85eeff", shape = 1, alpha = 0.5) +
    geom_hline(yintercept = 1, color = "white", linetype = 3) +
    # scale_fill_viridis_c(option = 16, name = "expected\nprevalence") +
    scale_fill_gradient(low = "black", high = "white", name = "expected\nprevalence") +
    # labs(y = bquote(frac(beta,delta*eta)), x = bquote(sum(p[i]))) + 
    # labs(y = bquote(frac(R[0], eta)), x = bquote(sum(p[i]))) +
    labs(y = bquote(R[0]/eta), x = bquote(sum(p[i]))) +
    ggtitle("c. CMEs in parameter space") +
    theme_classic() +
    theme(text = element_text(size = 12, family = "serif"),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1))

cmeforplot <- "c1"
insetp <- data.frame(R0 = ret[[cmeforplot]]) |> 
    ggplot(aes(x = R0)) +
    geom_density() +
    labs(x = bquote(R[0]/eta), y = "") +
    # labs(x = bquote(frac(R[0], eta)), y = "") +
    # labs(x = bquote(frac(beta,delta*eta))) + 
    ggtitle("b. Posterior density") +
    theme_classic() +
    theme(text = element_text(size = 9, family = "serif"),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1))

testdat <- create_dat(cmeforplot)
post <- c()
for (i in 1:1000) {
    print(i)
    # estR0 <- mean(ret$c4)
    zstar = 1 - 1 / (testdat$X * ret[[cmeforplot]][i])
    # zstar = 1 - 1 / (testdat$X * estR0)
    smax = log(exp(20 * zstar) + exp(20 * 0.01) ) / 20
    post <- c(post, 
              rnbinom(length(smax), mu = testdat$M * testdat$xcov * smax,
                      size = testdat$M * testdat$xcov ^ 2 / (testdat$xvar - testdat$xcov))
    )
}
realR <- data.frame(table(testdat$R)) |> 
    mutate(x = as.numeric(Var1), y = log(Freq)) |> 
    filter(y > 0) |> select(x, y) |> mutate(class = "data")
postR <- data.frame(table(post) / 1000) |> 
    mutate(x = as.numeric(post), y = log(Freq)) |> 
    filter(y > 0) |> select(x, y) |> mutate(class = "simulation")

mainp <- rbind(realR, postR) |> 
    mutate(class = factor(class, levels = c("simulation", "data"))) |> 
    arrange(class) |> 
    ggplot(aes(x = x, y = y, fill = class)) +
    geom_col(position = "identity", alpha = 0.7) +
    scale_fill_manual(values = c("black", "#85eeff"), name = "") +
    labs(x = "# reads per subject per bin", y = "log\ncount") +
    ggtitle(paste0("a. Fitting ", cmeforplot)) +
    guides(fill = guide_legend(reverse = T)) +
    theme_classic() +
    theme(text = element_text(size = 12, family = "serif"),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1))

(mainp + inset_element(insetp, 0.5, 0.5, 1, 1)) / p1



# lik <- function(data, par) {
#     -(sum(with(data, dnbinom(x = R,
#                             mu = M * xcov * pmax(1 - 1 / (X * par), 1e-12),
#                             size = M * xcov ^ 2 / (xvar - xcov),
#                             log = T))) - par)
# }

# optimize <- function(cme, M) {
#     print(which(cmes == cme))
#     optdf <- create_df(cme) |> select(xcov, xvar, M, X, R)
#     record <- rep(0, M + 1)
#     record[1] <- optim(par = 0,
#                  lower = 0,
#                  upper = 10000,
#                  fn = lik,
#                  data = optdf,
#                  method = "Brent")$par
#     for (i in 1:M) {
#         record[i + 1] <- optim(par = 0,
#                      lower = 0,
#                      upper = 10000,
#                      fn = lik,
#                      data = optdf[sample(1:nrow(optdf), nrow(optdf), replace = T), ], 
#                      method = "Brent")$par
#     }
#     record
# }

lik <- function(data, par) {
    tol <- 1e-12
    # alpha <- 20
    
    aids <- unique(data$aid)
    
    smax <- pmax(tol, 1 - 1 / (data$X * par))
    
    # zstar <- 1 - 1 / (data$X * par)
    # smax <- matrixStats::logSumExp(c(alpha * zstar, alpha * 0.01)) / alpha
    
    phi <- with(data, M * xcov ^ 2 / (xvar - xcov))
    mu <- with(data, M * xcov * smax)
    
    data$prob <- 1 - (phi / (phi + mu)) ^ phi
    
    # loglik <- -par
    loglik <- 0
    for (i in 1:length(aids)) {
        # print(i)
        # print(loglik)
        pp <- data$prob[data$aid == aids[i]]
        pp[pp < tol] <- tol
        pp[pp > 1 - tol] <- 1 - tol
        loglik <- loglik + log(max(dpoisbinom(sum(data$R[data$aid == aids[i]] > 0),
                                      pp = pp,
                                      log_d = F), tol))
    }
    
    -loglik
}

optimize <- function(cme, M) {
    print(which(cmes == cme))
    optdf <- create_df(cme)
    
    optim(par = 0,
          lower = 1.427928,
          upper = 100,
          fn = lik,
          data = optdf,
          method = "Brent")
    
    record <- rep(0, M + 1)
    record[1] <- optim(par = 0,
                       lower = 0,
                       upper = 20,
                       fn = lik,
                       data = optdf,
                       method = "Brent")$par
    for (i in 1:M) {
        record[i + 1] <- optim(par = 0,
                               lower = 0,
                               upper = 10,
                               fn = lik,
                               data = optdf[sample(1:nrow(optdf), nrow(optdf), replace = T), ], 
                               method = "Brent")$par
    }
    record
}

M <- 100
# out <- data.frame(sapply(cmes, \(x) optimize(x, M)))
# write_csv(out, "optim3.csv")
# out <- read_csv("optim3.csv")

ests <- data.frame(est = t(out[1, ]))
ests$cme <- rownames(ests)
plotout <- out[2:nrow(out), ] |> 
    pivot_longer(everything(), names_to = "cme") |> 
    left_join(Ps, by = "cme") |> 
    left_join(ests, by = "cme") |> 
    mutate(delta = value - est) |> 
    group_by(cme) |> 
    summarize(l = quantile(delta, 0.025),
              u = quantile(delta, 0.975),
              est = est[1],
              Ps = Ps[1]) |> 
    mutate(ymin = pmax(est - u, 0), ymax = est - l)

plotout |> 
    arrange(-est) |> mutate(idx = row_number()) |> 
    ggplot(aes(x = idx, y = est)) +
    geom_point(alpha = 0.8, shape = 1) +
    geom_hline(yintercept = 1, linetype = 2, color = "red") +
    geom_linerange(aes(ymin = ymin, ymax = ymax), alpha = 0.25) +
    theme_classic()

plotout |>
    ggplot(aes(x = Ps, y = est)) +
    geom_point() +
    theme_classic()

bootest <- apply(out[2:nrow(out), ], 2, median)
tester <- "c1"
plout <- plotsim(tester, bootest[tester])
plout[[1]]

testdat <- create_dat(tester)
sims <- dosim(testdat$X, testdat, bootest[tester])
hist(testdat$R)
hist(sims[, 1])







cme <- "c8"
tol <- 1e-9
# tol <- exp(-500)
data <- create_df(cme)
aids <- unique(data$aid)
phi <- with(data, M * xcov ^ 2 / (xvar - xcov))
pars <- seq(from = 0, to = 2.6, length.out = 1000)
Rdat <- sapply(1:42, \(i) sum(data$R[data$aid == aids[i]] > 0))
idcs <- sapply(1:42, \(i) data$aid == aids[i])
P <- data.frame(aid = data$aid, X = data$X) |>
    distinct() |> pull(X) |> mean()
lik <- function(par, tol, M, xcov, X, phi, Rdat, idcs) {
    # X <- data$X
    # xcov <- data$xcov
    # M <- data$M
    # par <- 2
    zstar <- pmax(0, 1 - 1 / (X * par))
    mu <- M * xcov * zstar
    prob <- 1 - (phi / (phi + mu)) ^ phi

    l <- length(prob)
    sum(apply(rbind(idcs, 1:42), 2, \(col) {
        # col <- rbind(idcs, 1:42)[, 3]
        log(max(dpoisbinom(Rdat[col[l + 1]], pp = prob[as.logical(col[1:l])]), tol))
        }))
}
lik3 <- function(par, tol, M, xcov, X, phi, Rdat, idcs) {
    # X <- data$X
    # xcov <- data$xcov
    # M <- data$M
    # par <- 1.8
    zstar <- pmax(0, 1 - 1 / (X * par))
    mu <- M * xcov * zstar
    prob <- 1 - (phi / (phi + mu)) ^ phi
    prods <- apply(idcs, 2, \(col) {prod(1 - prob[col])})
    dpoisbinom(Rdat, pp = prods, log = T)
}

lik2 <- function(par, tol, M, xcov, X, phi, R) {
    zstar <- pmax(tol, 1 - 1 / (X * par))
    sum(dnbinom(data$R, mu = M * xcov * zstar, size = phi, log = T))
}
sim <- function(par, tol, M, xcov, X, phi, Rdat, idcs) {
    zstar <- pmax(0, 1 - 1 / (X * par))
    mu <- M * xcov * zstar
    prob <- 1 - (phi / (phi + mu)) ^ phi
    
    l <- length(prob)
    apply(rbind(idcs, 1:42), 2, \(col) {
        # col <- rbind(idcs, 1:42)[, 3]
        rpoisbinom(1, pp = prob[as.logical(col[1:l])])
    })
}
sim2 <- function(par, tol, M, xcov, X, phi, R) {
    zstar <- pmax(tol, 1 - 1 / (X * par))
    mu <- M * xcov * zstar
    rnbinom(length(R), mu = mu, size = phi)
}

# ret <- sapply(pars, \(par) lik(par, tol, data$M, data$xcov, data$X, phi, Rdat, idcs))
# data.frame(R0 = pars, loglik = ret) |> 
#     ggplot(aes(x = R0, y = loglik)) +
#     geom_point() +
#     geom_vline(xintercept = pars[which.max(ret)], color = "red", linetype = 2) +
#     theme_classic()
# 
# table(Rdat) |> plot()
# table(sim(pars[which.max(ret)], tol, data$M, data$xcov, data$X, phi, Rdat, idcs)) |> plot()

# ret2 <- sapply(pars, \(par) lik2(par, tol, data$M, data$xcov, data$X, phi, data$R))
# data.frame(R0 = pars, loglik = ret2) |> 
#     ggplot(aes(x = R0, y = loglik)) +
#     geom_point() +
#     geom_vline(xintercept = pars[which.max(ret2)], color = "red", linetype = 2) +
#     theme_classic()
# 
# (table(data$R)) |> plot()
# (table(sim2(pars[which.max(ret2)], tol, data$M, data$xcov, data$X, phi, data$R))) |> plot()

ret3 <- sapply(pars, \(par) lik3(par, tol, data$M, data$xcov, data$X, phi, sum(Rdat == 0), idcs))
data.frame(R0 = pars, loglik = ret3) |> 
    ggplot(aes(x = R0, y = loglik)) +
    geom_point() +
    geom_vline(xintercept = pars[which.max(ret3)], color = "red", linetype = 2) +
    theme_classic()
table(Rdat) |> plot()
table(sim(pars[which.max(ret3)], tol, data$M, data$xcov, data$X, phi, Rdat, idcs)) |> plot()


fit_pi <- function(cme) {
    print(which(cme == cmes))
    dat <- create_dat(cme)
    dat$aid <- NULL
    dat$bin <- NULL
    dat$X <- NULL
    dat$thresh <- NULL
    
    fit <- mod$sample(
        data = dat,
        chains = 4,
        iter_warmup = 2000,
        iter_sampling = 2000,
        parallel_chains = 4)
    
    fit$draws(format = "df") |> pull(f)
    # data.frame(pi = fit$draws(format = "df") |> pull(f)) |>
    #     ggplot(aes(x = pi)) +
    #     geom_density() + theme_classic()
}

mod <- cmdstan_model("pi.stan")
ret <- data.frame(sapply(cmes, fit_pi))
write_csv(ret, "pifit.csv")

ret <- read_csv("pifit.csv")
c4pi <- mean(ret$c4)

data.frame(pi <- colMeans(ret)) |> 
    ggplot(aes(x = pi)) +
    geom_histogram(color = "white", fill = "black", bins = 100) +
    theme_classic() + theme(text = element_text(size = 12, family = "mono"))

xbar <- x |> group_by(host_bin) |> summarize(xbar = mean(x))
unvals <- sapply(unique(x$host_bin), \(b) unique(round(x$x[x$host_bin == b], 3)))

mod <- cmdstan_model("Z.stan")
dat <- create_dat(cme)
dat$s <- as.numeric(table(dat$aid))
dat$aid <- NULL
dat$bin <- NULL
dat$X <- NULL
dat$thresh <- NULL
dat$naid <- 42
fit <- mod$sample(
    data = dat,
    chains = 4,
    iter_warmup = 2000,
    iter_sampling = 2000,
    parallel_chains = 4)
Z <- as.numeric(fit$draws(format = "df") |> select(starts_with("Z")) |> colMeans())

testcme <- "c4"
testdf <- create_df(testcme)
eta <- 76
R0 <- 0.7
f <- mean(ret[[testcme]])
xbarvec <- xbar$xbar
names(xbarvec) <- xbar$host_bin
retvec <- c()
for (iteraid in unique(testdf$aid)) {
    # iteraid <- "EBZ"
    aiddf <- testdf |> filter(aid == iteraid)
    slimx <- x |> filter(aid == iteraid)
    
    pT <- aiddf$X[1]
    
    a = R0
    b = kappa + eta - R0 * pT
    c = -kappa * pT * f
    EY = (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a)
    
    p1s <- xbarvec[aiddf$host_bin]
    samplez <- rep(0, length(p1s))
    for (i in 1:length(p1s)) {
        # i <- 1
        Ex = p1s[i]
        Ex2 = p1s[i] * (kappa * p1s[i] + 1) / (kappa + 1)
        Varx = Ex2 - Ex ^ 2
        Ey = p1s[i] * (R0 * EY + kappa * f) / (R0 * EY + kappa + eta)
        Exy = (Ey * kappa * p1s[i] + p1s[i] ^ 2 * kappa * f + Ex2 * R0 * EY) / (2 * kappa + eta + R0 * EY)
        Ey2 = (Ey * (kappa * p1s[i] * f + 1) + R0 * EY * Exy) / (kappa + eta + R0 * EY + 1)
        Vary = Ey2 - Ey^2
        Covxy = Exy - Ex*Ey
        
        tmpx <- slimx$x[slimx$host_bin == names(p1s[i])]
        condEy <- Ey + (Covxy / Varx) * (tmpx - Ex)
        condVary <- Vary - Covxy ^ 2 / Varx
        
        # samples <- rnorm(1000, condEy, condVary)
        # samples <- samples[samples < tmpx & samples > 0 & samples < 1]
        # samples <- samples[samples < tmpx]
        # samplez[i] <- mean(samples / tmpx)
        samplez[i] <- condEy
    }
    l <- nrow(aiddf)
    retvec <- c(retvec, 
                sum(rnbinom(l, mu = aiddf$M * aiddf$xcov * samplez,
                            size = aiddf$M * aiddf$xcov ^ 2 / (aiddf$xvar - aiddf$xcov)) > 0))
}
plot(table(retvec))


