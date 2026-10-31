library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(patchwork)

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
         aid = all$aid,
         bin = all$host_bin)
}

fit_R0 <- function(cme) {
    dat <- create_dat(cme)
    dat$aid <- NULL
    fit <- mod$sample(
        data = dat,
        chains = 4,
        parallel_chains = 4,
        adapt_delta = 0.99)
    
    # mcmc_trace(R0_fit$draws(format = "df"))
    fit$draws(format = "df") |> pull(R0)
}

dosim <- function(simXi, testdat, testR0) {
    zstar = 1 - 1 / (simXi * testR0)
    smax = log(exp(100 * zstar) + exp(100 * 0.001) ) / 100
    
    data.frame(sapply(1:10, \(x) rnbinom(length(simXi), mu = testdat$M * testdat$xcov * smax,
            size = testdat$M * testdat$xcov ^ 2 / (testdat$xvar - testdat$xcov)))) |> 
        mutate(aid = testdat$aid)
}

plotsim <- function(cme) {
    # cme <- "c4"
    estR0 <- mean(ret[[cme]])
    testdat <- create_dat(cme)
    
    P <- data.frame(aid = testdat$aid, X = testdat$X) |> 
        distinct() |> pull(X) |> mean()
    
    simX <- matrix(rbeta(42 * 10, kappa * P, kappa * (1 - P)), nrow = 42, ncol = 10)
    simX <- simX |> data.frame() |> mutate(aid = unique(testdat$aid)) |> 
        right_join(data.frame(aid = testdat$aid), by = "aid")
    
    sims <- bind_rows(apply(simX[, 1:10], 2, \(x) dosim(x, testdat, testR0))) |> mutate(id1 = rep(1:10, each = n() / 10))

    realtab <- data.frame(aid = testdat$aid, real = testdat$R) |> 
        group_by(aid) |> summarize(real = sum(real >= 1)) |> 
        pull(real) |> table()
    realtab <- realtab / sum(realtab)
    
    simdf <- sims |> pivot_longer(cols = starts_with("X"), names_to = "id2", values_to = "sim") |> 
        mutate(id2 = as.numeric(str_remove(id2, "X"))) |> 
        group_by(aid, id1, id2) |> summarize(sim = sum(sim >= 1)) |> 
        ungroup() |> mutate(sim = factor(sim, levels = 0:max(sim)))
    simtab <- rep(0, length(levels(simdf$sim)))
    names(simtab) <- levels(simdf$sim)
    for (i in 1:10) {
        for (j in 1:10) {
            tab <- simdf |> filter(id1 == i, id2 == j) |> pull(sim) |> table()
            simtab <- simtab + tab
        }
    }
    simtab <- simtab / sum(simtab)

    rl <- as.numeric(names(simtab[simtab > 0.1]))
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
        theme_classic() + xlab("posterior")
    
    list(p1 / p2 + plot_annotation(title = cme) &
             theme(text = element_text(size = 12, family = "mono")),
        estR0, 1 / P)
}

mod <- cmdstan_model("R0.stan")

# ret <- data.frame(sapply(cmes, fit_R0))
# write_csv(ret, "ret.csv")
ret <- read_csv("ret.csv")

data.frame(t(apply(ret, 2, \(x) c(mean(x), var(x))))) |> 
    arrange(-X1) |> mutate(idx = row_number()) |> 
    mutate(ymin = X1 - sqrt(X2), ymax = X1 + sqrt(X2)) |> 
    mutate(ymin = ifelse(ymin < 0, 0, ymin)) |> 
    ggplot(aes(x = idx, y = X1)) +
    geom_hline(yintercept = 1, linetype = 4, color = "red") +
    geom_point(shape = 3, alpha = 0.5) +
    geom_linerange(aes(ymin = ymin, ymax  =ymax),
                   alpha = 0.2) +
    labs(x = "cme", y = "R0") +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))

tester <- "c4"
out <- plotsim(tester)
out[[1]]
print(paste0("R0 = ", round(out[[2]], 2)))
print(paste0("1 / pT = ", round(out[[3]], 2)))
