library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

setwd("~/Desktop/CME")

hosttax <- read.table("cme_v2/hosts.txt", sep = "\t", header = T)

create_dat <- function(cme) {
    filt <- assoc |> filter(cme_id == cme)
    bins <- unique(filt$host_bin)
    # otherbins <- setdiff(unique(hosttax$bin[hosttax$family %in% unique(hosttax$family[hosttax$bin %in% bins])]), bins)
    # otherfilt <- assoc |> filter(host_bin %in% otherbins) |> 
    #     mutate(cme_id = "c4", R = 0, M = 1)
    # bins <- c(bins, otherbins)
    # filt <- rbind(filt, otherfilt)
    
    M <- filt |> group_by(host_bin) |>
        summarize(M = max(M)) |> ungroup() 
    
    xs <- x |> filter(host_bin %in% bins)
    
    xcov <- assoc |> select(aid, host_bin, mu) |> 
        distinct() |> filter(host_bin %in% bins) |> 
        pivot_wider(id_cols = aid, names_from = host_bin, values_from = mu) |> 
        pivot_longer(-aid, names_to = "host_bin", values_to = "xcov")
    xvar <- assoc |> select(aid, host_bin, var) |> 
        distinct() |> filter(host_bin %in% bins) |> 
        pivot_wider(id_cols = aid, names_from = host_bin, values_from = var) |> 
        pivot_longer(-aid, names_to = "host_bin", values_to = "xvar")
    
    R <- filt |> group_by(aid, host_bin) |> summarize(R = sum(R)) |> ungroup() |> 
        pivot_wider(id_cols = aid, names_from = host_bin, values_from = R) |> 
        pivot_longer(-aid, names_to = "host_bin", values_to = "R")
    
    p <- colMeans(abundance[, colnames(abundance) %in% bins])
    p <- data.frame(host_bin = names(p), p = p)
    
    all <- full_join(xcov, xvar, by = c("aid", "host_bin")) |> 
        full_join(M, by = "host_bin") |> 
        full_join(R, by = c("aid", "host_bin")) |> 
        full_join(x, by = c("aid", "host_bin")) |> 
        full_join(p, by = c("host_bin")) |> 
        mutate(R = ifelse(is.na(R), 0, R)) |> 
        filter(!is.na(xcov)) |> filter(xcov > 0)
    
    all <- all |> 
        arrange(host_bin) |> 
        mutate(host_bin = as.numeric(as.factor(host_bin)))
    
    p <- unique(all$p)
    
    all <- all |> arrange(aid) |> mutate(aid = as.numeric(as.factor(aid)))
    
    # all <- all |> filter(aid %in% 1:10)
    
    X <- all |> group_by(aid) |> summarize(X = sum(x)) |> pull(X)
    
    list(n = nrow(all),
         nbin = length(p),
         naid = length(unique(all$aid)),
         bin = all$host_bin,
         p = p,
         aid = all$aid,
         Y = Z * X,
         x = all$x,
         mu = all$M * all$xcov,
         phi = all$M * all$xcov ^ 2 / (all$xvar - all$xcov),
         X = X,
         # R = all$R,
         kappa = kappa,
         P = sum(p),
         f = c4pi,
         stat = all |> group_by(aid) |> summarize(stat = sum(R > 0)) |> pull(stat),
         s = all |> group_by(aid) |> summarize(n = n()) |> pull(n)
         )
}

fit_R0 <- function(cme) {
    cme <- "c4"
    print(which(cme == cmes))
    dat <- create_dat(cme)

    fit <- mod$sample(
        data = dat,
        chains = 4,
        iter_warmup = 1000,
        iter_sampling = 1000,
        parallel_chains = 4, 
        fixed_param = T)
        # init = list(
        #     list(R0 = 5, a = 10),
        #     list(R0 = 4, a = 15),
        #     list(R0 = 3, a = 20),
        #     list(R0 = 2, a = 25)))
        # init = list(
        #     list(R0 = 2, a = 15, eta = 0.033),
        #     list(R0 = 3, a = 20, eta = 0.067),
        #     list(R0 = 4, a = 25, eta = 0.1),
        #     list(R0 = 5, a = 30, eta = 0.133)))
    
    mcmc_trace(fit$draws(format = "df", variables = "R0")) + theme_classic()
    mcmc_trace(fit$draws(format = "df", variables = "eta")) + theme_classic()
    mcmc_trace(fit$draws(format = "df", variables = "a")) + theme_classic()
    
    # data.frame(fit$draws(format = "df")) |> select(R0, eta) |>
    #     pivot_longer(everything(), names_to = "par", values_to = "posterior") |> 
    #     mutate(prior = c(rgamma(4000, shape = 100, rate = 1),
    #                      rgamma(4000, shape = 1, rate = 1))[rep(1:4000, each = 2) + rep(c(0, 4000), 4000)]) |> 
    #     pivot_longer(-par) |> 
    #     ggplot(aes(x = value, fill = name)) +
    #     geom_density(color = "white", alpha = 0.8) + 
    #     scale_fill_manual(values = c("#4287f5", "#ea9afc"), name = "") +
    #     facet_wrap(~par, scales = "free") +
    #     theme_classic()
    
    realtab <- data.frame(R = dat$R, aid = dat$aid) |> 
        group_by(aid) |> summarize(X = sum(R > 0)) |>
        mutate(X = factor(X, levels = 0:22)) |> 
        pull(X) |> table()
    plot(realtab / sum(realtab))
    
    statsimtab <- fit$draws(format = "df") |> 
        select(starts_with("stat_sim[")) |>
        t() |> as.data.frame() |> 
        apply(1, \(x) table(factor(x, levels = 0:22))) |> 
        rowSums()
    plot(statsimtab / sum(statsimtab), type = "h")
    
    data.frame(data = as.numeric(realtab / sum(realtab)),
               sim = as.numeric(statsimtab / sum(statsimtab)),
               bin = 0:22) |> 
        pivot_longer(-bin) |> 
        mutate(name = factor(name, levels = c("data", "sim"))) |> 
        arrange(desc(name)) |> 
        ggplot(aes(x = bin, y = value, fill = name)) + 
        geom_col(position = "identity", alpha = 1) +
        scale_fill_manual(name = "", values = c("#4287f5", "black")) +
        theme_classic() +
        theme(text = element_text(size = 12, family = "mono"))
    
    realR <- table(dat$R)
    realR <- realR[realR > 1]
    plot(log(realR / min(realR)))
    
    rsims <- fit$draws(format = "df") |> 
        select(starts_with("R_sim[")) |>
        t() |> as.data.frame()
    simtabs <- apply(rsims, 2, \(x) {
        table(factor(x, levels = 0:36))
    })
    plot(log(rowSums(simtabs) / min(rowSums(simtabs))), type = "h")
    
    realtab <- data.frame(R = dat$R, aid = dat$aid) |> 
        group_by(aid) |> summarize(X = sum(R > 0)) |>
        mutate(X = factor(X, levels = 0:22)) |> 
        pull(X) |> table()
    plot(realtab / sum(realtab))
    
    transfrsims <- rsims |> cbind(aid = dat$aid) |>
        group_by(aid) |> summarize_at(vars(starts_with("V")), \(x) sum(x > 0)) |>
        select(-aid) |> as.matrix()
    transsimtabs <- apply(transfrsims, 2, \(x) {
        table(factor(x, levels = 0:22))
        })
    plot(rowSums(transsimtabs) / sum(transsimtabs), type = "h")
}

mod <- cmdstan_model("mvn.stan")

ret <- data.frame(sapply("c4", fit_R0))
write_csv(ret, "ret.csv")
ret <- read_csv("ret.csv")



# sims
nR0 <- 1
naid <- 100
keep <- which(dat$p > 0.005)
nbin <- length(keep)
R0s <- c(0, 2 ^ (0:10))[nR0]
sim <- data.frame(matrix(0, ncol = naid + 3, nrow = nbin * nR0 * 2))
colnames(sim)[(naid + 1):(naid + 3)] <- c("state", "bin", "R0")
for (i in 1:nR0) {
    idxx <- (1:nbin) + nbin * 2 * (i - 1)
    idxy <- ((nbin + 1):(2 * nbin)) + nbin * 2 * (i - 1)

    sim[idxx, 1:naid] <- t(read_csv(paste0("simulations/x_", i, ".csv"), col_names = F))
    sim$state[idxx] <- "x"
    sim$bin[idxx] <- keep
    sim$R0[idxx] <- R0s[i]
    
    sim[idxy, 1:naid] <- t(read_csv(paste0("simulations/y_", i, ".csv"), col_names = F))
    sim$state[idxy] <- "y"
    sim$bin[idxy] <- keep
    sim$R0[idxy]<- R0s[i]
}
simproc <- sim |> select(-R0) |> 
    pivot_longer(starts_with("X"), names_to = "aid", values_to = "val") |> 
    mutate(aid = as.numeric(str_remove(aid, "X"))) |>
    # pivot_wider(id_cols = c(aid, bin, R0), names_from = "state", values_from = "val") |> 
    pivot_wider(id_cols = c(aid, bin), names_from = "state", values_from = "val") |> 
    filter(x > min(dat$x))
datdata <- data.frame(aid = dat$aid, bin = dat$bin, mu = dat$mu, phi = dat$phi) |> 
    distinct() |> filter(bin %in% keep)
simproc <- simproc |> left_join(datdata, by = c("bin", "aid")) |> na.omit()
simproc$z <- simproc$y / simproc$x
simproc$z <- 1 / (1 + exp(-20 * (simproc$z - 0.5)))
simproc$R <- with(simproc, rnbinom(nrow(simproc), mu = mu * z, size = phi))
simstat <- simproc |> 
    # group_by(aid, R0) |> 
    group_by(aid) |> 
    summarize(stat = sum(R > 0)) |> ungroup()
# simstat |> filter(R0 == 1) |> pull(stat) |> table() |> plot()
simstat |> pull(stat) |> table() |> plot()
