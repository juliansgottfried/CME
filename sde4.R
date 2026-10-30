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
         aid = all$aid)
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

plotsim <- function(cme, ret) {
    # cme <- "c4"
    estR0 <- mean(ret[[cme]])
    testdat <- create_dat(cme)
    
    zstar = 1 - 1 / (testdat$X * testR0)
    smax = log(exp(100 * zstar) + exp(100 * 0.001) )/ 100
    # R <- rnbinom(testdat$n, mu = testdat$M * testdat$xcov * smax, 
    #              size = testdat$M * testdat$xcov ^ 2 / (testdat$xvar - testdat$xcov))
    
    R <- round(testdat$M * testdat$xcov * smax)
    
    forplot <- data.frame(aid = testdat$aid, sim = R, real = testdat$R) |> 
        group_by(aid) |> summarize(sim = sum(sim > 0), real = sum(real > 0))
    simtab <- forplot$sim |> table()
    realtab <- forplot$real |> table()
    
    rl <- max(c(forplot$sim, forplot$real))
    
    p1 <- data.frame(count = as.numeric(names(realtab)), freq = as.numeric(realtab / sum(realtab))) |> 
        ggplot(aes(x = count, y = freq)) +
        geom_col(width = 1, color = "black", fill = "white") +
        ylim(0, 1) + xlim(-0.5, rl + 0.5) +
        theme_classic() + xlab("data") +
        theme(text = element_text(size = 12, family = "mono"))
    p2 <- data.frame(count = as.numeric(names(simtab)), freq = as.numeric(simtab / sum(simtab))) |> 
        ggplot(aes(x = count, y = freq)) +
        geom_col(width = 1, color = "black", fill = "white") +
        ylim(0, 1) + xlim(-0.5, rl + 0.5) +
        theme_classic() + xlab("posterior mean")
    
    P <- data.frame(aid = testdat$aid, X = testdat$X) |> 
        group_by(aid) |> slice_head(n = 1) |> 
        pull(X) |> mean()
    
    list(p1 / p2 + plot_annotation(title = cme) &
             theme(text = element_text(size = 12, family = "mono")),
        estR0,
        1 / P)
}

mod <- cmdstan_model("R0.stan")

tester <- "c4"
ret <- data.frame(sapply(tester, fit_R0))
# write_csv(ret, "ret.csv")
# ret <- read_csv("ret.csv")

out <- plotsim(tester, ret)
out[[1]]
print(paste0("R0 = ", round(out[[2]], 2)))
print(paste0("1/pT = ", round(out[[3]], 2)))

