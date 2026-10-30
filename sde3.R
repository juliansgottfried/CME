library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

setwd("~/Desktop/CME")

assoc <- read.table("cme_v2/associations_pre.txt", sep = "\t", header = T)
O2 <- assoc |> select(aid, host_bin, cme_id) |> distinct() |> 
    group_by(aid, cme_id) |> summarize(O2 = n()) |> ungroup()

cmes <- unique(assoc$cme_id)

hosts <- read.table("cme_v2/hosts.txt", sep = "\t", header = T) |> 
    select(aid, host_bin = bin, family)
assoc <- left_join(assoc, hosts, by = c("aid", "host_bin"))

abun <- read.table("cme_v2/host_abundance.txt", sep = "\t", header = T) |> 
    rename(host_bin = bin)
abun <- left_join(abun, hosts, by = c("aid", "host_bin"))
O1 <- abun |> filter(early > 0) |> select(aid, host_bin, family) |> distinct()

abun <- abun |> select(aid, host_bin, family, x = early) |> mutate(x = x / 100) |> 
    pivot_wider(id_cols = aid, names_from = "host_bin", values_from = "x")
abun[is.na(abun)] <- 0
minval <- min(abun[, 2:ncol(abun)][abun[, 2:ncol(abun)] > 0])
abun[, 2:ncol(abun)] <- abun[, 2:ncol(abun)] + minval / 10
rows <- rowSums(abun[, 2:ncol(abun)])
abun[, 2:ncol(abun)] <- apply(abun[, 2:ncol(abun)], 2, \(x) x / rows)
abun <- abun |> pivot_longer(cols = -aid, names_to = "host_bin", values_to = "x")
abun <- left_join(abun, hosts, by = c("aid", "host_bin"))

O1s <- NULL
abuns <- NULL
for (i in 1:length(cmes)) {
    print(i)
    fams <- assoc %>% filter(cme_id == cmes[i]) %>% pull(family) %>% unique
    
    O1s <- rbind(O1s, O1 |> filter(family %in% fams) |> 
                     group_by(aid) |> summarize(O1 = n()) |> 
                     mutate(cme_id = cmes[i]))
    
    abuns <- rbind(abuns, abun |> filter(family %in% fams) |> group_by(aid) |> 
        summarize(X = sum(x)) |> mutate(cme_id = cmes[i]))
}


O1s <- O1s |> 
    pivot_wider(id_cols = aid, names_from = "cme_id", values_from = "O1") |> 
    arrange(aid) |> select(-aid)
O1s[is.na(O1s)] <- 0
O2 <- O2 |> 
    pivot_wider(id_cols = aid, names_from = "cme_id", values_from = "O2") |> 
    arrange(aid) |> select(-aid)
O2[is.na(O2)] <- 0
abuns <- abuns |>
    pivot_wider(id_cols = aid, names_from = "cme_id", values_from = "X") |> 
    arrange(aid) |> select(-aid)
abuns[is.na(abuns)] <- 0
abun <- abun |>
    pivot_wider(id_cols = aid, names_from = "host_bin", values_from = "x") |> 
    arrange(aid) |> select(-aid)

O1s <- O1s[, order(colnames(O1s))]
O2 <- O2[, order(colnames(O2))]
abuns <- abuns[, order(colnames(abuns))]
abun <- abun[, order(colnames(abun))]

O1s[O1s == 0] <- 1
abuns[abuns == 0] <- 1e-10

kappa_dat <- list(n_sub = nrow(abun),
                  n_sp = ncol(abun),
                  x = as.matrix(abun),
                  p = as.numeric(colMeans(abun)))

kappa_mod <- cmdstan_model("kappa.stan")
kappa_fit <- kappa_mod$sample(
    data = kappa_dat,
    chains = 4,
    parallel_chains = 4)

kappa <- kappa_fit$draws(format = "df") |> pull(kappa) |> mean()

R0_mod <- cmdstan_model("R0.stan")

fit_R0 <- function(cme) {
    R0_dat <- list(n = nrow(O1s),
                   O1 = as.numeric(O1s[[cme]]),
                   O2 = as.numeric(O2[[cme]]),
                   X = as.numeric(abuns[[cme]]))
    R0_fit <- R0_mod$sample(
        data = R0_dat,
        chains = 4,
        parallel_chains = 4,
        adapt_delta = 0.99)
    
    # mcmc_trace(R0_fit$draws(format = "df"))
    
    R0_fit$draws(format = "df") |> pull(R0)
}

ret <- sapply(cmes, fit_R0)
ret <- data.frame(ret)
write_csv(ret, "ret.csv")
ret <- read_csv("ret.csv")

ret |> pivot_longer(everything()) |> 
    group_by(name) |> summarize(mean = mean(value), sd = sd(value)) |> 
    arrange(-mean) |> mutate(cme = row_number()) |> 
    ggplot(aes(x = cme, y = mean)) +
    geom_point(alpha = 0.5, shape = 3) +
    geom_linerange(aes(ymin = mean - sd, ymax = mean + sd),
                   alpha = 0.1) +
    geom_hline(yintercept = 1, linetype = 2, color = "red") +
    ylim(0, NA) + ylab("R0") +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))

tester <- "c4"
dat <- list(n = nrow(O1s),
            O1 = as.numeric(O1s[[tester]]),
            O2 = as.numeric(O2[[tester]]),
            X = as.numeric(abuns[[tester]]))

R0 <- mean(ret[[tester]])
ystar = 1 - 1 / (dat$X * R0)
soft_max = log(exp(100 * ystar) + exp(100 * 0.001)) / 100
avg <- floor(dat$O1 * soft_max)
sd <- floor(sqrt(dat$O1 * soft_max * (1 - soft_max)))
coll <- c()
for (i in 1:dat$n) {
    coll <- c(coll, max(avg[i] - sd[i], 0):(avg[i] + sd[i]))
}
sim <- table(coll)
real <- table(dat$O2)
data.frame(count = as.numeric(names(real)), freq = as.numeric(real / sum(real))) |> 
    ggplot(aes(x = count, y = freq)) +
    geom_col(width = 1, color = "black", fill = "white") +
    ylim(0, 1) + xlim(-0.5, 30) +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))
data.frame(count = as.numeric(names(sim)), freq = as.numeric(sim / sum(sim))) |> 
    ggplot(aes(x = count, y = freq)) +
    geom_col(width = 1, color = "black", fill = "white") +
    ylim(0, 1) + xlim(-0.5, 30) +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))

R0
1 / colMeans(abuns)[[tester]]





assoc <- read.table("cme_v2/associations_pre.txt", sep = "\t", header = T)
assoc <- assoc |> group_by(aid, host_bin, cme_id, host_segment) |> 
    slice_head(n = 1) |> ungroup(host_segment) |> 
    summarize(R = sum(reads), M = n(), mu = xcov[1], var = xvar[1]) |> 
    ungroup() |> mutate(var = ifelse(var <= mu, mu + 1e-10, var)) |> 
    arrange(aid, cme_id)

M <- assoc |> group_by(host_bin, cme_id) |> summarize(M = median(M)) |> 
    ungroup() |> pivot_wider(id_cols = host_bin, names_from = cme_id, values_from = M) |> 
    arrange(host_bin) |> select(-host_bin)

xcov <- assoc |> select(aid, host_bin, mu, var) |> distinct() |> 
    pivot_wider(id_cols = aid, names_from = host_bin, values_from = mu) |> 
    arrange(aid) |> select(-aid)
xcov <- xcov[, order(colnames(xcov))]

xvar <- assoc |> select(aid, host_bin, mu, var) |> distinct() |> 
    pivot_wider(id_cols = aid, names_from = host_bin, values_from = var) |> 
    arrange(aid) |> select(-aid)
xvar <- xvar[, order(colnames(xvar))]

makeR <- function(cme) {
    R <- assoc |> mutate(R = ifelse(cme_id == cme, R, 0)) |> 
        group_by(aid, host_bin) |> summarize(R = sum(R)) |> ungroup() |> 
        pivot_wider(id_cols = aid, names_from = host_bin, values_from = R) |> 
        arrange(aid) |> select(-aid)
    R[is.na(R)] <- 0
    R <- R[, order(colnames(R))]
    R
}

X <- abuns |> arrange(cme_id, aid) |> select(-aid)

cme <- "c4"

dat_pre <- list(n_aid = nrow(xcov),
            n_bin = ncol(xcov),
            M = M[[cme]],
            xcov = xcov,
            xvar = xvar,
            X = X |> filter(cme_id == cme) |> pull(X),
            R = makeR(cme))

M <- rep(dat_pre$M, each = dat_pre$n_aid)
xcov <- stack(dat_pre$xcov)$values
xvar <- stack(dat_pre$xvar)$values
X <- rep(dat$X, times = dat_pre$n_bin)
R <- stack(dat_pre$R)$values

keep <- which(!is.na(xcov) & !is.na(M) & xcov > 0 & M > 0)
Mfilt <- M[keep]
xcovfilt <- xcov[keep]
xvarfilt <- xvar[keep]
Xfilt <- X[keep]
Rfilt <- R[keep]

dat <- list(n = length(Mfilt),
            M = Mfilt, 
            xcov = xcovfilt,
            xvar = xvarfilt,
            X = Xfilt,
            R = Rfilt)

R0_mod <- cmdstan_model("R0.stan")

R0_fit <- R0_mod$sample(
    data = dat,
    chains = 4,
    parallel_chains = 4,
    adapt_delta = 0.9)

R0 <- R0_fit$draws(format = "df") |> pull(R0) |> mean()

zstar = 1 - 1 / (dat$X * R0)
smax = log(exp(100 * zstar) + exp(100 * 0.001) )/ 100

R <- rnbinom(dat$n, mu = dat$M * dat$xcov * smax, size = dat$M * dat$xcov ^ 2 / (dat$xvar - dat$xcov));
sim <- table(R) / sum(R)
data.frame(count = as.numeric(names(sim)), freq = as.numeric(sim)) |> 
    ggplot(aes(x = count, y = freq)) +
    geom_col(width = 1, color = "black", fill = "white") +
    ylim(0, 1) +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))
