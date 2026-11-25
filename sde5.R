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
        summarize(M = median(M)) |> ungroup() 
    
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
         # Y = Z * X,
         x = all$x,
         mu = all$M * all$xcov,
         phi = all$M * all$xcov ^ 2 / (all$xvar - all$xcov),
         X = X,
         R = all$R,
         kappa = kappa,
         P = sum(p),
         f = c4pi,
         stat = all |> group_by(aid) |> summarize(stat = sum(R > 0)) |> pull(stat),
         s = all |> group_by(aid) |> summarize(n = n()) |> pull(n),
         T0 = sum(all |> group_by(aid) |> summarize(stat = sum(R > 0)) |> pull(stat) == 0),
         T1 = sum(all |> group_by(aid) |> summarize(stat = sum(R > 0)) |> pull(stat) == 1)
         )
}

fit_R0 <- function(cme) {
    cme <- "c4"
    print(which(cme == cmes))
    dat <- create_dat(cme)

    fit <- mod$sample(
        data = dat[dat != "R"],
        chains = 4,
        iter_warmup = 1000,
        iter_sampling = 1000,
        parallel_chains = 4)
        # init = list(
        #     list(scR0 = 0.5, leta = 0.5),
        #     list(scR0 = 1, leta = 1),
        #     list(scR0 = 1.5, leta = 1.5),
        #     list(scR0 = 1.75, leta = 2)))
        # fixed_param = T)
    
    mcmc_trace(fit$draws(format = "df", variables = "scR0")) + theme_classic()
    mcmc_trace(fit$draws(format = "df", variables = "leta")) + theme_classic()

    GGally::ggpairs(fit$draws(format = "df") |> select(scR0, leta))
    
    realtab <- data.frame(R = dat$R, aid = dat$aid) |> 
        group_by(aid) |> summarize(X = sum(R > 0)) |>
        mutate(X = factor(X, levels = 0:22)) |> 
        pull(X) |> table()
    realdf <- data.frame(name = 0:22, m = as.numeric(realtab / sum(realtab))) |> 
        mutate(u = m, l = m, class = "data")
    
    statsimtab <- fit$draws(format = "df") |> 
        select(starts_with("stat_sim[")) |>
        apply(1, \(x) {
            tab <- table(factor(x, levels = 0:22))
            tab / sum(tab)
            }) |> t() |> data.frame()
    
    statsimtab |> 
        pivot_longer(everything()) |> 
        mutate(name = as.numeric(str_remove(name, "X"))) |>
        group_by(name) |> 
        summarize(m = mean(value),
                  u = quantile(value, 0.975),
                  l = quantile(value, 0.025)) |> 
        mutate(class = "fit") |> 
        rbind(realdf) |>
        ggplot(aes(x = name, y = m, color = class, shape = class)) +
        geom_linerange(aes(ymin = l, ymax = u)) +
        geom_point() +
        scale_color_manual(values = c("#36cfe3", "black"), name = "") +
        scale_shape_manual(values = c(19, 3), name = "") +
        labs(x = "# infected bins", y = "frequency") +
        ggtitle("c4") +
        theme_classic() +
        theme(text = element_text(size = 12, family = "mono"),
              axis.title.y = element_text(angle = 0, vjust = 0.5))
    
    fit$draws(format = "df") |> pull(scR0) |> mean()
}

mod <- cmdstan_model("mvnstar.stan")

ret <- data.frame(sapply("c4", fit_R0))
write_csv(ret, "ret.csv")
ret <- read_csv("ret.csv")

# sims
R0range <- 1:1
naid <- 100
keep <- which(dat$p > 0.005)
nbin <- length(keep)
R0s <- 1
sim <- data.frame(matrix(0, nrow = nbin * length(R0range) * 2, ncol = naid + 3))
colnames(sim)[(naid + 1):(naid + 3)] <- c("state", "bin", "R0")
for (i in 1:length(R0range)) {
    idxx <- (1:nbin) + nbin * 2 * (i - 1)
    idxy <- ((nbin + 1):(2 * nbin)) + nbin * 2 * (i - 1)

    sim[idxx, 1:naid] <- t(read_csv(paste0("simulations/x_", R0range[i], ".csv"), col_names = F))
    sim$state[idxx] <- "x"
    sim$bin[idxx] <- keep
    sim$R0[idxx] <- R0s[R0range[i]]
    
    sim[idxy, 1:naid] <- t(read_csv(paste0("simulations/y_", R0range[i], ".csv"), col_names = F))
    sim$state[idxy] <- "y"
    sim$bin[idxy] <- keep
    sim$R0[idxy]<- R0s[R0range[i]]
}
simproc <- sim |> select(-R0) |> 
    pivot_longer(starts_with("X"), names_to = "aid", values_to = "val") |> 
    mutate(aid = as.numeric(str_remove(aid, "X"))) |>
    # pivot_wider(id_cols = c(aid, bin, R0), names_from = "state", values_from = "val") |> 
    pivot_wider(id_cols = c(aid, bin), names_from = "state", values_from = "val") |> 
    mutate(z = y / x)
datdata <- data.frame(aid = dat$aid, bin = dat$bin, mu = dat$mu, phi = dat$phi) |> 
    distinct() |> filter(bin %in% keep)
simproc <- simproc |> left_join(datdata, by = c("bin", "aid")) |> na.omit()
simproc$z[simproc$z < 1e-16] = 1e-16
simproc$z[simproc$z > 1-1e-16] = 1-1e-16
simproc$R <- with(simproc, rnbinom(nrow(simproc), mu = mu * z, size = phi))
simstat <- simproc |> 
    # group_by(aid, R0) |> 
    group_by(aid) |> 
    summarize(stat = sum(R > 0)) |> ungroup()
# simstat |> filter(R0 == 1) |> pull(stat) |> table() |> plot()
simstat |> pull(stat) |> table() |> plot()
