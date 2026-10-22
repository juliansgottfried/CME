library(tidyverse)

setwd("~/Desktop/CME")

outx <- data.table::fread("outx_100.csv", header = F)
outy <- data.table::fread("outy_100.csv", header = F)

outx_df <- data.frame(outx)
outy_df <- data.frame(outy)

outx |> sample_n(100000) |>
    select(V1) |> 
    mutate(t = row_number()) |> 
    pivot_longer(-t, names_to = "sp", values_to = "x") |> 
    mutate(sp = as.numeric(str_remove(sp, "X"))) |> 
    ggplot(aes(x = t, y = x, group = sp)) +
    geom_line(alpha = 0.5) +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))
outy |> sample_n(1000) |>
    select(X100) |> 
    mutate(t = row_number()) |> 
    pivot_longer(-t, names_to = "sp", values_to = "y") |> 
    mutate(sp = as.numeric(str_remove(sp, "X"))) |> 
    ggplot(aes(x = t, y = y)) +
    geom_line(alpha = 0.5) +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))

outz <- outy / outx
outz <- outz[apply(outz, 1, \(x) sum(is.na(x) | x > 10000) == 0), ]
mat <- cor(outz)
mat |> data.frame() |> 
    mutate(sp1 = 1:10) |> 
    pivot_longer(-sp1) |> 
    mutate(sp2 = as.numeric(str_remove(name, "X"))) |> 
    ggplot(aes(x = sp1, y = sp2, fill = value)) +
    geom_tile() +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))
    
idx <- 1
cbind(outx[, idx], outy[, idx]) |>
    `colnames<-`(c("x","y")) |>
    mutate(z = y / x) |>
    sample_n(10000) |>
    ggplot(aes(x = x, y = z)) +
    geom_point(shape = 3, alpha = 0.5) +
    labs(y = "log[y/x]") +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))

sampled <- sample(1:nrow(outx), 1000, replace = F)
flattenx <- unlist(outx[sampled, ])
flatteny <- unlist(outy[sampled, ])

xc <- 0
cbind(flattenx, flatteny) |> data.frame() |>
    mutate(sp = rep(1:10, each = length(sampled))) |>
    `colnames<-`(c("x", "y", "sp")) |>
    mutate(thresh = ifelse(x > xc, y / (x - xc), y / x)) |>
    ggplot(aes(x = x, y = y)) +
    geom_point(shape = 12, alpha = 1) +
    labs(x = "N/K", y = "I/K") +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))





library(cmdstanr)
library(posterior)

hostx <- read.table("cipro_abundance/host_abundance_sr.txt", sep = "\t", header = T)
taxon <- read.table("cipro_abundance/hosts.txt", sep = "\t", header = T)
key <-  read.table("cipro_abundance/short_column_key.txt", sep = "\t", header = T)
canons <- key |> select(day, canon) |> distinct() |> 
    filter(day < 0) |> pull(canon) |> unique()
xdat <- hostx[, c("aid", "bin", canons)] |> 
    left_join(taxon |> select(aid, bin, species), by = c("aid", "bin")) |> 
    pivot_longer(all_of(canons), names_to = "time", values_to = "x") |> 
    mutate(x = x / 100) |> 
    group_by(aid, time, species) |> 
    summarize(x = sum(x)) |> mutate(x = x / sum(x)) |> 
    ungroup() |> group_by(species) |> 
    mutate(xbar = sum(x) / (42 * length(canons))) |> ungroup() |> 
    filter(xbar > 0) |> mutate(x = x + min(xbar) / 1000) |> 
    group_by(species) |> 
    mutate(xbar2 = sum(x) / (42 * length(canons))) |> ungroup() |> 
    mutate(sp = as.numeric(factor(species))) |> 
    select(sp, species, x, xbar, xbar2)

xdat |>
    ggplot(aes(x = xbar)) +
    geom_histogram(color = "white", fill = "black", bins = 100) +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))

dat <- list(n_obs = nrow(xdat),
            x = xdat$x,
            xbar = xdat$xbar2)

mod <- cmdstan_model("beta.stan")

fit <- mod$sample(
    data = dat,
    chains = 4,
    parallel_chains = 4,
    refresh = 500)

fit$draws(format = "df") |> select(kappa) |> 
    ggplot(aes(x = kappa)) +
    geom_density() +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))

kappa <- fit$draws(format = "df") |> pull(kappa) |> mean()

xdat |> select(species, xbar) |> 
    filter(species %in% c4sp) |> distinct() |> 
    select(xbar) |> 
    write_csv("xbar.csv", col_names = F)
