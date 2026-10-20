library(tidyverse)

setwd("~/Desktop/CME")

outx <- read_csv("outx_100.csv", col_names = F)
outy <- read_csv("outy_100.csv", col_names = F)

idx <- 1
cbind(outx[, idx], outy[, idx]) |> 
    `colnames<-`(c("x","y")) |> 
    mutate(z = y / x) |> 
    sample_n(10000) |> 
    ggplot(aes(x = x, y = log(z))) +
    geom_point(shape = 3, alpha = 0.5) +
    labs(y = "log[y/x]") +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))

sampled <- sample(1:nrow(outx), 1000, replace = F)
flattenx <- unlist(outx[sampled, ])
flatteny <- unlist(outy[sampled, ])

cbind(flattenx, flatteny) |> data.frame() |> 
    mutate(sp = rep(1:100, each = length(sampled))) |> 
    `colnames<-`(c("x", "y", "sp")) |> 
    mutate(z = y / x) |> 
    mutate(class = floor((sp - 1) / 10)) |> 
    ggplot(aes(x = x, y = log(z))) +
    geom_point(shape = 12, alpha = 0.1) +
    labs(x = "N/K", y = "log[I/N]") +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))



