library(tidyverse)
library(patchwork)

setwd("~/Desktop/CME")

df <- read_csv("outz_100.csv", col_names = F)

df[1000:nrow(df), ] |> 
    sample_n(1000) |> 
    pivot_longer(everything(), values_to = "z", names_to = "sp") |> 
    mutate(sp = as.numeric(str_remove(sp, "X"))) |> 
    ggplot(aes(x = z, color = sp, group = sp)) +
    geom_density() +
    scale_color_viridis_c(option = 6, alpha = 0.5) +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))

df[1000:nrow(df), ] |> 
    sample_n(10000) |> 
    pivot_longer(everything(), values_to = "z", names_to = "sp") |> 
    mutate(sp = as.numeric(str_remove(sp, "X"))) |> 
    group_by(sp) |> summarize(p0 = mean(z < 0.001),
                              mean = mean(z),
                              med = median(z),
                              max = max(z),
                              sd = sd(z)) |> 
    ggplot(aes(x = sp, y = p0)) +
    geom_point() +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))

cormat <- cor(df)
cormat |> data.frame() |>
    mutate(a = 1:100) |>
    pivot_longer(-a, names_to = "b", values_to = "cor") |>
    mutate(b = as.numeric(str_remove(b, "X"))) |>
    ggplot(aes(x = a, y = b, fill = cor)) +
    geom_tile() +
    scale_fill_viridis_c(option = 12) +
    theme_classic()

plotter <- function(i, j) {
    round(log(df[1000:nrow(df), c(i, j)]) / 10, 2) |>
        `colnames<-`(c("X1", "X2")) |>
        group_by(X1, X2) |>
        summarize(n = n()) |> ungroup() |>
        mutate(f = n / sum(n)) |>
        ggplot(aes(x = X1, y = X2, fill = log(f))) +
        geom_tile() +
        scale_fill_viridis_c(option = 12) +
        xlim(-1.5, 0) + ylim(-1.5, 0) +
        labs(x = paste0("X", i), y = paste0("X", j)) +
        theme_classic() +
        theme(text = element_text(size = 12, family = "mono"))
}
p1 <- plotter(1, 2)
p2 <- plotter(1, 100)
p3 <- plotter(99, 100)
p1 / p2 / p3

df[1000:nrow(df), ] |> 
    sample_n(1000) |> 
    pivot_longer(everything(), values_to = "z") |> 
    ggplot(aes(x = z)) + geom_density() +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))
