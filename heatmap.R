kappa <- 117.4605
n <- 300
P <- seq(from = 0.01, to = 1, length.out = n)
R0 <- seq(from = 0.01, to = 20, length.out = n)
record <- data.frame(P = rep(P, each = n), 
                     R0 = rep(R0, n),
                     mean = 0, p0 = 0)
idx <- 1
for (i in 1:n) {
    print(i)
    sims <- rbeta(10000, kappa * P[i], kappa * (1 - P[i]))
    for (j in 1:n) {
        zstar <- pmax(1 - 1 / (sims * R0[j]), 0)
        record$p0[idx] <- mean(zstar == 0)
        if (record$p0[idx] == 1) { record$mean[idx] <- 0
        } else nonzeromean <- record$mean[idx] <- mean(zstar[zstar > 0])
        idx <- idx + 1
    }
}
record |>
    ggplot(aes(x = P, y = R0)) +
    geom_tile(aes(fill = mean)) +
    geom_hline(yintercept = 1, color = "white", linetype = 3) +
    scale_fill_viridis_c(option = 16, name = "expected\nprevalence") +
    labs(y = bquote(R[0]/eta), x = bquote(sum(p[i]))) +
    ggtitle("parameter space") +
    theme_classic() +
    theme(text = element_text(size = 12, family = "serif"),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1))
 

Nt <- 10^13
kappa <- 117.4605
m <- kappa / Nt
eta <- m
R0 <- 5 * m
pi <- 0.1
p <- 0.01
pT <- 0.57
xbet <- rbeta(10000, kappa * p, kappa * (1 - p))
xnor <- truncnorm::rtruncnorm(10000, 
                              mean = p, 
                              sd = sqrt(p * (1 - p) / (kappa + 1)),
                              a = 0, b = 1)
data.frame(beta = xbet, normal = xnor) |> 
    arrange(normal) |> mutate(idx = row_number()) |> 
    pivot_longer(-idx) |> 
    ggplot(aes(x = value, fill = name)) +
    geom_histogram(color = "white") +
    theme_classic()


a <- R0
b <- m + eta - R0 * pT
c <- -m * pT * pi
EY <- (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a)

Ey <- p * (m * pi + R0 * EY) / (m + eta + R0 * EY)
Vary <- Ey * (1 - Ey)
ynor <- truncnorm::rtruncnorm(10000, 
                      mean = Ey, 
                      sd = sqrt(Vary),
                      a = 0, b = 1)
data.frame(x = xnor, y = ynor) |> 
    # pivot_longer(everything()) |> 
    ggplot(aes(x = x, y = y)) +
    geom_point() +
    theme_classic()
