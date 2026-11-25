runner <- function(R0, eta) {
    # R0 <- 1
    # eta <- 1000
    
    # a <- 1000

    AXY <- matrix(0, nrow = 3, ncol = 3)
    AxYy <- matrix(0, nrow = 6, ncol = 6)
    
    sxYy <- matrix(0, nrow = dat$nbin, ncol = 3)
    
    bxYy <- matrix(0, nrow = dat$nbin, ncol = 6)
    
    SigmaxYy <- matrix(0, dat$nbin, 6)
    
    Sigma22 <- matrix(0, 2, 2)
    
    z <- rep(0, dat$n)

    EY <- with(dat, (-(kappa + eta - R0 * P) + sqrt((kappa + eta - R0 * P) ^ 2 - 4 * R0 * (-kappa * P * f))) / (2 * R0))
    Ey <- with(dat, p * (kappa * f + R0 * EY) / (kappa + eta + R0 * EY))
    
    AXY[1, 1] = -2 * dat$kappa
    AXY[2, 1] = R0 * EY
    AXY[2, 2] = with(dat, -2 * kappa - eta + R0 * (P - 2 * EY))
    AXY[3, 2] = 2 * R0 * EY
    AXY[3, 3] = with(dat, -2 * kappa - 2 * eta + 2 * R0 * (P - 2 * EY))
    
    AxYy[1, 1] = -2 * dat$kappa;
    AxYy[2, 2] = with(dat, -2 * kappa - eta + R0 * (P - 2 * EY))
    AxYy[3, 1] = R0 * EY
    AxYy[3, 3] = -2 * dat$kappa - eta - R0 * EY
    AxYy[4, 4] = with(dat, -2 * kappa - 2 * eta + 2 * R0 * (P - 2 * EY))
    AxYy[5, 2] = R0 * EY
    AxYy[5, 5] = -dat$kappa - eta - R0 * EY
    AxYy[6, 3] = 2 * R0 * EY
    AxYy[6, 6] = with(dat, -2 * kappa - 2 * eta - 2 * R0 * EY)
    
    sXY <- c(with(dat, sqrt(P * (1 - P))), 
             sqrt(EY * (1 - EY)))
    
    sxYy[, 1] = with(dat, sqrt(p * (1 - p)))
    sxYy[, 2] = rep(sqrt(EY * (1 - EY)), dat$nbin)
    sxYy[, 3] = sqrt(Ey * (1 - Ey))
    
    bXY <- -c(sXY[1] ^ 2,
             sXY[1] * sXY[2],
             sXY[2] ^ 2)
    
    bxYy[, 1] = -sxYy[, 1] ^ 2
    bxYy[, 2] = -sxYy[, 1] * sxYy[, 2]
    bxYy[, 3] = -sxYy[, 1] * sxYy[, 3]
    bxYy[, 4] = -sxYy[, 2] ^ 2
    bxYy[, 5] = -sxYy[, 2] * sxYy[, 3]
    bxYy[, 6] = -sxYy[, 3] ^ 2
    
    SigmaXY <- solve(AXY) %*% bXY
    
    for (i in 1:dat$nbin) {
        AxYy[3, 2] = R0 * (dat$p[i] - Ey[i])
        AxYy[5, 4] = AxYy[3, 2]
        AxYy[6, 5] = 2 * AxYy[3, 2]
        SigmaxYy[i, ] = solve(AxYy) %*% bxYy[i, ]
    }
    
    # Ycondmu <- with(dat, EY + (SigmaXY[2] / SigmaXY[1]) * (X - P))
    # for (i in 1:dat$naid) {
    #     if (Ycondmu[i] < 1e-12) Ycondmu[i] = 1e-12
    #     if (Ycondmu[i] > 1-1e-12) Ycondmu[i] = 1-1e-12
    # }
    # Ycondvar <- SigmaXY[3] - (SigmaXY[2] / SigmaXY[1]) * SigmaXY[2]
    
    Sigma22[2, 2] <- SigmaXY[3]
    for (i in 1:dat$n) {
        Sigma22[1, 1] = SigmaxYy[dat$bin[i], 1]
        Sigma22[2, 1] = SigmaxYy[dat$bin[i], 2]
        Sigma22[1, 2] = SigmaxYy[dat$bin[i], 2]
        Sigma12 = with(dat, c(SigmaxYy[bin[i], 3], SigmaxYy[bin[i], 5]))
        # ycondmu = with(dat, Ey[bin[i]] + (Sigma12 %*% solve(Sigma22)) %*% c(x[i] - p[bin[i]], 0))
        ycondmu = with(dat, Ey[bin[i]] + (Sigma12 %*% solve(Sigma22))[1] * (x[i] - p[bin[i]]))
        z[i] = ycondmu / dat$x[i]
        # z[i] = -matrixStats::logSumExp(c(-matrixStats::logSumExp(c(a * z[i], a * 1e-16)), -a)) / a
        if (z[i] < 1e-16) { z[i] = 1e-16
        } else if (z[i] > 1-1e-16) z[i] = 1-1e-16
    }

    thetac <- rep(0, dat$n)
    thetac <- with(dat, 1 - exp(phi * (log(phi) - log(mu * z + phi))))

    ret <- rep(0, dat$naid)
    ret2 <- rep(0, dat$naid)
    ret3 <- rep(0, dat$naid)
    pos = 1
    
    for (i in 1:dat$naid) {
        hold <- poisbinom::rpoisbinom(1000, thetac[pos:(pos + dat$s[i] - 1)])
        ret[i] = mean(hold == 0)
        ret2[i] = mean(hold == 1)
        ret3[i] = mean(hold == 2)
        pos = pos + dat$s[i]
    }
    c(mean(ret), mean(ret2), mean(ret3), mean(z))
    
    # for (i in 1:dat$naid) {
    #     ret[i] = poisbinom::rpoisbinom(1, thetac[pos:(pos + dat$s[i] - 1)])
    #     pos = pos + dat$s[i]
    # }
    # plot(table(ret))
}

scR0 = 10000 / 75000
leta = log10(75000)
runner(scR0 * (10 ^ leta), 10 ^ leta)

target0 <- round(0.45238095, 1)
target1 <- round(0.19047619, 1)
# R0s <- seq(from = 100, to = 100000, length.out = 40)
# etas <- seq(from = 100, to = 100000, length.out = 40)
scR0s <- seq(from = 0.0001, to = 1/dat$P, length.out = 100)
letas <- 5
grid <- expand.grid(scR0s, letas)
ret <- cbind(grid, t(apply(grid, 1, \(x) runner(x[1] * (10 ^ x[2]), 10 ^ x[2]))))
colnames(ret) <- c("R0", "eta", "p0", "p1", "p2", "z")
# ret$viz <- round(ret$p0, 1) >= target0 & round(ret$p1, 1) == target1
ret |> 
    ggplot(aes(x = eta, y = R0, fill = 100 * p0)) +
    geom_tile() +
    scale_fill_viridis_c(option = 15, name = "%0") +
    scale_alpha_manual(values = c(0.75, 1), name = "match") +
    geom_abline(slope = 1/dat$P, color = "white", 
                linetype = 3, linewidth = 1) +
    guides(alpha = "none") +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))
ret |> mutate(R0star = R0 / eta) |> 
    ggplot(aes(x = R0star, y = p0, group = eta, color = eta)) +
    geom_line() +
    theme_classic()
# ret |> mutate(R0star = R0 / eta) |> 
#     ggplot(aes(x = log(R0star), y = p1, group = eta, color = eta)) +
#     geom_line() +
#     theme_classic()
# ret |> mutate(R0star = R0 / eta) |> 
#     ggplot(aes(x = log(R0star), y = p2, group = eta, color = eta)) +
#     geom_line() +
#     theme_classic()
# ret |> filter(!is.na(var) & !is.infinite(var)) |> 
#     group_by(R0, f) |> summarize(var = mean(var)) |> 
#     ggplot(aes(x = R0, y = f, fill = var)) +
#     geom_tile() +
#     theme_classic()
# ret |> filter(!is.na(var) & !is.infinite(var)) |> 
#     group_by(eta, f) |> summarize(var = mean(var)) |> 
#     ggplot(aes(x = eta, y = f, fill = var)) +
#     geom_tile() +
#     theme_classic()
ret |> 
    filter(eta == etas[1]) |> 
    ggplot(aes(x = R0, y = plot)) +
    geom_line() +
    theme_classic()
ret |> 
    filter(eta == etas[4]) |> 
    ggplot(aes(x = R0, y = plot)) +
    geom_line() +
    theme_classic()
ret |> 
    filter(eta == etas[10]) |> 
    ggplot(aes(x = R0, y = plot)) +
    geom_line() +
    theme_classic()
ret |> 
    filter(eta == etas[25]) |> 
    ggplot(aes(x = R0, y = plot)) +
    geom_line() +
    theme_classic()

juliasim <- data.frame(matrix(0, nrow = 30, ncol = 2))
for (i in 1:30) {
    tmpz <- read_csv(paste0("simulations_etaverysmall/z_", i, ".csv"), col_names = F)
    tmpz <- c(as.matrix(tmpz))
    # tmpz[tmpz < 0] <- 0
    # tmpz[tmpz > 1] <- 1
    tmpz <- tmpz[!is.na(tmpz)]
    juliasim[i, ] <- c(R0s[i], mean(tmpz))
}
colnames(juliasim) <- c("R0", "plot")
juliasim |> 
    ggplot(aes(x = R0, y = plot)) +
    geom_line() +
    theme_classic()
