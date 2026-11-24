runner <- function(f, R0, eta) {
    a <- 100
    
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
    AXY[3, 3] = -with(dat, 2 * kappa - 2 * eta + 2 * R0 * (P - 2 * EY))
    
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
    
    sxYy[, 1] = rep(sqrt(EY * (1 - EY)), dat$nbin)
    sxYy[, 2] = with(dat, sqrt(p * (1 - p)))
    sxYy[, 3] = sqrt(Ey * (1 - Ey))
    
    bXY <- c(-sXY[1] ^ 2,
             -sXY[1] * sXY[2],
             -sXY[2] ^ 2)
    
    bxYy[, 1] = -sxYy[, 1] ^ 2;
    bxYy[, 2] = -sxYy[, 1] * sxYy[, 2];
    bxYy[, 3] = -sxYy[, 1] * sxYy[, 3];
    bxYy[, 4] = -sxYy[, 2] ^ 2;
    bxYy[, 5] = -sxYy[, 2] * sxYy[, 3];
    bxYy[, 6] = -sxYy[, 3] ^ 2;
    
    SigmaXY <- solve(AXY) %*% bXY
    
    for (i in 1:dat$nbin) {
        AxYy[3, 2] = R0 * (dat$p[i] - Ey[i]);
        AxYy[5, 4] = AxYy[3, 2];
        AxYy[6, 5] = 2 * AxYy[3, 2];
        SigmaxYy[i, ] <- solve(AxYy) %*% bxYy[i, ]
    }
    
    Ycondmu <- with(dat, EY + (SigmaXY[2] %*% solve(SigmaXY[1])) %*% (X - P))
    Ycondvar <- SigmaXY[3] - (SigmaXY[2] %*% solve(SigmaXY[1])) %*% SigmaXY[2]
    
    Sigma22[2, 2] <- SigmaXY[3]
    for (i in 1:dat$n) {
        Sigma22[1, 1] = SigmaxYy[dat$bin[i], 1]
        Sigma22[2, 1] = SigmaxYy[dat$bin[i], 2]
        Sigma22[1, 2] = SigmaxYy[dat$bin[i], 2]
        Sigma12 = with(dat, c(SigmaxYy[bin[i], 3], SigmaxYy[bin[i], 5]))
    
        ycondmu = with(dat, Ey[bin[i]] + (Sigma12 %*% solve(Sigma22)) %*% c(x[i] - p[bin[i]], Y[aid[i]] - EY))
        z[i] = ycondmu / dat$x[i]
        z[i] = -matrixStats::logSumExp(c(-matrixStats::logSumExp(c(a * z[i], a * 1e-16)), -a)) / a
    }
    
    thetac <- rep(0, dat$n)
    thetac <- with(dat, 1 - exp(phi * (log(phi) - log(mu * z + phi))))
    
    ret <- rep(0, dat$naid)
    pos = 1
    
    for (i in 1:dat$naid) {
        ret[i] = mean(poisbinom::rpoisbinom(1000, thetac[pos:(pos + dat$s[i] - 1)]) == 0)
        pos = pos + dat$s[i]
    }
    mean(ret)
    
    # for (i in 1:dat$naid) {
    #     ret[i] = poisbinom::rpoisbinom(1, thetac[pos:(pos + dat$s[i] - 1)])
    #     pos = pos + dat$s[i]
    # }
    # plot(table(ret))
    
    # mean(z)
}

target <- 0.45238095
f <- 0.045
R0s <- seq(from = 0.0001, to = 1000000, length.out = 30)
etas <- seq(from = 0.0001, to = 1000000, length.out = 30)
grid <- expand.grid(f, R0s, etas)
ret <- cbind(grid, (apply(grid, 1, \(x) runner(x[1], x[2], x[3]))))
colnames(ret) <- c("f", "R0", "eta", "plot")
ret |> filter(!is.na(plot) & !is.infinite(plot)) |>
    group_by(eta, R0) |> summarize(plot = mean(plot)) |>
    ggplot(aes(x = eta, y = R0, fill = plot, alpha = round(plot, 1) > round(target, 1))) +
    geom_tile() +
    scale_fill_viridis_c(option = 15) +
    scale_alpha_manual(values = c(0.1, 1), name = "match") +
    theme_classic()

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
