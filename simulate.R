runner <- function(scR0, P) {
    eta <- 10 ^ 5
    R0 <- scR0 * eta
    
    AXY <- matrix(0, nrow = 3, ncol = 3)
    AxYy <- matrix(0, nrow = 6, ncol = 6)
    # Axy <- matrix(0, nrow = 3, ncol = 3)
    
    sxYy <- matrix(0, nrow = dat$nbin, ncol = 3)
    # sxy <- matrix(0, nrow = dat$nbin, ncol = 2)
    
    bxYy <- matrix(0, nrow = dat$nbin, ncol = 6)
    # bxy <- matrix(0, nrow = dat$nbin, ncol = 3)
    
    SigmaxYy <- matrix(0, dat$nbin, 6)
    Sigmaxy <- matrix(0, dat$nbin, 3)
    
    Sigma22 <- matrix(0, 2, 2)
    
    z <- rep(0, dat$n)

    EY <- (-(dat$kappa + eta - R0 * P) + sqrt((dat$kappa + eta - R0 * P) ^ 2 - 4 * R0 * (-dat$kappa * P * dat$f))) / (2 * R0)
    Ey <- dat$p * (dat$kappa * dat$f + R0 * EY) / (dat$kappa + eta + R0 * EY)
    
    AXY[1, 1] = -2 * dat$kappa
    AXY[2, 1] = R0 * EY
    AXY[2, 2] = -2 * dat$kappa - eta + R0 * (P - 2 * EY)
    AXY[3, 2] = 2 * R0 * EY
    AXY[3, 3] = -2 * dat$kappa - 2 * eta + 2 * R0 * (P - 2 * EY)
    
    AxYy[1, 1] = -2 * dat$kappa;
    AxYy[2, 2] = -2 * dat$kappa - eta + R0 * (P - 2 * EY)
    AxYy[3, 1] = R0 * EY
    AxYy[3, 3] = -2 * dat$kappa - eta - R0 * EY
    AxYy[4, 4] = -2 * dat$kappa - 2 * eta + 2 * R0 * (P - 2 * EY)
    AxYy[5, 2] = R0 * EY
    AxYy[5, 5] = -dat$kappa - eta - R0 * EY
    AxYy[6, 3] = 2 * R0 * EY
    AxYy[6, 6] = -2 * dat$kappa - 2 * eta - 2 * R0 * EY
    
    # Axy[1, 1] = -2 * dat$kappa;
    # Axy[2, 1] = R0 * EY
    # Axy[2, 2] = -2 * dat$kappa - eta - R0 * EY
    # Axy[3, 2] = 2 * R0 * EY
    # Axy[3, 3] = -2 * dat$kappa - 2 * eta - 2 * R0 * EY
    
    sXY <- c(sqrt(P * (1 - P)), 
             sqrt(EY * (1 - EY)))
    
    sxYy[, 1] = sqrt(dat$p * (1 - dat$p))
    sxYy[, 2] = rep(sqrt(EY * (1 - EY)), dat$nbin)
    sxYy[, 3] = sqrt(Ey * (1 - Ey))
    
    # sxy[, 1] = sqrt(dat$p * (1 - dat$p))
    # sxy[, 2] = sqrt(Ey * (1 - Ey))
    
    bXY <- -c(sXY[1] ^ 2,
             sXY[1] * sXY[2],
             sXY[2] ^ 2)
    
    bxYy[, 1] = -sxYy[, 1] ^ 2
    bxYy[, 2] = -sxYy[, 1] * sxYy[, 2]
    bxYy[, 3] = -sxYy[, 1] * sxYy[, 3]
    bxYy[, 4] = -sxYy[, 2] ^ 2
    bxYy[, 5] = -sxYy[, 2] * sxYy[, 3]
    bxYy[, 6] = -sxYy[, 3] ^ 2
    
    # bxy[, 1] = -sxy[, 1] ^ 2
    # bxy[, 2] = -sxy[, 1] * sxy[, 2]
    # bxy[, 3] = -sxy[, 2] ^ 2

    SigmaXY <- solve(AXY) %*% bXY
    
    for (i in 1:dat$nbin) {
        AxYy[3, 2] = R0 * (dat$p[i] - Ey[i])
        AxYy[5, 4] = AxYy[3, 2]
        AxYy[6, 5] = 2 * AxYy[3, 2]
        SigmaxYy[i, ] = solve(AxYy) %*% bxYy[i, ]
        # Sigmaxy[i, ] = solve(Axy) %*% bxy[i, ]
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
        Sigma12 = c(SigmaxYy[dat$bin[i], 3], SigmaxYy[dat$bin[i], 5])
        # ycondmu = with(dat, Ey[bin[i]] + (Sigma12 %*% solve(Sigma22)) %*% c(x[i] - p[bin[i]], 0))
        ycondmu = Ey[dat$bin[i]] + (Sigma12 %*% solve(Sigma22))[1] * (dat$x[i] - dat$p[dat$bin[i]])
        # ycondmu = Ey[dat$bin[i]] + (Sigmaxy[dat$bin[i], 2] / Sigmaxy[dat$bin[i], 1]) * (dat$x[i] - dat$p[dat$bin[i]])
        z[i] = ycondmu / dat$x[i]
        # z[i] = -matrixStats::logSumExp(c(-matrixStats::logSumExp(c(a * z[i], a * 1e-16)), -a)) / a
        if (z[i] < 1e-16) { z[i] = 1e-16
        } else if (z[i] > 1-1e-16) z[i] = 1-1e-16
    }

    # thetac <- rep(0, dat$n)
    # thetac <- 1 - exp(dat$phi * (log(dat$phi) - log(dat$mu * z + dat$phi)))

    thetal <- dat$phi * (log(dat$phi) - log(dat$mu * z + dat$phi))
    
    ret <- rep(0, dat$naid)
    # ret2 <- rep(0, dat$naid)
    # ret3 <- rep(0, dat$naid)
    # pos = 1
    
    for (i in 1:dat$naid) {
        # hold <- poisbinom::rpoisbinom(10000, thetac[pos:(pos + dat$s[i] - 1)])
        # ret[i] = mean(hold == 0)
        # ret2[i] = mean(hold == 1)
        # ret3[i] = mean(hold == 2)
        ret[i] = sum(thetal[pos:(pos + dat$s[i] - 1)])
        pos = pos + dat$s[i]
    }
    # c(mean(ret), mean(ret2), mean(ret3), mean(z))
    
    BayesRTMB::log_sum_exp(ret) - log(length(ret))
    
    # for (i in 1:dat$naid) {
    #     ret[i] = poisbinom::rpoisbinom(1, thetac[pos:(pos + dat$s[i] - 1)])
    #     pos = pos + dat$s[i]
    # }
    # plot(table(ret))
}

alldiag <- list.files("standiagnostics", full.names = T)
allfit <- list.files("stanfits", full.names = T)
gather <- NULL
for (i in 1:length(allfiles)) {
    print(i)
    if (sum(read_csv(allfiles[i])$num_divergent) == 0) {
        tmpdf <- read_csv(allfit[i])
        gather <- rbind(gather, 
                        cbind(c(mean(tmpdf$scR0), var(tmpdf$scR0), 
                                quantile(tmpdf$scR0, c(0.025, 0.5, 0.975))),
              c(mean(tmpdf$leta), var(tmpdf$scR0), 
                quantile(tmpdf$leta, c(0.025, 0.5, 0.975)))) |> 
            data.frame() |> 
            mutate(stat = c("mean", "var", "l", "m", "u"),
                   idx = i))
    }
}
colnames(gather)[1:2] <- c("scR0", "leta")
gather$cme <- str_remove(str_remove(allfit, "stanfits/"), ".csv")[gather$idx]
pimeansdf <- data.frame(pimeans)
pimeansdf$cme <- names(pimeans)

plotgather <- gather |> 
    pivot_wider(id_cols = cme, names_from = stat, values_from = scR0) |> 
    # filter(sqrt(var) / mean < 0.5) |> 
    left_join(pimeansdf, by = "cme") |> 
    select(cme, pimeans, mean, var, l, m, u) |> 
    filter(sqrt(var) / mean < 0.15)

scR0s <- exp(seq(from = log(min(plotgather$l)) - 0.1, to = -(log(min(plotgather$pimeans)) - 0.1), length.out = 40))
fs <-  exp(seq(from = log(min(plotgather$pimeans)) - 0.1, to = -(log(min(plotgather$l)) - 0.1), length.out = 40))
grid <- expand.grid(scR0s, fs)
# ret <- cbind(grid, t(apply(grid, 1, \(x) runner(x[1], x[2]))))
ret <- cbind(grid, apply(grid, 1, \(x) runner(x[1], x[2])))
# colnames(ret) <- c("scR0", "f", "p0", "p1", "p2", "z")
colnames(ret) <- c("scR0", "f", "theta")

ret |> 
    ggplot(aes(x = log(f), y = log(scR0), fill = BayesRTMB::log1m_exp(theta))) +
    geom_tile() +
    geom_abline(slope = -1, color = "black",
                linetype = 3, linewidth = 0.5) +
    geom_pointrange(data = plotgather,
               aes(x = log(pimeans), y = log(m), ymin = log(l), ymax = log(u),
                   color = cme != "c4"),
               fill = "white", 
               alpha = 0.6, shape = 10, stroke = 0.8, linewidth = 0.15) +
    scale_color_manual(values = c("red", "#35adbd"),
                       labels = c("c4", "other"),
                       name = "cme") +
    scale_fill_gradient(low = "white", high = "black", name = "log prob\ninfection") +
    guides(fill = guide_colorbar(order = 1), 
           color = guide_legend(order = 2)) +
    labs(x = bquote("log"~pi), y = bquote("log"~R[0]/eta)) +
    theme_classic() +
    coord_fixed(ratio = diff(range(log(fs))) / diff(range(log(scR0s)))) +
    theme(text = element_text(size = 12, family = "mono"),
           axis.title.y = element_text(angle = 0, vjust = 0.5))

# ret |> mutate(R0star = R0 / eta) |> 
#     ggplot(aes(x = R0star, y = p0, group = eta, color = eta)) +
#     geom_line() +
#     theme_classic()
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
