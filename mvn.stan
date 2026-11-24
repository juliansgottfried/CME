data {
    int n; // number of observations
    int nbin; // number of host bins
    int naid; // number of human subjects
    
    array[n] int bin; // bin id
    array[n] int aid; // aid aid

    vector[nbin] p; // bin source abundance
    real P; // total source abundance

    vector[n] x; // abundance per bin per aid
    vector[naid] X; // total abundance per aid
    vector[naid] Y; // total infected abundance per aid

    vector[n] mu; // read average per bin per aid
    vector[n] phi; // read overdispersion per bin per aid
    array[n] int R; // read counts per bin per aid
    
    array[naid] int stat; // number of non-zero bins per aid
    array[naid] int s; // number of bins per aid
    
    real kappa; // compound community parameter
    real f; // fraction of source community infected
}

parameters {
    // real<lower=0,upper=1> f; // fraction of source community infected
    // real<lower=0> R0; // NT * R0, scaled basic reproductive number
    // real<lower=0> eta; // NT * eta, scaled chance of segregation error
}

model {
    // priors
    // f ~ beta(1, 99);
    // R0 ~ gamma(1, 1);
    // eta ~ gamma(100, 1);
    // 
    // real a = 100;
    // 
    // real EY; // expected infected abundance per aid
    // vector[nbin] Ey; // expected infected abundance per bin per aid
    // 
    // matrix[3, 3] AXY = rep_matrix(0, 3, 3); // matrix for solving for XY covariance
    // matrix[6, 6] AxYy = rep_matrix(0, 6, 6); // matrix for solving for xYy covariance
    // 
    // vector[2] sXY; // diffusion values for X and Y
    // matrix[nbin, 3] sxYy; // diffusion values for x, Y, and y
    // 
    // vector[3] bXY; // vector for solving for XY covariance
    // matrix[nbin, 6] bxYy; // vector for solving for xYy covariance
    // 
    // vector[3] SigmaXY; // XY covariance values
    // array[nbin] vector[6] SigmaxYy; // xYy covariance values
    // 
    // matrix[2, 2] Sigma22; // xY covariance matrix, for inversion
    // row_vector[2] Sigma12; // xYy covariance, i.e. without variances
    // 
    // vector[naid] Ycondmu; // conditional infected abundance per aid, given X
    // real Ycondvar; // conditional variance of infected abundance per aid, given X
    // 
    // real ycondmu; // conditional infected abundance per bin per aid, given x, Y
    // real zmu; // expected infected fraction per bin per aid
    // vector[n] z; // vector of all fractions
    // 
    // EY = (-(kappa + eta - R0 * P) + sqrt((kappa + eta - R0 * P) ^ 2 - 4 * R0 * (-kappa * P * f))) / (2 * R0);
    // Ey = p * (kappa * f + R0 * EY) / (kappa + eta + R0 * EY);
    // 
    // AXY[1, 1] = -2 * kappa;
    // AXY[2, 1] = R0 * EY;
    // AXY[2, 2] = -2 * kappa - eta + R0 * (P - 2 * EY);
    // AXY[3, 2] = 2 * R0 * EY;
    // AXY[3, 3] = -2 * kappa - 2 * eta + 2 * R0 * (P - 2 * EY);
    // 
    // AxYy[1, 1] = -2 * kappa;
    // AxYy[2, 2] = -2 * kappa - eta + R0 * (P - 2 * EY);
    // AxYy[3, 1] = R0 * EY;
    // AxYy[3, 3] = -2 * kappa - eta - R0 * EY;
    // AxYy[4, 4] = -2 * kappa - 2 * eta + 2 * R0 * (P - 2 * EY);
    // AxYy[5, 2] = R0 * EY;
    // AxYy[5, 5] = -kappa - eta - R0 * EY;
    // AxYy[6, 3] = 2 * R0 * EY;
    // AxYy[6, 6] = -2 * kappa - 2 * eta - 2 * R0 * EY;
    // 
    // sXY[1] = sqrt(P * (1 - P));
    // sXY[2] = sqrt(EY * (1 - EY));
    // 
    // sxYy[, 1] = rep_vector(sqrt(EY * (1 - EY)), nbin);
    // sxYy[, 2] = sqrt(p .* (1 - p));
    // sxYy[, 3] = sqrt(Ey .* (1 - Ey));
    // 
    // bXY[1] = -sXY[1] ^ 2;
    // bXY[2] = -sXY[1] * sXY[2];
    // bXY[3] = -sXY[2] ^ 2;
    // 
    // bxYy[, 1] = -sxYy[, 1] ^ 2;
    // bxYy[, 2] = -sxYy[, 1] .* sxYy[, 2];
    // bxYy[, 3] = -sxYy[, 1] .* sxYy[, 3];
    // bxYy[, 4] = -sxYy[, 2] ^ 2;
    // bxYy[, 5] = -sxYy[, 2] .* sxYy[, 3];
    // bxYy[, 6] = -sxYy[, 3] ^ 2;
    // 
    // SigmaXY = mdivide_left_tri_low(AXY, bXY);
    // 
    // for (i in 1:nbin) {
    //     AxYy[3, 2] = R0 * (p[i] - Ey[i]);
    //     AxYy[5, 4] = AxYy[3, 2];
    //     AxYy[6, 5] = 2 * AxYy[3, 2];
    //     SigmaxYy[i] = mdivide_left_tri_low(AxYy, bxYy[i]');
    // }
    // 
    // Ycondmu = EY + (SigmaXY[2] / SigmaXY[1]) .* (X - P);
    // Ycondvar = SigmaXY[3] - (SigmaXY[2] / SigmaXY[1]) * SigmaXY[2];
    // 
    // Sigma22[2, 2] = SigmaXY[3];
    // for (i in 1:n) {
    //     Sigma22[1, 1] = SigmaxYy[bin[i]][1];
    //     Sigma22[2, 1] = SigmaxYy[bin[i]][2];
    //     Sigma22[1, 2] = SigmaxYy[bin[i]][2];
    //     Sigma12 = [SigmaxYy[bin[i]][3], SigmaxYy[bin[i]][5]];
    // 
    //     ycondmu = Ey[bin[i]] + (Sigma12 / Sigma22) * [x[i] - p[bin[i]], Y[aid[i]] - EY]';
    //     z[i] = ycondmu / x[i];
    //     z[i] = -log_sum_exp(-log_sum_exp(a * z[i], a * 1e-7), -a) / a;
    //     if (z[i] < 1e-7) { z[i] = 1e-7;
    //     } else if (z[i] > 1-1e-7) z[i] = 1-1e-7;
    // }
    // 
    // Y ~ normal(Ycondmu, Ycondvar);
    // // R ~ neg_binomial_2(mu .* z, phi);
    // 
    // // vector[n] ltheta;
    // // ltheta = phi .* (log(phi) - log(mu .* z + phi));
    // vector[n] thetac;
    // thetac = 1 - exp(phi .* (log(phi) - log(mu .* z + phi)));
    // for (i in 1:n) {
    //     if (thetac[i] < 1e-7) { thetac[i] = 1e-7;
    //     } else if (thetac[i] > 1-1e-7) thetac[i] = 1-1e-7;
    // }
    // 
    // int pos;
    // pos = 1;
    // for (i in 1:naid) {
    //     stat[i] ~ poisson_binomial(segment(thetac, pos, s[i]));
    //     // if (stat[i] == 0) {
    //     //     target += sum(segment(ltheta, pos, s[i]));
    //     // } else {
    //     //     target += log1m_exp(sum(segment(ltheta, pos, s[i])));
    //     // }
    //     pos = pos + s[i];
    // }
}

generated quantities {
    // real f = beta_rng(1, 99);
    real R0 = gamma_rng(1, 1);
    real eta = gamma_rng(100, 1);

    real a = 100;

    real EY; // expected infected abundance per aid
    vector[nbin] Ey; // expected infected abundance per bin per aid

    matrix[3, 3] AXY = rep_matrix(0, 3, 3); // matrix for solving for XY covariance
    matrix[6, 6] AxYy = rep_matrix(0, 6, 6); // matrix for solving for xYy covariance

    vector[2] sXY; // diffusion values for X and Y
    matrix[nbin, 3] sxYy; // diffusion values for x, Y, and y

    vector[3] bXY; // vector for solving for XY covariance
    matrix[nbin, 6] bxYy; // vector for solving for xYy covariance

    vector[3] SigmaXY; // XY covariance values
    array[nbin] vector[6] SigmaxYy; // xYy covariance values

    matrix[2, 2] Sigma22; // xY covariance matrix, for inversion
    row_vector[2] Sigma12; // xYy covariance, i.e. without variances

    vector[naid] Ycondmu; // conditional infected abundance per aid, given X
    real Ycondvar; // conditional variance of infected abundance per aid, given X

    real ycondmu; // conditional infected abundance per bin per aid, given x, Y
    real zmu; // expected infected fraction per bin per aid
    vector[n] z; // vector of all fractions

    EY = (-(kappa + eta - R0 * P) + sqrt((kappa + eta - R0 * P) ^ 2 - 4 * R0 * (-kappa * P * f))) / (2 * R0);
    Ey = p * (kappa * f + R0 * EY) / (kappa + eta + R0 * EY);

    AXY[1, 1] = -2 * kappa;
    AXY[2, 1] = R0 * EY;
    AXY[2, 2] = -2 * kappa - eta + R0 * (P - 2 * EY);
    AXY[3, 2] = 2 * R0 * EY;
    AXY[3, 3] = -2 * kappa - 2 * eta + 2 * R0 * (P - 2 * EY);

    AxYy[1, 1] = -2 * kappa;
    AxYy[2, 2] = -2 * kappa - eta + R0 * (P - 2 * EY);
    AxYy[3, 1] = R0 * EY;
    AxYy[3, 3] = -2 * kappa - eta - R0 * EY;
    AxYy[4, 4] = -2 * kappa - 2 * eta + 2 * R0 * (P - 2 * EY);
    AxYy[5, 2] = R0 * EY;
    AxYy[5, 5] = -kappa - eta - R0 * EY;
    AxYy[6, 3] = 2 * R0 * EY;
    AxYy[6, 6] = -2 * kappa - 2 * eta - 2 * R0 * EY;

    sXY[1] = sqrt(P * (1 - P));
    sXY[2] = sqrt(EY * (1 - EY));

    sxYy[, 1] = rep_vector(sqrt(EY * (1 - EY)), nbin);
    sxYy[, 2] = sqrt(p .* (1 - p));
    sxYy[, 3] = sqrt(Ey .* (1 - Ey));

    bXY[1] = -sXY[1] ^ 2;
    bXY[2] = -sXY[1] * sXY[2];
    bXY[3] = -sXY[2] ^ 2;

    bxYy[, 1] = -sxYy[, 1] ^ 2;
    bxYy[, 2] = -sxYy[, 1] .* sxYy[, 2];
    bxYy[, 3] = -sxYy[, 1] .* sxYy[, 3];
    bxYy[, 4] = -sxYy[, 2] ^ 2;
    bxYy[, 5] = -sxYy[, 2] .* sxYy[, 3];
    bxYy[, 6] = -sxYy[, 3] ^ 2;

    SigmaXY = mdivide_left_tri_low(AXY, bXY);

    for (i in 1:nbin) {
        AxYy[3, 2] = R0 * (p[i] - Ey[i]);
        AxYy[5, 4] = AxYy[3, 2];
        AxYy[6, 5] = 2 * AxYy[3, 2];
        SigmaxYy[i] = mdivide_left_tri_low(AxYy, bxYy[i]');
    }

    Ycondmu = EY + (SigmaXY[2] / SigmaXY[1]) .* (X - P);
    Ycondvar = SigmaXY[3] - (SigmaXY[2] / SigmaXY[1]) * SigmaXY[2];

    Sigma22[2, 2] = SigmaXY[3];
    for (i in 1:n) {
        Sigma22[1, 1] = SigmaxYy[bin[i]][1];
        Sigma22[2, 1] = SigmaxYy[bin[i]][2];
        Sigma22[1, 2] = SigmaxYy[bin[i]][2];
        Sigma12 = [SigmaxYy[bin[i]][3], SigmaxYy[bin[i]][5]];

        ycondmu = Ey[bin[i]] + (Sigma12 / Sigma22) * [x[i] - p[bin[i]], Y[aid[i]] - EY]';
        z[i] = ycondmu / x[i];
        z[i] = -log_sum_exp(-log_sum_exp(a * z[i], a * 1e-16), -a) / a;
        if (z[i] < 1e-16) { z[i] = 1e-16;
        } else if (z[i] > 1-1e-16) z[i] = 1-1e-16;
    }

    // array[n] int R_sim;
    // R_sim = neg_binomial_2_rng(mu .* z, phi);

    vector[n] thetac;
    thetac = 1 - exp(phi .* (log(phi) - log(mu .* z + phi)));
    for (i in 1:n) {
        if (thetac[i] < 1e-16) { thetac[i] = 1e-16;
        } else if (thetac[i] > 1-1e-16) thetac[i] = 1-1e-16;
    }

    array[naid] int stat_sim;
    int pos;
    pos = 1;
    for (i in 1:naid) {
        stat_sim[i] = poisson_binomial_rng(segment(thetac, pos, s[i]));
        pos = pos + s[i];
    }
}
