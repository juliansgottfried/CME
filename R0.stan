data {
    int n;
    vector[n] M;
    vector[n] xcov;
    vector[n] xvar;
    vector[n] X;
    array[n] int R;
    real thresh;
}

parameters {
    real<lower=1.427928> R0;
}

model {
    real zstar;
    vector[n] smax;
    
    R0 ~ exponential(0.5);
    
    for (i in 1:n) {
        zstar = 1 - 1 / (X[i] * R0);
        // smax[i] = log_sum_exp(100 * zstar, 100 * 0.001) / 100;
        if (zstar > thresh) {
            R[i] ~ neg_binomial_2(M[i] * xcov[i] * zstar, M[i] * xcov[i] ^ 2 / (xvar[i] - xcov[i]));
        } else {
            R[i] ~ neg_binomial_2(M[i] * xcov[i] * thresh, M[i] * xcov[i] ^ 2 / (xvar[i] - xcov[i]));
        }
    }
}
