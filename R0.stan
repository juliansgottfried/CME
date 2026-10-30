data {
    int n;
    vector[n] M;
    vector[n] xcov;
    vector[n] xvar;
    vector[n] X;
    array[n] int R;
}

parameters {
    real<lower=0> R0;
}

model {
    real zstar;
    vector[n] smax;
    
    R0 ~ exponential(1);
    
    for (i in 1:n) {
        zstar = 1 - 1 / (X[i] * R0);
        smax[i] = log_sum_exp(100 * zstar, 100 * 0.001) / 100;
    }
    
    R ~ neg_binomial_2(M .* xcov .* smax, M .* xcov ^ 2 ./ (xvar - xcov));
}
