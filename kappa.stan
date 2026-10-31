data {
    int n_sub;
    int n_bin;
    array[n_sub] row_vector[n_bin] x;
    simplex[n_bin] p;
}

parameters {
    real<lower=0> kappa;
}

model {
    kappa ~ lognormal(4, 1);
    for (i in 1:n_sub) {
        x[i] ~ beta_proportion(p, kappa);
    }
}
