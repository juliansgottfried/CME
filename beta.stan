data {
    int n_obs;
    vector<lower=0,upper=1>[n_obs] x;
    vector<lower=0,upper=1>[n_obs] xbar;
}

parameters {
    real<lower=0> kappa;
}

model {
    kappa ~ lognormal(4, 1);
    x ~ beta(kappa * xbar, kappa * (1 - xbar));
}
