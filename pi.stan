data {
    int n;
    vector[n] M;
    vector[n] xcov;
    vector[n] xvar;
    array[n] int R;
}

parameters {
    real<lower=0,upper=1> f;
}

model {
    f ~ beta(1, 1);
    R ~ neg_binomial_2(M .* xcov * f, M .* xcov .^ 2 ./ (xvar - xcov));
}
