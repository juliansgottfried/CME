data {
    int n;
    int naid;
    array[naid] int s;
    vector[n] M;
    vector[n] xcov;
    vector[n] xvar;
    array[n] int R;
}

parameters {
    vector<lower=0,upper=1>[naid] Z;
}

model {
    Z ~ beta(1, 1);
    int pos;
    pos = 1;
    for (i in 1:naid) {
        segment(R, pos, s[i]) ~ neg_binomial_2(segment(M, pos, s[i]) .* segment(xcov, pos, s[i]) * Z[i], segment(M, pos, s[i]) .* segment(xcov, pos, s[i]) .^ 2 ./ (segment(xvar, pos, s[i]) - segment(xcov, pos, s[i])));
        pos = pos + s[i];
  }
}
