# library(tidyverse)

N <- 10 # 100 / 1000
S0 <- N
I0 <- 0

dS <- 0.05 # 10 / 50
dI <- 0.0005 # 0.005 / 0.05
beta <- 0.001 / N # 0.01 / 0.01
gamma <- 0.001 # 0.01 / 0.01

G <- 100000
dt <- 1
J <- 100

S <- rep(0, G)
I <- rep(0, G)
S[1] <- S0
I[1] <- I0
for (i in 2:G) {
    muSI <- rbinom(1, floor(S[i - 1]), 1 - exp(-beta * I[i - 1] * dt))
    S[i] <- min(max(S[i - 1] + dS * dt - muSI - rbinom(1, floor(S[i - 1]), 1 - exp(-gamma * dt)), 0), N)
    I[i] <- min(max(I[i - 1] + dI * dt + muSI - rbinom(1, floor(I[i - 1]), 1 - exp(-gamma * dt)), 0), N)
}
data.frame(t = c(1:G, 1:G), val = c(S, I), state = c(rep("S", G), rep("I", G))) %>%
    mutate(state = factor(state, levels = c("S", "I"))) %>%
    ggplot(aes(x = t, y = val)) +
    geom_line() +
    facet_wrap(~ state) +
    theme_classic()

data.frame(obs = as.numeric(rbinom(G, 100, I / (1000 + I)) > 1)) %>% 
    ggplot(aes(x = 1:G, y = obs)) +
    geom_point() + 
    theme_classic()


store[j] <- mean(I[(G - 100):G])

store <- rep(0, J)
for (j in 1:J) {
    S <- rep(0, G)
    I <- rep(0, G)
    S[1] <- S0
    I[1] <- I0
    for (i in 2:G) {
        muSI <- rbinom(1, floor(S[i - 1]), 1 - exp(-beta * I[i - 1] * dt))
        S[i] <- min(max(S[i - 1] + dS * dt - muSI - rbinom(1, floor(S[i - 1]), 1 - exp(-gamma * dt)), 0), N)
        I[i] <- min(max(I[i - 1] + dI * dt + muSI - rbinom(1, floor(I[i - 1]), 1 - exp(-gamma * dt)), 0), N)
    }
    store[j] <- mean(I[(G - 100):G])
}
data.frame(rep = 1:J, val = store) %>%
    ggplot(aes(x = val)) +
    geom_density() +
    theme_classic()
