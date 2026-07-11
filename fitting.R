library(tidyverse)
setwd("~/Desktop/CME")

birth <- function(pX, mu, dt) rpois(1, pX * mu * dt)
death <- function(X, dt) rbinom(1, X, 1 - exp(-dt))
FOI <- function(beta, S, I, dt) {
    if (S == 0 && I == 0) return(0)
    rbinom(1, S, 1 - exp(-beta * I / (S + I) * dt))
}

singlerun <- function(S0, I0, beta, mu, p, G, dt) {
    S <- rep(0, G)
    I <- rep(0, G)
    S[1] <- S0
    I[1] <- I0
    for (i in 2:G) {
        muSI <- FOI(beta, S[i - 1], I[i - 1], dt)
        
        S[i] <- S[i - 1] - muSI
        I[i] <- I[i - 1] + muSI
        
        S[i] <- S[i] + birth(1 - p, mu, dt) - death(S[i], dt)
        I[i] <- I[i] + birth(p, mu, dt) - death(I[i], dt)
    }
    c(S, I)
}

plotsingle <- function(result, G) {
    data.frame(time = c(1:G, 1:G), count = result, state = c(rep("S", G), rep("I", G))) %>%
        mutate(state = factor(state, levels = c("S", "I"))) %>%
        ggplot(aes(x = time, y = count)) +
        geom_line() +
        facet_wrap(~state) +
        theme_classic()
}

beta <- seq(from=0,by=0.01,to=0.99)
mu <- seq(from=1,by=10,to=1000)
p <- seq(from=0.0001,by=0.01,to=1)

counts <- NULL
for (i in 1:100) {
    print(i)
    countstmp <- data.frame(t(read.table(paste0("JLresults/counts_",i,".csv"), sep=",")))
    # parstmp <- data.frame(t(read.table(paste0("writepars/pars_",i,".csv"), sep=",")))
    
    countstmp$beta <- beta[i]
    countstmp$mu <- 0
    countstmp$p <- 0
    for (j in 1:100) {
        for (k in 1:100) {
            idx <- 100 * (j - 1) + k
            countstmp$mu[idx] <- mu[j]
            countstmp$p[idx] <- p[k]
        }
    }
    counts <- rbind(counts, countstmp)
}
write_csv(counts, "JLresults/aggregated.csv")
counts <- read_csv("JLresults/aggregated.csv")

rowsums <- rowSums(counts %>% select(-c(beta, mu, p)))
counts <- counts %>% mutate_at(paste0("X", 1:100), \(x) x / rowsums)

counts[200, ] %>% 
    pivot_longer(starts_with("X"), names_to = "bin", values_to = "proportion") %>% 
    mutate(bin = as.numeric(str_remove(bin, "X")) - 1) %>%
    ggplot(aes(x = bin, y = proportion))+
    geom_point()+
    theme_classic()

finalbin <- counts %>% select(X21:X100) %>% rowSums
counts <- counts %>% select(-c(X21:X100)) %>% add_column(.before = "beta", X21 = finalbin)

cmemat <- read.table("eitan_data/host_cme_matrix.txt", header = T, sep = "\t")
cmemat$subject <- str_split(cmemat$host_id, "_", simplify = T)[, 1]
cmedat <- cmemat %>% 
    select(-host_id) %>% 
    pivot_longer(cols = -subject, names_to = "cme") %>% 
    group_by(cme, subject) %>% 
    summarize(count = sum(value)) %>% 
    ungroup

pseudo = 0.01
pars <- counts %>% select(beta:p)
cmes <- unique(cmedat$cme)
fits <- data.frame(cme = cmes, idx = 0, loglik = -Inf)
for (i in 1:length(cmes)) {
    print(i)
    tmpdat <- cmedat %>% filter(cme == cmes[i]) %>% pull(count)
    tmpdat <- table(tmpdat + 1)
    subset <- counts[, as.numeric(names(tmpdat))]
    loglik <- apply(subset, 1, \(x) sum(log(x + pseudo) * tmpdat))
    fits[i, ] <- c(cmes[i], which.max(loglik), max(loglik, na.rm = T))
}
write_csv(fits, "fits.csv")

pars[fits$idx, ] %>% 
    mutate(muI = mu * p, beta = beta + 0.001) %>%
    select(muI, beta) %>%
    ggplot(aes(x = log(muI), y = log(beta))) +
    geom_density_2d_filled() +
    annotate("text", x = -7, y = -1.8, label = "(0.0001, 0.11)", color = "white", family = "mono") + 
    annotate("text", x = -4.5, y = -0.8, label = "(0.052, 0.64)", color = "white", family = "mono") +
    annotate("text", x = -5.9, y = -4.2, label = "(0.015, 0.019)", color = "white", family = "mono") +
    annotate("text", x = -2.54, y = -6.67, label = "labels untransformed", color = "white", family = "mono") + 
    scale_fill_viridis_d(option = 14) +
    guides(fill="none") +
    ylab("log beta") + xlab("log muI") +
    theme_classic() +
    ggtitle("parameter fits across CMEs") +
    theme(
        text=element_text(size=12,family="mono"),
        axis.text.x=element_text(vjust=0),
        axis.title.y=element_text(vjust=3),
        axis.title.x=element_text(vjust=-1),
        legend.title=element_text(vjust=4),
        plot.margin=margin(r=15,t=15,l=15,b=15))

pars[fits$idx, ] %>% 
    filter(mu > 1 & p > 0.0001) %>% 
    mutate(muI = mu * p) %>%
    select(muI, beta) %>%
    ggplot(aes(x = muI, y = beta)) +
    geom_density_2d_filled() +
    scale_fill_viridis_d(option = 2)+ 
    guides(fill="none") +
    ylim(0, 1) + xlim(0, 1) +
    ylab("\u03B2") + xlab("Infected influx") +
    theme_classic()
pars[fits$idx, ] %>% 
    filter((mu == 1 | p == 0.0001) & beta > 0) %>%
    ggplot(aes(x = beta)) +
    geom_density() +
    xlim(0, 1) +
    xlab("\u03B2") +
    theme_classic()

tester <- "c4"
cmedat %>%
    filter(cme == tester) %>% 
    ggplot(aes(x = count)) +
    geom_bar() +
    xlim(-1,22)+
    xlab("bin")+
    theme_classic()
counts[fits$idx[fits$cme == tester], 1:21] %>% 
    pivot_longer(everything(), names_to = "bin", values_to="freq") %>% 
    mutate(bin=as.numeric(str_remove(bin,"X"))-1) %>% 
    ggplot(aes(x=bin, y=freq))+
    geom_col()+
    xlim(-1,22)+
    theme_classic()
testpars <- pars[fits$idx[fits$cme == tester], ]
testpars
plotsingle(singlerun(50, 0,
                     testpars$beta,
                     testpars$mu, 
                     testpars$p,
                     25000, 0.05), 25000)

pars[fits$idx, ]$mu %>% range
