library(tidyverse)
setwd("~/Desktop/CME")

singlerun <- function(alpha, f, R0, G, S, I) {
    Ss <- matrix(0, nrow = 168, ncol = G)
    Is <- matrix(0, nrow = 168, ncol = G)
    t <- rep(0, 168)
    g <- 0
    
    while (g < G) {
        pop = S + I
        if (pop > 0 & pop <= 168) {
            if (g >= t[pop]) {
                t[pop] <- t[pop] + 1
                Ss[pop, t[pop]] <- S
                Is[pop, t[pop]] <- I
            }
        }
        if (S + I == 0) { foi = 0
        } else {foi <- R0 * S * I / (S + I)}
        rate <- alpha + foi + S + I
        g <- g - 1 / rate * log(1 - runif(1))
        pick <- runif(1)
        if (pick < alpha * (1 - f) / rate) { S <- S + 1
        } else if (pick < alpha / rate) { I <- I + 1
        } else if (pick < (alpha + foi) / rate) { 
            S <- S - 1
            I <- I + 1
        } else if (pick < (alpha + foi + S) / rate) { S <- S - 1
        } else {I <- I - 1}
    }
    list(Ss, Is)
}

cmemat <- read.table("eitan_data/host_cme_matrix.txt", header = T, sep = "\t")
cmes <- colnames(cmemat %>% select(-host_id))
cmemat$subject <- str_split(cmemat$host_id, "_", simplify = T)[, 1]
cmemat$genome <- str_split(cmemat$host_id, "_", simplify = T)[, 2]
hostdat <-  read.table("eitan_data/hosts.txt", header = T, sep = "\t")
cmemat <- left_join(cmemat, hostdat %>% select(species, genus, family, host_id), by = "host_id")
cmetab <- read.table("eitan_data/cme_table.txt", header = T, sep = "\t")

cmedat <- cmemat %>% 
    group_by(subject, species) %>% 
    summarize_at(cmes, \(x) sum(x) != 0) %>% 
    select(-species) %>% 
    pivot_longer(cols = -subject, names_to = "cme") %>% 
    group_by(cme, subject) %>% 
    summarize(count = sum(value), .groups = "drop")

alpha <- rep(0, length(cmes))
for (i in 1:length(cmes)) {
    fams <- cmemat %>% filter(.[[cmes[i]]]) %>% pull(family) %>% unique
    alpha[i] <- cmemat %>% 
        filter(family %in% fams) %>% 
        group_by(subject) %>% 
        select(subject, species) %>%
        distinct %>% 
        summarize(n = n()) %>% 
        pull(n) %>% 
        mean
}
alpha <- data.frame(alpha = round(alpha, 0), cme = cmes)
alphalist <- sort(unique(alpha$alpha))
alphaidx <- list(tmp1 = alphalist[1:10],
                 tmp2 = alphalist[11:20],
                 tmp3 = alphalist[21:30],
                 tmp4 = alphalist[31:40],
                 tmp5 = alphalist[41:45])

fits <- read_csv("conditional_fits.csv")

NI <- data.frame(matrix(0, N * length(cmes), 3))
for (i in 1:length(cmes)) {
    print(i)
    fams <- cmemat %>% filter(.[[cmes[i]]]) %>% pull(family) %>% unique
    Nstats <- cmemat %>% 
        filter(family %in% fams) %>% 
        group_by(subject) %>% 
        select(subject, species) %>%
        distinct %>% 
        summarize(N = n())
    Istats <- cmemat %>% filter(.[[cmes[i]]]) %>% group_by(subject, species) %>% 
        summarize(n = n()) %>% summarize(I = n())
    statistics <- left_join(Nstats, Istats, by = "subject")
    statistics <- rbind(statistics,
        data.frame(subject = unique(hostdat$aid), N = 0, I = NA)[!(unique(hostdat$aid) %in% statistics$subject), ])
    statistics[is.na(statistics)] <- 0
    NI[N * (i - 1) + (1:N), ] <- statistics %>% select(N, I) %>% mutate(cme = cmes[i])
}

fits %>% 
    ggplot(aes(x = R0)) +
    geom_histogram() +
    theme_classic()

tester <- "c4"
parfit <- fits[fits$cme == tester, 2:4]
counts <- read_csv(paste0("alpharesults/alpha_",
                          which(alphalist ==  min(alpha$alpha[alpha$cme == tester], 94)),
                          ".csv"))
counts <- counts %>% filter(f == parfit$f, R0 == parfit$R0)
gatherdat <- colSums(counts[NI[NI$X3 == tester, 1], ])
gatherdat <- gatherdat / sum(gatherdat)

cmedat %>%
    filter(cme == tester) %>%
    pull(count) %>% 
    factor(levels = 0:(O-1)) %>% 
    table %>% 
    data.frame %>% 
    rename(count = '.', freq = Freq) %>% 
    mutate(count = as.numeric(count) - 1, 
           freq = freq / sum(freq)) %>% 
    ggplot(aes(x=count, y=freq)) +
    geom_col(width = 1, color = "black", fill = "white")+
    ylim(0,1)+xlim(-1,51)+
    theme_classic() +
    ggtitle(tester) +
    theme(text=element_text(size=12,family="mono"),
          legend.position = "bottom")

data.frame(count = 0:(O - 1), freq = gatherdat[1:O]) %>% 
    ggplot(aes(x=count, y=freq)) +
    geom_col(width = 1, color = "black", fill = "white")+
    ylim(0,1)+xlim(-1,51)+
    theme_classic() +
    ggtitle("best fit") +
    theme(text=element_text(size=12,family="mono"),
          legend.position = "bottom")


testalpha <- 69
testcounts <- read_csv(paste0("alpharesults/alpha_", which(alphalist == testalpha),".csv"))
testcounts[, 1:O] <- testcounts[, 1:O] + 1 / 10
rowsums <- rowSums(testcounts[, 1:O])
testcounts[, 1:O] <- testcounts[, 1:O] %>% 
    mutate_all(\(x) x / rowsums)

condN <- rpois(N, testalpha)

samplepars <- rbind(testcounts %>% filter(R0 == 0, f <= 0.1) %>% select(f, R0) %>% sample_n(100),
      testcounts %>% filter(R0 > 0, f <= 0.1) %>% select(f, R0) %>% sample_n(100))
# samplepars <- testcounts %>% select(f, R0) %>% sample_n(100)
testret <- matrix(0, nrow=200, ncol=3)
for (i in 1:200) {
    print(i)
    sim <- apply((testcounts %>% filter(f == samplepars[i, ]$f, R0 == samplepars[i, ]$R0))[condN, 1:O],
                 1, \(x)  which.max(rmultinom(1, 1, x)) - 1)
    statistics <- data.frame(N = condN, I = sim)
        
    factor <- 0
    loglik <- rep(0, 100*100)
    loglik0 <- rep(0, 100)
    for (k in unique(statistics$N)) {
        slimdata <- statistics %>% filter(N == k)
        slimI <- as.numeric(unname(table(factor(slimdata$I, levels = 0:(O-1)))))
        factor <- factor + log(factorial(sum(slimI))) - sum(log(factorial(slimI)))
        
        slimfreq <- testcounts %>% filter(N == k)
        slimfreq0 <- testcounts %>% filter(N == k, R0 == 0)
        
        loglik <- loglik + apply(slimfreq[, 1:O], 1, \(x) sum(slimI * log(x)))
        loglik0 <- loglik0 + apply(slimfreq0[, 1:O], 1, \(x) sum(slimI * log(x)))
    }
        
    testret[i, 3] <- pchisq(2 * abs(max(loglik) - max(loglik0)), 
                           1, lower.tail = F, log.p = F) / 2
        
    # lambda <- mean(rep(0:(O - 1), slimI))
    # exp <- N * dpois(0:(lO-1), lambda) + 0.0001
    # expected <- N * exp / sum(exp)
    # goodness[i, j] <- pchisq(sum((slimI - expected) ^ 2 / expected), 
    #                       O - 2, lower.tail = F, log.p = F)
        
    idx <- which.max(loglik)

    testret[i, 1] <- slimfreq$f[idx]
    testret[i, 2] <- slimfreq$R0[idx]
}

testcounts <- testcounts %>% 
    summarize(.by = c(starts_with("X"), N),
              alpha = alpha[1],
              f = f[1],
              R0 = R0[1],
              identifiable = n() == 1)

testcounts %>% ggplot(aes(x = f, y = R0, fill = identifiable))+
    geom_tile()+
    scale_fill_manual(values = c("black","white")) +
    theme_classic()

data.frame(cbind(testret, samplepars)) %>%
    ggplot(aes(x = R0, y = X2, color = f)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, color = "red", linetype = 2) +
    ylab("estimate") +
    theme_classic()
data.frame(cbind(testret, samplepars)) %>%
    ggplot(aes(x = f, y = X1, color = R0)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, color = "red", linetype = 2) +
    ylab("estimate") +
    theme_classic()

data.frame(cbind(testret, samplepars)) %>%
    mutate(is0 = ifelse(R0 == 0, "R0 = 0", "R0 > 0")) %>% 
    ggplot(aes(x = X3)) +
    geom_histogram() +
    geom_histogram(color = "black", fill = "white") + 
    facet_wrap(~is0) +
    xlab("p-value") +
    theme_classic() +
    theme(text = element_text(size = 12, family = "mono"))

# data.frame(cbind(testret, samplepars)) %>%
#     mutate(is0 = ifelse(R0 == 0, "R0 = 0", "R0 > 0")) %>% 
#     ggplot(aes(x = value)) +
#     geom_histogram() +
#     facet_wrap(~is0) +
#     theme_classic()
