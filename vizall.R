library(tidyverse)
setwd("~/Desktop/CME")

fits1 <- read_csv("fits_me_2.csv") %>% mutate(fit = "nb")
fits2 <- read_csv("fits_marc_2.csv") %>% mutate(fit = "gp")
fits3 <- read_csv("fits_marc_42.csv") %>% mutate(fit = "p")

fits <- rbind(fits1, fits2, fits3) %>% 
    mutate(fit = factor(fit, levels = c("p", "gp", "nb"))) %>% 
    group_by(fit) %>% 
    mutate(identified = !(R0l == min(R0l) & R0u == max(R0u)) &
               !(pil == min(pil) & piu == max(piu))) %>% 
    ungroup %>% mutate(identifiedfac = factor(identified, levels = c(T, F)))
fits <- left_join(fits, NI %>% group_by(cme) %>% summarize(N = mean(N)), by = "cme") %>% 
    mutate(R0m_eff = R0m * N, R0u_eff = R0u * N, R0l_eff = R0l * N) %>% 
    mutate(thresh = R0m_eff >= 1)

ordered <- fits %>% filter(fit == "p") %>% arrange(-R0m) %>% 
    mutate(idx = row_number()) %>% select(cme, idx)
fitsR0 <- left_join(fits, ordered, by = "cme")

fitsR0 %>%
    ggplot(aes(x = idx, y = R0m, color = identifiedfac)) +
    geom_linerange(aes(ymin = R0l, ymax = R0u), alpha = 0.15) +
    geom_point(alpha = 0.8, shape = 3) +
    scale_color_manual(name = "", labels = c("identifiable", "non"),
                       values = c("black", "#abc4c1")) +
    theme_classic() +
    ylab("R0") + xlab("cme") +
    facet_wrap(~fit) +
    theme(text = element_text(size = 12, family = "mono"))

ordered <- fits %>% filter(identified) %>% filter(fit == "p") %>% 
    arrange(-R0m) %>%  mutate(idx = row_number()) %>% select(cme, idx)
fitsR0_id <- left_join(fits %>% filter(identified), ordered, by = "cme")

fitsR0_id %>%
    ggplot(aes(x = idx, y = R0m, color = !thresh)) +
    geom_linerange(aes(ymin = R0l, ymax = R0u), alpha = 0.15) +
    geom_point(alpha = 0.8, shape = 3) +
    scale_color_manual(name = "R0 * E[N]", labels = c(">= 1", " < 1"),
                       values = c("black", "#abc4c1")) +
    theme_classic() +
    ylab("R0") + xlab("cme") +
    facet_wrap(~fit) +
    theme(text = element_text(size = 12, family = "mono"))

ordered <- fits %>% filter(fit == "p") %>% arrange(-pim) %>% 
    mutate(idx = row_number()) %>% select(cme, idx)
fitspi <- left_join(fits, ordered, by = "cme")

fitspi %>%
    ggplot(aes(x = idx, y = pim, color = identifiedfac)) +
    geom_linerange(aes(ymin = R0l, ymax = piu), alpha = 0.15) +
    geom_point(alpha = 0.8, shape = 3) +
    scale_color_manual(name = "", labels = c("identifiable", "non"),
                       values = c("black", "#abc4c1")) +
    theme_classic() +
    ylab("pi") + xlab("cme") +
    facet_wrap(~fit) +
    theme(text = element_text(size = 12, family = "mono"))

fits %>% filter(identified) %>% 
    ggplot(aes(x = R0m, y = pim)) +
    geom_point(shape = 1) +
    theme_classic() +
    xlab("R0") + ylab("pi") +
    facet_wrap(~fit) +
    theme(text = element_text(size = 12, family = "mono"))

fits %>% filter(identified) %>% 
    ggplot(aes(x = R0m, y = loglik)) +
    geom_point(shape = 1) +
    theme_classic() +
    xlab("R0") + 
    facet_wrap(~fit) +
    theme(text = element_text(size = 12, family = "mono"))
