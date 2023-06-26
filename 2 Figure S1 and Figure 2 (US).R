library(tidyverse)
library(httr)
library(jsonlite)
library(lubridate)
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(lemon)
library(tidyr)
library(covidcast)
library(utils)
library(summarytools)
library(data.table)
library(ggrepel)
library(npreg)

# 1 Data Processing -------------------------------------------------------

plot_us_vac <- new_us_vac %>%
  filter(Week >= as.Date("2021-02-05", "%Y-%m-%d") & Week <= as.Date("2021-05-15", "%Y-%m-%d")) %>%
  arrange(Study)

plot_us_vac <- plot_us_vac %>%
  mutate(
    #imprecision of benchmark
    pct_vaccinated_act_less_10 = pct_vaccinated_act * 0.9,
    pct_vaccinated_act_less_5 = pct_vaccinated_act * 0.95,
    pct_vaccinated_act_plus_5 = pct_vaccinated_act * 1.05,
    pct_vaccinated_act_plus_10 = pct_vaccinated_act * 1.1,
    
    #error
    error = pct_vaccinated - pct_vaccinated_act,
    error_less_10 = pct_vaccinated - pct_vaccinated_act_less_10,
    error_less_5 = pct_vaccinated - pct_vaccinated_act_less_5,
    error_plus_5 = pct_vaccinated - pct_vaccinated_act_plus_5,
    error_plus_10 = pct_vaccinated - pct_vaccinated_act_plus_10,
    
    #standard deviation
    sd_G = sqrt(pct_vaccinated_act * (1 - pct_vaccinated_act)),
    sd_G_less_10 = sqrt(pct_vaccinated_act_less_10 * (1 - pct_vaccinated_act_less_10)),
    sd_G_less_5 = sqrt(pct_vaccinated_act_less_5 * (1 - pct_vaccinated_act_less_5)),
    sd_G_plus_5 = sqrt(pct_vaccinated_act_plus_5 * (1 - pct_vaccinated_act_plus_5)),
    sd_G_plus_10 = sqrt(pct_vaccinated_act_plus_10 * (1 - pct_vaccinated_act_plus_10)),
    
    #ddc
    ddc = error / (sqrt((pop_total - n) / n) * sd_G),
    ddc_less_10 = error_less_10 / (sqrt((pop_total - n) / n) * sd_G_less_10),
    ddc_less_5 = error_less_5 / (sqrt((pop_total - n) / n) * sd_G_less_5),
    ddc_plus_5 = error_plus_5 / (sqrt((pop_total - n) / n) * sd_G_plus_5),
    ddc_plus_10 = error_plus_10 / (sqrt((pop_total - n) / n) * sd_G_plus_10),
    
    #drop out odds
    DO = (pop_total - n) / n,
    DO_sqrt = sqrt(DO),
    
    #effective sample size
    n_eff_star = (sd_G / error)^2,
    n_eff_star_cap = ifelse(n_eff_star > n, n, n_eff_star),
    n_eff_star_less_10 = (sd_G_less_10 / error_less_10)^2,
    n_eff_star_cap_less_10 = ifelse(n_eff_star_less_10 > n, n, n_eff_star_less_10),
    n_eff_star_less_5 = (sd_G_less_5 / error_less_5)^2,
    n_eff_star_cap_less_5 = ifelse(n_eff_star_less_5 > n, n, n_eff_star_less_5),
    n_eff_star_plus_5 = (sd_G_plus_5 / error_plus_5)^2,
    n_eff_star_cap_plus_5 = ifelse(n_eff_star_plus_5 > n, n, n_eff_star_plus_5),
    n_eff_star_plus_10 = (sd_G_plus_10 / error_plus_10)^2,
    n_eff_star_cap_plus_10 = ifelse(n_eff_star_plus_10 > n, n, n_eff_star_plus_10),
    
    #standard errors, MoEs and CIs
    se_samp = sqrt(pct_vaccinated * (1 - pct_vaccinated) / n),
    MoE_samp = 2 * se_samp,
    ci_2.5_samp = pct_vaccinated - MoE_samp,
    ci_97.5_samp = pct_vaccinated + MoE_samp,
    
    #reduction percentage
    pct_red = 1 - (n_eff_star_cap / n),
    pct_red_less_10 = 1 - (n_eff_star_cap_less_10 / n),
    pct_red_less_5 = 1 - (n_eff_star_cap_less_5 / n),
    pct_red_plus_5 = 1 - (n_eff_star_cap_plus_5 / n),
    pct_red_plus_10 = 1 - (n_eff_star_cap_plus_10 / n)
  )

# 2 Extended Figure 1 -----------------------------------------------------

df <- plot_us_vac

fig_us <- list()

fig_us[["EDA"]] <- ggplot(df, aes(x = Week, y = pct_vaccinated)) +
  geom_line(aes(col = Study, linetype = Study), size=0.5) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  geom_ribbon(data = filter(df3, Study == "axios"), 
              aes(ymin = ci_2.5_samp, ymax = ci_97.5_samp), 
              alpha=0.3, fill = "#F8766D") +
  geom_ribbon(data = filter(df3, Study == "fb_aggregated"), 
              aes(ymin = ci_2.5_samp, ymax = ci_97.5_samp), 
              alpha=0.3, fill = "#619CFF") +
  scale_color_manual(breaks = c("axios", "fb_aggregated", "cdc"),
                     values = c("axios" = "#F8766D", 
                                "fb_aggregated" = "#619CFF", 
                                "cdc" = "black"),
                     labels = c("axios" = "Axios-Ipsos (n=1,000)", 
                                "fb_aggregated" = "CTIS (n=250,000)", 
                                "cdc" = "CDC (benchmark)"),
                     drop = FALSE) +
  scale_linetype_manual(breaks = c("axios", "fb_aggregated", "cdc"),
                        values = c("axios"= "solid", 
                                   "fb_aggregated" = "solid", 
                                   "cdc" = "dashed"),
                        labels = c("axios" = "Axios-Ipsos (n=1,000)", 
                                   "fb_aggregated" = "CTIS (n=250,000)", 
                                   "cdc" = "CDC (benchmark)"),
                        drop = FALSE) +
  labs(x = "", y = "Vaccinated(%)") +
  theme(legend.position="none")

fig_us[["Error"]] <- ggplot(df, aes(x = Week, y = error)) +
  geom_line(aes(col = Study), size=0.5) + 
  ylim(-0.1, 0.5) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  geom_ribbon(data = filter(df, Study == "fb_aggregated"), aes(ymin = error_less_5, ymax = error_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "fb_aggregated"), aes(ymin = error_less_10, ymax = error_plus_10), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "axios"), aes(ymin = error_less_5, ymax = error_plus_5), 
              alpha=0.3, fill = "#F8766D") +
  geom_ribbon(data = filter(df, Study == "axios"), aes(ymin = error_less_10, ymax = error_plus_10), 
              alpha=0.3, fill = "#F8766D") +
  scale_color_manual(values = c("axios"= "#F8766D", 
                                "fb_aggregated" = "#619CFF", 
                                "cdc" = "black"),
                     labels = c("Axios-Ipsos", "Facebook", "CDC")) +
  labs(x = "", y = "Estimation Error") +
  theme(legend.position="none")

fig_us[["sdG"]] <- ggplot(data = filter(df, Study == "fb_aggregated"), aes(x = Week, y = sd_G)) +
  geom_line(size=0.5, linetype = "dashed", col = "black") + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  ylim(0.3, 0.5) +
  labs(x = "", y = expression(atop("Problem Difficulty", sigma[Y])))

fig_us[["DO"]] <- ggplot(df, aes(x = Week, y = DO_sqrt)) +
  geom_line(aes(col = Study), size=0.5) + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  scale_color_manual(values = c("axios"= "#F8766D", 
                                "fb_aggregated" = "#619CFF", 
                                "cdc" = "black"),
                     labels = c("Axios-Ipsos", "Facebook", "CDC")) +
  ylim(0, 1000) +
  labs(x = "", y = expression(atop("Data Quantity", sqrt((N - n) / n)))) + 
  theme(legend.position="none")

fig_us[["ddc"]] <- ggplot(df, aes(x = Week, y = ddc)) +
  geom_line(aes(col = Study), size=0.5) + 
  ylim(-0.005, 0.015) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  geom_ribbon(data = filter(df, Study == "fb_aggregated"), aes(ymin = ddc_less_5, ymax = ddc_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "fb_aggregated"), aes(ymin = ddc_less_10, ymax = ddc_plus_10), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "axios"), aes(ymin = ddc_less_5, ymax = ddc_plus_5), 
              alpha=0.3, fill = "#F8766D") +
  geom_ribbon(data = filter(df, Study == "axios"), aes(ymin = ddc_less_10, ymax = ddc_plus_10), 
              alpha=0.3, fill = "#F8766D") +
  scale_color_manual(values = c("axios"= "#F8766D", 
                                "fb_aggregated" = "#619CFF", 
                                "cdc" = "black"),
                     labels = c("Axios-Ipsos", "Facebook", "CDC")) +
  labs(x = "", y = expression(atop("Data Defect Correlation", paste(hat(rho)[paste("Y", ",", "R")])))) + 
  theme(legend.position="none")

fig_us[["eff"]] <- df %>%
  ggplot(aes(x = Week, y = n_eff_star_cap)) +
  scale_y_continuous(trans='log10',
                     labels = scales::comma) +
  geom_line(aes(col = Study), size=0.5) + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  geom_ribbon(data = filter(df, Study == "fb_aggregated"), 
              aes(ymin = n_eff_lower_5, ymax = n_eff_higher_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "axios"), 
              aes(ymin = n_eff_lower_5, ymax = n_eff_higher_5), 
              alpha=0.3, fill = "#F8766D") +
  scale_color_manual(values = c("axios"= "#F8766D", 
                                "fb_aggregated" = "#619CFF", 
                                "cdc" = "black"),
                     labels = c("Axios-Ipsos", "Facebook", "CDC")) +
  labs(x = "", y = expression(n[eff])) + 
  theme(legend.position="none")

ggarrange((fig_us[["EDA"]] +
             font("ylab", size = 9)),
          (fig_us[["Error"]] +
             font("ylab", size = 9)),
          (fig_us[["sdG"]] +
             font("ylab", size = 9)),
          (fig_us[["DO"]] +
             font("ylab", size = 9)),
          (fig_us[["ddc"]] +
             font("ylab", size = 9)),
          (fig_us[["eff"]] +
             font("ylab", size = 9)),
          nrow = 2, ncol = 3,
          hjust = -0.5,
          font.label = list(size = 10),
          labels =  c("A", "B", "C", "D", "E", "F"),
          align = "hv",
          common.legend = TRUE)

# 3 Figure 2 --------------------------------------------------------------

# 3.1 Data Processing -----------------------------------------------------

temp <- plot_us_vac %>%
  filter(Study == "fb_aggregated")

# 3.1.1 Successive Difference ---------------------------------------------

temp$difference_act = c(0, (temp$pct_vaccinated_act[-1] - temp$pct_vaccinated_act[-nrow(temp)]))

temp$difference = c(0, (temp$pct_vaccinated[-1] - temp$pct_vaccinated[-nrow(temp)]))

#error
temp$difference_error = temp$difference - temp$difference_act
temp$difference_error_less_10 = temp$difference - temp$difference_act*0.9
temp$difference_error_less_5 = temp$difference - temp$difference_act*0.95
temp$difference_error_plus_5 = temp$difference - temp$difference_act*1.05
temp$difference_error_plus_10 = temp$difference - temp$difference_act*1.1

#standard deviation
temp$difference_sd = c(0, (temp$sd_G[-1]^2 + temp$sd_G[-nrow(temp)]^2))
temp$difference_sd_less_10 = c(0, (temp$sd_G_less_10[-1]^2 + temp$sd_G_less_10[-nrow(temp)]^2))
temp$difference_sd_less_5 = c(0, (temp$sd_G_less_5[-1]^2 + temp$sd_G_less_5[-nrow(temp)]^2))
temp$difference_sd_plus_5 = c(0, (temp$sd_G_plus_5[-1]^2 + temp$sd_G_plus_5[-nrow(temp)]^2))
temp$difference_sd_plus_10 = c(0, (temp$sd_G_plus_10[-1]^2 + temp$sd_G_plus_10[-nrow(temp)]^2))

#effective sample size
index <- (temp$difference_sd[-1])/(temp$difference_error[-1]^2)
max <- ifelse(temp$n[-1] > temp$n[-nrow(temp)], temp$n[-1], temp$n[-nrow(temp)])
temp$difference_eff = c(0, ifelse(index > max, max, index))

index <- (temp$difference_sd_less_10[-1])/(temp$difference_error_less_10[-1]^2)
max <- ifelse(temp$n[-1] > temp$n[-nrow(temp)], temp$n[-1], temp$n[-nrow(temp)])
temp$difference_eff_less_10 = c(0, ifelse(index > max, max, index))

index <- (temp$difference_sd_less_5[-1])/(temp$difference_error_less_5[-1]^2)
max <- ifelse(temp$n[-1] > temp$n[-nrow(temp)], temp$n[-1], temp$n[-nrow(temp)])
temp$difference_eff_less_5 = c(0, ifelse(index > max, max, index))

index <- (temp$difference_sd_plus_5[-1])/(temp$difference_error_plus_5[-1]^2)
max <- ifelse(temp$n[-1] > temp$n[-nrow(temp)], temp$n[-1], temp$n[-nrow(temp)])
temp$difference_eff_plus_5 = c(0, ifelse(index > max, max, index))

index <- (temp$difference_sd_plus_10[-1])/(temp$difference_error_plus_10[-1]^2)
max <- ifelse(temp$n[-1] > temp$n[-nrow(temp)], temp$n[-1], temp$n[-nrow(temp)])
temp$difference_eff_plus_10 = c(0, ifelse(index > max, max, index))

#reduction percentage
temp$difference_red = c(0, (1 - temp$difference_eff[-1]/max))
temp$difference_red_less_10 = c(0, (1 - temp$difference_eff_less_10[-1]/max))
temp$difference_red_less_5 = c(0, (1 - temp$difference_eff_less_5[-1]/max))
temp$difference_red_plus_5 = c(0, (1 - temp$difference_eff_plus_5[-1]/max))
temp$difference_red_plus_10 = c(0, (1 - temp$difference_eff_plus_10[-1]/max))

# 3.1.2 Relative Successive Difference ------------------------------------

temp$relative_act = c(1, (temp$pct_vaccinated_act[-1] - temp$pct_vaccinated_act[-nrow(temp)])/temp$pct_vaccinated_act[-nrow(temp)])
temp$relative = c(1, (temp$pct_vaccinated[-1] - temp$pct_vaccinated[-nrow(temp)])/temp$pct_vaccinated[-nrow(temp)])

#error
temp$relative_error = temp$relative - temp$relative_act

temp$ratio_act = c(1, (temp$pct_vaccinated_act[-1])/temp$pct_vaccinated_act[-nrow(temp)])

#standard deviation
temp$relative_sd = c(0, (temp$ratio_act[-1])^2 * ((temp$sd_G[-1]^2)/(temp$pct_vaccinated_act[-1]^2) + (temp$sd_G[-nrow(temp)]^2)/(temp$pct_vaccinated_act[-nrow(temp)]^2)))

temp$relative_sd_less_10 = c(0, (temp$ratio_act[-1])^2 * ((temp$sd_G_less_10[-1]^2)/(temp$pct_vaccinated_act[-1]^2) + (temp$sd_G_less_10[-nrow(temp)]^2)/(temp$pct_vaccinated_act[-nrow(temp)]^2))) * (1/0.9) * (1/0.9)

temp$relative_sd_less_5 = c(0, (temp$ratio_act[-1])^2 * ((temp$sd_G_less_5[-1]^2)/(temp$pct_vaccinated_act[-1]^2) + (temp$sd_G_less_5[-nrow(temp)]^2)/(temp$pct_vaccinated_act[-nrow(temp)]^2))) * (1/0.95) * (1/0.95)

temp$relative_sd_plus_5 = c(0, (temp$ratio_act[-1])^2 * ((temp$sd_G_plus_5[-1]^2)/(temp$pct_vaccinated_act[-1]^2) + (temp$sd_G_plus_5[-nrow(temp)]^2)/(temp$pct_vaccinated_act[-nrow(temp)]^2))) * (1/1.05) * (1/1.05)

temp$relative_sd_plus_10 = c(0, (temp$ratio_act[-1])^2 * ((temp$sd_G_plus_10[-1]^2)/(temp$pct_vaccinated_act[-1]^2) + (temp$sd_G_plus_10[-nrow(temp)]^2)/(temp$pct_vaccinated_act[-nrow(temp)]^2))) * (1/1.1) * (1/1.1)

#effective sample size
index <- (temp$relative_sd[-1])/(temp$relative_error[-1]^2)
max <- ifelse(temp$n[-1] > temp$n[-nrow(temp)], temp$n[-1], temp$n[-nrow(temp)])
temp$relative_eff = c(0, ifelse(index > max, max, index))

index <- (temp$relative_sd_less_10[-1])/(temp$relative_error[-1]^2)
max <- ifelse(temp$n[-1] > temp$n[-nrow(temp)], temp$n[-1], temp$n[-nrow(temp)])
temp$relative_eff_less_10 = c(0, ifelse(index > max, max, index))

index <- (temp$relative_sd_less_5[-1])/(temp$relative_error[-1]^2)
max <- ifelse(temp$n[-1] > temp$n[-nrow(temp)], temp$n[-1], temp$n[-nrow(temp)])
temp$relative_eff_less_5 = c(0, ifelse(index > max, max, index))

index <- (temp$relative_sd_plus_5[-1])/(temp$relative_error[-1]^2)
max <- ifelse(temp$n[-1] > temp$n[-nrow(temp)], temp$n[-1], temp$n[-nrow(temp)])
temp$relative_eff_plus_5 = c(0, ifelse(index > max, max, index))

index <- (temp$relative_sd_plus_10[-1])/(temp$relative_error[-1]^2)
max <- ifelse(temp$n[-1] > temp$n[-nrow(temp)], temp$n[-1], temp$n[-nrow(temp)])
temp$relative_eff_plus_10 = c(0, ifelse(index > max, max, index))

#reduction percentage
temp$relative_red = c(0, (1 - temp$relative_eff[-1]/max))
temp$relative_red_less_10 = c(0, (1 - temp$relative_eff_less_10[-1]/max))
temp$relative_red_less_5 = c(0, (1 - temp$relative_eff_less_5[-1]/max))
temp$relative_red_plus_5 = c(0, (1 - temp$relative_eff_plus_5[-1]/max))
temp$relative_red_plus_10 = c(0, (1 - temp$relative_eff_plus_10[-1]/max))

df <- rbind(temp[-1, ])

# 3.1.3 Renaming ----------------------------------------------------------

df1 <- df %>%
  select(Week, Study, pct_vaccinated, pct_vaccinated_act,
         n_eff_star_cap, n_eff_star_cap_less_10, n_eff_star_cap_less_5, 
         n_eff_star_cap_plus_5, n_eff_star_cap_plus_10,
         pct_red, pct_red_less_10, pct_red_less_5, pct_red_plus_5, pct_red_plus_10)

df2 <- df %>%
  select(Week, difference, difference_act,
         difference_eff, difference_eff_less_10, difference_eff_less_5,
         difference_eff_plus_5, difference_eff_plus_10,
         difference_red, difference_red_less_10, difference_red_less_5, 
         difference_red_plus_5, difference_red_plus_10) %>%
  rename(pct_vaccinated = difference,
         pct_vaccinated_act = difference_act,
         n_eff_star_cap = difference_eff,
         n_eff_star_cap_less_10 = difference_eff_less_10,
         n_eff_star_cap_less_5 = difference_eff_less_5,
         n_eff_star_cap_plus_5 = difference_eff_plus_5,
         n_eff_star_cap_plus_10 = difference_eff_plus_10,
         pct_red = difference_red,
         pct_red_less_10 = difference_red_less_10,
         pct_red_less_5 = difference_red_less_5,
         pct_red_plus_5 = difference_red_plus_5,
         pct_red_plus_10 = difference_red_plus_10) %>%
  mutate(Study = "Difference")

df3 <- df %>%
  select(Week, relative, relative_act,
         relative_eff, relative_eff_less_10, relative_eff_less_5,
         relative_eff_plus_5, relative_eff_plus_10,
         relative_red, relative_red_less_10, relative_red_less_5, 
         relative_red_plus_5, relative_red_plus_10) %>%
  rename(pct_vaccinated = relative,
         pct_vaccinated_act = relative_act,
         n_eff_star_cap = relative_eff,
         n_eff_star_cap_less_10 = relative_eff_less_10,
         n_eff_star_cap_less_5 = relative_eff_less_5,
         n_eff_star_cap_plus_5 = relative_eff_plus_5,
         n_eff_star_cap_plus_10 = relative_eff_plus_10,
         pct_red = relative_red,
         pct_red_less_10 = relative_red_less_10,
         pct_red_less_5 = relative_red_less_5,
         pct_red_plus_5 = relative_red_plus_5,
         pct_red_plus_10 = relative_red_plus_10) %>%
  mutate(Study = "Relative")

temp1 <- rbind(df1, df2, df3) %>%
  mutate(n_eff_lower_5 = pmin(n_eff_star_cap, n_eff_star_cap_less_5, n_eff_star_cap_plus_5),
         n_eff_higher_5 = pmax(n_eff_star_cap, n_eff_star_cap_less_5, n_eff_star_cap_plus_5),
         n_eff_lower_10 = pmin(n_eff_star_cap, n_eff_star_cap_less_10, n_eff_star_cap_plus_10),
         n_eff_higher_10 = pmax(n_eff_star_cap, n_eff_star_cap_less_10, n_eff_star_cap_plus_10),
         
         pct_red_lower_10 =  pmin(pct_red, pct_red_less_10, pct_red_plus_10),
         pct_red_higher_10 =  pmax(pct_red, pct_red_less_10, pct_red_plus_10),
         pct_red_lower_5 =  pmin(pct_red, pct_red_less_5, pct_red_plus_5),
         pct_red_higher_5 =  pmax(pct_red, pct_red_less_5, pct_red_plus_5))

dfa <- df1 %>%
  select(Week, Study, pct_vaccinated)

dfb <- df1 %>%
  select(Week, pct_vaccinated_act) %>%
  rename(pct_vaccinated = pct_vaccinated_act) %>%
  mutate(Study = "cdc")

dfc <- df2 %>%
  select(Week, Study, pct_vaccinated)

dfd <- df2 %>%
  select(Week, pct_vaccinated_act) %>%
  rename(pct_vaccinated = pct_vaccinated_act) %>%
  mutate(Study = "Difference_cdc")

dfe <- df3 %>%
  select(Week, Study, pct_vaccinated)

dff <- df3 %>%
  select(Week, pct_vaccinated_act) %>%
  rename(pct_vaccinated = pct_vaccinated_act) %>%
  mutate(Study = "Relative_cdc")

temp2 <- rbind(dfa, dfb, dfc, dfd, dfe, dff)

# 3.2 Plotting ------------------------------------------------------------

ggplot(temp2, aes(x = Week, y = pct_vaccinated)) +
  geom_line(aes(linetype = Study, col = Study), size=0.5) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"), axis.text = element_text(size = 6)) +
  scale_color_manual(breaks = c("fb_aggregated",
                                "Difference",
                                "Relative",
                                "cdc",
                                "Difference_cdc",
                                "Relative_cdc"),
                     values  = c("fb_aggregated" = "#619CFF",
                                 "Difference" = "#619CFF",
                                 "Relative" = "#619CFF",
                                 "cdc" = "black",
                                 "Difference_cdc" = "black",
                                 "Relative_cdc" = "black"),
                     labels  = c("fb_aggregated" = "CTIS (per week)",
                                 "Difference" = "CTIS (weekly difference)",
                                 "Relative" = "CTIS (relative weekly difference)",
                                 "cdc" = "CDC (per week)",
                                 "Difference_cdc" = "CDC (weekly difference)",
                                 "Relative_cdc" = "CDC (relative weekly difference)")) +  
  scale_linetype_manual(breaks = c("fb_aggregated",
                                   "Difference",
                                   "Relative",
                                   "cdc",
                                   "Difference_cdc",
                                   "Relative_cdc"),
                        values  = c("fb_aggregated" = "solid",
                                    "Difference" = "longdash",
                                    "Relative" = "dotted",
                                    "cdc" = "solid",
                                    "Difference_cdc" = "longdash",
                                    "Relative_cdc" = "dotted"),
                        labels  = c("fb_aggregated" = "CTIS (per week)",
                                    "Difference" = "CTIS (weekly difference)",
                                    "Relative" = "CTIS (relative weekly difference)",
                                    "cdc" = "CDC (per week)",
                                    "Difference_cdc" = "CDC (weekly difference)",
                                    "Relative_cdc" = "CDC (relative weekly difference)")) +
  labs(x = "", y = "Vaccinated(%)")

temp1 %>%
  filter(Study == "fb_aggregated") %>%
  ggplot(aes(x = Week, y = n_eff_star_cap)) +
  scale_y_continuous(trans='log10',
                     limits = c(1, 250000),
                     labels = scales::comma) +
  geom_line(linetype = "solid", size = 0.5, col = "#619CFF") +
  geom_ribbon(aes(ymin = n_eff_lower_5, ymax = n_eff_higher_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_line(data = filter(temp1, Study == "Difference"),
            aes(y = n_eff_star_cap), size = 0.5, linetype = "longdash", col = "#619CFF") +
  geom_ribbon(data = filter(temp1, Study == "Difference"),
              aes(ymin = n_eff_lower_5, ymax = n_eff_higher_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_line(data = filter(temp1, Study == "Relative"),
            aes(y = n_eff_star_cap), size = 0.5, linetype = "dotted", col = "#619CFF") +
  geom_ribbon(data = filter(temp1, Study == "Relative"),
              aes(ymin = n_eff_lower_5, ymax = n_eff_higher_5), 
              alpha=0.3, fill = "#619CFF") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6),
        legend.position = "top") +
  labs(x = "", y = expression(n[eff]))
