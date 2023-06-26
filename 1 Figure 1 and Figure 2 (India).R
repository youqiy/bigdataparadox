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
library(cowplot)
library(gridExtra)
library(grid)


# 1 Data Processing -------------------------------------------------------

plot_ind_vac <- new_ind_vac %>%
  filter(Week >= as.Date("2021-05-16", "%Y-%m-%d") & Week <= as.Date("2021-09-18", "%Y-%m-%d")) %>%
  arrange(Study)

plot_ind_vac <- plot_ind_vac %>%
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

# 2 Figure 1 --------------------------------------------------------------

df <- plot_ind_vac

fig_ind <- list()

fig_ind[["EDA"]] <- ggplot(df, aes(x = Week, y = pct_vaccinated)) +
  geom_line(aes(col = Study, linetype = Study), size = 0.5) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"), axis.text = element_text(size = 6)) +
  geom_ribbon(data = filter(df3, Study == "cvoter"), 
              aes(ymin = ci_2.5_samp, ymax = ci_97.5_samp), 
              alpha=0.3, fill = "#F8766D") +
  geom_ribbon(data = filter(df3, Study == "fb_aggregated"), 
              aes(ymin = ci_2.5_samp, ymax = ci_97.5_samp), 
              alpha=0.3, fill = "#619CFF") +
  scale_color_manual(values = c("cvoter"= "#F8766D", 
                                "fb_aggregated" = "#619CFF", 
                                "cowin" = "black"),
                     labels = c("cvoter" = "CVoter (n=2,700)", 
                                "fb_aggregated" = "CTIS (n=25,000)", 
                                "cowin" = "CoWIN (benchmark)")) +  
  scale_linetype_manual(values = c("cvoter"= "solid", 
                                   "fb_aggregated" = "solid", 
                                   "cowin" = "dashed"),
                        labels = c("cvoter" = "CVoter (n=2,700)", 
                                   "fb_aggregated" = "CTIS (n=25,000)", 
                                   "cowin" = "CoWIN (benchmark)")) +
  labs(x = "", y = "Vaccinated(%)")

fig_ind[["Error"]] <- ggplot(df, aes(x = Week, y = error)) +
  geom_line(aes(col = Study, linetype = Study), size = 0.5) + 
  ylim(-0.1, 0.5) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  geom_ribbon(data = filter(df, Study == "fb_aggregated"), 
              aes(ymin = error_less_5, ymax = error_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "fb_aggregated"), aes(ymin = error_less_10, ymax = error_plus_10), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "cvoter"), 
              aes(ymin = error_less_5, ymax = error_plus_5), 
              alpha=0.3, fill = "#F8766D") +
  geom_ribbon(data = filter(df, Study == "cvoter"), 
              aes(ymin = error_less_10, ymax = error_plus_10), 
              alpha=0.3, fill = "#F8766D") +
  scale_color_manual(limits = c("cvoter", "fb_aggregated", "cowin"),
                     values = c("cvoter"= "#F8766D", 
                                "fb_aggregated" = "#619CFF",
                                "cowin" = "black"),
                     labels = c("cvoter"= "CVoter (n=2,700)", 
                                "fb_aggregated" = "CTIS (n=25,000)",
                                "cowin"  = "CoWIN (benchmark)"),
                     drop = FALSE) +
  scale_linetype_manual(limits = c("cvoter", "fb_aggregated", "cowin"),
                        values = c("cvoter"= "solid", 
                                   "fb_aggregated" = "solid",
                                   "cowin" = "dashed"),
                        labels = c("cvoter"= "CVoter (n=2,700)", 
                                   "fb_aggregated" = "CTIS (n=25,000)",
                                   "cowin"  = "CoWIN (benchmark)"),
                        drop = FALSE) +
  labs(x = "", y = "Estimation Error")

fig_ind[["sdG"]] <- ggplot(data = filter(df, Study == "fb_aggregated"), aes(x = Week, y = sd_G)) +
  geom_line(size = 0.5, linetype = "dashed", col = "black") + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  ylim(0.3, 0.5) +
  labs(x = "", y = expression(paste("Problem Difficulty", " ", sigma[Y])))

fig_ind[["DO"]] <- ggplot(df, aes(x = Week, y = DO_sqrt)) +
  geom_line(aes(col = Study), size=0.5) + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  scale_color_manual(limits = c("cvoter", "fb_aggregated", "cowin"),
                     values = c("cvoter"= "#F8766D", 
                                "fb_aggregated" = "#619CFF",
                                "cowin" = "black"),
                     drop = FALSE) + 
  ylim(0, 1000) +
  labs(x = "", y = expression(paste("Data Quantity", " ", sqrt((N - n) / n))))

fig_ind[["ddc"]] <- ggplot(df, aes(x = Week, y = ddc)) +
  geom_line(aes(col = Study), size=0.5) + 
  ylim(-0.005, 0.015) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  geom_ribbon(data = filter(df, Study == "fb_aggregated"), 
              aes(ymin = ddc_less_5, ymax = ddc_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "fb_aggregated"), 
              aes(ymin = ddc_less_10, ymax = ddc_plus_10), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "cvoter"), 
              aes(ymin = ddc_less_5, ymax = ddc_plus_5), 
              alpha=0.3, fill = "#F8766D") +
  geom_ribbon(data = filter(df, Study == "cvoter"), 
              aes(ymin = ddc_less_10, ymax = ddc_plus_10), 
              alpha=0.3, fill = "#F8766D") +
  scale_color_manual(values = c("cvoter"= "#F8766D", 
                                "fb_aggregated" = "#619CFF", 
                                "cowin" = "black"),
                     labels = c("CVoter", "Facebook", "CoWIN")) +
  labs(x = "", y = expression(paste("Data Defect Correlation", " ", paste(hat(rho)[paste("Y", ",", "R")]))))

fig_ind[["eff"]] <- df %>%
  ggplot(aes(x = Week, y = n_eff_star_cap)) +
  geom_line(aes(col = Study), size=0.5) + 
  ylim(0, 100) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  geom_ribbon(data = filter(df, Study == "fb_aggregated"), 
              aes(ymin = n_eff_star_less_5, ymax = n_eff_star_cap_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "cvoter"), 
              aes(ymin = n_eff_star_cap_less_5, ymax = n_eff_star_cap_plus_5), 
              alpha=0.3, fill = "#F8766D") +
  scale_color_manual(values = c("cvoter"= "#F8766D", 
                                "fb_aggregated" = "#619CFF", 
                                "cowin" = "black"),
                     labels = c("CVoter", "Facebook", "CoWIN")) +
  labs(x = "", y = expression(n[eff]))

fig_ind[["red"]] <- df %>%
  ggplot(aes(x = Week, y = pct_red)) +
  geom_line(aes(col = Study), size=0.5) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0.97, 1)) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  geom_ribbon(data = filter(df, Study == "fb_aggregated"), 
              aes(ymin = pct_red_less_5, ymax = pct_red_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "cvoter"), 
              aes(ymin = pct_red_less_5, ymax = pct_red_plus_5), 
              alpha=0.3, fill = "#F8766D") +
  scale_color_manual(values = c("cvoter"= "#F8766D", 
                                "fb_aggregated" = "#619CFF")) +
  labs(x = "", y = expression((n - n[eff])/n))

ggarrange((fig_ind[["EDA"]] +
             font("ylab", size = 9)),
          ggarrange((fig_ind[["Error"]] +
                       font("ylab", size = 9) +
                       theme(legend.position = 'none')),
                    (fig_ind[["sdG"]] +
                       font("ylab", size = 9) +
                       theme(legend.position = 'none')),
                    (fig_ind[["DO"]] +
                       font("ylab", size = 9) +
                       theme(legend.position = 'none')),
                    (fig_ind[["ddc"]] +
                       font("ylab", size = 9) +
                       theme(legend.position = 'none')),
                    (fig_ind[["eff"]] +
                       font("ylab", size = 9) +
                       theme(legend.position = 'none')),
                    (fig_ind[["red"]] + 
                       font("ylab", size = 9) +
                       theme(legend.position = 'none')),
                    nrow = 2, ncol = 3,
                    align = "hv",
                    font.label = list(size = 10),
                    labels =  c("B", "C", "D", "E", "F", "G"),
                    hjust = -1,
                    common.legend = FALSE),
          labels = "A",
          nrow = 1,
          ncol = 2,
          font.label = list(size = 10),
          hjust = -1,
          common.legend = TRUE)


# 3 Figure 2 --------------------------------------------------------------

# 3.1 Data Processing -----------------------------------------------------

temp <- plot_ind_vac %>%
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
  mutate(Study = "cowin")

dfc <- df2 %>%
  select(Week, Study, pct_vaccinated)

dfd <- df2 %>%
  select(Week, pct_vaccinated_act) %>%
  rename(pct_vaccinated = pct_vaccinated_act) %>%
  mutate(Study = "Difference_cowin")

dfe <- df3 %>%
  select(Week, Study, pct_vaccinated)

dff <- df3 %>%
  select(Week, pct_vaccinated_act) %>%
  rename(pct_vaccinated = pct_vaccinated_act) %>%
  mutate(Study = "Relative_cowin")

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
                                "cowin",
                                "Difference_cowin",
                                "Relative_cowin"),
                     values  = c("fb_aggregated" = "#619CFF",
                                 "Difference" = "#619CFF",
                                 "Relative" = "#619CFF",
                                 "cowin" = "black",
                                 "Difference_cowin" = "black",
                                 "Relative_cowin" = "black"),
                     labels  = c("fb_aggregated" = "CTIS (per week)",
                                 "Difference" = "CTIS (weekly difference)",
                                 "Relative" = "CTIS (relative weekly difference)",
                                 "cowin" = "CoWIN (per week)",
                                 "Difference_cowin" = "CoWIN (weekly difference)",
                                 "Relative_cowin" = "CoWIN (relative weekly difference)")) +  
  scale_linetype_manual(breaks = c("fb_aggregated",
                                   "Difference",
                                   "Relative",
                                   "cowin",
                                   "Difference_cowin",
                                   "Relative_cowin"),
                        values  = c("fb_aggregated" = "solid",
                                    "Difference" = "longdash",
                                    "Relative" = "dotted",
                                    "cowin" = "solid",
                                    "Difference_cowin" = "longdash",
                                    "Relative_cowin" = "dotted"),
                        labels  = c("fb_aggregated" = "CTIS (per week)",
                                    "Difference" = "CTIS (weekly difference)",
                                    "Relative" = "CTIS (relative weekly difference)",
                                    "cowin" = "CoWIN (per week)",
                                    "Difference_cowin" = "CoWIN (weekly difference)",
                                    "Relative_cowin" = "CoWIN (relative weekly difference)")) +
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
