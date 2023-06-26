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
library(betareg)


# 1 Data Processing -------------------------------------------------------

# 1.1 Axios-Ipsos ---------------------------------------------------------

#eda
ggplot(data = us_axios_acpt, aes(y = pct_hesitant_raw, x = pct_vaccinated)) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(y = "vaccine hesitancy", x = "vaccine uptake")

#beta regression
lm1 <- betareg(pct_hesitant_raw ~ pct_vaccinated, data = us_axios_acpt)
summary(lm1)

temp <- as.data.frame(us_axios_acpt$pct_vaccinated)
colnames(temp) <- "pct_vaccinated"
us_axios_acpt$pct_hesitant_fitted <- predict(lm1, newdata = temp)

temp <- as.data.frame(us_axios_acpt$pct_vaccinated_act)
colnames(temp) <- "pct_vaccinated"
us_axios_acpt$pct_hesitant_predicted <- predict(lm1, newdata = temp)

us_axios_acpt <- us_axios_acpt %>% 
  mutate(N = 255200373, 
         pct_hesitant = n/N*pct_hesitant_raw + pct_hesitant_predicted - n/N*pct_hesitant_fitted)

# 1.2 CTIS ----------------------------------------------------------------

#beta regression
temp <- as.data.frame(us_fb_acpt$pct_vaccinated)
colnames(temp) <- "pct_vaccinated"
us_fb_acpt$pct_hesitant_fitted <- predict(lm1, newdata = temp)

temp <- as.data.frame(us_fb_acpt$pct_vaccinated_act)
colnames(temp) <- "pct_vaccinated"
us_fb_acpt$pct_hesitant_predicted <- predict(lm1, newdata = temp)

us_fb_acpt <- us_fb_acpt %>%
  mutate(N = 255200373,
         pct_hesitant = n/N*pct_hesitant_raw + pct_hesitant_predicted - n/N*pct_hesitant_fitted)

#eda
df <- us_fb_acpt
ggplot(df, aes(x = Date, y = pct_hesitant)) +
  geom_line(col = "#619CFF", size=0.5) +
  scale_y_continuous(labels = scales::percent) + 
  geom_line(aes(y = pct_hesitant_raw), size = 0.5, linetype = "dashed", col = "#619CFF") +
  geom_line(data = us_axios_acpt, aes(y = pct_hesitant), size = 0.5, col = "#F8766D") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  labs(x = "", y = "% Hesitancy")

# 1.3 Error Decomposition -------------------------------------------------

temp <- us_axios_acpt %>%
  select(Date, pct_hesitant) %>%
  rename(pct_hesitant_act = pct_hesitant)

temp1 <- us_fb_acpt %>%
  select(Date, pct_hesitant, n, N) %>%
  rename(pop_total = N) %>%
  mutate(Study = "Adjusted")

temp2 <- us_fb_acpt %>%
  select(Date, pct_hesitant_raw, n, N) %>%
  rename(pop_total = N,
         pct_hesitant = pct_hesitant_raw) %>%
  mutate(Study = "Raw")

new_fb_acpt <- rbind(temp1, temp2)

us_acpt <- merge(temp, new_fb_acpt, by = "Date") %>%
  arrange(Study) %>%
  mutate(
    country = "US",
    
    #imprecision of benchmark
    pct_hesitant_act_less_10 = pct_hesitant_act * 0.9,
    pct_hesitant_act_less_5 = pct_hesitant_act * 0.95,
    pct_hesitant_act_plus_5 = pct_hesitant_act * 1.05,
    pct_hesitant_act_plus_10 = pct_hesitant_act * 1.1,
    
    #error
    error = pct_hesitant - pct_hesitant_act,
    error_less_10 = pct_hesitant - pct_hesitant_act_less_10,
    error_less_5 = pct_hesitant - pct_hesitant_act_less_5,
    error_plus_5 = pct_hesitant - pct_hesitant_act_plus_5,
    error_plus_10 = pct_hesitant - pct_hesitant_act_plus_10,
    
    #standard deviation
    sd_G = sqrt(pct_hesitant_act * (1 - pct_hesitant_act)),
    sd_G_less_10 = sqrt(pct_hesitant_act_less_10 * (1 - pct_hesitant_act_less_10)),
    sd_G_less_5 = sqrt(pct_hesitant_act_less_5 * (1 - pct_hesitant_act_less_5)),
    sd_G_plus_5 = sqrt(pct_hesitant_act_plus_5 * (1 - pct_hesitant_act_plus_5)),
    sd_G_plus_10 = sqrt(pct_hesitant_act_plus_10 * (1 - pct_hesitant_act_plus_10)),
    
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
    
    n_eff_lower_5 = pmin(n_eff_star_cap, n_eff_star_cap_less_5, n_eff_star_cap_plus_5),
    n_eff_higher_5 = pmax(n_eff_star_cap, n_eff_star_cap_less_5, n_eff_star_cap_plus_5),
    n_eff_lower_10 = pmin(n_eff_star_cap, n_eff_star_cap_less_10, n_eff_star_cap_plus_10),
    n_eff_higher_10 = pmax(n_eff_star_cap, n_eff_star_cap_less_10, n_eff_star_cap_plus_10),
    
    #standard errors, MoEs and CIs
    se_samp = sqrt(pct_hesitant * (1 - pct_hesitant) / n),
    MoE_samp = 2 * se_samp,
    ci_2.5_samp = pct_hesitant - MoE_samp,
    ci_97.5_samp = pct_hesitant + MoE_samp
  )

# 2 Figure 4 --------------------------------------------------------------

df <- us_acpt

fig_acpt <- list()

fig_acpt[["PCA"]] <- ggplot(df, aes(x = Date, y = pct_hesitant)) +
  geom_line(aes(col = Study, linetype = Study), size=0.5) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  geom_ribbon(data = filter(df, Study == "Raw"), aes(ymin = ci_2.5_samp, ymax = ci_97.5_samp), 
              alpha=0.3, fill = "#9BBB59") +
  geom_ribbon(data = filter(df, Study == "Adjusted"), aes(ymin = ci_2.5_samp, ymax = ci_97.5_samp), 
              alpha=0.3, fill = "#619CFF") +
  geom_line(data = filter(df, Study == "Adjusted"), aes(y = pct_hesitant_act), size=0.5, col = "#F8766D", linetype = "dashed") +
  scale_color_manual(breaks = c("axios",
                                "Adjusted",
                                "Raw"),
                     limits = c("axios",
                                "Adjusted",
                                "Raw"),
                     values = c("axios"= "#F8766D", 
                                "Adjusted" = "#619CFF", 
                                "Raw" = "#9BBB59"),
                     labels = c("axios" = "model-assisted Axios-Ipsos", 
                                "Adjusted" = "model-assisted CTIS", 
                                "Raw" = "raw CTIS")) +  
  scale_linetype_manual(breaks = c("axios",
                                   "Adjusted",
                                   "Raw"),
                        limits = c("axios",
                                   "Adjusted",
                                   "Raw"),
                        values = c("axios"= "dashed", 
                                   "Adjusted" = "solid", 
                                   "Raw" = "solid"),
                        labels = c("axios" = "model-assisted Axios-Ipsos", 
                                   "Adjusted" = "model-assisted CTIS", 
                                   "Raw" = "raw CTIS")) +  
  labs(x = "", y = "Hesitancy(%)") +
  theme(legend.position="none")

fig_acpt[["Error"]] <- ggplot(df, aes(x = Date, y = error)) +
  geom_line(aes(col = Study), size=0.5) + 
  scale_y_continuous(limits = c(-0.2, 0.05)) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  scale_color_manual(values = c("axios"= "#F8766D", 
                                "Adjusted" = "#619CFF", 
                                "Raw" = "#9BBB59"),
                     labels = c("Axios-Ipsos", "Facebook adjusted", "Facebook raw")) + 
  labs(x = "", y = "Estimation Error") +
  geom_ribbon(data = filter(df, Study == "Adjusted"), aes(ymin = error_less_5, ymax = error_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "Adjusted"), aes(ymin = error_less_10, ymax = error_plus_10), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "Raw"), aes(ymin = error_less_5, ymax = error_plus_5), 
              alpha=0.3, fill = "#9BBB59") +
  geom_ribbon(data = filter(df, Study == "Raw"), aes(ymin = error_less_10, ymax = error_plus_10), 
              alpha=0.3, fill = "#9BBB59") +
  theme(legend.position="none")

fig_acpt[["sdG"]] <- ggplot(data = filter(df, Study == "Adjusted"), aes(x = Date, y = sd_G)) +
  geom_line(size=0.5, linetype = "dashed", col = "#F8766D") + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  ylim(0.3, 0.5) +
  labs(x = "", y = expression(atop("Problem Difficulty", sigma[Y])))

fig_acpt[["DO"]] <- ggplot(filter(df, Study == "Raw"), aes(x = Date, y = DO_sqrt)) +
  geom_line(col = "#619CFF", size=0.5) + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6))  +
  ylim(0, 1000) +
  labs(x = "", y = expression(atop("Data Quantity", sqrt((N - n) / n)))) + 
  theme(legend.position="none")

fig_acpt[["ddc"]] <- ggplot(df, aes(x = Date, y = ddc)) +
  geom_line(aes(col = Study), size=0.5) + 
  ylim(-0.015, 0.005) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  scale_color_manual(values = c("axios"= "#F8766D", 
                                "Adjusted" = "#619CFF", 
                                "Raw" = "#9BBB59"),
                     labels = c("Axios-Ipsos", "Facebook adjusted", "Facebook raw")) + 
  geom_ribbon(data = filter(df, Study == "Adjusted"), aes(ymin = ddc_less_5, ymax = ddc_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "Adjusted"), aes(ymin = ddc_less_10, ymax = ddc_plus_10), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "Raw"), aes(ymin = ddc_less_5, ymax = ddc_plus_5), 
              alpha=0.3, fill = "#9BBB59") +
  geom_ribbon(data = filter(df, Study == "Raw"), aes(ymin = ddc_less_10, ymax = ddc_plus_10), 
              alpha=0.3, fill = "#9BBB59") +
  labs(x = "", y = expression(atop("Data Defect Correlation", paste(hat(rho)[paste("Y", ",", "R")])))) + 
  theme(legend.position="none")

fig_acpt[["eff"]] <- df %>%
  ggplot(aes(x = Date, y = n_eff_star_cap)) +
  geom_line(aes(col = Study), size=0.5) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  scale_y_continuous(trans='log10',
                     labels = scales::comma) +
  scale_color_manual(values = c("axios"= "#F8766D", 
                                "Adjusted" = "#619CFF", 
                                "Raw" = "#9BBB59"),
                     labels = c("Axios-Ipsos", "Facebook adjusted", "Facebook raw")) + 
  geom_ribbon(data = filter(df, Study == "Adjusted"), 
              aes(ymin = n_eff_lower_5,
                  ymax = n_eff_higher_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = filter(df, Study == "Raw"), 
              aes(ymin = n_eff_lower_5, ymax = n_eff_higher_5), 
              alpha=0.3, fill = "#9BBB59") +
  labs(x = "", y = expression(n[eff])) +
  theme(legend.position="none")

ggarrange((fig_acpt[["PCA"]] +
             font("ylab", size = 9)),
          (fig_acpt[["Error"]] +
             font("ylab", size = 9)),
          (fig_acpt[["sdG"]] +
             font("ylab", size = 9)),
          (fig_acpt[["DO"]] +
             font("ylab", size = 9)),
          (fig_acpt[["ddc"]] +
             font("ylab", size = 9)),
          (fig_acpt[["eff"]] +
             font("ylab", size = 9)),
          nrow = 2, ncol = 3,
          hjust = -0.5,
          font.label = list(size = 10),
          labels =  c("A", "B", "C", "D", "E", "F"),
          align = "hv",
          common.legend = TRUE)

