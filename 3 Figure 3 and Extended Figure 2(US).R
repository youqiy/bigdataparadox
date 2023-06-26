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
library(stringr) 
library(aspline)


# 1 Data Processing -------------------------------------------------------

temp1 <- new_us_vac %>%
  filter(Gender == "Male")

temp2 <- new_us_vac %>%
  filter(Gender == "Female")

df <- merge(x = temp1, y = temp2, by = "Week") %>%
  mutate(pop_total = pop_total.x + pop_total.y,
         n = n.x + n.y,
         pct_difference = pct_vaccinated.y - pct_vaccinated.x,
         pct_difference_act = pct_vaccinated_act.y - pct_vaccinated_act.x,
         pct_vaccinated = (pct_vaccinated.x * n.x)/n + (pct_vaccinated.y * n.y)/n,
         pct_vaccinated_act = (pct_vaccinated_act.x * pop_total.x)/pop_total + (pct_vaccinated_act.y * pop_total.y)/pop_total,
         pct_gender_act = pop_total.y/pop_total,
         pct_combined_act = (pct_vaccinated_act.y * pop_total.y)/pop_total) %>%
  select(Week, pop_total, n, pct_difference, pct_difference_act, pct_vaccinated, pct_vaccinated_act, pct_gender_act, pct_combined_act)

df <- df %>%
  mutate(
    ##gender difference
    #error
    error_difference = pct_difference - pct_difference_act,
    error_difference_less_10 = pct_difference - pct_difference_act*0.9,
    error_difference_less_5 = pct_difference - pct_difference_act*0.95,
    error_difference_plus_5 = pct_difference - pct_difference_act*1.05,
    error_difference_plus_10 = pct_difference - pct_difference_act*1.1,
    
    #drop out odds 
    DO_sqrt = sqrt((pop_total - n) / n),
         
    #variance
    var_difference = pct_vaccinated_act*(1-pct_vaccinated_act) + 4*pct_gender_act*(1-pct_gender_act) + 4*(pct_combined_act-pct_vaccinated_act*pct_gender_act),
    var_difference_less_10 = pct_vaccinated_act*0.9*(1-pct_vaccinated_act*0.9) + 4*pct_gender_act*(1-pct_gender_act) + 4*(pct_combined_act*0.9-pct_vaccinated_act*0.9*pct_gender_act),
    var_difference_less_5 = pct_vaccinated_act*0.95*(1-pct_vaccinated_act*0.95) + 4*pct_gender_act*(1-pct_gender_act) + 4*(pct_combined_act*0.95-pct_vaccinated_act*0.95*pct_gender_act),
    var_difference_plus_5 = pct_vaccinated_act*1.05*(1-pct_vaccinated_act*1.05) + 4*pct_gender_act*(1-pct_gender_act) + 4*(pct_combined_act*1.05-pct_vaccinated_act*1.05*pct_gender_act),
    var_difference_plus_10 = pct_vaccinated_act*1.1*(1-pct_vaccinated_act*1.1) + 4*pct_gender_act*(1-pct_gender_act) + 4*(pct_combined_act*1.1-pct_vaccinated_act*1.1*pct_gender_act),
         
    #standard deviation    
    sdG_difference = sqrt(var_difference),
    sdG_difference_less_10 = sqrt(var_difference_less_10),
    sdG_difference_less_5 = sqrt(var_difference_less_5),
    sdG_difference_plus_5 = sqrt(var_difference_plus_5),
    sdG_difference_plus_10 = sqrt(var_difference_plus_10),
         
    #ddc   
    ddc_difference = error_difference / (DO_sqrt * sdG_difference),
    ddc_difference_less_10 = error_difference_less_10 / (DO_sqrt * sdG_difference_less_10),
    ddc_difference_less_5 = error_difference_less_5 / (DO_sqrt * sdG_difference_less_5),
    ddc_difference_plus_5 = error_difference_plus_5 / (DO_sqrt * sdG_difference_plus_5),
    ddc_difference_plus_10 = error_difference_plus_10 / (DO_sqrt * sdG_difference_plus_10),
         
    #effective sample size
    neff_difference = (sdG_difference / error_difference)^2,
    neff_cap_difference = ifelse(neff_difference > n, n, neff_difference),
    neff_difference = (sdG_difference_less_10 / error_difference_less_10)^2,
    neff_cap_difference_less_10 = ifelse(neff_difference > n, n, neff_difference),
    neff_difference = (sdG_difference_less_5 / error_difference_less_5)^2,
    neff_cap_difference_less_5 = ifelse(neff_difference > n, n, neff_difference),
    neff_difference = (sdG_difference_plus_5 / error_difference_plus_5)^2,
    neff_cap_difference_plus_5 = ifelse(neff_difference > n, n, neff_difference),
    neff_difference = (sdG_difference_plus_10 / error_difference_plus_10)^2,
    neff_cap_difference_plus_10 = ifelse(neff_difference > n, n, neff_difference),
         
    #reduction percentange   
    red_difference = 1 - neff_cap_difference/n,
    red_difference_less_10 = 1 - neff_cap_difference_less_10/n,
    red_difference_less_5 = 1 - neff_cap_difference_less_5/n,
    red_difference_plus_5 = 1 - neff_cap_difference_plus_5/n,
    red_difference_plus_10 = 1 - neff_cap_difference_plus_10/n,
    
    ##overall estimate
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

# 2 Figure 3 --------------------------------------------------------------

df <- df %>%
  filter(Week >= as.Date("2021-02-14", "%Y-%m-%d") & Week <= as.Date("2021-05-15", "%Y-%m-%d")) %>%
  mutate(Study = "Difference")

fig <- list()

fig[["A"]] <- ggplot(df, aes(x = Week, y = error_difference)) +
  geom_line(aes(linetype = Study, col = Study), size=0.5) +
  geom_line(aes(y = error), size=0.5, col="#619CFF") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  labs(x = "", y = "Estimation Error") +
  scale_color_manual(limits = c("All",
                                "Difference",
                                "All_cdc",
                                "Difference_cdc"),
                     breaks = c("All",
                                "Difference",
                                "All_cdc",
                                "Difference_cdc"),
                     values  = c("All" = "#619CFF",
                                 "All_cdc" = "black",
                                 "Difference" = "#619CFF",
                                 "Difference_cdc" = "black"),
                     labels  = c("All" = "CTIS (overall rate)",
                                 "All_cdc" = "CDC (overall rate)",
                                 "Difference" = "CTIS (gender difference)",
                                 "Difference_cdc" = "CDC (gender difference)")) + 
  scale_linetype_manual(limits = c("All",
                                   "Difference",
                                   "All_cdc",
                                   "Difference_cdc"),
                        breaks = c("All",
                                   "Difference",
                                   "All_cdc",
                                   "Difference_cdc"),
                        values  = c("All" = "solid",
                                    "All_cdc" = "solid",
                                    "Difference" = "dotdash",
                                    "Difference_cdc" = "dotdash"),
                        labels  = c("All" = "CTIS (overall rate)",
                                    "All_cdc" = "CDC (overall rate)",
                                    "Difference" = "CTIS (gender difference)",
                                    "Difference_cdc" = "CDC (gender difference)")) +
  geom_ribbon(data = df, 
              aes(ymin = error_less_5, ymax = error_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = df, 
              aes(ymin = error_less_10, ymax = error_plus_10), 
              alpha=0.3, fill = "#619CFF") +  
  geom_ribbon(data = df, 
              aes(ymin = error_difference_less_5, ymax = error_difference_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = df, 
              aes(ymin = error_difference_less_10, ymax = error_difference_plus_10), 
              alpha=0.3, fill = "#619CFF") +
  theme(legend.position="none")

fig[["B"]] <- ggplot(df, aes(x = Week, y = sdG_difference)) +
  geom_line(linetype = "dotdash", col = "black", size=0.5) +
  geom_line(aes(y = sd_G), size=0.5, col="black") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  labs(x = "", y = expression(atop("Problem Difficulty", sigma[Y]))) +
  theme(legend.position="none")

fig[["C"]] <- ggplot(df, aes(x = Week, y = DO_sqrt)) +
  geom_line(size=0.5, col="#619CFF") +
  scale_y_continuous(limits = c(0, 100)) + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  labs(x = "", y = expression(atop("Data Quantity", sqrt((N - n) / n)))) +
  theme(legend.position="none")

fig[["D"]] <- ggplot(df, aes(x = Week, y = ddc_difference)) +
  geom_line(linetype = "dotdash", size=0.5, col="#619CFF") +
  geom_line(aes(y = ddc), size=0.5, col="#619CFF") +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  labs(x = "", y = expression(atop("Data Defect Correlation", paste(hat(rho)[paste("Y", ",", "R")])))) + 
  theme(legend.position="none")

fig[["E"]] <- ggplot(df, aes(x = Week, y = neff_cap_difference)) +
  geom_line(linetype = "dotdash", size=0.5, col="#619CFF") +
  geom_line(aes(y = n_eff_star_cap), size=0.5, col="#619CFF") +
  scale_y_continuous(trans='log10',
                     labels = scales::comma) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  labs(x = "", y = expression(n[eff])) +
  geom_ribbon(data = df, aes(ymin = n_eff_star_less_5, ymax = n_eff_star_cap_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(data = df, aes(ymin = neff_cap_difference_less_5, ymax = neff_cap_difference_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  theme(legend.position="none")

fig[["F"]] <- ggplot(df, aes(x = Week, y = 1 - red_difference)) +
  geom_line(linetype = "dotdash", size=0.5, col="#619CFF") +
  geom_line(aes(y = 1 - pct_red), size=0.5, col="#619CFF") +
  scale_y_continuous(trans='log10',
                     labels = scales::percent_format(accuracy = 0.01)) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  labs(x = "", y = expression(n[eff]/n)) +
  geom_ribbon(aes(ymin = 1 - pct_red_less_5, ymax = 1 - pct_red_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  geom_ribbon(aes(ymin = 1 - red_difference_less_5, ymax = 1 - red_difference_plus_5), 
              alpha=0.3, fill = "#619CFF") +
  theme(legend.position="none")

ggarrange((fig[["A"]] +
             font("ylab", size = 9)),
          (fig[["B"]] +
             font("ylab", size = 9)),
          (fig[["C"]] +
             font("ylab", size = 9)),
          (fig[["D"]] +
             font("ylab", size = 9)),
          (fig[["E"]] +
             font("ylab", size = 9)),
          (fig[["F"]] +
             font("ylab", size = 9)),
          nrow = 2, ncol = 3,
          hjust = -1,
          font.label = list(size = 8),
          labels =  c("A", "B", "C", "D", "E", "F"),
          align = "hv",
          common.legend = TRUE)

# 3 Extended Figure 2 -----------------------------------------------------

ggplot(temp_plot2, aes(x = Week, y = pct_vaccinated)) +
  geom_line(aes(linetype = Study, col = Study), size=0.5) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  geom_line(data = filter(temp_plot1, Study == "Male"),
            col = "#AECCFF", linetype = "solid", size=0.5) +  
  geom_segment(aes(x = as.Date("2021-02-14", "%Y-%m-%d"),
                   xend = as.Date("2021-03-28", "%Y-%m-%d"),
                   y = median(filter(temp_plot1, Study == "Male")$pct_vaccinated),
                   yend = median(filter(temp_plot1, Study == "Male")$pct_vaccinated)), 
               color="#AECCFF", linetype = "dashed", size=0.5) + 
  geom_line(data = filter(temp_plot1, Study == "Female"),
            col = "#156CFF", linetype = "solid", size=0.5) + 
  geom_segment(aes(x = as.Date("2021-02-14", "%Y-%m-%d"),
                   xend = as.Date("2021-03-28", "%Y-%m-%d"),
                   y = median(filter(temp_plot1, Study == "Female")$pct_vaccinated),
                   yend = median(filter(temp_plot1, Study == "Female")$pct_vaccinated)), 
               color="#156CFF", linetype = "dashed", size=0.5) + 
  geom_line(data = filter(temp_plot1, Study == "Male_cdc"),
            col = "#7A7A7A", linetype = "solid", size=0.5) + 
  geom_segment(aes(x = as.Date("2021-02-14", "%Y-%m-%d"),
                   xend = as.Date("2021-03-28", "%Y-%m-%d"),
                   y = median(filter(temp_plot1, Study == "Male_cdc")$pct_vaccinated),
                   yend = median(filter(temp_plot1, Study == "Male_cdc")$pct_vaccinated)), 
               color="#7A7A7A", linetype = "dashed", size=0.5) + 
  geom_line(data = filter(temp_plot1, Study == "Female_cdc"),
            col = "#333333", linetype = "solid", size=0.5) +
  geom_segment(aes(x = as.Date("2021-02-14", "%Y-%m-%d"),
                   xend = as.Date("2021-03-28", "%Y-%m-%d"),
                   y = median(filter(temp_plot1, Study == "Female_cdc")$pct_vaccinated),
                   yend = median(filter(temp_plot1, Study == "Female_cdc")$pct_vaccinated)), 
               color="#333333", linetype = "dashed", size=0.5) + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(size = 6)) +
  scale_color_manual(limits = c("Male",
                                "Female",
                                "Difference",
                                "Male_cdc",
                                "Female_cdc",
                                "Difference_cdc"),
                     breaks = c("Male",
                                "Female",
                                "Difference",
                                "Male_cdc",
                                "Female_cdc",
                                "Difference_cdc"),
                     values  = c("Male" = "#AECCFF",
                                 "Female" = "#156CFF",
                                 "Difference" = "#619CFF",                                
                                 "Male_cdc" = "#7A7A7A",
                                 "Female_cdc" = "#333333",
                                 "Difference_cdc" = "black"),
                     labels  = c("Male" = "CTIS (male)",
                                 "Female" = "CTIS (female)",
                                 "Difference" = "CTIS (gender difference)",
                                 "Male_cdc" = "CDC (male)",
                                 "Female_cdc" = "CDC (female)",
                                 "Difference_cdc" = "CDC (gender difference)")) +   
  scale_linetype_manual(limits = c("Male",
                                   "Female",
                                   "Difference",
                                   "Male_cdc",
                                   "Female_cdc",
                                   "Difference_cdc"),
                        breaks = c("Male",
                                   "Female",
                                   "Difference",
                                   "Male_cdc",
                                   "Female_cdc",
                                   "Difference_cdc"),
                        values  = c("Male" = "solid",
                                    "Female" = "solid",
                                    "Difference" = "dotdash",                                
                                    "Male_cdc" = "solid",
                                    "Female_cdc" = "solid",
                                    "Difference_cdc" = "dotdash"),
                        labels  = c("Male" = "CTIS (male)",
                                    "Female" = "CTIS (female)",
                                    "Difference" = "CTIS (gender difference)",
                                    "Male_cdc" = "CDC (male)",
                                    "Female_cdc" = "CDC (female)",
                                    "Difference_cdc" = "CDC (gender difference)")) +
  labs(x = "", y = "Vaccinated(%)")
