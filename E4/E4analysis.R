library(dplyr)
library(ggplot2)
library(tidyr)
library(rlang)
library(purrr)
library(gridExtra)
library(tinytex)
library(stringr)

load("stats_matrix_e4.RData")

model_type_labels <- c(
  "adcat" = "AC",
  "propodds" = "CL",
  "stoprat" = "SR"
)

stats_by_group$Sampling_Method <- factor(stats_by_group$Sampling_Method,
                                         levels = c("SRS", "Best-SRS",
                                                    "Z-Bal", "GA-Z-Bal",
                                                    "Y-Bal", "GA-Y-Bal",
                                                    "YZ-Bal", "GA-YZ-Bal"))
# Web Figure 17
p1 <- stats_by_group %>%
  filter(G_vec == "G1", Experiment == 4, Sampling_Method %in% c("SRS",
                                                                "Z-Bal",
                                                                "Y-Bal",
                                                                "YZ-Bal"), BETA1 == 25) %>%
  group_by(Sampling_Method, Analysis_Method, Sample_Size) %>%
  ggplot(aes(x = Sample_Size, y = Mean_Absolute_Relative_Bias, shape = Analysis_Method, color = Analysis_Method, linetype = "Incomplete Data")) +
  geom_line(linewidth = 1.2) +
  geom_line(aes(y = Mean_Absolute_Relative_Bias_Complete, linetype = "Complete Data"), linewidth = 0.5, alpha = 0.4) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_bw() +
  facet_grid(Model_Type ~ Sampling_Method, labeller = labeller(Model_Type = model_type_labels)) + 
  labs(
    title = "Mean Absolute Relative Bias by Sample Size (Baseline)- Exp.4-B25",  
    x = "Phase 2 Sample Size", 
    y = "Mean Absolute Relative Bias",  
    shape = "Analysis Method",  
    color = "Analysis Method",  
    linetype = "Data Completeness"
  ) +
  scale_x_continuous(
    breaks = c(500, 1500, 2500, 3500, 4500)
  ) +
  theme(
    axis.text.x = element_text(angle = 80, hjust = 1), legend.position = "bottom",
    legend.box = "vertical"
  )


######################################################
# Web Figure 18

p1.1 <- stats_by_group %>%
  filter(G_vec == "G1", Experiment == 4, Sampling_Method %in% c("Best-SRS",
                                                                "GA-Z-Bal",
                                                                "GA-Y-Bal",
                                                                "GA-YZ-Bal"), BETA1 == 25) %>%
  group_by(Sampling_Method, Analysis_Method, Sample_Size) %>%
  ggplot(aes(x = Sample_Size, y = Mean_Absolute_Relative_Bias, shape = Analysis_Method, color = Analysis_Method, linetype = "Incomplete Data")) +
  geom_line(linewidth = 1.2) +
  geom_line(aes(y = Mean_Absolute_Relative_Bias_Complete, linetype = "Complete Data"), linewidth = 0.5, alpha = 0.4) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_bw() +
  facet_grid(Model_Type ~ Sampling_Method, labeller = labeller(Model_Type = model_type_labels)) + 
  labs(
    title = "Mean Absolute Relative Bias by Sample Size (Optimal)- Exp.4-B25",  
    x = "Phase 2 Sample Size", 
    y = "Mean Absolute Relative Bias",  
    shape = "Analysis Method",  
    color = "Analysis Method",  
    linetype = "Data Completeness"
  ) +
  scale_x_continuous(
    breaks = c(500, 1500, 2500, 3500, 4500)
  ) +
  theme(
    axis.text.x = element_text(angle = 80, hjust = 1), legend.position = "bottom",
    legend.box = "vertical"
  )


# Web Figure 19

p1.3 <- stats_by_group %>%
  filter(G_vec == "G1", Experiment == 4, Sampling_Method %in% c("Best-SRS",
                                                                "GA-Z-Bal",
                                                                "GA-Y-Bal",
                                                                "GA-YZ-Bal"),
         Analysis_Method == "EM", BETA1 == 25) %>%
  group_by(Sampling_Method, Sample_Size) %>%
  ggplot(aes(x = Sample_Size, y = Mean_Absolute_Relative_Bias, linetype = "Incomplete Data")) +
  geom_line(linewidth = 1.2) +
  geom_line(aes(y = Mean_Absolute_Relative_Bias_Complete, linetype = "Complete Data"), linewidth = 0.5, alpha = 0.4) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_bw() +
  facet_grid(Model_Type ~ Sampling_Method, labeller = labeller(Model_Type = model_type_labels)) + 
  labs(
    title = "Mean Absolute Relative Bias by Sample Size (Optimal-EM only)- Exp.4-B25",  
    x = "Phase 2 Sample Size", 
    y = "Mean Absolute Relative Bias",  
    shape = "Analysis Method",  
    color = "Analysis Method",  
    linetype = "Data Completeness"
  ) +
  scale_x_continuous(
    breaks = c(500, 1500, 2500, 3500, 4500)
  ) +
  theme(
    axis.text.x = element_text(angle = 80, hjust = 1), legend.position = "bottom",
    legend.box = "vertical"
  )





###########################################################
# Web Figure 20

p3 <- stats_by_group %>%
  filter(G_vec == "G1", Experiment == 4, Sampling_Method != "Complete_Data", BETA1 == 25) %>%
  group_by(Sampling_Method, Analysis_Method, Sample_Size) %>%
  ggplot(aes(x = Sample_Size, y = Relative_Efficiency, shape = Analysis_Method, color = Analysis_Method, linetype = Analysis_Method)) +
  geom_line(linewidth=1) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  theme_bw ()  +
  facet_grid(Model_Type ~ Sampling_Method, labeller = labeller(Model_Type = model_type_labels))  + 
  labs(
    title = "MSE Ratio - Exp.4-B25",  
    x = "Phase 2 Sample Size", 
    y = "MSE Ratio",  
    shape = "Analysis Method",  
    color = "Analysis Method",  
    linetype = "Analysis Method"  
  )  +
  scale_x_continuous(
    breaks = c(500,  1500,  2500,  3500,  4500)  ) +
  theme(
    axis.text.x = element_text(angle = 80, hjust = 1), legend.position = "bottom",
    legend.box = "vertical")



################################################################################
################################################################################


############# Power & Type 1 Error #############

library(dplyr)
library(purrr)
library(stats)

calculate_power_type1_error <- function(data, alpha, statistic) {
  critical_value <- qchisq(1 - alpha, df = 1)
  
  statistic <- sym(statistic)
  
 
  
  clopper_pearson_ci <- function(successes, n, alpha = 0.05) {
    lower <- qbeta(alpha / 2, successes, n - successes + 1)
    upper <- qbeta(1 - alpha / 2, successes + 1, n - successes)
    return(c(lower = lower, upper = upper))
  }
  
  
  type1_errors <- data %>%
    filter(G_vec == "G0", Sample_Size %in% c(500,1000,1500,2000,2500)) %>%
    group_by(Sample_Size, Experiment, Sampling_Strategy, Method, Model_Type, BETA1) %>%
    summarise(
      Trials = n(),
      Successes = sum(!!statistic > critical_value, na.rm = TRUE),
      Type1_Error_Rate = Successes / Trials,
      .groups = 'drop'
    ) %>%
    mutate(
      CI = map2(Successes, Trials, clopper_pearson_ci),
      CI_Lower = map_dbl(CI, "lower"),
      CI_Upper = map_dbl(CI, "upper")
    ) %>%
    dplyr::select(-CI)
  
  power <- data %>%
    filter(G_vec == "G1", Sample_Size %in% c(500,1000,1500,2000,2500)) %>%
    group_by(Sample_Size, Experiment, Sampling_Strategy, Method, Model_Type, BETA1) %>%
    summarise(
      Trials = n(),
      Successes = sum(!!statistic > critical_value, na.rm = TRUE),
      Power = Successes / Trials,
      .groups = 'drop'
    ) %>%
    mutate(
      CI = map2(Successes, Trials, clopper_pearson_ci),
      CI_Lower = map_dbl(CI, "lower"),
      CI_Upper = map_dbl(CI, "upper")
    ) %>%
    dplyr::select(-CI)
  
  coverage <- data %>%
    filter(Sample_Size %in% c(500,1000,1500,2000,2500)) %>%
    group_by(Sample_Size, Experiment, Sampling_Strategy, Method, Model_Type, BETA1) %>%
    summarise(
      Trials = n(),
      Successes = sum((True_Beta1 >= beta1 - 1.96 * sqrt(var_beta1)) & (True_Beta1 <= beta1 + 1.96 * sqrt(var_beta1)), na.rm = TRUE),
      Coverage = Successes / Trials,
      .groups = 'drop'
    ) %>%
    mutate(
      CI = map2(Successes, Trials, clopper_pearson_ci),
      Coverage_CI_Lower = map_dbl(CI, "lower"),
      Coverage_CI_Upper = map_dbl(CI, "upper")
    ) %>%
    dplyr::select(-CI)
  
  # Combine Type 1 Error, Power, and Coverage into a single result
  combined_results <- full_join(type1_errors, power, by = c("Sample_Size", "Sampling_Strategy", "Experiment", "Method", "Model_Type", "BETA1")) %>%
    full_join(coverage, by = c("Sample_Size", "Sampling_Strategy", "Experiment", "Method", "Model_Type", "BETA1"))
  
  return(combined_results)
}

power_type1_error_results_0.05_LR <- calculate_power_type1_error(all_combined, 0.05, "LR")
power_type1_error_results_0.05_LR$Test_Statistic <- "LR"

power_type1_error_results_0.05_W <- calculate_power_type1_error(all_combined, 0.05, "W")
power_type1_error_results_0.05_W$Test_Statistic <- "Wald"

power_type1_error_results_0.05_S <- calculate_power_type1_error(all_combined, 0.05, "S")
power_type1_error_results_0.05_S$Test_Statistic <- "Score"

power_type1_error_results_0.05_all_tests <- rbind(power_type1_error_results_0.05_LR,power_type1_error_results_0.05_W,power_type1_error_results_0.05_S)


power_type1_error_results_0.05_all_tests <- power_type1_error_results_0.05_all_tests %>%
  mutate(
    Sampling_Method = case_when(
      Method == "Balanced-EM" ~ "Y-Bal",
      Method == "Balanced-Naive" ~ "Y-Bal",
      Method == "SRS-EM" ~ "SRS",
      Method == "SRS-Naive" ~ "SRS",
      Method == "YZ-Balanced-EM" ~ "YZ-Bal",
      Method == "YZ-Balanced-Naive" ~ "YZ-Bal",
      Method == "Z-Balanced-EM" ~ "Z-Bal",
      Method == "Z-Balanced-Naive" ~ "Z-Bal",
      Method == "Fit-SRS-Naive" ~ "Best-SRS",
      Method == "Fit_SRS-EM" ~ "Best-SRS",
      Method == "GA-Y-EM" ~ "GA-Y-Bal",
      Method == "GA-Y-Naive" ~ "GA-Y-Bal",
      Method == "GA-YZ-EM" ~ "GA-YZ-Bal",
      Method == "GA-YZ-Naive" ~ "GA-YZ-Bal",
      Method == "GA-Z-EM" ~ "GA-Z-Bal",
      Method == "GA-Z-Naive" ~ "GA-Z-Bal",
      Method == "Complete" ~ "Complete_Data",
      TRUE ~ NA_character_
    ),
    Analysis_Method = case_when(
      str_detect(Method, "EM") ~ "EM",
      str_detect(Method, "Naive") ~ "Naive",
      Method == "Complete" ~ "Complete",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    Sampling_Method = factor(Sampling_Method, levels = c("SRS", "Best-SRS",
                                                         "Z-Bal", "GA-Z-Bal",
                                                         "Y-Bal", "GA-Y-Bal",
                                                         "YZ-Bal", "GA-YZ-Bal", "Complete_Data"))
  )


# Web Figure 21
Power_Effect_E1_n500 <- power_type1_error_results_0.05_all_tests %>%
  filter(Test_Statistic == "LR", Experiment == 4, Sample_Size == 500, Analysis_Method != "Complete", Sampling_Method != "Complete") %>%
  group_by(Sampling_Method, Analysis_Method) %>%
  ggplot(aes(x = BETA1, y = Power, color = Analysis_Method)) +
  geom_line() + 
  geom_ribbon(aes(ymin = CI_Lower.y, ymax = CI_Upper.y, fill = Analysis_Method), alpha = 0.3, color = NA) +
  facet_grid(Model_Type ~ Sampling_Method, labeller = labeller(Model_Type = model_type_labels)) +
  labs(
    title = "Power vs Effect Size - Exp.4 (n=500)",
    x = "Effect Size",
    y = "Power",
    color = "Analysis Method",
    fill = "Analysis Method"
  )  +
  scale_x_continuous(
    breaks = c(15, 20, 25, 30, 35),
    labels = c("B15", "B20", "B25", "B30", "B35")
  ) +
  ylim(0.15, 1.02) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 80, hjust = 1), legend.position = "bottom",
    legend.box = "vertical")




