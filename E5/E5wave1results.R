load("w1analysis_dat.RData")



library(ggplot2)

ac_jp <- ggplot(adcat_table_bias_df, aes(x=w1samplesize, y=SumMeanAbsBias, color=Method)) +
  geom_line(size=1.2) +
  geom_point(size=2) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(
    title="Adjacent-Category: Sum Mean\nAbsolute Bias for Joint Probabilities",
    x="Wave 1 Sample Size",
    y="Sum of Mean Absolute Bias"
  ) 

ac_bes <- ggplot(adcat_g_bias_df, aes(x=w1samplesize, y=MeanAbsBias, color=Method)) +
  geom_line(size=1.2) +
  geom_point(size=2) +
  
  geom_hline(yintercept=0, linetype="dashed") +
  labs(
    title="Adjacent-Category: Mean\nAbsolute Relative Bias for Beta_1 Estimates",
    x="Wave 1 Sample Size",
    y="Mean Absolute Bias"
  ) 
















###############################
###############################


cl_jp <- ggplot(propodds_table_bias_df, aes(x=w1samplesize, y=SumMeanAbsBias, color=Method)) +
  geom_line(size=1.2) +
  geom_point(size=2) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(
    title="Cumulative-Logit: Sum Mean\nAbsolute Bias for Joint Probabilities",
    x="Wave 1 Sample Size",
    y="Sum of Mean Absolute Bias"
  ) 

cl_bes <- ggplot(propodds_g_bias_df, aes(x=w1samplesize, y=MeanAbsBias, color=Method)) +
  geom_line(size=1.2) +
  geom_point(size=2) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(
    title="Cumulative-Logit: Mean Absolute\n Relative Bias for Beta_1 Estimates",
    x="Wave 1 Sample Size",
    y="Mean Absolute Bias"
  ) 














###############################################
###############################################














sr_jp <- ggplot(stoprat_table_bias_df, aes(x=w1samplesize, y=SumMeanAbsBias, color=Method)) +
  geom_line(size=1.2) +
  geom_point(size=2) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(
    title="Stopping-Ratio: Sum Mean\nAbsolute Bias for Joint Probabilities",
    x="Wave 1 Sample Size",
    y="Sum of Mean Absolute Bias"
  ) 

sr_bes <- ggplot(stoprat_g_bias_df, aes(x=w1samplesize, y=MeanAbsBias, color=Method)) +
  geom_line(size=1.2) +
  geom_point(size=2) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(
    title="Stopping-Ratio: Mean Absolute\n Relative Bias for Beta_1 Estimates",
    x="Wave 1 Sample Size",
    y="Mean Absolute Bias"
  ) 







################################################################################


table_bias_df <- rbind(adcat_table_bias_df, propodds_table_bias_df, stoprat_table_bias_df)

g_bias_df <- rbind(adcat_g_bias_df, propodds_g_bias_df, stoprat_g_bias_df)

library(dplyr)


g_bias_df <- g_bias_df %>%
  mutate(Method=recode(Method,
                         "YBal"="Y-Bal",
                         "YZBal"="YZ-Bal",
                         "ZBal"="Z-Bal"))

table_bias_df <- table_bias_df %>%
  mutate(Method=recode(Method,
                         "YBal"="Y-Bal",
                         "YZBal"="YZ-Bal",
                         "ZBal"="Z-Bal"))

################################################################################




table_plot <- ggplot(table_bias_df, aes(x=w1samplesize, y=SumMeanAbsBias, color=Method)) +
  geom_line(aes(linetype=Method), size=0.8) +
  #geom_point(size=1.2) +
  facet_grid(. ~ Model) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(
    title="A) Mean Absolute Bias for Joint Probabilities",
    x="Wave 1 Sample Size",
    y="Mean Absolute Bias"
  ) + 
  theme_bw() +
  theme(strip.text=element_text(size=12, face="bold"),
        legend.position="none")




g_bias_plot <- ggplot(g_bias_df, aes(x=w1samplesize, y=MeanAbsBias, color=Method)) +
  geom_line(aes(linetype=Method), size=0.8) +
  #geom_point(size=1.2) +
  facet_grid(. ~ Model) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(
    title=bquote("B) Mean Absolute Relative Bias for " ~ beta[1]),
    x="Wave 1 Sample Size",
    y="Mean Absolute Bias"
  ) +
  theme_bw() +
  theme(
    strip.text=element_text(size=12, face="bold"),
    legend.position="bottom",              
    legend.title=element_text(size=11),
    legend.text=element_text(size=10),
    legend.box="horizontal"                
  )

library(gridExtra)

# Paper's Figure 1
grid.arrange(table_plot, g_bias_plot, ncol=1)