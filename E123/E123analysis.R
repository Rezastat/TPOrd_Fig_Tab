library(dplyr)
library(ggplot2)
library(tidyr)
library(rlang)
library(purrr)
library(gridExtra)
library(tinytex)
library(stringr)


# the dataset can be provided by the corresponding author upon reasonable request
load("stats_matrix_E123.RData")


E1_CL <- stats_by_group %>%
  filter(G_vec=="G1",Experiment==1,Sampling_Method != "Complete_Data",BETA1==25,Model_Type=="propodds",Sample_Size %in% c(200,500,1500)) %>%
  select(Sample_Size,Analysis_Method,Sampling_Method,Mean_Absolute_Relative_Bias)


library(dplyr)
library(tidyr)

# row order
methods <- c("SRS","Best-SRS","Y-Bal","GA-Y-Bal",
             "Z-Bal","GA-Z-Bal","YZ-Bal","GA-YZ-Bal")
analyses <- c("Naive","EM")

# Table 1 of the paper
eff <- E1_CL %>%
  filter(as.character(Sampling_Method) %in% methods) %>%
  mutate(
    Sampling_Method=factor(as.character(Sampling_Method),levels=methods),
    Analysis_Method=factor(Analysis_Method,levels=analyses)
  ) %>%
  arrange(Sampling_Method,Analysis_Method) %>%
  pivot_wider(
    id_cols  =c(Sampling_Method,Analysis_Method),
    names_from=Sample_Size,
    values_from=Mean_Absolute_Relative_Bias
  )

eff




model_type_labels <- c(
  "adcat"="AC",
  "propodds"="CL",
  "stoprat"="SR"
)
stats_by_group$Sampling_Method <- factor(stats_by_group$Sampling_Method,
                                         levels=c("SRS","Best-SRS",
                                                    "Z-Bal","GA-Z-Bal",
                                                    "Y-Bal","GA-Y-Bal",
                                                    "YZ-Bal","GA-YZ-Bal"))
# Web Figure 2
p1 <- stats_by_group %>%
  filter(G_vec=="G1",Experiment==1,Sampling_Method != "Complete_Data",BETA1==25) %>%
  group_by(Sampling_Method,Analysis_Method,Sample_Size) %>%
  ggplot(aes(x=Sample_Size,y=Mean_Absolute_Relative_Bias,shape=Analysis_Method,color=Analysis_Method,linetype="Incomplete Data")) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=Mean_Absolute_Relative_Bias_Complete,linetype="Complete Data"),linewidth=0.5,alpha=0.4) +
  geom_hline(yintercept=0,color="black",linetype="dashed") +
  theme_grey() +
  facet_grid(Model_Type ~ Sampling_Method,labeller=labeller(Model_Type=model_type_labels)) + 
  labs(
    title="Mean Absolute Relative Bias by Sample Size - Exp.1-B25", 
    x="Phase 2 Sample Size",
    y="Mean Absolute Relative Bias", 
    shape="Analysis Method", 
    color="Analysis Method", 
    linetype="Data Completeness"
  ) +
  theme_bw() +
  scale_x_continuous(
    breaks=c(500,1500,2500,3500,4500)
  ) +
  theme(
    axis.text.x=element_text(angle=80,hjust=1),legend.position="bottom",
    legend.box="vertical"
  ) 





###########################################################
# Web Figure 3
p3 <- stats_by_group %>%
  filter(G_vec=="G1",Experiment==1,Sampling_Method != "Complete_Data",BETA1==25) %>%
  group_by(Sampling_Method,Analysis_Method,Sample_Size) %>%
  ggplot(aes(x=Sample_Size,y=Relative_Efficiency,shape=Analysis_Method,color=Analysis_Method,linetype=Analysis_Method)) +
  geom_line(linewidth=1) +
  geom_hline(yintercept=1,color="black",linetype="dashed") +
  theme_grey ()  +
  facet_grid(Model_Type ~ Sampling_Method,labeller=labeller(Model_Type=model_type_labels))  
  labs(
    title="MSE Ratio - Exp.1-B25", 
    x="Phase 2 Sample Size",
    y="MSE Ratio", 
    shape="Analysis Method", 
    color="Analysis Method", 
    linetype="Analysis Method"  
  )  +
  theme_bw() +
  scale_x_continuous(
    breaks=c(500, 1500, 2500, 3500, 4500)  ) +
  theme(
    axis.text.x=element_text(angle=80,hjust=1),legend.position="bottom",
    legend.box="vertical")


###########################################################
# Web Figure 6
p4 <- stats_by_group %>%
  filter(G_vec=="G1",Experiment==2,Sampling_Method %in% c("Y-Bal","SRS","YZ-Bal","Z-Bal"),BETA1==25) %>%
  group_by(Sampling_Method,Analysis_Method,Sample_Size) %>%
  ggplot(aes(x=Sample_Size,y=Mean_Absolute_Relative_Bias,shape=Analysis_Method,color=Analysis_Method,linetype="Incomplete Data")) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=Mean_Absolute_Relative_Bias_Complete,linetype="Complete Data"),linewidth=0.5,alpha=0.4) +
  geom_hline(yintercept=0,color="black",linetype="dashed") +
  theme_grey() +
  facet_grid(Model_Type ~ Sampling_Method,labeller=labeller(Model_Type=model_type_labels)) + 
  labs(
    title="Mean Absolute Relative Bias by Sample Size - Exp.2-B25", 
    x="Phase 2 Sample Size",
    y="Mean Absolute Relative Bias", 
    shape="Analysis Method", 
    color="Analysis Method", 
    linetype="Data Completeness"
  ) +
  theme_bw() +
  scale_x_continuous(
    breaks=c(500,1500,2500,3500,4500)
  ) +
  theme(
    axis.text.x=element_text(angle=80,hjust=1),legend.position="bottom",
    legend.box="vertical"
  )


###########################################################
# Web Figure 7
p6 <- stats_by_group %>%
  filter(G_vec=="G1",Experiment==2,Sampling_Method %in% c("Y-Bal","SRS","YZ-Bal","Z-Bal"),BETA1==25) %>%
  group_by(Sampling_Method,Analysis_Method,Sample_Size) %>%
  ggplot(aes(x=Sample_Size,y=Relative_Efficiency,shape=Analysis_Method,color=Analysis_Method,linetype=Analysis_Method)) +
  geom_line(linewidth=1) +
  geom_hline(yintercept=1,color="black",linetype="dashed") +
  theme_grey ()  +
  facet_grid(Model_Type ~ Sampling_Method,labeller=labeller(Model_Type=model_type_labels))  +
  labs(
    title="MSE Ratio - Exp.2-B25", 
    x="Phase 2 Sample Size",
    y="MSE Ratio", 
    shape="Analysis Method", 
    color="Analysis Method", 
    linetype="Analysis Method"  
  )  +
  theme_bw() +
  scale_x_continuous(
    breaks=c(500, 1500, 2500, 3500, 4500)  ) +
  theme(
    axis.text.x=element_text(angle=80,hjust=1),legend.position="bottom",
    legend.box="vertical")

#############################################################
# Web Figure 9
p7 <- stats_by_group %>%
  filter(G_vec=="G1",Experiment==3,Sampling_Method %in% c("Y-Bal","SRS","YZ-Bal","Z-Bal"),BETA1==25) %>%
  group_by(Sampling_Method,Analysis_Method,Sample_Size) %>%
  ggplot(aes(x=Sample_Size,y=Mean_Absolute_Relative_Bias,shape=Analysis_Method,color=Analysis_Method,linetype="Incomplete Data")) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=Mean_Absolute_Relative_Bias_Complete,linetype="Complete Data"),linewidth=0.5,alpha=0.4) +
  geom_hline(yintercept=0,color="black",linetype="dashed") +
  theme_grey() +
  facet_grid(Model_Type ~ Sampling_Method,labeller=labeller(Model_Type=model_type_labels)) + 
  labs(
    title="Mean Absolute Relative Bias by Sample Size - Exp.3-B25", 
    x="Phase 2 Sample Size",
    y="Mean Absolute Relative Bias", 
    shape="Analysis Method", 
    color="Analysis Method", 
    linetype="Data Completeness"
  ) +
  theme_bw() +
  scale_x_continuous(
    breaks=c(500,1500,2500,3500,4500)
  ) +
  theme(
    axis.text.x=element_text(angle=80,hjust=1),legend.position="bottom",
    legend.box="vertical"
  )


######################################################

# Web Figure 10
p9 <- stats_by_group %>%
  filter(G_vec=="G1",Experiment==3,Sampling_Method %in% c("Y-Bal","SRS","YZ-Bal","Z-Bal"),BETA1==25) %>%
  group_by(Sampling_Method,Analysis_Method,Sample_Size) %>%
  ggplot(aes(x=Sample_Size,y=Relative_Efficiency,shape=Analysis_Method,color=Analysis_Method,linetype=Analysis_Method)) +
  geom_line(linewidth=1) +
  geom_hline(yintercept=1,color="black",linetype="dashed") +
  theme_grey ()  +
  facet_grid(Model_Type ~ Sampling_Method,labeller=labeller(Model_Type=model_type_labels))  +
  labs(
    title="MSE Ratio - Exp.3-B25", 
    x="Phase 2 Sample Size",
    y="MSE Ratio", 
    shape="Analysis Method", 
    color="Analysis Method", 
    linetype="Analysis Method"  
  )  +
  theme_bw() +
  scale_x_continuous(
    breaks=c(500, 1500, 2500, 3500, 4500)  ) +
  theme(
    axis.text.x=element_text(angle=80,hjust=1),legend.position="bottom",
    legend.box="vertical")

############################################################
############################################################
################################################################################


library(dplyr)
library(purrr)
library(stats)

calculate_power_type1_error <- function(data,alpha,statistic) {
  critical_value <- qchisq(1 - alpha,df=1)
  
  statistic <- sym(statistic)
  
  
  clopper_pearson_ci <- function(successes,n,alpha=0.05) {
    lower <- qbeta(alpha / 2,successes,n - successes + 1)
    upper <- qbeta(1 - alpha / 2,successes + 1,n - successes)
    return(c(lower=lower,upper=upper))
  }
  
  
  type1_errors <- data %>%
    filter(G_vec=="G0",Sample_Size %in% c(500,1000,1500,2000,2500)) %>%
    group_by(Sample_Size,Experiment,Sampling_Strategy,Method,Model_Type,BETA1) %>%
    summarise(
      Trials=n(),
      Successes=sum(!!statistic > critical_value,na.rm=TRUE),
      Type1_Error_Rate=Successes / Trials,
      .groups='drop'
    ) %>%
    mutate(
      CI=map2(Successes,Trials,clopper_pearson_ci),
      CI_Lower=map_dbl(CI,"lower"),
      CI_Upper=map_dbl(CI,"upper")
    ) %>%
    dplyr::select(-CI)
  
  power <- data %>%
    filter(G_vec=="G1",Sample_Size %in% c(500,1000,1500,2000,2500)) %>%
    group_by(Sample_Size,Experiment,Sampling_Strategy,Method,Model_Type,BETA1) %>%
    summarise(
      Trials=n(),
      Successes=sum(!!statistic > critical_value,na.rm=TRUE),
      Power=Successes / Trials,
      .groups='drop'
    ) %>%
    mutate(
      CI=map2(Successes,Trials,clopper_pearson_ci),
      CI_Lower=map_dbl(CI,"lower"),
      CI_Upper=map_dbl(CI,"upper")
    ) %>%
    dplyr::select(-CI)
  
  coverage <- data %>%
    filter(Sample_Size %in% c(500,1000,1500,2000,2500)) %>%
    group_by(Sample_Size,Experiment,Sampling_Strategy,Method,Model_Type,BETA1) %>%
    summarise(
      Trials=n(),
      Successes=sum((True_Beta1 >= beta1 - 1.96 * sqrt(var_beta1)) & (True_Beta1 <= beta1 + 1.96 * sqrt(var_beta1)),na.rm=TRUE),
      Coverage=Successes / Trials,
      .groups='drop'
    ) %>%
    mutate(
      CI=map2(Successes,Trials,clopper_pearson_ci),
      Coverage_CI_Lower=map_dbl(CI,"lower"),
      Coverage_CI_Upper=map_dbl(CI,"upper")
    ) %>%
    dplyr::select(-CI)
  
  # Combine Type 1 Error,Power,and Coverage into a single result
  combined_results <- full_join(type1_errors,power,by=c("Sample_Size","Sampling_Strategy","Experiment","Method","Model_Type","BETA1")) %>%
    full_join(coverage,by=c("Sample_Size","Sampling_Strategy","Experiment","Method","Model_Type","BETA1"))
  
  return(combined_results)
}

power_type1_error_results_0.05_LR <- calculate_power_type1_error(all_combined,0.05,"LR")
power_type1_error_results_0.05_LR$Test_Statistic <- "LR"

power_type1_error_results_0.05_W <- calculate_power_type1_error(all_combined,0.05,"W")
power_type1_error_results_0.05_W$Test_Statistic <- "Wald"

power_type1_error_results_0.05_S <- calculate_power_type1_error(all_combined,0.05,"S")
power_type1_error_results_0.05_S$Test_Statistic <- "Score"

power_type1_error_results_0.05_all_tests <- rbind(power_type1_error_results_0.05_LR,power_type1_error_results_0.05_W,power_type1_error_results_0.05_S)


power_type1_error_results_0.05_all_tests <- power_type1_error_results_0.05_all_tests %>%
  mutate(
    Sampling_Method=case_when(
      Method=="Balanced-EM" ~ "Y-Bal",
      Method=="Balanced-Naive" ~ "Y-Bal",
      Method=="SRS-EM" ~ "SRS",
      Method=="SRS-Naive" ~ "SRS",
      Method=="YZ-Balanced-EM" ~ "YZ-Bal",
      Method=="YZ-Balanced-Naive" ~ "YZ-Bal",
      Method=="Z-Balanced-EM" ~ "Z-Bal",
      Method=="Z-Balanced-Naive" ~ "Z-Bal",
      Method=="Fit-SRS-Naive" ~ "Best-SRS",
      Method=="Fit_SRS-EM" ~ "Best-SRS",
      Method=="GA-Y-EM" ~ "GA-Y-Bal",
      Method=="GA-Y-Naive" ~ "GA-Y-Bal",
      Method=="GA-YZ-EM" ~ "GA-YZ-Bal",
      Method=="GA-YZ-Naive" ~ "GA-YZ-Bal",
      Method=="GA-Z-EM" ~ "GA-Z-Bal",
      Method=="GA-Z-Naive" ~ "GA-Z-Bal",
      Method=="Complete" ~ "Complete_Data",
      TRUE ~ NA_character_
    ),
    Analysis_Method=case_when(
      str_detect(Method,"EM") ~ "EM",
      str_detect(Method,"Naive") ~ "Naive",
      Method=="Complete" ~ "Complete",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    Sampling_Method=factor(Sampling_Method,levels=c("SRS","Z-Bal","Y-Bal","YZ-Bal","Best-SRS","GA-Z-Bal","GA-Y-Bal","GA-YZ-Bal","Complete_Data"))
  )


################################################################################


the_order <- c("SRS","Best-SRS","Z-Bal","GA-Z-Bal","Y-Bal","GA-Y-Bal","YZ-Bal","GA-YZ-Bal")

power_type1_error_results_0.05_all_tests <- power_type1_error_results_0.05_all_tests %>%
  mutate(Sampling_Method=factor(Sampling_Method,levels=the_order))

# Web Figure 4
Power_Effect_E1_n500 <- power_type1_error_results_0.05_all_tests %>%
  filter(Test_Statistic=="LR",Experiment==1,Sample_Size==500,Analysis_Method != "Complete") %>%
  group_by(Sampling_Method,Analysis_Method) %>%
  ggplot(aes(x=BETA1,y=Power,color=Analysis_Method)) +
  geom_line() + 
  geom_ribbon(aes(ymin=CI_Lower.y,ymax=CI_Upper.y,fill=Analysis_Method),alpha=0.3,color=NA) +
  facet_grid(Model_Type ~ Sampling_Method,labeller=labeller(Model_Type=model_type_labels)) +
  labs(
    title="Power vs Effect Size - Exp.1 (n=500)",
    x="Effect Size",
    y="Power",
    color="Analysis Method",
    fill="Analysis Method"
  )  +
  theme_bw() +
  scale_x_continuous(
    breaks=c(15,20,25,30,35),
    labels=c("B15","B20","B25","B30","B35")
  ) +
  ylim(0.15,1.02) +
  theme(
    axis.text.x=element_text(angle=80,hjust=1),
    legend.position="bottom",
    legend.box="vertical"
  )


# Web Figure 8
Power_Effect_E2_n500 <-power_type1_error_results_0.05_all_tests %>%
  filter(Test_Statistic=="LR",Experiment==2,Sample_Size==500,Analysis_Method != "Complete") %>%
  group_by(Sampling_Method,Analysis_Method) %>%
  ggplot(aes(x=BETA1,y=Power,color=Analysis_Method)) +
  geom_line() + 
  geom_ribbon(aes(ymin=CI_Lower.y,ymax=CI_Upper.y,fill=Analysis_Method),alpha=0.3,color=NA) +
  facet_grid(Model_Type ~ Sampling_Method,labeller=labeller(Model_Type=model_type_labels)) +
  labs(
    title="Power vs Effect Size - Exp.2 (n=500)",
    x="Effect Size",
    y="Power",
    color="Analysis Method",
    fill="Analysis Method"
  )  +
  theme_bw() +
  scale_x_continuous(
    breaks=c(15,20,25,30,35),
    labels=c("B15","B20","B25","B30","B35")
  ) +
  ylim(0.15,1.02) +
  theme(
    axis.text.x=element_text(angle=80,hjust=1),legend.position="bottom",
    legend.box="vertical")


# Web Figure 11
Power_Effect_E3_n500 <-power_type1_error_results_0.05_all_tests %>%
  filter(Test_Statistic=="LR",Experiment==3,Sample_Size==500,Analysis_Method != "Complete",Sampling_Method %in% c("SRS","Y-Bal","YZ-Bal","Z-Bal")) %>%
  group_by(Sampling_Method,Analysis_Method) %>%
  ggplot(aes(x=BETA1,y=Power,color=Analysis_Method)) +
  geom_line() + 
  geom_ribbon(aes(ymin=CI_Lower.y,ymax=CI_Upper.y,fill=Analysis_Method),alpha=0.3,color=NA) +
  facet_grid(Model_Type ~ Sampling_Method,labeller=labeller(Model_Type=model_type_labels)) +
  labs(
    title="Power vs Effect Size - Exp.3 (n=500)",
    x="Effect Size",
    y="Power",
    color="Analysis Method",
    fill="Analysis Method"
  ) +
  theme_bw() +  
  scale_x_continuous(
    breaks=c(15,20,25,30,35),
    labels=c("B15","B20","B25","B30","B35")
  ) +
  ylim(0.15,1.02) +
  theme(
    axis.text.x=element_text(angle=80,hjust=1),legend.position="bottom",
    legend.box="vertical")


################################################################################
################################################################################




calculate_type1_error <- function(data,alpha,statistic) {
  critical_value <- qchisq(1 - alpha,df=1)
  statistic <- sym(statistic)
  
  clopper_pearson_ci <- function(successes,n,alpha=0.05) {
    lower <- qbeta(alpha / 2,successes,n - successes + 1)
    upper <- qbeta(1 - alpha / 2,successes + 1,n - successes)
    return(c(lower=lower,upper=upper))
  }
  
  type1_errors <- data %>%
    filter(G_vec=="G0",Sample_Size %in% c(500,1000,1500,2000,2500)) %>%
    # Pool across effect sizes by NOT grouping by BETA1.
    group_by(Sample_Size,Experiment,Sampling_Strategy,Method,Model_Type) %>%
    summarise(
      Trials=n(),
      Successes=sum(!!statistic > critical_value,na.rm=TRUE),
      Type1_Error_Rate=Successes / Trials,
      .groups='drop'
    ) %>%
    mutate(
      CI=map2(Successes,Trials,clopper_pearson_ci),
      CI_Lower=map_dbl(CI,"lower"),
      CI_Upper=map_dbl(CI,"upper")
    ) %>%
    select(-CI)
  
  return(type1_errors)
}




# For G0: Calculate the pooled type I error rate.
type1_error_results_LR <- calculate_type1_error(all_combined,0.05,"LR")
type1_error_results_LR$Test_Statistic <- "LR" 

type1_error_results_score <- calculate_type1_error(all_combined,0.05,"S")
type1_error_results_score$Test_Statistic <- "Score"  

type1_error_results_Wald <- calculate_type1_error(all_combined,0.05,"W")
type1_error_results_Wald$Test_Statistic <- "Wald"  

type1_error_results <- rbind(type1_error_results_LR,type1_error_results_score,type1_error_results_Wald)




type1_error_results <- type1_error_results %>%
  mutate(
    Sampling_Method=case_when(
      Method=="Balanced-EM" ~ "Y-Bal",
      Method=="Balanced-Naive" ~ "Y-Bal",
      Method=="SRS-EM" ~ "SRS",
      Method=="SRS-Naive" ~ "SRS",
      Method=="YZ-Balanced-EM" ~ "YZ-Bal",
      Method=="YZ-Balanced-Naive" ~ "YZ-Bal",
      Method=="Z-Balanced-EM" ~ "Z-Bal",
      Method=="Z-Balanced-Naive" ~ "Z-Bal",
      Method=="Fit-SRS-Naive" ~ "Best-SRS",
      Method=="Fit_SRS-EM" ~ "Best-SRS",
      Method=="GA-Y-EM" ~ "GA-Y-Bal",
      Method=="GA-Y-Naive" ~ "GA-Y-Bal",
      Method=="GA-YZ-EM" ~ "GA-YZ-Bal",
      Method=="GA-YZ-Naive" ~ "GA-YZ-Bal",
      Method=="GA-Z-EM" ~ "GA-Z-Bal",
      Method=="GA-Z-Naive" ~ "GA-Z-Bal",
      Method=="Complete" ~ "Complete_Data",
      TRUE ~ NA_character_
    ),
    Analysis_Method=case_when(
      str_detect(Method,"EM") ~ "EM",
      str_detect(Method,"Naive") ~ "Naive",
      Method=="Complete" ~ "Complete",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    Sampling_Method=factor(Sampling_Method,levels=c("SRS","Best-SRS","Z-Bal","GA-Z-Bal","Y-Bal","GA-Y-Bal","YZ-Bal" ,"GA-YZ-Bal","Complete_Data"))
  )



library(ggplot2)
library(dplyr)


type1_error_results <- type1_error_results %>%
  mutate(Model_Type=recode(Model_Type,
                             "propodds"="CL",
                             "adcat"="AC",
                             "stoprat"="SR"))

plot_data_baseline <- type1_error_results %>%
  filter(Analysis_Method != "Complete",
         !Sampling_Method %in% c("Best-SRS","GA-Z-Bal","GA-Y-Bal","GA-YZ-Bal"))

plot_data_baseline_e1 <- type1_error_results %>%
  filter(Analysis_Method != "Complete",
         !Sampling_Method %in% c("Best-SRS","GA-Z-Bal","GA-Y-Bal","GA-YZ-Bal"),
         Experiment==1)

plot_data_baseline_e2 <- type1_error_results %>%
  filter(Analysis_Method != "Complete",
         !Sampling_Method %in% c("Best-SRS","GA-Z-Bal","GA-Y-Bal","GA-YZ-Bal"),
         Experiment==2)

plot_data_baseline_e3 <- type1_error_results %>%
  filter(Analysis_Method != "Complete",
         !Sampling_Method %in% c("Best-SRS","GA-Z-Bal","GA-Y-Bal","GA-YZ-Bal"),
         Experiment==3)



#############################################################
# Web Figure 12
baseline_t1e_e1 <- ggplot(plot_data_baseline_e1,aes(x=Sample_Size,y=Type1_Error_Rate,color=Analysis_Method,fill=Analysis_Method)) +
  
  geom_ribbon(aes(ymin=CI_Lower,ymax=CI_Upper),alpha=0.2,color=NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept=0.05,linetype="dashed",color="black") +
  facet_grid(rows=vars(Model_Type),
             cols=vars(Sampling_Method,Test_Statistic)) +
  labs(title="Baseline Type I Error Rates (Experiment 1)",
       x="Sample Size",
       y="Type I Error Rate",
       color="Analysis Method",
       fill="Analysis Method") +
  theme_bw() +
  theme(
    axis.text.x=element_text(angle=70,hjust=1),
    legend.position="bottom",
    legend.box="horizontal",
    legend.title=element_text(size=11),
    legend.text=element_text(size=10),
    plot.title=element_text(size=14,face="bold",hjust=0.5)  
  )

#######################################################
# Web Figure 13
baseline_t1e_e2 <- ggplot(plot_data_baseline_e2,aes(x=Sample_Size,y=Type1_Error_Rate,color=Analysis_Method,fill=Analysis_Method)) +
  
  geom_ribbon(aes(ymin=CI_Lower,ymax=CI_Upper),alpha=0.2,color=NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept=0.05,linetype="dashed",color="black") +
  facet_grid(rows=vars(Model_Type),
             cols=vars(Sampling_Method,Test_Statistic)) +
  labs(title="Baseline Type I Error Rates (Experiment 2)",
       x="Sample Size",
       y="Type I Error Rate",
       color="Analysis Method",
       fill="Analysis Method") +
  theme_bw() +
  theme(
    axis.text.x=element_text(angle=70,hjust=1),
    legend.position="bottom",
    legend.box="horizontal",
    legend.title=element_text(size=11),
    legend.text=element_text(size=10),
    plot.title=element_text(size=14,face="bold",hjust=0.5)  
  )

#######################################################
# Web Figure 14
baseline_t1e_e3 <- ggplot(plot_data_baseline_e3,aes(x=Sample_Size,y=Type1_Error_Rate,color=Analysis_Method,fill=Analysis_Method)) +
  
  geom_ribbon(aes(ymin=CI_Lower,ymax=CI_Upper),alpha=0.2,color=NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept=0.05,linetype="dashed",color="black") +
  facet_grid(rows=vars(Model_Type),
             cols=vars(Sampling_Method,Test_Statistic)) +
  labs(title="Baseline Type I Error Rates (Experiment 3)",
       x="Sample Size",
       y="Type I Error Rate",
       color="Analysis Method",
       fill="Analysis Method") +
  theme_bw() +
  theme(
    axis.text.x=element_text(angle=70,hjust=1),
    legend.position="bottom",
    legend.box="horizontal",
    legend.title=element_text(size=11),
    legend.text=element_text(size=10),
    plot.title=element_text(size=14,face="bold",hjust=0.5)  
  )

#######################################################
# Web Figure 15

plot_data_optimal <- type1_error_results %>%
  filter(Analysis_Method != "Complete",Experiment==1,
         Sampling_Method %in% c("Best-SRS","GA-Z-Bal","GA-Y-Bal","GA-YZ-Bal"))

optimal_t1e <- ggplot(plot_data_optimal,aes(x=Sample_Size,y=Type1_Error_Rate,color=Analysis_Method,fill=Analysis_Method)) +
  
  geom_ribbon(aes(ymin=CI_Lower,ymax=CI_Upper),alpha=0.2,color=NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept=0.05,linetype="dashed",color="black") +
  facet_grid(rows=vars(Model_Type),
             cols=vars(Sampling_Method,Test_Statistic)) +
  labs(title="Single-Wave Optimal Sampling Type I Error Rates (Experiment 1)",
       x="Sample Size",
       y="Type I Error Rate",
       color="Analysis Method",
       fill="Analysis Method") +
  theme_bw() +
  theme(
    axis.text.x=element_text(angle=70,hjust=1),
    legend.position="bottom",
    legend.box="horizontal",
    legend.title=element_text(size=11),
    legend.text=element_text(size=10),
    plot.title=element_text(size=14,face="bold",hjust=0.5)  
  )


#############################################################
#############################################################


# Web figure 5
Coverage_Effect_E1_n500 <- power_type1_error_results_0.05_all_tests %>%
  filter(Test_Statistic=="Wald",Experiment==1,Sample_Size==500,
         Analysis_Method != "Complete") %>%
  ggplot(aes(x=factor(BETA1),y=Coverage)) +
  geom_hline(yintercept=0.95,linetype="dashed") +
  geom_point(aes(color=Analysis_Method),position=position_dodge(width=0.7),size=1.2) +  
  geom_errorbar(aes(ymin=Coverage_CI_Lower,ymax=Coverage_CI_Upper,color=Analysis_Method),
                position=position_dodge(width=0.7),width=0.2,size=1) +  
  facet_grid(Model_Type ~ Sampling_Method,labeller=labeller(Model_Type=model_type_labels)) +
  labs(
    title="Coverage vs Effect Size - Exp.1 (n=500)",
    x="Effect Size",
    y="Coverage",
    color="Analysis Method"
  ) +
  scale_x_discrete(labels=c("15"="B15","20"="B20","25"="B25","30"="B30","35"="B35")) +  
  ylim(0.8,0.98) +  
  theme_bw() +
  theme(
    axis.text.x=element_text(angle=80,hjust=1),
    legend.position="bottom",
    legend.box="vertical"
  )

Coverage_Effect_E1_n500
