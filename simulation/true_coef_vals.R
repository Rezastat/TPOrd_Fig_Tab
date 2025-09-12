library(TPOrd)
library(ggplot2)
library(reshape2)
library(dplyr)
library(parallel)



array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset=1))


m <- as.numeric(Sys.getenv("ITERATIONS", "10")) # number of iterations
num_cores=as.numeric(Sys.getenv("NUM_CORES", "5")) # number of cores we want to use




target_cor <- 0.75
fractions_cat1 <- c(0.35,0.25,0.2)
fractions_cat2 <- c(0.30,0.25,0.25,0.1)
N <- 5000
cor_YZ_break <- 0.85
num_categories <- 5
n2 <- 1000
model_type="Stopping_Ratio"
Beta0 <- c(0.6,0.2,-0.7,0.4)
Beta1 <- -0.25





tru_simulation <- function(iteration) {
  my_seed <-  iteration  + 1000 * array_id
  set.seed(my_seed)
  
  dat_sim <- sim_categorical_data (Beta0=Beta0 ,# intercept
                                   Beta1=Beta1, # effect size 
                                   fractions_cat1 =fractions_cat1,
                                   fractions_cat2 =fractions_cat2,
                                   target_cor=target_cor,
                                   N=N, # phase 1  sample size
                                   n2=n2, # phase 2 sample size
                                   cor_YZ_break=cor_YZ_break,
                                   num_categories=num_categories,
                                   model_type=model_type) 
  
  mod_propodds <- vglm(Y ~ G1, family=propodds(reverse=FALSE), data=dat_sim)
  mod_acat <- vglm(Y ~ G1, family=acat(reverse=TRUE, parallel=TRUE), data=dat_sim)
  mod_sratio <- vglm(Y ~ G1, family=sratio(reverse=FALSE, parallel=TRUE), data=dat_sim)

  # get coefficients
  coefs_propodds <- coef(mod_propodds)
  coefs_acat <- coef(mod_acat)
  coefs_sratio <- coef(mod_sratio)

  
  # probability matrices
  prop_odds_pmat <- calc_pvals(coef(mod_propodds)[1:(num_categories-1)], coef(mod_propodds)[num_categories], "Proportional_Odds", num_categories, 1:(length(fractions_cat1) + 1))
  adjac_cat_pmat <- calc_pvals(coef(mod_acat)[1:(num_categories-1)], coef(mod_acat)[num_categories], "Adjacent_Category", num_categories, 1:(length(fractions_cat1) + 1))
  stop_ratio_pmat <- calc_pvals(coef(mod_sratio)[1:(num_categories-1)], coef(mod_sratio)[num_categories], "Stopping_Ratio", num_categories, 1:(length(fractions_cat1) + 1))

  return(list(seed=my_seed, coefs_propodds=coefs_propodds, coefs_acat=coefs_acat, coefs_sratio=coefs_sratio, prop_odds_pmat=prop_odds_pmat, adjac_cat_pmat=adjac_cat_pmat, stop_ratio_pmat=stop_ratio_pmat)) #, stereotype_pmat=stereotype_pmat))
}

cl <- makeCluster(num_cores)

clusterExport(cl, varlist=c(
  "array_id", "m", "num_cores",
  "target_cor", "fractions_cat1", "fractions_cat2",
  "N", "cor_YZ_break", "num_categories",
  "n2", "model_type", "Beta0", "Beta1", "tru_simulation"
  
))
clusterEvalQ(cl, {
  library(TPOrd)
})

total.time <- system.time({results <- parLapply(cl, 1:m, tru_simulation)})[3]

stopCluster(cl)


############

#  coefficients for each model type
coefs_propodds <- do.call(rbind, lapply(results, `[[`, "coefs_propodds"))
coefs_acat <- do.call(rbind, lapply(results, `[[`, "coefs_acat"))
coefs_sratio <- do.call(rbind, lapply(results, `[[`, "coefs_sratio"))


seeds <- sapply(results, `[[`, "seed")

#  probability matrices for each model type
prop_odds_pmat_combined <- do.call(rbind, lapply(results, `[[`, "prop_odds_pmat"))
adjac_cat_pmat_combined <- do.call(rbind, lapply(results, `[[`, "adjac_cat_pmat"))
stop_ratio_pmat_combined <- do.call(rbind, lapply(results, `[[`, "stop_ratio_pmat"))

all_pmat <- rbind(prop_odds_pmat_combined,adjac_cat_pmat_combined,stop_ratio_pmat_combined) # ,stereotype_pmat_combined)
model_type_vector <- c(rep("Proportional_Odds", m*(length(fractions_cat1) + 1)), rep("Adjacent_Category", m*(length(fractions_cat1) + 1)), rep("Stopping_Ratio", m*(length(fractions_cat1) + 1)))
G_vec <- rep(1:(length(fractions_cat1) + 1), m*3)
iteration_vec <- rep(rep(1:m,each=(length(fractions_cat1) + 1)), 3)
pmat_df <- data.frame(iteration =iteration_vec, model_type=model_type_vector, G=G_vec, p1=all_pmat[,1], p2=all_pmat[,2], p3=all_pmat[,3],p4=all_pmat[,4], p5=all_pmat[,5])


save(seeds, total.time, coefs_propodds, coefs_acat, coefs_sratio, pmat_df,
     file=paste0("model_stoprat_800_", array_id, ".RData"))