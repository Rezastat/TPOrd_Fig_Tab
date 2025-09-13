##################################
## 2 wave 2 phase ##########
##################################


library(parallel)
library(TPOrd)

job_name <- Sys.getenv("JOB_NAME","Experiment5_adcat2wdata_R5678_BETA1_25")

Sys.setenv(OMP_NUM_THREADS="1",
           MKL_NUM_THREADS="1",
           OPENBLAS_NUM_THREADS="1",
           BLAS_NUM_THREADS="1",
           LAPACK_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1",
           NUMEXPR_NUM_THREADS="1")

Sys.getenv("MKL_NUM_THREADS","1")  
Sys.getenv("MKL_DYNAMIC","FALSE")      


array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID",unset=1))

m <- as.numeric(Sys.getenv("ITERATIONS","3"))
num_cores=as.numeric(Sys.getenv("NUM_CORES","3"))


N <- 5000 # Phase 1 size
cor_YZ_break <- 0.75
num_categories <- 5
n2_vector <- 1000 # phase 2 sample size. can be a vector.
model_type <-  "Adjacent_Category"
Beta0 <- c(-1.007578022,-0.510263682,2.252233149,2.27519740)
Beta1 <- -0.1238386604



joint_dist <- matrix(
  c(0.259860227,0.056747538,0.017921970,0.007525758,0.007703977,
    0.014901894,0.172798864,0.051742235,0.005044697,0.005158712,
    0.012423864,0.009817045,0.167688447,0.006222538,0.003882576,
    0.013159848,0.010351326,0.012348864,0.081296780,0.083402841),
  nrow=4,byrow=TRUE
)


run_simulation <- function(params) {
  iteration <- params$iter
  n2 <- params$n2

  fam <- switch(
    model_type,
    "Stopping_Ratio"    =sratio(reverse=FALSE,parallel=TRUE),
    "Proportional_Odds" =propodds(reverse=FALSE),
    "Adjacent_Category" =acat(reverse=TRUE,parallel=TRUE),
    stop("Invalid 'model_type' specified.")
  )

  Sys.setenv(OMP_NUM_THREADS="1",
             MKL_NUM_THREADS="1",
             OPENBLAS_NUM_THREADS="1",
             BLAS_NUM_THREADS="1",
             LAPACK_NUM_THREADS="1",
             VECLIB_MAXIMUM_THREADS="1",
             NUMEXPR_NUM_THREADS="1")


  my_seed <- iteration + 1000 * array_id
  set.seed(my_seed)






  dat_sim_w2_init <- sim_categorical_data_joint (Beta0=Beta0 ,# intercept
                                         Beta1=Beta1,# effect size of the sequencing variant
                                         joint_prob_G1_Z=joint_dist,
                                         N=N,# phase 1  sample size
                                         n2=n2,# phase 2 sample size
                                         cor_YZ_break=cor_YZ_break,
                                         num_categories=num_categories,
                                         model_type=model_type)


  n2.W1 <- floor(n2/2)
  n2.W2 <- n2 - n2.W1
  #######################################################
  #######################################################

  dat_sim_w2_init$YZ <- interaction(dat_sim_w2_init$Y,dat_sim_w2_init$Z)
  dat_sim_w2_init <- sample_dataframe(dat_sim_w2_init,n2.W1)
  dat_sim_w2_init$fZ <- factor(dat_sim_w2_init$Z)

  # R1:SRS,R2:Y-Bal,R3:Z-Bal,R4:YZ-Bal

  #######################################################
  pGZ_YBal_W1 <- table(dat_sim_w2_init[dat_sim_w2_init$R2==1,"G1"],dat_sim_w2_init[dat_sim_w2_init$R2==1,"Z"])/n2.W1


  pGZ_comp <- table(dat_sim_w2_init$G1,dat_sim_w2_init$Z)/N


  #################################################
  #################################################


  data1.R2.G1.W1 <- dat_sim_w2_init[dat_sim_w2_init$R2==1,c("Y","G1","Z","fZ")] # phase 2 data
  names(data1.R2.G1.W1)[names(data1.R2.G1.W1)=="G1"] <- "G"
  data0.R2.G1.W1 <- dat_sim_w2_init[dat_sim_w2_init$R2==0,c("Y","Z","fZ")] # phase 2 data complement
  resHa.YBal.G1.W1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R2.G1.W1,data1.R2.G1.W1,start.values=NULL,verbose=FALSE,n_second=n2.W1,model_type=model_type,num_categories=num_categories,N=N,iteration=iteration,array_id=array_id)

  fitcoef <- resHa.YBal.G1.W1$theta

  p_gz0 <- convertMeanFrequenciesToTable(pGZ_YBal_W1)


  ####################################################
  dat_sim_w2_init_w1 <- dat_sim_w2_init %>%
    filter(R2 == 1)

  dat_sim_w2_init_w1$R1 <- 1
  dat_sim_w2_init_w1$R3 <- 1
  dat_sim_w2_init_w1$R4 <- 1
  dat_sim_w2_init_w1$wave <- "w1"



  dat_sim_w2_init_not_w1 <- dat_sim_w2_init %>%
    filter(R2 != 1)

  dat_sim_w2_init_not_w1$wave <- "nw1"

  dat_sim_w2_init_w2_init <- rbind(dat_sim_w2_init_w1,dat_sim_w2_init_not_w1) %>%
    select(wave,Y,G1,Z,fZ,G0,YZ)


  dat_list <- sample_dataframe_balanced_matrix(dat_sim_w2_init_not_w1,n2.W2,1000)

  dat_list <- lapply(dat_list,function(mat) {

    rbind(matrix(1,nrow=n2.W1,ncol=ncol(mat)),mat)
  })

  ### The formula under the null is the same regardless
  Kind <- 1
  formula_Ha <- Y~G+fZ
  miscov_formula <- ~G
  formula_Ho <- Y~fZ
  auxvar=~Z
  strataformula_YZ=~YZ
  strataformula_Y=~Y
  strataformula_Z=~Z

  optMethod <- "Par-spec" # c("Par-spec","A-opt","D-opt")[1]

  dat=dat_sim_w2_init_w2_init[ ,c("Y","Z","fZ","YZ")]
  dat$YZ <- as.numeric(dat$YZ)



  ### Obtain GA design ----
  IM <- TPOrd:::obsIM_R_ord(formula=formula_Ha,
                            miscov=miscov_formula,auxvar=~Z,
                            model_type=model_type,dat,
                            beta=fitcoef,p_gz=p_gz0,num_categories=num_categories)




  indices_matrix_apply <- function(binary_matrix) {
    num_ones <- sum(binary_matrix[,1])  # Assuming all columns have the same number of `1`s
    apply(binary_matrix,2,function(col) which(col == 1))
  }

  pop1a_Y <- indices_matrix_apply(dat_list$R2)
  pop1a_YZ <- indices_matrix_apply(dat_list$R3)
  pop1a_Z <- indices_matrix_apply(dat_list$R4)


  fit1_Y <- apply(pop1a_Y,2,function(r){
    return(TPOrd:::fitnessTP_ord(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
  })

  fit1_YZ <- apply(pop1a_YZ,2,function(r){
    return(TPOrd:::fitnessTP_ord(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
  })

  fit1_Z <- apply(pop1a_Z,2,function(r){
    return(TPOrd:::fitnessTP_ord(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
  })





  pop2a <- sapply(1:1000,function(x){
    sa <- sample((n2.W1 + 1) : N,n2.W2)
  })


  pop2a <- rbind(matrix(rep(1:n2.W1,1000),nrow=n2.W1,ncol=1000),pop2a)



  fit2 <- apply(pop2a,2,function(r){
    return(TPOrd:::fitnessTP_ord(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
  })



  ### Put initial population (mix between LM and best candidates from SRS)
  pop1_Y <- pop1a_Y[,order(fit1_Y)[1:40]]
  pop1_YZ <- pop1a_YZ[,order(fit1_YZ)[1:40]]
  pop1_Z <- pop1a_Z[,order(fit1_Z)[1:40]]

  pop2 <- pop2a[,order(fit2)[1:20]]



  ### It may take a few minutes to run

  GA.sol1_Y <- TPOrd:::optimTP_GA_ord_w2(ncores=1,
                                      formula=formula_Ha,
                                      miscov=miscov_formula,
                                      auxvar=~Z,
                                      model_type=model_type,
                                      n=n2.W2,
                                      data=dat,
                                      beta=fitcoef,p_gz=p_gz0,
                                      ga.popsize=60,
                                      ga.propelit=0.5,
                                      ga.proptourney=0.9,
                                      ga.ngen=600,
                                      ga.mutrate=0.006,
                                      ga.initpop=t(cbind(pop1_Y,pop2)),
                                      optimMeasure=optMethod,K.idx=Kind,seed=1,constant_sample=n2.W1)


  GA.sol1_YZ <- TPOrd:::optimTP_GA_ord_w2(ncores=1,
                                          formula=formula_Ha,
                                          miscov=miscov_formula,
                                          auxvar=~Z,
                                          model_type=model_type,
                                          n=n2.W2,
                                          data=dat,
                                          beta=fitcoef,p_gz=p_gz0,
                                          ga.popsize=60,
                                          ga.propelit=0.5,
                                          ga.proptourney=0.9,
                                          ga.ngen=600,
                                          ga.mutrate=0.006,
                                          ga.initpop=t(cbind(pop1_YZ,pop2)),
                                          optimMeasure=optMethod,K.idx=Kind,seed=1,constant_sample=n2.W1)


  GA.sol1_Z <- TPOrd:::optimTP_GA_ord_w2(ncores=1,
                                         formula=formula_Ha,
                                         miscov=miscov_formula,
                                         auxvar=~Z,
                                         model_type=model_type,
                                         n=n2.W2,
                                         data=dat,
                                         beta=fitcoef,p_gz=p_gz0,
                                         ga.popsize=60,
                                         ga.propelit=0.5,
                                         ga.proptourney=0.9,
                                         ga.ngen=600,
                                         ga.mutrate=0.006,
                                         ga.initpop=t(cbind(pop1_Z,pop2)),
                                         optimMeasure=optMethod,K.idx=Kind,seed=1,constant_sample=n2.W1)


  ### Obtain the indicators
  dat_sim_w2_init_w2_init$R5 <- 1*(1:N %in% pop2a[,which(fit2 == min(fit2) )])
  dat_sim_w2_init_w2_init$R6 <- rep(0,N); dat_sim_w2_init_w2_init$R6[GA.sol1_Y$bestsol] <- 1
  dat_sim_w2_init_w2_init$R7 <- rep(0,N); dat_sim_w2_init_w2_init$R7[GA.sol1_YZ$bestsol] <- 1
  dat_sim_w2_init_w2_init$R8 <- rep(0,N); dat_sim_w2_init_w2_init$R8[GA.sol1_Z$bestsol] <- 1



  com_table <- table(dat_sim_w2_init$Y,dat_sim_w2_init$G1)
  R5_table <- table(dat_sim_w2_init[dat_sim_w2_init$R5==1,]$Y,dat_sim_w2_init[dat_sim_w2_init$R5==1,]$G1)
  R6_table <- table(dat_sim_w2_init[dat_sim_w2_init$R6==1,]$Y,dat_sim_w2_init[dat_sim_w2_init$R6==1,]$G1)
  R7_table <- table(dat_sim_w2_init[dat_sim_w2_init$R7==1,]$Y,dat_sim_w2_init[dat_sim_w2_init$R7==1,]$G1)
  R8_table <- table(dat_sim_w2_init[dat_sim_w2_init$R8==1,]$Y,dat_sim_w2_init[dat_sim_w2_init$R8==1,]$G1)


  stop_dat <- list(iteration= iteration,
                   seed=my_seed,
                   sample_freq=list(com_table=com_table, R5_table=R5_table,R6_table=R6_table,R7_table=R7_table,R8_table=R8_table),
                   dataset=dat_sim_w2_init_w2_init)


  return(stop_dat)
}

closeAllConnections()


cl <- makeCluster(num_cores)

clusterEvalQ(cl,{
  library(TPOrd)
})

clusterExport(cl,varlist=c(
  "array_id","m","num_cores",
  "N","cor_YZ_break","num_categories",
  "model_type","Beta0","Beta1",
  "run_simulation","joint_dist"
))


all_results <- list()

# simulations for each n2 value
for (n2 in n2_vector) {
  tasks <- lapply(1:m,function(iter) list(iter=iter,n2=n2))

  total.time <- system.time({
    results <- parLapply(cl,tasks,run_simulation)
  })[3]

  # Save results
  save(total.time,results,file=sprintf("%s_arrayid_%d_n2_%d_cores_%d_iterations_%d.RData",job_name,array_id,n2,num_cores,m))
}

stopCluster(cl)
closeAllConnections()








