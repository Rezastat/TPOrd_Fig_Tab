library(ggplot2)
library(dplyr)
library(tibble)
library(VGAM)
library(MASS)
library(TPOrd)

# The biopsy dataset
data("biopsy")
df <- biopsy

# column names
colnames(df) <- c("ID","ClumpThickness","CellSize","CellShape",
                  "MarginalAdhesion","SingleEpithelialCellSize","BareNuclei",
                  "BlandChromatin","NormalNucleoli","Mitoses","Class")

# Removing the id and also rows with missing values
df <- df[,-1]  
df <- na.omit(df)

df$Class <- as.numeric(df$Class) 

# The break points for CellSize and CellShape that we define to reduce the categories
cell_size_breaks <- c(0,4,6,8,10)
CellShape_breaks <- c(0,5,8,10)
clump_thickness_breaks <-  c(0,4,6,8,10)

# Categorize CellSize and CellShape and ClumpThickness
df$CellSizeCategory <- cut(df$CellSize,breaks=cell_size_breaks,labels=c(1,2,3,4),include.lowest=TRUE)
df$CellShapeCategory <- cut(df$CellShape ,breaks=CellShape_breaks,labels=c(1,2,3),include.lowest=TRUE)
df$ClumpCategory <- cut(df$ClumpThickness,breaks=clump_thickness_breaks,labels=c(1,2,3,4),include.lowest=TRUE)


df$CellSizeCategory <- as.numeric(df$CellSizeCategory)
df$CellShapeCategory <- as.numeric(df$CellShapeCategory)
df$ClumpCategory <- as.numeric(df$ClumpCategory) 

# The dataset with the outcome variable and the two chosen covariates
dat_sim <- df %>% 
  dplyr::mutate(Y=ClumpCategory,G=CellShapeCategory,Z=CellSizeCategory) %>%
  dplyr::select(Y,G,Z)

# The correlation 
correlation_matrix <- cor(dat_sim,use="pairwise.complete.obs")
print(correlation_matrix)

# Original (raw) data distributions
raw_clump_thickness_plot <- ggplot(df,aes(x=ClumpThickness)) +
  geom_histogram(binwidth=1,fill="skyblue",color="black") +
  labs(title="C) Raw Distribution of Clump Thickness",x="Clump Thickness",y="Count") +
  theme_minimal()

raw_cell_size_plot <- ggplot(df,aes(x=CellSize)) +
  geom_histogram(binwidth=1,fill="orange",color="black") +
  labs(title="B) Raw Distribution of Cell Size",x="Cell Size",y="Count") +
  theme_minimal()

raw_cell_shape_plot <- ggplot(df,aes(x=CellShape)) +
  geom_histogram(binwidth=1,fill="green",color="black") +
  labs(title="A) Raw Distribution of Cell Shape",x="Cell Shape",y="Count") +
  theme_minimal()

# New categorized data distributions
clump_thickness_plot <- ggplot(df,aes(x=ClumpCategory)) +
  geom_bar(fill="skyblue",color="black") +
  labs(title="F) Distribution of Clump Thickness Categories",x="Clump Thickness Category (Y)",y="Count") +
  theme_minimal()

cell_size_plot <- ggplot(df,aes(x=CellSizeCategory)) +
  geom_bar(fill="orange",color="black") +
  labs(title="E) Distribution of Cell Size Categories",x="Cell Size Category (Z)",y="Count") +
  theme_minimal()

cell_shape_plot <- ggplot(df,aes(x=CellShapeCategory)) +
  geom_bar(fill="green",color="black") +
  labs(title="D) Distribution of Cell Shape Categories",x="Cell Shape Category (G)",y="Count") +
  theme_minimal()

################################################################################


set.seed(1001)





##############################################
n2=70
model_type="Adjacent_Category"
N=length(dat_sim$Y)
fam <- acat(reverse =T,parallel=T)
num_categories=length(unique(dat_sim$Y))
############################################



dat_sim <- sample_dataframe(dat_sim,n2)


dat_sim$YZ <- interaction(dat_sim$Y,dat_sim$Z)

dat_list <- TPOrd:::sample_dataframe_balanced_matrix(dat_sim,n2,1000)
dat_sim$fZ <- factor(dat_sim$Z)
dat=dat_sim[ ,c("Y","Z","fZ","YZ")]
dat$YZ <- as.numeric(dat$YZ)


fit= vglm(Y ~ fZ,family=fam,data=dat)
fitcoef=c(fit@coefficients[1: (num_categories-1)],0,fit@coefficients[num_categories: length(fit@coefficients)])





pZ <- table(dat$Z)[1:(length(table(dat$Z))-1)]/N ## pZ

pG <- table(dat_sim$G)[1:(length(table(dat_sim$G))-1)]/N


Kind <- 1
formula_Ha <- Y~G+fZ
miscov_formula <- ~G



resultsgz <- simulateGZ(N,pG,pZ,0.81,num_sim=1000)

p_gz0 <- convertMeanFrequenciesToTable(resultsgz$MeanFrequencies)


### The formula under the null is the same regardless
formula_Ho <- Y~fZ
auxvar=~Z
strataformula_YZ=~YZ
strataformula_Y=~Y
strataformula_Z=~Z

optMethod <- "Par-spec" # c("Par-spec","A-opt","D-opt")[1]

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
  sa <- sample(N,n2)
})
fit2 <- apply(pop2a,2,function(r){
  return(TPOrd:::fitnessTP_ord(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
})



### Put initial population (mix between best outcome dependent candidates and best candidates from SRS)
pop1_Y <- pop1a_Y[,order(fit1_Y)[1:40]]
pop1_YZ <- pop1a_YZ[,order(fit1_YZ)[1:40]]
pop1_Z <- pop1a_Z[,order(fit1_Z)[1:40]]

pop2 <- pop2a[,order(fit2)[1:20]]




GA.sol1_Y <- TPOrd:::optimTP_GA_ord(ncores=1,
                                    formula=formula_Ha,
                                    miscov=miscov_formula,
                                    auxvar=~Z,
                                    model_type=model_type,
                                    n=n2,
                                    dat,
                                    beta=fitcoef,p_gz=p_gz0,
                                    ga.popsize=60,
                                    ga.propelit=0.9,
                                    ga.proptourney=0.9,
                                    ga.ngen=5000,# can change to 500
                                    ga.mutrate=0.001,
                                    ga.initpop=t(cbind(pop1_Y,pop2)),
                                    optimMeasure=optMethod,K.idx=Kind,seed=1)


GA.sol1_YZ <- TPOrd:::optimTP_GA_ord(ncores=1,
                                     formula=formula_Ha,
                                     miscov=miscov_formula,
                                     auxvar=~Z,
                                     model_type=model_type,
                                     n=n2,
                                     dat,
                                     beta=fitcoef,p_gz=p_gz0,
                                     ga.popsize=60,
                                     ga.propelit=0.9,
                                     ga.proptourney=0.9,
                                     ga.ngen=5000,# can change to 500
                                     ga.mutrate=0.001,
                                     ga.initpop=t(cbind(pop1_YZ,pop2)),
                                     optimMeasure=optMethod,K.idx=Kind,seed=1)


GA.sol1_Z <- TPOrd:::optimTP_GA_ord(ncores=1,
                                    formula=formula_Ha,
                                    miscov=miscov_formula,
                                    auxvar=~Z,
                                    model_type=model_type,
                                    n=n2,
                                    dat,
                                    beta=fitcoef,p_gz=p_gz0,
                                    ga.popsize=60,
                                    ga.propelit=0.9,
                                    ga.proptourney=0.9,
                                    ga.ngen=5000,# can change to 500
                                    ga.mutrate=0.001,
                                    ga.initpop= t(cbind(pop1_Z,pop2)),
                                    optimMeasure=optMethod,K.idx=Kind,seed=1)


### Obtain the indicators
dat_sim$R5 <- 1*(1:N %in% pop2a[,which(fit2 == min(fit2) )])
dat_sim$R6 <- rep(0,N); dat_sim$R6[GA.sol1_Y$bestsol] <- 1
dat_sim$R7 <- rep(0,N); dat_sim$R7[GA.sol1_YZ$bestsol] <- 1
dat_sim$R8 <- rep(0,N); dat_sim$R8[GA.sol1_Z$bestsol] <- 1
library(tibble)
dat_sim <- dat_sim %>% 
  rowid_to_column(var="id")


dat_sim$Model <- "AC"

save(dat_sim,file="biopsy_adcat_dat_sim.RData")

################################################################################

load("biopsy_adcat_dat_sim.RData")



dat_sim$Y <- as.integer(dat_sim$Y)
dat_sim$G <- as.integer(dat_sim$G)
dat_sim$Z <- as.integer(dat_sim$Z)

vglmfit.com <- vglm(factor(Y,ordered=T) ~ G+fZ,family=fam,data=dat_sim)
vglmfit.com.Ho <- vglm(factor(Y,ordered=T) ~ fZ,family=fam, data=dat_sim)
out_com0 <- c(coef(vglmfit.com)[num_categories],diag(vcov(vglmfit.com))[num_categories])
Wcom <- out_com0[1]^2/out_com0[2]
LRcom <- as.numeric(2*(logLik(vglmfit.com) - logLik(vglmfit.com.Ho)))
Scom <- Scom <- (score.stat(vglmfit.com,orig.SE=F)^2)[1]
out_com <- c( out_com0,Wcom,Scom,LRcom)
names(out_com)<-c("beta1","var_beta1","W","S","LR")

data1.R1<- dat_sim[dat_sim$R1==1,c("Y","G","Z","fZ")]
data1.R2<- dat_sim[dat_sim$R2==1,c("Y","G","Z","fZ")]
data1.R3<- dat_sim[dat_sim$R3==1,c("Y","G","Z","fZ")]
data1.R4<- dat_sim[dat_sim$R4==1,c("Y","G","Z","fZ")]
data1.R5<- dat_sim[dat_sim$R5==1,c("Y","G","Z","fZ")]
data1.R6<- dat_sim[dat_sim$R6==1,c("Y","G","Z","fZ")]
data1.R7<- dat_sim[dat_sim$R7==1,c("Y","G","Z","fZ")]
data1.R8<- dat_sim[dat_sim$R8==1,c("Y","G","Z","fZ")]

data0.R1<- dat_sim[dat_sim$R1==0,c("Y","Z","fZ")]
data0.R2<- dat_sim[dat_sim$R2==0,c("Y","Z","fZ")]
data0.R3<- dat_sim[dat_sim$R3==0,c("Y","Z","fZ")]
data0.R4<- dat_sim[dat_sim$R4==0,c("Y","Z","fZ")]
data0.R5<- dat_sim[dat_sim$R5==0,c("Y","Z","fZ")]
data0.R6<- dat_sim[dat_sim$R6==0,c("Y","Z","fZ")]
data0.R7<- dat_sim[dat_sim$R7==0,c("Y","Z","fZ")]
data0.R8<- dat_sim[dat_sim$R8==0,c("Y","Z","fZ")]


resHo.R1 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R1,data1.R1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
resHa.R1 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R1,data1.R1,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
index_of_G <- which(names(resHa.R1$theta) == "G")
res.R1 <-c( resHa.R1$theta[index_of_G],resHa.R1$var_theta[index_of_G],resHa.R1$Wobs,resHo.R1$Sobs,2*(resHa.R1$ll-resHo.R1$ll) )


resHo.R2 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R2,data1.R2,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
resHa.R2 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R2,data1.R2,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
index_of_G <- which(names(resHa.R2$theta) == "G")
res.R2 <-c( resHa.R2$theta[index_of_G],resHa.R2$var_theta[index_of_G],resHa.R2$Wobs,resHo.R2$Sobs,2*(resHa.R2$ll-resHo.R2$ll) )

resHo.R3 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R3,data1.R3,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
resHa.R3 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R3,data1.R3,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
index_of_G <- which(names(resHa.R3$theta) == "G")
res.R3 <-c( resHa.R3$theta[index_of_G],resHa.R3$var_theta[index_of_G],resHa.R3$Wobs,resHo.R3$Sobs,2*(resHa.R3$ll-resHo.R3$ll) )

resHo.R4 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R4,data1.R4,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
resHa.R4 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R4,data1.R4,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
index_of_G <- which(names(resHa.R4$theta) == "G")
res.R4 <-c( resHa.R4$theta[index_of_G],resHa.R4$var_theta[index_of_G],resHa.R4$Wobs,resHo.R4$Sobs,2*(resHa.R4$ll-resHo.R4$ll) )

resHo.R5 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R5,data1.R5,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
resHa.R5 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R5,data1.R5,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
index_of_G <- which(names(resHa.R5$theta) == "G")
res.R5 <-c( resHa.R5$theta[index_of_G],resHa.R5$var_theta[index_of_G],resHa.R5$Wobs,resHo.R5$Sobs,2*(resHa.R5$ll-resHo.R5$ll) )

resHo.R6 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R6,data1.R6,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
resHa.R6 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R6,data1.R6,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
index_of_G <- which(names(resHa.R6$theta) == "G")
res.R6 <-c( resHa.R6$theta[index_of_G],resHa.R6$var_theta[index_of_G],resHa.R6$Wobs,resHo.R6$Sobs,2*(resHa.R6$ll-resHo.R6$ll) )

resHo.R7 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R7,data1.R7,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
resHa.R7 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R7,data1.R7,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
index_of_G <- which(names(resHa.R7$theta) == "G")
res.R7 <-c( resHa.R7$theta[index_of_G],resHa.R7$var_theta[index_of_G],resHa.R7$Wobs,resHo.R7$Sobs,2*(resHa.R7$ll-resHo.R7$ll) )

resHo.R8 <- twoPhaseSPML_ord(formula=Y~fZ,miscov=~G,auxvar=~Z,family=fam,data0.R8,data1.R8,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
resHa.R8 <- twoPhaseSPML_ord(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=fam,data0.R8,data1.R8,start.values=NULL,verbose=FALSE,n_second=n2,model_type=model_type,num_categories=num_categories,N=N)
index_of_G <- which(names(resHa.R8$theta) == "G")
res.R8 <-c( resHa.R8$theta[index_of_G],resHa.R8$var_theta[index_of_G],resHa.R8$Wobs,resHo.R8$Sobs,2*(resHa.R8$ll-resHo.R8$ll) )

###############################################################

dat_nai.R1 <- dat_sim[dat_sim$R1==1,c("Y","G","Z","fZ")] # phase 2 data alone
ordinalfit.nai.R1 <- vglm(factor(Y,ordered=T) ~ G+fZ,family=fam, data=dat_nai.R1)
ordinalfit.nai.Ho.R1 <- vglm(factor(Y,ordered=T) ~ fZ,family=fam, data=dat_nai.R1)
out_nai0.R1 <- c(coef(ordinalfit.nai.R1)[num_categories],diag(vcov(ordinalfit.nai.R1))[num_categories])
Wnai.R1 <- out_nai0.R1[1]^2/out_nai0.R1[2]
LRnai.R1 <- as.numeric(2*(logLik(ordinalfit.nai.R1) - logLik(ordinalfit.nai.Ho.R1)))
Snai.R1 <- (score.stat(ordinalfit.nai.R1,orig.SE=T)^2)[1]
out_nai.R1 <- c( out_nai0.R1,Wnai.R1,Snai.R1,LRnai.R1)

dat_nai.R2 <- dat_sim[dat_sim$R2==1,c("Y","G","Z","fZ")] # phase 2 data alone
ordinalfit.nai.R2 <- vglm(factor(Y,ordered=T) ~ G+fZ,family=fam, data=dat_nai.R2)
ordinalfit.nai.Ho.R2 <- vglm(factor(Y,ordered=T) ~ fZ,family=fam, data=dat_nai.R2)
out_nai0.R2 <- c(coef(ordinalfit.nai.R2)[num_categories],diag(vcov(ordinalfit.nai.R2))[num_categories])
Wnai.R2 <- out_nai0.R2[1]^2/out_nai0.R2[2]
LRnai.R2 <- as.numeric(2*(logLik(ordinalfit.nai.R2) - logLik(ordinalfit.nai.Ho.R2)))
Snai.R2 <- (score.stat(ordinalfit.nai.R2,orig.SE=T)^2)[1]
out_nai.R2 <- c( out_nai0.R2,Wnai.R2,Snai.R2,LRnai.R2)

dat_nai.R3 <- dat_sim[dat_sim$R3==1,c("Y","G","Z","fZ")] # phase 2 data alone
ordinalfit.nai.R3 <- vglm(factor(Y,ordered=T) ~ G+fZ,family=fam, data=dat_nai.R3)
ordinalfit.nai.Ho.R3 <- vglm(factor(Y,ordered=T) ~ fZ,family=fam, data=dat_nai.R3)
out_nai0.R3 <- c(coef(ordinalfit.nai.R3)[num_categories],diag(vcov(ordinalfit.nai.R3))[num_categories])
Wnai.R3 <- out_nai0.R3[1]^2/out_nai0.R3[2]
LRnai.R3 <- as.numeric(2*(logLik(ordinalfit.nai.R3) - logLik(ordinalfit.nai.Ho.R3)))
Snai.R3 <- (score.stat(ordinalfit.nai.R3,orig.SE=T)^2)[1]
out_nai.R3 <- c( out_nai0.R3,Wnai.R3,Snai.R3,LRnai.R3)

dat_nai.R4 <- dat_sim[dat_sim$R4==1,c("Y","G","Z","fZ")] # phase 2 data alone
ordinalfit.nai.R4 <- vglm(factor(Y,ordered=T) ~ G+fZ,family=fam, data=dat_nai.R4)
ordinalfit.nai.Ho.R4 <- vglm(factor(Y,ordered=T) ~ fZ,family=fam, data=dat_nai.R4)
out_nai0.R4 <- c(coef(ordinalfit.nai.R4)[num_categories],diag(vcov(ordinalfit.nai.R4))[num_categories])
Wnai.R4 <- out_nai0.R4[1]^2/out_nai0.R4[2]
LRnai.R4 <- as.numeric(2*(logLik(ordinalfit.nai.R4) - logLik(ordinalfit.nai.Ho.R4)))
Snai.R4 <- (score.stat(ordinalfit.nai.R4,orig.SE=T)^2)[1]
out_nai.R4 <- c( out_nai0.R4,Wnai.R4,Snai.R4,LRnai.R4)

dat_nai.R5 <- dat_sim[dat_sim$R5==1,c("Y","G","Z","fZ")] # phase 2 data alone
ordinalfit.nai.R5 <- vglm(factor(Y,ordered=T) ~ G+fZ,family=fam, data=dat_nai.R5)
ordinalfit.nai.Ho.R5 <- vglm(factor(Y,ordered=T) ~ fZ,family=fam, data=dat_nai.R5)
out_nai0.R5 <- c(coef(ordinalfit.nai.R5)[num_categories],diag(vcov(ordinalfit.nai.R5))[num_categories])
Wnai.R5 <- out_nai0.R5[1]^2/out_nai0.R5[2]
LRnai.R5 <- as.numeric(2*(logLik(ordinalfit.nai.R5) - logLik(ordinalfit.nai.Ho.R5)))
Snai.R5 <- (score.stat(ordinalfit.nai.R5,orig.SE=T)^2)[1]
out_nai.R5 <- c( out_nai0.R5,Wnai.R5,Snai.R5,LRnai.R5)

dat_nai.R6 <- dat_sim[dat_sim$R6==1,c("Y","G","Z","fZ")] # phase 2 data alone
ordinalfit.nai.R6 <- vglm(factor(Y,ordered=T) ~ G+fZ,family=fam, data=dat_nai.R6)
ordinalfit.nai.Ho.R6 <- vglm(factor(Y,ordered=T) ~ fZ,family=fam, data=dat_nai.R6)
out_nai0.R6 <- c(coef(ordinalfit.nai.R6)[num_categories-1],diag(vcov(ordinalfit.nai.R6))[num_categories-1])
Wnai.R6 <- out_nai0.R6[1]^2/out_nai0.R6[2]
LRnai.R6 <- as.numeric(2*(logLik(ordinalfit.nai.R6) - logLik(ordinalfit.nai.Ho.R6)))
Snai.R6 <- (score.stat(ordinalfit.nai.R6,orig.SE=T)^2)[1]
out_nai.R6 <- c( out_nai0.R6,Wnai.R6,Snai.R6,LRnai.R6)

dat_nai.R7 <- dat_sim[dat_sim$R7==1,c("Y","G","Z","fZ")] # phase 2 data alone
ordinalfit.nai.R7 <- vglm(factor(Y,ordered=T) ~ G+fZ,family=fam, data=dat_nai.R7)
ordinalfit.nai.Ho.R7 <- vglm(factor(Y,ordered=T) ~ fZ,family=fam, data=dat_nai.R7)
out_nai0.R7 <- c(coef(ordinalfit.nai.R7)[num_categories-1],diag(vcov(ordinalfit.nai.R7))[num_categories-1])
Wnai.R7 <- out_nai0.R7[1]^2/out_nai0.R7[2]
LRnai.R7 <- as.numeric(2*(logLik(ordinalfit.nai.R7) - logLik(ordinalfit.nai.Ho.R7)))
Snai.R7 <- (score.stat(ordinalfit.nai.R7,orig.SE=T)^2)[1]
out_nai.R7 <- c( out_nai0.R7,Wnai.R7,Snai.R7,LRnai.R7)

dat_nai.R8 <- dat_sim[dat_sim$R8==1,c("Y","G","Z","fZ")] # phase 2 data alone
ordinalfit.nai.R8 <- vglm(factor(Y,ordered=T) ~ G+fZ,family=fam, data=dat_nai.R8)
ordinalfit.nai.Ho.R8 <- vglm(factor(Y,ordered=T) ~ fZ,family=fam, data=dat_nai.R8)
out_nai0.R8 <- c(coef(ordinalfit.nai.R8)[num_categories-1],diag(vcov(ordinalfit.nai.R8))[num_categories-1])
Wnai.R8 <- out_nai0.R8[1]^2/out_nai0.R8[2]
LRnai.R8 <- as.numeric(2*(logLik(ordinalfit.nai.R8) - logLik(ordinalfit.nai.Ho.R8)))
Snai.R8 <- (score.stat(ordinalfit.nai.R8,orig.SE=T)^2)[1]
out_nai.R8 <- c( out_nai0.R8,Wnai.R8,Snai.R8,LRnai.R8)

beta1_table_case <- rbind(out_com,out_nai.R1,out_nai.R3,out_nai.R2,out_nai.R4,res.R1,res.R3,res.R2,res.R4,out_nai.R5,out_nai.R8,out_nai.R6,out_nai.R7,res.R5,res.R8,res.R6,res.R7)
beta1_table_case[,"var_beta1"] <- sqrt(beta1_table_case[,"var_beta1"])
beta1_table_case <- round(beta1_table_case,2)


sampling_method <- c("Complete Data",rep(c("SRS","Z-Bal","Y-Bal","YZ-Bal"),2),rep(c("Best-SRS","GA-Z-Bal","GA-Y-Bal","GA-YZ-Bal"),2))
analysis_method <- c("Complete Data",rep(c("Naive","EM","Naive","EM"),each=4))
Case_Study_Res <- data.frame(cbind(sampling_method,analysis_method,beta1_table_case))
colnames(Case_Study_Res) <- c("Sampling Method","Analysis Method","Beta1","SE","Wald","Score","LRT")
rownames(Case_Study_Res) <- NULL


Case_Study_Res$Beta1 <- as.numeric(Case_Study_Res$Beta1)
Case_Study_Res$SE <- as.numeric(Case_Study_Res$SE)
Case_Study_Res$Wald <- as.numeric(Case_Study_Res$Wald)
Case_Study_Res$Score <- as.numeric(Case_Study_Res$Score)
Case_Study_Res$LRT <- as.numeric(Case_Study_Res$LRT)

Case_Study_Res <- Case_Study_Res %>%
  mutate(LowerCI=Beta1 - 1.96 * SE,
         UpperCI=Beta1 + 1.96 * SE)

Case_Study_Res$`Sampling Method` <- factor(Case_Study_Res$`Sampling Method`,
                                           levels=c("Complete Data","SRS","Z-Bal","Y-Bal","YZ-Bal",
                                                      "Best-SRS","GA-Z-Bal","GA-Y-Bal","GA-YZ-Bal"))
Case_Study_Res$`Analysis Method` <- factor(Case_Study_Res$`Analysis Method`,
                                           levels=c("Naive","EM","Complete Data"))

ordered_Case_Study_Res <- Case_Study_Res %>%
  arrange(`Sampling Method`,`Analysis Method`)

ordered_Case_Study_Res$`Analysis Method`[1] <- "Complete Data"

complete_data_CI <- ordered_Case_Study_Res %>%
  filter(`Sampling Method` == "Complete Data") %>%
  dplyr::select(LowerCI,UpperCI)

# Plot for single wave
case_estimate_plot <- ggplot(ordered_Case_Study_Res,aes(x=`Sampling Method`,y=Beta1,color=`Analysis Method`,shape=`Analysis Method`)) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=complete_data_CI$LowerCI,ymax=complete_data_CI$UpperCI),
            fill="gray90",color=NA,alpha=0.5) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=LowerCI,ymax=UpperCI),width=0.2,position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0,linetype="dashed",color="black") +
  labs(title="Estimates with 95% Confidence Intervals with AC",x="Sampling Method",y=expression(hat(beta)[1])) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=20,hjust=1),legend.position="bottom") +
  scale_color_manual(values=c("Complete Data"="black","Naive"="red","EM"="blue"))





