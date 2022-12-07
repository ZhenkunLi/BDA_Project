#set local working directory, where files will be imported/read/written
setwd("L:/Priv/Cin/ORD/NERL_NRMRL_JointModel")

#set path to local folder with R packages
my.path <- "C:/Users/rmartin/Desktop/R/win-library/3.2" #Roy
.libPaths(c(my.path,.libPaths())) #Roy

set.seed(9527)
#import packages to be used
library(gridExtra)
library(bayesplot)
library(qgraph)
library(coda)
library(rstan)
library(ggplot2)
library(loo)
library(posterior)

#set stan options for parallel processing of HMC chains
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Examine priors
#Functions for student t distribution
qt_ls <- function(prob, df, mu, sigma) qt(prob, df)*sigma + mu
dt_ls <- function(x, df, mu, sigma) 1/sigma * dt((x - mu)/sigma, df)
pt_ls <- function(x, df, mu, sigma) pt((x - mu)/sigma, df)
rt_ls <- function(n, df, mu, sigma) rt(n,df)*sigma + mu

#Gelman 2006 suggests abs(beta) < 5 for coefficients
#Find equivalent for probit
plogis(5)
qnorm(plogis(5)) #So, abs(beta) in probit < approx 2.47

#
qlogis(10e-9)
qnorm(plogis(qlogis(10e-9))) #So, abs(beta0) in probit < approx 5.6

#Gelman suggests T(df=1, mu=0, sigma=2.5) for logistic regression prior on coefficients
pt_ls(5, df=1, mu=0, sigma=2.5) #about 70% of density over values -5 to 5
pt_ls(qnorm(plogis(5)), df=1, mu=0, sigma=1.25) #equivalent cauchy for probit

#Gelman suggests T(df=1, mu=0, sigma=10) for prior on interecept or success prob
#for average case to be between 10e-9 and 1-10e-9.
pt_ls(-qlogis(10e-9), df=1, mu=0, sigma=10) #about 85% of density over values less than 20
pt_ls(-qnorm(plogis(qlogis(10e-9))), df=1, mu=0, sigma=3) #about 85% of density over values less than 6

#plot priors
plot(x=seq(-10,10,0.1),dt_ls(x=seq(-10,10,0.1), df=1,  mu=0, sigma=3),
     typ='l',
     xlim=c(-10,10),
     xlab=expression(prior: beta[0]),
     ylab="density")

plot(x=seq(-5,5,0.1),dt_ls(x=seq(-5,5,0.1), df=1,  mu=0, sigma=1.25),
     typ='l',
     xlim=c(-5,5),
     xlab=expression(prior: beta[k]),
     ylab="density")


#Function to scale and center predictors based on Gelman et al. 2006
MyStd <- function(x){ (x-mean(x))/(sd(x)*2)} #apply(x,2,MyStd)

#Functions for summarizing posteriors#
######################################

#compute posterior quantiles by catchment for a species from 3d array (e.g., z_predict, z_new)
post.table <- function(x, species=1) {
  n <- dim(x)[3]
  mat <- matrix(NA, n, 5)
  tmp <- cbind(apply(as.mcmc(x[, species, ]), 2, quantile, probs = 0.05),
               apply(as.mcmc(x[, species, ]), 2, quantile, probs = 0.25),
               apply(as.mcmc(x[, species, ]), 2, median),
               apply(as.mcmc(x[, species, ]), 2, quantile, probs = 0.75),
               apply(as.mcmc(x[, species, ]), 2, quantile, probs = 0.90))
  
  mat <- data.frame(tmp)
  colnames(mat) <- c("L90","L50","Median","U50","U90")
  return(mat)
}

#compute posterior quantiles by catchment from 2d matrix 
post.table.b <- function(x, rd=1) {
  n <- dim(x)[rd]
  mat <- matrix(NA, n, 5)
  tmp <- cbind(apply(as.mcmc(x), rd, quantile, probs = 0.05),
               apply(as.mcmc(x), rd, quantile, probs = 0.25),
               apply(as.mcmc(x), rd, median),
               apply(as.mcmc(x), rd, quantile, probs = 0.75),
               apply(as.mcmc(x), rd, quantile, probs = 0.90))
  
  mat <- data.frame(tmp)
  colnames(mat) <- c("L90","L50","Median","U50","U90")
  return(mat)
}

#from 2d array for beta0
post.table.2d <- function(x, rd=2, rn=FALSE) {
  n <- dim(x)[rd]
  mat <- matrix(NA, n, 5)
  tmp <- cbind(apply(as.mcmc(x), 2, quantile, probs = 0.05),
               apply(as.mcmc(x), 2, quantile, probs = 0.25),
               apply(as.mcmc(x), 2, median),
               apply(as.mcmc(x), 2, quantile, probs = 0.75),
               apply(as.mcmc(x), 2, quantile, probs = 0.95))
  for (i in 1:n){
    for (j in 1:5){
      mat[i,j] <- tmp[i,j]
    }
  }
  mat <- data.frame(mat)
  colnames(mat) <- c("L90","L50","Median","U50","U90")
  rownames(mat) <- ifelse(rn==FALSE,rep(1:rd,1),rn)
  mat <- mat[with(mat,order(Median)),]
  return(mat)
}

#from 3d array from beta (e.g. i=2000 sims, j=3 parameters, k=5 groups)
post.table.3d <- function(x, rd=2, vd=3, v=1, 
                          prob_1 = 0.95, prob_2 = 0.8,
                          point_est='median',
                          rn=FALSE){
  n <- dim(x)[rd]
  p <- dim(x)[vd]
  mat <- matrix(NA, n, 5)
  tmp <- cbind(apply(as.mcmc(x[, , v]), 2, quantile, prob = 0.05),
               apply(as.mcmc(x[, , v]), 2, quantile, prob = 0.25),
               apply(x[, ,v ], 2, median),
               apply(as.mcmc(x[, , v]), 2, quantile, prob = 0.75),
               apply(as.mcmc(x[, , v]), 2, quantile, prob = 0.95))
  for (i in 1:n){
    for (j in 1:5){
      mat[i,j] <- tmp[i,j]
    }
  }
  mat <- data.frame(mat)
  colnames(mat) <- c("L90","L50","Median","U50","U90")
  rownames(mat) <- ifelse(rn==FALSE,rep(1:rd,1),rn)
  mat <- mat[with(mat,order(Median)),]
  return(mat)
}

#################################################
#####Import data from working directory files#####
#################################################
fish <- read.csv(file="data/fish_data_OEPA_2012.csv", header=T)
streamCat <- read.csv(file="data/NHD_Plus_StreamCat.csv", header=T)
water <- read.csv(file="data/OEPA_WATER_2012.csv",header=T)

#merge files to working dataset
fishCat <- merge(fish, streamCat, by="COMID", all.x=TRUE)
fishWater <- merge(fishCat,water,by="STORET",all.x=TRUE)
modelData <- fishWater#[!is.na(fishWater$NO2),]

#Convert TP data to binary
#based on a management threshold (Table 9; Miltner 2010)
modelData$P_exc <- ifelse(modelData$P>=(0.1),1,0)

#convert fish abundance to presence-absence (binary; 0 or 1)
y <- as.matrix(ifelse(modelData[,c(6:82,146)]>0,1,0))

#response variable for TP
Y <- modelData[,146]

#export .csv file of working dataset to working directory
#for viewing and QA checking
write.table(modelData, 
            file="model_input_data.csv", 
            sep=",", 
            row.names = FALSE)

#set inputs for Stan model
#number of observations in model
N <- nrow(modelData)
#number of new data points to predict
#the model. For mapping predictions
#across sampled and unsampled catchments.
NN <- nrow(streamCat) #number of observations for predictions 

#number of simulated observations for counterfactual plots
M <- 100

#create matrix of predictors for Stan
x <- modelData[,c(85, 97, 100)]

#check for strong correlations (i.e., colinearities)
round(cor(x),1)

#Apply standardization to predictors
x <- as.matrix(cbind(MyStd(log(x[,1])),
                     MyStd(x[,2]),
                     MyStd(x[,3])
))

#Apply same standardization to variables used
#to predict to unsampled catchments
x_predict <-  as.matrix(cbind(MyStd(log(streamCat[,4])),
                              MyStd(streamCat[,16]),
                              MyStd(streamCat[,19])
))

#Check values and verify standardization
apply(x,2,sd)
apply(x_predict,2,sd)

#Assign simple variable names to predictor inputs
#logAREA = log(cumulative drainage area for catchment)
#CROP = percent cumulative crop cover for catchment
#SEPR = normalized density of septics for catchment
colnames(x) <- c("logAREA","CROP","SEP")

colnames(x_predict) <- c("logAREA","CROP","SEP")

#write csv for mapping
write.table(cbind(streamCat[,c(1, 4, 16, 19)],x_predict), 
            file="standardized_predictors.csv", 
            sep=",", 
            row.names = FALSE)


#######################
####END import data####
#######################

##########################
#####Data Description#####
##########################
library(corrplot)
descData <- read.csv("standardized_predictors.csv", header = TRUE)
summary(Y)
corrplot(cor(modelData[,c(85, 97, 100, 112:115)]))


idx2 <- sample(1:dim(modelData)[1], 15)
N2 <- length(idx2)

idx1 <- setdiff(1:dim(modelData)[1], idx2)
N1 <- length(idx1)

trials <- rep(1,N1+N2)

# Hier ----------------------------------------------------------------------

######################################
####STAN Hier probit model fitting####
######################################


HierfitDSO <- stan_model(file = "models/Hier_model_TP.stan")

#Data List for Stan input for real data
HierdataList <- list(y=y,
                     Y=Y,
                     x=x, 
                     X=x_predict,
                     x_m_1=as.matrix(cbind(seq(min(x[,1]),max(x[,1]),length.out=M),rep(0, M), rep(0, M))),
                     x_m_2=as.matrix(cbind(rep(0, M), seq(min(x[,2]),max(x[,2]),length.out=M), rep(0, M))),
                     x_m_3=as.matrix(cbind(rep(0, M), rep(0, M), seq(min(x[,3]),max(x[,3]),length.out=M))),
                     K=ncol(x), 
                     D=ncol(y), 
                     N1=N1,
                     N2=N2,
                     NN=NN,
                     M=M,
                     idx1=idx1,
                     idx2=idx2,
                     trials=trials,
                     x_data=x,
                     x_predict=x_predict)#real data 


#Run Stan model with input data and create model object
Hierfit <- sampling(object=HierfitDSO,
                    data=HierdataList,
                    chains=4,
                    iter=15000,
                    cores=4,
                    thin=1,
                    init=0,
                    seed=0,
                    control = list( #add options to improve sampling (divergent transitions)
                      adapt_delta=0.98, #default=0.8
                      max_treedepth =15 #default= 10
                    )
)


#print outputs in R
print(Hierfit, digits=2)

#Print model outputs to text file
options(max.print=999999)
sink(file="Hier_EFLMR_cooc_Final.txt")
print(Hierfit,digits=2)
sink(NULL)

Hier_f <- summary(Hierfit, pars = c("f"), probs = c(0.025, 0.5, 0.975), digits_summary = 4)$summary
Hier_y <- summary(Hierfit, pars = c("y_predict"), probs = c(0.025, 0.5, 0.975), digits_summary = 4)$summary
Hier_finv <- summary(Hierfit, pars = c("f_invlogit"), probs = c(0.025, 0.5, 0.975), digits_summary = 4)$summary
Hier_elpd <- summary(Hierfit, pars = c("log_y_new"), probs = c(0.025, 0.5, 0.975), digits_summary = 4)$summary
Hier_time <- sum(get_elapsed_time(Hierfit)[2])/(Hierfit@sim$iter - Hierfit@sim$warmup)

# Plots
dev.new()
par(mfrow=c(2,2))

hist(as.numeric(Hier_y[,1]>0.5), main="", xlab="Hier_y")
hist(y, main="", xlab="y")

plot(Hier_y[idx2,1], Hier_elpd[,1], col=4, cex=0.6, xlab="Hier12_y", ylab="Hier12_elpd"); abline(a=0, b=1, col="gray")

table(as.factor(as.numeric(Hier_y[,1]>0.5)), Y)

table(as.factor(as.numeric(Hier_y[idx2,1]>0.5)), Y[idx2])

loo1 <- loo(Hierfit, save_psis = TRUE)


# GP ----------------------------------------------------------------------

####################################
####STAN GP probit model fitting####
####################################

#Create dso file for Stan 
GPfitDSO <- stan_model(file = "models/GP.stan")

#Data List for Stan input for real data
GPdataList <- list(y=y,
                   Y=Y,
                    x=x, 
                    X=x_predict,
                    x_m_1=as.matrix(cbind(seq(min(x[,1]),max(x[,1]),length.out=M),rep(0, M), rep(0, M))),
                    x_m_2=as.matrix(cbind(rep(0, M), seq(min(x[,2]),max(x[,2]),length.out=M), rep(0, M))),
                    x_m_3=as.matrix(cbind(rep(0, M), rep(0, M), seq(min(x[,3]),max(x[,3]),length.out=M))),
                    K=ncol(x), 
                    D=ncol(y), 
                    N1=N1,
                    N2=N2,
                    NN=NN,
                    M=M,
                    idx1=idx1,
                    idx2=idx2,
                    trials=trials,
                    x_data=x,
                    x_predict=x_predict)#real data 


#Run Stan model with input data and create model object
GPfit <- sampling(object=GPfitDSO,
                   data=GPdataList,
                   chains=4,
                   iter=1000,
                   cores=4,
                   thin=1,
                   init=0,
                   seed=0,
                   control = list( #add options to improve sampling (divergent transitions)
                     adapt_delta=0.98, #default=0.8
                     max_treedepth =15 #default= 10
                   )
)


#print outputs in R
print(GPfit, digits=2)

#Print model outputs to text file
options(max.print=999999)
sink(file="GP_EFLMR_cooc_Final.txt")
print(GPfit,digits=2)
sink(NULL)

GP_f <- summary(GPfit, pars = c("f"), probs = c(0.025, 0.5, 0.975), digits_summary = 4)$summary
GP_y <- summary(GPfit, pars = c("y_predict"), probs = c(0.025, 0.5, 0.975), digits_summary = 4)$summary
GP_finv <- summary(GPfit, pars = c("f_invlogit"), probs = c(0.025, 0.5, 0.975), digits_summary = 4)$summary
GP_elpd <- summary(GPfit, pars = c("log_y_new"), probs = c(0.025, 0.5, 0.975), digits_summary = 4)$summary
GP_time <- sum(get_elapsed_time(GPfit)[2])/(GPfit@sim$iter - GPfit@sim$warmup)

#Extract outputs to model object for post processing in R
gp <- extract(GPfit)

#remove original stanfit model and dso objects from environment to free space
rm(GPfit)
rm(GPfitDSO)

# Plots
dev.new()
par(mfrow=c(2,2))

# Posterior predictive checking
dev.off()
pdf("post_check.pdf", width = 4, height = 4)
ppc_stat(Y, gp$f_invlogit)
dev.off()

# Posterior predictive checking-histplot 
dev.off()
pdf("hist_check.pdf", width = 5, height = 5)
ppc_hist(Y, gp$y_predict[1:8, ])
dev.off()


ppc_stat(Y, gp$f_invlogit)


hist(as.numeric(GP_y[,1]>0.5), main="", xlab="GP_y")
hist(y, main="", xlab="y")

plot(GP_y[idx2,1], GP_elpd[,1], col=4, cex=0.6, xlab="gp12_y", ylab="gp12_elpd"); abline(a=0, b=1, col="gray")

table(as.factor(as.numeric(GP_y[,1]>0.5)), Y)

table(as.factor(as.numeric(GP_y[idx2,1]>0.5)), Y[idx2])

loo_GP <- loo(GPfit, save_psis = TRUE)

