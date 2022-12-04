#set local working directory, where files will be imported/read/written
setwd("L:/Priv/Cin/ORD/NERL_NRMRL_JointModel")

#set path to local folder with R packages
my.path <- "C:/Users/rmartin/Desktop/R/win-library/3.2" #Roy
.libPaths(c(my.path,.libPaths())) #Roy

#import packages to be used
library(gridExtra)
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

# Plots
dev.new()
par(mfrow=c(2,2))

hist(as.numeric(GP_y[,1]>0.5), main="", xlab="GP_y")
hist(y, main="", xlab="y")

plot(GP_y[idx2,1], GP_elpd[,1], col=4, cex=0.6, xlab="gp12_y", ylab="gp12_elpd"); abline(a=0, b=1, col="gray")

table(as.factor(as.numeric(GP_y[,1]>0.5)), Y)

table(as.factor(as.numeric(GP_y[idx2,1]>0.5)), Y[idx2])

loo2 <- loo(GPfit, save_psis = TRUE)

# compare
plot(loo1,diagnostic = c("k", "n_eff"),main = "Hierarchical model")
plot(loo2,diagnostic = c("k", "n_eff"),main = "GP model")

comp = loo_compare(list("Hier_model" = loo1, "GP_model" = loo2))
comp

######################################
#####STAN MV probit model fitting#####
######################################

#Create dso file for Stan 
MVPfitDSO <- stan_model(file = "models/MV.stan")

#Data List for Stan input for real data
MVPdataList <- list(y=y[,c(1,29,32,34,42, #darters
                           10,15,51,54,56,57,64,69,72, #minnows
                           6,23,40,62,76, #suckers
                           14,67,77, #catfishes
                           33,36,49,60,65, #sunfishes
                           78)], #TP
                    x=x, 
                    D=ncol(y[,c(1,29,32,34,42, #darters
                                10,15,51,54,56,57,64,69,72, #minnows
                                6,23,40,62,76, #suckers
                                14,67,77, #catfishes
                                33,36,49,60,65, #sunfishes
                                78)]), 
                    x_m_1=as.matrix(cbind(seq(min(x[,1]),max(x[,1]),length.out=M),rep(0, M), rep(0, M))),
                    x_m_2=as.matrix(cbind(rep(0, M), seq(min(x[,2]),max(x[,2]),length.out=M), rep(0, M))),
                    x_m_3=as.matrix(cbind(rep(0, M), rep(0, M), seq(min(x[,3]),max(x[,3]),length.out=M))),
                    K=ncol(x), 
                    N=N,
                    NN=NN,
                    M=M,
                    x_data=x,
                    x_predict=x_predict)#real data 


#Run Stan model with input data and create model object
MVPfit <- sampling(object=MVPfitDSO,
                   data=MVPdataList,
                   chains=4,
                   iter=1000,
                   cores=4,
                   thin=1,
                   init=0,
                   control = list( #add options to improve sampling (divergent transitions)
                     adapt_delta=0.98, #default=0.8
                     max_treedepth =15 #default= 10
                   )
)


#print outputs in R
print(MVPfit, pars=c("beta0","beta","Omega","lp__"), digits=2)

#Print model outputs to text file
options(max.print=999999)
sink(file="MVP_EFLMR_cooc_Final.txt")
print(MVPfit,digits=2)
sink(NULL)

#Extract outputs to model object for post processing in R
la <- extract(MVPfit)

#remove original stanfit model and dso objects from environment to free space
rm(MVPfit)
rm(MVPfitDSO)

#post-checks

#check prevalence
dim(la$z_new)
count_check <- matrix(rep(NA, 2000*28), nrow=28, ncol=2000)
for(i in 1:28){
  for(j in 1:2000){
    count_check[i, j] <- sum(rbinom(n=83, size=1,prob=pnorm(la$z_new[j,i,])))
  }
}

dev.off()
tiff("species_prevalence_check.tiff", width = 8.5, height = 8.5, units = 'in', res = 800, compression = 'lzw')
par(mfrow=c(7,4),mar=c(2,2,1,1))
for(p in 1:28){
  hist(count_check[p,],breaks=20, xlim=c(0,70), xlab='count',cex.main=1,col='gray',
       main=colnames(MVPdataList$y)[p])
  abline(v=colSums(MVPdataList$y)[p], col='red',lwd=2)
}
dev.off()

#posterior intervals and measured values for each response
ppc_table <- list()

for(i in 1:28){ #28 responses
  ppc_table[[i]] <- data.frame(apply(post.table(la$z_new, species=i), 2, pnorm)) #TP is 28th column in MVPdataList$y
  ppc_table[[i]]$Measured <- MVPdataList$y[,i]
  ppc_table[[i]] <- ppc_table[[i]][order(ppc_table[[i]]$Median),]
  ppc_table[[i]]$Site <- seq(1,nrow(ppc_table[[i]]),1)
}

plot_ppc <- list()
for(i in 1:28){
  plot_ppc[[i]] <- ggplot(ppc_table[[i]],aes(x=Site,y=Median))
  plot_ppc[[i]] <- plot_ppc[[i]] + geom_linerange(ymin=ppc_table[[i]]$L90, ymax=ppc_table[[i]]$U90, size=0.5, color="gray")
  plot_ppc[[i]] <- plot_ppc[[i]] + geom_linerange(ymin=ppc_table[[i]]$L50, ymax=ppc_table[[i]]$U50, size=1, color='dark gray')
  plot_ppc[[i]] <- plot_ppc[[i]] + geom_point(aes(x=Site,y=Measured), size=1)
  plot_ppc[[i]] <- plot_ppc[[i]] + labs(title=colnames(MVPdataList$y)[i], x='site',y=expression(p))
  plot_ppc[[i]] <- plot_ppc[[i]] + theme(text = element_text(size=10, angle=360))
  plot_ppc[[i]] <- plot_ppc[[i]] + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
  plot_ppc[[i]] <- plot_ppc[[i]] + expand_limits(y=c(0, 1))
  plot_ppc[[i]] <- plot_ppc[[i]] + theme_bw()
}

dev.off()
tiff("PPC.tiff", width = 8.5, height = 8.25, units = 'in', res = 800, compression = 'lzw')
do.call("grid.arrange", c(plot_ppc, ncol=4))
dev.off()


#########################
#####Parameter plots#####
#########################

#Assign species to groups for plot legend (real data)
speciesGroups <- data.frame(cbind(colnames(MVPdataList$y),
                                  c(rep("Darters",5),
                                    rep("Minnows",9),
                                    rep("Suckers",5),
                                    rep("Catfishes",3),
                                    rep("Sunfishes",5),
                                    "TP"),
                                  c(rep("darkolivegreen3",5),
                                    rep("cadetblue3",9),
                                    rep("plum3",5),
                                    rep("tan",3),
                                    rep("orange2",5),
                                    "indianred3")))


colnames(speciesGroups) <- c("Species","Group","Color")
speciesGroups

#Make plots for each beta parameter
betaHPD <- data.frame(post.table.3d(la$beta,rd=2,vd=3,v=1,rn=colnames(MVPdataList$y)))

betaHPD$Species <- factor(rownames(betaHPD),levels=rownames(betaHPD))
betaHPD$Color <- speciesGroups[match(betaHPD$Species,speciesGroups$Species),]$Color

p1 <- ggplot(betaHPD,aes(x=Species,y=Median)) + geom_point(shape=108, color=betaHPD$Color, size=2)
p1 <- p1 + geom_linerange(ymin=betaHPD$L90, ymax=betaHPD$U90, color=betaHPD$Color, size=0.5) +coord_flip()
p1 <- p1 + geom_linerange(ymin=betaHPD$L50, ymax=betaHPD$U50, color=betaHPD$Color, size=1)
p1 <- p1 + geom_hline(yintercept=0,linetype=1,size=0.5,color='darkgray') 
p1 <- p1 + expand_limits(y=c(-3,15))
p1 <- p1 + labs(x="",y=expression(beta[log(AREA)]))
p1 <- p1 + theme_bw()
p1 <- p1 + theme(text = element_text(size=8))
p1


betaHPD <- data.frame(post.table.3d(la$beta,rd=2,vd=3,v=2,rn=colnames(MVPdataList$y)))

betaHPD$Species <- factor(rownames(betaHPD),levels=rownames(betaHPD))
betaHPD$Color <- speciesGroups[match(betaHPD$Species,speciesGroups$Species),]$Color

p2 <- ggplot(betaHPD,aes(x=Species,y=Median)) + geom_point(shape=108, color=betaHPD$Color, size=2)
p2 <- p2 + geom_linerange(ymin=betaHPD$L90, ymax=betaHPD$U90, color=betaHPD$Color, size=0.5) +coord_flip()
p2 <- p2 + geom_linerange(ymin=betaHPD$L50, ymax=betaHPD$U50, color=betaHPD$Color, size=1)
p2 <- p2 + geom_hline(yintercept=0,linetype=1,size=0.5,color='darkgray') 
p2 <- p2 + expand_limits(y=c(-2,5))
p2 <- p2 + labs(x="",y=expression(beta[CROP]))
p2 <- p2 + theme_bw()
p2 <- p2 + theme(text = element_text(size=8))
p2


betaHPD <- data.frame(post.table.3d(la$beta,rd=2,vd=3,v=3,rn=colnames(MVPdataList$y)))

betaHPD$Species <- factor(rownames(betaHPD),levels=rownames(betaHPD))
betaHPD$Color <- speciesGroups[match(betaHPD$Species,speciesGroups$Species),]$Color

p3 <- ggplot(betaHPD,aes(x=Species,y=Median)) + geom_point(shape=108, color=betaHPD$Color, size=2)
p3 <- p3 + geom_linerange(ymin=betaHPD$L90, ymax=betaHPD$U90, color=betaHPD$Color, size=0.5) +coord_flip()
p3 <- p3 + geom_linerange(ymin=betaHPD$L50, ymax=betaHPD$U50, color=betaHPD$Color, size=1)
p3 <- p3 + geom_hline(yintercept=0,linetype=1,size=0.5,color='darkgray') 
p3 <- p3 + expand_limits(y=c(-3,3))
p3 <- p3 + labs(x="",y=expression(beta[SEP]))
p3 <- p3 + theme_bw()
p3 <- p3 + theme(text = element_text(size=8))
p3


beta0HPD <- data.frame(post.table.2d(la$beta0, rn=colnames(MVPdataList$y)))

beta0HPD$Species <- factor(rownames(beta0HPD),levels=rownames(beta0HPD))
beta0HPD$Color <- speciesGroups[match(beta0HPD$Species,speciesGroups$Species),]$Color
beta0HPD$Group <- speciesGroups[match(beta0HPD$Species,speciesGroups$Species),]$Group

pI <- ggplot(beta0HPD,aes(x=Species,y=Median)) + geom_point(shape=108, color=beta0HPD$Color, size=2)
pI <- pI + geom_linerange(ymin=beta0HPD$L90, ymax=beta0HPD$U90, color=beta0HPD$Color, size=0.5) +coord_flip()
pI <- pI + geom_linerange(ymin=beta0HPD$L50, ymax=beta0HPD$U50, color=beta0HPD$Color, size=1)
pI <- pI + geom_hline(yintercept=0,linetype=1,size=0.5,color='darkgray') 
pI <- pI + expand_limits(y=c(-3,2))
pI <- pI + labs(x="",y=expression(beta[0]))
pI <- pI + theme_bw()
pI <- pI + scale_color_discrete(guide=FALSE)
pI <- pI + theme(text = element_text(size=8))
pI

#Extract Legend
#Fake plot
pFake <- ggplot(beta0HPD,aes(x=Species,y=Median)) + geom_point(shape=108, color=beta0HPD$Color, size=2)
pFake <- pFake + geom_errorbar(data=beta0HPD, aes(ymin=L90, ymax=U90, color=Group), size=1)
pFake <- pFake + scale_color_manual(name="Family", values=c("tan","darkolivegreen3","cadetblue3","plum3","orange2","indianred3"))
pFake + theme(legend.position="top", legend.direction="vertical")+ coord_flip()
pFake + theme(text = element_text(size=8))

#legend function
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

my_legend <- g_legend(pFake)

#Creat figure of plot
dev.off()
tiff("Betas.tiff", width = 8.5, height = 5, units = 'in', res = 800, compression = 'lzw')
grid.arrange(arrangeGrob(pI, p1, p2,
                         p3), my_legend,
             ncol=2,
             widths=c(3.2,0.8),
             respect=FALSE)
dev.off()

#Correlation graph#
###################
#load qgraph package
library(qgraph)

#Correlation quantiles function for Omega

post.test <- function(mat, probs = c(0.05, 0.95)) {
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  tmp <- matrix(NA, nrow=prod(nrow(mat))-nrow(mat),ncol=2)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- quantile(la$Omega[, i, j], probs = probs)
      p.mat[i, j] <- p.mat[j, i] <- prod(sign(tmp))<=0
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp[2]
    }
    
  }
  return(list(as.matrix(p.mat <- ifelse(p.mat==0,1,0)), as.matrix(lowCI.mat), as.matrix(uppCI.mat)))
}


#creat array with mode of correlations for each response
corOmega <- matrix(NA, MVPdataList$D, MVPdataList$D)
for (d in 1:MVPdataList$D){
  corOmega[d,] <- apply(la$Omega[,d,], 2, median)
}

colnames(corOmega) <- colnames(MVPdataList$y)
rownames(corOmega) <- colnames(MVPdataList$y)

#Calculate HPD for correlations
PrOmega <- post.test(la$Omega, prob= c(0.05, 0.90))

#Create figure from plot
qpNoWQ <-qgraph(corOmega*(PrOmega[[1]]),
                mar=c(1.5,1.5,1.5,0.5),
                graph='cor',
                vsize=6,
                esize=6,
                labels=colnames(MVPdataList$y),
                groups= as.factor(speciesGroups$Group),
                edge.labels=TRUE,
                edge.label.cex=0.6,
                posCol= "darkblue",
                negCol= "darkred",
                color=c("tan","darkolivegreen3","cadetblue3","plum3","orange2","indianred3"),
                shape='ellipse',
                layout='circular',
                overlay=FALSE,
                legend=TRUE,
                legend.cex=0.5)

dev.off()
tiff("Omega_Pr90.tiff", width = 8.5, height = 6, units = 'in', res = 800, compression = 'lzw')
qgraph(corOmega*(PrOmega[[1]]),
       mar=c(1.5,1.5,1.5,0.5),
       graph='cor',
       vsize=6,
       esize=6,
       labels=colnames(MVPdataList$y),
       groups= as.factor(speciesGroups$Group),
       edge.labels=TRUE,
       edge.label.cex=0.6,
       posCol= "darkblue",
       negCol= "darkred",
       color=c("tan","darkolivegreen3","cadetblue3","plum3","orange2","indianred3"),
       shape='ellipse',
       layout='circular',
       overlay=FALSE,
       legend=TRUE,
       legend.cex=0.5)
dev.off()
############################
####END parameter plots#####
############################

#marginal effect plots

#Code to create plots in ggplot2
#counterfactual plot: new data along x for marginal probs
mar_plot <- function(y_data=ifelse(MVPdataList$y==1,1,0)[,10], #pick species in data
                     x_data=MVPdataList$x[,1], #pick covariate in data
                     y_m=la$z_m_1[,10,],#pick species in posterior predictive dist
                     x_m=MVPdataList$x_m_1[,1], #pick covariate in "new data"
                     title="Silver shiner",
                     title_x="Std(logAREA)",
                     title_y="Pr(y=1)"){
  plotDat <- data.frame(cbind(LCI1=pnorm(apply(y_m, 2, quantile, prob=0.25)),
                              UCI1=pnorm(apply(y_m, 2, quantile, prob=0.75)),
                              LCI2=pnorm(apply(y_m, 2, quantile, prob=0.05)),
                              UCI2=pnorm(apply(y_m, 2, quantile, prob=0.95))
  ))
  
  addDat <- data.frame(cbind(y_data, x_data))
  
  p <- ggplot(plotDat,aes(x_m)) + 
    #geom_point(aes(y=med), size=1, col="black") +
    geom_linerange(aes(ymin=LCI1, ymax=UCI1), size=1.5, alpha=0.2) + 
    geom_linerange(aes(ymin=LCI2, ymax=UCI2), size=1, col="gray20", alpha=0.2)
  
  p <- p + geom_point(data=addDat, aes(y=y_data, x=x_data)) 
  p <- p + labs(x=title_x, y=title_y) + ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
  p
}


#cooc from joint for silverjaw minnow and TP along counterfactuals
library(mvtnorm)

temp1a <- matrix(NA, 2000, 100)#Silverjaw minnow & TP along logAREA
for(i in 1:2000){# 100 points along hypothetical gradient
  for (j in 1:100){ # 2000 draws
    temp1a[i,j] <- pmvnorm(lower=c(0,0), 
                           upper=c(Inf,Inf), 
                           mean=la$z_m_1[i,c(11,28),j],
                           corr=matrix(c(1,la$Omega[i,11,28],la$Omega[i,28,11],1),2,2))[1]
  }
}

temp1b <- matrix(NA, 2000, 100)#Silverjaw minnow & TP along CROP
for(i in 1:2000){# 100 points along hypothetical gradient
  for (j in 1:100){ # 2000 draws
    temp1b[i,j] <- pmvnorm(lower=c(0,0), 
                           upper=c(Inf,Inf), 
                           mean=la$z_m_2[i,c(11,28),j],
                           corr=matrix(c(1,la$Omega[i,11,28],la$Omega[i,28,11],1),2,2))[1]
  }
}

temp1c <- matrix(NA, 2000, 100)#Silverjaw minnow & TP along SEPTIC
for(i in 1:2000){# 100 points along hypothetical gradient
  for (j in 1:100){ # 2000 draws
    temp1c[i,j] <- pmvnorm(lower=c(0,0), 
                           upper=c(Inf,Inf), 
                           mean=la$z_m_3[i,c(11,28),j],
                           corr=matrix(c(1,la$Omega[i,11,28],la$Omega[i,28,11],1),2,2))[1]
  }
}

#counterfactual plot: new data along x for cooc 
joint_plot <- function(y_data=ifelse(MVPdataList$y==1,1,0)[,10], #pick species in data
                       x_data=MVPdataList$x[,1], #pick covariate in data
                       y_m=la$z_m_1[,10,],#pick species in posterior predictive dist
                       x_m=MVPdataList$x_m_1[,1], #pick covariate in "new data"
                       title="Silver shiner",
                       title_x="Std(logAREA)",
                       title_y="Pr(y=1)"){
  plotDat <- data.frame(cbind(LCI1=apply(y_m, 2, quantile, prob=0.25),
                              UCI1=apply(y_m, 2, quantile, prob=0.75),
                              LCI2=apply(y_m, 2, quantile, prob=0.05),
                              UCI2=apply(y_m, 2, quantile, prob=0.95)
  ))
  
  addDat <- data.frame(cbind(y_data, x_data))
  
  p <- ggplot(plotDat,aes(x_m)) + 
    #geom_point(aes(y=med), size=1, col="black") +
    geom_linerange(aes(ymin=LCI1, ymax=UCI1), size=1.5, alpha=0.2) + 
    geom_linerange(aes(ymin=LCI2, ymax=UCI2), size=1, col="gray20", alpha=0.2)
  
  p <- p + geom_point(data=addDat, aes(y=y_data, x=x_data)) 
  p <- p + labs(x=title_x, y=title_y) + ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
  p
}

#make plots
mp11 <- mar_plot(y_data=y[,57],#pick data species 
                 x_data=MVPdataList$x[,1], #pick data covariate
                 y_m=la$z_m_1[,11,], #match post pred z for marginal to species
                 x_m=MVPdataList$x_m_1[,1], #match grid covariate
                 title="Silverjaw minnow",
                 title_x="Std(logAREA)",
                 title_y=expression(P(y[new]==1)))

mp12 <- mar_plot(y_data=y[,57],#pick data species 
                 x_data=MVPdataList$x[,2], #pick data covariate
                 y_m=la$z_m_2[,11,], #match post pred z for marginal to species
                 x_m=MVPdataList$x_m_2[,2], #match grid covariate
                 title="",
                 title_x="Std(CROP)",
                 title_y=expression(P(y[new]==1)))

mp13 <- mar_plot(y_data=y[,57],#pick data species 
                 x_data=MVPdataList$x[,3], #pick data covariate
                 y_m=la$z_m_3[,11,], #match post pred z for marginal to species
                 x_m=MVPdataList$x_m_3[,3], #match grid covariate
                 title="",
                 title_x="Std(SEP)",
                 title_y=expression(P(y[new]==1)))

mp21 <- mar_plot(y_data=y[,78],#pick data species 
                 x_data=MVPdataList$x[,1], #pick data covariate
                 y_m=la$z_m_1[,28,], #match post pred z for marginal to species
                 x_m=MVPdataList$x_m_1[,1], #match grid covariate
                 title="TP",
                 title_x="Std(logAREA)",
                 title_y="")

mp22 <- mar_plot(y_data=y[,78],#pick data species 
                 x_data=MVPdataList$x[,2], #pick data covariate
                 y_m=la$z_m_2[,28,], #match post pred z for marginal to species
                 x_m=MVPdataList$x_m_2[,2], #match grid covariate
                 title="",
                 title_x="Std(CROP)",
                 title_y="")

mp23 <- mar_plot(y_data=y[,78],#pick data species 
                 x_data=MVPdataList$x[,3], #pick data covariate
                 y_m=la$z_m_3[,28,], #match post pred z for marginal to species
                 x_m=MVPdataList$x_m_3[,3], #match grid covariate
                 title="",
                 title_x="Std(SEP)",
                 title_y="")

jp1 <- joint_plot(y_data=ifelse(rowSums(y[,c(57,78)])>1, 1, 0),#pick data species 
                  x_data=MVPdataList$x[,1], #pick data covariate
                  y_m=temp1a, #match post pred z for marginal to species
                  x_m=MVPdataList$x_m_1[,1], #match grid covariate
                  title="Co-occurrence",
                  title_x="Std(logAREA)",
                  title_y="")

jp2 <- joint_plot(y_data=ifelse(rowSums(y[,c(57,78)])>1, 1, 0),#pick data species 
                  x_data=MVPdataList$x[,2], #pick data covariate
                  y_m=temp1b, #match post pred z for marginal to species
                  x_m=MVPdataList$x_m_2[,2], #match grid covariate
                  title="",
                  title_x="Std(CROP)",
                  title_y="")

jp3 <- joint_plot(y_data=ifelse(rowSums(y[,c(57,78)])>1, 1, 0),#pick data species 
                  x_data=MVPdataList$x[,3], #pick data covariate
                  y_m=temp1c, #match post pred z for marginal to species
                  x_m=MVPdataList$x_m_3[,3], #match grid covariate
                  title="",
                  title_x="Std(SEP)",
                  title_y="")

#create figure of counterfactual plots
dev.off()
tiff("Counterfactual_Plots.tiff", width = 8.5, height = 8.5, units = 'in', res = 800, compression = 'lzw')
grid.arrange(mp11, mp21, jp1, mp12, mp22, jp2, mp13, mp23, jp3, ncol=3)
dev.off()


#####################################
#######Predict to all catchments#####
#####################################

#expected marginal prob of presence/absence (median)
predict_catchments_marginal <- matrix(NA, nrow=MVPdataList$NN,ncol=MVPdataList$D)
for(nn in 1:MVPdataList$NN){
  for(d in 1:MVPdataList$D){
    predict_catchments_marginal[nn,d] <- pnorm(median(la$z_predict[,d,nn]))
  }
}

predict_catchments_marginal <- data.frame(predict_catchments_marginal)
colnames(predict_catchments_marginal) <- colnames(MVPdataList$y)
predict_catchments_marginal$COMID <- streamCat$COMID
head(predict_catchments_marginal)

#write predictions to .csv file
write.table(predict_catchments_marginal,file="predict_catchments_marginal.csv",sep=",",row.names=F)

#uncertainty in marginal predictions (proportion posterior where p < 0.1, > 0.5, > 0.9)

predict_catchments_marginal_unc01 <- data.frame(matrix(NA, nrow=MVPdataList$NN,ncol=MVPdataList$D))
predict_catchments_marginal_unc50 <- data.frame(matrix(NA, nrow=MVPdataList$NN,ncol=MVPdataList$D))
predict_catchments_marginal_unc90 <- data.frame(matrix(NA, nrow=MVPdataList$NN,ncol=MVPdataList$D))
for(d in 1:MVPdataList$D){
  predict_catchments_marginal_unc01[,d] <- apply(pnorm(la$z_predict[,d,]), 2, function(x) sum(x < 0.1)/dim(la$z_predict)[1])
  predict_catchments_marginal_unc50[,d] <- apply(pnorm(la$z_predict[,d,]), 2, function(x) sum(x > 0.5)/dim(la$z_predict)[1])
  predict_catchments_marginal_unc90[,d] <- apply(pnorm(la$z_predict[,d,]), 2, function(x) sum(x > 0.9)/dim(la$z_predict)[1])
}

colnames(predict_catchments_marginal_unc01) <- colnames(MVPdataList$y)
colnames(predict_catchments_marginal_unc50) <- colnames(MVPdataList$y)
colnames(predict_catchments_marginal_unc90) <- colnames(MVPdataList$y)

predict_catchments_marginal_unc01$COMID <- streamCat$COMID
predict_catchments_marginal_unc50$COMID <- streamCat$COMID
predict_catchments_marginal_unc90$COMID <- streamCat$COMID

write.table(predict_catchments_marginal_unc01,file="predict_catchments_marginal_unc01.csv",sep=",",row.names=F)
write.table(predict_catchments_marginal_unc50,file="predict_catchments_marginal_unc50.csv",sep=",",row.names=F)
write.table(predict_catchments_marginal_unc90,file="predict_catchments_marginal_unc90.csv",sep=",",row.names=F)



#expected cooccurrence (posterior median) all species from marginal (i.e., cor=0)
predict_catchments_marginal_cooc <- matrix(NA, nrow=MVPdataList$NN,ncol=MVPdataList$D-1)
for (nn in 1:MVPdataList$NN){
  for(d in 1:(MVPdataList$D-1)){
    predict_catchments_marginal_cooc[nn,d] <- median(pnorm(la$z_predict[,d,nn]) * pnorm(la$z_predict[,28,nn]))
  }
}

predict_catchments_marginal_cooc <- data.frame(predict_catchments_marginal_cooc)
colnames(predict_catchments_marginal_cooc) <- colnames(MVPdataList$y[,-28])
predict_catchments_marginal_cooc$COMID <- streamCat$COMID
head(predict_catchments_marginal_cooc)

#write to .csv file
write.table(predict_catchments_marginal_cooc,file="predict_catchments_marginal_cooc.csv",sep=",",row.names=F)

#uncertainty in cooccurrence all species

predict_catchments_marg_cooc_unc01 <- matrix(NA, nrow=MVPdataList$NN,ncol=MVPdataList$D-1)
predict_catchments_marg_cooc_unc50 <- matrix(NA, nrow=MVPdataList$NN,ncol=MVPdataList$D-1)
predict_catchments_marg_cooc_unc90 <- matrix(NA, nrow=MVPdataList$NN,ncol=MVPdataList$D-1)

for(d in 1:(MVPdataList$D-1)){
  predict_catchments_marg_cooc_unc01[,d] <- apply(pnorm(la$z_predict[,d,]) * pnorm(la$z_predict[,28,]), 2, function(x) sum(x < 0.1)/dim(la$z_predict)[1])
  predict_catchments_marg_cooc_unc50[,d] <- apply(pnorm(la$z_predict[,d,]) * pnorm(la$z_predict[,28,]), 2, function(x) sum(x > 0.5)/dim(la$z_predict)[1])
  predict_catchments_marg_cooc_unc90[,d] <- apply(pnorm(la$z_predict[,d,]) * pnorm(la$z_predict[,28,]), 2, function(x) sum(x > 0.9)/dim(la$z_predict)[1])
}

colnames(predict_catchments_marg_cooc_unc01) <- colnames(MVPdataList$y[,-28])
colnames(predict_catchments_marg_cooc_unc50) <- colnames(MVPdataList$y[,-28])
colnames(predict_catchments_marg_cooc_unc90) <- colnames(MVPdataList$y[,-28])

predict_catchments_marg_cooc_unc01 <- data.frame(predict_catchments_marg_cooc_unc01)
predict_catchments_marg_cooc_unc50 <- data.frame(predict_catchments_marg_cooc_unc50)
predict_catchments_marg_cooc_unc90 <- data.frame(predict_catchments_marg_cooc_unc90)

predict_catchments_marg_cooc_unc01$COMID <- streamCat$COMID
predict_catchments_marg_cooc_unc50$COMID <- streamCat$COMID
predict_catchments_marg_cooc_unc90$COMID <- streamCat$COMID

#write to .csv file
write.table(predict_catchments_marg_cooc_unc01,file="predict_catchments_marg_cooc_unc01.csv",sep=",",row.names=F)
write.table(predict_catchments_marg_cooc_unc50,file="predict_catchments_marg_cooc_unc50.csv",sep=",",row.names=F)
write.table(predict_catchments_marg_cooc_unc90,file="predict_catchments_marg_cooc_unc90.csv",sep=",",row.names=F)


#cooc from joint for silverjaw minnow and TP at all catchments
temp2 <- matrix(NA, 815, 2000)#Silverjaw minnow & TP
for(j in 1:815){ #815 catchments
  for (i in 1:2000){ #2000 draws
    temp3[j,i] <- pmvnorm(lower=c(0,0), 
                          upper=c(Inf,Inf), 
                          mean=la$z_predict[i,c(11,28),j],
                          corr=matrix(c(1,la$Omega[i,11,28],la$Omega[i,28,11],1),2,2))[1]
  }
}

#temp3 <- matrix(NA, 815, 2000)#Silver shiner & TP (not run)
#for(j in 1:815){
#  for (i in 1:2000){
#    temp3[j,i] <- pmvnorm(lower=c(0,0), 
#                          upper=c(Inf,Inf), 
#                          mean=la$z_predict[i,c(10,28),j],
#                          corr=matrix(c(1,la$Omega[i,10,28],la$Omega[i,28,10],1),2,2))[1]
#  }
#}

#temp4 <- matrix(NA, 815, 2000)#Silverjaw, Silver, & TP (not run)
#for(j in 1:815){
#  for (i in 1:2000){
#    temp4[j,i] <- pmvnorm(lower=c(0,0,0), 
#                          upper=c(Inf,Inf,Inf), 
#                          mean=la$z_predict[i,c(11,10,28),j],
#                          corr=matrix(c(1,la$Omega[i,11,10],la$Omega[i,11,28],
#                                        la$Omega[i,11,10],1,la$Omega[i,10,28],
#                                        la$Omega[i,11,28],la$Omega[i,10,28],1),
#                                      3,3,byrow=T))[1]
#  }
#}

#temp5 <- matrix(NA, 815, 2000)#Silverjaw & TP OR Silver & TP (not run)
#for(j in 1:815){
#  for (i in 1:2000){
#    temp5[j, i] <- (temp3[j, i] + temp4[j,i]) - temp5[j,i]
#  }
#}


predict_catchments_sjm_cooc <- apply(temp2,1,median)

predict_catchments_sjm_cooc <- data.frame(predict_catchments_sjm_cooc)
colnames(predict_catchments_sjm_cooc) <- "sens_sjm_cooc"
predict_catchments_sjm_cooc$COMID <- streamCat$COMID
head(predict_catchments_sjm_cooc)

write.table(predict_catchments_sjm_cooc,file="predict_catchments_sjm_cooc.csv",sep=",",row.names=F)


predict_catchments_sjm_cooc_unc01 <- apply(temp2, 1, function(x) sum(x < 0.1)/2000)

predict_catchments_sjm_cooc_unc01 <- data.frame(predict_catchments_sjm_cooc_unc01)
colnames(predict_catchments_sjm_cooc_unc01) <- "sjm_cooc_unc01"
predict_catchments_sjm_cooc_unc01$COMID <- streamCat$COMID
head(predict_catchments_sjm_cooc_unc01)


predict_catchments_sjm_cooc_unc50 <- apply(temp2, 1, function(x) sum(x > 0.5)/2000)

predict_catchments_sjm_cooc_unc50 <- data.frame(predict_catchments_sjm_cooc_unc50)
colnames(predict_catchments_sjm_cooc_unc50) <- "sjm_cooc_unc50"
predict_catchments_sjm_cooc_unc50$COMID <- streamCat$COMID
head(predict_catchments_sjm_cooc_unc50)


predict_catchments_sjm_cooc_unc90 <- apply(temp2, 1, function(x) sum(x > 0.9)/2000)

predict_catchments_sjm_cooc_unc90 <- data.frame(predict_catchments_sjm_cooc_unc90)
colnames(predict_catchments_sjm_cooc_unc90) <- "sjm_cooc_unc90"
predict_catchments_sjm_cooc_unc90$COMID <- streamCat$COMID
head(predict_catchments_sjm_cooc_unc90)

#write to .csv file
write.table(predict_catchments_sjm_cooc_unc01,file="predict_catchments_sjm_cooc_unc01.csv",sep=",",row.names=F)
write.table(predict_catchments_sjm_cooc_unc50,file="predict_catchments_sjm_cooc_unc50.csv",sep=",",row.names=F)
write.table(predict_catchments_sjm_cooc_unc90,file="predict_catchments_sjm_cooc_unc90.csv",sep=",",row.names=F)


#Posterior predictive summaries for all catchments for Silverjaw, TP, and cooc
#For example of predictive uncertainty

#TP

SJM_predict_table <- data.frame(apply(post.table(la$z_predict, species=11), 2, pnorm)) #Silverjam is 11th column in MVPdataList$y
TP_predict_table <- data.frame(apply(post.table(la$z_predict, species=28), 2, pnorm)) #TP is the 28th column in MVPdataList$y
COOC_predict_table <- data.frame(post.table.b(temp3)) #cooc of sjm and tp

SJM_predict_table <- SJM_predict_table[order(SJM_predict_table$Median),]
TP_predict_table <- TP_predict_table[order(TP_predict_table$Median),]
COOC_predict_table <- COOC_predict_table[order(COOC_predict_table$Median),]

SJM_predict_table$Site <- seq(1,nrow(SJM_predict_table),1)

TP_predict_table$Site <- seq(1,nrow(TP_predict_table),1)

COOC_predict_table$Site <- seq(1,nrow(COOC_predict_table),1)


#plots
plot_pp_sjm <- ggplot(SJM_predict_table,aes(x=Site,y=Median)) #+ geom_point(shape=3, size=1)
plot_pp_sjm <- plot_pp_sjm + geom_linerange(ymin=SJM_predict_table$L90, ymax=SJM_predict_table$U90, size=0.5, color="gray")
plot_pp_sjm <- plot_pp_sjm + geom_linerange(ymin=SJM_predict_table$L50, ymax=SJM_predict_table$U50, size=1, color='dark gray')
plot_pp_sjm <- plot_pp_sjm + theme_bw()
plot_pp_sjm <- plot_pp_sjm + labs(x='',y=expression(p), title="Silverjaw minnow")
plot_pp_sjm <- plot_pp_sjm + theme(text = element_text(size=12, angle=360))
plot_pp_sjm <- plot_pp_sjm + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
plot_pp_sjm <- plot_pp_sjm + expand_limits(y=c(0, 1))


plot_pp_tp <- ggplot(TP_predict_table,aes(x=Site,y=Median)) #+ geom_point(shape=3, size=1)
plot_pp_tp <- plot_pp_tp + geom_linerange(ymin=TP_predict_table$L90, ymax=TP_predict_table$U90, size=0.5, color="gray")
plot_pp_tp <- plot_pp_tp + geom_linerange(ymin=TP_predict_table$L50, ymax=TP_predict_table$U50, size=1, color='dark gray')
plot_pp_tp <- plot_pp_tp + theme_bw()
plot_pp_tp <- plot_pp_tp + labs(x='',y=expression(p), title="TP exceedance")
plot_pp_tp <- plot_pp_tp + theme(text = element_text(size=12, angle=360))
plot_pp_tp <- plot_pp_tp + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
plot_pp_tp <- plot_pp_tp + expand_limits(y=c(0, 1))


plot_pp_cooc <- ggplot(COOC_predict_table,aes(x=Site,y=Median)) #+ geom_point(shape=3, size=1)
plot_pp_cooc <- plot_pp_cooc + geom_linerange(ymin=COOC_predict_table$L90, ymax=COOC_predict_table$U90, size=0.5, color="gray")
plot_pp_cooc <- plot_pp_cooc + geom_linerange(ymin=COOC_predict_table$L50, ymax=COOC_predict_table$U50, size=1, color='dark gray')
plot_pp_cooc <- plot_pp_cooc + theme_bw()
plot_pp_cooc <- plot_pp_cooc + labs(x='Catchment',y=expression(p), title="Co-occurrence")
plot_pp_cooc <- plot_pp_cooc + theme(text = element_text(size=12, angle=360))
plot_pp_cooc <- plot_pp_cooc + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
plot_pp_cooc <- plot_pp_cooc + expand_limits(y=c(0, 1))


tiff("p_uncertainty_catchments.tiff", width = 9, height = 6.5, units = 'in', res = 800, compression = 'lzw')
grid.arrange(arrangeGrob(plot_pp_sjm, plot_pp_tp, plot_pp_cooc), ncol=1)
dev.off()

