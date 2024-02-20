library(runjags)
library(coda)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grafify)
library(grid)
library(gridExtra)
# load functions - we'll need the 'split mcmc' and 'wide-to-stacked' functions
functions_to_load <- list.files(
  "./R/functions/",
  full.names = TRUE
)
for(fn in functions_to_load){
  source(fn)
}

# read in the results and site-level covariates
cooc_covs <- read.csv("data/cooc_covs.csv")
canopy_mod <- readRDS("./results/canopy_model.RDS")

# convert the model file to an mcmc object; use MCMCvis to get summary statistics
mc <- as.mcmc(canopy_mod)
library(MCMCvis)
MCMCvis::MCMCsummary(mc,
                     params = c("b1_mu", "b1_sd", "b1"),       # this is just a subset of parameters
                     probs = c(0.025, 0.25, 0.5, 0.75, 0.975), 
                     round = 2)

# convert the mcmc object to a matrix
mc <- as.matrix(mc)
# sub-sample the mcmc matrix a bit as we don't really need to make predictions with 
# all 60K samples
set.seed(554)
mc_sub <- mc[sample(1:nrow(mc), 10000), ]
# and use split_mcmc (thanks @MFidino !)
mc <- split_mcmc(mc_sub)
rm(mc_sub)
# check out the dimensions of one of the list elements. for this example (the slope term for my 
# occupancy covariate), it should have a number of rows equal to the number of MCMC samples 
# (10000) & a number of columns equal to the number of cities sampled.
dim(mc$a1)
# [1] 10000  9
# the among-city terms (e.g., a0_mu) will have a 1 as the 2nd term
# the species-specific terms (p0) will have dimensions equal to the number of MCMC samples,
# the number of species sampled, then the number of cities

# in the next fiew big chunks of code, I am putting together predicted occupancy probabilities for
# each city. I will explain in detail this first example, then the sections are all fairly
# similar to one another
#### Austin, TX ####
# generate a sequence of impervious surface values for the city in question
range(cooc_covs$gmcCanopy[1:26])
# [1] -19.04090  46.17025
# I need to choose some 'pretty' numbers based on that
predV <- seq(-19, 47, length.out = 200)

# Generating a matrix that will contain values for predicting. This includes our sequence of
# impervious surface values, and adding a row of 1's for the intercept
pred_mat <- cbind(
  1,
  scale(predV)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  scale(predV),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_psi <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)

# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)
# fill it in. No species is easy, it's just a 1.
pred_psi[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi[,,2] <- exp(
  cbind(mc$a0[,1],mc$a1[,1]) %*% pred_mat
)
# Species 2 (gray squirrel)
pred_psi[,,3] <- exp(
  cbind(mc$b0[,1],mc$b1[,1]) %*% pred_mat
)
# Species 3 (red squirrel)
pred_psi[,,4] <- exp(
  cbind(mc$c0[,1],mc$c1[,1]) %*% pred_mat
)
# species 1 & 2
pred_psi[,,5] <- exp(
  cbind(mc$a0[,1],mc$a1[,1]) %*% pred_mat + 
  cbind(mc$b0[,1],mc$b1[,1]) %*% pred_mat + 
  cbind(mc$d0[,1],mc$d1[,1]) %*% pred_mat
)
# species 1 & 3
pred_psi[,,6] <- exp(
  cbind(mc$a0[,1],mc$a1[,1]) %*% pred_mat + 
  cbind(mc$c0[,1],mc$c1[,1]) %*% pred_mat + 
  cbind(mc$e0[,1],mc$e1[,1]) %*% pred_mat
)
# species 2 & 3
pred_psi[,,7] <- exp(
  cbind(mc$b0[,1],mc$b1[,1]) %*% pred_mat + 
  cbind(mc$c0[,1],mc$c1[,1]) %*% pred_mat + 
  cbind(mc$g0[,1],mc$g1[,1]) %*% pred_mat
)
# all species together
pred_psi[,,8] <- exp(
  cbind(mc$a0[,1],mc$a1[,1]) %*% pred_mat + 
    cbind(mc$b0[,1],mc$b1[,1]) %*% pred_mat + 
    cbind(mc$c0[,1],mc$c1[,1]) %*% pred_mat + 
    cbind(mc$d0[,1],mc$d1[,1]) %*% pred_mat +
    cbind(mc$e0[,1],mc$e1[,1]) %*% pred_mat +
    cbind(mc$g0[,1],mc$g1[,1]) %*% pred_mat
)
# convert to probability
prob_psi <- array(
  NA,
  dim = dim(pred_psi)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi[,i,] <- sweep(
    pred_psi[,i,],
    1,
    rowSums(pred_psi[,i,]),
    FUN = "/"
  )
}
# now we'll do the same thing, but this time including the autologistic term, phi
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept + autologistic term (remember to use 
#                           average values for all covariates except your variable of interest)
# create the matrix of the coefficients from the model
pred_psi_auto <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat_auto),
    8
  )
)
# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)

# fill it in. No species is easy, it's just a 1.
pred_psi_auto[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi_auto[,,2] <- exp(
  cbind(mc$a0[,1],mc$a1[,1],mc$phi[,1,1]) %*% pred_mat_auto
)
# Species 2 (gray squirrel)
pred_psi_auto[,,3] <- exp(
  cbind(mc$b0[,1],mc$b1[,1],mc$phi[,2,1]) %*% pred_mat_auto
)
# Species 3 (red squirrel)
pred_psi_auto[,,4] <- exp(
  cbind(mc$c0[,1],mc$c1[,1],mc$phi[,3,1]) %*% pred_mat_auto
)
# species 1 & 2
pred_psi_auto[,,5] <- exp(
  cbind(mc$a0[,1],mc$a1[,1],mc$phi[,1,1]) %*% pred_mat_auto + 
    cbind(mc$b0[,1],mc$b1[,1],mc$phi[,2,1]) %*% pred_mat_auto + 
    cbind(mc$d0[,1],mc$d1[,1]) %*% pred_mat
)
# species 1 & 3
pred_psi_auto[,,6] <- exp(
  cbind(mc$a0[,1],mc$a1[,1],mc$phi[,1,1]) %*% pred_mat_auto + 
    cbind(mc$c0[,1],mc$c1[,1],mc$phi[,3,1]) %*% pred_mat_auto + 
    cbind(mc$e0[,1],mc$e1[,1]) %*% pred_mat
)
# species 2 & 3
pred_psi_auto[,,7] <- exp(
  cbind(mc$b0[,1],mc$b1[,1],mc$phi[,2,1]) %*% pred_mat_auto + 
    cbind(mc$c0[,1],mc$c1[,1],mc$phi[,3,1]) %*% pred_mat_auto + 
    cbind(mc$g0[,1],mc$g1[,1]) %*% pred_mat
)
# all species together
pred_psi_auto[,,8] <- exp(
  cbind(mc$a0[,1],mc$a1[,1],mc$phi[,1,1]) %*% pred_mat_auto + 
    cbind(mc$b0[,1],mc$b1[,1],mc$phi[,2,1]) %*% pred_mat_auto + 
    cbind(mc$c0[,1],mc$c1[,1],mc$phi[,3,1]) %*% pred_mat_auto + 
    cbind(mc$d0[,1],mc$d1[,1]) %*% pred_mat +
    cbind(mc$e0[,1],mc$e1[,1]) %*% pred_mat +
    cbind(mc$g0[,1],mc$g1[,1]) %*% pred_mat
)
# convert to probability
prob_psi_auto <- array(
  NA,
  dim = dim(pred_psi_auto)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi_auto[,i,] <- sweep(
    pred_psi_auto[,i,],
    1,
    rowSums(pred_psi_auto[,i,]),
    FUN = "/"
  )
}

# calculating average psi
trueProb <- prob_psi / (
  prob_psi + (1 - prob_psi_auto)
)

AUTX <- apply(
  trueProb,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)
#### Chicago, IL ####
# generate a sequence of impervious surface values
range(cooc_covs$gmcCanopy[27:130])
# [1] -9.54655 47.73140
# I need to choose some 'pretty' numbers based on that
predV <- seq(-10, 48, length.out = 200)
# Generating a matrix that will contain values for predicting. This includes our sequence of
# impervious surface values, and adding a row of 1's for the intercept
pred_mat <- cbind(
  1,
  scale(predV)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  scale(predV),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_psi <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)

# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)
# fill it in. No species is easy, it's just a 1.
pred_psi[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi[,,2] <- exp(
  cbind(mc$a0[,2],mc$a1[,2]) %*% pred_mat
)
# Species 2 (gray squirrel)
pred_psi[,,3] <- exp(
  cbind(mc$b0[,2],mc$b1[,2]) %*% pred_mat
)
# Species 3 (red squirrel)
pred_psi[,,4] <- exp(
  cbind(mc$c0[,2],mc$c1[,2]) %*% pred_mat
)
# species 1 & 2
pred_psi[,,5] <- exp(
  cbind(mc$a0[,2],mc$a1[,2]) %*% pred_mat + 
    cbind(mc$b0[,2],mc$b1[,2]) %*% pred_mat + 
    cbind(mc$d0[,2],mc$d1[,2]) %*% pred_mat
)
# species 1 & 3
pred_psi[,,6] <- exp(
  cbind(mc$a0[,2],mc$a1[,2]) %*% pred_mat + 
    cbind(mc$c0[,2],mc$c1[,2]) %*% pred_mat + 
    cbind(mc$e0[,2],mc$e1[,2]) %*% pred_mat
)
# species 2 & 3
pred_psi[,,7] <- exp(
  cbind(mc$b0[,2],mc$b1[,2]) %*% pred_mat + 
    cbind(mc$c0[,2],mc$c1[,2]) %*% pred_mat + 
    cbind(mc$g0[,2],mc$g1[,2]) %*% pred_mat
)
# all species together
pred_psi[,,8] <- exp(
  cbind(mc$a0[,2],mc$a1[,2]) %*% pred_mat + 
    cbind(mc$b0[,2],mc$b1[,2]) %*% pred_mat + 
    cbind(mc$c0[,2],mc$c1[,2]) %*% pred_mat + 
    cbind(mc$d0[,2],mc$d1[,2]) %*% pred_mat +
    cbind(mc$e0[,2],mc$e1[,2]) %*% pred_mat +
    cbind(mc$g0[,2],mc$g1[,2]) %*% pred_mat
)
# convert to probability
prob_psi <- array(
  NA,
  dim = dim(pred_psi)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi[,i,] <- sweep(
    pred_psi[,i,],
    1,
    rowSums(pred_psi[,i,]),
    FUN = "/"
  )
}

# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept + autologistic term (remember to use 
#                           average values for all covariates except your variable of interest)
# create the matrix of the coefficients from the model
pred_psi_auto <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat_auto),
    8
  )
)
# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)

# fill it in. No species is easy, it's just a 1.
pred_psi_auto[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi_auto[,,2] <- exp(
  cbind(mc$a0[,2],mc$a1[,2],mc$phi[,1,2]) %*% pred_mat_auto
)
# Species 2 (gray squirrel)
pred_psi_auto[,,3] <- exp(
  cbind(mc$b0[,2],mc$b1[,2],mc$phi[,2,2]) %*% pred_mat_auto
)
# Species 3 (red squirrel)
pred_psi_auto[,,4] <- exp(
  cbind(mc$c0[,2],mc$c1[,2],mc$phi[,3,2]) %*% pred_mat_auto
)
# species 1 & 2
pred_psi_auto[,,5] <- exp(
  cbind(mc$a0[,2],mc$a1[,2],mc$phi[,1,2]) %*% pred_mat_auto + 
    cbind(mc$b0[,2],mc$b1[,2],mc$phi[,2,2]) %*% pred_mat_auto + 
    cbind(mc$d0[,2],mc$d1[,2]) %*% pred_mat
)
# species 1 & 3
pred_psi_auto[,,6] <- exp(
  cbind(mc$a0[,2],mc$a1[,2],mc$phi[,1,2]) %*% pred_mat_auto + 
    cbind(mc$c0[,2],mc$c1[,2],mc$phi[,3,2]) %*% pred_mat_auto + 
    cbind(mc$e0[,2],mc$e1[,2]) %*% pred_mat
)
# species 2 & 3
pred_psi_auto[,,7] <- exp(
  cbind(mc$b0[,2],mc$b1[,2],mc$phi[,2,2]) %*% pred_mat_auto + 
    cbind(mc$c0[,2],mc$c1[,2],mc$phi[,3,2]) %*% pred_mat_auto + 
    cbind(mc$g0[,2],mc$g1[,2]) %*% pred_mat
)
# all species together
pred_psi_auto[,,8] <- exp(
  cbind(mc$a0[,2],mc$a1[,2],mc$phi[,1,2]) %*% pred_mat_auto + 
    cbind(mc$b0[,2],mc$b1[,2],mc$phi[,2,2]) %*% pred_mat_auto + 
    cbind(mc$c0[,2],mc$c1[,2],mc$phi[,3,2]) %*% pred_mat_auto + 
    cbind(mc$d0[,2],mc$d1[,2]) %*% pred_mat +
    cbind(mc$e0[,2],mc$e1[,2]) %*% pred_mat +
    cbind(mc$g0[,2],mc$g1[,2]) %*% pred_mat
)
# convert to probability
prob_psi_auto <- array(
  NA,
  dim = dim(pred_psi_auto)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi_auto[,i,] <- sweep(
    pred_psi_auto[,i,],
    1,
    rowSums(pred_psi_auto[,i,]),
    FUN = "/"
  )
}

# calculating average psi
trueProb <- prob_psi / (
  prob_psi + (1 - prob_psi_auto)
)

CHIL <- apply(
  trueProb,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)
#### Denver, CO ####
# generate a sequence of impervious surface values
range(cooc_covs$gmcCanopy[131:169])
# [1] -3.746837 15.052608
# I need to choose some 'pretty' numbers based on that
predV <- seq(-4, 16, length.out = 200)
# Generating a matrix that will contain values for predicting. This includes our sequence of
# impervious surface values, and adding a row of 1's for the intercept
pred_mat <- cbind(
  1,
  scale(predV)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  scale(predV),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_psi <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)

# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)
# fill it in. No species is easy, it's just a 1.
pred_psi[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi[,,2] <- exp(
  cbind(mc$a0[,3],mc$a1[,3]) %*% pred_mat
)
# Species 3 (gray squirrel)
pred_psi[,,3] <- exp(
  cbind(mc$b0[,3],mc$b1[,3]) %*% pred_mat
)
# Species 3 (red squirrel)
pred_psi[,,4] <- exp(
  cbind(mc$c0[,3],mc$c1[,3]) %*% pred_mat
)
# species 1 & 2
pred_psi[,,5] <- exp(
  cbind(mc$a0[,3],mc$a1[,3]) %*% pred_mat + 
    cbind(mc$b0[,3],mc$b1[,3]) %*% pred_mat + 
    cbind(mc$d0[,3],mc$d1[,3]) %*% pred_mat
)
# species 1 & 3
pred_psi[,,6] <- exp(
  cbind(mc$a0[,3],mc$a1[,3]) %*% pred_mat + 
    cbind(mc$c0[,3],mc$c1[,3]) %*% pred_mat + 
    cbind(mc$e0[,3],mc$e1[,3]) %*% pred_mat
)
# species 2 & 3
pred_psi[,,7] <- exp(
  cbind(mc$b0[,3],mc$b1[,3]) %*% pred_mat + 
    cbind(mc$c0[,3],mc$c1[,3]) %*% pred_mat + 
    cbind(mc$g0[,3],mc$g1[,3]) %*% pred_mat
)
# all species together
pred_psi[,,8] <- exp(
  cbind(mc$a0[,3],mc$a1[,3]) %*% pred_mat + 
    cbind(mc$b0[,3],mc$b1[,3]) %*% pred_mat + 
    cbind(mc$c0[,3],mc$c1[,3]) %*% pred_mat + 
    cbind(mc$d0[,3],mc$d1[,3]) %*% pred_mat +
    cbind(mc$e0[,3],mc$e1[,3]) %*% pred_mat +
    cbind(mc$g0[,3],mc$g1[,3]) %*% pred_mat
)
# convert to probability
prob_psi <- array(
  NA,
  dim = dim(pred_psi)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi[,i,] <- sweep(
    pred_psi[,i,],
    1,
    rowSums(pred_psi[,i,]),
    FUN = "/"
  )
}

# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept + autologistic term (remember to use 
#                           average values for all covariates except your variable of interest)
# create the matrix of the coefficients from the model
pred_psi_auto <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat_auto),
    8
  )
)
# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)

# fill it in. No species is easy, it's just a 1.
pred_psi_auto[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi_auto[,,2] <- exp(
  cbind(mc$a0[,3],mc$a1[,3],mc$phi[,1,3]) %*% pred_mat_auto
)
# Species 2 (gray squirrel)
pred_psi_auto[,,3] <- exp(
  cbind(mc$b0[,3],mc$b1[,3],mc$phi[,2,3]) %*% pred_mat_auto
)
# Species 3 (red squirrel)
pred_psi_auto[,,4] <- exp(
  cbind(mc$c0[,3],mc$c1[,3],mc$phi[,3,3]) %*% pred_mat_auto
)
# species 1 & 2
pred_psi_auto[,,5] <- exp(
  cbind(mc$a0[,3],mc$a1[,3],mc$phi[,1,3]) %*% pred_mat_auto + 
    cbind(mc$b0[,3],mc$b1[,3],mc$phi[,2,3]) %*% pred_mat_auto + 
    cbind(mc$d0[,3],mc$d1[,3]) %*% pred_mat
)
# species 1 & 3
pred_psi_auto[,,6] <- exp(
  cbind(mc$a0[,3],mc$a1[,3],mc$phi[,1,3]) %*% pred_mat_auto + 
    cbind(mc$c0[,3],mc$c1[,3],mc$phi[,3,3]) %*% pred_mat_auto + 
    cbind(mc$e0[,3],mc$e1[,3]) %*% pred_mat
)
# species 2 & 3
pred_psi_auto[,,7] <- exp(
  cbind(mc$b0[,3],mc$b1[,3],mc$phi[,2,3]) %*% pred_mat_auto + 
    cbind(mc$c0[,3],mc$c1[,3],mc$phi[,3,3]) %*% pred_mat_auto + 
    cbind(mc$g0[,3],mc$g1[,3]) %*% pred_mat
)
# all species together
pred_psi_auto[,,8] <- exp(
  cbind(mc$a0[,3],mc$a1[,3],mc$phi[,1,3]) %*% pred_mat_auto + 
    cbind(mc$b0[,3],mc$b1[,3],mc$phi[,2,3]) %*% pred_mat_auto + 
    cbind(mc$c0[,3],mc$c1[,3],mc$phi[,3,3]) %*% pred_mat_auto + 
    cbind(mc$d0[,3],mc$d1[,3]) %*% pred_mat +
    cbind(mc$e0[,3],mc$e1[,3]) %*% pred_mat +
    cbind(mc$g0[,3],mc$g1[,3]) %*% pred_mat
)
# convert to probability
prob_psi_auto <- array(
  NA,
  dim = dim(pred_psi_auto)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi_auto[,i,] <- sweep(
    pred_psi_auto[,i,],
    1,
    rowSums(pred_psi_auto[,i,]),
    FUN = "/"
  )
}

# calculating average psi
trueProb <- prob_psi / (
  prob_psi + (1 - prob_psi_auto)
)

DECO <- apply(
  trueProb,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)
#### Indianapolis, IN ####
# generate a sequence of impervious surface values
range(cooc_covs$gmcCanopy[170:210])
# [1] -10.47935  23.83548
# I need to choose some 'pretty' numbers based on that
predV <- seq(-11, 24, length.out = 200)
# Generating a matrix that will contain values for predicting. This includes our sequence of
# impervious surface values, and adding a row of 1's for the intercept
pred_mat <- cbind(
  1,
  scale(predV)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  scale(predV),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_psi <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)

# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)
# fill it in. No species is easy, it's just a 1.
pred_psi[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi[,,2] <- exp(
  cbind(mc$a0[,4],mc$a1[,4]) %*% pred_mat
)
# Species 3 (gray squirrel)
pred_psi[,,3] <- exp(
  cbind(mc$b0[,4],mc$b1[,4]) %*% pred_mat
)
# Species 3 (red squirrel)
pred_psi[,,4] <- exp(
  cbind(mc$c0[,4],mc$c1[,4]) %*% pred_mat
)
# species 1 & 2
pred_psi[,,5] <- exp(
  cbind(mc$a0[,4],mc$a1[,4]) %*% pred_mat + 
    cbind(mc$b0[,4],mc$b1[,4]) %*% pred_mat + 
    cbind(mc$d0[,4],mc$d1[,4]) %*% pred_mat
)
# species 1 & 3
pred_psi[,,6] <- exp(
  cbind(mc$a0[,4],mc$a1[,4]) %*% pred_mat + 
    cbind(mc$c0[,4],mc$c1[,4]) %*% pred_mat + 
    cbind(mc$e0[,4],mc$e1[,4]) %*% pred_mat
)
# species 2 & 3
pred_psi[,,7] <- exp(
  cbind(mc$b0[,4],mc$b1[,4]) %*% pred_mat + 
    cbind(mc$c0[,4],mc$c1[,4]) %*% pred_mat + 
    cbind(mc$g0[,4],mc$g1[,4]) %*% pred_mat
)
# all species together
pred_psi[,,8] <- exp(
  cbind(mc$a0[,4],mc$a1[,4]) %*% pred_mat + 
    cbind(mc$b0[,4],mc$b1[,4]) %*% pred_mat + 
    cbind(mc$c0[,4],mc$c1[,4]) %*% pred_mat + 
    cbind(mc$d0[,4],mc$d1[,4]) %*% pred_mat +
    cbind(mc$e0[,4],mc$e1[,4]) %*% pred_mat +
    cbind(mc$g0[,4],mc$g1[,4]) %*% pred_mat
)
# convert to probability
prob_psi <- array(
  NA,
  dim = dim(pred_psi)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi[,i,] <- sweep(
    pred_psi[,i,],
    1,
    rowSums(pred_psi[,i,]),
    FUN = "/"
  )
}

# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept + autologistic term (remember to use 
#                           average values for all covariates except your variable of interest)
# create the matrix of the coefficients from the model
pred_psi_auto <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)
# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)

# fill it in. No species is easy, it's just a 1.
pred_psi_auto[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi_auto[,,2] <- exp(
  cbind(mc$a0[,4],mc$a1[,4],mc$phi[,1,4]) %*% pred_mat_auto
)
# Species 2 (gray squirrel)
pred_psi_auto[,,3] <- exp(
  cbind(mc$b0[,4],mc$b1[,4],mc$phi[,2,4]) %*% pred_mat_auto
)
# Species 3 (red squirrel)
pred_psi_auto[,,4] <- exp(
  cbind(mc$c0[,4],mc$c1[,4],mc$phi[,3,4]) %*% pred_mat_auto
)
# species 1 & 2
pred_psi_auto[,,5] <- exp(
  cbind(mc$a0[,4],mc$a1[,4],mc$phi[,1,4]) %*% pred_mat_auto + 
    cbind(mc$b0[,4],mc$b1[,4],mc$phi[,2,4]) %*% pred_mat_auto + 
    cbind(mc$d0[,4],mc$d1[,4]) %*% pred_mat
)
# species 1 & 3
pred_psi_auto[,,6] <- exp(
  cbind(mc$a0[,4],mc$a1[,4],mc$phi[,1,4]) %*% pred_mat_auto + 
    cbind(mc$c0[,4],mc$c1[,4],mc$phi[,3,4]) %*% pred_mat_auto + 
    cbind(mc$e0[,4],mc$e1[,4]) %*% pred_mat
)
# species 2 & 3
pred_psi_auto[,,7] <- exp(
  cbind(mc$b0[,4],mc$b1[,4],mc$phi[,2,4]) %*% pred_mat_auto + 
    cbind(mc$c0[,4],mc$c1[,4],mc$phi[,3,4]) %*% pred_mat_auto + 
    cbind(mc$g0[,4],mc$g1[,4]) %*% pred_mat
)
# all species together
pred_psi_auto[,,8] <- exp(
  cbind(mc$a0[,4],mc$a1[,4],mc$phi[,1,4]) %*% pred_mat_auto + 
    cbind(mc$b0[,4],mc$b1[,4],mc$phi[,2,4]) %*% pred_mat_auto + 
    cbind(mc$c0[,4],mc$c1[,4],mc$phi[,3,4]) %*% pred_mat_auto + 
    cbind(mc$d0[,4],mc$d1[,4]) %*% pred_mat +
    cbind(mc$e0[,4],mc$e1[,4]) %*% pred_mat +
    cbind(mc$g0[,4],mc$g1[,4]) %*% pred_mat
)
# convert to probability
prob_psi_auto <- array(
  NA,
  dim = dim(pred_psi_auto)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi_auto[,i,] <- sweep(
    pred_psi_auto[,i,],
    1,
    rowSums(pred_psi_auto[,i,]),
    FUN = "/"
  )
}

# calculating average psi
trueProb <- prob_psi / (
  prob_psi + (1 - prob_psi_auto)
)

ININ <- apply(
  trueProb,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)
#### Iowa City, IA ####
# generate a sequence of impervious surface values
range(cooc_covs$gmcCanopy[211:249])
# [1] -8.91934 46.20841
# I need to choose some 'pretty' numbers based on that
predV <- seq(-9, 47, length.out = 200)
# Generating a matrix that will contain values for predicting. This includes our sequence of
# impervious surface values, and adding a row of 1's for the intercept
pred_mat <- cbind(
  1,
  scale(predV)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  scale(predV),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_psi <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)

# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)
# fill it in. No species is easy, it's just a 1.
pred_psi[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi[,,2] <- exp(
  cbind(mc$a0[,5],mc$a1[,5]) %*% pred_mat
)
# Species 3 (gray squirrel)
pred_psi[,,3] <- exp(
  cbind(mc$b0[,5],mc$b1[,5]) %*% pred_mat
)
# Species 3 (red squirrel)
pred_psi[,,4] <- exp(
  cbind(mc$c0[,5],mc$c1[,5]) %*% pred_mat
)
# species 1 & 2
pred_psi[,,5] <- exp(
  cbind(mc$a0[,5],mc$a1[,5]) %*% pred_mat + 
    cbind(mc$b0[,5],mc$b1[,5]) %*% pred_mat + 
    cbind(mc$d0[,5],mc$d1[,5]) %*% pred_mat
)
# species 1 & 3
pred_psi[,,6] <- exp(
  cbind(mc$a0[,5],mc$a1[,5]) %*% pred_mat + 
    cbind(mc$c0[,5],mc$c1[,5]) %*% pred_mat + 
    cbind(mc$e0[,5],mc$e1[,5]) %*% pred_mat
)
# species 2 & 3
pred_psi[,,7] <- exp(
  cbind(mc$b0[,5],mc$b1[,5]) %*% pred_mat + 
    cbind(mc$c0[,5],mc$c1[,5]) %*% pred_mat + 
    cbind(mc$g0[,5],mc$g1[,5]) %*% pred_mat
)
# all species together
pred_psi[,,8] <- exp(
  cbind(mc$a0[,5],mc$a1[,5]) %*% pred_mat + 
    cbind(mc$b0[,5],mc$b1[,5]) %*% pred_mat + 
    cbind(mc$c0[,5],mc$c1[,5]) %*% pred_mat + 
    cbind(mc$d0[,5],mc$d1[,5]) %*% pred_mat +
    cbind(mc$e0[,5],mc$e1[,5]) %*% pred_mat +
    cbind(mc$g0[,5],mc$g1[,5]) %*% pred_mat
)
# convert to probability
prob_psi <- array(
  NA,
  dim = dim(pred_psi)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi[,i,] <- sweep(
    pred_psi[,i,],
    1,
    rowSums(pred_psi[,i,]),
    FUN = "/"
  )
}

# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept + autologistic term (remember to use 
#                           average values for all covariates except your variable of interest)
# create the matrix of the coefficients from the model
pred_psi_auto <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)
# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)

# fill it in. No species is easy, it's just a 1.
pred_psi_auto[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi_auto[,,2] <- exp(
  cbind(mc$a0[,5],mc$a1[,5],mc$phi[,1,5]) %*% pred_mat_auto
)
# Species 2 (gray squirrel)
pred_psi_auto[,,3] <- exp(
  cbind(mc$b0[,5],mc$b1[,5],mc$phi[,2,5]) %*% pred_mat_auto
)
# Species 3 (red squirrel)
pred_psi_auto[,,4] <- exp(
  cbind(mc$c0[,5],mc$c1[,5],mc$phi[,3,5]) %*% pred_mat_auto
)
# species 1 & 2
pred_psi_auto[,,5] <- exp(
  cbind(mc$a0[,5],mc$a1[,5],mc$phi[,1,5]) %*% pred_mat_auto + 
    cbind(mc$b0[,5],mc$b1[,5],mc$phi[,2,5]) %*% pred_mat_auto + 
    cbind(mc$d0[,5],mc$d1[,5]) %*% pred_mat
)
# species 1 & 3
pred_psi_auto[,,6] <- exp(
  cbind(mc$a0[,5],mc$a1[,5],mc$phi[,1,5]) %*% pred_mat_auto + 
    cbind(mc$c0[,5],mc$c1[,5],mc$phi[,3,5]) %*% pred_mat_auto + 
    cbind(mc$e0[,5],mc$e1[,5]) %*% pred_mat
)
# species 2 & 3
pred_psi_auto[,,7] <- exp(
  cbind(mc$b0[,5],mc$b1[,5],mc$phi[,2,5]) %*% pred_mat_auto + 
    cbind(mc$c0[,5],mc$c1[,5],mc$phi[,3,5]) %*% pred_mat_auto + 
    cbind(mc$g0[,5],mc$g1[,5]) %*% pred_mat
)
# all species together
pred_psi_auto[,,8] <- exp(
  cbind(mc$a0[,5],mc$a1[,5],mc$phi[,1,5]) %*% pred_mat_auto + 
    cbind(mc$b0[,5],mc$b1[,5],mc$phi[,2,5]) %*% pred_mat_auto + 
    cbind(mc$c0[,5],mc$c1[,5],mc$phi[,3,5]) %*% pred_mat_auto + 
    cbind(mc$d0[,5],mc$d1[,5]) %*% pred_mat +
    cbind(mc$e0[,5],mc$e1[,5]) %*% pred_mat +
    cbind(mc$g0[,5],mc$g1[,5]) %*% pred_mat
)
# convert to probability
prob_psi_auto <- array(
  NA,
  dim = dim(pred_psi_auto)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi_auto[,i,] <- sweep(
    pred_psi_auto[,i,],
    1,
    rowSums(pred_psi_auto[,i,]),
    FUN = "/"
  )
}

# calculating average psi
trueProb <- prob_psi / (
  prob_psi + (1 - prob_psi_auto)
)

ICIA <- apply(
  trueProb,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)
#### Madison, WI ####
# generate a sequence of impervious surface values
range(cooc_covs$gmcCanopy[250:271])
# [1] -14.50888  33.78841
# I need to choose some 'pretty' numbers based on that
predV <- seq(-15, 34, length.out = 200)
# Generating a matrix that will contain values for predicting. This includes our sequence of
# impervious surface values, and adding a row of 1's for the intercept
pred_mat <- cbind(
  1,
  scale(predV)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  scale(predV),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_psi <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)

# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)
# fill it in. No species is easy, it's just a 1.
pred_psi[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi[,,2] <- exp(
  cbind(mc$a0[,6],mc$a1[,6]) %*% pred_mat
)
# Species 3 (gray squirrel)
pred_psi[,,3] <- exp(
  cbind(mc$b0[,6],mc$b1[,6]) %*% pred_mat
)
# Species 3 (red squirrel)
pred_psi[,,4] <- exp(
  cbind(mc$c0[,6],mc$c1[,6]) %*% pred_mat
)
# species 1 & 2
pred_psi[,,5] <- exp(
  cbind(mc$a0[,6],mc$a1[,6]) %*% pred_mat + 
    cbind(mc$b0[,6],mc$b1[,6]) %*% pred_mat + 
    cbind(mc$d0[,6],mc$d1[,6]) %*% pred_mat
)
# species 1 & 3
pred_psi[,,6] <- exp(
  cbind(mc$a0[,6],mc$a1[,6]) %*% pred_mat + 
    cbind(mc$c0[,6],mc$c1[,6]) %*% pred_mat + 
    cbind(mc$e0[,6],mc$e1[,6]) %*% pred_mat
)
# species 2 & 3
pred_psi[,,7] <- exp(
  cbind(mc$b0[,6],mc$b1[,6]) %*% pred_mat + 
    cbind(mc$c0[,6],mc$c1[,6]) %*% pred_mat + 
    cbind(mc$g0[,6],mc$g1[,6]) %*% pred_mat
)
# all species together
pred_psi[,,8] <- exp(
  cbind(mc$a0[,6],mc$a1[,6]) %*% pred_mat + 
    cbind(mc$b0[,6],mc$b1[,6]) %*% pred_mat + 
    cbind(mc$c0[,6],mc$c1[,6]) %*% pred_mat + 
    cbind(mc$d0[,6],mc$d1[,6]) %*% pred_mat +
    cbind(mc$e0[,6],mc$e1[,6]) %*% pred_mat +
    cbind(mc$g0[,6],mc$g1[,6]) %*% pred_mat
)
# convert to probability
prob_psi <- array(
  NA,
  dim = dim(pred_psi)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi[,i,] <- sweep(
    pred_psi[,i,],
    1,
    rowSums(pred_psi[,i,]),
    FUN = "/"
  )
}

# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept + autologistic term (remember to use 
#                           average values for all covariates except your variable of interest)
# create the matrix of the coefficients from the model
pred_psi_auto <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)
# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)

# fill it in. No species is easy, it's just a 1.
pred_psi_auto[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi_auto[,,2] <- exp(
  cbind(mc$a0[,6],mc$a1[,6],mc$phi[,1,6]) %*% pred_mat_auto
)
# Species 2 (gray squirrel)
pred_psi_auto[,,3] <- exp(
  cbind(mc$b0[,6],mc$b1[,6],mc$phi[,2,6]) %*% pred_mat_auto
)
# Species 3 (red squirrel)
pred_psi_auto[,,4] <- exp(
  cbind(mc$c0[,6],mc$c1[,6],mc$phi[,3,6]) %*% pred_mat_auto
)
# species 1 & 2
pred_psi_auto[,,5] <- exp(
  cbind(mc$a0[,6],mc$a1[,6],mc$phi[,1,6]) %*% pred_mat_auto + 
    cbind(mc$b0[,6],mc$b1[,6],mc$phi[,2,6]) %*% pred_mat_auto + 
    cbind(mc$d0[,6],mc$d1[,6]) %*% pred_mat
)
# species 1 & 3
pred_psi_auto[,,6] <- exp(
  cbind(mc$a0[,6],mc$a1[,6],mc$phi[,1,6]) %*% pred_mat_auto + 
    cbind(mc$c0[,6],mc$c1[,6],mc$phi[,3,6]) %*% pred_mat_auto + 
    cbind(mc$e0[,6],mc$e1[,6]) %*% pred_mat
)
# species 2 & 3
pred_psi_auto[,,7] <- exp(
  cbind(mc$b0[,6],mc$b1[,6],mc$phi[,2,6]) %*% pred_mat_auto + 
    cbind(mc$c0[,6],mc$c1[,6],mc$phi[,3,6]) %*% pred_mat_auto + 
    cbind(mc$g0[,6],mc$g1[,6]) %*% pred_mat
)
# all species together
pred_psi_auto[,,8] <- exp(
  cbind(mc$a0[,6],mc$a1[,6],mc$phi[,1,6]) %*% pred_mat_auto + 
    cbind(mc$b0[,6],mc$b1[,6],mc$phi[,2,6]) %*% pred_mat_auto + 
    cbind(mc$c0[,6],mc$c1[,6],mc$phi[,3,6]) %*% pred_mat_auto + 
    cbind(mc$d0[,6],mc$d1[,6]) %*% pred_mat +
    cbind(mc$e0[,6],mc$e1[,6]) %*% pred_mat +
    cbind(mc$g0[,6],mc$g1[,6]) %*% pred_mat
)
# convert to probability
prob_psi_auto <- array(
  NA,
  dim = dim(pred_psi_auto)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi_auto[,i,] <- sweep(
    pred_psi_auto[,i,],
    1,
    rowSums(pred_psi_auto[,i,]),
    FUN = "/"
  )
}

# calculating average psi
trueProb <- prob_psi / (
  prob_psi + (1 - prob_psi_auto)
)

MAWI <- apply(
  trueProb,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)
#### Metro LA ####
# generate a sequence of impervious surface values
range(cooc_covs$gmcCanopy[272:314])
# [1] -3.575935 19.745199
# I need to choose some 'pretty' numbers based on that
predV <- seq(-4, 20, length.out = 200)
# Generating a matrix that will contain values for predicting. This includes our sequence of
# impervious surface values, and adding a row of 1's for the intercept
pred_mat <- cbind(
  1,
  scale(predV)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  scale(predV),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_psi <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)

# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)
# fill it in. No species is easy, it's just a 1.
pred_psi[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi[,,2] <- exp(
  cbind(mc$a0[,7],mc$a1[,7]) %*% pred_mat
)
# Species 3 (gray squirrel)
pred_psi[,,3] <- exp(
  cbind(mc$b0[,7],mc$b1[,7]) %*% pred_mat
)
# Species 3 (red squirrel)
pred_psi[,,4] <- exp(
  cbind(mc$c0[,7],mc$c1[,7]) %*% pred_mat
)
# species 1 & 2
pred_psi[,,5] <- exp(
  cbind(mc$a0[,7],mc$a1[,7]) %*% pred_mat + 
    cbind(mc$b0[,7],mc$b1[,7]) %*% pred_mat + 
    cbind(mc$d0[,7],mc$d1[,7]) %*% pred_mat
)
# species 1 & 3
pred_psi[,,6] <- exp(
  cbind(mc$a0[,7],mc$a1[,7]) %*% pred_mat + 
    cbind(mc$c0[,7],mc$c1[,7]) %*% pred_mat + 
    cbind(mc$e0[,7],mc$e1[,7]) %*% pred_mat
)
# species 2 & 3
pred_psi[,,7] <- exp(
  cbind(mc$b0[,7],mc$b1[,7]) %*% pred_mat + 
    cbind(mc$c0[,7],mc$c1[,7]) %*% pred_mat + 
    cbind(mc$g0[,7],mc$g1[,7]) %*% pred_mat
)
# all species together
pred_psi[,,8] <- exp(
  cbind(mc$a0[,7],mc$a1[,7]) %*% pred_mat + 
    cbind(mc$b0[,7],mc$b1[,7]) %*% pred_mat + 
    cbind(mc$c0[,7],mc$c1[,7]) %*% pred_mat + 
    cbind(mc$d0[,7],mc$d1[,7]) %*% pred_mat +
    cbind(mc$e0[,7],mc$e1[,7]) %*% pred_mat +
    cbind(mc$g0[,7],mc$g1[,7]) %*% pred_mat
)
# convert to probability
prob_psi <- array(
  NA,
  dim = dim(pred_psi)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi[,i,] <- sweep(
    pred_psi[,i,],
    1,
    rowSums(pred_psi[,i,]),
    FUN = "/"
  )
}

# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept + autologistic term (remember to use 
#                           average values for all covariates except your variable of interest)
# create the matrix of the coefficients from the model
pred_psi_auto <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)
# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)

# fill it in. No species is easy, it's just a 1.
pred_psi_auto[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi_auto[,,2] <- exp(
  cbind(mc$a0[,7],mc$a1[,7],mc$phi[,1,7]) %*% pred_mat_auto
)
# Species 2 (gray squirrel)
pred_psi_auto[,,3] <- exp(
  cbind(mc$b0[,7],mc$b1[,7],mc$phi[,2,7]) %*% pred_mat_auto
)
# Species 3 (red squirrel)
pred_psi_auto[,,4] <- exp(
  cbind(mc$c0[,7],mc$c1[,7],mc$phi[,3,7]) %*% pred_mat_auto
)
# species 1 & 2
pred_psi_auto[,,5] <- exp(
  cbind(mc$a0[,7],mc$a1[,7],mc$phi[,1,7]) %*% pred_mat_auto + 
    cbind(mc$b0[,7],mc$b1[,7],mc$phi[,2,7]) %*% pred_mat_auto + 
    cbind(mc$d0[,7],mc$d1[,7]) %*% pred_mat
)
# species 1 & 3
pred_psi_auto[,,6] <- exp(
  cbind(mc$a0[,7],mc$a1[,7],mc$phi[,1,7]) %*% pred_mat_auto + 
    cbind(mc$c0[,7],mc$c1[,7],mc$phi[,3,7]) %*% pred_mat_auto + 
    cbind(mc$e0[,7],mc$e1[,7]) %*% pred_mat
)
# species 2 & 3
pred_psi_auto[,,7] <- exp(
  cbind(mc$b0[,7],mc$b1[,7],mc$phi[,2,7]) %*% pred_mat_auto + 
    cbind(mc$c0[,7],mc$c1[,7],mc$phi[,3,7]) %*% pred_mat_auto + 
    cbind(mc$g0[,7],mc$g1[,7]) %*% pred_mat
)
# all species together
pred_psi_auto[,,8] <- exp(
  cbind(mc$a0[,7],mc$a1[,7],mc$phi[,1,7]) %*% pred_mat_auto + 
    cbind(mc$b0[,7],mc$b1[,7],mc$phi[,2,7]) %*% pred_mat_auto + 
    cbind(mc$c0[,7],mc$c1[,7],mc$phi[,3,7]) %*% pred_mat_auto + 
    cbind(mc$d0[,7],mc$d1[,7]) %*% pred_mat +
    cbind(mc$e0[,7],mc$e1[,7]) %*% pred_mat +
    cbind(mc$g0[,7],mc$g1[,7]) %*% pred_mat
)
# convert to probability
prob_psi_auto <- array(
  NA,
  dim = dim(pred_psi_auto)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi_auto[,i,] <- sweep(
    pred_psi_auto[,i,],
    1,
    rowSums(pred_psi_auto[,i,]),
    FUN = "/"
  )
}

# calculating average psi
trueProb <- prob_psi / (
  prob_psi + (1 - prob_psi_auto)
)

MELA <- apply(
  trueProb,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)
#### Rochester, NY ####
# generate a sequence of impervious surface values
range(cooc_covs$gmcCanopy[315:337])
# [1] -24.86679  48.39943
# I need to choose some 'pretty' numbers based on that
predV <- seq(-25, 49, length.out = 200)
# Generating a matrix that will contain values for predicting. This includes our sequence of
# impervious surface values, and adding a row of 1's for the intercept
pred_mat <- cbind(
  1,
  scale(predV)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  scale(predV),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_psi <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)

# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)
# fill it in. No species is easy, it's just a 1.
pred_psi[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi[,,2] <- exp(
  cbind(mc$a0[,8],mc$a1[,8]) %*% pred_mat
)
# Species 3 (gray squirrel)
pred_psi[,,3] <- exp(
  cbind(mc$b0[,8],mc$b1[,8]) %*% pred_mat
)
# Species 3 (red squirrel)
pred_psi[,,4] <- exp(
  cbind(mc$c0[,8],mc$c1[,8]) %*% pred_mat
)
# species 1 & 2
pred_psi[,,5] <- exp(
  cbind(mc$a0[,8],mc$a1[,8]) %*% pred_mat + 
    cbind(mc$b0[,8],mc$b1[,8]) %*% pred_mat + 
    cbind(mc$d0[,8],mc$d1[,8]) %*% pred_mat
)
# species 1 & 3
pred_psi[,,6] <- exp(
  cbind(mc$a0[,8],mc$a1[,8]) %*% pred_mat + 
    cbind(mc$c0[,8],mc$c1[,8]) %*% pred_mat + 
    cbind(mc$e0[,8],mc$e1[,8]) %*% pred_mat
)
# species 2 & 3
pred_psi[,,7] <- exp(
  cbind(mc$b0[,8],mc$b1[,8]) %*% pred_mat + 
    cbind(mc$c0[,8],mc$c1[,8]) %*% pred_mat + 
    cbind(mc$g0[,8],mc$g1[,8]) %*% pred_mat
)
# all species together
pred_psi[,,8] <- exp(
  cbind(mc$a0[,8],mc$a1[,8]) %*% pred_mat + 
    cbind(mc$b0[,8],mc$b1[,8]) %*% pred_mat + 
    cbind(mc$c0[,8],mc$c1[,8]) %*% pred_mat + 
    cbind(mc$d0[,8],mc$d1[,8]) %*% pred_mat +
    cbind(mc$e0[,8],mc$e1[,8]) %*% pred_mat +
    cbind(mc$g0[,8],mc$g1[,8]) %*% pred_mat
)
# convert to probability
prob_psi <- array(
  NA,
  dim = dim(pred_psi)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi[,i,] <- sweep(
    pred_psi[,i,],
    1,
    rowSums(pred_psi[,i,]),
    FUN = "/"
  )
}

# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept + autologistic term (remember to use 
#                           average values for all covariates except your variable of interest)
# create the matrix of the coefficients from the model
pred_psi_auto <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)
# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)

# fill it in. No species is easy, it's just a 1.
pred_psi_auto[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi_auto[,,2] <- exp(
  cbind(mc$a0[,8],mc$a1[,8],mc$phi[,1,8]) %*% pred_mat_auto
)
# Species 2 (gray squirrel)
pred_psi_auto[,,3] <- exp(
  cbind(mc$b0[,8],mc$b1[,8],mc$phi[,2,8]) %*% pred_mat_auto
)
# Species 3 (red squirrel)
pred_psi_auto[,,4] <- exp(
  cbind(mc$c0[,8],mc$c1[,8],mc$phi[,3,8]) %*% pred_mat_auto
)
# species 1 & 2
pred_psi_auto[,,5] <- exp(
  cbind(mc$a0[,8],mc$a1[,8],mc$phi[,1,8]) %*% pred_mat_auto + 
    cbind(mc$b0[,8],mc$b1[,8],mc$phi[,2,8]) %*% pred_mat_auto + 
    cbind(mc$d0[,8],mc$d1[,8]) %*% pred_mat
)
# species 1 & 3
pred_psi_auto[,,6] <- exp(
  cbind(mc$a0[,8],mc$a1[,8],mc$phi[,1,8]) %*% pred_mat_auto + 
    cbind(mc$c0[,8],mc$c1[,8],mc$phi[,3,8]) %*% pred_mat_auto + 
    cbind(mc$e0[,8],mc$e1[,8]) %*% pred_mat
)
# species 2 & 3
pred_psi_auto[,,7] <- exp(
  cbind(mc$b0[,8],mc$b1[,8],mc$phi[,2,8]) %*% pred_mat_auto + 
    cbind(mc$c0[,8],mc$c1[,8],mc$phi[,3,8]) %*% pred_mat_auto + 
    cbind(mc$g0[,8],mc$g1[,8]) %*% pred_mat
)
# all species together
pred_psi_auto[,,8] <- exp(
  cbind(mc$a0[,8],mc$a1[,8],mc$phi[,1,8]) %*% pred_mat_auto + 
    cbind(mc$b0[,8],mc$b1[,8],mc$phi[,2,8]) %*% pred_mat_auto + 
    cbind(mc$c0[,8],mc$c1[,8],mc$phi[,3,8]) %*% pred_mat_auto + 
    cbind(mc$d0[,8],mc$d1[,8]) %*% pred_mat +
    cbind(mc$e0[,8],mc$e1[,8]) %*% pred_mat +
    cbind(mc$g0[,8],mc$g1[,8]) %*% pred_mat
)
# convert to probability
prob_psi_auto <- array(
  NA,
  dim = dim(pred_psi_auto)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi_auto[,i,] <- sweep(
    pred_psi_auto[,i,],
    1,
    rowSums(pred_psi_auto[,i,]),
    FUN = "/"
  )
}

# calculating average psi
trueProb <- prob_psi / (
  prob_psi + (1 - prob_psi_auto)
)

RONY <- apply(
  trueProb,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)
#### Wilmington, DE ####
# generate a sequence of impervious surface values
range(cooc_covs$gmcCanopy[338:366])
# [1] -24.56086  57.18273
# I need to choose some 'pretty' numbers based on that
predV <- seq(-25, 58, length.out = 200)
# Generating a matrix that will contain values for predicting. This includes our sequence of
# impervious surface values, and adding a row of 1's for the intercept
pred_mat <- cbind(
  1,
  scale(predV)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  scale(predV),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_psi <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)

# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)
# fill it in. No species is easy, it's just a 1.
pred_psi[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi[,,2] <- exp(
  cbind(mc$a0[,9],mc$a1[,9]) %*% pred_mat
)
# Species 3 (gray squirrel)
pred_psi[,,3] <- exp(
  cbind(mc$b0[,9],mc$b1[,9]) %*% pred_mat
)
# Species 3 (red squirrel)
pred_psi[,,4] <- exp(
  cbind(mc$c0[,9],mc$c1[,9]) %*% pred_mat
)
# species 1 & 2
pred_psi[,,5] <- exp(
  cbind(mc$a0[,9],mc$a1[,9]) %*% pred_mat + 
    cbind(mc$b0[,9],mc$b1[,9]) %*% pred_mat + 
    cbind(mc$d0[,9],mc$d1[,9]) %*% pred_mat
)
# species 1 & 3
pred_psi[,,6] <- exp(
  cbind(mc$a0[,9],mc$a1[,9]) %*% pred_mat + 
    cbind(mc$c0[,9],mc$c1[,9]) %*% pred_mat + 
    cbind(mc$e0[,9],mc$e1[,9]) %*% pred_mat
)
# species 2 & 3
pred_psi[,,7] <- exp(
  cbind(mc$b0[,9],mc$b1[,9]) %*% pred_mat + 
    cbind(mc$c0[,9],mc$c1[,9]) %*% pred_mat + 
    cbind(mc$g0[,9],mc$g1[,9]) %*% pred_mat
)
# all species together
pred_psi[,,8] <- exp(
  cbind(mc$a0[,9],mc$a1[,9]) %*% pred_mat + 
    cbind(mc$b0[,9],mc$b1[,9]) %*% pred_mat + 
    cbind(mc$c0[,9],mc$c1[,9]) %*% pred_mat + 
    cbind(mc$d0[,9],mc$d1[,9]) %*% pred_mat +
    cbind(mc$e0[,9],mc$e1[,9]) %*% pred_mat +
    cbind(mc$g0[,9],mc$g1[,9]) %*% pred_mat
)
# convert to probability
prob_psi <- array(
  NA,
  dim = dim(pred_psi)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi[,i,] <- sweep(
    pred_psi[,i,],
    1,
    rowSums(pred_psi[,i,]),
    FUN = "/"
  )
}

# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept + autologistic term (remember to use 
#                           average values for all covariates except your variable of interest)
# create the matrix of the coefficients from the model
pred_psi_auto <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)
# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)

# fill it in. No species is easy, it's just a 1.
pred_psi_auto[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi_auto[,,2] <- exp(
  cbind(mc$a0[,9],mc$a1[,9],mc$phi[,1,9]) %*% pred_mat_auto
)
# Species 2 (gray squirrel)
pred_psi_auto[,,3] <- exp(
  cbind(mc$b0[,9],mc$b1[,9],mc$phi[,2,9]) %*% pred_mat_auto
)
# Species 3 (red squirrel)
pred_psi_auto[,,4] <- exp(
  cbind(mc$c0[,9],mc$c1[,9],mc$phi[,3,9]) %*% pred_mat_auto
)
# species 1 & 2
pred_psi_auto[,,5] <- exp(
  cbind(mc$a0[,9],mc$a1[,9],mc$phi[,1,9]) %*% pred_mat_auto + 
    cbind(mc$b0[,9],mc$b1[,9],mc$phi[,2,9]) %*% pred_mat_auto + 
    cbind(mc$d0[,9],mc$d1[,9]) %*% pred_mat
)
# species 1 & 3
pred_psi_auto[,,6] <- exp(
  cbind(mc$a0[,9],mc$a1[,9],mc$phi[,1,9]) %*% pred_mat_auto + 
    cbind(mc$c0[,9],mc$c1[,9],mc$phi[,3,9]) %*% pred_mat_auto + 
    cbind(mc$e0[,9],mc$e1[,9]) %*% pred_mat
)
# species 2 & 3
pred_psi_auto[,,7] <- exp(
  cbind(mc$b0[,9],mc$b1[,9],mc$phi[,2,9]) %*% pred_mat_auto + 
    cbind(mc$c0[,9],mc$c1[,9],mc$phi[,3,9]) %*% pred_mat_auto + 
    cbind(mc$g0[,9],mc$g1[,9]) %*% pred_mat
)
# all species together
pred_psi_auto[,,8] <- exp(
  cbind(mc$a0[,9],mc$a1[,9],mc$phi[,1,9]) %*% pred_mat_auto + 
    cbind(mc$b0[,9],mc$b1[,9],mc$phi[,2,9]) %*% pred_mat_auto + 
    cbind(mc$c0[,9],mc$c1[,9],mc$phi[,3,9]) %*% pred_mat_auto + 
    cbind(mc$d0[,9],mc$d1[,9]) %*% pred_mat +
    cbind(mc$e0[,9],mc$e1[,9]) %*% pred_mat +
    cbind(mc$g0[,9],mc$g1[,9]) %*% pred_mat
)
# convert to probability
prob_psi_auto <- array(
  NA,
  dim = dim(pred_psi_auto)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi_auto[,i,] <- sweep(
    pred_psi_auto[,i,],
    1,
    rowSums(pred_psi_auto[,i,]),
    FUN = "/"
  )
}

# calculating average psi
trueProb <- prob_psi / (
  prob_psi + (1 - prob_psi_auto)
)

WIDE <- apply(
  trueProb,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)
#### Among-city Average ####
# generate a sequence of impervious surface values that includes minimum and maximum across
# all cities surveyed
predV <- seq(-25, 58, length.out = 200)
# Generating a matrix that will contain values for predicting. This includes our sequence of
# impervious surface values, and adding a row of 1's for the intercept
pred_mat <- cbind(
  1,
  scale(predV)
)
pred_mat <- t(pred_mat)
# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept (remember to use average values for all
#                           covariates except your variable of interest)

# Generating a matrix that will contain values for predicting, including the autologistic term. 
# This includes our sequence of impervious surface values, adding a row of 1's for the 
# intercept & autologistic term
pred_mat_auto <- cbind(
  1,
  scale(predV),
  1
)
pred_mat_auto <- t(pred_mat_auto)
# create the matrix of the coefficients from the model
pred_psi <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)

# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)
# note: because this is the among-city average, these include the "mu" terms, NOT the city-
# specific terms
# fill it in. No species is easy, it's just a 1.
pred_psi[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi[,,2] <- exp(
  cbind(mc$a0_mu[,1],mc$a1_mu[,1]) %*% pred_mat
)
# Species 3 (gray squirrel)
pred_psi[,,3] <- exp(
  cbind(mc$b0_mu[,1],mc$b1_mu[,1]) %*% pred_mat
)
# Species 3 (red squirrel)
pred_psi[,,4] <- exp(
  cbind(mc$c0_mu[,1],mc$c1_mu[,1]) %*% pred_mat
)
# species 1 & 2
pred_psi[,,5] <- exp(
  cbind(mc$a0_mu[,1],mc$a1_mu[,1]) %*% pred_mat + 
    cbind(mc$b0_mu[,1],mc$b1_mu[,1]) %*% pred_mat + 
    cbind(mc$d0_mu[,1],mc$d1_mu[,1]) %*% pred_mat
)
# species 1 & 3
pred_psi[,,6] <- exp(
  cbind(mc$a0_mu[,1],mc$a1_mu[,1]) %*% pred_mat + 
    cbind(mc$c0_mu[,1],mc$c1_mu[,1]) %*% pred_mat + 
    cbind(mc$e0_mu[,1],mc$e1_mu[,1]) %*% pred_mat
)
# species 2 & 3
pred_psi[,,7] <- exp(
  cbind(mc$b0_mu[,1],mc$b1_mu[,1]) %*% pred_mat + 
    cbind(mc$c0_mu[,1],mc$c1_mu[,1]) %*% pred_mat + 
    cbind(mc$g0_mu[,1],mc$g1_mu[,1]) %*% pred_mat
)
# all species together
pred_psi[,,8] <- exp(
  cbind(mc$a0_mu[,1],mc$a1_mu[,1]) %*% pred_mat + 
    cbind(mc$b0_mu[,1],mc$b1_mu[,1]) %*% pred_mat + 
    cbind(mc$c0_mu[,1],mc$c1_mu[,1]) %*% pred_mat + 
    cbind(mc$d0_mu[,1],mc$d1_mu[,1]) %*% pred_mat +
    cbind(mc$e0_mu[,1],mc$e1_mu[,1]) %*% pred_mat +
    cbind(mc$g0_mu[,1],mc$g1_mu[,1]) %*% pred_mat
)
# convert to probability
prob_psi <- array(
  NA,
  dim = dim(pred_psi)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi[,i,] <- sweep(
    pred_psi[,i,],
    1,
    rowSums(pred_psi[,i,]),
    FUN = "/"
  )
}

# matrix should be long; ncol = length of your predictor sequence, 
#                        nrow = # covariates + intercept + autologistic term (remember to use 
#                           average values for all covariates except your variable of interest)
# create the matrix of the coefficients from the model
pred_psi_auto <- array(
  NA,
  dim = c(
    nrow(mc$a0),
    ncol(pred_mat),
    8
  )
)
# dimensions should be: number of MCMC interactions by length of sequence of predictor variable
# by number of states (8 states in a 3-species model)

# fill it in. No species is easy, it's just a 1.
pred_psi_auto[,,1] <- 1
# Species 1 (fox squirrel)
pred_psi_auto[,,2] <- exp(
  cbind(mc$a0_mu[,1],mc$a1_mu[,1],mc$phi_mu[,1]) %*% pred_mat_auto
)
# Species 2 (gray squirrel)
pred_psi_auto[,,3] <- exp(
  cbind(mc$b0_mu[,1],mc$b1_mu[,1],mc$phi_mu[,2]) %*% pred_mat_auto
)
# Species 3 (red squirrel)
pred_psi_auto[,,4] <- exp(
  cbind(mc$c0_mu[,1],mc$c1_mu[,1],mc$phi_mu[,3]) %*% pred_mat_auto
)
# species 1 & 2
pred_psi_auto[,,5] <- exp(
  cbind(mc$a0_mu[,1],mc$a1_mu[,1],mc$phi_mu[,1]) %*% pred_mat_auto + 
    cbind(mc$b0_mu[,1],mc$b1_mu[,1],mc$phi_mu[,2]) %*% pred_mat_auto + 
    cbind(mc$d0_mu[,1],mc$d1_mu[,1]) %*% pred_mat
)
# species 1 & 3
pred_psi_auto[,,6] <- exp(
  cbind(mc$a0_mu[,1],mc$a1_mu[,1],mc$phi_mu[,1]) %*% pred_mat_auto + 
    cbind(mc$c0_mu[,1],mc$c1_mu[,1],mc$phi_mu[,3]) %*% pred_mat_auto + 
    cbind(mc$e0_mu[,1],mc$e1_mu[,1]) %*% pred_mat
)
# species 2 & 3
pred_psi_auto[,,7] <- exp(
  cbind(mc$b0_mu[,1],mc$b1_mu[,1],mc$phi_mu[,2]) %*% pred_mat_auto + 
    cbind(mc$c0_mu[,1],mc$c1_mu[,1],mc$phi_mu[,3]) %*% pred_mat_auto + 
    cbind(mc$g0_mu[,1],mc$g1_mu[,1]) %*% pred_mat
)
# all species together
pred_psi_auto[,,8] <- exp(
  cbind(mc$a0_mu[,1],mc$a1_mu[,1],mc$phi_mu[,1]) %*% pred_mat_auto + 
    cbind(mc$b0_mu[,1],mc$b1_mu[,1],mc$phi_mu[,2]) %*% pred_mat_auto + 
    cbind(mc$c0_mu[,1],mc$c1_mu[,1],mc$phi_mu[,3]) %*% pred_mat_auto + 
    cbind(mc$d0_mu[,1],mc$d1_mu[,1]) %*% pred_mat +
    cbind(mc$e0_mu[,1],mc$e1_mu[,1]) %*% pred_mat +
    cbind(mc$g0_mu[,1],mc$g1_mu[,1]) %*% pred_mat
)
# convert to probability
prob_psi_auto <- array(
  NA,
  dim = dim(pred_psi_auto)
)
pb <- txtProgressBar(max = ncol(pred_mat))
for(i in 1:ncol(pred_mat)){
  setTxtProgressBar(pb, i)
  prob_psi_auto[,i,] <- sweep(
    pred_psi_auto[,i,],
    1,
    rowSums(pred_psi_auto[,i,]),
    FUN = "/"
  )
}

# calculating average psi
trueProb <- prob_psi / (
  prob_psi + (1 - prob_psi_auto)
)

ACA <- apply(
  trueProb,
  c(2,3),
  quantile,
  probs = c(0.025,0.25,0.5,0.75,0.975)
)

#### Graphing ####
# Fox squirrel marginal occupancy ####
# making the dataframe to pass to ggplot
# marginal occupancy = probability of all occupancy states including the speces of interest
# combined (e.g., fox only + fox&gray + fox&red + allSpecies)
FoxSQ <- rbind(data.frame(group = "Austin, TX\n(27%)",
                          imperv = seq(-19, 47, length.out = 200),
                          psi = (AUTX[3,,2] + AUTX[3,,5] + AUTX[3,,6] + AUTX[3,,8]),
                          upper = (AUTX[4,,2] + AUTX[4,,5] + AUTX[4,,6] + AUTX[4,,8]),
                          lower = (AUTX[2,,2] + AUTX[2,,5] + AUTX[2,,6] + AUTX[2,,8])),
               data.frame(group = "Chicago, IL\n(10%)",
                          imperv = seq(-10, 48, length.out = 200),
                          psi = (CHIL[3,,2] + CHIL[3,,5] + CHIL[3,,6] + CHIL[3,,8]),
                          upper = (CHIL[4,,2] + CHIL[4,,5] + CHIL[4,,6] + CHIL[4,,8]),
                          lower = (CHIL[2,,2] + CHIL[2,,5] + CHIL[2,,6] + CHIL[2,,8])),
               data.frame(group = "Denver, CO\n(4%)",
                          imperv = seq(-4, 16, length.out = 200),
                          psi = (DECO[3,,2] + DECO[3,,5] + DECO[3,,6] + DECO[3,,8]),
                          upper = (DECO[4,,2] + DECO[4,,5] + DECO[4,,6] + DECO[4,,8]),
                          lower = (DECO[2,,2] + DECO[2,,5] + DECO[2,,6] + DECO[2,,8])),
               data.frame(group = "Indianapolis, IN\n(11%)",
                          imperv = seq(-11, 24, length.out = 200),
                          psi = (ININ[3,,2] + ININ[3,,5] + ININ[3,,6] + ININ[3,,8]),
                          upper = (ININ[4,,2] + ININ[4,,5] + ININ[4,,6] + ININ[4,,8]),
                          lower = (ININ[2,,2] + ININ[2,,5] + ININ[2,,6] + ININ[2,,8])),
               data.frame(group = "Iowa City, IA\n(9%)",
                          imperv = seq(-9, 47, length.out = 200),
                          psi = (ICIA[3,,2] + ICIA[3,,5] + ICIA[3,,6] + ICIA[3,,8]),
                          upper = (ICIA[4,,2] + ICIA[4,,5] + ICIA[4,,6] + ICIA[4,,8]),
                          lower = (ICIA[2,,2] + ICIA[2,,5] + ICIA[2,,6] + ICIA[2,,8])),
               data.frame(group = "Los Angeles Metro,\nCA (4%)",
                          imperv = seq(-4, 20, length.out = 200),
                          psi = (MELA[3,,2] + MELA[3,,5] + MELA[3,,6] + MELA[3,,8]),
                          upper = (MELA[4,,2] + MELA[4,,5] + MELA[4,,6] + MELA[4,,8]),
                          lower = (MELA[2,,2] + MELA[2,,5] + MELA[2,,6] + MELA[2,,8])),
               data.frame(group = "Madison, WI\n(15%)",
                          imperv = 0,
                          psi = 0,
                          upper = 0,
                          lower = 0),
               data.frame(group = "Rochester, NY\n(28%)",
                          imperv = 0,
                          psi = 0,
                          upper = 0,
                          lower = 0),
               data.frame(group = "Wilmington, DE\n(25%)",
                          imperv = 0,
                          psi = 0,
                          upper = 0,
                          lower = 0))
FoxSQ$upper[FoxSQ$upper > 1] <- 1
FoxSQ$lower[FoxSQ$lower < 0] <- 0
FoxSQ$group <- factor(FoxSQ$group, levels = c("Los Angeles Metro,\nCA (4%)", 
                                              "Denver, CO\n(4%)",
                                              "Iowa City, IA\n(9%)",
                                              "Chicago, IL\n(10%)",
                                              "Indianapolis, IN\n(11%)",
                                              "Madison, WI\n(15%)",
                                              "Wilmington, DE\n(25%)",
                                              "Austin, TX\n(27%)",
                                              "Rochester, NY\n(28%)"
))
# saving the plot as its own object to make the multi-part figure
p1 <- ggplot(FoxSQ, aes(x=imperv, y=psi)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.3, show.legend=FALSE) +
  scale_fill_grafify(palette = "safe") +
  geom_line(aes(color=group), linewidth = 1.05, show.legend = FALSE) +
  scale_color_grafify(palette = "safe") +
  labs(x="", y="", color="City") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_blank(),
        axis.text = element_text(size=10))

# Gray squirrel marginal occupancy ####
GraySQ <- rbind(data.frame(group = "Austin, TX\n(27%)",
                           imperv = 0,
                           psi = 0,
                           upper = 0,
                           lower = 0),
                data.frame(group = "Chicago, IL\n(10%)",
                           imperv = seq(-10, 48, length.out = 200),
                          psi = (CHIL[3,,3] + CHIL[3,,5] + CHIL[3,,7] + CHIL[3,,8]),
                          upper = (CHIL[4,,3] + CHIL[4,,5] + CHIL[4,,7] + CHIL[4,,8]),
                          lower = (CHIL[2,,3] + CHIL[2,,5] + CHIL[2,,7] + CHIL[2,,8])),
                data.frame(group = "Denver, CO\n(4%)",
                           imperv = 0,
                           psi = 0,
                           upper = 0,
                           lower = 0),
                data.frame(group = "Indianapolis, IN\n(11%)",
                           imperv = seq(-11, 24, length.out = 200),
                           psi = (ININ[3,,3] + ININ[3,,5] + ININ[3,,7] + ININ[3,,8]),
                           upper = (ININ[4,,3] + ININ[4,,5] + ININ[4,,7] + ININ[4,,8]),
                           lower = (ININ[2,,3] + ININ[2,,5] + ININ[2,,7] + ININ[2,,8])),
               data.frame(group = "Iowa City, IA\n(9%)",
                          imperv = seq(-9, 47, length.out = 200),
                          psi = (ICIA[3,,3] + ICIA[3,,5] + ICIA[3,,7] + ICIA[3,,8]),
                          upper = (ICIA[4,,3] + ICIA[4,,5] + ICIA[4,,7] + ICIA[4,,8]),
                          lower = (ICIA[2,,3] + ICIA[2,,5] + ICIA[2,,7] + ICIA[2,,8])),
               data.frame(group = "Los Angeles Metro,\nCA (4%)",
                          imperv = 0,
                          psi = 0,
                          upper = 0,
                          lower = 0),
               data.frame(group = "Madison, WI\n(15%)",
                          imperv = seq(-15, 34, length.out = 200),
                          psi = (MAWI[3,,3] + MAWI[3,,5] + MAWI[3,,7] + MAWI[3,,8]),
                          upper = (MAWI[4,,3] + MAWI[4,,5] + MAWI[4,,7] + MAWI[4,,8]),
                          lower = (MAWI[2,,3] + MAWI[2,,5] + MAWI[2,,7] + MAWI[2,,8])),
               data.frame(group = "Rochester, NY\n(28%)",
                          imperv = seq(-25, 49, length.out = 200),
                          psi = (RONY[3,,3] + RONY[3,,5] + RONY[3,,7] + RONY[3,,8]),
                          upper = (RONY[4,,3] + RONY[4,,5] + RONY[4,,7] + RONY[4,,8]),
                          lower = (RONY[2,,3] + RONY[2,,5] + RONY[2,,7] + RONY[2,,8])),
               data.frame(group = "Wilmington, DE\n(25%)",
                          imperv = seq(-25, 58, length.out = 200),
                          psi = (WIDE[3,,3] + WIDE[3,,5] + WIDE[3,,7] + WIDE[3,,8]),
                          upper = (WIDE[4,,3] + WIDE[4,,5] + WIDE[4,,7] + WIDE[4,,8]),
                          lower = (WIDE[2,,3] + WIDE[2,,5] + WIDE[2,,7] + WIDE[2,,8])))
GraySQ$upper[GraySQ$upper > 1] <- 1
GraySQ$group <- factor(GraySQ$group, levels = c("Los Angeles Metro,\nCA (4%)", 
                                              "Denver, CO\n(4%)",
                                              "Iowa City, IA\n(9%)",
                                              "Chicago, IL\n(10%)",
                                              "Indianapolis, IN\n(11%)",
                                              "Madison, WI\n(15%)",
                                              "Wilmington, DE\n(25%)",
                                              "Austin, TX\n(27%)",
                                              "Rochester, NY\n(28%)"
))
p2 <- ggplot(GraySQ, aes(x=imperv, y=psi)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.3, show.legend=FALSE) +
  scale_fill_grafify(palette = "safe") +
  geom_line(aes(color=group), linewidth = 1.05, show.legend = FALSE) +
  scale_color_grafify(palette = "safe") +
  labs(x="", y="", color="City") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=10), axis.text.y = element_blank(),
        axis.title=element_blank())

# Red squirrel marginal occupancy ####
RedSQ <- rbind(data.frame(group = "Austin, TX\n(27%)",
                          imperv = 0,
                          psi = 0,
                          upper = 0,
                          lower = 0),
               data.frame(group = "Chicago, IL\n(10%)",
                          imperv = 0,
                          psi = 0,
                          upper = 0,
                          lower = 0),
               data.frame(group = "Denver, CO\n(4%)",
                          imperv = 0,
                          psi = 0,
                          upper = 0,
                          lower = 0),
               data.frame(group = "Indianapolis, IN\n(11%)",
                          imperv = seq(-11, 24, length.out = 200),
                          psi = (ININ[3,,4] + ININ[3,,6] + ININ[3,,7] + ININ[3,,8]),
                          upper = (ININ[4,,4] + ININ[4,,6] + ININ[4,,7] + ININ[4,,8]),
                          lower = (ININ[2,,4] + ININ[2,,6] + ININ[2,,7] + ININ[2,,8])),
               data.frame(group = "Iowa City, IA\n(9%)",
                          imperv = 0,
                          psi = 0,
                          upper = 0,
                          lower = 0),
               data.frame(group = "Los Angeles Metro,\nCA (4%)",
                          imperv = 0,
                          psi = 0,
                          upper = 0,
                          lower = 0),
               data.frame(group = "Madison, WI\n(15%)",
                          imperv = 0,
                          psi = 0,
                          upper = 0,
                          lower = 0),
               data.frame(group = "Rochester, NY\n(28%)",
                          imperv = seq(-25, 49, length.out = 200),
                           psi = (RONY[3,,4] + RONY[3,,6] + RONY[3,,7] + RONY[3,,8]),
                           upper = (RONY[4,,4] + RONY[4,,6] + RONY[4,,7] + RONY[4,,8]),
                           lower = (RONY[2,,4] + RONY[2,,6] + RONY[2,,7] + RONY[2,,8])),
               data.frame(group = "Wilmington, DE\n(25%)",
                          imperv = seq(-25, 58, length.out = 200),
                          psi = (WIDE[3,,4] + WIDE[3,,6] + WIDE[3,,7] + WIDE[3,,8]),
                          upper = (WIDE[4,,4] + WIDE[4,,6] + WIDE[4,,7] + WIDE[4,,8]),
                          lower = (WIDE[2,,4] + WIDE[2,,6] + WIDE[2,,7] + WIDE[2,,8]))
)
RedSQ$upper[RedSQ$upper > 1] <- 1
RedSQ$group <- factor(RedSQ$group, levels = c("Los Angeles Metro,\nCA (4%)", 
                                              "Denver, CO\n(4%)",
                                              "Iowa City, IA\n(9%)",
                                              "Chicago, IL\n(10%)",
                                              "Indianapolis, IN\n(11%)",
                                              "Madison, WI\n(15%)",
                                              "Wilmington, DE\n(25%)",
                                              "Austin, TX\n(27%)",
                                              "Rochester, NY\n(28%)"
                                              ))

p3 <- ggplot(RedSQ, aes(x=imperv, y=psi)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=.3, show.legend=FALSE) +
  scale_fill_grafify(palette = "safe") +
  geom_line(aes(color=group), linewidth = 1.05, show.legend = FALSE) +
  scale_color_grafify(palette = "safe") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x="", y="", color="City") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=10), axis.text.y = element_blank(),
        axis.title=element_blank())
# this 'p4' object pulls out the legend as its own object (for more info, see the 'get_legend'
# function of the 'ggpubr' package)
p4 <- ggplot(RedSQ, aes(x=imperv, y=psi))+
  geom_line(aes(color=group), linewidth = 1.05) +
  scale_color_grafify(palette = "safe") +
  theme_bw() +
  theme(legend.position = c(0.5,0.5), legend.title = element_blank(), legend.background=element_blank(),
        legend.text = element_text(size=6))

leg <- get_legend(p4)
leg <- as_ggplot(leg)
# the code below makes & saves Figure 2
jpeg("./results/MarginalOccupancy_all.jpeg", width = 6, height = 4, units = 'in', res = 300)
gridExtra::grid.arrange(p1,p2,p3,leg, ncol=4, widths=c(1.5,1.25,1.25,1),
                        bottom=textGrob("Mean-centered Canopy Cover (%)", gp=gpar(fontsize=12)),
                        left="Marginal Occupancy")
dev.off()

#### Fox and Gray conditional occupancy ####
# the formatting of this data.frame is a little different; here we are subsetting
# posterior distributions to their 95% credible interval
FG <- rbind(data.frame(city="Chicago,\nIL",
                       dist = subset(mc$d0[,2], mc$d0[,2] > -2.17 & mc$d0[,2] < 0.74),
                       var = "Intercept"),
            data.frame(city="Chicago,\nIL",
                       dist = subset(mc$d1[,2], mc$d1[,2] > -1.14 & mc$d1[,2] < 0.39),
                       var = "Canopy"),
            data.frame(city="Indianapolis,\nIN",
                       dist = subset(mc$d0[,4], mc$d0[,4] > -5.12 & mc$d0[,4] < -2.50),
                       var = "Intercept"),
            data.frame(city="Indianapolis,\nIN",
                       dist = subset(mc$d1[,4], mc$d1[,4] > -1.16 & mc$d1[,4] < 0.95),
                       var = "Canopy"),
            data.frame(city="Iowa City,\nIA",
                       dist = subset(mc$d0[,5], mc$d0[,5] > -0.52 & mc$d0[,5] < 3.22),
                       var = "Intercept"),
            data.frame(city="Iowa City,\nIA",
                       dist = subset(mc$d1[,5], mc$d1[,5] > -1.20 & mc$d1[,5] < 0.67),
                       var = "Canopy"),
            data.frame(city="Among City\nAverage",
                       dist = subset(mc$d0_mu[,1], mc$d0_mu[,1] > -1.75 & mc$d0_mu[,1] < 1.40),
                       var = "Intercept"),
            data.frame(city="Among City\nAverage",
                       dist = subset(mc$d1_mu[,1], mc$d1_mu[,1] > -0.68 & mc$d1_mu[,1] < 0.87),
                       var = "Canopy")
)

library(ggridges)
p1 <- ggplot(FG, aes(x=dist, y=city, fill=city, color=city)) +
  geom_vline(aes(xintercept = 0), linetype=2, color="gray") +
  stat_density_ridges(rel_min_height = 0.01, scale = 1.1, quantile_lines = TRUE, quantiles = 2, alpha=0.5) +
  scale_fill_manual(values=c("#000000","#117733","#332288","#DDCC77")) +
  scale_color_manual(values = c("#000000","#117733","#332288","#DDCC77")) +
  xlim(-6, 9) +
  scale_y_discrete(expand = expansion(mult=c(0.01, 0.3))) +
  facet_wrap(~factor(var, levels=c('Intercept','Canopy')), ncol=1, strip.position = "top") +
  theme_ridges() +
  theme(axis.title=element_blank(), axis.text.x=element_text(size=8), axis.text.y=element_text(size=6),
        strip.background = element_blank(), strip.text=element_text(size=12,hjust=0.01,margin=margin(1,1,2,1,"pt")),
        legend.position = "none")

#### Gray and Red conditional occupancy ####
GR <- rbind(data.frame(city="Indianapolis,\nIN",
                       dist = subset(mc$g0[,4], mc$g0[,4] > 1.29 & mc$g0[,4] < 5.14),
                       var = "Intercept"),
            data.frame(city="Indianapolis,\nIN",
                       dist = subset(mc$g1[,4], mc$g1[,4] > -0.72 & mc$g1[,4] < 1.76),
                       var = "Canopy"),
            data.frame(city="Rochester,\nNY",
                       dist = subset(mc$g0[,8], mc$g0[,8] > -5.94 & mc$g0[,8] < 8.55),
                       var = "Intercept"),
            data.frame(city="Rochester,\nNY",
                       dist = subset(mc$g1[,8], mc$g1[,8] > -0.85 & mc$g1[,8] < 1.97),
                       var = "Canopy"),
            data.frame(city="Wilmington,\nDE",
                       dist = subset(mc$g0[,9], mc$g0[,9] > -5.33 & mc$g0[,9] < 4.29),
                       var = "Intercept"),
            data.frame(city="Wilmington,\nDE",
                       dist = subset(mc$g1[,9], mc$g1[,9] > -1.75 & mc$g1[,9] < 1.54),
                       var = "Canopy"),
            data.frame(city="Among City\nAverage",
                       dist = subset(mc$g0_mu[,1], mc$g0_mu[,1] > -1.95 & mc$g0_mu[,1] < 3.17),
                       var = "Intercept"),
            data.frame(city="Among City\nAverage",
                       dist = subset(mc$g1_mu[,1], mc$g1_mu[,1] > -0.61 & mc$g1_mu[,1] < 1.23),
                       var = "Canopy")
)

p2 <- ggplot(GR, aes(x=dist, y=city, fill=city, color=city)) +
  geom_vline(aes(xintercept = 0), linetype=2, color="gray") +
  stat_density_ridges(rel_min_height = 0.01, scale = 1.1, quantile_lines = TRUE, quantiles = 2, alpha=0.5) +
  scale_fill_manual(values=c("#000000","#332288","#882255","#44AA99")) +
  scale_color_manual(values = c("#000000","#332288","#882255","#44AA99")) +
  xlim(-6, 9) +
  scale_y_discrete(expand = expansion(mult=c(0.01, 0.3))) +
  facet_wrap(~factor(var, levels=c('Intercept','Canopy')), ncol=1, strip.position = "top") +
  theme_ridges() +
  theme(axis.title=element_blank(), axis.text.x=element_text(size=8), axis.text.y=element_text(size=6),
        strip.background = element_blank(), strip.text=element_text(size=12,hjust=0.01,margin=margin(1,1,2,1,"pt")),
        legend.position = "none")

#### Fox and Red conditional occupancy ####
FR <- rbind(data.frame(city="Indianapolis,\nIN",
                       dist = subset(mc$e0[,4], mc$e0[,4] > 3.22 & mc$e0[,4] < 8.55),
                       var = "Intercept"),
            data.frame(city="Indianapolis,\nIN",
                       dist = subset(mc$e1[,4], mc$e1[,4] > -1.89 & mc$e1[,4] < 0.91),
                       var = "Canopy"),
            data.frame(city="Among City\nAverage",
                       dist = subset(mc$e0_mu[,1], mc$e0_mu[,1] > -4.15 & mc$e0_mu[,1] < 1.58),
                       var = "Intercept"),
            data.frame(city="Among City\nAverage",
                       dist = subset(mc$e1_mu[,1], mc$e1_mu[,1] > -1.37 & mc$e1_mu[,1] < 0.98),
                       var = "Canopy")
            )

p3 <- ggplot(FR, aes(x=dist, y=city, fill=city, color=city)) +
  geom_vline(aes(xintercept = 0), linetype=2, color="gray") +
  stat_density_ridges(rel_min_height = 0.01, scale = 1.1, quantile_lines = TRUE, quantiles = 2, alpha=0.5) +
  scale_fill_manual(values=c("#000000","#332288")) +
  scale_color_manual(values = c("#000000","#332288")) +
  xlim(-6, 9) +
  scale_y_discrete(expand = expansion(mult=c(0.01, 0.3))) +
  facet_wrap(~factor(var, levels=c('Intercept','Canopy')), ncol=1, strip.position = "top") +
  theme_ridges() +
  theme(axis.title=element_blank(), axis.text.x=element_text(size=8), axis.text.y=element_text(size=6),
        strip.background = element_blank(), strip.text=element_text(size=12,hjust=0.01,margin=margin(1,1,2,1,"pt")),
        legend.position = "none")
# the code below plots & saves Figure 3
jpeg("./results/CondDists.jpeg", width = 6, height = 4, units = 'in', res = 300)
gridExtra::grid.arrange(p1,p2,p3, ncol=3,
                        bottom=textGrob("Effect Size (logit scale)", gp=gpar(fontsize=12))
                        )
dev.off()

#### Figure 4 - Iowa City Example ####
# this code plots expected vs. model-estimated co-occurrence probabilites
FG_pred <- rbind(data.frame(canopy = seq(-9, 47, length.out = 200),
                            psi = ((ICIA[3,,2] + ICIA[3,,5] + ICIA[3,,6] + ICIA[3,,8])*(ICIA[3,,3] + ICIA[3,,5] + ICIA[3,,7] + ICIA[3,,8])),
                            upper = ((ICIA[5,,2] + ICIA[5,,5] + ICIA[5,,6] + ICIA[5,,8])*(ICIA[5,,3] + ICIA[5,,5] + ICIA[5,,7] + ICIA[5,,8])),
                            lower = ((ICIA[1,,2] + ICIA[1,,5] + ICIA[1,,6] + ICIA[1,,8])*(ICIA[1,,3] + ICIA[1,,5] + ICIA[1,,7] + ICIA[1,,8])),
                            state = "Expected Co-occurrence"),
              data.frame(canopy = seq(-9, 47, length.out = 200),
                          psi = (ICIA[3,,5] + ICIA[3,,8]) / (ICIA[3,,5] + ICIA[3,,6] + ICIA[3,,2] + ICIA[3,,8]),
                          lower = (ICIA[5,,5] + ICIA[5,,8]) / (ICIA[5,,5] + ICIA[5,,6] + ICIA[5,,2] + ICIA[5,,8]),
                          upper = (ICIA[1,,5] + ICIA[1,,8]) / (ICIA[1,,5] + ICIA[1,,6] + ICIA[1,,2] + ICIA[1,,8]),
                          state = "Actual Co-occurrence")
)

p1<-ggplot(FG_pred, aes(x=canopy, y=psi)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=state), alpha=.3, show.legend=FALSE) +
  scale_fill_manual(values=c("gray30","gray60"), guide = "none") +
  geom_line(aes(color=state, linetype=state), linewidth=1.05, show.legend = FALSE) +
  scale_color_manual(values=c("gray30","gray60"), guide = "none") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(title="", x="", y="", linetype="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=10))

# Indianapolis Example ####
GR_pred <- rbind(data.frame(canopy = seq(-11, 24, length.out = 200),
                          psi = ((ININ[3,,4] + ININ[3,,6] + ININ[3,,7] + ININ[3,,8])*(ININ[3,,3] + ININ[3,,5] + ININ[3,,7] + ININ[3,,8])),
                          upper = ((ININ[5,,4] + ININ[5,,6] + ININ[5,,7] + ININ[5,,8])*(ININ[5,,3] + ININ[5,,5] + ININ[5,,7] + ININ[5,,8])),
                          lower = ((ININ[1,,4] + ININ[1,,6] + ININ[1,,7] + ININ[1,,8])*(ININ[1,,3] + ININ[1,,5] + ININ[1,,7] + ININ[1,,8])),
                          state = "Expected\nCo-occurrence"),
               data.frame(canopy = seq(-11, 24, length.out = 200),
                          psi = (ININ[3,,7] + ININ[3,,8]) / (ININ[3,,3] + ININ[3,,5] + ININ[3,,7] + ININ[3,,8]),
                          upper = (ININ[5,,7] + ININ[5,,8]) / (ININ[5,,3] + ININ[5,,5] + ININ[5,,7] + ININ[5,,8]),
                          lower = (ININ[1,,7] + ININ[1,,8]) / (ININ[1,,3] + ININ[1,,5] + ININ[1,,7] + ININ[1,,8]),
                          state = "Actual\nCo-occurrence")
               )
GR_pred$upper[GR_pred$upper > 1] <- 1

p2<-ggplot(GR_pred, aes(x=canopy, y=psi)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=state), alpha=.3, show.legend=FALSE) +
  scale_fill_manual(values=c("gray30","gray60"), guide = "none") +
  geom_line(aes(color=state, linetype=state), linewidth=1.05, show.legend=FALSE) +
  scale_color_manual(values=c("gray30","gray60"), guide = "none") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(title="", x="", y="", linetype="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=10), axis.text.y=element_blank(),
        legend.spacing.y=unit(1, 'cm')) +
  guides(linetype = guide_legend(byrow = TRUE))

p3 <- ggplot(GR_pred, aes(x=canopy, y=psi)) +
  geom_line(aes(linetype=state), linewidth = 1.05) +
  theme_bw() +
  theme(legend.position = c(0.5,0.5), legend.title = element_blank(), legend.background=element_blank(),
        legend.text = element_text(size=8), legend.spacing.y=unit(1, 'cm')) +
  guides(linetype = guide_legend(byrow = TRUE))

leg <- get_legend(p3)
leg <- as_ggplot(leg)
jpeg("./results/CoOccupancy.jpeg", width = 6, height = 4, units = 'in', res = 300)
gridExtra::grid.arrange(p1,p2,leg, ncol=3, widths=c(1.75,1.6,1),
                        bottom=textGrob("Mean-centered Canopy Cover (%)", gp=gpar(fontsize=12)),
                        left="Co-occurrence Probability")
dev.off()
