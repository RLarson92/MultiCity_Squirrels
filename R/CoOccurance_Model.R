# load housekeeping & utility functions
functions_to_load <- list.files(
  "./R/functions/",
  full.names = TRUE
)
for(fn in functions_to_load){
  source(fn)
}
library(tidyr)

# read in the detection/non-detection data for each species. I tend to use the same abbreviation for the species
# names throughout: 'FS' = fox squirrel; 'GS' = gray squirrel; 'RS' = red squirrel
FS_occu <- read.csv("./data/FS.csv", stringsAsFactors = FALSE)
fox_dets <- FS_occu[,3:38] # pulling out the columns that have all the 1/0/NA's
GS_occu <- read.csv("./data/GS.csv", stringsAsFactors = FALSE)
gray_dets <- GS_occu[,3:38]
RS_occu <- read.csv("./data/RS.csv", stringsAsFactors = FALSE)
red_dets <- RS_occu[,3:38]

# take the three seperate squirrel data sets and stack the seasons on top of each other
FS_stack <- wide_to_stacked(fox_dets, 4, 1)
GS_stack <- wide_to_stacked(gray_dets, 4, 2)
RS_stack <- wide_to_stacked(red_dets, 4, 3)
# bind them together in order by species
y <- rbind(FS_stack, GS_stack, RS_stack)
# stack the replicates (weeks) on top of each other and order the datasheet. Data should be in order by species, then season, then city, then site
y_long <- gather(y, key = "Week", value = "Y", week1:week4, factor_key = TRUE)
y_long <- y_long[order(
  y_long[,"Species"], 
  y_long[,"Season"],
  y_long[,"City"],
  y_long[,"Site"]),]
# converting some text-based columns to numeric values
y_long$Week <- as.numeric(y_long$Week)
y_long$Site <- as.numeric(as.factor(y_long$Site))
y_long$City <- as.numeric(as.factor(y_long$City))

# create a blank array that has the correct dimensions, in this case 4 
# (species by site by week by season)
my_array <- array(
  NA,
  dim = c(
    max(y_long$Species),
    max(y_long$Site),
    max(y_long$Week),
    max(y_long$Season)
  )
)
# fill in the array with the information from the long-form data.frame
for(i in 1:nrow(y_long)){
  my_array[
    y_long$Species[i],
    y_long$Site[i],
    y_long$Week[i],
    y_long$Season[i]
  ] <- y_long$Y[i]
}

# reading in detection-level covariates (# active camera-days)
obs_covs <- read.csv("./data/cooc_obsCovs.csv", stringsAsFactors = FALSE)
# this data.frame needs to be in stacked format and ordered, the same as
# the detection histories, so let's do that
FS_obs <- wide_to_stacked(obs_covs[,3:38], 4, 1) 
GS_obs <- wide_to_stacked(obs_covs[,3:38], 4, 2)
RS_obs <- wide_to_stacked(obs_covs[,3:38], 4, 3)
obs <- rbind(FS_obs, GS_obs, RS_obs)
obs_long <- gather(obs, key = "Week", value = "J", week1:week4, factor_key = TRUE)
obs_long <- obs_long[order(
  obs_long[,"Species"],
  obs_long[,"Season"],
  obs_long[,"City"],
  obs_long[,"Site"]),]
obs_long$Week <- as.numeric(obs_long$Week)
obs_long$Site <- as.numeric(as.factor(obs_long$Site))

# make the array and fill it in
obsArray <- array(data = NA, dim = c(max(obs_long$Species),
                                     max(obs_long$Site),
                                     max(obs_long$Week),
                                     max(obs_long$Season)
                                     )
                  )
for(i in 1:nrow(obs_long)){
  obsArray[
    obs_long$Species[i],
    obs_long$Site[i],
    obs_long$Week[i],
    obs_long$Season[i]
  ] <- obs_long$J[i]
}

# we will handle occupancy covariates below, under the 'JAGS' heading, as they
# vary from model to model

# last but not least, making the array of possible occupancy states
# (e.g., 1/1/1 = all species present; 0/1/0 = only gray squirrel present, etc.)
Xcat <- matrix(c(1, 1, 1,
                 1, 0, 0,
                 0, 1, 0,
                 0, 0, 1,
                 1, 1, 0,
                 1, 0, 1,
                 0, 1, 1,
                 0, 0, 0), ncol = 3, byrow = TRUE)

library(runjags)
##### JAGS #####
#### Fit the Model (PCA model) ####
# read in site-level covariates
site_covariates <- read.csv("./data/cooc_covs.csv", stringsAsFactors = FALSE)
continuous_covariates <- site_covariates[,c(8,11,14)] # pulling out just the covariates, sans the city/site names
library(psych)
PCA <- prcomp(continuous_covariates, center = TRUE, scale. = TRUE)
summary(PCA)
library(ggfortify)
biplot(PCA)
# PC1: more +values = more impervious, more -values = more canopy+grass
# PC2: more +values = more canopy, more -values = more grass
covs_for_model <- PCA$x[,1:2]
covMatrix <- as.matrix(covs_for_model)

# Data list for model
data_list <- list(
  nspec = max(y_long$Species),
  ncity = max(y_long$City),
  nsite = max(y_long$Site),
  PC1 = covMatrix[,1],
  PC2 = covMatrix[,2],
  nrep = max(y_long$Week),
  trapDays = obsArray, #needs to be a 4-dimensional array
  y = my_array, #needs to be a 4-dimensional array
  nseason = max(y_long$Season),
  Xcat = Xcat,
  city_vec = as.numeric(factor(FS_occu$City))
)

# initial starting values for each chain
PCA_inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      x = matrix(1,
                 ncol = data_list$nseason,
                 nrow = data_list$nsite),
      p0_mu = rnorm(data_list$nspec),
      p0_tau = rgamma(data_list$nspec,1,1),
      p0 = matrix(1,
                  ncol = data_list$ncity,
                  nrow = data_list$nspec), 
      p1_mu = rnorm(data_list$nspec),
      p1_tau = rgamma(data_list$nspec,1,1),
      p1 = matrix(1,
                  ncol = data_list$ncity,
                  nrow = data_list$nspec), 
      phi_mu = rnorm(data_list$nspec),
      phi_tau = rgamma(data_list$nspec,1,1),
      phi = matrix(1,
                   ncol = data_list$ncity,
                   nrow = data_list$nspec), 
      a0_mu = rnorm(1),
      a0_tau = rgamma(1,1,1),
      a0 = rnorm(data_list$ncity),
      a1_mu = rnorm(1),
      a1_tau = rgamma(1,1,1),
      a1 = rnorm(data_list$ncity),
      a2_mu = rnorm(1),
      a2_tau = rgamma(1,1,1),
      a2 = rnorm(data_list$ncity),
      b0_mu = rnorm(1),
      b0_tau = rgamma(1,1,1),
      b0 = rnorm(data_list$ncity),
      b1_mu = rnorm(1),
      b1_tau = rgamma(1,1,1),
      b1 = rnorm(data_list$ncity),
      b2_mu = rnorm(1),
      b2_tau = rgamma(1,1,1),
      b2 = rnorm(data_list$ncity),
      c0_mu = rnorm(1),
      c0_tau = rgamma(1,1,1),
      c0 = rnorm(data_list$ncity),
      c1_mu = rnorm(1),
      c1_tau = rgamma(1,1,1),
      c1 = rnorm(data_list$ncity),
      c2_mu = rnorm(1),
      c2_tau = rgamma(1,1,1),
      c2 = rnorm(data_list$ncity),
      d0_mu = rnorm(1),
      d0_tau = rgamma(1,1,1),
      d0 = rnorm(data_list$ncity), 
      d1_mu = rnorm(1),
      d1_tau = rgamma(1,1,1),
      d1 = rnorm(data_list$ncity),
      d2_mu = rnorm(1),
      d2_tau = rgamma(1,1,1),
      d2 = rnorm(data_list$ncity),
      e0_mu = rnorm(1),
      e0_tau = rgamma(1,1,1),
      e0 = rnorm(data_list$ncity), 
      e1_mu = rnorm(1),
      e1_tau = rgamma(1,1,1),
      e1 = rnorm(data_list$ncity), 
      e2_mu = rnorm(1),
      e2_tau = rgamma(1,1,1),
      e2 = rnorm(data_list$ncity), 
      g0_mu = rnorm(1),
      g0_tau = rgamma(1,1,1),
      g0 = rnorm(data_list$ncity), 
      g1_mu = rnorm(1),
      g1_tau = rgamma(1,1,1),
      g1 = rnorm(data_list$ncity),
      g2_mu = rnorm(1),
      g2_tau = rgamma(1,1,1),
      g2 = rnorm(data_list$ncity)
    )
  }
  return(
    switch(
      chain,
      "1" = gen_list(chain),
      "2" = gen_list(chain),
      "3" = gen_list(chain),
      "4" = gen_list(chain),
      "5" = gen_list(chain),
      "6" = gen_list(chain),
      "7" = gen_list(chain),
      "8" = gen_list(chain)
    )
  )
}

# Fit the model
my_start <- Sys.time()
PCA_mod <- runjags::run.jags(
  model = "./JAGS/multiSpec_PCA_model.R",
  monitor = c("a0_mu", "a0_sd", "a1_mu", "a1_sd", "a2_mu", "a2_sd", "a0", "a1", "a2",
              "b0_mu", "b0_sd", "b1_mu", "b1_sd", "b2_mu", "b2_sd", "b0", "b1", "b2",
              "c0_mu", "c0_sd", "c1_mu", "c1_sd", "c2_mu", "c2_sd", "c0", "c1", "c2",
              "d0_mu", "d0_sd", "d1_mu", "d1_sd", "d2_mu", "d2_sd", "d0", "d1", "d2",
              "e0_mu", "e0_sd", "e1_mu", "e1_sd", "e2_mu", "e2_sd", "e0", "e1", "e2",
              "g0_mu", "g0_sd", "g1_mu", "g1_sd", "g2_mu", "g2_sd", "g0", "g1", "g2",
              "p0_mu", "p0_sd", "p1_mu", "p1_sd", "p0", "p1",
              "phi_mu", "phi_sd", "phi",
              "x"),
  data = data_list,
  n.chains = 3,
  inits = PCA_inits,
  burnin = 50000,
  sample = 20000,
  adapt = 500,
  modules = "glm",
  thin = 2,
  method = "parallel"
)
my_end <- Sys.time()

saveRDS(PCA_mod, "./results/PCA_model.RDS")

#### Imperv Only Model ####
# pulling out the continuous site covariates & scaling them
site_covariates <- read.csv("./data/cooc_covs.csv", stringsAsFactors = FALSE)
continuous_covariates <- site_covariates[,12]
continuous_covariates <- scale(continuous_covariates)
covs_for_model <- cbind.data.frame("imperv"=continuous_covariates[,1],
                                   "city"=site_covariates[,2],
                                   "site"=site_covariates[,1]
)
covs_for_model$city <- as.numeric(as.factor(covs_for_model$city))
covs_for_model$site <- as.numeric(as.factor(covs_for_model$site))
covMatrix <- as.matrix(covs_for_model)

# Data list for model
data_list <- list(
  nspec = max(y_long$Species),
  ncity = max(y_long$City),
  nsite = max(y_long$Site),
  imperv = covMatrix[,1],
  nrep = max(y_long$Week),
  trapDays = obsArray, #needs to be a 4-dimensional array
  y = my_array, #needs to be a 4-dimensional array
  nseason = max(y_long$Season),
  Xcat = Xcat,
  city_vec = as.numeric(factor(FS_occu$City))
)

# Fit the model
my_start <- Sys.time()
imperv_mod <- runjags::run.jags(
  model = "./JAGS/multiSpec_model_ranef.R",
  monitor = c("a0_mu", "a0_sd", "a1_mu", "a1_sd", "a0", "a1",
              "b0_mu", "b0_sd", "b1_mu", "b1_sd", "b0", "b1",
              "c0_mu", "c0_sd", "c1_mu", "c1_sd", "c0", "c1",
              "d0_mu", "d0_sd", "d1_mu", "d1_sd", "d0", "d1",
              "e0_mu", "e0_sd", "e1_mu", "e1_sd", "e0", "e1",
              "g0_mu", "g0_sd", "g1_mu", "g1_sd", "g0", "g1",
              "p0_mu", "p0_sd", "p1_mu", "p1_sd", "p0", "p1",
              "phi_mu", "phi_sd", "phi", "x"),
  data = data_list,
  n.chains = 3,
  inits = my_inits,
  burnin = 50000,
  sample = 20000,
  adapt = 500,
  modules = "glm",
  thin = 2,
  method = "parallel"
)
my_end <- Sys.time()

saveRDS(imperv_mod, "./results/imperv_model.RDS")
#### Canopy Only Model ####
# pulling out the continuous site covariates & scaling them
site_covariates <- read.csv("./data/cooc_covs.csv", stringsAsFactors = FALSE)
continuous_covariates <- site_covariates[,9]
continuous_covariates <- scale(continuous_covariates)
cov(continuous_covariates)
covs_for_model <- cbind.data.frame("canopy"=continuous_covariates[,1],
                                   "city"=site_covariates[,2],
                                   "site"=site_covariates[,1]
)
covs_for_model$city <- as.numeric(as.factor(covs_for_model$city))
covs_for_model$site <- as.numeric(as.factor(covs_for_model$site))
covMatrix <- as.matrix(covs_for_model)

# Data list for model
data_list <- list(
  nspec = max(y_long$Species),
  ncity = max(y_long$City),
  nsite = max(y_long$Site),
  canopy = covMatrix[,1],
  nrep = max(y_long$Week),
  trapDays = obsArray, #needs to be a 4-dimensional array
  y = my_array, #needs to be a 4-dimensional array
  nseason = max(y_long$Season),
  Xcat = Xcat,
  city_vec = as.numeric(factor(FS_occu$City))
)

# Fit the model
my_start <- Sys.time()
canopy_mod <- runjags::run.jags(
  model = "./JAGS/multiSpec_Canopy.R",
  monitor = c("a0_mu", "a0_sd", "a1_mu", "a1_sd", "a0", "a1",
              "b0_mu", "b0_sd", "b1_mu", "b1_sd", "b0", "b1",
              "c0_mu", "c0_sd", "c1_mu", "c1_sd", "c0", "c1",
              "d0_mu", "d0_sd", "d1_mu", "d1_sd", "d0", "d1",
              "e0_mu", "e0_sd", "e1_mu", "e1_sd", "e0", "e1",
              "g0_mu", "g0_sd", "g1_mu", "g1_sd", "g0", "g1",
              "p0_mu", "p0_sd", "p1_mu", "p1_sd", "p0", "p1",
              "phi_mu", "phi_sd", "phi","x"),
  data = data_list,
  n.chains = 3,
  inits = my_inits,
  burnin = 50000,
  sample = 20000,
  adapt = 500,
  modules = "glm",
  thin = 2,
  method = "parallel"
)
my_end <- Sys.time()
saveRDS(canopy_mod, "./results/canopy_model.RDS")
#### Grass Only Model #####
# pulling out the continuous site covariates & scaling them
site_covariates <- read.csv("./data/cooc_covs.csv", stringsAsFactors = FALSE)
continuous_covariates <- site_covariates[,15]
continuous_covariates <- scale(continuous_covariates)
cov(continuous_covariates)
covs_for_model <- cbind.data.frame("grass"=continuous_covariates[,1],
                                   "city"=site_covariates[,2],
                                   "site"=site_covariates[,1]
)
covs_for_model$city <- as.numeric(as.factor(covs_for_model$city))
covs_for_model$site <- as.numeric(as.factor(covs_for_model$site))
covMatrix <- as.matrix(covs_for_model)

# Data list for model
data_list <- list(
  nspec = max(y_long$Species),
  ncity = max(y_long$City),
  nsite = max(y_long$Site),
  grass = covMatrix[,1],
  nrep = max(y_long$Week),
  trapDays = obsArray, #needs to be a 4-dimensional array
  y = my_array, #needs to be a 4-dimensional array
  nseason = max(y_long$Season),
  Xcat = Xcat,
  city_vec = as.numeric(factor(FS_occu$City))
)

# Fit the model
my_start <- Sys.time()
grass_mod <- runjags::run.jags(
  model = "./JAGS/multiSpec_grass.R",
  monitor = c("a0_mu", "a0_sd", "a1_mu", "a1_sd", "a0", "a1",
              "b0_mu", "b0_sd", "b1_mu", "b1_sd", "b0", "b1",
              "c0_mu", "c0_sd", "c1_mu", "c1_sd", "c0", "c1",
              "d0_mu", "d0_sd", "d1_mu", "d1_sd", "d0", "d1",
              "e0_mu", "e0_sd", "e1_mu", "e1_sd", "e0", "e1",
              "g0_mu", "g0_sd", "g1_mu", "g1_sd", "g0", "g1",
              "p0_mu", "p0_sd", "p1_mu", "p1_sd", "p0", "p1",
              "phi_mu", "phi_sd", "phi","x"),
  data = data_list,
  n.chains = 3,
  inits = my_inits,
  burnin = 50000,
  sample = 20000,
  adapt = 500,
  modules = "glm",
  thin = 2,
  method = "parallel"
)
my_end <- Sys.time()
saveRDS(grass_mod, "./results/grass_model.RDS")
##### Null Model ####
# Data list for model
data_list <- list(
  nspec = max(y_long$Species),
  ncity = max(y_long$City),
  nsite = max(y_long$Site),
  nrep = max(y_long$Week),
  trapDays = obsArray, #needs to be a 4-dimensional array
  y = my_array, #needs to be a 4-dimensional array
  nseason = max(y_long$Season),
  Xcat = Xcat,
  city_vec = as.numeric(factor(FS_occu$City))
)
# initial starting values for each chain
null_inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      x = matrix(1,
                 ncol = data_list$nseason,
                 nrow = data_list$nsite),
      p0_mu = rnorm(data_list$nspec),
      p0_tau = rgamma(data_list$nspec,1,1),
      p0 = matrix(1,
                  ncol = data_list$ncity,
                  nrow = data_list$nspec), 
      p1_mu = rnorm(data_list$nspec),
      p1_tau = rgamma(data_list$nspec,1,1),
      p1 = matrix(1,
                  ncol = data_list$ncity,
                  nrow = data_list$nspec), 
      phi_mu = rnorm(data_list$nspec),
      phi_tau = rgamma(data_list$nspec,1,1),
      phi = matrix(1,
                   ncol = data_list$ncity,
                   nrow = data_list$nspec), 
      a0_mu = rnorm(1),
      a0_tau = rgamma(1,1,1),
      a0 = rnorm(data_list$ncity),
      b0_mu = rnorm(1),
      b0_tau = rgamma(1,1,1),
      b0 = rnorm(data_list$ncity),
      c0_mu = rnorm(1),
      c0_tau = rgamma(1,1,1),
      c0 = rnorm(data_list$ncity),
      d0_mu = rnorm(1),
      d0_tau = rgamma(1,1,1),
      d0 = rnorm(data_list$ncity), 
      e0_mu = rnorm(1),
      e0_tau = rgamma(1,1,1),
      e0 = rnorm(data_list$ncity), 
      g0_mu = rnorm(1),
      g0_tau = rgamma(1,1,1),
      g0 = rnorm(data_list$ncity) 
    )
  }
  return(
    switch(
      chain,
      "1" = gen_list(chain),
      "2" = gen_list(chain),
      "3" = gen_list(chain),
      "4" = gen_list(chain),
      "5" = gen_list(chain),
      "6" = gen_list(chain),
      "7" = gen_list(chain),
      "8" = gen_list(chain)
    )
  )
}

# Fit the model
my_start <- Sys.time()
null_mod <- runjags::run.jags(
  model = "./JAGS/multiSpec_Null.R",
  monitor = c("a0_mu", "a0_sd", "a0",
              "b0_mu", "b0_sd", "b0",
              "c0_mu", "c0_sd", "c0",
              "d0_mu", "d0_sd", "d0",
              "e0_mu", "e0_sd", "e0",
              "g0_mu", "g0_sd", "g0",
              "p0_mu", "p0_sd", "p1_mu", "p1_sd", "p0", "p1",
              "phi_mu", "phi_sd", "phi","x"),
  data = data_list,
  n.chains = 3,
  inits = null_inits,
  burnin = 50000,
  sample = 20000,
  adapt = 500,
  modules = "glm",
  thin = 2,
  method = "parallel"
)
my_end <- Sys.time()
saveRDS(null_mod, "./results/null_model.RDS")

#####
library(coda)
data_list <- readRDS("./data/data_list.RDS")
imperv_mod <- readRDS("./results/imperv_model.RDS")
mc <- do.call("rbind", imperv_mod$mcmc)
calculate_cpo(mm = mc, data_list = data_list)

canopy_mod <- readRDS("./results/canopy_model.RDS")
mc <- do.call("rbind", canopy_mod$mcmc)
calculate_cpo(mm = mc, data_list = data_list)

grass_mod <- readRDS("./results/grass_model.RDS")
mc <- do.call("rbind", grass_mod$mcmc)
calculate_cpo(mm = mc, data_list = data_list)

PCA_mod <- readRDS("./results/PCA_model.RDS")
mc <- do.call("rbind", PCA_mod$mcmc)
calculate_cpo(mm = mc, data_list = data_list)


# DIC stuff, but not sure it works for this model
canopy_mod <- readRDS("./results/canopy_model.RDS")
grass_mod <- readRDS("./results/grass_model.RDS")
PCA_mod <- readRDS("./results/PCA_model.RDS")
null_mod <- readRDS("./results/null_model.RDS")
runjags::extract(imperv_mod, 'dic')
runjags::extract(canopy_mod, 'dic')
runjags::extract(grass_mod, 'dic')
runjags::extract(PCA_mod, 'dic')
runjags::extract(null_mod, 'dic')
