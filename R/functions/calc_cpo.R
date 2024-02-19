###############################
#
# calc_cpo
#
# Written by Mason Fidino
#
###############################

# This function calculates the cpo for each datapoint and creates
#  the summary statistic -log(sum(cpo)) for each model. For clarity,
#  I've written this to loop through each step in the MCMC chain.
#  This slows down the code considerably. This function could be made faster
#  by applying matrix algebra to the whole mcmc chain instead,
#  but requires some manipulations to the transition matrices to do so.

# This was also set up for our own analysis where we only have coyote
#  influence the detection rates of the smaller species (and not the other way
#  around). To make it for all of the species the rho detection matrix would
#  need to be updated as well as the number of rows in the matrix l.

# Arguments for this function
# --------------------------- #
# mm = model matrix of jags output as a matrix
#  e.g., mm <- as.matrix(jags_output, chains = TRUE) via coda package
#  columns of this matrix represent different parameters while rows are
#  posterior simulations

# data_list = the same data_list supplied to JAGS to fit the DCOM
#   This was generated in fit_softmax_model.R

calculate_cpo <- function(mm = NULL, data_list = NULL){
  
  # ensure mm is a matrix
  if(!is.matrix(mm)) {
    stop("The object mm must be a matrix")
  }
  nmcmc <- nrow(mm)
  
  # Assigning all the objects in data_list so I don't have to
  # index them the whole time.
  for(i in 1:length(data_list)){
    assign(names(data_list)[i], data_list[[i]])
  }
  
  psi_cov <- cbind(
    1,
    PC1,
    PC2
  )
  p_cov <- cbind(
    1,
    trapDays
  )
  lik_mat <- matrix(0, ncol = nseason * nsite, nrow = nmcmc)
  mm <- split_mcmc(mm)
  # This is the likelihood matrix, which stores the likelihood of each
  # observations based off the parameters in the model for each step
  # of the mcmc chain.


  # get the observed state 
  state_observed <- array(
    NA,
    dim = c(nsite, nseason, nrep)
  )
  max_state <- matrix(
    NA,
    ncol = nseason,
    nrow = nsite
  )
  
  state_idx <- function(xx){
    if(length(grep("NA",xx)) == 1){
      return(NA)
    }
    switch(xx,
      "1-1-1" = 1,
      "1-0-0" = 2,
      "0-1-0" = 3,
      "0-0-1" = 4,
      "1-1-0" = 5,
      "1-0-1" = 6,
      "0-1-1" = 7,
      "0-0-0" = 8,
      "NA-NA-NA" = NA
    )
  }
  for(i in 1:nsite){
    for(k in 1:nseason){
      rep_vals <- rep(NA, 4)
      for(j in 1:nrep){
        rep_vals[j] <- state_idx(paste0(y[,i,j,k], collapse = "-"))
      }
      if(all(is.na(rep_vals))){
        next
      }
      state_observed[i,k,] <- rep_vals
    }
  }
  # determine max state (most species) observed across each survey
  for(i in 1:nsite){
    for(k in 1:nseason){
      tmp_max <- apply(
        y[,i,,k],
        1,
        max
      )
      if(all(is.na(tmp_max))){
        next
      }
      max_state[i,k] <- state_idx(paste0(tmp_max, collapse = "-"))
    }
  }
  
  pb <- txtProgressBar(max = nmcmc)
  for(o in 1:nmcmc){ # going through each step of the mcmc chain
    setTxtProgressBar(pb, o)
    # occupancy
    a <- cbind(mm$a0[o,], mm$a1[o,], mm$a2[o,])
    b <- cbind(mm$b0[o,], mm$b1[o,], mm$b2[o,])
    c <- cbind(mm$c0[o,], mm$c1[o,], mm$c2[o,])
    d <- cbind(mm$d0[o,], mm$d1[o,], mm$d2[o,])
    e <- cbind(mm$e0[o,], mm$e1[o,], mm$e2[o,])
    g <- cbind(mm$g0[o,], mm$g1[o,], mm$g2[o,])
    # detection
    p <- array(
      NA,
      dim = c(nspec, ncity, 2)
    )
    p[,,1] <- mm$p0[o,,]
    p[,,2] <- mm$p1[o,,]
    # estimated community state at time t and site j\
    x <- mm$x[o,,]
    
    # occupancy linear predictor
    fox  <- t(a %*% t(psi_cov))[cbind(1:nsite, city_vec)]
    # occupancy linear predictor
    gray  <- t(b %*% t(psi_cov))[cbind(1:nsite, city_vec)]
    # occupancy linear predictor
    red  <- t(c %*% t(psi_cov))[cbind(1:nsite, city_vec)]
    # occupancy linear predictor
    foxGray <- t(d %*% t(psi_cov))[cbind(1:nsite, city_vec)]
    # occupancy linear predictor
    foxRed  <- t(e %*% t(psi_cov))[cbind(1:nsite, city_vec)]
    # occupancy linear predictor
    grayRed  <- t(g %*% t(psi_cov))[cbind(1:nsite, city_vec)]
    
    phi <- mm$phi[o,,]
    
    # This is the numerator of the softmax function for each of the 8 community
    # states a site can be in.
    psi <- array(0, dim = c(nsite, 8, nseason))
    psi[,1,1] <- exp( fox + gray + red + foxGray + foxRed + grayRed ) #------|ABC
    psi[,2,1] <- exp( fox ) #------------------------------------------------|A
    psi[,3,1] <- exp( gray ) #-----------------------------------------------|B
    psi[,4,1] <- exp( red ) #------------------------------------------------|C
    psi[,5,1] <- exp( fox + gray + foxGray ) #-------------------------------|AB
    psi[,6,1] <- exp( fox + red + foxRed ) #---------------------------------|BC
    psi[,7,1] <- exp( gray + red + grayRed ) #-------------------------------|AC
    psi[,8,1] <- 1 #---------------------------------------------------------|U
    # fill in the remaining seasons, with the added
    #  autologistic term as needed.
    for(j in 2:nseason){
      tmp_x <- x[,j-1]
      tmp <- Xcat[tmp_x,]
      tmp <- tmp %*% phi
      tmp <- tmp[cbind(1:nsite, city_vec)]
      
      psi[,1,j] <- exp( 
        fox + gray + red + foxGray + foxRed + grayRed +
          ifelse(tmp_x == 1, tmp, 0) 
      ) #------|ABC
      psi[,2,j] <- exp(
        fox + ifelse(tmp_x == 2, tmp, 0)
      )
      psi[,3,j] <- exp(
        gray + ifelse(tmp_x == 3, tmp, 0)
      )
      psi[,4,j] <- exp(
        red + ifelse(tmp_x == 4, tmp, 0)
      )
      psi[,5,j] <- exp(
        fox + gray + foxGray +  ifelse(tmp_x == 5, tmp, 0)
      )
      psi[,6,j] <- exp(
        fox + red + foxRed + ifelse(tmp_x == 6, tmp, 0)
      )
      psi[,7,j] <- exp(
        gray + red + grayRed + ifelse(tmp_x == 7, tmp, 0)
      )
      psi[,8,j] <- 1
    }
    # convert to probability
    for(i in 1:nsite){
      for(j in 1:nseason){
        psi[i,,j] <- psi[i,,j] / sum(psi[i,,j])
      }
    }
    
    # calculate survey specific detection probability for each species
    # species x site x sample week x season
    rho <- array(
      NA,
      dim = dim(trapDays)
    )
    for(i in 1:nspec){
      for(j in 1:nrep){
        for(k in 1:nseason){
          # assuming here p0 is mcmc x species x city
          #  and trapDays is species x site x rep x season
          rho[i,1:nsite,j,k] <-  mm$p0[o,i,city_vec] + mm$p1[o,i,city_vec] * trapDays[i,,j,k]
        }
      }
    }
    # convert to probability
    rho <- plogis(rho)
    # now calculate the probability of detecting each state.
    rdm <- array(
      NA,
      dim = c(nsite, nrep, nseason, 8,8)
    )
    for(j in 1:nrep){
      for(k in 1:nseason){
        # assuming here p0 is mcmc x species x city
        #  and trapDays is species x site x rep x season
        # True state = 1
        rdm[,j,k,1,1] <- rho[1,,j,k] * rho[2,,j,k] * rho[3,,j,k]
        rdm[,j,k,2,1] <- rho[1,,j,k] * (1-rho[2,,j,k]) *(1- rho[3,,j,k])
        rdm[,j,k,3,1] <- (1-rho[1,,j,k]) * rho[2,,j,k] * (1-rho[3,,j,k])
        rdm[,j,k,4,1] <- (1-rho[1,,j,k]) * (1-rho[2,,j,k]) * rho[3,,j,k]
        rdm[,j,k,5,1] <- rho[1,,j,k] * rho[2,,j,k] * (1-rho[3,,j,k])
        rdm[,j,k,6,1] <- rho[1,,j,k] * (1-rho[2,,j,k]) * rho[3,,j,k]
        rdm[,j,k,7,1] <- (1-rho[1,,j,k]) * rho[2,,j,k] * rho[3,,j,k]
        rdm[,j,k,8,1] <- (1-rho[1,,j,k]) * (1-rho[2,,j,k]) * (1-rho[3,,j,k])
        # True state = 2
        rdm[,j,k,1,2] <- 0
        rdm[,j,k,2,2] <- rho[1,,j,k]
        rdm[,j,k,3,2] <- 0
        rdm[,j,k,4,2] <- 0
        rdm[,j,k,5,2] <- 0
        rdm[,j,k,6,2] <- 0
        rdm[,j,k,7,2] <- 0
        rdm[,j,k,8,2] <- (1-rho[1,,j,k])
        # True state = 3
        rdm[,j,k,1,3] <- 0
        rdm[,j,k,2,3] <- 0
        rdm[,j,k,3,3] <- rho[2,,j,k]
        rdm[,j,k,4,3] <- 0
        rdm[,j,k,5,3] <- 0
        rdm[,j,k,6,3] <- 0
        rdm[,j,k,7,3] <- 0
        rdm[,j,k,8,3] <- (1-rho[2,,j,k])
        # True state = 4
        rdm[,j,k,1,4] <- 0
        rdm[,j,k,2,4] <- 0
        rdm[,j,k,3,4] <- 0
        rdm[,j,k,4,4] <- rho[3,,j,k]
        rdm[,j,k,5,4] <- 0
        rdm[,j,k,6,4] <- 0
        rdm[,j,k,7,4] <- 0
        rdm[,j,k,8,4] <- (1-rho[3,,j,k])
        # True state = 5
        rdm[,j,k,1,5] <- 0
        rdm[,j,k,2,5] <- rho[1,,j,k] * (1 - rho[2,,j,k])
        rdm[,j,k,3,5] <- rho[2,,j,k] * (1 - rho[1,,j,k])
        rdm[,j,k,4,5] <- 0
        rdm[,j,k,5,5] <- rho[1,,j,k] * rho[2,,j,k]
        rdm[,j,k,6,5] <- 0
        rdm[,j,k,7,5] <- 0
        rdm[,j,k,8,5] <- (1-rho[1,,j,k]) * (1 - rho[2,,j,k])
        # True state = 6
        rdm[,j,k,1,6] <- 0
        rdm[,j,k,2,6] <- rho[1,,j,k] * (1 - rho[3,,j,k])
        rdm[,j,k,3,6] <- 0
        rdm[,j,k,4,6] <- (1-rho[1,,j,k]) * rho[3,,j,k]
        rdm[,j,k,5,6] <- 0
        rdm[,j,k,6,6] <- rho[1,,j,k] * rho[3,,j,k]
        rdm[,j,k,7,6] <- 0
        rdm[,j,k,8,6] <- (1-rho[1,,j,k]) * (1 - rho[3,,j,k])
        # True state = 7
        rdm[,j,k,1,7] <- 0
        rdm[,j,k,2,7] <- 0
        rdm[,j,k,3,7] <- rho[2,,j,k] * (1 - rho[3,,j,k])
        rdm[,j,k,4,7] <- (1-rho[2,,j,k]) *  rho[3,,j,k]
        rdm[,j,k,5,7] <- 0
        rdm[,j,k,6,7] <- rho[2,,j,k] * rho[3,,j,k]
        rdm[,j,k,7,7] <- 0
        rdm[,j,k,8,7] <- (1-rho[2,,j,k]) * (1 - rho[3,,j,k])
        # True state = 8
        rdm[,j,k,1,8] <- 0
        rdm[,j,k,2,8] <- 0
        rdm[,j,k,3,8] <- 0
        rdm[,j,k,4,8] <- 0
        rdm[,j,k,5,8] <- 0
        rdm[,j,k,6,8] <- 0
        rdm[,j,k,7,8] <- 0
        rdm[,j,k,8,8] <- 1
        
      }
    }
    
    # okay, so given the observed data, it is possible to be in a variety of
    #  states. As such, we need to incorporate those into the likelihood.
    lik <- matrix(0, ncol = nseason, nrow = nsite)
    for(i in 1:nsite){
      for(j in 1:nseason){
      my_obs <- state_observed[i,j,]
      max_obs <- max_state[i,j]
      if(is.na(max_obs)){
        next
      }
      if(max_obs == 1){
        # if state 8, that is all you could be in.
        tmp_rhos <- rdm[i, , j, my_obs , max_obs]
        # just need the diagonal
        tmp_rhos <- diag(tmp_rhos)
        lik[i,j] <- psi[i,max_obs, j] * prod(tmp_rhos)
        
      }
      if(max_obs == 2){
        # if state 2, you could be in state 1,2,5, or 6.
        
        # state 1
        lik1 <- psi[i,1, j] * prod(diag(rdm[i,,j,my_obs, 1]))
        # state 2
        lik2 <- psi[i,2, j] * prod(diag(rdm[i,,j,my_obs, 2]))
        # state 5
        lik3 <- psi[i,5, j] * prod(diag(rdm[i,,j,my_obs, 5]))
        # state 6
        lik4 <- psi[i,6, j] * prod(diag(rdm[i,,j,my_obs, 6]))
        
        # likelihood is the sum of these
        lik[i,j] <- lik1 + lik2 + lik3 + lik4
      }
      if(max_obs == 3){
        # if state 3, you could be in state 1,2,5, or 7.
        
        # state 1
        lik1 <- psi[i,1, j] * prod(diag(rdm[i,,j,my_obs, 1]))
        # state 2
        lik2 <- psi[i,2, j] * prod(diag(rdm[i,,j,my_obs, 2]))
        # state 5
        lik3 <- psi[i,5, j] * prod(diag(rdm[i,,j,my_obs, 5]))
        # state 7
        lik4 <- psi[i,7, j] * prod(diag(rdm[i,,j,my_obs, 7]))
        
        # likelihood is the sum of these
        lik[i,j] <- lik1 + lik2 + lik3 + lik4
        
      }
      if(max_obs == 4){
        # if state 4, you could be in state 1,2,6, or 7.
        
        # state 1
        lik1 <- psi[i,1, j] * prod(diag(rdm[i,,j,my_obs, 1]))
        # state 2
        lik2 <- psi[i,2, j] * prod(diag(rdm[i,,j,my_obs, 2]))
        # state 5
        lik3 <- psi[i,6, j] * prod(diag(rdm[i,,j,my_obs, 6]))
        # state 6
        lik4 <- psi[i,7, j] * prod(diag(rdm[i,,j,my_obs, 7]))
        
        # likelihood is the sum of these
        lik[i,j] <- lik1 + lik2 + lik3 + lik4
      }
      if(max_obs == 5){
        # if state 5, you could be in state 1 or 5.
        
        # state 1
        lik1 <- psi[i,1, j] * prod(diag(rdm[i,,j,my_obs, 1]))
        # state 5
        lik2 <- psi[i,5, j] * prod(diag(rdm[i,,j,my_obs, 5]))
        
        # likelihood is the sum of these
        lik[i,j] <- lik1 + lik2
      }
      if(max_obs == 6){
        # if state 6, you could be in state 1 or 6.
        
        # state 1
        lik1 <- psi[i,1, j] * prod(diag(rdm[i,,j,my_obs, 1]))
        # state 6
        lik2 <- psi[i,6, j] * prod(diag(rdm[i,,j,my_obs, 6]))
        
        # likelihood is the sum of these
        lik[i,j] <- lik1 + lik2
        
        
      }
      if(max_obs == 7){
        # if state 7, you could be in state 1 or 7.
        
        # state 1
        lik1 <- psi[i,1, j] * prod(diag(rdm[i,,j,my_obs, 1]))
        # state 7
        lik2 <- psi[i,7, j] * prod(diag(rdm[i,,j,my_obs, 7]))
        
        # likelihood is the sum of these
        lik[i,j] <- lik1 + lik2
        
      }
      if(max_obs == 8){
        # if state 8, you were either not detected for any of the other states
        #  or you simply were not there.
        # state 1
        lik1 <- psi[i,1, j] * prod(rdm[i,,j,8, 1])
        # state 2
        lik2 <- psi[i,2, j] * prod(rdm[i,,j,8, 2])
        # state 3
        lik3 <- psi[i,3, j] * prod(rdm[i,,j,8, 3])
        # state 4
        lik4 <- psi[i,4, j] * prod(rdm[i,,j,8, 4])
        # state 5
        lik5 <- psi[i,5, j] * prod(rdm[i,,j,8, 5])
        # state 6
        lik6 <- psi[i,6, j] * prod(rdm[i,,j,8, 6])
        # state 7
        lik7 <- psi[i,7, j] * prod(rdm[i,,j,8, 7])
        # state 8
        lik8 <- psi[i,8,j]
        lik[i,j] <- lik1 + lik2 + lik3 + lik4 + lik5 +
          lik6 + lik7 + lik8
      }
      
      }
    }
    lik_mat[o,] <- as.numeric(lik) # put lik into likelihood matrix
  }
  
  # these are sites by seasons that sampling did not happen
  # we do not want them to contribute towards the likelihood
  to_zero <- rep(0, ncol(lik_mat))
  
  for(i in 1:length(to_zero)){
    to_zero[i] <- sum(lik_mat[,i]==0)
  }
  # these are all unobserved sites
  to_go <- which(to_zero == nmcmc)
  # remove them from the likelihood matrix
  lik_mat <- lik_mat[,-to_go]
  
  # calculate cpo
  CPO <--sum(log(nmcmc /(apply(1/lik_mat, 2, sum))))
  # return cpo
  return(CPO)
}
