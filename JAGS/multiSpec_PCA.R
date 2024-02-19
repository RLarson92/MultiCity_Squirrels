# For this model, the occupancy state matrix Xcat should be of the following form, which matches
# the order presented in Appendix A:
# Xcat <- matrix(c(1, 1, 1,
#                  1, 0, 0,
#                  0, 1, 0,
#                  0, 0, 1,
#                  1, 1, 0,
#                  1, 0, 1,
#                  0, 1, 1,
#                  0, 0, 0), ncol = 3, byrow = TRUE)

model{
  ## PRIORS
  # latent state priors
    a0_mu ~ dnorm(0,0.25)
    a0_tau ~ dgamma(1,1)
    a0_sd <- 1 / sqrt(a0_tau)
    a1_mu ~ dnorm(0,0.25)
    a1_tau ~ dgamma(1,1)
    a1_sd <- 1 / sqrt(a1_tau)
    a2_mu ~ dnorm(0,0.25)
    a2_tau ~ dgamma(1,1)
    a2_sd <- 1 / sqrt(a1_tau)
    # gray squirrel
    b0_mu ~ dnorm(0,0.25)
    b0_tau ~ dgamma(1,1)
    b0_sd <- 1 / sqrt(b0_tau)
    b1_mu ~ dnorm(0,0.25)
    b1_tau ~ dgamma(1,1)
    b1_sd <- 1 / sqrt(b1_tau)
    b2_mu ~ dnorm(0,0.25)
    b2_tau ~ dgamma(1,1)
    b2_sd <- 1 / sqrt(b1_tau)
    # red squirrel
    c0_mu ~ dnorm(0,0.25)
    c0_tau ~ dgamma(1,1)
    c0_sd <- 1 / sqrt(c0_tau)
    c1_mu ~ dnorm(0,0.25)
    c1_tau ~ dgamma(1,1)
    c1_sd <- 1 / sqrt(c1_tau)
    c2_mu ~ dnorm(0,0.25)
    c2_tau ~ dgamma(1,1)
    c2_sd <- 1 / sqrt(c1_tau)
    # fox/gray
    d0_mu ~ dnorm(0,0.25)
    d0_tau ~ dgamma(1,1)
    d0_sd <- 1 / sqrt(d0_tau)
    d1_mu ~ dnorm(0,0.25)
    d1_tau ~ dgamma(1,1)
    d1_sd <- 1 / sqrt(d1_tau)
    d2_mu ~ dnorm(0,0.25)
    d2_tau ~ dgamma(1,1)
    d2_sd <- 1 / sqrt(d1_tau)
    # fox/red
    e0_mu ~ dnorm(0,0.25)
    e0_tau ~ dgamma(1,1)
    e0_sd <- 1 / sqrt(e0_tau)
    e1_mu ~ dnorm(0,0.25)
    e1_tau ~ dgamma(1,1)
    e1_sd <- 1 / sqrt(e1_tau)
    e2_mu ~ dnorm(0,0.25)
    e2_tau ~ dgamma(1,1)
    e2_sd <- 1 / sqrt(e1_tau)
    # gray/red
    g0_mu ~ dnorm(0,0.25)
    g0_tau ~ dgamma(1,1)
    g0_sd <- 1 / sqrt(g0_tau)
    g1_mu ~ dnorm(0,0.25)
    g1_tau ~ dgamma(1,1)
    g1_sd <- 1 / sqrt(g1_tau)
    g2_mu ~ dnorm(0,0.25)
    g2_tau ~ dgamma(1,1)
    g2_sd <- 1 / sqrt(e1_tau)
    
    for(z in 1:nspec){
      p0_mu[z] ~ dnorm(0, 0.25)
      p0_tau[z] ~ dgamma(1,1)
      p0_sd[z] <- 1/sqrt(p0_tau[z])
      p1_mu[z] ~ dnorm(0, 0.25)
      p1_tau[z] ~ dgamma(1,1)
      p1_sd[z] <- 1/sqrt(p1_tau[z])
      phi_mu[z] ~ dnorm(0, 0.25)
      phi_tau[z] ~ dgamma(1,1)
      phi_sd[z] <- 1/sqrt(phi_tau[z])
    }
    
  for(m in 1:ncity){ 
    # fox squirrel
    a0[m] ~ dnorm(a0_mu, a0_tau)
    a1[m] ~ dnorm(a1_mu, a1_tau)
    a2[m] ~ dnorm(a2_mu, a2_tau)
    # gray squirrel
    b0[m] ~ dnorm(b0_mu, b0_tau)
    b1[m] ~ dnorm(b1_mu, b1_tau)
    b2[m] ~ dnorm(b2_mu, b2_tau)
    # red squirrel
    c0[m] ~ dnorm(c0_mu, c0_tau)
    c1[m] ~ dnorm(c1_mu, c1_tau)
    c2[m] ~ dnorm(c2_mu, c2_tau)
    # fox/gray
    d0[m] ~ dnorm(d0_mu, d0_tau)
    d1[m] ~ dnorm(d1_mu, d1_tau)
    d2[m] ~ dnorm(d2_mu, d2_tau)
    # fox/red
    e0[m] ~ dnorm(e0_mu, e0_tau)
    e1[m] ~ dnorm(e1_mu, e1_tau)
    e2[m] ~ dnorm(e2_mu, e2_tau)
    # gray/red
    g0[m] ~ dnorm(g0_mu, g0_tau)
    g1[m] ~ dnorm(g1_mu, g1_tau)
    g2[m] ~ dnorm(g2_mu, g2_tau)
    
    # detection priors
    for(z in 1:nspec){
      p0[z,m] ~ dnorm(p0_mu[z],p0_tau[z])
      p1[z,m] ~ dnorm(p1_mu[z],p1_tau[z])
      phi[z,m] ~ dnorm(phi_mu[z],phi_tau[z])
    }
  }
  
  ## MODELS
  # loop over sites
  for(j in 1:nsite){
    ### SEASON 1
    ## STATE MODEL - multivariate categorical
    # natural parameters
    # 1 = fox, 2 = gray, 3 = red
    f1[j,1] <- a0[city_vec[j]] + a1[city_vec[j]]*PC1[j] + a2[city_vec[j]]*PC2[j]
    f2[j,1] <- b0[city_vec[j]] + b1[city_vec[j]]*PC1[j] + b2[city_vec[j]]*PC2[j]
    f3[j,1] <- c0[city_vec[j]] + c1[city_vec[j]]*PC1[j] + c2[city_vec[j]]*PC2[j]
    f12[j,1] <- d0[city_vec[j]] + d1[city_vec[j]]*PC1[j] + d2[city_vec[j]]*PC2[j]
    f13[j,1] <- e0[city_vec[j]] + e1[city_vec[j]]*PC1[j] + e2[city_vec[j]]*PC2[j]
    f23[j,1] <- g0[city_vec[j]] + g1[city_vec[j]]*PC1[j] + g2[city_vec[j]]*PC2[j]
    # Psi gives latent state 'category' for each site
    # for Psi[i,j,k], i = latent occupancy category (e.g., 100, 101, 010), j = individual site, 
    # k = individual season
    # the math below calculates the odds of each latent occupancy category for each iteration
    # all species present
    Psi[1,j,1] <- exp(f1[j,1] + f2[j,1] + f3[j,1] + f12[j,1] + f13[j,1] + f23[j,1])
    # only fox squirrels
    Psi[2,j,1] <- exp(f1[j,1])
    # only gray squirrels
    Psi[3,j,1] <- exp(f2[j,1])
    # only red squirrels
    Psi[4,j,1] <- exp(f3[j,1])
    # gray & fox squirrels
    Psi[5,j,1] <- exp(f1[j,1] + f2[j,1] + f12[j,1])
    # fox & red squirrels
    Psi[6,j,1] <- exp(f1[j,1] + f3[j,1] + f13[j,1])
    # gray & red squirrels
    Psi[7,j,1] <- exp(f2[j,1] + f3[j,1] + f23[j,1])
    # no squirrels present
    Psi[8,j,1] <- 1
    # model latent occupancy state as categorical random variable derived from categories in Psi
    x[j,1] ~ dcat(Psi[,j,1])
    
    ## OBSERVATION MODEL
    # loop over species
    for(i in 1:nspec){
      # loop over weeks
      for(k in 1:nrep){
        # observation modeled as function of trap days (# days camera was active per week per season)
        logit(p[i,j,k,1]) <- p0[i,city_vec[j]] + p1[i,city_vec[j]]*trapDays[i,j,k,1]
        # occupancy modeled as state * observation
        # note: Xcat is a matrix of 2^n rows and n columns, filled w/ 0 & 1, where n = # of species
        # (see above, y is the 0/1/NA array of detections, indexed by species, site, week, season)
        y[i,j,k,1] ~ dbern(Xcat[x[j,1],i]*p[i,j,k,1])
      }
    }
  
    ### SEASON >1
    ## STATE MODEL - multivariate categorical
    # natural parameters
    # for season >1, single species natural parameters have an additional parameter w/ coefficient
    # phi that is multiplied by that site's occupancy state from the previous season
    for(t in 2:nseason){
      f1[j,t] <- a0[city_vec[j]] + a1[city_vec[j]]*PC1[j] + a2[city_vec[j]]*PC2[j] + phi[1,city_vec[j]]*Xcat[x[j,t-1],1]
      f2[j,t] <- b0[city_vec[j]] + b1[city_vec[j]]*PC1[j] + b2[city_vec[j]]*PC2[j] + phi[2,city_vec[j]]*Xcat[x[j,t-1],2]
      f3[j,t] <- c0[city_vec[j]] + c1[city_vec[j]]*PC1[j] + c2[city_vec[j]]*PC2[j] + phi[3,city_vec[j]]*Xcat[x[j,t-1],3]
      f12[j,t] <- d0[city_vec[j]] + d1[city_vec[j]]*PC1[j] + d2[city_vec[j]]*PC2[j]
      f13[j,t] <- e0[city_vec[j]] + e1[city_vec[j]]*PC1[j] + e2[city_vec[j]]*PC2[j]
      f23[j,t] <- g0[city_vec[j]] + g1[city_vec[j]]*PC1[j] + g2[city_vec[j]]*PC2[j]
      Psi[1,j,t] <- exp(f1[j,t] + f2[j,t] + f3[j,t] + f12[j,t] + f13[j,t] + f23[j,t])
      Psi[2,j,t] <- exp(f1[j,t])
      Psi[3,j,t] <- exp(f2[j,t])
      Psi[4,j,t] <- exp(f3[j,t])
      Psi[5,j,t] <- exp(f1[j,t] + f2[j,t] + f12[j,t])
      Psi[6,j,t] <- exp(f1[j,t] + f3[j,t] + f13[j,t])
      Psi[7,j,t] <- exp(f2[j,t] + f3[j,t] + f23[j,t])
      Psi[8,j,t] <- 1
      # model latent occupancy state as categorical random variable derived from categories in Psi
      x[j,t] ~ dcat(Psi[,j,t])
      
      ## OBSERVATION MODEL
      # loop over species
      for(i in 1:nspec){
        # loop over weeks
        for(k in 1:nrep){
          logit(p[i,j,k,t]) <- p0[i,city_vec[j]] + p1[i,city_vec[j]]*trapDays[i,j,k,t]
          y[i,j,k,t] ~ dbern(Xcat[x[j,t],i]*p[i,j,k,t])
        }
      }
    }
  }
}
