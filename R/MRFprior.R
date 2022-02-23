## ------------------------------------------------------------------------
## Calculate normalizing constant for MRF prior for specific edge (i, j)
## ------------------------------------------------------------------------
B = function(K){
  return(data.matrix(expand.grid(replicate(K, 0:1, simplify = FALSE))))
}

mrf_C <- function(Theta, nui, B){
  ## Input
  ## Theta: graph similarity
  ## nu: edge specific prior
  ##
  mrf_C = 0

  for(i in 1:nrow(B)){
    g_ij = B[i, ,drop = F]
    mrf_C = mrf_C + exp(nui*sum(g_ij) + as.numeric(g_ij%*%Theta%*%t(g_ij))) # C
  }
  return(mrf_C)
}


## ------------------------------------------------------------------------
## update relatedness
## ------------------------------------------------------------------------
update_relatedness <- function(adj, Theta, nu, alpha, beta, my_w, B){
  # Input:
  #   adj: l x k edge inclusion
  #   Theta: initial value for graph similarity matrix
  #   nu: p vector of initial values for nu
  #   alpha: shape parameter for Theta slab
  #   beta: rate parameter for Theta slab
  #   a: first parameter for prior on nu
  #   b: second parameter for prior on nu
  #   my_w: Bernoulli prior parameter for graph similarity indicators
  # Output:
  #   Theta: updated Theta
  #   accep_theta: K x K matrix of acceptance of theta_kkprime (within-model move)
  #   accep_gamma: K x K matrix of acceptance of between-model moves
  #   within_model: within_model update 0/1

  K = ncol(Theta) # number of graph
  p = length(nu) # number of variables

  # Proposal parameters for MH steps
  alpha_prop = 1
  beta_prop = 1

  accep_gamma = accep_theta = within_model = matrix(0, K, K)

  for(k in 1:K){
    for(kprime in (k+1):K){
      ## --------------------
      ## Between model move
      ## --------------------
      if(Theta[k,kprime] == 0){
        theta_prop = rgamma(1, shape = alpha_prop, rate = beta_prop) # check which parameters alpha_prop, beta_prop
      }else{
        theta_prop = 0
      }
      Theta_prop = Theta
      Theta_prop[k, kprime] = Theta_prop[kprime, k] = theta_prop

      # log likelihood
      sum_over_edges = 0
      for(l in 1:p){ # regression p number of variables
        sum_over_edges = sum_over_edges + log(mrf_C(Theta, nu[l], B)) +
          2 * (theta_prop - Theta[k, kprime]) * adj[l, k] * adj[l, kprime] -
          log(mrf_C(Theta_prop, nu[l], B))
      }

      # MH ratio
      if(theta_prop == 0){ # 1 -> 0
        log_ar = alpha_prop*log(beta_prop) - lgamma(alpha_prop) +
          lgamma(alpha) - alpha*log(beta) -
          (alpha - alpha_prop)*log(Theta[k, kprime]) +
          (beta - beta_prop)*Theta[k, kprime] + sum_over_edges +
          log(1-my_w) - log(my_w)

      }else{ # 0 -> 1
        log_ar = alpha*log(beta) - lgamma(alpha) +
          lgamma(alpha_prop) -  alpha_prop*log(beta_prop) + # check the sign
          (alpha-alpha_prop)*log(theta_prop) -
          (beta - beta_prop)*theta_prop + sum_over_edges +
          log(my_w) - log(1-my_w)
      }

      # Accept or reject
      if(log_ar > log(runif(1))){
        Theta[k, kprime] = Theta[kprime, k] = theta_prop
        accep_gamma[k, kprime] = 1
      }

      ## --------------------
      ## within model move
      ## --------------------
      if(Theta[k, kprime] != 0){
        within_model[k, kprime] = 1

        theta_prop   = rgamma(1, shape = alpha_prop, rate = beta_prop)
        Theta_prop   = Theta
        Theta_prop[k, kprime] = Theta_prop[kprime, k] = theta_prop

        # log likelihood
        sum_over_edges = 0
        for(l in 1:p){ # regression p number of variables
          sum_over_edges = sum_over_edges + log(mrf_C(Theta, nu[l], B)) +
            2 * (theta_prop - Theta[k, kprime]) * adj[l, k] * adj[l, kprime] -
            log(mrf_C(Theta_prop, nu[l], B))
          }

        # MH ratio
        log_ar = (alpha - alpha_prop)*(log(theta_prop) - log(Theta[k, kprime])) +
          (beta- beta_prop)*(Theta[k, kprime] - theta_prop) + sum_over_edges

        # Accept or reject
        if(log_ar > log(runif(1))){
          Theta[k, kprime] = Theta[kprime, k] = theta_prop
          accep_theta[k, kprime] = 1
        }
      }
    }
  }

  return(list(Theta = Theta, accep_theta = accep_theta,
              accep_gamma = accep_gamma, within_model = within_model))
}


update_nu <- function(nu, adj, Theta, a, b, B){
  # Input:
  #   adj: l x k edge inclusion
  #   Theta: initial value for graph similarity matrix
  #   nu: p vector of initial values for nu
  #   a: first parameter for prior on nu
  #   b: second parameter for prior on nu
  #Output
  #   nu: updated nu
  #   accep_nu: 0/1 acceptance of nu

  a_prop = 2
  b_prop = 4

  accep_nu = vector("numeric", length(nu))

  for(l in 1:p){
   qu = rgamma(1, shape = a_prop, rate = b_prop)
   nu_prop = log(qu) - log(1-qu)

   log_ar = (nu_prop - nu[l])*log(a - a_prop + sum(adj[l,])) +
     log(mrf_C(Theta, nu[l],B)) +
     (a+b-a_prop-b_prop)*(log(1+exp(nu[l]))-log(1+exp(nu_prop))) -
     log(mrf_C(Theta, nu_prop, B))

   if(log_ar > log(runif(1))){
     nu[l]   = nu_prop
     accep_nu[l] = 1
     }
  }

  return(list = c(nu = nu, accep_nu = accep_nu))
}

## ------------------------------------------------------------------------
## update graph
## ------------------------------------------------------------------------

update_adj <- function(nu, adj, Theta){

  K       = nrow(Theta)
  p       = length(nu)
  prob    = matrix(NA, nrow = p, ncol = K)
  adj_new = adj
  for(l in 1:p){
  for(k in 1:K){
      tmp1 = sum(Theta[k,-k]* adj[l,-k])
      cat(sprintf("tmp1%f\n", tmp1))
      w =  exp(nu[l] + 2* tmp1)
      cat(sprintf("w%f\n", w))
      prob[l, k] = w/(1+w)
      cat(sprintf("prob[l, k]%d %d %f\n", l,k, prob[l, k]))

      adj_new[l,k] = prob[l, k] > runif(1)
    }
  }
  return(list(prob = prob, adj = adj_new))
}
