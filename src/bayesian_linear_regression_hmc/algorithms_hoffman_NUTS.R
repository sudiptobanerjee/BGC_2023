
# Algorithms in the Hoffman & Gelman (2014)


## Algorithm 1: HMC
##
## theta: model variable; locaiton
## r: auxiliary momentum variable
## epsilon: step size
## L: number of steps
## LL_log: log likelihood
## grad_LL_log: gradient of the log likelihood with respect to theta
## 0.5*(r dot r): kinetic energy:
## log joint density: negative energy

log_joint_density <- function(theta, r, LL_log){
  return(LL_log(theta) - sum(r^2)/2)
}

### Leapfrog integrator: update position theta and momentum r
Leapfrog_HMC <- function (theta, r, epsilon, LL_log, grad_LL_log){
  r_tilde = r + (epsilon/2) * grad_LL_log(theta,a = 2, b = 0.1)
  theta_tilde = theta + epsilon * r_tilde
  r_tilde = r_tilde + (epsilon / 2 ) * grad_LL_log(theta_tilde, a = 2, b = 0.1)
  return(list(theta_tilde = theta_tilde, r_tilde = r_tilde))
}

HMCOne <- function (current_theta, epsilon, L, LL_log, grad_LL_log){
  #theta=theta0
  #resample momentum from standard multivariate normal; Gibbs smapling uodate
  current_r = rnorm(length(current_theta),0,1)
  theta_p = current_theta
  r_p = current_r
  theta = current_theta
  
  # apply L leapfrog updates to position theta and momentum r
  # generate a proposal pair (theta_tilde, r_tilde)
  for (i in 1:L){
    #i=1
    list1 <- Leapfrog_HMC(theta_p, r_p, epsilon, LL_log, grad_LL_log)
    theta_p <- list1$theta_tilde
    r_p <- list1$r_tilde
  }
  
  #accept or reject the proposal: MH
  log_prop <- log_joint_density(theta_p, r_p, LL_log)
  log_current <- log_joint_density(current_theta, current_r, LL_log)
  accept_prob <- min(1,exp(log_prop - log_current))
  
  if (accept_prob == 0|is.nan(accept_prob)| is.na(accept_prob)){
    accept_prob <- runif(1,0,0.0001)}
  
  if (accept_prob > runif(1)){theta = theta_p}
  
  return(theta)
}

#M: Number of samples
HMC <- function (theta_0, epsilon, L, LL_log, grad_LL_log, M){
  
  theta_list = matrix(c(0), nrow = M, ncol = length(theta_0))
  theta_list[1,] <- theta_0
  
  for (m in c(2:M)){
    #m=4
    current_theta <- theta_list[m-1,]
    theta_list[m,] <- HMCOne(current_theta, epsilon, L, LL_log, grad_LL_log) 
  }
  
  theta_list[]
  
}

## Algorithm 2: Naive No-U-Turn Sampler

Leapfrog <- function (theta, r, epsilon, grad_LL_log){
  r_tilde = r + (epsilon/2) * grad_LL_log(theta,a = 2, b = 0.1)
  theta_tilde = theta + epsilon * r_tilde
  r_tilde = r_tilde + (epsilon / 2 ) * grad_LL_log(theta_tilde, a = 2, b = 0.1)
  return(list(theta_tilde = theta_tilde, r_tilde = r_tilde))
}

BuildTree <- function (theta, r, u, v, j, epsilon, sigma_max = 1000, LL_log, grad_LL_log){
  if (j == 0){
    
    list1 <- Leapfrog (theta, r, epsilon = v * epsilon, grad_LL_log)
    theta_p = list1$theta_tilde
    r_p = list1$r_tilde
    
    log_prop = log_joint_density(theta_p, r_p, LL_log)
    mar_prob = exp(log_prop)
    if (mar_prob <= 0| is.nan(mar_prob)| is.na(mar_prob)){mar_prob = runif(1,0,0.001)}
    
    if (u <= mar_prob){C_p = list(theta=theta_p, r=r_p)}else{
      C_p = list(theta=theta, r=r)
    }
    
    s_criteria = log(u) - sigma_max # recommend: sigma_max = 1000
    s_p = ifelse(log_prop > s_criteria, 1, 0)
    
    return(list(theta_b = theta_p, r_b = r_p, theta_f = theta_p, r_f = r_p, C_p = C_p, s_p = s_p))}
  else{
    #j=1
    list2 <- BuildTree(theta, r, u, v, j-1, epsilon, sigma_max = 1000, LL_log, grad_LL_log)
    theta_b = list2$theta_b
    r_b = list2$r_b
    theta_f = list2$theta_f
    r_f = list2$r_f
    C_p = list2$C_p
    s_p = list2$s_p
    
    if (v == -1){
      list3 <- BuildTree(theta_b, r_b, u, v, j-1, epsilon, sigma_max = 1000, LL_log, grad_LL_log)
      theta_b = list3$theta_b
      r_b = list3$r_b}else{
        list3 <- BuildTree(theta_f, r_f, u, v, j-1, epsilon, sigma_max = 1000, LL_log, grad_LL_log)
        theta_f = list3$theta_f
        r_f = list3$r_f
      }
    
    C_p2 = list3$C_p
    s_p2 = list3$s_p
    
    ip1<- sum((theta_f - theta_b)* r_b)
    sc1 = ifelse(ip1 >= 0, 1, 0)
    ip2<- sum((theta_f - theta_b)* r_f)
    sc2 = ifelse(ip2 >= 0, 1, 0)
    
    s_p =  s_p2 *sc1 * sc2
  

    C_p$theta <- unique(rbind(C_p$theta, C_p2$theta))
    C_p$r <- unique(rbind(C_p$r, C_p2$r))
    
    return(list(theta_b = theta_b, r_b = r_b, theta_f = theta_f, r_f = r_f,
                C_p = C_p, s_p = s_p))
  }
}

## u: slice variable
## C: candidate position-momentum states
## B: set of allposiiton-momentum states that leapfrog trace out
NaiveNUTOne <- function (current_theta, epsilon, LL_log, grad_LL_log, sigma_max = 1000){
  
  current_r <-rnorm(length(current_theta),0,1)
  
  log_current <- log_joint_density(current_theta, current_r, LL_log)
  joint_current <- exp(log_current)
  if (joint_current == 0| is.nan(joint_current)| is.na(joint_current)){
    u <- 0 }else{
    u <- runif(1, 0, joint_current)}
  
  theta_b = current_theta
  theta_f = current_theta
  r_b = current_r
  r_f = current_r
  j = 0 # current height of tree
  C = list(theta = current_theta, r = current_r)
  s = 1
  
  while (s == 1) {
    v_j = sample(c(-1,1),1)
    
    if (v_j == -1){
      listb <- BuildTree(theta_b, r_b, u, v = v_j, j, epsilon, sigma_max = 1000, LL_log, grad_LL_log)
      theta_b<- listb$theta_b
      r_b <- listb$r_b
    }else{
      listb <- BuildTree(theta_f, r_f, u, v= v_j, j, epsilon, sigma_max = 1000, LL_log, grad_LL_log)
      theta_f<- listb$theta_f
      r_f <- listb$r_f
    }
    
    C_p <- listb$C_p
    s_p <- listb$s_p
    
    if (s_p == 1 |!is.nan(s_p)){
      C$theta <- rbind(C$theta,C_p$theta)
      C$r <- rbind(C$r, C_p$r)
    }
    
    ip1<- sum((theta_f - theta_b)* r_b)
    sc1 = ifelse(ip1 >= 0, 1, 0)
    ip2<- sum((theta_f - theta_b)* r_f)
    sc2 = ifelse(ip2 >= 0, 1, 0)
    
    s = s_p * sc1 * sc2
    
    if(is.na(s)){s = 0}
    
    j = j + 1
  }
  
  theta = C$theta[sample(1:nrow(C$theta),1),]
  return(list(theta=theta))
}

NaiveNUT <- function(theta_0, epsilon, LL_log, grad_LL_log, M, sigma_max = 1000){
  theta_list = matrix(c(0), nrow = M, ncol = length(theta_0))
  theta_list[1,] <- theta_0
  
  for (m in c(2:M)){
    #m=2
    current_theta <- theta_list[m-1,]
    theta_list[m,] <- NaiveNUTOne(current_theta, epsilon, LL_log, grad_LL_log, sigma_max)$theta
  }
  
  theta_list[]
}


## Algorithm 3: Efficient No-U-Turn Sampler

BuildTree2 <- function (theta, r, u, v, j, epsilon, sigma_max = 1000, LL_log, grad_LL_log){
  if (j == 0){
    
    list1 <- Leapfrog (theta, r, epsilon = v * epsilon, grad_LL_log)
    theta_p = list1$theta_tilde
    r_p = list1$r_tilde
    
    log_prop = log_joint_density(theta_p, r_p, LL_log)
    n_p = ifelse(u <= exp(log_prop),1,0)
    
    s_criteria = log(u) - sigma_max # recommend: sigma_max = 1000
    s_p = ifelse(log_prop > s_criteria, 1, 0)
    
    return(list(theta_b = theta_p, r_b = r_p, theta_f = theta_p, r_f = r_p, 
                theta_p = theta_p, n_p = n_p, s_p = s_p))}
  else{
    #j=1
    list2 <- BuildTree2(theta, r, u, v, j-1, epsilon, sigma_max = 1000, LL_log, grad_LL_log)
    theta_b = list2$theta_b
    r_b = list2$r_b
    theta_f = list2$theta_f
    r_f = list2$r_f
    theta_p = list2$theta_p
    n_p = list2$n_p
    s_p = list2$s_p
    
    if (s_p == 1){
      if (v == -1){
        list3 <- BuildTree2(theta_b, r_b, u, v, j-1, epsilon, sigma_max = 1000, LL_log, grad_LL_log)
        theta_b = list3$theta_b
        r_b = list3$r_b}
      else{
        list3 <- BuildTree2(list2$theta_f, list2$r_f, u, v, j-1, epsilon, sigma_max = 1000, LL_log, grad_LL_log)
        theta_f = list3$theta_f
        r_f = list3$r_f
      }
      
      theta_p2 = list3$theta_p
      n_p2 = list3$n_p
      s_p2 = list3$s_p
    
      n_cri <- n_p2/(n_p+n_p2)
      
      if (runif(1) < n_cri ){
        theta_p = theta_p2
      }
      
      ip1<- sum((theta_f - theta_b)* r_b)
      sc1 = ifelse(ip1 >= 0, 1, 0)
      ip2<- sum((theta_f - theta_b)* r_f)         
      sc2 = ifelse(ip2 >= 0, 1, 0)
      
      s_p = s_p2 *sc1 * sc2
      n_p = n_p + n_p2
      
      #if(is.na(s_p)| is.nan(s_p)){s_p = 0}
    }
    
    return(list(theta_b = theta_b, r_b = r_b, theta_f = theta_f, r_f = r_f,
                theta_p = theta_p, n_p = n_p, s_p = s_p))
  }
}

EffNUTOne <- function (current_theta, epsilon, LL_log, grad_LL_log, sigma_max = 1000){
  
  current_r <-rnorm(length(current_theta),0,1)
  
  log_current <- log_joint_density(current_theta, current_r, LL_log)
  mar_prob_current <- exp(log_current)
  
  if (mar_prob_current == 0| is.nan(mar_prob_current) |is.na(mar_prob_current)){
    u <- 0}else{
    u <- runif(1, 0, mar_prob_current)}
  
  theta_b = current_theta
  theta_f = current_theta
  r_b = current_r
  r_f = current_r
  j = 0
  #theta = current_theta
  n = 1
  s = 1
  
  while (s == 1) {
    v_j = sample(c(-1,1),1)
    
    if (v_j == -1){
      listb <- BuildTree2(theta_b, r_b, u, v = v_j, j, epsilon, sigma_max = 1000, LL_log, grad_LL_log)
      theta_b<- listb$theta_b
      r_b <- listb$r_b
    }else{
      listb <- BuildTree2(theta_f, r_f, u, v= v_j, j, epsilon, sigma_max = 1000, LL_log, grad_LL_log)
      theta_f<- listb$theta_f
      r_f <- listb$r_f
    }
    
    n_p <- listb$n_p
    s_p <- listb$s_p
    theta_p <- listb$theta_p
    
    #if (is.nan(n_p) | is.na(n_p)){n_p <- 0}
    #if (is.nan(s_p) | is.na(s_p)){s_p <- 0}
    
    #alpha = min(1, n_p/n)
    #if (s_p == 1 && runif(1) < alpha){return(theta_p)}else{return(current_theta)}
    
    if (s_p == 1){alpha = min(1, n_p/n)
      if (runif(1) < alpha){return(theta_p)}else{return(current_theta)}}
    
    n = n + n_p
    ip1<- sum((theta_f - theta_b)* r_b)
    sc1 = ifelse(ip1 >= 0, 1, 0)
    ip2<- sum((theta_f - theta_b)* r_f)
    sc2 = ifelse(ip2 >= 0, 1, 0)
    
    s = s_p * sc1 * sc2
    
    j = j + 1
  }
  
}

EffNUT <- function(theta_0, epsilon, LL_log, grad_LL_log, M, sigma_max = 1000){
  theta_list = matrix(c(0), nrow = M, ncol = length(theta_0))
  theta_list[1,] <- theta_0
  
  for (m in c(2:M)){
    #m=3
    current_theta <- theta_list[m-1,]
    theta_list[m,] <- EffNUTOne(current_theta, epsilon, LL_log, grad_LL_log, sigma_max)
  }
  
  theta_list[]
}
# Algorithm 4: Heuristic for choosing an initial value of epsilon

FindReasonableEpsilon <- function(theta, LL_log, grad_LL_log){
  epsilon = 1 #default 1
  r = rnorm(length(theta))
  
  leap_list <- Leapfrog(theta, r, epsilon, grad_LL_log =  grad_LL_log)
  theta_p <- leap_list$theta_tilde
  r_p <- leap_list$r_tilde
  
  log_current <- log_joint_density(theta, r, LL_log)
  log_prop <- log_joint_density(theta_p, r_p, LL_log)
  prob <- exp(log_prop - log_current)
  if (prob <= 0| is.nan(prob) | is.na(prob)){prob <- runif(1,0,0.0001)}
  
  a = 2 * ifelse(prob > 0.5, 1, 0) - 1
  j = 0
  
  while ((prob^a > 2^(-a)) & (epsilon > 0)){
    epsilon <- 2^a * epsilon
    
    lista <- Leapfrog(theta, r, epsilon, grad_LL_log)
    theta_p <- lista$theta_tilde
    r_p <- lista$r_tilde
   
    log_current <- log_joint_density(theta, r, LL_log)
    log_prop <- log_joint_density(theta_p, r_p, LL_log)
    prob <- exp(log_prop - log_current)
    if (prob <= 0| is.nan(prob) | is.na(prob)){prob <- runif(1,0,0.0001)}
    
    j= j+1
  }
  
  return(epsilon)
}

# Algorithm 5: Hamilton Monte Carlo with Dual Averaging

HMCDualAvgOne <- function(m, current_theta, current_epsilon, current_epsilon_bar, current_H_bar, 
                          delta, lambda, LL_log, M_adapt, grad_LL_log, epsilon_bar_list,
                          mu, gamma, t_0, K){
  
  r_0 = rnorm(length(current_theta),0,1)
  theta = current_theta
  theta_t = current_theta
  r_t =r_0
  
  if (current_epsilon <= 0.2) {L_m = max(1, lambda/runif(1,0.2,0.5))} else{
    L_m = max(1, round(lambda/current_epsilon))
  }
  
  for (i in c(1:L_m)){
    list_o <- Leapfrog(theta_t, r_t, current_epsilon, grad_LL_log= grad_LL_log)
    theta_t <- list_o$theta_tilde
    r_t <- list_o$r_tilde
    }
  
  prob = exp(LL_log(theta_t) - sum(r_t^2)/2-LL_log(current_theta) + sum(r_0^2)/2)
  if (prob == 0 | is.nan(prob) | is.na(prob)){prob <- runif(1,0.01,0.2)}
  
  a_cri = min(1,prob)
  
  if (runif(1) < a_cri){theta = theta_t}
  
  if (m <= M_adapt){
    H_bar = (1 - 1/(m+t_0))*current_H_bar + 1/(m+t_0)*(delta - a_cri)
    epsilon = exp(mu -sqrt(m)/ gamma * H_bar)
    epsilon_bar = exp(m^(-K)*log(epsilon)+ (1-m^(-K))*log(current_epsilon_bar))}else{
    epsilon = epsilon_bar_list[M_adapt,]
    H_bar = (1 - 1/(m+t_0))*current_H_bar + 1/(m+t_0)*(delta - a_cri)
    epsilon_bar = exp(m^(-K)*log(epsilon)+ (1-m^(-K))*log(current_epsilon_bar))
    }  
  
  return(list(theta = theta, epsilon = epsilon, 
              epsilon_bar = epsilon_bar, H_bar = H_bar))
  
}

# delta: target mean acceptance prob
# lambda: target simulation length ~ epsilon*L
# M_adapt: number of iterations after which to stop the adaption

HMCDualAvg <- function(theta_0, delta, lambda, LL_log, grad_LL_log, M, M_adapt){
  #theta
  theta_list <- matrix(c(0), nrow = M, ncol = length(theta_0))
  theta_list[1,] <- theta_0
  #epsilon
  epsilon0 <- FindReasonableEpsilon(theta_0, LL_log, grad_LL_log)
  epsilon_list <- matrix(c(0), nrow = M, ncol = length(epsilon0))
  epsilon_list[1,] <- epsilon0
  #epsilon_bar 
  epsilon_bar_0 <- 1
  epsilon_bar_list <- matrix(c(0), nrow = M, ncol = length(epsilon_bar_0))
  epsilon_bar_list[1,] <- epsilon_bar_0
  #H_bar
  H_bar_0 <- 0
  H_bar_list <- matrix(c(0), nrow = M, ncol = length(H_bar_0))
  H_bar_list[1,] <- H_bar_0
  
  mu = log(10*epsilon0)
  gamma = 0.05
  t_0 = 10
  K = 0.75
  
  for (m in c(2:M)){
    #m=2
    current_theta <- theta_list[m-1,]
    current_epsilon <- epsilon_list[m-1,]
    current_epsilon_bar <- epsilon_bar_list[m-1,]
    current_H_bar <- H_bar_list[m-1,]
    
    output <- HMCDualAvgOne(m, current_theta, current_epsilon, current_epsilon_bar, current_H_bar, 
                            delta, lambda, LL_log, M_adapt, grad_LL_log, epsilon_bar_list,
                            mu, gamma, t_0, K)
    
    theta_list[m,] <- output$theta
    epsilon_list[m,] <- output$epsilon
    epsilon_bar_list[m,] <- output$epsilon_bar
    H_bar_list[m,] <- output$H_bar
  }
  
  theta_list[]
}

#Algorithm 6: No-U-Turn Sampler with Dual Averaging

BuildTree3 <- function (theta, r, u, v, j, epsilon, sigma_max = 1000, 
                        LL_log, grad_LL_log, theta0, r0){
  if (j == 0){
  
    list1 <- Leapfrog (theta, r, epsilon = v * epsilon, grad_LL_log) # proposal
    theta_p <- list1$theta_tilde
    r_p <- list1$r_tilde
    
    log_dens_prop = log_joint_density(theta_p, r_p,LL_log)
    n_p = ifelse(log(u) <= log_dens_prop, 1, 0) # slice?
    s_p = ifelse(log(u) < sigma_max + log_dens_prop , 1, 0)# accurate simulation?
    
    log_dens_current = log_joint_density(theta0, r0, LL_log)
    ra_dens = exp(log_dens_prop - log_dens_current)
    
    a0 = min(1, ra_dens)
    
    if (is.nan(a0))stop()
    if (is.nan(n_p) | is.na(n_p)){n_p <- 0}
    if (is.nan(s_p) | is.na(s_p)){s_p <- 0}
    
    return(list(theta_b = theta_p, r_b = r_p, theta_f = theta_p, r_f = r_p, 
                theta_p = theta_p, n_p = n_p, s_p = s_p, a_p = a0, n_a_p = 1))}
  else{
    #j=1
    list2 <- BuildTree3(theta, r, u, v, j-1, epsilon, sigma_max = 1000, 
                        LL_log, grad_LL_log, theta0, r0)
    theta_b = list2$theta_b
    r_b = list2$r_b
    theta_f = list2$theta_f
    r_f = list2$r_f
    theta_p = list2$theta_p
    n_p = list2$n_p
    s_p = list2$s_p
    a_p = list2$a_p
    n_a_p = list2$n_a_p
    
    if (s_p == 1){
      if (v == -1){
        list3 <- BuildTree3(theta_b, r_b, u, v, j-1, epsilon,  sigma_max = 1000, 
                            LL_log, grad_LL_log, theta0, r0)
        theta_b = list3$theta_b
        r_b = list3$r_b}
      else{
        list3 <- BuildTree3(theta_f, r_f, u, v, j-1, epsilon,  sigma_max = 1000, 
                            LL_log, grad_LL_log, theta0, r0)
        theta_f = list3$theta_f
        r_f = list3$r_f
      }
      
      theta_p2 = list3$theta_p
      n_p2 = list3$n_p
      s_p2 = list3$s_p
      a_p2 = list3$a_p
      n_a_p2 = list3$n_a_p
      
      n_s = n_p + n_p2
      n_cri <- n_p2/n_s
      
      if (n_s != 0 && runif(1) < n_cri){theta_p = theta_p2}
      
      a_p = a_p + a_p2
      n_a_p = n_a_p + n_a_p2
      
      ip1<- sum((theta_f - theta_b)* r_b)
      sc1 = ifelse(ip1 >= 0, 1, 0)
      ip2<- sum((theta_f - theta_b)* r_f)
      sc2 = ifelse(ip2 >= 0, 1, 0)
      
      s_p = s_p2 *sc1 * sc2
      n_p = n_p + n_p2
    }
    
   # if (is.nan(n_p) | is.na(n_p)){n_p <- 0}
    #if (is.nan(s_p) | is.na(s_p)){s_p <- 0}
    
    return(list(theta_b = theta_b, r_b = r_b, theta_f = theta_f, r_f = r_f,
                theta_p = theta_p, n_p = n_p, s_p = s_p, a_p = a_p, n_a_p = n_a_p))
  }
}


NUTDualAvgOne <- function (m, current_theta, current_epsilon, current_epsilon_bar, current_H_bar, 
                           delta, LL_log, M_adapt, grad_LL_log, 
                           mu, gamma, t_0, K, epsilon_bar_list){
  
  current_r <-rnorm(length(current_theta),0,1)
  
  current_joint_dens <- exp(log_joint_density(current_theta, current_r,LL_log))
  if (is.nan(current_joint_dens)|is.na(current_joint_dens)|current_joint_dens <= 0){
    u <- runif(1,0, 0.00001)}else{
    u <- runif(1,0,current_joint_dens)}
  
  theta_b = current_theta
  theta_f = current_theta
  r_b = current_r
  r_f = current_r
  theta = current_theta
  epilson_bar = current_epsilon_bar
  H_bar = current_H_bar
  epsilon = current_epsilon
  
  j = 0
  n = 1
  s = 1
  
  while (s == 1){
    v_j = sample(c(-1,1),1)
    
    if (v_j == -1){
      listb <- BuildTree3(theta_b, r_b, u, v = v_j, j, epsilon, sigma_max = 1000, 
                          LL_log, grad_LL_log, current_theta, current_r)
      theta_b<- listb$theta_b
      r_b <- listb$r_b
    }else{
      listb <- BuildTree3(theta_f, r_f, u, v = v_j, j, epsilon, sigma_max = 1000, 
                          LL_log, grad_LL_log, current_theta, current_r)
      theta_f<- listb$theta_f
      r_f <- listb$r_f
    }
    
    theta_p <- listb$theta_p
    n_p <- listb$n_p
    s_p <- listb$s_p
    a_p <- listb$a_p
    n_a_p <- listb$n_a_p
    
   #if (is.na(s_p) | is.nan(s_p)){s_p = 0}
    g = min(1, n_p/n)
    
    if (s_p == 1 && runif(1) < g){theta <- theta_p}

    n = n + n_p
    
    ip1<- sum((theta_f - theta_b)* r_b)
    sc1 = ifelse(ip1 >= 0, 1, 0)
    ip2<- sum((theta_f - theta_b)* r_f)
    sc2 = ifelse(ip2 >= 0, 1, 0)
    s = s_p * sc1 * sc2
   
    if (is.nan(s) | is.na(s)){s <- 0}
    
    j = j + 1
    
  }
  
  if (m <= M_adapt){
    H_bar = (1 - 1/(m+t_0))*current_H_bar + 1/(m+t_0)*(delta - a_p/n_a_p)
    epsilon = exp(mu -sqrt(m)/ gamma * H_bar)
    epsilon_bar = exp(m^(-K)*log(epsilon)+ (1-m^(-K))*log(current_epsilon_bar))
    }else{
      H_bar = (1 - 1/(m+t_0))*current_H_bar + 1/(m+t_0)*(delta - a_p/n_a_p)
      epsilon = epsilon_bar_list[M_adapt,]
      epsilon_bar = exp(m^(-K)*log(epsilon)+ (1-m^(-K))*log(current_epsilon_bar))
    }  
  
  return(list(theta = theta, epsilon = epsilon, epsilon_bar = epsilon_bar,
              H_bar = H_bar))
  
}

NUTDualAvg <- function(theta_0, delta, LL_log, grad_LL_log, M, M_adapt){
  #theta
  theta_list <- matrix(c(0), nrow = M, ncol = length(theta_0))
  theta_list[1,] <- theta_0
  #epsilon
  epsilon0 <- FindReasonableEpsilon(theta_0, LL_log, grad_LL_log)
  epsilon_list <- matrix(c(0), nrow = M, ncol = length(epsilon0))
  epsilon_list[1,] <- epsilon0
  #epsilon_bar 
  epsilon_bar_0 <- 1
  epsilon_bar_list <- matrix(c(0), nrow = M, ncol = length(epsilon_bar_0))
  epsilon_bar_list[1,] <- epsilon_bar_0
  #H_bar
  H_bar_0 <- 0
  H_bar_list <- matrix(c(0), nrow = M, ncol = length(H_bar_0))
  H_bar_list[1,] <- H_bar_0
  
  mu = log(10*epsilon0)
  gamma = 0.05 #default 0.05
  t_0 = 10 # >0
  K = 0.75 #(0.5,1]
  
  for (m in c(2:M)){
    #m=2
    current_theta <- theta_list[m-1,]
    current_epsilon <- epsilon_list[m-1,]
    current_epsilon_bar <- epsilon_bar_list[m-1,]
    current_H_bar <- H_bar_list[m-1,]
    
    output <- NUTDualAvgOne(m, current_theta, current_epsilon, current_epsilon_bar, current_H_bar, 
                            delta, LL_log, M_adapt, grad_LL_log, 
                            mu, gamma, t_0, K, epsilon_bar_list)
    theta_list[m,] <- output$theta
    epsilon_list[m,] <- output$epsilon
    epsilon_bar_list[m,] <- output$epsilon_bar
    H_bar_list[m,] <- output$H_bar
  }
  
  theta_list[]
}

##########try
#########Algorithm 6-2: No-U-Turn Sampler with Dual Averaging (considering my scale)

log_joint_density2 <- function(theta, r, LL_log, myscale){
  return(LL_log(theta) - sum(r^2 / diag(myscale))/2)
}

Leapfrog2 <- function (theta, r, epsilon, grad_LL_log, myscale){
  r_tilde = r + (epsilon/2) * grad_LL_log(theta,a = 2, b = 0.1)
  theta_tilde = theta + epsilon * r_tilde / diag(myscale)
  r_tilde = r_tilde + (epsilon / 2 ) * grad_LL_log(theta_tilde, a = 2, b = 0.1)
  return(list(theta_tilde = theta_tilde, r_tilde = r_tilde))
}

FindReasonableEpsilon2 <- function(theta, LL_log, grad_LL_log, myscale){
  epsilon = 1 #default 1
  r = rnorm(length(theta))
  
  leap_list <- Leapfrog2(theta, r, epsilon, grad_LL_log =  grad_LL_log, myscale)
  theta_p <- leap_list$theta_tilde
  r_p <- leap_list$r_tilde
  
  log_current <- log_joint_density2(theta, r, LL_log,myscale)
  log_prop <- log_joint_density2(theta_p, r_p, LL_log, myscale)
  prob <- exp(log_prop - log_current)
  if (prob <= 0| is.nan(prob) | is.na(prob)){prob <- runif(1,0,0.0001)}
  
  a = 2 * ifelse(prob > 0.5, 1, 0) - 1
  j = 0
  
  while ((prob^a > 2^(-a)) & (epsilon > 0)){
    epsilon <- 2^a * epsilon
    
    lista <- Leapfrog2(theta, r, epsilon, grad_LL_log, myscale)
    theta_p <- lista$theta_tilde
    r_p <- lista$r_tilde
    
    log_current <- log_joint_density2(theta, r, LL_log, myscale)
    log_prop <- log_joint_density2(theta_p, r_p, LL_log, myscale)
    prob <- exp(log_prop - log_current)
    if (prob <= 0| is.nan(prob) | is.na(prob)){prob <- runif(1,0,0.0001)}
    
    j= j+1
  }
  
  return(epsilon)
}

BuildTree4 <- function (theta, r, u, v, j, epsilon, sigma_max = 1000, 
                        LL_log, grad_LL_log, theta0, r0, myscale){
  if (j == 0){
    
    list1 <- Leapfrog2 (theta, r, epsilon = v * epsilon, grad_LL_log, myscale) # proposal
    theta_p <- list1$theta_tilde
    r_p <- list1$r_tilde
    
    log_dens_prop = log_joint_density2(theta_p, r_p, LL_log, myscale)
    n_p = ifelse(log(u) <= log_dens_prop, 1, 0) # slice?
    s_p = ifelse(log(u) < sigma_max + log_dens_prop , 1, 0)# accurate simulation?
    
    log_dens_current = log_joint_density2(theta0, r0, LL_log, myscale)
    ra_dens = exp(log_dens_prop - log_dens_current)
    
    a0 = min(1, ra_dens)
    
    if (is.nan(a0))stop()
    if (is.nan(n_p) | is.na(n_p)){n_p <- 0}
    if (is.nan(s_p) | is.na(s_p)){s_p <- 0}
    
    return(list(theta_b = theta_p, r_b = r_p, theta_f = theta_p, r_f = r_p, 
                theta_p = theta_p, n_p = n_p, s_p = s_p, a_p = a0, n_a_p = 1))}
  else{
    #j=1
    list2 <- BuildTree4(theta, r, u, v, j-1, epsilon, sigma_max = 1000, 
                        LL_log, grad_LL_log, theta0, r0, myscale)
    theta_b = list2$theta_b
    r_b = list2$r_b
    theta_f = list2$theta_f
    r_f = list2$r_f
    theta_p = list2$theta_p
    n_p = list2$n_p
    s_p = list2$s_p
    a_p = list2$a_p
    n_a_p = list2$n_a_p
    
    if (s_p == 1){
      if (v == -1){
        list3 <- BuildTree4(theta_b, r_b, u, v, j-1, epsilon,  sigma_max = 1000, 
                            LL_log, grad_LL_log, theta0, r0, myscale)
        theta_b = list3$theta_b
        r_b = list3$r_b}
      else{
        list3 <- BuildTree4(theta_f, r_f, u, v, j-1, epsilon,  sigma_max = 1000, 
                            LL_log, grad_LL_log, theta0, r0, myscale)
        theta_f = list3$theta_f
        r_f = list3$r_f
      }
      
      theta_p2 = list3$theta_p
      n_p2 = list3$n_p
      s_p2 = list3$s_p
      a_p2 = list3$a_p
      n_a_p2 = list3$n_a_p
      
      n_s = n_p + n_p2
      n_cri <- n_p2/n_s
      
      if (n_s != 0 && runif(1) < n_cri){theta_p = theta_p2}
      
      a_p = a_p + a_p2
      n_a_p = n_a_p + n_a_p2
      
      ip1<- sum((theta_f - theta_b)* r_b)
      sc1 = ifelse(ip1 >= 0, 1, 0)
      ip2<- sum((theta_f - theta_b)* r_f)
      sc2 = ifelse(ip2 >= 0, 1, 0)
      
      s_p = s_p2 *sc1 * sc2
      n_p = n_p + n_p2
    }
    
    # if (is.nan(n_p) | is.na(n_p)){n_p <- 0}
    #if (is.nan(s_p) | is.na(s_p)){s_p <- 0}
    
    return(list(theta_b = theta_b, r_b = r_b, theta_f = theta_f, r_f = r_f,
                theta_p = theta_p, n_p = n_p, s_p = s_p, a_p = a_p, n_a_p = n_a_p))
  }
}


NUTDualAvgOne2 <- function (m, current_theta, current_epsilon, current_epsilon_bar, current_H_bar, 
                           delta, LL_log, M_adapt, grad_LL_log, 
                           mu, gamma, t_0, K, epsilon_bar_list, myscale){
  
  current_r <-rnorm(length(current_theta),0, sqrt(diag(myscale)))
  
  current_joint_dens <- exp(log_joint_density2(current_theta, current_r,LL_log, myscale))
  if (is.nan(current_joint_dens)|is.na(current_joint_dens)|current_joint_dens <= 0){
    u <- runif(1,0, 0.00001)}else{
      u <- runif(1,0,current_joint_dens)}
  
  theta_b = current_theta
  theta_f = current_theta
  r_b = current_r
  r_f = current_r
  theta = current_theta
  epilson_bar = current_epsilon_bar
  H_bar = current_H_bar
  epsilon = current_epsilon
  
  j = 0
  n = 1
  s = 1
  
  while (s == 1){
    v_j = sample(c(-1,1),1)
    
    if (v_j == -1){
      listb <- BuildTree4(theta_b, r_b, u, v = v_j, j, epsilon, sigma_max = 1000, 
                          LL_log, grad_LL_log, current_theta, current_r,  myscale)
      theta_b<- listb$theta_b
      r_b <- listb$r_b
    }else{
      listb <- BuildTree4(theta_f, r_f, u, v = v_j, j, epsilon, sigma_max = 1000, 
                          LL_log, grad_LL_log, current_theta, current_r,  myscale)
      theta_f<- listb$theta_f
      r_f <- listb$r_f
    }
    
    theta_p <- listb$theta_p
    n_p <- listb$n_p
    s_p <- listb$s_p
    a_p <- listb$a_p
    n_a_p <- listb$n_a_p
    
    #if (is.na(s_p) | is.nan(s_p)){s_p = 0}
    g = min(1, n_p/n)
    
    if (s_p == 1 && runif(1) < g){theta <- theta_p}
    
    n = n + n_p
    
    ip1<- sum((theta_f - theta_b)* r_b)
    sc1 = ifelse(ip1 >= 0, 1, 0)
    ip2<- sum((theta_f - theta_b)* r_f)
    sc2 = ifelse(ip2 >= 0, 1, 0)
    s = s_p * sc1 * sc2
    
    if (is.nan(s) | is.na(s)){s <- 0}
    
    j = j + 1
    
  }
  
  if (m <= M_adapt){
    H_bar = (1 - 1/(m+t_0))*current_H_bar + 1/(m+t_0)*(delta - a_p/n_a_p)
    epsilon = exp(mu -sqrt(m)/ gamma * H_bar)
    epsilon_bar = exp(m^(-K)*log(epsilon)+ (1-m^(-K))*log(current_epsilon_bar))
  }else{
    H_bar = (1 - 1/(m+t_0))*current_H_bar + 1/(m+t_0)*(delta - a_p/n_a_p)
    epsilon = epsilon_bar_list[M_adapt,]
    epsilon_bar = exp(m^(-K)*log(epsilon)+ (1-m^(-K))*log(current_epsilon_bar))
  }  
  
  return(list(theta = theta, epsilon = epsilon, epsilon_bar = epsilon_bar,
              H_bar = H_bar))
  
}

NUTDualAvg2 <- function(theta_0, delta, LL_log, grad_LL_log, M, M_adapt, myscale){
  #theta
  theta_list <- matrix(c(0), nrow = M, ncol = length(theta_0))
  theta_list[1,] <- theta_0
  #epsilon
  epsilon0 <- FindReasonableEpsilon2(theta_0, LL_log, grad_LL_log, myscale)
  epsilon_list <- matrix(c(0), nrow = M, ncol = length(epsilon0))
  epsilon_list[1,] <- epsilon0
  #epsilon_bar 
  epsilon_bar_0 <- 1
  epsilon_bar_list <- matrix(c(0), nrow = M, ncol = length(epsilon_bar_0))
  epsilon_bar_list[1,] <- epsilon_bar_0
  #H_bar
  H_bar_0 <- 0
  H_bar_list <- matrix(c(0), nrow = M, ncol = length(H_bar_0))
  H_bar_list[1,] <- H_bar_0
  
  mu = log(10*epsilon0)
  gamma = 0.05 #default 0.05
  t_0 = 10 # >0
  K = 0.75 #(0.5,1]
  
  for (m in c(2:M)){
    #m=2
    current_theta <- theta_list[m-1,]
    current_epsilon <- epsilon_list[m-1,]
    current_epsilon_bar <- epsilon_bar_list[m-1,]
    current_H_bar <- H_bar_list[m-1,]
    
    output <- NUTDualAvgOne2(m, current_theta, current_epsilon, current_epsilon_bar, current_H_bar, 
                            delta, LL_log, M_adapt, grad_LL_log, 
                            mu, gamma, t_0, K, epsilon_bar_list, myscale)
    theta_list[m,] <- output$theta
    epsilon_list[m,] <- output$epsilon
    epsilon_bar_list[m,] <- output$epsilon_bar
    H_bar_list[m,] <- output$H_bar
  }
  
  theta_list[]
}