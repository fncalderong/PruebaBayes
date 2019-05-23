##### Bayesian Hidden Markov Model for count time series #####

##### A.1 The functions #####

# A.1.1 Transforming natural parameters to working
pois.HMM.pn2pw <- function(m, lambda, gamma, delta = NULL, stationary = TRUE){
  tlambda <- log(lambda)
  if(m == 1) return(tlambda)
  foo <- log(gamma / diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  if(stationary) {tdelta <- NULL}
  else{tdelta <- log(delta[-1]/delta[1])}
  parvect <- c(tlambda, tgamma, tdelta)
  return(parvect)
}

# A.1.2 Transforming working parameters to natural
pois.HMM.pw2pn <- function(m, parvect, stationary = TRUE ){
  lambda <- exp(parvect[1:m])
  gamma <- diag(m)
  if(m == 1) return(list(lambda = lambda , gamma = gamma , delta = 1))
  gamma[!gamma] <- exp(parvect[(m + 1):(m * m)])
  gamma <- gamma / apply(gamma ,1, sum)
  if(stationary){delta <- solve(t(diag(m) - gamma +1), rep(1,m))}
  else{foo <-c(1, exp(parvect[(m*m +1) :(m * m + m - 1)]))
  delta <- foo / sum(foo)}
  return(list(lambda = lambda , gamma = gamma , delta = delta))
}

# A.1.3 Computing minus the log-likelihood from the working parameters
pois.HMM.mllk <- function(parvect, x, m, stationary = TRUE,...){
  if(m == 1) return(-sum(dpois(x, exp(parvect), log = TRUE)))
  n <- length(x)
  pn <- pois.HMM.pw2pn(m, parvect, stationary = stationary )
  foo <- pn$delta * dpois(x[1] , pn$lambda)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  for(i in 2:n){
    if(!is.na(x[i])){P <- dpois(x[i], pn$lambda)}
    else{P <- rep(1,m)}
    foo <- foo %*% pn$gamma * P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  mllk <- -lscale
  return(mllk)
}

# A.1.4 Computing the MLEs, given starting values for the natural parameters
pois.HMM.mle <- function (x, m, lambda0, gamma0, delta0 = NULL, stationary = TRUE,...){
  parvect0 <- pois.HMM.pn2pw(m, lambda0, gamma0, delta0, stationary = stationary)
  mod <- nlm(pois.HMM.mllk, parvect0, x = x, m = m, stationary = stationary )
  pn <- pois.HMM.pw2pn(m = m, mod$estimate, stationary = stationary)
  mllk <- -mod$minimum
  np <- length(parvect0)
  AIC <- -2*mllk + 2*np
  n <- sum(!is.na(x))
  BIC <- -2*mllk + np*log(n)
  list(m = m, lambda = pn$lambda, gamma = pn$gamma, delta = pn$delta, code = mod$code,
       mllk = mllk, AIC = AIC, BIC = BIC)
}

# A.1.5 Generating a sample
pois.HMM.generate_sample <- function(ns , mod){
  mvect <- 1:mod$m
  state <- numeric(ns)
  state [1] <- sample(mvect, 1, prob = mod$delta)
  for(i in 2: ns) state[i] <- sample(mvect, 1, prob = mod$gamma[state[i - 1], ])
  x <- rpois(ns, lambda = mod$lambda[state])
  return(x)
}

# A.1.6 Global decoding by the Viterbi algorithm
pois.HMM.viterbi <- function(x, mod){
  n <- length(x)
  xi <- matrix(0, n, mod$m)
  foo <- mod$delta * dpois(x[1], mod$lambda)
  xi[1 ,] <- foo / sum(foo)
  for(i in 2:n){
    foo <- apply(xi[i - 1, ] * mod$gamma, 2, max) * dpois(x[i], mod$lambda)
    xi[i, ] <- foo / sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  for(i in (n - 1):1) iv[i] <- which.max(mod$gamma[ ,iv[i + 1]] * xi[i, ])
  return (iv)
}

# A.1.7 Computing log(forward probabilities)
pois.HMM.lforward <- function(x, mod){
  n <- length(x)
  lalpha <- matrix(NA, n, mod$m)
  foo <- mod$delta * dpois(x[1], mod$lambda)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[1, ] <- lscale + log(foo)
  for(i in 2:n){
    foo <- foo %*% mod$gamma * dpois(x[i], mod$lambda)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[i, ] <- log(foo) + lscale
  }
  return(lalpha)
}

# A.1.8 Computing log(backward probabilities)
pois.HMM.lbackward <- function(x, mod){
  n <- length(x)
  m <- mod$m
  lbeta <- matrix(NA, n, m)
  lbeta[n, ] <- rep(0, m)
  foo <- rep(1/m, m)
  lscale <- log(m)
  for(i in (n - 1):1){
    foo <- mod$gamma %*% (dpois(x[i + 1] , mod$lambda) * foo)
    lbeta[i, ] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo / sumfoo
    lscale <- lscale + log(sumfoo)
  }
  return(lbeta)
}

# A.1.9 Conditional probabilities
# Conditional probability that observation at time t equals xc, given all observations other than that at time t.
# Note : xc is a vector and the result (dxc) is a matrix .
pois.HMM.conditional <- function(xc, x, mod){
  n <- length(x)
  m <- mod$m
  nxc <- length(xc)
  dxc <- matrix(NA, nrow = nxc, ncol = n)
  Px <- matrix (NA, nrow = m, ncol = nxc)
  for(j in 1:nxc) Px[, j] <- dpois(xc[j], mod$lambda)
  la <- t(pois.HMM.lforward(x, mod))
  lb <- t(pois.HMM.lbackward(x, mod))
  la <- cbind(log(mod$delta), la)
  lafact <- apply(la, 2, max)
  lbfact <- apply(lb, 2, max)
  for(i in 1:n){
    foo <- (exp(la[, i] - lafact[i]) %*% mod$gamma) * exp(lb[, i] - lbfact[i])
    foo <- foo / sum(foo)
    dxc[, i] <- foo %*% Px
  }
  return(dxc)
}

# A.1.10 Pseudo-residuals
pois.HMM.pseudo_residuals <- function(x, mod){
  n <- length(x)
  cdists <- pois.HMM.conditional(xc = 0:max(x), x, mod)
  cumdists <- rbind(rep(0, n), apply(cdists, 2,cumsum))
  ulo <- uhi <- rep(NA ,n)
  for(i in 1:n){
    ulo[i] <- cumdists[x[i] + 1 ,i]
    uhi[i] <- cumdists[x[i] + 2 ,i]
  }
  umi <- 0.5 *(ulo + uhi)
  npsr <- qnorm(rbind(ulo, umi, uhi))
  return(npsr)
}

# A.1.11 State probabilities
pois.HMM.state_probs <- function(x, mod){
  n <- length(x)
  la <- pois.HMM.lforward(x, mod)
  lb <- pois.HMM.lbackward(x, mod)
  c <- max(la[n, ])
  llk <- c + log(sum(exp(la[n, ] - c)))
  stateprobs <- matrix(NA, nrow = n, ncol = mod$m)
  for(i in 1:n) stateprobs[i, ] <- exp (la[i, ] + lb[i, ] - llk )
  return(stateprobs)
}

# A.1.12 State prediction
# Note that state output `statepreds ' is a matrix even if h =1.
pois.HMM.state_prediction <- function(h = 1, x, mod){
  n <- length(x)
  la <- t(pois.HMM.lforward(x, mod))
  c <- max(la[, n])
  llk <- c + log(sum(exp(la[, n] - c)))
  statepreds <- matrix(NA, ncol = h, nrow = mod$m)
  foo <- exp(la[, n] - llk)
  for(i in 1:h){
    foo <- foo %*% mod$gamma
    statepreds[, i] <- foo
  }
  return(t(statepreds))
}

# A.1.13 Local decoding
pois.HMM.local_decoding <- function(x, mod){
  n <- length(x)
  stateprobs <- pois.HMM.state_probs(x, mod)
  ild <- rep(NA, n)
  for (i in 1:n) ild[i] <- which.max(stateprobs[i, ])
  ild
}

# A.1.14 Forecast probabilities
# Note that the output `dxf ' is a matrix .
pois.HMM.forecast <- function(xf, h=1, x, mod){
  n <- length(x)
  nxf <- length(xf)
  dxf <- matrix(0, nrow = h, ncol = nxf)
  foo <- mod$delta * dpois (x[1], mod$lambda)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  for(i in 2:n){
    foo <- foo %*% mod$gamma * dpois(x[i], mod$lambda)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  for(i in 1:h){
    foo <- foo %*% mod$gamma
    for (j in 1:mod$m) dxf[i, ] <- dxf[i, ] + foo[j] * dpois(xf, mod$lambda[j])
  }
  return(dxf)
}

# A.1.15 Stationary distribution
stadist <- function(modelo){
  #m <- ncol(gamma)
  #delta <- rep(1,m)%*%solve(diag(m) - gamma + matrix(1,m,m))
  G =  modelo$gamma
  delta = 1/sum(G[!diag(diag(G))])*G[!diag(diag(G))]
  return(delta)
}

# A.1.16 Poisson HMM moments
pois.HMM.moments <- function(modelo, lag.max=10) {
  require(expm)
  lambda <- modelo$lambda
  Gamma <- modelo$gamma
  delta <- modelo$delta
  EX <- delta%*%lambda
  VARX <- EX + delta%*%diag(lambda)%*%lambda - EX^2
  COR <- rep(NA, lag.max)
  for(i in 1:lag.max){
    nume <- (delta%*%diag(lambda)%*%expm::`%^%`(x = Gamma, k = i)%*%lambda - EX^2)
    deno <- (delta%*%diag(lambda)%*%lambda+EX-EX^2)
    COR[i] <-  nume / deno  
  }
  return(list(mu = EX, sigma = VARX, rho = COR))
}

##### PHMM and ZIPHMM #####

# A.1.17 Poisson Hidden Markov Model
bayes.PHMM <- function(y, m = 2, ...){
  standata <- list(N=length(y), m=m, y=y)
  out <- rstan::sampling(stanmodels$PHMM, data = standata, ...)
  return(out)
}

# A.1.18 ZIP-HMM (Zero inflated State == 1)
bayes.ZIPHMM1 <- function(y, m = 2, ...){
  standata <- list(N=length(y), m=m, y=y)
  mod <- ifelse(m==2,stanmodels$ZIPHMM2_INF_1,stanmodels$ZIPHMM_INF_1)
  out <- rstan::sampling(mod, data = standata, 
                     pars = c("theta","lambda","Gamma","lp__"),...)
  return(out)
}

# A.1.19 ZIP-HMM (Zero inflated all State)
bayes.ZIPHMM <- function(y, m = 2, ...){
  standata <- list(N=length(y), m=m, y=y)
  out <- rstan::sampling(stanmodels$ZIPHMM, data = standata, 
                     pars = c("theta","lambda","Gamma","lp__"),...)
  return(out)
}

#####