LDBayes <- function(x, y, ID, W, K, group, thin, burnin) {
  # Input:
  # x: n*p size matrix of genotype data
  # y: n dimension vector of phenotype data
  # ID: The line ID (note that one line corresponds to one genotype, but may correspond to multiple phenotypes)
  # W: the design matrix for environments (e.g. years or locations or year&location combinations)
  # K: Number of recorded MCMC samples (after doing thinning and droping burnin)
  # thin: for example thin=10, record every 10th of MCMC samples to reduce serial correlation among samples
  # burnin: Discard first number of MCMC samples to gurantee the samples converge to the target posterior discussion
  # group: Indices for marker clusters
  # optim_on: boolean whether to turn on code optimisations

  # Call R package mvtnorm, for generating random numbers from a multivariate normal distribution
  require(mvtnorm)

  # Function to monitor the iterations.
  mod <- function(x, m) {
    t1 <- floor(x / m)
    return(x - t1 * m)
  }



  ########################

  n <- dim(x)[1]

  p <- dim(x)[2]

  G <- unique(group)

  m <- length(G)

  # Precalculation, needed in MCMC interations
  Ssquare <- colSums(x^2)

  # total number of MCMC iterations
  iter <- K * thin + burnin

  # timing results
  timing <- rep(0, iter)

  # Allocate intitial values for all the parameters of interest, and specify hyperparameters
  b0 <- rep(0, K)
  b0u <- 0

  uid <- unique(ID)
  N <- length(uid)

  bw <- matrix(0, nrow = n, ncol = K)
  bwu <- matrix(0, nrow = n, ncol = 1)


  bw2 <- matrix(0, nrow = dim(W)[2], ncol = K)
  bwu2 <- matrix(0, nrow = dim(W)[2], ncol = 1)


  bu <- rep(0, p)
  b <- matrix(0, nrow = p, ncol = K)


  pzu <- 0.5
  pz <- rep(0, K)

  p0 <- 10

  pi <- 0.05



  Ru <- RBu <- rep(0, p)
  R <- RB <- matrix(0, nrow = p, ncol = K)

  ###############################################
  R2 <- 0.5

  MSx <- sum(apply(x, 2, var))

  alpha0 <- 5 / 2

  beta0 <- var(y) * (1 - R2) * (2 * alpha0 + 2) / 2


  alphab <- 5 / 2

  betab <- var(y) * R2 * (2 * alphab + 2) / (2 * MSx * pi)



  alphaw <- 1 / 10

  betaw <- 1 / 10


  ###############################################

  # Residual variance
  tau0u <- beta0 / (1 + alpha0)
  tau0 <- rep(1, K)

  # Marker variance
  tauu <- rep(1, p)
  tau <- matrix(1, nrow = p, ncol = K)

  # Experiment variance
  tauwu <- 1
  tauw <- rep(1, K)


  # Experiment variance
  tauwu2 <- 1
  tauw2 <- rep(1, K)


  # residuals
  r <- y - b0u - x %*% bu - bwu - W %*% bwu2

  count <- 0

  # RD
  WW <- t(W) %*% W
  #

  # Gibbs sampling steps (number of iterations: iter)
  for (i in 1:iter) {
    # print(i)

    b0old <- b0u
    z <- matrix(1, nrow = n, ncol = 1)
    mu <- (t(z) %*% z)^(-1) %*% t(z) %*% (r + z %*% b0old)
    sigma <- sqrt(1 / n * tau0u)
    b0u <- rnorm(1, mean = mu, sd = sigma)
    r <- r + z %*% (b0old - b0u)

    #######################
    bwold2 <- bwu2
    tau0u <- as.numeric(tau0u)
    COV <- solve(WW / tau0u + 1 / tauwu2)
    mu <- 1 / tau0u * COV %*% t(W) %*% (r + W %*% bwold2)
    bwu2 <- t(rmvnorm(1, mean = mu, sigma = COV))
    r <- r + W %*% (bwold2 - bwu2)

    #print("First bit")
    #print(end_time-start_time)

    ###########################
    bwold <- bwu

    start_time <- Sys.time()
    for (k in 1:N) {
      id <- which(ID == uid[k])
      W2 <- diag(length(id))
    
      bwuu <- bwold[id]
      rw <- r[id]
      # not sure why you would t(X) %*% X where X is the identity:
      # COV <- solve(W2 / tau0u + 1 / tauwu * W2)
      COV <- W2/(1/tau0u + 1/tauwu)
      # again no need to t(W2):
      mu <- 1 / tau0u * COV %*% t(W2) %*% (rw + W2 %*% bwuu)
      # is a transpose really necessary here? Should be fast though.
      bwuu <- t(rmvnorm(1, mean = mu, sigma = COV))
      bwu[id] <- bwuu
    }

    end_time <- Sys.time()
    #print("k loop:")
    #print(end_time-start_time)

    r <- r + (bwold - bwu)
    ###########################

    # change var to varvec because var() is a core R function
    # Should probably do the same with beta
    varvec <- tau0u / (Ssquare + tau0u / tauu)

    

sigma <- sqrt(varvec)    #update regression parameters beta    

        for (j in 1:p) {
        buoldj <- bu[j]
        but <- sum(r * x[, j])
        mu <- ((bu[j] * Ssquare[j] + but) / tau0u) * varvec[j]
        u1 <- pzu * sqrt(1 / tauu[j]) * sqrt(varvec[j]) * exp(mu^2 / (2 * varvec[j]))

        if (u1 < 10^100) {
          RBu[j] <- u1 / (1 - pzu + u1)
          Ru[j] <- rbinom(1, 1, prob = RBu[j])
        } else {
          Ru[j] <- 1
          RBu[j] <- 1
        }

        if (Ru[j] == 1) {
          bu[j] <- rnorm(1, mean = mu, sd = sigma[j])
        } else {
          bu[j] <- 0
        }

        r <- r + x[, j] * (buoldj - bu[j])
      }
    
    #Update inconlusion probablity
    pzu <- rbeta(1, sum(Ru) + p0 * pi, p - sum(Ru) + p0 * (1 - pi))
    #Update residual variance
    alpha <- n / 2 + alpha0
    beta <- t(r) %*% r / 2 + beta0
    tau0u <- 1 / rgamma(1, shape = alpha, scale = 1 / beta)
    
    #Update variance of regression parameters
    alpha <- 1 / 2 * tapply(Ru, group, sum) + alphab
     beta <- 1 / 2 * tapply(bu, group, function(x){sum(x^2)} ) + betab
      tauu_grouped <- 1 / rgamma(m, shape = alpha, scale = 1 / beta)
     tauu_grouped <- setNames(tauu_grouped,G)
     tauu <- tauu_grouped[group]
    
 ###########################################

    alpha <- 1 / 2 * n + alphaw
    beta <- 1 / 2 * sum(bwu^2) + betaw
    tauwu <- 1 / rgamma(1, shape = alpha, scale = 1 / beta)

    alpha <- 1 / 2 * dim(W)[2] + alphaw
    beta <- 1 / 2 * sum(bwu2^2) + betaw
    tauwu2 <- 1 / rgamma(1, shape = alpha, scale = 1 / beta)


    if (i > burnin & mod(i, thin) == 0) {
      count <- count + 1

      b0[count] <- b0u
      b[, count] <- bu
      bw[, count] <- bwu
      bw2[, count] <- bwu2
      R[, count] <- Ru
      RB[, count] <- RBu
      pz[count] <- pzu
      tau[, count] <- tauu
      tau0[count] <- tau0u
      tauw[count] <- tauwu
      tauw2[count] <- tauwu2
    }
    end_time <- Sys.time()
    #print("Last bit")
    #print(end_time-start_time)


    if (mod(i, 10) == 0) {
      # Print on the screen some message
      cat(paste0("iteration: ", i, "\n"))
    }
  }

  # print("Overall timing for that loop:")
  # print(mean(timing))


  # Output (MCMC samples):
  # b0: intercept
  # bw: environmental specific intercept
  # bw2: Random effect
  # b: Regression parameters (additive genetic effects)
  # R: indicator variable tells whether a marker (variable) should be included in the model
  # pz: Parameter for the prior of R (inclusion probability)
  # tau: Variance component for regression parameters
  # tau0: Variance component for the residual error
  # tauw: Variance component for the environmental specific intercept
  # tauw2: Variance component for the random effect

  return(
    list(
      b0 = b0,
      bw = bw,
      bw2 = bw2,
      b = b,
      R = R,
      pz = pz,
      tau = tau,
      tau0 = tau0,
      tauw = tauw,
      tauw2 = tauw2
    )
  )
}
