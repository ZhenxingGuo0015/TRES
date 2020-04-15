####################################################
## Some functions for EM algorithm
####################################################
require(matrixStats)

## EM algorithm for two component normal mixture
mix.2norm <- function(Y, pi=0.5, nmaxiter=100, TOL=1e-10*length(Y),
                      var.eq=TRUE) {
  n <- length(Y)
  sorty <- sort(Y)

  ## get initial values using pi0=0.5
  n0 <- round(n*pi)
  sd1 <- sd2 <- sd(Y)
  mu1 <- mean(sorty[1:n0])
  mu2 <- mean(sorty[(n0+1):n])
  tmp <- weights <- matrix(0, length(Y), 2)
  loglik.old <- 0

  for(i in 1:nmaxiter) {
    ## E-step: calculate responsibilities
    weights[,1] <- dnorm(Y, mu1, sd1)
    weights[,2] <- dnorm(Y, mu2, sd2)
    tmp <- sweep(weights, 2, c(1-pi, pi), "*")
    ll<- rowSums(tmp)
    gamma <- sweep(tmp, 1, ll, "/")
    loglik <- sum(log(ll))
    if(abs(loglik-loglik.old) < TOL)
      break
    loglik.old <- loglik

    s <- sum(gamma[,2])
    ## M-step: update mu's and sd's
    mu1 <- sum(gamma[,1]*Y)/(n-s)
    mu2 <- sum(gamma[,2]*Y)/s
    if(var.eq) {
      sd1 <- sd2 <- sqrt((sum(gamma[,1]*(Y-mu1)^2)+sum(gamma[,2]*(Y-mu2)^2))/n)
    }
    else {
      sd1 <- sqrt(sum(gamma[,1]*(Y-mu1)^2)/(n-s))
      sd2 <- sqrt(sum(gamma[,2]*(Y-mu2)^2)/s)
    }
    pi <- s/n
    cat(i, "loglik=", loglik,": mu1=", mu1, ", mu2=",
        mu2, "sd1=", sd1, "sd2=",sd2,"\n")
  }

  ## return a list
  list(mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2, pi=pi)
}

## EM algorithm for N component normal mixture
mix.Nnorm <- function(Y, N, pi, nmaxiter=100, TOL=1e-10*length(Y),
                      var.eq=TRUE) {
  n <- length(Y)
  sorty <- sort(Y)

  if(missing(pi))
      pi <- rep(1/N, N)
  lpi <- log(pi)

  ## get initial values
  nn <- c(1, cumsum(round(n*pi)))
  sds <- mus <- rep(0, N)
  for(i in 1:N) {
      ii <- nn[i]:nn[i+1]
      sds[i] <- sd(sorty)
      mus[i] <- mean(sorty[ii])
  }
  tmp <- weights <- matrix(0, length(Y), N)
  loglik <- 0

  mus.new <- mus
  sds.new <- sds
  pi.new <- pi

  ## start EM iteration
  for(i in 1:nmaxiter) {
      ## E-step: calculate responsibilities
      for(j in 1:N)
          weights[,j] <- dnorm(Y, mus[j], sds[j], log=TRUE)

      tmp <- sweep(weights, 2, lpi, "+")
      ll <- apply(tmp, 1, Rsumlog)
      loglik.new <- sum(ll)
      gamma <- exp(sweep(tmp, 1, ll, "-")) ## posterior probability
      s <- colSums(gamma)

      ## M-step: update mu's and sd's
      mus.new <- colSums(gamma*Y) / s
      if(var.eq) {
          tmp <- 0
          for(j in 1:N)
              tmp <- tmp + sum(gamma[,j]*(Y-mus[j])^2)
          sds.new[] <- sqrt(tmp/n)
      }
      else {
          for(j in 1:N)
              sds.new[j] <- sqrt(sum(gamma[,j]*(Y-mus[j])^2)/s[j])
      }
      pi.new <- s/n
      pi.new[pi.new==0] <- 0.001
      pi.new[pi.new==1] <- 1-0.001

      ## check convergence - I'll just use likelihood
      if(abs(loglik.new-loglik) < TOL)
          break
      loglik <- loglik.new
      mus <- mus.new
      sds <- sds.new
      pi <- pi.new
      lpi <- log(pi)

      cat(i, "loglik=", loglik,": mus=", mus, ", sds=", sds, "\n")
  }

  ## return a list
  list(mus=mus, sds=sds, pi=pi)
}

## EM algorithm for N component normal mixture.
## Data are higher dimensional: each observation is a multivariate normal
## Input Y is a matrix of MxK (M dimension, K observations)
## mus and sds are of dimension MxN: each dimension has its own mean and variance.
## Assume there's no correlation among M dimension.
mix.MVNnorm <- function(Y, N, pi, nmaxiter=100, TOL=0.001) {
    M <- nrow(Y)
    K <- ncol(Y)

    if(missing(pi))
        pi <- rep(1/N, N)
    lpi <- log(pi)

    ## get initial values - this is important. I'll do based on k-means
    result.kmeans = kmeans(t(Y), N)
    mus = sds = matrix(0,nrow=M, ncol = N)
    for(i in 1:max(result.kmeans$cluster)) {
        ii = result.kmeans$cluster == i
        mus[,i] = rowMeans(Y[,ii])
        sds[,i] = sqrt(rowVars(Y[,ii]))
    }
    p = as.numeric(table(result.kmeans$cluster) / length(result.kmeans$cluster))

    tmp <- weights <- matrix(0, nrow=K, ncol=N)
    loglik <- 0
    mus.new <- mus
    sds.new <- sds
    pi.new <- pi

    ## start EM iteration
    for(i in 1:nmaxiter) {
        ## E-step: calculate responsibilities
        for(j in 1:N)
            weights[,j] <- colSums(dnorm(Y, mus[,j], sds[,j], log=TRUE))

        tmp <- sweep(weights, 2, lpi, "+")
        ll <- apply(tmp, 1, Rsumlog)
        loglik.new <- sum(ll)
        gamma <- exp(sweep(tmp, 1, ll, "-")) ## posterior probability
        s <- colSums(gamma)

        if(any(is.nan(gamma)))
            browser()

        ## M-step: update mu's and sd's
        tmp = Y %*% gamma
        mus.new = sweep(tmp, 2, colSums(gamma), FUN="/")

        for(j in 1:N) {
            rss = (Y - mus[,j])^2
            sds.new[,j] = colSums(t(rss) * gamma[,j]) / sum(gamma[,j])
        }
        sds.new[sds.new<1e-5] = 1e-5 ## can't have 0 SD

        pi.new <- s/K
        pi.new[pi.new==0] <- 0.001
        pi.new[pi.new==1] <- 1-0.001

        ## check convergence - I'll just use likelihood
        if(abs(loglik.new-loglik) < TOL)
            break
        loglik <- loglik.new
        mus <- mus.new
        sds <- sds.new
        pi <- pi.new
        lpi <- log(pi)

        cat(i, "loglik=", loglik, "\n")
    }

    ## return a list
    list(mus=mus, sds=sds, pi=pi, gamma=gamma)
}

