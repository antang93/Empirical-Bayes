
ebpred <- function(y, X, XX, alpha, gam, sig2, log.f, M) {

  n <- nrow(X)
  p <- ncol(X)
  d <- nrow(XX)
  if(!exists("lars")) library(lars)
  o.lars <- lars(X, y, normalize=FALSE, intercept=FALSE, use.Gram=FALSE)
  cv <- cv.lars(X, y, plot.it=FALSE, se=FALSE, normalize=FALSE, intercept=FALSE, use.Gram=FALSE)
  b.lasso <- coef(o.lars, s=cv$index[which.min(cv$cv)], mode="fraction")
  S <- as.numeric(b.lasso != 0)
  if(is.null(sig2)) {
 
    z <- as.numeric(y - X[, S > 0] %*% b.lasso[S > 0])
    sig2 <- sum(z**2) / max(n - sum(S), 1)

  }
  v <- 1 + alpha / gam / sig2
  sq.v <- sqrt(v)
  B <- round(0.2 * M)
  SS <- matrix(0, nrow=B + M, ncol=p)
  ss <- numeric(B + M)
  ii <- integer(B + M)
  YY <- newname <- matrix(0, nrow=M, ncol=d)
  SS[1,] <- S
  ss[1] <- s <- sum(S)
  ii[1] <- i <- iuS <- nuS <- 1
  out <- list()
  out[[1]] <- get.lm.stuff(S, y, X)
  lprior <- log.f(0:n) - lchoose(p, 0:n)
  lpost <- lprior[s] - alpha * out[[1]]$sse / 2 / sig2 - s * log(v) / 2
  for(m in 1:(B + M)) {
    
    S.new <- rprop(S, n)
    s.new <- sum(S.new)
    i.new <- compare.to.rows(as.matrix(SS[iuS,]), S.new)
    if(i.new == 0) o.new <- get.lm.stuff(S.new, y, X) else o.new <- out[[i.new]]
    lpost.new <- lprior[s.new] - alpha * o.new$sse / 2 / sig2 - s.new * log(v) / 2
    if(runif(1) <= exp(lpost.new - lpost)) {

      S <- S.new
      s <- s.new
      lpost <- lpost.new
      if(i.new == 0) {

        nuS <- nuS + 1
        i <- nuS
        iuS <- c(iuS, m)
        out[[nuS]] <- o.new

      } else i <- i.new

    }
    SS[m,] <- S
    ss[m] <- s
    ii[m] <- i
    if(m > B) {

      out.mean <- XX[,S > 0] %*% out[[i]]$b.hat
      newname[m - B,] <- XX[,S > 0] %*% out[[i]]$b.hat
      #YY[m - B,] <- out.mean + sq.v * XX[,S > 0] %*% out[[i]]$U %*% rnorm(s) + sqrt(sig2) * rnorm(d)
#^noise can be taken out if we're only looking at mean
    }

  }
  return(list(b.lasso=b.lasso, Yhat = newname))

}

rprop <- function(S, n) {

 s <- sum(S)
 if(s == n) { S[sample(which(S == 1), 1)] <- 0 }
 else if(s  == 1) { S[sample(which(S == 0), 1)] <- 1 }
 else {

   if(runif(1) <= 0.5) S[sample(which(S == 1), 1)] <- 0
   else S[sample(which(S == 0), 1)] <- 1

 }
 return(S)

}


compare.to.rows <- function(SS, S) {
  if (ncol(SS)==1){
    SS <- t(SS)
  }
  numrows <- nrow(SS)
  newS <- matrix(rep(S,numrows), nrow=numrows, byrow=TRUE)
  subtracted <- abs(newS-SS)
  o <- rowSums(subtracted)
  if(all(o > 0)) return(0) else return(which(o == 0)[1])
}


compare.to.rows0 <- function(SS, S) {
  #since this is vector of 0, 1s, 
  #is it better to convert to the number that's represented by the binary numbers?
  #as.binary(S)
  h <- function(v) sum(abs(v - S))
  o <- apply(SS, 1, h)
  if(all(o > 0)) return(0) else return(which(o == 0)[1])
}


get.lm.stuff <- function(S, y, X) {

  if(!exists("ginv")) library(MASS)
  X.S <- as.matrix(X[, S > 0])
  o <- lm.fit(X.S, y, singular.ok = FALSE)
  sse <- sum(o$residuals**2)
  b.hat <- o$coefficients
  #XtX.S <- t(X.S) %*% X.S
  #V <- ginv(XtX.S)
  #U <- chol(V)
  return(list(sse=sse, b.hat=b.hat))

}


get.mode <- function(x) {

  # x is an integer vector
  xt <- table(x)
  mode <- as.integer(names(sort(-xt)[1]))
  return(mode)

}

lassopred <- function(X, y, X.new){ ##lasso
  o.lars <- lars(X, y, normalize=FALSE, intercept=FALSE, use.Gram=FALSE)
  cv <- cv.lars(X, y, plot.it=FALSE, se=FALSE, normalize=FALSE, intercept=FALSE, use.Gram=FALSE)
  b.lasso <- coef(o.lars, s=cv$index[which.min(cv$cv)], mode="fraction")
  pred <- X.new %*% b.lasso
  return (pred)
}

RRpred <- function(X, y, X.new){ ##ridge regression
  o.rr <- lm.ridge(y ~ X)
  b.rr <- o.rr$coef
  pred <- X.new %*% b.rr
  return (pred)
}

PLSRpred <- function(X, y, X.new, p){ ##partial least squares regression
  o.plsr <- plsr(y ~ X, method = pls.options()$plsralg)
  b.plsr <- o.plsr$coefficients[1:p]
  pred <- X.new %*% b.plsr
  return (pred)
}

BLpred <- function(X, y, X.new){ ##bayesian lasso
  o.bl <- blasso(X, y, normalize=FALSE, icept=FALSE, verb = 0)
  b.bl <- o.bl$beta
  pred <- X.new %*% colMeans(b.bl) 
  return (pred)
}

HSpred <- function(X, y, X.new){ ##horseshoe
  #tau.estimate <- HS.MMLE(y, 1)
  o.hs <- horseshoe(y, X, method.tau = "halfCauchy", method.sigma = "fixed", Sigma2 = 1)
  b.hs <- o.hs$BetaHat
  pred <- X.new %*% b.hs
  return (pred)
}

alasso.pred <- function(X, y, X.new){ ##adaptive lasso
  o.alasso <- adalasso(X, y, use.Gram = FALSE, intercept=FALSE)
  b.alasso <- o.alasso$coefficients.adalasso
  pred <- X.new %*% b.alasso
}

BRpred <- function(X, y, X.new){ ##bayesian ridge
  o.br <- bridge(X, y, normalize=FALSE, icept=FALSE, verb = 0)
  b.br <- o.br$beta
  pred <- X.new %*% colMeans(b.br)
  return (pred)
}

PCRpred <- function(X, y, X.new, p){ ##principal components regression
  o.pcr <- pcr(y ~ X, method = pls.options()$pcralg)
  b.pcr <- o.pcr$coefficients[1:p]
  pred <- X.new %*% b.pcr
  return (pred)
}

# Example...

dcomplex <- function(x, n, p, a, b) -x * (log(b) + a * log(p)) + log(x <= n)

library('MASS')
library('monomvn')
library('pls')
library('horseshoe')
library('parcor')

ebpred.sim <- function(reps=100, n=70, p=100, beta=rep(1, 5), r=0.5, M=5000) {

  sig2 <- 1
  alpha <- 0.999
  gam <- 1 - alpha
  log.f <- function(x) dcomplex(x, n, p, 0.05, 1)
  s0 <- length(beta)
  g <- function(i, j) r**(abs(i - j))
  R <- outer(1:p, 1:p, g)
  e <- eigen(R)
  sqR <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
  time <- mspe <- mspe.l <- time.l <- mspe.rr <- time.rr <- mspe.hs <- time.hs <- 
    mspe.alasso <- time.alasso <- mspe.plsr <- time.plsr <- mspe.bl <- time.bl <- 
    mspe.br <- time.br <- mspe.pcr <- time.pcr <- numeric(reps)
  for(k in 1:reps)  {
    
    X <- matrix(rnorm(n * p), nrow=n, ncol=p) %*% sqR
    X.new <- matrix(rnorm(n * p), nrow=n, ncol=p) %*% sqR
    y <- as.numeric(X[, 1:s0] %*% beta) + sqrt(sig2) * rnorm(n)
    y.new <- as.numeric(X.new[, 1:s0] %*% beta) + sqrt(sig2) * rnorm(n)
    time[k] <- system.time(o <- ebpred(y, X, X.new, alpha, gam, NULL, log.f, M))[3] 
    Y.hat <- colMeans(o$Yhat)
    if(reps == 1) hist(Y.hat, freq=FALSE, breaks=25, col="gray", border="white", main="")
    mspe[k] <- mean((y.new - Y.hat)**2) #mspe for empirical bayes
    mspe.l[k] <- mean((y.new - X.new %*% o$b.lasso)**2) #mspe for regular lasso
    time.l[k] <- system.time(lassopred(X, y, X.new))[3] #time for lasso
    mspe.rr[k] <- mean((y.new - RRpred(X, y, X.new))**2) #mspe for ridge regression
    time.rr[k] <- system.time(RRpred(X, y, X.new))[3] #time for ridge regression
    sink('/dev/null')
    mspe.hs[k] <- mean((y.new - HSpred(X, y, X.new))**2) #mspe for hs
    time.hs[k] <- system.time(HSpred(X, y, X.new))[3] #time for hs
    sink()
    mspe.alasso[k] <- mean((y.new - alasso.pred(X, y, X.new))**2) #mspe for adaptive lasso
    time.alasso[k] <- system.time(alasso.pred(X, y, X.new))[3] #time for adaptive lasso
    mspe.plsr[k] <- mean((y.new - PLSRpred(X, y, X.new, p))**2) #mspe for plsr
    time.plsr[k] <- system.time(PLSRpred(X, y, X.new, p))[3] #time for plsr
    #mspe.bl[k] <- mean((y.new - BLpred(X, y, X.new))**2) #mspe for bl
    #time.bl[k] <- system.time(BLpred(X, y, X.new))[3] #time for bl
    #mspe.br[k] <- mean((y.new - BRpred(X, y, X.new))**2) #mspe for br
    #time.br[k] <- system.time(BRpred(X, y, X.new))[3] #time for br
    mspe.pcr[k] <- mean((y.new - PCRpred(X, y, X.new, p))**2) #mspe for pcr
    time.pcr[k] <- system.time(PCRpred(X, y, X.new, p))[3] #time for pcr
    print (k)
  }
  return(cbind(n, p, mspe.eb=mean(mspe, trim=.1), mspe.lasso=mean(mspe.l, trim=.1), 
               mspe.rr=mean(mspe.rr, trim=.1), mspe.hs=mean(mspe.hs, trim=.1),
               mspe.alasso=mean(mspe.alasso, trim=.1), mspe.plsr=mean(mspe.plsr, trim=.1), 
               mspe.pcr=mean(mspe.pcr, trim=.1), time.eb=mean(time), 
               time.lasso=mean(time.l), time.rr=mean(time.rr), time.hs=mean(time.hs), 
               time.alasso=mean(time.alasso), time.plsr=mean(time.plsr), 
               time.pcr=mean(time.pcr), 
               sd.eb=sd(mspe), sd.lasso=sd(mspe.l), PI.eb=quantile(Y.hat, probs=c(.025, .975)), 
               PI.lasso=quantile(X.new %*% o$b.lasso, probs=c(.025, .975))))
}
