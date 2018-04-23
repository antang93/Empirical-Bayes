
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
  YY <- matrix(0, nrow=M, ncol=d)
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
      YY[m - B,] <- out.mean + sq.v * XX[,S > 0] %*% out[[i]]$U %*% rnorm(s) + sqrt(sig2) * rnorm(d)

    }

  }
  return(list(Yhat=YY, b.lasso=b.lasso))

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

  h <- function(v) sum(abs(v - S))
  o <- apply(SS, 1, h)
  if(all(o > 0)) return(0) else return(which(o == 0)[1])

}


get.lm.stuff <- function(S, y, X) {

  if(!exists("ginv")) library(MASS)
  X.S <- as.matrix(X[, S > 0])
  o <- lm(y ~ X.S - 1)
  sse <- sum(o$residuals**2)
  b.hat <- o$coefficients
  XtX.S <- t(X.S) %*% X.S
  V <- ginv(XtX.S)
  U <- chol(V)
  return(list(sse=sse, b.hat=b.hat, U=U))

}


get.mode <- function(x) {

  # x is an integer vector
  xt <- table(x)
  mode <- as.integer(names(sort(-xt)[1]))
  return(mode)

}



# Example...

dcomplex <- function(x, n, p, a, b) -x * (log(b) + a * log(p)) + log(x <= n)


ebpred.sim <- function(reps=100, n=70, p=100, beta=rep(1, 5), M=5000) {

  sig2 <- 1
  alpha <- 0.999
  gam <- 1 - alpha
  log.f <- function(x) dcomplex(x, n, p, 0.05, 1)
  s0 <- length(beta)
  r <- 0.5
  g <- function(i, j) r**(abs(i - j))
  R <- outer(1:p, 1:p, g)
  e <- eigen(R)
  sqR <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
  time <- mspe <- mspe.l <- numeric(reps)
  for(k in 1:reps)  {
    
    X <- matrix(rnorm(n * p), nrow=n, ncol=p) %*% sqR
    X.new <- matrix(rnorm(n * p), nrow=n, ncol=p) %*% sqR
    y <- as.numeric(X[, 1:s0] %*% beta) + sqrt(sig2) * rnorm(n)
    y.new <- as.numeric(X.new[, 1:s0] %*% beta) + sqrt(sig2) * rnorm(n)
    time[k] <- system.time(o <- ebpred(y, X, X.new, alpha, gam, NULL, log.f, M))[3]
    Y.hat <- apply(o$Yhat, 2, mean)
    if(reps == 1) hist(Y.hat, freq=FALSE, breaks=25, col="gray", border="white", main="")
    mspe[k] <- mean((y.new - Y.hat)**2)
    mspe.l[k] <- mean((y.new - X.new %*% o$b.lasso)**2)

  }
  return(cbind(mspe.lasso=mean(mspe.l), mspe.eb=mean(mspe), time.eb=mean(time)))

}
