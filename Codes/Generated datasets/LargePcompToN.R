set.seed(123)
library(matrixStats)
library(quantreg)

n <- 1e3   # samples
max_iter <- 100

for (p in c(100, 200, 400, 800, 1200)) {
  cat("\n============================\n")
  cat("Running for p =", p, "\n")
  cat("============================\n")
  
  ## Generate data
  X <- matrix(rnorm(n * (p-1)), n, (p-1))
  X <- cbind(1, X)
  beta_true <- c(2, 1:(p-1))
  y <- as.vector(X %*% beta_true + rnorm(n, sd = 1))
  
  ## --------------------------
  ## LAD Coordinate Descent
  ## --------------------------
  a <- numeric(max_iter)
  b <- matrix(0, nrow = max_iter, ncol = p-1)
  a_curr <- 1
  b_curr <- rep(1, p-1)
  
  MAE.loss <- numeric(max_iter)
  MAE.params <- numeric(max_iter)
  
  for (i in 1:max_iter) {
    for (j in 1:(p-1)) {
      x.j <- X[, 1+j]
      X.rest <- X[, -c(1, 1+j)]
      z <- y - a_curr - X.rest %*% b_curr[-j]
      z.bar <- z / x.j
      b_curr[j] <- weightedMedian(z.bar, abs(x.j))
      b[i, ] <- b_curr
    }
    w <- y - X[, -1] %*% b_curr
    a_curr <- median(w)
    a[i] <- a_curr
    
    MAE.loss[i] <- mean(abs(y - a_curr - X[, -1] %*% b_curr))
    MAE.params[i] <- mean(abs(c(a_curr, b_curr) - beta_true))
  }
  
  res_cd <- list(
    beta = c(a_curr, b_curr),
    MAE.loss = MAE.loss,
    MAE.params = MAE.params
  )
  
  ## --------------------------
  ## QuantReg (LP-based)
  ## --------------------------
  df <- as.data.frame(cbind(y, X[, -1]))
  names(df) <- c("y", paste0("x", 1:(p-1)))
  
  rq_fit <- tryCatch(
    rq(y ~ ., data = df, method = "br"),
    error = function(e) NULL
  )
  
  if (!is.null(rq_fit)) {
    beta_rq <- coef(rq_fit)
    yhat_rq <- predict(rq_fit, newdata = df)
    MAE_rq <- mean(abs(y - yhat_rq))
    MAE_params_rq <- mean(abs(beta_rq - beta_true))
  } else {
    beta_rq <- rep(NA, p)
    MAE_rq <- NA
    MAE_params_rq <- NA
  }
  
  res_rq <- list(
    beta = beta_rq,
    MAE = MAE_rq,
    MAE.params = MAE_params_rq
  )
  
  ## --------------------------
  ## Save results
  ## --------------------------
  saveRDS(list(
    res_cd = res_cd,
    res_rq = res_rq
  ), file = paste0("Large_P_fit_results\results_p", p, ".rds"))
  
  cat("LAD-CD Final MAE:", tail(res_cd$MAE.loss, 1), "\n")
  cat("QuantReg MAE    :", res_rq$MAE, "\n")
}
