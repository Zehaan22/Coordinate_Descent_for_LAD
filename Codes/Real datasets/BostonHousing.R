## ===============================
## Coordinate Descent for LAD in R
## ===============================

# Load packages
library(MASS)       # for Boston dataset
library(matrixStats) # for weightedMedian

# ---------- Algorithm Implementation ----------
LAD_CD <- function(X, y, max_iter = 100, tol = 1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Initialize beta (zeros)
  beta <- rep(0, p)
  
  # Helper function: weighted median
  wmedian <- function(z, w) {
    matrixStats::weightedMedian(z, w = w)
  }
  
  mae.loss <- numeric(max_iter)
  for (iter in 1:max_iter) {
    beta_old <- beta
    
    # ---- Update intercept (first column is 1 for intercept)
    r0 <- y - X[ , -1, drop = FALSE] %*% beta[-1]
    beta[1] <- median(r0)
    
    # ---- Update slopes
    for (j in 2:p) {
      rj <- y - X[ , -j, drop = FALSE] %*% beta[-j]  # partial residual
      z <- rj / X[, j]
      w <- abs(X[, j])
      beta[j] <- wmedian(z, w)
    }
    
    ## computing the losses
    r <- y - X %*% beta
    mae.loss[iter] <- mean(abs(r))
    
    # Convergence check
    if (sum(abs(beta - beta_old)) < tol) break
  }
  return(list(
    beta = beta,
    mae.loss = mae.loss
  ))
}

# ---------- Boston Housing Example ----------
data(Boston)
y <- Boston$medv
X <- as.matrix(cbind(1, Boston[ , -which(names(Boston)=="medv")]))  # add intercept

# Fit LAD via coordinate descent
set.seed(123)

max_iter <- 400
beta_hat <- LAD_CD(X, y, max_iter = max_iter)$beta

print("Estimated coefficients:")
print(beta_hat)

# ---------- Compare against quantreg (LP-based LAD) ----------
library(quantreg)
rq_fit <- rq(medv ~ ., data = Boston, method = "br")
print("quantreg coefficients (check consistency):")
print(coef(rq_fit))

## Diff between the 2 methods
sum(abs(beta_hat - coef(rq_fit)))

y_pred_cd <- X %*% beta_hat
y_pred_rq <- predict(rq_fit, newdata = Boston)

mae_cd <- mean(abs(y - y_pred_cd))
mae_rq <- mean(abs(y - y_pred_rq))

c(LAD_CD = mae_cd, QuantReg = mae_rq)

system.time(LAD_CD(X, y, max_iter = max_iter))
system.time(rq(medv ~ ., data = Boston, method = "br"))


# ---------- Visualization ----------
y_hat <- X %*% beta_hat
plot(y, y_hat, main = "Boston Housing: LAD CD Fit vs. QuantReg Fit", 
     xlab = "True MEDV", ylab = "Predicted MEDV", pch = 20, col = "blue")
points(y, y_pred_rq, col = "green", pch = 17, cex = 0.7)
abline(0, 1, col = "red", lwd = 2)
legend("topleft", legend = c("LAD CD Predictions", "QuantReg Predictions"), 
       col = c("blue", "green"), pch = c(20, 17))

## MAE loss over iterations
lad_fit <- LAD_CD(X, y, max_iter = max_iter)
plot(lad_fit$mae.loss, type = "b", pch = 20, col = "blue",
     xlab = "Iteration", ylab = "MAE Loss", main = "MAE Loss over Iterations",
     ylim = c(min(lad_fit$mae.loss, mae_rq) * 0.9, max(lad_fit$mae.loss, mae_rq) * 1.1))
abline(h = mae_rq, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("LAD CD MAE", "QuantReg MAE"), 
       col = c("blue", "red"), pch = c(20, NA), lty = c(1, 2), lwd = c(1, 2))
