## ===============================
## LAD vs QuantReg on AirQuality
## ===============================

library(datasets)
library(quantreg)
library(matrixStats)

# ---------- LAD Coordinate Descent Implementation ----------
LAD_CD <- function(X, y, max_iter = 500, tol = 1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  
  beta <- rep(0, p)
  
  wmedian <- function(z, w) {
    matrixStats::weightedMedian(z, w = w)
  }
  
  loss_hist <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    beta_old <- beta
    
    # Intercept update
    r0 <- y - X[, -1, drop = FALSE] %*% beta[-1]
    beta[1] <- median(r0)
    
    # Slopes
    for (j in 2:p) {
      xj <- X[, j]
      nz <- which(xj != 0)  # avoid division by zero
      rj <- y - X[, -j, drop = FALSE] %*% beta[-j]
      z <- rj[nz] / xj[nz]
      w <- abs(xj[nz])
      beta[j] <- wmedian(z, w)
    }
    
    # Store loss
    y_pred <- X %*% beta
    loss_hist[iter] <- mean(abs(y - y_pred))
    
    # Convergence check
    if (sum(abs(beta - beta_old)) < tol) {
      loss_hist <- loss_hist[1:iter]
      break
    }
  }
  
  list(beta = beta, loss = loss_hist)
}

# ---------- Air Quality Data ----------
data("airquality")
df <- na.omit(airquality)   # remove missing rows
y <- df$Ozone
X <- as.matrix(cbind(1, df[, -which(names(df) == "Ozone")]))

# ---------- Fit Models ----------
set.seed(123)
t1 <- system.time(fit_cd <- LAD_CD(X, y, max_iter = 500))
beta_cd <- fit_cd$beta
loss_cd <- fit_cd$loss

t2 <- system.time(fit_rq <- rq(Ozone ~ ., data = df, method = "br"))
beta_rq <- coef(fit_rq)

# ---------- Predictions & MAE ----------
y_pred_cd <- X %*% beta_cd
y_pred_rq <- predict(fit_rq, newdata = df)

mae_cd <- mean(abs(y - y_pred_cd))
mae_rq <- mean(abs(y - y_pred_rq))
coef_dist <- sum(abs(beta_cd - beta_rq))

cat("MAE LAD_CD:", mae_cd, "\n")
cat("MAE QuantReg:", mae_rq, "\n")
cat("Coefficient distance:", coef_dist, "\n")
cat("Runtime LAD_CD:", t1["elapsed"], "s\n")
cat("Runtime QuantReg:", t2["elapsed"], "s\n")

# ---------- Plots ----------
# Pred vs True
plot(y, y_pred_cd, col="blue", pch=20, 
     xlab="True Ozone", ylab="Predicted Ozone",
     main="Air Quality: LAD CD vs QuantReg")
points(y, y_pred_rq, col="green3", pch=17)
abline(0,1,col="red", lwd=2)
legend("topleft", legend=c("LAD CD", "QuantReg"), 
       col=c("blue","green3"), pch=c(20,17))

# Loss convergence
plot(loss_cd, type="o", col="blue", pch=20,
     main="MAE Loss over Iterations",
     xlab="Iteration", ylab="MAE Loss",
     ylim=c(min(c(loss_cd, mae_rq)), max(loss_cd)))
abline(h=mae_rq, col="red", lty=2, lwd=2)
legend("topright", legend=c("LAD CD MAE", "QuantReg MAE"), 
       col=c("blue","red"), lty=c(1,2), pch=c(20,NA))
