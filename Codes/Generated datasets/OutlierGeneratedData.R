set.seed(123)
library(matrixStats)
library(quantreg)

### 1. Generate synthetic regression dataset with outliers
n <- 1e3   # samples
p <- 2     # features

X <- matrix(rnorm(n * (p-1)), n, (p-1))
X.outliers <- matrix(rnorm(n/20 * (p-1), sd = 5), n/20, (p-1))  # outliers
X <- rbind(X, X.outliers)
X <- cbind(1,X)
beta_true <- c(2,5)
n <- 1.2 * n 

rows <- sample(1:nrow(X), n/2, replace = FALSE)
X.normal <- X[rows, ]
X.out <- X[-rows, ]

y.normal <- X.normal %*% beta_true + rnorm(nrow(X.normal), sd = 5)
y.out <- X.out %*% beta_true + rnorm(nrow(X.out), sd = 25)

y <- c(y.normal, y.out)
X <- rbind(X.normal, X.out)
n <- nrow(X)

x <- X[,2]

## Visualize
par(mfrow=c(1,1))
plot(X[,2], y, 
     main = "Generated Data with Outliers",
     xlab = "X", ylab = "y",
     pch = 19, col = ifelse(1:n <= nrow(X.normal), "blue", "red"))
legend("bottomright", legend = c("Normal Data", "Outliers"), 
       col = c("blue", "red"), pch = 19)

### 2. LAD–CD Optimization
max_iter <- 100
a <- numeric(max_iter); b <- numeric(max_iter)
a_curr <- 20; b_curr <- 20
MAE.loss <- numeric(max_iter); MAE.params <- numeric(max_iter)

for(i in 1:max_iter){
  z <- y - a_curr 
  z.bar <- z/x
  b_curr <- weightedMedian(z.bar, abs(x))
  b[i] <- b_curr
  
  w <- y - b_curr * x
  a_curr <- median(w)
  a[i] <- a_curr
  
  MAE.loss[i] <- mean(abs(y - a_curr - b_curr*x))
  MAE.params[i] <- mean(abs(c(a_curr, b_curr) - beta_true))
}

cat("True params:", beta_true, "\n")
cat("LAD-CD params:", c(a_curr, b_curr), "\n")

### 3. OLS baseline
ols_fit <- lm(y ~ x)
ols_coef <- coef(ols_fit)
cat("OLS params:", ols_coef, "\n")

### 4. QuantReg baseline
rq_fit <- rq(y ~ x, method = "br")
rq_coef <- coef(rq_fit)
cat("QuantReg params:", rq_coef, "\n")

### 5. Visualise fits
pdf("Plots/EMLAD_ModelFit_Outliers.pdf", width=8, height=6)
par(mfrow=c(1,1))
plot(x, y, pch=19, col="black",
     main="Fitted Lines under Outliers",
     xlab="X", ylab="y")
abline(a=a_curr, b=b_curr, col="red", lwd=2)         # LAD-CD
abline(ols_coef, col="green3", lwd=2, lty=2)         # OLS
abline(rq_coef, col="purple", lwd=3, lty=4)          # QuantReg
legend("bottomright", legend=c("LAD-CD","OLS","QuantReg"),
       col=c("red","green3","purple"), lwd=2, lty=c(1,2,3))
dev.off()

### 6. Plot MAE curves (with OLS & QuantReg baselines)

# Prediction MAE for OLS & QuantReg
mae_ols <- mean(abs(y - predict(ols_fit)))
mae_rq  <- mean(abs(y - predict(rq_fit)))

# Parameter MAE for OLS & QuantReg
mae_params_ols <- mean(abs(ols_coef - beta_true))
mae_params_rq  <- mean(abs(rq_coef - beta_true))


pdf("Plots/MAE_params_and_loss_Outliers.pdf", width=8, height=10)
par(mfrow=c(2,1))
# Prediction MAE
plot(1:max_iter, MAE.loss, col="red", type="l", lwd=2,
     main="Prediction MAE across iterations",
     xlab="Iteration", ylab="Prediction MAE")
abline(h=mae_ols, col="green3", lwd=2, lty=2)
abline(h=mae_rq, col="purple", lwd=2, lty=3)
legend("topright", legend=c("LAD-CD","OLS","QuantReg"),
       col=c("red","green3","purple"), lty=c(1,2,3), lwd=2)

# Parameter MAE
plot(1:max_iter, MAE.params, col="blue", type="l", lwd=2,
     main="Parameter MAE across iterations",
     xlab="Iteration", ylab="Parameter MAE",
     ylim = c(0, max(MAE.params, mae_params_ols, mae_params_rq)*1.1))
abline(h=mae_params_ols, col="green3", lwd=2, lty=2)
abline(h=mae_params_rq, col="purple", lwd=2, lty=3)
legend("topright", legend=c("LAD-CD","OLS","QuantReg"),
       col=c("blue","green3","purple"), lty=c(1,2,3), lwd=2)
dev.off()

