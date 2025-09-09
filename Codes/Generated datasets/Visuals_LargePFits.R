### Visualising the results from LargePcompToN

# Load the data
data.p100  <- readRDS("Large_P_fit_results/results_p100.rds")
data.p200  <- readRDS("Large_P_fit_results/results_p200.rds")
data.p400  <- readRDS("Large_P_fit_results/results_p400.rds")
data.p800  <- readRDS("Large_P_fit_results/results_p800.rds")
data.p1200 <- readRDS("Large_P_fit_results/results_p1200.rds")

max_iter <- 100

## Helper function to plot LAD-CD convergence with QuantReg baseline
plot_with_baseline <- function(cd_loss, rq_val, col_cd, col_rq, 
                               main, ylab="MAE") {
  plot(1:length(cd_loss), cd_loss, type="l", col=col_cd, lwd=2,
       xlab="Iteration", ylab=ylab, main=main,
       ylim = c(0, max(cd_loss, na.rm=TRUE)*1.1))
  if (!is.na(rq_val)) {
    abline(h=rq_val, col=col_rq, lwd=2, lty=2)
    legend("topright", legend=c("LAD-CD", "QuantReg"), 
           col=c(col_cd, col_rq), lty=c(1,2), lwd=2)
  } else {
    legend("bottomright", legend=c("LAD-CD", 
                                "No result for RegQuant"), 
           col=c(col_cd, col_rq), lty=1, lwd=2)
  }
}

## p = 100, 200


pdf("Plots/LargeP100and200.pdf")
par(mfrow=c(2,2))
plot_with_baseline(data.p100$res_cd$MAE.loss, data.p100$res_rq$MAE, 
                   "blue", "darkgreen", "Prediction MAE Loss (p=100)")
plot_with_baseline(data.p200$res_cd$MAE.loss, data.p200$res_rq$MAE, 
                   "purple2", "darkgreen", "Prediction MAE Loss (p=200)")

plot_with_baseline(data.p100$res_cd$MAE.params, data.p100$res_rq$MAE.params, 
                   "red", "darkgreen", "Parameter MAE Loss (p=100)")
plot_with_baseline(data.p200$res_cd$MAE.params, data.p200$res_rq$MAE.params, 
                   "gold", "darkgreen", "Parameter MAE Loss (p=200)")
dev.off()

## p = 400, 800


pdf("Plots/LargeP400and800.pdf")
par(mfrow=c(2,2))
plot_with_baseline(data.p400$res_cd$MAE.loss, data.p400$res_rq$MAE, 
                   "blue", "darkgreen", "Prediction MAE Loss (p=400)")
plot_with_baseline(data.p800$res_cd$MAE.loss, data.p800$res_rq$MAE, 
                   "purple2", "darkgreen", "Prediction MAE Loss (p=800)")

plot_with_baseline(data.p400$res_cd$MAE.params, data.p400$res_rq$MAE.params, 
                   "red", "darkgreen", "Parameter MAE Loss (p=400)")
plot_with_baseline(data.p800$res_cd$MAE.params, data.p800$res_rq$MAE.params, 
                   "gold", "darkgreen", "Parameter MAE Loss (p=800)")
dev.off()

## p = 1200


pdf("Plots/LargeP1200.pdf")
par(mfrow=c(1,2))
plot_with_baseline(data.p1200$res_cd$MAE.loss, data.p1200$res_rq$MAE, 
                   "blue", "darkgreen", "Prediction MAE Loss (p=1200)")
plot_with_baseline(data.p1200$res_cd$MAE.params, data.p1200$res_rq$MAE.params, 
                   "red", "darkgreen", "Parameter MAE Loss (p=1200)")
dev.off()
