plotROC <-
function(piNet, ...){
  plot(piNet[, "FPR"], piNet[, "Recall"], type = "l",
       xlab = "FP rate", ylab = "TP rate", main = "ROC Curve",
       xlim = 0:1, ylim = 0:1, ...)
  lines(0:1, 0:1, col = "black")
}
