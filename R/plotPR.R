plotPR <-
function(piNet, ...){
    plot(piNet[, "Recall"], piNet[, "Precision"], type = "l",
         xlab = "recall", ylab = "precision", main = "PR Curve",
         xlim = 0:1, ylim = 0:1, ...)
  }
