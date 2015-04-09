performanceIndex <-
function(testNet, gsNet){
  ans <- validate(testNet, gsNet)
  
  nEdges <- ans$tp + ans$fp
  recall <- ans$tp/(ans$tp + ans$fn)
  fpr <- ans$fp/(ans$fp + ans$tn)
  precision <- ans$tp/(ans$tp + ans$fp)
  precision[is.na(precision)] <- 0
  accuracy <- (ans$tp + ans$tn)/(ans$tp + ans$fn + ans$fp + ans$tn)
  Fscore <- 2*(precision*recall)/(precision + recall)
  Fscore[is.na(Fscore)] <- 0
  
  ans <- cbind(ans[, 1], nEdges, ans[, -1], recall, fpr, precision,
               accuracy, Fscore)
  colnames(ans) <- c("Thresh", "nEdges", "TP", "FP", "FN", "TN",
                     "Recall", "FPR", "Precision", "Accuracy", "Fscore")
  
  return(ans)
}
