.beta_k <-
function(priorHyperParam, cellCounts){
    n <- sum(cellCounts)
    p <- range(as.numeric(unlist(dimnames(cellCounts))))
    p <- (p[2] - p[1] + 1)^2
    
    if(priorHyperParam == "Jeffreys") prior <- 1/2 else
      if(priorHyperParam == "BLUnif") prior <- 1 else
        if(priorHyperParam == "Perks") prior <- 1/p else
          if(priorHyperParam == "MiniMax") prior <- sqrt(n)/p else
            stop("Unknown Prior")
    
    ans <- list(prior = prior, n = n, p = p)
    return(ans)
  }
.diagId <-
function(n){
    if(n < 1) stop("n must be positive")
    if(n == 1) return(1)
    ans <- c(1, rep(0, (n-1)))
    for(i in 2:n) ans[i] <- ans[i-1] + i
    return(ans)
  }
.getIndex <-
function(idx){
    rrow <- ceiling((-1 + sqrt(8*idx + 1))/2)
    ccol <- idx - rrow*(rrow - 1)/2
    ans <- c(rrow, ccol)
    return(ans)
  }
.getTask <-
function(n){
    tasks <- matrix(0, nrow = (n*n + n)/2, ncol = 2)
    iN <- seq_len(n)
    tasks[, 1] <- rep(iN, sort(iN, decreasing = T))
    tasks[, 2] <- unlist(lapply(iN, function(x) iN[x:n]))
    colnames(tasks) <- paste("Var", 1:2)
    return(tasks)
  }
.shrinkageIntensity <-
function(hatThetaML = hatThetaML, n = n,
           shrinkageTarget = shrinkageTarget){
    ans <- 1 - sum(hatThetaML^2)
    ans <- ans/((n - 1)*sum((shrinkageTarget - hatThetaML)^2))
    
    if(ans > 1) ans <- 1
    if(ans < 0) ans <- 0
    #check when n = c(0, 1)
    return(ans)
  }
.startId <-
function(n){
    if(n < 1) stop("n must be positive")
    if(n == 1) return(1)
    ans <- c(1, rep(0, (n-1)))
    for(i in 2:n) ans[i] <- ans[i-1] + (i-1)
    return(ans)
  }
.thetaBayes <-
function(cellCounts, priorHyperParam = priorHyperParam){
    npprior <- .beta_k(priorHyperParam = priorHyperParam,
                      cellCounts = cellCounts)
    B <- npprior$p * npprior$prior
    
    thetak <- (cellCounts + npprior$prior)/(npprior$n + B)
    theta0 <- npprior$prior/(npprior$n + B)
    
    ans <- list(thetak = thetak, theta0 = theta0, p = npprior$p)
    return(ans)
  }
.thetaGT <-
function(cellCounts){
    ans <- .thetaML(cellCounts)
    
    n <- sum(cellCounts)
    m1 <- sum(cellCounts == 1)
    if(m1 == n) m1 <- m1 - 1 #avoid (1 - m1/n) = 0
    
    ans <- (1 - m1/n)*ans
    return(ans)
  }
.thetaML <-
function(cellCounts) cellCounts/sum(cellCounts)
.thetaShrink <-
function(cellCounts, shrinkageTarget = shrinkageTarget){
    hatThetaML <- .thetaML(cellCounts)
    n <- sum(cellCounts)
    
    if(is.null(shrinkageTarget)){
      cellCounts <- as.matrix(cellCounts)
      p <- nrow(cellCounts)*ncol(cellCounts)
      shrinkageTarget <- 1/p
    }
    
    lambda <- .shrinkageIntensity(hatThetaML = hatThetaML, n = n,
                                 shrinkageTarget = shrinkageTarget)
    
    ans <- lambda*shrinkageTarget + (1 - lambda)*hatThetaML
    return(ans)
  }
