#' Dual Likelihood-ratio based Confidence Interval for Causal Effects 
#'
#' This function calculates an asymptotic Confidence Interval for the Total Causal Effect based on the Dual LRT.
#' It assumes a Linear Structural Causal Model (d nodes) and additive Gaussian errors with partial equal variance
#' between the potential cause and effect.
#' It targets the Total Causal Effect of X (first column) on Y (second column) based on observational (sample size n).
#'
#' @param data n x d matrix, first column X values (potential cause), second column Y values (potential effect)
#' @param alpha 1-confidence level, defaults to 0.05
#' @return A vector containing the bounds of the confidence region. First entry indicates if no effect is plausible.
#' @export
confDualLRTpartial <- function(data, alpha = 0.05) {
  data <- data.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  
  # Compute inverse covariance matrix
  siginv <- solve((n - 1) / n * stats::cov(data, use = 'complete.obs'))
  
  # Approximate MLE by ordering conditional variances
  varperm <- which.min(diag(siginv))
  val <- diag(siginv)[varperm]
  for (i in 1:(d-1)){
    help <- rep(Inf,d)
    for (j in setdiff(1:d,varperm)){
      help[j] <- siginv[j,j]-siginv[j,varperm]%*%solve(siginv[varperm,varperm])%*%siginv[varperm,j]
    }
    varperm <- c(which.min(help),varperm)
    val <- c(help[varperm[1]],val)
  }
  
  
  # Approximate threshold for the likelihood ratio test
  tresh <- -log((val[which(varperm == 1)] + val[which(varperm == 2)]) * sqrt(prod(val[-which(varperm == 1 | varperm == 2)])))
  
  # Find plausible causal orderings (parallelized)
  ordertrue <- parallel::mclapply(0:(3^(d - 2) - 1), orderDualLRTpartial, tresh = tresh, d = d, n = n, siginv = siginv, ord = c(1, 2), crit = stats::qchisq(1 - alpha, 2), mc.cores = parallel::detectCores())
  ordertrue <- do.call(rbind, ordertrue)
  
  orderfalse <- parallel::mclapply(0:(3^(d - 2) - 1), orderDualLRTpartial, tresh = tresh, d = d, n = n, siginv = siginv, ord = c(2, 1), crit = stats::qchisq(1 - alpha, 1), mc.cores = parallel::detectCores())
  orderfalse <- do.call(rbind, orderfalse)
  
  # Find MLE and test for zero effect
  L1 <- max(c(ordertrue[, d + 1], orderfalse[, d + 1]))
  zeropossible <- if (!is.null(orderfalse)) (L1 - max(orderfalse[, d + 1]) <= stats::qchisq(1 - alpha, 1) / (2 * n)) else FALSE
  
  # Reduce set of plausible orderings using true MLE
  if (!is.null(ordertrue)) {
    ind <- (L1 - ordertrue[, d + 1]) <= (stats::qchisq(1 - alpha, 2) / (2 * n))
    ordertrue <- ordertrue[ind, , drop = FALSE]
  }
  
  # If no plausible orderings, return NA otherwise test within intervals for effect
  if (nrow(ordertrue) == 0) {
    return(c(zeropossible, NA, NA))
  } else {
    intervals <- intervalsDualLRTpartial(siginv, ordertrue[, 1:d, drop = FALSE], L1, alpha, d, n)
    return(c(zeropossible, min(intervals), max(intervals)))
  }
}
