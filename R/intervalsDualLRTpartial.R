#' Compute the confidence intervals for the causal effect
#'
#' This (helper) function computes the confidence region of the total causal effect for a set of causal orderings based on Dual LRT.
#'
#' @param siginv Inverse covariance matrix of data
#' @param perm Set of causal orderings
#' @param L1 Likelihood under alternative
#' @param alpha 1-confidence level
#' @param d Number of variables
#' @param n Number of observations
#' @return Vector of plausible causal effect intervals
#' @export
intervalsDualLRTpartial <- function(siginv, perm, L1, alpha, d, n) {
  critval <- stats::qchisq(1 - alpha, 2)
  returnintervals <- c()
  intervals <- c()
  
  for (i in seq_len(nrow(perm))) {
    ord <- perm[i, ]
    
    # Compute likelihood for the ordering starting from sink nodes
    val <- numeric(d)
    val[d] <- siginv[ord[d], ord[d]]
    for (j in seq_len(d - 1)) {
      val[d - j] <- siginv[ord[d - j], ord[d - j]] -
        siginv[ord[d - j], ord[(d - j + 1):d]] %*%
        solve(siginv[ord[(d - j + 1):d], ord[(d - j + 1):d]]) %*%
        siginv[ord[(d - j + 1):d], ord[d - j]]
    }
    
    # treshold for dual likelihood ratio test
    tresh <- exp((critval - 2 * n * L1) / (2 * n)) / sqrt(prod(val[-which(ord == 1 | ord == 2)])) - val[which(ord == 2)]
    
    # Find the set of descendants of potential cause
    de <- ord[d:(which(ord == 1) + 1)]
    de <- de[de != 2]
    
    # compute the confidence intervals
    if (length(de) == 0) {
      K <- sqrt(siginv[1, 2]^2 - siginv[2, 2] * (siginv[1, 1] - tresh))
    } else {
      S_de <- solve(siginv[de, de])
      adj_1 <- siginv[2, de] %*% S_de %*% siginv[de, 1]
      adj_2 <- siginv[2, de] %*% S_de %*% siginv[de, 2]
      K <- sqrt((siginv[1, 2] - adj_1)^2 - (siginv[2, 2] - adj_2) * (siginv[1, 1] - siginv[1, de] %*% S_de %*% siginv[1, de] - tresh))
    }
    
    L <- (ifelse(length(de) > 0, adj_1, 0) - siginv[1, 2] - K) / (siginv[2, 2] - ifelse(length(de) > 0, adj_2, 0))
    U <- (ifelse(length(de) > 0, adj_1, 0) - siginv[1, 2] + K) / (siginv[2, 2] - ifelse(length(de) > 0, adj_2, 0))
    intervals <- c(intervals, L, U)
  }
  
  # Merge overlapping intervals
  while (length(intervals) != 0) {
    pos <- which(intervals[c(TRUE, FALSE)] <= intervals[which.min(intervals) + 1]) * 2 - 1
    pos2 <- which(intervals[c(TRUE, FALSE)] <= max(intervals[pos + 1])) * 2 - 1
    
    while (length(pos) != length(pos2)) {
      pos <- pos2
      pos2 <- which(intervals[c(TRUE, FALSE)] <= max(intervals[pos + 1])) * 2 - 1
    }
    
    returnintervals <- c(returnintervals, min(intervals), max(intervals[pos + 1]))
    intervals <- intervals[-c(pos, pos + 1)]
  }
  
  return(returnintervals)
}
