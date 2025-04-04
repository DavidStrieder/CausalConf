#' Test if a causal ordering is plausible
#'
#' This (helper) function tests whether a given causal ordering is plausible based on the dual likelihood ratio test.
#' Note that a causal ordering corresponds to a partition of the variables into three groups, before, between and after the cause and effect of interest.
#'
#' @param index Integer representing the causal ordering.
#' @param d Number of variables in the system.
#' @param n Sample size.
#' @param siginv Inverse of the covariance matrix.
#' @param tresh Likelihood ratio threshold.
#' @param ord ordering of potential cause and potential effect
#' @param crit Chi-square critical value for hypothesis testing.
#' @return Vector of plausible causal orderings with corresponding likelihood value.
#' @export
orderDualLRTpartial<- function(index, d, n, siginv, tresh, ord, crit) {
  
  # Convert ordering into partition (before cause, between cause and effect, after effect)
  parentset <- unname(gtools::baseOf(index, base = 3, len = d - 2) + 1)
  pre <- (3:d)[parentset == 1]
  mid <- (3:d)[parentset == 2]
  aft <- (3:d)[parentset == 3]
  
  common_inv <- solve(siginv[c(mid, ord[2], aft), c(mid, ord[2], aft)])
  val <- siginv[ord[1], ord[1]] - siginv[ord[1], c(mid, ord[2], aft)] %*% common_inv %*% siginv[c(mid, ord[2], aft), ord[1]] + siginv[ord[2], ord[2]]
  
  if (length(aft) != 0) {
    val <- val - siginv[ord[2], aft] %*% solve(siginv[aft, aft]) %*% siginv[aft, ord[2]]
  }
  
  # Update the likelihood value based on the ordering
  update_val <- function(subset, base_set) {
    for (i in seq_along(subset)) {
      vord <- c(subset[-(1:i)], base_set)
      if (length(vord)!=0) {
        val <<- val * sqrt(siginv[subset[i], subset[i]]  - siginv[subset[i], vord] %*% solve(siginv[vord, vord]) %*% siginv[vord, subset[i]])
      } else {
        val <<- val * sqrt(siginv[subset[i], subset[i]])
      }
    }
  }
  
  if (length(pre) != 0) update_val(pre, c(ord, mid, aft))
  if (length(mid) != 0) update_val(mid, c(ord[2], aft))
  if (length(aft) != 0) update_val(aft, NULL)
  
  # test the ordering with dual likelihood ratio test
  lik <- -log(val)
  if ((tresh - lik) > (crit / (2 * n))) {
    return()
  }
  
  return(c(pre, ord[1], mid, ord[2], aft, lik))
}
