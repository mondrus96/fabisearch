#===========================================================================
# Finds significant splits using the permutation distribution

#' sign_splits
#' @description This function is embedded in the main FaBiSearch function. It serves to compare the loss for the refitted splits to those of the permuted splits
#' and determine whether the loss is significantly lower for the refitted splits.
#'
#' @param orig.splits Matrix of candidate change points (\code{T.split}) and the change in loss (\code{chg.loss})
#' @param refit.splits List of refitted split loss values
#' @param perm.distr List of permuted split loss values
#' @param alpha Cuttoff value for statistical inference (e.g., 0.05 or 0.01)
#' @param test.type Type of statistical test to run, "t-test" for simple t-test (using \code{alpha}, returns \code{TRUE/FALSE}), or "pval.t-test" to
#' return the probability value for the t-test (no cutoff)
#'
#' @return \code{final.splits}, a matrix with columns \code{T.split} for the position of the candidate change point, and \code{significant} for result of the two sided
#' \eqn{t}-test
#' @export
#'
#' @examples

sign_splits = function(orig.splits, refit.splits, perm.distr, alpha, test.type){

  # refit.splits  = loss for refitted split values
  # perm.distr    = permutation distribution
  # alpha         = level of significance for statistical inference
  # test.type     = type of statistical test to run -> "ks" for the Kolmogorov-Smirnov, and "wilcox" for the Wilcoxon test

  # Find the splits where the loss metric is reduced, define the final.splits dataframe
  reduced.splits = sort(orig.splits[orig.splits$chg.loss < 0, ]$T.split)
  final.splits   = c()

  for (ji in 1:length(refit.splits)){

    # Perform the statistical test
    if (test.type == "t-test"){
      stat.result = t.test(refit.splits[[ji]], perm.distr[[ji]], alternative = "less")
    } else if (test.type == "pval.t-test"){
      stat.result = t.test(refit.splits[[ji]], perm.distr[[ji]], alternative = "less")
    }

    # Determine if the refitted value of loss is outside of the bounds of the distribution
    if (test.type == "t-test"){
      final.splits = rbind(final.splits, data.frame(reduced.splits[ji], stat.result$p.value < alpha))
    } else if (test.type == "pval.t-test"){
      final.splits = rbind(final.splits, data.frame(reduced.splits[ji], stat.result$p.value))
    }
  }
  # Rename column
  colnames(final.splits) = c("T.split", "significant")

  return(final.splits)
}
