#===========================================================================
# Finds significant splits using the permutation distribution

#' sign_splits
#'
#' @param orig.splits
#' @param refit.splits
#' @param perm.distr
#' @param alpha
#' @param test.type
#'
#' @return
#' @export
#'
#' @examples

sign_splits = function(orig.splits, refit.splits, perm.distr, alpha, test.type){

  # refit.splits  = loss for refitted split values
  # perm.distr    = permutation distribution
  # alpha         = level of significance for statistical inference
  # test.type   = type of statistical test to run -> "ks" for the Kolmogorov-Smirnov, and "wilcox" for the Wilcoxon test

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
