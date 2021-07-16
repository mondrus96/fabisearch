#===========================================================================
# Finds significant splits using the permutation distribution

sign_splits = function(orig.splits, refit.splits, perm.distr, alpha){

  # orig.splits  = original splits
  # refit.splits = refitted splits
  # perm.distr   = permutation distribution for all splits
  # alpha        = level of significance

  # Find the splits where the loss metric is reduced, define the final.splits dataframe
  reduced.splits = sort(orig.splits[orig.splits$chg.loss < 0, ]$T.split)
  final.splits   = c()

  if(!is.null(alpha)){
    for (ji in 1:length(refit.splits)){

      # Perform the statistical test
      stat.result = t.test(refit.splits[[ji]], perm.distr[[ji]], alternative = "less")

      # Determine if the refitted value of loss is outside of the bounds of the distribution
      if (is.numeric(alpha)){
        final.splits = rbind(final.splits, data.frame(reduced.splits[ji], stat.result$p.value < alpha))
      } else if (!is.numeric(alpha)){
        final.splits = rbind(final.splits, data.frame(reduced.splits[ji], stat.result$p.value))
      }
    }
    # Rename column
    colnames(final.splits) = c("T", "stat_test")
    rownames(final.splits) = NULL
  } else if (is.null(alpha)){

    # Create the final.splits list() variable
    final.splits = list("refit" = NULL, "perm" = NULL)

    for (ji in 1:length(refit.splits)){

      # Save the results in two tables
      final.splits$refit = cbind(final.splits$refit, refit.splits[[ji]])
      final.splits$perm = cbind(final.splits$perm, perm.distr[[ji]])
    }

    # Rename columns and rows
    colnames(final.splits$refit) = colnames(final.splits$perm) = reduced.splits
    rownames(final.splits$refit) = rownames(final.splits$perm) = NULL
  }

  return(final.splits)
}
