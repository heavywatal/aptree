#' @rdname likelihood
#' @export
likelihood.test = function(tree, model=c("yule", "pda"), alternative=c("two.sided", "less", "greater")) {
  if (class(tree)[[1]] != "treeshape") {
    stop("invalid arguments")
  }
  model = match.arg(model)
  alternative = match.arg(alternative)
  # number of internal nodes
  n <- nrow(tree$merge)
  if (n < 4) {
    stop("This test cannot be computed for trees with less than 4 leaves (negative variance)")
  }
  stat = shape.statistic(tree, norm = model)
  cat("statistic = ")
  cat(stat, "\n")
  cat("p.value = ")
  if (alternative == "two.sided") {
    p.value <- 2 * stats::pnorm(abs(stat), lower.tail=FALSE)
    cat(p.value, "\n")
    cat("alternative hypothesis: the tree does not fit the model")
  }
  else if (alternative == "less") {
    p.value <- stats::pnorm(stat)
    cat(p.value, "\n")
    cat("alternative hypothesis: the tree is more balanced than predicted by the model")
  }
  else if (alternative == "greater") {
    p.value <- stats::pnorm(stat, lower.tail=FALSE)
    cat(p.value, "\n")
    cat("alternative hypothesis: the tree is less balanced than predicted by the model")
  }
  cat("\n")
  list(model = model, statistic = stat, p.value = p.value, alternative = alternative)
}
