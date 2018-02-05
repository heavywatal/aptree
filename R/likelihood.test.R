#' @rdname likelihood
#' @export
likelihood.test = function(tree, model=c("yule", "pda"), alternative=c("two.sided", "less", "greater")) {
  stopifnot(inherits(tree, "treeshape"))
  model = match.arg(model)
  alternative = match.arg(alternative)
  # number of internal nodes
  n <- nrow(tree$merge)
  if (n < 4) {
    stop("This test cannot be computed for trees with less than 4 leaves (negative variance)")
  }
  stat = shape.statistic(tree, norm = model)
  cat("statistic = ", stat, "\n", sep = "")
  p.value = switch(alternative,
    "two.sided" = 2 * stats::pnorm(abs(stat), lower.tail = FALSE),
    "less" = stats::pnorm(stat),
    "greater" = stats::pnorm(stat, lower.tail = FALSE)
  )
  cat("p.value = ", p.value, "\n", sep = "")
  c(statistic = stat, p.value = p.value)
}
