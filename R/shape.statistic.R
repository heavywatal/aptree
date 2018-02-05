#' @rdname statistic
#' @export
normalize_shape_statistic = function(x, tree, norm=c("pda", "yule")) {
  norm = match.arg(norm)
  n <- nrow(tree$merge)
  if (norm == "pda") {
    (x - 2.03 * (n + 1) + 3.545 * sqrt(n)) / sqrt(1.570 * n * log(n) - 5.674 * n + 3.602 * sqrt(n) + 14.915)
  } else {
    (x - 1.204 * (n + 1) + log(n) + 2) / sqrt(0.168 * (n + 1) - 0.71)
  }
}

#' @rdname statistic
#' @export
shape.statistic = function(tree, norm=NULL) {
  clades = smaller.clade.spectrum(tree)
  res = sum(log(clades[, 1] - 1))
  if (is.null(norm)) {
    res
  } else {
    normalize_shape_statistic(res, tree, norm)
  }
}
