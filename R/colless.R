#' @rdname colless
#' @export
normalize_colless = function(ICN, tree, norm=c("pda", "yule")) {
  norm = match.arg(norm)
  leaf_nb <- nrow(tree$merge) + 1
  if (norm == "pda") {
    ICN / (leaf_nb ^ 1.5)
  } else {
    EICN <- leaf_nb * log(leaf_nb) + (0.57721566 - 1 - log(2)) * leaf_nb
    (ICN - EICN) / leaf_nb
  }
}

#' @rdname colless
#' @export
colless = function(tree, norm=NULL) {
  clades = smaller.clade.spectrum(tree)
  ICN = sum(abs(clades[, 1] - 2 * clades[, 2]))
  if (is.null(norm)) {
    ICN
  } else {
    normalize_colless(ICN, tree, norm)
  }
}
