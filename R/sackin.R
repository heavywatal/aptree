#' @rdname sackin
#' @export
normalize_sackin = function(INS, tree, norm=c("pda", "yule")) {
  norm = match.arg(norm)
  leaf_nb <- nrow(tree$merge) + 1
  if (norm == "pda") {
    INS / (leaf_nb ^ 1.5)
  } else {
    EINS <- 2 * leaf_nb * sum(1 / 2:leaf_nb)
    (INS - EINS) / leaf_nb
  }
}

#' @rdname sackin
#' @export
sackin = function(tree, norm=NULL) {
  if (identical(tree, NULL)) {
    stop("invalid tree", "\n")
  }
  clades = smaller.clade.spectrum(tree)
  INS = sum(clades[, 1])
  if (is.null(norm)) {
    INS
  } else {
    normalize_sackin(INS, tree, norm)
  }
}
