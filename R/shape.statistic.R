#' @rdname statistic
#' @export
shape.statistic = function(tree, norm=NULL) {
  clades = smaller.clade.spectrum(tree)
  stat = sum(log(clades[, 1] - 1))
  normalize_shape_statistic(stat, tree, norm=norm)
}

#' @rdname statistic
#' @export
colless = function(tree, norm=NULL) {
  clades = smaller.clade.spectrum(tree)
  stat = sum(abs(clades[, 1] - 2 * clades[, 2]))
  normalize_colless(stat, tree, norm=norm)
}

#' @rdname statistic
#' @export
sackin = function(tree, norm=NULL) {
  clades = smaller.clade.spectrum(tree)
  stat = sum(clades[, 1])
  normalize_sackin(stat, tree, norm=norm)
}

#' @rdname statistic
#' @export
normalize_shape_statistic = function(x, tree, norm=c(NA, "pda", "yule")) {
  norm = match.arg(norm)
  if (is.na(norm)) return(x)
  n <- nrow(tree$merge)
  if (norm == "pda") {
    (x - 2.03 * (n + 1) + 3.545 * sqrt(n)) / sqrt(1.570 * n * log(n) - 5.674 * n + 3.602 * sqrt(n) + 14.915)
  } else {
    (x - 1.204 * (n + 1) + log(n) + 2) / sqrt(0.168 * (n + 1) - 0.71)
  }
}

#' @rdname statistic
#' @export
normalize_colless = function(x, tree, norm=c(NA, "pda", "yule")) {
  norm = match.arg(norm)
  if (is.na(norm)) return(x)
  leaf_nb <- nrow(tree$merge) + 1
  if (norm == "pda") {
    x / (leaf_nb ^ 1.5)
  } else {
    EICN <- leaf_nb * log(leaf_nb) + (0.57721566 - 1 - log(2)) * leaf_nb
    (x - EICN) / leaf_nb
  }
}

#' @rdname statistic
#' @export
normalize_sackin = function(x, tree, norm=c(NA, "pda", "yule")) {
  norm = match.arg(norm)
  if (is.na(norm)) return(x)
  leaf_nb <- nrow(tree$merge) + 1
  if (norm == "pda") {
    x / (leaf_nb ^ 1.5)
  } else {
    EINS <- 2 * leaf_nb * sum(1 / 2:leaf_nb)
    (x - EINS) / leaf_nb
  }
}

#' @rdname statistic
#' @export
normalize = function(x, tree, index=c("shape", "colless", "sackin"), norm=c(NA, "pda", "yule")) {
  index = match.arg(index)
  switch(index,
    "shape" = normalize_shape_statistic(x, tree, norm=norm),
    "colless" = normalize_colless(x, tree, norm=norm),
    "sackin" = normalize_sackin(x, tree, norm=norm)
  )
}

#' @rdname statistic
#' @export
calc_stat = function(tree, index=c("shape", "colless", "sackin"), norm=c(NA, "pda", "yule")) {
  index = match.arg(index)
  switch(index,
    "shape" = shape.statistic(tree, norm=norm),
    "colless" = colless(tree, norm=norm),
    "sackin" = sackin(tree, norm=norm)
  )
}

#' @rdname statistic
#' @export
calc_stat_df = function(tree, index=c("shape", "colless", "sackin")) {
  index = match.arg(index)
  x = switch(index,
    "shape" = shape.statistic(tree),
    "colless" = colless(tree),
    "sackin" = sackin(tree)
  )
  x_pda = normalize(x, tree, index=index, norm="pda")
  x_yule = normalize(x, tree, index=index, norm="yule")
  tibble::tibble(
    norm = c(NA, "pda", "yule"),
    stat = c(x, x_pda, x_yule)
  )
}

#' @rdname statistic
#' @export
calc_stat_all = function(tree, index=c("shape", "colless", "sackin")) {
  purrr::map_dfr(stats::setNames(,index), ~calc_stat_df(tree, .x), .id = "index")
}
