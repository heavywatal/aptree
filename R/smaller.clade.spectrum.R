#' @rdname smaller.clade.spectrum
#' @export
smaller.clade.spectrum <- function(tree) {
  merge = tree$merge
  nrow_merge = nrow(merge)
  spec <- matrix(0, nrow = nrow_merge, ncol = 2)
  spec[nrow_merge, ] = c(2, 1)
  for (node in 2:nrow_merge) {
    lc = ifelse(merge[node, 1] < 0, 1, spec[nrow_merge - merge[node, 1] + 1, 1])
    rc = ifelse(merge[node, 2] < 0, 1, spec[nrow_merge - merge[node, 2] + 1, 1])
    spec[nrow_merge - node + 1, ] = c(rc + lc, min(rc, lc))
  }
  spec
}

#' @rdname smaller.clade.spectrum
#' @export
smaller.clade.spectrum2 <- function(tree) {
  if (identical(tree, NULL)) {
    stop("invalid tree", "\n")
  }
  options(warn = -1)

  merge <- tree$merge
  number <- as.numeric(tree$names)
  number[is.na(number)] <- 1
  number <- as.numeric(number)
  spec <- matrix(0, nrow = nrow(merge), ncol = 2)
  for (node in 1:nrow(merge)) {
    if (merge[node, 1] < 0) {
      lc <- number[abs(merge[node, 1])]
    }
    else {
      lc <- spec[nrow(merge) - merge[node, 1] + 1, 1]
    }
    if (merge[node, 2] < 0) {
      rc <- number[abs(merge[node, 2])]
    }
    else {
      rc <- spec[nrow(merge) - merge[node, 2] + 1, 1]
    }
    spec[nrow(merge) - node + 1, ] <- c(rc + lc, min(rc, lc))
  }

  options(warn = 0)
  spec
}
