#' @export
aldous.test <- function(tree, xmin = 20, ...) {
  if (identical(tree, NULL)) {
    stop("invalid tree", "\n")
  }
  if (xmin < 0.1) {
    xmin = 0.1
  }
  clade = smaller.clade.spectrum(tree)
  if (xmin > clade[1, 1]) {
    stop("'xmin' greater than the size of the tree")
  }
  ## clade = clade[clade[, 1] > xmin, ]
  xlab = "Size of parent clade (log scale)"
  ylab = "Size of smaller daughter clade (log scale)"
  graphics::plot(
    clade, pch = 20, xlab = xlab, log = "xy", ylab = ylab,
    xlim = c(xmin, max(clade[1, ]) + 10), ...
  )
  if (xmin >= 10) {
    x = xmin + 5
    graphics::text(x = x, y = 3 / 2 + 0.1, label = "PDA model", col = 4)
    graphics::text(x = x, y = x / 4 - 1, label = "Yule model", col = 4)
  }
  else {
    x = 20
    graphics::text(x = x, y = 3 / 2 + 0.1, label = "PDA model", col = 4)
    graphics::text(x = x, y = x / 4 + 2, label = "Yule model", col = 4)
  }
  mod = quantreg::rq(clade[, 2] ~ clade[, 1])
  x = sort(clade[, 1])
  coef = stats::coefficients(mod)
  graphics::abline(coef, untf = TRUE, col = 3)
  graphics::abline(0, 1 / 4, untf = TRUE)
  graphics::abline(3 / 2, 0, untf = TRUE)
  graphics::abline(v = 30, lty = 3)
}
