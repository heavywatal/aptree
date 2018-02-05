colless3 = function(n, sample_func, p=1 / 3) {
  if (n <= 2) return(0)
  l <- sample_func(n, p)
  colless3(l, sample_func, p) + colless3(n - l, sample_func, p) + abs(n - 2 * l)
}

#' @rdname colless
#' @export
colless.test = function(tree, model=c("yule", "pda"), alternative=c("less", "greater"), n.mc=500) {
  stopifnot(inherits(tree, "treeshape"))
  model = match.arg(model)
  alternative = match.arg(alternative)
  tip.number.mc <- nrow(tree$merge) + 1
  if (model == "pda") {
    trees.mc <- rtreeshape(n = n.mc, tip.number = tip.number.mc, p = 0.3, model = model)
    lind.mc <- sapply(trees.mc, FUN = colless, norm = NULL)
  } else {
    lind.mc <- NULL
    sample_func = switch(model,
      "biased" = sample_biased,
      "yule" = sample_yule,
      "aldous" = sample_aldous
    )
    for (i in 1:n.mc) {
      lind.mc <- c(lind.mc, colless3(tip.number.mc, sample_func))
    }
  }
  res <- stats::ecdf(lind.mc)
  stat <- colless(tree, norm = NULL)
  if (alternative == "less") {
    p.value <- res(stat)
  } else {
    p.value <- 1 - res(stat)
  }
  stat_norm <- colless(tree, norm = model)
  cat("Statistic = ", stat, " (normalized: ", stat_norm, ")\n", sep = "")
  cat("p-value = ", p.value, "\n", sep = "")
  c(statistic = stat_norm, p.value = p.value)
}

sample_biased = function(n, p=1 / 3) {
  sample(n - 1, 1, .5 * stats::dbinom(0:(n - 2), size = (n - 2), prob = p) + .5 * stats::dbinom(0:(n - 2), size = (n - 2), prob = (1 - p)))
}

sample_yule = function(n, p=NULL) {
  sample(n - 1, 1)
}

sample_aldous = function(n, p=NULL) {
  sample(n - 1, 1, prob = (1 / (1:(n - 1)) / ((n - 1):1)))
}
