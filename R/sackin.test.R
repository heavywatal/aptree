sackin3 = function(n, sample_func, p=1 / 3) {
  if (n <= 2) return(0)
  l <- sample_func(n, p)
  sackin3(l, sample_func, p) + sackin3(n - l, sample_func, p) + n
}

#' @rdname sackin
#' @export
sackin.test = function(tree, model=c("yule", "pda"), alternative=c("less", "greater"), n.mc=500) {
  stopifnot(inherits(tree, "treeshape"))
  model = match.arg(model)
  alternative = match.arg(alternative)
  tip.number.mc <- nrow(tree$merge) + 1
  if (model == "pda") {
    trees.mc <- rtreeshape(n = n.mc, tip.number = tip.number.mc, model = model)
    lind.mc <- sapply(trees.mc, FUN = sackin, norm = NULL)
  } else {
    lind.mc <- NULL
    sample_func = switch(model,
      "biased" = sample_biased,
      "yule" = sample_yule,
      "aldous" = sample_aldous
    )
    for (i in 1:n.mc) {
      lind.mc <- c(lind.mc, sackin3(tip.number.mc, sample_func))
    }
  }
  res <- stats::ecdf(lind.mc)
  stat <- sackin(tree, norm = NULL)
  if (alternative == "less") {
    p.value <- res(stat)
  } else {
    p.value <- 1 - res(stat)
  }
  stat_norm <- sackin(tree, norm = model)
  cat("Statistic = ", stat, " (normalized: ", stat_norm, ")\n", sep = "")
  cat("p-value = ", p.value, "\n", sep = "")
  c(statistic = stat_norm, p.value = p.value)
}
