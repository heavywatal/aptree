#' @rdname index-test
#' @export
index.test = function(tree, index=c("colless", "sackin"), model=c("yule", "pda"), alternative=c("less", "greater"), n.mc=500) {
  stopifnot(inherits(tree, "treeshape"))
  index = match.arg(index)
  model = match.arg(model)
  alternative = match.arg(alternative)
  index_func = switch(index,
    "colless" = colless,
    "sackin" = sackin
  )
  tip.number.mc <- nrow(tree$merge) + 1
  if (model == "pda") {
    trees.mc <- rtreeshape(n = n.mc, tip.number = tip.number.mc, p = 0.3, model = model)
    lind.mc <- sapply(trees.mc, FUN = index_func)
  } else {
    sample_func = switch(model,
      "biased" = sample_biased,
      "yule" = sample_yule,
      "aldous" = sample_aldous
    )
    index_func3 = switch(index,
      "colless" = colless3,
      "sackin" = sackin3
    )
    lind.mc <- NULL
    for (i in 1:n.mc) {
      lind.mc <- c(lind.mc, index_func3(tip.number.mc, sample_func))
    }
  }
  res <- stats::ecdf(lind.mc)
  stat <- index_func(tree)
  if (alternative == "less") {
    p.value <- res(stat)
  } else {
    p.value <- 1 - res(stat)
  }
  stat_norm <- normalize(stat, tree, index = index, norm = model)
  c(index = stat, statistic = stat_norm, p.value = p.value)
}

colless3 = function(n, sample_func, p=1 / 3) {
  if (n <= 2) return(0)
  l <- sample_func(n, p)
  colless3(l, sample_func, p) + colless3(n - l, sample_func, p) + abs(n - 2 * l)
}

sackin3 = function(n, sample_func, p=1 / 3) {
  if (n <= 2) return(0)
  l <- sample_func(n, p)
  sackin3(l, sample_func, p) + sackin3(n - l, sample_func, p) + n
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

#' @rdname index-test
#' @export
colless.test = function(tree, model = c("yule", "pda"), alternative = c("less", "greater"), n.mc = 500) {
  index.test(tree, index = "colless", model = model, alternative = alternative, n.mc = n.mc)
}

#' @rdname index-test
#' @export
sackin.test = function(tree, model = c("yule", "pda"), alternative = c("less", "greater"), n.mc = 500) {
  index.test(tree, index = "sackin", model = model, alternative = alternative, n.mc = n.mc)
}
