sackin3 = function(n, model, p=1 / 3) {
  if (n <= 2) return(0)
  l <- switch(model,
    "biased" = sample(n - 1, 1, .5 * stats::dbinom(0:(n - 2), size = (n - 2), prob = p) + .5 * stats::dbinom(0:(n - 2), size = (n - 2), prob = (1 - p))),
    "yule" = sample(n - 1, 1),
    "aldous" = sample(n - 1, 1, prob = (1 / (1:(n - 1)) / ((n - 1):1)))
  )
  sackin3(l, model, p) + sackin3(n - l, model, p) + n
}

#' @rdname sackin
#' @export
sackin.test = function(tree, model=c("yule", "pda"), alternative=c("less", "greater"), n.mc=500) {
  if (class(tree) != "treeshape") {
    stop("invalid arguments")
  }
  model = match.arg(model)
  alternative = match.arg(alternative)
  tip.number.mc <- nrow(tree$merge) + 1
  if (model == "pda") {
    cat("Computing Monte Carlo estimates...")
    trees.mc <- rtreeshape(n = n.mc, tip.number = tip.number.mc, model = model)
    lind.mc <- sapply(trees.mc, FUN = sackin, norm = NULL)
  }
  else {
    lind.mc <- NULL
    for (i in 1:n.mc) {
      lind.mc <- c(lind.mc, sackin3(tip.number.mc, model))
    }
  }
  res <- stats::ecdf(lind.mc)
  cat("\n\n")
  stat <- sackin(tree, norm = NULL)
  stat_norm <- sackin(tree, norm = model)
  cat("	Test of the ")
  cat(model)
  cat(" hypothesis using the Sackin index ", "\n")
  cat("Statistic = ")
  cat(stat, "\n")
  cat("Standardized Statistic = ")
  cat(stat_norm, "\n")
  cat("p-value = ")
  if (alternative == "less") {
    p.value <- res(stat)
    cat(p.value, "\n")
    cat("alternative hypothesis: the tree is more balanced than predicted by the ")
  }
  else {
    p.value <- 1 - res(stat)
    cat(p.value, "\n")
    cat("alternative hypothesis: the tree is less balanced than predicted by the ")
  }
  cat(model)
  cat(" model", "\n")
  cat("Note : the p-value was computed using a Monte-Carlo method")
  cat("\n")
  list(model = model, statistic = stat, p.value = p.value, alternative = alternative)
}
