sackin3 = function(n, sample_func, p=1 / 3) {
  if (n <= 2) return(0)
  l <- sample_func(n, p)
  sackin3(l, sample_func, p) + sackin3(n - l, sample_func, p) + n
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
  list(model = model, statistic = stat_norm, p.value = p.value, alternative = alternative)
}
