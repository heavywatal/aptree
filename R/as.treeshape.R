#' @rdname treeshape
#' @export
as.treeshape = function(x, ...) UseMethod("as.treeshape")

randomize_polytomy = function(phy, bin, model, p) {
  for (i in nrow(bin):1) {
    new.nodes = bin[i, 2] - 2
    for (j in 1:nrow(phy$edge)) {
      if (phy$edge[j, 1] > bin[i, 1]) {
        phy$edge[j, 1] = phy$edge[j, 1] + new.nodes
      }
      if (phy$edge[j, 2] > bin[i, 1]) {
        phy$edge[j, 2] = phy$edge[j, 2] + new.nodes
      }
    }
    tmp = switch(model,
      "pda" = rpda(bin[i, 2]),
      "yule" = ryule(bin[i, 2]),
      "biased" = rbiased(bin[i, 2], p),
      "aldous" = raldous(bin[i, 2]),
      stop("invalid model: ", model)
    )
    tmp <- as.phylo.treeshape(tmp)
    # print(tmp$edge)
    for (j in 1:length(tmp$edge)) {
      if (tmp$edge[j] > bin[i, 2]) {
        tmp$edge[j] <- tmp$edge[j] + bin[i, 1] - bin[i, 2] - 1
      }
    }
    # print(tmp$edge)
    new.values = phy$edge[phy$edge[, 1] == bin[i, 1], 2]
    for (j in 1:bin[i, 2]) {
      tmp$edge[tmp$edge == j] = new.values[j]
    }
    # print(tmp$edge)
    current.line = new.nodes + 1
    m = matrix(nrow = nrow(phy$edge) + new.nodes, ncol = 2)
    m[1:new.nodes, ] = tmp$edge[1:new.nodes, ]
    for (line in 1:nrow(phy$edge)) {
      if (phy$edge[line, 1] == bin[i, 1]) {
        m[line + new.nodes, ] = tmp$edge[current.line, ]
        current.line = current.line + 1
      } else {
        m[line + new.nodes, ] = phy$edge[line, ]
      }
    }
    phy$edge = m
  }
  phy
}

#' @rdname treeshape
#' @export
as.treeshape.phylo = function(x, model=NULL, p=0.3, ...) {
  if (identical(model, NULL)) return(NULL)
  phy <- x
  bin = is.binary.phylo(phy)
  randomize = !identical(bin, TRUE)
  if (randomize) {
    phy = randomize_polytomy(phy, bin, model, p)
  }
  height <- (nrow(phy$edge) / 2)
  merge <- matrix(0, ncol = 2, nrow = height)
  tip.number = height + 1
  total <- 2 * tip.number

  for (i in nrow(phy$edge):1) {
    if (as.numeric(phy$edge[i, 2]) > tip.number) {
      tmp = total - phy$edge[i, 2]
    } else {
      tmp = -phy$edge[i, 2]
    }
    if (merge[total - phy$edge[i, 1], 2] == 0) {
      merge[total - phy$edge[i, 1], 2] = tmp
    } else {
      merge[total - phy$edge[i, 1], 1] = tmp
    }
  }
  res = treeshape(merge, phy$tip.label)
  if (randomize) {
    class(res) = c("treeshape", "randomized.treeshape")
  }
  res
}

#' @export
as.treeshape.treebalance <- function(x, ...) {
  tree = x
  height = nrow(tree$merge)
  merge = matrix(NA, height, 2)
  current.tip = -1
  for (node in 1:height) {
    new.node = height + 1 - node
    if (tree$merge[node, 2] == 1) {
      res1 = current.tip
      current.tip = current.tip - 1
    } else {
      res1 = new.node - 1
    }
    if ((tree$merge[node, 1] - tree$merge[node, 2]) == 1) {
      res2 = current.tip
      current.tip = current.tip - 1
    } else {
      res2 = new.node - tree$merge[node, 2]
    }
    merge[new.node, ] = c(res1, res2)
  }
  treeshape(merge)
}
