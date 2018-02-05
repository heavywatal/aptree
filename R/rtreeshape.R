#' @export
ryule <- function(n) {
  lst <- (-n):(-1)

  Data <- matrix(NA, nrow = n - 1, ncol = 2)
  X <- NULL
  S <- 0
  for (i in 1:(n - 1)) {
    Data[i, 1:2] <- xx <- sort(sample(lst, 2, replace = FALSE))
    lst <- c(lst[ (lst != xx[1]) & (lst != xx[2]) ], i)
  }
  treeshape(Data)
}

#' @export
raldous <- function(tip.number) {
  if (tip.number < 1 | tip.number != floor(tip.number)) {
    stop("tip.number must be an integer greater than 1")
  }
  rtreeshape(n = 1, tip.number = tip.number, model = "aldous")
}

#' @export
rbiased <- function(tip.number=20, p=0.3) {
  tree <- matrix(NA, nrow = tip.number - 1, ncol = 2)
  tree[tip.number - 1, ] <- c(-1, -2)
  current.tip <- -3
  tip <- c(-1, -2)
  tab.prob <- c(p, 1 - p)
  for (node in (tip.number - 2):1) {
    id <- sample(tip, 1, prob = tab.prob)
    tree[which(tree == id)] <- node
    tree[node, ] <- c(id, current.tip)
    aux = sample(c(1, 2), 1)
    if (aux == 1) {
      tab.prob <- c(tab.prob[tip != id], tab.prob[tip == id] * c(p, 1 - p))
    } else {
      tab.prob <- c(tab.prob[tip != id], tab.prob[tip == id] * c(1 - p, p))
    }
    tip <- c(tip[tip != id], id, current.tip)
    current.tip <- current.tip - 1
  }
  treeshape(tree)
}

rtreeshape2 <- function(tip.number, FUN) {
  if (tip.number < 2 | tip.number != floor(tip.number)) {
    stop("tip.number must be an integer greater than 2")
  }
  merge = matrix(NA, tip.number - 1, 2)

  merge[1, 1] = tip.number
  for (node in 1:(tip.number - 1)) {
    prob = 0:(merge[node, 1])
    for (i in 0:(merge[node, 1])) {
      prob[i + 1] = FUN((merge[node, 1]), prob[i + 1])
    }
    i = sample(0:(merge[node, 1]), 1, prob = prob)

    merge[node, 2] = i
    if (i != 1) {
      merge[node + 1, 1] = i
    }
    if ((merge[node, 1] - i) != 1) {
      merge[node + i] = (merge[node, 1] - i)
    }
  }
  # names=as.character(1:tip.number)
  res = treebalance(nodes = merge)
  res = as.treeshape(res)
  res
}

Qyule <- function(n, i) {
  if (i == 0 | i == n) {
    0
  }
  else {
    1
  }
}

Qaldous <- function(n, i) {
  if (i == 0 | i == n) {
    0
  }
  else {
    1 / (i * (n - i))
  }
}

#' @export
rtreeshape <- function(n, tip.number, p=0.3, model="", FUN="") {
  if (n == 0) return(NULL)
  if (class(tip.number) == "list") {
    tip.number = unlist(tip.number)
  }
  if (length(tip.number) > 1) {
    res = list()
    current = 1
    for (i in 1:length(tip.number)) {
      tmp = rtreeshape(n, tip.number[i], p, model, FUN)
      for (j in 1:n) {
        res[[current]] = tmp[[j]]
        current = current + 1
      }
    }
    return(res)
  } else {
    if (n != floor(n) | n < 0) {
      stop("n must be a positive integer")
    }
    if (tip.number != floor(tip.number) | tip.number < 0) {
      stop("tip.number must be a positive integer")
    }
    if (identical(FUN, "") == TRUE & model == "") {
      stop("at least one option")
    }
    if (identical(FUN, "") == FALSE & model != "") {
      stop("at most one option")
    }

    if (identical(FUN, "") == FALSE) {
      trees = list()
      for (i in 1:n) {
        tree <- rtreeshape2(tip.number, FUN)
        trees[[i]] <- tree
      }
      return(trees)
    }
    if (model == "pda") {
      trees <- list()
      for (i in 1:n) {
        trees[[i]] <- rpda(tip.number)
      }
      return(trees)
    }
    if (model == "yule") {
      trees = list()
      for (i in 1:n) {
        trees[[i]] <- ryule(tip.number)
      }
      return(trees)
    }
    if (model == "aldous") {
      trees = list()
      for (i in 1:n) {
        trees[[i]] <- rtreeshape2(tip.number, Qaldous)
      }
      return(trees)
    }
    if (model == "biased") {
      trees <- list()
      for (i in 1:n) {
        trees[[i]] <- rbiased(tip.number = tip.number, p = p)
      }
      return(trees)
    }
    stop("model incorrect")
  }
}
