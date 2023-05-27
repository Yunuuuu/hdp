# func to initialise classqq matrix
# creates a 1-column matrix of zeros (one for every mutation type)
# hh is the vector of parameters for the base distribution
newclass <- function(hh) {
  qq <- as.integer(rep(0L, length(hh)))
  dim(qq) <- c(length(hh), 1L)
  return(qq)
}

# func to add data from one sample to the classqq matrix
# qq is the classqq matrix, counts the number of data items of
#   each mutation type in each class
# cc is a numeric vector of class allocations for each data point in one sample
# ss is a numeric vector of data observations (the mutation types) in one sample
additems <- function(qq, cc, ss) {
  sc <- table(
    factor(ss, levels = seq_len(nrow(qq))),
    factor(cc, levels = seq_len(ncol(qq)))
  )
  qq <- qq + sc
  qq <- matrix(qq, nrow = dim(qq)[1L], ncol = dim(qq)[2L])
  return(qq)
}

# func to add new classes to a hdp data structure
qq_addclass <- function(hdp, newclass) {
  HELDOUT <- 0L
  ACTIVE <- 2L

  oldcl <- hdp@base@numclass
  numcl <- oldcl + newclass
  hdp@base@numclass <- as.integer(numcl)
  hdp@base@classqq <- cbind(hdp@base@classqq, matrix(0L,
    nrow = nrow(hdp@base@classqq), ncol = newclass
  ))

  for (jj in 1:hdp@numdp) {
    if (hdp@dpstate[jj] != HELDOUT) {
      hdp@dp[[jj]]@classnd[(oldcl + 1L):(numcl + 1L)] <- 0L
      hdp@dp[[jj]]@classnt[(oldcl + 1L):(numcl + 1L)] <- 0L
    }
    if (hdp@dpstate[jj] == ACTIVE) {
      hdp@dp[[jj]]@beta[(oldcl + 1L):(numcl + 1L)] <- hdp@dp[[jj]]@beta[
        oldcl + 1L
      ] *
        randstick(hdp@dp[[jj]]@alpha, newclass + 1L)
    }
  }
  return(hdp)
}

# func to randomly assign a number of tables
randnumtable <- function(weights, maxtable) {
  numtable <- rep(0L, length(maxtable))
  B <- unique(sort(maxtable))
  J <- match(maxtable, B)

  weights <- log(weights)

  for (ii in seq_along(B)) {
    maxtable <- B[ii]
    if (maxtable > 0L) {
      mm <- 1:maxtable
      stirnum <- stirling(maxtable)
      for (jj in which(J == ii)) {
        clik <- mm * weights[jj]
        clik <- cumsum(stirnum * exp(clik - max(clik)))
        numtable[jj] <- 1L + sum(runif(1L) * clik[maxtable] > clik)
      }
    }
  }
  return(numtable)
}


# func to return weights from a random stick breaking process
randstick <- function(alpha, numclass) {
  zz <- c(rbeta(numclass - 1L, 1L, alpha), 1L) # proportion of stick broken off
  beta <- zz * cumprod(c(1L, 1L - zz[1:(numclass - 1L)])) # amount of stick remaining
  return(beta)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Report if an argument is a specific class
#'
#' @keywords internal
#' @noRd
assert_class <- function(x, is_class, msg, cross_msg = "{.cls {class(x)}} object", null_ok = FALSE, arg = substitute(x), call = parent.frame()) {
  if (is_scalar_character(is_class)) {
    class <- is_class
    is_class <- function(x) {
      inherits(x, what = class)
    }
    if (missing(msg)) {
      msg <- "{.cls {class}} object"
    }
  }
  if (null_ok) {
    msg <- paste(msg, "or {.code NULL}", sep = " ")
  }
  msg <- sprintf("{.arg {arg}} must be a %s", msg)
  is_right_class <- is_class(x)
  # is_class sometimes return `TRUE` for`NULL`
  if (is.null(x) && !is_right_class) {
    if (!null_ok) {
      cli::cli_abort(c(msg,
        "x" = "You've supplied a {.code NULL}"
      ), call = call)
    }
  } else if (!is_right_class) {
    if (!is.null(cross_msg)) {
      cross_msg <- sprintf("You've supplied a ", cross_msg)
      msg <- c(msg, x = cross_msg)
    }
    cli::cli_abort(msg, call = call)
  }
}

#' Report if an argument has specific length
#' @keywords internal
#' @noRd
assert_length <- function(x, length, null_ok = FALSE, arg = substitute(x), call = parent.frame()) {
  length <- as.integer(length)
  if (length == 1L) {
    msg <- "{.field scalar} object"
  } else {
    msg <- "length {.val {length}} object"
  }
  if (null_ok) {
    msg <- paste(msg, "or {.code NULL}", sep = " ")
  }
  msg <- sprintf("{.arg {arg}} must be a %s", msg)
  is_right_length <- length(x) == length
  if (is.null(x) && !is_right_length) {
    if (!null_ok) {
      cli::cli_abort(c(msg,
        "x" = "You've supplied a {.code NULL}"
      ), call = call)
    }
  } else if (!is_right_length) {
    cli::cli_abort(c(msg,
      "x" = "You've supplied a length {.val {length(x)}} object"
    ), call = call)
  }
}

is_scalar <- function(x) {
  length(x) == 1L
}

is_scalar_numeric <- function(x) {
  is_scalar(x) && is.numeric(x)
}

is_scalar_character <- function(x) {
  is_scalar(x) && is.character(x)
}

is_scalar_logical <- function(x) {
  is_scalar(x) && is.logical(x)
}

is_round_integer <- function(x) {
  all(as.integer(x) == x, na.rm = TRUE)
}
