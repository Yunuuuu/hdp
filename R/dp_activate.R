#' Activate DP nodes
#'
#' Specifiy the number of starting clusters, and activate the DP nodes to be
#' included in the posterior sampling process (\code{\link{hdp_posterior}}).
#' When initialised, the DP nodes are 'heldout' (not available for posterior sampling).
#'
#' Note that this step can be slow and memory-intensive for very large datasets.
#'
#' @param hdp A hdpState object
#' @param dpindex Indices of the DPs to activate (include all parent DPs)
#' @param initcc Number of data clusters to start with
#'  (every data item is randomly assigned to a cluster to start with)
#' @param seed The (integer) seed that can be set to reproduce output. Default is a
#'  random seed from 1 -- 10^7, reported in the output.
#' @return A hdpState object with activated DPs and an initial random cluster allocation for each data item. See \code{\link{hdpState-class}}
#' @seealso \code{\link{hdp_init}}, \code{\link{hdp_addconparam}}, \code{\link{hdp_adddp}},
#'  \code{\link{hdp_setdata}}, \code{\link{hdp_posterior}}
#' @export
#' @examples
#' my_hdp <- hdp_init(ppindex = 0, cpindex = 1, hh = rep(1, 6), alphaa = rep(1, 3), alphab = rep(2, 3))
#' my_hdp <- hdp_adddp(my_hdp, 2, 1, 2)
#' my_hdp <- hdp_adddp(my_hdp, 10, c(rep(2, 5), rep(3, 5)), 3)
#' my_hdp <- hdp_setdata(my_hdp, 4:13, example_data_hdp)
#' # active all DPs and start with two data clusters
#' my_hdp <- dp_activate(my_hdp, 1:13, 2)
dp_activate <- function(hdp, dpindex, initcc, seed = sample(1:10^7, 1)) {
  # input checks
  assert_class(hdp, function(x) {
    is(x, "hdpState") && validObject(x)
  }, msg = "{.cls hdpState} object")

  assert_class(dpindex, function(x) {
    all(x > 0L) && is_round_integer(x) &&
      !anyDuplicated(dpindex) &&
      all(x <= hdp@numdp)
  }, msg = "positive integers no greater than {.code numdp(hdp)} with no duplicates")

  assert_class(initcc, function(x) {
    all(x > 0L) &&
      is_scalar(x) &&
      is_round_integer(x)
  }, msg = "scalar positive integer")

  assert_class(seed, function(x) {
    is_scalar(x) &&
      is_round_integer(x)
  }, msg = "scalar integer")

  # set seed
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  hdp@seed_activate <- as.integer(seed)
  hdp@initcc <- as.integer(initcc)

  ACTIVE <- 2L
  HELDOUT <- 0L

  # initialize numclass and classqq
  if (initcc > hdp@base@numclass) {
    hdp <- qq_addclass(hdp, initcc - hdp@base@numclass)
  }

  dpindex <- sort(dpindex)
  # initialize state of HDP
  for (kk in seq_along(dpindex)) {
    jj <- dpindex[kk]
    pp <- hdp@ppindex[jj]
    cp <- hdp@cpindex[jj]
    tt <- hdp@ttindex[jj]

    if (hdp@dpstate[jj] == ACTIVE) {
      stop("Trying to activate a DP that is already activated")
    }
    if (pp > 0L) {
      if (hdp@dpstate[pp] != ACTIVE) {
        stop("Ancestors of to be activated DPs has to be already activated")
      }
    }

    if (hdp@dpstate[jj] == HELDOUT) {
      hdp@dp[[jj]]@datacc <- as.integer(ceiling(runif(hdp@dp[[jj]]@numdata) *
        hdp@base@numclass))
      hdp@base@classqq <- additems(
        hdp@base@classqq, hdp@dp[[jj]]@datacc,
        hdp@dp[[jj]]@datass
      )
      hdp@dp[[jj]]@classnd <- tabulate(hdp@dp[[jj]]@datacc,
        nbins = hdp@base@numclass + 1L
      )
      hdp@dp[[jj]]@classnt <- 0L
    }

    if (pp == 0L) {
      hdp@dp[[jj]]@beta <- rep(1L, hdp@base@numclass + 1L) /
        (hdp@base@numclass + 1L)
    } else {
      hdp@dp[[jj]]@beta <- hdp@dp[[pp]]@beta
    }
    hdp@dp[[jj]]@alpha <- hdp@conparam[[cp]]@alpha
    hdp@dpstate[jj] <- ACTIVE
  }

  for (kk in rev(seq_along(dpindex))) {
    jj <- dpindex[kk]
    pp <- hdp@ppindex[jj]
    cp <- hdp@cpindex[jj]
    tt <- hdp@ttindex[jj]
    alpha <- hdp@dp[[jj]]@alpha

    if (pp == 0L) {
      hdp@dp[[jj]]@classnt <- as.integer(hdp@dp[[jj]]@classnd > 0L)
    } else {
      hdp@dp[[pp]]@classnd <- as.integer(hdp@dp[[pp]]@classnd -
        hdp@dp[[jj]]@classnt)
      hdp@dp[[jj]]@classnt <- as.integer(randnumtable(
        alpha * hdp@dp[[pp]]@beta,
        hdp@dp[[jj]]@classnd
      ))
      hdp@dp[[pp]]@classnd <- as.integer(hdp@dp[[pp]]@classnd +
        hdp@dp[[jj]]@classnt)
    }

    hdp@conparam[[cp]]@totalnd[tt] <- as.integer(sum(hdp@dp[[jj]]@classnd))
    hdp@conparam[[cp]]@totalnt[tt] <- as.integer(sum(hdp@dp[[jj]]@classnt))
  }

  # check validity and return
  if (!validObject(hdp)) warning("Not a valid hdpState object.")
  return(hdp)
}
