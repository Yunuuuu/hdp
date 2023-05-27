#' Freeze DP nodes
#'
#' Freezes previously active DP nodes. A frozen DP node is not included in posterior sampling,
#' but its statistics \emph{are} considered in the sampling of other active DPs.
#' This is useful for conditioning on a previous dataset. First, set up a HDP
#' for one dataset, run the posterior sampling chain, and then freeze all old nodes
#' (except the top DP). Add new DP nodes with new data and run a
#' second posterior sampling chain over the new nodes (\emph{given} the information in the frozen nodes).
#'
#' @param hdp A hdpState object
#' @param dpindex Indices of the DPs to freeze
#' @return A hdpState object with the specified DP nodes frozen. See \code{\link{hdpState-class}}
#' @seealso \code{\link{hdp_init}}, \code{\link{hdp_addconparam}}, \code{\link{hdp_adddp}},
#'  \code{\link{hdp_setdata}}, \code{\link{dp_activate}}, \code{\link{hdp_posterior}}
#' @export
dp_freeze <- function(hdp, dpindex) {
  # input checks
  assert_class(hdp, function(x) {
    is(x, "hdpState") && validObject(x)
  }, msg = "{.cls hdpState} object")
  assert_class(dpindex, function(x) {
    all(x > 0L) && is_round_integer(x) &&
      !anyDuplicated(dpindex) &&
      all(x <= hdp@numdp)
  }, msg = "positive integers no greater than {.code numdp(hdp)} with no duplicates")

  ACTIVE <- 2L
  FROZEN <- 1L

  dpindex <- sort(dpindex)

  for (kk in seq_along(dpindex)) {
    jj <- dpindex[kk]
    if (hdp@dpstate[jj] != ACTIVE) {
      stop("Can only freeze a DP that is activated")
    }

    hdp@dp[[jj]]@alpha <- numeric(0)
    hdp@dpstate[jj] <- FROZEN
  }

  # check validity and return
  if (!validObject(hdp)) warning("Not a valid hdpState object.")
  return(hdp)
}
