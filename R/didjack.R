#' Cluster Jackknife (CV3) for Callaway-Sant'Anna DID
#'
#' Computes CV3 (delete-one-cluster) standard errors and inference
#' for an aggregated DID object from \pkg{did}.
#'
#' @details Updated: 2025-07-18
#' @param ag An `AGGTEobj` from [did::aggte()] (e.g. `aggte(gt, type="simple")`).
#' @param level Confidence level as a proportion (default 0.95). Values >1 are treated as percentages.
#' @param cluster Optional name of an alternative clustering variable (string). If `NULL`
#'   the original `clustervars` used in `att_gt()` is used.
#'
#' @return An object of class `didjack` with elements:
#'   \item{ATT}{Overall ATT}
#'   \item{cv3se, cv3t, cv3p, cv3lci, cv3uci, cv3df}{Inference components}
#'   \item{reps}{Number of jackknife replicates}
#'   \item{singclust}{Indicator for single treated cluster case}
#'   \item{atts}{Vector of delete-one ATTs}
#'   \item{clustvar, orig_clust}{Cluster variable actually used / original one}
#' @examples
#' \dontrun{
#' gt <- did::att_gt(yname="y", tname="year", gname="g",
#'                   data=df, clustervars="state", panel=FALSE,
#'                   est_method="reg", bstrap=FALSE)
#' ag <- did::aggte(gt, type="simple")
#' res <- didjack(ag)
#' summary(res)
#' }
#' @importFrom utils flush.console
#' @export

didjack <- function(ag, level = 0.95, cluster = NULL) {

  ## 0. Input check
  if (!inherits(ag, "AGGTEobj"))
    stop("didjack: pass the AGGTEobj returned by did::aggte(), run did::aggte() before running didjack, or update did to latest version")

  ## 1. Recover original att_gt() call
  atcall <- ag[["DIDparams"]][["call"]]
  data   <- eval(atcall$data, envir = parent.frame())

  orig_clust <- as.character(atcall$clustervars)
  clust      <- if (is.null(cluster)) orig_clust else cluster
  if (!clust %in% names(data))
    stop("didjack: cluster variable '", clust, "' not found in data.")

  gname <- as.character(atcall$gname)

  ## 2. Cluster structure
  cl_i     <- data[[clust]]
  clist    <- sort(unique(cl_i))
  G        <- length(clist)
  ever     <- unique(cl_i[data[[gname]] != 0])
  singclust<- length(ever) < 2L   # single treated cluster?

  ## 3. Overall ATT
  ATT  <- ag$overall.att
  type <- ag$type

  ## 4. Template args for refits
  baseargs <- as.list(atcall)[-1]

  ## 5. Leave-one-cluster-out jackknife
  atts <- numeric(0)
  cat("Jackknife (CV3): ")
  for (c in clist) {
    if (singclust && c %in% ever) next
    cat("."); flush.console()
    subdat        <- data[cl_i != c, , drop = FALSE]
    args          <- baseargs
    args$data     <- subdat
    gt_sub <- suppressWarnings(do.call(did::att_gt, args))
    atts   <- c(atts, did::aggte(gt_sub, type = type)$overall.att)
  }
  cat("\n")

  reps <- length(atts)

  if (reps == 0)
    stop("didjack: no jackknife replicates (check cluster variable).")

  ## 6. CV3 variance & inference
  df    <- reps - 1L
  cv3se <- sqrt(df / reps * sum((atts - ATT)^2))
  cv3t  <- ATT / cv3se
  cv3p  <- 2 * stats::pt(-abs(cv3t), df)
  crit  <- stats::qt(1 - (1 - level)/2, df)
  cv3ci <- ATT + c(-1, 1) * crit * cv3se

  ## 7. Return
  res <- list(
    ATT       = ATT,
    cv3se     = cv3se,
    cv3t      = cv3t,
    cv3p      = cv3p,
    cv3lci    = cv3ci[1],
    cv3uci    = cv3ci[2],
    cv3df     = df,
    reps      = reps,
    singclust = singclust,
    atts      = atts,
    clustvar  = clust,
    orig_clust= orig_clust,
    call      = match.call()
  )
  class(res) <- "didjack"
  res
}

#' @export
summary.didjack <- function(object, digits = 4, ...) {
  with(object, {
    cat("Cluster Jackknife (CV3) for DID\n")
    if (!identical(clustvar, orig_clust))
      cat(sprintf("Clustering on: %s (override of %s)\n", clustvar, orig_clust))
    else
      cat(sprintf("Clustering on: %s\n", clustvar))
    cat(sprintf("ATT           : %.*f\n", digits, ATT))
    cat(sprintf("SE (CV3)      : %.*f\n", digits, cv3se))
    cat(sprintf("t-stat (df=%d) : %.*f\n", cv3df, digits, cv3t))
    cat(sprintf("p-value        : %.4f\n", cv3p))
    cat(sprintf("95%% CI        : [%.4f , %.4f]\n", cv3lci, cv3uci))
    cat(sprintf("Replicates     : %d\n", reps))
    if (singclust) cat("Note: single treated cluster - treated cluster skipped in jackknife.\n")
  })
  invisible(object)
}

#' @export
print.didjack <- function(x, ...) summary(x, ...)



