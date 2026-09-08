#' Convert metafor objects into MAIVE input data
#'
#' Builds the data frame that \code{\link{maive}()} and \code{\link{waive}()}
#' expect (columns \code{bs}, \code{sebs}, \code{Ns}, and optionally
#' \code{study_id}) from a \pkg{metafor} \code{escalc} data frame or an
#' \code{rma.uni} model fit.
#'
#' The conversion takes the standard error as \code{sqrt(vi)}, so metafor's
#' sampling variance cannot end up in the \code{sebs} slot by mistake. Sample
#' sizes are never inferred from the variances; doing so would reintroduce the
#' spurious precision that MAIVE corrects for. They are resolved, in order, from
#' the \code{ni} argument, an \code{ni} column (or the fit's \code{ni}), the
#' \code{ni} attribute that \code{escalc()} stamps on the effect size column,
#' and finally the exact two-group total \code{n1i + n2i}. The column and the
#' attribute are only used when they agree with \code{n1i + n2i} wherever both
#' are present: the attribute does not follow a reorder such as
#' \code{dplyr::arrange()} and \code{rma()} copies it into the fit, so a stale
#' vector of the right length falls through to the two-group total instead of
#' binding every sample size to the wrong estimate. When none of these is
#' available the function stops and names what it needs.
#'
#' For \code{rma.uni} fits the effect sizes, variances, sample sizes, and any
#' vector supplied through \code{ni} or \code{study_id} are taken through the
#' same \code{subset} and missing-value masks the fit applied, so the rows stay
#' aligned. Multivariate (\code{rma.mv}), GLMM (\code{rma.glmm}), and other
#' \code{rma} subclasses are refused rather than silently flattened, and so are
#' \code{metafor::trimfill()} fits: their imputed studies would enter MAIVE
#' as observed estimates, so pass the original \code{rma()} fit.
#'
#' @param x An \code{escalc} data frame (or any data frame with \code{yi} and
#'   \code{vi} columns) or an \code{rma.uni} fit.
#' @param ni Optional sample sizes: a numeric vector, or the name of a column in
#'   \code{x} (or in the data stored with the fit). Overrides any sample sizes
#'   carried by \code{x}.
#' @param study_id Optional study identifiers: a vector, or the name of a column
#'   in \code{x} (or in the data stored with the fit). Vectors may have either
#'   the length of the original data or the number of rows the fit kept.
#' @return A data frame with columns \code{bs}, \code{sebs}, \code{Ns}, and,
#'   when identifiers were supplied, \code{study_id}, ready for
#'   \code{\link{maive}()} or \code{\link{waive}()}.
#' @examplesIf requireNamespace("metafor", quietly = TRUE)
#' dat <- metafor::escalc(
#'   measure = "SMD",
#'   m1i = c(2.1, 2.5, 1.9, 2.4, 2.2, 2.6),
#'   sd1i = c(1.0, 1.2, 0.9, 1.1, 1.0, 1.3),
#'   n1i = c(40, 50, 30, 45, 60, 35),
#'   m2i = c(1.8, 2.0, 1.7, 2.1, 1.9, 2.2),
#'   sd2i = c(1.0, 1.1, 1.0, 1.0, 1.1, 1.2),
#'   n2i = c(42, 55, 33, 48, 58, 37)
#' )
#' maive_dat <- maive_from_metafor(dat)
#' head(maive_dat)
#'
#' fit <- metafor::rma(yi, vi, data = dat)
#' identical(maive_from_metafor(fit), maive_dat)
#'
#' result <- maive(maive_dat,
#'   method = 3, weight = 0, instrument = 1,
#'   studylevel = 0, SE = 0, AR = 0
#' )
#' @export
maive_from_metafor <- function(x, ni = NULL, study_id = NULL) {
  if (inherits(x, "rma")) {
    return(maive_from_rma(x, ni = ni, study_id = study_id))
  }
  if (is.data.frame(x)) {
    return(maive_from_escalc(x, ni = ni, study_id = study_id))
  }
  cli::cli_abort(
    "{.arg x} must be an {.cls escalc} data frame or an {.cls rma.uni} fit, not {.cls {class(x)}}.",
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
maive_from_escalc <- function(x, ni = NULL, study_id = NULL) {
  yi_name <- attr(x, "yi.names")
  vi_name <- attr(x, "vi.names")
  if (is.null(yi_name) || !yi_name[1] %in% names(x)) yi_name <- "yi"
  if (is.null(vi_name) || !vi_name[1] %in% names(x)) vi_name <- "vi"
  yi_name <- yi_name[1]
  vi_name <- vi_name[1]

  missing_cols <- setdiff(c(yi_name, vi_name), names(x))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.arg x} must carry effect sizes and sampling variances; missing column{?s}: {.val {missing_cols}}.",
      call. = FALSE
    )
  }

  yi <- x[[yi_name]]
  vi <- x[[vi_name]]
  n_rows <- nrow(x)

  # Column-name arguments resolve against the data frame itself
  ni_vec <- maive_metafor_resolve_column(ni, x, "ni", n_rows)
  study_vec <- maive_metafor_resolve_column(study_id, x, "study_id", n_rows)

  if (is.null(ni_vec)) {
    ni_vec <- maive_metafor_sample_sizes(
      ni_col = x[["ni"]],
      yi_attr = attr(yi, "ni"),
      n1i = x[["n1i"]],
      n2i = x[["n2i"]],
      n_rows = n_rows,
      source = "the escalc data"
    )
  }

  maive_metafor_assemble(
    yi = as.numeric(yi), vi = as.numeric(vi), ni = ni_vec, study_id = study_vec,
    drop_missing = TRUE
  )
}

#' @keywords internal
#' @noRd
maive_from_rma <- function(x, ni = NULL, study_id = NULL) {
  # trimfill fits inherit rma.uni, so refuse them before the generic gate: their
  # yi/vi/ni slots hold observed plus imputed studies, and the imputed rows
  # carry the sample size of the study they mirror, so nothing marks them
  if (inherits(x, "rma.uni.trimfill")) {
    cli::cli_abort(
      c(
        "{.fn metafor::trimfill} fits contain imputed studies and are not supported.",
        "i" = "MAIVE must be estimated on observed estimates. Pass the original {.fn metafor::rma} fit."
      ),
      call. = FALSE
    )
  }
  if (!inherits(x, "rma.uni")) {
    cli::cli_abort(
      c(
        "Only {.cls rma.uni} fits are supported; got {.cls {class(x)[1]}}.",
        "i" = "Multivariate, multilevel, and GLMM fits carry structure that a single effect and variance per row cannot represent. Fit an {.fn metafor::rma} model or convert the underlying {.fn metafor::escalc} data instead."
      ),
      call. = FALSE
    )
  }

  yi <- as.numeric(x$yi)
  vi <- as.numeric(x$vi)
  k <- length(yi)

  # Row masks the fit applied: subset (over the original rows), then NA removal
  subset_mask <- x$subset
  not_na <- x$not.na
  n_after_subset <- if (is.null(not_na)) k else length(not_na)
  n_original <- if (is.null(subset_mask)) n_after_subset else length(subset_mask)

  apply_masks <- function(v, name) {
    if (length(v) == k) {
      return(v)
    }
    if (length(v) == n_original) {
      if (!is.null(subset_mask)) v <- v[subset_mask]
      if (!is.null(not_na)) v <- v[not_na]
      if (length(v) == k) {
        return(v)
      }
    }
    cli::cli_abort(
      "{.arg {name}} has length {length(v)}; expected {k} (rows kept by the fit) or {n_original} (rows of the original data).",
      call. = FALSE
    )
  }

  fit_data <- x$data
  resolve <- function(arg, name) {
    if (is.null(arg)) {
      return(NULL)
    }
    if (is.character(arg) && length(arg) == 1L) {
      if (is.null(fit_data)) {
        cli::cli_abort(
          "{.arg {name}} names a column, but the fit does not store its data. Pass a vector instead.",
          call. = FALSE
        )
      }
      if (!arg %in% names(fit_data)) {
        cli::cli_abort("Column {.val {arg}} ({name}) is not in the data stored with the fit.", call. = FALSE)
      }
      return(apply_masks(fit_data[[arg]], name))
    }
    apply_masks(arg, name)
  }

  ni_vec <- resolve(ni, "ni")
  study_vec <- resolve(study_id, "study_id")

  if (is.null(ni_vec)) {
    n1i <- if (!is.null(fit_data) && !is.null(fit_data[["n1i"]])) apply_masks(fit_data[["n1i"]], "n1i") else NULL
    n2i <- if (!is.null(fit_data) && !is.null(fit_data[["n2i"]])) apply_masks(fit_data[["n2i"]], "n2i") else NULL
    ni_vec <- maive_metafor_sample_sizes(
      ni_col = x$ni,
      yi_attr = attr(x$yi, "ni"),
      n1i = n1i,
      n2i = n2i,
      n_rows = k,
      source = "the rma fit"
    )
  }

  # The fit already removed rows with missing yi or vi, so nothing to drop here
  maive_metafor_assemble(yi = yi, vi = vi, ni = ni_vec, study_id = study_vec, drop_missing = FALSE)
}

#' Resolve a column-name or vector argument against a data frame
#' @keywords internal
#' @noRd
maive_metafor_resolve_column <- function(arg, x, name, n_rows) {
  if (is.null(arg)) {
    return(NULL)
  }
  # A single string names a column when one exists; otherwise it is a value
  if (is.character(arg) && length(arg) == 1L && arg %in% names(x)) {
    return(x[[arg]])
  }
  if (length(arg) != n_rows) {
    cli::cli_abort(
      "{.arg {name}} has length {length(arg)}; expected {n_rows} to match the rows of {.arg x}.",
      call. = FALSE
    )
  }
  arg
}

#' Pick sample sizes from the sources metafor provides, in priority order
#'
#' An \code{ni} column (or the fit's \code{ni} slot) and the \code{ni}
#' attribute on the effect sizes are only accepted when they agree with
#' \code{n1i + n2i} wherever both are present. The attribute is a plain
#' attribute on the column, so a reorder such as \code{dplyr::arrange()} leaves
#' it at full length in the original order, and \code{rma()} copies the same
#' stale vector into \code{fit$ni}; length alone cannot tell these apart.
#' @keywords internal
#' @noRd
maive_metafor_sample_sizes <- function(ni_col, yi_attr, n1i, n2i, n_rows, source) {
  two_group <- NULL
  if (!is.null(n1i) && !is.null(n2i) && length(n1i) == n_rows && length(n2i) == n_rows) {
    two_group <- as.numeric(n1i) + as.numeric(n2i)
  }

  usable <- function(candidate) {
    if (is.null(candidate) || length(candidate) != n_rows) {
      return(FALSE)
    }
    if (is.null(two_group)) {
      return(TRUE)
    }
    candidate <- as.numeric(candidate)
    both <- !is.na(candidate) & !is.na(two_group)
    isTRUE(all.equal(candidate[both], two_group[both]))
  }

  if (usable(ni_col)) {
    return(as.numeric(ni_col))
  }
  if (usable(yi_attr)) {
    return(as.numeric(yi_attr))
  }
  if (!is.null(two_group)) {
    return(two_group)
  }
  cli::cli_abort(
    c(
      "No sample sizes found in {source}.",
      "i" = "MAIVE needs the per-estimate sample size and never infers it from the variance.",
      "i" = "Supply {.arg ni} (a vector or column name), or build the data with {.fn metafor::escalc} using {.arg ni} or {.arg n1i}/{.arg n2i} so the sample size is recorded."
    ),
    call. = FALSE
  )
}

#' Assemble the MAIVE data frame, dropping rows with missing effects or variances
#' @keywords internal
#' @noRd
maive_metafor_assemble <- function(yi, vi, ni, study_id, drop_missing) {
  if (any(vi < 0, na.rm = TRUE)) {
    cli::cli_abort("Sampling variances ({.arg vi}) must be non-negative.", call. = FALSE)
  }

  keep <- rep(TRUE, length(yi))
  if (drop_missing) {
    keep <- !is.na(yi) & !is.na(vi)
    n_dropped <- sum(!keep)
    if (n_dropped > 0L) {
      cli::cli_alert_info("Dropped {n_dropped} row{?s} with missing effect size or sampling variance.")
    }
  }

  out <- data.frame(
    bs = yi[keep],
    sebs = sqrt(vi[keep]),
    Ns = ni[keep],
    stringsAsFactors = FALSE
  )
  if (!is.null(study_id)) {
    out$study_id <- study_id[keep]
  }
  rownames(out) <- NULL
  out
}
