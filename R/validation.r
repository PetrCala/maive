#' Validate MAIVE Data Structure and Scientific Constraints
#'
#' Internal function that validates the input data frame for MAIVE analysis.
#' Checks data structure, numeric types, minimum sample size, and other
#' scientific requirements.
#'
#' @param dat Data frame with required columns: bs, sebs, Ns
#' @param studylevel Study-level option (0-3) the analysis will run with. The
#'   degrees-of-freedom rule for study identifiers only applies when study
#'   dummies are fitted (odd \code{studylevel}); with \code{NULL} (unknown)
#'   the rule is applied whenever a \code{study_id} column is present.
#' @return Cleaned data frame (invisible) with empty rows removed
#' @keywords internal
#' @noRd
validate_maive_data <- function(dat, studylevel = NULL) {
  # Data frame type check
  if (!is.data.frame(dat)) {
    cli::cli_abort("Input 'dat' must be a data frame.", call. = FALSE)
  }

  # Required columns (check before row count so we can clean properly)
  required_cols <- c("bs", "sebs", "Ns")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  # Remove empty rows
  original_nrow <- nrow(dat)
  dat <- maive_drop_empty_rows(dat)

  # Numeric types (strict validation, no silent coercion) after removing empty rows
  for (col in required_cols) {
    if (!is.numeric(dat[[col]])) {
      cli::cli_abort(
        sprintf("Column '%s' must be numeric. Convert your data before calling MAIVE.", col),
        call. = FALSE
      )
    }
    if (any(is.na(dat[[col]]))) {
      cli::cli_abort(
        sprintf("Column '%s' must not contain missing values.", col),
        call. = FALSE
      )
    }
    if (any(!is.finite(dat[[col]]))) {
      cli::cli_abort(
        sprintf("Column '%s' must contain finite values only.", col),
        call. = FALSE
      )
    }
  }

  # Minimum rows after cleaning (scientific requirement for reliable estimation)
  if (nrow(dat) < 4) {
    if (nrow(dat) < original_nrow) {
      # Some rows were removed, use specific message
      cli::cli_abort(
        "After removing empty rows, insufficient data remains (< 4 observations).",
        call. = FALSE
      )
    } else {
      # No rows were removed, use general message
      cli::cli_abort(
        "Insufficient data: MAIVE requires at least 4 observations for reliable estimation.",
        call. = FALSE
      )
    }
  }

  # Positive standard errors (on non-missing values)
  invalid_se <- !is.na(dat$sebs) & dat$sebs <= 0
  if (any(invalid_se)) {
    n_invalid <- sum(invalid_se)
    cli::cli_abort(
      sprintf(
        "Standard errors ('sebs') must be positive. Found %d non-positive value%s.",
        n_invalid,
        if (n_invalid == 1) "" else "s"
      ),
      call. = FALSE
    )
  }

  # Positive sample sizes (on non-missing values)
  invalid_n <- !is.na(dat$Ns) & dat$Ns <= 0
  if (any(invalid_n)) {
    n_invalid <- sum(invalid_n)
    cli::cli_abort(
      sprintf(
        "Sample sizes ('Ns') must be positive. Found %d non-positive value%s.",
        n_invalid,
        if (n_invalid == 1) "" else "s"
      ),
      call. = FALSE
    )
  }

  # Study ID degrees of freedom check. The rule protects the study dummies
  # (one regressor per study), so it does not apply when the identifier is
  # unused (studylevel = 0) or only drives clustering (studylevel = 2); a
  # positionally inferred column must not abort an analysis that ignores it.
  uses_dummies <- is.null(studylevel) || (as.integer(studylevel) %% 2L) == 1L
  if ("study_id" %in% names(dat) && uses_dummies) {
    n_studies <- length(unique(dat$study_id))
    min_required <- n_studies + 3
    if (nrow(dat) < min_required) {
      cli::cli_abort(
        sprintf(
          "Insufficient degrees of freedom: %d observations with %d unique studies requires at least %d rows.",
          nrow(dat),
          n_studies,
          min_required
        ),
        call. = FALSE
      )
    }
  }

  invisible(dat)
}

#' Drop completely empty rows from a data frame
#'
#' Rows where every column is NA carry no information and are removed before
#' column resolution, which rejects missing values in the mapped columns.
#'
#' @param dat Data frame
#' @return Data frame without completely empty rows (row names reset)
#' @keywords internal
#' @noRd
maive_drop_empty_rows <- function(dat) {
  if (ncol(dat) == 0L || nrow(dat) == 0L) {
    return(dat)
  }
  empty <- rowSums(is.na(dat)) == ncol(dat)
  n_removed <- sum(empty)
  if (n_removed > 0L) {
    dat <- dat[!empty, , drop = FALSE]
    rownames(dat) <- NULL
    cli::cli_alert_info(
      sprintf(
        "Removed %d completely empty row%s",
        n_removed,
        if (n_removed == 1L) "" else "s"
      )
    )
  }
  dat
}

#' Resolve column mappings for MAIVE inputs
#'
#' Allows users to specify custom column names for estimate, SE, Ns, and study_id.
#' Falls back to the default names (bs, sebs, Ns, study_id) when mappings are
#' not provided. When no study_id mapping is given and no column is named
#' study_id, the fourth column is used positionally with a warning naming the
#' column, provided it was not consumed by the estimate, se, or n mapping. When
#' it was, the single remaining column is used instead (still with the
#' warning); if several columns remain the caller is asked to name study_id
#' rather than have one guessed, unless studylevel is 0 and no identifier is
#' needed. Performs presence, numeric, non-missing, and finite validation.
#'
#' @param studylevel Study-level option (0-3). With 0 an ambiguous fallback is
#'   skipped silently instead of raising an error, since the identifier is not
#'   used. \code{NULL} (unknown) behaves like a non-zero value.
#' @keywords internal
#' @noRd
resolve_maive_columns <- function(dat, estimate = NULL, se = NULL, n = NULL, study_id = NULL, # nolint: object_name_linter.
                                  studylevel = NULL) {
  col_for <- function(arg, default_name) {
    if (is.null(arg) || is.na(arg) || identical(arg, "")) {
      return(default_name)
    }
    as.character(arg)
  }

  est_col <- col_for(estimate, "bs")
  se_col <- col_for(se, "sebs")
  n_col <- col_for(n, "Ns")

  required <- c(est_col, se_col, n_col)
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    cli::cli_abort(
      sprintf("Missing required columns: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  check_vector <- function(vec, name) {
    if (!is.numeric(vec)) {
      cli::cli_abort(sprintf("Column '%s' must be numeric.", name), call. = FALSE)
    }
    if (any(is.na(vec))) {
      cli::cli_abort(sprintf("Column '%s' must not contain missing values.", name), call. = FALSE)
    }
    if (any(!is.finite(vec))) {
      cli::cli_abort(sprintf("Column '%s' must contain finite values only.", name), call. = FALSE)
    }
  }

  bs <- dat[[est_col]]
  sebs <- dat[[se_col]]
  Ns <- dat[[n_col]]

  check_vector(bs, est_col)
  check_vector(sebs, se_col)
  check_vector(Ns, n_col)

  # Study ID handling: use the mapped name if provided; otherwise a column named
  # study_id; otherwise fall back to a positional column with a warning.
  if (!is.null(study_id) && !identical(study_id, "")) {
    study_col <- as.character(study_id)
    if (!study_col %in% names(dat)) {
      cli::cli_abort(sprintf("Column '%s' (study_id) is missing.", study_col), call. = FALSE)
    }
    studyid <- dat[[study_col]]
  } else if ("study_id" %in% names(dat)) {
    studyid <- dat[["study_id"]]
  } else {
    fallback_col <- maive_positional_study_column(dat, est_col, se_col, n_col, studylevel)
    if (!is.null(fallback_col)) {
      which_col <- if (identical(fallback_col, names(dat)[[4L]])) "the fourth column" else "the only unmapped column"
      cli::cli_warn(
        c(
          "No 'study_id' column found; using {which_col} ('{fallback_col}') as the study identifier.",
          "i" = "Pass study_id = \"{fallback_col}\" to confirm, or drop the column if it is not a study identifier."
        ),
        call. = FALSE
      )
      studyid <- dat[[fallback_col]]
    } else {
      studyid <- NULL
    }
  }

  # Validate study_id if present
  if (!is.null(studyid)) {
    if (any(is.na(studyid))) {
      cli::cli_abort("Column 'study_id' must not contain missing values.", call. = FALSE)
    }
    if (!is.atomic(studyid)) {
      cli::cli_abort("Column 'study_id' must be a vector.", call. = FALSE)
    }
    if (length(studyid) != length(bs)) {
      cli::cli_abort("Column 'study_id' must have the same length as the estimates.", call. = FALSE)
    }
  }

  # Rebuild a clean data.frame with standardized names for downstream use
  cleaned <- data.frame(
    bs = bs,
    sebs = sebs,
    Ns = Ns,
    stringsAsFactors = FALSE
  )
  if (!is.null(studyid)) {
    cleaned$study_id <- studyid
  }

  list(dat = cleaned)
}

#' Pick the positional study identifier column
#'
#' Keeps the documented fourth-column convention when that column is still
#' free, that is, not consumed by the \code{estimate}, \code{se}, or \code{n}
#' mapping. When the fourth column was consumed, the single remaining column is
#' used because it is unambiguous. When several remain, the function refuses to
#' guess: it aborts naming the candidates, unless \code{studylevel} is 0, in
#' which case no identifier is needed and \code{NULL} is returned silently.
#'
#' @return The column name, or \code{NULL} when there is no candidate
#' @keywords internal
#' @noRd
maive_positional_study_column <- function(dat, est_col, se_col, n_col, studylevel = NULL) {
  candidates <- setdiff(names(dat), c(est_col, se_col, n_col))
  fourth <- if (ncol(dat) >= 4L) names(dat)[[4L]] else NULL

  if (!is.null(fourth) && fourth %in% candidates) {
    return(fourth)
  }
  if (length(candidates) == 1L) {
    return(candidates[[1L]])
  }
  if (length(candidates) == 0L) {
    return(NULL)
  }
  if (!is.null(studylevel) && as.integer(studylevel) == 0L) {
    return(NULL)
  }
  cli::cli_abort(
    c(
      "No 'study_id' column found, and the fourth column ('{fourth}') is already mapped as {.arg estimate}, {.arg se}, or {.arg n}.",
      "i" = "Several columns could be the study identifier: {.val {candidates}}.",
      "i" = "Pass study_id = \"<column>\" to name it, or use studylevel = 0 if the data has no study structure."
    ),
    call. = FALSE
  )
}

#' Validate MAIVE Parameter Values
#'
#' Internal function that validates parameter values are within expected ranges.
#'
#' @param method Method choice (1=PET, 2=PEESE, 3=PET-PEESE, 4=EK)
#' @param weight Weighting scheme (0=equal, 1=standard, 2=adjusted, 3=study)
#' @param instrument Variance instrumentation (0=disabled, 1=enabled)
#' @param studylevel Study-level effects (0-3)
#' @param SE Standard error treatment (0-3)
#' @param AR Anderson-Rubin confidence interval (0=no, 1=yes)
#' @return TRUE invisibly if all validations pass
#' @keywords internal
#' @noRd
validate_maive_parameters <- function(method, weight, instrument, studylevel, SE, AR) {
  # Method: 1=PET, 2=PEESE, 3=PET-PEESE, 4=EK
  if (!method %in% c(1, 2, 3, 4)) {
    cli::cli_abort(
      "Parameter 'method' must be 1 (PET), 2 (PEESE), 3 (PET-PEESE), or 4 (EK).",
      call. = FALSE
    )
  }

  # Weight: 0=equal, 1=standard, 2=adjusted, 3=study
  if (!weight %in% c(0, 1, 2, 3)) {
    cli::cli_abort(
      "Parameter 'weight' must be 0 (equal), 1 (standard), 2 (adjusted), or 3 (study).",
      call. = FALSE
    )
  }

  # Binary/categorical parameters
  if (!instrument %in% c(0, 1)) {
    cli::cli_abort("Parameter 'instrument' must be 0 or 1.", call. = FALSE)
  }

  if (!studylevel %in% c(0, 1, 2, 3)) {
    cli::cli_abort("Parameter 'studylevel' must be 0, 1, 2, or 3.", call. = FALSE)
  }

  if (!SE %in% c(0, 1, 2, 3)) {
    cli::cli_abort("Parameter 'SE' must be 0, 1, 2, or 3.", call. = FALSE)
  }

  if (!AR %in% c(0, 1)) {
    cli::cli_abort("Parameter 'AR' must be 0 or 1.", call. = FALSE)
  }

  invisible(TRUE)
}

#' Normalize MAIVE options and apply cross-parameter guardrails
#'
#' Coerces parameter inputs to integers, maps string first-stage names, applies
#' AR/instrument guardrails, and returns a normalized option list for downstream
#' pipeline use. Relies on existing data/parameter validators.
#'
#' @param dat Data frame with required columns: bs, sebs, Ns (or mapped), study_id optional
#' @param estimate Optional column name to use instead of 'bs'
#' @param se Optional column name to use instead of 'sebs'
#' @param n Optional column name to use instead of 'Ns'
#' @param study_id Optional column name for study identifiers (a column named
#'   study_id is used when absent; otherwise the 4th column, with a warning,
#'   unless that column is mapped as estimate/se/n, in which case the single
#'   remaining column is used or, when several remain, an error asks for it)
#' @param method Method choice (1=PET, 2=PEESE, 3=PET-PEESE, 4=EK)
#' @param weight Weighting scheme (0=equal, 1=standard, 2=adjusted, 3=study)
#' @param instrument Variance instrumentation (0=disabled, 1=enabled)
#' @param studylevel Study-level effects (0-3)
#' @param SE Standard error treatment (0-3)
#' @param AR Anderson-Rubin confidence interval (0=no, 1=yes)
#' @param first_stage First-stage specification: 1/\"log\" (default) or 0/\"levels\"
#' @return List of normalized options and derived metadata
#' @keywords internal
#' @noRd
normalize_maive_options <- function(dat,
                                    method,
                                    weight,
                                    instrument,
                                    studylevel,
                                    SE,
                                    AR,
                                    first_stage = 1L,
                                    estimate = NULL,
                                    se = NULL,
                                    n = NULL,
                                    study_id = NULL) {
  dat <- as.data.frame(dat)
  # Resolve custom column names first, then validate the resolved frame.
  # Completely empty rows are dropped before resolution because resolution
  # rejects missing values in the mapped columns rather than cleaning them.
  scalar_int <- function(value, name) {
    if (length(value) != 1L || is.na(value)) {
      cli::cli_abort(sprintf("Parameter '%s' must be a single non-missing value.", name), call. = FALSE)
    }
    as.integer(value)
  }

  # studylevel is needed before the data is resolved: it decides whether an
  # ambiguous positional study_id fallback matters and whether the
  # degrees-of-freedom rule for study dummies applies
  studylevel <- scalar_int(studylevel, "studylevel")

  dat <- maive_drop_empty_rows(dat)
  resolved <- resolve_maive_columns(
    dat,
    estimate = estimate, se = se, n = n, study_id = study_id,
    studylevel = studylevel
  )
  dat <- validate_maive_data(resolved$dat, studylevel = studylevel)

  normalize_first_stage <- function(first_stage) {
    if (missing(first_stage) || is.null(first_stage)) {
      return(1L)
    }
    if (is.character(first_stage)) {
      match_idx <- match(tolower(first_stage), c("levels", "log"))
      if (is.na(match_idx)) {
        cli::cli_abort("Parameter 'first_stage' must be one of 'levels' or 'log'.", call. = FALSE)
      }
      return(as.integer(match_idx - 1L))
    }
    scalar_int(first_stage, "first_stage")
  }

  apply_guardrails <- function(method, weight, instrument, AR, dat) {
    # Disable AR when incompatible with method/weights/instrument choice
    if (method == 4L || weight == 1L || instrument == 0L) {
      AR <- 0L
    }

    # Disable instrumentation (and AR) if Ns has no variation
    if (instrument == 1L) {
      # Check minimum sample size for IV estimation
      min_n_iv <- 10L
      if (nrow(dat) < min_n_iv) {
        cli::cli_warn(
          "Sample size ({nrow(dat)}) is small for IV estimation. Results may be unreliable. Consider using instrument=0 for small samples.",
          call. = FALSE
        )
      }

      Ns <- dat$Ns
      unique_Ns <- unique(stats::na.omit(Ns))
      if (length(unique_Ns) < 2L) {
        cli::cli_warn("Ns has no variation; disabling variance instrumentation (instrument=0).")
        instrument <- 0L
        AR <- 0L
      }
    }
    list(instrument = instrument, AR = AR)
  }

  method <- scalar_int(method, "method")
  weight <- scalar_int(weight, "weight")
  instrument <- scalar_int(instrument, "instrument")
  SE <- scalar_int(SE, "SE")
  first_stage <- normalize_first_stage(first_stage)

  guardrails <- apply_guardrails(method, weight, instrument, AR, dat)
  instrument <- guardrails$instrument
  AR <- guardrails$AR

  if (!first_stage %in% c(0L, 1L)) {
    cli::cli_abort("Parameter 'first_stage' must be 0 (levels) or 1 (log).", call. = FALSE)
  }

  type_map <- c("CR0", "CR1", "CR2")
  type_choice <- if (SE == 3L) "CR0" else type_map[SE + 1L]
  first_stage_type <- c("levels", "log")[first_stage + 1L]

  validate_maive_parameters(method, weight, instrument, studylevel, SE, AR)

  list(
    dat = dat,
    method = method,
    weight = weight,
    instrument = instrument,
    studylevel = studylevel,
    SE = SE,
    AR = AR,
    type_choice = type_choice,
    alpha_s = 0.05,
    first_stage = first_stage,
    first_stage_type = first_stage_type
  )
}

#' @keywords internal
maive_prepare_data <- function(dat, studylevel) {
  bs <- dat$bs
  sebs <- dat$sebs
  Ns <- dat$Ns

  M <- length(bs)

  cluster <- studylevel %/% 2L
  dummy <- studylevel %% 2L

  if ("study_id" %in% names(dat)) {
    studyid <- dat[["study_id"]]
  } else {
    if (studylevel > 0L) {
      cli::cli_warn(
        c(
          "studylevel = {studylevel} has no effect because the data has no study identifier.",
          "i" = "Add a 'study_id' column or pass study_id = \"<column>\" to fit study dummies or cluster by study."
        ),
        call. = FALSE
      )
    }
    studyid <- seq_len(M)
    dummy <- 0L
    cluster <- 0L
  }

  D <- maive_center_dummy_matrix(studyid)
  g <- if (cluster == 0L) seq_len(M) else studyid
  dat$g <- g

  list(
    dat = dat,
    bs = bs,
    sebs = sebs,
    Ns = Ns,
    M = M,
    studyid = studyid,
    dummy = dummy,
    cluster = cluster,
    g = g,
    D = D
  )
}
