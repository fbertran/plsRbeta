#' Predict method for plsRbeta models
#'
#' This function predicts responses or component scores from a fitted
#' \code{"plsRbetamodel"} object.
#'
#' @param object an object of class \code{"plsRbetamodel"}
#' @param newdata optional new predictor data. When omitted, fitted values or
#' training scores are returned.
#' @param comps number of components to use.
#' @param type type of prediction. \code{"response"} returns predicted response
#' values and \code{"scores"} returns component scores.
#' @param methodNA method to use when \code{newdata} contains missing values.
#' \code{"adaptative"} uses the standard score computation for complete rows and
#' the missing-data score reconstruction otherwise. \code{"missingdata"} applies
#' the missing-data reconstruction to every row.
#' @param verbose should info messages be displayed ?
#' @param ... additional arguments passed to the underlying \code{predict}
#' method when a final model is available.
#' @return A numeric vector/matrix of predictions or a matrix of component
#' scores.
#' @examples
#' data("GasolineYield", package = "betareg")
#' modpls <- plsRbeta(yield ~ ., data = GasolineYield, nt = 2, modele = "pls-beta")
#' head(predict(modpls))
#' head(predict(modpls, newdata = GasolineYield[1:3, -1, drop = FALSE]))
#' head(predict(modpls, type = "scores"))
#' @export
predict.plsRbetamodel <- function(object, newdata, comps = object$computed_nt,
                                  type = c("response", "scores"),
                                  methodNA = "adaptative", verbose = TRUE, ...) {
  if (!inherits(object, "plsRbetamodel")) {
    stop("Primary argument must be a plsRbetamodel object")
  }

  type <- match.arg(type)
  methodNA <- match.arg(methodNA, c("adaptative", "missingdata"))

  if (length(comps) != 1L || !is.numeric(comps) || is.na(comps)) {
    stop("'comps' must be a single positive integer")
  }
  comps <- as.integer(comps)
  if (comps < 1L || comps > object$computed_nt) {
    stop("Cannot predict using more components than extracted.")
  }

  tt <- if (missing(newdata) || is.null(newdata)) {
    .training_scores_plsRbetamodel(object, comps = comps)
  } else {
    newdata.frame <- .predictors_plsRbetamodel(object, newdata = newdata)
    .scores_plsRbetamodel(
      object,
      newdata.frame = newdata.frame,
      comps = comps,
      methodNA = methodNA,
      verbose = verbose
    )
  }

  if (type == "scores") {
    return(tt[, seq_len(comps), drop = FALSE])
  }

  if (is.null(object$FinalModel)) {
    return(.predict_pls_response_plsRbetamodel(object, tt = tt, comps = comps))
  }

  if (inherits(object$FinalModel, "polr") && comps < object$computed_nt) {
    stop(
      "Prediction with fewer components is not implemented for ordinal plsRbetamodel objects; refit with the requested number of components."
    )
  }

  if (comps < object$computed_nt) {
    return(.predict_glm_or_beta_subset_plsRbetamodel(object, tt = tt, comps = comps))
  }

  if (missing(newdata) || is.null(newdata)) {
    if (inherits(object$FinalModel, "polr")) {
      return(stats::predict(object$FinalModel, type = "probs", ...))
    }
    return(stats::predict(object$FinalModel, type = "response", ...))
  }

  final_term <- .final_model_term_plsRbetamodel(object)
  if (is.null(final_term)) {
    stop("The stored final model does not expose a compatible prediction term.")
  }

  newdata.final <- data.frame(tt = I(tt[, seq_len(object$computed_nt), drop = FALSE]))
  names(newdata.final) <- final_term

  if (inherits(object$FinalModel, "polr")) {
    return(stats::predict(object$FinalModel, newdata = newdata.final, type = "probs", ...))
  }

  stats::predict(object$FinalModel, newdata = newdata.final, type = "response", ...)
}

.formula_plsRbetamodel <- function(object) {
  if (!is.null(object$call$object)) {
    if (inherits(object$call$object, "formula")) {
      return(object$call$object)
    }
    if (is.call(object$call$object) && identical(object$call$object[[1L]], as.name("~"))) {
      return(as.formula(object$call$object))
    }
  }
  if (!is.null(object$call$formula)) {
    if (inherits(object$call$formula, "formula")) {
      return(object$call$formula)
    }
    if (is.call(object$call$formula) && identical(object$call$formula[[1L]], as.name("~"))) {
      return(as.formula(object$call$formula))
    }
  }
  NULL
}

.predictors_plsRbetamodel <- function(object, newdata) {
  formula_object <- .formula_plsRbetamodel(object)

  if (is.null(formula_object)) {
    newdata.frame <- as.data.frame(newdata)
    if (!all(vapply(newdata.frame, is.numeric, logical(1)))) {
      stop("`newdata` must be numeric when the model is fitted with the default interface.")
    }
  } else if (inherits(object$FinalModel, "betareg")) {
    formula_object <- Formula::as.Formula(formula_object)
    mtX <- terms(formula_object, data = newdata, rhs = 1L)
    attr(mtX, "intercept") <- 0L
    mtX <- delete.response(mtX)
    mf <- model.frame(mtX, data = newdata, na.action = na.pass, drop.unused.levels = FALSE)
    newdata.frame <- if (!is.empty.model(mtX)) {
      as.data.frame(model.matrix(mtX, mf))
    } else {
      as.data.frame(matrix(, NROW(mf), 0L))
    }
  } else {
    mt <- terms(formula_object, data = newdata)
    attr(mt, "intercept") <- 0L
    mt <- delete.response(mt)
    mf <- model.frame(mt, data = newdata, na.action = na.pass, drop.unused.levels = FALSE)
    newdata.frame <- if (!is.empty.model(mt)) {
      as.data.frame(model.matrix(mt, mf))
    } else {
      as.data.frame(matrix(, NROW(mf), 0L))
    }
  }

  .align_predictors_plsRbetamodel(object, newdata.frame = newdata.frame)
}

.align_predictors_plsRbetamodel <- function(object, newdata.frame) {
  expected <- colnames(object$ExpliX)

  if (ncol(newdata.frame) == 0L) {
    stop("`newdata` must contain at least one predictor.")
  }

  if (is.null(colnames(newdata.frame))) {
    if (ncol(newdata.frame) != length(expected)) {
      stop("`newdata` does not have the same number of predictors as the fitted model.")
    }
    colnames(newdata.frame) <- expected
  }

  extra <- setdiff(colnames(newdata.frame), expected)
  if (length(extra) > 0L) {
    stop(
      sprintf(
        "`newdata` has incompatible predictor columns: %s",
        paste(extra, collapse = ", ")
      )
    )
  }

  missing_cols <- setdiff(expected, colnames(newdata.frame))
  if (length(missing_cols) > 0L) {
    zeros <- as.data.frame(matrix(0, nrow = nrow(newdata.frame), ncol = length(missing_cols)))
    colnames(zeros) <- missing_cols
    rownames(zeros) <- rownames(newdata.frame)
    newdata.frame <- cbind(newdata.frame, zeros)
  }

  newdata.frame[, expected, drop = FALSE]
}

.scores_plsRbetamodel <- function(object, newdata.frame, comps, methodNA, verbose) {
  if (any(apply(is.na(newdata.frame), MARGIN = 1, FUN = all))) {
    stop("One of the rows of newdata is completely filled with missing data")
  }

  centers <- attr(object$ExpliX, "scaled:center")
  scales <- attr(object$ExpliX, "scaled:scale")
  newdata.scaled <- sweep(sweep(as.matrix(newdata.frame), 2, centers, FUN = "-"), 2, scales, FUN = "/")
  newdataNA <- !is.na(newdata.scaled)
  newdata.scaledwotNA <- newdata.scaled
  newdata.scaledwotNA[!newdataNA] <- 0

  tt <- matrix(0, nrow = nrow(newdata.scaledwotNA), ncol = object$computed_nt)

  if (methodNA == "adaptative") {
    for (ii in seq_len(nrow(newdata.scaledwotNA))) {
      if (all(newdataNA[ii, ])) {
        tt[ii, seq_len(comps)] <- newdata.scaledwotNA[ii, ] %*% object$wwetoile[, seq_len(comps), drop = FALSE]
      } else {
        if (verbose) {
          cat("Missing value in row ", ii, ".\n", sep = "")
        }
        tt[ii, seq_len(comps)] <- t(
          solve(
            t(object$pp[newdataNA[ii, ], seq_len(comps), drop = FALSE]) %*%
              object$pp[newdataNA[ii, ], seq_len(comps), drop = FALSE]
          ) %*%
            t(object$pp[newdataNA[ii, ], seq_len(comps), drop = FALSE]) %*%
            (newdata.scaledwotNA[ii, ])[newdataNA[ii, ]]
        )
      }
    }
  } else {
    if (verbose) {
      cat("Prediction as if missing values in every row.\n")
    }
    for (ii in seq_len(nrow(newdata.scaledwotNA))) {
      tt[ii, seq_len(comps)] <- t(
        solve(
          t(object$pp[newdataNA[ii, ], seq_len(comps), drop = FALSE]) %*%
            object$pp[newdataNA[ii, ], seq_len(comps), drop = FALSE]
        ) %*%
          t(object$pp[newdataNA[ii, ], seq_len(comps), drop = FALSE]) %*%
          (newdata.scaledwotNA[ii, ])[newdataNA[ii, ]]
      )
    }
  }

  colnames(tt) <- paste("Comp_", seq_len(object$computed_nt), sep = "")
  rownames(tt) <- rownames(newdata.frame)
  tt
}

.training_scores_plsRbetamodel <- function(object, comps) {
  tt <- matrix(0, nrow = nrow(object$tt), ncol = object$computed_nt)
  tt[, seq_len(comps)] <- object$tt[, seq_len(comps), drop = FALSE]
  colnames(tt) <- colnames(object$tt)
  rownames(tt) <- rownames(object$tt)
  tt
}

.predict_pls_response_plsRbetamodel <- function(object, tt, comps) {
  drop(
    attr(object$RepY, "scaled:center") +
      attr(object$RepY, "scaled:scale") *
      tt[, seq_len(comps), drop = FALSE] %*%
      object$CoeffC[seq_len(comps)]
  )
}

.predict_glm_or_beta_subset_plsRbetamodel <- function(object, tt, comps) {
  coeffs <- object$CoeffCFull[seq_len(comps + 1L), comps]
  eta <- drop(coeffs[1L] + tt[, seq_len(comps), drop = FALSE] %*% coeffs[-1L])

  if (inherits(object$FinalModel, "betareg")) {
    return(object$FinalModel$link$mean$linkinv(eta))
  }

  if (inherits(object$FinalModel, "glm")) {
    return(object$FinalModel$family$linkinv(eta))
  }

  stop("Prediction with fewer components is not implemented for this fitted model.")
}

.final_model_term_plsRbetamodel <- function(object) {
  term_labels <- attr(terms(object$FinalModel), "term.labels")
  if (length(term_labels) != 1L) {
    return(NULL)
  }
  term_labels
}
