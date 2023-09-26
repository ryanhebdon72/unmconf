#' Generate synthetic data
#'
#' [unm_summary()] generates synthetic data for use in examples of models with
#' unmeasured confounders (can take up to two unmeasured confounders). Defaults to
#' one unmeasured confounder.
#'
#'
#' @param mod Output from [unm_glm()].
#' @param data The data `mod` was generated with.
#' @param quantiles A numeric vector of quantiles.
#' @return A `tibble`
#' @export
#' @name unm_summary
#' @export
#'
unm_summary <- function(mod, data, quantiles = c(.025, .975)) {

  param <- NULL; rm(param)
  true_value <- NULL; rm(true_value)
  sd <- NULL; rm(sd)

  summary <- summary(mod, quantiles = quantiles)

  stats <- summary$statistics |>
    as.data.frame() |>
    tibble::as_tibble(rownames = "param") |>
    janitor::clean_names() |>
    dplyr::select(param, mean, sd)

  cis <- summary$quantiles |>
    as.data.frame() |>
    tibble::as_tibble(rownames = "param")

  out <- stats |>
    dplyr::left_join(cis, by = "param")

  if (!missing(data)) {
    out <- out |>
      dplyr::left_join(
        tibble::enframe(attr(data, "param"), name = "param", value = "true_value"),
        by = "param"
      ) |>
      dplyr::relocate(true_value, .after = param)
  }

  out

}

