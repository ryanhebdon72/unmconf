#' Generate synthetic data
#'
#' [runm()] generates synthetic data for use of modeling with unmeasured
#' confounders. Defaults to the case of one unmeasured confounder present and
#' fixed parameter values. Can be customized. Currently set up to have at most
#' two unmeasured confounders to pair with [unm_glm()].
#'
#' @param n Number of observations. When `type = "int"`, `n` is a vector of
#'   length 1. When `type = "ext"`, `n` can either be a vector of length 1 or 2.
#'   For the case when `n` is of length 2, `n = (n_main, n_external)`, where
#'   `n_main` corresponds to the main study sample size and `n_external`
#'   corresponds to the external validation sample size. For the case when `n`
#'   is of length 1, `n` will be split evenly between main study and external
#'   validation observations, with the main study getting the additional
#'   observation when `n` is odd.
#' @param type Type of validation source. Can be `"int"` for internal validation
#'   or `"ext"` for external validation. Defaults to `"int"`.
#' @param missing_prop Proportion of missing values for internal validation
#'   scenario (i.e., when `type = "int"`).
#' @param response `"norm"`, `"bin"`, `"pois"`, or `"gam"`. Defaults to `"bin"`.
#' @param response_param Nuisance parameters for response type. For `"norm"`,
#'   the default standard deviation is 1. For `"gam"`, the default shape
#'   parameter is 2. For `"pois"`, an offset variable is added to the dataset
#'   that is uniformly distributed from 1 to 10.
#' @param response_model_coefs A named vector of coefficients to generate data
#'   from the response model. This must include an intercept (`"int" = `), a
#'   coefficient for each covariate specified, a coefficient for each unmeasured
#'   confounder, and a treatment coefficient (`"x" = `). The coefficients for
#'   the covariates and treatment will be denoted with `"beta[.]"` and the
#'   unmeasured confounders with `"lambda[.]"`.
#' @param treatment_model_coefs A named vector of coefficients to generate data
#'   from the treatment model. This must include an intercept (`"int" = `), a
#'   coefficient for each covariate specified, and a coefficient for each
#'   unmeasured confounder. The coefficients for the covariates and unmeasured
#'   confounders will be denoted with `"eta[.]"`.
#' @param covariate_fam_list A list of either `"norm"` or `"bin"`, where the
#'   length of the list matches the number of covariates in the model.
#' @param covariate_param_list A list of parameters for the respective
#'   distributions in `covariate_fam_list`, where the length of the list matches
#'   the length of `covariate_fam_list`.
#' @param unmeasured_fam_list A list of either `"norm"` or `"bin"`, where the
#'   length of the list matches the number of unmeasured confounders in the
#'   model. This can be at most a length of 2 to pair with [unm_glm()].
#' @param unmeasured_param_list A list of parameters for the respective
#'   distributions in `unmeasured_fam_list`, where the length of the list
#'   matches the length of `unmeasured_fam_list`.
#'
#' @return A `tibble`
#'
#' @name runm
#' @examples
#'
#' runm(100)
#' runm(n = 100, type = "int", missing_prop = .75)
#' runm(n = 100, type = "int", missing_prop = .75) |> attr("params")
#' runm(100, type = "int", response = "norm")
#' runm(100, type = "int", response = "norm") |> attr("params")
#' runm(100, type = "int", response = "norm", response_param = 3) |> attr("params")
#' runm(100, type = "int", response = "gam")
#' runm(100, type = "int", response = "gam", response_param = 5) |> attr("params")
#' runm(100, type = "int", missing_prop = .5, response = "pois")
#'
#' runm(n = 100, type = "ext")
#' runm(n = 100, type = "ext") |> attr("params")
#' runm(n = c(10, 10), type = "ext")
#' runm(100, type = "ext", response = "norm")
#' runm(100, type = "int", response = "norm", response_param = 3) |> attr("params")
#' runm(100, type = "ext", response = "gam")
#' runm(100, type = "ext", response = "pois")
#'
#' runm(
#'   n = 100,
#'   type = "int",
#'   missing_prop = .80,
#'   response = "norm",
#'   response_param = c("si_y" = 2),
#'   response_model_coefs = c("int" = -1, "z" = .4,
#'                            "u1" = .75, "u2" = .75, "x" = .75),
#'   treatment_model_coefs = c("int" = -1, "z" = .4,
#'                             "u1" = .75, "u2" = .75),
#'   covariate_fam_list = list("norm"),
#'   covariate_param_list = list(c(mean = 0, sd = 1)),
#'   unmeasured_fam_list = list("norm", "bin"),
#'   unmeasured_param_list = list(c(mean = 0, sd = 1), c(.3))
#' )
#'
#' runm(
#'   n = c(20, 30),
#'   type = "ext",
#'   response = "norm",
#'   response_param = c("si_y" = 2),
#'   response_model_coefs = c("int" = -1, "z1" = .4, "z2" = .5, "z3" = .4,
#'                            "u1" = .75, "u2" = .75, "x" = .75),
#'   treatment_model_coefs = c("int" = -1, "z1" = .4, "z2" = .5, "z3" = .4,
#'                             "u1" = .75, "u2" = .75),
#'   covariate_fam_list = list("norm", "bin", "norm"),
#'   covariate_param_list = list(c(mean = 0, sd = 1), c(.3), c(0, 2)),
#'   unmeasured_fam_list = list("norm", "bin"),
#'   unmeasured_param_list = list(c(mean = 0, sd = 1), c(.3))
#' )


#' @rdname runm
#' @export
runm <- function(
    n,
    type = "int",
    missing_prop = .80,
    response = "bin", response_param = NULL,
    response_model_coefs = c("int" = -1, "z1" = .5, "z2" = .5, "z3" = .5,
                             "u1" = .5, "x" = .5),
    treatment_model_coefs = c("int" = -1, "z1" = .5, "z2" = .5, "z3" = .5,
                              "u1" = .5),
    covariate_fam_list = list("norm", "bin", "norm"),
    covariate_param_list = list(c(mean = 0, sd = 1), prob = .3, c(0, 2)),
    unmeasured_fam_list = list("norm"),
    unmeasured_param_list = list(c(mean = 0, sd = 1))
) {
  if (type == "int") {

    runm_int(
      n,
      missing_prop,
      response, response_param,
      response_model_coefs,
      treatment_model_coefs,
      covariate_fam_list, covariate_param_list,
      unmeasured_fam_list, unmeasured_param_list
    )

  } else  if (type == "ext") {

    if (length(n) == 2) {
      n_main <- n[1]
      n_external <- n[2]
    } else if (length(n) == 1) {
      if (n %% 2 == 0) {
        n_main <- n_external <- n / 2
      } else {
        n_main <- floor(n / 2) + 1
        n_external <- floor(n / 2)
      }

    } else {
      stop(" externally validated data, `n` should be c(`n_main`, `n_external`).", call. = TRUE)
    }

    runm_ext(
      n_main, n_external,
      response, response_param,
      response_model_coefs,
      treatment_model_coefs,
      covariate_fam_list, covariate_param_list,
      unmeasured_fam_list, unmeasured_param_list
    )
  } else {
    stop('`type` should be `"int"` or `"ext"`.', call. = FALSE)
  }

}




#' @noRd
runm_full <- function(
    n,
    response, response_param,
    response_model_coefs,
    treatment_model_coefs,
    covariate_fam_list, covariate_param_list,
    unmeasured_fam_list, unmeasured_param_list,
    treatment = TRUE
) {

  if (length(covariate_fam_list) != length(covariate_param_list)) {
    stop("`covariate_fam_list` and `covariate_param_list` must have same length.")
  }

  if (length(unmeasured_fam_list) != length(unmeasured_param_list)) {
    stop("`unmeasured_fam_list` and `unmeasured_param_list` must have same length.")
  }

  if (any(treatment_model_coefs[-1] %in% "1.|Int|int|Intercept")) {
    stop("First element in `treatment_model_coefs` must be coefficient for
         intercept, named '1', 'Int', 'int', or 'Intercept'.")
  }

  if (any(names(response_model_coefs[-1]) %in% "1.|Int|int|Intercept")) {
    stop("First element in `response_model_coefs` must be coefficient for
         intercept, named '1', 'Int', 'int', or 'Intercept'.")
  }

  if (names(response_model_coefs[length(response_model_coefs)]) != "x") {
    stop("Last element in `response_model_coefs` must be coefficient for
         treatment, named 'x'.")
  }

  if (any(names(response_model_coefs[-length(response_model_coefs)]) !=
      names(treatment_model_coefs))) {
    stop("`response_model_coefs` and `treatment_model_coefs` must have same names
         attribute and be in the same order (excluding `x` in `response_model_coefs`).")
  }

  a <- covariate_fam_list
  b <- covariate_param_list
  (covariate_df <- as.data.frame(
    Map(function(dist, param) {
      if(dist == "norm") {
        rnorm(n, param[1], param[2])
      } else if(dist == "bin") {
        rbinom(n, 1, param)
      }
    }, a, b),
    col.names = sapply(1:length(a), function(.) {
      glue::glue("z", .)
    })
  )
  )

  a <- unmeasured_fam_list
  b <- unmeasured_param_list
  (unm_conf_df <- as.data.frame(
    Map(function(dist, param) {
      if(dist == "norm") {
        rnorm(n, param[1], param[2])
      } else if(dist == "bin") {
        rbinom(n, 1, param)
      }
    }, a, b),
    col.names = sapply(1:length(a), function(.) {
      glue::glue("u", .)
    })
  )
  )

  df_combined <- cbind(covariate_df, unm_conf_df)

  z_param_vec <- treatment_model_coefs
  (length(z_param_vec) - 1) == ncol(df_combined)
  names(df_combined) <- names(z_param_vec)[-which(grepl("1.|Int|int|Intercept",
                                                        names(z_param_vec)))]
  vec_names <- names(z_param_vec)
  names(z_param_vec) <- sapply(vec_names, function(.) {
    glue::glue("et_{.}")
  })
  #make treatment
  if(treatment == TRUE) {
    W <- model.matrix(~ ., data = df_combined)
    et <- z_param_vec
    p_et <- length(et) # = # non-confounder params in treatment model
    df_combined$x <- rbinom(n, 1, binomial()$linkinv(W %*% et))
  } else df_combined$x <- 0 # For the case of external validation data

  #make response
  X <- model.matrix(~ ., data = df_combined)


  (length(response_model_coefs) - 1) == ncol(df_combined)
  names(df_combined) <- names(response_model_coefs)[-1]
  vec_names <- names(response_model_coefs)
  names(response_model_coefs) <- sapply(vec_names, function(.) {
    if(grepl("z|nt|x", .)) {
      glue::glue("be_{.}")
    } else glue::glue("la_{.}")
  })
  be <- response_model_coefs
  p_be <- length(be) # = # non-confounder params in response model
  # p_la <- length(la) # = # confounder params in response model
  (th <- be)


  (if (is.null(response_param)) {
    nuis_param <- switch(response,
                         "norm" = c("si_y" = 1),
                         "gam" = c("al_y" = 2))
  } else {
    nuis_param <- switch(response,
                         "norm" = c(response_param),
                         "gam" = c(response_param)
    )
  })

  names(nuis_param) <- switch(response,
                              "norm" = "si_y",
                              "gam" = "al_y"
  )
  if (response == "pois") df_combined$t <- runif(n, 1, 10)

  invlink_y <- switch(response,
                      "norm" = function(x) gaussian()$linkinv(as.numeric(x)),
                      "bin" = function(x) binomial()$linkinv(as.numeric(x)),
                      "pois" = function(x) poisson()$linkinv(as.numeric(x)),
                      "gam" = function(x) Gamma(link = "log")$linkinv(as.numeric(x)),
                      stop()
  )

  resp <- switch(response,
                 "norm" = function(n) rnorm(n, invlink_y(X %*% th), nuis_param),
                 "bin" = function(n) rbinom(n, 1, invlink_y(X %*% th)),
                 "pois" = function(n) rpois(n, lambda = df_combined$t * invlink_y(X %*% th)),
                 "gam" = function(n) rgamma(n, shape = nuis_param, rate = nuis_param / invlink_y(X %*% th))
  )

  df_combined$y <- resp(n)


  # add metadata
  if(treatment == TRUE) {
    params <- drop_nulls(c(et, be, nuis_param))
  } else params <- drop_nulls(c(be, nuis_param))
  names(params) <- expand_labels(names(params))
  params

  df_combined <- tibble::as_tibble(df_combined)
  attr(df_combined, "params") <- params


  # return
  df_combined

}







#' @noRd
runm_int <- function(
    n,
    missing_prop,
    response, response_param,
    response_model_coefs,
    treatment_model_coefs,
    covariate_fam_list, covariate_param_list,
    unmeasured_fam_list, unmeasured_param_list
) {

  df <- runm_full(
    n,
    response, response_param,
    response_model_coefs,
    treatment_model_coefs,
    covariate_fam_list, covariate_param_list,
    unmeasured_fam_list, unmeasured_param_list
  )

  # remove unmeasured confounder(s)
  # ndcs_to_remove <- sample(1:n, ceiling(n*missing_prop))
  ndcs_to_remove <- (n-ceiling(n*missing_prop)+1):n
  df$u1[ndcs_to_remove] <- NA
  if("u2" %in% colnames(df)) df$u2[ndcs_to_remove] <- NA else NULL

  df

}




#' @noRd
runm_ext <- function(
    n_main, n_external,
    response, response_param,
    response_model_coefs,
    treatment_model_coefs,
    covariate_fam_list, covariate_param_list,
    unmeasured_fam_list, unmeasured_param_list
) {

  df_main <- runm_full(
    n_main,
    response, response_param,
    response_model_coefs,
    treatment_model_coefs,
    covariate_fam_list, covariate_param_list,
    unmeasured_fam_list, unmeasured_param_list
  )

  df_main$u1 <- NA
  df_main$u2 <- NA

  df_external <- runm_full(
    n_external,
    response, response_param,
    response_model_coefs,
    treatment_model_coefs,
    covariate_fam_list, covariate_param_list,
    unmeasured_fam_list, unmeasured_param_list,
    treatment = FALSE
  )


  rbind(df_main[names(df_external)], df_external)

}

