#' Generate synthetic data
#'
#' [runm()] generates customized synthetic data for use in examples of models
#' with unmeasured confounders. Currently set up to have at most 2 unmeasured
#' confounders to pair with [unm_glm()].
#'
#' @param n Number of observations
#' @param response `"norm"`, `"bin"`, `"pois"`, or `"gam"`
#' @param response_param Nuisance parameters for response type.
#' @param response_model_coefs A named vector of coefficients to generate data
#'   from the response model. This must include an intercept (`"int" = `), a
#'   coefficient for each covariate specified, a coefficient for each unmeasured
#'   confounder, and a treatment coefficient (`"x" = `).
#' @param treatment_model_coefs A named vector of coefficients to generate data
#'   from the treatment model. This must include an intercept (`"int" =`), a
#'   coefficient for each covariate specified, and a coefficient for each
#'   unmeasured confounder.
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
#' @param type `"int"` or `"ext"`
#' @param missing_prop Proportion of missing values. Only used when `type =
#'   "int"`
#' @param treatment Input documentation here.
#' @param n_main Input documentation here.
#' @param n_external Input documentation here.
#' @return A `tibble`
#' @name runm
#' @examples
#'
#' runm(
#'   n = 100,
#'   type = "int",
#'   missing_prop = .80,
#'   response = "norm",
#'   response_param = c("si_y" = 1),
#'   response_model_coefs = c("int" = -1, "z1" = .4, "z2" = .5, "z3" = .4,
#'                            "u1" = .75, "u2" = .75, "x" = .75),
#'   treatment_model_coefs = c("int" = -1, "z1" = .4, "z2" = .5, "z3" = .4,
#'                             "u1" = .75, "u2" = .75),
#'   covariate_fam_list = list("norm", "bin", "norm"),
#'   covariate_param_list = list(c(mean = 0, sd = 1), c(.3), c(0, 2)),
#'   unmeasured_fam_list = list("norm", "bin"),
#'   unmeasured_param_list = list(c(mean = 0, sd = 1), c(.3))
#' )
#'
#' runm(
#'   n = 100,
#'   type = "ext",
#'   response = "norm",
#'   response_param = c("si_y" = 1),
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
    response, response_param = NULL,
    response_model_coefs,
    treatment_model_coefs,
    covariate_fam_list, covariate_param_list,
    unmeasured_fam_list, unmeasured_param_list
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




#' @rdname runm
runm_full <- function(
    n,
    response, response_param,
    response_model_coefs,
    treatment_model_coefs,
    covariate_fam_list, covariate_param_list,
    unmeasured_fam_list, unmeasured_param_list,
    treatment = TRUE
) {

  a <- covariate_fam_list #list("norm", "bin", "norm")
  b <- covariate_param_list #list(c(mean = 0, sd = 1), c(.3), c(10, 2))
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


  a <- unmeasured_fam_list # list("norm", "bin")
  b <- unmeasured_param_list # list(c(mean = 0, sd = 1), c(.3))
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

  x_param_vec <- treatment_model_coefs # c("int" = -1, "z1" = .4, "z2" = .5, "z3" = .4, "u1" = .75, "u2" = .75)
  (length(x_param_vec) - 1) == ncol(df_combined)
  names(df_combined) == names(x_param_vec)[-1]
  vec_names <- names(x_param_vec)
  names(x_param_vec) <- sapply(vec_names, function(.) {
    glue::glue("et_{.}")
    # if(grepl("z|nt", .)) {
    #   glue::glue("ga_{.}")
    # } else glue::glue("ze_{.}")
  })
  #make treatment
  if(treatment == TRUE) {
    W <- model.matrix(~ ., data = df_combined)
    et <- x_param_vec
    p_et <- length(et) # = # non-confounder params in treatment model
    df_combined$x <- rbinom(n, 1, binomial()$linkinv(W %*% et))
  } else df_combined$x <- 0 # For the case of external validation data

  #make response
  X <- model.matrix(~ ., data = df_combined)

  (length(response_model_coefs) - 1) == ncol(df_combined)
  names(df_combined) == names(response_model_coefs)[-1]
  vec_names <- names(response_model_coefs)
  names(response_model_coefs) <- sapply(vec_names, function(.) {
    if(grepl("z|nt|x", .)) {
      glue::glue("be_{.}")
    } else glue::glue("la_{.}")
  })
  be <- response_model_coefs # c("int" = -1, "z1" = .4, "z2" = .5, "z3" = .4, "u1" = .75, "u2" = .75, "x" = .75)

  # c("be_1" = -1, "be_z" = .75, "la_u1" = .75, "la_u2" = .75, "be_x" = .75)
  p_be <- length(be) # = # non-confounder params in response model
  # p_la <- length(la) # = # confounder params in response model
  (th <- be)


  si_y <- if (response == "norm") c("si_y" = response_param) else NULL
  if (response == "pois") df$t <- runif(n, 1, 10)
  al_y <- if (response == "gam") c("al_y" = response_param) else NULL

  invlink_y <- switch(response,
                      "norm" = function(x) gaussian()$linkinv(as.numeric(x)),
                      "bin" = function(x) binomial()$linkinv(as.numeric(x)),
                      "pois" = function(x) poisson()$linkinv(as.numeric(x)),
                      "gam" = function(x) Gamma(link = "log")$linkinv(as.numeric(x))
  )

  resp <- switch(response,
                 "norm" = function(n) rnorm(n, invlink_y(X %*% th), si_y),
                 "bin" = function(n) rbinom(n, 1, invlink_y(X %*% th)),
                 "pois" = function(n) rpois(n, lambda = df$t * invlink_y(X %*% th)),
                 "gam" = function(n) rgamma(n, shape = al_y, rate = al_y / invlink_y(X %*% th))
  )

  df_combined$y <- resp(n)


  # add metadata
  if(treatment == TRUE) {
    params <- drop_nulls(c(et, be, si_y, al_y))
  } else params <- drop_nulls(c(be, si_y, al_y))
  names(params) <- expand_labels(names(params))
  params

  df_combined <- tibble::as_tibble(df_combined)
  attr(df_combined, "params") <- params


  # return
  df_combined

}







#' @rdname runm
runm_int <- function(
    n,
    response, response_param,
    response_model_coefs,
    treatment_model_coefs,
    covariate_fam_list, covariate_param_list,
    unmeasured_fam_list, unmeasured_param_list,
    missing_prop
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




#' @rdname runm
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

