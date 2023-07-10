#' Generate synthetic data
#'
#' [runm()] generates synthetic data for use in examples of models with
#' unmeasured confounders (can take up to two unmeasured confounders). Defaults to
#' one unmeasured confounder.
#'
#' @param n Number of observations
#' @param response `"norm"`, `"bin"`, `"pois"`, or `"gam"`
#' @param confounder1 `"norm"` or `"bin"`
#' @param confounder2 `"norm"` or `"bin"`
#' @param type `"int"` or `"ext"` (not implemented)
#' @param missing_prop Proportion of missing values. Only used when `type = "int"`
#' @return A `tibble`
#' @export
#' @name runm
#' @details
#'   \deqn{
#'   x  \sim \text{Ber}(\text{logit}(\gamma_0 + \gamma_z z + \zeta_{u_1} u_1 + \zeta_{u_2} u_2))\\
#'   y  \sim N(\beta_0 + \beta_{x} x + \beta_z z + \lambda_{u_1} u_1 + \lambda_{u_2} u_2, \sigma^2_y = 1),}
#'
#'   Where \eqn{(\beta_0 = -1, \beta_x = 0.75, \beta_z = 0.75, \lambda_{u_1} = 0.75, \lambda_{u_2} = 0.75,
#'    \gamma_0 = -1, \gamma_{z} = .4, \zeta_{u_1} = 0.75, \zeta_{u_2} = 0.75).}
#'
#' @examples
#'
#' # internally validated data
#' runm(20)
#' runm(20, missing_prop = .50)
#'
#' runm(20, response = "bin")
#' runm(20, response = "pois")
#' runm(20, response = "gam")
#' runm(20, response = "norm", confounder1 = "bin")
#' runm(20, response =  "bin", confounder1 = "bin")
#' runm(20, response = "pois", confounder1 = "bin", missing_prop = .7)
#' runm(20, response = "gam", confounder1 = "bin", missing_prop = .7)
#'
#' runm(20, response = "norm", confounder1 = "bin", confounder2 = "norm")
#' runm(20, response =  "bin", confounder1 = "bin", confounder2 = "norm")
#' runm(20, response = "pois", confounder1 = "bin", confounder2 = "bin", missing_prop = .7)
#' runm(20, response = "gam", confounder1 = "bin", confounder2 = "bin", missing_prop = .7)
#'
#' runm(c(15,5), type = "ext")
#' runm(c(15,5), response = "norm", confounder1 = "norm", type = "ext")
#' runm(c(15,5), response = "norm", confounder1 = "bin", confounder2 = "norm", type = "ext")


runm <- function(n, response = "norm",
                 confounder1 = "norm", confounder2 = NA,
                 type = "int", missing_prop = .80) {

  if (type == "int") {
    runm_int(n, response, confounder1, confounder2, missing_prop)
  } else  if (type == "ext") {
    if (length(n) == 2) {
      n_main <- n[1]
      n_external <- n[2]
    } else if (length(n) == 1) {
      if (n %% 2 == 0) {
        n_main <- n_external <- n/2
      } else {
        n_main <- floor(n/2) + 1
        n_external <- floor(n/2)
      }
    } else {
      stop(" externally validated data, `n` should be c(`n_main`, `n_external`).", call. = TRUE)
    }
    runm_ext(n_main, n_external, response, confounder1, confounder2)
  } else {
    stop('`type` should be `"int"` or `"ext"`.', call. = FALSE)
  }

}


runm_full <- function(n, response = "norm",
                      confounder1 = "norm", confounder2 = NA, treatment = TRUE) {

  df <- data.frame("x" = rnorm(n))

  rconf1 <- switch(confounder1,
                   "norm" = function(n) rnorm(n, 1, 1),
                   "bin" = function(n) rbinom(n, 1, .7)
  )

  rconf2 <- switch(confounder2,
                   "norm" = function(n) rnorm(n, 1, 1),
                   "bin" = function(n) rbinom(n, 1, .7)
  )
  #make confounders
  df$u1 <- rconf1(n)
  if(!is.na(confounder2)) df$u2 <- rconf2(n)

  #make treatment
  if(treatment == TRUE) {
    W <- model.matrix(~ ., data = df)
    if(is.na(confounder2)) {ga <- c("ga_1" = -1, "ga_x" = .4, "ga_u1" = .75)
    } else ga <- c("ga_1" = -1, "ga_x" = .4, "ga_u1" = .75, "ga_u2" = .75)
    p_ga <- length(ga) # = # non-confounder params in treatment model
    df$trt <- rbinom(n, 1, binomial()$linkinv(W %*% ga))
  } else df$trt <- 0 # For the case of external validation data

  #make response
  X <- model.matrix(~ ., data = df)
  be1 <- c("be_1" = -1, "be_x" = .75)
  be2 <- c("be_trt" = .75)
  p_be <- length(be1) + length(be2) # = # non-confounder params in response model
  if(is.na(confounder2)) {la <- c("la_u1" = .75)
  } else la <- c("la_u1" = .75, "la_u2" = .75)
  p_la <- length(la) # = # confounder params in response model
  (th <- c(be1, la, be2))


  si_y <- if (response == "norm") c("si_y" = 1) else NULL
  if (response == "pois") df$t <- runif(n, 1, 10)
  al_y <- if (response == "gam") c("al_y" = 2) else NULL

  invlink_y <- switch(response,
                      "norm" = function(x) gaussian()$linkinv(as.numeric(x)),
                      "bin" = function(x) binomial()$linkinv(as.numeric(x)),
                      "pois" = function(x) poisson()$linkinv(as.numeric(x)),
                      "gam" = function(x) Gamma(link = "log")$linkinv(as.numeric(x))
  )
  rresp <- switch(response,
                  "norm" = function(n) rnorm(n, invlink_y(X %*% th), si_y),
                  "bin" = function(n) rbinom(n, 1, invlink_y(X %*% th)),
                  "pois" = function(n) rpois(n, lambda = df$t * invlink_y(X %*% th)),
                  "gam" = function(n) rgamma(n, shape = al_y, rate = al_y / invlink_y(X %*% th))
  )
  df$y <- rresp(n)


  # add metadata
  if(treatment == TRUE) {
    params <- drop_nulls(c(ga, be1, be2, la, si_y, al_y))
  } else params <- drop_nulls(c(be1, be2, la, si_y, al_y))
  names(params) <- expand_labels(names(params))
  params

  df <- tibble::as_tibble(df)
  attr(df, "params") <- params


  # return
  df
}








runm_int <- function(n, response = "norm",
                     confounder1 = "norm", confounder2 = NA, missing_prop) {

  df <- runm_full(n, response, confounder1, confounder2)

  # remove unmeasured confounder(s)
  # ndcs_to_remove <- sample(1:n, ceiling(n*missing_prop))
  ndcs_to_remove <- (n-ceiling(n*missing_prop)+1):n
  df$u1[ndcs_to_remove] <- NA
  if("u2" %in% colnames(df)) df$u2[ndcs_to_remove] <- NA else NULL

  df

}





runm_ext <- function(n_main, n_external, response, confounder1, confounder2) {

  df_main <- runm_full(n_main, response, confounder1, confounder2)
  df_main$u1 <- NA
  df_main$u2 <- NA

  df_external <- runm_full(n_external, response, confounder1, confounder2,
                           treatment = FALSE)


  rbind(df_main[names(df_external)], df_external)
}




runm_full <- function(n, response = "norm",
                      confounder1 = "norm", confounder2 = NA, treatment = TRUE) {

  df <- data.frame("x" = rnorm(n))

  rconf1 <- switch(confounder1,
                   "norm" = function(n) rnorm(n, 1, 1),
                   "bin" = function(n) rbinom(n, 1, .7)
  )

  rconf2 <- switch(confounder2,
                   "norm" = function(n) rnorm(n, 1, 1),
                   "bin" = function(n) rbinom(n, 1, .7)
  )
  #make confounders
  df$u1 <- rconf1(n)
  if(!is.na(confounder2)) df$u2 <- rconf2(n)

  #make treatment
  if(treatment == TRUE) {
    W <- model.matrix(~ ., data = df)
    if(is.na(confounder2)) {ga <- c("ga_1" = -1, "ga_x" = .4, "ga_u1" = .75)
    } else ga <- c("ga_1" = -1, "ga_x" = .4, "ga_u1" = .75, "ga_u2" = .75)
    p_ga <- length(ga) # = # non-confounder params in treatment model
    df$trt <- rbinom(n, 1, binomial()$linkinv(W %*% ga))
  } else df$trt <- 0 # For the case of external validation data

  #make response
  X <- model.matrix(~ ., data = df)
  be1 <- c("be_1" = -1, "be_x" = .75)
  be2 <- c("be_trt" = .75)
  p_be <- length(be1) + length(be2) # = # non-confounder params in response model
  if(is.na(confounder2)) {la <- c("la_u1" = .75)
  } else la <- c("la_u1" = .75, "la_u2" = .75)
  p_la <- length(la) # = # confounder params in response model
  (th <- c(be1, la, be2))


  si_y <- if (response == "norm") c("si_y" = 1) else NULL
  if (response == "pois") df$t <- runif(n, 1, 10)
  al_y <- if (response == "gam") c("al_y" = 2) else NULL

  invlink_y <- switch(response,
                      "norm" = function(x) gaussian()$linkinv(as.numeric(x)),
                      "bin" = function(x) binomial()$linkinv(as.numeric(x)),
                      "pois" = function(x) poisson()$linkinv(as.numeric(x)),
                      "gam" = function(x) Gamma(link = "log")$linkinv(as.numeric(x))
  )
  rresp <- switch(response,
                  "norm" = function(n) rnorm(n, invlink_y(X %*% th), si_y),
                  "bin" = function(n) rbinom(n, 1, invlink_y(X %*% th)),
                  "pois" = function(n) rpois(n, lambda = df$t * invlink_y(X %*% th)),
                  "gam" = function(n) rgamma(n, shape = al_y, rate = al_y / invlink_y(X %*% th))
  )
  df$y <- rresp(n)


  # add metadata
  if(treatment == TRUE) {
    params <- drop_nulls(c(ga, be1, be2, la, si_y, al_y))
  } else params <- drop_nulls(c(be1, be2, la, si_y, al_y))
  names(params) <- expand_labels(names(params))
  params

  df <- tibble::as_tibble(df)
  attr(df, "params") <- params


  # return
  df
}


runm_extended <- function(n,
                          covariate_fam_list,
                          covariate_param_list,
                          unmeasured_fam_list,
                          unmeasured_param_list,
                          treatment_model_coefs,
                          response_model_coefs,
                          response,
                          response_param,
                          type = "int", missing_prop = .80) {
  if (type == "int") {
    runm_int_extended(n,
                      covariate_fam_list,
                      covariate_param_list,
                      unmeasured_fam_list,
                      unmeasured_param_list,
                      treatment_model_coefs,
                      response_model_coefs,
                      response,
                      response_param, missing_prop)
  } else  if (type == "ext") {
    if (length(n) == 2) {
      n_main <- n[1]
      n_external <- n[2]
    } else if (length(n) == 1) {
      if (n %% 2 == 0) {
        n_main <- n_external <- n/2
      } else {
        n_main <- floor(n/2) + 1
        n_external <- floor(n/2)
      }
    } else {
      stop(" externally validated data, `n` should be c(`n_main`, `n_external`).", call. = TRUE)
    }
    runm_ext_extended(n_main, n_external,
                      covariate_fam_list,
                      covariate_param_list,
                      unmeasured_fam_list,
                      unmeasured_param_list,
                      treatment_model_coefs,
                      response_model_coefs,
                      response,
                      response_param)
  } else {
    stop('`type` should be `"int"` or `"ext"`.', call. = FALSE)
  }

}





runm_full_extended <- function(n,
                               covariate_fam_list,
                               covariate_param_list,
                               unmeasured_fam_list,
                               unmeasured_param_list,
                               treatment_model_coefs,
                               response_model_coefs,
                               response,
                               response_param,
                               treatment = TRUE) {

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

  trt_param_vec <- treatment_model_coefs # c("int" = -1, "z1" = .4, "z2" = .5, "z3" = .4, "u1" = .75, "u2" = .75)
  (length(trt_param_vec) - 1) == ncol(df_combined)
  names(df_combined) == names(trt_param_vec)[-1]
  vec_names <- names(trt_param_vec)
  names(trt_param_vec) <- sapply(vec_names, function(.) {
    if(grepl("z|nt", .)) {
      glue::glue("ga_{.}")
    } else glue::glue("ze_{.}")
  })
  #make treatment
  if(treatment == TRUE) {
    W <- model.matrix(~ ., data = df_combined)
    ga <- trt_param_vec
    p_ga <- length(ga) # = # non-confounder params in treatment model
    df_combined$trt <- rbinom(n, 1, binomial()$linkinv(W %*% ga))
  } else df_combined$trt <- 0 # For the case of external validation data

  #make response
  X <- model.matrix(~ ., data = df_combined)

  (length(response_model_coefs) - 1) == ncol(df_combined)
  names(df_combined) == names(response_model_coefs)[-1]
  vec_names <- names(response_model_coefs)
  names(response_model_coefs) <- sapply(vec_names, function(.) {
    if(grepl("z|nt|trt", .)) {
      glue::glue("be_{.}")
    } else glue::glue("la_{.}")
  })
  be <- response_model_coefs # c("int" = -1, "z1" = .4, "z2" = .5, "z3" = .4, "u1" = .75, "u2" = .75, "trt" = .75)

  # c("be_1" = -1, "be_x" = .75, "la_u1" = .75, "la_u2" = .75, "be_trt" = .75)
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
  rresp <- switch(response,
                  "norm" = function(n) rnorm(n, invlink_y(X %*% th), si_y),
                  "bin" = function(n) rbinom(n, 1, invlink_y(X %*% th)),
                  "pois" = function(n) rpois(n, lambda = df$t * invlink_y(X %*% th)),
                  "gam" = function(n) rgamma(n, shape = al_y, rate = al_y / invlink_y(X %*% th))
  )
  df_combined$y <- rresp(n)


  # add metadata
  if(treatment == TRUE) {
    params <- drop_nulls(c(ga, be, si_y, al_y))
  } else params <- drop_nulls(c(be, si_y, al_y))
  names(params) <- expand_labels(names(params))
  params

  df_combined <- tibble::as_tibble(df_combined)
  attr(df_combined, "params") <- params


  # return
  df_combined

}








runm_int_extended <- function(n,
                              covariate_fam_list,
                              covariate_param_list,
                              unmeasured_fam_list,
                              unmeasured_param_list,
                              treatment_model_coefs,
                              response_model_coefs,
                              response,
                              response_param, missing_prop) {

  df <- runm_full_extended(n,
                           covariate_fam_list,
                           covariate_param_list,
                           unmeasured_fam_list,
                           unmeasured_param_list,
                           treatment_model_coefs,
                           response_model_coefs,
                           response,
                           response_param)

  # remove unmeasured confounder(s)
  # ndcs_to_remove <- sample(1:n, ceiling(n*missing_prop))
  ndcs_to_remove <- (n-ceiling(n*missing_prop)+1):n
  df$u1[ndcs_to_remove] <- NA
  if("u2" %in% colnames(df)) df$u2[ndcs_to_remove] <- NA else NULL

  df

}





runm_ext_extended <- function(n_main, n_external,
                              covariate_fam_list,
                              covariate_param_list,
                              unmeasured_fam_list,
                              unmeasured_param_list,
                              treatment_model_coefs,
                              response_model_coefs,
                              response,
                              response_param) {

  df_main <- runm_full_extended(n_main,
                                covariate_fam_list,
                                covariate_param_list,
                                unmeasured_fam_list,
                                unmeasured_param_list,
                                treatment_model_coefs,
                                response_model_coefs,
                                response,
                                response_param)
  df_main$u1 <- NA
  df_main$u2 <- NA

  df_external <- runm_full_extended(n_external,
                                    covariate_fam_list,
                                    covariate_param_list,
                                    unmeasured_fam_list,
                                    unmeasured_param_list,
                                    treatment_model_coefs,
                                    response_model_coefs,
                                    response,
                                    response_param,
                                    treatment = FALSE)


  rbind(df_main[names(df_external)], df_external)

}

