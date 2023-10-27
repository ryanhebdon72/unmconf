#' Fitting Multi-Staged Bayesian Regression Model with Unmeasured Confounders
#'
#' [unm_glm()] fits a multi-staged Bayesian regression model that accounts for
#' unmeasured confounders. Users can input model information into [unm_glm()] in
#' a similar manner as they would for the standard [stats::glm()] function,
#' providing arguments like `formula`, `family`, and `data`. Results are stored
#' as MCMC iterations.
#'
#' @param form1 The formula specification for the response model (stage I)
#' @param form2 The formula specification for the first unmeasured confounder
#'   model (stage II)
#' @param form3 The formula specification for the second unmeasured confounder
#'   model (stage III)
#' @param family1,family2,family3 The family object, communicating the types of
#'   models to be used for response (`form1`) and unmeasured confounder (`form2,
#'   form3`) models. See [stats::family()] for details
#' @param data The dataset containing all variables (this function currently
#'   only supports a single dataset containing internally validated data)
#' @param n.iter `n.iter` argument of [rjags::coda.samples()]
#' @param n.adapt `n.adapt` argument of [rjags::jags.model()]
#' @param thin `thin` argument of [rjags::coda.samples()]
#' @param n.chains `n.chains` argument of [rjags::jags.model()]
#' @param filename File name where to store jags code
#' @param quiet The `quiet` parameter of [rjags::jags.model()]. Defaults to
#'   `TRUE`, but you can change it on a per-session basis with
#'   `options(unm_quiet = FALSE)`.
#' @param progress.bar The `progress.bar` parameter of [rjags::update.jags()].
#'   Defaults to `"none"`, but you can change it on a per-session basis with
#'   `options(unm_progress.bar = "text")`.
#' @param code_only Should only the code be created?
#' @param default_prior The default prior to use on the regression coefficients.
#' @param priors Custom priors to use on regression coefficients, see examples.
#' @param
#'   response_nuisance_priors,confounder1_nuisance_priors,confounder2_nuisance_priors
#'   JAGS code for the nuisance priors on parameters in a JAGS model (see
#'   examples)
#' @param
#'   response_params_to_track,confounder1_params_to_track,confounder2_params_to_track
#'   Additional parameters to track when nuisance parameter priors are used (see
#'   examples)
#' @param ... Additional arguments to pass into [rjags::jags.model()], such as
#'   `inits`
#' @return (Invisibly) The output of [rjags::coda.samples()], an object of class
#'   `mcmc.list`, along with attributes `code` containing the jags code used and
#'   `file` containing the filename of the jags code.
#' @name unm_glm
#' @seealso [runm()], [rjags::dic.samples()]
#' @examples
#'
#' # ~~ One Unmeasured Confounder Examples (II-Stage Model) ~~
#'
#' # normal response, normal confounder model with internally validated data
#' (df <- runm(20, response = "norm"))
#' (unm_mod <- unm_glm(y ~ x + z1 + z2 + z3 + u1,
#'                     u1 ~ x + z1 + z2 + z3,
#'                     family1 = gaussian(),
#'                     family2 = gaussian(), data = df))
#' (unm_mod <- unm_glm(y ~ .,
#'                     u1 ~ . - y,
#'                     family1 = gaussian(),
#'                     family2 = gaussian(), data = df))
#' glm(y ~ x + z1 + z2 + z3, data = df)
#' coef(unm_mod)
#'
#' jags_code(unm_mod)
#' unm_glm(y ~ .,
#'         u1 ~ . - y,
#'         family1 = gaussian(),
#'         family2 = gaussian(), data = df, code_only = TRUE)
#'
#'
#'
#'
#' # a normal-normal model - external validation
#' (df <- runm(c(10, 10), type = "ext", response = "norm"))
#' (unm_mod <- unm_glm(y ~ x + z1 + z2 + z3 + u1,
#'                     u1 ~ x + z1 + z2 + z3,
#'                     family1 = gaussian(),
#'                     family2 = gaussian(), data = df))
#'
#'
#'
#' # setting custom priors
#' unm_glm(y ~ .,
#'         u1 ~ . - y,
#'         family1 = gaussian(),
#'         family2 = gaussian(), data = df, code_only = TRUE
#' )
#' unm_glm(y ~ .,
#'         u1 ~ . - y,
#'         family1 = gaussian(),
#'         family2 = gaussian(), data = df, code_only = FALSE,
#'         priors = c("lambda[u1]" = "dnorm(1, 10)"),
#'         response_nuisance_priors = "tau_{y} <- sigma_{y}^-2; sigma_{y} ~ dunif(0, 100)",
#'         response_params_to_track = "sigma_{y}",
#'         confounder1_nuisance_priors = "tau_{u1} <- sigma_{u1}^-2; sigma_{u1} ~ dunif(0, 100)",
#'         confounder1_params_to_track = "sigma_{u1}"
#' )
#'
#'
#'
#' # turn progress tracking on
#' options("unm_progress.bar" = "text")
#'
#'
#'
#' # more complex functional forms _for non-confounder predictors only_
#' # zero-intercept model
#' unm_glm(y ~ . - 1,
#'         u1 ~ . - y,
#'         family1 = gaussian(),
#'         family2 = gaussian(), data = df)
#' glm(y ~ . - 1, data = df)
#'
#' # polynomial model
#' unm_glm(y ~ x + poly(z1, 2) + u1,
#'         u1 ~ x + z1,
#'         family1 = gaussian(),
#'         family2 = gaussian(), data = df)
#' glm(y ~ x + poly(z1, 2), data = df)
#'
#' # interaction model
#' unm_glm(y ~ x*z1 + u1,
#'         u1 ~ x*z1,
#'         family1 = gaussian(),
#'         family2 = gaussian(), data = df)
#' glm(y ~ x*z1, data = df)
#'
#'
#'
#' # a binomial-binomial model
#' (df <- runm(50,
#'             missing_prop = .75,
#'             response = "bin",
#'             unmeasured_fam_list = list("bin"),
#'             unmeasured_param_list = list(.5)))
#' (unm_mod <- unm_glm(
#'   y ~ .,
#'   u1 ~ . - y,
#'   family1 = binomial(),
#'   family2 = binomial(),
#'   data = df
#' ))
#' glm(y ~ . - u1, family = binomial(), data = df)
#'
#'
#'
#' # a poisson-normal model
#' (df <- runm(25,
#'             response = "pois",
#'             response_model_coefs = c("int" = -1, "z" = .5, "u1" = .5, "x" = .5),
#'             treatment_model_coefs = c("int" = -1, "z" = .5, "u1" = .5),
#'             covariate_fam_list = list("norm"),
#'             covariate_param_list = list(c(mean = 0, sd = 1)),
#'             unmeasured_fam_list = list("norm"),
#'             unmeasured_param_list = list(c(0, 1))))
#'
#' (unm_mod <- unm_glm(
#'   y ~ x + z + u1 + offset(log(t)),
#'   u1 ~ x + z,
#'   family1 = poisson(),
#'   family2 = gaussian(),
#'   data = df
#' ))
#' glm(y ~ x + z + offset(log(t)), family = poisson(), data = df)
#'
#'
#'
#' # a poisson-binomial model
#' (df <- runm(25,
#'             response = "pois",
#'             response_model_coefs = c("int" = -1, "z" = .5, "u1" = .5, "x" = .5),
#'             treatment_model_coefs = c("int" = -1, "z" = .5, "u1" = .5),
#'             covariate_fam_list = list("norm"),
#'             covariate_param_list = list(c(mean = 0, sd = 1)),
#'             unmeasured_fam_list = list("bin"),
#'             unmeasured_param_list = list(.5)))
#' (unm_mod <- unm_glm(
#'   y ~ x + z + u1 + offset(log(t)), family1 = poisson(),
#'   u1 ~ x + z,                      family2 = binomial(),
#'   data = df
#' ))
#' glm(y ~ x + z + offset(log(t)), family = poisson(), data = df)
#'
#'
#'
#' # a gamma-normal model
#' (df <- runm(25,
#'             response = "gam",
#'             response_model_coefs = c("int" = -1, "z" = .5, "u1" = .5, "x" = .5),
#'             treatment_model_coefs = c("int" = -1, "z" = .5, "u1" = .5),
#'             covariate_fam_list = list("norm"),
#'             covariate_param_list = list(c(mean = 0, sd = 1)),
#'             unmeasured_fam_list = list("norm"),
#'             unmeasured_param_list = list(c(0, 1))))
#' (unm_mod <- unm_glm(
#'   y ~ x + z + u1, family1 = Gamma(),
#'   u1 ~ x + z, family2 = gaussian(),
#'   data = df
#' ))
#' glm(y ~ x + z, family = Gamma(link = "log"), data = df)
#'
#'
#'
#' # a gamma-binomial model
#' (df <- runm(25,
#'             response = "gam",
#'             response_model_coefs = c("int" = -1, "z" = .5, "u1" = .5, "x" = .5),
#'             treatment_model_coefs = c("int" = -1, "z" = .5, "u1" = .5),
#'             covariate_fam_list = list("norm"),
#'             covariate_param_list = list(c(mean = 0, sd = 1)),
#'             unmeasured_fam_list = list("bin"),
#'             unmeasured_param_list = list(.5)))
#' (unm_mod <- unm_glm(
#'   y ~ x + z + u1, family1 = Gamma(),
#'   u1 ~ x + z, family2 = binomial(),
#'   data = df
#' ))
#' glm(y ~ x + z, family = Gamma(link = "log"), data = df)
#' print(df, n = 25)
#'
#'
#'
#' # the output of unm_glm() is classed jags output
#' (df <- runm(20, response = "norm"))
#' (unm_mod <- unm_glm(y ~ ., u1 ~ . - y, family1 = gaussian(), family2 = gaussian(), data = df))
#' class(unm_mod)
#' jags_code(unm_mod)
#' unm_glm(y ~ ., u1 ~ . - y, data = df, code_only = TRUE)
#'
#'
#'
#'
#'
#' # visualizing output
#' library("ggplot2")
#' library("bayesplot"); bayesplot_theme_set(ggplot2::theme_minimal())
#' mcmc_hist(unm_mod, facet_args = list(labeller = label_parsed))
#' mcmc_hist(unm_mod)
#' mcmc_trace(unm_mod, facet_args = list(labeller = label_parsed))
#'
#' # more extensive visualization with the tidyverse
#' mcmc_intervals(unm_mod, prob = .90) +
#'   geom_point(
#'     aes(value, name), data = tibble::enframe(attr(df, "params")),
#'     color = "red", fill = "pink", size = 4, shape = 21
#'   )
#'
#'
#' library("dplyr")
#' library("tidyr")
#' unm_mod %>% as.matrix() %>% as_tibble() %>%
#'   pivot_longer(everything(), names_to = "var", values_to = "val") %>%
#'   ggplot(aes("0", val)) +
#'     geom_jitter() +
#'     geom_point(
#'       aes("0", value), data = tibble::enframe(attr(df, "params"), name = "var"),
#'       color = "red", fill = "pink", size = 4, shape = 21
#'     ) +
#'     coord_flip() +
#'     facet_grid(var ~ ., scales = "free_y", labeller = label_parsed) +
#'     theme_bw() +
#'     theme(
#'       axis.title = element_blank(),
#'       axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#'       strip.text.y = element_text(angle = 0)
#'     )
#'
#'
#' # getting draws out
#' (samps <- posterior::as_draws_df(unm_mod))
#' samps$`.chain`
#' samps$`.iteration`
#' samps$`.draw`
#'
#'
#'
#' # implementation is variable-name independent
#' (df <- runm(100, response = "norm"))
#' df$ht <- df$y
#' df$age <- df$u1
#' df$biom <- df$x
#' (unm_mod <- unm_glm(ht ~ x + biom + age,
#'                     age ~ x + biom,
#'                     family1 = gaussian(),
#'                     family2 = gaussian(), data = df))
#' jags_code(unm_mod)
#'
#' # ~~ Two Unmeasured Confounders Examples (III-Stage Model) ~~
#' # a normal-normal-normal model - internal validation
#' (df <- runm(50,
#'             missing_prop = .75,
#'             response = "norm",
#'             response_model_coefs = c("int" = -1, "z" = .5,
#'                                      "u1" = .5, "u2" = .5, "x" = .5),
#'             treatment_model_coefs = c("int" = -1, "z" = .5,
#'                                       "u1" = .5, "u2" = .5),
#'             covariate_fam_list = list("norm"),
#'             covariate_param_list = list(c(mean = 0, sd = 1)),
#'             unmeasured_fam_list = list("norm", "norm"),
#'             unmeasured_param_list = list(c(0, 1), c(0, 1))))
#' (unm_mod <- unm_glm(y ~ x + z + u1 + u2,
#'                     u1 ~ x + z + u2,
#'                     u2 ~ x + z,
#'                     family1 = gaussian(),
#'                     family2 = gaussian(),
#'                     family3 = gaussian(), data = df))
#' glm(y ~ x + z, data = df)
#' coef(unm_mod)
#'
#' unm_glm(y ~ x + z + u1 + u2,
#'         u1 ~ x + z + u2,
#'         u2 ~ x + z,
#'         family1 = gaussian(),
#'         family2 = gaussian(),
#'         family3 = gaussian(), data = df, code_only = TRUE)
#'
#'
#'
#' # a normal-normal-normal model - external validation
#' (df <- runm(c(20, 20),
#'             type = "ext",
#'             response = "norm",
#'             response_model_coefs = c("int" = -1, "z" = .5,
#'                                      "u1" = .5, "u2" = .5, "x" = .5),
#'             treatment_model_coefs = c("int" = -1, "z" = .5,
#'                                       "u1" = .5, "u2" = .5),
#'             covariate_fam_list = list("norm"),
#'             covariate_param_list = list(c(mean = 0, sd = 1)),
#'             unmeasured_fam_list = list("norm", "norm"),
#'             unmeasured_param_list = list(c(0, 1), c(0, 1))))
#' (unm_mod <- unm_glm(y ~ x + z + u1 + u2,
#'                     u1 ~ x + z + u2,
#'                     u2 ~ x + z,
#'                     family1 = gaussian(),
#'                     family2 = gaussian(),
#'                     family3 = gaussian(), data = df))
#'
#'
#'
#' # a binomial-binomial-binomial model - internal validation
#' (df <- runm(25,
#'             response = "bin",
#'             response_model_coefs = c("int" = -1, "z" = .5,
#'                                      "u1" = .5, "u2" = .5, "x" = .5),
#'             treatment_model_coefs = c("int" = -1, "z" = .5,
#'                                       "u1" = .5, "u2" = .5),
#'             covariate_fam_list = list("norm"),
#'             covariate_param_list = list(c(mean = 0, sd = 1)),
#'             unmeasured_fam_list = list("bin", "bin"),
#'             unmeasured_param_list = list(.5, .75)))
#' unm_glm(y ~ x + z + u1 + u2, family1 = binomial(),
#'         u1 ~ x + z + u2, family2 = binomial(),
#'         u2 ~ x + z, family3 = binomial(),
#'         data = df, code_only = TRUE)






#' @export
#' @rdname unm_glm

unm_glm <- function(
    form1, form2 = NA, form3 = NA,
    family1 = binomial(), family2 = NA, family3 = NA,
    data,
    n.iter = 2000, n.adapt = 1000, thin = 1, n.chains = 4,
    filename = tempfile(fileext = ".jags"),
    quiet = getOption("unm_quiet"),
    progress.bar = getOption("unm_progress.bar"),
    code_only = FALSE,
    default_prior = NULL,
    priors,
    response_nuisance_priors, response_params_to_track,
    confounder1_nuisance_priors, confounder1_params_to_track,
    confounder2_nuisance_priors, confounder2_params_to_track,
    ...
) {

  if (is.null(default_prior)) {
    default_prior <- switch(family1$family,
                            "gaussian" = "dnorm(0, .001)",
                            "binomial" = "dnorm(0, .1)",
                            "Gamma" = "dnorm(0, .001)",
                            "poisson" = "dnorm(0, .1)")
  }

  g <- glue::glue

  if (!inherits(family2, "family")) family2 <- list("family" = "none")
  if (!inherits(family3, "family")) family3 <- list("family" = "none")

  y <- deparse(form1[[2]]) # e.g. "y", character name of response var
  u1 <- if(inherits(form2, "formula")) deparse(form2[[2]]) else NULL # e.g. "u1", character name of confounding var
  u2 <- if(inherits(form3, "formula")) deparse(form3[[2]]) else NULL # e.g. "u2", character name of confounding var

  if(!is.null(u1)) {
    conf_piece <- "+ inprod(U[i,], lambda)"
  } else {
    conf_piece <- ""
  }

  response_model_code <- switch(family1$family,
                                "gaussian" = g("{y}[i] ~ dnorm(mu_{y}[i], tau_{y})
{'    '}mu_{y}[i] <- inprod(X[i,], beta) {conf_piece}"),
                                "binomial" = g("{y}[i] ~ dbern(p_{y}[i])
{'    '}logit(p_{y}[i]) <- inprod(X[i,], beta) {conf_piece}"),
                                # note: r parameterizes the gamma as shape al and rate la, with mean = al / la
                                #    jags parameterizes the gamma as shape  r and rate la, with mean =  r / la, the same
                                "Gamma" = g("{y}[i] ~ dgamma(alpha_{y}, d[i])
{'    '}d[i] <- alpha_{y} / mu_{y}[i]
{'    '}log(mu_{y}[i]) <- inprod(X[i,], beta) {conf_piece}"),
                                "poisson" = g("{y}[i] ~ dpois(mu_{y}[i])
{'    '}log(mu_{y}[i]) <- log_e[i] + inprod(X[i,], beta) {conf_piece}"))


  if (missing(response_nuisance_priors)) {

    response_nuisance_priors <- switch(family1$family,
                                       "gaussian" = g("tau_{y} ~ dt(0, 1, 3) I(0, )"),
                                       "binomial" = " ",
                                       "Gamma" = g("alpha_{y} ~ dgamma(.1, .1)"),
                                       "poisson" = " "
    )

    response_params_to_track <- switch(family1$family,
                                       "gaussian" = g("tau_{y}"),
                                       "binomial" = character(0),
                                       "Gamma" = g("alpha_{y}"),
                                       "poisson" = character(0)
    )

  } else {

    response_nuisance_priors <- g(response_nuisance_priors)
    response_params_to_track <- g(response_params_to_track)

  }

  if(!is.null(u2)) {
    conf2_piece <- "+ inprod(U[i, 2], zeta)"
  } else {
    conf2_piece <- ""
  }

  confounder1_model_code <- switch(family2$family,
                                   "gaussian" = g("U[i,1] ~ dnorm(mu_{u1}[i], tau_{u1})
                   {'    '}mu_{u1}[i] <- inprod(W[i,], gamma) {conf2_piece}"),
                                   "binomial" = g("U[i,1] ~ dbern(p_{u1}[i])
                   {'    '}logit(p_{u1}[i]) <- inprod(W[i,], gamma) {conf2_piece}"),
                                   "none" = "")

  if (missing(confounder1_nuisance_priors)) {

    confounder1_nuisance_priors <- switch(family2$family,
                                          # "gaussian" = g("tau_{u1} <- 1 / ( sigma_{u1} * sigma_{u1} ); sigma_{u1} ~ dunif(0, 100)"),
                                          "gaussian" = g("tau_{u1} ~ dt(0, 1, 3) I(0, )"),
                                          "binomial" = " ",
                                          "none" = ""
    )

    confounder1_params_to_track <- switch(family2$family,
                                          # "gaussian" = g("sigma_{u1}"),
                                          "gaussian" = g("tau_{u1}"),
                                          "binomial" = character(0),
                                          "none" = character(0)
    )

  } else {

    confounder1_nuisance_priors <- g(confounder1_nuisance_priors)
    confounder1_params_to_track <- g(confounder1_params_to_track)

  }

  confounder2_model_code <- switch(family3$family,
                                   "gaussian" = g("U[i,2] ~ dnorm(mu_{u2}[i], tau_{u2})
                 {'    '}mu_{u2}[i] <- inprod(V[i,], delta)"),
                                   "binomial" = g("U[i,2] ~ dbern(p_{u2}[i])
                 {'    '}logit(p_{u2}[i]) <- inprod(V[i,], delta)"),
                                   "none" = "")

  if (missing(confounder2_nuisance_priors)) {

    confounder2_nuisance_priors <- switch(family3$family,
                                          # "gaussian" = g("tau_{u1} <- 1 / ( sigma_{u1} * sigma_{u1} ); sigma_{u1} ~ dunif(0, 100)"),
                                          "gaussian" = g("tau_{u2} ~ dt(0, 1, 3) I(0, )"),
                                          "binomial" = " ",
                                          "none" = ""
    )

    confounder2_params_to_track <- switch(family3$family,
                                          # "gaussian" = g("sigma_{u1}"),
                                          "gaussian" = g("tau_{u2}"),
                                          "binomial" = character(0),
                                          "none" = character(0)
    )

  } else {

    confounder2_nuisance_priors <- g(confounder2_nuisance_priors)
    confounder2_params_to_track <- g(confounder2_params_to_track)

  }

  # make design matrices
  # https://stackoverflow.com/questions/5616210/model-matrix-with-na-action-null
  current_na_action <- getOption("na.action")
  on.exit(options(na.action = current_na_action), add = TRUE)
  options(na.action = "na.pass")
  (X <- model.matrix(form1, data = data))
  if(!is.null(u1)) (W <- model.matrix(form2, data = data)) else (W <- NULL)
  if(!is.null(u2)) (V <- model.matrix(form3, data = data)) else (V <- NULL)


  # partition response model matrix into observed and unobserved parts
  if (is.null(u1)) {
    (U <- NULL)
  } else if (is.null(u2)) {
    (U <- X[, c(u1), drop = FALSE])
  } else {
    (U <- X[, c(u1, u2), drop = FALSE])
  }
  (X <- X[, setdiff(colnames(X), colnames(U)), drop = FALSE])

  (U2 <- W[, c(u2), drop = FALSE])
  (W <- W[, setdiff(colnames(W), c(u2)), drop = FALSE])

  # determine numbers of parameters
  p_be <- ncol(X)  # = # non-confounder params in response model
  p_la <- ncol(U)  # = # confounder params in response model
  p_ga <- ncol(W)  # = # non-confounder params in unmeasured1 model
  p_ze <- ncol(U2) # = # confounder params in unmeasured1 model
  p_de <- ncol(V)  # = # non-confounder params in unmeasured2 model

  # make conversion from, e.g., beta[2] to beta[x]
  X_vars <- colnames(X) # corresponds to beta
  U_vars <- colnames(U) # corresponds to lambda
  W_vars <- colnames(W) # corresponds to gamma
  U2_vars <- colnames(U2) # corresponds to zeta
  V_vars <- colnames(V) # corresponds to delta

  pretty_X_vars <- gsub("\\(Intercept\\)", "1", X_vars)
  pretty_W_vars <- gsub("\\(Intercept\\)", "1", W_vars)
  pretty_V_vars <- gsub("\\(Intercept\\)", "1", V_vars)

  # make priors
  jags_coefs <-
    if(is.null(u1)){
      c( g("beta[{1:p_be}]") )
    } else if(is.null(u2)){
      c( g("beta[{1:p_be}]"), g("lambda[{1:p_la}]"), g("gamma[{1:p_ga}]") )
    } else {
      c( g("beta[{1:p_be}]"), g("lambda[{1:p_la}]"), g("gamma[{1:p_ga}]"),
         g("zeta[{1:p_ze}]"), g("delta[{1:p_de}]") )
    }
  real_coefs <- c( g("beta[{X_vars}]"), g("lambda[{U_vars}]"), g("gamma[{W_vars}]"),
                   g("zeta[{U2_vars}]"), g("delta[{V_vars}]") )
  real_coefs <- gsub("\\(Intercept\\)", "1", real_coefs)

  # make function to convert coefficient names, e.g. gamma[x] -> gamma[2]
  real_to_jags_coefs <- function(x) {
    dict <- structure(jags_coefs, names = real_coefs)
    x[x %in% real_coefs] <- unname(dict[x])[x %in% real_coefs]
    x
  }
  # real_to_jags_coefs(c("beta[1]", "gamma[x]", "a", "lambda[u]"))
  # [1] "beta[1]" "gamma[2]" "a" "lambda[1]"

  jags_to_real_coefs <- function(x) {
    dict <- structure(real_coefs, names = jags_coefs)
    x[x %in% jags_coefs] <- unname(dict[x])[x %in%jags_coefs]
    x
  }
  # jags_to_real_coefs(c("beta[1]", "gamma[2]", "a", "lambda[1]"))
  # [1] "beta[1]" "gamma[x]" "a" "lambda[u]"

  jags_priors <- rep(default_prior, length(jags_coefs))
  names(jags_priors) <- jags_coefs
  if (!missing(priors)) jags_priors[real_to_jags_coefs(names(priors))] <- priors


  # make jags priors code
  jags_priors <- g("{names(jags_priors)} ~ {jags_priors} \t\t# = {real_coefs}")


  # make jags data list
  jd <- list(
    "n" = nrow(X), "X" = X, "U" = U, "W" = W, "V" = V,
    "log_e" = model.offset(model.frame(form1, data = data))
  )
  jd[[y]] <- data[[y]]

  # this wipes the log_e if not present
  jd <- drop_nulls(jd)


  # initialize chain - commented out because need to make compatible with
  # user prior specification
  # inits <- replicate(
  #   n.chains,
  #   list("beta" = rnorm(p_be), "lambda" = rnorm(p_la), "gamma" = rnorm(p_ga)),
  #   simplify = FALSE
  # )


  # write the jags file
  code <- g("
  model {
    # models
    for (i in 1:n) {
      {{response_model_code}}

      {{confounder1_model_code}}
      {{confounder2_model_code}}
    }

    # priors
    {{paste(jags_priors, collapse = '\n    ')}}
    {{response_nuisance_priors}}
    {{confounder1_nuisance_priors}}
    {{confounder2_nuisance_priors}}
  }
  ", .open = "{{", .close = "}}")
  if (code_only) {
    cat(code)
    cat("\n")
    return(invisible())
  }
  cat(code, file = filename)

  # compile chain and adapt
  jm <- jags.model(filename, data = jd,
                   n.adapt = n.adapt, n.chains = n.chains,
                   quiet = TRUE, ...
  )

  params_of_interest <- unique(sub("\\[.+\\]", "", real_coefs))

  # sample
  samps <- coda.samples(
    jm,
    c(
      params_of_interest,
      response_params_to_track,
      confounder1_params_to_track,
      confounder2_params_to_track
    ),
    n.iter = n.iter, thin = thin, progress.bar = progress.bar
  )


  # fix names - does jags simplify lambda[1] to just lambda?
  names <- dimnames(samps[[1]])[[2]]
  names[names == "beta"] <- "beta[1]"
  names[names == "lambda"] <- "lambda[1]"
  names[names == "gamma"] <- "gamma[1]"
  names[names == "delta"] <- "delta[1]"
  names[names == "zeta"] <- "zeta[1]"
  names <- jags_to_real_coefs(names)


  # format nuisance parameters
  names <- gsub("(\\w+)_(\\w+)", "\\1[\\2]", names)
  # cbind(dimnames(samps[[1]])[[2]], names)

  # change names in samples output
  samps <- lapply(samps, function(mcmc) {
    dimnames(mcmc) <- list(NULL, names)
    mcmc
  })
  class(samps) <- c("unm_int", "unm_mod", "mcmc.list")



  # add metadata
  attr(samps, "file") <- filename
  attr(samps, "code") <- code
  attr(samps, "form1") <- form1
  attr(samps, "family1") <- family1
  attr(samps, "form2") <- form2
  attr(samps, "family2") <- family2
  attr(samps, "form3") <- form3
  attr(samps, "family3") <- family3
  attr(samps, "call") <- match.call()
  attr(samps, "jm") <- jm
  attr(samps, "n.iter") <- n.iter
  attr(samps, "n.adapt") <- n.adapt
  attr(samps, "thin") <- thin
  attr(samps, "n.chains") <- n.chains


  # return
  samps

}




#' @param mod The output of [unmconf::unm_glm()]
#'
#' @export
#' @rdname unm_glm
jags_code <- function(mod) {
  stopifnot(inherits(mod, "mcmc.list"), !is.null(attr(mod, "code")))
  cat(attr(mod, "code"))
  invisible(mod)
}


#' @param x Object to be printed
#' @param digits Number of digits to round to; defaults to 3
#' @param print_call Should the call be printed? Defaults to `TRUE`, but can be
#'   turned off with `options("unm_print_call" = FALSE)`
#' @export
#' @rdname unm_glm
print.unm_int <- function(x, digits = 3, ...,
                          print_call = getOption("unm_print_call")) { #

  # print the call
  call <- attr(x, 'call')
  n <- length(call)

  if (print_call) {
    cat("Call:\n")
    cat("", call[[1]], "(\n", sep = "")
    for (k in 2:(n-1)) {
      cat("  ", glue::glue("{names(call)[k]} = {deparse(call[[k]])},"), "\n")
    }
    cat("  ", glue::glue("{names(call)[n]} = {deparse(call[[n]])}"), "\n")
    cat(")\n")
    cat("\n")
  }


  # print the summaries
  summary <- summary(x, quantiles = numeric())
  ests <- summary$statistic[,1]
  one_to_intercept <- function(x) { x[x == "1"] <- "(Intercept)"; x }

  cat("Response model coefficients:\n")
  response_ests <- ests[grepl("(beta|lambda)", names(ests))]
  names(response_ests) <- one_to_intercept(ez_extract_subset(names(response_ests)))
  print(round(response_ests, digits = digits))
  cat("\n")

  cat("Confounder model coefficients 1:\n")
  confounder1_ests <- ests[grepl("gamma|zeta", names(ests))]
  names(confounder1_ests) <- one_to_intercept(ez_extract_subset(names(confounder1_ests)))
  print(round(confounder1_ests, digits = digits))
  cat("\n")

  cat("Confounder model coefficients 2:\n")
  confounder2_ests <- ests[grepl("delta", names(ests))]
  names(confounder2_ests) <- one_to_intercept(ez_extract_subset(names(confounder2_ests)))
  print(round(confounder2_ests, digits = digits))
  cat("\n")

  other_ests <- ests[!grepl("(beta|lambda|gamma|delta|zeta)", names(ests))]
  if (length(other_ests) > 0) {
    cat("Other estimates:\n")
    print(round(other_ests, digits = digits))
  }
  # cat("\n")
  #
  # mod_name <- eval(deparse(substitute(x)), envir = parent.frame())
  # cat(g("To compute the deviance use `unm_dic({mod_name})`"))

  invisible(x)

}



#' @param object Model object for which the coefficients are desired
#' @export
#' @rdname unm_glm
coef.unm_int <- function(object, ...){
  summary <- summary(object, quantiles = numeric())
  summary$statistic[,1]
}

