#' Convert to Greek expressions
#'
#' Convert to Greek expressions for plotting
#'
#' @param labs A character vector of greek symbols of the form `ga_x` and `be_1`.
#' @param s A character vector of Greek short hand codes, e.g. `"si"`.
#' @param mod Output from [unm_glm()].
#' @return A character vector.
#' @name helpers
#' @examples
#'
#' labs <- c("ga_1", "ga_treatment", "ga_x", "be_1",
#'   "be_treatment", "be_x", "la_u", "al_y", "si")
#' expand_labels(labs)
#'
#'


#' @rdname helpers
#' @export
expand_labels <- function(labs) {
  # e.g. "ga_x" -> "gamma[x]"
  vapply(
    strsplit(labs, "_"),
    function(s) {
      if (length(s) == 1L) return(greek_expander(s))
      sprintf("%s[%s]", greek_expander(s[1]), s[2])
    },
    character(1)
  )
}




#' @rdname helpers
#' @export
greek_expander <- function(s) {
  c(
    "si" = "sigma",
    "ga" = "gamma",
    "de" = "delta",
    "ze" = "zeta",
    "la" = "lambda",
    "be" = "beta",
    "al" = "alpha",
    "et" = "eta"
  )[s]
}


#' @rdname helpers
#' @export
make_greek_coefs <- function(mod) {
  structure(
    lapply(
      mod,
      function(mcmc) {
        attr(mcmc, "dimnames") <- list(NULL, expand_labels(attr(mcmc, "dimnames")[[2]]))
        mcmc
      }
    ),
    class = attr(mod, "class"),
    file = attr(mod, "file"),
    code = attr(mod, "code")
  )
}



ez_trunc <- function(string, width, ellipsis = "...") {
  too_long <- !is.na(string) & nchar(string) > width
  width... <- width - nchar(ellipsis)
  string[too_long] <- paste0(substr(string[too_long], 1, width...), ellipsis)
  string
}
# ez_trunc("Hello how are you today?", 10)



ez_distill <- function(string, with = " ") gsub("\\s+", with, string)



ez_trim <- function(string) {
  string <- gsub("^\\s+", "", string)
  gsub("\\s+$", "", string)
}
# ez_trim("      hello! ")
# ez_trim(c("      hello! ", "      hello! "))



ez_extract_subset <- function(string) {
  starts <- regexpr("\\[", string) + 1
  stops  <- regexpr("\\]", string) - 1
  substr(string, start = starts, stop = stops)
}




#' @param x Character vector.
#' @rdname helpers
#' @export
drop_nulls <- function(x) {
  x[vapply(x, function(.) !is.null(.), logical(1))]
}
# drop_nulls(c(1, NULL, 4))
# drop_nulls(list(1, NULL, 4))
