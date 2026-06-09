## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~ Preference and switch-point formulas (shared) ~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Sourced by RunMe.R before visualize_stan_model.R and compare_models.R so that
## both consumers use a single definition of the preference and switch-point math.
##
##   model_family(model_name)               "Logistic" or "vanLeeuwen"
##   preference(x, w, q, model_name)        drift preference at log-ratio x
##   log_switch_point(y, w, q, model_name)  log-ratio at which preference = y
##
## Model family is resolved from the model name (the same string used to locate
## stan_model_<name>.stan and model_output_<name>.RDS). Names starting with
## "Logistic" use the Yodzis logistic preference form; names starting with
## "vanLeeuwen" use the vanLeeuwen form. Add new families here when introducing
## additional model variants.
##
## All numeric arguments may be vectors; expressions are elementwise arithmetic.

model_family <- function(model_name) {
  if (grepl("^Logistic",   model_name)) return("Logistic")
  if (grepl("^vanLeeuwen", model_name)) return("vanLeeuwen")
  stop("Unknown preference family for model '", model_name,
       "'. Expected the name to start with 'Logistic' or 'vanLeeuwen'; ",
       "extend model_family() in preference_helpers.R to register new families.",
       call. = FALSE)
}

preference <- function(x, w, q, model_name) {
  if (model_family(model_name) == "Logistic") {
    1 - 1 / (1 + exp(w + q * x))
  } else {
    1 - (1 + exp((q - 4) + x)) /
        (1 + exp(log(2) + (q - 4) + x) + exp(q + 2 * x))
  }
}

log_switch_point <- function(y, w, q, model_name) {
  if (model_family(model_name) == "Logistic") {
    (log(-y / (y - 1)) - w) / q
  } else {
    log(- (2 * y) /
        (exp(q - 4) * (2 * y - 1) -
         exp(q) * sqrt(exp(-2 * q) *
                       (exp(2 * (q - 4)) * (1 - 2 * y)^2 -
                        4 * exp(q) * (y - 1) * y))))
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
