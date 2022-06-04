#' Estimate a causal effect in a convenience sample
#'
#' Propensity adjusted causal effects are estimated for convenience samples and
#' sampling bias is adjusted for. Sampling weights are estimated for the
#' convenience sample using an auxiliary dataset. The weights are estimated with
#' logistic regression and stepwise forward selection. The estimated sampling weights
#' are then used to estimate propensity scores. The propensity scores are estimated
#' using logistic regression and stepwise selection with quadratic terms and two-way
#' interactions. The causal effect is estimated with a propensity adjusted GLM that
#' uses the estimated sampling weights. The estimated causal effect is returned with a
#' standard error that is estimated using a simultaneous estimating equation that accounts
#'  for uncertainty from the weight estimating and propensity score estimation process.
#'  More details can be found in our forthcoming paper Propensity Scores in
#'  Convenience Samples by Bernstein Morgan et al.
#'
#' @param convSamp
#' @param SRS
#' @param sampwt_vars
#' @param PS_vars
#' @param treatment_var
#' @param response_var
#' @param outcome_family
#'
#' @return
#' @export
#'
#' @examples
convPS = function(convSamp, SRS,
                  sampwt_vars, PS_vars,
                  treatment_var, response_var,
                  outcome_family = stats::gaussian){


  #### check function inputs ####
  # check that the vars in  sampwt_vars, PS_vars, treatment_var, response_var
  # are in the correct data frames
  if(sum(sampwt_vars %in% colnames(convSamp)) != length(sampwt_vars)){
    stop("Data set up incorrectly. Make sure the variables in sampwt_vars are in the convSamp data frame")
  }
  if(sum(sampwt_vars %in% colnames(SRS)) != length(sampwt_vars)){
    stop("Data set up incorrectly. Make sure the variables in sampwt_vars are in the SRS data frame")
  }
  if(sum(PS_vars %in% colnames(convSamp)) != length(PS_vars)){
    stop("Data set up incorrectly. Make sure the variables in PS_vars are in the convSamp data frame")
  }
  if(sum(treatment_var %in% colnames(convSamp)) != length(treatment_var)){
    stop("Data set up incorrectly. Make sure the treatment_var variable is in the convSamp data frame")
  }
  if(sum(response_var %in% colnames(convSamp)) != length(response_var)){
    stop("Data set up incorrectly. Make sure the response_var variable is in the convSamp data frame")
  }
  # check that the treatment variable (treatment_var) is binary
  if(sum(unique(convSamp[,treatment_var]) == c(1,0)) == 2){
    stop("Treatment variable (treatment_var) should be a binary variable where 1 denotes the treatment group and 0 denotes the control group.")
  }

  ##### Estimate sampling weights #####

  # format data for estweight function and create biased sample indicator
  r.samp = data.frame(SRS[,sampwt_vars], biased = 0)
  b.samp = data.frame(convSamp[,sampwt_vars], biased = 1)
  names(r.samp)[1:length(sampwt_vars)] = sampwt_vars
  names(b.samp)[1:length(sampwt_vars)] = sampwt_vars

  # concatenate 2 datasets
  Xcomb = rbind(b.samp, r.samp)
  Xfit = data.frame(ID = 1:nrow(Xcomb), Xcomb)

  # estimate sampling weights
  weights_and_fit = estweight(data = Xfit, weight_model = "logistic",
                              return_model_fit = TRUE)
  convSamp$htweight = weights_and_fit$weights$htweight
  estwt_fit = weights_and_fit$estwt_fit

  #### Estimate propensity scores ####
  if(length(PS_vars)==1){
    fs = is.factor(convSamp[,PS_vars])
    if(is.factor(convSamp[,PS_vars])){
      fact_vars_prop = PS_vars
      cont_vars = NULL
    } else{
      fact_vars_prop = NULL
      cont_vars_prop = PS_vars
    }
  } else{
    fs = sapply(convSamp[,PS_vars],is.factor)
    fact_vars_prop = PS_vars[fs]
    cont_vars_prop = PS_vars[!fs]
  }

  form_min_prop = stats::as.formula(paste(treatment_var," ~ 1"))

  interactions = paste0(
    "(",paste(c(cont_vars_prop, fact_vars_prop),collapse = "+"),
    ")^2")
  form_max_prop = stats::as.formula(paste0(treatment_var,"~",
                                           paste0(interactions, "+",
                                                  paste0(paste0("I(",cont_vars_prop,"^2)"),collapse = "+"))))

  # fit_min = do.call("svyglm", args = list(form_min_prop,
  #                                         design = survey::svydesign(ids = ~0, weights = convSamp$htweight, data = convSamp),
  #                                         family = stats::quasibinomial))

  fit_min = survey::svyglm(form_min_prop,
                           design = survey::svydesign(ids = ~0,
                                                      weights = convSamp$htweight,
                                                      data = convSamp),
                           family = stats::quasibinomial)

  attach(convSamp) # using dAIC with svyglm() doesn't pass data the same way glm() does
  forward = stats::step(fit_min,scope=list(lower=form_min_prop,upper=form_max_prop),
                        direction="forward", trace = 0)
  detach(convSamp)

  # fit selected model
  estprop_form = stats::formula(forward)
  estprop = survey::svyglm(estprop_form,
                           design = survey::svydesign(ids = ~0, weights = convSamp$htweight, data = convSamp),
                           family = stats::quasibinomial)

  convSamp$propscore = fitted(estprop, "response")


  #### Estimate causal effect ####
  outcome = svyglm(paste0(response_var,"~",treatment_var,"+ propscore"),
                   design = survey::svydesign(ids = ~0, weights = convSamp$htweight, data = convSamp),
                   family = outcome_family)
  coef.save = outcome$coefficients

  # proposed SE estimate
  se.prop.save = treatmentSE(estwt_fit = estwt_fit,
                             estprop_fit = estprop,
                             fit_outcome = outcome,
                             biased = Xfit$biased,
                             outcome_family = outcome_family)

  out = list(Treatment_effect_est = coef.save[2],
             Treatment_effect_SE = se.prop.save[2])

  return(out)
}
