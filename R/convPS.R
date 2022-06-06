#' Estimate a causal effect in a convenience sample
#'
#' Propensity adjusted causal effects are estimated for convenience samples and
#' sampling bias is adjusted for, with standard errors that accounts for uncertainy
#' in estimating sampling weights and propensity scores.
#'
#' @param convSamp Convenience sample stored as a data frame. It should contain variables for
#' estimating the sampling weight, for estimating the propensity score, the treatment variable, and the response variable.
#' Note that all factors in the dataset should either be binary or of class factor.
#' @param repSamp Representative sample stored as a data frame. It should contain variables for
#' estimating the sampling weight. The variables for estimating the sampling weight should be named
#' the same as they are in the convenience sample (convSamp).
#' Note that all factors in the dataset should either be binary or of class factor.
#' @param sampwt_vars A vector of column names used for estimating the sampling weights. These variables
#' should be present in both the convenience sample (convSamp) and the representative sample (repSamp).
#' @param PS_vars A vector of column names used for estimating the sampling weights. These variables
#' should be present in the convenience sample (convSamp).
#' @param treatment_var The column name for the treatment variable. This variable
#' should be present in the convenience sample (convSamp).
#' @param response_var The column name for the response variable. This variable
#' should be present in the convenience sample (convSamp).
#' @param outcome_family Error distribution and link function for the response variable. Defaults to a gaussian distribution with an identity link. (See
#' \code{family} for details of family function)
#'
#' @return \code{convPS} returns a list containing the estimated treatment effect
#' (Treatment_effect_est) and corresponding standard error (Treatment_effect_SE)
#' that accounts for uncertainty from estimating the sampling weights and propensity scores
#' @export
#'
#' @details Sampling weights are estimated for the
#' convenience sample using an auxiliary dataset. The weights are estimated with
#' logistic regression and stepwise forward selection. The estimated sampling weights
#' are then used to estimate propensity scores. The propensity scores are estimated
#' using logistic regression and stepwise selection with quadratic terms and two-way
#' interactions. The causal effect is estimated with a propensity adjusted GLM that
#' uses the estimated sampling weights. The estimated causal effect is returned with a
#' standard error that is estimated using a simultaneous estimating equation that accounts
#'  for uncertainty from the weight estimating and propensity score estimation process.
#'
#'  @references “Causal Inference in Convenience Samples.”
#'  Bernstein, OM; Vegetabile, BG; Grill, JD; Gillen, DL. Forthcoming.
#'
#' @examples
#'# Simulated a convenience sample and a representative sample
#'require(MASS)
#'expit = function(x){exp(x)/(1+exp(x))}
#'n.pop = 10000; n = 2000
#'
#'Sigma = matrix(c(1,.9,.9,1),nrow = 2)
#'cont = mvrnorm(n.pop, mu = c(0,0), Sigma = Sigma)
#'K.full = ifelse(cont[,1] > 0, 1, 0)
#'x1.full = ifelse(cont[,2] > 0, 1, 0)
#'prob.sample = .6*K.full + .2
#'
#'# get convenince sample
#'CSamp = sample(1:n.pop, size = n, prob = prob.sample)
#'x1.c = x1.full[CSamp]
#'K.c = K.full[CSamp]
#'
#'# get simple random sample (SRS)
#'SRS = sample(1:n.pop, size = n)
#'x1.s = x1.full[SRS]
#'K.s = K.full[SRS]
#'
#'# get remaining covariates for biased sample
#'x2.c = rnorm(n,0,2)
#'x3.c = rnorm(n,0,2)
#'pi.c = expit((log(1.3)*x2.c + log(.4)*x3.c)*(1-K.c) +
#' (log(2)*x2.c + log(1.5)*x3.c)*(K.c))
#' T.c = rbinom(n,1,pi.c)
#' x4.c = rnorm(n,0,1)
#' x5.c = rnorm(n,0,1)
#' y.c = rnorm(n, mean = (0 + 1*T.c + 3*T.c*K.c + 1.5*x2.c - 2*x3.c^2 - 1*x4.c + 1.5*x5.c^3), sd = 1)
#' convSamp = data.frame(x1 = x1.c, x2 = x2.c, x3 = x3.c,
#'  x4 = x4.c, x5 = x5.c, Tx = T.c, y = y.c)
#'  repSamp = data.frame(x1 = x1.s)
#'
#'convPS(convSamp = convSamp, repSamp = repSamp,
#'sampwt_vars = "x1", PS_vars = paste0("x",1:4),
#'treatment_var = "Tx", response_var = "y")
#'
#' @importFrom survey svyglm
#'
convPS = function(convSamp, repSamp,
                  sampwt_vars, PS_vars,
                  treatment_var, response_var,
                  outcome_family = stats::gaussian){


  #### check function inputs ####
  # check that the vars in  sampwt_vars, PS_vars, treatment_var, response_var
  # are in the correct data frames
  if(sum(sampwt_vars %in% colnames(convSamp)) != length(sampwt_vars)){
    stop("Data set up incorrectly. Make sure the variables in sampwt_vars are in the convSamp data frame")
  }
  if(sum(sampwt_vars %in% colnames(repSamp)) != length(sampwt_vars)){
    stop("Data set up incorrectly. Make sure the variables in sampwt_vars are in the repSamp data frame")
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
  if(sum(sort(unique(convSamp[,treatment_var])) == c(0,1)) != 2){
    stop("Treatment variable (treatment_var) should be a binary variable where 1 denotes the treatment group and 0 denotes the control group.")
  }

  ##### Estimate sampling weights #####

  # format data for estweight function and create biased sample indicator
  r.samp = data.frame(repSamp[,sampwt_vars], biased = 0)
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

  convSamp$propscore = stats::fitted(estprop, "response")


  #### Estimate causal effect ####
  outcome = survey::svyglm(paste0(response_var,"~",treatment_var,"+ propscore"),
                   design = survey::svydesign(ids = ~0, weights = convSamp$htweight, data = convSamp),
                   family = outcome_family)
  coef.save = outcome$coefficients

  # proposed SE estimate
  se.prop.save = treatmentSE(estwt_fit = estwt_fit,
                             estprop_fit = estprop,
                             fit_outcome = outcome,
                             biased = Xfit$biased,
                             outcome_family = outcome_family)

  out = list(Treatment_effect_est = coef.save[2][[1]],
             Treatment_effect_SE = se.prop.save[2])

  return(out)
}
