#' Estimate propensity weights for biased samples
#'
#' Estimate propensity weights for biased samples using information from a representative sample.
#'
#' @param data Combined dataset including the convenience sample and representative sample with a subject
#' identifier, \code{ID}, in the first column, any covariates needed for matching, and an indicator, \code{biased},
#' that is 1 for subjects in the biased sample in the last column. Note that all factors in the dataset
#' should either  be binary or of class factor
#' @param weight_model Which propensity weight estimation model should be used. Either \code{"logistic"},
#' \code{"randomForest"}, \code{"CBPS"}, or \code{"entbal"}
#' @param incl_wts_for_rep A logical value indicating whether propensity weights for
#' observations from the representative sample should be included in the output.
#'
#' @return \code{estweight} returns a dataframe containing IDs and the corresponding estimated
#' propensity weights
#' @export
#'
#' @examples
#' #' data("mtcars")
#' repsample = mtcars
#' n = nrow(repsample)
#' expit = function(x) {exp(x) / (1 + exp(x))}
#'
#' # Calculate probability of being oversampled
#' repsample$sampprob = expit(.01*(repsample$am*4 + repsample$carb*3
#' + repsample$drat*.9 -repsample$mpg*repsample$disp*.05 + .002*repsample$hp^2 + 80))
#'
#' # draw biased and representative samples
#' b.samp = repsample[sample(1:n, 500, prob = repsample$sampprob, replace = TRUE), ]
#' r.samp = repsample[sample(1:n, 500, replace = TRUE), ]
#'
#' # Create indicator of biased sample membership
#' b.samp$biased = 1; r.samp$biased = 0
#'
#'# Format data to pass to function
#'Xcomb = data.frame(ID = 1:(1000), rbind(b.samp, r.samp))
#'Xfit = Xcomb[,c("ID", colnames(Xcomb)[c(2:8,10:12)], "biased")]
#'
#'# Estimate propensity weights
#'estweight(data = Xfit, weight_model = "logistic")
estweight = function(data, weight_model = "logistic",
                     incl_wts_for_rep = FALSE){

  # Check inputs
  if(!(weight_model %in% c("logistic", "randomForest", "CBPS","entbal"))){
    warning("Invalid weight_model input. Will use default logistic weight estimation model")
    weight_model = "logistic"
  }

  cols = colnames(data)
  if(sum(cols%in%c("ID","biased"))!=2){
    stop("Data set up incorrectly. Make sure the data includes columns named 'ID' and 'biased'")
  }

  if(!is.logical(incl_wts_for_rep)){
    stop("incl_wts_for_rep needs to be a logical value")
  }

  # Identify covariate names and which are factors
  cov = cols[!(cols%in%c("ID","biased"))]
  if(length(cov)==1){
    fs = is.factor(data[,cov])
    if(is.factor(data[,cov])){
      fact_vars = cov
      cont_vars = NULL
    } else{
      fact_vars = NULL
      cont_vars = cov
    }
  } else{
    fs = sapply(data[,cov],is.factor)
    fact_vars = cov[fs]
    cont_vars = cov[!fs]
  }

  ## Estimate selection weights

  if(weight_model == "logistic"){

    # Use forward-selection AIC to select a model
    # Model scope includes second order terms
    form_min = stats::as.formula("biased ~ 1")
    form_max = stats::as.formula(paste0("biased ~",paste0(c(cov, paste0("I(",cont_vars,"^2)")),collapse = "+")))
    fit_min = stats::glm(formula = form_min, family = stats::binomial, data = data)
    forward = stats::step(fit_min,scope=list(lower=form_min,upper=form_max),
                          direction="forward", trace = 0)

    # fit selected model
    estwt_form = stats::formula(forward)
    estwt_fit = stats::glm(formula = estwt_form, family = stats::binomial, data = data)
    prob_bias = stats::fitted(estwt_fit, "response")
    htweight_unnorm = (1-prob_bias)/prob_bias
    htweight = htweight_unnorm/sum(htweight_unnorm)

  } else if(weight_model == "randomForest"){

    form_max = stats::as.formula(paste0("factor(biased) ~",paste0(cov,collapse = "+")))

    fit_rf = randomForest::randomForest(formula = form_max, data = data, type = "classification")

    prob_bias1 = stats::predict(fit_rf,type =  "prob")[,"1"]
    prob_bias = ifelse(prob_bias1==0,.01,ifelse(prob_bias1==1,.99,prob_bias1)) # deal with est wts of 0 or 1
    htweight_unnorm = (1-prob_bias)/prob_bias
    htweight = htweight_unnorm/sum(htweight_unnorm)
  } else if(weight_model == "CBPS"){

    form_max = stats::as.formula(paste0("factor(biased) ~",paste0(c(fact_vars, paste0("poly(",cont_vars,",2)")),collapse = "+")))
    fit_cbps = CBPS::CBPS(formula = form_max, data = data, ATT = 2) # Representative sample should be the 'treatment' group
    prob_bias = 1 - fit_cbps$fitted.values # want prob c2c (not nhanes)
    htweight_unnorm = (1-prob_bias)/prob_bias
    htweight = htweight_unnorm/sum(htweight_unnorm)
  } else if(weight_model == "entbal"){

    form_max = stats::as.formula(paste0("factor(biased) ~",paste0(cov,collapse = "+")))

    eb_pars_att = list(exp_type = 'binary',
                       estimand = 'ATT',
                       n_moments = 3,
                       optim_method = 'L-BFGS-B',
                       verbose = T,
                       opt_constraints = c(-250,250),
                       bal_tol = 1e-8,
                       max_iters = 1000,
                       which_z = 0) # this will match covariates to those in NHANES-REP
    estwts_eb_att = entbal::entbal(form_max, data = data,
                                   eb_pars = eb_pars_att)

    htweight = estwts_eb_att$wts
  } else{
    stop("No sampling weight estimation model selected!")
  }

  ## create dataframe for outcome model
  weights = data.frame(ID = data$ID,htweight)
  if(!incl_wts_for_rep){
    weights = weights[data$biased == 1,]
  }

  return(weights)

}
