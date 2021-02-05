convGLM = function(data, outcome_formula, response,
                   weight_model = "logistic",
                   outcome_family = stats::gaussian){

  # data:            should be the concatonated combination of the biased sample and representative
  #                  sample with an ID in the first column and an indicator, biased, that is 1
  #                  for subjects in the biased sample in the last column
  #                  factors should be of class factor or binary
  # outcome_formula: and object of class "formula" (or one that can be coerced to that class)
  #                  as in the glm function
  # weight_model:    which sampling weight estimation model should be used. Either "logistic",
  #                  "randomForest", "CBPS", or "entbal"
  # outcome_family:  Error distribution and link function to be used in the outcome model (See
  #                  family for details of family function)

  # Check weight_model selection
  if(!(weight_model %in% c("logistic", "randomForest", "CBPS","entbal"))){
    warning("Invalid weight_model input. Will use default logistic weight estimation model")
    weight_model = "logistic"
  }

  # check data set up
  cols = colnames(data)
  if(sum(cols%in%c("ID","biased"))!=2){
    stop("Data set up incorrectly. Make sure the data includes columns named 'ID' and 'biased'")
  }
  if((sum(colnames(response)%in%c("ID","y"))!=2) | (ncol(response)!=2)){
    stop("Response set up incorrectly. Make sure it includes exactly 2 columns named 'ID' and 'y'")
  }

  # check all subjects in y are in data with biased == 1
  t = data[data$biased==1,"ID"]
  if(sum(response$ID %in% t)!=sum(data$biased)){
    stop("Data or response set up incorrectly. Make sure all subjects in response are also in data with biased == 1.")
  }

  # Identify covariate names and which are factors
  cov = cols[!(cols%in%c("ID","biased"))]
  fs = sapply(data[,cov],is.factor)
  fact_vars = cov[fs]
  cont_vars = cov[!fs]

  # Combine data and response
  data_use = merge(data, response, by = "ID", all = TRUE)

  ## Estimate selection weights

  if(weight_model == "logistic"){

    # Use forward-selection AIC to select a model
    # Model scope includes second order terms
    form_min = stats::as.formula("biased ~ 1")
    form_max = stats::as.formula(paste0("biased ~",paste0(c(cov, paste0("I(",cont_vars,"^2)")),collapse = "+")))
    fit_min = stats::glm(formula = form_min, family = stats::binomial, data = data_use)
    forward = stats::step(fit_min,scope=list(lower=form_min,upper=form_max),
                   direction="forward", trace = 0)

    # fit selected model
    estwt_form = stats::formula(forward)
    estwt_fit = stats::glm(formula = estwt_form, family = stats::binomial, data = data_use)
    prob_bias = stats::fitted(estwt_fit, "response")
    htweight_unnorm = (1-prob_bias)/prob_bias
    htweight = htweight_unnorm/sum(htweight_unnorm)

  } else if(weight_model == "randomForest"){

    form_max = stats::as.formula(paste0("factor(biased) ~",paste0(cov,collapse = "+")))

    fit_rf = randomForest::randomForest(formula = form_max, data = data_use, type = "classification")

    prob_bias1 = stats::predict(fit_rf,type =  "prob")[,"1"]
    prob_bias = ifelse(prob_bias1==0,.01,ifelse(prob_bias1==1,.99,prob_bias1)) # deal with est wts of 0 or 1
    htweight_unnorm = (1-prob_bias)/prob_bias
    htweight = htweight_unnorm/sum(htweight_unnorm)
  } else if(weight_model == "CBPS"){

    form_max = stats::as.formula(paste0("factor(biased) ~",paste0(c(fact_vars, paste0("poly(",cont_vars,",2)")),collapse = "+")))
    fit_cbps = CBPS::CBPS(formula = form_max, data = data_use, ATT = 2) # Representative sample should be the 'treatment' group
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
    estwts_eb_att = entbal::entbal(form_max, data = data_use,
                           eb_pars = eb_pars_att)

    htweight = estwts_eb_att$wts
  } else{
    stop("No sampling weight estimation model selected!")
  }

  ## create dataframe for outcome model
  data_outcome = subset(cbind(data_use,htweight), biased == 1)

  ## Fit outcome model
  fit_outcome = survey::svyglm(outcome_formula,
                       design = survey::svydesign(ids = ~0, weights = data_outcome$htweight, data = data_outcome),
                       family = outcome_family)
  coef = fit_outcome$coefficients

  ## Get analytic variance estimate
  if(weight_model == "logistic"){

    analyticSE = logistic_var(estwt_fit, fit_outcome, biased = data_use$biased)

  } else{
    analyticSE = summary(fit_outcome)$coefficients[,2]
  }

  return(data.frame(coef = coef,
                    analyticSE = analyticSE))


}
