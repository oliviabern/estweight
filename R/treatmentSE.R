#' Standard error estimation for treatment effects estimated with convenience samples
#'
#' Estimate standard errors with a simultaneous estimating equation approach for
#' propensity adjusted estimates of a treatment effect when sampling weights are estimated
#' using an auxiliary sample and the sampling weights are used to estimate the propensity score
#' and the treatment effect
#'
#' @param estwt_fit A \code{glm} object from the sampling weight estimation
#' @param estprop_fit A \code{svyglm} object from the propensity score estimation
#' @param fit_outcome A \code{svyglm} object from the treatment effect estimation
#' @param biased An indicator of biased sample membership
#' @param outcome_family Error distribution and link function to be used in the outcome model
#'  (See \code{family} for details of family function). Defaults to stats::gaussian
#'
#' @return Standard error estimates for treatment effect estimates that account for
#' uncertainty in sampling-weight estimation, propensity score estimation, and
#' estimating the treatment effect
#' @export
#'
treatmentSE = function(estwt_fit, estprop_fit, fit_outcome, biased, outcome_family){

  # check that the propensity score model was fit with svyglm
  if(estprop_fit$call[1]!="svyglm()"){
    stop("The propensity score model (estprop_fit) must be fit with svyglm()")
  }
  # check that the outcome model was fit with svyglm
  if(estprop_fit$call[1]!="svyglm()"){
    stop("The outcome model (fit_outcome) must be fit with svyglm()")
  }

  ## Get design matrices
  # sampling weight estimation
  if (is.matrix(estwt_fit$x)){
    X1=estwt_fit$x
  } else {
    mf=stats::model.frame(estwt_fit)
    X1=stats::model.matrix(stats::terms(estwt_fit),mf)
  }
  # propensity estimation
  if (is.matrix(estprop_fit$x)){
    X2=estprop_fit$x
  } else {
    mf=stats::model.frame(estprop_fit)
    X2=stats::model.matrix(stats::terms(estprop_fit),mf)
  }
  # outcome model
  if (is.matrix(fit_outcome$x)){
    X3=fit_outcome$x
  } else {
    mf=stats::model.frame(fit_outcome)
    X3=stats::model.matrix(stats::terms(fit_outcome),mf)
  }

  # get n parameters
  K = ncol(X1) # number of parameters in weight est. model
  L = ncol(X2) # number of parameters in prop. est. model


  ## I_TT
  I_TT = solve(summary(estwt_fit)$cov.unscaled)

  ## I_SS
  e = stats::fitted(estprop_fit)
  w = 1/estprop_fit$survey.design$prob
  w = w/sum(w)*length(w)
  I_SS = t(X2)%*%(w*diag(e*(1-e))%*%X2)

  ## I_UU
  mu = stats::fitted(fit_outcome)
  V_mu = (stats::residuals(fit_outcome,"response")/stats::residuals(fit_outcome,"pearson"))^2 #
  dmudeta = stats::residuals(fit_outcome,"response")/stats::residuals(fit_outcome,"working")
  I_UU = t(X3)%*%(diag(w/V_mu*(dmudeta)^{2})%*%X3)

  # # these match...
  # vcov(fit_outcome)
  # solve(I_UU)%*%t(U_score)%*%U_score%*%solve(I_UU)

  w = w/sum(w) # now we wants weights to sum to 1

  ##  I_ST
  # subset X1 down to subjects from biased sample
  X1.conv = X1[biased==1,]#X1[as.logical(biased),] # updated 1/11/22 OB
  a_e = stats::residuals(estprop_fit, type = 'response')
  I_ST = t(X2)%*%(diag(w*(a_e))%*%X1.conv)

  ## I_US
  beta2 = fit_outcome$coefficients["propscore"]
  y_mu = stats::residuals(fit_outcome, type = "response")

  # Get d/dmu (V(mu))

  if(outcome_family()$family == "gaussian"){
    dVdmu = rep(0, length(mu))
  } else if(outcome_family()$family %in% c("binomial", "quasibinomial")){
    dVdmu = (1-2*mu)
  } else if(outcome_family()$family %in% c("poisson", "quasipoisson")){
    dVdmu = rep(1, length(mu))
  } else{
    stop("outcome_family must be gaussian, binomial, quasibinomial, poisson, or quasipoisson")
  }

  # get d/dmu (dmu/deta)
  if(outcome_family()$link == "identity"){
    ddmu_dmudeta = rep(0, length(mu))
  } else if(outcome_family()$link == "logit"){
    ddmu_dmudeta = (1-2*mu)
  } else if(outcome_family()$link == "log"){
    ddmu_dmudeta = rep(1, length(mu))
  } else{
    stop("link function must be identity, log, or logit")
  }

  one_nc = matrix(1,ncol = 1, nrow = nrow(X3)) # column matrix
  D_US1 = diag(w*e*(1-e)*beta2*(
    -1/V_mu * dmudeta^2
    -y_mu/V_mu^2 * dVdmu * dmudeta^2
    + y_mu/V_mu * ddmu_dmudeta * dmudeta))
  D_US2 = diag(w*e*(1-e)*(y_mu/V_mu*dmudeta))
  I_US = -t(X3)%*%D_US1%*%X2 - matrix(c(0,0,1),ncol = 1)%*%t(one_nc)%*%D_US2%*%X2
  #I_US = matrix(0,nrow = nrow(I_US), ncol = ncol(I_US))

  # bernoulli test
  # test = (t(X3)%*%diag(w*e*(1-e)*beta2*mu*(1-mu))%*%X2 + matrix(c(0,0,1),ncol = 1)%*%t(one_nc)%*%diag(-w*e*(1-e)*y_mu)%*%X2)

  ## I_UT
  D_UT1 = diag(w/V_mu*dmudeta*(y_mu - beta2*a_e*(
    -dmudeta*(1+y_mu/V_mu*dVdmu)
    + y_mu*ddmu_dmudeta
  )))
  D_UT2 = diag(-w/V_mu*dmudeta*y_mu*a_e)
  I_UT = t(X3)%*%D_UT1%*%X1.conv + matrix(c(0,0,1),ncol = 1)%*%t(one_nc)%*%D_UT2%*%X1.conv
  # I_UT = matrix(0,nrow = nrow(I_UT), ncol = ncol(I_UT))

  ## create I
  I = rbind(cbind(I_TT, matrix(0, nrow = K, ncol = L + 3)),
            cbind(I_ST, I_SS, matrix(0, nrow = L, ncol = 3)),
            cbind(I_UT, I_US, I_UU))

  ## T
  T_score_all = stats::residuals(estwt_fit,"working")*estwt_fit$weights*X1
  # subset to subjects from biased sample
  T_score = T_score_all[biased==1,]#[as.logical(biased),] # updated 1/11/22 OB

  ## S
  S_score = stats::residuals(fit_outcome,"working")*fit_outcome$weights*X2

  ## U
  U_score = stats::residuals(fit_outcome,"working")*fit_outcome$weights*X3

  score = cbind(T_score, S_score, U_score)
  scorescoreT = t(score)%*%score
  SR = sqrt(diag(solve(I)%*%scorescoreT%*%solve(I))[(K+L+1):(K+L+3)])
  SR
  return(SR)
}
