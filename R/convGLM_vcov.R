#' Variance covariance matrix for propensity-weighted generalized linear models with logistic weight estimation
#'
#' Estimate covariance matrix with a simultaneous estimating equation approach for a propensity-weighted generalized linear model
#' when weights are estimated via logistic regression
#'
#' @param estwt_fit A \code{glm} object from the propensity weight estimation
#' @param fit_outcome A \code{glm} object from the final outcome model
#' @param biased An indicator of biased sample membership
#'
#' @return Covariance matrix estimate that accounts for uncertainty in propensity-weight estimation
#' for model coefficients from \code{fit_outcome}
#' @export
#'
convGLM_vcov = function(estwt_fit, fit_outcome, biased){
  ## I_TT
  I_TT = solve(summary(estwt_fit)$cov.unscaled)

  ## I_UU
  mu = stats::fitted(fit_outcome)
  V_mu = (stats::residuals(fit_outcome,"response")/stats::residuals(fit_outcome,"pearson"))^2 #
  dmu_deta = stats::residuals(fit_outcome,"response")/stats::residuals(fit_outcome,"working")
  X = stats::model.matrix(fit_outcome)
  pi = 1/fit_outcome$survey.design$prob
  pi = pi/sum(pi)*length(pi)
  I_UU = t(X)%*%(pi*diag(1/V_mu*(dmu_deta)^{2})%*%X)

  ##  I_UT
  # subset W down to subjects from biased sample
  W = stats::model.matrix(estwt_fit)[as.logical(biased),]
  y_mu = stats::residuals(fit_outcome, type = 'response')
  I_UT = t(X)%*%(pi*diag(dmu_deta*y_mu/V_mu)%*%W)

  ## T
  if (is.matrix(estwt_fit$x)){
    xmat=estwt_fit$x
  } else {
    mf=stats::model.frame(estwt_fit)
    xmat=stats::model.matrix(stats::terms(estwt_fit),mf)
  }
  T_score_all = stats::residuals(estwt_fit,"working")*estwt_fit$weights*xmat
  # subset to subjects from biased sample
  T_score = T_score_all[as.logical(biased),]

  ## U
  if (is.matrix(fit_outcome$x)){
    xmat=fit_outcome$x
  } else {
    mf=stats::model.frame(fit_outcome)
    xmat=stats::model.matrix(stats::terms(fit_outcome),mf)
  }
  U_score = (stats::residuals(fit_outcome,"working")*fit_outcome$weights*xmat)


  # var/cov matrix
  vcov <- solve(I_UU)%*%t(U_score)%*%(U_score)%*%solve(I_UU) -
                    solve(I_UU)%*%I_UT%*%solve(I_TT)%*%(t(T_score))%*%((U_score))%*%solve(I_UU)

  return(vcov)
}
