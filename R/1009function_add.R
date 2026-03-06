 # Requires: CVXR
# Assumes: test.data$post_Z is an n x p numeric matrix (no NAs)
#          test.data$Y      is an n x J matrix/data.frame of 0/1 outcomes (no NAs)
#          fit.prep$pred.mat is n x p; fit.prep$fit.sep.coef is compatible with p
# Note: Like your original, we *fix* the intercept at the MLE of the aggregated model.

group_DRO <- function(fit.prep, test.data, est.alpha, solver = c("ECOS", "SCS")) {
  solver <- match.arg(solver)

  Z   <- as.matrix(test.data$post_Z)            
  n0  <- nrow(Z)
  p   <- ncol(Z)
  J   <- nrow(est.alpha)
  B = t(est.alpha)
  P <- as.matrix(fit.prep$pred.mat)    
  
  tGamma = t(fit.prep$pred.mat%*%B) %*% (fit.prep$pred.mat %*%B) / n0
  a <- Variable(J, name = "a")          
  obj <- quad_form(a, tGamma)
  constraints <- list(a >= 0,sum_entries(a) == 1)

  prob   <- Problem(Minimize(obj), constraints)
  result <- tryCatch(solve(prob, solver = solver),
                     error = function(e) {
                       # fall back to the other solver if the first fails
                       alt <- setdiff(c("ECOS","SCS"), solver)[1]
                       message("Primary solver failed (", solver, "), retrying with ", alt, "...\n", e$message)
                       solve(prob, solver = alt)
                     })

  a.opt  <- as.numeric(result$getValue(a))
  v.opt  <- as.numeric(B %*% a.opt)

  # ---------- (5) Same post-processing as your original ----------
  tau     <- as.vector(fit.prep$fit.sep.coef %*% v.opt)
  pred.te <- as.vector(P %*% v.opt)

  list(
    tau        = tau,
    v.opt      = v.opt,       # p-vector inside conv hull of {coef(Y_j ~ Z)}
    a.opt      = a.opt,       # simplex weights (convex-combination coefficients)
    pred.te    = pred.te,
    B          = B,        
    status     = result$status,
    value      = result$value
  )
}















###### New fit_maximin #####
select_delta_new_noGEO_los <- function(Z, est_O, Y,  est_gamma, quan  = 1, exclude_null = F, upper_b = 0.9){
  # Z is est.Z, O and Y is observed. different from generate_gamma, here we used only observed data, instead of population level empirical expectation
  nz = ncol(Z)
  #O_glm = glm(est_O ~ Z, family = "binomial")
  #logLik_mle = as.numeric(logLik(O_glm))
  #logLik_null = as.numeric(logLik(glm(est_O ~ 1, family = "binomial")))
  
  logLik_mle = loss_est2(est_gamma,b0 = 0, X = Z, y = est_O)
  logLik_null =loss_est2(rep(0,nz),b0 = 0, X = Z, y = est_O)
  
  delta.null = logLik_mle - logLik_null
  # inside this function sample size not change so no over n.sub
  # binary outcome
  get_delta_one <- function(y){
    outcome = as.numeric(y == 2) # this prob is generated from g^(-1)(-simul.alpha%*%true.Z), but estimated with g^(-1)(est.alpha%*%est.Z), IRT give estimation opposite with true, so est.Z and true.Z are on same direction
    glm1 = glm(outcome ~ Z, family = "binomial") # use est.Z, then est.gamma corresponding to est.alpha, right combo for external outcomes
    #linear_pre = Z%*%glm1$coefficients[1:nz+1]
    #fit <- glm(est_O ~ 1 + offset(linear_pre), family = "binomial")
    # return(tail(coef(glm1), nz)) # only return gamma, not intercept= 
    #return(as.numeric(logLik(fit)))# return loglik
    alpha_y = as.vector(glm1$coefficients[1:nz+1])
    logLik_y = loss_est2(alpha_y,b0 = 0, X = Z, y = est_O)
    return(logLik_y)# return loglik
  }
  
  logLik_obs.outcome = apply(Y, 2, get_delta_one)
  deviance_obs.outcome = (logLik_mle - logLik_obs.outcome)/delta.null
  
  if (exclude_null){
    delta_select = min(quantile(deviance_obs.outcome, probs = quan, na.rm = TRUE),upper_b)
  } else{
    delta_select = quantile(deviance_obs.outcome, probs = quan, na.rm = TRUE)
  }
  return(delta_select)
}

loss_est <- function(b, b0, X, y){
  # not loss, it's loglik, the larger the better
  l1 = sum((X %*% b+ b0) *y) - sum(CVXR::logistic(X %*% b + b0)) # must use CVXR version logistic()  !!!, o.w. it cannot identify convex
  # l1 = -sum(X %*% b *(1- y) + b0) - sum(log(1 + exp(-X %*% b - b0))) # for testing ML  
  return(l1)
}

softplus_los <- function(z) {
  z <- as.numeric(z)              
  pmax(z, 0) + log1p(exp(-abs(z)))
  
}

loss_est2 <- function(b, b0, X, y){
  eta <- as.vector(X %*% b + b0)
  ll  <- sum(y * eta) - sum(sapply(eta, function(x)softplus_los(x)))
  return(ll)
}

fit_maximin_new_noGEO_los <- function(scale.delta, fit.prep, est_GEO, est_O, test.data,solver = c("MOSEK","ECOS","SCS")) {
  n0 = nrow(test.data)
  nz = ncol(test.data$post_Z)
  # scale.delta represents the scaled delta w.r.t delta.null
  
  #new.data = tibble(response = est_O, post_Z = test.data$post_Z)
  #res_glm <- glm(response ~ post_Z, data = new.data, family = "binomial") 
  #logLik_mle = as.numeric(logLik(res_glm))
  #logLik_null = as.numeric(logLik(glm(response ~ 1, data = new.data, family = "binomial")))
  #loss <- function(v){loss_logistic(v, b0 = as.numeric(res_glm$coefficients[1]), X = as.matrix(new.data$post_Z), y = as.numeric(new.data$response))}
  
  logLik_mle = loss_est2(est_GEO, b0 = 0, X = as.matrix(test.data$post_Z), y = as.numeric(est_O))
  logLik_null =loss_est2(rep(0,nz), b0 = 0, X = as.matrix(test.data$post_Z), y = as.numeric(est_O))
  loss <- function(v){loss_est(v, b0 = 0, X = as.matrix(data_est$post_Z), y = as.numeric(est_O))}
  
  delta = scale.delta * (logLik_mle - logLik_null)
  
  # note that loss() is defined inside of this function, and your X, y should not has missing value(because you loss_logistic() didn't consider missing value), o.w., CVXR would have error
  
  
  delta = scale.delta * (logLik_mle - logLik_null)
  # note that loss() is defined inside of this function, and your X, y should not has missing value(because you loss_logistic() didn't consider missing value), o.w., CVXR would have error
  
  tGamma = t(fit.prep$pred.mat) %*% fit.prep$pred.mat / n0
  v = Variable(nz, name = "v")
  obj = quad_form(v, tGamma)
  constraints = list(loss(v) >= logLik_mle - delta) # (cgi_clm$logLik-delta)# is that because some operation in loss() can not be recognize by CVRX, try add operation step by step?
  prob = Problem(Minimize(obj), constraints)
  result = solve(prob)
  v.opt = result$getValue(v)
  
  tau = fit.prep$fit.sep.coef %*% v.opt
  pred.te = as.vector(fit.prep$pred.mat %*% v.opt)
  return(list(tau = tau, v.opt = v.opt, pred.te = pred.te, logLik_mle = logLik_mle, deviance = logLik_mle - logLik_null))
} # note test.data contains response and post_Z

