##### Without GEO ######

geo_center <- function(L_chol, alphas,solver = c("MOSEK","ECOS","SCS")) {
  A <- as.matrix(alphas); K <- ncol(A); J <- nrow(A)
  L <- as.matrix(L_chol)
  gamma <- Variable(K); tvar <- Variable(1)
  constr <- lapply(1:J, function(j) {
    # || L (gamma - alpha_j) ||_2 <= t
    norm( L %*% (gamma - A[j,]), "2") <= tvar
  })
  prob <- Problem(Minimize(tvar), constr)
  sol <- tryCatch(solve(prob, solver = match.arg(solver)),
                  error = function(e) solve(prob)) 
  #sol = solve(prob, solver = "MOSEK") #"ECOS") 
  list(gamma = as.vector(sol$getValue(gamma)),
       worst_value = as.numeric(sol$value),
       status = sol$status)
}

# fit maximin with given scaled delta
fit_maximin_noGEO <- function(scale.delta, fit.prep, est.gammaGEO, Cov_estZ, test.data,solver = c("MOSEK","ECOS","SCS")) {
  # scale.delta represents the scaled delta w.r.t delta.null
  # note that loss() is defined inside of this function, and your X, y should not has missing value(because you loss_logistic() didn't consider missing value), o.w., CVXR would have error
  est.gammaGEO = matrix(est.gammaGEO,ncol=1)
  L <- chol(Cov_estZ)
  n0 = nrow(test.data)
  nz = ncol(test.data$post_Z)
  
  
  tGamma = t(fit.prep$pred.mat) %*% fit.prep$pred.mat / n0
  v = Variable(nz, name = "v")
  obj = quad_form(v, tGamma)
  constraints <- list(norm(L %*% (v - est.gammaGEO), "2") <= scale.delta)  
  prob = Problem(Minimize(obj), constraints)
  result = tryCatch(solve(prob, solver = solver),
                    error = function(e) solve(prob))
  v.opt = result$getValue(v)
  
  tau = fit.prep$fit.sep.coef %*% v.opt
  pred.te = as.vector(fit.prep$pred.mat %*% v.opt)
  
  return(list(tau = tau, v.opt = v.opt, pred.te = pred.te))
} # note test.data contains response and post_Z




kl_fun <- function(a,b){
  b*(log(b)-log(a)) + (1-b) * (log(1-b)-log(1-a))
}


kl_mean <- function(test.data, a, b, b0=0) {
  eta1 <- 1/(1+exp(-test.data$post_Z %*% a + b0))
  eta2 <- 1/(1+exp(-test.data$post_Z %*% b + b0))
  kl_vec <- kl_fun(eta1,eta2) 
  return(mean(kl_vec))
}


###### New fit_maximin #####
select_delta_new_noGEO <- function(Z, est_O, Y,  est_gamma, quan  = 1, exclude_null = F, upper_b = 0.9){
  # Z is est.Z, O and Y is observed. different from generate_gamma, here we used only observed data, instead of population level empirical expectation
  nz = ncol(Z)
  O_glm = glm(est_O ~ Z, family = "binomial")
  logLik_mle = as.numeric(logLik(O_glm))
  logLik_null = as.numeric(logLik(glm(est_O ~ 1, family = "binomial")))
  
  #logLik_mle = loss_est2(est_gamma,b0 = 0, X = Z, y = est_O)
  #logLik_null =loss_est2(rep(0,nz),b0 = 0, X = Z, y = est_O)
  
  delta.null = logLik_mle - logLik_null
   # inside this function sample size not change so no over n.sub
  # binary outcome
  get_delta_one <- function(y){
    outcome = as.numeric(y == 2) # this prob is generated from g^(-1)(-simul.alpha%*%true.Z), but estimated with g^(-1)(est.alpha%*%est.Z), IRT give estimation opposite with true, so est.Z and true.Z are on same direction
    glm1 = glm(outcome ~ Z, family = "binomial") # use est.Z, then est.gamma corresponding to est.alpha, right combo for external outcomes
    linear_pre = Z%*%glm1$coefficients[1:nz+1]
    fit <- glm(est_O ~ 1 + offset(linear_pre), family = "binomial")
    # return(tail(coef(glm1), nz)) # only return gamma, not intercept= 
    return(as.numeric(logLik(fit)))# return loglik
    #alpha_y = as.vector(glm1$coefficients[1:nz+1])
    #logLik_y = loss_est2(alpha_y,b0 = 0, X = Z, y = est_O)
    #return(logLik_y)# return loglik
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

  
fit_maximin_new_noGEO <- function(scale.delta, fit.prep, est_GEO, est_O, test.data,solver = c("MOSEK","ECOS","SCS")) {
  n0 = nrow(test.data)
  nz = ncol(test.data$post_Z)
  # scale.delta represents the scaled delta w.r.t delta.null
  new.data = tibble(response = est_O, post_Z = test.data$post_Z)
  res_glm <- glm(response ~ post_Z, data = new.data, family = "binomial") 
  logLik_mle = as.numeric(logLik(res_glm))
  
  logLik_null = as.numeric(logLik(glm(response ~ 1, data = new.data, family = "binomial")))
  delta = scale.delta * (logLik_mle - logLik_null)
  
  # note that loss() is defined inside of this function, and your X, y should not has missing value(because you loss_logistic() didn't consider missing value), o.w., CVXR would have error
  loss <- function(v){loss_logistic(v, b0 = as.numeric(res_glm$coefficients[1]), X = as.matrix(new.data$post_Z), y = as.numeric(new.data$response))}
  #logLik_mle = loss_est2(est_GEO, b0 = 0, X = as.matrix(test.data$post_Z), y = as.numeric(est_O))
  #logLik_null =loss_est2(rep(0,nz), b0 = 0, X = as.matrix(test.data$post_Z), y = as.numeric(est_O))
  #loss <- function(v){loss_est(v, b0 = 0, X = as.matrix(data_est$post_Z), y = as.numeric(est_O))}
  
  
  
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





##### sampling from MVN to S^{p-1} #####
r_acg <- function(n, Sigma) {
  p <- nrow(Sigma)
  if (ncol(Sigma) != p) stop("Sigma must be p x p.")
  # chol gives upper-triangular R so that t(R) %*% R = Sigma
  R <- chol(Sigma)
  Z <- matrix(rnorm(n * p), n, p)        # i.i.d. N(0,1)
  Y <- Z %*% R                           # N(0, Sigma)
  X <- Y / sqrt(rowSums(Y^2))            # radial projection to unit sphere
  return(X)                                       # n x p, rows on S^{p-1}
}



### vMF sampling for gamma_Y's  ###

sample_t <- function(n, p, kappa) {
  if (kappa <= 0) {              # uniform case
    # t has density ∝ (1 - t^2)^{ν - 1/2}; draw via Beta, then symmetrize
    nu <- p/2 - 1
    u  <- rbeta(n, 0.5, nu + 0.5)
    return(2*u - 1)
  }
  b  <- (-2 * kappa + sqrt(4 * kappa^2 + (p - 1)^2)) / (p - 1)
  x0 <- (1 - b) / (1 + b)
  c  <- kappa * x0 + (p - 1) * log(1 - x0^2)
  
  tvals <- numeric(n)
  a <- (p - 1) / 2
  for (i in seq_len(n)) {
    repeat {
      z <- rbeta(1, a, a)
      w <- (1 - (1 + b) * z) / (1 - (1 - b) * z)  # candidate t
      if (kappa * w + (p - 1) * log(1 - x0 * w) - c >= log(runif(1))) {
        tvals[i] <- w
        break
      }
    }
  }
  return(tvals)
}

# Draw v uniformly on the (p-2)-sphere in the tangent space of mu
tangent_unit_vectors <- function(n, mu) {
  mu <- mu / sqrt(sum(mu^2))
  p  <- length(mu)
  G  <- matrix(rnorm(n * p), n, p)
  # project each row onto mu^\perp
  proj <- as.numeric(G %*% mu)
  V <- G - proj %o% mu
  return(V / sqrt(rowSums(V^2)))
}

# Main sampler: n draws from vMF(mu, kappa) using x = t mu + sqrt(1 - t^2) v
rvmf_rt <- function(n, mu, kappa) {
  mu <- as.numeric(mu)
  p  <- length(mu)
  if (p < 2) stop("Dimension p must be >= 2.")
  if (all(mu == 0)) stop("mu must be nonzero.")
  mu <- mu / sqrt(sum(mu^2))
  
  t  <- sample_t(n, p, kappa)               # radial component in [-1, 1]
  V  <- tangent_unit_vectors(n, mu)         # unit vectors in mu^\perp
  X  <- t %o% mu + (sqrt(1 - t^2) %o% rep(1, p)) * V
  return(X)
}


# 1. function to simulate data ------------

# g^(-1) 
## in case of Inf appears, use this for g^(-1)
inverse_link <- function(log_odds) {
  if (is.finite(log_odds)) # cause add this line, this function can only apply to single element. if want to apply to vector/matrix, use apply()
    return(1 / (1 + exp(-log_odds))) # if log_odds = Inf, this would be 1; 
  else if (log_odds == Inf){
    return(1)
  } else if (log_odds == -Inf){
    return(0)
  }
}

# f(x)
f_w <- function(w, Lambda) { # Lambda: nw*nz
  return(w %*% Lambda) # return to the row vector of E[Z}: 1*nz
}

# modify Y ~ alpha%*%Z, no negative sign later, if so remember to modify get_delta, cause that used Y
generate_FAS_0930 <- function(nsub, ny, cat.y, nz, nx, Lambda, tau, alpha, gamma_GEO, sigma_GEO = 1, centerZ = T, ny_cont = NULL, sigma = NULL, alpha_cont = NULL, beta_cont = NULL){
  # Lambda should be Lambda 0 first, then by Lambda1, note the generated X will automatically last col = 1 as intercept
  # gamma_GEO should be of length nz+1, the last one is for intercept
  A <- rbinom(nsub, 1, 0.5) # treatment
  X <- cbind(2*matrix(rnorm((nx - 1)*nsub), nrow = nsub, ncol = nx - 1), 1)  # covariate on z, nsub*nw; last column is for intercept!!!
  lambda_indices = A + 1
  Z = t(sapply(1:nsub, function(ix) {f_w(X[ix,], Lambda[[lambda_indices[ix]]])})) + matrix(rnorm(nsub*nz), nrow = nsub, ncol = nz) # control random error scale
  # Z.test = matrix(NA, nrow = nsub, ncol = nz) # for test
  # for (ix in 1:nsub) {
  #   Z.test[ix,] <- f_w(W[ix,], Lambda[[A[ix] + 1]]) #+ rnorm(nz) # W should be in a larger scale than N(0,1) random error, o.w, Lambda estimation not good
  # }
  
  if (centerZ == T){
    Z <- scale(Z, center = T, scale = F) # default center Z before generating O
  }
  # generate binary GEO
  logits <- cbind(Z, 1) %*% gamma_GEO 
  probabilities <- 1 / (1 + exp(-logits))
  O <- rbinom(n = nsub, size = 1, prob = probabilities) # primary outcome 
  
  
  # # generate continuous GEO  
  # generate_oneO <- function(z){
  #   o = stats::rnorm(1, mean = gamma_GEO[nz + 1] + gamma_GEO[1:nz] %*% z, sd = sigma_GEO)
  # }
  # O <- apply(Z, 1, generate_oneO)
  
  # calculate prob to generate y, for one subject
  generate_oneY <- function(z){
    Gamma = apply(sweep(tau, MARGIN = 1, STATS = -alpha %*% z, FUN = "+"), c(1,2), inverse_link)
    Pi = Gamma - cbind(rep(0, length(cat.y)), Gamma[,1:(max(cat.y) - 1)])
    y = apply(Pi, 1, function(row) sample(x = 1:length(row), size = 1, prob = row))
  }
  generate_oneY_cont <- function(z){ # note later we will use mapply, so input z is a vector, not a matrix; so the mean formula coding will be a little different from g_function_cont()
    y_cont = stats::rnorm(ny_cont, mean = as.numeric(-alpha_cont %*% z ), sd = sigma)
  }
  
  if (is.null(ny_cont)){
    # Y <- t(mapply(generate_oneY, z = split(Z, row(Z)), x = split(X, row(X)))) # post treatment item response, nsub*ny
    Y = t(apply(Z, 1, generate_oneY))
    return(list(X = X, Z = Z, Y = Y, A = A, O = O))
  } else {
    # Y <- t(mapply(generate_oneY, z = split(Z, row(Z)), x = split(X, row(X)))) # post treatment item response, nsub*ny
    # Y_cont <- t(mapply(generate_oneY_cont, z = split(Z, row(Z)), x = split(X, row(X)))) # cluster has error?? Error in stats::rnorm(ny_cont, mean = as.numeric(-alpha_cont %*% z + beta_cont %*%  : invalid arguments
    Y = t(apply(Z, 1, generate_oneY))
    Y_cont = t(apply(Z, 1, generate_oneY_cont))
    return(list(X = X, Z = Z, Y = Y, A = A, Y_cont = Y_cont, O = O))
  }
}

# modified based on 240701maximin_function.R, but refer to new simulated data generation frame
loss_logistic <- function(b, b0, X, y){
  # not loss, it's loglik, the larger the better
  l1 = -sum(X[y <= 0, ] %*% b + b0) - sum(CVXR::logistic(-X %*% b - b0)) # must use CVXR version logistic()  !!!, o.w. it cannot identify convex
  # l1 = -sum(X[y == 0, ] %*% b + b0) - sum(log(1 + exp(-X %*% b - b0))) # for testing ML  
  return(l1)
}

# calculate empirical Xi, and TE coef for each Z
fit_prep <- function(x, Lambda1, Lambda0) { #add Lambda1, Lambda0 as arguments input, since we are not using get_TE() for step 1 # fit model for individual factor
  # x default has intercept term
  fit.sep.coef = Lambda1 - Lambda0 # nx*nz
  pred.mat = x %*% fit.sep.coef # nsub*nz
  return(list(pred.mat = pred.mat, fit.sep.coef = fit.sep.coef))
}

# fit maximin with given scaled delta
fit_maximin <- function(scale.delta, fit.prep, test.data) {
  # scale.delta represents the scaled delta w.r.t delta.null
  res_glm <- glm(response ~ post_Z, data = test.data, family = "binomial") 
  logLik_mle = as.numeric(logLik(res_glm))
  
  logLik_null = as.numeric(logLik(glm(response ~ 1, data = test.data, family = "binomial")))
  delta = scale.delta * (logLik_mle - logLik_null)
  
  # note that loss() is defined inside of this function, and your X, y should not has missing value(because you loss_logistic() didn't consider missing value), o.w., CVXR would have error
  loss <- function(v){loss_logistic(v, b0 = as.numeric(res_glm$coefficients[1]), X = as.matrix(test.data$post_Z), y = as.numeric(test.data$response))}
  n0 = nrow(test.data)
  nz = ncol(test.data$post_Z)
  
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

# currently can use true best delta, and leave the part for self selected delta later

# fit baseline method with S
fit_baseline <- function(test.data, fit.new = FALSE, new.data = NULL){ # note test.data contains S and A and X     
  test.data1 = test.data %>% mutate(A = 1) 
  test.data0 = test.data %>% mutate(A = 0) 
  test.lm1 = glm(response ~ A + X + A:X, data = test.data, family = "binomial")
  if (fit.new == FALSE){
    TE = predict(test.lm1, newdata = test.data1) - predict(test.lm1, newdata = test.data0) # >0 means A=1 more likely to have better S response
  }else{
    new.data1 = new.data %>% mutate(A = 1)
    new.data0 = new.data %>% mutate(A = 0)
    TE = predict(test.lm1, newdata = new.data1) - predict(test.lm1, newdata = new.data0)
  }
  return(list(TE=TE, mod = test.lm1))
}

# not used in new simulation setting
generate_U <- function(Z, gamma_mat){ # modify min_r later to allow lower distance for gamma
  # binary test outcomes
  generate_oneU <- function(gamma, Z){# generate one external outcome for all subjects
    probabilities <- sapply(Z %*% gamma + 0, inverse_link)
    U = rbinom(nsub, 1, probabilities)
    return(U)
  }
  U <- apply(gamma_mat, 1, generate_oneU, Z = Z)
  
  return(list(Z = Z, gamma_mat = gamma_mat, U = U))
}

# get reference treatment effect for outcome U for evaluation purpose
ref_TE <- function(U, data){ # exchange the U to the first argument (compared to 240620) for ease of apply
  tmp.data =  tibble(U = U, X = data$X, A = data$A) # data contains X and A
  tmp.data1 = tmp.data %>% mutate(A = 1) 
  tmp.data0 = tmp.data %>% mutate(A = 0) 
  test.lm1 = glm(U ~ A + X + A:X, data = tmp.data, family = "binomial")
  TE = predict(test.lm1, newdata = tmp.data1) - predict(test.lm1, newdata = tmp.data0) # >0 means A=1 more likely to have better S response
  return(TE)
}
# get reference treatment effect parameters for outcome U with true parameter
ref_TE_true <- function(gamma_one, Lambda1, Lambda0){
  return(as.vector((Lambda1-Lambda0) %*% gamma_one)) # true lambda and true gamma, corresponding to true TE on factorized outcome
}


# for use of calculate loss
loss_logistic2 <- function(b, b0, X, y){ # for normally calculate loglik loss
  # modify to avoid Inf then NAN
  linear_pred = X %*% b + b0
  prob = sapply(linear_pred, inverse_link)
  l1 = sum(ifelse(y == 1, log(prob), log(1 - prob)))
  return(l1)
} # not really loss, just loglikelihood # for calculating in general to debug. (log_logistic is for CVXR)

# selected delta based on Y, not tested
select_delta <- function(Z, O, Y, quan  = 1){
  # Z is est.Z, O and Y is observed. different from generate_gamma, here we used only observed data, instead of population level empirical expectation
  nz = ncol(Z)
  O_glm = glm(O ~ Z, family = "binomial")
  logLik_mle = as.numeric(logLik(O_glm ))
  logLik_null = as.numeric(logLik(glm(O ~ 1, family = "binomial")))
  
  delta.null = logLik_mle - logLik_null # inside this function sample size not change so no over n.sub
  # binary outcome
  get_delta_one <- function(y){
    outcome = as.numeric(y == 1) # this prob is generated from g^(-1)(-simul.alpha%*%true.Z), but estimated with g^(-1)(est.alpha%*%est.Z), IRT give estimation opposite with true, so est.Z and true.Z are on same direction
    glm1 = glm(outcome ~ Z, family = "binomial") # use est.Z, then est.gamma corresponding to est.alpha, right combo for external outcomes
    # return(tail(coef(glm1), nz)) # only return gamma, not intercept
    return(as.numeric(logLik(glm1)) )# return loglik
  }
  # gamma.Y = -t(apply(Y, 2, get_delta_one)) # take negative,since fitted GLM use as.numeric(y == 1), and y take 1 or 2 (generated through graded response model)
  # loss2 <- function(v){loss_logistic2(v, b0 = as.numeric(O_glm$coefficients[1]), X = as.matrix(Z), y = O)}
  # logLik_obs.outcome = apply(gamma.Y, 1, function(gamma){loss2(gamma)})
  logLik_obs.outcome = apply(Y, 2, get_delta_one)
  deviance_obs.outcome = (logLik_mle - logLik_obs.outcome)/delta.null
  return(quantile(deviance_obs.outcome, probs = quan, na.rm = TRUE))
}


select_delta_new <- function(Z, O, Y,  quan  = 1, exclude_null = F, upper_b = 0.9){
  # Z is est.Z, O and Y is observed. different from generate_gamma, here we used only observed data, instead of population level empirical expectation
  nz = ncol(Z)
  O_glm = glm(O ~ Z, family = "binomial")
  logLik_mle = as.numeric(logLik(O_glm ))
  logLik_null = as.numeric(logLik(glm(O ~ 1, family = "binomial")))
  
  delta.null = logLik_mle - logLik_null # inside this function sample size not change so no over n.sub
  # binary outcome
  get_delta_one <- function(y){
    outcome = as.numeric(y == 2) # this prob is generated from g^(-1)(-simul.alpha%*%true.Z), but estimated with g^(-1)(est.alpha%*%est.Z), IRT give estimation opposite with true, so est.Z and true.Z are on same direction
    glm1 = glm(outcome ~ Z, family = "binomial") # use est.Z, then est.gamma corresponding to est.alpha, right combo for external outcomes
    linear_pre = Z%*%glm1$coefficients[1:nz+1]
    fit <- glm(O ~ 1 + offset(linear_pre), family = "binomial")
    # return(tail(coef(glm1), nz)) # only return gamma, not intercept= 
    return(as.numeric(logLik(fit)))# return loglik
  }
  # gamma.Y = -t(apply(Y, 2, get_delta_one)) # take negative,since fitted GLM use as.numeric(y == 1), and y take 1 or 2 (generated through graded response model)
  # loss2 <- function(v){loss_logistic2(v, b0 = as.numeric(O_glm$coefficients[1]), X = as.matrix(Z), y = O)}
  # logLik_obs.outcome = apply(gamma.Y, 1, function(gamma){loss2(gamma)})
  logLik_obs.outcome = apply(Y, 2, get_delta_one)
  deviance_obs.outcome = (logLik_mle - logLik_obs.outcome)/delta.null
  if (exclude_null){
    delta_select = min(quantile(deviance_obs.outcome, probs = quan, na.rm = TRUE),upper_b)
  } else{
    delta_select = quantile(deviance_obs.outcome, probs = quan, na.rm = TRUE)
  }
  return(delta_select)
}

select_delta_seq <- function(Z, O, Y,  delta_list, quan  = 1, exclude_null = F, upper_b = 0.9){
  # Z is est.Z, O and Y is observed. different from generate_gamma, here we used only observed data, instead of population level empirical expectation
  nz = ncol(Z)
  O_glm = glm(O ~ Z, family = "binomial")
  logLik_mle = as.numeric(logLik(O_glm ))
  logLik_null = as.numeric(logLik(glm(O ~ 1, family = "binomial")))
  
  delta.null = logLik_mle - logLik_null # inside this function sample size not change so no over n.sub
  # binary outcome
  get_delta_one <- function(y){
    outcome = as.numeric(y == 2) # this prob is generated from g^(-1)(-simul.alpha%*%true.Z), but estimated with g^(-1)(est.alpha%*%est.Z), IRT give estimation opposite with true, so est.Z and true.Z are on same direction
    glm1 = glm(outcome ~ Z, family = "binomial") # use est.Z, then est.gamma corresponding to est.alpha, right combo for external outcomes
    linear_pre = Z%*%glm1$coefficients[1:nz+1]
    fit <- glm(O ~ 1 + offset(linear_pre), family = "binomial")
    # return(tail(coef(glm1), nz)) # only return gamma, not intercept= 
    return(as.numeric(logLik(fit)))# return loglik
  }
  # gamma.Y = -t(apply(Y, 2, get_delta_one)) # take negative,since fitted GLM use as.numeric(y == 1), and y take 1 or 2 (generated through graded response model)
  # loss2 <- function(v){loss_logistic2(v, b0 = as.numeric(O_glm$coefficients[1]), X = as.matrix(Z), y = O)}
  # logLik_obs.outcome = apply(gamma.Y, 1, function(gamma){loss2(gamma)})
  logLik_obs.outcome = apply(Y, 2, get_delta_one)
  deviance_obs.outcome = (logLik_mle - logLik_obs.outcome)/delta.null
  delta_select = quantile(deviance_obs.outcome, probs = quan, na.rm = TRUE)
  if (exclude_null){
    opt = sapply(delta_list*delta_select, function(x)min(x,upper_b))
  } else{
    opt = delta_list*delta_select
  }
  return(opt)
}



select_delta_record <- function(Z, O, Y,  quan  = 1, upper_b = 0.9){
  # Z is est.Z, O and Y is observed. different from generate_gamma, here we used only observed data, instead of population level empirical expectation
  nz = ncol(Z)
  O_glm = glm(O ~ Z, family = "binomial")
  logLik_mle = as.numeric(logLik(O_glm ))
  logLik_null = as.numeric(logLik(glm(O ~ 1, family = "binomial")))
  
  delta.null = logLik_mle - logLik_null # inside this function sample size not change so no over n.sub
  # binary outcome
  get_delta_one <- function(y){
    outcome = as.numeric(y == 2) # this prob is generated from g^(-1)(-simul.alpha%*%true.Z), but estimated with g^(-1)(est.alpha%*%est.Z), IRT give estimation opposite with true, so est.Z and true.Z are on same direction
    glm1 = glm(outcome ~ Z, family = "binomial") # use est.Z, then est.gamma corresponding to est.alpha, right combo for external outcomes
    linear_pre = Z%*%glm1$coefficients[1:nz+1]
    fit <- glm(O ~ 1 + offset(linear_pre), family = "binomial")
    # return(tail(coef(glm1), nz)) # only return gamma, not intercept= 
    return(as.numeric(logLik(fit)))# return loglik
  }
  # gamma.Y = -t(apply(Y, 2, get_delta_one)) # take negative,since fitted GLM use as.numeric(y == 1), and y take 1 or 2 (generated through graded response model)
  # loss2 <- function(v){loss_logistic2(v, b0 = as.numeric(O_glm$coefficients[1]), X = as.matrix(Z), y = O)}
  # logLik_obs.outcome = apply(gamma.Y, 1, function(gamma){loss2(gamma)})
  logLik_obs.outcome = apply(Y, 2, get_delta_one)
  deviance_obs.outcome = (logLik_mle - logLik_obs.outcome)
  return(list(delta_select = quantile(deviance_obs.outcome, probs = quan, na.rm = TRUE),
              delta_null = delta.null))
}



getLambda <- function(nx,nz,beta){
  lam_0 = c()
  lam_1 = c()
  for (idz in 1:nz){
    begin_num=(idz-1)*2*nx+1
    lam_0=cbind(lam_0, beta[c((begin_num+1):(begin_num+nx-1),begin_num)])
    lam_1=cbind(lam_1, beta[c((begin_num+1+nx):(begin_num+2*nx-1),begin_num+nx)]+beta[c((begin_num+1):(begin_num+nx-1),begin_num)])
  }
  return(list(lam_0,lam_1))
}
