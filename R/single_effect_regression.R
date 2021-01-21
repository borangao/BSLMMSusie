single_effect_regression =function(Y, X, U, D,UtXD1,UtXD1_square,sigma_b,sigma_k,sigma2 = 1, prior_weights = NULL,check_null_threshold = 0) {
    Xty = crossprod(X,Y)
    Uty = crossprod(U,Y)
    n = length(Y)
    XtHy=Xty/sigma2-tcrossprod(UtXD1,t(Uty))
    S_l<-1/(1/sigma_b+(n-1)/sigma2-UtXD1_square)
    sigma2_k = (n-1)/sigma2-UtXD1_square
    ####First Update sigma_b^2 Use EB
    if (is.null(prior_weights))
    prior_weights = rep(1/ncol(X),ncol(X))
    
    sigma_b = optimize_prior_variance(XtHy,sigma_b,sigma2_k,prior_weights,
                                alpha = NULL,post_mean2 = NULL,
                                check_null_threshold = check_null_threshold)
    sigma_b
    ##Given Updated Sigma_B, update alpha and mu mu^2
    S_l<-1/(1/sigma_b+sigma2_k)
    lbf = 0.5*(XtHy*S_l*XtHy)+0.5*log(S_l)-0.5*log(sigma_b)
    lbf[is.infinite(S_l)] = 0 
    maxlbf = max(lbf)
    w = exp(lbf - maxlbf) # w = BF/BFmax
    w_weighted = w * prior_weights
    weighted_sum_w = sum(w_weighted)
    alpha = w_weighted / weighted_sum_w
    
    
    post_var = S_l # Posterior variance.
    post_mean = S_l * XtHy
    post_mean2 = post_var + post_mean^2 # Second moment.
    
    ###BF for single effect model
    lbf_model = maxlbf + log(weighted_sum_w)
    loglik = lbf_model + ((-log(2)-log(pi)-log(sigma2))*n-sum(log(1+sigma_k*D))-1/sigma2*sum(Y^2)+1/sigma2*sum(Uty^2*sigma_k*D/(1+sigma_k*D)))/2
       return(list(alpha = alpha,mu = post_mean,mu2 = post_mean2,lbf = lbf,
                lbf_model = lbf_model,sigma_b = sigma_b,loglik = loglik))
  }
loglik = function (sigma_b, XtHy, sigma2_k, prior_weights) {
  S_l<-1/(1/sigma_b+sigma2_k)
 # S_l = sigma_b/(sigma2_k*sigma_b+1)
  lbf = 0.5*(XtHy*S_l*XtHy)+0.5*(log(S_l)-log(sigma_b))
  lbf[is.infinite(S_l)] = 0
  maxlbf = max(lbf)
  w = exp(lbf - maxlbf) # w = BF/BFmax
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w
  return(log(weighted_sum_w) + maxlbf)
}

neg.loglik.logscale = function(log_sigma_b, XtHy, sigma2_k,prior_weights)
  -loglik(exp(log_sigma_b),XtHy, sigma2_k,prior_weights)

optimize_prior_variance = function (XtHy, sigma_b,sigma2_k, prior_weights,
                                    alpha = NULL, post_mean2 = NULL,
                                    V_init = NULL, check_null_threshold = 0) {
      V = V_init
      lV = optim(par = log(max(c((XtHy*sigma2_k)^2-sigma2_k,1),na.rm = TRUE)),
                 fn = neg.loglik.logscale,XtHy = XtHy,sigma2_k = sigma2_k,
                 prior_weights = prior_weights,method = "Brent",lower = -30,
                 upper = 15)$par
      V = exp(lV)
   #  if (loglik(V,XtHy,sigma2_k,prior_weights)<=0)
   #   V = 0
      return(V)

  }
  




