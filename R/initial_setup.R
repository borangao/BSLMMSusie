init_setup = function (n, p, L, scaled_prior_variance, residual_variance,kinship_variance,
                       prior_weights, null_weight, varY, standardize) {
  if (!is.numeric(scaled_prior_variance) || scaled_prior_variance < 0)
    stop("Scaled prior variance should be positive number")
  if (scaled_prior_variance > 1 && standardize)
    stop("Scaled prior variance should be no greater than 1 when ",
         "standardize = TRUE")
  if(is.null(residual_variance))
    residual_variance = varY
  if(is.null(prior_weights))
    prior_weights = rep(1/p,p)
  else
    prior_weights = prior_weights / sum(prior_weights)
  if(length(prior_weights) != p)
    stop("Prior weights must have length p")
  if (p < L)
    L = p
  s = list(alpha  = matrix(1/p,nrow = L,ncol = p),
           mu     = matrix(0,nrow = L,ncol = p),
           mu2    = matrix(0,nrow = L,ncol = p),
           Xr     = rep(0,n),
           KL     = rep(as.numeric(NA),L),
           lbf    = rep(as.numeric(NA),L),
           lbf_variable = matrix(as.numeric(NA),L,p),
           sigma2 = residual_variance,
           sigma_b = scaled_prior_variance*varY,
           sigma_k = kinship_variance,
           pi     = prior_weights)
  if (is.null(null_weight))
    s$null_index = 0
  else
    s$null_index = p
  class(s) = "susie"
  return(s)
}

initial_sigma_k<-function(X,Y){
  z = crossprod(X,Y)/sqrt(length(Y))
  ldscore = apply((t(X)%*%X/(length(Y)-1))^2,1,sum)
  sigma_k_initial = mean(z^2-1)/mean(ldscore)
  if(sigma_k_initial < 0 )
    sigma_k_initial = 0
  if(sigma_k_initial>1)
    sigma_k_initial = 1
  return(sigma_k_initial)
  
}