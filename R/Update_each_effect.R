# @title update each effect once
# @param X an n by p matrix of regressor variables
# @param Y an n vector of response variable
# @param s a SuSiE fit
# @param estimate_prior_variance boolean indicating whether to
#   estimate prior variance
# @param check_null_threshold float a threshold on the log scale to
#   compare likelihood between current estimate and zero the null
update_each_effect = function (X, Y, s,D,U,UtX_square, UtXD1,UtXD1_square, X_square,estimate_prior_variance = FALSE,check_null_threshold)
  {
  # Repeat for each effect to update.
  L = nrow(s$alpha)
  if (L > 0)
    for (l in 1:L) {
      
      # Remove lth effect from fitted values.
      s$Xr = s$Xr - tcrossprod(X,t(s$alpha[l,] * s$mu[l,]))
      
      # Compute residuals.
      R = Y - s$Xr
      
      res = single_effect_regression(R,X,U,D, UtXD1,UtXD1_square,s$sigma_b[l],s$sigma_k,s$sigma2,s$pi,check_null_threshold)
      
      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$sigma_b[l]      = res$sigma_b
      s$lbf[l]    = res$lbf_model
      s$lbf_variable[l,] = res$lbf
      s$KL[l]     = -res$loglik + SER_posterior_e_loglik(X,R,s$sigma2,s$sigma_k,res$alpha * res$mu,res$alpha * res$mu2,D,U,UtX_square,X_square)
      
      s$Xr = s$Xr + tcrossprod(X,t(s$alpha[l,] * s$mu[l,]))
    }
  return(s)
}


