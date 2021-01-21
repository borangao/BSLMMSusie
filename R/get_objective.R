####ELBO function
get_objective = function (X, Y, s,D, U,UtX_square,X_square) {
  return(Eloglik(X,Y,s,D, U,UtX_square,X_square) - sum(s$KL))
}

# Expected loglikelihood for a susie fit.
Eloglik = function (X, Y, s, D, U,UtX_square,X_square) {
  n = nrow(X)
  rl_bar = Y - s$Xr
  rl_bar_U = crossprod(U,rl_bar)
  D_vec_1 = s$sigma_k*D/(1+s$sigma_k*D)
  var_b = colSums((s$alpha * s$mu2)-(s$alpha * s$mu)^2) 
  loglik = -(1/2) *n* (log(2*pi*s$sigma2))-1/2*sum(log(1+D*s$sigma_k))-(sum(rl_bar^2)-sum(rl_bar_U^2*D_vec_1) -sum(tcrossprod(UtX_square,t(var_b))*D_vec_1)+sum(tcrossprod(X_square,t(var_b))))/2/s$sigma2
  
  return(loglik)
}

SER_posterior_e_loglik = function (X, Y, s2, s3,Eb, Eb2,D,U,UtX_square,X_square) {
  n = length(Y)
  rl_bar = Y - tcrossprod(X, t(Eb))
  rl_bar_U = crossprod(U,rl_bar)
  var_b = Eb2-(Eb)^2 
  D_vec_1 = s3*D/(1+s3*D)
  loglik = -(1/2) *n* (log(2*pi*s2))-1/2*sum(log(1+D*s3))-(sum(rl_bar^2)-sum(rl_bar_U^2*D_vec_1) +sum(tcrossprod(UtX_square,t(var_b))*D_vec_1)-sum(tcrossprod(X_square,t(var_b))))/2/s2
  return(loglik)
}



