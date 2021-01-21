estimate_residual_variance<-function(X,Y,s,D,U,UtX,UtX_square,X_square){
  n = length(Y)
  rl = Y - s$Xr
  rl_U = crossprod(U,rl)
  rl_U_2 = rl_U^2
  D_sigma = 1/(1+s$sigma_k*D)
  D_sigma_2 = 1-D_sigma
  var_b = colSums((s$alpha * s$mu2)-(s$alpha * s$mu)^2) 
  sigma2 = (sum(tcrossprod(X_square,t(var_b)))-sum(tcrossprod(UtX_square,t(var_b))*D_sigma)+sum(rl^2)-sum(rl_U_2*D_sigma_2))/n

 # sum(diag(t(rl_U)%*%diag(D_sigma)%*%(rl_U)))
 # D_sigma_2 = s$sigma_k*D/(1+s$sigma_k*D)
  
  D_sigma_k_first_order_1 = D_sigma*D
  D_sigma_k_first_order_2 = D_sigma_k_first_order_1*D_sigma
  sigma_k_first_order = (sum(rl_U_2*D_sigma_k_first_order_2)+sum(tcrossprod(UtX_square,t(var_b))*D_sigma_k_first_order_2))/2/s$sigma2-sum(D_sigma_k_first_order_1)/2

  D_sigma_k_second_order_1 = D_sigma_k_first_order_2*D
  D_sigma_k_second_order_2 = D_sigma_k_second_order_1*D_sigma
  
  sigma_k_second_order= sum(D_sigma_k_second_order_1)/2-(sum(rl_U_2* D_sigma_k_second_order_2)+sum(tcrossprod(UtX_square,t(var_b))*D_sigma_k_second_order_2))/s$sigma2
  
  sigma_k = s$sigma_k-sigma_k_first_order/sigma_k_second_order
   # D_sigma_k_second_order = 
   # sigma2 = (sum(diag(X%*%diag(var_b)%*%t(X))*D_sigma)+sum(diag(t(rl_U)%*%diag(D_sigma)%*%(rl_U))))/dim(X)[2]
 # Xr_L = tcrossprod(X,s$alpha * s$mu)
 # Xr_L_square = Xr_L^2 
 # postb2 = s$alpha * s$mu2
 # Xr_L_second_all<-c()
 # for(L in 1:10){
 #   Xr_L_second = diag(X%*%diag(postb2[L,])%*%t(X))
 #   Xr_L_second_all<-rbind(Xr_L_second_all,Xr_L_second )
 # }
 #sum(t(Xr_L_second_all)-Xr_L_square)
 # Xr_L_square = Xr_L^2 
 # sum(rl^2)
  
 # D_sigma_k =(D-2*s$sigma_k*D*D/(1+s$sigma_k*D)+(s$sigma_k)^2*D*D*D/(1+s$sigma_k*D)/(1+s$sigma_k*D) )/s$sigma2
  
 # gradient = 0.5*sum(rl_U^2*D_sigma_k)+0.5*sum(diag(UtX%*%diag(var_b)%*%t(UtX))*D_sigma_k)-0.5*sum(diag(U%*%diag(D-s$sigma_k*D/(1+s$sigma_k*D))%*%t(U)))
 # s$sigma_k = s$sigma_k
  return(list(sigma2 = sigma2,sigma_k =sigma_k))
  
}
