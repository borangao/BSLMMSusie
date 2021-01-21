set_X_attributes = function(X, center = TRUE, scale = TRUE) {
  
  # if X is a trend filtering matrix
  if (!is.null(attr(X,"matrix.type"))) {
    order = attr(X,"order")
    n = ncol(X)
    
    # Set three attributes for X.
    attr(X,"scaled:center") = compute_tf_cm(order,n)
    attr(X,"scaled:scale") = compute_tf_csd(order,n)
    attr(X,"d") = compute_tf_d(order,n,attr(X,"scaled:center"),
                               attr(X,"scaled:scale"),scale,center)
    if (!center)
      attr(X,"scaled:center") = rep(0,n)
    if (!scale)
      attr(X,"scaled:scale") = rep(1,n)
  } else {
    
    # If X is either a dense or sparse ordinary matrix.
    # Get column means.
    cm = colMeans(X,na.rm = TRUE)
    
    # Get column standard deviations.
    csd = compute_colSds(X)
    
    # Set sd = 1 when the column has variance 0.
    csd[csd == 0] = 1
    if (!center)
      cm = rep(0,length = length(cm))
    if (!scale) 
      csd = rep(1,length = length(cm))
    X.std = (t(X) - cm)/csd
    
    # Set three attributes for X.
    attr(X,"d") = rowSums(X.std * X.std)
    attr(X,"scaled:center") = cm
    attr(X,"scaled:scale") = csd
  }
  return(X)
}

eigen_X = function(X) {
  p = dim(X)[2]
 # X = scale(X)
  XtX = t(X)%*%X/p
  eigen_decomposition = eigen(XtX,TRUE)
  eigen_vec = eigen_decomposition$vectors
  U = X%*%(sweep(eigen_vec,2,FUN="/",sqrt(eigen_decomposition$values)*sqrt(p)))
  UtX = t(U)%*%X
  UtX_square = apply(UtX,2,function(x)x^2)
  D = eigen_decomposition$values
  return(eigen_X_list = list(D = D,U = U, UtX = UtX, UtX_square = UtX_square))
}




  
