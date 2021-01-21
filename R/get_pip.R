susie_get_pip = function (res, prune_by_cs = FALSE, prior_tol = 1e-9) {
  
  if (inherits(res,"susie")) {
    
    # Drop null weight columns.
    if (res$null_index > 0)
      res$alpha = res$alpha[,-res$null_index,drop=FALSE]
    
    # Drop the single-effects with estimated prior of zero.
    if (is.numeric(res$V))
      include_idx = which(res$V > prior_tol)
    else
      include_idx = 1:nrow(res$alpha)
    
    # Only consider variables in reported CS.
    # This is not what we do in the SuSiE paper.
    # So by default prune_by_cs = FALSE means we do not run the
    # following code.
    if (!is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = intersect(include_idx,res$sets$cs_index)
    if (is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = numeric(0)
    
    # now extract relevant rows from alpha matrix
    if (length(include_idx) > 0)
      res = res$alpha[include_idx,,drop = FALSE]
    else
      res = matrix(0,1,ncol(res$alpha))
  }
  
  return(as.vector(1 - apply(1 - res,2,prod)))
}


