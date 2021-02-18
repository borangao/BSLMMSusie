#' @title BSLMM extension of Sum of Single Effects (SuSiE) Regression
#'
#' @description Performs Bayesian multiple linear regression of Y on
#'   X; that is, this function fits the regression model \eqn{Y = \sum_l
#'   X b_{l=1}^L + z + e}, where elements of e are \emph{i.i.d.} normal with
#'   zero mean and variance \code{residual_variance}, and
#'   \eqn{\sum_{l=1}^L b_l} is a vector of length p representing the
#'   effects to be estimated. The \dQuote{susie assumption} is that each
#'   \eqn{b_l} has exactly one non-zero element. The prior on the
#'   non-zero element is normal with zero mean and variance \code{var(Y)
#'   * scaled_prior_variance}. The model is fitted using the
#'   \dQuote{Iterative Bayesian Stepwise Selection} (IBSS) algorithm.
#'   See also \code{\link{susie_trendfilter}} for applying susie to
#'   non-parametric regression, particularly changepoint problems.
#'
#' @details \code{susie_suff_stat} performs sum of single-effect
#' linear regression with summary statistics. The required summary
#' data are either: \code{bhat}, \code{shat}, the p by p symmetric,
#' positive semidefinite correlation (or covariance) matrix \code{R},
#' the sample size \code{n}, and the variance of y; or the p by p
#' matrix \eqn{X'X}, the p-vector \eqn{X'y}, the sum of squares
#' \eqn{y'y}, and the sample size \code{n}. The summary statistics
#' should come from the same individuals. Both the columns of X and
#' the vector y should be centered to have mean zero before computing
#' these summary statistics; you may also want to scale each column of
#' X and y to have variance 1 (see examples).
#'
#' \code{susie_rss} performs sum of single-effect linear regression
#' with z scores; all posterior calculations are for z-scores. This
#' function fits the regression model \eqn{z = \sum_l R*b_l + e},
#' where e is \eqn{N(0,residual_var*R)} and \eqn{\sum_l b_l} is a
#' p-vector of effects to be estimated. The required summary data are
#' the p by p correlation matrix, \code{R}, and the p-vector
#' \code{z}. The summary stats should come from the same individuals
#' (samples).
#'
#'
#'
#'
#' susie_auto is an attempt to automate reliable running of susie even
#' on hard problems. It implements a three-stage strategy for each L:
#' first, fit susie with very small residual error; next, estimate
#' residual error; finally, estimate the prior variance. If the last
#' step estimates some prior variances to be zero, stop. Otherwise,
#' double L, and repeat. Initial runs are performed with relaxed
#' tolerance; the final run is performed using the default susie
#' tolerance.
#'
#' @param X An n by p matrix of covariates.
#'
#' @param Y The observed responses, a vector of length n.
#'
#' @param L Number of components (nonzero coefficients) in the susie
#'   regression model. If L is larger than the number of covariates, p,
#'   L is set to p.
#'
#' @param scaled_prior_variance The scaled prior variance. This is
#'   either a scalar or a vector of length \code{L}. The prior variance
#'   of each non-zero element of b is set to \code{var(Y) *
#'   scaled_prior_variance}. If \code{estimate_prior_variance = TRUE},
#'   this provides initial estimates of the prior variances.
#'
#' @param residual_variance Variance of the residual. If
#'   \code{estimate_residual_variance = TRUE}, this value provides the
#'   initial estimate of the residual variance. By default, it is
#'   \code{var(Y)}.
#'
#' @param prior_weights A vector of length p, in which each entry
#'   gives the prior probability that corresponding column of X has a
#'   nonzero effect on the outcome, Y.
#'
#' @param null_weight Prior probability of no effect (a number between
#'   0 and 1, and cannot be exactly 1).
#'
#' @param standardize If \code{standardize = TRUE}, standardize the
#'   columns of X (or XtX and Xty) to unit variance prior to
#'   fitting. Note that \code{scaled_prior_variance} specifies the prior
#'   on the coefficients of X \emph{after} standardization (if it is
#'   performed). If you do not standardize, you may need to think more
#'   carefully about specifying \code{scaled_prior_variance}. Whatever
#'   your choice, the coefficients returned by \code{coef} are given for
#'   \code{X} on the original input scale. Any column of \code{X} that
#'   has zero variance is not standardized.
#'
#' @param intercept If \code{intercept = TRUE}, the intercept is
#'   fitted; it \code{intercept = FALSE}, the intercept is set to
#'   zero. Setting \code{intercept = FALSE} is generally not
#'   recommended.
#'
#' @param estimate_residual_variance If
#'   \code{estimate_residual_variance = TRUE}, the residual variance is
#'   estimated, using \code{residual_variance} as an initial value. If
#'   \code{estimate_residual_variance = FALSE}, the residual variance is
#'   fixed to the value supplied by \code{residual_variance}.
#'
#' @param estimate_prior_variance If \code{estimate_prior_variance =
#'   TRUE}, the prior variance is estimated (this is a separate
#'   parameter for each of the L effects). If provided,
#'   \code{scaled_prior_variance} is then used as an initial value for
#'   the optimization. When \code{estimate_prior_variance = FALSE}, the
#'   prior variance for each of the L effects is determined by the
#'   value supplied to \code{scaled_prior_variance}.
#'
#' @param estimate_prior_method The method used for estimating prior
#'   variance. When \code{estimate_prior_method = "simple"} is used, the
#'   likelihood at the specified prior variance is compared to the
#'   likelihood at a variance of zero, and the setting with the larger
#'   likelihood is retained.
#'
#' @param check_null_threshold When the prior variance is estimated,
#'   compare the estimate with the null, and set the prior variance to
#'   zero unless the log-likelihood using the estimate is larger by this
#'   threshold amount. For example, if you set
#'   \code{check_null_threshold = 0.1}, this will "nudge" the estimate
#'   towards zero when the difference in log-likelihoods is small. A
#'   note of caution that setting this to a value greater than zero may
#'   lead the IBSS fitting procedure to occasionally decrease the ELBO.
#'
#' @param prior_tol When the prior variance is estimated, compare the
#'   estimated value to \code{prior_tol} at the end of the computation,
#'   and exclude a single effect from PIP computation if the estimated
#'   prior variance is smaller than this tolerance value.
#'
#' @param residual_variance_upperbound Upper limit on the estimated
#'   residual variance. It is only relevant when
#'   \code{estimate_residual_variance = TRUE}.
#'
#' @param s_init A previous susie fit with which to initialize.
#'
#' @param coverage A number between 0 and 1 specifying the
#'   \dQuote{coverage} of the estimated confidence sets.
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param compute_univariate_zscore If \code{compute_univariate_zscore
#'   = TRUE}, the univariate regression z-scores are outputted for each
#'   variable.
#'
#' @param na.rm Drop any missing values in Y from both X and Y.
#'
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param tol A small, non-negative number specifying the convergence
#'   tolerance for the IBSS fitting procedure. The fitting procedure
#'   will halt when the difference in the variational lower bound, or
#'   \dQuote{ELBO} (the objective function to be maximized), is
#'   less than \code{tol}.
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#'   and a summary of the optimization settings, are printed to the
#'   console.
#'
#' @param track_fit If \code{track_fit = TRUE}, \code{trace}
#'   is also returned containing detailed information about the
#'   estimates at each iteration of the IBSS fitting procedure.
#'
#' @param residual_variance_lowerbound Lower limit on the estimated
#'   residual variance. It is only relevant when
#'   \code{estimate_residual_variance = TRUE}.
#'
#' @return A \code{"susie"} object with some or all of the following
#'   elements:
#'
#' \item{alpha}{An L by p matrix of posterior inclusion probabilites.}
#'
#' \item{mu}{An L by p matrix of posterior means, conditional on
#'   inclusion.}
#'
#' \item{mu2}{An L by p matrix of posterior second moments,
#'   conditional on inclusion.}
#'
#' \item{Xr}{A vector of length n, equal to \code{X \%*\% colSums(alpha
#'   * mu)}.}
#'
#' \item{lbf}{log-Bayes Factor for each single effect.}
#'
#' \item{lbf_variable}{log-Bayes Factor for each variable and single effect.}
#'
#' \item{intercept}{Intercept (fixed or estimated).}
#'
#' \item{sigma2}{Residual variance (fixed or estimated).}
#'
#' \item{V}{Prior variance of the non-zero elements of b, equal to
#'   \code{scaled_prior_variance * var(Y)}.}
#'
#' \item{elbo}{The value of the variational lower bound, or
#'   \dQuote{ELBO} (objective function to be maximized), achieved at
#'   each iteration of the IBSS fitting procedure.}
#'
#' \item{fitted}{Vector of length n containing the fitted values of
#'   the outcome.}
#'
#' \item{sets}{Credible sets estimated from model fit; see
#'   \code{\link{susie_get_cs}} for details.}
#'
#' \item{pip}{A vector of length p giving the (marginal) posterior
#'   inclusion probabilities for all p covariates.}
#'
#' \item{z}{A vector of univariate z-scores.}
#'
#' \item{niter}{Number of IBSS iterations that were performed.}
#'
#' \item{converged}{\code{TRUE} or \code{FALSE} indicating whether
#'   the IBSS converged to a solution within the chosen tolerance
#'   level.}
#'
#' \code{susie_suff_stat} returns also outputs:
#'
#' \item{XtXr}{A p-vector of \code{t(X)} times the fitted values,
#'   \code{X \%*\% colSums(alpha*mu)}.}
#'
#' \code{susie_rss} also outputs:
#'
#' \item{Rr}{An p-vector of \code{t(X)} times fitted values, \code{X
#'   \%*\% colSums(alpha*mu)}.}
#'
#' @export
susie_BSLMM <- function (X,Y,L = min(10,ncol(X)),
                   scaled_prior_variance = 0.2,
                   scaled_prior_kinship_variance = 0,
                   residual_variance = NULL,
                   prior_weights = NULL,
                   null_weight = NULL,
                   standardize = TRUE,
                   intercept = TRUE,
                   estimate_residual_variance = TRUE,
                   estimate_prior_variance = TRUE,
                   check_null_threshold = 0,
                   prior_tol = 1e-9,
                   residual_variance_upperbound = Inf,
                   s_init = NULL,
                   coverage = 0.95,
                   min_abs_corr = 0.5,
                   compute_univariate_zscore = FALSE,
                   na.rm = FALSE,
                   max_iter = 100,
                   tol = 1e-3,
                   verbose = FALSE,
                   track_fit = FALSE,
                   residual_variance_lowerbound = var(drop(Y))/1e4) {



  if (any(is.na(X)))
    stop("Input X must not contain missing values")
  if (any(is.na(Y))) {
    if (na.rm) {
      samples_kept = which(!is.na(Y))
      Y = Y[samples_kept]
      X = X[samples_kept,]
    } else
      stop("Input Y must not contain missing values")
  }

  # Check input Y.
  p = ncol(X)
  n = nrow(X)
  mean_y = mean(Y)



  eigen_X_list = eigen_X(X)
  D = eigen_X_list$D
  U = eigen_X_list$U
  UtX = eigen_X_list$UtX
  UtX_square = eigen_X_list$UtX_square
  X_square = X^2


  prior_weights = 1/p

  residual_variance =var(drop(Y))


  s = list(alpha  = matrix(1/p,nrow = L,ncol = p),
           mu     = matrix(0,nrow = L,ncol = p),
           mu2    = matrix(0,nrow = L,ncol = p),
           Xr     = rep(0,n),
           KL     = rep(as.numeric(NA),L),
           lbf    = rep(as.numeric(NA),L),
           lbf_variable = matrix(as.numeric(NA),L,p),
           sigma2 = residual_variance,
           sigma_b = scaled_prior_variance*residual_variance,
           sigma_k = scaled_prior_kinship_variance,
           pi     = prior_weights,
           null_index = 0)
  class(s) = "susie"


  elbo = rep(as.numeric(NA),max_iter + 1)
  elbo[1] = -Inf;
  tracking = list()
  n_iter = 0

  for (i in 1:max_iter) {
    D1 = s$sigma_k*D/(1+s$sigma_k*D)/s$sigma2
    UtXD1 = sweep(t(UtX),2,FUN="*",D1)
    UtXD1_square = apply(sweep(UtX_square,1,FUN="*",D1),2,sum)

    if (track_fit)
      tracking[[i]] = susie_slim(s)

    s = update_each_effect(X,Y,s,D,U,UtX_square,UtXD1,UtXD1_square,X_square,estimate_prior_variance = FALSE,check_null_threshold=0)
    elbo[i+1] = get_objective(X,Y,s,D, U,UtX_square,X_square)
    if (abs(elbo[i+1] - elbo[i]) < tol) {
      s$converged = TRUE
      break
    }


    variance_estimate =  estimate_residual_variance(X,Y,s,D,U,UtX,UtX_square,X_square)
    s$sigma2 = variance_estimate$sigma2
    #s$sigma_k = 0
    s$sigma_k = variance_estimate$sigma_k
    if (s$sigma_k < 0)
      s$sigma_k = 0
    if (s$sigma_k > 1)
      s$sigma_k = 1
    n_iter = n_iter + 1

  }
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_cs(s,coverage = coverage,X = X,
                          min_abs_corr = min_abs_corr)
    s$pip = susie_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }
  if (track_fit)
    s$trace = tracking
  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    s$converged = FALSE
  }

  s$beta =  colSums(s$alpha*s$mu)
  s$niter = n_iter
  s$ELBO = elbo[2:i+1]
  return(s)
}


