#' Variable selection using group LASSO (gglasso, grpreg that use distinct mathematical approaches), sparse group LASSO, glmnet, stability selection with error control
#'
#' @param X input a matrix of exposure variables for selection
#' @param Y input a matrix of outcome variables
#' @param group input the grouping structure. Group inputs should be integers
#'
#' @return Returns a list of lists that contain beta coefficient for each exposure. Each list correspond to a outcome variable provided in the input
#' @export var_selec
#'
#' @import gglasso
#' @import grpreg
#' @import sparsegl
#' @import glmnet
#' @import stabs
#' @import pacman
#'
#' @example
#' \donttest{
#' Define the predictor variables
#' X <- as.matrix(dat[c(1:32)])
#' Y <- as.matrix(dat[c(33:35)])
#'Group index for X variables
#' group<- as.integer(c(rep(1,times=3), rep(2,times=6), rep(3, times=3),
#'                      rep(4,times=4), rep(5, times=3), rep(6,times=2),
#'                      rep(7,times=6), rep(8, times=2), rep(9,times=3)))
#'res<- var_selec(X, Y, group)
#' }
var_selec <- function(X, Y, group) {
  results <- list()
  for (i in 1:3) {
    y <- Y[, i]
    # Group graphical Lasso
    gr_cv <- cv.gglasso(X, y, group=group, loss="ls", pred.loss="L2",  nfolds=10)
    gr_min_beta <- coef(gr_cv, s = gr_cv$lambda.min)[-1]

    # Group Lasso
    grpp_cv <- cv.grpreg(X, y, group = group, penalty="grLasso",seed=5678,nfolds = 10)
    grpp_min_beta <- coef(grpp_cv, s = grpp_cv$lambda.min)[-1]

    #Sparse lasso
    sparse_cv<- cv.sparsegl(X, y, group = group, family = "gaussian", nfolds = 10)
    sparse_min_beta<- coef(sparse_cv, s= sparse_cv$lambda.min)[-1]


    #Stability selection with error control - input cross-validated lambda.min from cv-glmnet
    stab_lambda_min <- cv.glmnet(X, y, nfolds=10)$lambda.min
    stab_maxCoef <- stabsel(X, y, fitfun = glmnet.lasso_maxCoef, args.fitfun = list(lambda = stab_lambda_min), cutoff = 0.75, PFER = 1)
    stab_maxCoef_selec<- stab_maxCoef$max

    # Store results in list
    results[[paste0("outcome", i)]] <- as.data.frame(list(gr_lasso = gr_min_beta,
                                                          grpp_lasso = grpp_min_beta,
                                                          sparse_lasso = sparse_min_beta,
                                                          stab_cv_glmnet = stab_maxCoef_selec))
  }
  # Return list of results
  return(results)
}
