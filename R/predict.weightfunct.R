#' Predicted Values for 'weightfunct' Objects
#'
#' This function calculates predicted conditional means and their corresponding standard errors for objects of class weightfunct.
#' @param object an object of class weightfunct
#' @param values a vector or matrix specifying the values of the moderator variables for which predicted values should be calculated; defaults to \code{NULL}
#' @param ... other arguments
#' @keywords weightr
#' @details \code{predict(object)} requires that the user specify a vector or matrix of predictor values. Without specifying values, the function will not work.
#'
#' For models including \code{y} number of moderator variables, users should set \code{values} equal to a \code{k} x \code{y} matrix, where \code{k} is the number of rows of data (i.e., "new" studies). In the example code, for example, there are 3 moderator variables and one row of data, so \code{values} is a 1 x 3 matrix. The intercept is incldued by default.
#'
#' Note that \code{weightfunct} handles categorical moderators automatically. To include them here, the appropriate contrast (dummy) variables must be explicitly specified. The \code{contrasts} function can help to understand the contrast matrix for a given factor.
#' @export
#' @return The function returns a list containing the following components: \code{unadjusted}, \code{adjusted}, and \code{values}. The \code{values} section simply prints the \code{values} matrix for verification. The \code{unadjusted} and \code{adjusted} sections print the conditional means for each row of new data, unadjusted and adjusted for publication bias (respectively), and their standard errors.
#' @examples
#' \dontrun{
#' test <- weightfunct(effect, v, mods=~x1 + x2 + x3, steps)
#' values <- matrix(c(0,1,0),ncol=3) # An arbitrary set of 3 dummy-coded moderators
#' predict(test, values)
#' }

predict.weightfunct <- function(object, values = NULL, ...){
  if (!inherits(object, "weightfunct"))
    stop("Argument 'object' must be an object of class \"weightfunct\".")
  if(is.null(object$mods)){
    stop("Error: You cannot specify values for moderators that do not exist.")
  }
  if(!(is.vector(values) || inherits(values, "matrix"))){
    stop(paste0("Argument 'values' should be a vector or matrix, but is of class '", class(values), "'."))
  }
  if(object$npred == 1L){
    k.new <- length(values)
    X.new <- cbind(c(values))
  }else{
    if(is.vector(values) || nrow(values) == 1L){ # user gives one vector or one row
      k.new <- 1
      X.new <- rbind(values)
    }else{ # user gives multiple rows
      k.new <- nrow(values)
      X.new <- cbind(values)
    }
  }
  if(inherits(X.new[1,1], "character")){
    stop("Argument 'values' should only contain numeric variables.")
  }
  if(ncol(X.new) != object$npred){
    stop("Number of values does not match number of moderator variables.")
  }
  if(object$fe){
    X.new <- cbind(int=rep(1,k.new), X.new)
  }else{
    X.new <- cbind(vc=rep(0,k.new), int=rep(1,k.new), X.new)
  }
  params <- matrix(object[1][[1]]$par, ncol=1, byrow=TRUE)
  pred <- apply(X.new, 1, function(x) {x %*% params}) # conditional means, unadjusted
  hess <- object[1][[1]]$hessian
  inv.hess <- solve(hess)
  se <- sqrt(apply(X.new, 1, function(x) {matrix(x,nrow=1) %*%
      tcrossprod(inv.hess, matrix(x,nrow=1))})) # standard errors for conditional means, unadj
  unadj <- data.frame(cbind(pred),cbind(se)) # make it look pretty
  params_adj <- matrix(object[2][[1]]$par, ncol=1, byrow=TRUE)
  weights.new <- matrix(rep(rep(0,(object$nsteps-1)),k.new),ncol=(object$nsteps-1),nrow=k.new) # adding placeholders for weights
  X.new_adj <- cbind(X.new,weights.new)
  pred_adj <- apply(X.new_adj, 1, function(x) {x %*% params_adj}) # conditional means, adjusted
  if(is.null(object$weights)){ # if this is V & H, NOT V & W, calculate SEs
    hess_adj <- object[2][[1]]$hessian
    inv.hess_adj <- solve(hess_adj)
    se_adj <- sqrt(apply(X.new_adj, 1, function(x) {matrix(x,nrow=1) %*%
        tcrossprod(inv.hess_adj, matrix(x,nrow=1))})) # standard errors for conditional means, adjusted
    adj <- data.frame(cbind(pred_adj),cbind(se_adj)) # make it look pretty
    colnames(adj) <- c("pred", "se")
  }else{ # otherwise, don't
    se_adj <- NULL
    adj <- data.frame(pred=cbind(pred_adj))
    colnames(adj) <- c("pred")
  }
  return(list(unadjusted = unadj, adjusted = adj,
              values = values)) # print this stuff out as a list
}
