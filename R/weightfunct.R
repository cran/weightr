#' Estimate the Vevea and Hedges (1995) Weight-Function Model
#'
#' This function allows the user to estimate the Vevea and Hedges (1995) weight-function model for publication bias.
#' @param effect a vector of meta-analytic effect sizes.
#' @param v a vector of meta-analytic sampling variances; needs to match up with the vector of effects, such that the first element in the vector of effect sizes goes with the first element in the vector of sampling variances, and so on.
#' @param pval defaults to \code{NULL}. A vector containing observed p-values for the corresponding effect sizes. If not provided, p-values are calculated.
#' @param steps a vector of p-value cutpoints. The default only distinguishes between significant and non-significant effects (p < 0.05).
#' @param mods defaults to \code{NULL}. A formula specifying the linear model.
#' @param weights defaults to \code{FALSE}. A vector of prespecified weights for p-value cutpoints to estimate the Vevea and Woods (2005) model.
#' @param fe defaults to \code{FALSE}. Indicates whether to estimate a fixed-effect model.
#' @param table defaults to \code{FALSE}. Indicates whether to print a table of the p-value intervals specified and the number of effect sizes per interval.
#' @importFrom stats model.matrix optim pchisq pnorm qnorm model.frame na.action na.omit
#' @keywords weightr
#' @details This function allows meta-analysts to estimate both the
#' weight-function model for publication bias that was originally published in
#' Vevea and Hedges (1995) and the modified version presented in Vevea and Woods
#' (2005). Users can estimate both of these models with and without predictors and
#' in random-effects or fixed-effect situations. The function does not currently
#' accommodate models without an intercept.
#'
#' The Vevea and Hedges (1995) weight-function model is a tool for modeling publication
#' bias using weighted distribution theory. The model first estimates an unadjusted
#' fixed-, random-, or mixed-effects model, where the observed effect sizes are
#' assumed to be normally distributed as a function of predictors. This unadjusted
#' model is no different from the traditional meta-analytic model. Next, the Vevea
#' and Hedges (1995) weight-function model estimates an adjusted model that includes
#' not only the original mean model, fixed-, random-, or mixed-effects, but a series
#' of weights for any pre-specified p-value intervals of interest. This produces mean,
#' variance component, and covariate estimates adjusted for publication bias, as well
#' as weights that reflect the likelihood of observing effect sizes in each specified
#' interval.
#'
#' It is important to remember that the weight for each
#' estimated p-value interval must be interpreted relative to the first interval,
#' the weight for which is fixed to 1 so that the model is identified. In other
#' words, a weight of 2 for an interval indicates that effect sizes in that p-value
#' interval are about twice as likely to be observed as those in the first interval.
#' Finally, it is also important to remember that the model uses p-value cutpoints
#' corresponding to one-tailed p-values. This allows flexibility in the selection
#' function, which does not have to be symmetric for effects in the opposite direction;
#' a two-tailed p-value of 0.05 can therefore be represented as p < .025 or p > .975.
#'
#' After both the unadjusted and adjusted meta-analytic models are estimated, a
#' likelihood-ratio test compares the two. The degrees of freedom for this test are
#' equal to the number of weights being estimated. If the likelihood-ratio test is
#' significant, this indicates that the adjusted model is a better fit for the data,
#' and that publication bias may be a concern.
#'
#' To estimate a large number of weights for p-value intervals, the Vevea and Hedges
#' (1995) model works best with large meta-analytic datasets. It may have trouble
#' converging and yield unreliable parameter estimates if researchers, for instance,
#' specify a p-value interval that contains no observed effect sizes. However,
#' meta-analysts with small datasets are still likely to be interested in assessing
#' publication bias, and need tools for doing so. Vevea and Woods (2005)
#' attempted to solve this problem by adapting the Vevea and Hedges (1995) model to
#' estimate fewer parameters. The meta-analyst can specify p-value cutpoints,
#' as before, and specify corresponding fixed weights for those cutpoints. Then the
#' model is estimated. For the adjusted model, only the variance component and mean
#' model parameters are estimated, and they are adjusted relative to the fixed weights.
#' For example, weights of 1 for each p-value interval specified describes a situation
#' where there is absolutely no publication bias, in which the adjusted estimates are
#' identical to the unadjusted estimates. By specifying weights that depart from 1 over various p-value intervals, meta-analysts can
#' examine how various one-tailed or two-tailed selection patterns would alter their
#' effect size estimates. If changing the pattern of weights drastically changes
#' the estimated mean, this is evidence that the data may be vulnerable to
#' publication bias.
#'
#' For more information, consult the papers listed in the References section here.
#' Also, feel free to email the maintainer of \code{weightr} at kcoburn@ucmerced.edu.
#'
#' @export
#' @return The function returns a list of three lists. The first contains the following components: \code{output_unadj}, \code{output_adj}, \code{steps}, \code{mods}, \code{weights}, \code{fe}, \code{table}, \code{effect}, \code{v}, \code{npred}, \code{nsteps}, \code{k}, \code{p}, \code{removed}, and \code{XX}. \code{output_unadj} and \code{output_adj} return the results of the unadjusted and adjusted models, respectively, including Hessian matrices and model parameters. The other elements of this list are the arguments from \code{weightfunct}, as described above.
#'
#' The second list contains the following: \code{unadj_est}, \code{unadj_se}, \code{adj_est}, \code{adj_se}, \code{z_unadj}, \code{z_adj}, \code{p_unadj}, \code{p_adj}, \code{ci.lb_unadj}, \code{ci.ub_unadj}, \code{ci.lb_adj}, and \code{ci.ub_adj}. These are vectors of, respectively, the unadjusted and adjusted parameter estimates, standard errors, z-statistics, p-values, and 95% confidence interval boundaries.
#'
#' The third list contains information pertaining to heterogeneity tests: \code{QE}, \code{QEp}, \code{QM}, \code{QMp} (for the unadjusted model), and \code{QM2} and \code{QMp2} (for the adjusted model). These are the Q-values for tests of overall or excess heterogeneity and tests of moderators, along with their p-values. If no moderators are specified, the QM values will be \code{NA}.
#'
#' @references Coburn, K. M. & Vevea, J. L. (2015). Publication bias as a function
#' of study characteristics. Psychological Methods, 20(3), 310.
#'
#' Vevea, J. L. & Hedges, L. V. (1995). A general linear model for
#' estimating effect size in the presence of publication bias. Psychometrika, 60(3),
#' 419-435.
#'
#' Vevea, J. L. & Woods, C. M. (2005). Publication bias in research synthesis:
#' Sensitivity analysis using a priori weight functions. Psychological Methods, 10(4),
#' 428-443.
#' @examples
#' \dontrun{
#' # Uses the default p-value cutpoints of 0.05 and 1:
#'
#' weightfunct(effect, v)
#'
#' # Estimating a fixed-effect model, again with the default cutpoints:
#'
#' weightfunct(effect, v, fe=TRUE)
#'
#' # Specifying cutpoints:
#'
#' weightfunct(effect, v, steps=c(0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50, 1.00))
#'
#' # Including a linear model, where moderators are denoted as 'mod1' and mod2':
#'
#' weightfunct(effect, v, mods=~mod1+mod2)
#'
#' # Specifying cutpoints and weights to estimate Vevea and Woods (2005):
#'
#' weightfunct(effect, v, steps=c(0.01, 0.05, 0.50, 1.00), weights=c(1, .9, .7, .5))
#'
#' # Specifying cutpoints and weights while including a linear model:
#'
#' weightfunct(effect, v, mods=~mod1+mod2, steps=c(0.05, 0.10, 0.50, 1.00), weights=c(1, .9, .8, .5))
#' }

weightfunct <- function(effect, v, steps=c(0.025,1.00), mods=NULL,
                        weights=NULL, fe=FALSE, table=FALSE, pval=NULL){

  ## Function calculating the negative log-likelihood of the unadjusted
  ## meta-analytic model ##

  neglike_unadj <- function(pars) {
    if(fe == FALSE){
      vc = pars[1]
      beta = pars[2:(npred+2)]

      mn = XX%*%beta
      eta = sqrt(v + vc)
      b = 1/2 * sum(((effect-mn)/eta)^2)
      c = sum(log(eta))
    }
    else{
      beta = pars[1:(npred+1)]
      mn = XX%*%beta
      eta = sqrt(v)
      b = 1/2 * sum(((effect-mn)/eta)^2)
      c = sum(log(eta))
    }
    return(b+c)
  }

  ## Function calculating the negative log-likelihood of the adjusted
  ## meta-analytic model ##

  neglike_adj <- function(pars) {
    if(fe == FALSE){
      vc = pars[1]
      beta = pars[2:(npred+2)]
      if(is.null(weights) == FALSE){
        w = weights
      }
      else{
        w = c(1,pars[(npred+3):( (nsteps - 2) + (npred+3) )])
      }
      contrib = log(w[wt])

      mn = XX%*%beta
      a = sum(contrib)
      eta = sqrt(v + vc)
    }
    else{
      beta = pars[1:(npred+1)]
      if(is.null(weights) == FALSE){
        w = weights
      }
      else{
        w = c(1,pars[(npred+2):( (nsteps - 2) + (npred+2) )])
      }
      contrib = log(w[wt])

      mn = XX%*%beta
      a = sum(contrib)
      eta = sqrt(v)
    }
    b = 1/2 * sum(((effect-mn)/eta)^2)
    c = sum(log(eta))
    Bij <- matrix(rep(0,number*nsteps),nrow=number,ncol=nsteps)
    bi = -si * qnorm(steps[1])
    Bij[,1] = 1-pnorm((bi-mn)/eta)
    if(nsteps > 2){
      for(j in 2:(length(steps)-1)) {
        bi = -si * qnorm(steps[j])
        bilast = -si * qnorm(steps[j-1])
        Bij[,j] = pnorm((bilast-mn)/eta) - pnorm((bi-mn)/eta)
      }
    }
    bilast = -si * qnorm(steps[length(steps)-1])
    Bij[,length(steps)] = pnorm((bilast-mn)/eta)

    swbij = 0
    for(j in 1:length(steps)) swbij = swbij + w[j]*Bij[,j]
    d = sum(log(swbij))

    return(-a + b + c + d)
  }

  if(is.null(mods)){
    npred <- 0
    data <- data.frame(effect, v)
  }
  else{
    if(typeof(mods)=="language"){
        XX <- model.matrix(mods, model.frame(mods, na.action='na.pass'))
        npred <- dim(XX)[2]-1
        data <- data.frame(effect, v, XX)
    }else{
      stop('Moderators must be entered as "mods= ~ x1 + x2"')
    }
  }
  if(any(is.na(data))){
    data <- na.omit(data)
    removed <- as.numeric(na.action(data))
  }
  else{
    removed <- NULL
  }
  effect <- data[,1]
  v <- data[,2]
  if(npred == 0){
    XX <- cbind(rep(1,length(effect)))
  }
  else{
    XX <- as.matrix(data[,(3:(npred+3))])
  }

  if(length(effect)!=length(v)){
    stop('Your vector of effect sizes and your vector of sampling variances are not the same length. Please check your data.')
  }

  if(identical(effect,v)){
    stop('Your vector of effect sizes is exactly the same as your vector of sampling variances. Please check your data.')
  }

  if(min(v) < 0){
    stop('Sampling variances cannot be negative. Please check your data.')
  }

  si <- sqrt(v)
  ################# To add: If users specify p-values, trap any missing p-values and replace them
  ### with manually calculated ones.
  if(is.null(pval)){
    p <- 1-pnorm(effect/sqrt(v))
  }
  else{
    p <- pval
  }
  if(max(steps)!=1){
    steps <- c(steps,1)
  }
  if(max(steps) > 1 | min(steps) < 0){
    stop('p-value steps must be bounded by 0 and 1.')
  }
  if(length(unique(steps)) != length(steps)){
    stop('Two or more p-value cutpoints are identical.')
  }
  if(is.null(weights)){
    steps <- sort(steps)
  }
  if(is.null(weights) == FALSE){
    if(min(weights) < 0){
      stop('Weights for p-value intervals cannot be negative.')
    }
    if(length(weights)!=length(steps)){
      stop('The number of weights does not match the number of p-value intervals created.')
    }
    new <- cbind(steps, weights)
    steps <- new[order(steps),1]
    weights <- new[order(steps),2]
  }

  number <- length(effect)
  nsteps <- length(steps)

  wt <- rep(1,number)
  for(i in 1:number) {
    if(p[i] <= steps[1]) wt[i] = 1
    for(j in 2:nsteps) {
      if (steps[j-1] < p[i] && p[i] <= steps[j]) wt[i] = j
    }
    if(  p[i] > steps[nsteps-1] ) wt[i] = nsteps
  }

  pvalues <- as.numeric(table(intervaltally(p, steps)))

  if(sum(table(intervaltally(p,steps)) == 0) >= 1){
    warning('At least one of the p-value intervals contains no effect sizes, leading to estimation problems. Consider re-specifying the cutpoints.')
  }

  if(sum( table(intervaltally(p,steps)) > 0 & table(intervaltally(p, steps)) <= 3) >= 1){
    warning('At least one of the p-value intervals contains three or fewer effect sizes, which may lead to estimation problems. Consider re-specifying the cutpoints.')
  }

  if(is.null(mods)){

    pars <- c(mean(v)/4, mean(effect), rep(0,(nsteps-1)))

    if(fe == FALSE){

      output_unadj <- optim(par=pars[1:2],fn=neglike_unadj,lower=c(0,-Inf),method="L-BFGS-B",hessian=TRUE)

      output_adj <- optim(par=pars,fn=neglike_adj,lower=c(0,-Inf, rep(0.01,(nsteps-1))),method="L-BFGS-B",hessian=TRUE)

    }

    if(fe == TRUE){

      output_unadj <- optim(par=pars[2],fn=neglike_unadj,lower=c(-Inf),method="L-BFGS-B",hessian=TRUE)

      output_adj <- optim(par=pars[2:length(pars)],fn=neglike_adj,lower=c(-Inf, rep(0.01,(nsteps-1))),method="L-BFGS-B",hessian=TRUE)

    }

  }

  else{

    if(typeof(mods)=="language"){

      pars <- c(mean(v)/4, mean(effect), rep(0, npred), rep(1, (nsteps - 1)))

      if(fe == FALSE){

        output_unadj <- optim(par=pars[1:(npred+2)],fn=neglike_unadj,lower=c(0,rep(-Inf, (npred+1))),method="L-BFGS-B",hessian=TRUE)

        output_adj <- optim(par=pars,fn=neglike_adj,lower=c(0,rep(-Inf, (npred+1)),rep(0.01,(nsteps-1))),method="L-BFGS-B",hessian=TRUE)

      }

      if(fe == TRUE){

        output_unadj <- optim(par=pars[2:(npred+2)],fn=neglike_unadj,lower=c(rep(-Inf, (npred+1))),method="L-BFGS-B",hessian=TRUE)

        output_adj <- optim(par=pars[2:length(pars)],fn=neglike_adj,lower=c(rep(-Inf, (npred+1)),rep(0.01,(nsteps-1))),method="L-BFGS-B",hessian=TRUE)

      }
    }
  }

  unadj_est <- c(output_unadj$par)
  unadj_se <- c(sqrt(diag(solve(output_unadj$hessian))))
  adj_est <- c(output_adj$par)
  if(is.null(weights)){
    adj_se <- c(sqrt(diag(solve(output_adj$hessian))))
  }else{
    adj_se <- c(NULL)
  }

  ## This is a test for heterogeneity when there are no moderators specified
  ## and a test for residual heterogeneity when there ARE. The difference
  ## can be indicated conditionally in the print statement. It's the
  ## same for both the unadjusted and adjusted models as it doesn't involve
  ## any actual parameter estimates.

  w.FE <- diag((1/v), nrow=number, ncol=number)
  stXWX <- tcrossprod(qr.solve( (sqrt(w.FE) %*% XX), diag(number) ))
  P <- w.FE - w.FE %*% XX %*% stXWX %*% crossprod(XX, w.FE)
  QE <- max(0, c(crossprod(as.matrix(effect), P) %*% as.matrix(effect)))
  QEp <- pchisq(QE, df=(number), lower.tail=FALSE)

  ## This is the test for heterogeneity in the unadjusted model when
  ## there are moderators present.
  ## There's one of these tests for each of the adjusted and unadjusted models
  ## because the conditional mean estimates are adjusted.
  if(npred > 0){
    if(fe==FALSE){
      beta_vect <- c(output_unadj$par[2:(npred+2)])
      w.M <- 1/(v+output_unadj$par[1])
      beta_vect2 <- c(output_adj$par[2:(npred+2)])
      w.M2 <- 1/(v+output_adj$par[1])
    }else{
      beta_vect <- c(output_unadj$par[1:(npred+1)])
      w.M <- 1/v
      beta_vect2 <- c(output_adj$par[1:(npred+1)])
      w.M2 <- 1/v
    }
    w.M <- diag(w.M, nrow=number, ncol=number)
    w.M2 <- diag(w.M2, nrow=number, ncol=number)
    swx <- sqrt(w.M) %*% XX
    swx2 <- sqrt(w.M2) %*% XX
    res.qrs <- qr.solve(swx, diag(number))
    res.qrs2 <- qr.solve(swx2, diag(number))
    vb <- tcrossprod(res.qrs)
    vb2 <- tcrossprod(res.qrs2)
    QM <- as.vector(t(beta_vect[2:(npred+1)]) %*% chol2inv(chol(vb[(2:(npred+1)),(2:(npred+1))])) %*% beta_vect[2:(npred+1)])
    QMp <- pchisq(QM, df=(length(beta_vect)-1), lower.tail=FALSE)
    QM2 <- as.vector(t(beta_vect2[2:(npred+1)]) %*% chol2inv(chol(vb2[(2:(npred+1)),(2:(npred+1))])) %*% beta_vect2[2:(npred+1)])
    QMp2 <- pchisq(QM2, df=(length(beta_vect2)-1), lower.tail=FALSE)
  }else{
    ## If there's no moderators, QM doesn't even get printed, so
    ## it's replaced with NA.
    QM <- NA
    QMp <- NA
    QM2 <- NA
    QMp2 <- NA
  }

    results <- c(list(output_unadj=output_unadj,
                      output_adj=output_adj,
                      steps=steps,
                      mods=mods,
                      weights=weights,
                      fe=fe,
                      table=table,
                      effect=effect,
                      v=v,
                      npred=npred,
                      nsteps=nsteps,
                      k=number,
                      p=p,
                      removed=removed,
                      XX=XX),
                 list(unadj_est=cbind(unadj_est),
                      unadj_se=cbind(unadj_se),
                      adj_est=cbind(adj_est),
                      adj_se=cbind(adj_se),
                      z_unadj=cbind(unadj_est/unadj_se),
                      z_adj=cbind(adj_est/adj_se),
                      p_unadj=cbind(2*pnorm(-abs(unadj_est/unadj_se))),
                      p_adj=cbind(2*pnorm(-abs(adj_est/adj_se))),
                      ci.lb_unadj=cbind(unadj_est - qnorm(0.975) * unadj_se),
                      ci.ub_unadj=cbind(unadj_est + qnorm(0.975) * unadj_se),
                      ci.lb_adj=cbind(adj_est - qnorm(0.975) * adj_se),
                      ci.ub_adj=cbind(output_adj$par + qnorm(0.975) * adj_se)
                      ),
                 list(
                   QE = QE,
                   QEp = QEp,
                   QM = QM,
                   QMp = QMp,
                   QM2 = QM2,
                   QMp2 = QMp2
                 )
                 )

  class(results) <- c("weightfunct")
  return(results)

}
