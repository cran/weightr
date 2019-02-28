#' Create a Funnel Plot
#'
#' This function allows you to create a funnel plot using a vector of effect sizes and a vector of their corresponding sampling variances.
#' @param effect a vector of meta-analytic effect sizes
#' @param v a vector of sampling variances
#' @param type \code{v} for variance or \code{se} for standard error; defaults to standard error
#' @keywords weightr
#' @details This funnel plot, by default, plots the effect sizes on the y-axis and the measure of study size (either variance or standard error) on the x-axis. If no asymmetry is present, the plot should resemble a horizontal funnel.
#'
#' Users can choose either standard error (default) or sampling variance as a measure of study size. The choice is mostly arbitrary. In both cases, however, \code{v} must be a vector of variances, the same as that required by \code{weightfunct}. The conversion to standard error is automatic.
#' @export
#' @examples
#' \dontrun{
#' # Funnel plot using standard error (default):
#' funnel(effect, v)
#' # Funnel plot using sampling variance:
#' funnel(effect, v, type='v')
#' }

funnel <- function(effect, v, type='se'){
  se <- sqrt(v)
  xinc <- 0.05*(max(se) - min(se))
  yinc <- 0.05*(max(effect) - min(effect))
  xmin <- min(se) - xinc
  xmax <- max(se) + xinc
  ymin <- min(effect) - yinc
  ymax <- max(effect) + yinc
  plot(se, effect,
       xlab="Standard Error", ylab="Effect Size",
       xlim=c(xmin,xmax),
       ylim=c(ymin,ymax))
  box(lwd=1.2)
  if(type=='v'){
    xinc <- 0.05*(max(v) - min(v))
    yinc <- 0.05*(max(effect) - min(effect))
    xmin <- min(v) - xinc
    xmax <- max(v) + xinc
    ymin <- min(effect) - yinc
    ymax <- max(effect) + yinc
    plot(v, effect,
         xlab="Sampling Variance", ylab="Effect Size",
         xlim=c(xmin,xmax),
         ylim=c(ymin,ymax))
    box(lwd=1.2)
  }
}
