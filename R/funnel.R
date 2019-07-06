#' Create a Funnel Plot
#'
#' This function allows you to create a funnel plot using a vector of effect sizes and a vector of their corresponding sampling variances.
#' @param effect a vector of meta-analytic effect sizes
#' @param v a vector of sampling variances
#' @param type \code{v} for variance or \code{se} for standard error; defaults to standard error
#' @param flip \code{FALSE} (default) to plot effect sizes on the vertical axis (for a horizontal funnel); \code{TRUE} to plot them on the horizontal axis (for a vertical funnel)
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
#' # Vertical funnel plot using standard error:
#' funnel(effect, v, flip=TRUE)
#' }

funnel <- function(effect, v, type='se', flip=FALSE){
  x <- sqrt(v)
  if(type=='v'){
    x <- v
  }
  xinc <- 0.05*(max(x) - min(x))
  yinc <- 0.05*(max(effect) - min(effect))
  xmin <- min(x) - xinc
  xmax <- max(x) + xinc
  ymin <- min(effect) - yinc
  ymax <- max(effect) + yinc
  if(type=='se'){
    if(flip == FALSE){
      plot(x, effect,
           xlab="Standard Error", ylab="Effect Size",
           xlim=c(xmin,xmax),
           ylim=c(ymin,ymax))
    }
    if(flip == TRUE){
      plot(effect, x,
           xlab="Effect Size", ylab="Standard Error",
           xlim=c(ymin,ymax),
           ylim=c(xmin,xmax))
    }
  }
  if(type=='v'){
    if(flip == FALSE){
      plot(x, effect,
           xlab="Sampling Variance", ylab="Effect Size",
           xlim=c(xmin,xmax),
           ylim=c(ymin,ymax))
    }
    if(flip == TRUE){
      plot(effect, x,
           xlab="Effect Size", ylab="Sampling Variance",
           xlim=c(ymin,ymax),
           ylim=c(xmin,xmax))
    }
  }
  box(lwd=1.2)
}
