#' Create a Density Plot
#'
#' This function allows you to create a plot displaying the unadjusted and adjusted densities of the specified model. Note that you must first specify a model using \code{weightfunct}.
#' @param x an object of class weightfunct
#' @param ... other arguments
#' @importFrom ggplot2 ggplot aes geom_line theme_bw geom_hline labs xlab ylab theme element_text element_blank scale_x_continuous scale_y_continuous
#' @importFrom stats median
#' @importFrom graphics box plot
#' @importFrom scales pretty_breaks
#' @keywords weightr
#' @details This function produces an approximate graphical illustration of the estimated unweighted and weighted densities. The unweighted density is represented by a dashed line and the weighted density by a solid line. For the unweighted density, the effect sizes are assumed to be normally distributed, with a mean equal to their unadjusted mean and a variance equal to their unadjusted variance component plus their individual sampling variances. This plot is an approximation because it is necessary to use a fixed sampling variance; here, we fix the sampling variance to the median of the distribution of sampling variances.
#'
#' For the adjusted density, the expected density for effect sizes within each specified p-value interval is multiplied by the estimated weight for the corresponding interval. Greater density in an interval then represents a greater likelihood of effect-size survival. (Remember, of course, that the weight for the first interval is fixed to one, and other intervals should be interpreted relative to it.) Each discontinuity in the solid line, therefore, represents a p-value cutpoint.
#'
#' Users may wonder why the adjusted density, or the solid line, sometimes falls outside of the unadjusted density, or the dashed line. In answer, recall that the mean and variance of the adjusted density also differ. Based on the severity of this difference, the adjusted density may fall outside of its unadjusted counterpart.
#' @export
#' @examples
#' \dontrun{
#' test <- weightfunct(effect, v, steps)
#' density(test)
#' }

density <- function(x, ...){
  if (!inherits(x, "weightfunct")){
    stop("Argument 'x' must be an object of class \"weightfunct\".")}

  unadj_est <- cbind(c(x[[1]]$par[1:2]))
  adj_est <- cbind(c(x[[2]]$par[1:2]))
  weights <- cbind(c(x[[2]]$par[(x$npred+3):(length(x$steps)+(x$npred+1))]))
  vc1 <- unadj_est[1]
  mu1 <- unadj_est[2]
  vc2 <- adj_est[1]
  mu2 <- adj_est[2]
  cuts <- x$steps
  x_low_lim <- min(x$effect) - 2
  x_up_lim <- max(x$effect) + 2

  xfull <- seq(x_low_lim,x_up_lim,.01)

  vi <- median(x$v)

  fx <- ( 1/(sqrt(2*pi*(vi + vc1))) ) * exp( -1/2*( (xfull - mu1)^2 / (vi + vc1) ) )
  yfull <- fx
  A0 <- sum(rep(.01,length(xfull))*yfull)

  fx2 <- ( 1/(sqrt(2*pi*(vi + vc1))) ) * exp( -1/2*( (xfull - mu1)^2 / (vi + vc1) ) )

  testlist <- -1 * qnorm(x$steps, 0, sqrt(vi + vc2))

  testxfull <- findInterval(xfull,sort(testlist))
  xlist <- split(xfull, testxfull)

  ylist <- split(fx2, testxfull)
  weights2 <- rev(c(1, weights))
  testyfull <- mapply("*", ylist, weights2)

  A1 <- sum(rep(.01,length(unlist(xlist)))*unlist(testyfull))

  test <- data.frame(rep(xfull, 2), c((yfull/A0), (as.numeric(unlist(testyfull))/A1)),
                     c(rep("Unadjusted", length((yfull/A0))), rep("Adjusted",
                                                             length(as.numeric(unlist(testyfull))))))
  test2 <- subset(test, (test[,2] > 0.001))
  colnames(test2) <- c("xval", "density", "model")
  ggplot(data = test2, aes(x = test2$xval, y = test2$density, group = test2$model)) +
    geom_line(aes(linetype = test2$model)) +
    theme_bw() +
    geom_hline(yintercept = 0) +
    labs(title = "Adjusted and Unadjusted Densities", linetype = "Model") +
    xlab("Sample Effect Size") + ylab("Density") +
    theme(text = element_text(size=17), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

}
