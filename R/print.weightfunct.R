#' Print Model Results
#'
#' This function allows you to print the model results.
#' @param x an object of class weightfunct
#' @param ... other arguments
#' @keywords weightr
#' @export
#' @importFrom stats model.matrix optim pchisq pnorm qnorm
#' @examples
#' \dontrun{
#' print(weightfunct(d,v))
#' }
print.weightfunct <- function(x, ...){
  if (!inherits(x, "weightfunct"))
    stop("Argument 'x' must be an object of class \"weightfunct\".")

  ####### Unadjusted model ########

    cat("\n")
    cat("Unadjusted Model (k = ", x$k, "):", sep="")
    cat("\n\n")
    # Heterogeneity estimates
    if(x$fe == FALSE){
      cat("tau^2 (estimated amount of total heterogeneity): ", formatC(round(x$unadj_est[1], digits = 4), digits = 4, format = "f"), " (SE = ", formatC(round(x$unadj_se[1], digits = 4),digits = 4, format = "f"), ")", sep="")
      cat("\n")
      cat("tau (square root of estimated tau^2 value): ", formatC(round(sqrt(x$unadj_est[1]), digits = 4),digits = 4, format = "f"))
      cat("\n\n")
      if(x$npred==0){
        cat("Test for Heterogeneity:")
      }
      if(x$npred > 0){
        cat("Test for Residual Heterogeneity:")
      }
      cat("\n")
      cat("Q(df = ", (x$k-x$npred-1), ") = ", formatC(round(x$QE, digits=4), digits=4, format = "f"), ", p-val = ", x$QEp, sep="")
      cat("\n\n")
    }else{
      cat("Test for Residual Heterogeneity:")
      cat("\n")
      cat("Q(df = ", (x$k-x$npred-1), ") = ", formatC(round(x$QE, digits=4), digits=4, format = "f"), ", p-val = ", x$QEp, sep="")
      cat("\n\n")
    }
    if(x$npred > 0){
        cat("Test of Moderators (coefficients 2:", (x$npred + 1), "):", sep="")
        cat("\n")
        cat("QM(df = ", x$npred, ") = ", formatC(round(x$QM, digits=4), digits=4, format = "f"), ", p-val = ", x$QMp, sep="")
        cat("\n\n")
    }
    cat("Model Results:")
    cat("\n\n")
    if(x$fe == FALSE){
      unadj_est <- cbind(x[[1]]$par[2:(x$npred+2)])
      unadj_se <- cbind(sqrt(diag(solve(x[[1]]$hessian)))[2:(x$npred+2)])
    }
    if(x$fe == TRUE){
      unadj_est <- cbind(x[[1]]$par[1:(x$npred+1)])
      unadj_se <- cbind(sqrt(diag(solve(x[[1]]$hessian)))[1:(x$npred+1)])
    }
    z_stat <- unadj_est/unadj_se
    p_val <- (2*pnorm(-abs(z_stat)))
    ci.lb <- unadj_est - qnorm(0.975) * unadj_se
    ci.ub <- unadj_est + qnorm(0.975) * unadj_se
    res.table <- data.frame(matrix(c(unadj_est, unadj_se, z_stat, p_val, ci.lb, ci.ub), nrow=(x$npred+1), byrow=F),stringsAsFactors=FALSE)
    rowlabels <- rep(0, (x$npred+1))
    rowlabels[1] <- "Intercept"
    if(x$npred > 0){
      for(i in 2:(x$npred+1)){
        rowlabels[i] <- paste(c(colnames(x$XX)[i]))
      }
    }
    row.names(res.table) <- c(rowlabels)
    colnames(res.table) <- c("estimate","std.error","z-stat","p-val","ci.lb","ci.ub")
    res.table[,4] <- format.pval(res.table[,4])
    res.table[,c(1,2,3,5,6)] <- format(res.table[,c(1,2,3,5,6)], digits=4)
    print.data.frame(res.table)

    ####### Adjusted model ########

      cat("\n")
      cat("Adjusted Model (k = ", x$k, "):", sep="")
      cat("\n\n")
      if(x$fe == FALSE){
        if(is.null(x$weights)){
          if(is.nan(suppressWarnings(sqrt(diag(solve(x[[2]]$hessian)))[1]))){
            warning('The adjusted variance component is so close to zero that a border condition prevents a meaningful iterative solution. As long as other model estimates are still \nreasonable, the results are identical to those from a fixed-effect analysis.')
          }
          suppressWarnings(
          cat("tau^2 (estimated amount of total heterogeneity): ", formatC(round(x$adj_est[1], digits = 4),digits = 4, format = "f"), " (SE = ", formatC(round(x$adj_se[1], digits = 4),digits = 4, format = "f"), ")", sep="")
          )
        }
        if(is.null(x$weights) == FALSE){
          cat("tau^2 (estimated amount of total heterogeneity): ", formatC(round(x$adj_est[1], digits = 4),digits = 4, format = "f"), " (SE = ", "---", ")", sep="")
        }
        cat("\n")
        cat("tau (square root of estimated tau^2 value): ", formatC(round(sqrt(x$adj_est[1]), digits = 4),digits = 4, format = "f"))
        cat("\n\n")
        if(x$npred==0){
          cat("Test for Heterogeneity:")
        }
        if(x$npred > 0){
          cat("Test for Residual Heterogeneity:")
        }
        cat("\n")
        cat("Q(df = ", (x$k-x$npred-1), ") = ", formatC(round(x$QE, digits=4), digits=4, format = "f"), ", p-val = ", x$QEp, sep="")
        cat("\n\n")
      }else{
        cat("Test for Residual Heterogeneity:")
        cat("\n")
        cat("Q(df = ", (x$k-x$npred-1), ") = ", formatC(round(x$QE, digits=4), digits=4, format = "f"), ", p-val = ", x$QEp, sep="")
        cat("\n\n")
      }

      if(x$npred > 0){
        if(is.null(x$weights)){
          cat("Test of Moderators (coefficients 2:", (x$npred + 1), "):", sep="")
          cat("\n")
          cat("QM(df = ", x$npred, ") = ", formatC(round(x$QM2, digits=4), digits=4, format = "f"), ", p-val = ", x$QMp2, sep="")
          cat("\n\n")
        }else{

        }
      }

      cat("Model Results:")
      cat("\n\n")

      if(x$fe == FALSE){
        if(is.null(x$weights)){
          adj_int_est <- cbind(x$adj_est[2:( (x$nsteps - 1) + (x$npred+2) )])
          adj_int_se <- cbind(x$adj_se[2:( (x$nsteps - 1) + (x$npred+2) )])
        }
        else{
          adj_int_est <- cbind(c(
            round(x$adj_est[2:( (x$npred+2) )], digits=4),
            sprintf('%.4f', x$weights[2:length(x$weights)])
            ))
          adj_int_se <- cbind(rep("---", length(x[[2]]$par[2:length(x[[2]]$par)])))
        }
      }
      if(x$fe == TRUE){
        if(is.null(x$weights)){
          adj_int_est <- cbind(x$adj_est[1:( (x$nsteps - 1) + (x$npred+1) )])
          adj_int_se <- cbind(x$adj_se[1:( (x$nsteps - 1) + (x$npred+1) )])
        }
        else{
          adj_int_est <- cbind(c(
            round(x$adj_est[1:( (x$npred+1) )], digits=4),
            sprintf('%.4f', x$weights[2:length(x$weights)])
            ))
          adj_int_se <- cbind(rep("---", length(x[[2]]$par[1:length(x[[2]]$par)])))
        }
      }

      if(is.null(x$weights)){
        z_stat_int <- adj_int_est/adj_int_se
        p_val_int <- (2*pnorm(-abs(z_stat_int)))
        ci.lb_int <- adj_int_est - qnorm(0.975) * adj_int_se
        ci.ub_int <- adj_int_est + qnorm(0.975) * adj_int_se
      }else{
        if(x$fe == FALSE){
          length_a <- length(x[[2]]$par[2:length(x[[2]]$par)])
          z_stat_int <- rep("---", length_a)
          p_val_int <- rep("---", length_a)
          ci.lb_int <- rep("---", length_a)
          ci.ub_int <- rep("---", length_a)
        }else{
          length_aF <- length(x[[2]]$par[1:length(x[[2]]$par)])
          z_stat_int <- rep("---", length_aF)
          p_val_int <- rep("---", length_aF)
          ci.lb_int <- rep("---", length_aF)
          ci.ub_int <- rep("---", length_aF)
        }
      }
      res.table <- data.frame(matrix(c(adj_int_est, adj_int_se, z_stat_int, p_val_int, ci.lb_int, ci.ub_int), nrow=(x$npred+1+(x$nsteps-1)), byrow=F),stringsAsFactors=FALSE)

      rowlabels1 <- rep(0, (x$npred+1))
      rowlabels1[1] <- "Intercept"
      if(x$npred > 0){
        for(i in 2:length(rowlabels1)){
          rowlabels1[i] <- paste(c(colnames(x$XX)[i]))
        }
      }
      rowlabels2 <- rep(0, (x$nsteps-1))
      for(i in 1:(length(rowlabels2))){
        rowlabels2[i] <- paste(c(x$steps[i], "< p <", x$steps[i + 1]), collapse=" ")
      }
      row.names(res.table) <- c(rowlabels1,rowlabels2)
      colnames(res.table) <- c("estimate","std.error","z-stat","p-val","ci.lb","ci.ub")
      if(is.null(x$weights)){
        res.table[,"p-val"] <- format.pval(res.table[,"p-val"])
      }
      res.table[,c(1,2,3,5,6)] <- format(res.table[,c(1,2,3,5,6)], digits=4)
      print.data.frame(res.table)

      ####### LRT ########

        if(is.null(x$weights)){
          cat("\n")
          cat("Likelihood Ratio Test:")
          cat("\n")
          df <- length(x[[2]]$par) - length(x[[1]]$par)
          lrchisq <- 2*(abs(x[[1]]$value - x[[2]]$value))
          pvalue <- 1-pchisq(lrchisq,df)
          cat("X^2(df = ", df, ") = ", lrchisq, ", p-val = ", format.pval(pvalue), sep="")
        }else{
          cat("\n")
          cat("Note: The symbol --- appears because the user has specified weights,\nchoosing to use the Vevea and Woods model, which does not estimate \nweights for p-value intervals and therefore cannot produce meaningful \nstandard errors. The likelihood ratio test is also not interpretable.")
        }

      ####### Interval table ########

        if(x$table == TRUE){
          pvalues <- as.numeric(table(intervaltally(x$p, x$steps)))
          cat("\n")
          cat("\n")
          cat("Number of Effect Sizes per Interval:")
          cat("\n")
          cat("\n")
          format(print.data.frame(sampletable(x$p, pvalues, x$steps)))
        }

      if(is.null(x$removed)==FALSE){
        cat("\n")
        cat("There were ", length(x$removed), "cases removed from your dataset due to the presence of missing data. To view the row numbers of these cases, use the attribute '$removed'.")
      }

    invisible()
}
