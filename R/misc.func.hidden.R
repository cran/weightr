## Note to User:
## This file contains code for functions that are useful in calculating
## or formatting model output. They operate "behind the scenes."

## Function to tally the number of effect sizes per p-value interval ##

intervaltally <- function(p, steps) {
  p1 <- cut(p, breaks=c(-Inf,steps), labels=steps)
  return(p1) }

## Function to format the above tally into a table and add labels ##

sampletable <- function(p, pvalues, steps){
  nsteps <- length(steps)
  results <- matrix(nrow=length(pvalues),ncol=1)
  results[,1] <- pvalues
  rowlabels <- c(0, length(results[,1]))
  rowlabels[1] <- paste(c("p-values <", steps[1]), collapse="")
  for(i in 2:nsteps){
    rowlabels[i] <- paste(c(steps[i - 1], "< p-values <", steps[i]), collapse=" ")
  }
  resultsb <- data.frame(results, row.names=c(rowlabels))
  colnames(resultsb) <- c("Frequency")
  return(resultsb)
}

##
