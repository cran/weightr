#' Studies From Smith, Glass, and Miller's (1980) Meta-Analysis of Psychotherapy Outcomes
#'
#' An arbitrary subset of 74 studies from a meta-analysis assessing the effectiveness of psychotherapy. Contains two moderator variables.
#'
#' @docType data
#'
#' @usage dat.smith
#'
#' @format A data frame containing the following columns:
#' \describe{
#'  \item{\code{es}}{standardized mean difference effect sizes}
#'  \item{\code{v}}{corresponding sampling variances}
#'  \item{\code{age}}{continuous moderator representing clients' average age in years}
#'  \item{\code{diagnosis}}{categorical moderator representing disorder for which clients were treated; 1 = complex phobia, 2 = simple phobia, 3 = other}
#'  }
#'
#' @keywords datasets
#'
#' @details This dataset consists of an arbitrarily selected subset of 74 studies assessing the effectiveness of psychotherapy. Smith, Glass, and Miller (1980) published a meta-analysis designed to explore the current state of knowledge about psychotherapy effectiveness. Their original meta-analysis contains more than 1,700 effect sizes from 475 studies with multiple moderators and outcome measures. This subset is vastly simplified and intended solely for the purpose of demonstration.
#'
#' @references Smith, M. L., Glass, G. V., & Miller, T. I. (1980). Meta-analysis of psychotherapy. American Psychologist, 41, 165-180.
#'
#' @source Smith, M. L., Glass, G. V., & Miller, T. I. (1980). Meta-analysis of psychotherapy. American Psychologist, 41, 165-180.
#'
#' @examples
#' \dontrun{
#' dat.smith
#' effect <- dat.smith$es
#' v <- dat.smith$v
#' weightfunct(effect, v)
#' }
"dat.smith"
