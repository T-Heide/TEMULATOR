#' 150 synthetic sequencing datasets from MOBSTER
#'
#' A list of 150 synthetic sequencing data ('temulator_result_object')
#' used in Caravagna et al. (2020).
#'
#' @format A names list of 150 temulator_result_objects.
#' 
"mobster_simulations"


#' Summary statistics of synthetic sequencing datasets used in MOBSTER
#
#' @format A data frame with 1755 rows and 8 variables:
#' \describe{
#'   \item{subclone_start}{insertion time of the subclone}
#'   \item{subclone_birthrate}{birth rate of the subclone}
#'   \item{n}{simulation number (9 per parameter combination)}
#'   \item{seed}{seed used for the simulation}
#'   \item{n_mutations_before_insertions}{number of mutations that occured prior to mutation}
#'   \item{subclone_fraction}{fraction of all cells that subclone made up}
#'   \item{assigned_id}{id used as identifier in the paper}
#' }
"mobster_summary_of_simulations"