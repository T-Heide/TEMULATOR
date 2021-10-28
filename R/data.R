#' 150 synthetic sequencing datasets from MOBSTER
#'
#' A list of 150 synthetic sequencing data ('temulator_result_object')
#' used in Caravagna et al. (2020).
#'
#' @format A names list of 150 temulator_result_objects.
#' @usage data(mobster_simulations)
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
#' @usage data(mobster_summary_of_simulations)
"mobster_summary_of_simulations"


#' Simulated sequencing data of a sweeping clone at various time points
#'
#' Details can be found in the "Time series data" vignette. 
#'
#' @format A data frame with 14850 rows and 7 variables:
#' \describe{
#'   \item{clone}{Cell clone the mutation is part of.}
#'   \item{alt}{Simulated number of variant reads.}
#'   \item{depth}{Simulated coverage of mutated site.}
#'   \item{id}{Mutation identifier.}
#'   \item{vaf}{Simulated VAF of the mutation.}
#'   \item{n_reactions}{Total number of reactions at present time.}
#'   \item{time}{Gillespie time.}
#' }
#' @usage data(samples_selection_sweep)
"samples_selection_sweep"


#' Clone frequencies of a sweeping clone at various time points
#' 
#' Details can be found in the "Time series data" vignette. 
#' 
#' @format A data frame with 50 rows and 7 variables:
#' \describe{
#'   \item{reactions}{Total number of reactions at present time.}
#'   \item{time}{Gillespie time.}
#'   \item{cells}{Total number of cells.}
#'   \item{clone1}{Number of cells in the ancestral clone.}
#'   \item{clone2}{Number of cells in the subclone.}
#'   \item{seed}{Random seed used to generate the simulation.}
#'   \item{start_time}{Time point (reactions) at which the subclone was added.}
#' }
#' @usage data(clone_f_selection_sweep)
"clone_f_selection_sweep"


#' Example of a time series dataset for a simulation used in the MOBSTER paper.
#
#' @format A data frame with 1100 rows and 5 variables:
#' \describe{
#'   \item{reactions}{number of reactions at time point}
#'   \item{t}{gillespie time at time point}
#'   \item{cells_a}{number of cells for the ancestral clone  at time point}
#'   \item{cells_sc}{number of cells for the (selected) subclone  at time point}
#'   \item{cells}{total number of cells at time point}
#' }
#' @usage data(time_series_data_example)
"time_series_data_example"

