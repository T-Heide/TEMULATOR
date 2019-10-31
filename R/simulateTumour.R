#' Main function used to produce a simulation
#'
#' @param birthrates Birth rate(s) of subpopulations
#' @param deathrates Death rate(s) of subpopulations
#' @param mutation_rates Mutation rate(s) per division
#' @param clone_start_times Timepoint (in number of reaction) a subpopulation is introduced
#' @param fathers Subpopulation (0 based) a subpopulations originate from
#' @param simulation_end_time Number of reactions performed
#' @param seed Simulation seed
#' @param number_clonal_mutations Number of clonal variants
#' @param purity Assumed purity of the tumour
#' @param min_vaf Minimum VAF to report
#' @param depth Sequencing depth
#' @param depth_model Number specificing the distribution of the sequencing depth to use (1:  poisson, 2: overdispersed beta binomial (default), 3: fixed) 
#' @param verbose Print progress?
#'
#' @return List containing i) properties of the subclones ('clone_parameters'), ii) the main simulation parameters ('simulation_parameters'), iii) the sequencing simulation parameters ('sequencing_parameters'), iv) the final population structure ('cell_numbers'), and v) the simulated mutation data.
#' @export
#'
#' @examples 
#' simulateTumour()
simulateTumour =
  function(birthrates=c(1.0, 1.0),
           deathrates=c(0.2, 0.2),
           mutation_rates=c(16, 16),
           clone_start_times=c(0, 256),
           fathers=c(0,0),
           simulation_end_time=1048576,
           seed=42,
           number_clonal_mutations=100,
           purity=1.0,
           min_vaf=0.01,
           depth=100,
           depth_model=2,
           verbose=FALSE) {


  # Create structures containing simulation parameters:
  clone_params = list(birthrates=birthrates,
                      deathrates=deathrates,
                      mutation_rates=mutation_rates,
                      clone_start_times=clone_start_times,
                      fathers=fathers)

  sim_params = c(simulation_end_time=simulation_end_time,
                 seed=seed)

  seq_params = c(number_clonal_mutations=number_clonal_mutations,
                 purity=purity,
                 min_vaf=min_vaf,
                 depth=depth,
                 depth_model=depth_model)


  # Test inputs:
  if (length(unique(sapply(clone_params, length))) != 1) {
    stop("Clone parameters must be of equal length.\n")
  }

  # ToDo:... additional tests of inputs needed!


  # Call simulation:
  params = c(clone_params,  as.list(sim_params), as.list(seq_params),  verbose=verbose)
  sim_results = do.call(SimulateTumor, params)


  # Modify objects for return:

  clone_params = as.data.frame(clone_params)
  rownames(clone_params) = paste0("clone_", seq_len(NROW(clone_params)))

  params = list(clone_parameters=clone_params,
                simulation_parameters=sim_params,
                sequencing_parameters=seq_params)

  names(sim_results$cell_numbers) = rownames(clone_params)

  sim_results$mutation_data = tibble::as_tibble(sim_results$mutation_data)

  sim_results = c(params, sim_results)

  invisible(sim_results)
}
