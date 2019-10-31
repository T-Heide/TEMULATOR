#' simulateTumour
#'
#' @param birthrates
#' @param deathrates
#' @param mutation_rates
#' @param clone_start_times
#' @param fathers
#' @param simulation_end_time
#' @param seed
#' @param number_clonal_mutations
#' @param purity
#' @param min_vaf
#' @param depth
#' @param depth_model Number specificing the depth model to use (1:  poisson distributed depth model, 2: overdispersed beta binomial (default), 3: fixed depth) 
#'
#' @return
#' @export
#'
#' @examples
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
