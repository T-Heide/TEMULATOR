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
#' @param depth_model Number specificing the distribution of the sequencing depth to use (0: no sequencing, 1:  poisson, 2: overdispersed beta binomial (default), 3: fixed) 
#' @param verbose Print progress?
#' @param subset_fractions Optional, Numeric vector specificing fraction variants are subset to.
#'
#' @return A temulator_result_object.
#' @export
#'
#' @examples 
#' simulateTumour(simulation_end_time=1000)
#' simulateTumour(subset_fractions=c("WES"=0.03), simulation_end_time=1000) 
#' simulateTumour(depth=c(30,50,100), simulation_end_time=1000) # three different sequencing depths  
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
           verbose=FALSE,
           subset_fractions=numeric()) {

  stopifnot(is.vector(subset_fractions))
  stopifnot(is.numeric(subset_fractions))
  stopifnot(all(subset_fractions > 0 & subset_fractions < 1))

  n_seq_parms = list(purity, min_vaf, depth, depth_model) %>% sapply(length)
  stopifnot(all(n_seq_parms == max(n_seq_parms) | n_seq_parms == 1))

  # Put parameters into structures, that are later returned
  sim_params = 
    c(
      simulation_end_time=simulation_end_time,
      seed=seed
     )
    
  clone_params = 
    data.frame(
      birthrate=birthrates, 
      deathrate=deathrates, 
      mutationrate=mutation_rates,
      start_time=clone_start_times,
      father=fathers,
      row.names=paste0("clone_", seq_along(mutation_rates))
    )
  
  seq_params =
    data.frame(
      number_clonal_mutations=number_clonal_mutations,
      purity=purity,
      min_vaf=min_vaf,
      depth=depth,
      depth_model=depth_model
     ) %>% dplyr::filter(depth_model != 0)
  
  
  # Create simulation object
  simulation = 
    new(
      TEMULATOR_object,
        clone_params, 
        simulation_end_time, 
        number_clonal_mutations,
        seed
    )
  
  # run the simulation, this might take a while
  success = simulation$run(verbose)
  if (!success) stop("Simulation failed")

  
  # Print results if verbose
  if (verbose) simulation$print()
  if (verbose) simulation$print_cell_types()
  
  
  n_cells = simulation$cell_counts; names(n_cells) = rownames(clone_params)
  sim_data = c(time=simulation$simulation_time, reactions=simulation$n_reactions)
  
  # Sample the results:
  mutation_data = list()
  for (i in seq_len(NROW(seq_params))) {
    
    mutation_data[[i]] =
      simulation$sample(
        seq_params$min_vaf[i],
        seq_params$purity[i],
        seq_params$depth[i],
        seq_params$depth_model[i]
      ) %>% as_tibble()
      
  }
  
  rm(simulation)
  
  
  # for consistency with old version change colnames:
  colnames(clone_params) = c("birthrates",
                             "deathrates",
                             "mutation_rates",
                             "clone_start_times",
                             "fathers")
  
  
  result_object = 
    list(
      clone_parameters = clone_params,
      simulation_parameters = sim_params,
      sequencing_parameters = seq_params,
      cell_numbers = n_cells,
      simulation_data = sim_data,
      mutation_data = mutation_data
    )
  
  class(result_object) = "temulator_result_object"
  
  
  
  # optionally subset the variant set
  if (length(subset_fractions)) {
    
    for (i in seq_along(mutation_data)) {
      
      mut_q = # mutation quantile (hash of variant id, scaled to (0,1])
        as.character(result_object$mutation_data[[i]]$id) %>% 
        sapply(digest::digest, algo="xxhash32", raw=TRUE, seed=seed) %>% 
        (function(x) as.numeric(paste0("0x", x)) / (2^32))
      
      mutation_subsets = 
        lapply(subset_fractions, function(x) { 
          result_object$mutation_data[[i]][mut_q < x,]
        })
      
      if (is.null(names(mutation_subsets))) { # set names if non given
        names(mutation_subsets) = as.character(subset_fractions)
      }
      
      result_object[["mutation_subsets"]][[i]] = mutation_subsets
      
    }
  }
  
  
  # revert to old structure if only one sample taken
  if (NROW(result_object$sequencing_parameters) == 1) {
    
    result_object$sequencing_parameters = 
      unlist(result_object$sequencing_parameters)
    
    result_object$mutation_data = 
      result_object$mutation_data[[1]]
  }
  
  invisible(result_object)
}

