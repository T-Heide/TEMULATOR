
tested_parameter_sets = 
  list(
    "f69c2d6fb3f0383bf4d3254af47ae59d3414c091"=
      list(seed=1, simulation_end_time=1e4),
    
    "622e40df34a87644e9912660ad35414f0822bc63"=
      list(seed=2, simulation_end_time=1e4),
    
    "6aff59db789ac36422afb7d18bf703c659da2370"=
      list(seed=2, simulation_end_time=1e5),
    
    "29379e8bf3ebbf5fe0b872e085600b7adea97ce8"=
      list(seed=2, simulation_end_time=1e4,
           birthrates=c(1,1.2), 
           deathrates=c(0,0), 
           mutation_rates=c(16,16),
           clone_start_times=c(0,100),
           fathers=c(0,0)),
    
    "d15c9cbc59b05f2411d1a19e742669c5f8b06e39"=
      list(seed=2, simulation_end_time=1e4,
           birthrates=c(1,1.2), 
           deathrates=c(0,0.1), 
           mutation_rates=c(12,16),
           clone_start_times=c(0,100),
           fathers=c(0,0))
  )


tested_simulations = 
  tested_parameter_sets %>%
  lapply(do.call, what=simulateTumour)



test_that("object structure", {

  test_object_structure = 
    function(x) {
      
      # returned structure is a named list
      # containing the different parameters and results:
      
      expect_identical(class(x), "list")
      
      expect_named(x)
      
      expect_identical(
        names(x),
        c(
          "clone_parameters",
          "simulation_parameters",
          "sequencing_parameters",
          "cell_numbers",
          "simulation_data",
          "mutation_data"
        )
      )
      
      
      # 1) list element clone parameters
      
      n_clones = NROW(x$clone_parameters)
      
      expect_identical( # rownames are sequentialy numbered clone ids 
        rownames(x$clone_parameters), 
        paste0("clone_", seq_len(NROW(x$clone_parameters)))
      )
      
      expect_identical( # names parameters for each clone
        colnames(x$clone_parameters),
        c(
          "birthrates",
          "deathrates",
          "mutation_rates",
          "clone_start_times",
          "fathers"
        )
      )
      
      
      # 2) list element simulation parameters
      
      expect_vector(x$simulation_parameters)
      
      expect_identical(
        names(x$simulation_parameters), 
        c(
          "simulation_end_time", 
          "seed"
          )
      )
      
      
      # 3) list element sequencing parameters
      
      expect_vector(x$sequencing_parameters)
      
      expect_identical(
        names(x$sequencing_parameters),
        c(
          "number_clonal_mutations",
          "purity",
          "min_vaf",
          "depth",
          "depth_model"
        )
      )
      
      
      # 4) list element cell numbers
      
      expect_vector(x$cell_numbers)
      
      expect_length(x$cell_numbers, NROW(x$clone_parameters))
      
      expect_identical(
        names(x$cell_numbers),
        rownames(x$clone_parameters)
      )
        
      
      # 5) list element simulation data
      
      expect_vector(x$simulation_data)
      
      expect_identical(
        names(x$simulation_data),
        c(
          "time",
          "reactions"
        )
      )
      
      # basic mutation properties:
      expect_true(all(x$simulation_data>0))
      
      
      # 6) list element mutation data
      
      expect_true(is.data.frame(x$mutation_data))
      
      expect_identical(
        colnames(x$mutation_data),
        c(
          "clone",
          "alt",
          "depth",
          "id"
        )
      )
      
      # basic mutation properties:
      expect_true(all(x$mutation_data$alt>0))
      expect_true(all(x$mutation_data$depth>0))
      expect_true(all(x$mutation_data$alt<=x$mutation_data$depth))
      
    }
  
  lapply(tested_simulations, test_object_structure)
  
  expect_invisible(simulateTumour())
  
})


test_that("reproducibility of simulations", {
  
  test_reproducibility = 
    function(params) {
      
      sim_results =
        lapply(1:4, function(i) # 4 repeated simulations
          do.call(simulateTumour, params)) %>%
        
        lapply(function(x) {
          x$mutation_data$id = NULL
          x$mutation_data$clone = NULL
          return(x)
        }) # drop mutation ids these differ ...
      
      
      for (i in seq_along(sim_results)) { # compare all against the first 
        expect_true(all.equal(sim_results[[1]], sim_results[[i]]))
      }
      
    }
  
  calculate_repoducible_checksum = 
    function(x) {
      x$mutation_data$id = NULL
      x$mutation_data$clone = NULL
      return(digest::sha1(x))
    }
  
  lapply(tested_parameter_sets, test_reproducibility)
  
  checksums = sapply(tested_simulations, calculate_repoducible_checksum)
  expect_equal(as.character(checksums), names(checksums))
  
})


test_that("seeding works", {
  
  test_seeding = 
    function(params) {
      
      if ("seed" %in% names(params)) {
        params[["seed"]] = NULL 
      }
      
      sim_results =
        lapply(1:4, function(i) # 4 repeated simulations
          do.call(simulateTumour, c(params, seed=i))) %>%
        
        lapply(function(x) {
          x$mutation_data$id = NULL
          x$mutation_data$clone = NULL
          return(x)
        }) # drop mutation ids these differ ...
      
      
      # compare all against the first 
      compared_elements = c("cell_numbers","mutation_data")
      ref = sim_results[[1]]$cell_numbers[compared_elements]
      for (i in seq_along(sim_results)[-1]) { 
        expect_false(isTRUE(all.equal(ref, sim_results[[i]][compared_elements])))
      }
      
    }
  
  lapply(tested_parameter_sets, test_seeding)
})



test_that("results are reasonable", {
  
  params = list(
    birthrates = 1,
    deathrates = 0,
    fathers = 0,
    mutation_rates = 10,
    clone_start_times = 0,
    depth_model = 3
  )
  
  for (n in c(10,100,500,1000,5000,10000)) {
    sim = do.call(simulateTumour, c(params, list(simulation_end_time=n)))
    expect_identical(sim$cell_numbers, c("clone_1"=n+1))
  }
  
})
