
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
  lapply(tested_simulations, validate_temulator_result_object)
  expect_invisible(simulateTumour(simulation_end_time=1e4))
})


test_that("reproducibility of simulations", {
  
  calculate_repoducible_checksum = 
    function(x) {
      class(x) = "list" # revert to the old class type
      x$mutation_data$id = NULL
      x$mutation_data$clone = NULL
      return(digest::sha1(x))
    }
  
  
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



test_that("simulations can be interupted", {
  
  params = 
    data.frame(
      birthrate=c(1,1.2),
      deathrate=c(0.2,0.2),
      mutationrate=c(1,1),
      start_time=c(0,50),
      father=c(0,0)
    )
  
  
  n_reactions = c(60, 100, 600, 10000)
  seed = 1
  n_clonal = 100
  
  sim_multistep = new(TEMULATOR_object, params, 0, 9, seed)
  
  for (i in seq_along(n_reactions)) {
    
    sim_multistep$end_time = n_reactions[i]    
    sim_multistep$run(FALSE)
    
    sim_at_once = new(TEMULATOR_object, params, n_reactions[i], n_clonal, seed)
    sim_at_once$run(FALSE)
    
    expect_equal(sim_multistep$cell_counts, sim_at_once$cell_counts)
  }
 
})


test_that("sample function", {
  
  params = 
    data.frame(
      birthrate=c(1,1.2),
      deathrate=c(0.2,0.2),
      mutationrate=c(1,1),
      start_time=c(0,50),
      father=c(0,0)
    )
  
  n_reactions = 1e5
  seed = 1
  n_clonal = 100
  
  sim = new(TEMULATOR_object, params, n_reactions, n_clonal, seed)
  sim$run(FALSE)

  # Invalid parameter min_vaf
  expect_error(sim$sample(-0.1, 1, 100, 1))
  expect_error(sim$sample(1.0, 1, 100, 1))
  expect_error(sim$sample(1.1, 1, 100, 1))
  expect_identical(class(sim$sample(0.5, 1, 100, 1)), "data.frame")
  
  
  # Invalid parameter purity
  expect_error(sim$sample(0.1, 1.1, 100, 1))
  expect_error(sim$sample(0.1, -0.1, 100, 1))
  expect_error(sim$sample(0.1, -0.0, 100, 1))
  expect_identical(class(sim$sample(0.1, 1.0, 100, 1)), "data.frame")
  expect_identical(class(sim$sample(0.1, 0.1, 100, 1)), "data.frame")
  
  
  # Invalid parameter depth
  expect_error(sim$sample(0.1, 1.0, -1.0, 1))
  expect_error(sim$sample(0.1, 1.0, 0.0, 1))
  expect_identical(class(sim$sample(0.1, 1.0, 100, 1)), "data.frame")

  
  # Invalid parameter purity
  expect_error(sim$sample(0.1, 1.0, 100, 0))
  expect_error(sim$sample(0.1, 1.0, 100, 100))
  expect_identical(class(sim$sample(0.1, 1.0, 100, 1)), "data.frame")
  expect_identical(class(sim$sample(0.1, 1.0, 100, 2)), "data.frame")
  expect_identical(class(sim$sample(0.1, 1.0, 100, 3)), "data.frame")
  
  rm(sim)
})
