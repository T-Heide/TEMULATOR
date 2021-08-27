
assign_mutation_label = function(d) {
  
  stopifnot(is.data.frame(d))
  stopifnot(c("clone","id") %in% colnames(d))
  
  # group of mutations that co-occure (i.e. true frequency clusters)
  f_group = base::strtoi(paste0("0x", gsub("X.*$", "", d$id)))
  f_group = f_group - min(f_group)
  
  
  # group of cell types in which a mutation occures (i.e. clone clusters)
  # 0 means multiple (i.e. selected)
  c_group = d$clone
  c_group = if_else(c_group == 0, 1, c_group - min(c_group[c_group != 0]) + 2)
  c_group[f_group == 0] = 0
  
  # f_ids of clonal cluster for each sub linage
  #f_group_clonal = tapply(f_group, c_group, min)
  
  label = as.character(c_group)
}



#' Validates a TEMULATOR result object
#'
#' @param x object of class 'temulator_result_object'.
#'
#' @return TRUE if successfull. Throws an error on failure.
validate_temulator_result_object = 
  function(x) {
    
    # returned structure is a named list
    # containing the different parameters and results:
    
    testthat::expect_identical(class(x), "temulator_result_object")
    
    testthat::expect_named(x)
    
    testthat::expect_identical(
      names(x)[1:6],
      c(
        "clone_parameters",
        "simulation_parameters",
        "sequencing_parameters",
        "cell_numbers",
        "simulation_data",
        "mutation_data"
      )
    )
    
    if (length(x) > 7) {
      testthat::expect_identical(names(x)[7], "mutation_subsets")
    }
    
    
    # 1) list element clone parameters
    
    n_clones = NROW(x$clone_parameters)
    
    testthat::expect_identical( # rownames are sequentialy numbered clone ids 
      rownames(x$clone_parameters), 
      paste0("clone_", seq_len(NROW(x$clone_parameters)))
    )
    
    testthat::expect_identical( # names parameters for each clone
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
    
    testthat::expect_vector(x$simulation_parameters)
    
    testthat::expect_identical(
      names(x$simulation_parameters), 
      c(
        "simulation_end_time", 
        "seed"
      )
    )
    
    
    # 3) list element sequencing parameters
    
    testthat::expect_vector(x$sequencing_parameters)
    
    testthat::expect_identical(
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
    
    testthat::expect_vector(x$cell_numbers)
    
    testthat::expect_length(x$cell_numbers, NROW(x$clone_parameters))
    
    testthat::expect_identical(
      names(x$cell_numbers),
      rownames(x$clone_parameters)
    )
    
    
    # 5) list element simulation data
    
    testthat::expect_vector(x$simulation_data)
    
    testthat::expect_identical(
      names(x$simulation_data),
      c(
        "time",
        "reactions"
      )
    )
    
    # basic mutation properties:
    testthat::expect_true(all(x$simulation_data>0, na.rm = TRUE))
    
    
    # 6) list element mutation data
    
    check_mdata = function(d) {
      
      # data frame
      testthat::expect_true(is.data.frame(d))
      testthat::expect_identical(colnames(d), c("clone","alt","depth","id","ccf"))

      # basic mutation properties:
      testthat::expect_true(all(d$alt > 0))
      testthat::expect_true(all(d$depth > 0))
      testthat::expect_true(all(d$alt <= d$depth))
      testthat::expect_true(all(between(d$ccf, 0, 1)))

    }
    
    if (!is.data.frame(x$sequencing_parameters)) { # old structure, one sequencing dataset
      check_mdata(x$mutation_data)
    } else {
      lapply(x$mutation_data, check_mdata)
    }
      
      
    if ("mutation_subsets" %in% names(x)) { # optional subsets of sequencing datasets
      for (i in seq_along(x$mutation_subsets)) {
        if (!is.data.frame(x$sequencing_parameters)) { # old structure, one sequencing dataset
          check_mdata(x$mutation_subsets[[i]])
        } else {
          lapply(x$mutation_subsets[[i]], check_mdata)
        }
      }
    }
    
  }



#' Print method for TEMULATOR result object
#'
#' @param x object of class 'temulator_result_object'.
#' @param ... unused.

#'
print.temulator_result_object = function(x, ...) {
  
  cat ("[ TEMULATOR result object ]\n\n")
  cat("> Clones:\n\n")
  
  for (i in seq_len(NROW(x$clone_parameters))) {
    
    cat("  Clone #", i, ":\n", sep="")
    cat("    - Birth rate:", x$clone_parameters$birthrates[i], "\n")
    cat("    - Death rate:", x$clone_parameters$deathrates[i],  "\n")
    cat("    - Mutation rate:", x$clone_parameters$mutation_rates[i], "\n")
    
    if (i > 1) {
      cat("    - Start time:", x$clone_parameters$clone_start_times[i],  "\n")
      cat("    - Accestor: Clone #", x$clone_parameters$fathers[i], "\n", sep="")
    }
    
    rn = rownames(x$clone_parameters)[i]
    if (rn %in% names(x$cell_numbers)) {
      frac = round(x$cell_numbers[rn] / sum(x$cell_numbers) * 100)
      cat("    => ", x$cell_numbers[i], " cells (", frac, "%)\n", sep="")
    }
    
    cat("\n")
  }
  
  cat("\n")
  cat("> Other:\n\n")
  cat(" - Seed:", x$simulation_parameters["seed"], "\n")
  cat(" -> Gillespie time:", x$simulation_parameters["time"], "\n")
  cat(" -> Total reactions:", x$simulation_data["reactions"], "\n")
  cat("\n\n")
  
  if (is.data.frame(x$sequencing_parameters)) {
    n_datasets = NROW(x$sequencing_parameters)
  } else {
    n_datasets = as.numeric(!is.null(x$mutation_data))
  }
  
  
  if (n_datasets) {
    cat("> Sequencing data:\n\n")
    
    if ("mutation_subsets" %in% names(x)) {
      n_subsets = length(x$mutation_subsets[[1]])
    } else {
      n_subsets = 0
    }

    for (j in seq_len(n_datasets)) {
      
      idx = (j - 1) * (n_subsets + 1) + 1
      cat("", j, ")\n\n", sep="")
      cat("  # Clonal mutations:", x$sequencing_parameters[["number_clonal_mutations"]][j], "\n")
      cat("  # Purity:", x$sequencing_parameters[["purity"]][j], "\n")
      cat("  # Depth:", x$sequencing_parameters[["depth"]][j], "\n")
      dm = as.character(x$sequencing_parameters[["depth_model"]][j])
      dm_label = c("1"="poisson distributed", "2"="overdispersed beta binomial", "3"="constant depth")[dm]
      cat("  # Depth model:", dm_label, "\n")
      cat("  # VAF cutoff:", x$sequencing_parameters[["min_vaf"]][j], "\n\n")
      print(get_sequencing_data(x, idx))
      cat("\n")
      cat("-> Call get_sequencing_data(x, idx=", idx, ") to retrieve these.\n\n", sep="")
      cat("\n")
      
      if ("mutation_subsets" %in% names(x)) {
        
        cat("=> Also", n_subsets, "subset(s):\n")
        for (k in seq_len(n_subsets)) {
          cat("     ", k, ") ", names(x$mutation_subsets[[j]])[k]," (idx=", idx + k, ")\n", sep="")
        }
        cat("\n")
        cat("-> Call get_sequencing_data(x, idx=i) to retrieve these.\n\n")
      }
    }
    
  }
 

  

  invisible(NULL)
}

#' Plot method for TEMULATOR result object
#'
#' @param x object of class 'temulator_result_object'.#' 
#' @param idx optional, index of the mutation dataset to plot (see output of print(x)).
#' @param quite optional, logical indicating if result should be plotted.
#' @param ... additional parameters passed to get_sequencing_data and assign_mutation_label.
#'
#' @return A ggplot object.
plot.temulator_result_object = function(x, idx=1, quite=FALSE, ...) {
  
  stopifnot(is.logical(quite))
  stopifnot(length(quite) == 1)
  
  # mutation data
  m_data = get_sequencing_data(x, idx=idx, ...)
  
  
  # index of row in sequencing_parameters
  if ("mutation_subsets" %in% names(x)) {
    n_subsets = length(x$mutation_subsets[[1]])
  } else {
    n_subsets = 0
  }
  idx_s = (idx-1) %/% (n_subsets+1) + 1
  
  
  # cell counts
  n_c = x$cell_numbers
  n_c_sum = sum(n_c)
  
  # clone data
  c_data = x$clone_parameters
  c_data = c_data[names(which(n_c > 0)),] # drop those without cells
  n_clones = NROW(c_data)
  
  
  # labels
  selection = sum(!duplicated(c_data[,c("birthrates","deathrates","mutation_rates")])) > 1
  sim_type = if_else(selection, "Non-neutral", "Neutral")
  main_label = paste0("TEMULATOR - ", sim_type, " simulation")
  
  sublabel_per_clone = paste0(round(n_c / n_c_sum, 2) * 100, "% C", seq_along(n_c), collapse=", ")
  sublabel = paste0(n_c_sum, " cells (", sublabel_per_clone, ")")
  
  
  ## clone parameter labels: 
  mr = paste0("(", paste0(c_data$mutation_rates, collapse=","), ")")
  br = paste0("(", paste0(c_data$birthrates, collapse=","), ")")
  dr = paste0("(", paste0(c_data$deathrates, collapse=","), ")")
  t = paste0("(", paste0(c_data$clone_start_times, collapse=","), ")")
  t_end = x$simulation_parameters["simulation_end_time"]
  C = x$sequencing_parameters[["depth"]][idx_s]
  p = x$sequencing_parameters[["purity"]][idx_s]
  min_vaf = x$sequencing_parameters[["min_vaf"]][idx_s]
  N_clonal = x$sequencing_parameters[["number_clonal_mutations"]][idx_s]
  caption = bquote(lambda==.(br)*","~mu==.(dr)*","~m==.(mr)*","~t==.(t)*";"~t[end]==.(t_end)*","~bar(C)==.(C)*","~N[clonal]==.(N_clonal))
  
    
  # Positions of subclonal clusters
  if (n_clones > 1) {
    cl_frac = c(clonal=1, n_c[-1] / n_c_sum)
  } else {
    cl_frac = c(clonal=1)
  }
  cl_locations = cl_frac * x$sequencing_parameters[["purity"]][idx_s] / 2
  cl_data = data.frame(clone=names(cl_locations), vaf=cl_locations)
  
  
  plot = 
    m_data %>% 
    ggplot(aes(x=vaf, fill=factor(label))) + 
      geom_histogram(breaks=seq(from=min_vaf, 1, by = 0.01)) + 
      geom_vline(data=cl_data, aes(xintercept=vaf), linetype=2, color="gray20") + 
      scale_fill_brewer(palette = "Set1") + 
      xlab("VAF") + 
      ylab("Number of mutations") + 
      labs(fill="") +
      ggtitle(label=main_label, subtitle=sublabel) +
      labs(caption = caption) + 
      theme(plot.caption=element_text(hjust=0)) + 
      xlim(0, max(m_data$vaf))

  
  # plot the figure and return it invisibliy
  if (!quite) plot(plot)
  
  invisible(plot)
}


get_sequencing_data = 
  function(x, ...) UseMethod("get_sequencing_data")

#' Gets sequencing data from a TEMULATOR result object
#'
#' @param x object of class 'temulator_result_object'.
#' @param idx optional index of the mutation dataset to return (see output of print(x)).
#' @param ... additional parameters passed to assign_mutation_label. 
#' 
#' @return tibble object
get_sequencing_data.temulator_result_object = 
  function(x, idx=1, ...) {
    
    if ("mutation_subsets" %in% names(x)) {
      n_subsets = length(x$mutation_subsets[[1]])
    } else {
      n_subsets = 0
    }
    
    if (!is.data.frame(x$sequencing_parameters)) { # old structure, one sequencing dataset
      x$mutation_data = list(x$mutation_data)      # convert to new structure
      if ("mutation_subsets" %in% names(x)) { 
        for (i in seq_along(x$mutation_subsets)) {
          x$mutation_subsets = list(x$mutation_subsets)
        }
      }
    }
    
    
    n_datasets = length(x$mutation_data)
    max_idx = (n_subsets + 1) * n_datasets
  
    if (!is.numeric(idx)) stop("Index not numeric.")
    if (!length(idx) == 1) stop("Index not of length 0.")
    if (!length(idx) <= max_idx & idx > 0) stop("Index out of range.")
    
    idx_a = (idx-1) %/% (n_subsets+1) + 1
    idx_b = (idx-1) %% (n_subsets+1)
    
    if (idx_b == 0) {
      mdata = x$mutation_data[[idx_a]]
    } else {
      mdata = x$mutation_subsets[[idx_a]][[idx_b]]
    }
      
    mdata %>% 
      mutate(vaf=alt/depth) %>% 
      mutate(label=assign_mutation_label(., ...)) %>% 
      select(alt, depth, vaf, id, label)
  }


get_clone_frequency = 
  function(x) UseMethod("get_clone_frequency")


#' Gets clone fractions from a TEMULATOR result object
#'
#' @param x object of class 'temulator_result_object'.
#'
#' @return a named numeric vector
get_clone_frequency.temulator_result_object =
  function(x) {
    x$cell_numbers/sum(x$cell_numbers)
  }
  