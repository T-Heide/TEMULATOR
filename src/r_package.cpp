#ifdef _DEBUG_
#define D(x) x
#else
#define D(x)
#endif

#include <stdlib.h>
#include "extern_global_variables.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <time.h>
#include "Phylogeny.h"
#include "Universe.h"
#include "Cell.h"
#include "CellType.h"
// [[Rcpp::depends(BH)]]
#include "Rcpp.h"

boost::random::mt19937_64 rng;         // produces randomness out of thin air

// [[Rcpp::export]]
Rcpp::List SimulateTumor(std::vector< double > birthrates,
                         std::vector< double > deathrates,
                         std::vector< double > mutation_rates,
                         std::vector< unsigned int > clone_start_times,
                         std::vector< int > fathers,
                         int simulation_end_time,
                         int seed,
                         int number_clonal_mutations,
                         double purity,
                         double min_vaf,
                         int depth,
                         int depth_model,
                         bool verbose)
{

  SimulationParameterSet sim_params(mutation_rates,
                                    birthrates,
                                    deathrates,
                                    clone_start_times,
                                    fathers,
                                    simulation_end_time,
                                    seed,
                                    number_clonal_mutations,
                                    purity,
                                    depth,
                                    min_vaf,
                                    depth_model,
                                    "");

  Universe universe(sim_params);
  universe.RunSimulation(verbose);
  if (verbose) universe.PrintCellTypes();

  // Sampling:
  std::vector <int> clone_mutation;
  std::vector <int> alt_mutation;
  std::vector <int> depth_mutation;
  std::vector <std::string> id_mutations;

  universe.Sample(clone_mutation, alt_mutation, depth_mutation, id_mutations);


  // Copy vectors to output R-Vector:
  Rcpp::IntegerVector out_clone(clone_mutation.size());
  Rcpp::IntegerVector out_alt(alt_mutation.size());
  Rcpp::IntegerVector out_depth(depth_mutation.size());
  Rcpp::CharacterVector out_ids(id_mutations.size());

  for(std::vector<double>::size_type i = 0; i < clone_mutation.size(); i++) {
    out_clone[i] = clone_mutation[i];
    out_alt[i] = alt_mutation[i];
    out_depth[i] = depth_mutation[i];
    out_ids[i] = id_mutations[i];
  }


  // Get total cell number:
  std::vector <unsigned int> cells_per_clone = universe.CellCounts();
  Rcpp::NumericVector out_cells_per_clone(cells_per_clone.size());
  for (int i = 0; i < cells_per_clone.size(); i++) {// for each clone seperatly:
    out_cells_per_clone[i] += cells_per_clone[i];
  }

  // Simulation parameters:
  Rcpp::NumericVector sim_data = 
    Rcpp::NumericVector::create(
      Rcpp::Named("time", universe.Time()),
      Rcpp::Named("reactions", universe.Reactions())
    );
  

  // Return list of results:
  return Rcpp::List::create(Rcpp::Named("cell_numbers") = out_cells_per_clone,
                            Rcpp::Named("simulation_data") = sim_data,
                            Rcpp::Named("mutation_data") =
                                Rcpp::DataFrame::create(
                                    Rcpp::Named("clone") = out_clone,
                                    Rcpp::Named("alt") = out_alt,
                                    Rcpp::Named("depth") = out_depth,
                                    Rcpp::Named("id") = out_ids));
}
