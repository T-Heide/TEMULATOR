// [[Rcpp::depends(BH)]]
#include "Rcpp.h"
#include "Universe.h"
#include "extern_global_variables.h"

boost::random::mt19937_64 rng;


//' @name TEMULATOR_object
//' @title A TEMULATER simulation
//' @description Type the name of the class to see its methods
//' @field print Print the object
//' @field print_cell_types Print contained cell types
//' @field run Start or restart the simulation \itemize{
//' \item Parameter: verbose - bool, Flag indicating if progress should be printed.
//' \item Returns: bool, Indicating success.
//' }
//' @field sample Sample mutations from the population \itemize{
//' \item Parameter: min_vaf - double, Minimum VAF reported.
//' \item Parameter: purity - double, Assumed purity of the tumour.
//' \item Parameter: depth - double, Sequencing depth.
//' \item Parameter: depth_model - int,  Distribution of the sequencing depth (1: poisson, 2: overdispersed beta binomial (default), 3: fixed).
//' \item Returns: data.frame containing mutation information.
//' }
RCPP_MODULE(temulator_module) {
  using namespace Rcpp;
  
  class_<Universe>("TEMULATOR_object")
    
  .constructor<Rcpp::DataFrame, unsigned int, unsigned int, int>()
    
  .method("print", &Universe::Print)
  .method("print_cell_types", &Universe::PrintCellTypes)
  .method("show", &Universe::Print)
    
  .method("run", &Universe::RunSimulation)
  .method("sample", &Universe::SampleRcpp)
    
  .property("cell_counts", &Universe::CellCounts)
  .property("simulation_time", &Universe::Time)
  .property("n_reactions", &Universe::Reactions)
  .property("end_time", &Universe::get_SimulationEndTime, &Universe::set_SimulationEndTime)
    
  ;
}
