/*
    Copyright (C) 2018 Timon Heide (timon.heide@icr.ac.uk)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "extern_global_variables.h"
#include "CellType.h"
#include "Cell.h"
#include "Universe.h"
#include "Phylogeny.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <array>
#include <vector>
#include <fstream>
#include <cmath>
#include <boost/random/uniform_real.hpp>
#include "Rcpp.h"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#define TICKSOURCE generation

// Universe ////////////////////////////////////////////////////////////////////


Universe::Universe(
    Rcpp::DataFrame clone_params,
    unsigned int end_time,
    unsigned int clonal_mutations,
    int seed
    )
  : 
  mNumberOfClones(0),
  mNumberOfCells(0), 
  mNumberOfReactions(0),
  mTime(0.0), 
  
  mSimulationEndTime(end_time),
  mNextClone(0),
  mSeed(seed),
  
  mMutationrates(Rcpp::as<std::vector<double> >(clone_params["mutationrate"])),
  mBirthrates(Rcpp::as<std::vector<double> >(clone_params["birthrate"])),
  mDeathrates(Rcpp::as<std::vector<double> >(clone_params["deathrate"])),
  mCloneStartTimes(Rcpp::as<std::vector<unsigned int> >(clone_params["start_time"])),
  mFathers(Rcpp::as<std::vector<unsigned int> >(clone_params["father"])),
  mClonalMutations(clonal_mutations)
{
  
  mNumberOfClones = static_cast <int> (mMutationrates.size());
  
  
  // Basic checks of the input:
  double last_start_time = 0.0;
  
  for (int i = 0; i < mNumberOfClones; i++) {
    
    // Ensure that the clone start times are ordered:
    if (mCloneStartTimes[i] < last_start_time) {
      Rcpp::stop("Defined clone start times have to be sorted!");
    }
    
    last_start_time = mCloneStartTimes[i];
    
    // Ensure that fathers are defined:
    if (i == 0 && mFathers[i] != 0) {
      Rcpp::stop("Fathers of first clone has to be 0.");
    }
    
    if (i > 0 && mFathers[i] >= i) {
      Rcpp::stop("Fathers of clones have to be introduced first!");
    }
  }
}


// Destructor:
Universe::~Universe(){

  // Remove all CellTypes:
  for(std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++) {
    delete mpTypes[i++];
  }

  // Remove the complete phylogeny:
  for(std::vector<PhylogenyRoot*>::size_type i = 0; i < mpPhylogenies.size();) {
    delete mpPhylogenies[i++];
  }
}


// Simulate tumour:
bool Universe::RunSimulation(bool verbose)
{

  // Create reusable pointers to handled cells, types
  // and other variables used for simulations:
  Cell* pCell;
  CellType* pType;
  long double dt = 0.0L;
  int action = -1;

  

  if (mNumberOfReactions == 0) { // Not started simulation
    
    if (verbose) Rcpp::Rcout << "Starting stimulation:" << std::endl;
    
    // Set random seed
    rng.seed(mSeed);
    
    
    // Create and insert a cell of the first type:
    pType = new CellType(mBirthrates[0], mDeathrates[0], mMutationrates[0]);
    pCell = new Cell(pType);
    this->InsertCell(pCell);
    pCell->MutateCellFixedNumber(mClonalMutations);
    
    
    // universe time to start time of fist clone:
    // Debug messages:
    /*if (verbose) {
      Rcpp::Rcout << std::endl;
      Rcpp::Rcout << "########## Comment #############" << std::endl;
      Rcpp::Rcout << " New clone: " << mNextClone << std::endl;
      Rcpp::Rcout << " Time: " << mNumberOfReactions << std::endl;
      Rcpp::Rcout << "################################" << std::endl;
      Rcpp::Rcout << std::endl;
    }*/
    
    mNextClone++;
    
  } else {
    if (verbose) Rcpp::Rcout << "Restarting stimulation:" << std::endl;
  }
  

  // Run till limit of universe has been reached:
  Progress p(mSimulationEndTime-mNumberOfReactions, verbose);
  
  while (mNumberOfReactions < mSimulationEndTime) {
    
    mNumberOfReactions++;
    
    // Check if new clones need to be introduced in each cycle:
    for (; mNextClone < mNumberOfClones &&
           mCloneStartTimes[mNextClone] <= mNumberOfReactions;
           mNextClone++)
    {
      // Introduce a new clone:
      int c_father = mFathers[mNextClone];
      pType = new CellType(mBirthrates[mNextClone], mDeathrates[mNextClone], mMutationrates[mNextClone]);
      mpTypes[c_father]->RandomMember()->Type(pType);

      // Debug messages:
      /*if (verbose) {
        Rcpp::Rcout << std::endl;
        Rcpp::Rcout << "########## Comment #############" << std::endl;
        Rcpp::Rcout << " New clone: " << mNextClone << std::endl;
        Rcpp::Rcout << " Time: " << mNumberOfReactions << std::endl;
        Rcpp::Rcout << "################################" << std::endl;
        Rcpp::Rcout << std::endl;
      }*/
      
    } // end of for loop for introduction of new clones

    // Select reaction type and members:
    this->NextReaction(&dt, &action)->RandomMember()->DoAction(&action);
    this->IncrementTimeBy(dt);
    
    p.increment(); 
    if (Progress::check_abort()) {
      return false;
    }
    
  } // stop running after reaching limit

  return true;
}


// Getter functions:
double Universe::Time() const {return mTime;}

unsigned int Universe::Reactions() const {return mNumberOfReactions;}

CellType* Universe::NextReaction(long double* r_delta_time, int* action) const{
  // Samples and returns the next reaction that occures.

  // Exit if the universe contains no types:
  if (mpTypes.size() == 0) {
    Rcpp::stop("In Universe::NextReaction: No type to choose.");
  }

  // Keep track of minum index and dt:
  double delta_time_minimum;
  int index_minimum;

  boost::uniform_real<double> runi_dist(0.0, 1.0);

  // Calculate dt for all types in the universe:
  for(std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++) {

    long double runi = runi_dist(rng);
    long double c_birthrate = mpTypes[i]->Birthrate();
    long double c_number = mpTypes[i]->NumMembers();
    long double delta_time_current = (log(1) - log(runi)) / (c_number * c_birthrate);

    // Keep track of minum index and dt:
    if (i == 0 || delta_time_minimum > delta_time_current) {
      index_minimum = i;
      delta_time_minimum = delta_time_current;
    }
  }

  // Return results:
  *r_delta_time = delta_time_minimum;
  *action = 1; // action divide

  return(mpTypes[index_minimum]);
}


Rcpp::DataFrame Universe::SampleRcpp(double min_vaf, double purity, double depth, int depth_model) const {
  
  // Check input:
  if (min_vaf>=1 || min_vaf<0) {
    Rcpp::stop("Argument min_vaf should be in the interval [0,1).");
  }
  
  if (purity>1 || purity<=0) {
    Rcpp::stop("Argument purity should be in the interval (0,1].");
  }
  
  if (depth<=0) {
    Rcpp::stop("Argument depth should be in the interval (0,Inf].");
  }
  
  
  // Sampling:
  std::vector <int> clone_mutation;
  std::vector <int> alt_mutation;
  std::vector <int> depth_mutation;
  std::vector <std::string> id_mutations;
  
  this->Sample(clone_mutation, alt_mutation, depth_mutation, id_mutations, min_vaf, purity, depth, depth_model);
  
  
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
  
  
  return
    Rcpp::DataFrame::create(
        Rcpp::Named("clone") = out_clone,
        Rcpp::Named("alt") = out_alt,
        Rcpp::Named("depth") = out_depth,
        Rcpp::Named("id") = out_ids
    );
  ;
}


bool Universe::Sample(std::vector <int> &clone,
                      std::vector <int> &alt,
                      std::vector <int> &depth,
                      std::vector <std::string> &ids,
                      double min_vaf, 
                      double purity, 
                      double avg_depth, 
                      int depth_model) const
{

  // Determine total number of cells:
  int ncells = 0;
  for (std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++){
    ncells += mpTypes[i]->NumMembers();
  }

  // Traverse the Phylogenies:
  int ncells2 = 0;
  for (std::vector<CellType*>::size_type i = 0; i < mpPhylogenies.size(); i++){
    ncells2 += mpPhylogenies[i]->Root()->SampleNode(clone, alt, depth, ids,
                                                    min_vaf, 
                                                    purity, 
                                                    avg_depth, 
                                                    depth_model,
                                                    ncells);
  }

  if (ncells != ncells2) {
    Rcpp::Rcerr << "sum_cells != total_cells: " << ncells << " vs " << ncells2 << std::endl;
  }

  return true;
}


std::vector <unsigned int> Universe::CellCounts() const {

  // Get total count:
  std::vector <unsigned int> cells_in_clone;
  for(std::vector<CellType *>::size_type i = 0; i < mpTypes.size(); i++) {
    cells_in_clone.push_back(mpTypes[i]->NumMembers());
  }

  return(cells_in_clone);
}

unsigned int Universe::get_SimulationEndTime() const {return mSimulationEndTime;}


// Setter functions:
void Universe::IncrementTimeBy(long double Delta){ mTime += Delta; }

bool Universe::InsertCell(Cell* pCell) {
  return InsertCell(pCell, true);
}

bool Universe::InsertCell(Cell* pCell, bool is_new_lineage) {

  // Insert type into universe
  if (is_new_lineage) {
    this->RegisterType(pCell->Type());
  }

  // Set universe variable in cell:
  pCell->AssociatedUniverse(this);

  // Update universe:
  mNumberOfCells++;

  // Register new phylogeny:
  if (pCell->AssociatedNode() == 0) {
    PhylogenyRoot* pNewPhylo = new PhylogenyRoot(pCell);
    mpPhylogenies.push_back(pNewPhylo);
  }

  return true;
}

void Universe::RegisterType(CellType* p_new_type){
  std::vector<CellType *>::size_type i = 0;

  while(i < mpTypes.size()) {
    if (mpTypes[i] == p_new_type)
      break;
    i++;
  }

  if (i == mpTypes.size()) { // reached end.
    mpTypes.push_back(p_new_type);
  }
}


void Universe::set_SimulationEndTime(unsigned int time){
  
  if (time <= mNumberOfReactions) {
    Rcpp::stop("New time has to be larger than current time.");
  }
  
  mSimulationEndTime = time;
  
}




// Output Functions:
void Universe::Print() const {
  
  // Print arguments:
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "########## Options #############" << std::endl;
  // clone level data
  for (int i = 0; i < mNumberOfClones; i++) { // for each clone seperatly:
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "  Mutation rate: " << mMutationrates[i] << std::endl;
    Rcpp::Rcout << "  Birth rate: " << mBirthrates[i] << std::endl;
    Rcpp::Rcout << "  Death rate: " << mDeathrates[i] << std::endl;
    Rcpp::Rcout << "  Clone start time: " << mCloneStartTimes[i] << std::endl;
    Rcpp::Rcout << "  Father: " << mFathers[i] << std::endl;
  } // end printing of clone level data
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "  Simulation end time: " << mSimulationEndTime << std::endl;
  Rcpp::Rcout << "  Clonal mutations: " << mClonalMutations << "\n";
  Rcpp::Rcout << "  Random seed: " << mSeed << std::endl;
  Rcpp::Rcout << "################################" << std::endl;
  Rcpp::Rcout << std::endl;
}

void Universe::PrintCellTypes() const {

  // Get total count:
  unsigned long sum_cells = 0;
  for(std::vector<CellType *>::size_type i = 0; i < mpTypes.size(); i++) {
    sum_cells += mpTypes[i]->NumMembers();
  }

  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "########## Result ##############" << std::endl;
  Rcpp::Rcout << "Total number of cells: " << sum_cells << std::endl;

  for(std::vector<CellType *>::size_type i = 0; i < mpTypes.size(); i++) {
    int cells = mpTypes[i]->NumMembers();
    Rcpp::Rcout << "  Clone " << i << ": " << cells
              << " (" << std::setprecision(3) << cells*100.0/sum_cells << "%)"
              << std::endl;
  }

  Rcpp::Rcout << "################################" << std::endl;
  Rcpp::Rcout << std::endl;
}

