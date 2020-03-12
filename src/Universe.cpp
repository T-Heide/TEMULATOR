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

#define TICKSOURCE generation

// Universe ////////////////////////////////////////////////////////////////////


// Constructor:
Universe::Universe(SimulationParameterSet params)
  : mTime(0.0), mNumberOfCells(0), mNumberOfReactions(0), mParameters(params) {}


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

  if (verbose) { // Print arguments:
    mParameters.Print();
  } // end of argument printing:


  // Set random seed
  rng.seed(mParameters.Seed());


  // Basic checks of the input:
  double last_start_time = 0.0;
  for (int i = 0; i < mParameters.NumberOfClones(); i++) {
    // Ensure that the clone start times are ordered:
    if (mParameters.CloneStartTimes()[i] < last_start_time) {
      Rcpp::Rcerr << "Defined clone start times have to be sorted!" << std::endl;
      return 0;
    }
    last_start_time = mParameters.CloneStartTimes()[i];

    // Ensure that fathers are defined:
    if ( (i == 0 && mParameters.Fathers()[i] != 0) ||
         (i > 0 && mParameters.Fathers()[i] >= i)
       )
    {
      Rcpp::Rcerr << "Fathers of first clone has to be 0 and" << std::endl;
      Rcpp::Rcerr << "Fathers of clones have to be introduced first!";
      Rcpp::Rcerr << std::endl;
      return 0;
    }
  }


  // Create reusable pointers to handled cells, types
  // and other variables used for simulations:
  Cell* pCell;
  CellType* pType;
  long int generation = 0;
  long double dt = 0.0L;
  int action = -1;


  // Create all cell types and put them into a vector:
  std::vector <CellType*> pvCelltypes; // vector of all types
  for (int i = 0; i < mParameters.NumberOfClones(); i++) {

    // New cell type:
    pType = new CellType(mParameters.Birthrates()[i],
                         mParameters.Deathrates()[i],
                         mParameters.Mutationrates()[i]);

    // Insert cell types:
    pvCelltypes.push_back(pType);
  }


  // Create and insert a cell of the first type:
  pCell = new Cell(pvCelltypes[0]);
  this->InsertCell(pCell);
  pCell->MutateCellFixedNumber(mParameters.NumberClonalMutations());


  // universe time to start time of fist clone:
  int next_clone = 1; // count number of next type to introduce:
  // Debug messages:
  if (verbose) {
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "########## Comment #############" << std::endl;
    Rcpp::Rcout << " New clone: " << next_clone - 1 << std::endl;
    Rcpp::Rcout << " Time: " << generation << std::endl;
    Rcpp::Rcout << "################################" << std::endl;
    Rcpp::Rcout << std::endl;
  }


  // Run till limit of universe has been reached:
  while (++generation <= mParameters.SimulationEndTime()) {

    // Check if new clones need to be introduced in each cycle:
    for (; next_clone < mParameters.NumberOfClones() &&
           mParameters.CloneStartTimes()[next_clone] <= generation;
           next_clone++)
    {
      // Introduce a new clone:
      int c_father = mParameters.Fathers()[next_clone];
      pvCelltypes[c_father]->RandomMember()->Type(pvCelltypes[next_clone]);

      // Debug messages:
      if (verbose) {
        Rcpp::Rcout << std::endl;
        Rcpp::Rcout << "########## Comment #############" << std::endl;
        Rcpp::Rcout << " New clone: " << next_clone << std::endl;
        Rcpp::Rcout << " Time: " << generation << std::endl;
        Rcpp::Rcout << "################################" << std::endl;
        Rcpp::Rcout << std::endl;
      }
    } // end of for loop for introduction of new clones

    // Select reaction type and members:
    this->NextReaction(&dt, &action)->RandomMember()->DoAction(&action);
    this->IncrementTimeBy(dt);
    mNumberOfReactions++;
  } // stop running after reaching limit

  return true;
}

bool Universe::RunSimulation() {
  return RunSimulation(false);
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

bool Universe::Sample(std::vector <int> &clone , std::vector <double>& vaf) const {

  // Determine total number of cells:
  int ncells = 0;
  for (std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++){
    ncells += mpTypes[i]->NumMembers();
  }

  // Traverse the Phylogenies:
  for (std::vector<CellType*>::size_type i = 0; i < mpPhylogenies.size(); i++){
    mpPhylogenies[i]->Root()->SampleNode(clone, vaf,
                                         mParameters.SequencingMinVaf(),
                                         mParameters.TumourPurity(),
                                         mParameters.SequencingDepth(),
                                         mParameters.SequencingDepthModel(),
                                         ncells);
  }

  return true;
}



bool Universe::Sample(std::vector <int> &clone,
                      std::vector <int> &alt,
                      std::vector <int> &depth,
                      std::vector <std::string> &ids)
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
                                                    mParameters.SequencingMinVaf(),
                                                    mParameters.TumourPurity(),
                                                    mParameters.SequencingDepth(),
                                                    mParameters.SequencingDepthModel(),
                                                    ncells);
  }

  if (ncells != ncells2) {
    Rcpp::Rcerr << "sum_cells != total_cells: " << ncells << " vs " << ncells2 << std::endl;
  }

  return true;
}


bool Universe::SampleToFile(std::string out_file) const {

  std::ofstream ostream(out_file);
  if (ostream.is_open()) {
    ostream << "VAF\tALT\tDP\tCLONE\tTRUE_VAF\tTRUE_CLUSTER" << std::endl;

    // Determine total number of cells:
    int ncells = 0;
    for (std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++){
      ncells += mpTypes[i]->NumMembers();
    }

    // Traverse the Phylogenies:
    for (std::vector<CellType*>::size_type i = 0; i <mpPhylogenies.size(); i++){
      mpPhylogenies[i]->Root()->SampleNode(ostream,
                                           mParameters.SequencingMinVaf(),
                                           mParameters.TumourPurity(),
                                           mParameters.SequencingDepth(),
                                           mParameters.SequencingDepthModel(),
                                           ncells);
    }

  } else {
    Rcpp::Rcout << "Error." << std::endl;
    Rcpp::Rcout << "  Unable to open sequencing output file:" << std::endl;
    Rcpp::Rcout << "    " << out_file << std::endl;
    return false;
  }
  ostream.close();
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


// Output Functions:
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

bool Universe::WriteCellCountsToFile(std::string out_file) const {
  unsigned long sum_cells = 0;
  char delim = '\t';

  std::ofstream outstream(out_file);
  if (outstream.is_open()) {

    // Write a header:
    outstream << "total";
    for(std::vector<CellType *>::size_type i = 0; i < mpTypes.size(); i++) {
      outstream << delim << "clone" << i;
      sum_cells += mpTypes[i]->NumMembers();
    }
    outstream << std::endl;

    // Write the cell counts to the next line:
    outstream << sum_cells;
    for(std::vector<CellType *>::size_type i = 0; i < mpTypes.size(); i++) {
      outstream << delim << mpTypes[i]->NumMembers();
    }
    outstream << std::endl;

    // Write the expected peak position assuming CN == 2 to the next line:
    outstream << sum_cells * 0.5 / sum_cells;
    for(std::vector<CellType *>::size_type i = 0; i < mpTypes.size(); i++) {
      outstream << delim << std::setprecision(5)
                << mpTypes[i]->NumMembers() * 0.5 / sum_cells;
    }
    outstream << std::endl;

  } else {

    Rcpp::Rcout << "Error." << std::endl;
    Rcpp::Rcout << "  Unable to open cell count output file:" << std::endl;
    Rcpp::Rcout << "    " << out_file << std::endl;
    return false;

  }
  outstream.close();
  return true;
}
