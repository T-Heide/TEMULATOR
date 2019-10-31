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


#include "SimulationParameterSet.h"
#include <iostream>
#include <fstream>
#include "Rcpp.h"


// Constructor:
SimulationParameterSet::SimulationParameterSet(
    std::vector <double> mutationrates,
    std::vector <double> birthrates,
    std::vector <double> deathrates,
    std::vector <unsigned int> clone_start_times,
    std::vector <int> fathers,
    unsigned int simulation_end_time,
    int seed,
    int number_clonal_mutations,
    double tumour_purity,
    double sequencing_depth,
    double minimum_vaf,
    int depth_model,
    std::string output_prefix
  ) :
    mMutationrates(mutationrates),
    mBirthrates(birthrates),
    mDeathrates(deathrates),
    mCloneStartTimes(clone_start_times),
    mFathers(fathers),
    mSimulationEndTime(simulation_end_time),
    mSeed(seed),
    mNumberClonalMutations(number_clonal_mutations),
    mTumourPurity(tumour_purity),
    mSequencingDepth(sequencing_depth),
    mMinSequencingVaf(minimum_vaf),
    mDepthModel(depth_model),
    mOutputPrefix(output_prefix)
  {
    if (  mutationrates.size() != birthrates.size() ||
          mutationrates.size() != deathrates.size() ||
          mutationrates.size() != clone_start_times.size() ||
          mutationrates.size() != fathers.size()
        )
      {
        Rcpp::Rcerr << "Error during creation of SimulationParameterSet:\n";
        Rcpp::Rcerr << "  Parameter vectors not of same length.\n";
        Rcpp::Rcerr << "    Mutation rates: " << mutationrates.size() << "\n";
        Rcpp::Rcerr << "    Birth rates: " << birthrates.size() << "\n";
        Rcpp::Rcerr << "    Death rates: " << deathrates.size() << "\n";
        Rcpp::Rcerr << "    Start times: " << clone_start_times.size() << "\n";
        Rcpp::Rcerr << "    Father types: " << fathers.size() << "\n";
        Rcpp::Rcerr << std::endl;
        Rcpp::stop("Parameter vectors not of same length.");
      }
      mNumberOfClones = static_cast <int> (mutationrates.size());
  }


// Getters:
int SimulationParameterSet::NumberOfClones() const {
  return mNumberOfClones;
}

std::vector <double> SimulationParameterSet::Mutationrates() const {
  return mMutationrates;
}

std::vector <double> SimulationParameterSet::Birthrates() const {
  return mBirthrates;
}

std::vector <double> SimulationParameterSet::Deathrates() const {
  return mDeathrates;
}

std::vector <unsigned int> SimulationParameterSet::CloneStartTimes() const {
  return mCloneStartTimes;
}

std::vector <int> SimulationParameterSet::Fathers() const {
  return mFathers;
}

unsigned int SimulationParameterSet::SimulationEndTime() const {
  return mSimulationEndTime;
}

int SimulationParameterSet::Seed() const {
  return mSeed;
}

double SimulationParameterSet::TumourPurity() const {
  return mTumourPurity;
}

int SimulationParameterSet::NumberClonalMutations() const {
  return mNumberClonalMutations;
}

double SimulationParameterSet::SequencingDepth() const {
  return mSequencingDepth;
}

double SimulationParameterSet::SequencingMinVaf() const {
  return mMinSequencingVaf;
}

int SimulationParameterSet::SequencingDepthModel() const {
  return mDepthModel;
}

std::string SimulationParameterSet::OutputPrefix() const {
  return mOutputPrefix;
}


// Output functions:
void SimulationParameterSet::Print() const {
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
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "  Clonal mutations: " << mNumberClonalMutations << "\n";
  Rcpp::Rcout << "  Tumour cellularity: " << mTumourPurity << "\n";
  Rcpp::Rcout << "  VAF cutoff: " << mMinSequencingVaf << "\n";
  Rcpp::Rcout << "  Sequencing depth: " << mSequencingDepth << "\n";
  Rcpp::Rcout << "  Sequencing depth model: " << mDepthModel << "\n";
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "  Random seed: " << mSeed << std::endl;
  Rcpp::Rcout << "  Output prefix: " << mOutputPrefix << std::endl;
  Rcpp::Rcout << "################################" << std::endl;
  Rcpp::Rcout << std::endl;
}

bool SimulationParameterSet::WriteParamsToFile(std::string param_output_file) const {

  std::ofstream output_stream(param_output_file);

  if (output_stream.is_open()) {
    // Then other singlular params:
    output_stream << "# Seed: " << mSeed << std::endl;
    output_stream << "# N clones: " << mNumberOfClones << std::endl;
    output_stream << "# N clonal: " << mNumberClonalMutations << std::endl;
    output_stream << "# Tumour cellularity: " << mTumourPurity << std::endl;
    output_stream << "# Simulation end: " << mSimulationEndTime << std::endl;
    output_stream << "# Sequencing depth: " << mSequencingDepth << std::endl;
    output_stream << "# Sequencing depth model: " << mDepthModel << std::endl;
    output_stream << "# LOD sequencing: " << mMinSequencingVaf << std::endl;
    output_stream << "# Output prefix: " << mOutputPrefix << std::endl;

    // Now the header:
    char delim = '\t'; // TSV seperated.
    output_stream << "type_number" << delim
                  <<  "mutation_rate" << delim
                  <<  "birth_rate" << delim
                  <<  "death_rates" << delim
                  <<  "start_time" << delim
                  <<  "father" << std::endl;

    // And now clone types line by line:
    for (int i = 0; i < mNumberOfClones; i++) {
      output_stream << i << delim // Column by column ...
                    << mMutationrates[i] << delim
                    << mBirthrates[i] << delim
                    << mDeathrates[i] << delim
                    << mCloneStartTimes[i] << delim
                    << mFathers[i] << std::endl;
    }
  } else { // failed to open out_file
    Rcpp::Rcout << "Error." << std::endl;
    Rcpp::Rcout << "  Unable to open sim. parameter output file:" << std::endl;
    Rcpp::Rcout << "    " << param_output_file << std::endl;
    return false;
  }

  output_stream.close();
  return true;
}


bool SimulationParameterSet::WriteAdjacencyMatrixToFile(std::string out_file) const {
  char delim = '\t';

  std::ofstream outstream(out_file);
  if (outstream.is_open()) {
    for(std::vector<int>::size_type i = 0; i < mFathers.size(); i++) {
      for(std::vector<int>::size_type j = 0; j < mFathers.size(); j++) {
        if (mFathers[j] == i && j != i) {
          outstream << 1 << delim;
        } else {
          outstream << 0 << delim;
        }
      }
      outstream << std::endl;
    }
  } else {
    Rcpp::Rcout << "Error." << std::endl;
    Rcpp::Rcout << "  Unable to open adjcaceny matrix output file:" << std::endl;
    Rcpp::Rcout << "    " << out_file << std::endl;
    return false;
  }

  outstream.close();
  return true;
}
