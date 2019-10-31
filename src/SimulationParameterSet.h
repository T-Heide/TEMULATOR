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


#ifndef SIMULATIONARAMETERSET_H
#define SIMULATIONARAMETERSET_H

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>
#include <string>

class SimulationParameterSet {
  public:

    // Cell type related:
    int mNumberOfClones;
    std::vector <double> mMutationrates;
    std::vector <double> mBirthrates;
    std::vector <double> mDeathrates;
    std::vector <unsigned int> mCloneStartTimes;
    std::vector <int> mFathers;
    unsigned int mSimulationEndTime;
    int mSeed;

    // Sequencing related:
    int mNumberClonalMutations;
    double mTumourPurity;
    double mSequencingDepth;
    double mMinSequencingVaf;
    int mDepthModel;

    // Output related:
    std::string mOutputPrefix;

  public:
    // Constructor:
    SimulationParameterSet(std::vector <double>,       // mutation rates
                           std::vector <double>,       // birth rates
                           std::vector <double>,       // death rates
                           std::vector <unsigned int>, // clone start times
                           std::vector <int>,          // fathers
                           unsigned int,               // simulation end time
                           int,                        // seed
                           int,                        // number of clonal mutations
                           double,                     // tumour purity
                           double,                     // sequencing depth
                           double,                     // minimum vaf
                           int,                        // depth model
                           std::string);               // output prefix


    // Getters:
    int NumberOfClones() const;
    std::vector <double> Mutationrates() const;
    std::vector <double> Birthrates() const;
    std::vector <double> Deathrates() const;
    std::vector <unsigned int> CloneStartTimes() const;
    std::vector <int> Fathers() const;
    unsigned int SimulationEndTime() const;
    int Seed() const;

    int NumberClonalMutations() const;
    double TumourPurity() const;
    double SequencingDepth() const;
    double SequencingMinVaf() const;
    int SequencingDepthModel() const;

    std::string OutputPrefix() const;


    // Output functions:
    void Print() const;
    bool WriteParamsToFile(std::string) const;
    bool WriteAdjacencyMatrixToFile(std::string) const;
};

#endif // SIMULATIONARAMETERSET_H
