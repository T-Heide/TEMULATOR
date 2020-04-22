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


#ifndef UNIVERSE_H
#define UNIVERSE_H

// Forward declerations: ///////////////////////////////////////////////////////
class Cell;
class CellType;
class PhylogenyRoot;
class Shape;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>
#include "Rcpp.h"
#include <boost/random.hpp>

// Universe ////////////////////////////////////////////////////////////////////

class Universe {
  
   // summaries
    unsigned int mNumberOfClones;
    unsigned int mNumberOfCells;
    unsigned int mNumberOfReactions;
    long double mTime;

    // other parameters
    unsigned int mSimulationEndTime;
    unsigned int mNextClone;
    int mSeed;
    
    // clone parameters
    std::vector <double> mMutationrates;
    std::vector <double> mBirthrates;
    std::vector <double> mDeathrates;
    std::vector <unsigned int> mCloneStartTimes;
    std::vector <unsigned int> mFathers;
    unsigned int mClonalMutations;
    
    // other objects
    std::vector<CellType*> mpTypes;
    std::vector<PhylogenyRoot*> mpPhylogenies;
    boost::random::mt19937_64 mRngState;
    
  public:
    // Constructor:
    Universe(
      Rcpp::DataFrame,
      unsigned int,
      unsigned int,
      int);
    
    // Destructor:
    ~Universe();

    // Simulate tumour:
    bool RunSimulation(bool);

    // Getter functions:
    double Time() const;
    unsigned int Reactions() const;
    std::vector <unsigned int> CellCounts() const;
    
    class CellType* NextReaction(long double*, int*) const;
    bool Sample(std::vector <int>&, std::vector <int>&, std::vector <int> &, std::vector <std::string>&, double, double, double, int) const;
    Rcpp::DataFrame SampleRcpp(double, double, double, int) const;
    Rcpp::DataFrame SampleSeededRcpp(double, double, double, int, int) const;
    Rcpp::DataFrame Sample(double, double, double, int) const;
    unsigned int get_SimulationEndTime() const;
      
    
    // Setter functions:
    void IncrementTimeBy(long double);
    bool InsertCell(Cell*);
    bool InsertCell(Cell*, bool);
    void RegisterType(CellType *);
    void set_SimulationEndTime(unsigned int);
    

    // Output Functions:
    void Print() const;
    void PrintCellTypes() const;

};

#endif // UNIVERSE_H
