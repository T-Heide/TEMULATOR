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
#include "SimulationParameterSet.h"

// Universe ////////////////////////////////////////////////////////////////////

class Universe {
    long double mTime;
    unsigned int mNumberOfCells;
    unsigned int mNumberOfReactions;
    std::vector<CellType*> mpTypes;
    std::vector<PhylogenyRoot*> mpPhylogenies;
    SimulationParameterSet mParameters;

  public:
    // Constructor:
    Universe(SimulationParameterSet);

    // Destructor:
    ~Universe();

    // Simulate tumour:
    bool RunSimulation(bool);
    bool RunSimulation();

    // Getter functions:
    double Time() const;
    unsigned int Reactions() const;
    
    class CellType* NextReaction(long double*, int*) const;
    bool Sample(std::vector <int>&, std::vector <double>&) const;
    bool Sample(std::vector <int>&, std::vector <int>&, std::vector <int> &, std::vector <std::string>&);
    bool SampleToFile(std::string) const;
    std::vector <unsigned int> CellCounts() const;

    // Setter functions:
    void IncrementTimeBy(long double);
    bool InsertCell(Cell*);
    bool InsertCell(Cell*, bool);
    void RegisterType(CellType *);

    // Output Functions:
    void PrintCellTypes() const;
    bool WriteCellCountsToFile(std::string) const;

    //void TypesToCsvFile(std::string, int tum_id);
    //void PutNodesToStream(PhylogenyNode*, std::ofstream&);
    //void PhylogeniesToFile(std::string, int tum_id);

};

#endif // UNIVERSE_H
