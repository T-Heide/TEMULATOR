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

#ifndef CELL_H
#define CELL_H

// Forward declerations: ///////////////////////////////////////////////////////
class CellType;
class Universe;
class PhylogenyNode;

//class Phylogeny_Node;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>
#include <array>
#include "CellType.h"

// Cell ////////////////////////////////////////////////////////////////////////

class Cell {
    // Identifiers:
    const unsigned long mId;
    static unsigned long msNextId;

    // Location related:
    Universe* mpUniverse;
    int mTypeIndex;

    // Properties:
    CellType* mpType;
    PhylogenyNode* mpNode;

  public:
    // Constructors:
    Cell ();
    Cell (CellType*);

    // Destructor:
    ~Cell ();

    // Getter functions:
    unsigned long Id() const;
    CellType* Type() const;
    bool AsProgenitorDies() const;
    PhylogenyNode* AssociatedNode() const;
    Universe* AssociatedUniverse() const;
    int TypeIndex() const;

    // Setter functions:
    void Type(CellType*);
    void AssociatedNode(PhylogenyNode*);
    void AssociatedUniverse(Universe*);
    void TypeIndex(int);

    // Other functions:
    void MutateCell();
    void MutateCell(double);
    void MutateCellFixedNumber(int);

    void DoAction(int*);
    void Divide();
};

#endif // CELL_H
