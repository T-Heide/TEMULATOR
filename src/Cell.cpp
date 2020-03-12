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
#include "Cell.h"
#include "CellType.h"
#include "Phylogeny.h"
#include "Universe.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <iostream>
#include <array>
#include <vector>
#include <string>
#include "Rcpp.h"



// Statics:
unsigned long Cell::msNextId = 1;


// Constructors:
Cell::Cell ()
  : mId(msNextId++),
    mpUniverse(0),
    mpType(0),
    mpNode(0)
  {}

Cell::Cell (CellType* pType)
  : mId(msNextId++),
    mpUniverse(0),
    mpType(pType),
    mpNode(0)
  {
    mpType->RegisterMember(this);
  }


// Destructor:
Cell::~Cell(){ // Proper deletion of a cell.

  mpType->DeregisterMember(this);      // Let the associated type forget.

  if (mpNode != 0) {
    mpNode->AssociatedCell(0); // Unlink cell from the associated node.

    // Trim back the tree:
    if (mpNode->LeftNode() == 0 && mpNode->RightNode() == 0) { // trim
      PhylogenyNode* cNode = mpNode->UpNode(); // current node.
      PhylogenyNode* lNode = mpNode;           // last node.

      while (cNode != 0 && // didn't reach root of tree.
             (cNode->LeftNode() == 0 || cNode->RightNode() == 0) && // not branching at this level
             cNode->AssociatedCell() == 0) // node does not hold any cell
      {
        lNode = cNode;
        cNode = cNode->UpNode();
      }

      if (cNode != 0) { // Don't delete the root node ...
        delete lNode;
      }
    }
    mpNode = 0; // mpNode is gonne/should not be deleted ...
  }
  // Removal of cell complete.
}


// Getter functions:
unsigned long Cell::Id() const {return(mId);};
CellType* Cell::Type() const {return mpType;};
PhylogenyNode* Cell::AssociatedNode() const {return mpNode;}; // Assoc. node
Universe* Cell::AssociatedUniverse() const {return mpUniverse;}; // Assoc. node
int Cell::TypeIndex() const {return mTypeIndex;};

bool Cell::AsProgenitorDies() const { // Random realization to die with prob alpha.
  if (mpType->NumMembers() == 1) {
    return false;
  } else {
    double alpha = mpType->Alpha();
    boost::random::bernoulli_distribution<> bern_alpha(alpha);
    bool res = bern_alpha(rng);
    return res;
  }
}



// Setter functions:
void Cell::Type(CellType* pNewType) {
  mpUniverse->RegisterType(pNewType);
  mpType->DeregisterMember(this);
  pNewType->RegisterMember(this);
  mpType = pNewType;

  if (mpNode != 0) { // If the cell has a associated node, member of a universe,
    mpNode->TypeId(pNewType->Id());
  }

  // Update the TypeId of all ancestral nodes to 0 (undertermined/mix type):
  PhylogenyNode *pCurrentNode = mpNode;
  while ((pCurrentNode = pCurrentNode->UpNode()) != 0) {
    pCurrentNode->TypeId(0);
  }
};

void Cell::AssociatedNode(PhylogenyNode* new_node) { mpNode = new_node; };
void Cell::AssociatedUniverse(Universe* new_universe) {
  mpUniverse = new_universe;
}

void Cell::TypeIndex(int newIndex) {mTypeIndex = newIndex;};


// Other functions:
void Cell::MutateCell() {
  Cell::MutateCell(mpType->Mu());
}

void Cell::MutateCell(double mu){
  // Sample number of new muts from poisson:
  boost::random::poisson_distribution<int> dist_mutations(mu);
  int new_muts = dist_mutations(rng);

  // Append muts to node:
  mpNode->AddNewMutations(new_muts);
}

void Cell::MutateCellFixedNumber(int new_muts){
  // Append muts to node:
  mpNode->AddNewMutations(new_muts);
}



void Cell::Print() const {
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "############ Cell ##############" << std::endl;
  Rcpp::Rcout << "   ID: " << mId << std::endl;
  Rcpp::Rcout << "   Location:" << std::endl;
  Rcpp::Rcout << "       Universe: " << mpUniverse << std::endl;
  Rcpp::Rcout << "   Type: " << mpType->Id() << std::endl;
  Rcpp::Rcout << "       Birth rate: " << mpType->Birthrate() << std::endl;
  Rcpp::Rcout << "###############################" << std::endl;
}

void Cell::DoAction(int *action) {
  switch(*action) {
    case 1: // Action 'Divide'
      this->Divide();
      break;
    default:
      Rcpp::stop("Unknow action.\n");
      break;
  }
}

void Cell::Divide() {

  // Cells that were not introduced into a universe can't divide!
  if (mpUniverse == 0) {
    Rcpp::Rcerr << "Cells not introduced into universe can't divide!" << std::endl;
    return;
  }

  if(this->AsProgenitorDies()){
    delete this;
  } else {

    // Store current associated node, then branch to the left and mutate:
    PhylogenyNode* old_node = this->mpNode;
    PhylogenyNode* new_left_node = new PhylogenyNode(this, old_node);
    old_node->LeftNode(new_left_node);
    this->MutateCell();

    // Create daughter, insert on branch to the right of old node and mutate:
    Cell* pDaughter = new Cell(mpType);
    PhylogenyNode* new_right_node = new PhylogenyNode(pDaughter, old_node);
    old_node->RightNode(new_right_node);
    mpUniverse->InsertCell(pDaughter, false);
    pDaughter->MutateCell();
  }
}
