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


#define DEPTH_OVERDISPERSION 0.08


#include "extern_global_variables.h"
#include "Phylogeny.h"
#include "Cell.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/beta_distribution.hpp>
#include "Rcpp.h"


// PhylogenyRoot ///////////////////////////////////////////////////////////////

// Statics:
unsigned int PhylogenyRoot::msNextId = 0;

// Constructor:
PhylogenyRoot::PhylogenyRoot(Cell* pCell) : mId(msNextId++) {
  mpRoot = new PhylogenyNode(pCell);
}

PhylogenyNode* PhylogenyRoot::Root(){ return mpRoot; }

//Destructor
PhylogenyRoot::~PhylogenyRoot() {
  delete mpRoot;
}

// PhylogenyNode ///////////////////////////////////////////////////////////////

// Statics:
unsigned long PhylogenyNode::msNextId = 0;

// Constructors
PhylogenyNode::PhylogenyNode(Cell* pCell) //
  : mId(msNextId++),
    mpCell(pCell),
    mpUp(0),
    mpLeft(0),
    mpRight(0),
    mGeneration(1),
    mNumMutsGeneration(0)
  {
    if (pCell->AssociatedNode() != 0) {
      pCell->AssociatedNode()->AssociatedCell(0);
    }
    pCell->AssociatedNode(this);

    mTypeId = pCell->Type() == 0 ? 0 : pCell->Type()->Id();
  }

PhylogenyNode::PhylogenyNode(Cell* pCell, PhylogenyNode* pUp)
  : mId(msNextId++),
  mpCell(pCell),
  mpUp(pUp),
  mpLeft(0),
  mpRight(0),
  mNumMutsGeneration(0)
  {
    mGeneration = pUp->Generation() + 1;
    if (pCell->AssociatedNode() != 0) {
      pCell->AssociatedNode()->AssociatedCell(0);
    }
    pCell->AssociatedNode(this);

    mTypeId = pCell->Type() == 0 ? 0 : pCell->Type()->Id();
  }

//Destructor
PhylogenyNode::~PhylogenyNode() {

  // Delete left node:
  if (mpLeft != 0) {
    delete mpLeft;
    mpLeft = 0;
  }

  // Delete right node:
  if (mpRight != 0) {
    delete mpRight;
    mpRight = 0;
  }

  // Delete associated cell:
  if (mpCell != 0) {
    mpCell->AssociatedNode(0); // unlink cell first
    delete mpCell;
    mpCell = 0;
  }

  // Clean pointers in the up node:
  if (mpUp != 0) {
    if (mpUp->LeftNode() == this) {
      mpUp->LeftNode(0);
    } else if (mpUp->RightNode() == this) {
      mpUp->RightNode(0);
    }
  }
}



// Getters:
unsigned long PhylogenyNode::Id() const {return mId;}
unsigned int PhylogenyNode::Generation() const {return mGeneration;}
unsigned int PhylogenyNode::NumMutations() const {return mNumMutsGeneration;}
Cell* PhylogenyNode::AssociatedCell() const {return mpCell;}
PhylogenyNode* PhylogenyNode::UpNode() const {return mpUp;}
PhylogenyNode* PhylogenyNode::LeftNode() const {return mpLeft;}
PhylogenyNode* PhylogenyNode::RightNode() const {return mpRight;}
int PhylogenyNode::TypeId() const {return mTypeId;};


void simulateNGS(double vaf, int dp, double& sim_vaf, int& sim_alt, int& sim_dp,
                 int depth_model)
{
  switch (depth_model) {

    case 1: {// poisson distributed depth model
      boost::random::poisson_distribution<int> dist_depth1(dp);
      sim_dp = dist_depth1(rng);
      break;
    }

    case 2: {// overdispersed beta binomial (default).
      double mu = 0.6;
      double rho = DEPTH_OVERDISPERSION;
      double sh1 = mu * (1.0 / rho - 1.0); // scale param 1 (beta)
      double sh2 = (mu - 1.0) * (rho - 1.0) / rho; // scale param 2 (beta)

      // beta component
      boost::random::beta_distribution <double> beta_component_depth2(sh1, sh2);
      double cbd = beta_component_depth2(rng);

      // binomial component
      boost::random::binomial_distribution<int> dist_depth2(dp / mu, cbd);
      sim_dp = dist_depth2(rng);
      break;
    }

    case 3: { // fixed depth
      sim_dp = dp;
      break;
    }

    default: {
      Rcpp::stop("Error. Undefined model type '"+std::to_string(depth_model)+"!\n");
    }
  }

  // binomial sampling
  boost::random::binomial_distribution<int> binom_alt(sim_dp, vaf);
  sim_alt = binom_alt(rng);
  sim_vaf = sim_alt * 1.0 / sim_dp;

  return;
}


int PhylogenyNode::SampleNode(std::vector <int>& clone,
                              std::vector <int>& alt,
                              std::vector <int>& depth,
                              std::vector <std::string>& ids,
                              double minVAF,
                              double purity,
                              int dp,
                              int dp_model,
                              int total_cells)
{
  int n_cells = 0;

  if (mpLeft != 0) {
    n_cells += mpLeft->SampleNode(clone, alt, depth, ids,
                                  minVAF, purity, dp, dp_model,
                                  total_cells);
  }

  if (mpRight != 0) {
    n_cells += mpRight->SampleNode(clone, alt, depth, ids,
                                  minVAF, purity, dp, dp_model,
                                  total_cells);
  }

  if (mpCell != 0) {
    n_cells += 1;
  }

  if (n_cells != 0) {
    double exp_vaf = 0.5 * n_cells / total_cells * purity;
    double sim_vaf = 0.0;
    int sim_alt = 0;
    int sim_dp = 0;

    // Convert node hex:
    std::stringstream out_stream;
    out_stream << std::hex << mId;
    std::string hex_node_id = out_stream.str();

    for (unsigned int i = 0; i < mNumMutsGeneration; i++) {
      // Convert mutation id to a two compound hex:
      out_stream.str("");
      out_stream << std::hex << i;
      std::string mutation_id = hex_node_id + "X" + out_stream.str();

      simulateNGS(exp_vaf, dp, sim_vaf, sim_alt, sim_dp, dp_model);
      if (sim_vaf >= minVAF) {
        clone.push_back(this->TypeId());
        alt.push_back(sim_alt);
        depth.push_back(sim_dp);
        ids.push_back(mutation_id);
      }
    }
  }

  return n_cells;
}

int PhylogenyNode::SampleNode(std::vector <int>&clone,
                              std::vector <double>& vaf,
                              double minVAF, double purity,
                              int dp, int dp_model,
                              int total_cells) const
{
  int n_cells = 0;

  if (mpLeft != 0) {
    n_cells += mpLeft->SampleNode(clone, vaf, minVAF, purity, dp, dp_model,
                                  total_cells);
  }

  if (mpRight != 0) {
    n_cells += mpRight->SampleNode(clone, vaf, minVAF, purity, dp, dp_model,
                                   total_cells);
  }

  if (mpCell != 0) {
    n_cells += 1;
  }

  if (n_cells != 0) {
    double exp_vaf = 0.5 * n_cells / total_cells * purity;
    double sim_vaf = 0.0;
    int sim_alt = 0;
    int sim_dp = 0;

    for (unsigned int i = 0; i < mNumMutsGeneration; i++) {
      simulateNGS(exp_vaf, dp, sim_vaf, sim_alt, sim_dp, dp_model);
      if (sim_vaf >= minVAF) {
        vaf.push_back(sim_vaf);
        clone.push_back(this->TypeId());
      }
    }
  }

  return n_cells;
}



int PhylogenyNode::SampleNode(std::ofstream& ostream,
                              double minVAF, double purity,
                              int dp, int dp_model,
                              int total_cells) const
{
  int n_cells = 0;

  if (mpLeft != 0) {
    n_cells += mpLeft->SampleNode(ostream, minVAF, purity, dp, dp_model,
                                  total_cells);
  }

  if (mpRight != 0) {
    n_cells += mpRight->SampleNode(ostream, minVAF, purity, dp, dp_model,
                                   total_cells);
  }

  if (mpCell != 0) {
    n_cells += 1;
  }

  if (n_cells != 0) {
    double exp_vaf = 0.5 * n_cells / total_cells * purity;
    double sim_vaf = 0.0;
    int sim_alt = 0;
    int sim_dp = 0;

    for (unsigned int i = 0; i < mNumMutsGeneration; i++) {
      simulateNGS(exp_vaf, dp, sim_vaf, sim_alt, sim_dp, dp_model);
      if (sim_vaf >= minVAF) {
        std::stringstream ss;
        ss << std::scientific << std::setprecision(3) << sim_vaf << "\t";  // simulated VAF
        ss << sim_alt << "\t" << sim_dp << "\t" << this->TypeId() << "\t"; // simulated depth
        ss << std::scientific << std::setprecision(5) << exp_vaf << "\t";  // expected VAF
        ss << std::hex << mId;                                             // node id/true cluster
        ostream << ss.str() << std::endl;
      }
    }
  }

  return n_cells;
}


// Setters:
void PhylogenyNode::AddNewMutations(unsigned int num_new_mutations) {
  mNumMutsGeneration += num_new_mutations;
}

void PhylogenyNode::AssociatedCell(Cell* pCell) { mpCell = pCell; }

void PhylogenyNode::UpNode(PhylogenyNode* pNode) {mpUp = pNode;};
void PhylogenyNode::LeftNode(PhylogenyNode* pNode) {mpLeft = pNode;};
void PhylogenyNode::RightNode(PhylogenyNode* pNode) {mpRight = pNode;};
void PhylogenyNode::TypeId(int newTypeId) {mTypeId = newTypeId;};



// Output:
void PhylogenyNode::Print() const {
  unsigned int up = (mpUp == 0) ? 0 : mpUp->Id();
  unsigned int left = (mpLeft == 0) ? 0 : mpLeft->Id();
  unsigned int right = (mpRight == 0) ? 0 : mpRight->Id();

  Rcpp::Rcout << "###### Phylogenetic Node ######" << std::endl;
  Rcpp::Rcout << "   ID: " << mId << std::endl;
  Rcpp::Rcout << "   Up: " << up << std::endl;
  Rcpp::Rcout << "   Left: " << left << std::endl;
  Rcpp::Rcout << "   Right: " << right << std::endl;
  Rcpp::Rcout << "   Generation: " << mGeneration << std::endl;
  Rcpp::Rcout << "   Mutations: " << mNumMutsGeneration << std::endl;
  if (mpCell != 0) { Rcpp::Rcout << "   Cell: " << mpCell << std::endl; }
  Rcpp::Rcout << "###############################" << std::endl;
}

void PhylogenyNode::PrintAncestry() const {

  if (mpLeft == 0 && mpRight == 0) { // If this is a leaf
    std::string cell_id = (mpCell == 0) ? "NA" : std::to_string(mpCell->Id());
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "###### Ancestry of cell #######" << std::endl;
    Rcpp::Rcout << "  Cell ID: " << cell_id << std::endl;
    Rcpp::Rcout << std::endl;
  }

  this->Print();

  if (mpUp != 0) { // If next is not root
    Rcpp::Rcout << "              | " << std::endl;
    Rcpp::Rcout << "              | " << std::endl;
    mpUp->PrintAncestry();
  } else {
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "###############################" << std::endl;
  }
}
