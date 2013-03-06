/* -*- c++ -*- */
#ifndef BANKLENNARDJONESPARAMETERTABLE_H
#define BANKLENNARDJONESPARAMETERTABLE_H

#include <protomol/topology/LennardJonesParameterTable.h>
#include <vector>

namespace ProtoMol {
  //________________________________________ BankLennardJonesParameterTable
  /**
   * Keeps a bank of different LennardJonesParameterTable's, iSG
   */

  struct BankLennardJonesParameterTable {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    BankLennardJonesParameterTable() {};

    /// function which creates the correct # of parameter tables
    void createTables(unsigned int numComp, unsigned int tableSize) {
      // compute the # of needed tables
      numTables = numComp * numComp;
      numComps = numComp;

      // allocate memory for the parameter tables
      lennardJonesBank.resize(numTables);

      // properly size each individual parameter table
      for (unsigned int i = 0; i < numTables; i++)
        lennardJonesBank[i].resize(tableSize);
    }

    /// determine the proper index # into to the bank
    const LennardJonesParameters &operator()(int i, int j, int type1,
                                             int type2) const {
      // **TIM** look over notes and add a comment explaining why bankIndex
      // is always = i*NumComps+j

      // get the LJ parameters for this pair
      return lennardJonesBank[i * numComps + j](type1, type2);
    }

    /// function that stores a set of LJ parameters for one pair of atoms
    void set(int i,
             int j,
             int type1,
             int type2,
             const LennardJonesParameters &paramsij) {
      lennardJonesBank[i * numComps + j].set(type1, type2, paramsij);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    /// The bank of necessary LJ parameter tables
    std::vector<LennardJonesParameterTable> lennardJonesBank;
    unsigned int numTables;
    unsigned int numComps;
  };
  //________________________________________ INLINES
}
#endif /* BANKLENNARDJONESPARAMETERS_H */
