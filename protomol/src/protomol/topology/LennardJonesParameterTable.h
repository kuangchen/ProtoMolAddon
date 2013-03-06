/* -*- c++ -*- */
#ifndef LENNARDJONESPARAMETERTABLE_H
#define LENNARDJONESPARAMETERTABLE_H

#include <protomol/type/Real.h>
#include <protomol/topology/LennardJonesParameters.h>

namespace ProtoMol {
  //________________________________________ LennardJonesParameterTable

  /**
   * This table contains the Lennard-Jones parameters for all atom type pairs.
   * It is guaranteed to be symmetric.
   */
  class LennardJonesParameterTable {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    LennardJonesParameterTable();
    ~LennardJonesParameterTable();


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LennardJonesParameterTable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void resize(int count);
    ///< Set the number of atom types whose information is stored in the table.

    void set(int type1, int type2, const LennardJonesParameters &params);
    // Set the parameters for one pair of atom types.

    const LennardJonesParameters &operator()(int type1, int type2) const {
      return myData[type1 * myCurrentSize + type2];
    }
    ///< Get the parameters for one pair of atom types.

    int size() const {return myCurrentSize;}
    ///< Get the current size of the table - 5/25/01 TMC

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    int myCurrentSize;
    LennardJonesParameters *myData;
    ///< Data contains a complete matrix, which is symmetric.
  };
  //________________________________________ INLINES
}
#endif /* not LENNARDJONESPARAMETERTABLE_H */
