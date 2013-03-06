/* -*- c++ -*- */
#ifndef COULOMBSCPISMPARAMETERTABLE_H
#define COULOMBSCPISMPARAMETERTABLE_H

#include <protomol/type/Real.h>
#include <protomol/topology/CoulombSCPISMParameters.h>
#include <protomol/io/SCPISMReader.h>

#include <string>
#include <map>

namespace ProtoMol {
  //________________________________________ CoulombSCPISMParameterTable

  /**
   * This table contains the Lennard-Jones parameters for all atom type pairs.
   * It is guaranteed to be symmetric.
   */
  struct CoulombSCPISMParameterTable {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    //    CoulombSCPISMParameterTable() {}
    //    ~CoulombSCPISMParameterTable() {}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombSCPISMParameterTable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void populateTable();

    void displayTable();

    //    void resize(int count);
    ///< Set the number of atom types whose information is stored in the table.

    //    void set(int type1, const CoulombSCPISMParameters &params);
    // Set the parameters for one pair of atom types.

    //      int size() const {return myCurrentSize;}
    ///< Get the current size of the table - 5/25/01 TMC

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    std::map<std::string, CoulombSCPISMParameters> myData;
  };
  //________________________________________ INLINES
}
#endif /* not COULOMBSCPISMPARAMETERTABLE_H */
