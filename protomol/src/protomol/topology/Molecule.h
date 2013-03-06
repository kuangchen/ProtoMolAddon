/*  -*- c++ -*-  */
#ifndef MOLECULE_H
#define MOLECULE_H

#include <protomol/type/Real.h>
#include <protomol/type/Vector3D.h>
#include <protomol/type/SimpleTypes.h>
#include <string>
#include <vector>
#include <set>

namespace ProtoMol {
  //________________________________________ Molecule

  /**
     This class defines the information for one molecule.  It contains the mass
     of the molecule, its center of mass position and momentum, and a list of
     the atoms on the molecule.
   */
  struct Molecule {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Molecule() : mass(0.0),
      position(0.0, 0.0,
               0.0),
      momentum(0.0, 0.0, 0.0), water(false), type(0), newtype(0), lambda(0.0) {}


    size_t size() {return atoms.size();}
    ///< Return number of atoms of this molecule

    int operator[](int i) {return atoms[i];}
    ///< Return atom number

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Real mass;
    ///< The mass of the molecule.

    Vector3D position;
    ///<  The xyz coordinate of the molecule's center of mass.

    Vector3D momentum;
    ///<  The xyz momentum of the molecule's center of mass.

    std::vector<int> atoms;
    ///<  List of the ID#s of the atoms on the molecule.

    std::vector<PairInt> pairs;
    ///<  List of pairs (bond, angle, dihedral and improper)

    bool water;
    ///< If molecule is a water molecule


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // iSG/oSG molecule parts are listed below
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    std::string name;
    ///< name of the molecule

    unsigned int type, newtype;
    ///< Type of the molecule...needed for ISG simulations
    ///< type is the molecule's type in the Lambda = 0 state;
    ///< newtype is the molecule's type in the Lambda = 1 state;

    Real lambda;
    ///< Parameter controlling the identity (or type) of the molecule.

    std::vector<int> bondList;
    std::vector<int> angleList;
    std::vector<int> dihedralList;
    std::vector<int> improperList;
    ///< list of the index #s of each bond, angle, dihedral,
    ///< and improper on this molecule...may not need this because of pairs 
    ///< above!!!!
  };
}
#endif /* MOLECULE_H */

