/*  -*- c++ -*-  */
#ifndef PSF_H
#define PSF_H

#include <string>
#include <vector>

#include <protomol/type/Real.h>

namespace ProtoMol {
  //_________________________________________________________________PSF
  /**
   * PSF container holding the topology of the system
   */
  class PSF {
  public:
    //___________________________________________________________________PsfAtom
    /**
     * This class holds data for a basic .psf atom.  The data stored includes:
     * atom number, segment identifier, residue sequence, residue name, atom 
     * name, atom2 name, charge and mass.
     */
    struct Atom {
      Atom() : number(0), seg_id(""), residue_sequence(0), residue_name(""),
        atom_name(""), atom_type(""), charge(0.0), mass(0.0), identity(0) {}
      Atom(int a,
           std::string b,
           int c,
           std::string d,
           std::string e,
           std::string f,
           Real g,
           Real h,
           unsigned int i) :
        number(a), seg_id(b), residue_sequence(c), residue_name(d),
        atom_name(e), atom_type(f), charge(g), mass(h), identity(i) {}

      int number;               ///< atom number
      std::string seg_id;       ///< segment identifier
      int residue_sequence;     ///< residue sequence
      std::string residue_name; ///< residue name
      std::string atom_name;    ///< atom name
      std::string atom_type;    ///< atom type
      Real charge;              ///< charge [e]
      Real mass;                ///< mass [amu]
      unsigned int identity;    ///< atom's identity (for iSGMD)
    };


    //___________________________________________________________________Bond
    /**
     * This structure holds data for a basic .psf bond.  The data stored
     * includes the bond number and the numbers of the two atoms involved.
     */
    struct Bond {
      Bond() : number(0), atom1(0), atom2(0) {}
      Bond(int a, int b, int c) : number(a), atom1(b), atom2(c) {}

      int number;     ///< bond number
      int atom1;      ///< bonded atom 1 number
      int atom2;      ///< bonded atom 2 number
    };


    //___________________________________________________________________Angle
    /**
     * This structure holds data for a basic .psf angle.  The data stored
     * includes the angle number and the numbers of the three atoms involved.
     */
    struct Angle {
      Angle() : number(0), atom1(0), atom2(0), atom3(0) {}
      Angle(int a, int b, int c,
            int d) : number(a), atom1(b), atom2(c), atom3(d) {}

      int number;     ///< angle number
      int atom1;      ///< angle atom 1 number
      int atom2;      ///< angle atom 2 number
      int atom3;      ///< angle atom 3 number
    };


    //__________________________________________________________________Dihedral
    /**
     * This structure holds data for a basic .psf dihedral.  The data stored
     * includes the dihedral number and the numbers of the four atoms involved.
     */
    struct Dihedral {
      Dihedral() : number(0), atom1(0), atom2(0), atom3(0), atom4(0) {}
      Dihedral(int a, int b, int c, int d,
               int e) : number(a), atom1(b), atom2(c), atom3(d), atom4(e) {}

      int number;     ///< dihedral number
      int atom1;      ///< dihedral atom 1 number
      int atom2;      ///< dihedral atom 2 number
      int atom3;      ///< dihedral atom 3 number
      int atom4;      ///< dihedral atom 4 number
    };

    //__________________________________________________________________Improper
    /**
     * This structure holds data for a basic .psf improper.  The data stored
     * includes the improper number and the numbers of the four atoms involved.
     */
    struct Improper {
      Improper() : number(0), atom1(0), atom2(0), atom3(0), atom4(0) {}
      Improper(int a, int b, int c, int d,
               int e) : number(a), atom1(b), atom2(c), atom3(d), atom4(e) {}

      int number;     ///< improper number
      int atom1;      ///< improper atom 1 number
      int atom2;      ///< improper atom 2 number
      int atom3;      ///< improper atom 3 number
      int atom4;      ///< improper atom 4 number
    };

    struct RBDihedral {

      RBDihedral() : number(0), atom1(0), atom2(0), atom3(0), atom4(0) {}
      RBDihedral(int a, int b, int c, int d,
               int e) : number(a), atom1(b), atom2(c), atom3(d), atom4(e) {}
      

     int number;     ///< dihedral number
     int atom1;      ///< dihedral atom 1 number
     int atom2;      ///< dihedral atom 2 number
     int atom3;      ///< dihedral atom 3 number
     int atom4;      ///< dihedral atom 4 number
      
    };

    //_____________________________________________________________________Donor
    /**
     * This structure holds data for a basic .psf donor.  The data stored
     * includes the donor number and the numbers of the two donors involved
     */
    struct Donor {
      Donor() : number(0), atom1(0), atom2(0) {}
      Donor(int a, int b, int c) : number(a), atom1(b), atom2(c) {}

      int number;     ///< donor number
      int atom1;      ///< donor atom 1 number
      int atom2;      ///< donor atom 2 number
    };

    //__________________________________________________________________Acceptor
    /**
     * This structure holds data for a basic .psf acceptor.  The data stored
     * includes the acceptor number and the numbers of the two acceptors.
     */
    struct Acceptor {
      Acceptor() : number(0), atom1(0), atom2(0) {}
      Acceptor(int a, int b, int c) : number(a), atom1(b), atom2(c) {}

      int number;     ///< acceptor number
      int atom1;      ///< acceptor atom 1 number
      int atom2;      ///< acceptor atom 2 number
    };

    //_________________________________________________________________Nonbonded
    /**
     * This structure holds data for a basic .psf nonbonded.  The data stored
     * includes the number of the nonbonded structure as well as the number of
     * the nonbonded atom.
     */
    struct Nonbonded {
      Nonbonded() : number(0), atom1(0) {}
      Nonbonded(int a, int b) : number(a), atom1(b) {}

      int number;     ///< nonbonded structure number
      int atom1;      ///< nonbonded atom number
    };

    //______________________________________________________________________Ngrp
    /**
     * This structure holds data for a basic .psf ngrp.  The data stored
     * includes the number of the ngrp as well as the numbers of the three
     * atoms involved.
     */
    struct Ngrp {
      Ngrp() : number(0), atom1(0), atom2(0), atom3(0) {}
      Ngrp(int a, int b, int c,
           int d) : number(a), atom1(b), atom2(c), atom3(d) {}

      int number;     ///< ngrp number
      int atom1;      ///< ngrp atom 1 number
      int atom2;      ///< ngrp atom 2 number
      int atom3;      ///< ngrp atom 3 number
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class PSF
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void clear();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::vector<PSF::Atom> atoms;
    std::vector<PSF::Bond> bonds;
    std::vector<PSF::Angle> angles;
    std::vector<PSF::Dihedral> dihedrals;
    std::vector<PSF::Improper> impropers;
    std::vector<PSF::Donor> donors;
    std::vector<PSF::Acceptor> acceptors;
    std::vector<PSF::Nonbonded> nonbondeds;
    std::vector<PSF::Ngrp> ngrp;

    //Amber style dihedrals
    std::vector<PSF::RBDihedral> rb_dihedrals;

  };
}

#endif /* PSF_H */
