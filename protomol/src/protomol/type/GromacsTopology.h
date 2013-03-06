/*  -*- c++ -*-  */
#ifndef GROMACSTOPOLOGY_H
#define GROMACSTOPOLOGY_H

#include <string>
#include <vector>

#include <protomol/type/Real.h>


namespace ProtoMol {

  //container holding GROMACS topology information

  class GromacsTopology {

   public:
    
     //data structures 
 
     struct Atom {

       Atom() : number(0), atom_type(""),residue_number(0),residue_name(""),
          atom_name(""),cgnr(0),charge(0),mass(0) {}

       Atom(int b,
          std::string a_type,
          int res_no,
          std::string res_name,
          std::string a_name,
          int c_gnr,
          Real c,
          Real m):
        number(b), atom_type(a_type), residue_number(res_no), residue_name(res_name),
        atom_name(a_name), cgnr(c_gnr), charge(c), mass(m) {}
   
       int number;
       std::string atom_type;
       int residue_number;
       std::string residue_name;
       std::string atom_name;
       int cgnr;
       Real charge;
       Real mass;

     };

     struct Bond {

       Bond() : number(0), atom1(0), atom2(0), funct(0){}

       Bond(int n, int a1, int a2, int f):
           number(n), atom1(a1), atom2(a2), funct(f) {}

       int number;
       int atom1;
       int atom2;
       int funct;

     };

     struct Pair {

       Pair() : number(0), atom1(0), atom2(0), funct(0) {}

       Pair(int n, int a1, int a2, int f) : number(n), atom1(a1), atom2(a2), funct(f) {}

       int number;
       int atom1;
       int atom2;
       int funct;

     };

     struct Angle {

       Angle() : number(0), atom1(0), atom2(0), atom3(0), funct(0) {}

       Angle(int n, int a1, int a2, int a3, int f) : number(n), atom1(a1), atom2(a2), atom3(a3), funct(f) {}

       int number;
       int atom1;
       int atom2;
       int atom3;
       int funct;

     };

     struct RBDihedral {

       RBDihedral() : number(0), atom1(0), atom2(0), atom3(0), atom4(0), funct(0) {}

       RBDihedral(int n, int a1, int a2, int a3, int a4, int f) :
            number(n), atom1(a1), atom2(a2), atom3(a3), atom4(a4), funct(f) {}

       int number;
       int atom1;
       int atom2;
       int atom3;
       int atom4;
       int funct;
     };

     struct Proper_Dihedral {

       Proper_Dihedral() : number(0), atom1(0), atom2(0), atom3(0), atom4(0), funct(0) {}

       Proper_Dihedral(int n, int a1, int a2, int a3, int a4, int f) :
            number(n), atom1(a1), atom2(a2), atom3(a3), atom4(a4), funct(f) {}

       int number;
       int atom1;
       int atom2;
       int atom3;
       int atom4;
       int funct;

     };

     //methods
     void clear();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::vector<GromacsTopology::Atom> atoms;
    std::vector<GromacsTopology::Bond> bonds;
    std::vector<GromacsTopology::Pair> pairs;
    std::vector<GromacsTopology::Angle> angles;
    std::vector<GromacsTopology::RBDihedral> rb_dihedrals;
    std::vector<GromacsTopology::Proper_Dihedral> dihedrals;     

  };

}

#endif /* GROMACSTOPOLOGYREADER_H */
