/*  -*- c++ -*-  */
#ifndef GROMACSPARAMETERS_H
#define GROMACSPARAMETERS_H

#include <string>
#include <vector>
#include <map>

#include <protomol/type/Real.h>

namespace ProtoMol {

  class GromacsParameters {

    public:
      //"defaults" record
      //This record tells us how to deal with exclusions and scaling
      //factors for LJ, Coulomb etc.
      //Check Gromacs manual file format section for details

      struct Default {

        Default() : nonbonded_function(0), combination_rule(0),
              gen_pairs(0), fudgeLJ(0), fudgeQQ(0){}

        Default(
          int n_f, int c_r, int gen_p, Real f_LJ, Real f_Q) :
           nonbonded_function(n_f), combination_rule(c_r),
           gen_pairs(gen_p), fudgeLJ(f_LJ), fudgeQQ(f_Q){} 
 
        int nonbonded_function;
        int combination_rule;
        int gen_pairs;
        Real fudgeLJ;
        Real fudgeQQ;
      };

      struct BondType {

        BondType() {}

        BondType(int n, std::string typea, std::string typeb, int f, Real fc, Real dist) :
              number(n), atom_typeA(typea), atom_typeB(typeb), funct(f), forceConstant(fc), distance(dist){}
       
        int number;  
        std::string atom_typeA;
        std::string atom_typeB;
        int funct;
        Real forceConstant; /* kb in the parameter file */
        Real distance; /* b0 in the parameter file */

      };

      struct AngleType {

        AngleType() {}

        AngleType(int n, std::string typea, std::string typeb, std::string typec,
            int f, Real fc, Real dist):
          number(n), atom_typeA(typea), atom_typeB(typeb), atom_typeC(typec),
          funct(f), forceConstant(fc), distance(dist){}

        int number;
        std::string atom_typeA;
        std::string atom_typeB;
        std::string atom_typeC;
        int funct;
        Real forceConstant; /* cth0 in the parameter file */
        Real distance; /* equilibrium angle th0 in parameter file */

      };


      struct Dihedral_base {

        Dihedral_base() {}
        Dihedral_base(int n, std::string typea, std::string typeb, std::string typec, std::string typed, int f,
           Real ph, Real k, int p, Real cc0, Real cc1, Real cc2, Real cc3, Real cc4, Real cc5):
           number(n), atom_typeA(typea), atom_typeB(typeb), atom_typeC(typec), atom_typeD(typed), funct(f),
           phase(std::vector<Real>(1,ph)), 
           kd(std::vector<Real>(1,k)), 
           pn(std::vector<int>(1,p)), 
           C0(cc0), C1(cc1), C2(cc2), C3(cc3), C4(cc4), C5(cc5) {}

        Dihedral_base(int n, std::string typea, std::string typeb, std::string typec, std::string typed, int f,
           std::vector<Real> ph, std::vector<Real> k, std::vector<int> p, Real cc0, Real cc1, 
           Real cc2, Real cc3, Real cc4, Real cc5) :
           number(n), atom_typeA(typea), atom_typeB(typeb), atom_typeC(typec), atom_typeD(typed), funct(f),
           phase(ph), kd(k), pn(p), C0(cc0), C1(cc1), C2(cc2), C3(cc3), C4(cc4), C5(cc5) {}

        int number;
        std::string atom_typeA;
        std::string atom_typeB;
        std::string atom_typeC;
        std::string atom_typeD;
        int funct;
        int multiplicity;
        std::vector<Real> phase;
        std::vector<Real> kd;
        std::vector<int> pn;
        Real C0;
        Real C1;
        Real C2;
        Real C3;
        Real C4;
        Real C5;

     };


     //Nonbonded LJ parameters
     struct Atom_base {

        Atom_base() {}
        Atom_base(int n, std::string b_type, Real m, Real c,
              std::string p_type, Real s, Real e) :
        number(n), bond_type(b_type), mass(m), charge(c),
        particle_type(p_type), sigma(s), epsilon(e) {}
     
        int number;
        std::string bond_type;
        Real mass;
        Real charge;
        std::string particle_type;
        Real sigma;
        Real epsilon;
   
     };


     typedef std::map<std::string, Atom_base> AtomType;

     //Generalized Born parameters
     struct GBParameters {

       GBParameters() {}
       GBParameters(std::string type, Real r , Real ig, Real ia, Real id,
          Real sg, Real sa, Real sd, Real GBd, Real ai) :
         atom_type_identifier(type), radius(r), igamma(ig), ialpha(ia),
         idelta(id), sgamma(sg), salpha(sa), sdelta(sd), GBdistcorr(GBd),
          a_i(ai) {}
     
       
       std::string atom_type_identifier;
       Real  radius;
       Real  igamma;
       Real  ialpha;
       Real  idelta;
       Real  sgamma;
       Real  salpha;
       Real  sdelta;
       Real  GBdistcorr; 
       Real  a_i;

     };

     //methods
     void clear();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::vector<GromacsParameters::Default> defaults;
    std::vector<GromacsParameters::BondType> bondTypes;
    std::vector<GromacsParameters::AngleType> angleTypes; 
    std::vector<GromacsParameters::Dihedral_base> dihedralTypes;
    AtomType atomTypes; 

    std::vector<GromacsParameters::GBParameters> gb_parameters;

  };

}

#endif /* GROMACSPARAMETERS_H */
