/* -*- c++ -*- */
#ifndef COULOMBSCPISM_H
#define COULOMBSCPISM_H

#include <protomol/type/Real.h>
#include <protomol/topology/AtomType.h>

namespace ProtoMol {
  //________________________________________ CoulombSCPISMParameters

  /// The Coulomb SCPISM for an atom type
  struct CoulombSCPISMParameters {
  public:

    CoulombSCPISMParameters(std::string s = "", Real al = 0.0, Real hb = 0.5, Real r = 0.5,
                            Real rc = 0.37, Real ga = 0.52,
                            Hbonded is = NO)
      : atom_name(s), alpha_i(al), hbond_factor(hb), R_iw(r),
        r_cov(rc), gamma_i(ga), isHbonded(is) {}
    
    void set(std::string s = "", Real al = 0.0, Real hb = 0.5, Real r = 0.5, 
             Real rc = 0.37, Real ga = 0.0052, Hbonded is = NO) {
      atom_name = s;
      alpha_i = al;
      hbond_factor = hb;
      R_iw = r;
      r_cov = rc;
      gamma_i = ga;
      isHbonded = is;
      //find SCPISM paramiters by base atom type
      char typ;
      if(s != "HE" && s!= "NE" && s != "CAL")
        typ = (s.c_str())[0];
      else
        typ = ' ';
      switch(typ){
        case 'C': A_i = 31.0; B_i = 7.92; C_i = 0.32;
                  break;
        case 'O': A_i = 40.0; B_i = 8.52; C_i = 0.29;
                  break;
        case 'N': A_i = 60.0; B_i = 9.66; C_i = 0.22;
                  break;
      case 'P': A_i = 60.0; B_i = 9.66; C_i = 0.22; // JAIP 9/29/09
	          break;
        case 'S': A_i = 60.0; B_i = 9.10; C_i = 0.22;
                  break;
        case 'H': if(isHbonded == PH){
                    A_i = 1.20; B_i = 3.68; C_i = 0.80;
                  }else{
                    A_i = 17.00; B_i = 9.00; C_i = 0.50;
                  }
                  break;
        default: A_i = 0.0; B_i = 0.0; C_i = 0.0;
                  break;
      }

    }

    std::string atom_name; //text name of atom
    Real alpha_i; // Alpha_i controls slope of D(r) around atom type i
    Real hbond_factor; // hbond_factor (Polar H) * hbond_factor (PA) controls
                       // increment of Born radius to correct h bonding strength
    Real R_iw; // extension of Born radius R_iw to obtain R_ip
    Real r_cov; // Covalent radius (only values for C,N,O,S,H)
    Real gamma_i;  // hydrophobic energy term
    Hbonded isHbonded; // whether atom is involved in H-bonding
    Real A_i; // A_i in (2) on SCPISM.DOC
    Real B_i; // B_i in (2) on SCPISM.DOC
    Real C_i; // C_i in (2) on SCPISM.DOC
  };
}
#endif /* not COULOMBSCPISM_H */
