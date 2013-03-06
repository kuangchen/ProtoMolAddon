/*  -*- c++ -*-  */
#ifndef PAR_H
#define PAR_H

#include <string>
#include <vector>
#include <map>

#include <protomol/type/Real.h>
#include <protomol/base/Report.h>

namespace ProtoMol {
  //_________________________________________________________________PAR
  /**
   * Container class for Charmm19/28/XPLOR parameters @n
   *
   * NB:
   * - angles are always kept in degrees
   * - simga's for Nonbonded are asummed to be Charmm28
   */
  class PAR {
  public:
    //______________________________________________________________________Bond
    /**
     * This structure holds data for a bond, including the bond number, the two
     * atoms involved, the force constant and the distance.
     */
    struct Bond {
      Bond() {}
      Bond(int a, std::string b, std::string c, Real d,
           Real e) : number(a), atom1(b), atom2(c), forceConstant(d),
        distance(e) {}


      int number;         ///< bond number
      std::string atom1;  ///< atom 1 number
      std::string atom2;  ///< atom 2 number
      Real forceConstant; ///< force constant
      Real distance;      ///< distance

      friend Report::MyStreamer &operator<<(Report::MyStreamer &OS,
                                            const Bond &p);
    };

    //___________________________________________________________________Angle
    /**
     * This structure holds data for an angle.  There is the angle number, the
     * three atoms involved, the force constant, the actual value of the angle,
     * and the Urey-Bradley constants if they exist.  The ub_flag will be 1 if
     *  they exist and 0 otherwise.
     */
    struct Angle {
      Angle() {}
      Angle(int a, std::string b, std::string c, std::string d, Real e,
            Real f, bool g, Real h, Real i) :
        number(a), atom1(b), atom2(c), atom3(d), forceConstant(e), angleval(f),
        ub_flag(g), k_ub(h), r_ub(i) {}


      int number;         ///< angle number
      std::string atom1;  ///< atom 1 number
      std::string atom2;  ///< atom 2 number
      std::string atom3;  ///< atom 3 number
      Real forceConstant; ///< force constant
      Real angleval;      ///< angle value
      /**
       * Urey-Bradley flag - '1' if there are Urey-Bradley constants following
       * If '0', ignore the next two data members
       */
      bool ub_flag;
      Real k_ub;      ///< Urey-Bradley force constant
      Real r_ub;      ///< Urey-Bradley radius

      friend Report::MyStreamer &operator<<(Report::MyStreamer &OS,
                                            const Angle &p);
    };


    //_______________________________________________________________Dihedral
    /** This structure holds data for a dihedral, consisting of the number, the
     * four atoms involved, the multiplicity (default = 1), and the force
     * constant, the periodicity, and the phase shift
     */
    struct Dihedral {
      Dihedral() {}
      Dihedral(int a, std::string b, std::string c, std::string d,
               std::string e, Real f, int g, Real h) :
        number(a), atom1(b), atom2(c), atom3(d), atom4(e), multiplicity(1),
        forceConstant(std::vector<Real>(1, f)),
        periodicity(std::vector<int>(1, g)),
        phaseShift(std::vector<Real>(1, h)) {}

      Dihedral(int a, std::string b, std::string c, std::string d,
               std::string e, int f, std::vector<Real> g, std::vector<int> h,
               std::vector<Real> i) :
        number(a), atom1(b), atom2(c), atom3(d), atom4(e), multiplicity(f),
        forceConstant(g), periodicity(h), phaseShift(i) {}


      int number;         ///< dihedral number
      std::string atom1;  ///< atom 1 number
      std::string atom2;  ///< atom 2 number
      std::string atom3;  ///< atom 3 number
      std::string atom4;  ///< atom 4 number
      int multiplicity;   ///< multiplicity
      std::vector<Real> forceConstant;  ///< force constant
      std::vector<int> periodicity;     ///< periodicity
      std::vector<Real> phaseShift;     ///< phase shift

      friend Report::MyStreamer &operator<<(Report::MyStreamer &OS,
                                            const Dihedral &p);
    };


    //__________________________________________________________________Improper
    /** This structure holds data for an improper.  The data held is the same 
     * as that for a dihedral - the number, four atoms involved, the force
     * constant, the periodicity, and the phase shift
     */
    struct Improper {
      Improper() {}
      Improper(int a, std::string b, std::string c, std::string d,
               std::string e, Real f, int g, Real h) :
        number(a), atom1(b), atom2(c), atom3(d), atom4(e), forceConstant(f),
        periodicity(g), phaseShift(h) {}

      int number;         ///< improper number
      std::string atom1;  ///< atom 1 number
      std::string atom2;  ///< atom 2 number
      std::string atom3;  ///< atom 3 number
      std::string atom4;  ///< atom 4 number
      Real forceConstant; ///< force constant
      int periodicity;    ///< periodicity
      Real phaseShift;    ///< phase shift

      friend Report::MyStreamer &operator<<(Report::MyStreamer &OS,
                                            const Improper &p);
    };

    struct RBDihedral {

      RBDihedral() {}
      RBDihedral(int a, std::string b, std::string c, std::string d,
           std::string e, Real cc0, Real cc1, Real cc2, Real cc3, Real cc4, Real cc5):
           number(a), atom1(b), atom2(c), atom3(d), atom4(e),
           C0(cc0), C1(cc1), C2(cc2), C3(cc3), C4(cc4), C5(cc5) {} 

      int number;         ///< dihedral number
      std::string atom1;  ///< atom 1 number
      std::string atom2;  ///< atom 2 number
      std::string atom3;  ///< atom 3 number
      std::string atom4;  ///< atom 4 number
      Real C0;
      Real C1;
      Real C2;
      Real C3;
      Real C4;
      Real C5;

    };

    struct GB_gromacs {

      GB_gromacs() {}
      GB_gromacs(std::string a,Real r , Real ig, Real id, Real ia, Real sg, Real sa, Real sd,
               Real GBd, Real ai) : atom_name(a),
         radius(r),igamma(ig),ialpha(ia),idelta(id),sgamma(sg),salpha(sa),
         sdelta(sd), GBdistcorr(GBd), a_i(ai){} 

       std::string atom_name;
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

    //_________________________________________________________________Nonbonded
    /// This structure holds data for a nonbonded - including
    struct Nonbonded {
      Nonbonded() {}
      Nonbonded(int a, std::string b, Real c, Real d, Real e, bool f, bool g,
                Real h, Real i, Real j, bool k) :
        number(a), atom(b), polarizability(c), epsilon(d), sigma(e),
        negative(f), vdw(g), polarizability2(h), epsilon14(i), sigma14(j),
        negative2(k) {}

      int number;          ///< nonbonded number
      std::string atom;    ///< atom number
      Real polarizability; ///< polarizability or ignore (see description of
                           ///<negative below), default to zero
      Real epsilon;        ///< well depth or number of effective electrons
                           ///<(see description of negative below)
      Real sigma;          ///< minimum radius divided by 2
      /**
       * flag for if the second term is negative - if so, second_term = epsilon
       * or well-depth and the first term is ignored, otherwise second_term =
       * number of effective electrons and the first term is the polarizability
       * default to true
       */
      bool negative;
      /**
       *	 default to true - likely there will be epsilon 1:4 and sigma 1:4,
       * flag to see if there is a second set, indicating VDW parameters
       */
      bool vdw;
      Real polarizability2; ///< VDW parameter polarizability, default to zero
      Real epsilon14;       ///< VDW parameter well depth or number of effective
                            ///<electrons (see above)
      Real sigma14;         ///< VDW parameter minimum radius divided by 2
      bool negative2;       ///< flag for a negative VDW paramenter second term
                            ///< (see above), default to true


      static const Real SIGMA_CHARMM19_TO_CHARMM28;
      static const Real SIGMA_CHARMM28_TO_CHARMM19;

      friend Report::MyStreamer &operator<<(Report::MyStreamer &OS,
                                            const Nonbonded &p);
    };

    //_________________________________________________________________Nbfix
    /**
     * This structure holds data for atom pairs with modifiable VDW
     *  interactions. Data includes the number, two atoms, epsilon, sigma,
     *  epsilon 14, and sigma 14.
     */
    struct Nbfix {
      Nbfix() {}
      Nbfix(int a_, std::string b_, std::string c, Real d, Real e, 
            Real f, Real g)
        : number(a_), atom1(b_), atom2(c), a(d), b(e), a14(f), b14(g) {}

      int number;         ///< nbfi number
      std::string atom1;  ///< atom 1 number
      std::string atom2;  ///< atom 2 number
      Real a;             ///< epsilon
      Real b;             ///< sigma
      Real a14;           ///< epsilon 14
      Real b14;           ///< sigma 14

      friend Report::MyStreamer &operator<<(Report::MyStreamer &OS,
                                            const Nbfix &p);
    };
    //_________________________________________________________________Hbond
    /// This structure holds data for hydrogen bonds, including the well depth
    /// and the minimum radius
    struct Hbond {
      Hbond() {}
      Hbond(int a, std::string b, std::string c, Real d,
            Real e) : number(a), atom1(b), atom2(c), emin(d), rmin(e) {}

      int number;        ///< hydrogen bond number
      std::string atom1; ///< atom 1 number
      std::string atom2; ///< atom 2 number
      Real emin;         ///< well depth
      Real rmin;         ///< minimum radius

      friend Report::MyStreamer &operator<<(Report::MyStreamer &OS,
                                            const Hbond &p);
    };

    //_________________________________________________________________PAR

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// Supported Charmm/XPLOR type
    enum CharmmTypeEnum {
      UNDEFINED,
      CHARMM28,
      CHARMM19
    };


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class PAR
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void clear();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PAR container
    std::vector<PAR::Bond> bonds;
    std::vector<PAR::Angle> angles;
    std::vector<PAR::Dihedral> dihedrals;
    std::vector<PAR::Improper> impropers;
    std::vector<PAR::Nonbonded> nonbondeds;
    std::vector<PAR::Nbfix> nbfixs;
    std::vector<PAR::Hbond> hbonds;

    //for GROMACS
    std::vector<PAR::RBDihedral> rb_dihedrals;
    std::map<std::string,PAR::GB_gromacs> gb_parameters;

    Real fudgeLJ;
    Real fudgeQQ;

  };
}
#endif /* PAR_H */
