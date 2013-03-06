#ifndef PORTGROMACSPARAMETERS_H
#define PORTGROMACSPARAMETERS_H

#include <iostream>
#include <protomol/type/Real.h>

//#define PROPER_DIHEDRAL 1
//#define RYCKERT_BELLEMAN_DIHEDRAL 3

//
//Porting Gromacs parameter into PSF and PAR structure (with required
//unit conversions etc.). We can then use any existing intergrator in
//ProtoMol with the GROMACS/AMBER force fields.
//

namespace ProtoMol {

  class PSF;
  class PAR;
  class GromacsTopology;
  class GromacsParameters;


  class PortGromacsParameters {

    private:
      static const int PROPER_DIHEDRAL = 1;
      static const int RYCKERT_BELLEMAN_DIHEDRAL = 3;

    public:
      PortGromacsParameters();
       ~PortGromacsParameters() {}

     bool Read_Basic_Gromacs_Parameters(PSF &psf, PAR &par, GromacsTopology &gTopo, GromacsParameters &gParams, 
            std::string filename, std::string pathname);

     bool Read_Gromacs_GB_Parameters(std::string filename);

     void Port_Parameters();

     void Port_GB_Parameters();

   private:
     void PortAtoms();
     void PortBonds();
     void PortAngles();
     void PortDihedrals();
     void PortImpropers();
     void PortPARBonds();
     void PortPARAngles();
     void PortPARNonbonded();
     void PortGromacsDihedrals();

    private:
      PSF *myPSF;
      PAR *myPAR;
      GromacsTopology *myGromacsTopo;
      GromacsParameters *myGromacsParam;

      

  };


}

#endif /* PORTGROMACSPARAMETERS_H */
