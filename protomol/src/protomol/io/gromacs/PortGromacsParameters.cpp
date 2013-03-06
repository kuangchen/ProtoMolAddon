#include <protomol/io/gromacs/PortGromacsParameters.h>

#include <protomol/type/PSF.h>
#include <protomol/type/PAR.h>
#include <protomol/type/GromacsTopology.h>
#include <protomol/type/GromacsParameters.h>
#include <protomol/io/gromacs/GromacsTopologyReader.h>
#include <protomol/io/gromacs/GromacsParameterFileReader.h>
#include <protomol/io/gromacs/GromacsBondedParameterFileReader.h>
#include <protomol/io/gromacs/GromacsNonbondedParameterFileReader.h>
#include <protomol/io/gromacs/GromacsGBParameterFileReader.h>

#include <protomol/base/PMConstants.h>

#include <iterator>
#include <algorithm>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

typedef map<string, GromacsParameters::Atom_base>::iterator atomtype_iterator;

typedef vector<GromacsParameters::Dihedral_base>::iterator dihedral_iterator;

PortGromacsParameters::PortGromacsParameters() {

  myPSF = NULL;
  myPAR = NULL;
  myGromacsTopo = NULL;
  myGromacsParam = NULL;

}

bool PortGromacsParameters::Read_Basic_Gromacs_Parameters(PSF &psf, PAR &par, GromacsTopology &gTopo, 
   GromacsParameters &gParams, string filename, string pathname) {

  myPSF = &psf;
  myPAR = &par;
  myGromacsTopo = &gTopo;
  myGromacsParam = &gParams;

  GromacsTopologyReader myReader(filename); 

  if (!myReader.read(*myGromacsTopo)) {
    report << error <<"Error reading GROMACS topology file";
    return false;
  }

  string parameter_filename = myReader.GetParameterFilename();

  parameter_filename = pathname+"/"+parameter_filename;

  GromacsParameterFileReader gParamReader(parameter_filename);

  if (!gParamReader.read(*myGromacsParam)) {
     report << error <<"Error reading top level GROMACS force field parameter file"<<endr;
     return false;
  }

  string bonded_parameter_file = gParamReader.GetBondedParameterFile();

  string nonbonded_parameter_file = gParamReader.GetNonbondedParameterFile();

  bonded_parameter_file = pathname+"/"+bonded_parameter_file;

  nonbonded_parameter_file = pathname+"/"+nonbonded_parameter_file;


  report << debug(2)  << bonded_parameter_file<<endl;
  report << debug(2)  << nonbonded_parameter_file<<endl;

  GromacsBondedParameterFileReader gBondedParamReader(bonded_parameter_file);

  GromacsNonbondedParameterFileReader gNonbondedParamReader(nonbonded_parameter_file);

  if (!gBondedParamReader.read(gParams)) {
    report << error <<"Error reading bonded parameter file"<<endr;
    return false;
  }

  if (!gNonbondedParamReader.read(gParams)) {
    report << error <<"Error reading nonbonded parameter file"<<endr;
    return false;
  }

    
  Port_Parameters();

  if (myGromacsParam->defaults.size() == 0) {
     report << plain <<"Cant find parameters for modifying the 1-4 interactions with"
                <<" GROMACS/AMBER force field"<<endr;
  }else {
     par.fudgeQQ = myGromacsParam->defaults[0].fudgeQQ;
     par.fudgeLJ = myGromacsParam->defaults[0].fudgeLJ;
  }

  return true;

}

bool PortGromacsParameters::Read_Gromacs_GB_Parameters(string filename) {

  if (( myPAR == NULL) || (myGromacsParam == NULL)) {
     report << error <<"PAR and GromacsParameter data structures not available"<<endr;
     return false;
  }

  GromacsGBParameterFileReader gb_reader(filename);

  if (!gb_reader.read(*myGromacsParam)) {
     report << error <<"Error reading Generalized Born parameter file"<<endr;
     return false;
  }

  Port_GB_Parameters();

  return true;
}

void PortGromacsParameters::Port_GB_Parameters() {

  unsigned int num_records = myGromacsParam->gb_parameters.size();

  for (unsigned int i=0;i<num_records;i++) {
     PAR::GB_gromacs gb;

     string s = myGromacsParam->gb_parameters[i].atom_type_identifier;
     string atom_type_name = myGromacsParam->atomTypes[s].bond_type;

     gb.atom_name = atom_type_name;
     gb.radius = myGromacsParam->gb_parameters[i].radius;
     gb.igamma = myGromacsParam->gb_parameters[i].igamma;
     gb.ialpha = myGromacsParam->gb_parameters[i].ialpha;
     gb.idelta = myGromacsParam->gb_parameters[i].idelta;
     gb.sgamma = myGromacsParam->gb_parameters[i].sgamma;
     gb.salpha = myGromacsParam->gb_parameters[i].salpha;
     gb.sdelta = myGromacsParam->gb_parameters[i].sdelta;
     gb.GBdistcorr = myGromacsParam->gb_parameters[i].GBdistcorr;
     gb.a_i = myGromacsParam->gb_parameters[i].a_i;

     //pair<string,PAR::GB_gromacs> myPair(atom_type_name,gb);

    report << debug(2) <<atom_type_name << " "<<gb.radius<<endr;

     myPAR->gb_parameters[atom_type_name] = gb;
  }
}

void PortGromacsParameters::Port_Parameters() {

  PortAtoms();
  PortBonds();
  PortAngles();
  PortDihedrals();
  PortImpropers();
  PortPARBonds();
  PortPARAngles();
  PortPARNonbonded();
  PortGromacsDihedrals();

}


//What about unit conversion?
void PortGromacsParameters::PortAtoms() {

  //total no of atoms.
  unsigned int num_atoms = myGromacsTopo->atoms.size();

  for (unsigned int i=0;i<num_atoms;i++) {
     PSF::Atom myAtom;

     myAtom.number = i+1;
     myAtom.residue_sequence = myGromacsTopo->atoms[i].residue_number;
     myAtom.residue_name = myGromacsTopo->atoms[i].residue_name;
     myAtom.atom_name = myGromacsTopo->atoms[i].atom_name;
     
     /* why do they call this bond type in ffAMBER parameter file? */
     myAtom.atom_type = myGromacsParam->atomTypes[myGromacsTopo->atoms[i].atom_type].bond_type; 
     myAtom.charge = myGromacsTopo->atoms[i].charge;
     myAtom.mass = myGromacsTopo->atoms[i].mass;
     
     myPSF->atoms.push_back(myAtom); 

  }    

  report << debug(2)  <<" PSF atoms size "<<myPSF->atoms.size()<<", "<<num_atoms<<endl;

  //for(unsigned int i=0;i<num_atoms;i++) {
  
}

void PortGromacsParameters::PortBonds() {

  unsigned int num_bonds = myGromacsTopo->bonds.size();

  for (unsigned int i=0;i<num_bonds;i++) {

     PSF::Bond myBond;

     myBond.number = i+1;
     myBond.atom1 = myGromacsTopo->bonds[i].atom1;
     myBond.atom2 = myGromacsTopo->bonds[i].atom2;

     myPSF->bonds.push_back(myBond);
  }

  report << debug(2)  <<" PSF bonds size "<<myPSF->bonds.size()<<", "<<num_bonds<<endl;

}

void PortGromacsParameters::PortAngles() {

  unsigned int num_angles = myGromacsTopo->angles.size();

  for(unsigned int i=0;i<num_angles;i++) {
     PSF::Angle myAngle;

     myAngle.number = i+1;
     myAngle.atom1 = myGromacsTopo->angles[i].atom1;
     myAngle.atom2 = myGromacsTopo->angles[i].atom2;
     myAngle.atom3 = myGromacsTopo->angles[i].atom3;

     myPSF->angles.push_back(myAngle);
  }

  report << debug(2)  <<" PSF angle size "<<myPSF->angles.size()<<", "<<num_angles<<endl;
}

void PortGromacsParameters::PortDihedrals() {

  unsigned int num_dihedrals = myGromacsTopo->rb_dihedrals.size();


  for(unsigned int i=0;i<num_dihedrals;i++) {

     PSF::RBDihedral myD;

     myD.number = i+1;
     myD.atom1 = myGromacsTopo->rb_dihedrals[i].atom1;
     myD.atom2 = myGromacsTopo->rb_dihedrals[i].atom2;
     myD.atom3 = myGromacsTopo->rb_dihedrals[i].atom3;
     myD.atom4 = myGromacsTopo->rb_dihedrals[i].atom4;

     myPSF->rb_dihedrals.push_back(myD);
  }
  report << debug(2)  <<" PSF dihedral size "<<myPSF->rb_dihedrals.size()<<", "<<num_dihedrals<<endl;

}

void PortGromacsParameters::PortImpropers() {

  unsigned int num_impropers = myGromacsTopo->dihedrals.size();

  for(unsigned int i=0;i<num_impropers;i++) {
     PSF::Dihedral myD;

     myD.number = i+1;
     myD.atom1 = myGromacsTopo->dihedrals[i].atom1;
     myD.atom2 = myGromacsTopo->dihedrals[i].atom2;
     myD.atom3 = myGromacsTopo->dihedrals[i].atom3;
     myD.atom4 = myGromacsTopo->dihedrals[i].atom4;

     myPSF->dihedrals.push_back(myD);
  }

}

void PortGromacsParameters::PortPARBonds() {

  unsigned int k = myGromacsParam->bondTypes.size();

  PAR::Bond par_bond;

  for(unsigned int i=0;i<k;i++) {
     PAR::Bond par_bond;

     par_bond.number = i+1;
     par_bond.atom1 = myGromacsParam->bondTypes[i].atom_typeA;
     par_bond.atom2 = myGromacsParam->bondTypes[i].atom_typeB;

     par_bond.forceConstant = myGromacsParam->bondTypes[i].forceConstant
                                * Constant::KJ_KCAL * Constant::INV_NM_ANGSTROM * Constant::INV_NM_ANGSTROM * 0.5;
     par_bond.distance = myGromacsParam->bondTypes[i].distance * Constant::NM_ANGSTROM;
     myPAR->bonds.push_back(par_bond);
  }

}

void PortGromacsParameters::PortPARAngles() {

  unsigned int k = myGromacsParam->angleTypes.size();


  for(unsigned int i=0;i<k;i++) {

     PAR::Angle par_angle;
     par_angle.number = i+1;

     par_angle.atom1 = myGromacsParam->angleTypes[i].atom_typeA;
     par_angle.atom2 = myGromacsParam->angleTypes[i].atom_typeB;
     par_angle.atom3 = myGromacsParam->angleTypes[i].atom_typeC;
     par_angle.forceConstant = myGromacsParam->angleTypes[i].forceConstant*Constant::KJ_KCAL*0.5;
     //par_angle.forceConstant = myGromacsParam->angleTypes[i].forceConstant;
     par_angle.angleval = myGromacsParam->angleTypes[i].distance;
     par_angle.ub_flag = 0; /*means ignore ub parameters */
  
    myPAR->angles.push_back(par_angle);

  }

}

void PortGromacsParameters::PortPARNonbonded() {

   unsigned int i=0;


   report << debug(2)  <<"Nonbonded interactions..."<<endl;

   for(atomtype_iterator itr = myGromacsParam->atomTypes.begin(); itr !=  myGromacsParam->atomTypes.end();itr++) {
      PAR::Nonbonded par_nb;
      //sequence no?
      par_nb.number = i;
      i++;

      //atom type?
      //par_nb.atom = (*itr).first;
      //first = amberxx_xx...thats not the atom type.
      par_nb.atom = (*itr).second.bond_type;

      //epsilon
      //par_nb.epsilon = (*itr).second.epsilon*Constant::NM_ANGSTROM;
      par_nb.epsilon = (*itr).second.epsilon*Constant::KJ_KCAL;
       
      //sigma
      //par_nb.sigma = (*itr).second.sigma*Constant::KJ_KCAL;
      par_nb.sigma = (*itr).second.sigma*Constant::NM_ANGSTROM;

      report << debug(2)  <<"Atom "<<par_nb.atom<<", epsilon "<<par_nb.epsilon<<", "<<", sigma "<<par_nb.sigma<<endl;

      myPAR->nonbondeds.push_back(par_nb);
   }      

}

void PortGromacsParameters::PortGromacsDihedrals() {

  for(dihedral_iterator d = myGromacsParam->dihedralTypes.begin() ; d != myGromacsParam->dihedralTypes.end() ; d++) {

     int f = (*d).funct;

     if (f == RYCKERT_BELLEMAN_DIHEDRAL) {
        PAR::RBDihedral rb_dih;
        
         rb_dih.number = myPAR->rb_dihedrals.size()+1; 
         rb_dih.atom1 = (*d).atom_typeA;
         rb_dih.atom2 = (*d).atom_typeB;
         rb_dih.atom3 = (*d).atom_typeC;
         rb_dih.atom4 = (*d).atom_typeD;

         rb_dih.C0 = (*d).C0*Constant::KJ_KCAL;
         rb_dih.C1 = (*d).C1*Constant::KJ_KCAL;
         rb_dih.C2 = (*d).C2*Constant::KJ_KCAL;
         rb_dih.C3 = (*d).C3*Constant::KJ_KCAL;
         rb_dih.C4 = (*d).C4*Constant::KJ_KCAL;
         rb_dih.C5 = (*d).C5*Constant::KJ_KCAL;

         //check non zero force
         //if(rb_dih.C0 != 0 || rb_dih.C1 != 0 || rb_dih.C2 != 0 || 
         //   rb_dih.C3 != 0 || rb_dih.C4 != 0 || rb_dih.C5 != 0 ) {

         myPAR->rb_dihedrals.push_back(rb_dih);
         //}

     }else if (f == PROPER_DIHEDRAL) {

         PAR::Dihedral dih;

         dih.number = myPAR->dihedrals.size()+1;
         dih.atom1 = (*d).atom_typeA;
         dih.atom2 = (*d).atom_typeB;
         dih.atom3 = (*d).atom_typeC;
         dih.atom4 = (*d).atom_typeD;

        dih.multiplicity = 1;

        //dih.forceConstant.push_back((*d).kd*Constant::KJ_KCAL);
        //dih.forceConstant.resize((*d).kd.size());
        //copy((*d).kd.begin(),(*d).kd.end(),dih.forceConstant.begin());
        copy((*d).kd.begin(),(*d).kd.end(),back_insert_iterator<vector<Real> >(dih.forceConstant));

        for(unsigned int kk = 0 ;kk<dih.forceConstant.size();kk++) dih.forceConstant[kk] *= Constant::KJ_KCAL;

        //dih.periodicity.push_back((int)(*d).pn);
        //copy((*d).pn.begin(),(*d).pn.end(),dih.periodicity.begin());
        copy((*d).pn.begin(),(*d).pn.end(),back_insert_iterator<vector<int> >(dih.periodicity));

        //dih.phaseShift.push_back((*d).phase);
        //copy((*d).phase.begin(),(*d).phase.end(),dih.phaseShift.begin());
        copy((*d).phase.begin(),(*d).phase.end(),back_insert_iterator<vector<Real> >(dih.phaseShift));
       
        myPAR->dihedrals.push_back(dih);

     } else {
        // f = 2
     }

  }

}


