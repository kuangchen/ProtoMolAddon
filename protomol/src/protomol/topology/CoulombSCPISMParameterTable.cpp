#include <protomol/topology/CoulombSCPISMParameters.h>
#include <protomol/topology/CoulombSCPISMParameterTable.h>
#include <protomol/base/Report.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

void CoulombSCPISMParameterTable::populateTable() {
  CoulombSCPISMParameters temp;

//  if(!scpismFile.empty()) {

//    SCPISMReader mReader(scpismFile);

//    if(!mReader.read( myData ))
//      report << error << "Invalid SCPISM parameter file " << scpismFile << "." << endr;

//  } else {

  // Order:
  // name, alpha_i, hbond_factor, R_iw, r_cov, gamma_i, isHbonded

  std::string atom_name;

  atom_name = "H"; 
  temp.set(atom_name, 0.4906, 0.6259, 0.5000, 0.3700, 0.0052, PH);
  myData[atom_name] = temp; // polar H
  atom_name = "HC"; 
  temp.set(atom_name, 0.5300, 9.6746, 0.5000, 0.3700, 0.0052, PH);
  myData[atom_name] = temp; // N-ter H
  atom_name = "HA"; 
  temp.set(atom_name, 0.5300, 0.5000, 0.5000, 0.3700, 0.0052, NO);
  myData[atom_name] = temp; // nonpolar H
  atom_name = "HT"; 
  temp.set(atom_name, 0.4906, 2.3800, 0.5000, 0.3700, 0.0052, PH);
  myData[atom_name] = temp; // TIPS3P Water H
  atom_name = "HP"; 
  temp.set(atom_name, 0.5274, 0.5000, 0.5000, 0.3700, 0.0052, NO);
  myData[atom_name] = temp; // aromatic H
  atom_name = "HB"; 
  temp.set(atom_name, 0.4858, 0.5000, 0.5000, 0.3700, 0.0052, NO);
  myData[atom_name] = temp; // backbone H
  atom_name = "HR1"; 
  temp.set(atom_name, 0.4651, 0.5000, 0.5000, 0.3700, 0.0052, NO);
  myData[atom_name] = temp; //his he1, (+) his HG,HD2
  atom_name = "HR2"; 
  temp.set(atom_name, 0.4580, 0.5000, 0.5000, 0.3700, 0.0052, NO);
  myData[atom_name] = temp; //(+) his HE1
  atom_name = "HR3"; 
  temp.set(atom_name, 0.4893, 0.5000, 0.5000, 0.3700, 0.0052, NO);
  myData[atom_name] = temp; //neutral his HG, HD2
  atom_name = "HS"; 
  temp.set(atom_name, 0.4859, 0.5000, 0.5000, 0.3700, 0.0052, NO);
  myData[atom_name] = temp; // thiol H
  atom_name = "HA1"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.3700, 0.0052, NO);
  myData[atom_name] = temp; // for alkene; RHC=CR, wildcard Rvdw
  atom_name = "HA2"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.3700, 0.0052, NO);
  myData[atom_name] = temp; // for alkene; H2C=CR, wildcard Rvdw
  atom_name = "C"; 
  temp.set(atom_name, 0.5298, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; // polar C
  atom_name = "CA"; 
  temp.set(atom_name, 0.5192, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; // aromatic C
  atom_name = "CT1"; 
  temp.set(atom_name, 0.5209, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //aliphatic sp3 C for CH
  atom_name = "CT2"; 
  temp.set(atom_name, 0.5298, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //aliphatic sp3 C for CH2
  atom_name = "CT3"; 
  temp.set(atom_name, 0.5082, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //aliphatic sp3 C for CH3
  atom_name = "CPH1"; 
  temp.set(atom_name, 0.4538, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //his CG and CD2 carbons
  atom_name = "CPH2"; 
  temp.set(atom_name, 0.4974, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //his CE1 carbon
  atom_name = "CPT"; 
  temp.set(atom_name, 0.5217, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //trp C between rings
  atom_name = "CY"; 
  temp.set(atom_name, 0.5240, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //TRP C in pyrrole ring
  atom_name = "CP1"; 
  temp.set(atom_name, 0.4763, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //tetrahedral C (proline CA)
  atom_name = "CP2"; 
  temp.set(atom_name, 0.4599, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //tetrahedral C (proline CB/CG)
  atom_name = "CP3"; 
  temp.set(atom_name, 0.4627, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //tetrahedral C (proline CD)
  atom_name = "CC"; 
  temp.set(atom_name, 0.4970, 0.5000, 0.5000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //carbonyl C for asn,asp,gln,glu
  atom_name = "CD"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //carbonyl C for none amides, asp,glu,cter
  atom_name = "CPA"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //heme alpha-C
  atom_name = "CPB"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //heme beta-C
  atom_name = "CPM"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //heme meso-C
  atom_name = "CM"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //heme CO carbon
  atom_name = "CS"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //thiolate carbon
  atom_name = "CE1"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //for alkene; RHC=CR
  atom_name = "CE2"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO);
  myData[atom_name] = temp; //for alkene; H2C=CR
  atom_name = "N"; 
  temp.set(atom_name, 0.4894, 0.5000, 0.5000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; // proline N
  atom_name = "NR1"; 
  temp.set(atom_name, 0.4598, 0.5000, 0.5000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; // neutral his protonated ring N
  atom_name = "NR2"; 
  temp.set(atom_name, 0.4839, 19.4824, 0.6956,0.7400, 0.0052, PA);
  myData[atom_name] = temp; // neutral his unprotonated ring N
  atom_name = "NR3"; 
  temp.set(atom_name, 0.4742, 0.5000, 0.5000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; // charged his  ring N
  atom_name = "NH1"; 
  temp.set(atom_name, 0.5103, 0.5000, 0.5000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; // peptide N
  atom_name = "NH2"; 
  temp.set(atom_name, 0.5045, 0.5000, 0.5000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; // amide N
  atom_name = "NH3"; 
  temp.set(atom_name, 0.5300, 0.5000, 0.5000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; // ammonium N
  atom_name = "NC2"; 
  temp.set(atom_name, 0.5299, 0.5000, 0.5000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; // guanidinium N
  atom_name = "NY"; 
  temp.set(atom_name, 0.5284, 0.5000, 0.5000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; // TRP N in pyrrole ring
  atom_name = "NP"; 
  temp.set(atom_name, 0.4894, 0.5000, 0.5000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; // Proline ring NH2+ (N-terminal)
  atom_name = "NPH"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; //heme pyrrole N
  atom_name = "O"; 
  temp.set(atom_name, 0.4810, 4.0085, 0.5000, 0.7400, 0.0052, PA);
  myData[atom_name] = temp; // carbonyl oxygen
  atom_name = "OB"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; //carbonyl oxygen in acetic acid
  atom_name = "OC"; 
  temp.set(atom_name, 0.6000, 8.1838, 0.5000, 0.7400, 0.0052, PA);
  myData[atom_name] = temp; // carboxlate oxygen
  atom_name = "OH1"; 
  temp.set(atom_name, 0.4843, 3.3433, 0.5000, 0.7400, 0.0052, PA);
  myData[atom_name] = temp; // hydroxyl oxygen
  atom_name = "OS"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; // ester oxygen
  atom_name = "OT"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; // TIPS3P Water O
  atom_name = "OM"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.7400, 0.0052, NO);
  myData[atom_name] = temp; // heme CO/O2 oxygen
  atom_name = "S"; 
  temp.set(atom_name, 0.5054, 5.0037, 0.5000, 1.0400, 0.0052, PA);
  myData[atom_name] = temp; // sulphur
  atom_name = "SM"; 
  temp.set(atom_name, 0.5054, 0.5000, 0.5000, 1.0400, 0.0052, NO);
  myData[atom_name] = temp; // sulphur C-S-S-C type
  atom_name = "SS"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 1.0400, 0.0052, NO);
  myData[atom_name] = temp; // thiolate sulphur
  atom_name = "HE"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, NO);
  myData[atom_name] = temp; // helium
  atom_name = "NE"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, NO);
  myData[atom_name] = temp; // neon
  atom_name = "CAL"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, NO);
  myData[atom_name] = temp; // calcium 2+
  atom_name = "ZN"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, NO);
  myData[atom_name] = temp; // zinc (II) cation
  atom_name = "FE"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, NO);
  myData[atom_name] = temp; // heme iron 56
  atom_name = "DUM"; 
  temp.set(atom_name, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, NO);
  myData[atom_name] = temp; // dummy atom

}

void CoulombSCPISMParameterTable::displayTable() {
  report << debug(2) << "Coulomb SCPISM Parameter Table." << endl;
  
  map<string, CoulombSCPISMParameters>::const_iterator it;
  for (it = myData.begin(); it != myData.end(); it++) {
    report << debug(2) << it->first << "\t";
    report << it->second.alpha_i << "\t";
    report << it->second.hbond_factor << "\t";
    report << it->second.R_iw << "\t";
    report << it->second.r_cov << "\t";
    report << it->second.gamma_i << "\t";
    report << it->second.isHbonded << "\t";
    report << it->second.A_i << "\t";
    report << it->second.B_i << "\t";
    report << it->second.C_i << "\t";
    report << endl;
  }
}

