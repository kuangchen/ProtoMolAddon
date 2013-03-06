#include <protomol/io/gromacs/GromacsBondedParameterFileReader.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>

using namespace ProtoMol::Report;
using namespace ProtoMol;
using namespace std;

GromacsBondedParameterFileReader::GromacsBondedParameterFileReader():
  Reader(NULL), gParams(NULL) {}

GromacsBondedParameterFileReader::GromacsBondedParameterFileReader(const string &filename):
  Reader(filename), gParams(NULL) {}

GromacsBondedParameterFileReader::~GromacsBondedParameterFileReader(){

}

bool GromacsBondedParameterFileReader::tryFormat() {

  if (!open()) return false;
  //cant find any header string which will determine file type.
  return true;

}

bool GromacsBondedParameterFileReader::read(GromacsParameters &gP) {
  gParams = &gP;
  return (read());
 
}

bool GromacsBondedParameterFileReader::read() {

   if (!tryFormat())
    return false;
  if (!open())
    return false;

  string param_head;

  while (!file.eof()) {

    param_head = getline();

     if (equalStartNocase(";",param_head)) {
        //report << plain <<param_head<<endr;
        continue;
     }else if (equalStartNocase("[",param_head)){

        stringstream ss(param_head);
        string token;
        ss >> token;
        ss >> token;

        if (equal("bondtypes",token)) {
           read_bondtypes();
        } 

        if (equal("angletypes",token)) {
           read_angletypes();
        }
        
        if (equal("dihedraltypes",token)) {
           read_dihedrals_and_impropers();
        }
     }

  }

  return true;    

}

void GromacsBondedParameterFileReader::read_bondtypes(){

   string line;
 
   line = getline();

   // get rid of comment lines starting with ';'
   while (equalStartNocase(";",line)) {
       line = getline();
   }
 
   while ((!line.empty()) && (!file.eof())) {
      FillGromacsParameterBondTypes(line);
      line = getline();
   }

   if (file.eof()) report <<error<<"End of file reached while reading bondtype information."<<endr;

}

void GromacsBondedParameterFileReader::read_angletypes(){
 
   string line;
 
   line = getline();
 
   // get rid of comment lines starting with ';'
   while (equalStartNocase(";",line)) line = getline();
 
   while ((!line.empty()) && (!file.eof())) {
      FillGromacsParameterAngleTypes(line);
      line = getline();
   }
 
   if (file.eof()) report <<error<<"End of file reached while reading bondtype information."<<endr;
 
}

void GromacsBondedParameterFileReader::read_dihedrals_and_impropers(){
 
   string line;
 
   line = getline();
 
   // get rid of comment lines starting with ';'
   while (equalStartNocase(";",line)) line = getline();
 
   while ((!line.empty()) && (!file.eof())) {
      FillGromacsParameterDihedralTypes(line);
      line = getline();
   }
 
   //if (file.eof()) report <<error<<"End of file reached while reading bondtype information."<<endr;
 
}


void GromacsBondedParameterFileReader::FillGromacsParameterBondTypes(string line) {

   GromacsParameters::BondType myBondType;

   string s;
   stringstream ss(line);
 
   //report << line <<endr;
   //ss << removeBeginEndBlanks(line);

   myBondType.number = gParams->bondTypes.size();

   ss >> s;
   myBondType.atom_typeA = s;

   ss >> s;
   myBondType.atom_typeB = s;

   ss >> s;
   myBondType.funct = toInt(s); 

   ss >> s;
   myBondType.distance = toReal(s);

   ss >> s;
   myBondType.forceConstant = toReal(s);


   report << debug(2) <<"BondType in Gromacs paramete file: "<<myBondType.number<<" "<<myBondType.atom_typeA<<" "<<myBondType.atom_typeB<<" "<< myBondType.funct<<" "<<myBondType.forceConstant<<" "<<myBondType.distance<<endr;

   gParams->bondTypes.push_back(myBondType);

}

 
void GromacsBondedParameterFileReader::FillGromacsParameterAngleTypes(string line) {

   GromacsParameters::AngleType myAngleType;
 
   string s;
    stringstream ss;
 
   //report << line <<endr;
   ss << removeBeginEndBlanks(line);

   myAngleType.number = gParams->angleTypes.size();

   ss >> s;
   myAngleType.atom_typeA = s;

   ss >> s;
   myAngleType.atom_typeB = s;

   ss >> s;
   myAngleType.atom_typeC = s;

   ss >> s;
   myAngleType.funct = toInt(s);

   ss >> s;
   myAngleType.distance = toReal(s);

   ss >> s;
   myAngleType.forceConstant = toReal(s);
  
   report << debug(2) <<"AngleType in Gromacs paramete file: "<<myAngleType.number<<" "<<myAngleType.atom_typeA<<" "<<myAngleType.atom_typeB<<" "<< myAngleType.atom_typeC<<" "<<myAngleType.funct<<" "<<myAngleType.forceConstant<<" "<<myAngleType.distance<<endr;
 
   gParams->angleTypes.push_back(myAngleType);
}

void GromacsBondedParameterFileReader::FillGromacsParameterDihedralTypes(string line) {

   GromacsParameters::Dihedral_base myD;

   stringstream ss;
   string s;

   //report << line <<endr;
   ss << removeBeginEndBlanks(line);

   myD.number = gParams->dihedralTypes.size();

   ss >> s;
   myD.atom_typeA = s;

   ss >> s;
   myD.atom_typeB = s;

   ss >> s;
   myD.atom_typeC = s;

   ss >> s;
   myD.atom_typeD = s;

   ss >> s;
   myD.funct = toInt(s);

   if (myD.funct == 1) {

      myD.multiplicity = 1;

      ss >> s;
      //myD.phase = toReal(s);
      myD.phase.push_back(toReal(s));

      ss >> s;
      //myD.kd = toReal(s);
      myD.kd.push_back(toReal(s));

      ss >> s;
      //myD.pn = toReal(s);
      //myD.pn.push_back(toReal(s));
      myD.pn.push_back(toInt(s));

#if 0
      report << debug(2) <<"Dihedral type : funct "<<myD.funct<<" "<<myD.number<<" "<<myD.atom_typeA<<" "<<myD.atom_typeB<<" "<<myD.atom_typeC<<" "<<myD.atom_typeD<<" "<<myD.phase<<" "<<myD.kd<<" "<<myD.pn<<endr;
#endif


   }else if (myD.funct == 3) {

      ss >> s;
      myD.C0 = toReal(s);

      ss >> s;
      myD.C1 = toReal(s);

      ss >> s;
      myD.C2 = toReal(s);

      ss >> s;
      myD.C3 = toReal(s);

      ss >> s;
      myD.C4 = toReal(s);

      ss >> s;
      myD.C5 = toReal(s);

      report << debug(2) <<"Dihedral type : funct "<<myD.funct<<" "<<myD.number<<" "<<myD.atom_typeA<<" "<<myD.atom_typeB<<" "<<myD.atom_typeC<<" "<<myD.atom_typeD<<" "<<myD.C0<<" "<<myD.C1<<" "<<myD.C2<<" "<<myD.C3<<" "<<myD.C4<<" "<<myD.C5<<endr;


   }
   
   gParams->dihedralTypes.push_back(myD);
}
