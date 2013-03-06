#include <protomol/io/gromacs/GromacsNonbondedParameterFileReader.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>


using namespace ProtoMol::Report;
using namespace ProtoMol;
using namespace std;

GromacsNonbondedParameterFileReader::GromacsNonbondedParameterFileReader():
  Reader(NULL), gParams(NULL) {}

GromacsNonbondedParameterFileReader::GromacsNonbondedParameterFileReader(const string &filename):
  Reader(filename), gParams(NULL) {}

GromacsNonbondedParameterFileReader::~GromacsNonbondedParameterFileReader(){

}

bool GromacsNonbondedParameterFileReader::tryFormat() {


  if (!open()) return false;
  //cant find any header string which will determine file type.
  return true;

}

bool GromacsNonbondedParameterFileReader::read(GromacsParameters &gP) {

  gParams = &gP;
  return (read());

}

bool GromacsNonbondedParameterFileReader::read() {

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

        if (equal("atomtypes",token)) {
           read_atomtypes();
        }

    }

  }

  return true;
}

void GromacsNonbondedParameterFileReader::read_atomtypes() {

   string line;

   line = getline();

   // get rid of comment lines starting with ';'
   while (equalStartNocase(";",line)) line = getline();

   while ((!line.empty()) && (!file.eof())) {
      FillGromacsParameterAtomTypes(line);
      line = getline();
   }

   //if (file.eof()) report <<error<<"End of file reached while reading atomtype information."<<endr;

}

void GromacsNonbondedParameterFileReader::FillGromacsParameterAtomTypes(string line) {

  if (equalStartNocase(";",line)) return;

  GromacsParameters::Atom_base atom_base;

   string s;
    stringstream ss;

   //report << line <<endr;
   ss << removeBeginEndBlanks(line);

   string str;

   ss >> str;

   ss >> s;
   atom_base.bond_type = s;

   ss >> s;
   atom_base.mass = toReal(s);

   ss >> s;
   atom_base.charge = toReal(s);

   ss >> s;
   atom_base.particle_type = s;

   ss >> s;
   atom_base.sigma = toReal(s);

   ss >> s;
   atom_base.epsilon = toReal(s);

   pair<string,GromacsParameters::Atom_base> myPair(str,atom_base);

   report << debug(2) <<" Gromacs nonbonded parameters : atom_type "<<str<<" "<<atom_base.bond_type<<" "<<atom_base.mass<<" "<<atom_base.charge<<" "<<atom_base.particle_type<<" "<<atom_base.sigma<<" "<<atom_base.epsilon<<endr;

   gParams->atomTypes.insert(myPair);

}
