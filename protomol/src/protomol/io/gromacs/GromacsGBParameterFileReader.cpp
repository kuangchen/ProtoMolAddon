#include <protomol/io/gromacs/GromacsGBParameterFileReader.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>
 
 
using namespace ProtoMol::Report;
using namespace ProtoMol;
using namespace std;

GromacsGBParameterFileReader::GromacsGBParameterFileReader():
   Reader(NULL), gParams(NULL) {}

GromacsGBParameterFileReader::GromacsGBParameterFileReader(const string &filename):
   Reader(filename), gParams(NULL) {}

GromacsGBParameterFileReader::~GromacsGBParameterFileReader() {}

bool GromacsGBParameterFileReader::tryFormat() {

  if (!open()) return false;
  //cant find any header string which will determine file type.
  return true;
}

bool GromacsGBParameterFileReader::read(GromacsParameters &gP) {

  gParams = &gP;
 
  return (read());
}

bool GromacsGBParameterFileReader::read() {

   if (!tryFormat())
    return false;
  if (!open())
    return false;
 
  string line;
 
  while (!file.eof()) {

    line = getline();

    if ((equalStartNocase("#",line)) || (equalStartNocase("@",line))) {
        continue;
    }else {
        //Fill in data structure
       if (!line.empty()) {
        if (!read_gromacs_GB_Parameters(line)) return false;
       }
            
    }

  }

  report <<"No of GB records = "<<gParams->gb_parameters.size()<<endr;

  return true;
}

bool GromacsGBParameterFileReader::read_gromacs_GB_Parameters(string line) {

   GromacsParameters::GBParameters gb_param;

   string s;
   stringstream ss;

   ss << line;

   ss >> s;

   gb_param.atom_type_identifier = s;

   ss >> s;
   gb_param.radius = toReal(s);

   ss >> s;
   gb_param.igamma =  toReal(s);

   ss >> s;
   gb_param.ialpha =  toReal(s);

   ss >> s;
   gb_param.idelta =  toReal(s);

   ss >> s;
   gb_param.sgamma =  toReal(s);

   ss >> s;
   gb_param.salpha =  toReal(s);

   ss >> s;
   gb_param.sdelta =  toReal(s);

   ss >> s;
   gb_param.GBdistcorr = toReal(s);

   ss >> s;
   gb_param.a_i =  toReal(s);

   report << debug(2) <<gb_param.atom_type_identifier<<" "<<gb_param.radius<<" "<<gb_param.igamma<<" "<<gb_param.ialpha<<" "<<gb_param.idelta<<" "<<gb_param.sgamma<<" "<<gb_param.salpha<<" "<< gb_param.sdelta<<" "<<gb_param.GBdistcorr<<" "<<gb_param.a_i<<endr;

   gParams->gb_parameters.push_back(gb_param);

   return true;
}
