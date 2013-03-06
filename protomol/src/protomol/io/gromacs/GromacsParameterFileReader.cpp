#include <protomol/io/gromacs/GromacsParameterFileReader.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>


using namespace ProtoMol::Report;
using namespace ProtoMol;
using namespace std;

GromacsParameterFileReader::GromacsParameterFileReader():
   Reader(NULL), gParams(NULL) {}

GromacsParameterFileReader::GromacsParameterFileReader(const string &filename):
   Reader(filename), gParams(NULL) {}

GromacsParameterFileReader::~GromacsParameterFileReader() {}


bool GromacsParameterFileReader::tryFormat() {

  if (!open()) return false;
  //cant find any header string which will determine file type.
  return true;

}

bool GromacsParameterFileReader::read(GromacsParameters &gP) {

  gParams = &gP;

  return (read());

}

bool GromacsParameterFileReader::read() {

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
     }else if(equalStartNocase("#include \"", param_head)) {
        string filename = ParseInclude(param_head);
        if (equalEnd("bon.itp",filename)) bondedParameterFile = filename;
        else if (equalEnd("nb.itp",filename)) nonbondedParameterFile = filename;
     }else if (equalStartNocase("[",param_head)){
        stringstream ss(param_head);
        string token;
        ss >> token;
        ss >> token;

        if (equal("defaults",token)) {
           read_defaults_information();
        }
      }
  } 

  return true;
}

void GromacsParameterFileReader::read_defaults_information() {

   string line;

   line = getline();

   // get rid of comment lines starting with ';'
   while (equalStartNocase(";",line)) line = getline();

   while ((!line.empty()) && (!file.eof())) {
      FillGromacsParameterDefaults(line);
      line = getline();
   }

   if (file.eof()) report <<error<<"End of file reached while reading defaults info"<<endr;

}

void GromacsParameterFileReader::FillGromacsParameterDefaults(string line) {

   GromacsParameters::Default myD;

   string s;
    stringstream ss;

   //report << line <<endr;
   ss << line;

   ss >> s;
   myD.nonbonded_function = toInt(s);

   ss >> s;
   myD.combination_rule = toInt(s);

   ss >> s;
   if (equal(s,"yes")) myD.gen_pairs = 1;
   else myD.gen_pairs = 0;

   ss >> s;
   myD.fudgeLJ = toReal(s);

   ss >> s;
   myD.fudgeQQ = toReal(s);

   report << debug(2) <<"Defaults in Gromacs paramete file: "<<myD.nonbonded_function<<" "<<myD.combination_rule<<" "<<myD.gen_pairs<<" "<<myD.fudgeLJ<<" "<<myD.fudgeQQ<<endr;

   gParams->defaults.push_back(myD);

}

//
//If there is an include directive, this function will parse it and
// return the included filename.
//
string GromacsParameterFileReader::ParseInclude(string s) {
 
  //report<<plain<<"String in ParseInclude "<<s<<endr;
  const vector<string> &vec_s = splitString(s);
  if (vec_s.size() != 2) report <<error<<"GromacsTopologyReader::Error reading include directive"<<endr;
 
  //report << plain <<vec_s[1]<<endr;
 
  string t(vec_s[1]);
 
  t.erase(0,1);
 
  int s_t = t.length();
  t.erase(s_t-1,1);
  report << debug(2) <<"ParseInclude : t "<<t<<endr;
  return t;
 
}

