#include <protomol/io/gromacs/GromacsTopologyReader.h>

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____GromacsTopologyReader

GromacsTopologyReader::GromacsTopologyReader() :
   Reader(), gTopo(NULL) {}

GromacsTopologyReader::GromacsTopologyReader(const std::string& filename):
   Reader(filename), gTopo(NULL){}

GromacsTopologyReader::~GromacsTopologyReader() {

}
bool GromacsTopologyReader::tryFormat() {

   if (!open()) return false;
   string gromacs_topo_head;
   file >> gromacs_topo_head;

   return file.good() && equalStartNocase(";", gromacs_topo_head);

}

bool GromacsTopologyReader::read(GromacsTopology &gT) {

  gTopo = &gT;

  return (read());

}

bool GromacsTopologyReader::read() {

   if (!tryFormat())
    return false;
  if (!open())
    return false;

  //process header
  string topo_head;
  //topo_head = getline();

  while (!file.eof()) {

    topo_head = getline();
    if (equalStartNocase(";",topo_head)) {
       //report << plain <<topo_head<<endr;
       continue;
    }else if (equalStartNocase("#include \"", topo_head)) {
     string parm_filename = ParseInclude(topo_head);
     parameter_files.push_back(parm_filename);
    }else if (equalStartNocase("[",topo_head)) {
      stringstream ss(topo_head);
      string token;
      ss >> token;
      ss >> token;

      report<< debug(2) <<"First token "<<token<<endr;
     
      if (equal("atoms",token)) {
           read_atom_information();
           report <<debug(2)<<"Done reading atoms from topology"<<endr;
      }

      if (equal("bonds",token)) {
          read_bond_information();
           report <<debug(2)<<"Done reading bonds from topology"<<endr;
      }

      if (equal("angles",token)) {
         read_angle_information();
         report <<debug(2)<<"Done reading angles from topology"<<endr;
      }

      if (equal("dihedrals",token)) {
         read_dihedral_information();
         report <<debug(2)<<"Done reading dihedrals from topology"<<endr;
      }

    }
 
  }

#if 0
  for (unsigned int ii = 0; ii<parameter_files.size();ii++) report << debug(2) <<parameter_files[ii]<<endl; 
#endif
     
  return true;

}



//
//If there is an include directive, this function will parse it and
// return the included filename.
//
string GromacsTopologyReader::ParseInclude(string s) {

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

void GromacsTopologyReader::read_atom_information() {

  string line;

  line = getline();

  // get rid of comment lines starting with ';'
  while (equalStartNocase(";",line)) line = getline();

  while ((!line.empty()) && (!file.eof())) {
     FillGromacsTopologyAtom(line);
     line = getline();
  }

  if (file.eof()) report <<error<<"End of file reached while reading atom info"<<endr;

}

void GromacsTopologyReader::read_bond_information() {

   string line = getline();
  
  // get rid of comment lines starting with ';'
  while (equalStartNocase(";",line)) line = getline();


  while ((!line.empty()) && (!file.eof())) {
     FillGromacsTopologyBond(line);
     line = getline();
  }

}

void GromacsTopologyReader::read_angle_information(){

  string line = getline();

  // get rid of comment lines starting with ';'
  while (equalStartNocase(";",line)) line = getline();

  while ((!line.empty()) && (!file.eof())) {
    FillGromacsTopologyAngle(line);
    line = getline();
  }

}

void GromacsTopologyReader::read_dihedral_information() {

  string line = getline();

 // get rid of comment lines starting with ';'
  while (equalStartNocase(";",line)) line = getline();


  while ((!line.empty()) && (!file.eof())) {
    FillGromacsTopologyDihedrals(line);
    line = getline();
  }
 
} 
   
void GromacsTopologyReader::FillGromacsTopologyAtom(string line) {

    GromacsTopology::Atom myAtom;
    string s;
     stringstream ss;
     //report << line <<endr;
     ss << line;
     ss >> s;
     myAtom.number = toInt(s);
     ss >> s; /* atom type ...need to read parameter file */
  
     myAtom.atom_type = s;    

     ss >> s; /* residue sequence*/
     myAtom.residue_number = toInt(s);

     ss >> s; /* residue name */
     myAtom.residue_name = s;

     ss >> s; /* atom name */
     myAtom.atom_name = s;

     ss >> s ; /* cgnr */
     myAtom.cgnr = toInt(s);

     ss >> s ; /* charge */

     myAtom.charge = toReal(s);

     ss >> s; /* mass */

     myAtom.mass = toReal(s);

     report << debug(2) << "Atom information in GromacsTopology::Atom : "<<myAtom.number<<" "<<myAtom.atom_type<<" "<<myAtom.residue_number<<" "<<myAtom.residue_name<<" "<<myAtom.atom_name<<" "<<myAtom.cgnr<<" "<<myAtom.charge<<" "<<myAtom.mass<<endl;

     gTopo->atoms.push_back(myAtom);
}

void GromacsTopologyReader::FillGromacsTopologyBond(string line) {

   GromacsTopology::Bond myBond;

   string s;
   stringstream ss;

   //report << plain << line <<endr;
   ss << removeBeginEndBlanks(line);
   ss >> s;

   myBond.number = gTopo->bonds.size();
   
   myBond.atom1 = toInt(s);

   ss >> s;

   myBond.atom2 = toInt(s);

   ss >> s;

   myBond.funct = toInt(s);

   gTopo->bonds.push_back(myBond);

   report << debug(2) <<"Bond information in GromacsTopology::Bond : "<<myBond.number<<" "<<myBond.atom1<<" "<<myBond.atom2<<" "<<myBond.funct<<endr;
  
}

void GromacsTopologyReader::FillGromacsTopologyAngle(string line) {

  GromacsTopology::Angle myAngle;

   string s;
   stringstream ss;

   //report << plain << line <<endr;

   ss << removeBeginEndBlanks(line);

   ss >> s;

   myAngle.number = gTopo->angles.size();

   myAngle.atom1 = toInt(s);

   ss >> s;

   myAngle.atom2 = toInt(s);

   ss >> s;

   myAngle.atom3 = toInt(s);

   ss >> s;

   myAngle.funct = toInt(s);

   gTopo->angles.push_back(myAngle);

   report << debug(2) <<"Angle information in GromacsTopology::Angle : "<<myAngle.number<<" "<<myAngle.atom1<<" "<<myAngle.atom2<<" "<<myAngle.atom3<<" "<<myAngle.funct<<endr;  

}

void GromacsTopologyReader::FillGromacsTopologyDihedrals(string line) {

  string s1,s2,s3,s4,s5;
  stringstream ss;

  //report << plain << line <<endr;

  ss << removeBeginEndBlanks(line);

  ss >> s1;
  ss >> s2;
  ss >> s3;
  ss >> s4;
  ss >> s5;

  int rb_dih_num = gTopo->rb_dihedrals.size();

  int dih_num = gTopo->dihedrals.size();

  if (toInt(s5) == 3) {
    //Ryckert-Belleman
    GromacsTopology::RBDihedral g_rb_dih;
    g_rb_dih.number = rb_dih_num;
    g_rb_dih.atom1 = toInt(s1);
    g_rb_dih.atom2 = toInt(s2);
    g_rb_dih.atom3 = toInt(s3);
    g_rb_dih.atom4 = toInt(s4);
    g_rb_dih.funct = toInt(s5);
  
    gTopo->rb_dihedrals.push_back(g_rb_dih);

    report << debug(2) <<"Improper information in GromacsTopology::Improper "<<g_rb_dih.number<<" "<<g_rb_dih.atom1<<" "<<g_rb_dih.atom2<<" "<<g_rb_dih.atom3<<" "<<g_rb_dih.atom4<<" "<<g_rb_dih.funct<<endr;
  }else if (toInt(s5) == 1){
    //a dihedral
    GromacsTopology::Proper_Dihedral g_dih;
    g_dih.number = dih_num;
    g_dih.atom1 = toInt(s1);
    g_dih.atom2 = toInt(s2);
    g_dih.atom3 = toInt(s3);
    g_dih.atom4 = toInt(s4);
    g_dih.funct = toInt(s5);

    gTopo->dihedrals.push_back(g_dih);

    report << debug(2) <<"Dihedral information in GromacsTopology::Dihedral "<<g_dih.number<<" "<<g_dih.atom1<<" "<<g_dih.atom2<<" "<<g_dih.atom3<<" "<<g_dih.atom4<<" "<<g_dih.funct<<endr;

  }

}

string GromacsTopologyReader::GetParameterFilename() {

  //find the include file which starts with ff
  string fname_start("ff");
  
  for(unsigned int ii=0;ii<parameter_files.size();ii++) {
     if (equalStart(fname_start,parameter_files[ii])) return parameter_files[ii];
  }

  //if control reaches here then it means parameter file is not found.
  //return an empty string
  string s("");
  return s;
} 
