#include <protomol/io/PARReader.h>

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <algorithm>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____PARReader

PARReader::PARReader(PAR::CharmmTypeEnum charmmType) :
  Reader(), myPAR(NULL), myCharmmType(charmmType),
  myCharmmTypeDetected(PAR::UNDEFINED) {}

PARReader::PARReader(const string &filename,
                     PAR::CharmmTypeEnum charmmType) :
  Reader(filename), myPAR(NULL), myCharmmType(charmmType),
  myCharmmTypeDetected(PAR::UNDEFINED) {}

PARReader::PARReader(const char *filename, PAR::CharmmTypeEnum charmmType) :
  Reader(string(filename)), myPAR(NULL), myCharmmType(charmmType),
  myCharmmTypeDetected(PAR::UNDEFINED) {}

PARReader::~PARReader() {
  if (myPAR != NULL)
    delete myPAR;
}

bool PARReader::openWith(const string &filename,
                         PAR::CharmmTypeEnum charmmType) {
  myCharmmType = charmmType;
  myCharmmTypeDetected = PAR::UNDEFINED;
  return open(filename);
}

bool PARReader::openWith(PAR::CharmmTypeEnum charmmType) {
  myCharmmType = charmmType;
  myCharmmTypeDetected = PAR::UNDEFINED;
  return open();
}

PAR *PARReader::orphanPAR() {
  PAR *tmp = myPAR;
  myPAR = NULL;
  return tmp;
}

bool PARReader::tryFormat() {
  if (!open())
    return false;
  while (!file.eof()) {
    string str;
    str = getline();
    stringstream ss(str);
    str.resize(find(str.begin(), str.end(), '!') - str.begin());
    str.resize(find(str.begin(), str.end(), '*') - str.begin());
    str.resize(find(str.begin(), str.end(), '{') - str.begin());
    str.resize(find(str.begin(), str.end(), '}') - str.begin());
    string i;
    ss >> i;
    if (!i.empty() &&
        (isKeywordCharmm19(i) ||
         isKeywordCharmm28(i)) && !equalStartNocase("END", i)) {
      close();
      return true;
    }
  }

  file.setstate(ios::failbit);
  close();
  return false;
}

bool PARReader::read() {
  if (myPAR == NULL)
    myPAR = new PAR();
  return read(*myPAR);
}

bool PARReader::read(PAR &par) {
  if (!tryFormat())
    return false;
  if (!open())
    return false;
  par.clear();

  vector<vector<string> > input;   // Stripped PAR file input
  vector<string> signatures;       // 'd' for Real or int, 'w' for a word
  vector<int> lines;               // line number in PAR file
  int charmm19 = 0;
  int charmm28 = 0;
  int count = 0;                   // PAR file line counter

  int comment = 0;
  while (!file.eof()) {
    vector<string> data;
    string numbers = "";
    string str(removeBeginEndBlanks(getline()));
    ++count;

    // Remove {} comments
    if (!str.empty() && str[0] == '*')
      continue;

    if (find(str.begin(), str.end(), '}') == str.end() &&
        find(str.begin(), str.end(), '{') == str.end()) {
      if (comment > 0)
        continue;
      else {
        string tmp = "";
        for (unsigned int i = 0; i < str.size(); ++i) {
          if (str[i] == '{')
            ++comment;
          else if (str[i] == '}') {
            --comment;
            tmp += " ";
            // Froce to close comments if } found at the end of a line
            if (comment > 0) {
              stringstream ss(str.substr(i, str.size() - i));
              string rest;
              ss >> rest;
              if (rest == "}" && ss.eof()) {
                comment = 0;
                report << warning << "The comments at line " << count <<
                " in PAR file are not properly closed." << endr;
              }
            }
            if (comment < 0) {
              report << warning << "The comments at line " << count <<
              " in PAR file has more closing \'}\'." << endr;
              comment = 0;
            }
          } else if (comment == 0)
            tmp += str[i];
        }

        str = tmp;
      }
    }

    // Remove ! comments
    str.resize(find(str.begin(), str.end(), '!') - str.begin());

    stringstream ss(str);
    while (!ss.eof()) {
      string i;
      ss >> i;
      // Skip * comments
      if (i[0] == '*')
        break;
      if (!i.empty()) {
        data.push_back(i);
        numbers += (isReal(i) ? "d" : "w");
      }
    }

    // Remove REMARK and SET
    if (!data.empty() &&
        (equalStartNocase("REMARK",
           data[0]) || equalStartNocase("SET", data[0]))) {
      data.clear();
      numbers = "";
    }

    // Store line
    if (!data.empty()) {
      input.push_back(data);
      signatures.push_back(numbers);
      lines.push_back(count);

      // Statistics for format selection
      if (numbers.size() == 1 && isKeywordCharmm28(data[0]))
        ++charmm28;
      else if ((numbers.size() > 1 && isKeywordCharmm19(data[0])) ||
               find(numbers.begin(), numbers.end(), 'w') == numbers.end())
        ++charmm19;
      else if (numbers.size() > 1 && !isKeywordCharmm19(data[0]) &&
               !isKeywordCharmm28(data[0]))
        ++charmm28;
    }
  }

  // Selection of format
  myCharmmTypeDetected = (charmm19 > charmm28 ? PAR::CHARMM19 : PAR::CHARMM28);
  bool isCharmm19 =
    (myCharmmType == PAR::UNDEFINED ? charmm19 >
     charmm28 : (myCharmmType == PAR::CHARMM19));

  // Parse
  PARRecordTypeEnum type = UNDEFINED;
  vector<string>::const_iterator sig = signatures.begin();
  vector<int>::const_iterator line = lines.begin();

  for (vector<vector<string> >::const_iterator i = input.begin();
       i != input.end();
       ++i, ++sig, ++line) {
    vector<string>::const_iterator j = i->begin();
    string s(*sig);
    if (isCharmm19) {   // Charmm19
      if (equalStartNocase("AEXP", (*j)) ||
          equalStartNocase("REXP", (*j)) ||
          equalStartNocase("HAEX", (*j)) ||
          equalStartNocase("AAEX", (*j)) ||
          equalStartNocase("NBOND", (*j)) ||
          equalStartNocase("CUTNB", (*j)) ||
          equalStartNocase("END", (*j)) ||
          equalStartNocase("CTONN", (*j)) ||
          equalStartNocase("EPS", (*j)) ||
          equalStartNocase("VSWI", (*j)) ||
          equalStartNocase("NBXM", (*j)) ||
          equalStartNocase("INHI", (*j)))
        continue;
      if (equalStartNocase("bond", (*j)))
        type = BOND;
      else if (equalStartNocase("angl", (*j)))
        type = ANGLE;
      else if (equalStartNocase("dihe", (*j)))
        type = DIHEDRAL;
      else if (equalStartNocase("impr", (*j)))
        type = IMPROPER;
      else if (equalStartNocase("nonb", (*j)))
        type = NONBONDED;
      else if (equalStartNocase("nbfi", (*j)))
        type = NBFIX;
      else if (equalStartNocase("hbon", (*j)))
        type = HBOND;
      else {
        report << warning << "Unknown parameter \'" << (*j) <<
        "\' in X-Plor parameter file at line " << (*line) << "." << endr;
        continue;
      }
      // Increment array pointer and remove first signature since the first
      // entry is a Charmm19 keyword
      ++j;
      s = s.substr(1, s.size() - 1);
    } else {   // Charmm28
      if (equalStartNocase("bond", (*j))) {
        type = BOND;
        continue;
      } else if (equalStartNocase("angl",
                   (*j)) || equalStartNocase("thet", (*j))) {
        type = ANGLE;
        continue;
      } else if (equalStartNocase("dihe",
                   (*j)) || equalStartNocase("phi", (*j))) {
        type = DIHEDRAL;
        continue;
      } else if (equalStartNocase("impr",
                   (*j)) || equalStartNocase("imph", (*j))) {
        type = IMPROPER;
        continue;
      } else if (equalStartNocase("nonb",
                   (*j)) || equalStartNocase("nbon", (*j))) {
        type = NONBONDED;
        continue;
      } else if (equalStartNocase("nbfi", (*j))) {
        type = NBFIX;
        continue;
      } else if (equalStartNocase("hbon", (*j))) {
        type = HBOND;
        continue;
      } else if (equalStartNocase("end",
                   (*j)) || equalStartNocase("nbxm", (*j)) ||
                 equalStartNocase("grou",
                   (*j)) || equalStartNocase("cdie", (*j)) ||
                 equalStartNocase("shif",
                   (*j)) || equalStartNocase("vgro", (*j)) ||
                 equalStartNocase("vdis",
                   (*j)) || equalStartNocase("vswi", (*j)) ||
                 equalStartNocase("cutn",
                   (*j)) || equalStartNocase("ctof", (*j)) ||
                 equalStartNocase("cton",
                   (*j)) || equalStartNocase("eps", (*j)) ||
                 equalStartNocase("e14f",
                   (*j)) || equalStartNocase("wmin", (*j)) ||
                 equalStartNocase("aexp",
                   (*j)) || equalStartNocase("rexp", (*j)) ||
                 equalStartNocase("haex",
                   (*j)) || equalStartNocase("aaex", (*j)) ||
                 equalStartNocase("noac",
                   (*j)) || equalStartNocase("hbno", (*j)) ||
                 equalStartNocase("cuth",
                   (*j)) || equalStartNocase("ctof", (*j)) ||
                 equalStartNocase("cton",
                   (*j)) || equalStartNocase("cuth", (*j)) ||
                 equalStartNocase("ctof",
                   (*j)) || equalStartNocase("cton", (*j))) {
        if (type == NONBONDED || type == NBFIX || type == HBOND)
          continue;
        report << warning << "Unknown parameter \'" << (*j) <<
        "\' in Charmm28 parameter file at line " << (*line) << "." << endr;
        continue;
      }
    }

    // The definition as one string
    string definition;
    for (vector<string>::const_iterator l = j; l != i->end(); ++l)
      definition += (definition.empty() ? "" : " ") + (*l);

    // Add the type to the corresponding container ...
    switch (type) {
    case BOND: {
      if (s == "wwdd")
        par.bonds.push_back(PAR::Bond(par.bonds.size() + 1, j[0], j[1],
            toReal(j[2]), toReal(j[3])));
      else
        report << warning << "Unknown bond definition \'" << definition <<
        "\' (" << s << ") in PAR file at line " << (*line) << "." << endr;
      break;
    }

    case ANGLE: {
      if (s == "wwwdd")
        par.angles.push_back(PAR::Angle(par.angles.size() + 1, j[0], j[1],
            j[2], toReal(j[3]), toReal(j[4]), false, 0.0, 0.0));
      else if (s == "wwwdddd")
        par.angles.push_back(PAR::Angle(par.angles.size() + 1, j[0], j[1],
                                        j[2], toReal(j[3]), toReal(j[4]),
                                        true, toReal(j[5]), toReal(j[6])));
      else if (s == "wwwddwdd")
        par.angles.push_back(PAR::Angle(par.angles.size() + 1, j[0], j[1],
                                        j[2], toReal(j[3]), toReal(j[4]), true,
                                        toReal(j[6]), toReal(j[7])));
      else
        report << warning << "Unknown angle definition \'" << definition <<
        "\' (" << s << ") in PAR file at line " << (*line) << "." << endr;
      break;
    }

    case DIHEDRAL: {
      if (s == "wwwwddd") {
        // Charmm19 with multiplicity == 1
        // Charmm28 with multiplicity >= 1
        string atom1(j[0]);
        string atom2(j[1]);
        string atom3(j[2]);
        string atom4(j[3]);
        int l = par.dihedrals.size() - 1;
        if (l >= 0 && !isCharmm19 &&
            ((par.dihedrals[l].atom1 == atom1 && par.dihedrals[l].atom2 ==
              atom2 &&
              par.dihedrals[l].atom3 == atom3 && par.dihedrals[l].atom4 ==
              atom4) ||
             (par.dihedrals[l].atom1 == atom4 && par.dihedrals[l].atom2 ==
              atom3 &&
              par.dihedrals[l].atom3 == atom2 && par.dihedrals[l].atom4 ==
              atom1))) {
          // Charmm28 and previous is the same, increase multiplicity
          par.dihedrals[l].multiplicity++;
          par.dihedrals[l].forceConstant.push_back(toReal(j[4]));
          par.dihedrals[l].periodicity.push_back(toInt(j[5]));
          par.dihedrals[l].phaseShift.push_back(toReal(j[6]));
        } else
          par.dihedrals.push_back(PAR::Dihedral(par.dihedrals.size() + 1,
              j[0], j[1], j[2], j[3], toReal(j[4]), toInt(j[5]), toReal(j[6])));
      } else if (s == "wwwwwdddd") {
        // Charmm19 with multiplicity > 1
        // Add the dihedral ...
        PAR::Dihedral dihedral(
          par.dihedrals.size() + 1, j[0], j[1], j[2], j[3], toReal(
            j[6]), toInt(j[7]), toReal(j[8]));
        par.dihedrals.push_back(dihedral);
        int multiplicity = toInt(j[5]);
        // ... and the multiplicity
        if (multiplicity > 1)
          for (int m = 1; m < multiplicity; ++m) {
            // Next line
            ++i;
            ++sig;
            ++line;
            if (i == input.end()) {
              report << warning <<
              "Uncomplete dihedral definition in PAR file at line " <<
              (*line) << "." << endr;
              break;
            }
            if (!equalStart("ddd", (*sig))) {
              report << warning <<
              "Wrong multiple dihedral definition in PAR file at line " <<
              (*line) << "." << endr;
              --i;
              --sig;
              --line;
              break;
            }
            j = i->begin();
            par.dihedrals[par.dihedrals.size() - 1].multiplicity++;
            par.dihedrals[par.dihedrals.size() - 1].forceConstant.push_back(
              toReal(j[0]));
            par.dihedrals[par.dihedrals.size() - 1].periodicity.push_back(
              toInt(j[1]));
            par.dihedrals[par.dihedrals.size() - 1].phaseShift.push_back(
              toReal(j[2]));
          }

      } else
        report << warning << "Unknown dihedral definition \'" <<
        definition << "\' (" << s << ") in PAR file at line " << (*line) <<
        "." << endr;
      break;
    }

    case IMPROPER: {
      if (s == "wwwwddd")
        par.impropers.push_back(PAR::Improper(par.impropers.size() + 1, j[0],
            j[1], j[2], j[3], toReal(j[4]), toInt(j[5]), toReal(j[6])));
      else
        report << warning << "Unknown improper definition \'" <<
        definition << "\' (" << s << ") in PAR file at line " << (*line) <<
        "." << endr;
      break;
    }

    case NONBONDED: {
      if ((!isCharmm19 && s != "wdddddd" &&
           s != "wddd") || (isCharmm19 && s != "wdddd")) {
        report << warning << "Unknown nonbonded definition \'" <<
        definition << "\' (" << s << ") in PAR file at line " << (*line) <<
        "." << endr;
        break;
      }
      PAR::Nonbonded nonbonded;

      nonbonded.number = par.nonbondeds.size() + 1;
      nonbonded.atom = j[0];
      nonbonded.polarizability = (isCharmm19 ? 0.0 : toReal(j[1]));
      nonbonded.epsilon = toReal(j[1 + (isCharmm19 ? 0 : 1)]);
      nonbonded.sigma =
        toReal(j[2 +
                 (isCharmm19 ? 0 : 1)]) *
        (isCharmm19 ? PAR::Nonbonded::SIGMA_CHARMM19_TO_CHARMM28 : 1.0);
      nonbonded.polarizability2 = (s.size() >= 7 ? toReal(j[4]) :
                                   0.0);
      nonbonded.epsilon14 =
        (s.size() >=
         5 ? toReal(j[3 + (isCharmm19 ? 0 : 2)]) :
         nonbonded.epsilon);
      nonbonded.sigma14 =
        (s.size() >=
         5 ? toReal(j[4 +
                      (isCharmm19 ? 0 : 2)]) :
         nonbonded.sigma) *
        (isCharmm19 ? PAR::Nonbonded::SIGMA_CHARMM19_TO_CHARMM28 : 1.0);

      nonbonded.negative = (nonbonded.epsilon < 0.0);
      nonbonded.vdw = (s.size() >= 5);
      nonbonded.negative2 = (nonbonded.epsilon14 < 0.0);
      par.nonbondeds.push_back(nonbonded);
      break;
    }

    case NBFIX: {
      if (s == "wwdddd")
        par.nbfixs.push_back(PAR::Nbfix(par.nbfixs.size() + 1, j[0], j[1],
            toReal(j[2]), toReal(j[3]), toReal(j[4]), toReal(j[5])));
      else if (s == "wwdd")
        par.nbfixs.push_back(PAR::Nbfix(par.nbfixs.size() + 1, j[0], j[1],
            toReal(j[2]), toReal(j[3]), toReal(j[2]), toReal(j[3])));
      else
        report << warning << "Unknown nbfix definition \'" << definition <<
        "\' (" << s << ") in PAR file at line " << (*line) << "." << endr;

      break;
    }

    case HBOND: {
      if (s == "wwdd")
        par.hbonds.push_back(PAR::Hbond(par.hbonds.size() + 1, j[0], j[1],
            toReal(j[2]), toReal(j[3])));
      else
        report << warning << "Unknown hbond definition \'" << definition <<
        "\' (" << s << ") in PAR file at line " << (*line) << "." << endr;
      break;
    }

    case UNDEFINED:
    default: {
      report << warning << "Unknown definition \'" << definition << "\' (" <<
      s << ") in PAR file at line " << (*line) << "." << endr;
      break;
    }
    }
  }

  close();
  return !file.fail();
}

//____isKeywordCharmm19
bool PARReader::isKeywordCharmm19(const string &word) {
  return equalStartNocase("bond", word) || equalStartNocase("angl", word) ||
         equalStartNocase("dihe", word) || equalStartNocase("impr", word) ||
         equalStartNocase("nonb", word) || equalStartNocase("nbfi", word) ||
         equalStartNocase("hbon", word);
}

//____isKeywordCharmm28
bool PARReader::isKeywordCharmm28(const string &word) {
  return equalStartNocase("BOND", word) || equalStartNocase("ANGL", word) ||
         equalStartNocase("THET", word) || equalStartNocase("DIHE", word) ||
         equalStartNocase("PHI", word) || equalStartNocase("IMPH", word) ||
         equalStartNocase("IMPROPER",
    word) || equalStartNocase("NBOND", word) ||
         equalStartNocase("NONBONDED", word) || equalStartNocase("NBFIX",
    word) ||
         equalStartNocase("HBOND", word) || equalStartNocase("END", word);
}

namespace ProtoMol {
  PARReader &operator>>(PARReader &parReader, PAR &par) {
    parReader.read(par);
    return parReader;
  }
}

