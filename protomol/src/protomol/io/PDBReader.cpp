#include <protomol/io/PDBReader.h>

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>

//____#define DEBUG_PDB

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____PDBReader

PDBReader::PDBReader() :
  Reader(), myCoords(NULL), myAtoms(NULL), myTers(NULL) {}

PDBReader::PDBReader(const string &filename) :
  Reader(filename), myCoords(NULL), myAtoms(NULL), myTers(NULL) {}

PDBReader::~PDBReader() {
  if (myCoords != NULL)
    delete myCoords;
  if (myAtoms != NULL)
    delete myAtoms;
  if (myTers != NULL)
    delete myTers;
}

bool PDBReader::tryFormat() {
  if (!open())
    return false;
  do {
    string record, str;
    record = getline();
    stringstream ss(record);
    ss >> str;
    if ("ATOM" == str) {
      close();
      return true;
    } else if ("HETATM" == str) {
      close();
      return true;
    }
  } while (!file.eof());

  file.setstate(ios::failbit);
  close();
  return false;
}

bool PDBReader::read() {
  if (myCoords == NULL)
    myCoords = new Vector3DBlock();
  if (myAtoms == NULL)
    myAtoms = new vector<PDB::Atom>();
  if (myTers == NULL)
    myTers = new vector<PDB::Ter>();
  return read(*myCoords, *myAtoms, *myTers);
}

bool PDBReader::read(PDB &pdb) {
  return read(pdb.coords, pdb.atoms, pdb.ters);
}

bool PDBReader::read(Vector3DBlock &coords, vector<PDB::Atom> &atoms,
                     vector<PDB::Ter> &ters) {
  if (!tryFormat())
    return false;
  if (!open())
    return false;
  coords.clear();
  atoms.clear();
  ters.clear();
  std::vector<Vector3D> tmpdata;
  // Now we want to read data in until the record name is "END", then stop.
  int big = 0;
  int toBig = 0;
  comment = "";
  do {
    string record(getline());
    if (equalStart("END", record))
      break;
    if (equalStart("REMARK", record) || record.empty()) {
      record = removeBeginEndBlanks(record);
      if (!record.empty())
        comment += (comment.empty() ? "" : "\n") + record;
      continue;
    }
    if (equalStart("ATOM", record) || equalStart("HETATM", record)) {
      record = getRightFill(record, 80);
      string str =
        removeBeginEndBlanks(record.substr(PDB::Atom::S_RES_SEQ,
            PDB::Atom::L_RES_SEQ));
      int resSeq = toInt(str);
      if (!isInt(str)) {
        if (str.size() == 4 &&
            isInt(string(&str[1],
                &str[4])) && str[0] >= 'A' && str[0] <= 'Z') {
          resSeq =
            (str[0] - 'A') * 1000 + 10000 + toInt(string(&str[1], &str[4]));
          ++big;
        } else {
          resSeq = -1;
          ++toBig;
        }
      }

      atoms.push_back(PDB::Atom(removeBeginEndBlanks(record.substr(PDB::Atom::
                S_RECORD_NAME, PDB::Atom::L_RECORD_NAME)),
          toInt(record.substr(PDB::Atom::S_SERIAL, PDB::Atom::L_SERIAL)),
          removeBeginEndBlanks(record.substr(PDB::Atom::S_ATOM_NAME,
              PDB::Atom::L_ATOM_NAME)),
          removeBeginEndBlanks(record.substr(PDB::Atom::S_ALT_LOC,
              PDB::Atom::L_ALT_LOC)),
          removeBeginEndBlanks(record.substr(PDB::Atom::S_RES_NAME,
              PDB::Atom::L_RES_NAME + 1)),   // enable to read TIP3
          removeBeginEndBlanks(record.substr(PDB::Atom::S_CHAIN_ID,
              PDB::Atom::L_CHAIN_ID)),
          resSeq,
          removeBeginEndBlanks(record.substr(PDB::Atom::S_I_CODE,
              PDB::Atom::L_I_CODE)),
          toReal(record.substr(PDB::Atom::S_OCCUP, PDB::Atom::L_OCCUP)),
          toReal(record.substr(PDB::Atom::S_TEMP_FACT,
                               PDB::Atom::L_TEMP_FACT)),
          removeBeginEndBlanks(record.substr(PDB::Atom::S_SEG_ID,
              PDB::Atom::L_SEG_ID)),
          removeBeginEndBlanks(record.substr(PDB::Atom::S_ELEMENT_SYMBOL,
              PDB::Atom::L_ELEMENT_SYMBOL)),
          removeBeginEndBlanks(record.substr(PDB::Atom::S_CHARGE,
              PDB::Atom::L_CHARGE)),
          0));
      tmpdata.push_back(Vector3D(toReal(record.substr(PDB::Atom::S_X,
              PDB::Atom::L_X)),
          toReal(record.substr(PDB::Atom::S_Y, PDB::Atom::L_Y)),
          toReal(record.substr(PDB::Atom::S_Z, PDB::Atom::L_Z))));
    } else if (equalStart("TER", record)) {
      record = getRightFill(record, 27);
      string str =
        removeBeginEndBlanks(record.substr(PDB::Ter::S_RES_SEQ,
            PDB::Ter::L_RES_SEQ));
      int resSeq = toInt(str);
      if (!isInt(str)) {
        if (str.size() == 4 &&
            isInt(string(&str[1],
                &str[4])) && str[0] >= 'A' && str[0] <= 'Z') {
          resSeq =
            (str[0] - 'A') * 1000 + 10000 + toInt(string(&str[1], &str[4]));
          ++big;
        } else {
          resSeq = -1;
          ++toBig;
        }
      }

      ters.push_back(PDB::Ter(removeBeginEndBlanks(record.substr(PDB::Ter::
                S_RECORD_NAME, PDB::Ter::L_RECORD_NAME)),
          toInt(record.substr(PDB::Ter::S_SERIAL, PDB::Ter::L_SERIAL)),
          removeBeginEndBlanks(record.substr(PDB::Ter::S_RES_NAME,
              PDB::Ter::L_RES_NAME + 1)),  // enable to read TIP3
          removeBeginEndBlanks(record.substr(PDB::Ter::S_CHAIN_ID,
              PDB::Ter::L_CHAIN_ID)),
          resSeq,
          removeBeginEndBlanks(record.substr(PDB::Ter::S_I_CODE,
              PDB::Ter::L_I_CODE))));
    } else
      report << recoverable << "[PDB::read] Record unknown:\'" <<
      removeBeginEndBlanks(record) << "\'." << endr;
  } while (!(file.eof()));

  coords.resize(tmpdata.size());
  for (unsigned int i = 0; i < tmpdata.size(); i++)
    coords[i] = tmpdata[i];
  //for (unsigned int i = 0; i < coords.size(); i++)
  //  if (coords[i].mode) cout << coords[i].mode << endl;
  if (big > 0)
    report << hint << "[PDB::read] Found " << big <<
    " X-Plor residue number(s) starting with a character." << endr;
  if (toBig > 0)
    report << recoverable << "[PDB::read] Found " << toBig <<
    " non interger/X-Plor residue number(s)." << endr;
  close();
  return !file.fail();
}

Vector3DBlock *PDBReader::orphanCoords() {
  Vector3DBlock *tmp = myCoords;
  myCoords = NULL;
  return tmp;
}

vector<PDB::Atom> *PDBReader::orphanAtoms() {
  vector<PDB::Atom> *tmp = myAtoms;
  myAtoms = NULL;
  return tmp;
}

vector<PDB::Ter> *PDBReader::orphanTers() {
  vector<PDB::Ter> *tmp = myTers;
  myTers = NULL;
  return tmp;
}

PDB PDBReader::getPDB() const {
  PDB res;
  if (myCoords != NULL)
    res.coords = (*myCoords);
  if (myAtoms != NULL)
    res.atoms = (*myAtoms);
  if (myTers != NULL)
    res.ters = (*myTers);
  return res;
}

namespace ProtoMol {
  PDBReader &operator>>(PDBReader &pdbReader, PDB &pdb) {
    pdbReader.read(pdb.coords, pdb.atoms, pdb.ters);
    return pdbReader;
  }
  
  PDBReader &operator>>(PDBReader &pdbReader, Vector3DBlock &coords) {
    if (pdbReader.myAtoms == NULL)
      pdbReader.myAtoms = new vector<PDB::Atom>();
    if (pdbReader.myTers == NULL)
      pdbReader.myTers = new vector<PDB::Ter>();
    pdbReader.read(coords, *pdbReader.myAtoms, *pdbReader.myTers);
    return pdbReader;
  }

  PDBReader &operator>>(PDBReader &pdbReader, XYZ &xyz) {
    if (pdbReader.myAtoms == NULL)
      pdbReader.myAtoms = new vector<PDB::Atom>();
    if (pdbReader.myTers == NULL)
      pdbReader.myTers = new vector<PDB::Ter>();
    if (pdbReader.read(xyz.coords, *pdbReader.myAtoms, *pdbReader.myTers)) {
      xyz.names.resize(xyz.coords.size());
      for (unsigned int i = 0; i < xyz.coords.size(); ++i)
        xyz.names[i] = (*pdbReader.myAtoms)[i].elementName;
    }
    return pdbReader;
  }
}
