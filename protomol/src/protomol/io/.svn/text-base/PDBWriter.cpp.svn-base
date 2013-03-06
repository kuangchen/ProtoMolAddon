#include <protomol/io/PDBWriter.h>

#include <protomol/type/PDB.h>
#include <protomol/base/Report.h>
#include <protomol/base/StringUtilities.h>

#include <iomanip>
#include <map>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____PDBWriter

PDBWriter::PDBWriter() :
  Writer() {}

PDBWriter::PDBWriter(const string &filename) :
  Writer(filename) {}

bool PDBWriter::write(const Vector3DBlock &coords, const PDB &pdb) {
  return write(coords, pdb.atoms, pdb.ters);
}

bool PDBWriter::write(const PDB &pdb) {
  return write(pdb.coords, pdb.atoms, pdb.ters);
}

bool PDBWriter::write(const Vector3DBlock &coords,
                      const vector<PDB::Atom> &atoms) {
  return write(coords, atoms, vector<PDB::Ter>());
}

bool PDBWriter::write(const Vector3DBlock &coords,
                      const vector<PDB::Atom> &atoms,
                      const vector<PDB::Ter> &ters) {
  if (!open())
    return false;

  const unsigned int count = atoms.size();
  if (coords.size() != count)
    report << error << "[PDBWriter::write]"
           << " Size of coordinates(" << coords.size() << ") and atoms("
           << atoms.size() << ") differ." << endr;

  map<int, int> tersMap;
  for (unsigned int i = 0; i < ters.size(); ++i)
    tersMap[ters[i].elementNum - 1] = i;

  int big = 0;
  int toBig = 0;
  string last = "";
  for (unsigned int i = 0; i < count; ++i) {
    string line(80, ' ');
    const PDB::Atom &a(atoms[i]);
    const Vector3D &c(coords[i]);
    string resSeq = toString(a.residueNum);
    if (a.residueNum >= 10000) {
      if (a.residueNum < 36000) {
        resSeq = string(1, 'A' + (a.residueNum - 10000) / 1000) + string(
          toString(a.residueNum % 1000));
        ++big;
      } else {
        ++toBig;
        resSeq = "-1";
      }
    }
    line.replace(PDB::Atom::S_RECORD_NAME, PDB::Atom::L_RECORD_NAME,
                 getRightFill(a.elementType,
                              PDB::Atom::L_RECORD_NAME));
    line.replace(PDB::Atom::S_SERIAL,
                 PDB::Atom::L_SERIAL,
                 getLeftFill(toString(a.elementNum), PDB::Atom::L_SERIAL));
    line.replace(PDB::Atom::S_ATOM_NAME,
                 PDB::Atom::L_ATOM_NAME,
                 getRightFill(a.elementName, PDB::Atom::L_ATOM_NAME));
    line.replace(PDB::Atom::S_ALT_LOC,
                 PDB::Atom::L_ALT_LOC,
                 getRightFill(a.altLoc, PDB::Atom::L_ALT_LOC));
    line.replace(PDB::Atom::S_RES_NAME,
                 PDB::Atom::L_RES_NAME,
                 // enable to write TIP3
                 getRightFill(a.residueName, PDB::Atom::L_RES_NAME + 1));   
    line.replace(PDB::Atom::S_CHAIN_ID,
                 PDB::Atom::L_CHAIN_ID,
                 getRightFill(a.chainID, PDB::Atom::L_CHAIN_ID));
    line.replace(PDB::Atom::S_RES_SEQ,
                 PDB::Atom::L_RES_SEQ,
                 getLeftFill(resSeq, PDB::Atom::L_RES_SEQ));
    line.replace(PDB::Atom::S_I_CODE,
                 PDB::Atom::L_I_CODE,
                 getRightFill(a.insertionCode, PDB::Atom::L_I_CODE));
    line.replace(PDB::Atom::S_X + 1,
                 PDB::Atom::L_X - 1,
                 getLeftFill(toString(c.c[0]), PDB::Atom::L_X - 1));
    line.replace(PDB::Atom::S_Y + 1,
                 PDB::Atom::L_Y - 1,
                 getLeftFill(toString(c.c[1]), PDB::Atom::L_Y - 1));
    line.replace(PDB::Atom::S_Z + 1,
                 PDB::Atom::L_Z - 1,
                 getLeftFill(toString(c.c[2]), PDB::Atom::L_Z - 1));
    line.replace(PDB::Atom::S_OCCUP,
                 PDB::Atom::L_OCCUP,
                 getLeftFill(toString(a.occupancy), PDB::Atom::L_OCCUP));
    line.replace(PDB::Atom::S_TEMP_FACT,
                 PDB::Atom::L_TEMP_FACT,
                 getLeftFill(toString(a.tempFactor), PDB::Atom::L_TEMP_FACT));
    line.replace(PDB::Atom::S_SEG_ID,
                 PDB::Atom::L_SEG_ID,
                 getLeftFill(a.segID, PDB::Atom::L_SEG_ID));
    line.replace(PDB::Atom::S_ELEMENT_SYMBOL,
                 PDB::Atom::L_ELEMENT_SYMBOL,
                 getLeftFill(a.symbol, PDB::Atom::L_ELEMENT_SYMBOL));
    line.replace(PDB::Atom::S_CHARGE,
                 PDB::Atom::L_CHARGE,
                 getLeftFill(a.charge, PDB::Atom::L_CHARGE));
    file << line << endl;
    if (tersMap.find(a.elementNum) != tersMap.end()) {
      string line(80, ' ');
      const PDB::Ter &t(ters[tersMap[a.elementNum]]);
      string resSeq = toString(t.residueNum);
      if (t.residueNum >= 10000) {
        if (t.residueNum < 36000) {
          resSeq = string(1, 'A' + (t.residueNum - 10000) / 1000) + string(
            toString(t.residueNum % 1000));
          ++big;
        } else {
          ++toBig;
          resSeq = "-1";
        }
      }
      line.replace(PDB::Ter::S_RECORD_NAME, PDB::Ter::L_RECORD_NAME,
                   getRightFill(t.elementType,
                                PDB::Ter::L_RECORD_NAME));
      line.replace(PDB::Ter::S_SERIAL,
                   PDB::Ter::L_SERIAL,
                   getLeftFill(toString(t.elementNum), PDB::Ter::L_SERIAL));
      line.replace(PDB::Ter::S_RES_NAME,
                   PDB::Ter::L_RES_NAME,
                   // enable to write TIP3
                   getRightFill(t.residueName, PDB::Ter::L_RES_NAME + 1));   
      line.replace(PDB::Ter::S_CHAIN_ID,
                   PDB::Ter::L_CHAIN_ID,
                   getRightFill(t.chainID, PDB::Ter::L_CHAIN_ID));
      line.replace(PDB::Ter::S_RES_SEQ,
                   PDB::Ter::L_RES_SEQ,
                   getLeftFill(resSeq, PDB::Ter::L_RES_SEQ));
      line.replace(PDB::Ter::S_I_CODE,
                   PDB::Ter::L_I_CODE,
                   getRightFill(t.insertionCode, PDB::Ter::L_I_CODE));
      file << line << endl;
    }
  }

  file << "END" << endl;
  if (big > 0)
    report << hint << "[PDB::write] Wrote " << big <<
    " X-Plor residue number(s) starting with a character." << endr;
  if (toBig > 0)
    report << recoverable << "[PDB::write] Wrote " << toBig <<
    " non interger/X-Plor residue number(s)." << endr;
  close();
  return !file.fail();
}

namespace ProtoMol {
  PDBWriter &operator<<(PDBWriter &pdbWriter, const PDB &pdb) {
    pdbWriter.write(pdb);
    return pdbWriter;
  }
}
