#include <protomol/type/PDB.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____PDB
PDB::Atom::Atom() :
  elementType(""), elementNum(0), elementName(""), altLoc(""), residueName(""),
  chainID(""), residueNum(0), insertionCode(""), occupancy(0.0),
  tempFactor(0.0), segID(""), symbol(""), charge(""), hvyAtomGrpsize(0) {}

PDB::Atom::Atom(string etype, int anum, string ename, string altloc,
                string rname, string chain, int rnum, string insertion,
                Real occ, Real tf, string segname, string symname, string c,
                int ha) :
  elementType(etype), elementNum(anum), elementName(ename), altLoc(altloc),
  residueName(rname), chainID(chain), residueNum(rnum),
  insertionCode(insertion), occupancy(occ), tempFactor(tf), segID(segname),
  symbol(symname), charge(c), hvyAtomGrpsize(ha) {}

PDB::Ter::Ter() :
  elementType(""), elementNum(0), residueName(""), chainID(""), residueNum(0),
  insertionCode("") {}

PDB::Ter::Ter(string etype, int anum, string rname, string chain, int rnum,
              string insertion) :
  elementType(etype), elementNum(anum), residueName(rname), chainID(chain),
  residueNum(rnum), insertionCode(insertion) {}

void PDB::clear() {
  coords.clear();
  atoms.clear();
  ters.clear();
}

//____ Atom
//____
//____ COLUMNS   DATA TYPE     FIELD       DEFINITION
//____ ------------------------------------------------------------------------
//____  1 -  6   Record name   "ATOM  "
//____  7 - 11   Integer       serial      Atom serial number.
//____ 13 - 16   Atom          name        Atom name.
//____ 17 - 17   Character     altLoc      Alternate location indicator.
//____ 18 - 20   Residue name  resName     Residue name.
//____ 22 - 22   Character     chainID     Chain identifier.
//____ 23 - 26   Integer       resSeq      Residue sequence number.
//____ 27 - 27   AChar         iCode       Code for insertion of residues.
//____ 31 - 38   Real(8.3)     x           Orthogonal coordinates for X in
//____                                     Angstroms.
//____ 39 - 46   Real(8.3)     y           Orthogonal coordinates for Y in
//____                                     Angstroms.
//____ 47 - 54   Real(8.3)     z           Orthogonal coordinates for Z in
//____                                     Angstroms.
//____ 55 - 60   Real(6.2)     occupancy   Occupancy.
//____ 61 - 66   Real(6.2)     tempFactor  Temperature factor.
//____ 73 - 76   LString(4)    segID       Segment identifier, left-justified.
//____ 77 - 78   LString(2)    element     Element symbol, right-justified.
//____ 79 - 80   LString(2)    charge      Charge on the atom.
//____ NOTE: The PDB says the length of the residue name is only 3 characters
//____  whereas XPLOR allows 4 character names.  We choose 4 for compatability
//____  with both systems (since we never change the length, we you give us is
//____  what we use)

namespace ProtoMol {
  MyStreamer &operator<<(MyStreamer &OS, const PDB::Atom &p) {
    OS << p.elementType << "," << p.elementNum << "," << p.elementName <<
      "," << p.altLoc << "," << p.residueName << "," << p.chainID << "," <<
      p.residueNum << "," << p.insertionCode << "," << p.occupancy << "," <<
      p.tempFactor << "," << p.segID << "," << p.symbol << "," << p.charge;
    return OS;
  }
}

