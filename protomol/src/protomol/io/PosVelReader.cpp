#include <protomol/io/PosVelReader.h>

#include <protomol/type/Vector.h>
#include <protomol/base/SystemUtilities.h>

#include <protomol/io/XYZBinReader.h>
#include <protomol/io/XYZReader.h>
#include <protomol/io/PDBReader.h>

using namespace std;
using namespace ProtoMol;
//____ File

PosVelReader::PosVelReader() :
  filename(""), myOk(true), myType(PosVelReaderType::UNDEFINED) {}

PosVelReader::PosVelReader(const string &filename) :
  filename(filename), myOk(SystemUtilities::isAccessible(filename)),
  myType(PosVelReaderType::UNDEFINED) {}

PosVelReader::operator void*() const {
  return !myOk ? 0 : const_cast<PosVelReader *>(this);
}

bool PosVelReader::operator!() const {
  return !myOk;
}

void PosVelReader::setFilename(const string &filename) {
  this->filename = filename;
  myOk = true;
  myType = PosVelReaderType::UNDEFINED;
}

bool PosVelReader::open() {
  myOk = SystemUtilities::isAccessible(filename);
  return myOk;
}

bool PosVelReader::open(const string &filename) {
  setFilename(filename);
  return open();
}

bool PosVelReader::tryFormat(PosVelReaderType::Enum type) {
  if (type == PosVelReaderType::PDB) {
    PDBReader reader(filename);
    return reader.tryFormat();
  } else if (type == PosVelReaderType::XYZ) {
    XYZReader reader(filename);
    return reader.tryFormat();
  } else if (type == PosVelReaderType::XYZBIN) {
    XYZBinReader reader(filename);
    return reader.tryFormat();
  } else
    return SystemUtilities::isAccessible(filename);
}

PosVelReaderType PosVelReader::getType() const {
  return myType;
}

namespace ProtoMol {
  PosVelReader &operator>>(PosVelReader &posReader, PDB &pdb) {
    posReader.myType = PosVelReaderType::UNDEFINED;

    // PDB
    PDBReader pdbReader(posReader.filename);
    posReader.myOk = pdbReader.tryFormat();
    if (posReader.myOk) {
      posReader.myOk = (pdbReader >> pdb ? true : false);
      posReader.myType = PosVelReaderType::PDB;
    }

    return posReader;
  }

  PosVelReader &operator>>(PosVelReader &posReader, XYZ &xyz) {
    posReader.myType = PosVelReaderType::UNDEFINED;

    // XYZ
    XYZReader xyzReader(posReader.filename);
    posReader.myOk = xyzReader.tryFormat();
    if (posReader.myOk) {
      posReader.myOk = (xyzReader >> xyz ? true : false);
      posReader.myType = PosVelReaderType::XYZ;
    }

    // XYZ binary
    if (!posReader.myOk) {
      XYZBinReader xyzBinReader(posReader.filename);
      if (xyzBinReader.tryFormat()) {
        posReader.myOk = (xyzBinReader >> xyz ? true : false);
        posReader.myType = PosVelReaderType::XYZBIN;
      }
    }

    // PDB
    if (!posReader.myOk) {
      PDBReader pdbReader(posReader.filename);
      if (pdbReader.tryFormat()) {
        posReader.myOk = (pdbReader >> xyz ? true : false);
        posReader.myType = PosVelReaderType::PDB;
      }
    }

    return posReader;
  }

  PosVelReader &operator>>(PosVelReader &posReader, Vector3DBlock &coords) {
    posReader.myType = PosVelReaderType::UNDEFINED;
    
    // XYZ
    XYZReader xyzReader(posReader.filename);
    posReader.myOk = xyzReader.tryFormat();
    if (posReader.myOk) {
      posReader.myOk = (xyzReader >> coords ? true : false);
      posReader.myType = PosVelReaderType::XYZ;
    }
    
    // XYZ binary
    if (!posReader.myOk) {
      XYZBinReader xyzBinReader(posReader.filename);
      if (xyzBinReader.tryFormat()) {
        posReader.myOk = (xyzBinReader >> coords ? true : false);
        posReader.myType = PosVelReaderType::XYZBIN;
      }
    }
    
    // PDB
    if (!posReader.myOk) {
      PDBReader pdbReader(posReader.filename);
      if (pdbReader.tryFormat()) {
        posReader.myOk = (pdbReader >> coords ? true : false);
        posReader.myType = PosVelReaderType::PDB;
      }
    }
    
    return posReader;
  }
}

