/*  -*- c++ -*-  */
#ifndef POS_VEL_READER_H
#define POS_VEL_READER_H

#include <protomol/io/PosVelReaderType.h>
#include <protomol/type/XYZ.h>
#include <protomol/type/PDB.h>

namespace ProtoMol {
  class Configuration;

  //________________________________________________________ PosVelReader
  class PosVelReader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    PosVelReader();
    explicit PosVelReader(const std::string &filename);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class PosVelReader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void setFilename(const std::string &filename);

    bool open();
    bool open(const std::string &filename);
    bool tryFormat(PosVelReaderType::Enum type);

    operator void*() const;
    bool operator!() const;
    // enable expression testing

    PosVelReaderType getType() const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend PosVelReader &operator>>(PosVelReader &posReader, PDB &pdb);
    friend PosVelReader &operator>>(PosVelReader &posReader, XYZ &xyz);
    friend PosVelReader &operator>>(PosVelReader &posReader,
                                    Vector3DBlock &coords);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    std::string filename;
    bool myOk;
    PosVelReaderType myType;
  };
}
#endif
