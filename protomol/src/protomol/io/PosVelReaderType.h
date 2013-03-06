/* -*- c++ -*- */
#ifndef POS_VEL_READERTYPE_H
#define POS_VEL_READERTYPE_H

#include <protomol/type/AbstractEnumType.h>

namespace ProtoMol {
  //_____________________________________________________ PosVelReaderEnum

  class PosVelReaderEnum {
  public:
    virtual ~PosVelReaderEnum() {}
    enum Enum {
      FIRST = 0,       // Used internally only
      UNDEFINED = 0,  // PosVelReader returned when no string matches
      PDB,
      XYZ,
      XYZBIN,
      LAST              // Used internally only
    };

  protected:
    static const std::string str[];
  };

  //_____________________________________________________ PosVelReaderType
  typedef AbstractEnumType<PosVelReaderEnum> PosVelReaderType;
}
#endif /* POS_VEL_READERTYPE_H */
