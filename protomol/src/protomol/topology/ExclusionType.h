/*  -*- c++ -*-  */
#ifndef EXCLUSIONTYPE_H
#define EXCLUSIONTYPE_H

#include <protomol/type/AbstractEnumType.h>

namespace ProtoMol {
  //________________________________________ ExclusionEnum
  /**
   * Exclusion types for intra-molecular interactions
   */
  class ExclusionEnum {
  public:
    virtual ~ExclusionEnum() {}
    enum Enum {
      FIRST = 0, // Used internally only
      UNDEFINED = 0,  // Value returned when no string matches
      NONE,         ///< no exclusions at all
      ONE2,         ///< exclude 1. neighbors
      ONE3,         ///< exclude 1. and 2. neighbors
      ONE4,         ///< exclude 1., 2. and 3. neighbors
      ONE4MODIFIED, ///< exclude 1. and 2. and modify 3. neighbors
      LAST // Used internally only
    };
    static const std::string str[];
  };


  //________________________________________ ExclusionType

  typedef AbstractEnumType<ExclusionEnum> ExclusionType;
}
#endif
