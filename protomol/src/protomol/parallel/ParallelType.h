/*  -*- c++ -*-  */
#ifndef PARALLELTYPE_H
#define PARALLELTYPE_H

#include <protomol/type/AbstractEnumType.h>

namespace ProtoMol {
  //____ ParallelEnum

  /**
   * The types of parallel strategies.
   *
   */
  class ParallelEnum  {
  public:
    virtual ~ParallelEnum() {}

    enum Enum {
      FIRST = 0,       // Only internal purpose
      UNDEFINED = 0,   // Value returned when no string matches
      STATIC,          ///< Static load balancing
      DYNAMIC,         ///< Dynamic load balancing, similar to master-slave but
                       ///< with working master
      MASTERSLAVE,     ///< Master-slave
      LAST             // Only internal purpose
    };
    static const std::string str[];
  };

  //____ ParallelType

  typedef AbstractEnumType<ParallelEnum> ParallelType;
}
#endif //  PARALLELTYPE_H
