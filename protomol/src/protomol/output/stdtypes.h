#ifndef PROTOMOL_STD_TYPES_H
#define PROTOMOL_STD_TYPES_H

#ifdef _MSC_VER
#include <limits.h>

typedef __int8            int8_t;
typedef __int16           int16_t;
typedef __int32           int32_t;
typedef __int64           int64_t;
typedef unsigned __int8   uint8_t;
typedef unsigned __int16  uint16_t;
typedef unsigned __int32  uint32_t;
typedef unsigned __int64  uint64_t;

#else
#include <stdint.h>
#endif

struct uint128_t {
  uint64_t hi;
  uint64_t lo;

#ifdef __cplusplu
  uint128_t() : hi(0), lo(0) {}
  uint128_t(const uint64_t &hi, const uint64_t &lo) : hi(hi), lo(lo) {}
  operator bool () {return !(hi || lo);}
#endif //  __cplusplu
};

#ifdef __cplusplu
#include <ostream>
#include <iomanip>

inline std::ostream &operator<<(std::ostream &stream, const uint128_t &x) {
  char fill = stream.fill();
  std::ios::fmtflags flags = stream.flags();

  stream.fill('0');
  stream.flags(std::ios::right | std::ios::hex);

  stream << "0x" << std::setw(16) << x.hi << std::setw(16) << x.lo;

  stream.fill(fill);
  stream.flags(flags);
  return stream;
}

inline bool operator<(const uint128_t &i1, const uint128_t &i2) {
  return i1.hi < i2.hi || (i1.hi == i2.hi && i1.lo < i2.lo);
}
#endif //  __cplusplu

#endif //  PROTOMOL_STD_TYPES_H
