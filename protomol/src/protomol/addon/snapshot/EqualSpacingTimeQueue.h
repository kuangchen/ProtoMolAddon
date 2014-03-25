#ifndef _EQUAL_SPACING_TIME_QUEUE_H
#define _EQUAL_SPACING_TIME_QUEUE_H

#include <iosfwd>

using std::istream;

namespace ProtoMolAddon {
  namespace Snapshot {
  
    struct EqualSpacingTimeQueue {
      double t0;
      double dt;
      size_t n;
      size_t current;

      EqualSpacingTimeQueue(double t0=0, double dt=0, size_t n=0);
      double PopFront();
      bool IsDue(double time);

      friend istream& operator>> (istream &is, EqualSpacingTimeQueue &spec); 
    };
  }
}

#endif
