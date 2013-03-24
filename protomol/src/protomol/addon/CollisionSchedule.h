#ifndef _COLLISION_SCHEDULE_H
#define _COLLISION_SCHEDULE_H

#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

namespace ProtoMolAddon {

  class Collision { 
  private:
    int atomId; 
    double time; 
  
  public:
    Collision(int atomId=0, double time=0);

    double GetTime() const;
    double GetAtomId() const;
    bool operator<(Collision other) const;
  };

  class CollisionSchedule: public vector<Collision> {
  public:
    void Generate(double start, double end, int atom_count, gsl_rng *r, double freq);
  };

}
#endif
