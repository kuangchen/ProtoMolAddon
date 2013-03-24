#include <algorithm>
#include <protomol/addon/CollisionSchedule.h>


using namespace ProtoMolAddon;
using namespace std;

double Collision::GetTime() const {
  return time;
}

double Collision::GetAtomId() const {
  return atomId;
}

Collision::Collision(int atomId, double time) :
  atomId(atomId),
  time(time) 
{}

bool Collision::operator<(Collision other) const {
  return time < other.time;
}

void CollisionSchedule::Generate(double start, double end, int atom_count, gsl_rng *r, double freq)  
{
  for (int n=0; n<atom_count; n++) 
    for (double t=start; t<end; t += gsl_ran_exponential(r, 1.0/freq))
      push_back(Collision(n, t));

  sort(this->begin(), this->end());

  
}

