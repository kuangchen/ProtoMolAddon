#ifndef __BUFFER_GAS_H_
#define __BUFFER_GAS_H_

#include <protomol/addon/LuaState.h>
#include <protomol/type/Vector3D.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <string>
#include <protomol/ProtoMolApp.h>

using namespace ProtoMol;
using namespace std;
using namespace ProtoMolAddon::Lua;

namespace ProtoMolAddon {

  namespace Collision {
    struct Event {
      unsigned int atomID;
      double time;
      
    public:
      Event(unsigned int atomID=0, double time=0) : 
  	atomID(atomID), 
  	time(time) {}
      static bool Compare(const Event p1, const Event p2);
    };

    class BufferGas {
    private:
      double _mass;
      double _temperature;
      double _freq;
      double _vn;

      vector<Event> _schedule;
      int _nextEvent;

      const gsl_rng_type *_T;
      gsl_rng *_r;

    public:
      BufferGas();
      BufferGas(LuaState& L);
      ~BufferGas();

      
      void collide(double mass, Vector3D& pos, Vector3D& vel);
      vector<Event> createCollisionSchedule(double start, double end, int atomCount);

      bool isCollisionFinished() const;
      void collide(ProtoMolApp *app);
      void scheduleCollision(double start, double end, int numAtom);
      double getNextCollisionEventTime() const;
      
    };

  }
  
}


#endif
