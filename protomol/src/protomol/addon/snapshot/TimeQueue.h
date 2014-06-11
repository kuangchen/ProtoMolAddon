#ifndef _TIME_QUEUE_H
#define _TIME_QUEUE_H

#include <deque>
#include <iostream>
#include <stdexcept>

namespace ProtoMolAddon {
  namespace Snapshot {

    class TimeQueue {
    private:
      std::deque<double> q_unsaved;
      std::deque<double> q_saved;
      size_t current;
      
    public:
      TimeQueue() :
	q_unsaved(), q_saved(), current(0) {}

      TimeQueue(size_t n, double t0, double dt) :
	q_unsaved(),	
	current(0) {

	for (size_t i=0; i<n; i++) q_unsaved.push_back(t0+dt*i);
      }

      void Pop() { 
	q_saved.push_back(q_unsaved.front());
	q_unsaved.pop_front(); 
      }
      
      size_t Size() const { return q_unsaved.size() + q_saved.size(); }
      
      bool IsDue(double now) const { 
	if (q_unsaved.empty()) return false;
	else 
	  return q_unsaved.front() < now; 
      }
      
      const std::deque<double>& SavedTime() const { return q_saved; }

      friend std::istream& operator>> (std::istream &is, TimeQueue &tq) {
	double t0, dt;
	size_t n;
	
	is >> n >> t0 >> dt;
	tq = TimeQueue(n, t0, dt);

	return is;
      }      
    };
  }
}

#endif

