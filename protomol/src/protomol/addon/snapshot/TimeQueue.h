#ifndef _TIME_QUEUE_H
#define _TIME_QUEUE_H

#include <queue>
#include <iostream>

namespace ProtoMolAddon {
  namespace Snapshot {

    class TimeQueue: public std::queue<double> {
    private:
      size_t init_len;
      
    public:
      TimeQueue(): init_len(0) {};
      TimeQueue(size_t n, double t0, double dt);
      size_t get_init_len() const { return init_len; }
      friend std::istream& operator>> (std::istream &is, TimeQueue &tq);
    };
  }
}
    

    
    // class TimeQueue: public std::vector<double> {
    // private:
    //   typedef std::vector<double>::const_iterator iter;
    //   iter curr;
      
    // public:
    //   TimeQueue() {}
      
    //   TimeQueue(size_t n, double t0, double dt) {

    // 	for (size_t i=0; i<n; i++) push_back(t0+dt*i);
    // 	curr = begin();
    //   }

    //   void step_curr() { curr++; }

    //   const_iterator current() const { return curr; }
    //   size_t saved_size() const { return curr - begin(); }

    //   friend std::istream& operator>> (std::istream &is, TimeQueue &tq) {
    // 	double t0, dt;
    // 	size_t n;
	
    // 	is >> n >> t0 >> dt;
    // 	tq = TimeQueue(n, t0, dt);

    //   return is;
    // }
#endif

