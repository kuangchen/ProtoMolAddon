/*  -*- c++ -*-  */
#ifndef TIMERSTATISTIC_H
#define TIMERSTATISTIC_H

#include <string>

#include <protomol/base/Timer.h>
#include <protomol/base/Report.h>

namespace ProtoMol {
  //____________________________________________________________ TimerStatistic
  /**
   * Simple timer statistic array implemented with Timer class and static's.@n
   *
   * TimerStatistic::timer[TimerStatistic::RUN].start();@n
   * TimerStatistic::timer[TimerStatistic::RUN].stop();@n
   * TimerStatistic::timer[TimerStatistic::RUN].reset();@n
   */
  class TimerStatistic {
public:
    enum Enum {
      FIRST = 0,       // Only internal purpose
      WALL = 0,
      RUN,
      INTEGRATOR,
      FORCES,
      COMMUNICATION,
      IDLE,
      LAST             // Only internal purpose
    };
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
    //  TimerStatistic();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class TimerStatistic
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
    static void setParallel(bool b) {myIsParallel = b;}
    static bool isParallel(void)    {return myIsParallel;}
    friend Report::MyStreamer &operator<<(Report::MyStreamer &os,
                                          const TimerStatistic &);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
    static Timer timer[];
private:
    static bool myIsParallel;
  };

  Report::MyStreamer &operator<<(Report::MyStreamer &os,
                                 const TimerStatistic &);
}
#endif /* TIMERSTATISTIC_H */
