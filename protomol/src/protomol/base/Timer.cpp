/*
 * Timer.C
 *
 * Copyright (C) 1994, Christoph Streit <streit@iam.unibe.ch>
 *                     Stephan Amann <amann@iam.unibe.ch>
 *                     University of Berne, Switzerland
 *
 * All rights reserved.
 *
 * This software may be freely copied, modified, and redistributed
 * provided that this copyright notice is preserved on all copies.
 *
 * You may not distribute this software, in whole or in part, as part of
 * any commercial product without the express consent of the authors.
 *
 * There is no warranty or other guarantee of fitness of this software
 * for any purpose.  It is provided solely "as is".
 *
 */

#ifdef SVR4
#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/param.h>
#elif defined (_WIN32)
#include <time.h>
#include <windows.h>  // mmsystem.h cannot exist without this.
#include <mmsystem.h> // timeGetTime()
#else
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <protomol/base/Timer.h>
#include <iostream>

using namespace ProtoMol::Report;
using namespace ProtoMol;
using namespace std;

//____ TimeRep

TimeRep::TimeRep() :
  myRealTime(0.0), myUserTime(0.0), mySystemTime(0.0)
{}

TimeRep::TimeRep(double realTime, double userTime, double sysTime) :
  myRealTime(realTime), myUserTime(userTime), mySystemTime(sysTime)
{}

TimeRep TimeRep::operator+(const TimeRep &time) const {
  return TimeRep(myRealTime + time.myRealTime,
    myUserTime + time.myUserTime,
    mySystemTime + time.mySystemTime);
}

TimeRep TimeRep::operator-(const TimeRep &time) const {
  return TimeRep(myRealTime - time.myRealTime,
    myUserTime - time.myUserTime,
    mySystemTime - time.mySystemTime);
}

TimeRep &TimeRep::operator+=(const TimeRep &time) {
  *this = *this + time;
  return *this;
}

TimeRep &TimeRep::operator-=(const TimeRep &time) {
  *this = *this - time;
  return *this;
}

MyStreamer &ProtoMol::operator<<(MyStreamer &os, const TimeRep &time) {
  os.setf(ios::showpoint | ios::fixed);
  os.precision(5);
  os << time.myRealTime << "[s] real, ";
  os.precision(5);
  os << time.myUserTime << "[s] user, ";
  os.precision(5);
  os << time.mySystemTime << "[s] sys";
  return os;
}

void TimeRep::set(double realTime, double userTime, double sysTime) {
  myRealTime = realTime;
  myUserTime = userTime;
  mySystemTime = sysTime;
}

void TimeRep::reset() {
  myRealTime = 0.0;
  myUserTime = 0.0;
  mySystemTime = 0.0;
}

//____ Timer

Timer::Timer() :
  myRunningFlag(false) {
  reset();
}

void Timer::start() {
  if (myRunningFlag)
    return;

  myRunningFlag = true;
  myStartTime = getCurrentTime();
  myLastLapTime = myStartTime;
  myLapTime.reset();
}

void Timer::stop() {
  if (!myRunningFlag)
    return;

  myTotalTime = myTotalTime + (getCurrentTime() - myStartTime);
  myRunningFlag = false;
}

TimeRep Timer::getActualTime() const {
  if (!myRunningFlag)
    return myTotalTime;
  return myTotalTime + (getCurrentTime() - myStartTime);
}

TimeRep Timer::lap() {
  if (!myRunningFlag)
    return TimeRep();

  TimeRep t = getCurrentTime();
  myLapTime = t - myLastLapTime;
  myLastLapTime = t;
  return myLapTime;
}

void Timer::reset() {
  myRunningFlag = false;

  myStartTime.reset();
  myTotalTime.reset();
  myLastLapTime.reset();
  myLapTime.reset();
}

TimeRep Timer::getTime() const {
  if (myRunningFlag)
    return myTotalTime + (getCurrentTime() - myStartTime);
  else
    return myTotalTime;
}

TimeRep Timer::getLapTime() const {
  if (!myRunningFlag)
    return TimeRep();

  return myLapTime;
}

Timer &Timer::operator+=(const TimeRep &time) {
  if (!myRunningFlag)
    myTotalTime += time;
  else
    myStartTime -= time;

  return *this;
}

Timer &Timer::operator-=(const TimeRep &time) {
  if (!myRunningFlag)
    myTotalTime -= time;
  else
    myStartTime += time;

  return *this;
}

MyStreamer &ProtoMol::operator<<(MyStreamer &os, const Timer &timer) {
  os << timer.getTime();
  return os;
}

TimeRep Timer::getCurrentTime() {
  double realTime, userTime, sysTime;

#ifdef SVR4
  struct timeval tv;

  gettimeofday(&tv, NULL);
  realTime = (double)tv.tv_sec + (double)tv.tv_usec / 1000000.;

  struct tms time;

  times(&time);
  userTime = (double)time.tms_utime / (double)HZ;
  sysTime = (double)time.tms_stime / (double)HZ;

#elif defined (WIN32)
#define EPOCHFILETIME (116444736000000000LL)
  FILETIME ft;
  LARGE_INTEGER li;
  __int64 t;
  GetSystemTimeAsFileTime(&ft);
  li.LowPart = ft.dwLowDateTime;
  li.HighPart = ft.dwHighDateTime;
  t = li.QuadPart;
  t -= EPOCHFILETIME;
  t /= 10;
  realTime =
    (double)((long)(t / 1000000)) + (double)((long)(t % 1000000)) * .000001;
  userTime = clock() / 1000.0;       // don't know how to measure user sys time
  sysTime = 0;                  // this is not a unix system.
#else
  struct timeval tv;

  gettimeofday(&tv, NULL);
  realTime = (double)tv.tv_sec + (double)tv.tv_usec * .000001;

  struct rusage usage;

  getrusage(RUSAGE_SELF, &usage);
  userTime = (double)usage.ru_utime.tv_sec +
             (double)usage.ru_utime.tv_usec * .000001;
  sysTime = (double)usage.ru_stime.tv_sec +
            (double)usage.ru_stime.tv_usec * .000001;
#endif

  return TimeRep(realTime, userTime, sysTime);
}

