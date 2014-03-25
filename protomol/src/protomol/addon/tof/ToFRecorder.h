#ifndef __TOF_RECORDER_H
#define __TOF_RECORDER_H

#include <protomol/ProtoMolApp.h>
#include <protomol/addon/tof/AtomRecord.h>
#include <protomol/addon/Constants.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>
#include <string>
#include <vector>
#include <memory>
#include <iterator>

using std::endl;
using std::shared_ptr;
using std::string;
using std::vector;
using std::ostream;
using std::cout;
using std::ostream_iterator;
using namespace ProtoMol;
using namespace ProtoMolAddon::Constant;

namespace ProtoMolAddon {
  
  namespace ToF {

    template <class Detector>
      class ToFRecorder {

    private:
      vector< shared_ptr<AtomRecord> > record_list;
      Detector detector;

    public:
      ToFRecorder(const typename Detector::init_type& det_init= typename Detector::init_type());

      void Initialize(const ProtoMolApp *app);
      void UpdateRecord(const ProtoMolApp *app);

      template <class D>
      friend ostream& operator<< (ostream &os, const ToFRecorder<D> &detector);

      friend class OutputToFDetector;
    };

    template<class Detector>
      ToFRecorder<Detector>::ToFRecorder(const typename Detector::init_type &det_init):
	record_list(), detector(det_init)
    {}

    template <class Detector>
      void ToFRecorder<Detector>::Initialize(const ProtoMolApp *app) {
      for (size_t i=0; i<app->positions.size(); i++) 
	record_list.push_back(shared_ptr<AtomRecord>(new AtomRecord(app->topology->atoms.at(i).name, 
								    app->positions[i], 
								    app->velocities[i])));
    }
  
    template <class Detector>
      void ToFRecorder<Detector>::UpdateRecord(const ProtoMolApp *app) {
      double now = app->topology->time * TIME_CONV;

      for (auto &r: record_list) 
	detector.Process(r.get(), now);
    }
  
    
    template <class Detector>
    ostream& operator<< (ostream &os, const ToFRecorder<Detector> &recorder) {
      for (auto &r: recorder.record_list) 
	os << *r;

      return os;
    }

  }
}

#endif
 
