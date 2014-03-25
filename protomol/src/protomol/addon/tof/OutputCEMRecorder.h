#ifndef _OUTPUT_CEM_RECORDER_H
#define _OUTPUT_CEM_RECORDER_H

#include <protomol/addon/tof/ToFRecorder.h>
#include <protomol/addon/tof/CEMDetector.h>
#include <protomol/output/Output.h>
#include <vector>
#include <string>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMolAddon;

namespace ProtoMolAddon {
  
  namespace ToF {

    class OutputCEMRecorder : public Output {
    private:
      ToFRecorder<CEMDetector> recorder;

      Real detector_pos;
      string output_filename;

    public:
      OutputCEMRecorder(const string &output_filename, double detector_pos);
      OutputCEMRecorder();
      ~OutputCEMRecorder();
  
      Output *doMake(const vector<Value> &values) const;
      void doInitialize();
      void doRun(int step); 
      void doFinalize(int step);
  
    public:
      static const string keyword;
      void getParameters(vector<Parameter> &parameter) const;
      string getIdNoAlias() const {return keyword;};
    
    };

  }
}

#endif
