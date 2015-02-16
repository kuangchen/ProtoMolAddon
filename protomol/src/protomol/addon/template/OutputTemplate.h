#ifndef _OUTPUT_TEMPLATE_H
#define _OUTPUT_TEMPLATE_H

#include <vector>
#include <string>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/module/MainModule.h>

namespace ProtoMolAddon {
  namespace Template {

    using namespace ProtoMol;
    
    template<class Recorder>
    class GenericOutput: public Output {

    private:
      std::string fname;
      Recorder recorder;
      
    private:
      virtual GenericOutput *doMake(const std::vector<Value> &values) const {
	return new GenericOutput(values[0]);
      }

    public:
      GenericOutput(): Output(-1) {}
      GenericOutput(const std::string &fname):
	Output(1),
	fname(fname),
	recorder(fname) {}
      
      ~GenericOutput() {} // ?????
            
      void doInitialize() { recorder.Initialize(app); }
      
      void doRun(int step) {
	recorder.Update(app->topology->time * ToSI::time);
      }
      
      void doFinalize(int step) { recorder.Finalize(); }

      void getParameters(std::vector<Parameter> &parameter) const {
	parameter.push_back(Parameter(getId(),
				      Value(fname, ConstraintValueType::NotEmpty())));
      }

    public:
      static const std::string keyword;
      
      virtual std::string getKeyword() const { return keyword; }
      virtual std::string getIdNoAlias() const { return keyword; }

    };
    
    template <class Recorder>
    const std::string GenericOutput<Recorder>::keyword(Recorder::GetName()); 
  }
}

#endif
