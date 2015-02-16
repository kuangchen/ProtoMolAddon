#ifndef _OUTPUT_SNAPSHOT_H
#define _OUTPUT_SNAPSHOT_H

#include <protomol/ProtoMolApp.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>

#include <protomol/output/Output.h>
#include <protomol/addon/snapshot/Snapshot.h>
#include <protomol/addon/Constants.h>
#include <string>

namespace ProtoMolAddon {
  namespace Snapshot {

    using namespace ProtoMol;
    
    template<class Storage>
    class OutputSnapshot : public Output {
    private:
      std::string fname;
      Snapshot<Storage> ss;

    public:
      OutputSnapshot(const std::string& fname) :
	Output(1),
	fname(fname),
	ss(fname) {
      }

      OutputSnapshot() : Output(-1) {}

      ~OutputSnapshot() {}
  
      Output *doMake(const vector<Value> &values) const {
	return new OutputSnapshot(values[0]);
      }


      void doInitialize() {
	ss.Initialize(app);
      }

      void doRun(int step) {
	ss.Run( app->topology->time * Constant::ToSI::time);
      }

      void doFinalize(int step) {
	ss.Finalize();
      }
  
    public:
      void getParameters(vector<Parameter> &parameter) const {
	parameter.push_back(Parameter(getId(), Value(fname, ConstraintValueType::NotEmpty())));
      }

      std::string getIdNoAlias() const {return "Snapshot" + Storage::GetName();};
    };
  }
}

#endif
