#ifndef OUTPUT_ION_SNAPSHOT_H
#define OUTPUT_ION_SNAPSHOT_H

extern "C" {
#include <lua5.1/lua.h>
#include <lua5.1/lauxlib.h>
#include <lua5.1/lualib.h>
#include <hdf5.h>
}

#include <protomol/output/Output.h>
#include <protomol/output/LinearPaulTrap.h>
#include <protomol/output/LuaState.h>
#include <vector>
#include <string>

using namespace std;
using LinearPaulTrap::Lqt;
using namespace LuaState;

namespace ProtoMol {

  class OutputIonSnapshot : public Output {
  private:
    string filename;
    vector<hid_t> fileId;
    vector<int> timeMark;

    string outputDir;
    int atomCount;
    int nextSnapshot;
    Real timeMin;
    Real timeMax;
    Real timeStep;
    
    Real h;
    Lqt trap;
    LuaState::LuaState L;

    int GenerateSnapshot(int snapshotIndex);

  public:
    OutputIonSnapshot(const string& filename);
    OutputIonSnapshot();
    ~OutputIonSnapshot();
  
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

#endif
