#ifndef OUTPUT_ION_SNAPSHOT_H
#define OUTPUT_ION_SNAPSHOT_H

#ifdef ARCHLINUX

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#else
extern "C" {
#include <lua5.1/lua.h>
#include <lua5.1/lauxlib.h>
#include <lua5.1/lualib.h>
}


#endif


#include <protomol/output/Output.h>
#include <protomol/addon/LinearPaulTrap.h>
#include <protomol/addon/LuaConfigReader.h>
#include <protomol/addon/Snapshot.h>
#include <vector>
#include <string>

using namespace std;
using namespace ProtoMolAddon::LinearPaulTrap;
using namespace ProtoMolAddon::Lua;
using namespace ProtoMolAddon::IonSnapshot;

namespace ProtoMol {

  class OutputIonSnapshot : public Output {
  private:
    string filename;
    vector<Snapshot> ssList;
    vector<double> ssStart;
    int numFrame;
    vector<int> timeMark;

    string outputDir;
    int numAtom;
    int nextSnapshot;
    
    Real h;
    Lqt trap;
    LuaConfigReader reader;

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
