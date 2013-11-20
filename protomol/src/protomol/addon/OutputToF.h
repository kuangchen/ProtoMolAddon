#ifndef OUTPUT_TOF_H
#define OUTPUT_TOF_H

extern "C" {
#include <lua5.1/lua.h>
#include <lua5.1/lauxlib.h>
#include <lua5.1/lualib.h>
}

#include <protomol/addon/ToFDetector.h>
#include <protomol/output/Output.h>
#include <protomol/addon/LuaConfigReader.h>
#include <vector>
#include <string>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMolAddon;
using namespace ProtoMolAddon::Lua;

namespace ProtoMolAddon {

  class OutputToF : public Output {
  private:
    LuaConfigReader reader;
    ToFDetector detector;
    string output_filename;

  public:
    OutputToF(const string& filename);
    OutputToF();
    ~OutputToF();
  
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
