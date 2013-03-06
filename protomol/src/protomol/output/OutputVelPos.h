#ifndef OUTPUT_VEL_POS_H
#define OUTPUT_VEL_POS_H

extern "C" {
#include <lua5.1/lua.h>
#include <lua5.1/lauxlib.h>
#include <lua5.1/lualib.h>
}

#include <protomol/output/Output.h>

namespace ProtoMol {
  using namespace std;

  class OutputVelPos : public Output{
  private:
    typedef struct task{
      double start;
      unsigned int progress;
      ofstream* file;
    }task;
    
    string filename;

    vector<task> task_array;
    unsigned task_count;
    double duration;
    int step;
    int freq;
    string dir;

    double sim_step;
    unsigned int atom_count;
    unsigned int step_count_per_duration;

  private:
    int CreateOutputFile();
    void WriteHeader();
    int Close();
    bool TaskRunning(int i, double time);

  public:
    OutputVelPos(const string& filename);
    OutputVelPos();
    ~OutputVelPos();
    void SetAtomCount(int n);
    void SetSimStep(double s);
    
  
    Output *doMake(const vector<Value> &values) const;
    void doInitialize();
    void doRun(int step);
    void doFinalize(int step);
  
  public:
    const string keyword;
    void getParameters(vector<Parameter> &parameter) const;
    string getIdNoAlias() const {return keyword;};

  };
}

#endif
