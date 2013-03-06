#include "OutputVelPos.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <sstream>
#include <protomol/ProtoMolApp.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/io/LuaConfigReader.h>

extern int errno;
using namespace ProtoMol;

using namespace std;
using namespace Util;

int OutputVelPos::CreateOutputFile()
{
  int status = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if ( status )
    if ( errno != EEXIST )
      throw (status);

  for ( int i = 0; i < task_count; i++ ){
    ostringstream ss;
    ss << dir << "/time_" << task_array[i].start << ".txt";

    task_array[i].file = new ofstream(ss.str().c_str());
    ss.flush();
  }
  
  return status;
}

OutputVelPos::OutputVelPos(): Output(-1)
{}

OutputVelPos::~OutputVelPos()
{}

OutputVelPos::OutputVelPos(const std::string& filename):Output(1)
{
  vector<double> v(0);

  std::cout << "opening " << filename << "\n";
  lua_State *L = luaL_newstate();
  luaL_openlibs(L);
    
  if (luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0))
    luaL_error(L, "cannot read file: %s", lua_tostring(L, -1));

  ExtractDouble(L, &duration, "duration",  "duration should be a number\n");
  ExtractInteger(L, &freq, "freq",    "freq should be a number\n");
  ExtractArray(L, v, "time_array", "time_array should be an array");
  ExtractString(L,dir, "dir", "dir should be a string");
  lua_close(L);

  task_count = v.size();
  task_array.resize(task_count);

  for ( unsigned int i = 0; i < task_count; i++)
    task_array[i].start = v[i];

  CreateOutputFile();  
}

void OutputVelPos::SetAtomCount(int n)
{
  atom_count = n;
}

void OutputVelPos::SetSimStep(double s)
{
  sim_step = s;
  step_count_per_duration = int(duration / s);
}

int OutputVelPos::Close()
{
  for ( int i = 0; i < task_count; i++ )
    task_array[i].file->close();
}

bool OutputVelPos::TaskRunning(int i, double time)
{
  return (time >= task_array[i].start) && (time < task_array[i].start + duration);
}

void OutputVelPos::WriteHeader()
{
  for ( int i = 0 ; i < task_count ; i++ )
    (*task_array[i].file) << " # Duration = " << duration 
			  << " , Sim step size = " << sim_step 
			  << " , Step counter per duration = " << step_count_per_duration << "\n";
    
}

Output* OutputVelPos::doMake(const std::vector<Value> &values) const
{
  return new OutputVelPos(values[0]);
}

void OutputVelPos::getParameters(std::vector<Parameter> &parameter) const
{
  parameter.push_back(Parameter(getId(), Value(filename, ConstraintValueType::NotEmpty())));
}

     
void OutputVelPos::doInitialize()
{
  SetAtomCount(app->topology->atoms.size());
  SetSimStep(app->integrator->getTimestep() / Constant::SI::TIME_FS);
  WriteHeader();
}

void OutputVelPos::doRun(int step)
{
  double current_time = app->topology->time / Constant::SI::TIME_FS;
  for ( unsigned int task_index = 0; task_index < task_count; task_index++ )
    if ( TaskRunning(task_index, current_time) ){
      if ( ( task_array[task_index].progress++ % freq ) == 0 ){
	(*task_array[task_index].file) << "t = "<< scientific << current_time << "\n";
	for ( unsigned int atom_index = 0; atom_index < atom_count ; atom_index++ ){
	  Vector3D pos(app->positions[atom_index] * 1e-10);
	  Vector3D vel(app->velocities[atom_index] * 1e-10 * ProtoMol::Constant::SI::TIME_FS / ProtoMol::Constant::TIMEFACTOR);
	  
	  (*task_array[task_index].file) << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\t";
	  (*task_array[task_index].file) << vel[0] << "\t" << vel[1] << "\t" << vel[2] << "\n";
	}
	(*task_array[task_index].file) << "-----------------------------\n";  
      }
    }
    
    else
      task_array[task_index].progress = 0;
}


void OutputVelPos::doFinalize(int step)
{
  Close();
}

