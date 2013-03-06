
#include <protomol/ProtoMolApp.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/output/OutputIonTemperatureLS.h>
#include <protomol/base/PMConstants.h>
#include <iostream>
#include <iterator>
#include <sys/types.h>
#include <sys/stat.h>
#include <omp.h>
#include <errno.h>
#include <sstream>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Constant;

extern int errno;

const string OutputIonTemperatureLS::keyword("IonTemperatureLSFile");

Output * OutputIonTemperatureLS::doMake(const std::vector<Value> &values) const
{
  return new OutputIonTemperatureLS(values[0], values[1], values[2], values[3]);
}

OutputIonTemperatureLS::OutputIonTemperatureLS(): Output(-1)
{}

OutputIonTemperatureLS::OutputIonTemperatureLS(const std::string &output_dir, const std::string &input_filename, double average_time, unsigned int output_freq):
  Output(1), input_filename(input_filename), output_dir(output_dir), average_time(average_time), output_freq(output_freq)
{
  input_file.open(input_filename.c_str());
  std::copy( std::istream_iterator<double>(input_file), std::istream_iterator<double>(), std::back_inserter(output_time) );
  std::sort( output_time.begin(), output_time.end() );
  input_file.close();

  time_count = output_time.size();
  output_counter.resize(time_count);

  int status;
  status = mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (status)
    if ( errno != EEXIST )
      throw(status);
  
  output_file.resize(0);
  
  for ( int i=0; i<time_count; i++ )
    {
      ostringstream ss;
      ss << output_dir << "/time_" << output_time[i] << ".txt";
      output_file.push_back(new ofstream(ss.str().c_str()));
      ss.flush();
    }
}

void OutputIonTemperatureLS::getParameters(std::vector<Parameter> &parameter) const
{
  parameter.push_back(Parameter(getId(), Value(output_dir, ConstraintValueType::NotEmpty())));
  parameter.push_back(Parameter("-input_filename", Value(input_filename)));
  parameter.push_back(Parameter("-average_time", Value(average_time)));
  parameter.push_back(Parameter("-output_freq", Value(output_freq)));
}


void OutputIonTemperatureLS::doInitialize()
{
  step_count_per_average_time = int(average_time / (app->integrator->getTimestep() / Constant::SI::TIME_FS));

  for (int i=0; i<time_count; i++)
      (*output_file[i]) << " # Average Time = " << average_time
		     << ", Simulation Step Size = " << app->integrator->getTimestep() / Constant::SI::TIME_FS 
			<< ", Step per average time = " << step_count_per_average_time << "\n";

  atom_count = app -> topology -> atoms.size();
}

void OutputIonTemperatureLS::doRun(int step)
{
  double current_time = app->topology->time / Constant::SI::TIME_FS;
  for ( unsigned int time_index = 0; time_index < time_count; time_index++ )
    if ( current_time >= output_time[time_index] && current_time < output_time[time_index] + average_time )
	{
	  if ( ( output_counter[time_index]++ % output_freq ) == 0 )
	    {
	      (*output_file[time_index]) << "t = " << scientific << current_time << "\n";
	      for ( unsigned int atom_index = 0; atom_index < atom_count; atom_index++ )
		{
		  Vector3D position(app->positions[atom_index] * 1e-10);
		  Vector3D velocity(app->velocities[atom_index] * 1e-10 * ProtoMol::Constant::SI::TIME_FS / ProtoMol::Constant::TIMEFACTOR);
	      
		  (*output_file[time_index]) << position[0] << "\t" << position[1] << "\t" << position[2] << "\t"; 
		  (*output_file[time_index]) << velocity[0] << "\t" << velocity[1] << "\t" << velocity[2] << "\n"; 
		}
	      (*output_file[time_index]) << "-----------------------------\n";  
	    }
	}
    else
      output_counter[time_index] == 0;
}


void OutputIonTemperatureLS::doFinalize(int step)
{
  for ( unsigned int i=0; i<time_count; i++ )
    (*output_file[i]).close();
}

