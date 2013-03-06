#include <protomol/output/OutputIonSnapshot.h>
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

const string OutputIonSnapshot::keyword("IonSnapshot");

OutputIonSnapshot::OutputIonSnapshot(): Output(-1)
{}

OutputIonSnapshot::~OutputIonSnapshot() {}

OutputIonSnapshot::OutputIonSnapshot(const std::string& filename):
  Output(1),
  numAtom(0),
  nextSnapshot(0),
  outputDir(""),
  L(filename) {

  trap = Lqt(L);
  outputDir = L.get<string>("dir");  

  int status = mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (status)
    if (errno != EEXIST)
      throw (status);

  numFrame = L.get<int>("ss.num_frame");
  ssStart = L.get< vector<double> >("ss.start");

}

Output* OutputIonSnapshot::doMake(const std::vector<Value> &values) const
{
  return new OutputIonSnapshot(values[0]);
}

void OutputIonSnapshot::getParameters(std::vector<Parameter> &parameter) const
{
  parameter.push_back(Parameter(getId(), Value(filename, ConstraintValueType::NotEmpty())));
}

     
void OutputIonSnapshot::doInitialize()
{
  numAtom = app->topology->atoms.size();
  h = getOutputFreq() * app->integrator->getTimestep() / Constant::SI::TIME_FS;

  int fps = static_cast<int>(1.0/h+0.5);
  
  ostringstream s;
  s.unsetf(ios::floatfield);
  s.precision(5);

  for (int i=0; i<ssStart.size(); i++) {
    int n = static_cast<int> (ssStart[i]/h + 0.5);
    ssStart[i] = h*n;

    SnapshotConfig config;
    config.start = ssStart[i];
    config.fps = fps;
    config.numFrame = numFrame;
    config.numAtom = numAtom;

    s << scientific << outputDir << "/time_" << ssStart[i] << ".hd5";
    //ssList.push_back(Snapshot(s.str(), ssStart[i], numFrame, numAtom));
    ssList.push_back(Snapshot(s.str(), config));
    s.flush();
    s.str("");
  }

  for (int i=0; i<ssStart.size(); i++) 
    ssList[i].AddHeader();
    
}

void OutputIonSnapshot::doRun(int step)
{
  Real now = app->topology->time /Constant::SI::TIME_FS;
  Real conversion = 1e-10 * ProtoMol::Constant::SI::TIME_FS / ProtoMol::Constant::TIMEFACTOR;

  bool active = false;
  for (int i=0; i<ssList.size(); i++) {
    active |= ssList[i].IsActive(now);
    if (active) break;
  }

  if (active) {
    Real *data = new double[numAtom * 14];

    vector<double> totEnergy, secEnergy;

    for (int i=0; i<numAtom; i++) {
      trap.GetEnergy(app->positions[i]*1e-10, app->velocities[i]*conversion, now, totEnergy, secEnergy);

      for (int j=0; j<3; j++) {
	// Positions
	data[i*14+j] = app->positions[i][j] * 1e-10;

	// Velocities
	data[i*14+j+3] = app->velocities[i][j] * conversion;
      
	// Total Energy and Secular Energy
	data[i*14+j+6] = totEnergy[j];
	data[i*14+j+9] = secEnergy[j];
      
	// Average Total
	data[i*14+12] = (data[i*14+6] + data[i*14+7]  + data[i*14+8])/3;
	data[i*14+13] = (data[i*14+9] + data[i*14+10] + data[i*14+11])/3;
      }
    }
    
    for (int i=0; i<ssList.size(); i++)
      ssList[i].AddFrame(now, data);

    delete []data;
  }
}

void OutputIonSnapshot::doFinalize(int step) {
}

