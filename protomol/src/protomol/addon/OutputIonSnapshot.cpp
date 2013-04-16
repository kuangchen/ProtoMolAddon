#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <sstream>
#include <protomol/ProtoMolApp.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <protomol/addon/OutputIonSnapshot.h>
#include <string>

extern int errno;
using namespace ProtoMol;
using namespace std;
using namespace ProtoMolAddon::Constant;

const string OutputIonSnapshot::keyword("IonSnapshot");

OutputIonSnapshot::OutputIonSnapshot(): Output(-1)
{}

OutputIonSnapshot::~OutputIonSnapshot() {}

OutputIonSnapshot::OutputIonSnapshot(const string& filename):
  Output(1),
  outputDir(""),
  numAtom(0),
  nextSnapshot(0),
  reader(filename) {

  trap = Lqt(reader);
  outputDir = reader.GetValue<string>("dir");  

  int status = mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (status)
    if (errno != EEXIST)
      throw (status);

  numFrame = reader.GetValue<int>("ss.num_frame");
  ssStart = reader.GetValue< vector<double> >("ss.start");

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
  
  int modifiedNumFrame = 0;
  if (numFrame<0) {
    vector<double> w = trap.GetSecularFrequency();
    vector<double>::iterator it = min_element(w.begin(), w.end());
    modifiedNumFrame = 2 * static_cast<int>(fps / (*it));
    numFrame = modifiedNumFrame;
  }  

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
    config.numFrameperMM = static_cast<int> (fps/trap.GetFrequency());

    s << scientific << outputDir << "/snapshot_" << i << ".hd5";
    ssList.push_back(Snapshot(s.str(), config));
    s.flush();
    s.str("");
  }

  for (vector<Snapshot>::iterator iter=ssList.begin(); iter!=ssList.end(); iter++) 
    iter->AddHeader();
    
}

void OutputIonSnapshot::doRun(int step)
{
  Real now = app->topology->time * TIME_CONV;
  bool active = false;

  for (Snapshot &s : ssList) {
    active |= s.IsActive(now);
    if (active) break;
  }

  if (active) {
    Real *data = new double[numAtom * 15];

    vector<double> totEnergy, secEnergy;

    Vector3D vel, pos;
    for (int i=0; i<numAtom; i++) {
      vel = app->velocities[i] * VELOCITY_CONV;
      pos = app->positions[i] * POSITION_CONV;

      trap.GetEnergy(pos, vel, now, totEnergy, secEnergy);
      double *offset = data+i*15;

      std::copy(pos.c, pos.c+3, offset);
      std::copy(vel.c, vel.c+3, offset+3);
      std::copy(totEnergy.begin(), totEnergy.end(), offset+6);
      std::copy(secEnergy.begin(), secEnergy.end(), offset+9);

      // for (int j=0; j<3; j++) {
      // 	// Positions
      // 	data[i*15+j] = pos[j];

      // 	// Velocities
      // 	data[i*15+j+3] = vel[j];
      
      // 	// Total Energy and Secular Energy
      // 	data[i*15+j+6] = totEnergy[j];
      // 	data[i*15+j+9] = secEnergy[j];
      // }      
	// Average Total
      data[i*15+12] = (data[i*15+6] + data[i*15+7]  + data[i*15+8])/3;
      data[i*15+13] = (data[i*15+9] + data[i*15+10] + data[i*15+11])/3;
      data[i*15+14] = app->topology->atoms[i].scaledCharge / Constant::SQRTCOULOMBCONSTANT;
    }
    
    for (Snapshot &s : ssList) 
	s.AddFrame(now, data);

    delete []data;
  }
}

void OutputIonSnapshot::doFinalize(int step) {
}

