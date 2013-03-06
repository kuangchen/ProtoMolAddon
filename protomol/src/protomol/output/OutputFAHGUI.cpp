//  Updated for standalone GUI
#if defined (HAVE_GUI) || defined (HAVE_LIBFAH)

#include <protomol/output/OutputFAHGUI.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/module/MainModule.h>
#include <protomol/ProtoMolApp.h>

#ifdef HAVE_LIBFAH
#include <fah/core/GUIServer.h>
#endif

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

#ifdef HAVE_LIBFAH
using namespace FAH;
#endif

#ifdef HAVE_LIBFAH
const string OutputFAHGUI::keyword("FAHGUI");
#else
const string OutputFAHGUI::keyword("Gui");
#endif


OutputFAHGUI::OutputFAHGUI() :
  name("ProtoMol"), port(0), portRange(0), timeout(0), pause(0),
  server(0) {}


OutputFAHGUI::OutputFAHGUI(const string &name, int freq, int port, int prange,
                           const string &projn, double timeout, bool pause) :
  Output(freq), name(name), port(port), portRange(prange),
  projName(projn), timeout(timeout), pause(pause), server(0) {

  // default pause timeout to 10 second
  if (pause && !timeout) timeout = 10;
}


void OutputFAHGUI::doInitialize() {
#ifdef HAVE_LIBFAH
  server = new GUIServer(projName.c_str(), app->topology->atoms.size(),
                         app->topology->bonds.size());
#else
  server = new GUIServer(projName.c_str(), app->topology->atoms.size(),
                         app->topology->bonds.size(), port, portRange);
  server->info.iterations = app->lastStep / 1000;
  server->info.frames = app->lastStep;
  server->current.frames_done = 0;
  // Mode number valid?
  if (app->eigenInfo.currentMode != -1)
    server->current.iterations_done = app->eigenInfo.currentMode;
  else server->current.iterations_done = 0;

  setAtoms();
  setBonds();
  // tart server now initial data set
  server->startServer();
#endif

  // timers for GuiTimeout
  guiTimer.reset();
  guiTimer.start();

}


void OutputFAHGUI::doRun(int step) {
  GUIServer::request_t request = server->getRequest();

  // GuiPause?
  if (pause && request == GUIServer::GS_NO_REQUEST){
    while ((request = server->getRequest()) == GUIServer::GS_NO_REQUEST){
      if ((guiTimer.getTime()).getRealTime() > timeout)
        THROW("GUI pause timed out.");
    }
  }

  switch (request) {
  case GUIServer::GS_META_REQUEST:
  case GUIServer::GS_COORD_REQUEST:
    guiTimer.reset(); // restart watchdog
    guiTimer.start();
    server->startUpdate();

    if (request == GUIServer::GS_META_REQUEST) {
      server->info.iterations = app->lastStep / 1000;
      server->info.frames = app->lastStep;
      setAtoms();
      setBonds();

    } else {
#ifdef HAVE_LIBFAH
      server->current.iterations_done = app->currentStep / 1000;
#endif
      server->current.frames_done = app->currentStep;
      server->current.energy = kineticEnergy(app->topology, &app->velocities);
      server->current.temperature =
        temperature(app->topology, &app->velocities);
      //  Mode number valid?
      if (app->eigenInfo.currentMode != -1)
        server->current.iterations_done = app->eigenInfo.currentMode;
      setCoords();
    }

    server->endUpdate();
    break;

  default: break;
  }

  if (timeout && (guiTimer.getTime()).getRealTime() > timeout)
    THROW("GUI communication timeout.");
}

void OutputFAHGUI::doFinalize(int step) {
  doRun(step);
  if (server) {
    delete server;
    server = 0;
  }
}


Output *OutputFAHGUI::doMake(const vector<Value> &values) const {
  return new OutputFAHGUI(values[0], values[1], values[2], values[3],
                          values[4], values[5], values[6]);
}


bool OutputFAHGUI::isIdDefined(const Configuration *config) const {
  return config->valid(getId());
}


void OutputFAHGUI::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(getId(), Value(name, ConstraintValueType::NotEmpty())));
  Output::getParameters(parameter);
  parameter.push_back
    (Parameter(keyword + "Port",
               Value(port, ConstraintValueType::Positive())));
  parameter.push_back
    (Parameter(keyword + "PortRange",
               Value(portRange, ConstraintValueType::Positive())));
  parameter.push_back
    (Parameter(keyword + "Proj",
               Value(projName, ConstraintValueType::NoConstraints())));
  parameter.push_back
    (Parameter(keyword + "Timeout",
               Value(timeout, ConstraintValueType::NotNegative())));
  parameter.push_back
    (Parameter(keyword + "Pause",
               Value(pause, ConstraintValueType::NoConstraints())));
}


bool OutputFAHGUI::adjustWithDefaultParameters(vector<Value> &values,
                                               const Configuration *config)
const {
  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  if (!values[0].valid()) values[0] = name;
  if (!values[2].valid()) values[2] = 52753;
  if (!values[3].valid()) values[3] = 1;
  if (!values[4].valid()) values[4] = "Protomol_3.0";
  if (!values[5].valid()) values[5] = 0;
  if (!values[6].valid()) values[6] = false;

  return checkParameters(values);
}


void OutputFAHGUI::setCoords() {
  Real x, y, z, sz;
  Vector3DBlock posMi;

  posMi.resize(app->positions.size());

  x = y = z = 0.0;
  for (unsigned int i = 0; i < app->positions.size(); i++) {
	posMi[i] = app->topology->minimalPosition(app->positions[i]);
    x += posMi[i].c[0];
    y += posMi[i].c[1];
    z += posMi[i].c[2];
  }

  sz = app->positions.size();
  x /= sz; y /= sz; z /= sz;
  for (unsigned int i = 0; i < app->positions.size(); i++) {
    server->xyz[i].x = posMi[i].c[0] - x;
    server->xyz[i].y = posMi[i].c[1] - y;
    server->xyz[i].z = posMi[i].c[2] - z;
  }
}


void OutputFAHGUI::setBonds() {
  for (unsigned int i = 0; i < app->topology->bonds.size(); i++ ) {
    //  FAH requires a < b
    if (app->topology->bonds[i].atom1 < app->topology->bonds[i].atom2) {
      server->bonds[i].a = app->topology->bonds[i].atom1;
      server->bonds[i].b = app->topology->bonds[i].atom2;

    } else {
      server->bonds[i].a = app->topology->bonds[i].atom2;
      server->bonds[i].b = app->topology->bonds[i].atom1;
    }
  }
}


void OutputFAHGUI::setAtoms() {
  float radius = 0;

  for (unsigned int i = 0; i < (app->topology->atoms).size(); i++) {
    //  Determine diamiter / set nam
    int atomNameLen = app->topology->atoms[i].name.length();
    for (int j=0;j<min(atomNameLen,4);j++)
      server->atoms[i].type[j] = (app->topology->atoms[i].name.c_str())[j];

    if (atomNameLen < 4)
        for (int j=atomNameLen;j<4;j++) server->atoms[i].type[j] = 0;

    switch (server->atoms[i].type[0]){
        case 'H':	radius = 1.2;
                    break;
        case 'C':	radius = 1.7;
                    break;
        case 'N':	radius = 1.55;
                    break;
        case 'O':	radius = 1.52;
                    break;
        case 'S':	radius = 1.85;
                    break;
        case 'P':	radius = 1.9;
                    break;
        default:	radius = 1.9;
                    break;
    }

    //  add charg
    server->atoms[i].charge = app->topology->atoms[i].scaledCharge;
    server->atoms[i].radius = radius / 2.0;
  }
}

#endif //  HAVE_GUI,HAVE_LIBFAH
