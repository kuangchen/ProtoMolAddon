#ifndef PROTOMOLAPP_H
#define PROTOMOLAPP_H

#include <protomol/factory/TopologyFactory.h>
#include <protomol/factory/ForceFactory.h>
#include <protomol/factory/IntegratorFactory.h>
#include <protomol/factory/OutputFactory.h>

#include <protomol/config/Configuration.h>
#include <protomol/config/CommandLine.h>
#include <protomol/output/OutputCache.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/type/EigenvectorInfo.h>
#include <protomol/type/PSF.h>
#include <protomol/type/PAR.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/CoulombSCPISMParameterTable.h>

#include <ostream>

namespace ProtoMol {
  class OutputCollection;
  class Integrator;
  class GenericTopology;
  class ModuleManager;

  class ProtoMolApp {
  public:
    ModuleManager *modManager;

    // Data
    Vector3DBlock positions;
    Vector3DBlock velocities;
    EigenvectorInfo eigenInfo;
    PSF psf;
    PAR par;
    ScalarStructure energies;
    CoulombSCPISMParameterTable *SCPISMParameters;

    // Factories
    TopologyFactory topologyFactory;
    ForceFactory forceFactory;
    IntegratorFactory integratorFactory;
    OutputFactory outputFactory;

    // Containers
    CommandLine cmdLine;
    Configuration config;
    mutable OutputCache outputCache;
    OutputCollection *outputs;
    Integrator *integrator;
    GenericTopology *topology;

    // Run
    int currentStep;
    int lastStep;

    ProtoMolApp() {}
    ProtoMolApp(ModuleManager *modManager);
    ~ProtoMolApp();

    static void splash(std::ostream &stream);
    void load(const std::string &configfile);
    bool load(int argc = 0, char *argv[] = 0);
    bool load(const std::vector<std::string> &args);
    void configure();
    void configure(const std::string &configfile);
    bool configure(int argc, char *argv[]);
    bool configure(const std::vector<std::string> &args);
    void build();
    void print(std::ostream &stream);
    bool step(unsigned inc = 0);
    void finalize();
  };
};
#endif // PROTOMOLAPP_H
