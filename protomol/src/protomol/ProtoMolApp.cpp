#include <protomol/ProtoMolApp.h>

#include <protomol/base/ModuleManager.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/TimerStatistic.h>
#include <protomol/base/Zap.h>
#include <protomol/base/Report.h>

#include <protomol/module/MainModule.h>
#include <protomol/module/IOModule.h>
#include <protomol/module/ConfigurationModule.h>

#include <protomol/type/String.h>

#include <protomol/config/CommandLine.h>
#include <protomol/config/Configuration.h>

#include <protomol/io/ConfigurationReader.h>

#include <protomol/factory/TopologyFactory.h>
#include <protomol/factory/OutputFactory.h>

#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/BuildTopology.h>
#include <protomol/topology/BuildTopologyFromTpr.h>
#include <protomol/topology/TopologyUtilities.h>

#include <protomol/output/OutputCollection.h>

#include <protomol/parallel/Parallel.h>

#include <iomanip>
#ifdef HAVE_PACKAGE_H
#include <protomol/package.h>
#endif

#ifdef HAVE_LIBFAH
#include <fah/core/Core.h>
#endif

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;


ProtoMolApp::ProtoMolApp(ModuleManager *modManager) :
  modManager(modManager), SCPISMParameters(0), cmdLine(&config), outputs(0),
  integrator(0), topology(0) {

  TimerStatistic::timer[TimerStatistic::WALL].start();

  modManager->init(this);

  topologyFactory.registerAllExemplarsConfiguration(&config);
  outputFactory.registerAllExemplarsConfiguration(&config);
}


ProtoMolApp::~ProtoMolApp() {}


void ProtoMolApp::splash(ostream &stream) {
	if( Parallel::iAmMaster() ){
		const int w = 16;
		stream
		<< headerRow("ProtoMol") << endl
		<< setw(w) << "Description: ";
		fillFormat(stream, "A rapid PROTOtyping MOLecular dynamics object-oriented "
				 "component based framework.", w, w);
		stream
#ifdef HAVE_PACKAGE_H
		<< setw(w) << "Version: " << PACKAGE_VERSION << endl
		<< setw(w) << "SVN revision: " << PACKAGE_REVISION << endl
		<< setw(w) << "Repository: " << PACKAGE_SOURCE_REPO << endl
		<< setw(w) << "Homepage: " << PACKAGE_HOMEPAGE << endl
		<< setw(w) << "Report bugs to: " << PACKAGE_BUGREPORT << endl
		<< setw(w) << "Compiler: " << PACKAGE_COMPILER << " "
		<< PACKAGE_COMPILER_VERSION << endl
		<< setw(w) << "Flags: " << PACKAGE_COMPILER_FLAGS << endl
		<< setw(w) << "Extra libs: " << PACKAGE_COMPILER_LIBS << endl
		<< setw(w) << "Built by: " << PACKAGE_BUILT_BY << endl
		<< setw(w) << "Build platform: " << PACKAGE_PLATFORM << endl
#endif // HAVE_PACKAGE_H
		<< setw(w) << "Build date: " <<  __DATE__ << ", " << __TIME__ << endl
#ifdef HAVE_LIBFAH
		<< setw(w) << "Checksumming: "
		<< "Enabled for Folding@Home file protection." << endl
#endif // HAVE_LIBFAH
#ifdef HAVE_PACKAGE_H
		<< setw(w) << "Please cite: ";
		fillFormat(stream, PACKAGE_CITE, w, w);
		stream
#endif // HAVE_PACKAGE_H
		<< PROTOMOL_HR << endl;
	}
}


void ProtoMolApp::load(const string &configfile) {
  vector<string> args;

  args.push_back("ProtoMol");
  args.push_back(configfile);
  load(args);
}


bool ProtoMolApp::load(int argc, char *argv[]) {
	Parallel::init( argc, argv );
  return load(vector<string>(argv, argv + argc));
}


bool ProtoMolApp::load(const vector<string> &args) {
  // Parse command line
  if (cmdLine.parse(args)) return false;

  return true;
}


void ProtoMolApp::configure() {
  // Read Config file
  if (config.valid(InputConfig::keyword))
    SystemUtilities::chdir(SystemUtilities::dirname
                           (config[InputConfig::keyword]));
  else THROW("Configuration file not set.");

  modManager->configure(this);
}


void ProtoMolApp::configure(const string &configfile) {
  load(configfile);
  configure();
}


bool ProtoMolApp::configure(int argc, char *argv[]) {
  if (!load(argc, argv)) return false;
  configure();
  return true;
}


bool ProtoMolApp::configure(const vector<string> &args) {
  if (!load(args)) return false;

  configure();
  
  return true;
}


void ProtoMolApp::build() {
  // TPR input for topology, positions and velocities?
  // Then check for Gromacs support
#if !defined(HAVE_GROMACS)
  if (config.valid(InputGromacsTprFile::keyword))
    THROWS("GROMACS support not available for '" <<
           config[InputGromacsTprFile::keyword].getString() << "'.");
#endif

  //floag for TPR
  bool GROMACRTPR(false);

  //test TPR file
  if( config.valid(InputGromacsTprFile::keyword) ) GROMACRTPR = true;

  // Read data if not TPR
  if( !GROMACRTPR ) modManager->read(this);

  // Build topology
  try {
    topology = topologyFactory.make(&config);
  } catch (const Exception &e) {

    if( !GROMACRTPR ){
      // Try to get some defaults with the postions known ...
      const GenericTopology *prototype =
        topologyFactory.find(config[GenericTopology::keyword].getString());

      if (prototype) {
        vector<Parameter> parameters = prototype->getDefaults(positions);

        for (unsigned int i = 0; i < parameters.size(); i++)
          if (!config.valid(parameters[i].keyword) &&
              parameters[i].value.valid()) {
            config.set(parameters[i].keyword, parameters[i].value);
            report << hint << parameters[i].keyword << " undefined, using "
                   << parameters[i].value.getString() << "." << endr;
          }

        topology = topologyFactory.make(&config);
      }
    }

    if (!topology) throw e;
  }

  if (!GROMACRTPR) {
    // Using SCPISM parameter? Flag or filename
    if (config[InputDoSCPISM::keyword] || SCPISMParameters) {

      if(config[InputDoSCPISM::keyword])
        topology->doSCPISM = config[InputDoSCPISM::keyword];

      if ((topology->doSCPISM < 1 || topology->doSCPISM > 3) &&
          !SCPISMParameters)
        THROW("doscpism should be between 1 and 3 or an input file should be "
              "used.");

      if  (SCPISMParameters) {
        if (!config[InputDoSCPISM::keyword]) topology->doSCPISM = 4;
      } else {
        SCPISMParameters = new CoulombSCPISMParameterTable;
        SCPISMParameters->populateTable();
      }

      report
        << "SCPISM: doSCPISM set to " << topology->doSCPISM << "." << endr;

      if (topology->doSCPISM == 3) {
        // Quartic switch parameters
        SCPISMParameters->myData["H"].hbond_factor = 0.4695;
        SCPISMParameters->myData["HC"].hbond_factor = 7.2560;
      }

      SCPISMParameters->displayTable();

    }

    // Using GBSA with openMM
    if (config[InputDoGBSAObc::keyword]) {
      topology->doGBSAOpenMM = 1;
      topology->obcType = config[InputDoGBSAObc::keyword];

      if (topology->obcType == 1) {
        topology->alphaObc = 0.8;
        topology->betaObc = 0;
        topology->gammaObc = 2.91;

      } else if (topology->obcType == 2) {
        //obctype = 2
        topology->alphaObc = 1.0;
        topology->betaObc = 0.8;
        topology->gammaObc = 4.85;
      }
      topology->dielecOffset = 0.09;

    }

    // Find force field type before building topology
    if (config.valid(InputGromacsTopo::keyword) &&
        config.valid(InputGromacsParamPath::keyword)) {
      topology->forceFieldFlag = GROMACS;
    }

    // Build the topology
    buildTopology(topology, psf, par, config[InputDihedralMultPSF::keyword],
                  SCPISMParameters);

  } else {
    //TPR/GROMACS if here
    topology->forceFieldFlag = GROMACS;//TPR;

    // Build the topology from the tpr file
    buildTopologyFromTpr( topology, positions, velocities,
                          config[InputGromacsTprFile::keyword].getString() );
  }

  // Register Forces
  modManager->registerForces(this);

  // Build the integrators and forces
  integrator =
    integratorFactory.make(config[InputIntegrator::keyword], &forceFactory);

  // Setup run parameters (used for GUI so required here)
  currentStep = config[InputFirststep::keyword];
  lastStep = currentStep + (int)config[InputNumsteps::keyword];

  // Create outputs
  // TODO if !Parallel::iAmMaster() turn off outputs
	if( !Parallel::iAmMaster() ){
		outputs = new OutputCollection;
	}else{
		if (config[InputOutput::keyword]){
			outputs = outputFactory.makeCollection(&config);
		}
	}
  
  // Post build processing
  modManager->postBuild(this);

  report << plain << "Actual start temperature : "
         << temperature(topology, &velocities) << "K" << endr;


  // Add Integrator Modifiers
  modManager->addModifiers(this);

  // Initialize
  energies.molecularVirial(config[InputMolVirialCalc::keyword]);
  energies.virial(config[InputVirialCalc::keyword]);
  report << plain << "Virial tensor : " << energies.virial() << endr;
  report << plain << "Molecular virial tensor : "
         << energies.molecularVirial() << endr;

  topology->time =
    (Real)config[InputFirststep::keyword] * integrator->getTimestep();

  integrator->initialize(this);
  outputs->initialize(this);
  outputCache.initialize(this);

  // Init cache
  //outputs->addToCache(pdbAtoms); // TODO fix this
  outputCache.add(psf);
  outputCache.add(par);

  // Print Factories
  if ((int)config[InputDebug::keyword] >= 5  &&
      (int)config[InputDebugLimit::keyword] <= 5)
    cout
      << headerRow("Factories")     << endl
      << headerRow("Configuration") << endl << config            << endl
      << headerRow("Topology")      << endl << topologyFactory   << endl
      << headerRow("Integrator")    << endl << integratorFactory << endl
      << headerRow("Force")         << endl << forceFactory      << endl
      << headerRow("Output")        << endl << outputFactory     << endl;

  // Clear all factories
  topologyFactory.unregisterAllExemplars();
  integratorFactory.unregisterAllExemplars();
  forceFactory.unregisterAllExemplars();
  outputFactory.unregisterAllExemplars();

  TimerStatistic::timer[TimerStatistic::RUN].reset();
  TimerStatistic::timer[TimerStatistic::INTEGRATOR].reset();
  TimerStatistic::timer[TimerStatistic::FORCES].reset();
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].reset();
  TimerStatistic::timer[TimerStatistic::IDLE].reset();
}


bool ProtoMolApp::step(unsigned inc) {
  if (currentStep >= lastStep) return false;

  TimerStatistic::timer[TimerStatistic::RUN].start();

  if (outputs->run(currentStep)) {
#ifdef HAVE_LIBFAH
    // Make sure we save the latest checksum information after writing data.
    if (FAH::Core::isActive()) FAH::Core::instance().checkpoint();
#endif
  }

  if (!inc) inc = outputs->getNext() - currentStep;
  inc = std::min(lastStep, (int)(currentStep + inc)) - currentStep;

  TimerStatistic::timer[TimerStatistic::INTEGRATOR].start();

  integrator->run(inc);

  TimerStatistic::timer[TimerStatistic::INTEGRATOR].stop();

  // moved here so that current step is valid in integrator
  currentStep += inc;

  TimerStatistic::timer[TimerStatistic::RUN].stop();

  return true;
}


void ProtoMolApp::finalize() {
  outputs->finalize(currentStep);

  // Clean up
  zap(topology);
  zap(integrator);
  zap(outputs);
  zap(SCPISMParameters);

  TimerStatistic::timer[TimerStatistic::WALL].stop();

	if( Parallel::iAmMaster() ){
		report << plain << "Timing: " << TimerStatistic() << "." << endr;
	}
	
	Parallel::finalize();
}


void ProtoMolApp::print(ostream &stream) {
	if( Parallel::iAmMaster() ){
		// Output
		stream << headerRow("Outputs") << endl;

		for (OutputCollection::const_iterator itr =
			 const_cast<const OutputCollection *>(outputs)->begin();
		   itr != const_cast<const OutputCollection*>(outputs)->end(); itr++) {

		stream << "Output " << (*itr)->getId();

		vector<Parameter> parameters;
		(*itr)->getParameters(parameters);
		for (unsigned int i = 0; i < parameters.size(); i++)
		  stream << " " << parameters[i].value.getString();
		stream << "." << endl;
		}

		if (!((bool)config[InputOutput::keyword]))
		stream << "All output suppressed!" << endl;


		// Integrator
		stream << headerRow("Integrator") << endl;
		vector<IntegratorDefinition> inter = integrator->getIntegratorDefinitionAll();
		stream  << InputIntegrator::keyword << " {" << endl;

		for (int i = inter.size() - 1; i >= 0; i--)
		stream << Constant::PRINTINDENT << "Level "
			   << i << " " << inter[i].print() << endl;

		stream << "}" << endl;


		// Topology
		stream << headerRow("Topology") << endl;
		stream << topology->print(&positions) << endl;

		stream << PROTOMOL_HR << endl;
	}
}
