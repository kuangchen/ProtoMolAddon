#include <protomol/module/IOModule.h>

#include <protomol/io/PosVelReader.h>
#include <protomol/io/PSFReader.h>
#include <protomol/io/PARReader.h>

//for GROMACS
#include <protomol/io/gromacs/PortGromacsParameters.h>
#include <protomol/type/GromacsTopology.h>
#include <protomol/type/GromacsParameters.h>
//
#include <protomol/config/Configuration.h>
#include <protomol/base/Exception.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>
#include <protomol/type/String.h>
#include <protomol/type/PDB.h>
#include <protomol/module/MainModule.h>
#include <protomol/io/SCPISMReader.h>
#include <protomol/topology/CoulombSCPISMParameterTable.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;

defineInputValue(InputVelocities, "velfile");
defineInputValueWithAliases(InputPositions, "posfile",
  ("coords")("coordinates"));
defineInputValueWithAliases(InputPSF, "psffile", ("structure"));
defineInputValueWithAliases(InputPAR, "parfile", ("parameters"));
defineInputValue(InputPDBScaling, "pdbScaling");
defineInputValue(InputDihedralMultPSF, "dihedralMultPSF");
defineInputValue(InputSCPISM, "scpismfile");

//for GROMACS
defineInputValue(InputGromacsTopo, "gromacstopologyfile");
defineInputValue(InputGromacsParamPath, "gromacsparameterpath");
defineInputValue(InputGromacsGBParameterFile, "gromacsGBparameterfile");

defineInputValue(InputGromacsTprFile, "gromacstprfile");

void IOModule::init(ProtoMolApp *app) {
  Configuration *config = &app->config;
  
  InputPositions::registerConfiguration(config);
  InputVelocities::registerConfiguration(config);
  InputPSF::registerConfiguration(config);
  InputPAR::registerConfiguration(config);
  InputPDBScaling::registerConfiguration(config);
  InputDihedralMultPSF::registerConfiguration(config);  
  InputSCPISM::registerConfiguration(config);

  //for GROMACS
  InputGromacsTopo::registerConfiguration(config);
  InputGromacsParamPath::registerConfiguration(config);
  InputGromacsGBParameterFile::registerConfiguration(config);

  InputGromacsTprFile::registerConfiguration(config);

}

void IOModule::read(ProtoMolApp *app) {
  Configuration &config = app->config;

  // Positions
  PosVelReader reader;
  if (!reader.open(config[InputPositions::keyword]))
    THROW(string("Can't open position file '") +
      config[InputPositions::keyword].getString() + "'.");

  if (reader.tryFormat(PosVelReaderType::PDB)) {
    PDB pdb;
	if (!(reader >> pdb)) {
      THROW(string("Could not parse PDB position file '") +
        config[InputPositions::keyword].getString() + "'.");
	}
    
    // Use swap function of Vector3DBlock
    app->positions.swap(pdb.coords);

    // Add to output cache
    app->outputCache.add(pdb.atoms);

  } else if (reader.tryFormat(PosVelReaderType::XYZ) ||
             reader.tryFormat(PosVelReaderType::XYZBIN)) {
	  XYZ xyz;
	  if (!(reader >> xyz)) {
		THROW(string("Could not parse position file '") +
		  config[InputPositions::keyword].getString() +
		  "'. Supported formats are : " +
		  PosVelReaderType::getPossibleValues(", ") + ".");
	  }

	  app->positions.swap(xyz.coords);

  } else THROW("Error reading positions file");

  report << plain << "Using " << reader.getType() << " posfile '"
         << config[InputPositions::keyword] << "' ("
         << app->positions.size() << ")." << endr;

  // Velocities
  if (config.valid(InputVelocities::keyword)) {
    if (!reader.open(config[InputVelocities::keyword]))
      THROW(string("Can't open velocity file '") +
        config[InputVelocities::keyword].getString() + "'.");

    if (!(reader >> app->velocities))
      THROW(string("Could not parse velocity file '") +
        config[InputVelocities::keyword].getString() +
        "'. Supported formats are : " +
        PosVelReaderType::getPossibleValues(", ") + ".");

    report << plain << "Using " << reader.getType() << " velfile '"
           << config[InputVelocities::keyword] << "' ("
           << app->velocities.size() << ")." << endr;

    if (reader.getType() == "PDB" && (bool)config[InputPDBScaling::keyword]) {
      for (unsigned int i = 0; i < app->velocities.size(); i++)
        app->velocities[i] /= PDBVELSCALINGFACTOR;

      report << plain << "PDB velocities scaled." << endr;
    }

  } else if (config.valid(InputTemperature::keyword)) {
    app->velocities.resize(app->positions.size());

    report << plain << "Using temperature "
           << config[InputTemperature::keyword] << "K for the velocities  ("
           << app->velocities.size() << ")." << endr;
    // Create velocities later, we need the topology for that ...

  } else THROW("Neither temperature nor velocity file specified.");

  // Gromacs/AMBER input files
  if (config.valid(InputGromacsTopo::keyword) && 
        config.valid(InputGromacsParamPath::keyword)){
     GromacsTopology gTopo;
     GromacsParameters gParams;
     PortGromacsParameters gromacs_port;
       
     if (!gromacs_port.Read_Basic_Gromacs_Parameters
         (app->psf,app->par, gTopo, gParams, 
          (const string)config[InputGromacsTopo::keyword],
          (const string)config[InputGromacsParamPath::keyword])){
       THROW(string("Cant read GROMACS parameters into PSF and PAR"));
     }

     //check if Generalized Born Parameter file has been passed
     if (config.valid(InputGromacsGBParameterFile::keyword)) {
         if (!gromacs_port.Read_Gromacs_GB_Parameters
             ((const string)config[InputGromacsGBParameterFile::keyword])){
           THROW(string("Cant read Gromacs parameters for Generalized Born"));
         }
     }

  } else {
    // PSF
    PSFReader psfReader;
    if (!psfReader.open(config[InputPSF::keyword]))
      THROWS("Can't open PSF file '"
             << (string)config[InputPSF::keyword] << "'.");

    if (!(psfReader >> app->psf))
      THROWS("Could not parse PSF file '"
             << (string)config[InputPSF::keyword] << "'.");

    report << plain << "Using PSF file '" << (string)config[InputPSF::keyword]
           << "' (" << app->psf.atoms.size() << ")." << endr;

    // PAR
    PARReader parReader;
    if (!parReader.open(config[InputPAR::keyword]))
      THROW(string("Can't open PAR file '") +
        config[InputPAR::keyword].getString() + "'.");

    if (!(parReader >> app->par))
      THROW(string("Could not parse PAR file '") +
        config[InputPAR::keyword].getString() + "'.");

    report << plain << "Using PAR file '" << config[InputPAR::keyword]
           << "', " << (parReader.getCharmmTypeDetected() != PAR::CHARMM28 ?
                        "old" : "new") << " charmm force field.";

    if (!config[InputDihedralMultPSF::keyword].valid())
      config[InputDihedralMultPSF::keyword] =
        (parReader.getCharmmTypeDetected() != PAR::CHARMM28);

    if (config[InputDihedralMultPSF::keyword])
      report << " Dihedral multiplictity defined by PSF.";
    report << endr;

  }

  //SCPISM
  if (config.valid(InputSCPISM::keyword)) {

    app->SCPISMParameters = new CoulombSCPISMParameterTable;

    SCPISMReader mReader(config[InputSCPISM::keyword]);

    if(!mReader.read( app->SCPISMParameters->myData ))
      report << error << "Invalid SCPISM parameter file "
             << config[InputSCPISM::keyword] << "." << endr;

    report << plain << "Using SCPISM file '" << config[InputSCPISM::keyword]
           << "'." << endr;


  }

  // Test input
  if (app->positions.size() != app->velocities.size() ||
      app->positions.size() != app->psf.atoms.size())
    THROWS("Positions, velocities and PSF input have different number "
           "of atoms. positions=" << app->positions.size()
           << " velocities=" << app->velocities.size()
           << " atoms=" << app->psf.atoms.size());
}
