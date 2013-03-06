#include "OutputEnergies.h"

#include <protomol/config/Configuration.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/output/OutputCache.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/ProtoMolApp.h>

#include <iomanip>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

const string OutputEnergies::keyword("allEnergiesFile");


OutputEnergies::OutputEnergies(const string &filename, int freq,
                               bool doMolTemp) :
  Output(freq), filename(filename), doMolecularTemperature(doMolTemp) {}


void OutputEnergies::doInitialize() {
  SystemUtilities::ensureDirectory(SystemUtilities::dirname(filename));

#ifdef HAVE_LIBFAH
  cb::File allEnergiesHeaderFile;
#else
  ofstream allEnergiesHeaderFile;
#endif

  allEnergiesHeaderFile.open((filename + ".header").c_str(), ios::out);
  if (!allEnergiesHeaderFile)
    report << error << " Can not open \'" << filename << ".header\' for "
           << getId() << "." << endr;

  allEnergiesHeaderFile
    << setw(14) << "Time(fs)" << " "
    << setw(14) << "E_potential" << " "
    << setw(14) << "E_kinetic" << " "
    << setw(14) << "E_total" << " "
    << setw(14) << "Temperature" << " "
    << setw(14) << "E_bond" << " "
    << setw(14) << "E_angle" << " "
    << setw(14) << "E_dihedral" << " "
    << setw(14) << "E_improper" << " "
    << setw(14) << "E_VdW" << " "
    << setw(14) << "E_coulomb" << " "
    << setw(14) << "E_other" << " "
    << setw(14) << "Volume(A^3)";

  if (app->energies.virial())
    allEnergiesHeaderFile  << " " << setw(14) << "Pressure(bar)";

  if (app->energies.molecularVirial())
    allEnergiesHeaderFile << " " << setw(14) << "Mol_Pres(bar)";

  if (doMolecularTemperature)
    allEnergiesHeaderFile << " " << setw(14) << "Mol_Temp(K)";

  allEnergiesHeaderFile << " " << setw(20) << "E_shadow" << endl;

  allEnergiesHeaderFile.close();

  file.open(filename.c_str(), ios::out);
  if (!file) THROWS("Failed to open all energies file '" << filename << "'");
}


void OutputEnergies::doRun(int) {
  file
    << resetiosflags(ios::showpoint | ios::fixed | ios::floatfield)
    << setw(14) << setprecision(2) << setiosflags(ios::showpoint | ios::fixed)
    << app->outputCache.getTime() << " "
    << resetiosflags(ios::showpoint | ios::fixed | ios::floatfield)
    << setiosflags(ios::floatfield) << setprecision(8)
    << setw(14) << app->outputCache.getPotentialEnergy() << " "
    << setw(14) << app->outputCache.getKineticEnergy() << " "
    << setw(14) << app->outputCache.getTotalEnergy() << " "
    << setw(14) << app->outputCache.getTemperature() << " "
    << setw(14) << app->energies[ScalarStructure::BOND] << " "
    << setw(14) << app->energies[ScalarStructure::ANGLE] << " "
    << setw(14) << app->energies[ScalarStructure::DIHEDRAL] << " "
    << setw(14) << app->energies[ScalarStructure::IMPROPER] << " "
    << setw(14) << app->energies[ScalarStructure::LENNARDJONES] << " "
    << setw(14) << app->energies[ScalarStructure::COULOMB] << " "
    << setw(14) << app->energies[ScalarStructure::OTHER] << " "
    << setw(14) << app->outputCache.getVolume();

  if (app->energies.virial())
    file << " " << setw(14) << app->outputCache.getPressure();

  if (app->energies.molecularVirial())
    file << " " << setw(14) << app->outputCache.getMolecularPressure();

  if (doMolecularTemperature)
    file << " " << setw(14) << app->outputCache.getMolecularTemperature();

  file
    << " " << setw(20) << setprecision(16) //   High precision needed.
    << app->energies[ScalarStructure::SHADOW]
    << endl;
}


void OutputEnergies::doFinalize(int step) {
  file.close();
}


Output *OutputEnergies::doMake(const vector<Value> &values) const {
  return new OutputEnergies(values[0], values[1], values[2]);
}


void OutputEnergies::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(keyword, Value(filename, ConstraintValueType::NotEmpty())));
  Output::getParameters(parameter);
  parameter.push_back(Parameter("molecularTemperature",
                                Value(doMolecularTemperature), false));
}
