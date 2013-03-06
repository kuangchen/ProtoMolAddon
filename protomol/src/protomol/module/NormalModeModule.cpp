#include <protomol/module/NormalModeModule.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/type/String.h>
#include <protomol/io/EigenvectorReader.h>
#include <protomol/io/EigenvectorTextReader.h>
#include <protomol/io/XYZReader.h>
#include <protomol/integrator/normal/NormalModeLangevin.h>
#include <protomol/integrator/normal/NormalModeLangLf.h>
#include <protomol/integrator/normal/NormalModeMinimizer.h>
#include <protomol/integrator/normal/NormalModeDiagonalize.h>
#include <protomol/integrator/normal/NormalModeMori.h>
#include <protomol/integrator/normal/NormalModeRelax.h>
#include <protomol/integrator/normal/NormalModeBrownian.h>
#include <protomol/integrator/normal/NormalModeDamping.h>
#include <protomol/integrator/normal/NormalModeQuadratic.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

defineInputValue(InputEigenVectors, "eigfile");
defineInputValue(InputEigTextFile, "eigtextfile");
defineInputValue(InputEigenValues, "eigvaluefile");

void NormalModeModule::init(ProtoMolApp *app) {
  InputEigenVectors::registerConfiguration(&app->config);
  InputEigTextFile::registerConfiguration(&app->config);
  InputEigenValues::registerConfiguration(&app->config);

  app->integratorFactory.registerExemplar(new NormalModeLangevin());
  app->integratorFactory.registerExemplar(new NormalModeLangLf());
  app->integratorFactory.registerExemplar(new NormalModeMinimizer());
  app->integratorFactory.registerExemplar(new NormalModeDiagonalize());
  app->integratorFactory.registerExemplar(new NormalModeMori());
  app->integratorFactory.registerExemplar(new NormalModeRelax());
  app->integratorFactory.registerExemplar(new NormalModeBrownian());
  app->integratorFactory.registerExemplar(new NormalModeDamping());
  app->integratorFactory.registerExemplar(new NormalModeQuadratic());
}

void NormalModeModule::read(ProtoMolApp *app) {
  Configuration &config = app->config;

  // Eigenvectors/value
  if (config.valid(InputEigTextFile::keyword)) {
    EigenvectorTextReader evTextReader;

    if (config.valid(InputEigTextFile::keyword)) {
      if (!evTextReader.open(config[InputEigTextFile::keyword]))
        THROWS("Can't open eigenvector file '"
               << config[InputEigTextFile::keyword].getString() << "'.");

      if (!(evTextReader >> app->eigenInfo)) {
        THROWS("Could not parse eigenvector file '"
               << (string)config[InputEigTextFile::keyword] << "'");

        if (app->eigenInfo.myEigenvectorLength != (double)app->positions.size())
          THROWS("Eigenvector length is wrong, should be "
                 << app->positions.size() << " got "
                 << app->eigenInfo.myEigenvectorLength << ".");

        if (app->eigenInfo.myNumEigenvectors < 1 ||
            app->eigenInfo.myNumEigenvectors > (double)app->positions.size())
          THROWS("Wrong number of eigenvectors ("
                 << app->eigenInfo.myNumEigenvectors << ").");
      }

      report << plain << "Using eigfile '" << config[InputEigTextFile::keyword]
             << "' (" << app->eigenInfo.myEigenvectorLength << ")." << endr;
    }

  } else if (config.valid(InputEigenVectors::keyword)) {
    EigenvectorReader evReader;

    if (!evReader.open(config[InputEigenVectors::keyword]))
      THROWS("Can't open eigenvector file '"
             << (string)config[InputEigenVectors::keyword] << "'.");

    if (!(evReader >> app->eigenInfo)) {
      if (app->eigenInfo.myEigenvectorLength != (double)app->positions.size())
        THROWS("Eigenvector length is wrong, should be "
               << app->positions.size() << " got "
               << app->eigenInfo.myEigenvectorLength << ".");

      if (app->eigenInfo.myNumEigenvectors < 1 ||
          app->eigenInfo.myNumEigenvectors > (double)app->positions.size())
        THROWS("Wrong number of eigenvectors ("
               << app->eigenInfo.myNumEigenvectors << ").");
    }

    report << plain << "Using eigfile '"
           << config[InputEigenVectors::keyword] << "' ("
           << app->eigenInfo.myEigenvectorLength << ")." << endr;

  }

  // Eigenvalues file
  if (config.valid(InputEigenValues::keyword)) {
    // eigenvalues
    XYZReader valReader;
    Vector3DBlock tempEVal;
    if (!valReader.open(config[InputEigenValues::keyword])){
      THROWS(string("Can't open eigenvalue file '") +
        config[InputEigenValues::keyword].getString() + "'.");
    }

    //if (!(valReader >> app->eigenInfo.myEigenvalues)){
    if (!(valReader >> tempEVal)){
        THROWS(string("Could not parse eigenvalue file '") +
          config[InputEigenValues::keyword].getString() + "'. ");
    }

    //copy data accross
    int evsize = tempEVal.size() * 3;
    for ( int i=0; i<evsize; i++ ){
      app->eigenInfo.myEigenvalues.push_back(tempEVal.c[i]);
    }

    report << plain << "Using eigvaluefile '"
           << config[InputEigenValues::keyword] << "' ("
           << app->eigenInfo.myEigenvalues.size() << ")." << endr;

  }
}

