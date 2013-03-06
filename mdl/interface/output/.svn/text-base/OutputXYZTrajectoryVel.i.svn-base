%module OutputXYZTrajectoryVel
%{
#include <protomol/output/OutputXYZTrajectoryVel.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/output/Output.h>
#include <protomol/output/OutputFile.h>
#include <protomol/output/OutputCache.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>

%include "std_string.i"
%include <protomol/output/Output.h>
%include <protomol/output/OutputFile.h>
%include <protomol/output/OutputCache.h>
%include <protomol/output/OutputXYZTrajectoryVel.h>

%extend ProtoMol::OutputXYZTrajectoryVel {
void uncache(ProtoMolApp* app) {
   app->outputCache.uncache();
}
};
