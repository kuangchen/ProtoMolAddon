%module OutputEnergies
%{
#include <protomol/output/OutputEnergies.h>
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
%include <protomol/output/OutputEnergies.h>

%extend ProtoMol::OutputEnergies {
void uncache(ProtoMolApp* app) {
   app->outputCache.uncache();
}
};
