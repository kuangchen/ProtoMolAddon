%module OutputFAHGUI
%{
#include <protomol/output/OutputFAHGUI.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/output/Output.h>
#include <protomol/output/OutputFile.h>
#include <protomol/output/OutputCache.h>
using namespace ProtoMol;
using namespace std;
using std::string;
%}


%include <protomol/type/Real.h>
%include "std_string.i"

%include <protomol/output/Output.h>
%include <protomol/output/OutputFile.h>
%include <protomol/output/OutputCache.h>
%include <protomol/output/OutputFAHGUI.h>

%extend ProtoMol::OutputFAHGUI {
void uncache(ProtoMolApp* app) {
   app->outputCache.uncache();
}
};
