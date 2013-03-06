#include <protomol/output/OutputFAHFile.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/module/MainModule.h>
#include <protomol/module/OutputModule.h>
#include <protomol/type/String.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;


const string OutputFAHFile::keyword("FAHFile");


OutputFAHFile::OutputFAHFile() {}


OutputFAHFile::OutputFAHFile(const string &filename, int freq) :
  Output(freq), filename(filename) {}


void OutputFAHFile::doInitialize() {}


void OutputFAHFile::doRun(int step) {
  open(filename.c_str(), ios::out | ios::trunc);

  if (is_open()) {
    string str;
    str = String(app->positions.size()) + '\t' + "ProtoMol" + "\n";
    file.write(str.c_str(), str.length());

    for (unsigned int i = 0; i < app->positions.size(); i++) {
      str = String(i + 1) + '\t' +
        app->topology->atomTypes[app->topology->atoms[i].type].name + '\t' +
        String(app->positions[i].c[0]) + '\t' +
        String(app->positions[i].c[1]) + '\t' +
        String(app->positions[i].c[2]) + '\t' +
        String(1) + "\n";

      file.write(str.c_str(), str.length());
    }

    close();
  }
}


void OutputFAHFile::doFinalize(int step) {}


Output *OutputFAHFile::doMake(const vector<Value> &values) const {
  return new OutputFAHFile(values[0], values[1]);
}


bool OutputFAHFile::isIdDefined(const Configuration *config) const {
  return config->valid(getId());
}


void OutputFAHFile::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(getId(), Value(filename, ConstraintValueType::NotEmpty())));
  Output::getParameters(parameter);
}
