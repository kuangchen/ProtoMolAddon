#include <protomol/integrator/IntegratorDefinition.h>
#include <protomol/force/Force.h>

using namespace std;
using namespace ProtoMol;

//____ IntegratorDefinition
string IntegratorDefinition::print() const {
  string res = integrator.id + " {";
  for (unsigned int i = 0; i < integrator.parameters.size(); ++i) {
    res += "\n      " + integrator.parameters[i].keyword + " " +
           integrator.parameters[i].value.getString();
    if (!integrator.parameters[i].text.empty())
      res += "\t # " + integrator.parameters[i].text;
  }

  for (unsigned int i = 0; i < forces.size(); ++i) {
    res += "\n    " + Force::scope + " " + forces[i].id;
    for (unsigned int j = 0; j < forces[i].parameters.size(); ++j) {
      if (!forces[i].parameters[j].keyword.empty())
        res += "\n      " + forces[i].parameters[j].keyword;
      res += " " + forces[i].parameters[j].value.getString();
      if (!forces[i].parameters[j].text.empty())
        res += "\t # " + forces[i].parameters[j].text;
    }
  }

  res += "\n  }";
  return res;
}

