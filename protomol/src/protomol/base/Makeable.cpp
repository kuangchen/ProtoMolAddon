#include "Makeable.h"

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol;


string MakeableBase::getId() const {
  if (!alias.empty()) return getAlias();
  else return getIdNoAlias();
}


string MakeableBase::getAlias() const {
  return alias;
}


string MakeableBase::setAlias(const string &id) {
  string tmp(alias);

  if (equalNocase(id, getIdNoAlias())) alias = "";
  else alias = id;

  return tmp;
}


void MakeableBase::assertParameters(const vector<Value> &values) const {
  string err;
  vector<Parameter> tmp;
  getParameters(tmp);

  if (tmp.size() != values.size())
    err += " Expected " + toString(tmp.size()) +
      " value(s), but got " + toString(values.size()) + ".";

  for (unsigned int i = 0; i < values.size(); ++i) {

    if (!values[i].valid()) 
      err += " Parameter " + toString(i) + " '" + tmp[i].keyword +
                "' " + tmp[i].value.getDefinitionTypeString() +
                " undefined/missing or with non-valid value '" +
	values[i].getString() + "'.";

    if (!values[i].equalType(tmp[i].value)) 
      err += " Expected type " + tmp[i].value.getDefinitionTypeString() +
        " for parameter " + toString(i) + " '" + tmp[i].keyword +
        "', but got " + values[i].getDefinitionTypeString() + ".";
  }

  if (!err.empty()) THROW(getId() + ":" + err);
}


bool MakeableBase::checkParameterTypes(const vector<Value> &values) const {
  vector<Parameter> tmp;
  getParameters(tmp);
  if (tmp.size() != values.size())
    return false;

  for (unsigned i = 0; i < values.size(); ++i)
    if (!values[i].equalType(tmp[i].value) && values[i].isDefined())
      return false;

  return true;
}


bool MakeableBase::checkParameters(const vector<Value> &values) const {
  vector<Parameter> tmp;
  getParameters(tmp);

  if (tmp.size() != values.size()) return false;

  for (unsigned i = 0; i < values.size(); ++i)
    if (!values[i].valid() || !values[i].equalType(tmp[i].value))
      return false;

  return true;
}


MakeableDefinition MakeableBase::getDefinition() const {
  vector<Parameter> parameters;
  getParameters(parameters);

  return MakeableDefinition(getId(), parameters);
}

