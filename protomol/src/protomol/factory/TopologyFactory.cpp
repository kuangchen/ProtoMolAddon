#include <protomol/factory/TopologyFactory.h>

#include <protomol/config/Configuration.h>
#include <protomol/factory/HelpTextFactory.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

//____ TopologyFactory
void TopologyFactory::
  registerAllExemplarsConfiguration(Configuration *config) const {
  for (const_iterator i = begin(); i != end(); ++i) {
    vector<Parameter> parameter;
    (*i)->getParameters(parameter);
    for (unsigned int i = 0; i < parameter.size(); i++)
      config->registerKeyword(parameter[i].keyword, parameter[i].value);
  }

  config->registerKeyword(GenericTopology::getKeyword(),
    Value(string(""), ConstraintValueType::NotEmpty()));
  cache = false;
}

GenericTopology *TopologyFactory::make(const Configuration *config) const {
  string id = config->get(GenericTopology::getKeyword()).getString();
  const GenericTopology *prototype = getPrototype(id);

  vector<Parameter> parameters;
  if (prototype) prototype->getParameters(parameters);

  return make(id, config->get(parameters));
}

GenericTopology *TopologyFactory::make(const string &id,
                                       const vector<Value> &values) const {
  const GenericTopology *prototype = getPrototype(id);

  if (!prototype)
    THROWS(" Could not find any match for '" << id << "' in "
           << GenericTopology::scope << "Factory.\nPossible topologies are:\n"
           << *this);

  // Make
  GenericTopology *newObj = prototype->make(values);
  if (!newObj) THROW("Could not make topology from prototype.");

  // Adjust external alias
  newObj->setAlias(id);
  return newObj;
}

void TopologyFactory::registerHelpText() const {
  for (exemplars_t::const_iterator i = exemplars.begin();
       i != exemplars.end(); ++i) {
    HelpText helpText;
    i->second->getParameters(helpText.parameters);
    helpText.id = i->second->getIdNoAlias();
    helpText.text = i->second->getText();
    helpText.scope = i->second->getScope();
    HelpTextFactory::registerExemplar(i->second->getId(), helpText);

    HelpText alias;
    alias.text = "alias for \'" + i->second->getId() + "\'";
    alias.scope = i->second->getScope();
    for (exemplars_t::const_iterator j = aliasExemplars.begin();
         j != aliasExemplars.end(); ++j)
      if (j->second->getIdNoAlias() == i->second->getIdNoAlias()) {
        alias.id = j->first;
        HelpTextFactory::registerExemplar(alias.id, alias);
      }
  }
}

