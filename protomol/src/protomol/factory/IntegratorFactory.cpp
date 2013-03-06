#include <protomol/factory/IntegratorFactory.h>

#include <protomol/integrator/STSIntegrator.h>
#include <protomol/integrator/MTSIntegrator.h>
#include <protomol/integrator/NonStandardIntegrator.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/factory/ForceFactory.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/base/Report.h>
#include <protomol/factory/HelpTextFactory.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

Integrator *IntegratorFactory::make(const string &definition,
                                    ForceFactory *forceFactory) const {
  string str;
  vector<IntegratorInput> integratorInput;
  stringstream ss(definition);

  while (ss >> str) {
    // Parse level and integrator type
    string levelStr, integrator, d;
    ss >> levelStr >> integrator >> d;
    if (!(equalNocase(str, "level") && isUInt(levelStr) &&
          !integrator.empty() && d == "{"))
      THROW("Integration definition mismatch, expecting \'level <number> "
            "<integrator> { ... }.");

    const Integrator *prototype = getPrototype(integrator);
    if (prototype == NULL)
      THROWS(" Could not find any match for '" << integrator << "' in "
             << Integrator::scope << "Factory. Possible integrators are:\n"
             << *this);

    // Read first integrator parameters and then force definitions
    string parameterStr, forceStr;
    while (ss >> str) {
      if (str == "}") break;

      if (equalNocase(str, "Force") || !forceStr.empty())
        forceStr += (forceStr.empty() ? "" : " ") + str;

      else parameterStr += (parameterStr.empty() ? "" : " ") + str;
    }

    parameterStr += " ";  // some compilers need this

    // Expand vector
    unsigned int level = toUInt(levelStr);
    if (integratorInput.size() <= level)
      integratorInput.resize(level + 1);

    // Return if already defined
    if (integratorInput[level].prototype != NULL)
      THROW(" Level " + toString(level) + " already defined with " +
        integratorInput[level].prototype->getId() + ".");

    // Parse integrator parameter
    vector<Parameter> parameters;
    prototype->getParameters(parameters);
    integratorInput[level].values.resize(parameters.size());

    stringstream ssp(parameterStr);
    bool foundLast = false;

    for (unsigned int i = 0; i < parameters.size(); ++i) {
      // parse the the parameters one by one
      string strp;
      integratorInput[level].values[i] = parameters[i].value;
      integratorInput[level].values[i].clear();

      bool found = false;
      bool retry = true;
      if (!(ssp)) {
        ssp.str(parameterStr);
        ssp.seekg(ios::beg);
        ssp.clear();
      }

      while (ssp >> strp || retry) {
        if (!(ssp) && retry) {
          ssp.str(parameterStr);
          ssp.seekg(ios::beg);
          ssp.clear();
          retry = false;
          strp = "";
          continue;
        }

        if (equalNocase(strp, parameters[i].keyword) &&
            !parameters[i].keyword.empty()) {
          ssp >> integratorInput[level].values[i];
          found = true;
          break;

        } else if (foundLast && parameters[i].keyword.empty()) {
          ssp.seekg((-1) * static_cast<int>(strp.size()), ios::cur);
          ssp.clear();
          ssp >> integratorInput[level].values[i];
          found = true;
          break;
        }
      }

      foundLast = found;

      // If still undefined take default value is available
      if (!found && parameters[i].defaultValue.valid())
        integratorInput[level].values[i].set(parameters[i].defaultValue);
    }

    try {
      prototype->assertParameters(integratorInput[level].values);
    } catch (const Exception &e) {
      THROW(string(" Level " + toString(level) + " " + e.getMessage()));
    }

    // Parse forces
    stringstream ssf(forceStr);
    Value force(ValueType::Force(""));
    while (ssf >> force)
      integratorInput[level].forces.push_back(force.getString());

    // Set prototype
    integratorInput[level].prototype = prototype;
  }

  // Check if we have a definition for each level
  string err;
  for (unsigned int i = 0; i < integratorInput.size(); ++i)
    if (!integratorInput[i].prototype)
      err += " " + toString(i);

  if (!err.empty())
    THROW("Missing integrator definitions of level(s): " + err + ".");

  // Check if the chain is ok
  for (unsigned int i = 0; i < integratorInput.size(); ++i) {
    const Integrator *prototype = integratorInput[i].prototype;
    if (!prototype) THROW("null prototype");

    if (dynamic_cast<const StandardIntegrator *>(prototype)) {
      if (!((i == 0 && dynamic_cast<const STSIntegrator *>(prototype)) ||
            (i > 0 && dynamic_cast<const MTSIntegrator *>(prototype)))) {

        if (i > 0)
          err += "\nIntegrator " + toString(prototype->getId()) +
            " at level " + toString(i) +
            " is a STS integrator, expected MTS.";
        else
          err += "\nIntegrator " + toString(prototype->getId()) +
            " at level " + toString(i) +
            " is a MTS integrator, expected STS.";
      }

    } else if (dynamic_cast<const NonStandardIntegrator *>(prototype)) {
      err += "\nNonStandardIntegrator (level " + toString(i) + " " +
        toString(prototype->getId()) + ") are not supported by " +
        Force::scope + "Factory yet.";

    } else
      err += "\n[IntegratorFactory::make] Found an integrator '" +
        prototype->getId() +
        "' neither a StandardIntegrator nor NonStandardIntegrator.";
  }

  if (!err.empty()) THROW(err);

  // Now make the integrator chain ... with all forces
  StandardIntegrator *integrator = NULL;
  for (unsigned int i = 0; i < integratorInput.size(); ++i) {
    ForceGroup *forceGroup = new ForceGroup();
    for (unsigned int j = 0; j < integratorInput[i].forces.size(); ++j) {
      Force *force = forceFactory->make(integratorInput[i].forces[j]);
      if (force == NULL) {
        delete forceGroup;
        if (integrator) delete integrator;
        return NULL;
      }
      forceGroup->addForce(force);
    }

    StandardIntegrator *newIntegrator = NULL;
    if (i > 0)
      newIntegrator =
        dynamic_cast<const MTSIntegrator *>(integratorInput[i].prototype)->
          make(integratorInput[i].values, forceGroup, integrator);
    else
      newIntegrator =
        dynamic_cast<const STSIntegrator *>(integratorInput[i].prototype)->
          make(integratorInput[i].values, forceGroup);

    if (newIntegrator == NULL) {
      delete forceGroup;
      if (integrator) delete integrator;
      return NULL;
    }
    integrator = newIntegrator;
  }

  return integrator;
}

void IntegratorFactory::registerHelpText() const {
  for (map<string, const Integrator *,
           ltstrNocase>::const_iterator i = exemplars.begin();
       i != exemplars.end();
       ++i) {
    HelpText helpText;
    i->second->getParameters(helpText.parameters);
    helpText.id = i->second->getIdNoAlias();
    helpText.text = i->second->getText();
    helpText.scope = i->second->getScope();
    HelpTextFactory::registerExemplar(i->second->getId(), helpText);

    HelpText alias;
    alias.text = "alias for \'" + i->second->getId() + "\'";
    alias.scope = i->second->getScope();
    for (map<string, const Integrator *,
             ltstrNocase>::const_iterator j = aliasExemplars.begin();
         j != aliasExemplars.end();
         ++j)
      if (j->second->getIdNoAlias() == i->second->getIdNoAlias()) {
        alias.id = j->first;
        HelpTextFactory::registerExemplar(alias.id, alias);
      }
  }
}

