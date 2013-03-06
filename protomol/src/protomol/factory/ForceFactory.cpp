#include <protomol/factory/ForceFactory.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/type/SimpleTypes.h>
#include <protomol/force/CompareForce.h>
#include <protomol/force/TimeForce.h>
#include <protomol/parallel/Parallel.h>
#include <algorithm>
#include <sstream>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

//____ ForceFactory
ForceFactory::ForceFactory() :
  lastCompareForce(NULL) {}

ForceFactory::~ForceFactory() {}

Force *ForceFactory::make(const string &idInput, vector<Value> values) const {
  // Just try .. and see if we can find it directly
  string id = normalizeString(idInput);
  const Force *prototype = getPrototype(id);
  
  // Check that paramter size matches values
  if (prototype != NULL) {
      vector<Parameter> tmp;
      prototype->getParameters(tmp);
      if (values.size() != tmp.size()) prototype = NULL;
  }

  if (prototype == NULL) {
    // Lookup tables
    if (!cache) updateCache();

    vector<string> splitId(splitForceString(id));

    // Find the force keyword(s)
    string keyword = splitId[0];
    forceTypes_t::const_iterator itr = forceTypes.find(keyword);

    if (itr == forceTypes.end()) {
      // sort force keyword(s) and try again
      vector<string> tmp(splitString(keyword));
      sort(tmp.begin(), tmp.end(), ltstrNocaseOp);
      keyword = mergeString(tmp);
      if (forceTypesSorted.find(keyword) != forceTypesSorted.end()) {
        keyword = forceTypesSorted[keyword];
        itr = forceTypes.find(keyword);
      }
    }

    string newId;
    if (itr != forceTypes.end()) {
      // retrieve all valid policies for the force keyword(s) from the string
      newId = keyword;
      for (unsigned int i = 1; i < splitId.size(); ++i)
        if (forceTypes[keyword].policy.find(splitId[i]) !=
            forceTypes[keyword].policy.end())
          newId += " " + (*forceTypes[keyword].policy.find(splitId[i]));

      prototype = getPrototype(newId);
    }

    if (prototype == NULL && itr != forceTypes.end()) {
      // ... sort policies and try again ...
      newId = sortForceString(newId);
      for (map<string, string, ltstrNocase>::const_iterator j =
             itr->second.policiesSorted.begin();
           j != itr->second.policiesSorted.end(); ++j)
        if (equalNocase(j->first, newId))
          newId = j->second;

      prototype = getPrototype(newId);
    }

    if (prototype == NULL && itr != forceTypes.end()) {
      // ... remove double policies and try again ...
      newId = uniqueForceString(newId);
      for (forceTypesSorted_t::const_iterator j =
             itr->second.policiesSorted.begin();
           j != itr->second.policiesSorted.end(); ++j)
        if (equalNocase(uniqueForceString(j->first), newId))
          newId = j->second;

      prototype = getPrototype(newId);
    }

    //  ... try with alias!
    if (prototype == NULL) {
      prototype = getPrototype(splitId[0]);
      if (prototype != NULL)
        newId = prototype->getId();
    }

    // ... ok final last call SK945 ... time and compare
    if (prototype == NULL && equalStartNocase(CompareForce::keyword, id)) {
      Force *actualForce =
        make(id.substr(CompareForce::keyword.size() + 1), values);
      if (actualForce == NULL)
        return NULL;
      CompareForce *compareForce = NULL;
      if (lastCompareForce == NULL) {
        compareForce = actualForce->makeCompareForce(actualForce, NULL);
        lastCompareForce = compareForce;
      } else {
        compareForce =
          actualForce->makeCompareForce(actualForce, lastCompareForce);
        report << hint << "Comparing " << compareForce->getIdNumber() / 2 <<
        " :\'"
               << lastCompareForce->getForceObject()->getId()
               << "\' with \'" << actualForce->getId() << "\'" << endr;
        lastCompareForce = NULL;

        if (Parallel::isParallel())
          report << hint << "Comparing forces in parallel environment "
                 << "not supported." << endr;
      }
      return compareForce;
    } else if (prototype == NULL &&
               equalStartNocase(TimeForce::keyword, idInput)) {
      Force *actualForce =
        make(id.substr(TimeForce::keyword.size() + 1), values);
      if (actualForce == NULL)
        return NULL;
      TimeForce *timeForce = actualForce->makeTimeForce(actualForce);
      report << hint << "Timing " << timeForce->getIdNumber() << " :\'"
             << actualForce->getId() << "\'" << endr;
      return timeForce;
    }

    string parametersString;
    if (prototype != NULL) {
      // ok, we found a match, retrieve the parameters to
      // be passed the make method ... now remove
      // all force keyword(s) and policies.
      vector<string> splitNewId(splitForceString(newId));
      set<string, ltstrNocase> unique;
      for (unsigned int i = 0; i < splitNewId.size(); ++i)
        unique.insert(splitNewId[i]);

      for (unsigned int i = 0; i < splitId.size(); ++i)
        if (unique.find(splitId[i]) == unique.end())
          parametersString += (parametersString.empty() ? "" : " ") +
                              splitId[i];

      parametersString += " ";  // some compiler need this

      vector<Parameter> parameters;
      prototype->getParameters(parameters);
      values.resize(parameters.size());

      stringstream ssp(parametersString);
      bool foundLast = false;

      for (unsigned int i = 0; i < parameters.size(); ++i) {
        // parse the the parameters one by one
        string strp;
        values[i] = parameters[i].value;
        values[i].clear();

        bool found = false;
        bool retry = true;
        if (!(ssp)) {
          ssp.str(parametersString);
          ssp.clear();
          ssp.seekg(ios::beg);
        }
        while (ssp >> strp || retry) {
          if (!(ssp) && retry) {
            ssp.str(parametersString);
            ssp.clear();
            ssp.seekg(ios::beg);
            retry = false;
            strp = "";
            continue;
          }
          if (equalNocase(strp, parameters[i].keyword) &&
              !parameters[i].keyword.empty()) {
            ssp >> values[i];
            found = true;
            break;
          } else if (foundLast && parameters[i].keyword.empty()) {
            ssp.seekg((-1) * static_cast<int>(strp.size()), ios::cur);
            ssp.clear();
            ssp >> values[i];
            found = true;
            break;
          }
        }

        foundLast = found;

        // If still undefined take default value is available
        if (!found && parameters[i].defaultValue.valid())
          values[i].set(parameters[i].defaultValue);
      }

      id = newId;
    } else if (itr != forceTypes.end()) {
      string err;

      err += " Could not find a complete match for \'" + idInput +
        "\' in " + Force::scope + "Factory of force \'" + keyword +
        "\'. Possible definitions are:";

      for (policy_t::const_iterator j = itr->second.policies.begin();
           j != itr->second.policies.end(); ++j)
        err += string("\n") + itr->first + string(" ") + (*j);

      THROW(err);

    } else {
      string err;

      err += " Could not find any match for \'" + idInput + "\' in " +
        Force::scope + "Factory. Possible forces are:";

      for (forceTypes_t::const_iterator i =
             forceTypes.begin(); i != forceTypes.end(); ++i)
        err += string("\n") + i->first;

      THROW(err);
    }
  }

  // Make
  Force *newObj = NULL;

  if (prototype != NULL) {
    // Check parameter list
    prototype->assertParameters(values);

    newObj = prototype->make(values);
    if (newObj != NULL) newObj->setAlias(id);
  }

  if (!newObj) THROW("Could not make force");

  return newObj;
}

ostream &ForceFactory::print(ostream &stream) const {
  if (!cache) updateCache();

  forceTypes_t::const_iterator i;
  for (i = forceTypes.begin(); i != forceTypes.end(); ++i) {
    stream << i->first << endl;

    policy_t::const_iterator j;
    for (j = i->second.policies.begin(); j != i->second.policies.end(); ++j) {
      if (!(*j).empty())
        stream << Constant::PRINTINDENT << (*j) << endl;

      string protoName = i->first + ((*j).empty() ? "" : " " + (*j));
      vector<Parameter> parameters;
      getPrototype(protoName)->getParameters(parameters);
      for (unsigned int k = 0; k < parameters.size(); k++)
        parameters[k].print(stream);
    }
  }

  stream << "Aliases:" << endl;
  printAliases(stream);

  return stream;
}

void ForceFactory::updateCache() const {
  forceTypes.clear();
  forceTypesSorted.clear();

  set<PairString> removes;

  for (pointers_t::const_iterator i = begin(); i != end(); ++i) {
    string id = (*i)->getId();
    string idSorted = sortForceString(id);
    vector<string> idSplit = splitForceString(id);
    string keyword = idSplit[0];

    vector<string> tmp(splitString(keyword));
    sort(tmp.begin(), tmp.end(), ltstrNocaseOp);
    string keywordSorted = mergeString(tmp);

    for (unsigned int j = 1; j < idSplit.size(); ++j)
      forceTypes[keyword].policy.insert(idSplit[j]);

    if (!equalNocase(keywordSorted, keyword))
      forceTypesSorted[keywordSorted] = keyword;

    tmp.resize(idSplit.size() - 1);
    copy(idSplit.begin() + 1, idSplit.end(), tmp.begin());
    forceTypes[keyword].policies.insert(mergeString(tmp));

    if (forceTypes[keyword].policiesSorted.find(idSorted) !=
        forceTypes[keyword].policiesSorted.end() &&
        !equalNocase(forceTypes[keyword].policiesSorted[idSorted], id))
      removes.insert(PairString(keyword, idSorted));
    else
      forceTypes[keyword].policiesSorted[idSorted] = id;
  }

  // Remove all sorted id entries which match when sorted
  for (set<PairString>::const_iterator i = removes.begin();
       i != removes.end();
       ++i)
    forceTypes[i->first].
      policiesSorted.erase(forceTypes[i->first].
        policiesSorted.find(i->second));

  cache = true;
}

vector<string> ForceFactory::splitForceString(const string &id) const {
  vector<string> res(1);
  vector<string> tmp(splitString(id));
  unsigned int i = 0;
  while (i < tmp.size() && tmp[i][0] != '-')
    res[0] += (res[0].empty() ? "" : " ") + tmp[i++];

  for (; i < tmp.size(); ++i)
    if (tmp[i][0] == '-' && tmp[i].size() > 1 && !isdigit(tmp[i][1]))
      res.push_back(tmp[i]);
    else
      res[res.size() - 1] += " " + tmp[i];

  return res;
}

vector<string> ForceFactory::splitForceStringSorted(const string &id) const {
  vector<string> res(splitForceString(id));

  vector<string> tmp(splitString(res[0]));
  sort(tmp.begin(), tmp.end(), ltstrNocaseOp);
  res[0] = mergeString(tmp);

  sort(res.begin() + 1, res.end(), ltstrNocaseOp);
  return res;
}

string ForceFactory::sortForceString(const string &id) const {
  return mergeString(splitForceStringSorted(id));
}

string ForceFactory::uniqueForceString(const string &id) const {
  vector<string> tmp(splitForceString(id));
  string res = tmp[0];
  set<string, ltstrNocase> unique;
  for (unsigned int i = 1; i < tmp.size(); ++i) {
    if (unique.find(tmp[i]) == unique.end())
      res += " " + tmp[i];
    unique.insert(tmp[i]);
  }

  return res;
}

