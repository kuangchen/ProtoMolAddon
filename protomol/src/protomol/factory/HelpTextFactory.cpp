#include <protomol/factory/HelpTextFactory.h>
#include <protomol/config/Configuration.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ HelpTextFactory
HelpTextFactory *HelpTextFactory::obj(NULL);

HelpTextFactory::HelpTextFactory() {}

HelpTextFactory::~HelpTextFactory() {}

HelpTextFactory::HelpTextFactory(const HelpTextFactory &) {}

HelpTextFactory &HelpTextFactory::operator=(const HelpTextFactory &) {
  return *this;
}

void HelpTextFactory::kill() {
  HelpTextFactory *p = obj;
  obj = NULL;
  p->~HelpTextFactory();
}

void HelpTextFactory::registerExemplar(const string &id,
                                       const HelpText &helpText) {
  HelpTextFactory::instance().doRegisterExemplar(id, helpText);
}

void HelpTextFactory::registerExemplars(const Configuration *config) {
  if (config != NULL)
    HelpTextFactory::instance().doRegisterExemplars(config);
}

bool HelpTextFactory::unregisterExemplar(const string &id) {
  return HelpTextFactory::instance().doUnregisterExemplar(id);
}

void HelpTextFactory::unregisterAllExemplars() {
  HelpTextFactory::instance().myExemplars.clear();
}

bool HelpTextFactory::empty() {
  return HelpTextFactory::instance().myExemplars.empty();
}

string HelpTextFactory::search(const string &id) {
  return HelpTextFactory::instance().doSearch(id);
}

string HelpTextFactory::keywords() {
  return HelpTextFactory::instance().doKeywords();
}

HelpTextFactory &HelpTextFactory::instance() {
  // We have to do it ourself ... M$ problem ...
  if (obj == NULL) {
    obj = new HelpTextFactory();
    atexit(kill);
  }
  return *obj;
}

void HelpTextFactory::doRegisterExemplar(const string &id,
                                         const HelpText &helpText) {
  myExemplars[id] = helpText;
}

void HelpTextFactory::doRegisterExemplars(const Configuration *config) {
  for (Configuration::const_iterator i = config->begin();
       i != config->end();
       ++i) {
    HelpText helpText;
    helpText.id = i->first;
    helpText.defaultValue = i->second;
    helpText.scope = "Input";
    helpText.text = config->getText(helpText.id);
    doRegisterExemplar(helpText.id, helpText);

    HelpText alias;
    alias.text = "alias for \'" + i->first + "\'";
    alias.scope = helpText.scope;
    vector<string> aliases = config->getAliases(i->first);
    for (unsigned int j = 0; j < aliases.size(); ++j)
      if (!equalNocase(i->first, aliases[j])) {
        alias.id = aliases[j];
        HelpTextFactory::registerExemplar(alias.id, alias);
      }
  }
}

bool HelpTextFactory::doUnregisterExemplar(const string &id) {
  iterator i = myExemplars.find(id);
  if (i != myExemplars.end()) {
    myExemplars.erase(i);
    return true;
  }
  return false;
}

string HelpTextFactory::doSearch(const string &id) const {
  string res;
  const_iterator i = myExemplars.find(id);
  if (i != myExemplars.end()) {
    HelpText txt = i->second;
    res = id;

    if (txt.defaultValue.isDefined() && txt.parameters.empty()) {
      res += " is " + txt.scope + ", ";
      if (txt.defaultValue.getConstraintType() >
          ConstraintValueType::NOCONSTRAINTS)
        res += txt.defaultValue.getConstraintTypeString() + " ";
      res += txt.defaultValue.getTypeString();
      if (txt.defaultValue.valid())
        res += " with default \'" + txt.defaultValue.getString() + "\'";
      res += (txt.text.empty() ? "" : ", ") + txt.text + ".";
    } else {
      res += " is " + txt.scope + (txt.text.empty() ? "" : ", ") + txt.text +
             ".";
      for (unsigned int i = 0; i < txt.parameters.size(); ++i) {
        if (!txt.parameters[i].keyword.empty())
          res += "\n" + std::string(Constant::PRINTINDENT) + getRightFill(
            txt.parameters[i].keyword,
            Constant::PRINTMAXWIDTH);
        else
          res += "\n" + std::string(Constant::PRINTINDENT) + getRightFill("",
            Constant::PRINTMAXWIDTH);
        if (txt.parameters[i].value.getConstraintType() >
            ConstraintValueType::NOCONSTRAINTS)
          res += txt.parameters[i].value.getConstraintTypeString() + " ";
        res += txt.parameters[i].value.getTypeString();
        if (txt.parameters[i].defaultValue.valid())
          res += " with default \'" +
                 txt.parameters[i].defaultValue.getString() + "\'";
        if (!txt.parameters[i].text.empty())
          res += ", " + txt.parameters[i].text;
      }
    }
  } else
    res = "Sorry, could not find any entry for \'" + id + "\'.";
  return res;
}

string HelpTextFactory::doKeywords() const {
  string res;
  for (const_iterator i = myExemplars.begin(); i != myExemplars.end(); ++i) {
    if (i != myExemplars.begin())
      res += "\n";
    res += Constant::PRINTINDENT + getRightFill(i->first,
      Constant::PRINTMAXWIDTH) + " : " + i->second.scope +
           (i->second.text.empty() ? "" : ", " + i->second.text);
  }

  return res;
}

