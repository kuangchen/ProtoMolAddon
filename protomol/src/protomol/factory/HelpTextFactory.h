/* -*- c++ -*- */
#ifndef HELPTEXTFACTORY_H
#define HELPTEXTFACTORY_H

#include <protomol/config/Parameter.h>
#include <vector>
#include <map>
#include <set>

namespace ProtoMol {
  class Configuration;

  //________________________________________ HelpTextFactory
  /**
   */

  struct HelpText {
    std::string id;
    Value defaultValue;
    std::string text;
    std::string scope;
    std::vector<Parameter> parameters;
  };

  class HelpTextFactory {
  private:
    typedef std::map<std::string, HelpText, ltstrNocase> HelpTextMapType;
    typedef HelpTextMapType::const_iterator const_iterator;
    typedef HelpTextMapType::iterator iterator;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    HelpTextFactory();
    ~HelpTextFactory();
    HelpTextFactory(const HelpTextFactory &);
    HelpTextFactory &operator=(const HelpTextFactory &);
  private:
    /// Call by atexit() to clean up.
    static void kill();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class HelpTextFactory
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    static void registerExemplar(const std::string &id,
                                 const HelpText &helpText);
    static void registerExemplars(const Configuration *config);
    static bool unregisterExemplar(const std::string &id);
    static void unregisterAllExemplars();
    static bool empty();
    static std::string search(const std::string &id);
    static std::string keywords();


  private:
    static HelpTextFactory &instance();
    void doRegisterExemplar(const std::string &id, const HelpText &helpText);
    void doRegisterExemplars(const Configuration *config);
    bool doUnregisterExemplar(const std::string &id);
    std::string doSearch(const std::string &id) const;
    std::string doKeywords() const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // private data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    static HelpTextFactory *obj;
    HelpTextMapType myExemplars;
  };
}
#endif /* HELPTEXTFACTORY_H */
