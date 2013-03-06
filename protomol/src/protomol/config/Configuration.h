/*  -*- c++ -*-  */
#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <protomol/config/Parameter.h>
#include <protomol/base/StringUtilities.h>
#include <ostream>
#include <vector>
#include <map>

namespace ProtoMol {
  
  //________________________________________________________ Configuration
  /**
   * Container (map) holding a system configuration, a keyword together
   * with the value and an optional help text.  The value uses Value, which
   * comes with additional information about the type and constraints.
   */
  class Configuration {
    typedef std::map<std::string, Value, ltstrNocase>       ValueMapType;
    typedef std::map<std::string, std::string, ltstrNocase> AliasMapType;
    typedef std::map<std::string, std::string, ltstrNocase> TextMapType;
    typedef ValueMapType::iterator iterator;
  public:
    typedef ValueMapType::const_iterator const_iterator;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Configuration
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// Register new keyword with its value
    void registerKeyword(const std::string &keyword, Value value);
    /// Register aliases for a given keyword
    void registerAliases(const std::string &keyword,
                         const std::vector<std::string> &aliases);
    /// Unregister a keyword
    void unregisterKeyword(const std::string &keyword);

    /// Test if the keyword exists
    bool empty(const std::string &keyword = std::string("")) const;
    /// Test if the keyword has a defined type
    bool defined(const std::string &keyword) const;
    /// Test if the keyword's value is valid
    bool valid(const std::string &keyword) const;

    // Set
    bool set(const std::string &keyword, Value);
    // Uses ValueTraits<>::init() if val == ""
    bool set(const std::string &keyword, const std::string &val);
    bool set(const std::string &keyword,
             const std::vector<std::vector<std::string> > &values);
    bool set(const std::vector<std::vector<std::string> > &values);
    bool setText(const std::string &keyword, const std::string &text);

    // Get
    Value get(const std::string &keyword) const;
    std::vector<Value> get(const std::vector<Parameter> &parameters) const;
    std::string getText(const std::string &keyword) const;
    std::vector<std::string> getAliases(const std::string & keyword) const;

    // Get and set; direct access, unsafe
    Value &operator[](const std::string &keyword);
    const Value &operator[](const std::string &keyword) const;

    // Global print & test
    std::ostream &print(std::ostream &stream) const;
    friend
    std::ostream &operator<<(std::ostream &stream, const Configuration &c);
    bool hasUndefinedKeywords() const;
    std::string printUndefinedKeywords() const;

    // Iterators
    const_iterator begin() const {return myValues.begin();}
    const_iterator end() const {return myValues.end();}
    const_iterator find(const std::string &keyword) const;

  private:
    iterator begin() {return myValues.begin();}
    iterator end() {return myValues.end();}
    iterator find(const std::string &keyword);

  private:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // private data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    ValueMapType myValues;
    AliasMapType myAliases;
    TextMapType myTexts;
  };

  inline std::ostream &
  operator<<(std::ostream &stream, const Configuration &c) {
    return c.print(stream);
  }
}
#endif
