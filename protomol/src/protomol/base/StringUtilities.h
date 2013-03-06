/*  -*- c++ -*-  */
#ifndef STRINGUTILITIES_H
#define STRINGUTILITIES_H

#include <string>
#include <cstring>

#ifdef HAVE_NO_SSTREAM
#include <protomol/base/sstream_local.h>
#else
#include <sstream>
#endif

#include <protomol/type/Vector3D.h>
#include <vector>
#include <ostream>

namespace ProtoMol {
  //__________________________________________________________________ uppercase
  std::string uppercase(const std::string &word);
  //__________________________________________________________________ lowercase
  std::string lowercase(const std::string &word);


  //_____________________________________________________________________ equal
  bool equal(const std::string &s1, const std::string &s2);
  //________________________________________________________________ equalNocase
  bool equalNocase(const std::string &s1, const std::string &s2);
  //_________________________________________________________________ equalBegin
  bool equalBegin(const std::string &s1, const std::string &s2);
  //___________________________________________________________ equalBeginNocase
  bool equalBeginNocase(const std::string &s1, const std::string &s2);
  //_________________________________________________________________ equalStart
  bool equalStart(const std::string &s1, const std::string &s2);
  //___________________________________________________________ equalStartNocase
  bool equalStartNocase(const std::string &s1, const std::string &s2);
  //___________________________________________________________________ equalEnd
  bool equalEnd(const std::string &s1, const std::string &s2);
  //_____________________________________________________________ equalEndNocase
  bool equalEndNocase(const std::string &s1, const std::string &s2);
  //_____________________________________________________________ equalTerminate
  bool equalTerminate(const std::string &s1, const std::string &s2);
  //_______________________________________________________ equalTerminateNocase
  bool equalTerminateNocase(const std::string &s1, const std::string &s2);



  //_____________________________________________________________________ isReal
  bool isReal(const std::string &word);
  //_____________________________________________________________________ toReal
  bool toReal(const std::string &word, Real &r);
  //_____________________________________________________________________ toReal
  Real toReal(const std::string &word);
  //_____________________________________________________________________ isInt
  bool isInt(const std::string &word);
  //_____________________________________________________________________ toInt
  bool toInt(const std::string &word, int &i);
  //_____________________________________________________________________ toInt
  int toInt(const std::string &word);
  //_____________________________________________________________________ isUInt
  bool isUInt(const std::string &word);
  //_____________________________________________________________________ toUInt
  bool toUInt(const std::string &word, unsigned int &i);
  //_____________________________________________________________________ toUInt
  unsigned int toUInt(const std::string &word);
  //_____________________________________________________________________ isBool
  bool isBool(const std::string &word);
  //_____________________________________________________________________ toBool
  bool toBool(const std::string &word, bool &b);
  //_____________________________________________________________________ toBool
  bool toBool(const std::string &word);
  //_________________________________________________________________ isVector3D
  bool isVector3D(const std::string &word);
  //_________________________________________________________________ toVector3D
  bool toVector3D(const std::string &word, Vector3D &c);
  //_________________________________________________________________ toVector3D
  Vector3D toVector3D(const std::string &word);
  //___________________________________________________________________ isVector
  bool isVector(const std::string &word);
  //___________________________________________________________________ toVector
  std::vector<Real> toVector(const std::string & word);
  //___________________________________________________________________ toVector
  bool toVector(const std::string &word, std::vector<Real> &c);
  //___________________________________________________________________ toString
  std::string toString(Real x);
  //___________________________________________________________________ toString
  std::string toString(Real x, unsigned int n, unsigned int m);
  //___________________________________________________________________ toString
  std::string toString(bool x);
  //___________________________________________________________________ toString
  std::string toString(const Vector3D &x);
  //___________________________________________________________________ toString
  std::string toString(const std::vector<Real> &x);
  //___________________________________________________________________ toString
  /// NB: Needed for symmetry reason!
  inline const std::string &toString(const std::string &x) {
    return x;
  }
  //_______________________________________________________________ toString !NB
  /// NB: Template to catch other types ...
  template<class T>
  inline std::string toString(T x) {
    // http://www.bespecific.com/dialog/becodetalk/archive/980405/0058.html
    std::stringstream ss;
    ss << x;
    return std::string(ss.str());
  }


  //____________________________________________________________________ isBlank
  bool isBlank(const std::string &word);
  //____________________________________________________________________ isblank
  bool isblankchar(char c);
  //________________________________________________________________ isPrintable
  bool isPrintable(const std::string &word);
  //____________________________________________________________ isprintablechar
  bool isprintablechar(char c);

  //___________________________________________________________________ getBegin
  std::string getBegin(const std::string &s, std::string::size_type n);
  //_____________________________________________________________________ getEnd
  std::string getEnd(const std::string &s, std::string::size_type n);

  //_______________________________________________________________ getRightFill
  std::string getRightFill(const std::string &s, std::string::size_type n);
  //________________________________________________________________ getLeftFill
  std::string getLeftFill(const std::string &s, std::string::size_type n);

  //_______________________________________________________ removeBeginEndBlanks
  std::string removeBeginEndBlanks(const std::string &s);

  //________________________________________________________________ ltstrNocase
  struct ltstrNocase {bool operator()(const std::string &s1,
                                      const std::string &s2) const;};
  //______________________________________________________________ ltstrNocaseOp
  bool ltstrNocaseOp(const std::string &s1, const std::string &s2);

  //______________________________________________________________ equalWildcard
  /**
   * Wildcard specifications:
   * - * : matches any string of characters (including none),
   * - \% : matches any single character,
   * - # : matches any string of digits (including none),
   * - + : matches any single digit.
   *
   * Return:
   * - 2 : match without wildcards
   * - 1 : match with wildcards
   * - 0 : no match at all
   */
  int equalWildcard(const std::string &wildcard, const std::string &name);

  //________________________________________________________________ splitString
  std::vector<std::string> splitString(const std::string & str);
  //________________________________________________________________ mergeString
  std::string mergeString(const std::vector<std::string> &str);
  //____________________________________________________________ normalizeString
  std::string normalizeString(const std::string &str);
  //_________________________________________________________________ headString
  std::string headString(const std::string &str);
  //_________________________________________________________________ tailString
  std::string tailString(const std::string &str);

  std::string headerRow(const std::string &title, unsigned int maxColumn = 80);
  void fillFormat(std::ostream &stream, const std::string &str,
                  unsigned int currentColumn = 0, unsigned int indent = 0,
                  unsigned int maxColumn = 80);

  template <typename T>
  std::string Append( const std::string& inData, T value ){
    std::ostringstream retStream;

    retStream << inData << value;

    return std::string( retStream.str() );
  }
}

#endif
