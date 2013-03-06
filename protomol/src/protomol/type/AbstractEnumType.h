/*  -*- c++ -*-  */
#ifndef ABSTRACTENUMTYPE_H
#define ABSTRACTENUMTYPE_H

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>


namespace ProtoMol {
  //________________________________________________________ AbstractEnumType
  /**
   * Provides a simple interface to implement a simple map between enum's
   * and string's. AbstractEnumType supports conversions and comparison
   * between string and
   * enum/int simplifying the parsing. It is array based (search is of O(N)).
   */

  template<typename T, bool noCase = true>
  class AbstractEnumType : public T {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
    typedef typename T::Enum Enum;
public:
    AbstractEnumType() : myType(T::UNDEFINED) {}
    AbstractEnumType(Enum n) : myType(n) {}
    AbstractEnumType(const std::string &s) : myType(getEnum(s)) {}
    AbstractEnumType(const char s[]) : myType(getEnum(std::string(s))) {}
    //virtual ~AbstractEnumType(){};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class AbstractEnumType
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
    std::string getString() const {return T::str[static_cast<int>(myType)];}
    Enum getEnum() const {return myType;}
    bool valid() const {return myType > T::FIRST && myType < T::LAST;}
    static bool valid(Enum e) {return e > T::FIRST && e < T::LAST;}
    static const std::string &getString(Enum e) {return T::str[static_cast<int>(
                                                                 e)];}
    static std::string getPossibleValues(std::string separator = " | ");
public:
    // Conversion operators
    operator Enum() const {return myType;}
    operator std::string() const {return getString(myType);}
    AbstractEnumType &operator=(Enum a) {myType = a; return *this;}
    // Comparison operators
    bool operator==(Enum a)  {return a == myType;}
    bool operator==(const AbstractEnumType &a)  {return a.myType == myType;}
    bool operator!=(Enum a)  {return a != myType;}
    bool operator!=(const AbstractEnumType &a)  {return a.myType != myType;}
private:
    // Forbidden
    // operator int() const;

private:
    static Enum getEnum(const std::string &s);
    static Enum getEnum(int n);
    static const std::string &getString(int n) {return T::str[static_cast<int>(
                                                                getEnum(n))];}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // friends of AbstractEnumType
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
    friend Report::MyStreamer &operator<<(Report::MyStreamer &OS,
                                          const AbstractEnumType &b) {return OS
     << b.getString();}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // private data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
private:
    Enum myType;
  };

  template<typename T, bool noCase>
  typename AbstractEnumType<T, noCase>::Enum AbstractEnumType<T,
                                                              noCase>::getEnum(
    const std::string &s) {
    for (int i = static_cast<int>(T::FIRST); i < static_cast<int>(T::LAST);
         ++i) {
      if (!noCase)
        if (equal(T::str[i], s))
          return static_cast<typename AbstractEnumType<T, noCase>::Enum>(i);
      if (noCase)
        if (equalNocase(T::str[i], s))
          return static_cast<typename AbstractEnumType<T, noCase>::Enum>(i);
    }

    return T::UNDEFINED;
  }

  template<typename T, bool noCase>
  typename AbstractEnumType<T, noCase>::Enum AbstractEnumType<T,
                                                              noCase>::getEnum(
    int n) {
    if (n < static_cast<int>(T::FIRST) || n >= static_cast<int>(T::LAST))
      return T::UNDEFINED;
    return static_cast<typename AbstractEnumType<T, noCase>::Enum>(n);
  }

  template<typename T, bool noCase>
  std::string AbstractEnumType<T, noCase>::getPossibleValues(
    std::string separator) {
    std::string tmp = "";
    for (int i = static_cast<int>(T::FIRST) + 1;
         i < static_cast<int>(T::LAST);
         ++i) {
      tmp += T::str[static_cast<int>(i)];
      if (i < static_cast<int>(T::LAST) - 1)
        tmp += separator;
    }

    return tmp;
  }
}
#endif //  ABSTRACTENUMTYPE_H
