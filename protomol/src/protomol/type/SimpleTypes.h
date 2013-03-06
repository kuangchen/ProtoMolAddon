/*  -*- c++ -*-  */
#ifndef SIMPLETYPES_H
#define SIMPLETYPES_H

#include <map>
#include <string>

namespace ProtoMol {
  /// Pair of std::string
  typedef std::pair<std::string, std::string> PairString;

  /// Pair of int
  typedef std::pair<int, int> PairInt;

  /// Pair of unsigned int
  typedef std::pair<unsigned int, unsigned int> PairUInt;

  /// Pair of sorted int, where first <= second
  struct PairIntSorted {
    PairIntSorted() : first(0), second(0) {}
    PairIntSorted(unsigned int a, unsigned int b) :
      first(std::min(a, b)), second(std::max(a, b)) {}
    bool operator<(const PairIntSorted &p) const {
      if (first < p.first)
        return true;
      else if (first > p.first)
        return false;
      else if (second < p.second)
        return true;
      return false;
    }
    bool operator==(const PairIntSorted &p) const {
      return first == p.first && second == p.second;
    }

    unsigned int first, second;
  };

  /// Triple of int
  struct TripleInt {
    TripleInt() : h(0), k(0), l(0) {}
    TripleInt(int a, int b, int c) : h(a), k(b), l(c) {}
    int h; int k; int l;
  };

  /// Place holder for any text string, avoiding implicit conversion
  /// and distinguishing a text from and any other string
  struct Text {
    Text() : text("") {}
    Text(const std::string &t) : text(t) {}
    Text(const char *t) : text(t) {}

    std::string text;
  };
}
#endif /* TYPESELECTION_H */
