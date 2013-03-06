/*  -*- c++ -*-  */
#ifndef VECTOR_H
#define VECTOR_H

#include <vector>

namespace ProtoMol {
  template<typename T>
  class Vector;
  //________________________________________________________ VectorPacker
  template<typename T, unsigned int N>
  class VectorPacker {
    std::vector<T> &myVector;

    VectorPacker<T, N>(std::vector<T> &v) : myVector(v) {}

    VectorPacker < T, N + 1 > operator()(T elm) {
      if (myVector.size() > N)
        myVector.resize(N);
      myVector.push_back(elm);
      return VectorPacker < T, N + 1 > (myVector);
    }
  public:

    operator std::vector < T> () {return myVector;};
    friend class VectorPacker < T, N - 1 >;
    friend class Vector<T>;
  };

  //________________________________________________________ Vector
  /**
   * Vector enables to initialize a STL std::vector with arbitrary number of
   * elements@n
   *
   * static const std::vector<std::string> @n
   * months(Vector<std::string>("Jan")("Feb")("Mar")("Apr")("May")
   *        ("Jun")("Jul")("Aug")("Sep")("oct")("Nov")("Dec"));
   */
  template<typename T>
  class Vector  {
    std::vector<T> myVector;

  public:
    Vector(void) {}
    explicit Vector(T elm) {
      myVector.push_back(elm);
    }

    operator std::vector < T> () {return myVector;};

    VectorPacker<T, 2> operator()(T elm) {
      if (myVector.size() > 1)
        myVector.resize(1);
      myVector.push_back(elm);
      return VectorPacker<T, 2>(myVector);
    }
  };
}
#endif
