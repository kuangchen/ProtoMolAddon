/*  -*- c++ -*-  */
#ifndef TYPESELECTION_H
#define TYPESELECTION_H

#include <protomol/type/Real.h>

namespace ProtoMol {
  namespace Private {
    template<bool cmp, class A, class B> struct SelectTypeHelper {
      typedef B type;
    };

    template<class A, class B> struct SelectTypeHelper<true, A, B> {
      typedef A type;
    };

    template<unsigned int size, class A, class B> struct SelectType {
      typedef typename SelectTypeHelper<sizeof(A) == size, A, B>::type type;
    };

    template<bool cmp, class A> struct SelectTypeCheckHelper {};

    template<class A> struct SelectTypeCheckHelper<true, A> {
      typedef A type;
    };

    template<class A, unsigned int size> struct SelectTypeCheck {
      typedef typename SelectTypeCheckHelper<sizeof(A) == size, A>::type type;
    };
  }

  /**
   * Enables to select the right type of an int or float
   * with a given sizeof, bails out if there is not adequate
   * type. @n
   *
   * Usage:@n
   * typedef TypeSelection::Int<4>::type int32;@n
   * typedef TypeSelection::Float<4>::type float32;@n
   */
  namespace TypeSelection {
    /**
     * Select the right type among short, int, long or ong long according the 
     * given sizeof.
     */
    template<unsigned int size> struct Int {
      typedef Private::SelectType<size, int, short> ST1;
      typedef Private::SelectType<size, long, long long> ST2;
      typedef Private::
      SelectType<size, typename ST1::type, typename ST2::type> ST3;
      typedef typename Private::
      SelectTypeCheck<typename ST3::type, size>::type type;
    };

    /**
     * Select the right type among float or double according the given sizeof.
     */
    template<unsigned int size> struct Float {
      typedef Private::SelectType<size, float, double> ST1;
      typedef Private::SelectType<size, typename ST1::type, Real> ST2;
      typedef typename Private::
      SelectTypeCheck<typename ST2::type, size>::type type;
    };
  }
}
#endif /* TYPESELECTION_H */
