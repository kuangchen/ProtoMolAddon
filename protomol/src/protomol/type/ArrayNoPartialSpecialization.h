/*  -*- c++ -*-  */
#ifndef ARRAY_NOPARTIALSPECIALIZATION_H
#define ARRAY_NOPARTIALSPECIALIZATION_H

#include <algorithm>
#include <vector>
#include <assert.h>

namespace ProtoMol {
  //****************************************************************************
  // Class Template for a Generic Resizable N Dimensional Array       (for N>=2)
  // By Giovanni Bavestrelli                 Copyright 1999 Giovanni Bavestrelli
  // Any feedback is welcome, you can contact me at gbavestrelli@yahoo.com
  //
  // This is the version for Microsoft Visual C++ 6.0, and perhaps for other
  // compilers which do not support partial specialization.
  // To make my classes work, I had to use a trick. I moved the RefArray classes
  // as nested classes withing class Array, removing template parameter T from
  // them, allowing me to use full specialization for class RefArray.
  // I don't think this is up to standard C++. If your compiler supports partial
  // specialization, be sure to use the other version of these classes.
  // Thanks to Andrei Alexandrescu for suggesting this trick, which makes my
  // classes work well with Visual C++, although the solution is not standard.
  //****************************************************************************


  //============================================================================
  //Classes for passing a typesafe vector of dimensions to the Array constructor
  //============================================================================

  template<unsigned int N>
  class ArraySize {
    std::vector<unsigned int> &Vector;

    ArraySize<N>(std::vector<unsigned int> &v)
      : Vector(v)
    {}

  public:

    ArraySize < N + 1 > operator()(unsigned int dim) {
      if (Vector.size() > N) Vector.resize(N);
      Vector.push_back(dim);
      return ArraySize < N + 1 > (Vector);
    }
    const std::vector<unsigned int> &Vect() const {
      assert(Vector.size() == N);
      return Vector;
    }

    friend class ArraySizes;
    friend class ArraySize < N - 1 >;
  };

  class ArraySizes {
    std::vector<unsigned int> Vector;

  public:

    explicit ArraySizes(unsigned int dim) {
      Vector.push_back(dim);
    }

    ArraySize<2> operator()(unsigned int dim) {
      if (Vector.size() > 1) Vector.resize(1);
      Vector.push_back(dim);
      return ArraySize<2>(Vector);
    }
  };

  //============================================================================
  // Class Template for a Generic Resizable N Dimensional Array
  //============================================================================

  template<typename T, unsigned int N>
  class Array {
  public:

    //--------------------------------------------------------------------------
    // Class Template for N Dimensional SubArrays within an Array
    //--------------------------------------------------------------------------

    template<unsigned int N>
    class RefArray {
    public:

      // STL-like types
      typedef T value_type;
      typedef T &reference;
      typedef const T &const_reference;
      typedef T *pointer;
      typedef const T *const_pointer;
      typedef T *iterator;
      typedef const T *const_iterator;
      typedef size_t size_type;
      typedef ptrdiff_t difference_type;

      // Give access to number of dimensions
      enum  {array_dims = N};

    private:

      const size_type *const m_pNDimensions; // Array dimensions
      const size_type *const m_pSubArrayLen; // SubArray dimensions

      T   *const m_pElements; // Point to SubArray with elements within Array

      RefArray<N>(T * pElements, const size_type * pNDimensions,
                  const size_type * pSubArrayLen)
        : m_pElements(pElements), m_pNDimensions(pNDimensions), m_pSubArrayLen(
          pSubArrayLen) {
        assert(m_pElements && m_pNDimensions && m_pSubArrayLen);
        assert(m_pNDimensions[0] > 0 && m_pSubArrayLen[0] > 0);
      }

    public:

      RefArray < N - 1 > operator[](size_type Index) {
        assert(m_pElements);
        assert(Index < m_pNDimensions[0]);
        return RefArray < N - 1 > (&m_pElements[Index * m_pSubArrayLen[0]],
                                   m_pNDimensions + 1, m_pSubArrayLen + 1);
      }

      const RefArray < N - 1 > operator[](size_type Index) const {
        assert(m_pElements);
        assert(Index < m_pNDimensions[0]);
        return RefArray < N - 1 > (&m_pElements[Index * m_pSubArrayLen[0]],
                                   m_pNDimensions + 1, m_pSubArrayLen + 1);
      }

      // Return STL-like iterators
      iterator begin()       {return m_pElements;}
      const_iterator begin() const {return m_pElements;}
      iterator end()         {return m_pElements + size();}
      const_iterator end()   const {return m_pElements + size();}

      // Return size of array
      size_type size() const {return m_pNDimensions[0] * m_pSubArrayLen[0];}

      // Return size of subdimensions
      size_type size(unsigned int Dim) const {
        assert(Dim >= 1 && Dim <= N);
        return m_pNDimensions[Dim - 1];
      }

      // Return number of dimensions
      unsigned int dimensions()  const {return N;}

    protected:

      // The following are protected mainly because they are not exception-safe
      // but the way they are used in the rest of the class is exception-safe

      // Copy the elements of another subarray on this one where possible
      // Where not possible, initialize them to a specified value Init
      void copy(const RefArray<N> &SA, const T &Init = T()) {
        size_type below = std::min(size(1), SA.size(1));
        size_type above = size(1);

        // Copy the elements we can copy
        for (size_type i = 0; i < below; ++i)
          (*this)[i].copy(SA[i], Init);

        // Reset the elements we can't copy
        for (size_type j = below; j < above; ++j)
          (*this)[j].initialize(Init);
      }

      // Reset all the elements
      void initialize(const T &Init = T()) {
        std::fill(begin(), end(), Init);
      }

      friend class Array<T, N>;
      friend class Array < T, N + 1 >;
      friend class RefArray < N + 1 >;
    };


    //--------------------------------------------------------------------------
    // Partial Specialization for Monodimensional SubArray within an Array
    //--------------------------------------------------------------------------

    template<>
    class RefArray<1> {
    public:

      // STL-like types
      typedef T value_type;
      typedef T &reference;
      typedef const T &const_reference;
      typedef T *pointer;
      typedef const T *const_pointer;
      typedef T *iterator;
      typedef const T *const_iterator;
      typedef size_t size_type;
      typedef ptrdiff_t difference_type;

      // Give access to number of dimensions
      enum  {array_dims = 1};

    private:

      const size_type *const m_pNDimensions; // Array dimension

      T   *const m_pElements; // Point to elements within Array

      RefArray<1>(T * pElements, const size_type * pNDimensions,
                  const size_type * pSubArrayLen)
        : m_pElements(pElements), m_pNDimensions(pNDimensions) {
        assert(m_pElements && m_pNDimensions && pSubArrayLen);
        // We found the elements
        assert(m_pNDimensions[0] > 0 && pSubArrayLen[0] == 1);
      }

    public:

      reference operator[](size_type Index) {
        assert(m_pElements);
        assert(Index < m_pNDimensions[0]);
        return m_pElements[Index];
      }

      const_reference operator[](size_type Index) const {
        assert(m_pElements);
        assert(Index < m_pNDimensions[0]);
        return m_pElements[Index];
      }

      // Return STL-like iterators
      iterator begin()       {return m_pElements;}
      const_iterator begin() const {return m_pElements;}
      iterator end()         {return m_pElements + size();}
      const_iterator end()   const {return m_pElements + size();}

      // Return size of array
      size_type size()  const {return m_pNDimensions[0];}

      // Return size of subdimensions
      size_type size(unsigned int Dim) const {
        assert(Dim == 1);
        return m_pNDimensions[0];
      }

      // Return number of dimensions
      unsigned int dimensions()  const {return 1;}

    protected:

      // The following are protected mainly because they are not exception-safe
      // but the way they are used in the rest of the class is exception-safe

      // Copy the elements of another subarray on this one where possible
      // Where not possible, initialize them to a specified value Init
      void copy(const RefArray<1> &SA, const T &Init = T()) {
        size_type below = std::min(size(1), SA.size(1));
        size_type above = size(1);

        // Copy the elements we can copy
        for (size_type i = 0; i < below; ++i)
          m_pElements[i] = SA.m_pElements[i];

        // Reset the elements we can't copy
        for (size_type j = below; j < above; ++j)
          m_pElements[j] = Init;
      }

      // Reset all the elements
      void initialize(const T &Init = T()) {
        std::fill(begin(), end(), Init);
      }

      friend class Array<T, 1>;
      friend class Array<T, 2>;
      friend class RefArray<2>;
    };


    //--------------------------------------------------------------------------
    // Class Template for a Generic Resizable N Dimensional Array
    //--------------------------------------------------------------------------

  public:

    // STL-like types
    typedef T value_type;
    typedef T &reference;
    typedef const T &const_reference;
    typedef T *pointer;
    typedef const T *const_pointer;
    typedef T *iterator;
    typedef const T *const_iterator;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    // Give access to number of dimensions
    enum  {array_dims = N};

  private:

    T *m_pArrayElements; // Pointer to actual array elements
    size_type m_nArrayElements; // Total number of array elements

    size_type m_NDimensions[N];  // Size of the N array dimensions
    size_type m_SubArrayLen[N];  // Size of each subarray

  public:

    // Default constructor
    Array<T, N>()
      : m_pArrayElements(NULL), m_nArrayElements(0) {
      std::fill(m_NDimensions, m_NDimensions + N, 0);
      std::fill(m_SubArrayLen, m_SubArrayLen + N, 0);
    }

    // This takes an array of N values representing the size of the N dimensions
    explicit Array<T, N>(const unsigned int *Dimensions, const T &Init = T())
      : m_pArrayElements(NULL), m_nArrayElements(0) {
      std::fill(m_NDimensions, m_NDimensions + N, 0);
      std::fill(m_SubArrayLen, m_SubArrayLen + N, 0);

      resize(Dimensions, Init);
    }

    // This takes an ArraySize object with the N dimensions
    explicit Array<T, N>(const ArraySize<N> &Dimensions, const T &Init = T())
      : m_pArrayElements(NULL), m_nArrayElements(0) {
      std::fill(m_NDimensions, m_NDimensions + N, 0);
      std::fill(m_SubArrayLen, m_SubArrayLen + N, 0);

      resize(Dimensions, Init);
    }

    // Copy constructor
    Array<T, N>(const Array<T, N> &A)
      : m_pArrayElements(NULL), m_nArrayElements(0) {
      std::fill(m_NDimensions, m_NDimensions + N, 0);
      std::fill(m_SubArrayLen, m_SubArrayLen + N, 0);

      Array<T, N> Temp;
      if (!A.empty() && Temp.resize(A.m_NDimensions))
        std::copy(A.begin(), A.end(), Temp.begin());
      swap(Temp);
    }

    // Destructor
    ~Array<T, N>()
    {
      delete[] m_pArrayElements;
    }

    // Indexing Array
    RefArray < N - 1 > operator[](size_type Index) {
      assert(m_pArrayElements);
      assert(Index < m_NDimensions[0]);
      return RefArray < N - 1 > (&m_pArrayElements[Index * m_SubArrayLen[0]],
                                 m_NDimensions + 1, m_SubArrayLen + 1);
    }

    // Indexing Constant Array
    const RefArray < N - 1 > operator[](size_type Index) const {
      assert(m_pArrayElements);
      assert(Index < m_NDimensions[0]);
      return RefArray < N - 1 > (&m_pArrayElements[Index * m_SubArrayLen[0]],
                                 m_NDimensions + 1, m_SubArrayLen + 1);
    }

    // Return RefArray referencing entire Array
    RefArray<N> GetRefArray() {
      assert(m_pArrayElements);
      return RefArray<N>(m_pArrayElements, m_NDimensions, m_SubArrayLen);
    }

    // Return constant RefArray referencing entire Array
    const RefArray<N> GetRefArray() const {
      assert(m_pArrayElements);
      return RefArray<N>(m_pArrayElements, m_NDimensions, m_SubArrayLen);
    }

    // Set the size of each array dimension
    //Visual C++ does not accept parameter defined so: const unsigned int (&)[N]
    // so I accepted a solution which is not type-safe: use it judiciously
    bool resize(const unsigned int *Dimensions,
                const T &Init = T(), bool PreserveElems = false) {
      assert(Dimensions);

      Array<T, N> Temp;

      // Calculate all the information you need to use the array
      Temp.m_nArrayElements = 1;
      for (int i = 0; i < N; ++i) {
        if (Dimensions[i] == 0)
          return false;   // Check that no dimension was zero
        Temp.m_nArrayElements *= Dimensions[i];
        Temp.m_NDimensions[i] = Dimensions[i];
        Temp.m_SubArrayLen[i] = 1;
        for (int k = N - 1; k > i; k--)
          Temp.m_SubArrayLen[i] *= Dimensions[k];
      }

      // Allocate new elements, let exception propagate
      Temp.m_pArrayElements = new T[Temp.m_nArrayElements];

      // Some compilers might not throw exception if allocation fails
      if (!Temp.m_pArrayElements)
        return false;

      // Copy the elements from the previous array if requested
      if (PreserveElems && !empty())
        Temp.copy(*this, Init);
      // Otherwise initialize them to the specified value
      else
        Temp.initialize(Init);

      // Now swap this object with the temporary
      swap(Temp);

      return true;
    }

    // resize accepting a fixed ArraySize, this solution is type-safe
    bool resize(const ArraySize<N> &Dimensions,
                const T &Init = T(), bool PreserveElems = false) {
      unsigned int Dims[N];
      std::copy(Dimensions.Vect().begin(), Dimensions.Vect().end(), Dims);
      return resize(Dims, Init, PreserveElems);
    }

    // Delete the complete Array
    void clear() {
      delete[] m_pArrayElements;
      m_pArrayElements = NULL;
      m_nArrayElements = 0;

      std::fill(m_NDimensions, m_NDimensions + N, 0);
      std::fill(m_SubArrayLen, m_SubArrayLen + N, 0);
    }

    // Assignment operator
    Array<T, N> &operator=(const Array<T, N> &A) {
      if (&A != this) { // For efficiency
        Array<T, N> Temp(A);
        swap(Temp);
      }
      return *this;
    }

    // Return STL-like iterators
    iterator begin()       {return m_pArrayElements;}
    const_iterator begin() const {return m_pArrayElements;}
    iterator end()         {return m_pArrayElements + m_nArrayElements;}
    const_iterator end()   const {return m_pArrayElements + m_nArrayElements;}

    // Some more STL-like size members
    size_type size()       const {return m_nArrayElements;}

    // Return the size of each dimension, 1 to N
    size_type size(unsigned int Dim) const {
      assert(Dim >= 1 && Dim <= N);
      return m_NDimensions[Dim - 1];
    }

    // Say if the array is empty
    bool empty()           const {return m_nArrayElements == 0;}

    // Return number of dimensions
    unsigned int dimensions()  const {return N;}

    // Swap this array with another, a'la STL
    void swap(Array<T, N> &A) {
      std::swap(m_pArrayElements, A.m_pArrayElements);
      std::swap(m_nArrayElements, A.m_nArrayElements);

      std::swap_ranges(m_NDimensions, m_NDimensions + N, A.m_NDimensions);
      std::swap_ranges(m_SubArrayLen, m_SubArrayLen + N, A.m_SubArrayLen);
    }

  protected:

    // The following are protected mainly because they are not exception-safe
    // but the way they are used in the rest of the class is exception-safe

    // Copy the elements of another array on this one where possible
    // Where not possible, initialize them to a specified value Init
    void copy(const Array<T, N> &A, const T &Init = T()) {
      size_type below = std::min(size(1), A.size(1));
      size_type above = size(1);

      // Copy the elements we can copy
      for (size_type i = 0; i < below; ++i)
        (*this)[i].copy(A[i], Init);

      // Reset the elements we can't copy
      for (size_type j = below; j < above; ++j)
        (*this)[j].initialize(Init);
    }

    // Initialize all the array elements
    void initialize(const T &Init = T()) {
      std::fill(begin(), end(), Init);
    }

    // Prefer non-member operator ==, but it needs to be a friend
    friend bool operator==(const Array<T, N> &A, const Array<T, N> &B);
  };


  // Test for equality between two arrays
  template<typename T, unsigned int N>
  bool operator==(const Array<T, N> &A, const Array<T, N> &B) {
    return std::equal(A.m_NDimensions, A.m_NDimensions + N, B.m_NDimensions) &&
           std::equal(A.begin(), A.end(), B.begin());
  }

  // Test for inequality between two arrays
  template<typename T, unsigned int N>
  bool operator!=(const Array<T, N> &A, const Array<T, N> &B) {
    return !(A == B);
  }


  /* The following don't work for Visual C++

     // Not implemented, meaningless to have 0 dimensions
     template <typename T>
     class Array<T,0>
     {
     };

     // Not implemented, use std::vector for one dimensional arrays
     template <typename T>
     class Array<T,1>
     {
     };

   */
}
#endif


