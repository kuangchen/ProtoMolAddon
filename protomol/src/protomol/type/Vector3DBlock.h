/*  -*- c++ -*-  */
#ifndef VECTOR3DBLOCK_H
#define VECTOR3DBLOCK_H

#include <vector>
#include <iostream>
using namespace std;
#include <protomol/type/Vector3D.h>
#include <protomol/base/Proxy.h>

namespace ProtoMol {
  //_____________________________________________________________ Vector3DBlock
  /**
   * Container holding a vector (array) of 3D coordinates/vectors
   */
  class Vector3DBlock : public Proxy {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types & enum's
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    enum {LIMIT = 30};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /* No change */
    Vector3DBlock()  : Proxy(), vec() {}

    explicit Vector3DBlock(unsigned int n) : Proxy() {
      c = new Real[3*n];
      for (unsigned int i = 0; i < n; i++)
        vec.push_back(Vector3DB(c+3*i));
    }

    Vector3DBlock(unsigned int n, const Vector3D &t) : Proxy() {
      c = new Real[3*n];
      for (unsigned int i = 0; i < n; i++)
        vec.push_back(Vector3DB(t, c+3*i));
    }

    Vector3DBlock(const Vector3DBlock &rhs) {
      this->intoAssign(rhs);
    }

    //Vector3DBlock(const std::vector<>) {

    Vector3DBlock &operator=(const Vector3DBlock &rhs) {
      this->intoAssign(rhs);
      return *this;
    }
    // NOTE: A destructor is not needed here, because when vec is deleted
    // each Vector3D will be deallocated?
    ~Vector3DBlock() {
      if (size() != 0) delete [] c;
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  New methods of class Vector3DBlock from std::vector
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    unsigned int size() const {return vec.size();}
    bool empty() const {return vec.empty();}
    Vector3DB& operator[](unsigned int n) {return vec[n];}
    const Vector3DB& operator[](unsigned int n) const {return vec[n];}
    void push_back(const Vector3D &t) {
      resize(size()+1, t);
    }
    //void pop_back() {return vec.pop_back();}
    void swap(Vector3DBlock &x) {
      vec.swap(x.vec);
      Real* tmp = c;
      c = x.c;
      x.c = tmp;
    }

    void clear() {vec.clear();}
    void resize(unsigned int n, Vector3D t = Vector3D(0.0, 0.0, 0.0)) {
      // Three cases.
      // Case 1: n is the same as the current size.
      if (n == size()) {
        // Do nothing.
        return;
      }
      // Case 2: n is less than the current size (we're removing data)
      else if (n < size()) {
        // 1. Create a smaller array.
        Real* newdata = new Real[n*3];
        // 2. Copy data from c
        for (unsigned int i = 0; i < 3*n; i++) {
          newdata[i] = c[i];
          // 2b. Also set internal vector3d pointers to new data
          if (i % 3 == 0)
            vec[i].c = newdata+3*i;
        }
        // 3. Delete old c
        delete c;
        // 4. Reset c to the new data
        c = newdata;
        // 5. Resize STL vector (will just remove entries at the end)
        vec.resize(n);
      }
      // Case 3: n is greater than the current size (must add data)
      else if (n > size()) {
        // 1. Create a bigger array
        Real* newdata = new Real[n*3];
        // 2. Copy data from c
        int currentsize = size();
        for (int i = 0; i < currentsize*3; i++) {
          newdata[i] = c[i];
          // 2b. Also set internal vector3d pointers to new data
          if (i % 3 == 0) {
            vec[i].c = newdata+3*i;
          }
        }
        // 3. Resize STL vector (will just add entries to the end)
        // 5. Delete old c
        if (currentsize != 0) delete c;
        // 6. Reset c to the new data
        c = newdata;
        // 4. Point internal vector3d pointers to new data
        for (unsigned int i = currentsize; i < n; i++) {
          vec.push_back(Vector3DB(t.c[0], t.c[1], t.c[2], c+3*i));
        }
      }
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Vector3DBlock
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// Clear (set to zero) each element of the vector.
    void zero(int n = -1) {
      for (unsigned int i = 0; i < size()*3; i++)
        c[i] = 0;
      if (n >= 0 && n != (int)size())
        resize(n, Vector3D(0,0,0));
    }


    Vector3DBlock &intoAssign(const Vector3DBlock &x) {
      if (size() != x.size()) resize(x.size());
      for (unsigned int i = 0; i < 3*size(); i++)
        c[i] = x.c[i];
      return *this;
    }

    /// Add another Vector3DBlock into this one.
    Vector3DBlock &intoAdd(const Vector3DBlock &x) {
      const unsigned int count = 3*size();
      for (unsigned int i = 0; i < count; ++i)
        c[i] += x.c[i];
      return *this;
    }

    Vector3DBlock& operator+=(const Vector3DBlock &x) {
      return this->intoAdd(x);
    }

    /// Subtract another Vector3DBlock from this one.
    Vector3DBlock &intoSubtract(const Vector3DBlock &x) {
      const unsigned int count = 3*size();

      for (unsigned int i = 0; i < count; ++i)
        c[i] -= x.c[i];
      return *this;
    }

    Vector3DBlock& operator-=(const Vector3DBlock &x) {
      return this->intoSubtract(x);
    }

    /// Add a scalar multiple of another Vector3DBlock into this one.
    Vector3DBlock &intoWeightedAdd(Real weight, const Vector3DBlock &x) {
      const unsigned int count = 3*size();
      for (unsigned int i = 0; i < count; ++i)
        c[i] += x.c[i] * weight;
      return *this;
    }

    /// Assign a scalar multiple of another Vector3DBlock into this one.
    Vector3DBlock &intoWeighted(Real weight, const Vector3DBlock &x) {
      const unsigned int count = 3*size();
      for (unsigned int i = 0; i < count; ++i)
        c[i] = x.c[i] * weight;
      return *this;
    }

    /// Subtract a scalar multiple of another Vector3DBlock from this one.
    Vector3DBlock &intoWeightedSubtract(Real weight, const Vector3DBlock &x) {
      const unsigned int count = 3*size();

      for (unsigned int i = 0; i < count; ++i)
        c[i] -= x.c[i] * weight;
      return *this;
    }

    Vector3DBlock operator+(const Vector3DBlock& rhs) {
      Vector3DBlock retval(size());

      for (unsigned int i = 0; i < 3*size(); ++i)
        retval.c[i] = c[i] + rhs.c[i];
      return retval;
    }

    Vector3DBlock operator-(const Vector3DBlock& rhs) {
      Vector3DBlock retval(size());
      for (unsigned int i = 0; i < 3*size(); ++i)
        retval.c[i] = c[i] - rhs.c[i];
      return retval;
    }

    Vector3DBlock operator*(Real scale) {
      Vector3DBlock retval(size());
      for (unsigned int i = 0; i < 3*size(); ++i)
        retval.c[i] = c[i] * scale;
      return retval;
    }

    Vector3DBlock operator/(Real scale) {
      Vector3DBlock retval(size());
      for (unsigned int i = 0; i < 3*size(); ++i)
        retval.c[i] = c[i] / scale;
      return retval;
    }

    /// Compute the boinding box over all elements
    void boundingbox(Vector3D &minbb, Vector3D &maxbb) const;

    /// Compute the sum over all elements
    Vector3D sum() const;

    /// Compute regression plane by SVD

    /// Streams
    friend std::ostream &operator<<(std::ostream &OS, const Vector3DBlock &vblock) {
      unsigned int blkSz = vblock.size();

      OS << blkSz;

      for (unsigned int i=0; i< blkSz; i++)
        OS << std::endl << vblock[i];

      return OS;
    }

    friend std::istream &operator>>(std::istream &OS, Vector3DBlock &vblock) {
      Vector3D coords;
      unsigned int blkSz;

      OS >> blkSz;

      vblock.resize( blkSz );
      for (unsigned int i=0; i< blkSz; i++) {
        OS >> vblock[i];
      }

      return OS;
    }


  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::vector<Vector3DB> vec;
    Real* c;
  };

}



#endif /* VECTOR3DBLOCK_H */
