/*  -*- c++ -*-  */
#ifndef LOOKUPTABLEBASE_H
#define LOOKUPTABLEBASE_H

#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/type/TypeSelection.h>

namespace ProtoMol {
  //____ LookupTable

  /**
   * Generic look up table covering the range [2^n0 <= from, to < 2^n1] with 
   * look up (hashed) values defined by TLookupValue and digit/bit precision 
   * PRE. lookupExp is the exponent of the hash value, yields also for the
   * range passed to the constructor, e.g., if you like to use a look up table
   * for the squared norm, then set lookupExp to 2. LookupTable does also
   * support switching functions. LookupTable is the base class of any
   * tabulated potential (e.g. CoulombForce). TLookupValue implements the
   * function to be hashed with arbitrary switching function.
   */
  template<class TLookupValue, unsigned int PRE, typename TReal>
  class LookupTable {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef's & Enum's
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    typedef typename TypeSelection::Int < sizeof(TReal) > ::type TInt;
    typedef union {TReal f; TInt i;} Conv;
  public:
    /// number of values per entry (0: value, 2: gradient, 2: 2. interp.,
    /// 3: 3. interp.)
    enum {NV = 4};
    /// number of bits representing the exponent
    enum {E = BinaryExponent < std::numeric_limits<TReal>::max_exponent -
              std::numeric_limits<TReal>::min_exponent > ::N};
    /// entries per look up up
    enum {STEP = TLookupValue::ENTRIES * NV};
    /// precision range/size of mantissa
    enum {N = (1 << PRE)};
    /// significant bits of mantissa
    enum {PRECISION = PRE};
    /// full digit mantissa mask
    enum {M = (sizeof(TReal) * 8 - 1 - E)};
    /// look up mask for hashing
    enum {SHIFT = (sizeof(TReal) * 8 - 1 - E - PRE)};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    LookupTable() :
      myTable(0), myFirst(0), myLast(0), mySize(0), myAlign(0), myMem(0) {}

    LookupTable(const LookupTable &rhs) {
      if (&rhs != this) {
        mySize = rhs.mySize;
        myFirst = rhs.myFirst;
        myLast = rhs.myLast;
        myAlign = rhs.myAlign;
        myMem = new TReal[mySize + myAlign / sizeof(TReal)];
        myTable = myMem;
        while ((long)(myTable) % myAlign)
          ++myTable;

        for (int i = 0; i < mySize; i++)
          myTable[i] = rhs.myTable[i];
      }
    }

    template<class TSwitchingFunction>
    LookupTable(TReal from, TReal to, int lookupExp, const TLookupValue &lvf,
                const TSwitchingFunction &swf, unsigned int align = 128) :
      myTable(0), myFirst(absindex(lowerboundExp2(from))),
      myLast(((expn(to) + 1 - expn(from)) * N) * STEP),
      mySize(((expn(to) + 1 - expn(from)) * N + 1) * STEP), myAlign(align),
      myMem(new TReal[((expn(to) + 1 - expn(from)) * N + 1) * STEP + align /
                      sizeof(TReal)]) {
      // align table
      myTable = myMem;
      while ((long)(myTable) % align)
        ++myTable;

      from = lowerboundExp2(from);
      // compute look up table
      Real errv = -1.0;
      Real errg = -1.0;
      for (int i = 0; i < mySize / STEP; i++) {
        // value/index
        Real base = from * (1 << (i / N));
        Real del = base / static_cast<TReal>(N);
        Real r = base + del * (i % N);
        Real r1;
        Real r2;

        // retrieve r and r^2
        if (lookupExp == -2) {
          r1 = 1.0 / sqrt(r);
          r2 = 1.0 / r;
        } else if (lookupExp == -1) {
          r1 = 1.0 / r;
          r2 = 1.0 / square(r);
        } else if (lookupExp == 1) {
          r1 = r;
          r2 = square(r);
        } else if (lookupExp == 2) {
          r1 = sqrt(r);
          r2 = r;
        } else {
          r1 = pow(r, 1.0 / static_cast<float>(lookupExp));
          r2 = pow(r, 2.0 / static_cast<float>(lookupExp));
        }

        // assign look up value(s) and clear 2. and 3. interpolation values
        Real v, d;
        swf(v, d, r2);
        lvf.assign(r, r1, lookupExp, v, d, &myTable[i * STEP]);
        for (unsigned int j = 0; j < static_cast<unsigned int>(STEP);
             j += NV) {
          myTable[2 + i * STEP + j] = 0.0;
          myTable[3 + i * STEP + j] = 0.0;
        }

        // interpolate 2. and 3. values
        if (i > 0) {
          Real x = from * (1 << ((i - 1) / N)) / static_cast<TReal>(N);
          TReal *t = &myTable[(i - 1) * STEP];
          for (unsigned int j = 0; j < static_cast<unsigned int>(STEP);
               j += NV) {
            Real v1 = t[0 + j];
            Real g1 = t[1 + j];
            Real v2 = t[0 + STEP + j];
            Real g2 = t[1 + STEP + j];
            // explicit formulas for v1 + g1 x + c x^2 + d x^3
            Real c = (3.0 * (v2 - v1) - x * (2.0 * g1 + g2)) / (x * x);
            Real d = (-2.0 * (v2 - v1) + x * (g1 + g2)) / (x * x * x);
            // since v2 - v1 is imprecise, we refine c and d numerically
            // important because we need accurate forces (more than energies!)
            for (unsigned int k = 0; k < 2; ++k) {
              Real dv = (v1 - v2) + ((d * x + c) * x + g1) * x;
              Real dg = (g1 - g2) + (3.0 * d * x + 2.0 * c) * x;
              c -= (3.0 * dv - x * dg) / (x * x);
              d -= (-2.0 * dv + x * dg) / (x * x * x);
            }

            t[2 + j] = c;
            t[3 + j] = d;
            Real ev = fabs(((t[3 + j] * x + t[2 + j]) * x + t[1 + j]) * x + 
                           t[0 + j] - t[0 + STEP + j]);
            Real eg = fabs((3.0 * t[3 + j] * x + 2.0 * t[2 + j]) * x +
                           t[1 + j] - t[1 + STEP + j]);
            errv = (errv < ev ? ev : errv);
            errg = (errg < eg ? eg : errg);
          }
        }
      }

      Report::report
        << Report::debug(2)
        << "Number of interpolation points:" << mySize / STEP << ","
        << (long)(myTable - myMem);

      if (errv > 0.0 || errg > 0.0)
        Report::report << ", error: energy=" << errv << ", force=" << errg;

      Report::report << Report::endr;
    }

    virtual ~LookupTable() {
      if (myMem != 0) delete[] myMem;
      myMem = 0;
    }

    LookupTable &operator=(LookupTable const &rhs) {
      if (&rhs != this) {
        mySize = rhs.mySize;
        myFirst = rhs.myFirst;
        myLast = rhs.myLast;
        myAlign = rhs.myAlign;
        if (myMem != 0) delete[] myMem;
        myMem = new TReal[mySize + myAlign / sizeof(TReal)];
        myTable = myMem;
        while ((long)(myTable) % myAlign)
          ++myTable;

        for (int i = 0; i < mySize; i++)
          myTable[i] = rhs.myTable[i];
      }
      return *this;
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LookupTable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    /// computes the index for the hash table and the delta r
    void index(Real r, int &i, Real &dt) const {
      Conv u;
      u.f = r;
      i = (u.i >> SHIFT) * STEP - myFirst;
      dt = 0.0;
      if (i < 0) i = 0;
      else if (i > myLast) i = myLast;
      else {
        u.i &= myMask;
        dt = r - u.f;
      }
    }

    /// interpolates the enery and force
    void interpolate(Real a, Real b, Real c, Real d, Real dt, Real &energy,
                     Real &force) const {
      energy = ((dt * d + c) * dt + b) * dt + a;
      force = -2.0 * ((3.0 * dt * d + 2.0 * c) * dt + b);
    }

  private:
    /// retrieve exponent from a
    static int expn(Real a) {
      Conv r;
      r.f = lowerboundExp2(a);
      return r.i >> M;
    }

    /// absolute index based on precision and exponent
    static int absindex(Real r) {
      Conv u;
      u.f = r;
      return (u.i >> SHIFT) * STEP;
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    TReal *myTable;
  private:
    int myFirst;
    int myLast;
    int mySize;
    unsigned int myAlign;
    TReal *myMem;
    static const TInt myMask;
  };

  template<class TLookupValue, unsigned int PRE, typename TReal>
  const typename LookupTable<TLookupValue, PRE, TReal>::TInt
  LookupTable<TLookupValue, PRE, TReal>::myMask = ~((1LL << SHIFT) - 1LL);
}
#endif /* LOOKUPTABLEBASE_H */
