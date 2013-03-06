#ifndef RANDOM_H
#define RANDOM_H

#include <protomol/type/Real.h>
#include <protomol/base/Singleton.h>

#include <vector>
#include <fstream>
#include <iostream>

namespace ProtoMol {

  //_____________________________________________________________ Random
  /**
   * Random number generation
   * "Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added
   * safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
   * values). Do not alter intermediate_rand between
   * successive deviates in a sequence. LESS_THAN_ONE should approximate the largest floating value that is
   * less than 1.
   */
  class Random : public Singleton<Random> {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Types & enum's
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    enum {SEED_DEFAULT = 1234};
    enum {MULTIPLIER = 16807};
    enum {MUDULUS = 2147483647};
    enum {SCHRAGE_Q = 127773};
    enum {SCHRAGE_R = 2836};
    enum {SHUFFLE_LEN_VAL = 32};
    enum {SHUFFLE_DIVISOR = (1+(MUDULUS-1)/SHUFFLE_LEN_VAL)};

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Constructors, destructors, assignment
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Random();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  New methods of class Random
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void seed( unsigned int seed );
    double rand();

    friend std::ostream &operator<<(std::ostream &stream, Random ob);
    friend std::istream &operator>>(std::istream &stream, Random &ob);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // My data members
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    int intermediate_rand;
    std::vector<int> shuffle_array;
    int shuffle_val;
    static const int SHUFFLE_LEN;
    static const double EPSILON, INVERSE_MODULUS, LESS_THAN_ONE;
  };
};

#endif // RANDOM_H
