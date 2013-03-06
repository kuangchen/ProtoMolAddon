#include <protomol/base/Random.h>
#include <protomol/base/Report.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

//____ Random

const int Random::SHUFFLE_LEN = SHUFFLE_LEN_VAL;
const double Random::EPSILON = 3.0e-16;
const double Random::INVERSE_MODULUS = 1.0 / MUDULUS;
const double Random::LESS_THAN_ONE = 1.0 - EPSILON;

Random::Random(){
  int j, k;

  shuffle_array.resize(SHUFFLE_LEN);
  shuffle_val = 0;
  intermediate_rand = SEED_DEFAULT;
  //Initialize.
  if (intermediate_rand < 1) intermediate_rand = 1;  //Be sure to prevent intermediate_rand = 0.
  for (j=SHUFFLE_LEN+7;j>=0;j--) {                   //Load the shuffle table (after 8 warm-ups).
    k = intermediate_rand / SCHRAGE_Q;
    intermediate_rand = MULTIPLIER * (intermediate_rand - k * SCHRAGE_Q) - SCHRAGE_R * k;
    if (intermediate_rand < 0) intermediate_rand += MUDULUS;
    if (j < SHUFFLE_LEN) shuffle_array[j] = intermediate_rand;
  }
  shuffle_val=shuffle_array[0];

}

void Random::seed( unsigned int seed ) {
  int j, k;

  shuffle_val = 0;
  intermediate_rand = seed;
  //Initialize.
  if (intermediate_rand < 1) { intermediate_rand = 1;  //Be sure to prevent intermediate_rand = 0.
    report << plain << "Warning: Seed must not be less than unity!" << endr;
  }
  for (j=SHUFFLE_LEN+7;j>=0;j--) {   //Load the shuffle table (after 8 warm-ups).
    k = intermediate_rand / SCHRAGE_Q;
    intermediate_rand = MULTIPLIER * (intermediate_rand - k * SCHRAGE_Q) - SCHRAGE_R * k;
    if (intermediate_rand < 0) intermediate_rand += MUDULUS;
    if (j < SHUFFLE_LEN) shuffle_array[j] = intermediate_rand;
  }
  shuffle_val = shuffle_array[0];
}

double Random::rand() {

  int j, k;
  double temp;

  k = intermediate_rand / SCHRAGE_Q;                          //Start here when not initializing.
  intermediate_rand = MULTIPLIER * (intermediate_rand - k * SCHRAGE_Q)
                        - SCHRAGE_R * k;                      //Compute intermediate_rand=(MULTIPLIER*intermediate_rand) % MUDULUS without overiflows by Scrage's method.
  if(intermediate_rand < 0) intermediate_rand += MUDULUS;     //flows by Schrage's method.
  j = shuffle_val / SHUFFLE_DIVISOR;                          //Will be in the range 0..SHUFFLE_LEN-1.
  shuffle_val = shuffle_array[j];                             //Output previously stored value and refill the
  shuffle_array[j] = intermediate_rand;                       //shuffle table.

  if ((temp = INVERSE_MODULUS*shuffle_val) > LESS_THAN_ONE)
    return LESS_THAN_ONE;   //Users expect (0,1)
  else return temp;

}

namespace ProtoMol {

  ostream & operator<<(ostream &stream, Random ob)
  {
    stream.precision(10);
    for (int i=0; i < ob.SHUFFLE_LEN; i++)
      stream << ob.shuffle_array[i] << ' ';

    stream << ob.shuffle_val << ' ' << ob.intermediate_rand;

    return stream;
  }

  istream & operator>>(istream &stream, Random &ob)
  {
    for (int i=0; i < ob.SHUFFLE_LEN; i++)
      stream >> ob.shuffle_array[i];

    stream >> ob.shuffle_val >> ob.intermediate_rand;

    return stream;
  }

}
