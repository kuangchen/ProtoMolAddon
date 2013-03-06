#include <protomol/base/PMConstants.h>

#ifdef HAVE_NO_LIMITS
#ifndef WIN32
#include <values.h>
#endif
#else /* HAVE_NO_LIMITS */
#include <limits>
#endif /* HAVE_NO_LIMITS */

#include <protomol/base/StringUtilities.h>

using namespace std;

namespace ProtoMol {
  namespace Constant {
#ifdef HAVE_NO_LIMITS
#ifdef USE_REAL_IS_DOUBLE
    const Real MAXREAL = MAXDOUBLE;
    const Real MINREAL = MINDOUBLE;
#endif
#ifdef USE_REAL_IS_FLOAT
    const Real MAXREAL = MAXFLOAT;
    const Real MINREAL = MINFLOAT;
#endif

    const Real REAL_INFINITY = 2.0 * MAXREAL;
    const Real REAL_NAN = REAL_INFINITY;
    const int MAX_INT = MAXINT;
    const int MAX_INT_2 = MAXINT / 2;

#else /* HAVE_NO_LIMITS */
    const Real MAXREAL = numeric_limits<Real>::max();
    const Real MINREAL = numeric_limits<Real>::min();
    const Real REAL_INFINITY =
      (numeric_limits<Real>::has_infinity ? numeric_limits<Real>::
         infinity() :
       2.0 * MAXREAL);
    const Real REAL_NAN =
      (numeric_limits<Real>::has_quiet_NaN ? numeric_limits<Real>::
         quiet_NaN() :
       REAL_INFINITY);
    const int MAX_INT = numeric_limits<int>::max();
    const int MAX_INT_2 = numeric_limits<int>::max() / 2;
#endif /* HAVE_NO_LIMITS */

    const Real EPSILON = 1.0e-14;
    const Real TINY = 1.0e-20;

    const Real TIMEFACTOR = 48.88821290839616; // TIMEUNIT is 1 / sqrt(4.184e-4)
    const Real INV_TIMEFACTOR = 0.02045482828087295; //  1 / TIMEFACTOR

    const Real PERIODIC_BOUNDARY_TOLERANCE = 3.0;

    const Real EPS_GOURAUD_THRESHOLD = 0.1;
    const Real EPS_SMOOTH_LINE_FACTOR = 0.06;

    const char* PROTOMOL_HR = "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+";
    const char* PRINTINDENT = "  ";
    const unsigned int PRINTMAXWIDTH = 30;

    const int FASTDELTAMAX = 32;

    const Real SQRTCOULOMBCONSTANT = 18.2226123264;
    // ENERGY_TO_SI/(LENGTH_TO_SI*LENGTH_TO_SI*LENGTH_TO_SI)*1e-5; // bar
    const Real PRESSUREFACTOR = 69478.0593635551;
    const Real BOLTZMANN = 0.001987191;
    const Real PDBVELSCALINGFACTOR = 20.45482706;

    // Constants for GROMACS-PROTOMOL conversion

    const Real ANGSTROM_NM = 0.1;
    const Real NM_ANGSTROM = 10.0;
    const Real INV_ANGSTROM_NM = 10.0;
    const Real INV_NM_ANGSTROM = 0.1;
    const Real KCAL_KJ = 4.184;
    const Real KJ_KCAL = 1.0 / 4.184;
    const Real FS_PS = 0.001;
    const Real PS_FS = 1000.0;

    namespace SI {
      const Real C = 299792458.0;    // [m/s]
      const Real COULOMB_FACTOR = C * C * 1e-7;       // [Vm/C]
      const Real ELECTRON_CHARGE = 1.6021892e-19;  // [C]
      const Real LENGTH_AA = 1e+10;          // [AA]
      const Real AVOGADRO = 6.022045e+23;   // [1/mol]
      const Real AMU = 1.6605655e-27;  // [kg]
      const Real KCAL = 1.0 / 4184.0;     // [J]
      const Real TIME_FS = 1e+15;          // [fs]
      const Real BOLTZMANN = 1.380662e-23;   // [J/K]
    }

    string print() {
      string res;
      res =
        "MAXREAL       = " + toString(MAXREAL) + "\n" +
        "MINREAL       = " + toString(MINREAL) + "\n" +
        "REAL_INFINITY = " + toString(REAL_INFINITY) + "\n" +
        "MAX_INT       = " + toString(MAX_INT) + "\n" +
        "MAX_INT_2     = " + toString(MAX_INT_2) + "\n" +
        "EPSILON       = " + toString(EPSILON) + "\n" +
        "TINY          = " + toString(TINY) + "\n";
      return res;
    }
  }
}
