/*  -*- c++ -*-  */
#ifndef FFTCOMPLEX_H
#define FFTCOMPLEX_H

  // Define the data type for FFT.
#ifndef _SGI_FFT_
  // Do not define if already defined by <fft.h>
  typedef struct {
    double re;
    double im;
  } zomplex;
#endif

#include <protomol/type/Real.h>

namespace ProtoMol {

  class FFTInternal;

  //_________________________________________________________________ FFTComplex
  /**
   * Wrapper class for complex 3D FFT. Library specific calls and  
   * initialization are implemented in FFTInternal.@n
   *
   * 1D and 2D: Overload method initialize()@n
   *
   * Run-time selection: Make FFTInternal abstract and the implementation
   * in local classes.
   */
  class FFTComplex {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    FFTComplex();
    ~FFTComplex();
  private:
    FFTComplex(const FFTComplex&);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class FFTComplex
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(int x, int y, int z, zomplex* a);
    void backward();
    void forward();

  public:
    /// initialize MPI with master id
    static void FFTComplexMPIInit(int master);
    /// finalize MPI
    static void FFTComplexMPIFinalize();    
    /// test if parallel environment
    static bool isParallel();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    FFTInternal* myFFTInternal;
  };
}
#endif /* FFTCOMPLEX_H */
