//
// FFT includes and forward declarations
//
#ifdef HAVE_FFT_SGI
#  include <protomol/parallel/FFTComplex.h> !!!
                   // We use the definition of zomplex in <fft.h>
#else
#  ifdef HAVE_FFT_FFTW3
#    include <fftw3.h>
#  else
#     ifdef HAVE_FFT_ESSL
#       define _ESV_COMPLEX_ 1
#       include <essl.h>
#     else
#       ifdef HAVE_FFT_FFTW2
#         include <fftw.h>
#       else
#         ifdef HAVE_FFT_FFTW2_MPI
#           include <fftw.h>
#           include <fftw_mpi.h>
#         else /* ZFFT */
#           include <protomol/parallel/FFTComplex.h>     // Defines zomplex
extern "C" {
  extern zomplex *zfftm1di(int m, zomplex *save);
  extern int zfftm1d(int sign, int m, int n, zomplex *array, 
                     int incI, int incJ, zomplex *save);
  
  extern zomplex *zfft2di( int n1, int n2, zomplex *save);
  extern int zfft2d(int sign, int n1, int n2, zomplex *array, 
                    int ld, zomplex *save);
  
  extern zomplex *zfft3di(int n1, int n2, int n3, zomplex *save);
  extern int zfft3d(int sign, int n1, int n2, int n3, zomplex *array, 
                    int ld1, int ld2, zomplex *save);
}
#         endif
#       endif
#     endif
#  endif
#endif

//
// Other includes
//
#include <protomol/parallel/FFTComplex.h>
#include <protomol/base/Timer.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/base/Report.h>
using namespace ProtoMol::Report;

namespace ProtoMol {

  //
  // Define the right implementation of FFTInternal
  //
#if defined(HAVE_FFT_SGI) || !defined(HAVE_FFT)
  //________________________________________________________________ FFTInternal
  //
  // complib.sgimath and ZFFT
  //
  class FFTInternal {

  public:
    FFTInternal():myNX(0),myNY(0),myNZ(0),myArray(NULL),myFFTCoeff(NULL){}
    ~FFTInternal();

    void initialize(int x, int y, int z,zomplex* a);
    void forward() {zfft3d( 1,myNZ,myNY,myNX,myArray,myNZ,myNY,myFFTCoeff);}
    void backward(){zfft3d(-1,myNZ,myNY,myNX,myArray,myNZ,myNY,myFFTCoeff);}
  private:
    FFTInternal(const FFTInternal&);

  private:
    int myNX;
    int myNY;
    int myNZ;
    zomplex* myArray;
    zomplex* myFFTCoeff;
  };

  FFTInternal::~FFTInternal(){
    if(myFFTCoeff != NULL)
      free(myFFTCoeff);    
  }

  void FFTInternal::initialize(int x, int y, int z,zomplex* a){
    if(x == myNX && y == myNY && z == myNZ && a == myArray)
      return;
    myNX = x;
    myNY = y;
    myNZ = z;
    myArray = a;
    if(myFFTCoeff != NULL)
      free(myFFTCoeff);
    myFFTCoeff = zfft3di(myNZ,myNY,myNX,NULL);    
  }

#else
#  ifdef HAVE_FFT_FFTW3
  //________________________________________________________________ FFTInternal
  //
  // FFTW3
  //
  class FFTInternal {

  public:
    FFTInternal():
      myNX(0),myNY(0),myNZ(0),myArray(NULL),myPlanForward(NULL),
      myPlanBackward(NULL){}
    ~FFTInternal();

    void initialize(int x, int y, int z,zomplex* a);
    void forward(){ fftw_execute(myPlanForward);}
    void backward(){fftw_execute(myPlanBackward);}

  private:
    FFTInternal(const FFTInternal&);

  private:
    int myNX;
    int myNY;
    int myNZ;
    fftw_complex* myArray;
    fftw_plan myPlanForward;
    fftw_plan myPlanBackward;
  };

  FFTInternal::~FFTInternal(){
    if(myPlanForward != NULL)
      fftw_destroy_plan(myPlanForward);
    if(myPlanBackward != NULL)
      fftw_destroy_plan(myPlanBackward);     
  }

  void FFTInternal::initialize(int x, int y, int z,zomplex* a){
    if(x == myNX && y == myNY && z == myNZ && (fftw_complex*)(a) == myArray)
      return;
    Timer t;
    t.start();
    myNX = x;
    myNY = y;
    myNZ = z;
    myArray = (fftw_complex*)(a);
    if(myPlanForward != NULL)
      fftw_destroy_plan(myPlanForward);
    if(myPlanBackward != NULL)
      fftw_destroy_plan(myPlanBackward);
    myPlanForward  = fftw_plan_dft_3d(myNX,myNY,myNZ,myArray,myArray,
                                      FFTW_FORWARD,  FFTW_MEASURE);
    myPlanBackward = fftw_plan_dft_3d(myNX,myNY,myNZ,myArray,myArray,
                                      FFTW_BACKWARD, FFTW_MEASURE);
    t.stop();
    report << hint <<"FFTW3 Initialization: "<<t.getTime()<<"."<<endr;
  }

#  else
#     ifdef HAVE_FFT_ESSL
  //________________________________________________________________ FFTInternal
  //
  // ESSL
  // 
  class FFTInternal {

  public:
    FFTInternal():myNX(0),myNY(0),myNZ(0),myArray(NULL),myAux(NULL){}
    ~FFTInternal();

    void initialize(int x, int y, int z,zomplex* a);
    void forward(){
      dcft3(myArray, myNZ, myNY*myNZ, myArray, myNZ, myNY*myNZ,
            myNZ, myNY, myNX,  1, 1.0, myAux, myNaux);
    }
    void backward(){
      dcft3(myArray, myNZ, myNY*myNZ, myArray, myNZ, myNY*myNZ, myNZ, myNY,
            myNX, -1, 1.0, myAux, myNaux);
    }
  private:
    int case1(int n1,int n2, int n3);
    int case2(int n1,int n2, int n3);
    int case3(int n1,int n2, int n3);
    FFTInternal(const FFTInternal&);

  private:
    int myNX;
    int myNY;
    int myNZ;
    zomplex* myArray;
    double*  myAux;
    int myNaux;
  };

  FFTInternal::~FFTInternal(){
    if(myAux != NULL)
      delete [] myAux;
  }

  void FFTInternal::initialize(int x, int y, int z,zomplex* a){
    if(x == myNX && y == myNY && z == myNZ && a == myArray)
      return;
    myNX = x; // n3
    myNY = y; // n2
    myNZ = z; // n1
    myArray = a;
  
    if(std::max(myNY, myNX) < 252){
      myNaux=case1(myNZ,myNY,myNX);
    }
    else if (myNX < 252) {
      myNaux=case2(myNZ,myNY,myNX);
    }
    else if (myNY < 252) {
      myNaux=case3(myNZ,myNY,myNX);
    }
    else {
      myNaux=std::max(case2(myNZ,myNY,myNX),case3(myNZ,myNY,myNX));
    }
    if(myAux != NULL)
      delete [] myAux;
    myAux = new double [myNaux];
    report << hint <<"ESSL Initialization: naux="<<myNaux<<"."<<endr;
  }

  int FFTInternal::case1(int n1,int n2, int n3){
    return static_cast<int>(60000 + (n1 <= 2048 ? 0 : 4.56*n1));
  }

  int FFTInternal::case2(int n1,int n2, int n3){
    return static_cast<int>(60000+(2*n2+256)*(std::min(64, n1)+4.56) +
                            (n1 <= 2048 ? 0 : 4.56*n1));
  }

  int FFTInternal::case3(int n1,int n2, int n3){
    return static_cast<int>(60000+(2*n3+256)*(std::min(64, n1*n2)+4.56) +
                            (n1 <= 2048 ? 0 : 4.56*n1));
  }
#     else
#       ifdef HAVE_FFT_FFTW2
  //________________________________________________________________ FFTInternal
  //
  // FFTW2
  //
  class FFTInternal {

  public:
    FFTInternal():
      myNX(0),myNY(0),myNZ(0),myArray(NULL),
      myPlanForward(NULL),myPlanBackward(NULL){}
    ~FFTInternal();

    void initialize(int x, int y, int z,zomplex* a);
    void forward(){ fftwnd_one(myPlanForward, myArray,NULL);}
    void backward(){fftwnd_one(myPlanBackward,myArray,NULL);}

  private:
    FFTInternal(const FFTInternal&);

  private:
    int myNX;
    int myNY;
    int myNZ;
    fftw_complex* myArray;
    fftwnd_plan myPlanForward;
    fftwnd_plan myPlanBackward;
  };

  FFTInternal::~FFTInternal(){
    if(myPlanForward != NULL)
      fftwnd_destroy_plan(myPlanForward);
    if(myPlanBackward != NULL)
      fftwnd_destroy_plan(myPlanBackward);     
  }

  void FFTInternal::initialize(int x, int y, int z,zomplex* a){
    if(x == myNX && y == myNY && z == myNZ && (fftw_complex*)(a) == myArray)
      return;
    Timer t;
    t.start();
    myNX = x;
    myNY = y;
    myNZ = z;
    myArray = (fftw_complex*)(a);
    if(myPlanForward != NULL)
      fftwnd_destroy_plan(myPlanForward);
    if(myPlanBackward != NULL)
      fftwnd_destroy_plan(myPlanBackward);
    myPlanForward  = fftw3d_create_plan(myNX,myNY,myNZ,FFTW_FORWARD,
                                        FFTW_MEASURE|FFTW_IN_PLACE);
    myPlanBackward = fftw3d_create_plan(myNX,myNY,myNZ,FFTW_BACKWARD,
                                        FFTW_MEASURE|FFTW_IN_PLACE);
    t.stop();
    report << hint <<"FFTW2 Initialization: "<<t.getTime()<<"."<<endr;
  }
#       else
#         ifdef HAVE_FFT_FFTW2_MPI
  //________________________________________________________________ FFTInternal
  //
  // FFTW2 MPI
  //

#ifdef USE_REAL_IS_FLOAT
#define MY_MPI_REAL MPI_FLOAT
#endif

#ifdef USE_REAL_IS_DOUBLE
#define MY_MPI_REAL MPI_DOUBLE
#endif


  class FFTInternal {

  public:
    FFTInternal():
      myNX(0),myNY(0),myNZ(0),myArray(NULL),myArrayTmp(NULL),
      myLocalComm(MPI_COMM_NULL),
      myPlanForward(NULL),myPlanBackward(NULL),data(NULL),work(NULL){}
    ~FFTInternal();

    void initialize(int x, int y, int z,zomplex* a);
    void forward();
    void backward();

  public:
    static void FFTInternalMPIInit(int master);
    static void FFTInternalMPIFinalize();

  private:
    FFTInternal(const FFTInternal&);

  private:
    int myNX;
    int myNY;
    int myNZ;
    fftw_complex* myArray;
    fftw_complex* myArrayTmp;
    // All possible available nodes, master exlcuded if master-slave
    static MPI_Comm myComm;
    int myNum;

    MPI_Comm myLocalComm;    // Nodes to be used by FFT
    fftwnd_mpi_plan myPlanForward;
    fftwnd_mpi_plan myPlanBackward;
    int local_nx;
    int local_x_start;
    int local_ny_after_transpose;
    int local_y_start_after_transpose;
    int total_local_size;
    fftw_complex* data;
    fftw_complex* work;
  };

  MPI_Comm  FFTInternal::myComm  = MPI_COMM_NULL;

  FFTInternal::~FFTInternal(){
    if(myComm  != MPI_COMM_NULL){
      if(myPlanForward != NULL)
	fftwnd_mpi_destroy_plan(myPlanForward);
      if(myPlanBackward != NULL)
	fftwnd_mpi_destroy_plan(myPlanBackward);     
      if(myLocalComm != MPI_COMM_NULL && myLocalComm !=
         MPI_COMM_WORLD && myLocalComm != myComm )
	MPI_Comm_free(&myLocalComm);    
    }
    if(work != NULL)
      delete [] work;
    if(data != NULL)
      delete [] data;
    if(myArrayTmp != NULL)
      delete [] myArrayTmp;
  }

  void FFTInternal::initialize(int x, int y, int z,zomplex* a){
    if(x == myNX && y == myNY && z == myNZ && (fftw_complex*)(a) == myArray)
      return;
    Timer t;
    t.start();
    myNX = x;
    myNY = y;
    myNZ = z;
    myArray = (fftw_complex*)a;

    // Clean
    if(myPlanForward != NULL){
      fftwnd_mpi_destroy_plan(myPlanForward);
      myPlanForward = NULL;
    }
    if(myPlanBackward != NULL){
      fftwnd_mpi_destroy_plan(myPlanBackward);
      myPlanBackward = NULL;
    }
    if(work != NULL){
      delete [] work;
      work = NULL;
    }
    if(data != NULL){
      delete [] data;
      data = NULL;
    }
    if(myArrayTmp != NULL){
      delete [] myArrayTmp;
      myArrayTmp = NULL;
    }
    if(myLocalComm != MPI_COMM_NULL &&
       myLocalComm != MPI_COMM_WORLD && myLocalComm != myComm )
      MPI_Comm_free(&myLocalComm);    


    // Select the right number of nodes
    int rank, size;
    MPI_Comm_rank(myComm,&rank);
    MPI_Comm_size(myComm,&size);
    myNum = std::min(myNX,size);
    if(myNum < size){
      // We need to remove some
      MPI_Group worldGroup = MPI_GROUP_NULL;
      MPI_Group slaveGroup = MPI_GROUP_NULL;
      int* excl = new int [size-myNum];
      for(int i=myNum;i<size;i++)
	excl[i-myNum] = i;
      MPI_Comm_group(myComm,&worldGroup);
      MPI_Group_excl(worldGroup,size-myNum,excl,&slaveGroup);
      MPI_Comm_create(myComm,slaveGroup,&myLocalComm);
      MPI_Group_free(&worldGroup);
      MPI_Group_free(&slaveGroup);
      delete [] excl;    
    }
    else {
      // Ok, use all available nodes
      myLocalComm = myComm;
    }
    if(size > 1)
      myArrayTmp = new  fftw_complex [myNX*myNY*myNZ];

    if(myLocalComm != MPI_COMM_NULL){
      myPlanForward  =
        fftw3d_mpi_create_plan(myLocalComm,myNX,myNY,myNZ,FFTW_FORWARD,
                               FFTW_MEASURE|FFTW_IN_PLACE);
      myPlanBackward =
        fftw3d_mpi_create_plan(myLocalComm,myNX,myNY,myNZ,FFTW_BACKWARD,
                               FFTW_MEASURE|FFTW_IN_PLACE);
      if(myNum > 1){
	fftwnd_mpi_local_sizes(myPlanForward, &local_nx, &local_x_start,
			       &local_ny_after_transpose,
			       &local_y_start_after_transpose,
			       &total_local_size);
	work = new fftw_complex [total_local_size];
	data = new fftw_complex [total_local_size];
	// Report if a node does nothing
	if(local_nx == 0){
	  int n;
	  MPI_Comm_rank(MPI_COMM_WORLD,&n);
	  report << allnodes << hint << "FFTW2 MPI: Node "<<n
             <<" does nothing."<<endr; 
	}
      }
    }
    t.stop();
    if(rank == 0){
      report << allnodes << hint <<"FFTW2 MPI Initialization with "
             <<myNum<<" node(s) of "<<size<<": "<<t.getTime()<<"."<<endr;
    }
  }

  void FFTInternal::forward(){ 
    if(myPlanForward){
      if(myNum <= 1){
	fftwnd_mpi(myPlanForward, 1,myArray,NULL,FFTW_NORMAL_ORDER);
      }
      else {
	for (int x = 0; x < local_nx; ++x)
	  for (int y = 0; y < myNY; ++y)
	    for (int z = 0; z < myNZ; ++z)
	      data[(x*myNY + y)*myNZ + z] =
            myArray[((x + local_x_start)*myNY +  y)*myNZ + z];

	fftwnd_mpi(myPlanForward, 1,data,work,FFTW_TRANSPOSED_ORDER);

	for(int i = 0;i<myNX*myNY*myNZ;i++){
	  myArrayTmp[i].re = 0.0;
	  myArrayTmp[i].im = 0.0;
	}
	for (int y = 0; y < local_ny_after_transpose; ++y)
	  for (int x = 0; x < myNX; ++x)
	    for (int z = 0; z < myNZ; ++z)
	      myArrayTmp[(x*myNY + local_y_start_after_transpose +  y)*myNZ + z] =
            data[(y*myNX + x) * myNZ + z];

	// Split up allreduce if not all nodes were used for FFT
	if(myComm != myLocalComm){
	  MPI_Reduce((Real*)myArrayTmp, (Real*)myArray,
                 2*myNX*myNY*myNZ, MY_MPI_REAL, MPI_SUM, 0, myLocalComm);
	}
	else {
	  MPI_Allreduce((Real*)myArrayTmp, (Real*)myArray,
                    2*myNX*myNY*myNZ, MY_MPI_REAL, MPI_SUM, myLocalComm);  
	}
      }
    }
    // If we excluded nodes we have to broad cast ...
    if(myComm != myLocalComm){
      MPI_Bcast((Real*)myArray,2*myNX*myNY*myNZ, MY_MPI_REAL,0, myComm);      
    }
  }

  void FFTInternal::backward(){
    if(myPlanBackward){
      if(myNum <= 1){
	fftwnd_mpi(myPlanBackward,1,myArray,NULL,FFTW_NORMAL_ORDER);
      }
      else {
	for (int x = 0; x < local_nx; ++x)
	  for (int y = 0; y < myNY; ++y)
	    for (int z = 0; z < myNZ; ++z)
	      data[(x*myNY + y)*myNZ + z] =
            myArray[((x + local_x_start)*myNY +  y)*myNZ + z];

	fftwnd_mpi(myPlanBackward, 1,data,work,FFTW_TRANSPOSED_ORDER);

	for(int i = 0;i<myNX*myNY*myNZ;i++){
	  myArrayTmp[i].re = 0.0;
	  myArrayTmp[i].im = 0.0;
	}
	for (int y = 0; y < local_ny_after_transpose; ++y)
	  for (int x = 0; x < myNX; ++x)
	    for (int z = 0; z < myNZ; ++z)
	      myArrayTmp[(x*myNY + local_y_start_after_transpose +  y)*myNZ + z] =
            data[(y*myNX + x) * myNZ + z];

	// Split up allreduce if not all nodes were used for FFT
	if(myComm != myLocalComm){
	  MPI_Reduce((Real*)myArrayTmp, (Real*)myArray, 2*myNX*myNY*myNZ,
                 MY_MPI_REAL, MPI_SUM, 0, myLocalComm);
	}
	else {
	  MPI_Allreduce((Real*)myArrayTmp, (Real*)myArray, 2*myNX*myNY*myNZ,
                    MY_MPI_REAL, MPI_SUM, myLocalComm);  
	}
      }
    }
    // If we excluded nodes we have to broad cast ...
    if(myComm != myLocalComm){
      MPI_Bcast((Real*)myArray,myNX*myNY*myNZ*2, MY_MPI_REAL,0 , myComm);      
    }
  }



  void FFTInternal::FFTInternalMPIInit(int master){
    if(myComm != MPI_COMM_NULL && myComm != MPI_COMM_WORLD){
      MPI_Comm_free(&myComm);
    }
    if(master < 0){
      // All nodes are slaves, no master
      myComm = MPI_COMM_WORLD;
    }
    else {
      // Create intracommunicator only with slaves
      MPI_Group worldGroup = MPI_GROUP_NULL;
      MPI_Group slaveGroup = MPI_GROUP_NULL;
      int excl[] = {master};

      MPI_Comm_group(MPI_COMM_WORLD,&worldGroup);
      MPI_Group_excl(worldGroup,1,excl,&slaveGroup);
      MPI_Comm_create(MPI_COMM_WORLD,slaveGroup,&myComm);
      MPI_Group_free(&worldGroup);
      MPI_Group_free(&slaveGroup);
    }
  }
  void FFTInternal::FFTInternalMPIFinalize(){
    if(myComm != MPI_COMM_NULL && myComm != MPI_COMM_WORLD){
      MPI_Comm_free(&myComm);
    }
    myComm = MPI_COMM_NULL;
  }
#         endif
#       endif
#     endif
#  endif
#endif


  //_______________________________________________________________ FFTComplex
  FFTComplex::FFTComplex():myFFTInternal(NULL){
    myFFTInternal = new FFTInternal();
  }

  FFTComplex::~FFTComplex(){
    delete myFFTInternal;
  }

  void FFTComplex::initialize(int x, int y, int z,zomplex* a){
    myFFTInternal->initialize(x,y,z,a);
  }

  void FFTComplex::backward(){
    myFFTInternal->backward();
  }

  void FFTComplex::forward(){
    myFFTInternal->forward();
  }

#if defined(HAVE_FFT_FFTW2_MPI)
  void FFTComplex::FFTComplexMPIInit(int master){
    FFTInternal::FFTInternalMPIInit(master);
  }

  void FFTComplex::FFTComplexMPIFinalize(){
    FFTInternal::FFTInternalMPIFinalize();
  }    

  bool FFTComplex::isParallel(){
    return true;
  }
#else
  void FFTComplex::FFTComplexMPIInit(int){}

  void FFTComplex::FFTComplexMPIFinalize(){}    

  bool FFTComplex::isParallel(){
    return false;
  }
#endif
}
