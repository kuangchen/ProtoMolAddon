/*  -*- c++ -*-  */
#ifndef GRID_H
#define GRID_H

#include <protomol/type/Array.h>
#include <protomol/type/Vector3D.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/parallel/FFTComplex.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>

namespace ProtoMol {
  //_________________________________________________________________ Grid
  /**
   * A simple Grid class using T as interpolation scheme, assuming periodic 
   * boundary conditions.
   */


  template<class TInterpolation>
  class Grid {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef & const
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    struct Int3D {int x; int y; int z;};
    /// 3d interpolation
    struct Interpolation3D {
      TInterpolation x;
      TInterpolation y;
      TInterpolation z;
      Interpolation3D(){};
      Interpolation3D(unsigned int order):x(order),y(order),z(order){}
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Grid();
    ~Grid();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Grid
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:


    void anterpolateCharge(Real q, const Vector3D& coord, unsigned int index);
    void fftBack(){myFFT.backward();}
    Real scalarSum(ScalarStructure* energies);
    Real scalarSum(ScalarStructure* energies, unsigned int block, unsigned int nBlocks);
    void fftForward(){myFFT.forward();}
    void interpolateForce(Real q, unsigned int index, Vector3D& force);

    void initialize(Real width, Real length, Real height, Real alpha,
		    unsigned int nx, unsigned int ny, unsigned int nz, 
		    unsigned int interOrder,
		    unsigned int atomCount);

    void clear();
    void getQ(Real*& begin, Real*& end) {begin=&(myQ.begin()->re);end=&(myQ.end()->re);}
    void print();
  private:
    void dftmod(unsigned int order, unsigned int n, Real* interpolation, Real* interpolationMod);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Array<zomplex,3> myQ;  
    Array<zomplex,3> myQTmp;  
    unsigned int myNX;
    unsigned int myNY;
    unsigned int myNZ;
    int myNXOffset;
    int myNYOffset;
    int myNZOffset;
    Real myWidth;
    Real myLength;
    Real myHeight;
    Real myWidthr;
    Real myLengthr;
    Real myHeightr;
    Real myV;
    Real myHX;
    Real myHY;
    Real myHZ;
    Real myHXr;
    Real myHYr;
    Real myHZr;
    Real myAlpha;
    std::vector<Int3D> myScaledParticleIntPositions;
    std::vector<Interpolation3D> myInterpolations;
    Real* myInerpolationModX;
    Real* myInerpolationModY;
    Real* myInerpolationModZ;
    Real* myExpX;
    Real* myExpY;
    Real* myExpZ;
    unsigned int myInterOrder;
    FFTComplex myFFT;
    unsigned int myAtomCount;
    Real myFac;
  };
  //______________________________________________________________________ INLINES
  template<class TInterpolation>
  inline void Grid<TInterpolation>::anterpolateCharge(Real q, const Vector3D& coord, unsigned int index)  {
    Real x = coord.c[0]*myHXr;
    Real y = coord.c[1]*myHYr;
    Real z = coord.c[2]*myHZr;
    while(x < 0.0) x += myNX;
    while(x >= myNX) x -= myNX;
    while(y < 0.0) y += myNY;
    while(y >= myNY) y -= myNY;
    while(z < 0.0) z += myNZ;
    while(z >= myNZ) z -= myNZ;
    int intX = (int)x;
    int intY = (int)y;
    int intZ = (int)z;
    int i0 = intX+myNXOffset;
    int j0 = intY+myNYOffset;
    int k0 = intZ+myNZOffset;
    myScaledParticleIntPositions[index].x = i0;
    myScaledParticleIntPositions[index].y = j0;
    myScaledParticleIntPositions[index].z = k0;
    myInterpolations[index].x.set(x-intX);
    myInterpolations[index].y.set(y-intY);
    myInterpolations[index].z.set(z-intZ);
    Real*  thetaX = myInterpolations[index].x.theta;
    Real*  thetaY = myInterpolations[index].y.theta;
    Real*  thetaZ = myInterpolations[index].z.theta;
    for(unsigned int i=0;i<myInterOrder;i++){
      Real a = q*thetaX[i];
      int i1 = (i+i0) % myNX; 
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<zomplex,3>::RefArray<2> rQX = myQ[i1];
#else
      RefArray<zomplex,2> rQX = myQ[i1];
#endif
      for(unsigned int j=0;j<myInterOrder;j++){
	Real ab = a*thetaY[j];
	int j1 = (j+j0) % myNY; 
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<zomplex,3>::RefArray<1> rQXY = rQX[j1];
#else
	RefArray<zomplex,1> rQXY = rQX[j1];
#endif
	for(unsigned int k=0;k<myInterOrder;k++){
	  int k1 = (k+k0) % myNZ; 
	  rQXY[k1].re += ab*thetaZ[k];
	}
      }
    }
  }



  template<class TInterpolation>
  inline void Grid<TInterpolation>::interpolateForce(Real q, unsigned int index, Vector3D& force){
    int i0 = myScaledParticleIntPositions[index].x;
    int j0 = myScaledParticleIntPositions[index].y;
    int k0 = myScaledParticleIntPositions[index].z;
    Real fx = 0.0;
    Real fy = 0.0;
    Real fz = 0.0;
    Real*  thetaX = myInterpolations[index].x.theta;
    Real*  thetaY = myInterpolations[index].y.theta;
    Real*  thetaZ = myInterpolations[index].z.theta;
    Real*  dThetaX = myInterpolations[index].x.dTheta;
    Real*  dThetaY = myInterpolations[index].y.dTheta;
    Real*  dThetaZ = myInterpolations[index].z.dTheta;
    for(unsigned int i=0;i<myInterOrder;i++){
      int i1 = (i+i0)%myNX; 
      for(unsigned int j=0;j<myInterOrder;j++){
	int j1 = (j+j0)%myNY;
	Real xij = dThetaX[i]*thetaY[j];
	Real yij = thetaX[i]*dThetaY[j];
	Real zij = thetaX[i]*thetaY[j];
	for(unsigned int k=0;k<myInterOrder;k++){
	  int k1 = (k+k0)%myNZ;
	  Real term = myQ[i1][j1][k1].re;
	  fx -= term*xij*thetaZ[k];
	  fy -= term*yij*thetaZ[k];
	  fz -= term*zij*dThetaZ[k];
	}
      }
    }
    force += Vector3D(fx*myHXr*q,fy*myHYr*q,fz*myHZr*q);
  }

  template<class TInterpolation>
  Grid<TInterpolation>::Grid(): 
    myNX(0), 
    myNY(0), 
    myNZ(0),
    myNXOffset(0), 
    myNYOffset(0), 
    myNZOffset(0),
    myWidth(0.0),
    myLength(0.0),
    myHeight(0.0),
    myWidthr(0.0),
    myLengthr(0.0),
    myHeightr(0.0),
    myV(0.0),
    myHX(0.0),
    myHY(0.0),
    myHZ(0.0),
    myHXr(0.0),
    myHYr(0.0),
    myHZr(0.0),
    myAlpha(0.0),
    myInerpolationModX(NULL),
    myInerpolationModY(NULL),
    myInerpolationModZ(NULL),
    myExpX(NULL),
    myExpY(NULL),
    myExpZ(NULL),
    myInterOrder(0),
    myAtomCount(0),
    myFac(0.0){
  }

  template<class TInterpolation>
  Grid<TInterpolation>::~Grid(){
    if(myInerpolationModX != NULL) delete [] myInerpolationModX;
    if(myInerpolationModY != NULL) delete [] myInerpolationModY;
    if(myInerpolationModZ != NULL) delete [] myInerpolationModZ;
    if(myExpX != NULL) delete [] myExpX;
    if(myExpY != NULL) delete [] myExpY;
    if(myExpZ != NULL) delete [] myExpZ;
  }



  template<class TInterpolation>
  Real Grid<TInterpolation>::scalarSum(ScalarStructure* energies){
    Real energy   = 0.0;
    Real virialxx = 0.0;
    Real virialxy = 0.0;
    Real virialxz = 0.0;
    Real virialyy = 0.0;
    Real virialyz = 0.0;
    Real virialzz = 0.0;
    bool doVirial = energies->virial();
    bool doMolVirial = energies->molecularVirial();
    Real piVr = 1.0/(M_PI*myV);
    int count = 0;
    for (unsigned int i = 0; i < myNX; i++){
      int i0 = i <= myNX/2 ? i : i-myNX;
      Real mi = i0*myWidthr;
      Real ex = myExpX[i]*piVr;
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<zomplex,3>::RefArray<2> rQX = myQ[i];
#else
      RefArray<zomplex,2> rQX = myQ[i];
#endif
      for (unsigned int j = 0 ; j < myNY; j++){
	int j0 = j <= myNY/2 ? j : j-myNY;
	Real interpolationModXY = myInerpolationModX[i]*myInerpolationModY[j];
	Real mj = j0*myLengthr;
	Real mij = mi*mi + mj*mj;
	Real exy = ex*myExpY[j];
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<zomplex,3>::RefArray<1> rQXY = rQX[j];
#else
	RefArray<zomplex,1> rQXY = rQX[j];
#endif
	for (unsigned int k = (i!=0 || j!=0 ? 0:1); k <= myNZ/2; k++){
	  count++;
	  int k0 = k <= myNZ/2 ? k : k-myNZ;
	  //Vector3D mHat(i0*myWidthr,j0*myLengthr,k0*myHeightr);
	  //Real mHatSquared = mHat.normSquared();
	  //Real theta = myInerpolationModX[i]*myInerpolationModY[j]*myInerpolationModZ[k]*exp(-fac*mHatSquared)/(mHatSquared*M_PI*myV);
	  Real mk = k0*myHeightr;
	  Real mHatSquared = mij+mk*mk;
	  Real theta = interpolationModXY*myInerpolationModZ[k]*exy*myExpZ[k]/mHatSquared;
        
	  // Energy
	  Real q = power<2>(rQXY[k].re)+power<2>(rQXY[k].im);
	  //Real q = power<2>(myQ[i][j][k].re)+power<2>(myQ[i][j][k].im);
	  Real e = q*theta;
	  Real v = 2.0*(1.0/mHatSquared + myFac);

	  // Symmetric 
	  if(k > 0 && ((k != myNZ/2) || (myNZ & 1))){
	    e *= 2.0;
	    zomplex& w = myQ[(myNX-i)%myNX][(myNY-j)%myNY][(myNZ-k)%myNZ];
	    w.re *= theta;
	    w.im *= theta;         
	  }

	  energy += e;

	  // Virial
	  if(doVirial){
	    virialxx += e *(1.0 - v * mi*mi);
	    virialxy -= e *(v * mi*mj);
	    virialxz -= e *(v * mi*mk);
	    virialyy += e *(1.0 - v * mj*mj);
	    virialyz -= e *(v * mj*mk);
	    virialzz += e *(1.0 - v * mk*mk);
	  }

	  rQXY[k].re *= theta;
	  rQXY[k].im *= theta;
	  //myQ[i][j][k].re *= theta;
	  //myQ[i][j][k].im *= theta;
	}
      }
    }
    // Just clear (0,0,0) since we did this not the nested loop.
    myQ[0][0][0].re = 0.0;
    myQ[0][0][0].im = 0.0;

    if(doVirial){
      (*energies)[ScalarStructure::VIRIALXX] += virialxx*0.5;
      (*energies)[ScalarStructure::VIRIALXY] += virialxy*0.5;
      (*energies)[ScalarStructure::VIRIALXZ] += virialxz*0.5;
      (*energies)[ScalarStructure::VIRIALYX] += virialxy*0.5;
      (*energies)[ScalarStructure::VIRIALYY] += virialyy*0.5;
      (*energies)[ScalarStructure::VIRIALYZ] += virialyz*0.5;
      (*energies)[ScalarStructure::VIRIALZX] += virialxz*0.5;
      (*energies)[ScalarStructure::VIRIALZY] += virialyz*0.5;
      (*energies)[ScalarStructure::VIRIALZZ] += virialzz*0.5;
    }

    // Molecular Virial
    if(doMolVirial) {
      (*energies)[ScalarStructure::MOLVIRIALXX] += virialxx*0.5;
      (*energies)[ScalarStructure::MOLVIRIALXY] += virialxy*0.5;
      (*energies)[ScalarStructure::MOLVIRIALXZ] += virialxz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALYX] += virialxy*0.5;
      (*energies)[ScalarStructure::MOLVIRIALYY] += virialyy*0.5;
      (*energies)[ScalarStructure::MOLVIRIALYZ] += virialyz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALZX] += virialxz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALZY] += virialyz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALZZ] += virialzz*0.5;
    }
    return energy*0.5;
  }

  template<class TInterpolation>
  Real Grid<TInterpolation>::scalarSum(ScalarStructure* energies, unsigned int block, unsigned int nBlocks){
    Real energy   = 0.0;
    Real virialxx = 0.0;
    Real virialxy = 0.0;
    Real virialxz = 0.0;
    Real virialyy = 0.0;
    Real virialyz = 0.0;
    Real virialzz = 0.0;
    bool doVirial = energies->virial();
    bool doMolVirial = energies->molecularVirial();
    Real piVr = 1.0/(M_PI*myV);

    myQTmp.resize(ArraySizes(myNX)(myNY)(myNZ));
    int m = myQ.size();
    zomplex *q = myQ.begin();
    zomplex *t = myQTmp.begin();
    for(int i=0;i<m;i++){
      t[i].re = q[i].re;
      t[i].im = q[i].im;      
      q[i].re = 0.0;
      q[i].im = 0.0;
    }

    int nx = myNX;
    int ny = myNY;
    int nz = myNZ/2+1;

    int nyz = ny*nz;
    int n  = nx*nyz;
    int sn = (n*block)/nBlocks + (block==0?1:0); // Add 1 to skip i,j,k == 0
    int en = (n*(block+1))/nBlocks - 1;

    int count = 0;
    int size = en-sn+1;
    if(size == 0)
      return 0.0;

    int k = sn % nz;
    int j = (sn / nz) % ny;
    int i = (sn / nyz);
    int ez = (en % nz)+1;
    int ey = ((en / nz) % ny)+1;
    int ex = (en / nyz)+1;
    if(j < ey-1)
      ez = nz;
    if(i < ex-1){
      ey = ny;
      ez = nz;
    }

    for (; i < ex; i++,j=0){
      int i0 = i <= static_cast<int>(myNX/2) ? i : i-myNX;
      Real mi = i0*myWidthr;
      Real ex = myExpX[i]*piVr;
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<zomplex,3>::RefArray<2> rQX  = myQ[i];
      Array<zomplex,3>::RefArray<2> rTQX = myQTmp[i];
#else
      RefArray<zomplex,2> rQX  = myQ[i];
      RefArray<zomplex,2> rTQX = myQTmp[i];
#endif
      for (; j < ey; j++,k=0){
	int j0 = j <= static_cast<int>(myNY/2) ? j : j-myNY;
	Real interpolationModXY = myInerpolationModX[i]*myInerpolationModY[j];
	Real mj = j0*myLengthr;
	Real mij = mi*mi + mj*mj;
	Real exy = ex*myExpY[j];
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<zomplex,3>::RefArray<1> rQXY  = rQX[j];
	Array<zomplex,3>::RefArray<1> rTQXY = rTQX[j];
#else
	RefArray<zomplex,1> rQXY  = rQX[j];
	RefArray<zomplex,1> rTQXY = rTQX[j];
#endif
	//for (unsigned int k = (i!=0 || j!=0 ? 0:1); k < myNZ; k++){
	for (; k < ez; k++){
	  int k0 = k <= static_cast<int>(myNZ/2) ? k : k-myNZ;
	  //Vector3D mHat(i0*myWidthr,j0*myLengthr,k0*myHeightr);
	  //Real mHatSquared = mHat.normSquared();
	  //Real theta = myInerpolationModX[i]*myInerpolationModY[j]*myInerpolationModZ[k]*exp(-fac*mHatSquared)/(mHatSquared*M_PI*myV);
	  Real mk = k0*myHeightr;
	  Real mHatSquared = mij+mk*mk;
	  Real theta = interpolationModXY*myInerpolationModZ[k]*exy*myExpZ[k]/mHatSquared;

	  // Energy
	  Real q = power<2>(rTQXY[k].re)+power<2>(rTQXY[k].im);
	  //Real q = power<2>(myQ[i][j][k].re)+power<2>(myQ[i][j][k].im);
	  Real e = q*theta;
	  Real v = 2.0*(1.0/mHatSquared + myFac);

	  // Symmetric 
	  if(k > 0 && ((k != static_cast<int>(myNZ/2)) || (myNZ & 1))){
	    e *= 2.0;
	    myQ[(myNX-i)%myNX][(myNY-j)%myNY][(myNZ-k)%myNZ].re = myQTmp[(myNX-i)%myNX][(myNY-j)%myNY][(myNZ-k)%myNZ].re * theta;
	    myQ[(myNX-i)%myNX][(myNY-j)%myNY][(myNZ-k)%myNZ].im = myQTmp[(myNX-i)%myNX][(myNY-j)%myNY][(myNZ-k)%myNZ].im * theta;         
	  }

	  energy += e;

	  // Virial
	  if(doVirial){
	    virialxx += e *(1.0 - v * mi*mi);
	    virialxy -= e *(v * mi*mj);
	    virialxz -= e *(v * mi*mk);
	    virialyy += e *(1.0 - v * mj*mj);
	    virialyz -= e *(v * mj*mk);
	    virialzz += e *(1.0 - v * mk*mk);
	  }
	  //
	  rQXY[k].re = rTQXY[k].re * theta;
	  rQXY[k].im = rTQXY[k].im * theta;
	  //myQ[i][j][k].re *= theta;
	  //myQ[i][j][k].im *= theta;

	  count++;
	  if(count >= size){
	    i = myNX;
	    j = myNY;
	    k = myNZ;
	  }
	}
      }
    }
    // Just clear (0,0,0) since we did this not the nested loop.
    myQ[0][0][0].re = 0.0;
    myQ[0][0][0].im = 0.0;

    if(doVirial){
      (*energies)[ScalarStructure::VIRIALXX] += virialxx*0.5;
      (*energies)[ScalarStructure::VIRIALXY] += virialxy*0.5;
      (*energies)[ScalarStructure::VIRIALXZ] += virialxz*0.5;
      (*energies)[ScalarStructure::VIRIALYX] += virialxy*0.5;
      (*energies)[ScalarStructure::VIRIALYY] += virialyy*0.5;
      (*energies)[ScalarStructure::VIRIALYZ] += virialyz*0.5;
      (*energies)[ScalarStructure::VIRIALZX] += virialxz*0.5;
      (*energies)[ScalarStructure::VIRIALZY] += virialyz*0.5;
      (*energies)[ScalarStructure::VIRIALZZ] += virialzz*0.5;
    }

    // Molecular Virial
    if(doMolVirial) {
      (*energies)[ScalarStructure::MOLVIRIALXX] += virialxx*0.5;
      (*energies)[ScalarStructure::MOLVIRIALXY] += virialxy*0.5;
      (*energies)[ScalarStructure::MOLVIRIALXZ] += virialxz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALYX] += virialxy*0.5;
      (*energies)[ScalarStructure::MOLVIRIALYY] += virialyy*0.5;
      (*energies)[ScalarStructure::MOLVIRIALYZ] += virialyz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALZX] += virialxz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALZY] += virialyz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALZZ] += virialzz*0.5;
    }
    return energy*0.5;
  }

  template<class TInterpolation>
  void Grid<TInterpolation>::initialize(Real width, Real length, Real height, Real alpha,
					unsigned int nx, unsigned int ny, unsigned int nz, 
					unsigned int interOrder,
					unsigned int atomCount){
    if(!myQ.resize(ArraySizes(nx)(ny)(nz)))
      report << error <<"[Grid<>::initialize] Could not allocate memory for Q["
	     <<nx<<"]["<<ny<<"]["<<nz<<"]."<<endr;
    myNX = nx; 
    myNY = ny;
    myNZ = nz;
    myNXOffset = -((int)interOrder-1)/2 + (int)nx;
    myNYOffset = -((int)interOrder-1)/2 + (int)ny;
    myNZOffset = -((int)interOrder-1)/2 + (int)nz;
    myWidth = width; 
    myLength = length; 
    myHeight = height;
    myWidthr = 1.0/width; 
    myLengthr = 1.0/length; 
    myHeightr = 1.0/height;
    myV = width*length*height;
    myHX = width/nx;
    myHY = length/ny;
    myHZ = height/nz;
    myHXr = nx/width;
    myHYr = ny/length;
    myHZr = nz/height;
    myAlpha = alpha;
    myInterOrder = interOrder;
    myAtomCount = atomCount;
    myFac =  M_PI*M_PI/(myAlpha*myAlpha);

    // Precompute exp(-pi^2m^2/alpha^2)
    if(myExpX != NULL) 
      delete [] myExpX;
    myExpX = new Real[nx];
    if(myExpY != NULL)
      delete [] myExpY;
    myExpY = new Real[ny];
    if(myExpZ != NULL) 
      delete [] myExpZ;
    myExpZ = new Real[nz];

    for (unsigned int i = 0; i < myNX; i++){
      int i0 = i <= myNX/2 ? i : i-myNX;
      myExpX[i] = exp(-myFac*power<2>(i0*myWidthr));
    }
    for (unsigned int j = 0 ; j < myNY; j++){
      int j0 = j <= myNY/2 ? j : j-myNY;
      myExpY[j] = exp(-myFac*power<2>(j0*myLengthr));
    }

    //for (unsigned int k = 0; k < myNZ; k++){
    for (unsigned int k = 0; k<= myNZ/2 ; k++){
      int k0 = k <= myNZ/2 ? k : k-myNZ;
      myExpZ[k] = exp(-myFac*power<2>(k0*myHeightr));
    }


    // Precompute the mod TInterpolation, B(m1,m2,m3)
    if(myInerpolationModX != NULL) 
      delete [] myInerpolationModX;
    myInerpolationModX = new Real[nx];
    if(myInerpolationModY != NULL)
      delete [] myInerpolationModY;
    myInerpolationModY = new Real[ny];
    if(myInerpolationModZ != NULL) 
      delete [] myInerpolationModZ;
    myInerpolationModZ = new Real[nz];

    TInterpolation interpolation = TInterpolation(myInterOrder,0.0);
    dftmod(myInterOrder,nx,interpolation.theta,myInerpolationModX);
    dftmod(myInterOrder,ny,interpolation.theta,myInerpolationModY);
    dftmod(myInterOrder,nz,interpolation.theta,myInerpolationModZ);
    //for(unsigned int i=0;i<nx;i++)
    //  report << plain <<"myInerpolationModX["<<i<<"]:"<<1.0/myInerpolationModX[i]<<endr;
    //for(unsigned int i=0;i<ny;i++)
    //  report << plain <<"myInerpolationModY["<<i<<"]:"<<1.0/myInerpolationModY[i]<<endr;
    //for(unsigned int i=0;i<nz;i++)
    //  report << plain <<"myInerpolationModZ["<<i<<"]:"<<1.0/myInerpolationModZ[i]<<endr;



    // Resize the vector data members
    myScaledParticleIntPositions.resize(myAtomCount);
    myInterpolations.resize(myAtomCount,Interpolation3D(myInterOrder));

    myFFT.initialize(myNX,myNY,myNZ,&myQ[0][0][0]);
  }

  template<class TInterpolation>
  void Grid<TInterpolation>::dftmod(unsigned int order, unsigned int n, Real* interpolation, Real* interpolationMod){
    for(unsigned int i=0;i<n;i++){
      Real sumCos = 0.0;
      Real sumSin = 0.0;
      for(unsigned int j=0;j<order;j++){
	Real x = M_PI*2.0*i*j/(Real)n;
	sumCos += interpolation[j]*cos(x);
	sumSin += interpolation[j]*sin(x);
      }
      interpolationMod[i] = 1.0/(sumCos*sumCos + sumSin*sumSin);
    }
  }

  template<class TInterpolation>
  void Grid<TInterpolation>::clear(){
    int n = myQ.size();
    zomplex *q = myQ.begin();
    for(int i=0;i<n;i++){
      q[i].re = 0.0;
      q[i].im = 0.0;
    }
  }

  template<class TInterpolation>
  void Grid<TInterpolation>::print(){
    Real q = 0.0;
    report << plain;
    for (unsigned int i = 0; i < myNX; i++){
      for (unsigned int j = 0 ; j < myNY; j++){
	report << "Q["<<i<<"]["<<j<<"][0-"<<myNZ-1<<"] : ";
	for (unsigned int k = 0; k < myNZ; k++){
	  report << "("<<myQ[i][j][k].re <<","<<myQ[i][j][k].im<<")";
	  q += myQ[i][j][k].re;
	}
	report <<std::endl;
      }
    }
    report <<"Sum Q :"<<q<<endr;
  }
}
#endif
