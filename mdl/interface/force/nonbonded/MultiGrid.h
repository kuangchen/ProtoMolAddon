/*  -*- c++ -*-  */
#ifndef MULTIGRID_H
#define MULTIGRID_H

#include <protomol/type/Array.h>
#include <protomol/type/Vector3D.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/base/Report.h>

//#define DEBUG_MULTIGRID
//#define DEBUG_MULTIGRID_TIMING

namespace ProtoMol {
  //_________________________________________________________________ MultiGrid
  static const Real BORDER_TOLERANCE = 0.0001;
  static const int BORDER = 2;


  /**
   * Multi grid algorithm based on TInterpolation scheme using TKernel
   * handling non-periodic and periodic boundary conditions.
   */
  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  class MultiGrid {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef & const
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    struct Int3D  {int x, y, z;};
    /// 3D interpolation
    struct Interpolation3D {
      TInterpolation x;
      TInterpolation y;
      TInterpolation z;
      Interpolation3D(){};
      Interpolation3D(unsigned int order):x(order),y(order),z(order){};
    };
  

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    MultiGrid();
#ifdef DEBUG_MULTIGRID_TIMING
    ~MultiGrid(){
      if(myCounterDirect+myCounterCorrection > 0)
	Report::report << allnodes << plain 
	       << "[MultiGrid] Smooth part: direct="
	       << myCounterDirect<<", correction="
	       << myCounterCorrection<<"."<<Report::endr;
    }
    long getCounter(){return (myCounterDirect + myCounterCorrection);}
#endif
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class MultiGrid
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    /// Interpolation weights for the position and  anterpolate the charge q onto grid level 0, q -> Q(0)
    void anterpolateCharge(Real q, const Vector3D& coord, unsigned int index);

    /// Anterpolates Q from level to level+1, Q(level) -> Q(level+1)
    void fineToCoarse(int level, unsigned int block, unsigned int nBlocks);
    void fineToCoarse(int level){fineToCoarse(level,0,1);}

    /// Computes the potential for the top level, Q(maxLevels-1) -> V(maxLevels-1)
    void direct(unsigned int block, unsigned int nBlocks);
    void direct(){direct(0,1);};

    /// Interpolates Q from level to level-1, Q(level) -> Q(level-1)
    void coarseToFine(int level, unsigned int block, unsigned int nBlocks);
    void coarseToFine(int level){coarseToFine(level,0,1);}

    /// Adds the correction term Q(level) -> V(level)
    void correction(int level, unsigned int block, unsigned int nBlocks);
    void correction(int level){correction(level,0,1);};

    /// Computes the energy of a given level, V(level) -> energy
    Real energy(int level, unsigned int block, unsigned int nBlocks);
    Real energy(int level){return energy(level,0,1);};

    /// Interpolates the force from grid, level 0, Q(0) -> force 
    void interpolateForce(Real q, unsigned int index, Vector3D& force);


    void initialize(unsigned int n, Real s, int levels,
		    int nx, int ny, int nz, 
		    int interOrder, int ratio, 
		    Vector3D min, Vector3D max, 
		    Vector3D h, Vector3D origin); 

    /** Update the length, width and height of the grid according min and max such
     * that anter-/interpolated works regardless order and boundary conditions.
     * In case of non-PBC it adds a border such that all points inside [min,max]
     * are correctly anter-/interpolated.
     */
    void updateSize(Vector3D min, Vector3D max);

    void getV(int level, Real*& begin, Real*& end) {begin=myV[level].begin();end=myV[level].end();}
    void getQ(int level, Real*& begin, Real*& end) {begin=myQ[level].begin();end=myQ[level].end();}
    void clear();

    Real sumV(int level);
    Real sumQ(int level);
    void print(int level){printQ(level);printV(level);}
    void printV(int level);
    void printQ(int level);
    void printConst();
  private:
    void precomputeG();
    void blocksCorrection(unsigned int nBlocks);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    std::vector<Array<Real,3> > myQ;  ///< Arrays of charges
    std::vector<Array<Real,3> > myV;  ///< Arrays of potentials
    Array<Real,3> myGDirect;     ///< G_direct for top level, G->V
    Array<Real,3> myGCorrection; ///< G_Correction for other levels, G->V   
    int  myGCorrDimX;            ///< Dimensions of G_correction
    int  myGCorrDimY;
    int  myGCorrDimZ;
    Vector3D myMin;           // 
    Vector3D myMax;           // 
    Vector3D myCMin;          // 
    Vector3D myCMax;          // 
    Vector3D myH;             // 
    Vector3D myOrigin;        // 
    int myNX;                    ///< x dimension of the grid level 0
    int myNY;                    ///< y dimension of the grid level 0
    int myNZ;                    ///< z dimension of the grid level 0
    int myNXOffset;              ///< x index offset between to grid levels
    int myNYOffset;              ///< y index offset between to grid levels
    int myNZOffset;              ///< z index offset between to grid levels
    Real myHX;                   ///< h_x at level 0
    Real myHY;                   ///< h_y at level 0
    Real myHZ;                   ///< h_z at level 0
    Real myHXr;                  ///< 1/h_x at level 0
    Real myHYr;                  ///< 1/h_y at level 0
    Real myHZr;                  ///< 1/h_z at level 0
    std::vector<Int3D> myScaledParticleIntPositions;  ///< Starting index for anter/-interpolation
    // for each particle
    std::vector<Interpolation3D> myInterpolations;    ///< Interpolation weights for each particle 
    std::vector<TInterpolation> myGridInterpolation;  ///< Interpolation between grids
    Real myS;                    ///< Softening distance
    Real myRS;                   ///< 1 / Softening distance
    int myInterOrder;            ///< Interpolation order
    int myRatio;                 ///<  Ratio between to levels, usually 2
    int myLevels;                ///<  Number of levels, 0,1,2,...,myLevels-1
    int mySigma;                 ///< 1 if interpolation weights ...,0 for w=0, 0 else
    std::vector<Int3D> myDim;    ///< Dimension of the grids
    std::vector<Real> myScale;   ///< Scaling factor for each level
    std::vector<std::vector<int> > myBlocksCorrection; ///< List of blocks for the correction
#ifdef DEBUG_MULTIGRID_TIMING
    long myCounterDirect;
    long myCounterCorrection;
#endif
  };
  //______________________________________________________________________ INLINES

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  inline void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::anterpolateCharge(Real q, const Vector3D& coord, unsigned int index)  {
    Real x = (coord.c[0]-myMin.c[0])*myHXr;
    Real y = (coord.c[1]-myMin.c[1])*myHYr;
    Real z = (coord.c[2]-myMin.c[2])*myHZr;
    if(pbcX) while(x < 0.0) x += myNX;
    if(pbcX) while(x >= myNX) x -= myNX;
    if(pbcY) while(y < 0.0) y += myNY;
    if(pbcY) while(y >= myNY) y -= myNY;
    if(pbcZ) while(z < 0.0) z += myNZ;
    if(pbcZ) while(z >= myNZ) z -= myNZ;
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
    Array<Real,3>& rQ = myQ[0];
    for(int i=0;i<myInterOrder;i++){
      Real a = q*thetaX[i];
      int i1 = i+i0;
      if(pbcX) i1 = i1 % myNX; 
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<Real,3>::RefArray<2> rQX = rQ[i1];
#else
      RefArray<Real,2> rQX = rQ[i1];
#endif
      for(int j=0;j<myInterOrder;j++){
	Real ab = a*thetaY[j];
	int j1 = j+j0;
	if(pbcY) j1 = j1 % myNY; 
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<Real,3>::RefArray<1> rQXY = rQX[j1];
#else
	RefArray<Real,1> rQXY = rQX[j1];
#endif
	for(int k=0;k<myInterOrder;k++){
	  int k1 = k+k0;
	  if(pbcZ) k1 = k1 % myNZ; 
	  rQXY[k1] += ab*thetaZ[k];
	}
      }
    }
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  inline void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::fineToCoarse(int level, unsigned int block, unsigned int nBlocks){
    int nx = myDim[level].x;
    int ny = myDim[level].y;
    int nz = myDim[level].z;
    int nx2 = myDim[level+1].x;
    int ny2 = myDim[level+1].y;
    int nz2 = myDim[level+1].z;
    Array<Real,3>& rQ0 = myQ[level];
    Array<Real,3>& rQ1 = myQ[level+1];

    int nyz = ny*nz;
    int n  = nx*nyz;
    int sn = (n*block)/nBlocks;
    int en = (n*(block+1))/nBlocks - 1;
    int count = 0;
    int size = en-sn+1;
    if(size == 0)
      return;

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
      int i1 = i/myRatio;
      Real*  thetaX = myGridInterpolation[i%myRatio].theta;
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<Real,3>::RefArray<2> rQ0X = rQ0[i];
#else
      RefArray<Real,2> rQ0X = rQ0[i];
#endif
      for (; j < ey; j++,k=0){
	int j1 = j/myRatio;
	Real*  thetaY = myGridInterpolation[j%myRatio].theta;
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<Real,3>::RefArray<1> rQ0XY = rQ0X[j];
#else
	RefArray<Real,1> rQ0XY = rQ0X[j];
#endif
	for (; k < ez; k++){
	  int k1 = k/myRatio;
	  Real*  thetaZ = myGridInterpolation[k%myRatio].theta;
	  Real q = rQ0XY[k];
	  for(int i0=0;i0<myInterOrder;i0++){
	    int i2 = i1+i0;
	    if(pbcX) i2 = (myNXOffset+i2) % nx2;
	    Real qx = q*thetaX[i0];
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	    Array<Real,3>::RefArray<2> rQ1X = rQ1[i2];
#else
	    RefArray<Real,2> rQ1X = rQ1[i2];
#endif
	    for(int j0=0;j0<myInterOrder;j0++){
	      int j2 = j1+j0;
	      if(pbcY) j2 = (myNYOffset+j2) % ny2; 
	      Real qxy = qx*thetaY[j0];
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	      Array<Real,3>::RefArray<1> rQ1XY = rQ1X[j2];
#else
	      RefArray<Real,1> rQ1XY = rQ1X[j2];
#endif
	      for(int k0=0;k0<myInterOrder;k0++){
		int k2 = k1+k0;
		if(pbcZ) k2 = (myNZOffset+k2) % nz2;
		rQ1XY[k2] += qxy*thetaZ[k0];
	      }
	    }
	  }
	  count++;
	  if(count >= size)
	    return;        
	}
      }
    }
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  inline void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::direct(unsigned int block, unsigned int nBlocks){
    int nx = myDim[myLevels-1].x;
    int ny = myDim[myLevels-1].y;
    int nz = myDim[myLevels-1].z;
    Array<Real,3>& rV = myV[myLevels-1];
    Array<Real,3>& rQ = myQ[myLevels-1];
    Real g0 = myGDirect[0][0][0];
    int nyz = ny*nz;
    int n  = nx*nyz;
    int sn = 0;
    int en = n -1;
    if(block > 0)
      sn = (int)((-1.0+2.0*n-sqrt(power<2>(1.0-2.0*n)-4.0*n*(n-1.0)*(block)/(Real)nBlocks))/2.0+0.5);
    if(block < nBlocks-1)
      en = (int)((-1.0+2.0*n-sqrt(power<2>(1.0-2.0*n)-4.0*n*(n-1.0)*(block+1.0)/(Real)nBlocks))/2.0  - 0.5);

    int count = 0;
    int size = en-sn+1;
    if(size == 0)
      return;

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
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<Real,3>::RefArray<2> rV1X  = rV[i];
      Array<Real,3>::RefArray<2> rQ1X  = rQ[i];
#else
      RefArray<Real,2> rV1X  = rV[i];
      RefArray<Real,2> rQ1X  = rQ[i];
#endif
      for (; j < ey; j++,k=0){
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<Real,3>::RefArray<1> rV1XY  = rV1X[j];
	Array<Real,3>::RefArray<1> rQ1XY  = rQ1X[j];
#else
	RefArray<Real,1> rV1XY  = rV1X[j];
	RefArray<Real,1> rQ1XY  = rQ1X[j];
#endif
	for (; k < ez; k++){
	  int l,m,n;
	  Real v = 0.0;
	  Real q = rQ1XY[k];
	  for (l = i, m = j, n = k+1; l < nx; l++, m=0){
	    int i0 = i-l;
	    if(pbcX){
	      i0 = (i0+nx)%nx;
	      i0 = std::min(i0,nx-i0);
	    }
	    else {
	      i0 = i0 < 0 ? -i0 : i0;
	    }
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	    Array<Real,3>::RefArray<2> rGDX = myGDirect[i0];
	    Array<Real,3>::RefArray<2> rVX  = rV[l];
	    Array<Real,3>::RefArray<2> rQX  = rQ[l];
#else
	    RefArray<Real,2> rGDX = myGDirect[i0];
	    RefArray<Real,2> rVX  = rV[l];
	    RefArray<Real,2> rQX  = rQ[l];
#endif
	    for (; m < ny; m++,n=0){
	      int j0 = j-m;
	      if(pbcX){
		j0 = (j0+ny)%ny;
		j0 = std::min(j0,ny-j0);
	      }
	      else {
		j0 = j0 < 0 ? -j0 : j0;
	      }
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	      Array<Real,3>::RefArray<1> rGDXY = rGDX[j0];
	      Array<Real,3>::RefArray<1> rVXY  = rVX[m];
	      Array<Real,3>::RefArray<1> rQXY  = rQX[m];
#else
	      RefArray<Real,1> rGDXY = rGDX[j0];
	      RefArray<Real,1> rVXY  = rVX[m];
	      RefArray<Real,1> rQXY  = rQX[m];
#endif
	      for (; n < nz; n++){
		int k0 = k-n;
		if(pbcZ){
		  k0 = (k0+nz)%nz;
		  k0 = std::min(k0,nz-k0);
		}
		else {
		  k0 = k0 < 0 ? -k0 : k0;
		}
		Real g = rGDXY[k0];
		v += rQXY[n]*g;
		rVXY[n] += q*g;
#ifdef DEBUG_MULTIGRID_TIMING
		myCounterDirect++;
#endif
	      }
	    }
	  }
	  rV1XY[k] += q*g0+v;
	  count++;
	  if(count >= size)
	    return;        
	}
      }
    }
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  inline void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::coarseToFine(int level, unsigned int block, unsigned int nBlocks){
    int nx = myDim[level-1].x;
    int ny = myDim[level-1].y;
    int nz = myDim[level-1].z;
    int nx2 = myDim[level].x;
    int ny2 = myDim[level].y;
    int nz2 = myDim[level].z;
    Array<Real,3>& rV0 = myV[level-1];
    Array<Real,3>& rV1 = myV[level];

    int nyz = ny*nz;
    int n  = nx*nyz;
    int sn = (n*block)/nBlocks;
    int en = (n*(block+1))/nBlocks - 1;
    int count = 0;
    int size = en-sn+1;
    if(size == 0)
      return;

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
      int i1 = i/myRatio;
      Real*  thetaX = myGridInterpolation[i%myRatio].theta;
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<Real,3>::RefArray<2> rV0X  = rV0[i];
#else
      RefArray<Real,2> rV0X  = rV0[i];
#endif
      for (; j < ey; j++,k=0){
	int j1 = j/myRatio;
	Real*  thetaY = myGridInterpolation[j%myRatio].theta;
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<Real,3>::RefArray<1> rV0XY  = rV0X[j];
#else
	RefArray<Real,1> rV0XY  = rV0X[j];
#endif
	for (; k < ez; k++){
	  int k1 = k/myRatio;
	  Real*  thetaZ = myGridInterpolation[k%myRatio].theta;
	  Real v = 0.0;
	  for(int i0=0;i0<myInterOrder;i0++){
	    int i2 = i1+i0;
	    if(pbcX) i2 = (myNXOffset+i2) % nx2;
	    Real x = thetaX[i0];
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	    Array<Real,3>::RefArray<2> rV1X  = rV1[i2];
#else
	    RefArray<Real,2> rV1X  = rV1[i2];
#endif
	    for(int j0=0;j0<myInterOrder;j0++){
	      int j2 = j1+j0;
	      if(pbcY) j2 = (myNYOffset+j2) % ny2; 
	      Real xy = x*thetaY[j0];
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	      Array<Real,3>::RefArray<1> rV1XY  = rV1X[j2];
#else
	      RefArray<Real,1> rV1XY  = rV1X[j2];
#endif
	      for(int k0=0;k0<myInterOrder;k0++){
		int k2 = k1+k0;
		if(pbcZ) k2 = (myNZOffset+k2) % nz2;
		v += rV1XY[k2]*xy*thetaZ[k0];
	      }
	    }
	  }
	  rV0XY[k] += v;
	  count++;
	  if(count >= size)
	    return;        
	}
      }
    }
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  inline void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::correction(int level, unsigned int block, unsigned int nBlocks){
    Real scale = 1.0/myScale[level];
    int nx = myDim[level].x;
    int ny = myDim[level].y;
    int nz = myDim[level].z;
    Array<Real,3>& rV = myV[level];
    Array<Real,3>& rQ = myQ[level];

    int nyz = ny*nz;
    int n  = nx*nyz;
    int sn = 0;
    int en = n - 1;
    if(nBlocks > 1){
      if(myBlocksCorrection.empty() || 
	 myBlocksCorrection[level].size() != nBlocks+1 ||
	 myBlocksCorrection[level][myBlocksCorrection[level].size()-1] !=  en +1)
	blocksCorrection(nBlocks);
      sn=myBlocksCorrection[level][block];
      en=myBlocksCorrection[level][block+1]-1;
    }
    int count = 0;
    int size = en-sn+1;
    if(size == 0)
      return;

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
      int hi_l = i+myGCorrDimX;
      int so_l = -myGCorrDimX+1; 
      int lo_l = i+so_l;
      if(!pbcX){ 
	hi_l = std::min(hi_l,nx);
	so_l+=-std::min(lo_l,0);
	lo_l = std::max(lo_l,0);
      }
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<Real,3>::RefArray<2> rVX  = rV[i];
#else
      RefArray<Real,2> rVX  = rV[i];
#endif
      for (; j < ey; j++,k=0){
	int hi_m = j+myGCorrDimY;
	int so_m = -myGCorrDimY+1; 
	int lo_m = j+so_m;
	if(!pbcY){
	  hi_m = std::min(hi_m,ny);
	  so_m+=-std::min(lo_m,0);
	  lo_m = std::max(lo_m,0);
	}
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<Real,3>::RefArray<1> rVXY  = rVX[j];
#else
	RefArray<Real,1> rVXY  = rVX[j];
#endif
	for (; k < ez; k++){
	  int hi_n = k+myGCorrDimZ;
	  int so_n = -myGCorrDimZ+1; 
	  int lo_n = k+so_n;
	  if(!pbcZ){
	    hi_n = std::min(hi_n,nz);
	    so_n+=-std::min(lo_n,0);
	    lo_n = std::max(lo_n,0);
	  }
	  Real v = 0.0;
	  for (int l = lo_l, l2 = so_l; l < hi_l; l++,l2++){
	    int l0 = l;
	    if(pbcX) l0 = (l0+nx)%nx;
	    int id = std::max(-l2,l2);
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	    Array<Real,3>::RefArray<2> rQX  = rQ[l0];
	    Array<Real,3>::RefArray<2> rGCX = myGCorrection[id];
#else
	    RefArray<Real,2> rQX  = rQ[l0];
	    RefArray<Real,2> rGCX = myGCorrection[id];
#endif
	    for (int m = lo_m, m2 = so_m; m < hi_m; m++,m2++){
	      int m0 = m;
	      if(pbcY) m0 = (m0+ny)%ny;	    
	      int jd = std::max(-m2,m2);
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	      Array<Real,3>::RefArray<1> rQXY  = rQX[m0];
	      Array<Real,3>::RefArray<1> rGCXY = rGCX[jd];
#else
	      RefArray<Real,1> rQXY  = rQX[m0];
	      RefArray<Real,1> rGCXY = rGCX[jd];
#endif
	      for (int n = lo_n, n2 = so_n; n < hi_n; n++,n2++){
		int n0 = n;
		if(pbcZ) n0 = (n0+nz)%nz;
		int kd = std::max(-n2,n2);
		v += rQXY[n0]*rGCXY[kd];
#ifdef DEBUG_MULTIGRID_TIMING
		myCounterCorrection++;
#endif
	      }
	    }
	  }
	  rVXY[k] += v*scale;
	  count++;
	  if(count >= size)
	    return;        
	}
      }
    }
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  inline Real MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::energy(int level, unsigned int block, unsigned int nBlocks){
    Real e = 0.0;
    int nx = myDim[level].x;
    int ny = myDim[level].y;
    int nz = myDim[level].z;

    int nyz = ny*nz;
    int n  = nx*nyz;
    int sn = (n*block)/nBlocks;
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

    Array<Real,3>& rV = myV[level];
    Array<Real,3>& rQ = myQ[level];
    for (; i < ex; i++,j=0){
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<Real,3>::RefArray<2> rVX  = rV[i];
      Array<Real,3>::RefArray<2> rQX  = rQ[i];
#else
      RefArray<Real,2> rVX  = rV[i];
      RefArray<Real,2> rQX  = rQ[i];
#endif
      for (; j < ey; j++,k=0){
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<Real,3>::RefArray<1> rVXY  = rVX[j];
	Array<Real,3>::RefArray<1> rQXY  = rQX[j];
#else
	RefArray<Real,1> rVXY  = rVX[j];
	RefArray<Real,1> rQXY  = rQX[j];
#endif
	for (; k < ez; k++){
	  e += rVXY[k] * rQXY[k];
	  count++;
	  if(count >= size)
	    return 0.5*e;        
	}
      }
    }
    //   Real *q = myQ[level].begin();
    //   Real *v = myV[level].begin();
    //   for(int i=sn;i<en+1;i++)
    //     e += v[i]*q[i];
    return 0.5*e;
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  inline Real MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::sumQ(int level){
    Real e = 0.0;
    Array<Real,3>& rQ = myQ[level];
    for (int i = 0; i < myDim[level].x; i++){
      for (int j = 0 ; j < myDim[level].y; j++){
	for (int k = 0; k < myDim[level].z; k++){
	  e += rQ[i][j][k];
	}
      }
    }
    return e;
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  inline Real MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::sumV(int level){
    Real e = 0.0;
    Array<Real,3>& rV = myV[level];
    for (int i = 0; i < myDim[level].x; i++){
      for (int j = 0 ; j < myDim[level].y; j++){
	for (int k = 0; k < myDim[level].z; k++){
	  e += rV[i][j][k];
	}
      }
    }
    return e;
  }


  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  inline void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::interpolateForce(Real q, unsigned int index, Vector3D& force){
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
    Array<Real,3>& rV = myV[0];
    for(int i=0;i<myInterOrder;i++){
      int i1 = i+i0;
      if(pbcX) i1 = i1 % myNX; 
      Real ax = dThetaX[i];
      Real ay = thetaX[i];
      Real az = thetaX[i];
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<Real,3>::RefArray<2> rVX  = rV[i1];
#else
      RefArray<Real,2> rVX  = rV[i1];
#endif
      for(int j=0;j<myInterOrder;j++){
	int j1 = j+j0;
	if(pbcY) j1 = j1 % myNY; 
	Real abx = ax * thetaY[j];
	Real aby = ay * dThetaY[j];
	Real abz = az * thetaY[j];
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<Real,3>::RefArray<1> rVXY  = rVX[j1];
#else
	RefArray<Real,1> rVXY  = rVX[j1];
#endif
	for(int k=0;k<myInterOrder;k++){
	  int k1 = k+k0;
	  if(pbcZ) k1 = k1 % myNZ; 
	  Real term = rVXY[k1];
	  fx -= term*abx*thetaZ[k];
	  fy -= term*aby*thetaZ[k];
	  fz -= term*abz*dThetaZ[k];
	}
      }
    }
    force += Vector3D(fx*q*myHXr,fy*q*myHYr,fz*q*myHZr);
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  inline void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::precomputeG(){
    if(!myGDirect.resize(ArraySizes(myDim[myLevels-1].x)(myDim[myLevels-1].y)(myDim[myLevels-1].z)))
      Report::report << Report::error 
	     <<"[MultiGrid<>::precomputeG] Could not allocate memory for GDirect["
	   
	     <<myDim[myLevels-1].x<<"]["<<myDim[myLevels-1].y<<"]["<<myDim[myLevels-1].z<<"]."
	     <<Report::endr;
    Real scale = 1.0/myScale[myLevels-1];
    for (int i = 0; i < myDim[myLevels-1].x; i++){
      for (int j = 0 ; j < myDim[myLevels-1].y; j++){
	for (int k = 0; k < myDim[myLevels-1].z; k++){
	  Real r = sqrt(power<2>(i*myHX)+power<2>(j*myHY)+power<2>(k*myHZ));
	  myGDirect[i][j][k] = TKernel::smoothKernel(r,myS,myRS)*scale;
	}
      }
    }

    // No correction if just one level specified.
    if(myLevels < 2 )
      return;
    myGCorrDimX = std::min((int)ceil(myRatio*myS/myHX),myNX);
    myGCorrDimY = std::min((int)ceil(myRatio*myS/myHY),myNY);
    myGCorrDimZ = std::min((int)ceil(myRatio*myS/myHZ),myNZ);
    if(!myGCorrection.resize(ArraySizes(myGCorrDimX)(myGCorrDimY)(myGCorrDimZ)))
      Report::report << Report::error <<"[MultiGrid<>::precomputeG] Could not allocate memory for GCorrection["
	     <<myGCorrDimX<<"]["<<myGCorrDimY<<"]["<<myGCorrDimZ<<"]."<<Report::endr;
    for (int i = 0; i < myGCorrDimX; i++){
      for (int j = 0 ; j < myGCorrDimY; j++){
	for (int k = 0; k < myGCorrDimZ; k++){
	  Real r = sqrt(power<2>(i*myHX)+power<2>(j*myHY)+power<2>(k*myHZ));
	  myGCorrection[i][j][k] = TKernel::smoothKernel(r,myS,myRS)-TKernel::smoothKernel(r,myS*myRatio,myRS/myRatio);
	}
      }
    }
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  inline void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::blocksCorrection(unsigned int nBlocks){

    if(myLevels < 2)
      return;
    myBlocksCorrection.resize(myLevels-1);

    long nxyz = (myGCorrDimX*2-1)*(myGCorrDimY*2-1)*(myGCorrDimZ*2-1);
  
    for(int level=0;level<myLevels-1;level++){
      int nx = myDim[level].x;
      int ny = myDim[level].y;
      int nz = myDim[level].z;
      long size = 0;
      for (int i = 0; i < nx; i++){
	int hi_l = i+myGCorrDimX;
	int so_l = -myGCorrDimX+1; 
	int lo_l = i+so_l;
	if(!pbcX){ 
	  hi_l = std::min(hi_l,nx);
	  lo_l = std::max(lo_l,0);
	}
	for (int j = 0 ; j < ny; j++){
	  int hi_m = j+myGCorrDimY;
	  int so_m = -myGCorrDimY+1; 
	  int lo_m = j+so_m;
	  if(!pbcY){
	    hi_m = std::min(hi_m,ny);
	    lo_m = std::max(lo_m,0);
	  }
	  for (int k = 0; k < nz; k++){
	    int hi_n = k+myGCorrDimZ;
	    int so_n = -myGCorrDimZ+1; 
	    int lo_n = k+so_n;
	    if(!pbcZ){
	      hi_n = std::min(hi_n,nz);
	      lo_n = std::max(lo_n,0);
	    }
	    for (int l = lo_l; l < hi_l; l++){
	      for (int m = lo_m; m < hi_m; m++){
		for (int n = lo_n; n < hi_n; n++){
		  size++;
		}
	      }
	    }
	    size -= nxyz;
	  }
	}
      }
      //Report::report << allnodesserial<< plain << level<<":"<<(Real)nx*(Real)ny*(Real)nz*(Real)nxyz-(Real)size<<","<<size <<Report::endr;

      std::vector<int> blocks;
      blocks.push_back(0);
      long count = 0;
      long block = 0;

      for (int i = 0; i < nx; i++){
	int hi_l = i+myGCorrDimX;
	int so_l = -myGCorrDimX+1; 
	int lo_l = i+so_l;
	if(!pbcX){ 
	  hi_l = std::min(hi_l,nx);
	  lo_l = std::max(lo_l,0);
	}
	for (int j = 0 ; j < ny; j++){
	  int hi_m = j+myGCorrDimY;
	  int so_m = -myGCorrDimY+1; 
	  int lo_m = j+so_m;
	  if(!pbcY){
	    hi_m = std::min(hi_m,ny);
	    lo_m = std::max(lo_m,0);
	  }
	  for (int k = 0; k < nz; k++){
	    int hi_n = k+myGCorrDimZ;
	    int so_n = -myGCorrDimZ+1; 
	    int lo_n = k+so_n;
	    if(!pbcZ){
	      hi_n = std::min(hi_n,nz);
	      lo_n = std::max(lo_n,0);
	    }
	    for (int l = lo_l; l < hi_l; l++){
	      for (int m = lo_m; m < hi_m; m++){
		for (int n = lo_n; n < hi_n; n++){
		  count++;
		}
	      }
	    }
	    count -= nxyz;
	    block++;
	    if((((Real)block*(Real)nxyz-(Real)count)*(Real)nBlocks)/((Real)nx*(Real)ny*(Real)nz*(Real)nxyz-(Real)size) >= (Real)blocks.size())
	      blocks.push_back(block);            
	  }
	}
      }
      block = nx*ny*nz;
      blocks[blocks.size()-1]=block;

      //Report::report << plain << level<<": ("<<blocks.size() <<") ";
      while(blocks.size()<= nBlocks)
	blocks.push_back(block);      
      myBlocksCorrection[level] = blocks;
      //Report::report <<myBlocksCorrection[level].size() <<" :";
      //for(int b=0;b<myBlocksCorrection[level].size();b++)
      //  Report::report << myBlocksCorrection[level][b]<<" ";
      //Report::report <<Report::endr;
    }
  }




  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::MultiGrid(): 
    myQ(0),
    myV(0),
    myGDirect(ArraySizes(0)(0)(0)),
    myGCorrDimX(0), 
    myGCorrDimY(0), 
    myGCorrDimZ(0),
    myMin(Vector3D(Constant::MAXREAL,Constant::MAXREAL,Constant::MAXREAL)),
    myMax(Vector3D(-Constant::MAXREAL,-Constant::MAXREAL,-Constant::MAXREAL)),
    myCMin(Vector3D(Constant::MAXREAL,Constant::MAXREAL,Constant::MAXREAL)),
    myCMax(Vector3D(-Constant::MAXREAL,-Constant::MAXREAL,-Constant::MAXREAL)),
    myH(Vector3D(0.0,0.0,0.0)),
    myOrigin(Vector3D(0.0,0.0,0.0)),
    myNX(0), 
    myNY(0), 
    myNZ(0),
    myNXOffset(0), 
    myNYOffset(0), 
    myNZOffset(0),
    myHX(0.0),
    myHY(0.0),
    myHZ(0.0),
    myHXr(0.0),
    myHYr(0.0),
    myHZr(0.0),
    myScaledParticleIntPositions(0),
    myInterpolations(0),
    myGridInterpolation(0),
    myS(0.0),
    myInterOrder(0),
    myRatio(0),
    myLevels(0),
    mySigma(0)
#ifdef DEBUG_MULTIGRID_TIMING
							       , 
    myCounterDirect(0),
    myCounterCorrection(0)
#endif
  {
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::initialize(unsigned int n, Real s, int levels,
								    int nx, int ny, int nz, 
								    int interOrder, int ratio, 
								    Vector3D min, Vector3D max,
								    Vector3D h, Vector3D origin){
  
    //Report::report << Report::debug << s<<","<<levels <<","<< nx<<","<<ny <<","<<nz<<","<<interOrder<<","<<ratio<<","<<min<<","<<max<<","<<h<<","<<origin<<Report::endr;

    //Report::report.reset();
    //Report::report.setf(std::ios::scientific);
    //Report::report.precision(18);
    // Frist do some fast checks on the input ...
    if(s <= 0.0)
      Report::report << Report::error << "[MultiGrid::initialize] s (="<<s<<") must be > 0."<<Report::endr;

    if(levels < 1)
      Report::report << Report::error << "[MultiGrid::initialize] levels (="<<levels<<") must be > 0."<<Report::endr;

    if(ratio < 2)
      Report::report << Report::error << "[MultiGrid::initialize] ratio (="
	     <<ratio<<") must be > 1."<<Report::endr;

    if(interOrder < 2 || interOrder % 2 != 0)
      Report::report << Report::error << "[MultiGrid::initialize] order (="
	     <<interOrder<<") must be > 1 and even."<<Report::endr;

    if((pbcX && nx < interOrder) || (pbcY && ny < interOrder) || (pbcZ && nz < interOrder)){
      Report::report << Report::error << "[MultiGrid::initialize]";
      if(pbcX) Report::report << " nx (=" <<nx<<") must be > order.";
      if(pbcY) Report::report << " ny (=" <<ny<<") must be > order.";
      if(pbcZ) Report::report << " nz (=" <<nz<<") must be > order.";
      Report::report<<Report::endr;
    }
    if((!pbcX && h.c[0] <= 0.0) || (!pbcX && h.c[0] <= 0.0) || (!pbcX && h.c[0] <= 0.0)){
      Report::report << Report::error << "[MultiGrid::initialize]";
      if(pbcX) Report::report << " h.x (=" <<h.c[0]<<") must be > 0.";
      if(pbcY) Report::report << " h.y (=" <<h.c[1]<<") must be > 0.";
      if(pbcZ) Report::report << " h.h (=" <<h.c[2]<<") must be > 0.";
      Report::report<<Report::endr;
    }    

    myScaledParticleIntPositions.resize(n);
    myInterpolations.resize(n,Interpolation3D(interOrder));

    myDim.resize(levels);
    myScale.resize(levels);
    myQ.resize(levels);
    myV.resize(levels);
    myBlocksCorrection.resize(0);

    // Tells us if the interpolation is 1 for one theta and 0 for the others (when w=0),
    // such that we can save one interpolation point for the next grid.
    // We have to keep track of that for the allocation of arrays since we
    // will need an extra point such that we can keep the interpolation loops
    // easy.
    // Hermitian interpolation has this nice property.
    mySigma = (TInterpolation::isSigma(interOrder)?1:0);

    myRatio = ratio;

    // The interpolation between the grids.
    myGridInterpolation.resize(myRatio,TInterpolation(interOrder));
    for(int i=0;i<myRatio;i++)
      myGridInterpolation[i].set((Real)i/(Real)myRatio);

    myS  = s;
    myRS = 1.0/s;
    myLevels = levels;
    myInterOrder = interOrder;
    myH = h;
    myOrigin = origin;

    // Toplevel
    myDim[levels-1].x = nx;
    myDim[levels-1].y = ny;
    myDim[levels-1].z = nz;


    for(int i=0;i<levels;i++){
      myScale[i] = power(ratio,i);
    }


    updateSize(min,max);


    Report::report << hint <<"MultiGrid: s="<<myS<<", ratio="<<myRatio<<", h="<<myHX
	   <<","<<myHY<<","<<myHZ<<", "<<TInterpolation::keyword
	   <<" "<<myInterOrder<<"th order "<<TKernel::keyword<<", levels=";
    for(int i=0;i<levels;i++)
      Report::report << "("<< myDim[i].x<<","<<myDim[i].y<<","<<myDim[i].z<<"),";
    Report::report<<" finest grid="<<myMin<<"-"<<myMax<<", Gc("<<myGCorrDimX<<","<<myGCorrDimY<<","<<myGCorrDimZ<<")."<<Report::endr;

  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::updateSize(Vector3D min, Vector3D max){

    //Report::report.reset();
    //Report::report.setf(std::ios::scientific);
    //Report::report.precision(4);
    Vector3D inMin(min),inMax(max);
    int nx=0;
    int ny=0;
    int nz=0;
    int nx1=0;
    int ny1=0;
    int nz1=0;
    int b = myInterOrder-mySigma-1;

    // Vacuum: Minimal number of grid points of the finest grid
    if(!pbcX){
      min.c[0] -= myH.c[0]*BORDER_TOLERANCE;
      max.c[0] += myH.c[0]*BORDER_TOLERANCE;    
      nx1 = (int)(ceil((max.c[0]-myOrigin.c[0])/myH.c[0])-floor((min.c[0]-myOrigin.c[0])/myH.c[0])+ b+1);
      nx = nx1+BORDER;
    }
    if(!pbcY){
      min.c[1] -= myH.c[1]*BORDER_TOLERANCE;
      max.c[1] += myH.c[1]*BORDER_TOLERANCE;    
      ny1 = (int)(ceil((max.c[1]-myOrigin.c[1])/myH.c[1])-floor((min.c[1]-myOrigin.c[1])/myH.c[1]) + b+1);
      ny = ny1+BORDER;
    }
    if(!pbcZ){
      min.c[2] -= myH.c[2]*BORDER_TOLERANCE;
      max.c[2] += myH.c[2]*BORDER_TOLERANCE;    
      nz1 = (int)(ceil((max.c[2]-myOrigin.c[2])/myH.c[2])-floor((min.c[2]-myOrigin.c[2])/myH.c[2])+ b+1);
      nz = nz1+BORDER;
    }
    // Vacuum: Compute number of grid points for the coarsest grid
    for(int i=1;i<myLevels;i++){
      if(!pbcX)
	nx = (int)ceil((Real)(nx-1)/(Real)(myRatio))+ b + 1;
      if(!pbcY)
	ny = (int)ceil((Real)(ny-1)/(Real)(myRatio))+ b + 1;
      if(!pbcZ)
	nz = (int)ceil((Real)(nz-1)/(Real)(myRatio))+ b + 1;
    }

    // Vacuum: No re-size if inside finest grid and
    //                       no major grid size changes
    if(!pbcX && !pbcY && !pbcZ){
      if(min.c[0] > myCMin.c[0] && min.c[1] > myCMin.c[1] && min.c[2] > myCMin.c[2] &&
	 max.c[0] < myCMax.c[0] && max.c[1] < myCMax.c[1] && max.c[2] < myCMax.c[2] &&
	 nx+BORDER>=myDim[myLevels-1].x && 
	 ny+BORDER>=myDim[myLevels-1].y && 
	 nz+BORDER>=myDim[myLevels-1].z)
	return;
    }

    if(pbcX)
      nx=myDim[myLevels-1].x;
    if(pbcY)
      ny=myDim[myLevels-1].y;
    if(pbcZ)
      nz=myDim[myLevels-1].z;

    bool update = false;

    if(pbcX || pbcY || pbcZ ||
       nx1 > myNX || ny1 > myNY || nz1 > myNZ || 
       fabs((float)(nx-myDim[myLevels-1].x))>BORDER ||
       fabs((float)(ny-myDim[myLevels-1].y))>BORDER ||
       fabs((float)(nz-myDim[myLevels-1].z))>BORDER){

      update = true;

      // Toplevel, coarsest grid
      if(!myQ[myLevels-1].resize(ArraySizes(nx+(pbcX?0:mySigma))(ny+(pbcY?0:mySigma))(nz+(pbcZ?0:mySigma))))
	Report::report << Report::error <<"[MultiGrid<>::initialize] Could not allocate memory for Q["
	       <<(nx+(pbcX?0:mySigma))<<"]["<<(ny+(pbcY?0:mySigma))<<"]["
	       <<(nz+(pbcZ?0:mySigma))<<"]."<<Report::endr;
      if(!myV[myLevels-1].resize(ArraySizes(nx+(pbcX?0:mySigma))(ny+(pbcY?0:mySigma))(nz+(pbcZ?0:mySigma))))
	Report::report << Report::error <<"[MultiGrid<>::initialize] Could not allocate memory for V["
	       <<(nx+(pbcX?0:mySigma))<<"]["<<(ny+(pbcY?0:mySigma))<<"]["
	       <<(nz+(pbcZ?0:mySigma))<<"]."<<Report::endr;
    
      myDim[myLevels-1].x = nx;
      myDim[myLevels-1].y = ny;
      myDim[myLevels-1].z = nz;

      // The other levels
      for(int i=myLevels-2;i>=0;i--){
	if(pbcX)
	  nx = nx * myRatio;
	else
	  nx = (nx-myInterOrder+mySigma) * myRatio + 1;
      
	if(pbcY)
	  ny = ny * myRatio;
	else
	  ny = (ny-myInterOrder+mySigma) * myRatio + 1;
      
	if(pbcZ)
	  nz = nz * myRatio;
	else
	  nz = (nz-myInterOrder+mySigma) * myRatio + 1;
      
	if(!myQ[i].resize(ArraySizes(nx+(pbcX?0:mySigma))(ny+(pbcY?0:mySigma))(nz+(pbcZ?0:mySigma))))
	  Report::report << Report::error <<"[MultiGrid<>::initialize] Could not allocate memory for Q["
		 <<(nx+(pbcX?0:mySigma))<<"]["<<(ny+(pbcY?0:mySigma))<<"]["
		 <<(nz+(pbcZ?0:mySigma))<<"]."<<Report::endr;
	if(!myV[i].resize(ArraySizes(nx+(pbcX?0:mySigma))(ny+(pbcY?0:mySigma))(nz+(pbcZ?0:mySigma))))
	  Report::report << Report::error <<"[MultiGrid<>::initialize] Could not allocate memory for V["
		 <<(nx+(pbcX?0:mySigma))<<"]["<<(ny+(pbcY?0:mySigma))<<"]["
		 <<(nz+(pbcZ?0:mySigma))<<"]."<<Report::endr;
	myDim[i].x = nx;
	myDim[i].y = ny;
	myDim[i].z = nz;
      }

      myNX = myDim[0].x; 
      myNY = myDim[0].y;
      myNZ = myDim[0].z;

      // Offset for mapping gridpoints between two levels
      myNXOffset = -(myInterOrder-1)/2;
      if(pbcX) myNXOffset += myNX;
      myNYOffset = -(myInterOrder-1)/2;
      if(pbcY) myNYOffset += myNY;
      myNZOffset = -(myInterOrder-1)/2;
      if(pbcZ) myNZOffset += myNZ;

      // Mesh size of the finest grid
      if(pbcX)
	myHX = (max.c[0] - min.c[0])/((Real)myNX);
      else
	myHX = myH.c[0];
      if(pbcY)
	myHY = (max.c[1] - min.c[1])/((Real)myNY);
      else
	myHY = myH.c[1];
      if(pbcZ)
	myHZ = (max.c[2] - min.c[2])/((Real)myNZ);    
      else
	myHZ = myH.c[2];    
    
      myHXr = 1.0/myHX;
      myHYr = 1.0/myHY;
      myHZr = 1.0/myHZ;
    }
  
    if(!pbcX){
      min.c[0] = floor((min.c[0]-myOrigin.c[0])/myH.c[0]-(myNX-nx1)/2)*myH.c[0]+myOrigin.c[0];    
      max.c[0] = min.c[0]+(myNX-1.0-b)*myH.c[0];
      myCMin.c[0] = min.c[0];
      myCMax.c[0] = max.c[0];
      min.c[0] -= myH.c[0]* (b/2);
      max.c[0] += myH.c[0]* (b-b/2);
    }
    if(!pbcY){
      min.c[1] = floor((min.c[1]-myOrigin.c[1])/myH.c[1]-(myNY-ny1)/2)*myH.c[1]+myOrigin.c[1];
      max.c[1] = min.c[1]+(myNY-1.0-b)*myH.c[1];
      myCMin.c[1] = min.c[1];
      myCMax.c[1] = max.c[1];
      min.c[1] -= myH.c[1]* (b/2);
      max.c[1] += myH.c[1]* (b-b/2);
    }
    if(!pbcZ){
      min.c[2] = floor((min.c[2]-myOrigin.c[2])/myH.c[2]-(myNZ-nz1)/2)*myH.c[2]+myOrigin.c[2];    
      max.c[2] = min.c[2]+(myNZ-1.0-b)*myH.c[2];
      myCMin.c[2] = min.c[2];
      myCMax.c[2] = max.c[2];
      min.c[2] -= myH.c[2]* (b/2);
      max.c[2] += myH.c[2]* (b-b/2);
    }
  
    if(!update && ((myMin-min).normSquared() > Constant::EPSILON || 
		   (myMax-max).normSquared() > Constant::EPSILON)){
      Report::report.reset();
      Report::report.precision(7);
      Report::report << Report::hint << "MultiGrid:"<<myMin<<"-"<<myMax<<" to "<<min<<"-"<<max<<Report::endr;  
    }

    myMin = min;
    myMax = max;
  
    if(update){
      myBlocksCorrection.resize(0);
      if(!myGDirect.empty()){
	Report::report.reset();
	Report::report.precision(7);
	Report::report << Report::hint << "MultiGrid: Finest grid("<<myNX<<","<<myNY<<","<<myNZ<<") was re-sized to "<<myMin<<"-"<<myMax<<" ("<<inMin<<"-"<<inMax<<")."<<Report::endr;    
      }
      precomputeG();    
    }
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  inline void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::clear(){
    for(int level=0;level<myLevels;level++){
      //     Array<Real,3>& rV = myV[level];
      //     Array<Real,3>& rQ = myQ[level];
      //     for (int i = 0; i < myDim[level].c[0]+(pbcX?0:mySigma); i++){
      //       for (int j = 0 ; j < myDim[level].c[1]+(pbcY?0:mySigma); j++){
      // 	for (int k = 0; k < myDim[level].c[2]+(pbcZ?0:mySigma); k++){
      // 	  rV[i][j][k] = 0.0;
      // 	  rQ[i][j][k] = 0.0;
      // 	}
      //       }
      //     }
      int n = myQ[level].size();
      Real *q = myQ[level].begin();
      Real *v = myV[level].begin();
      for(int i=0;i<n;i++){
	q[i] = 0.0;
	v[i] = 0.0;
      }
    }
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::printQ(int level){
    Report::report << plain << "Level "<<level<<":"<<Report::endr;
    Real q = 0.0;

    for (int i = 0; i < myDim[level].x; i++){
      for (int j = 0 ; j < myDim[level].y; j++){
	Report::report << plain << "Q["<<i<<"]["<<j<<"][0-"<<myDim[level].z-1<<"] : ";
	for (int k = 0; k < myDim[level].z; k++){
	  Report::report <<myQ[level][i][j][k]<<" ";
	  q += myQ[level][i][j][k];
	}
	Report::report <<Report::endr;
      }
    }
    Report::report << plain <<"Sum Q :"<<q<<Report::endr;
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::printV(int level){
    Report::report << plain << "Level "<<level<<":"<<Report::endr;
    Real v = 0.0;
    for (int i = 0; i < myDim[level].x; i++){
      for (int j = 0 ; j < myDim[level].y; j++){
	Report::report << plain << "V["<<i<<"]["<<j<<"][0-"<<myDim[level].z-1<<"] : ";
	for (int k = 0; k < myDim[level].z; k++){
	  Report::report <<myV[level][i][j][k]<<" ";
	  v += myV[level][i][j][k];
	}
	Report::report <<Report::endr;
      }
    }
    Report::report << plain <<"Sum V :"<<v<<Report::endr;
  }

  template<class TInterpolation, class TKernel, bool pbcX, bool pbcY, bool pbcZ>
  void MultiGrid<TInterpolation,TKernel,pbcX,pbcY,pbcZ>::printConst(){
    for (int i = 0; i < myDim[myLevels-1].x; i++){
      for (int j = 0 ; j < myDim[myLevels-1].y; j++){
	Report::report << plain << "Gd["<<i<<"]["<<j<<"][0-"<<myDim[myLevels-1].z-1<<"] : ";
	for (int k = 0; k < myDim[myLevels-1].z; k++){
	  Report::report <<myGDirect[i][j][k]<<" ";
	}
	Report::report <<Report::endr;
      }
    }

    // No correction if just one level specified.
    if(myLevels < 2 )
      return;
    for (int i = 0; i < myGCorrDimX; i++){
      for (int j = 0 ; j < myGCorrDimY; j++){
	Report::report << plain << "Gc["<<i<<"]["<<j<<"][0-"<<myGCorrDimZ-1<<"] : ";
	for (int k = 0; k < myGCorrDimZ; k++){
	  Report::report <<myGCorrection[i][j][k]<<" ";
	}
	Report::report <<Report::endr;
      }
    }
  }
}
#endif
