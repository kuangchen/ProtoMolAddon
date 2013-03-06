/* -*- c++ -*- */
#ifndef NONBONDEDFULLEWALDSYSTEMFORCE_H
#define NONBONDEDFULLEWALDSYSTEMFORCE_H

#include <protomol/force/system/SystemForce.h>
#include <protomol/parallel/Parallel.h>
#include <protomol/type/SimpleTypes.h>

//using namespace ProtoMol::Report;

//#define DEBUG_EWALD_TIMING
//#define DEBUG_EWALD_ENERGIES
//#define USE_EWALD_EXACT_ERF
//#define USE_EWALD_NO_SINCOS_TABLE

namespace ProtoMol {

  //_________________________________________________________________ NonbondedFullEwaldSystemForce
  
  template<class TBoundaryConditions, 
	   class TCellManager,
	   bool  real,
	   bool  reciprocal,
	   bool  correction,
	   class TSwitchingFunction>
  class NonbondedFullEwaldSystemForce: public SystemForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    typedef Topology<TBoundaryConditions, TCellManager> RealTopologyType;
    typedef typename RealTopologyType::Enumerator EnumeratorType;
    typedef CellPair CellPairType;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NonbondedFullEwaldSystemForce();
    NonbondedFullEwaldSystemForce(Real alpha, Real accuracy, Real expansionFactor);

    virtual ~NonbondedFullEwaldSystemForce();
  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NonbondedFullEwaldSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

  private:
    virtual Force* doMake(const std::vector<Value> &values) const;
    void initialize(const RealTopologyType* realTopo, const Vector3DBlock* positions);

    void realTerm(const RealTopologyType* realTopo,
		  const Vector3DBlock* positions, 
		  Vector3DBlock* forces, 
		  ScalarStructure* energies,
		  Real& realEnergy,
		  unsigned int n);

    void reciprocalTerm(const RealTopologyType* realTopo,
			const Vector3DBlock* positions, 
			Vector3DBlock* forces, 
			ScalarStructure* energies,
			Real& reciprocalEnergy,
			unsigned int from, unsigned int to);

    void correctionTerm(const RealTopologyType* realTopo,
			const Vector3DBlock* positions, 
			Vector3DBlock* forces, 
			ScalarStructure* energies,
			Real& intraMolecularEnergy,
			unsigned int from, unsigned int to);

    void surfaceDipoleTerm(const RealTopologyType* realTopo,
			   const Vector3DBlock* positions, 
			   Vector3DBlock* forces, 
			   ScalarStructure* energies,
			   Real& surfaceDipoleEnergy);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual void evaluate(const GenericTopology*, 
			  const Vector3DBlock*, 
			  Vector3DBlock*, 
			  ScalarStructure*);

    virtual void parallelEvaluate(const GenericTopology*, 
				  const Vector3DBlock*, 
				  Vector3DBlock*, 
				  ScalarStructure*);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual unsigned int numberOfBlocks(const GenericTopology*,const Vector3DBlock*);
    virtual std::string getKeyword() const{return "NonbondedFullEwald";}
    virtual void uncache(){myCached=false;};
  private:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual std::string getIdNoAlias() const;
    virtual void getParameters(std::vector<Parameter>& parameters) const;
    virtual unsigned int getParameterSize() const{return (TBoundaryConditions::PERIODIC ? 2:3);}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  private:
    bool myCached;

    std::vector<Vector3D> myK;  // Reciprocal lattice vector (h*2PI/Lx,k*2PI/Ly,l*2PI/Lz)
    std::vector<Real> myKSquared;  // Squared norm of the reciprocal lattice vector
#ifndef USE_EWALD_NO_SINCOS_TABLE
    std::vector<TripleInt> myKInt;      // Reciprocal lattice vector (h,k,l)
#endif

    Real myExpansionFactor;
    Real myTRate;             //
    Real myAccuracy;

    Real myAlpha;             // 
    bool myAlphaDefault;

    Real myAlphaSquared;
    Real myAlphaSquaredr;
    Real my2AlphaPI;
    Real myFac;

    Real myRc;                // Cutoff real term
    Real myRcSquared;         // Cutoff squared real term
    Real myKc;                // Cutoff reciprocal term

    int myHmax;
    int myKmax;
    int myLmax;
    int myHKLmax;

    Real* mySinCosA;          // Look up tables
    Real* myLastSinCos;
    Vector3D* mySinCos;

    Real myLX, myLY, myLZ;
    Real myLXr,myLYr,myLZr;
    Real myV,myVr;
    Vector3D myOrigin;

    Real myPointSelfEnergy;   // Precomputed energy terms
    Real myChargedSystemEnergy;

    PeriodicBoundaryConditions boundaryConditions;

    TSwitchingFunction switchingFunction;
    EnumeratorType enumerator;
    std::vector<Vector3D> myLattice;
    unsigned int myOldAtomCount; // Keep track of old values and hope
    int myOldHKLmax;             // we do not need to reallocate memory ...
#if defined(DEBUG_EWALD_TIMING)
    Timer myReal;
    Timer myReciprocal;
    Timer myIntra;
    Timer mySurface;
#endif
  };

  //______________________________________________________________________ INLINES
  template <class TBoundaryConditions, 
	    class TCellManager, 
	    bool  real, 
	    bool  reciprocal, 
	    bool  correction, 
	    class TSwitchingFunction>
  NonbondedFullEwaldSystemForce<TBoundaryConditions,
				TCellManager,
				real,
				reciprocal,
				correction,
				TSwitchingFunction>::NonbondedFullEwaldSystemForce() 
				  : SystemForce(),
				    myCached(false),
				    myExpansionFactor(0.0),
				    myTRate(0.0),
				    myAccuracy(0.0),
				    myAlpha(-1.0),
				    myAlphaDefault(true),
				    myRc(0.0),
				    myKc(0.0),
				    mySinCosA(NULL),
				    myLastSinCos(NULL),
				    mySinCos(NULL),
				    myV(-1.0),
				    myOldAtomCount(0),
				    myOldHKLmax(0){
#if defined(DEBUG_EWALD_TIMING)
    myReal.reset();
    myReciprocal.reset();
    myIntra.reset();
    mySurface.reset();
#endif
  }
  template <class TBoundaryConditions, 
	    class TCellManager, 
	    bool  real, 
	    bool  reciprocal, 
	    bool  correction, 
	    class TSwitchingFunction>
  NonbondedFullEwaldSystemForce<TBoundaryConditions,
				TCellManager,
				real,reciprocal,
				correction,
				TSwitchingFunction>::NonbondedFullEwaldSystemForce(Real alpha, Real accuracy, Real expansionFactor)
				  : SystemForce(),
				    myCached(false),
				    myExpansionFactor(expansionFactor),
				    myTRate(0.0),
				    myAccuracy(accuracy),
				    myAlpha(alpha),
				    myAlphaDefault(alpha<= 0.0),
				    myRc(0.0),
				    myKc(0.0),
				    mySinCosA(NULL),
				    myLastSinCos(NULL),
				    mySinCos(NULL),
				    myV(-1.0),
				    myOldAtomCount(0),
				    myOldHKLmax(0) {
#if defined(DEBUG_EWALD_TIMING)
    myReal.reset();
    myReciprocal.reset();
    myIntra.reset();
    mySurface.reset();
#endif
  }

  template <class TBoundaryConditions, 
	    class TCellManager, 
	    bool  real, 
	    bool  reciprocal, 
	    bool  correction, 
	    class TSwitchingFunction>
  NonbondedFullEwaldSystemForce<TBoundaryConditions,
				TCellManager,
				real,
				reciprocal,
				correction,
				TSwitchingFunction>::~NonbondedFullEwaldSystemForce(){
    if(reciprocal){
      delete [] mySinCosA;
      delete [] mySinCos;
      delete [] myLastSinCos;
    }
#if defined(DEBUG_EWALD_TIMING)
    if(boundaryConditions.getVolume() > Constnat::EPSILON){
      report.setf(std::ios::showpoint|std::ios::fixed);
      report << allnodes << plain <<"Timing ("
	     <<Parallel::getId()<<") Ewald:"
	     <<": real:" <<myReal.getTime().getProcessTime()<<"[s]"
	     <<", reciprocal:"<<myReciprocal.getTime().getProcessTime()<<"[s]"
	     <<", intra:"<<myIntra.getTime().getProcessTime()<<"[s]"
	     <<", dipole:"<<mySurface.getTime().getProcessTime()<<"[s]"
	     <<", Rc="<<myRc
	     <<", Kc="<<myKc
	     <<", Kn="<<myK.size()
	     <<", n="<<myLattice.size()
	     <<", alpha="<<myAlpha
	     <<", accuracy="<<myAccuracy<<"."<<endr;
    }
#endif
  }
  
  template <class TBoundaryConditions, 
	    class TCellManager, 
	    bool  real, 
	    bool  reciprocal, 
	    bool  correction, 
	    class TSwitchingFunction>
  void NonbondedFullEwaldSystemForce<TBoundaryConditions,
				     TCellManager,
				     real,
				     reciprocal,
				     correction,
				     TSwitchingFunction>::evaluate(const GenericTopology* topo, 
								   const Vector3DBlock* positions, 
								   Vector3DBlock* forces, 
								   ScalarStructure* energies) {
    
    const RealTopologyType* realTopo = (RealTopologyType*) topo;  
    // Initialize data members and precompute tables & short cuts
    if(!myCached)
      initialize(realTopo,positions);
    
    // Intra-molecular and surface diplol term
    Real intraMolecularEnergy = 0.0;
    Real surfaceDipoleEnergy = 0.0;
    Real pointSelfEnergy = 0.0;
    Real chargedSystemEnergy = 0.0;
    if(correction){
      correctionTerm(realTopo,positions,forces,energies,
		     intraMolecularEnergy,
		     0,realTopo->exclusions.getTable().size());
      pointSelfEnergy     = myPointSelfEnergy;
      chargedSystemEnergy = myChargedSystemEnergy;
      if(energies->virial()){
	(*energies)[ScalarStructure::VIRIALXX] += myChargedSystemEnergy;
	(*energies)[ScalarStructure::VIRIALYY] += myChargedSystemEnergy;
	(*energies)[ScalarStructure::VIRIALZZ] += myChargedSystemEnergy;
      }
      if(false)
	surfaceDipoleTerm(realTopo,positions,forces,energies,surfaceDipoleEnergy);
    }

    // Real-space term
    Real realEnergy = 0.0;
    if(real){
      realTopo->updateCellLists(positions);
      enumerator.initialize(realTopo, myRc);

      realTerm(realTopo,positions,forces,energies,realEnergy,realTopo->cellLists.size());
    }  

    // Reciprocal-space term
    Real reciprocalEnergy = 0.0;
    if(reciprocal)
      reciprocalTerm(realTopo,positions,forces,energies,reciprocalEnergy,
		     0,myK.size());

    // Sum of all energy terms
    // Sum of all energy terms
    Real e = 
      realEnergy+
      reciprocalEnergy+
      intraMolecularEnergy+
      surfaceDipoleEnergy+
      pointSelfEnergy+
      chargedSystemEnergy;

    (*energies)[ScalarStructure::COULOMB] += e;

#if defined(DEBUG_EWALD_ENERGIES)
    report.setf(std::ios::showpoint|std::ios::fixed);
    report << plain <<"Ewald: point="<<myPointSelfEnergy
	   <<", charged="<<myChargedSystemEnergy
	   <<", real="<<realEnergy
	   <<", reciprocal="<<reciprocalEnergy
	   <<", intra="<<intraMolecularEnergy
	   <<", surface="<<surfaceDipoleEnergy
	   <<", total="<<e 
	   <<endr;
#endif
  }

  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void NonbondedFullEwaldSystemForce<TBoundaryConditions,
				     TCellManager,
				     real,
				     reciprocal,
				     correction,
				     TSwitchingFunction>::parallelEvaluate(const GenericTopology* topo, 
									   const Vector3DBlock* positions, 
									   Vector3DBlock* forces, 
									   ScalarStructure* energies) {

    const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);  

    // Initialize data members and precompute tables & short cuts
    if(!myCached)
      initialize(realTopo,positions);
    
    // Intra-molecular and surface diplol term
    Real intraMolecularEnergy = 0.0;
    Real surfaceDipoleEnergy = 0.0;
    Real pointSelfEnergy = 0.0;
    Real chargedSystemEnergy = 0.0;
    if(correction){
      unsigned int n = realTopo->exclusions.getTable().size();
      if(n > 0){
	unsigned int count = std::min(n,static_cast<unsigned int>(Parallel::getAvailableNum()));    
	for(unsigned int i = 0;i<count;i++){
	  if(Parallel::next()){
	    unsigned int to = (n*(i+1))/count;
	    if(to > n)
	      to = n;
	    unsigned int from = (n*i)/count;
	    correctionTerm(realTopo,positions,forces,energies,
			   intraMolecularEnergy,from,to);
	  }
	}
      }
      if(Parallel::getAvailableId() == 0){
	pointSelfEnergy     = myPointSelfEnergy;
	chargedSystemEnergy = myChargedSystemEnergy;
	if(energies->virial()){
	  (*energies)[ScalarStructure::VIRIALXX] += myChargedSystemEnergy;
	  (*energies)[ScalarStructure::VIRIALYY] += myChargedSystemEnergy;
	  (*energies)[ScalarStructure::VIRIALZZ] += myChargedSystemEnergy;
	}
	if(false)
	  surfaceDipoleTerm(realTopo,positions,forces,energies,surfaceDipoleEnergy);
      }
    }

    // Real-space term
    Real realEnergy = 0.0;
    if(real){
      realTopo->updateCellLists(positions);
      enumerator.initialize(realTopo, myRc);
      unsigned int n = realTopo->cellLists.size();
      unsigned int count = Parallel::getNumberOfPackages(n);

      for(unsigned int i = 0;i<count;i++){
	unsigned int l = (n*(i+1))/count - (n*i)/count;
	if(Parallel::next())
	  realTerm(realTopo,positions,forces,energies,realEnergy,l);
	else
	  enumerator.nextNewPair(l);	
      }
    }  

    // Reciprocal-space term
    Real reciprocalEnergy = 0.0;
    if(reciprocal){   
      unsigned int count = std::min(static_cast<unsigned int>(myK.size()),static_cast<unsigned int>(Parallel::getAvailableNum()));

      for(unsigned int i = 0;i<count;i++)
	if(Parallel::next())	
	  reciprocalTerm(realTopo,positions,forces,energies,reciprocalEnergy,
			 (myK.size()*i)/count,(myK.size()*(i+1))/count);      

    }
    // Sum of all energy terms
    Real e = 
      realEnergy+
      reciprocalEnergy+
      intraMolecularEnergy+
      surfaceDipoleEnergy+
      pointSelfEnergy+
      chargedSystemEnergy;

    (*energies)[ScalarStructure::COULOMB] += e;

#if defined(DEBUG_EWALD_ENERGIES)
    report.setf(std::ios::showpoint|std::ios::fixed);
    report << allnodes << plain <<"Ewald: point="<<myPointSelfEnergy
	   <<", charged="<<myChargedSystemEnergy
	   <<", real="<<realEnergy
	   <<", reciprocal="<<reciprocalEnergy
	   <<", intra="<<intraMolecularEnergy
	   <<", surface="<<surfaceDipoleEnergy
	   <<", total="<<e 
	   <<endr;
#endif
  }

  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  unsigned int NonbondedFullEwaldSystemForce<TBoundaryConditions,
					     TCellManager,
					     real,
					     reciprocal,
					     correction,
					     TSwitchingFunction>::numberOfBlocks(const GenericTopology* topo,
										 const Vector3DBlock* positions){
    unsigned int n = 0;
    if(correction)
      n +=std::min(static_cast<int>(topo->exclusions.getTable().size()),static_cast<int>(Parallel::getAvailableNum()));

    if(reciprocal){
      if(!myCached)
	initialize(dynamic_cast<const RealTopologyType*>(topo),positions);

      n += std::min(static_cast<unsigned int>(myK.size()),static_cast<unsigned int>(Parallel::getAvailableNum()));
    }

    if(real){
      const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);
      realTopo->updateCellLists(positions);    
      n += Parallel::getNumberOfPackages(realTopo->cellLists.size());
    }

    return n;
  }


  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void NonbondedFullEwaldSystemForce<TBoundaryConditions,
				     TCellManager,
				     real,
				     reciprocal,
				     correction,
				     TSwitchingFunction>::initialize(const RealTopologyType* realTopo, const Vector3DBlock* positions) {
    bool dontHint = (myV >= 0.0 && TBoundaryConditions::PERIODIC);
    if(dontHint)
      Report::report << Report::donthint;

    if(TBoundaryConditions::VACUUM){
      // We are have a non-periodic case. Do an expansion of the simulation 
      // box by the expansion factor, if need
      realTopo->updateCellLists(positions);
      Vector3D d2 = (boundaryConditions.getMax()-boundaryConditions.getMin());
      Vector3D d0 = d2/myExpansionFactor;
      Vector3D d1 = realTopo->max-realTopo->min;
      if (fabs(1-d0.c[0]/d1.c[0]) > 0.1 || fabs(1-d0.c[1]/d1.c[1]) > 0.1 || fabs(1-d0.c[1]/d1.c[1]) > 0.1 || 
	  d2.c[0] <= d1.c[0] || d2.c[1] <= d1.c[1] || d2.c[2] <= d1.c[2] ){
	// The boundaries are to big or small, adjust
	Vector3D e = (realTopo->max-realTopo->min)*myExpansionFactor;
	Vector3D origin = realTopo->min*myExpansionFactor + e*0.5;
	boundaryConditions.set(Vector3D(e.c[0],0,0),
			       Vector3D(0,e.c[1],0),
			       Vector3D(0,0,e.c[2]),
			       origin);
	report << hint << "The boundaries for Ewald force evaluation were re-sized to "<<boundaryConditions.getMin()<<"-"<<boundaryConditions.getMax()<<"."<<endr;
      }
      else {
	// The particles are just inside +-10% our previous simulation box. No update needed
	myV = boundaryConditions.getVolume();
	return;      
      }
    }
    else {
      boundaryConditions.set(realTopo->boundaryConditions.e1(),
			     realTopo->boundaryConditions.e2(),
			     realTopo->boundaryConditions.e3(),
			     realTopo->boundaryConditions.origin());
    }

    if(!boundaryConditions.isOrthogonal())
      report << error << "[NonbondedFullEwaldSystemForce::initialize] Not orthogonal, aborting."<<endr;

    const unsigned int atomCount = realTopo->atoms.size();

      
    myTRate     = 5.5;     // From Moldy, rate between real and reciprocal
    //myAccuracy  = 0.00001; // From Moldy, accuracy
    //myAccuracy  = 1.e-6; // From NAMD, accuracy
  
    // Dimension of the simulation box
    myLX   = boundaryConditions.e1().c[0];
    myLY   = boundaryConditions.e2().c[1];
    myLZ   = boundaryConditions.e3().c[2];
    myV    = boundaryConditions.getVolume();
    myLXr  = boundaryConditions.e1r().c[0];
    myLYr  = boundaryConditions.e2r().c[1];
    myLZr  = boundaryConditions.e3r().c[2];
    myVr   = 1.0/myV;
    myOrigin = boundaryConditions.origin();

    // Short cuts
    if(myAlphaDefault)
      myAlpha         = sqrt(M_PI)*pow(myTRate*atomCount/(power<2>(myV)),1.0/6.0);
    myAlphaSquared  = myAlpha*myAlpha;
    myAlphaSquaredr = 1.0/myAlphaSquared;
    my2AlphaPI      = 2.0*myAlpha/sqrt(M_PI);
    Real p          = -log(myAccuracy);
    myRc            = sqrt(p)/myAlpha;
    myRcSquared     = myRc*myRc;
    myKc            = 2.0*myAlpha*sqrt(p);
    myFac           =  1.0/(4.0*myAlpha*myAlpha);

    switchingFunction=TSwitchingFunction(myRc);

    // Reciprocal part
    // Build the lattice k-vectors 2*PI(n0/Lx,n1/Ly,n2/Lz)
    // Maximum values of h, k, l  s.t. |k| < myKc
    myHmax = (int)floor(myKc/(2.0*M_PI)*myLX);
    myKmax = (int)floor(myKc/(2.0*M_PI)*myLY);
    myLmax = (int)floor(myKc/(2.0*M_PI)*myLZ);
    myHKLmax = std::max(2,std::max(myHmax,std::max(myKmax,myLmax))+1);
    myK.clear();
    myKSquared.clear();
#ifndef USE_EWALD_NO_SINCOS_TABLE
    myKInt.clear();
#endif
    int lastH = Constant::MAX_INT;
    int lastK = Constant::MAX_INT;
    int misses = 0;
    Real kcSquared = myKc*myKc;
    if(reciprocal){
      for(int h = 0; h <= myHmax; h++){
	Real kx = 2.0*M_PI*h*myLXr;
	for(int k = (h==0 ? 0 : -myKmax); k <= myKmax; k++){
	  Real ky = 2.0*M_PI*k*myLYr;
	  for(int l = (h==0 && k==0 ? 1 : -myLmax); l <= myLmax; l++){
	    Real kz = 2.0*M_PI*l*myLZr;
	    if(kx*kx + ky*ky + kz*kz < kcSquared){
	      myK.push_back(Vector3D(kx,ky,kz));
	      myKSquared.push_back(kx*kx + ky*ky + kz*kz);
#ifndef USE_EWALD_NO_SINCOS_TABLE
	      TripleInt tmp(h,k,l);
	      myKInt.push_back(tmp);
	      if(lastH != h || lastK != k)
		misses++;
	      lastH = h;
	      lastK = k;
#endif
	    }
	  }
	}
      }
    
      // Allocate memeory for our tables
      if(mySinCosA != NULL && (atomCount != myOldAtomCount)){
	delete [] mySinCosA;
	mySinCosA = NULL;
      }
      if(mySinCos != NULL && (atomCount != myOldAtomCount || myHKLmax != myOldHKLmax)){
	delete [] mySinCos;
	mySinCos = NULL;
      }
      if(myLastSinCos != NULL && (atomCount != myOldAtomCount)){
	delete [] myLastSinCos;
	myLastSinCos = NULL;
      }
      if(mySinCosA == NULL)
	mySinCosA    = new Real[2*atomCount];
      if(mySinCosA == NULL)
	report << error << "[NonbondedFullEwaldSystemForce::evaluate] Not enough memory, requesting "<<2*atomCount*sizeof(Real)<<" bytes."<<endr;
#ifndef USE_EWALD_NO_SINCOS_TABLE
      if(mySinCos == NULL)
	mySinCos     = new Vector3D[2*atomCount*myHKLmax];
      if(myLastSinCos == NULL)
	myLastSinCos = new Real[2*atomCount];
      if(mySinCos == NULL || myLastSinCos == NULL)
	report << error << "[NonbondedFullEwaldSystemForce::evaluate] Not enough memory, requesting "<<2*atomCount*myHKLmax *sizeof(Real)*3 + 2*atomCount *sizeof(Real)<<" bytes."<<endr;
      myOldHKLmax = myHKLmax;
#endif
      myOldAtomCount = atomCount;
    }
  
    //
    // Point self-energy
    //
    myPointSelfEnergy = 0.0;
    if(correction){
      Real q = 0.0;
      for(unsigned int i=0;i<atomCount;i++){
	q += realTopo->atoms[i].scaledCharge*realTopo->atoms[i].scaledCharge;
      }
      myPointSelfEnergy = -q*myAlpha/sqrt(M_PI);
    }
  
    //
    // Charged system energy
    //
    myChargedSystemEnergy = 0.0;
    if(correction){
      Real q = 0.0;
      for(unsigned int i=0;i<atomCount;i++){
	q += realTopo->atoms[i].scaledCharge;
      }
      if(fabs(q * 0.00268283) > 1.0e-5)
	myChargedSystemEnergy = -M_PI/(2.0*myV*myAlphaSquared)*q*q;
    }
  
    myLattice = boundaryConditions.buildLatticeVectors(myRc);
    myLattice.insert(myLattice.begin(),Vector3D(0,0,0));

    report << hint <<"Ewald";
#if defined(USE_EWALD_EXACT_ERF) && defined(USE_EWALD_NO_SINCOS_TABLE)
    report << "(Exact)";
#else
#ifdef USE_EWALD_EXACT_ERF
    report << "(Exact erf)";
#endif
#ifdef USE_EWALD_NO_SINCOS_TABLE
    report << "(No sincos table)";
#endif
#endif  
    report <<": alpha="<<toString(myAlpha)<<", V="<<myV<<", Rc="
	   <<toString(myRc)<<", Kc ("<<myK.size()<<")="<<toString(myKc)<<", n="<<myLattice.size()<<", accuracy="
	   <<myAccuracy<<", misses="
	   <<(Real)100*misses/(Real)(myK.empty() ? 1:myK.size()) <<"%."<<endr;

    // Now we have all pre-computed stuff, non-periodic case may need an update when 
    // the simulation box expands or shrinks to much (+-10%)
    if(TBoundaryConditions::PERIODIC)
      myCached = true;

    myV = boundaryConditions.getVolume();

    if(dontHint)
      Report::report << Report::dohint;
  }


  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void NonbondedFullEwaldSystemForce<TBoundaryConditions,
				     TCellManager,
				     real,
				     reciprocal,
				     correction,
				     TSwitchingFunction>::realTerm(const RealTopologyType* realTopo,
								   const Vector3DBlock* positions, 
								   Vector3DBlock* forces, 
								   ScalarStructure* energies,
								   Real& realEnergy,
								   unsigned int n) {
#if defined(DEBUG_EWALD_TIMING)
    myReal.start();
#endif

    // Real-space term
    CellPairType thisPair;
    bool doVirial = energies->virial();
    bool doMolVirial = energies->molecularVirial();
    unsigned int count = 0;
    for (; !enumerator.done(); enumerator.next()) {
      enumerator.get(thisPair);
      bool notSameCell = enumerator.notSameCell();
    
      if(!notSameCell){
	count++;
	if(count > n)
	  break;
      }

      for(int i=thisPair.first; i!=-1; i=realTopo->atoms[i].cellListNext){
	Real qi  = realTopo->atoms[i].scaledCharge;
	Vector3D ri((*positions)[i]),fi;
	int mi = realTopo->atoms[i].molecule;
	for(int j=(notSameCell ? thisPair.second:i); j!=-1; j=realTopo->atoms[j].cellListNext){
	  Vector3D rijMinimal(boundaryConditions.minimalDifference(ri,(*positions)[j]));
	  int mj = realTopo->atoms[j].molecule;
	  Real qj = realTopo->atoms[j].scaledCharge;
	  bool same = (mi==mj);
	  ExclusionClass excl = (same?realTopo->exclusions.check(i,j):EXCLUSION_NONE);
	  if(i == j)
	    excl = EXCLUSION_FULL;
	  for(unsigned int k=0;k<myLattice.size();k++){
	    // Check for an exclusion.
	    if (excl == EXCLUSION_FULL){ 
	      excl = EXCLUSION_NONE;
	      continue;
	    }
	    Vector3D rij(rijMinimal+myLattice[k]);
	    Real rSquared = rij.normSquared();
        
	    // Do switching function rough test.
	    if (rSquared>myRcSquared)
	      continue;
        
        
	    // Energy
	    Real r = sqrt(rSquared);
	    Real qq = qi*qj;
	    if (excl == EXCLUSION_MODIFIED) 
	      qq *= realTopo->coulombScalingFactor;
	    // Approximation Abramowitz & Stegun p299.
	    // Energy
#ifndef USE_EWALD_EXACT_ERF
	    Real rr = 1.0/r;
	    Real ar = myAlpha*r;
	    Real e = qq*exp(-ar*ar);
	    Real energy = poly5(ar)*e*rr;
	    Real force = ((energy+my2AlphaPI*e)*rr*rr);
#else	  
	    Real a = erfc(myAlpha*r)/r;
	    Real energy = qq*a;
	    Real force = qq*(a+my2AlphaPI*exp(-myAlphaSquared*rSquared))/rSquared;
#endif
	    // Calculate the switched force and energy.
	    Real switchingValue, switchingDeriv;
	    switchingFunction(switchingValue, switchingDeriv, rSquared);
	    // This has a - sign because the force is the negative of the 
	    // derivative of the energy (divided by the distance between the atoms).
	    force = force * switchingValue - energy * switchingDeriv;
	    energy = energy * switchingValue;
	    // Force, F_ij
	    realEnergy += energy; 
	    Vector3D fij(rij*force);
	    fi -= fij;
	    (*forces)[j] += fij;

	    // compute the vector between molecular centers of mass
	    if(!same && doMolVirial){
	      // Add to the atomic and molecular virials
	      energies->addVirial(fij,rij,realTopo->boundaryConditions.minimalDifference(realTopo->molecules[mi].position,
											  realTopo->molecules[mj].position));
	    }
	    else if(doVirial) {
	      energies->addVirial(fij,rij);
	    }
	    excl = EXCLUSION_NONE;
	  }
	}
	(*forces)[i] += fi;
      }
    }
#if defined(DEBUG_EWALD_TIMING)
    myReal.stop();
#endif
  }



  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void NonbondedFullEwaldSystemForce<TBoundaryConditions,
				     TCellManager,
				     real,
				     reciprocal,
				     correction,
				     TSwitchingFunction>::reciprocalTerm(const RealTopologyType* realTopo,
									 const Vector3DBlock* positions, 
									 Vector3DBlock* forces, 
									 ScalarStructure* energies,
									 Real& reciprocalEnergy,
									 unsigned int from, unsigned int to) {

#if defined(DEBUG_EWALD_TIMING)
    myReciprocal.start();
#endif
    const unsigned int atomCount = realTopo->atoms.size();
    

    Real energy = 0.0;
#ifndef USE_EWALD_NO_SINCOS_TABLE
    // Precompute/ cache cos/ sin (r*N*2*PI/L) for the lattice vectors in each dimension 
    for(unsigned int j=0;j<atomCount;j++){
      int l = 2*j*myHKLmax;
      Vector3D r(boundaryConditions.minimalPosition((*positions)[j]));
      // Multiply charge only with x-coord of each particle
      // since we use the add theorem 
      Real qi = realTopo->atoms[j].scaledCharge;
      Real x = r.c[0]*2.0*M_PI*myLXr;
      Real y = r.c[1]*2.0*M_PI*myLYr;
      Real z = r.c[2]*2.0*M_PI*myLZr;
      Real xsin = sin(x);
      Real ysin = sin(y);
      Real zsin = sin(z);
      Real xcos = cos(x);
      Real ycos = cos(y);
      Real zcos = cos(z);
      // The first two cos/ sin values
      // sin(r*0*2*PI/L)
      mySinCos[l  ].c[0] = 0.0;
      mySinCos[l  ].c[1] = 0.0;
      mySinCos[l  ].c[2] = 0.0;
      // cos(r*0*2*PI/L)
      mySinCos[l+1].c[0] = qi*1.0;
      mySinCos[l+1].c[1] = 1.0;
      mySinCos[l+1].c[2] = 1.0;
      // sin(r*1*2*PI/L)
      mySinCos[l+2].c[0] = qi*xsin;
      mySinCos[l+2].c[1] = ysin;
      mySinCos[l+2].c[2] = zsin;
      // cos(r*1*2*PI/L)
      mySinCos[l+3].c[0] = qi*xcos;
      mySinCos[l+3].c[1] = ycos;
      mySinCos[l+3].c[2] = zcos;
      
      // Using add theorem to compute sin(r*2*2*PI/L to r*myHKLmax*2*PI/L) 
      // and cos(r*2*2*PI/L to r*myHKLmax*2*PI/L)
      for(int i=4;i<2*myHKLmax;i+=2){
	mySinCos[l+i  ].c[0] = xsin*mySinCos[l+i-1].c[0]+xcos*mySinCos[l+i-2].c[0];
	mySinCos[l+i  ].c[1] = ysin*mySinCos[l+i-1].c[1]+ycos*mySinCos[l+i-2].c[1];
	mySinCos[l+i  ].c[2] = zsin*mySinCos[l+i-1].c[2]+zcos*mySinCos[l+i-2].c[2];
	mySinCos[l+i+1].c[0] = xcos*mySinCos[l+i-1].c[0]-xsin*mySinCos[l+i-2].c[0];
	mySinCos[l+i+1].c[1] = ycos*mySinCos[l+i-1].c[1]-ysin*mySinCos[l+i-2].c[1];
	mySinCos[l+i+1].c[2] = zcos*mySinCos[l+i-1].c[2]-zsin*mySinCos[l+i-2].c[2];
      }
    }
    
    //int nk = myK.size();
    int lastH = Constant::MAX_INT;
    int lastK = Constant::MAX_INT;
#endif
    // atomic virial
    Real virialxx = 0.0;
    Real virialxy = 0.0;
    Real virialxz = 0.0;
    Real virialyy = 0.0;
    Real virialyz = 0.0;
    Real virialzz = 0.0;

    // molecular virial
    bool doMolVirial = energies->molecularVirial();
    bool doVirial = energies->virial();

    for(unsigned int l=from;l<to;l++){
      // Energy
      Vector3D k(myK[l]);
      Real sumSin = 0.0;
      Real sumCos = 0.0;
      Real kSquared = myKSquared[l];

#ifndef USE_EWALD_NO_SINCOS_TABLE
      int indexH = myKInt[l].h;
      int indexK = myKInt[l].k;
      int indexKabs = abs(indexK);
      int indexKsign = 1;
      if(indexK < 0)
	indexKsign = -1;
      int indexL = myKInt[l].l;
      int indexLabs = abs(indexL);
      int indexLsign = 1;
      if(indexL < 0)
	indexLsign = -1;
    
      // Precompute and cache sin/ cos for h and k
      // using the precompute table of sin/ cos
      // Hit/miss rate vary from 6:1 to 12:1
      if(indexH != lastH || indexK != lastK){
	for(unsigned int i=0;i<atomCount;i++){
	  Real xsin =            mySinCos[i*myHKLmax*2+indexH*2  ].c[0];
	  Real xcos =            mySinCos[i*myHKLmax*2+indexH*2+1].c[0];
	  Real ysin = indexKsign*mySinCos[i*myHKLmax*2+indexKabs*2  ].c[1];
	  Real ycos =            mySinCos[i*myHKLmax*2+indexKabs*2+1].c[1];
	
	  myLastSinCos[i*2  ] = xsin*ycos + xcos*ysin;
	  myLastSinCos[i*2+1] = xcos*ycos - xsin*ysin;
	}
	lastH = indexH;
	lastK = indexK;
      }
#endif    
      for(unsigned int i=0;i<atomCount;i++){
#ifdef USE_EWALD_NO_SINCOS_TABLE
	Real qi = realTopo->atoms[i].scaledCharge;
	// It does not matter if coordinates are not in the minimal image since
	// they are multiplied by 2PI/l, which is a shift of 2PI of a. 
	Real a = k.dot(boundaryConditions.minimalPosition((*positions)[i]));
	//Real a = k.dot((*positions)[i]);
	Real sinA = qi*sin(a);
	Real cosA = qi*cos(a);
#else
      
	Real xysin = myLastSinCos[i*2  ];
	Real xycos = myLastSinCos[i*2+1];
      
	Real zsin = indexLsign*mySinCos[i*myHKLmax*2+indexLabs*2  ].c[2];
	Real zcos =            mySinCos[i*myHKLmax*2+indexLabs*2+1].c[2];
      
	Real sinA = xysin*zcos + xycos*zsin;
	Real cosA = xycos*zcos - xysin*zsin;
#endif
	mySinCosA[i*2  ] = sinA;
	mySinCosA[i*2+1] = cosA;
	sumSin += sinA;
	sumCos += cosA;
      }

      // Energy
      Real b = 1.0/kSquared*exp(-kSquared*myAlphaSquaredr/4.0);
      Real e = b*(sumSin*sumSin+sumCos*sumCos);
      energy += e;
    
      // Virial
      if(doVirial){
	Real c = 2.0*(1.0/kSquared+myFac);
	virialxx += e * (1.0-c*k.c[0] * k.c[0]);
	virialxy -= e * c * k.c[0] * k.c[1];
	virialxz -= e * c * k.c[0] * k.c[2];
	virialyy += e * (1.0-c*k.c[1] * k.c[1]);
	virialyz -= e * c * k.c[1] * k.c[2];
	virialzz += e * (1.0-c*k.c[2] * k.c[2]);
      }
      // Force, F_i 
      Real a = 8.0*M_PI*myVr*b;
      for(unsigned int i=0;i<atomCount;i++){

	// compute the force on atom i from the reciprocal space part
	Vector3D fi(k*(a*(mySinCosA[i*2  ]*sumCos - mySinCosA[i*2+1]*sumSin)));
	(*forces)[i] += fi;     

	// compute the reciprocal space contribution to the molecular virial
	// this expression is taken from Alejandre, Tildesley, and Chapela, J. Chem. Phys. 102 (11), 4574.
	if(doMolVirial){
	  // get the ID# of the molecule to which this atom belongs
	  int Mi = realTopo->atoms[i].molecule;
	  
	  // compute the vector from atom i to the center of mass of the molecule
	  Vector3D ria(boundaryConditions.minimalPosition((*positions)[i]));
	  Vector3D mri(realTopo->boundaryConditions.minimalDifference(ria,realTopo->molecules[Mi].position));

	  // compute the reciprocal space contribution to the molecular virial
	  // this expression is taken from Darden, et al. J. Chem. Phys. 103 (19), 8577.
	  energies->addMolVirial(fi,mri);
	}
      }
    }

    Real c = 4.0*M_PI*myVr;
    reciprocalEnergy += c*energy;

    // atomic virial
    if(doVirial){
      (*energies)[ScalarStructure::VIRIALXX] += c*virialxx;
      (*energies)[ScalarStructure::VIRIALXY] += c*virialxy;
      (*energies)[ScalarStructure::VIRIALXZ] += c*virialxz;
      (*energies)[ScalarStructure::VIRIALYX] += c*virialxy;
      (*energies)[ScalarStructure::VIRIALYY] += c*virialyy;
      (*energies)[ScalarStructure::VIRIALYZ] += c*virialyz;
      (*energies)[ScalarStructure::VIRIALZX] += c*virialxz;
      (*energies)[ScalarStructure::VIRIALZY] += c*virialyz;
      (*energies)[ScalarStructure::VIRIALZZ] += c*virialzz;
    }
    // molecular virial
    if(doMolVirial){
      (*energies)[ScalarStructure::MOLVIRIALXX] += c*virialxx;
      (*energies)[ScalarStructure::MOLVIRIALXY] += c*virialxy;
      (*energies)[ScalarStructure::MOLVIRIALXZ] += c*virialxz;
      (*energies)[ScalarStructure::MOLVIRIALYX] += c*virialxy;
      (*energies)[ScalarStructure::MOLVIRIALYY] += c*virialyy;
      (*energies)[ScalarStructure::MOLVIRIALYZ] += c*virialyz;
      (*energies)[ScalarStructure::MOLVIRIALZX] += c*virialxz;
      (*energies)[ScalarStructure::MOLVIRIALZY] += c*virialyz;
      (*energies)[ScalarStructure::MOLVIRIALZZ] += c*virialzz;
    }
#if defined(DEBUG_EWALD_TIMING)
    myReciprocal.stop();
#endif
  }


  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void NonbondedFullEwaldSystemForce<TBoundaryConditions,
				     TCellManager,
				     real,
				     reciprocal,
				     correction,
				     TSwitchingFunction>::correctionTerm(const RealTopologyType* realTopo,
									 const Vector3DBlock* positions, 
									 Vector3DBlock* forces, 
									 ScalarStructure* energies,
									 Real& intraMolecularEnergy,
									 unsigned int from, unsigned int to){
#if defined(DEBUG_EWALD_TIMING)
    myIntra.start();
#endif
    // Intra-molecular term
    bool doVirial = energies->virial();
    const std::vector<ExclusionPair>& exclusions = realTopo->exclusions.getTable();
    for(unsigned int i=from;i<to;i++){
      ExclusionPair excl = exclusions[i];
      Real rSquared;
      Vector3D rij(realTopo->boundaryConditions.minimalDifference((*positions)[excl.a1],(*positions)[excl.a2],rSquared));
      Real qq = realTopo->atoms[excl.a1].scaledCharge*realTopo->atoms[excl.a2].scaledCharge;
      if (excl.excl == EXCLUSION_MODIFIED)
	qq *= 1-realTopo->coulombScalingFactor;
      Real r  = sqrt(rSquared);
      Real rr = 1/r;
      Real e = erf(myAlpha*r)*rr;
      // Intra-molecular selv energy
      intraMolecularEnergy -= qq*e;
      // Intra-molecular selv force
      Vector3D fij(rij*(qq*(my2AlphaPI*exp(-myAlphaSquared*rSquared)-e)*rr*rr));
      (*forces)[excl.a1] -= fij;
      (*forces)[excl.a2] += fij;
      if(doVirial)
	energies->addVirial(fij,rij);
    }
#if defined(DEBUG_EWALD_TIMING)
    myIntra.stop();
#endif
  }

  // Surface diplole term
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void NonbondedFullEwaldSystemForce<TBoundaryConditions,
				     TCellManager,
				     real,
				     reciprocal,
				     correction,
				     TSwitchingFunction>::surfaceDipoleTerm(const RealTopologyType* realTopo,
									    const Vector3DBlock* positions, 
									    Vector3DBlock* forces, 
									    ScalarStructure* energies,
									    Real& surfaceDipoleEnergy){
#if defined(DEBUG_EWALD_TIMING)
    mySurface.start();
#endif
    const unsigned int atomCount = realTopo->atoms.size();
    Vector3D sum(0,0,0);
    for(unsigned int i=0;i<atomCount;i++)
      sum += boundaryConditions.minimalPosition((*positions)[i])*realTopo->atoms[i].scaledCharge;
    // Energy
    surfaceDipoleEnergy = 2.0/3.0*M_PI*myVr*sum.normSquared();
  
    Real virialxx = 0.0;
    Real virialxy = 0.0;
    Real virialxz = 0.0;
    Real virialyy = 0.0;
    Real virialyz = 0.0;
    Real virialzz = 0.0;
    bool doVirial = energies->virial();

    // Force, F_i and virial_i (not confirmed)
    sum *= 2.0/3.0*M_PI*myVr;
    for(unsigned int i=0;i<atomCount;i++){
      Vector3D force(sum*realTopo->atoms[i].scaledCharge);
      Vector3D ri(boundaryConditions.minimalPosition((*positions)[i]));
      (*forces)[i] += force;
      if(doVirial){
	virialxx += force.c[0]*ri.c[0];
	virialxy += force.c[0]*ri.c[1];
	virialxz += force.c[0]*ri.c[2];
	virialyy += force.c[1]*ri.c[1];
	virialyz += force.c[1]*ri.c[2];
	virialzz += force.c[2]*ri.c[2];
      }
    }
    if(doVirial){
      (*energies)[ScalarStructure::VIRIALXX] += virialxx;
      (*energies)[ScalarStructure::VIRIALXY] += virialxy;
      (*energies)[ScalarStructure::VIRIALXZ] += virialxz;
      (*energies)[ScalarStructure::VIRIALYX] += virialxy;
      (*energies)[ScalarStructure::VIRIALYY] += virialyy;
      (*energies)[ScalarStructure::VIRIALYZ] += virialyz;
      (*energies)[ScalarStructure::VIRIALZX] += virialxz;
      (*energies)[ScalarStructure::VIRIALZY] += virialyz;
      (*energies)[ScalarStructure::VIRIALZZ] += virialzz;
    }
#if defined(DEBUG_EWALD_TIMING)
    mySurface.stop();
#endif
  }


  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  std::string NonbondedFullEwaldSystemForce<TBoundaryConditions,
					    TCellManager,
					    real,
					    reciprocal,
					    correction,
					    TSwitchingFunction>::getIdNoAlias() const{
    return (" -algorithm "+ getKeyword() + 
	    std::string((real)       ? std::string(" -real")       : std::string("")) + 
	    std::string((reciprocal) ? std::string(" -reciprocal") : std::string("")) +
	    std::string((correction) ? std::string(" -correction") : std::string("")) +	  
	    std::string((TSwitchingFunction::getId() != CutoffSwitchingFunction::getId()) ? 
			std::string(std::string(" -switchingFunction " + TSwitchingFunction::getId())) : std::string("")));
  }

  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void NonbondedFullEwaldSystemForce<TBoundaryConditions,
				     TCellManager,
				     real,
				     reciprocal,
				     correction,
				     TSwitchingFunction>::getParameters(std::vector<Parameter>& parameters) const{
    parameters.push_back(Parameter("-alpha",Value(myAlpha),-1.0,Text("splitting")));    
    parameters.push_back(Parameter("-accuracy",Value(myAccuracy,ConstraintValueType::Positive()),0.00001));
    if(TBoundaryConditions::VACUUM)
      parameters.push_back(Parameter("-j",Value(myExpansionFactor,ConstraintValueType::Positive()),3.0));
  }

  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  Force*  NonbondedFullEwaldSystemForce<TBoundaryConditions,
					TCellManager,
					real,
					reciprocal,
					correction,
					TSwitchingFunction>::doMake(const std::vector<Value>& values) const{
    
    Real alpha           = values[0];
    Real accuracy        = values[1];
    Real expansionFactor = (TBoundaryConditions::VACUUM?(Real)values[2]:3.0);
    std::string err      = "";
    if(!values[0].valid())
      err +=" alpha \'"+values[0].getString()+"\' not valid.";

    if(!values[1].valid())
      err +=" accuracy \'"+values[1].getString()+"\' not valid.";

    if(TBoundaryConditions::VACUUM && !values[2].valid())
      err +=" expansionFactor \'"+values[2].getString()+"\' not valid.";
    else if(expansionFactor <= 1.0)
      err += getKeyword() + " simulation box expansion factor (="+toString(expansionFactor)+") > 1.0.";

    if(!err.empty()){
      err = " force "+getKeyword()+" :"+err;
      return NULL;
    }
    return (new NonbondedFullEwaldSystemForce(alpha, accuracy,expansionFactor));
  }
}
#endif /* NONBONDEDFULLEWALDSYSTEMFORCE_H */
