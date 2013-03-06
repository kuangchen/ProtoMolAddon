/* -*- c++ -*- */
#ifndef NONBONDEDPMEWALDSYSTEMFORCE_H
#define NONBONDEDPMEWALDSYSTEMFORCE_H

#include "Grid.h"
#include <protomol/force/system/SystemForce.h>
#include <protomol/base/TimerStatistic.h>
#include <protomol/parallel/Parallel.h>
#include <protomol/type/SimpleTypes.h>

using namespace ProtoMol::Report;

//#define DEBUG_PME_TIMING
//#define DEBUG_PME_ENERGIES
//#define USE_PME_EXACT_ERF
namespace ProtoMol {
  //_________________________________________________________________ NonbondedPMEwaldSystemForce
  
  template<class TBoundaryConditions, 
	   class TCellManager,
	   bool  real,
	   bool  reciprocal,
	   bool  correction,
	   class TInterpolation,
	   class TSwitchingFunction>
  class NonbondedPMEwaldSystemForce: public SystemForce {
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
    NonbondedPMEwaldSystemForce();
    NonbondedPMEwaldSystemForce(unsigned int nx, unsigned int ny, unsigned int nz, unsigned int order, Real cutoff, Real accuracy, Real alpha, Real expansionFactor);

    virtual ~NonbondedPMEwaldSystemForce();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NonbondedPMEwaldSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

  private:
    void initialize(const RealTopologyType* realTopo,
		    const Vector3DBlock* positions);

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
			Real& reciprocalEnergy);

    void reciprocalTermParallel(const RealTopologyType* realTopo,
				const Vector3DBlock* positions, 
				Vector3DBlock* forces, 
				ScalarStructure* energies,
				Real& reciprocalEnergy);

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
    virtual std::string getKeyword() const{return "PMEwald";}
    virtual void uncache(){myCached=false;};
  private:
    virtual Force* doMake(const std::vector<Value>& values) const;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual std::string getIdNoAlias() const;
    virtual void getParameters(std::vector<Parameter>& parameters) const;
    virtual unsigned int getParameterSize() const{return (TBoundaryConditions::PERIODIC ? 7:8);}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  private:
    bool myCached;
    Real myExpansionFactor;
    Real myTRate;             //
    Real myAccuracy;

    Real myAlpha;             // 
    bool myAlphaDefault;

    Real myAlphaSquared;
    Real myAlphaSquaredr;
    Real my2AlphaPI;

    Real myRc;                // Cutoff real term
    Real myRcSquared;         // Cutoff squared real term
    Real myKc;                // Cutoff reciprocal term

    Real myLX, myLY, myLZ;
    Real myLXr,myLYr,myLZr;
    Real myV,myVr;
    Vector3D myOrigin;
    unsigned int myNX;
    unsigned int myNY;
    unsigned int myNZ;
    Grid<TInterpolation> myGrid;
    unsigned int myInterOrder;

    Real myPointSelfEnergy;   // Precomputed energy terms
    Real myChargedSystemEnergy;

    PeriodicBoundaryConditions boundaryConditions;

    TSwitchingFunction switchingFunction;
    EnumeratorType enumerator;
    std::vector<Vector3D> myLattice;
#if defined(DEBUG_PME_TIMING)
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
	    class TInterpolation,
	    class TSwitchingFunction>
  NonbondedPMEwaldSystemForce<TBoundaryConditions,
			      TCellManager,
			      real,
			      reciprocal,
			      correction,
			      TInterpolation,
			      TSwitchingFunction>::NonbondedPMEwaldSystemForce() : SystemForce(),myCached(false),myExpansionFactor(3.0),myTRate(0.0),myAccuracy(0.0),myAlpha(-1.0),myAlphaDefault(true),myRc(0.0),myKc(0.0),myV(-1.0),myNX(0),myNY(0),myNZ(0),myInterOrder(0) {
#if defined(DEBUG_PME_TIMING)
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
	    class TInterpolation,
	    class TSwitchingFunction>
  NonbondedPMEwaldSystemForce<TBoundaryConditions,
			      TCellManager,
			      real,
			      reciprocal,
			      correction,
			      TInterpolation,
			      TSwitchingFunction>::NonbondedPMEwaldSystemForce(unsigned int nx, unsigned int ny, unsigned int nz, unsigned int order, Real cutoff, Real accuracy, Real alpha, Real expansionFactor) : SystemForce(),myCached(false),myExpansionFactor(expansionFactor),myTRate(0.0),myAccuracy(accuracy),myAlpha(alpha),myAlphaDefault(alpha<= 0.0),myRc(cutoff),myKc(0.0),myV(-1.0),myNX(nx),myNY(ny),myNZ(nz),myInterOrder(order){
    switchingFunction=TSwitchingFunction(myRc);
#if defined(DEBUG_PME_TIMING)
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
	    class TInterpolation,
	    class TSwitchingFunction>
  NonbondedPMEwaldSystemForce<TBoundaryConditions,
			      TCellManager,
			      real,
			      reciprocal,
			      correction,
			      TInterpolation,
			      TSwitchingFunction>::~NonbondedPMEwaldSystemForce(){
#if defined(DEBUG_PME_TIMING)
    if(boundaryConditions.getVolume() > Constant::EPSILON){
      report.setf(std::ios::showpoint|std::ios::fixed);
      report << allnodes << plain <<"Timing ("
             <<Parallel::getId()<<") PME:"
	     <<" real:" <<myReal.getTime().getProcessTime()<<"[s]"
	     <<", reciprocal:"<<myReciprocal.getTime().getProcessTime()<<"[s]"
	     <<", intra:"<<myIntra.getTime().getProcessTime()<<"[s]"
	     <<", dipole:"<<mySurface.getTime().getProcessTime()<<"[s]"
	     <<", Rc="<<myRc
	     <<", Kc="<<myKc
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
	    class TInterpolation,
	    class TSwitchingFunction>
  void NonbondedPMEwaldSystemForce<TBoundaryConditions,
				   TCellManager,
				   real,
				   reciprocal,
				   correction,
				   TInterpolation,
				   TSwitchingFunction>::evaluate(const GenericTopology* topo, 
								 const Vector3DBlock* positions, 
								 Vector3DBlock* forces, 
								 ScalarStructure* energies) {
    const RealTopologyType* realTopo = (const RealTopologyType *)topo;  

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

      realTerm(realTopo,positions,forces,energies,realEnergy,
	       realTopo->cellLists.size());
    }  

    // Reciprocal-space term
    Real reciprocalEnergy = 0.0;
    if(reciprocal)
      reciprocalTerm(realTopo,positions,forces,energies,reciprocalEnergy);

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

#if defined(DEBUG_PME_ENERGIES)
    report.setf(std::ios::showpoint|std::ios::fixed);
    report << plain <<"PME: point="<<myPointSelfEnergy
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
	    class TInterpolation,
	    class TSwitchingFunction>
  void NonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				   real,reciprocal,correction,TInterpolation,TSwitchingFunction>::parallelEvaluate(const GenericTopology* topo, 
														   const Vector3DBlock* positions, 
														   Vector3DBlock* forces, 
														   ScalarStructure* energies) {

    bool dontHint = (myV >= 0.0 && TBoundaryConditions::PERIODIC);
    if(dontHint)
      Report::report << Report::donthint;

    const RealTopologyType* realTopo = (const RealTopologyType *)topo;  

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
      reciprocalTermParallel(realTopo,positions,forces,energies,reciprocalEnergy);      
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

#if defined(DEBUG_PME_ENERGIES)
    report << allnodes << plain <<"PME: point="<<myPointSelfEnergy
	   <<", charged="<<myChargedSystemEnergy
	   <<", real="<<realEnergy
	   <<", reciprocal="<<reciprocalEnergy
	   <<", intra="<<intraMolecularEnergy
	   <<", surface="<<surfaceDipoleEnergy
	   <<", total="<<e 
	   <<endr;
#endif

    if(dontHint)
      Report::report << Report::dohint;

  }

  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  unsigned int NonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
					   real,reciprocal,correction,TInterpolation,TSwitchingFunction>::numberOfBlocks(const GenericTopology* topo,
															 const Vector3DBlock* positions){
    unsigned int n = 0;
    if(correction)
      n +=std::min(static_cast<int>(topo->exclusions.getTable().size()),static_cast<int>(Parallel::getAvailableNum()));

    if(reciprocal)
      ;

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
	    class TInterpolation,
	    class TSwitchingFunction>
  void NonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				   real,reciprocal,correction,TInterpolation,
				   TSwitchingFunction>::initialize(const RealTopologyType* realTopo, 
								   const Vector3DBlock* positions) {
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
	report << hint << "The boundaries for PME force evaluation were re-sized to "<<boundaryConditions.getMin()<<"-"<<boundaryConditions.getMax()<<"."<<endr;
      }
      else {
	// The particles are just inside +-10% our previous simulation box. No update needed
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
    if(myAlphaDefault){
      myAlpha = 1.0;
      //myRc = 6.5;
      while ( erfc(myAlpha*myRc)/myRc >= myAccuracy ) 
	myAlpha *= 2.0;

      Real low = 0.;
      Real high = myAlpha;
      for(unsigned int i=0; i<100; ++i){
	myAlpha = 0.5 * (low + high);
	if ( erfc(myAlpha*myRc)/myRc >= myAccuracy ) {
	  low = myAlpha;
	} else {
	  high = myAlpha;
	}
      }
      //myAlpha         = sqrt(M_PI)*pow(myTRate*atomCount/(power<2>(myV)),1.0/6.0);
    }
    Real p          = -log(myAccuracy);
    myKc            = 2.0*myAlpha*sqrt(p);
    //myRc            = sqrt(p)/myAlpha;
    myRcSquared     = myRc*myRc;


    myAlphaSquared  = myAlpha*myAlpha;
    myAlphaSquaredr = 1.0/myAlphaSquared;
    my2AlphaPI      = 2.0*myAlpha/sqrt(M_PI);
  
    // The grid
    if(reciprocal)
      myGrid.initialize(myLX,myLY,myLZ,myAlpha,myNX,myNY,myNZ,myInterOrder,atomCount);

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
    // Charged sytem energy
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
  
    report << hint <<"PME";
#ifdef USE_PME_EXACT_ERF
    report << "(Exact)";
#endif  
    report <<": alpha="<<toString(myAlpha)<<", V="<<myV<<", Rc="
	   <<toString(myRc)<<", Kc ="<<toString(myKc)<<", n="<<myLattice.size()<<", accuracy="
	   <<myAccuracy<<", interpolation="<<TInterpolation::keyword
	   <<", order="<<myInterOrder<<"."<<endr;

    // Now we have all pre-computed stuff
    if(TBoundaryConditions::PERIODIC)
      myCached = true;

    if(dontHint)
      Report::report << Report::dohint;
  }


  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void NonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				   real,reciprocal,correction,TInterpolation,
				   TSwitchingFunction>::realTerm(const RealTopologyType* realTopo,
								 const Vector3DBlock* positions, 
								 Vector3DBlock* forces, 
								 ScalarStructure* energies,
								 Real& realEnergy,
								 unsigned int n) {
#if defined(DEBUG_PME_TIMING)
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
#ifndef USE_PME_EXACT_ERF
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
	    else if(doVirial) 
	      energies->addVirial(fij,rij);
	  
	    excl = EXCLUSION_NONE;
	  }
	}
	(*forces)[i] += fi;
      }
    }
#if defined(DEBUG_PME_TIMING)
    myReal.stop();
#endif
  }



  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void NonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				   real,reciprocal,correction,TInterpolation,
				   TSwitchingFunction>::reciprocalTermParallel(const RealTopologyType* realTopo,
									       const Vector3DBlock* positions, 
									       Vector3DBlock* forces, 
									       ScalarStructure* energies,
									       Real& reciprocalEnergy) {

#if defined(DEBUG_PME_TIMING)
    myReciprocal.start();
#endif
    const unsigned int atomCount = realTopo->atoms.size();
    const unsigned int nBlocks = Parallel::getAvailableNum();
    const unsigned int block = Parallel::getAvailableId();
    const unsigned int i0 = (atomCount*block)/nBlocks;
    const unsigned int i1 = (atomCount*(block+1))/nBlocks;
    Real* begin = NULL;
    Real* end   = NULL;

    // Anterpolate
    myGrid.clear();
    for(unsigned int i=i0;i<i1;i++)
      myGrid.anterpolateCharge(realTopo->atoms[i].scaledCharge,(*positions)[i],i);
    myGrid.getQ(begin,end);
    TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
    Parallel::reduceSlaves(begin,end);
    TimerStatistic::timer[TimerStatistic::FORCES] -= TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();

    // FFT back
    if(FFTComplex::isParallel()){
      myGrid.fftBack();
    }
    else {
      if(block == 0)
	myGrid.fftBack();
      myGrid.getQ(begin,end);
      TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
      Parallel::bcastSlaves(begin,end);
      TimerStatistic::timer[TimerStatistic::FORCES] -= TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
    }
  
    // Sum
    Real energy = myGrid.scalarSum(energies,block,nBlocks);
    myGrid.getQ(begin,end);
    TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
    Parallel::reduceSlaves(begin,end);
    TimerStatistic::timer[TimerStatistic::FORCES] -= TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
    //  myGrid.print();

    // FFT forward
    if(FFTComplex::isParallel()){
      myGrid.fftForward();
    }
    else{
      if(block == 0)
	myGrid.fftForward();
      myGrid.getQ(begin,end);
      TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
      Parallel::bcastSlaves(begin,end);
      TimerStatistic::timer[TimerStatistic::FORCES] -= TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
    }

    // Interploate
    bool doMolVirial = energies->molecularVirial();
    for(unsigned int i=i0;i<i1;i++) {
      // compute the force on atom i from the reciprocal space part
      Vector3D fi;
      myGrid.interpolateForce(realTopo->atoms[i].scaledCharge,i,fi);
      (*forces)[i] += fi;

      if(doMolVirial){
	// compute the vector from atom i to the center of mass of the molecule
	Vector3D ria(boundaryConditions.minimalPosition((*positions)[i]));
	Vector3D mri(realTopo->boundaryConditions.minimalDifference(ria,realTopo->molecules[realTopo->atoms[i].molecule].position));

	// compute the reciprocal space contribution to the molecular virial
	// this expression is taken from Darden, et al. J. Chem. Phys. 103 (19), 8577.
	energies->addMolVirial(fi,mri);
      }
    } // end loop over atoms (i)

    reciprocalEnergy += energy;
      
#if defined(DEBUG_PME_TIMING)
    myReciprocal.stop();
#endif
  }
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void NonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				   real,reciprocal,correction,TInterpolation,
				   TSwitchingFunction>::reciprocalTerm(const RealTopologyType* realTopo,
								       const Vector3DBlock* positions, 
								       Vector3DBlock* forces, 
								       ScalarStructure* energies,
								       Real& reciprocalEnergy) {

#if defined(DEBUG_PME_TIMING)
    myReciprocal.start();
#endif
    const unsigned int atomCount = realTopo->atoms.size();
      
    myGrid.clear();
    for(unsigned int i=0;i<atomCount;i++)
      myGrid.anterpolateCharge(realTopo->atoms[i].scaledCharge,(*positions)[i],i);

    myGrid.fftBack();

    Real energy = myGrid.scalarSum(energies);
    //  myGrid.print();

    myGrid.fftForward();

    bool doMolVirial = energies->molecularVirial();
    for(unsigned int i=0;i<atomCount;i++) {

      // compute the force on atom i from the reciprocal space part
      Vector3D fi;
      myGrid.interpolateForce(realTopo->atoms[i].scaledCharge,i,fi);
      (*forces)[i] += fi;
      //myGrid.interpolateForce(realTopo->atoms[i].scaledCharge,i,(*forces)[i]);

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
  } // end loop over atoms

    reciprocalEnergy += energy;
      
#if defined(DEBUG_PME_TIMING)
    myReciprocal.stop();
#endif
  }


  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void NonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				   real,reciprocal,correction,TInterpolation,
				   TSwitchingFunction>::correctionTerm(const RealTopologyType* realTopo,
								       const Vector3DBlock* positions, 
								       Vector3DBlock* forces, 
								       ScalarStructure* energies,
								       Real& intraMolecularEnergy,
								       unsigned int from, unsigned int to){
#if defined(DEBUG_PME_TIMING)
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

#if defined(DEBUG_PME_TIMING)
    myIntra.stop();
#endif
  }

  // Surface diplole term
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void NonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				   real,reciprocal,correction,TInterpolation,
				   TSwitchingFunction>::surfaceDipoleTerm(const RealTopologyType* realTopo,
									  const Vector3DBlock* positions, 
									  Vector3DBlock* forces, 
									  ScalarStructure* energies,
									  Real& surfaceDipoleEnergy){
#if defined(DEBUG_PME_TIMING)
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
#if defined(DEBUG_PME_TIMING)
    mySurface.stop();
#endif
  }


  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  std::string NonbondedPMEwaldSystemForce<TBoundaryConditions, TCellManager,
					  real,reciprocal,correction,TInterpolation,TSwitchingFunction>::getIdNoAlias() const{
    return (CoulombForce::getId() + " -algorithm "+ getKeyword() + 
	    std::string((real)       ? std::string(" -real")       : std::string("")) + 
	    std::string((reciprocal) ? std::string(" -reciprocal") : std::string("")) +
	    std::string((correction) ? std::string(" -correction") : std::string("")) +
	    " -interpolation "+TInterpolation::getKeyword()+
	    std::string((TSwitchingFunction::getId() != CutoffSwitchingFunction::getId()) ? 
			std::string(std::string(" -switchingFunction " + TSwitchingFunction::getId())) : std::string("")));
  }

  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void NonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				   real,reciprocal,correction,TInterpolation,TSwitchingFunction>::getParameters(std::vector<Parameter>& parameters) const{
    parameters.push_back(Parameter("-gridsize",Value(myNX,ConstraintValueType::Positive())));
    parameters.push_back(Parameter("",Value(myNY,ConstraintValueType::Positive())));
    parameters.push_back(Parameter("",Value(myNZ,ConstraintValueType::Positive())));
    parameters.push_back(Parameter("-cutoff",Value(myRc,ConstraintValueType::Positive()),Text("Rc cutoff")));
    parameters.push_back(Parameter("-order",Value(myInterOrder,ConstraintValueType::Positive()),4,Text("interpolation")));
    parameters.push_back(Parameter("-accuracy",Value(myAccuracy,ConstraintValueType::Positive()),1.e-6));
    parameters.push_back(Parameter("-alpha",Value(myAlpha),-1.0,Text("splitting")));    
    if(TBoundaryConditions::VACUUM)
      parameters.push_back(Parameter("-j",Value(myExpansionFactor),3.0));
  }

  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  Force* NonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				     real,reciprocal,correction,TInterpolation,TSwitchingFunction>::doMake(const std::vector<Value>& values) const{
    int  nx              = values[0];
    int  ny              = values[1];
    int  nz              = values[2];
    Real cutoff          = values[3];
    int  order           = values[4];
    Real accuracy        = values[5];
    Real alpha           = values[6];
    Real expansionFactor = (TBoundaryConditions::VACUUM?(Real)values[7]:3.0);
    std::string err = "";
    if(TBoundaryConditions::VACUUM && !values[7].valid())
      err +=" expansionFactor \'"+values[2].getString()+"\' for "+getId()+" not valid.";
    if(expansionFactor <= 1.0)
      err += getKeyword() + " simulation box expansion factor (="+toString(expansionFactor)+") > 1.0.";

    if(order < 2 || (order % 2) != 0 || !values[4].valid())
      err += getKeyword() + " order (="+values[4].getString()+") > 1 and even.";

    if(!values[0].valid() || !values[1].valid() || !values[2].valid() || nx < order || ny < order || nz < order)
      err += getKeyword() + " force: "+values[4].getString()+" <= nx (="+values[0].getString()+"), "+values[4].getString()+" <= ny (="+values[1].getString()+"), "+values[4].getString()+" <= nz (="+values[2].getString()+").";

    if(cutoff <= 0.0 || !values[3].valid())
      err += getKeyword() + " cutoff (="+values[3].getString()+") > 0.";

    if(accuracy <= 0.0 || !values[5].valid())
      err += getKeyword() + " accuracy (="+values[5].getString()+") > 0.";

    if(!values[6].valid())
      err +=" alpha \'"+values[6].getString()+"\' not valid.";

    if(!err.empty()){
      err += " force "+getKeyword()+" :"+err;
      report << error << err << endr;
      return NULL;
    }

    return (new NonbondedPMEwaldSystemForce((unsigned int)nx,(unsigned int)ny,(unsigned int)nz,(unsigned int)order,cutoff,accuracy,alpha,expansionFactor));
  }
}
#endif /* NONBONDEDPMEWALDSYSTEMFORCE_H */
