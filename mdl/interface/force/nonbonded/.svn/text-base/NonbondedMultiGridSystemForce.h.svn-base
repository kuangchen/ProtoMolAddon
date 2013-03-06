/* -*- c++ -*- */
#ifndef NONBONDEDMULTIGRIDSYSTEMFORCE_H
#define NONBONDEDMULTIGRIDSYSTEMFORCE_H

/*#include "SystemForce.h"
#include "mathutilities.h"
#include "systemutilities.h"
#include "NonbondedMultiGridSystemForceBase.h"
#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
#include "TimerStatistic.h"*/

#include <protomol/parallel/Parallel.h>
#include <protomol/force/system/SystemForce.h>
#include <protomol/base/TimerStatistic.h>
#include <protomol/type/SimpleTypes.h>

using namespace ProtoMol::Report;

//#define DEBUG_MULTIGRID
//#define DEBUG_MULTIGRID_TIMING
// Same define's in "MultiGrid.h".
#include "MultiGrid.h"

namespace ProtoMol {
  //_________________________________________________________________ NonbondedMultiGridSystemForce
  
  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction,
	   bool smooth>
  class NonbondedMultiGridSystemForce: public SystemForce {
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
    NonbondedMultiGridSystemForce();
    NonbondedMultiGridSystemForce(unsigned int nx, unsigned int ny, unsigned int nz, unsigned int levels, 
				  Real s, unsigned int order, unsigned int ratio, Vector3D h, Vector3D origin );

    virtual ~NonbondedMultiGridSystemForce();
  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NonbondedMultiGridSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

  private:
    void initialize(const RealTopologyType* realTopo,
		    const Vector3DBlock* positions);
    void shortRangeTerm(const RealTopologyType* topo, 
			const Vector3DBlock* positions,   
			Vector3DBlock* forces, 
			Real& shortRange, 
			ScalarStructure* energies,
			unsigned int n);
    void longRangeTerm(const RealTopologyType* realTopo,
		       const Vector3DBlock* positions, 
		       Vector3DBlock* forces, 
		       Real& longRange,
		       ScalarStructure*);
    void correctionTerm(const RealTopologyType* realTopo,
			const Vector3DBlock* positions, 
			Vector3DBlock* forces, 
			Real& intraMolecularEnergy,
			ScalarStructure* energies,
			unsigned int from, unsigned int to);

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

    virtual std::string getKeyword() const{return "MultiGrid";}
    virtual void uncache(){myCached=false;};
  private:
    virtual Force* doMake(const std::vector<Value>& values) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual std::string getIdNoAlias() const;
    virtual void getParameters(std::vector<Parameter>& parameters) const;
    virtual unsigned int getParameterSize() const{return (1 + (TBoundaryConditions::PERIODIC && smooth?3:0) + (smooth?3:0) + (TBoundaryConditions::VACUUM && smooth?2:0));}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  private:
    bool myCached;
    Real myV;

    int myLevels;
    Real myS;
    Real myS2;
    Real myRS;
    unsigned int myNX;
    unsigned int myNY;
    unsigned int myNZ;
    MultiGrid<TInterpolation,TKernel,TBoundaryConditions::PERIODIC,TBoundaryConditions::PERIODIC,TBoundaryConditions::PERIODIC> myGrid;
    unsigned int myInterOrder;
    unsigned int myRatio;
    Vector3D myMin;
    Vector3D myMax;
    Vector3D myH;
    Vector3D myOrigin;
    Real myPointSelfEnergy;// Precomputed energy terms
    EnumeratorType enumerator;
#ifdef DEBUG_MULTIGRID_TIMING
    Timer myTimeInit;    
    Timer myTimeShort;   
    Timer myTimeLong;    
    Timer myTimeCorrection;
    Timer myTimeAnterpolate;
    std::vector<Timer> myTimeFineToCoarse;
    std::vector<Timer> myTimeCorrectionMG; 
    Timer myTimeDirect;
    std::vector<Timer> myTimeCoarseToFine;
    Timer myTimeInterpolate;
    long myCounterDirect;
    long myCounterCorrection;
#endif

  };

  //______________________________________________________________________ INLINES
  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  NonbondedMultiGridSystemForce<TBoundaryConditions,
				TCellManager,
				TInterpolation,
				TKernel,
				direct, 
				correction,
				smooth>::NonbondedMultiGridSystemForce() : 
				  SystemForce(),
				  myCached(false),
				  myV(-1.0),
				  myLevels(0),
				  myS(0.0),
				  myNX(0),
				  myNY(0),
				  myNZ(0),
				  myInterOrder(0),
				  myRatio(0),
				  myH(Vector3D(0,0,0)),
				  myOrigin(Vector3D(0,0,0)) {
#ifdef DEBUG_MULTIGRID_TIMING
    myTimeInit.reset();
    myTimeShort.reset();
    myTimeLong.reset();
    myTimeCorrection.reset();
    myTimeAnterpolate.reset();
    myTimeFineToCoarse.resize(0);
    myTimeCorrectionMG.resize(0); 
    myTimeDirect.reset();
    myTimeCoarseToFine.resize(0);
    myTimeInterpolate.reset();
    myCounterDirect = 0;
    myCounterCorrection = 0;
#endif
  }

  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  NonbondedMultiGridSystemForce<TBoundaryConditions,
				TCellManager,
				TInterpolation,
				TKernel,
				direct, 
				correction,
				smooth>::NonbondedMultiGridSystemForce(unsigned int nx, 
								       unsigned int ny, 
								       unsigned int nz, 
								       unsigned int levels, 
								       Real s, 
								       unsigned int order, 
								       unsigned int ratio, 
								       Vector3D h, 
								       Vector3D origin) : 
				  SystemForce(),
				  myCached(false),
				  myV(-1.0),
				  myLevels(levels),
				  myS(s),
				  myNX(nx),
				  myNY(ny),
				  myNZ(nz),
				  myInterOrder(order),
				  myRatio(ratio),
				  myH(h),
				  myOrigin(origin){
#ifdef DEBUG_MULTIGRID
    report << hint << "NonbondedMultiGridSystemForce:"<<(direct?"direct":"")<<" "<<(correction?"correction":"")<<" "<<(smooth?"smooth":"")<<" "<<endr;
#endif
#ifdef DEBUG_MULTIGRID_TIMING
    myTimeInit.reset();
    myTimeShort.reset();
    myTimeLong.reset();
    myTimeCorrection.reset();
    myTimeAnterpolate.reset();
    myTimeFineToCoarse.resize(levels,Timer());
    myTimeCorrectionMG.resize(levels,Timer()); 
    myTimeDirect.reset();
    myTimeCoarseToFine.resize(levels,Timer());
    myTimeInterpolate.reset();
    myCounterDirect = 0;
    myCounterCorrection = 0;
#endif
  }

  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  NonbondedMultiGridSystemForce<TBoundaryConditions,
				TCellManager,
				TInterpolation,
				TKernel,
				direct, 
				correction,
				smooth>::~NonbondedMultiGridSystemForce(){
#ifdef DEBUG_MULTIGRID_TIMING
    report.setf(std::ios::showpoint|std::ios::fixed);
    if((myTimeInit.getTime()+
	myTimeShort.getTime()+
	myTimeLong.getTime()+
	myTimeCorrection.getTime()).getProcessTime() > EPSILON){
      report << allnodes << plain <<"Timing ("
	     << Parallel::getId()<<") MultiGrid:"<<std::endl
	     <<"  Init:        "<< myTimeInit.getTime().getProcessTime()<<"[s]."<<std::endl
	     <<"  Short:       "<< myTimeShort.getTime().getProcessTime()<<"[s] ("<<myCounterDirect<<")."<<std::endl
	     <<"  Long:        "<<myTimeLong.getTime().getProcessTime()<<"[s] ("<<myGrid.getCounter()<<")."<<std::endl
	     <<"  Total:       "<<(myTimeInit.getTime()+myTimeShort.getTime()+myTimeLong.getTime()+myTimeCorrection.getTime()).getProcessTime()
	     <<" ("<<(myTimeShort.getTime()+myTimeLong.getTime()+myTimeCorrection.getTime()).getProcessTime()<<")[s]."<<std::endl<<std::endl
	     <<"  Correction:  "<< myTimeCorrection.getTime().getProcessTime()<<"[s] ("<<myCounterCorrection<<")."<<std::endl
	     <<"  Anterpolate: "<< myTimeAnterpolate.getTime().getProcessTime()<<"[s]."<<std::endl
	     <<"  Interpolate: "<< myTimeInterpolate.getTime().getProcessTime()<<"[s]."<<endr;
      for(int i=0;i<myLevels-1;i++){
	report << allnodes << plain << "  Level "<<i<<" :"
	       <<" FineToCoarse: "<< myTimeFineToCoarse[i].getTime().getProcessTime()<<"[s]"
	       <<" CoarseToFine: "<< myTimeCoarseToFine[i+1].getTime().getProcessTime()<<"[s]"
	       <<" Correction: "<< myTimeCorrectionMG[i].getTime().getProcessTime()<<"[s]."<<endr;
      }
      report << allnodes << plain << "  Level "<<myLevels-1<<" :"
	     <<" Direct:       "<< myTimeDirect.getTime().getProcessTime()<<"[s]."<<endr;
    }
#endif
  }

  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  void NonbondedMultiGridSystemForce<TBoundaryConditions,
				     TCellManager,
				     TInterpolation,
				     TKernel,
				     direct, 
				     correction,
				     smooth>::evaluate(const GenericTopology* topo, 
						       const Vector3DBlock* positions, 
						       Vector3DBlock* forces, 
						       ScalarStructure* energies) {

    const RealTopologyType* realTopo = (const RealTopologyType *)topo;  

    // Initialize data members and precompute tables & short cuts
#ifdef DEBUG_MULTIGRID_TIMING
    myTimeInit.start();
#endif
    if(!myCached)
      initialize(realTopo,positions);

    if(direct){
      realTopo->updateCellLists(positions);
      enumerator.initialize(realTopo, myS);
    }
    if(smooth){
      if(TBoundaryConditions::VACUUM){
	if(direct){
	  myGrid.updateSize(realTopo->min,realTopo->max);
	}
	else {
	  positions->boundingbox(myMin,myMax);
	  myGrid.updateSize(myMin,myMax);
	}
      }
    }

#ifdef DEBUG_MULTIGRID_TIMING
    myTimeInit.stop();
    myTimeLong.start();
#endif

    Real longRangeEnergy = 0.0;
    if(smooth)
      longRangeTerm(realTopo,positions,forces,longRangeEnergy,energies);

#ifdef DEBUG_MULTIGRID_TIMING
    myTimeLong.stop();
    myTimeShort.start();
#endif

    Real shortRangeEnergy = 0.0;
    if(direct)
      shortRangeTerm(realTopo,positions,forces,shortRangeEnergy,energies,
		     realTopo->cellLists.size());

#ifdef DEBUG_MULTIGRID_TIMING
    myTimeShort.stop();
    myTimeCorrection.start();
#endif

    Real intraMolecularEnergy = 0.0;
    Real pointSelfEnergy = 0.0;
    if(correction){
      correctionTerm(realTopo,positions,forces,intraMolecularEnergy,energies,
		     0,realTopo->exclusions.getTable().size());
      pointSelfEnergy = myPointSelfEnergy;
    }
#ifdef DEBUG_MULTIGRID_TIMING
    myTimeCorrection.stop();
#endif

    (*energies)[ScalarStructure::COULOMB] += 
      shortRangeEnergy+
      longRangeEnergy+
      intraMolecularEnergy+
      pointSelfEnergy;
#if defined(DEBUG_MULTIGRID)
    report << plain <<"MultiGrid : shortRangeEnergy="<<shortRangeEnergy<<", longRangeEnergy="<<longRangeEnergy<<", intraMolecularEnergy="<<intraMolecularEnergy<<", pointSelfEnergy="<<pointSelfEnergy<<", short+point="<<shortRangeEnergy+pointSelfEnergy<<", total="<<shortRangeEnergy+longRangeEnergy+intraMolecularEnergy+pointSelfEnergy<<"."<<endr;
#endif

  }


  //______________________________________________________________________ INLINES
  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  void NonbondedMultiGridSystemForce<TBoundaryConditions,
				     TCellManager,
				     TInterpolation,
				     TKernel,
				     direct, 
				     correction,
				     smooth>::parallelEvaluate(const GenericTopology* topo, 
							       const Vector3DBlock* positions, 
							       Vector3DBlock* forces, 
							       ScalarStructure* energies) {

    const RealTopologyType* realTopo = (const RealTopologyType *)topo;  

    //
    // Initialize data members and precompute tables & short cuts
    //
#ifdef DEBUG_MULTIGRID_TIMING
    myTimeInit.start();
#endif
    if(!myCached)
      initialize(realTopo,positions);

    if(direct){
      realTopo->updateCellLists(positions);
      enumerator.initialize(realTopo, myS);
    }
    if(smooth){
      if(TBoundaryConditions::VACUUM){
	if(direct){
	  myGrid.updateSize(realTopo->min,realTopo->max);
	}
	else {
	  positions->boundingbox(myMin,myMax);
	  myGrid.updateSize(myMin,myMax);
	}
      }
    }
#ifdef DEBUG_MULTIGRID_TIMING
    myTimeInit.stop();

    //
    // Long/smooth part
    //
    myTimeLong.start();
#endif
    Real longRangeEnergy = 0.0;
    if(smooth){
      longRangeTerm(realTopo,positions,forces,longRangeEnergy,energies);
    }
#ifdef DEBUG_MULTIGRID_TIMING
    myTimeLong.stop();


    //
    // Direct/short part
    //
    myTimeShort.start();
#endif
    Real shortRangeEnergy = 0.0;
    if(direct){
      unsigned int n = realTopo->cellLists.size();
      unsigned int count = Parallel::getNumberOfPackages(n);
    
      for(unsigned int i = 0;i<count;i++){
	unsigned int l = (n*(i+1))/count - (n*i)/count;
	if(Parallel::next())
	  shortRangeTerm(realTopo,positions,forces,shortRangeEnergy,energies,l);
	else
	  enumerator.nextNewPair(l);	
      }
    }
#ifdef DEBUG_MULTIGRID_TIMING
    myTimeShort.stop();

    //
    // Correction/exclusion and self point energy part
    //
    myTimeCorrection.start();
#endif
    Real intraMolecularEnergy = 0.0;
    Real pointSelfEnergy = 0.0;
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
	    correctionTerm(realTopo,positions,forces,intraMolecularEnergy,energies,from,to);
	  }
	}
      }
      if(Parallel::getAvailableId() == 0)
	pointSelfEnergy = myPointSelfEnergy;
    }
#ifdef DEBUG_MULTIGRID_TIMING
    myTimeCorrection.stop();
#endif

    //
    // Sum the enegries
    //
    (*energies)[ScalarStructure::COULOMB] += 
      shortRangeEnergy+
      longRangeEnergy+
      intraMolecularEnergy+
      pointSelfEnergy;
#if defined(DEBUG_MULTIGRID)
    report << allnodes << plain <<"MultiGrid : shortRangeEnergy="<<shortRangeEnergy<<", longRangeEnergy="<<longRangeEnergy<<", intraMolecularEnergy="<<intraMolecularEnergy<<", pointSelfEnergy="<<pointSelfEnergy<<", short+point="<<shortRangeEnergy+pointSelfEnergy<<", total="<<shortRangeEnergy+longRangeEnergy+intraMolecularEnergy+pointSelfEnergy<<"."<<endr;
#endif

  }

  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  void NonbondedMultiGridSystemForce<TBoundaryConditions,
				     TCellManager,
				     TInterpolation,
				     TKernel,
				     direct, 
				     correction,
				     smooth>::initialize(const RealTopologyType* realTopo, const Vector3DBlock* positions) {
  
    if(!realTopo->boundaryConditions.isOrthogonal())
      report << error << "[NonbondedMultiGridSystemForce::initialize] Not orthogonal, aborting."<<endr;

    bool dontHint = (myV >= 0.0 && TBoundaryConditions::PERIODIC);
    if(dontHint)
      Report::report << Report::donthint;
    myV = realTopo->getVolume(*positions);

    const unsigned int atomCount = realTopo->atoms.size();

    // Short cuts
    myS2 = myS*myS;
    myRS = 1.0/myS;
  
    // Dimension of the simulation box
    myMin = Vector3D(0.0,0.0,0.0);
    myMax = Vector3D(0.0,0.0,0.0);
    if(smooth){
      if(TBoundaryConditions::PERIODIC){
	myMin  = realTopo->boundaryConditions.getMin();
	myMax  = realTopo->boundaryConditions.getMax();
      }
      else {
	positions->boundingbox(myMin,myMax);
      }
      // The grid
      myGrid.initialize(positions->size(),myS,myLevels,myNX,myNY,myNZ,myInterOrder,myRatio,myMin,myMax,myH,myOrigin);
      //myGrid.printConst();
    }

    //
    // Point self-energy
    //
    myPointSelfEnergy = 0.0;
    if(correction){
      Real q = 0.0;

      for(unsigned int i=0;i<atomCount;i++)
	q += realTopo->atoms[i].scaledCharge*realTopo->atoms[i].scaledCharge;
    
      myPointSelfEnergy = -0.5*q*TKernel::smooth0(myS,myRS);
    }

    // Now we have all pre-computed stuff
    myCached = true;
    if(dontHint)
      Report::report << Report::dohint;
  }



  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  void NonbondedMultiGridSystemForce<TBoundaryConditions,
				     TCellManager,
				     TInterpolation,
				     TKernel,
				     direct, 
				     correction,
				     smooth>::longRangeTerm(const RealTopologyType* realTopo,
							    const Vector3DBlock* positions, 
							    Vector3DBlock* forces, 
							    Real& energy,
							    ScalarStructure* energies) {

    const unsigned int atomCount = realTopo->atoms.size();
    const unsigned int nBlocks = Parallel::getAvailableNum();
    const unsigned int block = Parallel::getAvailableId();
    const unsigned int i0 = (atomCount*block)/nBlocks;
    const unsigned int i1 = (atomCount*(block+1))/nBlocks;
    Real* begin = NULL;
    Real* end   = NULL;
    //report << allnodes <<debug << "("<<Parallel::getId()<<") : q= "<<block<<", p = "<<nBlocks<<" ["<<i0<<","<<i1<<"]."<<endr;


    myGrid.clear();

#ifdef DEBUG_MULTIGRID_TIMING
    myTimeAnterpolate.start();
#endif

    for(unsigned int i=i0;i<i1;i++)
      myGrid.anterpolateCharge(realTopo->atoms[i].scaledCharge,(*positions)[i],i);

#ifdef DEBUG_MULTIGRID_TIMING
    myTimeAnterpolate.stop();
#endif


    myGrid.getQ(0,begin,end);
    TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
    Parallel::reduceSlaves(begin,end);
    TimerStatistic::timer[TimerStatistic::FORCES] -= TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();

    //report << allnodes <<debug << "("<<Parallel::getId()<<") : 2 Q[10][10][10]= "<<begin[707]<<endr;
    //report << allnodes <<debug << "("<<Parallel::getId()<<") : "<<myGrid.sumQ(0)<<endr;
    //if(q != 0)
    //  return;
    //p = 1;
    //q = 0;
    //i0 = (q*atomCount)/p;
    //i1 = ((q+1)*atomCount)/p;


    //myGrid.printQ(0);
    for(int i=0;i<myLevels-1;i++){
#ifdef DEBUG_MULTIGRID_TIMING
      myTimeFineToCoarse[i].start();
#endif
      myGrid.fineToCoarse(i,block,nBlocks);
#ifdef DEBUG_MULTIGRID_TIMING
      myTimeFineToCoarse[i].stop();
#endif
      myGrid.getQ(i+1,begin,end);
      TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
      Parallel::reduceSlaves(begin,end);
      TimerStatistic::timer[TimerStatistic::FORCES] -= TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
#ifdef DEBUG_MULTIGRID_TIMING
      myTimeCorrectionMG[i].start();
#endif
      //myGrid.printQ(i+1);
      myGrid.correction(i,block,nBlocks);
      //myGrid.printV(i);
#ifdef DEBUG_MULTIGRID_TIMING
      myTimeCorrectionMG[i].stop();
#endif
    }

#ifdef DEBUG_MULTIGRID_TIMING
    myTimeDirect.start();
#endif
    myGrid.direct(block,nBlocks);
#ifdef DEBUG_MULTIGRID_TIMING
    myTimeDirect.stop();
#endif
    myGrid.getV(myLevels-1,begin,end);
    TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
    Parallel::reduceSlaves(begin,end);
    TimerStatistic::timer[TimerStatistic::FORCES] -= TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
    //myGrid.printV(myLevels-1);
    //Real e1 = myGrid.energy(myLevels-1);
    //report << plain << "Multigrid : level "<<myLevels-1<<": e= "<<e1<<endr;

    for(int i=myLevels-1; i>0;i--){
#ifdef DEBUG_MULTIGRID_TIMING
      myTimeCoarseToFine[i].start();
#endif
      myGrid.coarseToFine(i,block,nBlocks);
#ifdef DEBUG_MULTIGRID_TIMING
      myTimeCoarseToFine[i].stop();
#endif
      myGrid.getV(i-1,begin,end);
      TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
      Parallel::reduceSlaves(begin,end);
      TimerStatistic::timer[TimerStatistic::FORCES] -= TimerStatistic::timer[TimerStatistic::COMMUNICATION].lap();
      //myGrid.printV(i-1);
      //Real e1 = myGrid.energy(i-1);
      //report << plain << "Multigrid : level "<<i-1<<": e= "<<e1<<endr;
    }

#ifdef DEBUG_MULTIGRID_TIMING
    myTimeInterpolate.start();
#endif
    energy += myGrid.energy(0,block,nBlocks);

    bool doMolVirial = energies->molecularVirial();
    for(unsigned int i=i0;i<i1;i++){
      // compute the force on atom i from the reciprocal space part
      Vector3D fi;
      myGrid.interpolateForce(realTopo->atoms[i].scaledCharge,i,fi);
      (*forces)[i] += fi;

      if(doMolVirial){
	// compute the vector from atom i to the center of mass of the molecule
	Vector3D ria(realTopo->boundaryConditions.minimalPosition((*positions)[i]));
	Vector3D mri(realTopo->boundaryConditions.minimalDifference(ria,realTopo->molecules[realTopo->atoms[i].molecule].position));

	// compute the reciprocal space contribution to the molecular virial
	// this expression is taken from Darden, et al. J. Chem. Phys. 103 (19), 8577.
	energies->addMolVirial(fi,-mri);
      }
    }
#ifdef DEBUG_MULTIGRID_TIMING
    myTimeInterpolate.stop();
#endif
  }

  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  void NonbondedMultiGridSystemForce<TBoundaryConditions,
				     TCellManager,
				     TInterpolation,
				     TKernel,
				     direct, 
				     correction,
				     smooth>::shortRangeTerm(const RealTopologyType* realTopo, 
							     const Vector3DBlock* positions,
							     Vector3DBlock* forces, 
							     Real& shortRange,
							     ScalarStructure* energies, 
							     unsigned int n){

    CellPairType thisPair;
    unsigned int count = 0;
    bool doVirial = energies->virial();
    for (; !enumerator.done(); enumerator.next()) {
      enumerator.get(thisPair);
      bool notSameCell = enumerator.notSameCell();
    
      if(!notSameCell){
	count++;
	if(count > n)
	  break;
      }
      for(int i=thisPair.first; i!=-1; i=realTopo->atoms[i].cellListNext){
	Real qi = realTopo->atoms[i].scaledCharge;
	Vector3D ri((*positions)[i]);
	Vector3D fi(0.0,0.0,0.0);
	int mi = realTopo->atoms[i].molecule;
	for(int j=(notSameCell ? thisPair.second:realTopo->atoms[i].cellListNext); j!=-1; j=realTopo->atoms[j].cellListNext){
	  Real rSquared;
	  Vector3D rij(realTopo->boundaryConditions.minimalDifference(ri,(*positions)[j],rSquared));
	  if(rSquared > myS2)
	    continue;
	  int mj = realTopo->atoms[j].molecule;
	  bool same = (mi==mj);
	  ExclusionClass excl = (same?realTopo->exclusions.check(i,j):EXCLUSION_NONE);
	  if (excl == EXCLUSION_FULL)
	    continue;
#ifdef DEBUG_MULTIGRID_TIMING
	  myCounterDirect++;
#endif
	  Real qq  = qi*realTopo->atoms[j].scaledCharge;
	  Real r   = sqrt(rSquared);
	  Real rr = 1.0/r;
	  if (realTopo->coulombScalingFactor != 1.0 && excl == EXCLUSION_MODIFIED) 
	    qq *= realTopo->coulombScalingFactor;
	  Real force = (TKernel::dKernelR(rr) - TKernel::dSmooth(r,myS,myRS))*qq;
	  Vector3D fij(rij*rr*force);
	  fi           += fij;
	  (*forces)[j] -= fij;
	  if(doVirial)
	    energies->addVirial(fij,rij);
	  shortRange += (TKernel::kernelR(rr)-TKernel::smooth(r,myS,myRS))*qq; 
	}
	(*forces)[i] += fi;
      }
    }
  }

  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  void NonbondedMultiGridSystemForce<TBoundaryConditions,
				     TCellManager,
				     TInterpolation,
				     TKernel,
				     direct, 
				     correction,
				     smooth>::correctionTerm(const RealTopologyType* realTopo,
							     const Vector3DBlock* positions, 
							     Vector3DBlock* forces, 
							     Real& intraMolecularEnergy,
							     ScalarStructure* energies,
							     unsigned int from, unsigned int to){
    // Intra-molecular term
    const std::vector<ExclusionPair>& exclusions = realTopo->exclusions.getTable();
    bool doVirial = energies->virial();
    for(unsigned int i=from;i<to;i++){
      ExclusionPair excl = exclusions[i];
      Real rSquared;
      Vector3D rij(realTopo->boundaryConditions.minimalDifference((*positions)[excl.a1],(*positions)[excl.a2],rSquared));
      Real qq = realTopo->atoms[excl.a1].scaledCharge*realTopo->atoms[excl.a2].scaledCharge;
      if (excl.excl == EXCLUSION_MODIFIED)
	qq *= 1.0-realTopo->coulombScalingFactor;
      Real r     = sqrt(rSquared);
      intraMolecularEnergy -= TKernel::smoothKernel(r,myS,myRS)*qq;

      Real force = -TKernel::dSmoothKernel(r,myS,myRS)*qq;
      Vector3D fij(rij*force/r);
      (*forces)[excl.a1] += fij;
      (*forces)[excl.a2] -= fij;
      if(doVirial)
	energies->addVirial(fij,rij);
#ifdef DEBUG_MULTIGRID_TIMING
      myCounterCorrection++;
#endif
    }
  }


  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  unsigned int NonbondedMultiGridSystemForce<TBoundaryConditions,
					     TCellManager,
					     TInterpolation,
					     TKernel,
					     direct, 
					     correction,
					     smooth>::numberOfBlocks(const GenericTopology* topo, const Vector3DBlock* positions){
    const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);

    realTopo->updateCellLists(positions);

    unsigned int n = 0;
    
    // Direct part
    if(direct){
      n += Parallel::getNumberOfPackages(realTopo->cellLists.size());
    }
    // Correction term
    if(correction){
      n += std::min(static_cast<int>(realTopo->exclusions.getTable().size()),static_cast<int>(Parallel::getAvailableNum()));
    }
    // Smooth part
    if(smooth){
      //n++;
      // All slaves have to go through this ...
      ;
    }

    return n;
  }


  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  std::string NonbondedMultiGridSystemForce<TBoundaryConditions,
					    TCellManager,
					    TInterpolation,
					    TKernel,
					    direct, 
					    correction,
					    smooth>::getIdNoAlias() const{
    return (TKernel::getForceKeyword()+ 
	    " -algorithm " +getKeyword()+
	    " -interpolation "+TInterpolation::getKeyword()+
	    " -kernel "+TKernel::getKeyword() + 
	    (direct && (!smooth || !correction) ? std::string(" -direct")    : std::string("") ) +
	    (smooth && (!direct || !correction) ? std::string(" -smooth")    : std::string("") ) +
	    (correction && (!direct || !smooth) ? std::string(" -correction"): std::string("") ));
  }

  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  void NonbondedMultiGridSystemForce<TBoundaryConditions,
				     TCellManager,
				     TInterpolation,
				     TKernel,
				     direct, 
				     correction,
				     smooth>::getParameters(std::vector<Parameter>& parameters) const{
    parameters.push_back(Parameter("-s",Value(myS,ConstraintValueType::Positive()),Text("smoothing distance")));
    if(smooth){
      if(TBoundaryConditions::PERIODIC){
	parameters.push_back(Parameter("-toplevelgrid",Value(myNX,ConstraintValueType::Positive())));
	parameters.push_back(Parameter("",Value(myNY,ConstraintValueType::Positive())));
	parameters.push_back(Parameter("",Value(myNZ,ConstraintValueType::Positive())));
      }
      if(TBoundaryConditions::VACUUM){
	parameters.push_back(Parameter("-h",Value(myH,ConstraintValueType::NotZero())));
	parameters.push_back(Parameter("-origin",Value(myOrigin),Vector3D(0.0,0.0,0.0)));
      }
      parameters.push_back(Parameter("-levels",Value(myLevels,ConstraintValueType::Positive())));
      parameters.push_back(Parameter("-order",Value(myInterOrder,ConstraintValueType::Positive()),4,Text("interpolation")));
      parameters.push_back(Parameter("-ratio",Value(myRatio,ConstraintValueType::Positive()),2));
    }
  }

  template<class TBoundaryConditions, 
	   class TCellManager,
	   class TInterpolation,
	   class TKernel, 
	   bool direct, 
	   bool correction, 
	   bool smooth>
  Force* NonbondedMultiGridSystemForce<TBoundaryConditions,
				       TCellManager,
				       TInterpolation,
				       TKernel,
				       direct, 
				       correction,
				       smooth>::doMake(const std::vector<Value>& values) const{
    int i = 0;
    std::string err = "";
    Real s = values[i++];
    int nx = 0;
    int ny = 0;
    int nz = 0;
    Vector3D h(0.0,0.0,0.0);
    Vector3D origin(0.0,0.0,0.0);
    int  levels = 0;
    int  order = 0;
    int  ratio = 0;
    if(smooth){
      if(TBoundaryConditions::PERIODIC){
	nx = values[i++];
	ny = values[i++];
	nz = values[i++];
      }
      if(TBoundaryConditions::VACUUM){
	h = values[i++];
	origin = values[i++];
      }

      levels = values[i++];
      order  = values[i++];
      ratio  = values[i++];

      if(order < 2 || (order % 2) != 0)
	err += " order (="+toString(order)+") > 1 and even.";
      if(ratio < 2 )
	err += " ratio (="+toString(ratio)+") > 1.";
      if(levels <= 0){
	err += " levels (="+toString(levels)+") > 0.";

	if(TBoundaryConditions::PERIODIC){
	  if(levels == 1 && (nx < order || ny < order || nz < order))
	    err += " 1 level, "+toString(order)+" <= nx (="+toString(nx)+"), "+toString(order)+" <= ny (="+toString(ny)+"), "+toString(order)+" <= nz (="+toString(nz)+").";
	}
	if(TBoundaryConditions::VACUUM){
	  if(h.c[0] <= 0.0 || h.c[1] <= 0.0 || h.c[2] <= 0.0)
	    err += " hx (="+toString(h.c[0])+") > 0, hy (="+toString(h.c[1])+")) > 0, hz (="+toString(h.c[2])+") > 0..";
	}

      }
    }

    if(s <= 0.0)
      err += " s (="+toString(s)+") > 0.";
    
    if(!err.empty()){
      report << error << " force "+getKeyword()+" :"+err << endr;
      return NULL;
    }
    
    return (new NonbondedMultiGridSystemForce(nx,ny,nz,levels,s,order,ratio, h, origin));
  }
}
#endif /* NONBONDEDMULTIGRIDSYSTEMFORCE_H */
