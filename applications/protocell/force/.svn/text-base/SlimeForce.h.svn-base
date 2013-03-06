/* -*- c++ -*- */
#ifndef SLIMEFORCE_H
#define SLIMEFORCE_H

#include <protomol/force/system/SystemForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/parallel/Parallel.h>

#include <string>
#include <map>
#include <iostream>     // required for cout

#include <protomol/base/Report.h>
using namespace ProtoMol::Report;

namespace ProtoMol {
  //____ SlimeForce

  template<class TBoundaryConditions>
  class SlimeForce : public SystemForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    SlimeForce() : myForce(0.0), myReverseSteps(0), 
            myDwellSteps(0), setInitialData(true) {}
    SlimeForce( double frc, int rstep, int dstep ) : myForce(frc),
        myReverseSteps(rstep), myDwellSteps(dstep), 
        setInitialData(true) {
    }
    virtual ~SlimeForce() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class SlimeForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void calcInteraction(const GenericTopology *topo,
                  const TBoundaryConditions &boundary, const int atomIndex,
                  const Vector3DBlock *positions,
                  Vector3DBlock *forces,
                  ScalarStructure *energies);

    void initializeData( const GenericTopology *topo );
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void evaluate(const GenericTopology *topo,
                          const Vector3DBlock *positions,
                          Vector3DBlock *forces,
                          ScalarStructure *energies);

    virtual void parallelEvaluate(const GenericTopology *topo,
                                  const Vector3DBlock *positions,
                                  Vector3DBlock *forces,
                                  ScalarStructure *energies);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getKeyword() const {return "Slime";}

    virtual unsigned int numberOfBlocks(const GenericTopology *topo,
                                        const Vector3DBlock *pos);

  private:
    virtual Force *doMake(const std::vector<Value> &) const;// {
    //  return new SlimeForce();
    //}
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return getKeyword();}
    virtual void getParameters(std::vector<Parameter> &) const;// {}
    virtual unsigned int getParameterSize() const;//{return 2;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    double myForce;
    int myReverseSteps;
    int myDwellSteps;
    
    bool setInitialData;

    //per bacteria motion control
    struct motion{
        bool forward, dwell;
        int currentStep;
    };

    //bacteria data storage
    map< int, motion > bacteria;
    vector< int > bacteriaLookup;

  };

  //____ INLINES

  template<class TBoundaryConditions>
  inline void SlimeForce<TBoundaryConditions>::evaluate(
    const GenericTopology *topo, const Vector3DBlock *positions,
    Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      ((SemiGenericTopology<TBoundaryConditions> &)(*topo)).
        boundaryConditions;

    //initialize?
    if( setInitialData ){
        initializeData( topo );
        setInitialData = false;
    }

    //do updates
    if( myReverseSteps ){

        const int sz = bacteriaLookup.size();

        for( int i=0; i<sz; i++){
            
            const int bindx = bacteriaLookup[i];
            
            if(++bacteria[bindx].currentStep > myReverseSteps){ //increment reverse counter

                if(bacteria[bindx].currentStep > ( myReverseSteps + myDwellSteps )){
                    bacteria[bindx].dwell = false;
                    if(bacteria[bindx].forward) bacteria[bindx].forward = false;
                    else bacteria[bindx].forward = true;
                    bacteria[bindx].currentStep %= myReverseSteps + myDwellSteps;
                }else{
                    bacteria[bindx].dwell = true;
                }
            }            
        }
    }

    const unsigned int count = positions->size();

    for (unsigned int i = 0; i < count; i++){
        calcInteraction(topo, boundary, i, 
                positions, forces, energies);
    }
    
  }

  template<class TBoundaryConditions>
  inline void SlimeForce<TBoundaryConditions>::parallelEvaluate(
    const GenericTopology *topo, const Vector3DBlock *positions,
    Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions> &>(*topo)).
        boundaryConditions;

    //initialize?
    if( setInitialData ){
        initializeData( topo );
        setInitialData = false;
    }

    //do updates
    if( myReverseSteps ){

        const int sz = bacteriaLookup.size();

        for( int i=0; i<sz; i++){

            const int bindx = bacteriaLookup[i];

            if(++bacteria[bindx].currentStep > myReverseSteps){ //increment reverse counter

                if(bacteria[bindx].currentStep > ( myReverseSteps + myDwellSteps )){
                    bacteria[bindx].dwell = false;
                    if(bacteria[bindx].forward) bacteria[bindx].forward = false;
                    else bacteria[bindx].forward = true;
                    bacteria[bindx].currentStep %= myReverseSteps + myDwellSteps;
                }else{
                    bacteria[bindx].dwell = true;
                }
            }
        }
    }
    
    unsigned int n = positions->size();
    unsigned int count = numberOfBlocks(topo, positions);

    for (unsigned int i = 0; i < count; i++)
      if (Parallel::next()) {
        int to = (n * (i + 1)) / count;
        if (to > static_cast<int>(n))
          to = n;
        int from = (n * i) / count;
        for (int j = from; j < to; j++){
            calcInteraction(topo, boundary, i, 
                    positions, forces, energies);
        }
      }

  }

  template<class TBoundaryConditions>
  inline void SlimeForce<TBoundaryConditions>::calcInteraction(
    const GenericTopology *topo,
    const TBoundaryConditions &boundary, const int atomIndex,
    const Vector3DBlock *positions,
    Vector3DBlock *forces,
    ScalarStructure *energies) {

    //by convention CE?, GE? etc. for end cells
    //E2 is the tail, so add slime force
    if( (topo->atoms[atomIndex].name.c_str())[1] == 'E' ){

      //find index of bacteria map
      const int resindx = topo->atoms[atomIndex].residue_seq;

      //dwell? then return
      if( bacteria[resindx].dwell ){
        return;
      }


      //Head or Tail motility force (slime)
      if( 
          ( bacteria[resindx].forward && (topo->atoms[atomIndex].name.c_str())[2] == '2' && atomIndex > 0 )
            || ( !bacteria[resindx].forward && (topo->atoms[atomIndex].name.c_str())[2] == '1' && atomIndex  < (int)positions->size() - 1 )
               ){

          //Assume previous node is connected, hence atomIndex>0
          //get MINIMAL difference vector, removes PBC extents if necesary
          Vector3D diff;

          if(bacteria[resindx].forward) //CE2?
                diff = boundary.minimalDifference( (*positions)[atomIndex],
                                                            (*positions)[atomIndex - 1] );
          else //then CE1
                diff = boundary.minimalDifference( (*positions)[atomIndex],
                                                            (*positions)[atomIndex + 1] );

          //find norm of vector
          const Real norm = diff.norm();

          if( norm ){
            // Add to the total force.
            (*forces)[atomIndex] += diff * myForce / norm;

            // Add energy
            //(*energies)[ScalarStructure::OTHER] += 0.0;
          }
      }

      //Head or Tail social force (pilli)
      if(
          ( !bacteria[resindx].forward && (topo->atoms[atomIndex].name.c_str())[2] == '2' && atomIndex > 0 )
            || ( bacteria[resindx].forward && (topo->atoms[atomIndex].name.c_str())[2] == '1' && atomIndex  < (int)positions->size() - 1 )
               ){

      }
    }
    
    // Add virial
    /*if (energies->virial())
      energies->addVirial(force1, r12);*/
  }

  template<class TBoundaryConditions>
   unsigned int SlimeForce<TBoundaryConditions>::numberOfBlocks(
    const GenericTopology *topo, const Vector3DBlock *positions) {
    return std::min(Parallel::getAvailableNum(),
                    static_cast<int>(positions->size()));
  }

  template<class TBoundaryConditions>
   inline Force* SlimeForce<TBoundaryConditions>::doMake(const std::vector<Value> & values) const {
        return new SlimeForce(values[0], values[1], values[2]);

  }

  template<class TBoundaryConditions>
  inline void SlimeForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const {
      parameters.push_back(Parameter("-force", Value(myForce)));
      parameters.push_back(Parameter("-reversesteps", Value(myReverseSteps), 0));
      parameters.push_back(Parameter("-dwellsteps", Value(myDwellSteps), 0));
  }

  template<class TBoundaryConditions>
  inline unsigned int SlimeForce<TBoundaryConditions>::getParameterSize() const {
      return 3;
  }

  template<class TBoundaryConditions>
  inline void SlimeForce<TBoundaryConditions>::initializeData( const GenericTopology *topo ) {

    const unsigned int csize = topo->atoms.size();

    //loop to find cells/bacteria
    for (unsigned int i=0; i < csize; i++) {

        const unsigned int sz = bacteria.size();

        const int resindx = topo->atoms[i].residue_seq;

        bacteria[resindx];

        //did map increase in size? f so new entry
        if( bacteria.size() != sz ){

            //create motion struct
            motion mot;

            mot.dwell = false;

            if( myReverseSteps ){
                mot.currentStep = (int)( randomNumber() * (myReverseSteps + myDwellSteps) );

                //initial direction
                const Real rand = randomNumber();

                if( rand > 0.5 )
                    mot.forward = true;
                else
                    mot.forward = false;
            }else{
                mot.currentStep = 0;
                mot.forward = true;
            }

            //store it
            bacteria[resindx] = mot;

            //update lookup
            bacteriaLookup.push_back( resindx );

        }
    }

  }

}

#endif /* SLIMEFORCE_H */
