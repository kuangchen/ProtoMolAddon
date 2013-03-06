/* -*- c++ -*- */
#ifndef VESSELFORCE_H
#define VESSELFORCE_H

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/parallel/Parallel.h>

#include <string>

#include <protomol/base/Report.h>
using namespace ProtoMol::Report;

namespace ProtoMol {
  //____ VesselForce

  template<class TBoundaryConditions>
  class VesselForce : public ExtendedForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    VesselForce() : mySigma(0.0), myEpsilon(0.0), 
                    myDiameter(0.0), myGranularity(0.0) {}
    VesselForce( double sig, double eps, double dia, double gran ) : mySigma(sig),
        myEpsilon(eps), myDiameter(dia), myGranularity(gran) {}
    virtual ~VesselForce() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class VesselForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void calcInteraction(const GenericTopology *topo,
                  const TBoundaryConditions &boundary, const int atomIndex,
                  const Vector3DBlock *positions, const Vector3DBlock *velocities,
                  Vector3DBlock *forces,
                  ScalarStructure *energies);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void evaluate(const GenericTopology *topo,
                          const Vector3DBlock *positions,
                          const Vector3DBlock *velocities,
                          Vector3DBlock *forces,
                          ScalarStructure *energies);

    virtual void parallelEvaluate(const GenericTopology *topo,
                                  const Vector3DBlock *positions,
                                  const Vector3DBlock *velocities,
                                  Vector3DBlock *forces,
                                  ScalarStructure *energies);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getKeyword() const {return "Vessel";}

    virtual unsigned int numberOfBlocks(const GenericTopology *topo,
                                        const Vector3DBlock *pos);

  private:
    virtual Force *doMake(const std::vector<Value> &) const;// {
    //  return new VesselForce();
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
    double mySigma, myEpsilon, myDiameter, myGranularity;

  };

  //____ INLINES

  template<class TBoundaryConditions>
  inline void VesselForce<TBoundaryConditions>::evaluate(
    const GenericTopology *topo, const Vector3DBlock *positions,
    const Vector3DBlock *velocities,
    Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      ((SemiGenericTopology<TBoundaryConditions> &)(*topo)).
        boundaryConditions;

    const unsigned int count = positions->size();

    for (unsigned int i = 0; i < count; i++){
        calcInteraction(topo, boundary, i, 
                positions, velocities, forces, energies);
    }
    
  }


  template<class TBoundaryConditions>
  inline void VesselForce<TBoundaryConditions>::parallelEvaluate(
    const GenericTopology *topo, const Vector3DBlock *positions,
    const Vector3DBlock *velocities,
    Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions> &>(*topo)).
        boundaryConditions;

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
                    positions, velocities, forces, energies);
        }
      }

  }

  template<class TBoundaryConditions>
  inline void VesselForce<TBoundaryConditions>::calcInteraction(
    const GenericTopology *topo,
    const TBoundaryConditions &boundary, const int atomIndex,
    const Vector3DBlock *positions, const Vector3DBlock *velocities,
    Vector3DBlock *forces,
    ScalarStructure *energies) {

    //get atom type
    const int type = topo->atoms[atomIndex].type;
    
    //Get parameters
    Real sigma = topo->atomTypes[type].sigma;
    Real epsilon = topo->atomTypes[type].epsilon;

    //calculate factors
    Real r_ij = sigma + mySigma;
    Real e_ij = sqrt(epsilon * myEpsilon);
    Real A = power<12>(r_ij) * e_ij;
    Real B = 2 * power<6>(r_ij) * e_ij;

    //report << hint << "i " << atomIndex << " name " << topo->atoms[atomIndex].name <<  " A " << A << " B " << B << " sigma " << sigma << " epsilon "<< epsilon << endr;
    
    //Calculate SCE-Vessel force
    // Distance squared, diff vector
    Real distSquared(0.0);

    Vector3D diff(0.0,0.0,0.0);

    //get MINIMAL difference vector, removes PBC extents
    Vector3D diff1 = boundary.minimalDifference((*positions)[atomIndex], Vector3D(0.0,0.0,0.0) );

    //get smallest Y distance to cylinder
    Real ydiff = fabs(diff1[1] - myDiameter) < fabs(diff1[1] + myDiameter) ?
                diff1[1] - myDiameter : diff1[1] + myDiameter;

    //then Z-component smallest distance
    Real zdiff = fabs(diff1[2] - myDiameter) < fabs(diff1[2] + myDiameter) ?
                diff1[2] - myDiameter : diff1[2] + myDiameter;

    //find closest, Y or Z and change the diff vector (direction of force) to represent it
    //diff zero here
    //
    if(fabs(ydiff) < fabs(zdiff)){ 
      diff[1] = ydiff;
    }else{
      diff[2] = zdiff;
    }

    //figure out granularity of vessel, to represent cells

    //find last "cell" and next cell on vessel
    double last_vessel_cell = diff1[0] - floor( diff1[0] / myGranularity ) * myGranularity;

    double next_vessel_cell = myGranularity - last_vessel_cell;

    //find direction of SCE, chose trailing restraint
    //x velocity > 0?
    if( (*velocities)[atomIndex][0] > 0.0 ) {
        diff[0] = last_vessel_cell;
    }else{
        diff[0] = -next_vessel_cell;
    }

    //get distance to vessel cell
    distSquared = diff.normSquared();

    //report << hint << "Distance " << distSquared << endr;
    
    //Test for non zero distance squared
    if( distSquared ){
        Real rDistSquared = 1.0 / distSquared;

        // Calculate force on atom due to vessel.
        Real r6 = rDistSquared * rDistSquared * rDistSquared;
        Real r12 = r6 * r6;
        Real r6B = B * r6;
        Real r12A = A * r12;
        Real energy = r12A - r6B;
        Real force = 12.0 * r12A * rDistSquared - 6.0 * r6B * rDistSquared;

        //report << hint << "Force " << force << endr;

        //add it
        Vector3D fiv(diff * force);

        // Add to the total force.
        (*forces)[atomIndex] -= fiv;

        // Add energy
        (*energies)[ScalarStructure::LENNARDJONES] += energy; //or OTHER?
    }
    
    // Add virial
    /*if (energies->virial())
      energies->addVirial(force1, r12);*/
  }

  template<class TBoundaryConditions>
   unsigned int VesselForce<TBoundaryConditions>::numberOfBlocks(
    const GenericTopology *topo, const Vector3DBlock *positions) {
    return std::min(Parallel::getAvailableNum(),
                    static_cast<int>(positions->size()));
  }

  template<class TBoundaryConditions>
   inline Force* VesselForce<TBoundaryConditions>::doMake(const std::vector<Value> & values) const {
        return new VesselForce(values[0], values[1], values[2], values[3]);

  }

  template<class TBoundaryConditions>
  inline void VesselForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const {
      parameters.push_back(Parameter("-sigma", Value(mySigma)));
      //Parameter("-D", Value(D, ConstraintValueType::NotNegative()),
      //         80, Text("Bulk solvent dielectric"))
      parameters.push_back(Parameter("-epsilon", Value(myEpsilon)));
      parameters.push_back(Parameter("-diameter", Value(myDiameter)));
      parameters.push_back(Parameter("-granularity", Value(myGranularity)));
  }

  template<class TBoundaryConditions>
  inline unsigned int VesselForce<TBoundaryConditions>::getParameterSize() const {
      return 4;
  }

}

#endif /* VESSELFORCE_H */
