/* -*- c++ -*- */
#ifndef PLANEFORCE_H
#define PLANEFORCE_H

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/parallel/Parallel.h>

#include <string>
#include <cmath>

#include <protomol/base/Report.h>
using namespace ProtoMol::Report;

namespace ProtoMol {
  //____ PlaneForce

  template<class TBoundaryConditions>
  class PlaneForce : public ExtendedForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    PlaneForce() : mySigma(0.0), myEpsilon(0.0),
                    myNormX(0.0), myNormY(0.0), myNormZ(0.0),
                    myPointX(0.0), myPointY(0.0), myPointZ(0.0) {}
    PlaneForce( double sig, double eps, double nx, double ny, double nz,
                    double px, double py, double pz) : mySigma(sig),
                    myEpsilon(eps), myNormX(nx), myNormY(ny), myNormZ(nz),
                    myPointX(px), myPointY(py), myPointZ(pz) {

        //norm
        double norm = sqrt( myNormX * myNormX + myNormY * myNormY + myNormZ * myNormZ);
        
        if( norm > 0.0 ){
            myNormX /= norm;
            myNormY /= norm;
            myNormZ /= norm;
        }
        
        //constant
        d =  -myNormX * myPointX - myNormY * myPointY - myNormZ * myPointZ;
    }
    virtual ~PlaneForce() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class PlaneForce
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
    virtual std::string getKeyword() const {return "Plane";}

    virtual unsigned int numberOfBlocks(const GenericTopology *topo,
                                        const Vector3DBlock *pos);

  private:
    virtual Force *doMake(const std::vector<Value> &) const;// {
    //  return new PlaneForce();
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
    double mySigma, myEpsilon, myNormX, myNormY, myNormZ,
            myPointX, myPointY, myPointZ, d;

  };

  //____ INLINES

  template<class TBoundaryConditions>
  inline void PlaneForce<TBoundaryConditions>::evaluate(
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
  inline void PlaneForce<TBoundaryConditions>::parallelEvaluate(
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
  inline void PlaneForce<TBoundaryConditions>::calcInteraction(
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

    //get MINIMAL difference vector, removes PBC extents
    Vector3D diff1 = boundary.minimalDifference(Vector3D(0.0,0.0,0.0), (*positions)[atomIndex] );

    // Distance squared, diff vector
    Real distSquared( diff1.c[0] * myNormX +
                        diff1.c[1] * myNormY +
                            diff1.c[2] * myNormZ + d );

    distSquared *= distSquared;

    Vector3D diff(myNormX, myNormY, myNormZ);

    //report << hint << "Distance " << distSquared << " y " << diff1.c[1] << " real y " << (*positions)[atomIndex].c[1] << endr;
    
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
        (*forces)[atomIndex] += fiv;

        // Add energy
        (*energies)[ScalarStructure::LENNARDJONES] += energy; //or OTHER?
    }
    
    // Add virial
    /*if (energies->virial())
      energies->addVirial(force1, r12);*/
  }

  template<class TBoundaryConditions>
   unsigned int PlaneForce<TBoundaryConditions>::numberOfBlocks(
    const GenericTopology *topo, const Vector3DBlock *positions) {
    return std::min(Parallel::getAvailableNum(),
                    static_cast<int>(positions->size()));
  }

  template<class TBoundaryConditions>
   inline Force* PlaneForce<TBoundaryConditions>::doMake(const std::vector<Value> & values) const {
        return new PlaneForce(values[0], values[1], values[2], values[3],
                                values[4], values[5], values[6], values[7]);

  }

  template<class TBoundaryConditions>
  inline void PlaneForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const {
      parameters.push_back(Parameter("-sigma", Value(mySigma)));
      //Parameter("-D", Value(D, ConstraintValueType::NotNegative()),
      //         80, Text("Bulk solvent dielectric"))
      parameters.push_back(Parameter("-epsilon", Value(myEpsilon)) );
      parameters.push_back(Parameter("-normx", Value(myNormX), 0.0) );
      parameters.push_back(Parameter("-normy", Value(myNormY), 1.0) );
      parameters.push_back(Parameter("-normz", Value(myNormZ), 0.0) );
      parameters.push_back(Parameter("-pointx", Value(myPointX), 0.0) );
      parameters.push_back(Parameter("-pointy", Value(myPointY), 0.0) );
      parameters.push_back(Parameter("-pointz", Value(myPointZ), 0.0) );
  }

  template<class TBoundaryConditions>
  inline unsigned int PlaneForce<TBoundaryConditions>::getParameterSize() const {
      return 8;
  }

}

#endif /* PLANEFORCE_H */
