/* -*- c++ -*- */
#ifndef COULOMBTABLEFORCE_H
#define COULOMBTABLEFORCE_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/config/Parameter.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/force/table/LookupTable.h>
#include <string>

namespace ProtoMol {
  /**
   * Defines the function(s) to be tabulated.
   */
  class CoulombLookupValues {
  public:
    enum {ENTRIES = 1};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombLookupValues
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    template<typename TReal>
    void assign(Real r, Real r1, int ex, Real v, Real d, TReal *val) const {
      Real e = (ex == -1 ? r : power < -1 > (r1)); // value
      Real f = (ex == -3 ? r : power < -3 > (r1)); // gradient
      val[0] = e * v;
      val[1] = -0.5 * (f * v - e * d);
    }
  };

  //____ CoulombTableForce
  template<class TSwitchingFunction, unsigned int PRE, typename TReal = Real>
  class CoulombTableForce :
    public LookupTable<CoulombLookupValues, PRE, TReal> {
    
    typedef LookupTable<CoulombLookupValues, PRE, TReal> LookupTable_t;

  public:
    /// no need for reciprocal squared distance
    enum {DIST_R2 = 0};
    /// has a valid cutoff
    enum {CUTOFF = 1};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CoulombTableForce() : myCutoff(0.0), myCutoff2(0.0) {}

    /// Constructor with a cutoff switching function
    CoulombTableForce(TSwitchingFunction swf) :
      LookupTable_t
      (0.1, swf.cutoffSquared(), 2, CoulombLookupValues(), swf, 128),
      myCutoff(swf.cutoff()), myCutoff2(swf.cutoffSquared()),
      switchingFunction(swf) {}

    /// Constructor without a cutoff switching function
    CoulombTableForce(TSwitchingFunction swf, Real rc) :
      LookupTable_t(0.1, square(rc), 2, CoulombLookupValues(), swf, 128),
      myCutoff(rc), myCutoff2(rc * rc), switchingFunction(swf) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombTableForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force, Real distSquared,
                    Real /*rDistSquared*/, const Vector3D & /*diff*/,
                    const GenericTopology *topo, int atom1, int atom2,
                    ExclusionClass excl) const {
      Real q = topo->atoms[atom1].scaledCharge *
               topo->atoms[atom2].scaledCharge *
               ((topo->coulombScalingFactor != 1.0 && excl ==
                 EXCLUSION_MODIFIED) ? topo->coulombScalingFactor : 1.0);

      Real dt;
      int i;
      // get index and interpolation delta dt
      this->index(distSquared, i, dt);

      // compute the coefficients
      Real a = this->myTable[i + 0] * q;
      Real b = this->myTable[i + 1] * q;
      Real c = this->myTable[i + 2] * q;
      Real d = this->myTable[i + 3] * q;

      // interpolate energy and force
      this->interpolate(a, b, c, d, dt, energy, force);
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      (*energies)[ScalarStructure::COULOMB] += energy;
    }

    static Real getEnergy(const ScalarStructure *energies) {
      return (*energies)[ScalarStructure::COULOMB];
    }

    static std::string getKeyword() {return "CoulombTable";}
    
    // Parsing
    static std::string getId() {
      return getKeyword() +
        std::string((!TSwitchingFunction::USE) ? "" :
                    std::string(" -switchingFunction " +
                                TSwitchingFunction:: getId()));
    }

    static unsigned int getParameterSize() {
      return (TSwitchingFunction::CUTOFF ? 0 : 1) +
             TSwitchingFunction::getParameterSize();
    }

    void getParameters(std::vector<Parameter> &parameters) const {
      switchingFunction.getParameters(parameters);
      if (!TSwitchingFunction::CUTOFF)
        parameters.push_back
          (Parameter("-cutoff",
                     Value(myCutoff, ConstraintValueType::Positive()),
                     Text("cutoff for table look up")));
    }

    static CoulombTableForce make(const std::vector<Value> &values) {
      if (!TSwitchingFunction::CUTOFF) {
        std::vector<Value> sValues(values.begin(), values.end() - 1);
        return CoulombTableForce(TSwitchingFunction::make(sValues),
                                 values[values.size() - 1]);

      } else return CoulombTableForce(TSwitchingFunction::make(values));
    }

    Real cutoffSquared() const {return myCutoff2;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myCutoff;
    Real myCutoff2;
    TSwitchingFunction switchingFunction;
  };
}
#endif /* COULOMBTABLEFORCE_H */
