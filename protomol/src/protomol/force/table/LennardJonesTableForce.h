/* -*- c++ -*- */
#ifndef LENNARDJONESTABLEFORCE_H
#define LENNARDJONESTABLEFORCE_H

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
  class LennardJonesLookupValues {
  public:
    enum {ENTRIES = 2};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LennardJonesLookupValues
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    template<typename TReal>
    void assign(Real r, Real r1, int ex, Real v, Real d, TReal *val) const {
      Real e = (ex == -12 ? r : power < -12 > (r1)); // value
      Real f = 12.0 * (ex == -14 ? r : power < -14 > (r1)); // gradient
      val[0] = e * v;
      val[1] = -0.5 * (f * v - e * d);
      e = -(ex == -6 ? r : power < -6 > (r1)); // value
      f = -6.0 * (ex == -8 ? r : power < -8 > (r1)); // gradient
      val[4] = e * v;
      val[5] = -0.5 * (f * v - e * d);
    }
  };

  //____ LennardJonesTableForce
  template<class TSwitchingFunction, unsigned int PRE, typename TReal = Real>
  class LennardJonesTableForce :
    public LookupTable<LennardJonesLookupValues, PRE, TReal> {

    typedef LookupTable<LennardJonesLookupValues, PRE, TReal> LookupTable_t;
    
  public:
    enum {DIST_R2 = 0};
    enum {CUTOFF = 1};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    LennardJonesTableForce() : myCutoff(0.0), myCutoff2(0.0) {}

    LennardJonesTableForce(TSwitchingFunction swf)
      : LookupTable_t
        (0.1, swf.cutoffSquared(), 2, LennardJonesLookupValues(), swf, 128),
        myCutoff(swf.cutoff()), myCutoff2(swf.cutoffSquared()),
        switchingFunction(swf) {}

    LennardJonesTableForce(TSwitchingFunction swf, Real rc) :
      LookupTable_t
      (0.1, square(rc), 2,  LennardJonesLookupValues(), swf, 128),
      myCutoff(rc), myCutoff2(rc * rc), switchingFunction(swf) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LennardJonesTableForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force,
                    Real distSquared, Real /*rDistSquared*/,
                    const Vector3D & /*diff*/, const GenericTopology *topo,
                    int atom1, int atom2, ExclusionClass excl) const {
      const LennardJonesParameters &params =
        topo->lennardJonesParameters(topo->atoms[atom1].type,
                                     topo->atoms[atom2].type);

      Real A = (excl != EXCLUSION_MODIFIED ? params.A : params.A14);
      Real B = (excl != EXCLUSION_MODIFIED ? params.B : params.B14);

      Real dt;
      int i;
      this->index(distSquared, i, dt);

      Real a = this->myTable[i + 0] * A + this->myTable[i + 4] * B;
      Real b = this->myTable[i + 1] * A + this->myTable[i + 5] * B;
      Real c = this->myTable[i + 2] * A + this->myTable[i + 6] * B;
      Real d = this->myTable[i + 3] * A + this->myTable[i + 7] * B;

      this->interpolate(a, b, c, d, dt, energy, force);
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      (*energies)[ScalarStructure::LENNARDJONES] += energy;
    }

    static Real getEnergy(const ScalarStructure *energies) {
      return (*energies)[ScalarStructure::LENNARDJONES];
    }

    static std::string getKeyword() {return "LennardJonesTable";}

    // Parsing
    static std::string getId() {
      return getKeyword() +
        std::string(!TSwitchingFunction::USE ? "" :
                    std::string(" -switchingFunction " +
                                TSwitchingFunction::getId()));
    }

    void getParameters(std::vector<Parameter> &parameters) const {
      switchingFunction.getParameters(parameters);
      if (!TSwitchingFunction::CUTOFF)
        parameters.push_back
          (Parameter("-cutoff",
                     Value(myCutoff, ConstraintValueType::Positive()),
                     Text("cutoff for table look up")));
    }

    static unsigned int getParameterSize() {
      return TSwitchingFunction::getParameterSize() +
        (TSwitchingFunction::CUTOFF ? 0 : 1);
    }

    static LennardJonesTableForce make(const std::vector<Value> &values) {
      if (!TSwitchingFunction::CUTOFF) {
        std::vector<Value> sValues(values.begin(), values.end() - 1);
        return LennardJonesTableForce(TSwitchingFunction::make(sValues),
                                      values[values.size() - 1]);

      } else return LennardJonesTableForce(TSwitchingFunction::make(values));
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
#endif /* LENNARDJONESTABLEFORCE_H */
