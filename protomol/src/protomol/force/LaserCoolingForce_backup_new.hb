#ifndef LASERCOOLING_FORCE_H
#define LASERCOOLING_FORCE_H

extern "C" {
#include <lua5.1/lua.h>
#include <lua5.1/lauxlib.h>
#include <lua5.1/lualib.h>
}


#include <cmath>
#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>

#include <protomol/io/LuaConfigReader.h>
#include <omp.h>

namespace ProtoMol{

  using namespace Util;

  template<class TBoundaryConditions>
  class LaserCoolingForce: public ExtendedForce
  {
  public:
    LaserCoolingForce(): intensity_sat(1), intensity(1), detuning(1), natural_linewidth(1), wavenumber(1), t_off(1), sat_parameter(1)
    {
    };

    LaserCoolingForce(Real i_sat, Real i, Real d, Real lw, Real wn, Real t): 
      intensity_sat(i_sat), intensity(i), detuning(d), natural_linewidth(lw), wavenumber(wn), t_off(t), sat_parameter(i/i_sat)
	{
	  double hbar = 1.05457148e-34;
	  double reduced_detuning = detuning / natural_linewidth;
	  double a = (1 + sat_parameter + 4* reduced_detuning * reduced_detuning);
	  double conversion = 1/Constant::SI::AMU / Constant::SI::TIME_FS * Constant::TIMEFACTOR;
	  alpha = 8* hbar* wavenumber * wavenumber * 1e4 * reduced_detuning * sat_parameter / a / a; // See T.Loftus Thesis P25
	  alpha *= conversion;
	};

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

    virtual void getParameters(std::vector<Parameter>& ) const; 
    static const std::string keyword;
    virtual std::string getKeyword() const{return keyword;}
    virtual std::string getIdNoAlias() const{return keyword;}





  private:
    virtual Force* doMake(const std::vector<Value> &values) const
    {
      return new LaserCoolingForce(values[0], values[1], values[2], values[3], values[4], values[5]);
    };

    Real intensity_sat, intensity, detuning, natural_linewidth, wavenumber, t_off;
    
    Real sat_parameter, alpha;
  };
  
  template <class TBoundaryConditions>
  const std::string LaserCoolingForce<TBoundaryConditions>::keyword("LaserCooling");  

  template<class TBoundaryConditions>
  inline void LaserCoolingForce<TBoundaryConditions>::evaluate(const GenericTopology *topo,
							       const Vector3DBlock *positions,
							       const Vector3DBlock *velocities,
							       Vector3DBlock *forces,
							       ScalarStructure *energies)
  {
    if(topo->time/Constant::SI::TIME_FS < t_off)
      
#pragma omp parallel for schedule(dynamic, 3)
      for(unsigned int i=0;i<topo->atoms.size();i++)
	(*forces)[i] -= (*velocities)[i]*alpha;
  }
  

  template<class TBoundaryConditions>
  inline void LaserCoolingForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology *topo,
								       const Vector3DBlock *positions,
								       const Vector3DBlock *velocities,
								       Vector3DBlock *forces,
								       ScalarStructure *energies)
  {
    evaluate(topo, positions, velocities, forces, energies); // Not implemented right now
  }

  template<class TBoundaryConditions>
  inline void LaserCoolingForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const
  {
    parameters.push_back(Parameter("-intensity_sat", Value(intensity_sat)));
    parameters.push_back(Parameter("-intensity", Value(intensity)));
    parameters.push_back(Parameter("-detuning", Value(detuning)));
    parameters.push_back(Parameter("-natural_linewidth", Value(natural_linewidth)));
    parameters.push_back(Parameter("-wavenumber", Value(wavenumber)));
    parameters.push_back(Parameter("-t_off", Value(t_off)));
  }
}
#endif
