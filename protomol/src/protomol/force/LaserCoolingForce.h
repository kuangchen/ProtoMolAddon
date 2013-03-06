#ifndef LASERCOOLING_FORCE_H
#define LASERCOOLING_FORCE_H

extern "C" {
#include <lua5.1/lua.h>
#include <lua5.1/lauxlib.h>
#include <lua5.1/lualib.h>
}

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
      return new LaserCoolingForce(values[0]);
    };


  private:
    string def_filename;
    double intensity_sat, intensity, detuning, natural_linewidth, wavenumber, t_off, t_on;
    double sat_parameter, alpha;

  public:
    LaserCoolingForce(){}
    ~LaserCoolingForce(){}
    
    LaserCoolingForce(const string& filename)
    {
      std::cout << "opening " << filename << "\n";
      lua_State *L = luaL_newstate();
      luaL_openlibs(L);
    
      if (luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0))
	luaL_error(L, "cannot read file: %s", lua_tostring(L, -1));

      ExtractDouble(L, &intensity_sat, "intensity_sat", "intensity_sat should be a number\n");
      ExtractDouble(L, &intensity, "intensity", "intensity should be a number\n");
      ExtractDouble(L, &detuning, "detuning", "intensity should be a number\n");
      ExtractDouble(L, &natural_linewidth, "natural_linewidth"," natural_linewidth should be a number\n");
      ExtractDouble(L, &wavenumber, "wavenumber", "wavenumber should be a number\n");
      ExtractDouble(L, &t_on, "t_on", "t_on should be a number\n");
      ExtractDouble(L, &t_off, "t_off", "t_off should be a number\n");
      lua_close(L);
    

      sat_parameter = intensity / intensity_sat;

      double hbar = 1.05457148e-34;
      double reduced_detuning = detuning / natural_linewidth;
      double a = (1 + sat_parameter + 4* reduced_detuning * reduced_detuning);
      double conversion = 1/Constant::SI::AMU / Constant::SI::TIME_FS * Constant::TIMEFACTOR;
      alpha = 8* hbar* wavenumber * wavenumber * 1e4 * reduced_detuning * sat_parameter / a / a; // See T.Loftus Thesis P25
      alpha *= conversion;
    }
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
    if ( topo->time/Constant::SI::TIME_FS < t_off && topo->time/Constant::SI::TIME_FS > t_on )
#pragma omp parallel for schedule(dynamic, 3)
      for ( unsigned int i=0 ; i<topo->atoms.size() ; i++ )
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
    parameters.push_back(Parameter("-laser_cooling_def_filename", Value(def_filename)));
  }
}
#endif
