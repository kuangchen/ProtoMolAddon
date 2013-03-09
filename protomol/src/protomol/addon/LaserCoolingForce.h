#ifndef LASERCOOLING_FORCE_H
#define LASERCOOLING_FORCE_H

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>

#include <protomol/addon/LuaState.h>

using namespace ProtoMol;
using namespace ProtoMol::Constant;
using namespace ProtoMolAddon::Lua;

namespace ProtoMolAddon {

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
    virtual Force* doMake(const std::vector<Value> &values) const {
      return new LaserCoolingForce(values[0]);
    };


  private:
    LuaState _L;
    string def_filename;
    //double intensity_sat, intensity, detuning, natural_linewidth, wavenumber, t_off, t_on;

    double _Isat, _I, _delta, _gamma, _k, _t_off, _t_on;
    double _sat, _alpha;

    double sat_parameter, alpha;

  public:
    LaserCoolingForce(){}
    ~LaserCoolingForce(){}
    
    LaserCoolingForce(const string& filename):
      _L(filename)
    {
      _Isat = _L.get<double>("laser_cooling.Isat");
      _I = _L.get<double>("laser_cooling.I");
      _delta = _L.get<double>("laser_cooling.delta");
      _gamma = _L.get<double>("laser_cooling.gamma");
      _k = _L.get<double>("laser_cooling.k");
      _t_off = _L.get<double>("laser_cooling.t_off");
      _t_on = _L.get<double>("laser_cooling.t_on");

      _sat = _I / _Isat;

      double hbar = 1.05457148e-34;
      double reduced_delta = _delta / _gamma;
      double a = (1 + _sat + 4* reduced_delta * reduced_delta);
      double conversion = 1/SI::AMU / SI::TIME_FS * TIMEFACTOR;

      alpha = 8* hbar* _k * _k * 1e4 * _delta * _sat / a / a; // See T.Loftus Thesis P25
      alpha *= conversion;

      // std::cout << "opening " << filename << "\n";
      // lua_State *L = luaL_newstate();
      // luaL_openlibs(L);
    
      // if (luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0))
      // 	luaL_error(L, "cannot read file: %s", lua_tostring(L, -1));

      // ExtractDouble(L, &intensity_sat, "intensity_sat", "intensity_sat should be a number\n");
      // ExtractDouble(L, &intensity, "intensity", "intensity should be a number\n");
      // ExtractDouble(L, &detuning, "detuning", "intensity should be a number\n");
      // ExtractDouble(L, &natural_linewidth, "natural_linewidth"," natural_linewidth should be a number\n");
      // ExtractDouble(L, &wavenumber, "wavenumber", "wavenumber should be a number\n");
      // ExtractDouble(L, &t_on, "t_on", "t_on should be a number\n");
      // ExtractDouble(L, &t_off, "t_off", "t_off should be a number\n");
      // lua_close(L);
    
      // sat_parameter = intensity / intensity_sat;

      // double hbar = 1.05457148e-34;
      // double reduced_detuning = detuning / natural_linewidth;
      // double a = (1 + sat_parameter + 4* reduced_detuning * reduced_detuning);
      // double conversion = 1/Constant::SI::AMU / Constant::SI::TIME_FS * Constant::TIMEFACTOR;
      // alpha = 8* hbar* wavenumber * wavenumber * 1e4 * reduced_detuning * sat_parameter / a / a; // See T.Loftus Thesis P25
      // alpha *= conversion;
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
    if ( topo->time/SI::TIME_FS < _t_off && topo->time/SI::TIME_FS > _t_on )
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
