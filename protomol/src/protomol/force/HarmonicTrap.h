#ifndef HARMONIC_TRAP_H
#define HARMONIC_TRAP_H

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>
#include <protomol/io/LuaConfigReader.h>
#include <vector>

#include <omp.h>

namespace ProtoMol{
  using namespace Util;

  template<class TBoundaryConditions>
  class HarmonicTrapForce: public ExtendedForce
  {
  public:
  HarmonicTrapForce(): ht_def_filename(""), omega(0)
    {
      omega.resize(0);
    }

    HarmonicTrapForce(const std::string& def_filename): 
    ht_def_filename(def_filename),
      omega(3)
    {
      std:: cout << "opening " << def_filename <<"\n";
      lua_State *L = luaL_newstate();
      luaL_openlibs(L);
    
      if (luaL_loadfile(L, def_filename.c_str()) || lua_pcall(L, 0, 0, 0))
	luaL_error(L, "cannot run configuration file: %s",
		   lua_tostring(L, -1));

      lua_getglobal(L, "omega");
      
      if (lua_type(L, -1) != LUA_TTABLE)
	luaL_error(L, "omega is not an array");

      int len = lua_objlen(L, -1);

      if ( len <3 )
	luaL_error(L, "omega should contain three elements");
      
      for (int i = 1; i <= 3; i++){
	lua_rawgeti(L, -1, i);
	omega[i-1] = lua_tonumber(L, -1);
	lua_pop(L, 1);
      }

      lua_close(L);

    }
    
    virtual void evaluate(const GenericTopology* topo,
			  const Vector3DBlock* positions,
                          const Vector3DBlock *velocities,
			  Vector3DBlock* forces,
			  ScalarStructure* energies);

    virtual void parallelEvaluate(const GenericTopology* topo,
				  const Vector3DBlock* positions,
				  const Vector3DBlock *velocities,
				  Vector3DBlock* forces,
				  ScalarStructure* energies);

    virtual void getParameters(std::vector<Parameter>& ) const; 
    static const std::string keyword;
    virtual std::string getKeyword() const{return keyword;}
    virtual std::string getIdNoAlias() const{return keyword;}

  private:
    virtual Force* doMake(const std::vector<Value> &values) const
    {
      return new HarmonicTrapForce(values[0]);
    };

    std::string ht_def_filename;
    std::vector<double> omega;
  };
  
  template <class TBoundaryConditions>
  const std::string HarmonicTrapForce<TBoundaryConditions>::keyword("HarmonicTrapForce");  

  template<class TBoundaryConditions>
  inline void HarmonicTrapForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
						      const Vector3DBlock* positions,
						      const Vector3DBlock *velocities,
						      Vector3DBlock* forces,
						      ScalarStructure* energies)
  {
    double conversion = Constant::SI::KCAL * Constant::SI::AVOGADRO * 1e-10;

#pragma omp parallel for schedule(dynamic, 3)
    for(unsigned int i=0;i<topo->atoms.size();i++)
      {
	Vector3D f;
	Vector3D pos((*positions)[i]);

	for (int j = 0; j<3; j++)
	  f[j] = -pos[j] * 1e-10 * omega[j] * omega[j] * topo->atoms[0].scaledMass * 1.6605655e-27;

	(*forces)[i] += f * conversion;
      }
  }


  template<class TBoundaryConditions>
  inline void HarmonicTrapForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
							      const Vector3DBlock* positions,
							      const Vector3DBlock *velocities,
							      Vector3DBlock* forces,
							      ScalarStructure* energies)
  {
    evaluate(topo, positions, velocities, forces, energies); // Not implemented right now
  }

  template<class TBoundaryConditions>
  inline void HarmonicTrapForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const
  {
    parameters.push_back(Parameter("-ht_def_filename", Value(ht_def_filename)));
			 
  }
}
#endif
