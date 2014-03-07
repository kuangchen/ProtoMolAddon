#ifndef STRAYFIELD_H
#define STRAYFIELD_H

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


namespace ProtoMol{
  
  template<class TBoundaryConditions>
  class StrayFieldForce: public ExtendedForce
  {
    

  public:
    StrayFieldForce(): sf_def_filename(""), f(0)
    {}

    StrayFieldForce(const std::string& filename): sf_def_filename(filename), f(0)
    {
      std::cout << "opening " << sf_def_filename << "\n";
      lua_State *L = luaL_newstate();
      luaL_openlibs(L);
    
      if (luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0))
	luaL_error(L, "cannot read file: %s", lua_tostring(L, -1));

      f.resize(3);
      ExtractArray(L, f, "field", "field should be an array");
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
      return new StrayFieldForce(values[0]);
    };

    std::string sf_def_filename;
    std::vector<double> f;
  };
  
  template <class TBoundaryConditions>
  const std::string StrayFieldForce<TBoundaryConditions>::keyword("StrayField");  

  template<class TBoundaryConditions>
  inline void StrayFieldForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
						      const Vector3DBlock* positions,
						      const Vector3DBlock *velocities,
						      Vector3DBlock* forces,
						      ScalarStructure* energies)
  {
    double conversion = Constant::SI::KCAL * Constant::SI::AVOGADRO * 1e-10;

    for ( unsigned int i=0 ; i<topo->atoms.size() ; i++ )
      for ( unsigned int j=0 ; j<3 ; j++ )
	(*forces)[i][j] += f[j] * 1.60217646e-19 * conversion;
  }


  template<class TBoundaryConditions>
  inline void StrayFieldForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
							      const Vector3DBlock* positions,
							      const Vector3DBlock *velocities,
							      Vector3DBlock* forces,
							      ScalarStructure* energies)
  {
    evaluate(topo, positions, velocities, forces, energies); // Not implemented right now
  }

  template<class TBoundaryConditions>
  inline void StrayFieldForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const
  {
    parameters.push_back(Parameter("-sf_def_filename", Value(sf_def_filename)));
			 
  }
}
#endif
