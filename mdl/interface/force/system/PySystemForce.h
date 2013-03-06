#ifndef PYSYSTEMFORCE_H
#define PYSYSTEMFORCE_H


#include <protomol/force/system/SystemForce.h>

namespace ProtoMol {

  class PySystemForce : public SystemForce {

  public:
    PySystemForce() : SystemForce() {}
    PySystemForce(PyObject* obj) : myPythonForce(obj) {}
    void evaluate(const GenericTopology *topo,
		  const Vector3DBlock *positions,
		  Vector3DBlock *forces,
		  ScalarStructure *energies);
    
    // Just adding to override, these will never be invoked.
    virtual std::string getKeyword() const {return "PySystemForce";}
    virtual std::string getIdNoAlias() const {return "PySystemForce";}
    virtual void getParameters(std::vector<Parameter> &) const {}
    virtual Force *doMake(const std::vector<Value> &) const {
      return new PySystemForce();
    }
    
    PyObject* myPythonForce;
  };
  

  inline  void PySystemForce::evaluate(const GenericTopology *topo,
				       const Vector3DBlock *positions,
				       Vector3DBlock *forces,
				       ScalarStructure *energies) {
    if (PyObject_GetAttrString(myPythonForce, "eval") == NULL)
	cout << "WARNING: NULL function" << endl;
    PyEval_CallObject(PyObject_GetAttrString(myPythonForce, "eval"),
		      Py_BuildValue("()"));
  }
};

#endif
