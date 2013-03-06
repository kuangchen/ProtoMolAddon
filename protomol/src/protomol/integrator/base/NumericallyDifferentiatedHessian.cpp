#include <protomol/integrator/base/NumericallyDifferentiatedHessian.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ NumericallyDifferentiatedHessian

const string
NumericallyDifferentiatedHessian::keyword("NumericallyDifferentiatedHessian");

NumericallyDifferentiatedHessian::NumericallyDifferentiatedHessian() :
  STSIntegrator() {}

NumericallyDifferentiatedHessian::NumericallyDifferentiatedHessian(
  Real timestep, Real epsil, ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces), epsilon(epsil) {
  hsn.findForces(overloadedForces);     //find forces and parameters
}

NumericallyDifferentiatedHessian::~NumericallyDifferentiatedHessian() {}

void NumericallyDifferentiatedHessian::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
  //
  _N = app->positions.size();
  _3N = 3 * _N;
  hsn.initialData(_3N);
  hsn.clear();
  //
}

void NumericallyDifferentiatedHessian::run(int numTimesteps) {
  double *numHess;
  Real maxHessError = 0.0;

  if (numTimesteps < 1) return;

  preStepModify();
  calculateForces();
  //true for mass re-weight;
  hsn.evaluate(&app->positions, app->topology, false);
  numHess = new double[_3N * _3N];
  cout << "timesteps " << numTimesteps << endl;
  for (int k = 0; k < numTimesteps; k++) {
    maxHessError = 0.0;
    if (k != 0) epsilon /= 2.0;

    //Hessian
    double* f_plus_2h = new double[_3N];
    double* f_plus_h = new double[_3N];
    double* f_x = new double[_3N];
    double* f_h = new double[_3N];
    double* f_2h = new double[_3N];
    double* orig_pos = new double[_3N];
    Real tempErr;
    for (unsigned int i = 0;i < _3N; i++)
      {
        orig_pos[i] = app->positions[i / 3][i % 3];
      }
    for (unsigned int i = 0; i < _3N; i++)
      {
	app->positions[i / 3][i % 3] = orig_pos[i] + 2.0 * epsilon;
	calculateForces();
	for (unsigned int j = 0; j < _3N; j++)
	  {
	    f_plus_2h[j] = (*myForces)[j / 3][j % 3];
	  }
	      
	app->positions[i / 3][i % 3] = orig_pos[i] + epsilon;
	calculateForces();
	for (unsigned int j = 0; j < _3N; j++)
	  {
	    f_plus_h[j] = (*myForces)[j / 3][j % 3];
	  }

	app->positions[i / 3][i % 3] = orig_pos[i] - epsilon;
	calculateForces();
	for (unsigned int j = 0; j < _3N; j++)
	  {
	    f_h[j] = (*myForces)[j / 3][j % 3];
	  }

	app->positions[i / 3][i % 3] = orig_pos[i] - 2.0 * epsilon;
	calculateForces();
	for (unsigned int j = 0; j < _3N; j++)
	  {
	    f_2h[j] = (*myForces)[j / 3][j % 3];
	  }

	// restore original position
	app->positions[i / 3][i % 3] = orig_pos[i];
	calculateForces();
	for (unsigned int j = 0; j < _3N; j++)
	  {
	    f_x[j] = (*myForces)[j / 3][j % 3];
	  }

	// Five-point stencil -- error is order h^4
	for (unsigned int j = 0; j < _3N; j++)
	  {
	    numHess[i * _3N + j] = -1.0 *
          (8.0 * f_plus_h[j] - 8.0 * f_h[j] + f_2h[j] - f_plus_2h[j]) /
          (12.0 * epsilon);
	    tempErr = fabs(numHess[i * _3N + j] - hsn.hessM[i * _3N + j]);
	    if (tempErr > maxHessError)
	      {
		maxHessError = tempErr;
		cout.precision(10);
		cout << "numerical hessian " << numHess[i * _3N + j] << endl;
		cout << "analytical hessian " << hsn.hessM[i * _3N + j] << endl;
	      }
	  }
      }

    delete[] f_plus_2h;
    delete[] f_plus_h;
    delete[] f_x;
    delete[] f_h;
    delete[] f_2h;
    delete[] orig_pos;


    report.precision(10);
    report << debug(1)
           << "[NumericallyDifferentiatedHessian::run] Hessian error = "
           << maxHessError
           << ", epsilon = " << epsilon << endr;
  }

  report << hint << "[NumericallyDifferentiatedHessian::run] Hessian error = "
         << maxHessError << ", epsilon = " << epsilon << endr;

  //remove storage
  delete[] numHess;

  //
  postStepModify();
}

void NumericallyDifferentiatedHessian::
getParameters(vector<Parameter> &parameters) const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("epsilon", Value(epsilon, ConstraintValueType::Positive()), 1.0,
               Text("epsilon")));
}

STSIntegrator *NumericallyDifferentiatedHessian::
doMake(const vector<Value> &values, ForceGroup *fg) const {
  return new NumericallyDifferentiatedHessian(values[0], values[1], fg);
}

