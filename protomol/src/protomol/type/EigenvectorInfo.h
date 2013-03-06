/*  -*- c++ -*-  */
#ifndef EIGENVECTORINFO_H
#define EIGENVECTORINFO_H

#include <vector>

namespace ProtoMol {
	/**
	 * Container holding coordinates/Vector3D and names
	 */
	struct EigenvectorInfo {
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		EigenvectorInfo();
		EigenvectorInfo( unsigned int n, unsigned int m );

		~EigenvectorInfo();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class EigenvectorInfo
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		bool initializeEigenvectors();
		float *getFloatEigPointer();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// Eigenvector information
		unsigned int myEigenvectorLength;
		unsigned int myNumEigenvectors;
		unsigned int myNumUsedEigenvectors;
		double myMaxEigenvalue;
		double *myEigenvectors;

		//Current and original max subspace eigenvalue,
		//for adaptive timestep
		double myOrigCEigval, myNewCEigval;
		double myOrigTimestep;

		//re-diagonalization flag
		bool reDiagonalize;
		bool havePositionsChanged;
		
		bool OpenMMMinimize;
		unsigned int RediagonalizationCount;

		//OpenMM single precision interface
		float *mySingleEigs;
		bool myEigVecChanged;

		double myMinimumLimit;

		//Analytic integrator
		std::vector< double > myEigenvalues;
		int currentMode;
	};
}

#endif // EIGENVECTORINFO_H
