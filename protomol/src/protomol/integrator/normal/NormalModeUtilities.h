/* -*- c++ -*- */
#ifndef NORMALMODEUTILITIES_H
#define NORMALMODEUTILITIES_H

#include <protomol/type/Vector3DBlock.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/ProtoMolApp.h>

#include <protomol/type/EigenvectorInfo.h>
#include <protomol/integrator/Integrator.h>

namespace ProtoMol {
  class ScalarStructure;
    /**
    *
    * Specific NormalModeUtilities routines
    *
    */
    class NormalModeUtilities{

    public:
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Types and Enums
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        enum {COMPLIMENT_FORCES = 1};
        enum {GEN_COMP_NOISE    = 2};
        enum {NO_NM_FLAGS       = 0};

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Constructors, destructors, assignment
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
        NormalModeUtilities();
        NormalModeUtilities( int firstmode, int nummode, Real gamma, int seed, Real temperature );
        //
        virtual ~NormalModeUtilities();

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // New methods of class NormalModeUtilities
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    protected:  
        void genProjGauss(Vector3DBlock *gaussRandCoord, GenericTopology *myTopo);
        void genProjGaussC(Vector3DBlock *gaussRandCoord, Vector3DBlock *gaussRandCoordm, GenericTopology *myTopo);
        void nmlDrift(Vector3DBlock *myPositions, Vector3DBlock *myVelocities, Real dt, GenericTopology *myTopo);
        bool testRemainingModes();
        int diagHessian(double *eigVecO, double *eigValO, double *hsnhessM, int dim, int &numFound);
        void absSort(double *eigVec, double *eigVal, int *eigIndx, int dim);
        double calcRayleigh(double *rQ, double *boundRq, double *hsnhessM, int numv, double raylAverage);
        int minimizer(Real peLim, int numloop, bool simpM, bool reDiag, bool nonSubspace, int *forceCalc, Real *lastLambda, 
            ScalarStructure *myEnergies, Vector3DBlock *myPositions, GenericTopology *myTopo);
        virtual void utilityCalculateForces(){};
        void getInnerHess(double *eigVec, double *hess, double *innerHess);
        void getNewEigs(double *eigVec, double *origEigVec, double *innerEigVec);
    private:
    public:
        Vector3DBlock *subspaceForce(Vector3DBlock * force, Vector3DBlock * iPforce);
        Vector3DBlock *subspaceVelocity(Vector3DBlock * force, Vector3DBlock * iPforce);
        Vector3DBlock *nonSubspaceForce(Vector3DBlock * force, Vector3DBlock * iPforce);
        Vector3DBlock *nonSubspacePosition(Vector3DBlock * force, Vector3DBlock * iPforce);
        void subSpaceSift(Vector3DBlock *velocities, Vector3DBlock *forces);
      //double *vector3DBlockTOvect(Vector3DBlock* blkDat, double* vecDat);
      //Vector3DBlock *vectTOvector3DBlock(double* vecDat, Vector3DBlock* blkDat);
        void initialize(int sz, ProtoMolApp *app, Vector3DBlock *myForces, int nm_flags);
        virtual void forceProjection();
        void setIntegratorSetPointers(Integrator *integrator, EigenvectorInfo *eipt, bool eiValid);
        Vector3DBlock *cartSpaceProj(double *tmpC, Vector3DBlock * iPos, Vector3DBlock * ex0);
        double *modeSpaceProj(double *cPos, Vector3DBlock * iPos, Vector3DBlock * ex0);

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Friends of class NormalModeUtilities
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // private and public data members
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
        //pointers to eig pointers, used by all integrators in chain
        double **Q;
        double *eigValP;
        bool *eigVecChangedP;
        unsigned int numEig;
        //molecule constants
        int _N, _3N, _rfM;
        //linear arrays for blas
      double /**tmpFX,*/ *tmpC;
        //inputs for projection/integrator
        int firstMode;
        int numMode;
        Real myGamma;
        int mySeed;
        Real myTemp;
        //eig pointers, used if topmost integrator
        unsigned int numEigvectsu;
        //compliment randoms, compatibility with NormModeInt
        Vector3DBlock tempV3DBlk; 
        Vector3DBlock temp2V3DBlk;
        bool genCompNoise;
        //Random purtubation visible to outer Diag method
        Vector3DBlock gaussRandCoord1;
        //Last diagonalization point
        Vector3DBlock diagAt;
        //flag for diag
        bool newDiag;

    protected:
        double *invSqrtMass, *sqrtMass;
        int numSteps;
        Vector3DBlock gaussRandCoord2, posTemp;
        Vector3DBlock *myForcesP;
        bool complimentForces;

    };
}
#endif
