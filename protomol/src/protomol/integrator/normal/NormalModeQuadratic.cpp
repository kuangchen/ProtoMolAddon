#include <protomol/integrator/normal/NormalModeQuadratic.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

#include <protomol/base/Lapack.h>

using namespace std;

using namespace ProtoMol::Report;

using std::string;
using std::vector;


namespace ProtoMol
{
  //__________________________________________________ NormalModeQuadratic

  const string NormalModeQuadratic::keyword( "NormalModeQuadratic" );

  NormalModeQuadratic::NormalModeQuadratic() : STSIntegrator(), NormalModeUtilities()
  {    
    ex0 = NULL; cPos = NULL;
  }

  NormalModeQuadratic::NormalModeQuadratic( Real timestep,
      int firstmode, int nummode, Real temperature, int cyclest, bool scle, bool stepm,
      ForceGroup *overloadedForces )
      : STSIntegrator( timestep, overloadedForces ), NormalModeUtilities( firstmode, nummode, 91.0, 1234, temperature ),
      cycleSteps( cyclest ), fScale( scle ), stepModes( stepm )
  {
    ex0 = NULL; cPos = NULL;
  }

  NormalModeQuadratic::~NormalModeQuadratic()
  {
    if ( ex0 != NULL ) {
      delete ex0; ex0 = NULL;
    }

    if ( cPos != NULL ) {
      delete [] cPos; cPos = NULL;
    }
  }

  void NormalModeQuadratic::initialize( ProtoMolApp* app )
  {
    STSIntegrator::initialize( app );
    initializeForces();

    //if previous NM then it is diagonalize, so find pointer
    if ( top() != this ) {
      myPreviousNormalMode = dynamic_cast<NormalModeUtilities*>( myPreviousIntegrator );
    } else {
      myPreviousNormalMode = 0;
    }
    //NM initialization if OK //last for complimentary forces, no gen noise
    NormalModeUtilities::initialize( ( int )app->positions.size(), app, myForces, COMPLIMENT_FORCES );

    //modes
    cPos = new double[_3N];
    std::fill( cPos, cPos + _3N, 0.0 );

    //total steps
    numSteps = 0;
    currMode = firstMode - 1;
    if ( stepModes ) {
      app->eigenInfo.currentMode = firstMode;
    }

    //save initial positions
    ex0 = new Vector3DBlock;
    *ex0 = app->positions;

    total_time = 0.0;
  }

  void NormalModeQuadratic::run( int numTimesteps )
  {

    //check valid eigenvectors
    if ( *Q == NULL ){
      report << error << "No Eigenvectors for NormalMode integrator." << endr;
    }

    if ( app->eigenInfo.myEigenvalues.size() <
         (unsigned)firstMode + numMode - 1 ){
      report << error << "Insufficient Eigenvalues for Quadratic integrator (" << app->eigenInfo.myEigenvalues.size() << "." << endr;
    }

    if ( numTimesteps < 1 ){
      return;
    }

    //timestep
    Real h = getTimestep() * Constant::INV_TIMEFACTOR;

    //time calculated in forces! so fix here
    Real actTime = app->topology->time + numTimesteps * getTimestep();

    Real tempFrq;
    Real tempKt = sqrt( 2.0 * myTemp * Constant::BOLTZMANN );

    //find ex0 if re-diag present
    if ( myPreviousNormalMode && myPreviousNormalMode->newDiag ) {
      //reset ex0 and restart time for \omega t if new diagonalization
      *ex0 = myPreviousNormalMode->diagAt;

      total_time = 0;
    }

    //main loop
    for ( int i = 0; i < numTimesteps; i++ ) {
      preStepModify();

      numSteps++;
      if ( stepModes ) {

        //set mode
        if ( !( numSteps % cycleSteps ) ) {
          cPos[currMode++] = 0.0;
          app->positions = *ex0;

          if ( currMode >= firstMode + numMode - 1 ) {
            currMode = firstMode - 1;
          }
        }

        app->eigenInfo.currentMode = currMode + 1;

        //****Analytical mode integrator loop*****
        tempFrq = sqrt( fabs( app->eigenInfo.myEigenvalues[currMode] ) );
        cPos[currMode] = tempKt * sin( ( double )( numSteps % cycleSteps ) / ( double )cycleSteps * 2.0 * M_PI );

        cPos[currMode] /= fScale ? tempFrq : sqrt( tempFrq );

        subspaceProj( cPos, &app->positions );
      } else {
        //full dynamics here
        int startm = firstMode - 1;
        if ( startm < 7 ) {
          startm = 7;
        }

        for ( int i = startm; i < firstMode + numMode - 1;i++ ) {
          //****Analytical mode integrator loop*****
          tempFrq = sqrt( fabs( app->eigenInfo.myEigenvalues[i] ) );
          cPos[i] = tempKt * sin( total_time * tempFrq );

          cPos[i] /= fScale ? tempFrq : sqrt( tempFrq );
        }

        subspaceProj( cPos, &app->positions );
      }

      postStepModify();
      total_time += h;
    }

    //fix time
    app->topology->time = actTime;

  }

  //Project from subspace to 3D space
  Vector3DBlock* NormalModeQuadratic::subspaceProj( double *tmpC, Vector3DBlock * iPos )
  {
    // Transpose Q, LAPACK checks only first character N/V
    char transA = 'N';

    //sizes
    int m = _3N; int n = _rfM; int incxy = 1;

    //multiplyers, see Blas docs.
    double alpha = 1.0; double beta = 0.0;

    Lapack::dgemv( &transA, &m, &n, &alpha, ( *Q ), &m, tmpC, &incxy, &beta, iPos->c, &incxy );

    //add ex0
    for ( int i = 0; i < _N; i++ ) {
      (*iPos)[i] /= sqrt( app->topology->atoms[i].scaledMass );
      (*iPos)[i] += (*ex0)[i];
    }

    return iPos;
  }

  void NormalModeQuadratic::getParameters( vector<Parameter>& parameters ) const
  {
    STSIntegrator::getParameters( parameters );
    parameters.push_back( Parameter(
                              "firstmode",
                              Value( firstMode, ConstraintValueType::NoConstraints() ),
                              -1,
                              Text( "First mode to use in set" )
                            )
                        );

    parameters.push_back( Parameter(
                            "numbermodes",
                            Value( numMode, ConstraintValueType::NoConstraints() ),
                            -1,
                            Text( "Number of modes propagated" )
                          )
                        );

    parameters.push_back( Parameter(
                            "temperature",
                            Value( myTemp, ConstraintValueType::NotNegative() ),
                            300.0,
                            Text( "Simulation temperature" )
                          )
                        );

    parameters.push_back( Parameter(
                            "cycleSteps",
                            Value( cycleSteps, ConstraintValueType::NotNegative() ),
                            100,
                            Text( "Number of steps per mode cycle" )
                          )
                        );

    parameters.push_back( Parameter(
                            "frequScale",
                            Value( fScale, ConstraintValueType::NoConstraints() ),
                            true,
                            Text( "Scale for frequency" )
                          )
                        );

    parameters.push_back( Parameter(
                            "stepModes",
                            Value( stepModes, ConstraintValueType::NoConstraints() ),
                            true,
                            Text( "Step through modes?" )
                          )
                        );


  }

  STSIntegrator* NormalModeQuadratic::doMake( const vector<Value>& values, ForceGroup* fg )const
  {
    return new NormalModeQuadratic( values[0], values[1], values[2], values[3],
                                    values[4], values[5], values[6], fg         );
  }
}

