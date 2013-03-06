#include <protomol/integrator/openMM/NormalModeOpenMM.h>
#include <protomol/base/Report.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

#include <LTMD/Parameters.h>

using namespace ProtoMol::Report;

namespace ProtoMol {
	const string NormalModeOpenMM::keyword( "NormalModeOpenMM" );

	NormalModeOpenMM::NormalModeOpenMM() : OpenMMIntegrator(), NormalModeUtilities() {
	
	}

	NormalModeOpenMM::NormalModeOpenMM( const std::vector<Value>& base, const std::vector<Value>& params, ForceGroup *forces )
		: OpenMMIntegrator( base, forces ), NormalModeUtilities( params[0], params[1], base[2], base[3], base[1] ) {
		isLTMD = true;
		
		mModes = params[1];
		mResiduesPerBlock = params[2];
		mBlockDOF = params[3];
		mBlockDelta = params[4];
		mRediagonalizationFrequency = params[5];
		mMinimizationLimit = params[6];
		mBlockPlatform = params[7];
		shoudForceRediagOnMinFail = params[8];
		mSDelta = params[9];
		mProtomolDiagonalize = params[10];
		mRediagOnQuadratic = params[11];
	}

	NormalModeOpenMM::~NormalModeOpenMM() {

	}

	void NormalModeOpenMM::initialize( ProtoMolApp *app ) {
		report << plain << "OpenMM NML Vector information: Vector number " << app->eigenInfo.myNumEigenvectors << ", length " << app->eigenInfo.myEigenvectorLength << "." << endr;
		
		//NM initialization
		NormalModeUtilities::initialize( ( int )app->positions.size(), app,
										 myForces, NO_NM_FLAGS );

		//Set up minimum limit
		app->eigenInfo.myMinimumLimit = mMinimizationLimit;

		//Set number of eigenvectors in use
		app->eigenInfo.myNumUsedEigenvectors = _rfM;

		// Setup LTMD Parameters
		mLTMDParameters.blockDelta = mBlockDelta * Constant::ANGSTROM_NM;
		mLTMDParameters.sDelta = mSDelta * Constant::ANGSTROM_NM;
		mLTMDParameters.bdof = mBlockDOF;
		mLTMDParameters.res_per_block = mResiduesPerBlock;
		mLTMDParameters.modes = mModes;
		mLTMDParameters.rediagFreq = mRediagonalizationFrequency;
		mLTMDParameters.minLimit = mMinimizationLimit * Constant::KCAL_KJ;

		if( mProtomolDiagonalize ) {
			mLTMDParameters.ShouldProtoMolDiagonalize = true;
			std::cout << "Block Diagonalization: ProtoMol" << std::endl;
		} else {
			mLTMDParameters.ShouldProtoMolDiagonalize = false;
			std::cout << "Block Diagonalization: OpenMM" << std::endl;
		}

		if( shoudForceRediagOnMinFail ) {
			mLTMDParameters.ShouldForceRediagOnMinFail = true;
			std::cout << "Failure Rediagonalization: True" << std::endl;
		} else {
			mLTMDParameters.ShouldForceRediagOnMinFail = false;
			std::cout << "Failure Rediagonalization: False" << std::endl;
		}
		
		if( mRediagOnQuadratic ){
			mLTMDParameters.ShouldForceRediagOnQuadratic = true;
			std::cout << "Force Rediagonalization on Quadratic Minimization: True" << std::endl;
		}else{
			mLTMDParameters.ShouldForceRediagOnQuadratic = false;
			std::cout << "Force Rediagonalization on Quadratic Minimization: False" << std::endl;
		}
		
		if( !mProtomolDiagonalize ){
			switch( mBlockPlatform ){
				case 0: // Reference
					std::cout << "OpenMM Block Diagonalization Platform: Reference" << std::endl;
					mLTMDParameters.BlockDiagonalizePlatform = OpenMM::LTMD::Preference::Reference;
					break;
				case 1:	// OpenCL
					std::cout << "OpenMM Block Diagonalization Platform: OpenCL" << std::endl;
					mLTMDParameters.BlockDiagonalizePlatform = OpenMM::LTMD::Preference::OpenCL;
					break;
				case 2: // CUDA
					std::cout << "OpenMM Block Diagonalization Platform: CUDA" << std::endl;
					mLTMDParameters.BlockDiagonalizePlatform = OpenMM::LTMD::Preference::CUDA;
					break;
			}
		}
		
		int current_res = app->topology->atoms[0].residue_seq;
		int res_size = 0;
		for( int i = 0; i < app->topology->atoms.size(); i++ ) {
			if( app->topology->atoms[i].residue_seq != current_res ) {
				mLTMDParameters.residue_sizes.push_back( res_size );
				current_res = app->topology->atoms[i].residue_seq;
				res_size = 0;
			}
			res_size++;
		}
		mLTMDParameters.residue_sizes.push_back( res_size );
		
		if( mLTMDParameters.ShouldProtoMolDiagonalize ){
			//app->eigenInfo.OpenMMMinimize = true;
		}
		
		//initialize base
		OpenMMIntegrator::initialize( app );

		initializeForces();
	}

	typedef std::vector<OpenMM::Vec3> EigenVector;

	void NormalModeOpenMM::run( int numTimesteps ) {
		if( numTimesteps < 1 ) {
			return;
		}

		//check valid eigenvectors
		if( mProtomolDiagonalize && *Q == NULL ) {
			report << error << "No Eigenvectors for NormalMode integrator." << endr;
		}

		if( mProtomolDiagonalize && app->eigenInfo.myEigVecChanged && myPreviousIntegrator != NULL ) {
			OpenMM::LTMD::Integrator *integ = dynamic_cast<OpenMM::LTMD::Integrator *>( integrator );
			if( integ ) {
				const unsigned int count = app->eigenInfo.myNumUsedEigenvectors;
				const unsigned int length = app->eigenInfo.myEigenvectorLength * 3;

				std::vector<EigenVector> vectors( count );

				for( unsigned int i = 0; i < count; i++ ) {
					vectors[i].resize( length / 3 );

					for( unsigned int j = 0; j < length; j++ ) {
						vectors[i][j / 3][j % 3] = app->eigenInfo.myEigenvectors[i * length + j];
					}
				}

				integ->setProjectionVectors( vectors );
			}

			app->eigenInfo.myEigVecChanged = false;
		}
		
		if( mLTMDParameters.ShouldProtoMolDiagonalize && app->eigenInfo.havePositionsChanged ){
			const unsigned int sz = app->positions.size();
			
			std::vector<OpenMM::Vec3> positions;
			positions.reserve( sz );
			
			OpenMM::Vec3 openMMvecp;
			for( unsigned int i = 0; i < sz; ++i ) {
				for( int j = 0; j < 3; j++ ) {
					openMMvecp[j] = app->positions[i].c[j] * Constant::ANGSTROM_NM;
				}
				positions.push_back( openMMvecp );
			}
			
			context->setPositions( positions );
			
			app->eigenInfo.havePositionsChanged = false;
		}
		
		OpenMMIntegrator::run( numTimesteps );
		
		if( mLTMDParameters.ShouldProtoMolDiagonalize && mLTMDParameters.ShouldForceRediagOnMinFail ){
			OpenMM::LTMD::Integrator *ltmd = dynamic_cast<OpenMM::LTMD::Integrator*>( integrator );
			if( ltmd ){
				const unsigned int completed = ltmd->CompletedSteps();
				const unsigned int remaining = numTimesteps - completed;
				
				if( completed != numTimesteps ) {
					app->eigenInfo.reDiagonalize = true;
					
					//fix time as no forces calculated
					app->topology->time -= remaining * getTimestep();
					
					//Fix steps
					app->currentStep -= remaining;
					std::cout << "OpenMM Failed Minimization" << std::endl;
				}
			}
		}
	}

	void NormalModeOpenMM::getParameters( vector<Parameter>& parameters ) const {
		OpenMMIntegrator::getParameters( parameters );

		// Normal Mode Utilities Parameters
		parameters.push_back( Parameter( "firstmode", Value( firstMode, ConstraintValueType::NoConstraints() ), 1, Text( "First mode to use in set" ) ) );
		parameters.push_back( Parameter( "numbermodes", Value( mModes, ConstraintValueType::NoConstraints() ), 1, Text( "Number of modes propagated" ) ) );

		// LTMD OpenMM Parameters
		parameters.push_back( Parameter( "resPerBlock", Value( mResiduesPerBlock, ConstraintValueType::NotNegative() ), 1 ) );
		parameters.push_back( Parameter( "bdof", Value( mBlockDOF, ConstraintValueType::NotNegative() ), 12 ) );
		parameters.push_back( Parameter( "blockEpsilon", Value( mBlockDelta, ConstraintValueType::NotNegative() ), 1e-3 ) );
		parameters.push_back( Parameter( "rediagFreq", Value( mRediagonalizationFrequency, ConstraintValueType::NotNegative() ), 1000 ) );
		parameters.push_back( Parameter( "minimlim", Value( mMinimizationLimit, ConstraintValueType::NotNegative() ), 0.1, Text( "Minimizer target PE difference kcal mole^{-1}" ) ) );
		parameters.push_back( Parameter( "blockHessianPlatform", Value( mBlockPlatform, ConstraintValueType::NoConstraints() ), 0 ) );
		parameters.push_back( Parameter( "forceRediagOnMinFail", Value( shoudForceRediagOnMinFail, ConstraintValueType::NoConstraints() ), false ) );
		parameters.push_back( Parameter( "sEpsilon", Value( mSDelta, ConstraintValueType::NotNegative() ), 1e-3 ) );
		parameters.push_back( Parameter( "ProtomolDiag", Value( mProtomolDiagonalize, ConstraintValueType::NoConstraints() ), false ) );
		parameters.push_back( Parameter( "forceRediagOnQuadratic", Value( mRediagOnQuadratic, ConstraintValueType::NoConstraints() ), true ) );
	}

	STSIntegrator *NormalModeOpenMM::doMake( const vector<Value>& values, ForceGroup *fg ) const {
		const std::vector<Value> base( values.begin(), values.begin() + OpenMMIntegrator::getParameterSize() );
		const std::vector<Value> params( values.begin() + OpenMMIntegrator::getParameterSize(), values.end() );
		
		return ( STSIntegrator * ) new NormalModeOpenMM( base, params, fg );
	}

	unsigned int NormalModeOpenMM::getParameterSize() const {
		return OpenMMIntegrator::getParameterSize() + 11;
	}
}

