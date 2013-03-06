#include <protomol/integrator/normal/NormalModeBrownian.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

#include <protomol/integrator/normal/ModifierForceProjection.h>


using namespace std;

using namespace ProtoMol::Report;

using std::string;
using std::vector;


namespace ProtoMol {
	//__________________________________________________ NormalModeBrownian

	const string NormalModeBrownian::keyword( "NormalModeBrownian" );

	NormalModeBrownian::NormalModeBrownian() : STSIntegrator(), NormalModeUtilities() 
	{
		myWriter = NULL;
		myWriter2 = NULL;	

	}

	NormalModeBrownian::NormalModeBrownian(Real timestep, int firstmode, int nummode, Real gamma, int seed, Real temperature,  
		std::string avff, std::string inff,//####added avff, inff for diagnostics
		ForceGroup *overloadedForces) 
		: STSIntegrator(timestep,overloadedForces), NormalModeUtilities( firstmode, nummode, gamma, seed, temperature), 
		avForceFile(avff), inForceFile(inff) //####added avForceFile, inForceFile for diagnostics
	{
		myWriter = NULL;
		myWriter2 = NULL;
	}

	NormalModeBrownian::~NormalModeBrownian() 
	{  
		if(myWriter != NULL) delete myWriter;
		if(myWriter2 != NULL) delete myWriter2;
	}

	void NormalModeBrownian::initialize(ProtoMolApp* app){
			STSIntegrator::initialize(app);
			initializeForces();
			//NM initialization if OK
            //changed NO_NM_FLAGS to COMPLIMENT_FORCES so that this integrator uses the complement 'space' 
            //compared to the previous integrator, see NormalModeMinimizey
            NormalModeUtilities::initialize((int)app->positions.size(), app, myForces, COMPLIMENT_FORCES);
                                            //NO_NM_FLAGS); //last for non-complimentary forces
			//
			//initialize minimizer noise vars
			randStp = 0.0;
			//zero instantaneous and average force Vector3DBlock
			tempV3DBlk.resize(_N);
			temp2V3DBlk.resize(_N);
			//***********************************************************	
			//####diagnostics
			if(avForceFile != ""){
				myWriter = new XYZTrajectoryWriter();
				if(!myWriter->openWith(avForceFile, app->topology->atoms, app->topology->atomTypes))
					report << error << "Can't open output file '"<<avForceFile<<"'."<<endr;
			}

			if(inForceFile != ""){
				myWriter2 = new XYZTrajectoryWriter();
				if(!myWriter2->openWith(inForceFile, app->topology->atoms, app->topology->atomTypes))
					report << error << "Can't open output file '"<<inForceFile<<"'."<<endr;
			}

	}

	void NormalModeBrownian::run(int numTimesteps) {
		Real h = getTimestep() * Constant::INV_TIMEFACTOR;
		Real actTime;

		if( numTimesteps < 1 )
			return;

		//check valid eigenvectors
		if(*Q == NULL)
			report << error << "No Eigenvectors for NormalMode integrator."<<endr;
		//time calculated in forces! so fix here
		actTime = app->topology->time + numTimesteps * getTimestep();
		//main loop
		//zero average force
		tempV3DBlk.zero(_N);
		aveForceCount = 0;
		//
		for( int i = 0; i < numTimesteps; i++ ) {
			//****main loop*************************************
			//
			preStepModify();
			//      
			calculateForces();
			//
			for( unsigned int j = 0; j < app->positions.size(); ++j)
				app->positions[j]  += (*myForces)[j] * h / (app->topology->atoms[j].scaledMass * myGamma);
			//add random force
			genProjGauss(&gaussRandCoord1, app->topology);
			randStp = sqrt(2 * Constant::BOLTZMANN * myTemp * h / myGamma);
			app->positions.intoWeightedAdd(randStp,gaussRandCoord1);
			//
			//
			postStepModify();
		}	
		//fix average, and output
		if(aveForceCount){
			for( unsigned int i=0;i < app->positions.size(); i++) tempV3DBlk[i] /= (Real)aveForceCount;
			//####diagnostics
			if(avForceFile != "")
				myWriter->write(tempV3DBlk);

			//####

		}
		//fix time
		app->topology->time = actTime;
		//
	}  

	void NormalModeBrownian::getParameters(vector<Parameter>& parameters) const {
		STSIntegrator::getParameters(parameters);
		parameters.push_back(Parameter("firstmode",Value(firstMode,ConstraintValueType::NoConstraints()),-1,Text("First mode to use in set")));
		parameters.push_back(Parameter("numbermodes",Value(numMode,ConstraintValueType::NoConstraints()),-1,Text("Number of modes propagated")));
		parameters.push_back(Parameter("gamma",Value(myGamma*(1000 * Constant::INV_TIMEFACTOR),ConstraintValueType::NotNegative()),80.0,Text("Langevin Gamma")));
		parameters.push_back(Parameter("seed",Value(mySeed,ConstraintValueType::NotNegative()),1234,Text("Langevin random seed")));
		parameters.push_back(Parameter("temperature",Value(myTemp,ConstraintValueType::NotNegative()),300.0,Text("Langevin temperature")));
		//####diagnostics
		parameters.push_back(Parameter("avForceFile",Value(avForceFile,ConstraintValueType::NoConstraints()),std::string(""),Text("Average force filename")));
		parameters.push_back(Parameter("inForceFile",Value(inForceFile,ConstraintValueType::NoConstraints()),std::string(""),Text("Instantaneous force filename")));

		//####
	}

	STSIntegrator* NormalModeBrownian::doMake(const vector<Value>& values,ForceGroup* fg)const{
		return new NormalModeBrownian(values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7], //####last 2 for diagnostics
		fg);
	}

	void NormalModeBrownian::addModifierAfterInitialize(){
		adoptPostForceModifier(new ModifierForceProjection(this));
		STSIntegrator::addModifierAfterInitialize();
	}

	//override force projection for post force modifier, and find average in complement space
	void NormalModeBrownian::forceProjection(){
		unsigned int count = myForces->size();
		if((*Q) != NULL){
			// myForces has total forces
                       temp2V3DBlk.resize(count);
                       temp2V3DBlk.intoAssign(*myForces);
			// project myForces onto fast subspace 
			subspaceForce(myForces, myForces);
			// difference between old myForces stored in temp2V3DBlk and myForces
			// gives us instantaneous slow force
			for( unsigned int i=0;i < count; i++) temp2V3DBlk[i] -= (*myForces)[i];
			//###diagnostics
			if(inForceFile != "")
				myWriter2->write(temp2V3DBlk);

			// add slow force to running sum in tempV3DBlk[i] 
			for( unsigned int i=0;i < count; i++) tempV3DBlk[i] += temp2V3DBlk[i];
			aveForceCount++;
		}
	}


}

