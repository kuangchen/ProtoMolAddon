#include <protomol/integrator/normal/NormalModeDamping.h>
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


namespace ProtoMol {
    //__________________________________________________ NormalModeDamping

    const string NormalModeDamping::keyword( "NormalModeDamping" );

    NormalModeDamping::NormalModeDamping() : STSIntegrator(), NormalModeUtilities() 
    {
        myWriter = NULL;
        waterForces = NULL;
        velC = posC = NULL;

    }

    NormalModeDamping::NormalModeDamping(Real timestep, int firstmode, int nummode, //Real gamma, int seed, Real temperature,  
        std::string dcdf, std::string modeff, std::string atomff, Real temp,
        ForceGroup *overloadedForces) 
      : STSIntegrator(timestep,overloadedForces), NormalModeUtilities( firstmode, nummode, 91.0, 1234, 300.0), myTemperature(temp),
        myDCDFile(dcdf), modeForceFile(modeff), atomForceFile(atomff)
    {
        myWriter = NULL;
        waterForces = NULL;
        velC = posC = NULL;
    }

    NormalModeDamping::~NormalModeDamping() 
    {  
        if(myFile.is_open()){
            myFile.seekp(0, ios::beg);
            myFile << counter;
            myFile.close();
        }

        if(myWriter != NULL) delete myWriter;
        if(waterForces != NULL) delete waterForces;
        if(velC != NULL) delete velC;
        if(posC != NULL) delete posC;
    }

    void NormalModeDamping::initialize(ProtoMolApp* app){
            STSIntegrator::initialize(app);
            initializeForces();
            //NM initialization if OK
            NormalModeUtilities::initialize((int)app->positions.size(), app, myForces, NO_NM_FLAGS); //last for non-complimentary forces
            //
            //initialize minimizer noise vars
            randStp = 0.0;
            //zero instantaneous and average force Vector3DBlock
            tempV3DBlk.resize(_N);
            temp2V3DBlk.resize(_N);
            //***********************************************************	
            if(modeForceFile != ""){
                //clear file
                myFile.open(modeForceFile.c_str(),ofstream::out);
                myFile << "             " << getTimestep() << " " 
                            << firstMode-1+numMode << " " << sqrt(myTemperature * Constant::BOLTZMANN) << " " << Constant::INV_TIMEFACTOR << endl;
                //myFile.close();
                counter = 0;
            }
            //
            int dcdpos = myDCDFile.rfind(".dcd");
            string dcdx = myDCDFile; dcdx.insert(dcdpos,"Pos");
            string dcdf = myDCDFile; dcdf.insert(dcdpos,"Force");
            string dcdv = myDCDFile; dcdv.insert(dcdpos,"Vel");

            if(!myReaderX.open(dcdx)) report << error << "Cannot open position file "<<dcdx<<"."<<endr;
            //	
            if(!myReaderF.open(dcdf)) report << error << "Cannot open force file "<<dcdf<<"."<<endr;
            //
            if(!myReaderV.open(dcdv)) report << error << "Cannot open velocity file "<<dcdv<<"."<<endr;
            //
            if(atomForceFile != ""){
                //myWriter = new XYZTrajectoryWriter();
                //if(!myWriter->open(atomForceFile, app->topology->atoms, app->topology->atomTypes))
                //	report << error << "Can't open output file '"<<atomForceFile<<"'."<<endr;
            }
            //
            waterForces = new Vector3DBlock(app->positions.size());
            x0 = app->positions;
            tempV3DBlk = app->positions;
            //
            posC = new double[_3N];
            velC = new double[_3N];
            //

    }

    void NormalModeDamping::run(int numTimesteps) {

        if( numTimesteps < 1 )
            return;

        //check valid eigenvectors
        if(*Q == NULL)
            report << error << "No Eigenvectors for NormalMode integrator."<<endr;
        //main loop
        //zero average force
        tempV3DBlk.zero(_N);
        aveForceCount = 0;
        //
        Real t_fact = 0;
        if(myTemperature) t_fact = 1.0 / sqrt(myTemperature * Constant::BOLTZMANN);
        //
        for( int i = 0; i < numTimesteps; i++ ) {
            //****main loop*************************************
            //
            preStepModify();
            //if((*myReaderF >> *waterForces)){
            if((myReaderX >> app->positions) && (myReaderF >> *waterForces) && (myReaderV >> app->velocities)){
                subspaceVelocity(&app->velocities, &app->velocities);
                //report <<hint<<"[AutoCorrelatorOuter::run] in reads."<<endr;
                calculateForces();
                waterForces->intoSubtract(*myForces);
                //project forces due to waters into mode space
                modeProjector(waterForces, tmpC, true);
                //and positions
                app->positions.intoSubtract(tempV3DBlk); //tempV3DBlk gets updated during rediag!
                //myPositions->intoSubtract(x0);
                modeProjector(&app->positions, posC, false);
                app->positions.intoAdd(x0);
                //and Velocities
                modeProjector(&app->velocities, velC, false);
                //output per mode force
                if(modeForceFile != ""){
                    counter++;
                    //Output modes for analysis
                    //myFile.open(modeForceFile.c_str(),ofstream::app);
                    myFile.precision(10);
                    for(int ii=firstMode-1;ii<firstMode-1+numMode;ii++) 
                        myFile << ii+1 << " " << posC[ii] << " " << velC[ii] << " " << tmpC[ii]*t_fact << endl;
                    //close file
                    //myFile.close();
                }
                //per atom
                if(atomForceFile != "") myWriter->write(*waterForces);
                //
            }

            //   
            postStepModify();
        }	
        //
    }  

    void NormalModeDamping::modeProjector(Vector3DBlock *inp, double *tmpc, bool forcep){
        temp2V3DBlk = *inp;
        //calculate M^{1/2}(I-M^{1/2}\hat{Q}\hat{Q}^TM^{-1/2})M^{-1/2}f using BLAS
        //f'=M^{-1/2}*f
        if(forcep){ //force or pos/vel?
            for( int i=0; i < _3N; i++)
                        temp2V3DBlk.c[i] *= invSqrtMass[i/3];
        }else{
            for( int i=0; i < _3N; i++)
                        temp2V3DBlk.c[i] *= sqrtMass[i/3];
        }
        //c=hQ^T*M^{-1/2}*f
        char transA = 'T';							// Transpose Q, LAPACK checks only first character N/V
        int m = _3N; int n = _rfM; int incxy = 1;	//sizes
        double alpha = 1.0;	double beta = 0.0;		//multiplyers, see Blas docs.

        Lapack::dgemv(&transA, &m, &n, &alpha, (*Q), &m, temp2V3DBlk.c, &incxy, &beta, tmpc, &incxy);
    }

    void NormalModeDamping::getParameters(vector<Parameter>& parameters) const {
        STSIntegrator::getParameters(parameters);
        parameters.push_back(Parameter("firstmode",Value(firstMode,ConstraintValueType::NoConstraints()),-1,Text("First mode to use in set")));
        parameters.push_back(Parameter("numbermodes",Value(numMode,ConstraintValueType::NoConstraints()),-1,Text("Number of modes propagated")));
        parameters.push_back(Parameter("dcdFile",Value(myDCDFile,ConstraintValueType::NoConstraints())));
        parameters.push_back(Parameter("modeForceFile",Value(modeForceFile,ConstraintValueType::NoConstraints()),std::string(""),Text("Mode Force output filename")));
        parameters.push_back(Parameter("atomForceFile",Value(atomForceFile,ConstraintValueType::NoConstraints()),std::string(""),Text("Atom Force output filename")));
        parameters.push_back(Parameter("temperature",Value(myTemperature,ConstraintValueType::NotNegative())));
    }

    STSIntegrator* NormalModeDamping::doMake(const vector<Value>& values,ForceGroup* fg)const{
        return new NormalModeDamping(values[0],values[1],values[2],values[3],values[4],values[5],values[6],
        fg);
    }
}

