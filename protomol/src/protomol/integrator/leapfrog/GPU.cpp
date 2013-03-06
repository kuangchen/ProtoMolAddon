#include <protomol/integrator/leapfrog/GPU.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

#include <protomol/force/LennardJonesForce.h>
#include <protomol/switch/CnSwitchingFunction.h>
//#include <protomol/switch/C2SwitchingFunction.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ GPU

//Cn switch parameters
Real GPU::swcoef[][MAXEQNN] = {
  {10., -15., 6., 0, 0, 0, 0},
  {35., -84., 70., -20., 0, 0, 0},
  {126., -420., 540., -315., 70., 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {1716., -9009., 20020., -24024., 16380., -6006., 924.0}
};

Real GPU::dswcoef[][MAXEQNN] = {
  {30., -60., 30., 0, 0, 0, 0},
  {140., -420., 420., -140., 0, 0, 0},
  {630., -2520., 3780., -2520., 630., 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {12012., -72072., 180180., -240240., 180180., -72072., 12012.0}
};

Real GPU::d2swcoef[][MAXEQNN] = {
  {60., -180., 120., 0, 0, 0, 0},
  {420., -1680., 2100., -840., 0, 0, 0},
  {2520., -12600., 22680., -17640., 5040., 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {72072., -504504., 1441440., -2162160., 1801800., -792792., 144144.0}
};
//

const string GPU::keyword("GPU");

GPU::GPU() :
  STSIntegrator() {  
      pbcCell = hPbcCell = pbcOrig = mDiff = NULL;

}

GPU::GPU(Real timestep, bool emu, int diag, ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces), emulate(emu), diagnostics(diag){
    //Find forces and paramiters
    findForces(overloadedForces);
    //pbc
    pbcCell = new Real[3]; hPbcCell = new Real[3]; pbcOrig = new Real[3]; mDiff = new Real[3];	
    //
}

GPU::~GPU() {
    //output stats
    if((gpuTime.getTime()).getRealTime()>0.0){
        report.precision(5);
        report <<plain<<"GPU Timing: "<<(gpuTime.getTime()).getRealTime()<<"[s]."<<endl;
    }
    //Clean up
    if(pbcCell != NULL) delete [] pbcCell; 
    if(hPbcCell != NULL) delete [] hPbcCell; 
    if(pbcOrig != NULL) delete [] pbcOrig;
    if(mDiff != NULL) delete [] mDiff;
    //

}

void GPU::initialize(ProtoMolApp *app) {
    STSIntegrator::initialize(app);
    initializeForces();
    //timers/counters for diagnostics
    gpuTime.reset();
    //PBC?
    if(((SemiGenericTopology<PeriodicBoundaryConditions>*)app->topology)->boundaryConditions.PERIODIC){
        myPeriodic = true;

        Vector3D *ei, origin;
        ei = new Vector3D[3];

        ei[0] = ((SemiGenericTopology<PeriodicBoundaryConditions>*)app->topology)->boundaryConditions.e1();
        ei[1] = ((SemiGenericTopology<PeriodicBoundaryConditions>*)app->topology)->boundaryConditions.e2();
        ei[2] = ((SemiGenericTopology<PeriodicBoundaryConditions>*)app->topology)->boundaryConditions.e3();
        origin = ((SemiGenericTopology<PeriodicBoundaryConditions>*)app->topology)->boundaryConditions.origin();
        report.precision(5);
        report <<debug(1)<<"PBC: e1: "<<ei[0].c[0]<<", "<<ei[0].c[1]<<", "<<ei[0].c[2]<<", e2: "<<
                        ei[1].c[0]<<", "<<ei[1].c[1]<<", "<<ei[1].c[2]<<", e3: "<<
                            ei[2].c[0]<<", "<<ei[2].c[1]<<", "<<ei[2].c[2]<<", origin: "<<
                            origin.c[0]<<", "<<origin.c[1]<<", "<<origin.c[2]<<"."<<endl;
        //Orthogonal?
        if((ei[0].c[1] != 0.0 || ei[0].c[2] != 0.0 || ei[1].c[0] != 0.0 || ei[1].c[2] != 0.0 || ei[2].c[0] != 0.0 || ei[2].c[1] != 0.0))
            report << error << "PBC basis not orthogonal."<<endr;
        //into variables
		pbcMin2 = fabs(ei[0].c[0]) * 0.5;
        for(int i=0;i<3;i++){
            pbcCell[i] = ei[i].c[i];
            hPbcCell[i] = pbcCell[i] * 0.5;
			if(fabs(hPbcCell[i]) > pbcMin2) pbcMin2 = fabs(hPbcCell[i]);
            pbcOrig[i] = origin.c[i];
        }
		pbcMin2 *= pbcMin2;
        delete [] ei;
    }else{
        myPeriodic = false;
        report <<debug(1)<<"Vacuum BC."<<endl;
    }
    //Switch initialize
    if(lSwitch == 3){
        ordIdx = (int)lOrder - 2;
        myIRange[0] = pow(lSwitchon - lCutoff, -(lOrder + 1));
        for (int i = 1; i < lOrder + 1; i++)
            myIRange[i] = myIRange[i - 1] / (lSwitchon - lCutoff);
    }
    //
}


void GPU::run(int numTimesteps) {
  if (numTimesteps < 1)
    return;

  preStepModify();
  //start timer
  gpuTime.start();
  //do 'n' steps
  if(emulate) doStepsE(numTimesteps);
  else doSteps(numTimesteps);
  //stop timer
  gpuTime.stop();
  postStepModify();  

}

//####GPU code######################################################################
void GPU::doSteps(int numTimesteps){

}

//####Code to emulate GPU###########################################################
void GPU::doStepsE(int numTimesteps){

  doHalfKickdoDrift();
  gpuCalculateForces();
  for (int i = 1; i < numTimesteps; i++) {
    doKickdoDrift();
    gpuCalculateForces();
  }
  doHalfKick();	

}

void GPU::gpuCalculateForces(){
    Real rDistSquared, force, energy, deriv, value;

    //clear forces/energies
    for(unsigned int i=0;i<app->positions.size()*3;i++) myForces->c[i] = 0.0;
    app->energies.clear();
    //Lennard jones
    if (myLennardJones){
        for (unsigned int i = 0; i < app->topology->atoms.size(); i++){
            for (unsigned int j = i + 1; j < app->topology->atoms.size(); j++) {
                // find minimal distance
                Real distSquared = minimalDist(&(app->positions.c[i*3]), &(app->positions.c[j*3]), pbcCell, pbcOrig, mDiff);
                if((!lSwitch || distSquared < lCutoff*lCutoff) && distSquared > 0.0){
                    rDistSquared = 1.0 / distSquared;
                    //calculate force
                    const LennardJonesParameters & params = app->topology->lennardJonesParameters(
                                                                app->topology->atoms[i].type, app->topology->atoms[j].type);
                    Real A, B;
                    ExclusionClass excl = app->topology->exclusions.check(i, j);
                    if (excl != EXCLUSION_MODIFIED) {
                        A = params.A;
                        B = params.B;
                    } else {
                        A = params.A14;
                        B = params.B14;
                    }
                    // Fast LJ computation 
                    Real r6 = rDistSquared * rDistSquared * rDistSquared;
                    Real r12 = r6 * r6;
                    Real r6B = B * r6;
                    Real r12A = A * r12;
                    energy = r12A - r6B;
                    force = 12.0 * r12A * rDistSquared - 6.0 * r6B * rDistSquared;
                    //Switch  
                    deriv = 0.0; value = 1.0;
                    if(lSwitch == 3){
                        value = 0.0;
                        if (distSquared > lCutoff * lCutoff) value = 0.0;
                        else if (distSquared >= lSwitchon * lSwitchon) {
                            Real c[MAXEQNN + 1], swDiff, sqrtd = sqrt(distSquared);

                            swDiff = sqrtd - lCutoff;
                            c[0] = pow(swDiff, lOrder);
                            for (int i = 1; i < lOrder + 2; i++)
                                c[i] = c[i - 1] * swDiff;

                            for (int i = 0; i < lOrder + 1; i++) {
                                value += swcoef[ordIdx][i] * c[i + 1] * myIRange[i];
                                deriv += dswcoef[ordIdx][i] * c[i] * myIRange[i];
                            }
                            deriv /= sqrtd;
                        }else{
                            value = 1.0;
                        }
                    }
					//save for diagnostics
					Real old_force = force; Real old_energy = energy;
                    //accumulate to forces and energies
                    force = force * value - energy * deriv;
                    energy = energy * value;
                    app->energies[ScalarStructure::LENNARDJONES] += energy;
                    for(int k=0;k<3;k++){
                        myForces->c[i*3+k] -= mDiff[k] * force;
                        myForces->c[j*3+k] += mDiff[k] * force;
                    }
                    if(diagnostics == 1) outputDiagnostics(old_force, old_energy, value, deriv, distSquared, i, j);
                }
            }
        }
    }
    //Update potential for Shadow
    myPotEnergy = app->energies.potentialEnergy();
    //
}

Real GPU::minimalDist(Real *p1, Real *p2, Real *pbcCell, Real * pbcOrig, Real *mDiff){

	Real norm2 = 0.0;
	for(int i=0;i<3;i++){
		mDiff[i] = p2[i]-p1[i];
		norm2 += mDiff[i]*mDiff[i];
	}
	//
	if(!myPeriodic) return norm2;
	else{
		if(norm2 > pbcMin2){
			norm2 = 0.0;
			for(int i=0;i<3;i++){
				while(mDiff[i] > hPbcCell[i]) mDiff[i] -= pbcCell[i];
				while(mDiff[i] < -hPbcCell[i]) mDiff[i] += pbcCell[i];
				norm2 += mDiff[i]*mDiff[i];
			}
		}
	}
    return norm2;
}

void GPU::doHalfKickdoDrift() {
  if (anyPreDriftOrNextModify()) {
    doHalfKick();
    doDriftOrNextIntegrator();
  } else {
    Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    const unsigned int count = app->positions.size();

    //  Do a half kick on beta.
    updateBeta(0.5 * h);

    for (unsigned int i = 0; i < count; ++i) {
      app->velocities[i] += (*myForces)[i] * h * 0.5 /
                            app->topology->atoms[i].scaledMass;
    }
    
    app->positions += app->velocities*h;

    buildMolecularCenterOfMass(&app->positions, app->topology);
    buildMolecularMomentum(&app->velocities, app->topology);
    postDriftOrNextModify();
  }
}

void GPU::doKickdoDrift() {
  if (anyPreDriftOrNextModify() || anyPreStepModify() ||
      anyPostStepModify()) {
    if (anyPreStepModify() || anyPostStepModify()) {
      doHalfKick();
      postStepModify();
      preStepModify();
      doHalfKick();
    } else
      doKick();
    doDriftOrNextIntegrator();
  } else {
    Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    const unsigned int count = app->positions.size();

    updateBeta(h);

    for (unsigned int i = 0; i < count; ++i) {
      app->velocities[i] +=
        (*myForces)[i] * h / app->topology->atoms[i].scaledMass;
      //  app->positions[i] += app->velocities[i] * h;
    }

    app->positions += app->velocities*h;

    buildMolecularCenterOfMass(&app->positions, app->topology);
    buildMolecularMomentum(&app->velocities, app->topology);
    postDriftOrNextModify();
  }
}

void GPU::outputDiagnostics(Real force, Real energy, Real value, Real deriv, Real distSquared, int i, int j){

    ExclusionClass ec = app->topology->exclusions.check(i, j);
    LennardJonesForce hForce;
    Vector3D rij = app->topology->minimalDifference(app->positions[i], app->positions[j]);
    Real a = rij.normSquared();
    Real rawE = 0.0, rawF = 0.0;
    hForce(rawE, rawF, a, 1.0 / a, rij, app->topology, i, j, ec);
    report.precision(5);
    report <<hint<<"force= "<<force<<" energy= "<<energy<<" LF force= "<<rawF<<" LJ energy= "<<rawE<<endl;
    report <<hint<<"dist2= "<<distSquared<<" a= "<<a<<endl;
    report <<hint<<"Diff = "<<mDiff[0]<<", "<<mDiff[1]<<", "<<mDiff[2]<<endl;
    report <<hint<<"rij = "<<rij.c[0]<<", "<<rij.c[1]<<", "<<rij.c[2]<<endl;
    report <<hint<<"pos1= "<<app->positions.c[i*3]<<", "<<app->positions.c[i*3+1]<<", "<<app->positions.c[i*3+2]<<", "<<endl;
    report <<hint<<"pos2= "<<app->positions.c[j*3]<<", "<<app->positions.c[j*3+1]<<", "<<app->positions.c[j*3+2]<<", "<<endl;
    //C2SwitchingFunction c2sf(lSwitchon, lCutoff);
    CnSwitchingFunction cnsf(lSwitchon, lCutoff, lOrder, lSwitchoff);
    Real swtchV, swtchD;
    cnsf(swtchV, swtchD, a);
    report <<hint<<"swV= "<<value<<" swD= "<<deriv<<" CswV= "<<swtchV<<" CswD= "<<swtchD<<" order= "<<lOrder<<endl;
}

//  --------------------------------------------------------------------  //
//  This function is necessary to compute the shadow Hamiltonian and it   //
//  is integrator specific.  This version is written to work with LF.     //
//  Update beta: beta -= dt * ( q * F + 2 U )                             //
//  --------------------------------------------------------------------  //

void GPU::updateBeta(Real dt) {
  //  ----------------------------------------------------------------  //
  //  The shadow calculation is done in a postStep modifier.  If there  //
  //  aren't any, then obviously we don't need to do this calculation.  //
  //  It's possible that a different poststep modifier could make this  //
  //  execute, but no harm would be done ... only some extra cycles.    //
  //  ----------------------------------------------------------------  //

  if (!(anyPostStepModify() || top()->anyPostStepModify()))
    return;

  Real posDotF = 0.;

  for (unsigned int i = 0; i < app->positions.size(); i++)
    posDotF += app->positions[i].dot((*myForces)[i]);

  myBeta -= dt * (posDotF + 2. * myPotEnergy);
}

//####Find forces###################################################################
void GPU::findForces(ForceGroup *overloadedForces) {

  vector<Force *> ListForces = overloadedForces->getForces();
  //
  lCutoff = lSwitchon = lSwitch = cCutoff = cSwitchon = cSwitch = 0.0;
  lOrder = cOrder = lSwitchoff = cSwitchoff = 0.0;
  D = 78.0; S = 0.3; epsi = 1.0;
  myBond = myAngle = myCoulomb = myCoulombDielec = myCoulombSCPISM =
                                                     myCoulombBornRadii =
                                                       myLennardJones =
                                                         myDihedral =
                                                           myImproper = false;
  for (unsigned int i = 0; i < ListForces.size(); i++) {
    if (equalNocase(ListForces[i]->getId(), "Bond"))
      myBond = true;
    else if (equalNocase(ListForces[i]->getId(), "Angle"))
      myAngle = true;
    else if (equalStartNocase("CoulombDiElec", ListForces[i]->getId())) {
      myCoulombDielec = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          cCutoff = Fparam[j].value;
          if (cSwitch == 0) cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-switchon")) {
          cSwitchon = Fparam[j].value;
          if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
          else if (equalStartNocase("Cn", Fparam[j].text))
            cSwitch = 3;
          else cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-n")) {
          cOrder = Fparam[j].value;
          cSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
          cSwitchoff = Fparam[j].value;
          cSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-D"))
          D = Fparam[j].value;
        else if (equalNocase(Fparam[j].keyword, "-S"))
          S = Fparam[j].value;
        else if (equalNocase(Fparam[j].keyword, "-EPS"))
          epsi = Fparam[j].value;
      }
    } else if (equalStartNocase("CoulombSCPISM", ListForces[i]->getId())) {
      myCoulombSCPISM = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          cCutoff = Fparam[j].value;
          if (cSwitch == 0) cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-switchon")) {
          cSwitchon = Fparam[j].value;
          if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
          else if (equalStartNocase("Cn", Fparam[j].text)) cSwitch = 3;
          else cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-n")) {
          cOrder = Fparam[j].value;
          cSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
          cSwitchoff = Fparam[j].value;
          cSwitch = 3;
        }
      }
    } else if (equalStartNocase("CoulombBornRadii", ListForces[i]->getId())) {
      myCoulombBornRadii = true;
      swt = 1;         //default
    } else if (equalStartNocase("Coulomb", ListForces[i]->getId())) {
      myCoulomb = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          cCutoff = Fparam[j].value;
          if (cSwitch == 0) cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-switchon")) {
          cSwitchon = Fparam[j].value;
          if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
          else if (equalStartNocase("Cn", Fparam[j].text)) cSwitch = 3;
          else cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-n")) {
          cOrder = Fparam[j].value;
          cSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
          cSwitchoff = Fparam[j].value;
          cSwitch = 3;
        }
      }
    } else if (equalStartNocase("LennardJones", ListForces[i]->getId())) {
      myLennardJones = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          lCutoff = Fparam[j].value;
          if (lSwitch == 0) lSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-switchon")) {
          lSwitchon = Fparam[j].value;
          if (equalStartNocase("C2", Fparam[j].text)) lSwitch = 2;
          else if (equalStartNocase("Cn", Fparam[j].text)) lSwitch = 3;
          else lSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-n")) {
          lOrder = Fparam[j].value;
          lSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
          lSwitchoff = Fparam[j].value;
          lSwitch = 3;
        }
      }
    } else if (equalStartNocase("Dihedral", ListForces[i]->getId()))
      myDihedral = true;
    else if (equalStartNocase("Improper", ListForces[i]->getId()))
      myImproper = true;
  }
}

//####Standard integrator methods###################################################

void GPU::getParameters(vector<Parameter> &parameters)
const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back(
      Parameter("emulate",Value(emulate,ConstraintValueType::NoConstraints()),false,Text("Emulate GPU")));
  parameters.push_back(
      Parameter("diagnostics",Value(diagnostics,ConstraintValueType::NotNegative()),0,Text("Diagnostics?")));
}

STSIntegrator *GPU::doMake(const vector<Value> &values,
                                          ForceGroup *fg) const {
  return new GPU(values[0],values[1],values[2], fg);
}

