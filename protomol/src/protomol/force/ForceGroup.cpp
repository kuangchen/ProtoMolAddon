#include <protomol/force/ForceGroup.h>

#include <protomol/force/system/SystemForce.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/force/MollyForce.h>
#include <protomol/force/MetaForce.h>
#include <protomol/base/TimerStatistic.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/parallel/Parallel.h>

//____#define DEBUG_OUTSTANDING_MSG

#ifdef DEBUG_OUTSTANDING_MSG
#include <mpi.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

//____ ForceGroup
ForceGroup::ForceGroup() {}

ForceGroup::~ForceGroup() {
  for (list<SystemForce *>::iterator currentForce =
         mySystemForcesList.begin(); currentForce != mySystemForcesList.end();
       ++currentForce)
    delete (*currentForce);

  for (list<ExtendedForce *>::iterator currentForce =
         myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end();
       ++currentForce)
    delete (*currentForce);

  for (list<MollyForce *>::iterator currentForce = myMollyForcesList.begin();
       currentForce != myMollyForcesList.end();
       ++currentForce)
    delete (*currentForce);

  for (list<MetaForce *>::iterator currentForce = myMetaForcesList.begin();
       currentForce != myMetaForcesList.end();
       ++currentForce)
    delete (*currentForce);
}

void ForceGroup::evaluateSystemForces(ProtoMolApp *app,
                                      Vector3DBlock *forces) const {
	if (mySystemForcesList.empty()) return;
	app->topology->uncacheCellList();

	Parallel::distribute(&app->energies, forces);

	if (Parallel::isDynamic()) {
		// Collecting the number of blocks of each force.
		vector<int> blocks;
		list<SystemForce *>::const_iterator currentForce;
		for (currentForce = mySystemForcesList.begin(); currentForce != mySystemForcesList.end(); ++currentForce) {
		  blocks.push_back((*currentForce)->numberOfBlocks(app->topology, &app->positions));
		}

		Parallel::resetNext(blocks);
	}

	if (Parallel::iAmSlave()) {
		Parallel::resetNext();

		TimerStatistic::timer[TimerStatistic::FORCES].start();
		list<SystemForce *>::const_iterator currentForce;
		for (currentForce = mySystemForcesList.begin(); currentForce != mySystemForcesList.end(); ++currentForce){
			if (Parallel::isParallel()){
				(*currentForce)->parallelEvaluate(app->topology, &app->positions, forces, &app->energies);
				
				if( (*currentForce)->getId().find( "BornRadii" ) != std::string::npos ){
					/* Copy Radii
					std::vector<double> radii( app->topology->atoms.size() );
					for( unsigned int i = 0; i < app->topology->atoms.size(); i++ ){
						radii[i] = app->topology->atoms[i].mySCPISM_A->bornRadius;
					}
					
					std::vector<double> radiiResult( app->topology->atoms.size() );
					MPI_Allreduce(&radii[0], &radiiResult[0], app->topology->atoms.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
					
					for( unsigned int i = 0; i < app->topology->atoms.size(); i++ ){
						app->topology->atoms[i].mySCPISM_A->bornRadius = radiiResult[i];
					}*/
										
					/* Copy DS 
					std::vector<double> ds( app->topology->atoms.size() );
					for( unsigned int i = 0; i < app->topology->atoms.size(); i++ ){
						ds[i] = app->topology->atoms[i].mySCPISM_A->D_s;
					}
					
					std::vector<double> dsresult( app->topology->atoms.size() );
					MPI_Allreduce(&ds[0], &dsresult[0], app->topology->atoms.size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
					
					for( unsigned int i = 0; i < app->topology->atoms.size(); i++ ){
						app->topology->atoms[i].mySCPISM_A->D_s = dsresult[i];
					}*/
				}
			}else{
				(*currentForce)->evaluate(app->topology, &app->positions, forces, &app->energies);
			}
		}
		
		TimerStatistic::timer[TimerStatistic::FORCES].stop();
	}

  Parallel::reduce(&app->energies, forces);
}
//____ Evaluate all system forces in this group.

void ForceGroup::evaluateExtendedForces(ProtoMolApp *app, Vector3DBlock *forces) const {
  if (myExtendedForcesList.empty()) return;

  app->topology->uncacheCellList();

  Parallel::distribute(&app->energies, forces);

  if (Parallel::isDynamic()) {
    // Collecting the number of blocks of each force.
    vector<int> blocks;
    list<ExtendedForce *>::const_iterator currentForce;
    for (currentForce = myExtendedForcesList.begin();
         currentForce != myExtendedForcesList.end(); ++currentForce)
      blocks.push_back
        ((*currentForce)->numberOfBlocks(app->topology, &app->positions));

    Parallel::resetNext(blocks);
  }

  if (Parallel::iAmSlave()) {
    Parallel::resetNext();

  TimerStatistic::timer[TimerStatistic::FORCES].start();
  list<ExtendedForce *>::const_iterator currentForce;
  for (currentForce = myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end(); ++currentForce)
    if (Parallel::isParallel())
      (*currentForce)->parallelEvaluate(app->topology, &app->positions,
                                        &app->velocities, forces,
                                        &app->energies);
    else
      (*currentForce)->evaluate(app->topology, &app->positions,
                                &app->velocities, forces, &app->energies);

  TimerStatistic::timer[TimerStatistic::FORCES].stop();
}

#ifdef DEBUG_OUTSTANDING_MSG
  report << allnodes << plain << "Node " << Parallel::getId() << " done."
         << endr;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Status status;
  int test = 0;
  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test, &status);
  if (test != 0)
    report << plain << allnodes << "Node " << Parallel::getId() <<
      " outstanding msg from " << status.MPI_SOURCE << endr;
#endif

  Parallel::reduce(&app->energies, forces);
}

void ForceGroup::evaluateMollyForces(GenericTopology *topo,
                                     const Vector3DBlock *positions,
                                     vector<ReducedHessAngle> *angleFilter)
const {
  if (myMollyForcesList.empty())
    return;

  topo->uncacheCellList();

  TimerStatistic::timer[TimerStatistic::FORCES].start();

  for (list<MollyForce *>::const_iterator currentForce =
         myMollyForcesList.begin();
       currentForce != myMollyForcesList.end();
       ++currentForce)
    (*currentForce)->evaluate(topo, positions, angleFilter);

  TimerStatistic::timer[TimerStatistic::FORCES].stop();
}

void ForceGroup::addSystemForce(SystemForce *force) {
  if (force != NULL)
    mySystemForcesList.push_back(force);
}

void ForceGroup::addExtendedForce(ExtendedForce *force) {
  if (force != NULL)
    myExtendedForcesList.push_back(force);
}

void ForceGroup::addMollyForce(MollyForce *force) {
  if (force != NULL)
    myMollyForcesList.push_back(force);
}

void ForceGroup::addMetaForce(MetaForce *force) {
  if (force != NULL)
    myMetaForcesList.push_back(force);
}

void ForceGroup::addForce(Force *force) {
  force->addToForceGroup(this);
}

void ForceGroup::getDefinition(vector<MakeableDefinition> &forces) const {
  for (list<SystemForce *>::const_iterator currentForce =
         mySystemForcesList.begin();
       currentForce != mySystemForcesList.end();
       ++currentForce)
    forces.push_back((*currentForce)->getDefinition());

  for (list<ExtendedForce *>::const_iterator currentForce =
         myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end();
       ++currentForce)
    forces.push_back((*currentForce)->getDefinition());

  for (list<MollyForce *>::const_iterator currentForce =
         myMollyForcesList.begin();
       currentForce != myMollyForcesList.end();
       ++currentForce)
    forces.push_back((*currentForce)->getDefinition());

  for (list<MetaForce *>::const_iterator currentForce =
         myMetaForcesList.begin();
       currentForce != myMetaForcesList.end();
       ++currentForce)
    forces.push_back((*currentForce)->getDefinition());
}

void ForceGroup::uncache() {
  for (list<SystemForce *>::iterator currentForce = mySystemForcesList.begin();
       currentForce != mySystemForcesList.end();
       ++currentForce)
    (*currentForce)->uncache();

  for (list<ExtendedForce *>::iterator currentForce =
         myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end();
       ++currentForce)
    (*currentForce)->uncache();

  for (list<ExtendedForce *>::iterator currentForce =
         myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end();
       ++currentForce)
    (*currentForce)->uncache();

  for (list<MetaForce *>::iterator currentForce = myMetaForcesList.begin();
       currentForce != myMetaForcesList.end();
       ++currentForce)
    (*currentForce)->uncache();
}

vector<Force *> ForceGroup::getForces() const {
  vector<Force *> res;
  for (list<SystemForce *>::const_iterator currentForce =
         mySystemForcesList.begin();
       currentForce != mySystemForcesList.end();
       ++currentForce)
    res.push_back(*currentForce);

  for (list<ExtendedForce *>::const_iterator currentForce =
         myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end();
       ++currentForce)
    res.push_back(*currentForce);

  for (list<ExtendedForce *>::const_iterator currentForce =
         myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end();
       ++currentForce)
    res.push_back(*currentForce);

  for (list<MetaForce *>::const_iterator currentForce =
         myMetaForcesList.begin();
       currentForce != myMetaForcesList.end();
       ++currentForce)
    res.push_back(*currentForce);

  return res;
}

vector<Force *> ForceGroup::getDeepMetaForces() const {
  vector<Force *> res;
  for (list<MetaForce *>::const_iterator currentForce =
         myMetaForcesList.begin();
       currentForce != myMetaForcesList.end();
       ++currentForce)
    (*currentForce)->getDeepForces(res);

  return res;
}

