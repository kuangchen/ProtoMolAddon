#include <protomol/integrator/Integrator.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/modifier/Modifier.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/Zap.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

//____ Integrator

//____  Initialize static members.
const string Integrator::scope("Integrator");
Real Integrator::myBeta = 0.;

Integrator::Integrator() :
  myPotEnergy(0), app(0), myForces(0), myForcesToEvaluate(0),
  myForward(true), myOldForces(0) {}

Integrator::Integrator(ForceGroup *forceGroup) :
  myPotEnergy(0), app(0), myForces(new Vector3DBlock),
  myForcesToEvaluate(forceGroup), myForward(true),
  myOldForces(new Vector3DBlock) {}

Integrator::~Integrator() {
  zap(myForces);
  zap(myOldForces);
  zap(myForcesToEvaluate);

  for (iterator i = myListModifiers.begin();
       i != myListModifiers.end(); ++i)
    delete (*i);
}

void Integrator::initialize(ProtoMolApp *app) {
  this->app = app;

  myForces->zero(app->positions.size());
  myOldForces->zero(app->positions.size());

  buildMolecularCenterOfMass((&app->positions), app->topology);
  buildMolecularMomentum((&app->velocities), app->topology);

  // Initialize only external modifiers,
  // where internal modifiers will be added
  // and initialized at appropriated time
  deleteInternalModifiers();
  initializeModifiers();
}

Integrator *Integrator::top() {
  Integrator *i = this;
  for (; i->previous() != 0; i = i->previous()) ;

  return i;
}

const Integrator *Integrator::top() const {
  const Integrator *i = this;
  for (; i->previous() != 0; i = i->previous()) ;

  return i;
}

Integrator *Integrator::bottom() {
  Integrator *i = this;
  for (; i->next() != 0; i = i->next()) ;

  return i;
}

const Integrator *Integrator::bottom() const {
  const Integrator *i = this;
  for (; i->next() != 0; i = i->next()) ;

  return i;
}

int Integrator::level() const {
  int n = 0;
  for (const Integrator *i = this; i->next() != 0; i = i->next())
    n++;

  return n;
}

IntegratorDefinition Integrator::getIntegratorDefinition() const {
  IntegratorDefinition tmp;

  // Integrator definition
  tmp.integrator.id = this->getId();
  this->getParameters(tmp.integrator.parameters);

  // Force definitions
  if (myForcesToEvaluate != 0)
    myForcesToEvaluate->getDefinition(tmp.forces);

  return tmp;
}

vector<IntegratorDefinition> Integrator::getIntegratorDefinitionAll() const {
  vector<IntegratorDefinition> res;
  for (const Integrator *i = bottom(); i != 0; i = i->previous())
    res.push_back(i->getIntegratorDefinition());

  return res;
}

void Integrator::uncache() {
  for (Integrator *i = top(); i != 0; i = i->next()) {
    if (i->myForcesToEvaluate != 0)
      i->myForcesToEvaluate->uncache();
    i->doUncache();
  }
}

void Integrator::forward() {
  for (Integrator *i = top(); i != 0; i = i->next())
    i->myForward = true;
}

void Integrator::backward() {
  for (Integrator *i = top(); i != 0; i = i->next())
    i->myForward = false;
}

void Integrator::preStepModify() {
  report << debug(10) << "[Integrator::preStepModify] (" << (long)this << ") "
         << (app ? app->topology->time : 0.0) << endr;

  for (iterator i = myPreStepModifiers.begin();
       i != myPreStepModifiers.end(); ++i)
    (*i)->execute(this);
}

void Integrator::preDriftOrNextModify() {
  report << debug(10) << "[Integrator::preDriftOrNextModify] (" << (long)this
         << ") " << (app ? app->topology->time : 0.0) << endr;

  for (iterator i = myPreDriftOrNextModifiers.begin();
       i != myPreDriftOrNextModifiers.end(); ++i)
    (*i)->execute(this);
}

void Integrator::postDriftOrNextModify() {
  report << debug(10) << "[Integrator::postDriftOrNextModify] ("
         << (long)this << ") "
         << (app ? app->topology->time : 0.0) << endr;

  for (iterator i = myPostDriftOrNextModifiers.begin();
       i != myPostDriftOrNextModifiers.end();
       ++i)
    (*i)->execute(this);
}

void Integrator::preForceModify() {
  report << debug(10) << "[Integrator::preForceModify] (" << (long)this << ") "
         << (app ? app->topology->time : 0.0) << endr;

  for (iterator i = myPreForceModifiers.begin();
       i != myPreForceModifiers.end();
       ++i)
    (*i)->execute(this);
}

void Integrator::mediForceModify() {
  report << debug(10) << "[Integrator::mediForceModify] (" << (long)this << ") "
         << (app ? app->topology->time : 0.0) << endr;

  for (iterator i = myMediForceModifiers.begin();
       i != myMediForceModifiers.end();
       ++i)
    (*i)->execute(this);
}

void Integrator::postForceModify() {
  report << debug(10) << "[Integrator::postForceModify] (" << (long)this << ") "
         << (app ? app->topology->time : 0.0) << endr;

  for (iterator i = myPostForceModifiers.begin();
       i != myPostForceModifiers.end(); ++i) {
    (*i)->execute(this);
  }
}

void Integrator::postStepModify() {
  report << debug(10) << "[Integrator::postStepModify] (" << (long)this << ") "
         << (app ? app->topology->time : 0.0) << endr;

  for (iterator i = myPostStepModifiers.begin();
       i != myPostStepModifiers.end(); ++i)
    (*i)->execute(this);
}

void Integrator::adoptPreStepModifier(Modifier *modifier) {
  report << debug(10) << "[Integrator::adoptPreStepModifier] "
         << modifier->getId() << "(" << (long)modifier << ") "
         << (app ? app->topology->time : 0.0) << endr;

  if (app) modifier->initialize(app, myForces);
  myPreStepModifiers.insert(modifier);

  addModifier(modifier);
}

void Integrator::adoptPreDriftOrNextModifier(Modifier *modifier) {
  report << debug(10) << "[Integrator::adoptPreDriftOrNextModifier] "
         << modifier->getId() << "(" << (long)modifier << ") "
         << (app ? app->topology->time : 0.0) << endr;

  if (app) modifier->initialize(app, myForces);
  myPreDriftOrNextModifiers.insert(modifier);
  addModifier(modifier);
}

void Integrator::adoptPostDriftOrNextModifier(Modifier *modifier) {
  report << debug(10) << "[Integrator::adoptPostDriftOrNextModifier] "
         << modifier->getId() << "(" << (long)modifier << ") "
         << (app ? app->topology->time : 0.0) << endr;

  if (app) modifier->initialize(app, myForces);
  myPostDriftOrNextModifiers.insert(modifier);
  addModifier(modifier);
}

void Integrator::adoptPreForceModifier(Modifier *modifier) {
  report << debug(10) << "[Integrator::adoptPreForceModifier] "
         << modifier->getId() << "(" << (long)modifier << ") "
         << (app ? app->topology->time : 0.0) << endr;

  if (app) modifier->initialize(app, myForces);
  myPreForceModifiers.insert(modifier);
  addModifier(modifier);
}

void Integrator::adoptMediForceModifier(Modifier *modifier) {
  report << debug(10) << "[Integrator::adoptMediForceModifier] "
         << modifier->getId() << "(" << (long)modifier << ") "
         << (app ? app->topology->time : 0.0) << endr;

  if (app) modifier->initialize(app, myForces);
  myMediForceModifiers.insert(modifier);
  addModifier(modifier);
}

void Integrator::adoptPostForceModifier(Modifier *modifier) {
  report << debug(10) << "[Integrator::adoptPostForceModifier] "
         << modifier->getId() << "(" << (long)modifier << ") "
         << (app ? app->topology->time : 0.0) << endr;

  if (app) modifier->initialize(app, myForces);
  myPostForceModifiers.insert(modifier);
  addModifier(modifier);
}

void Integrator::adoptPostStepModifier(Modifier *modifier) {
  report << debug(10) << "[Integrator::adoptPostStepModifier] "
         << modifier->getId() << "(" << (long)modifier << ") "
         << (app ? app->topology->time : 0.0) << endr;

  if (app) modifier->initialize(app, myForces);
  myPostStepModifiers.insert(modifier);
  addModifier(modifier);
}

void Integrator::deleteInternalModifiers() {
  report << debug(10) << "[Integrator::deleteInternalModifiers]" << endr;

  for (iterator i = myPreStepModifiers.begin();
       i != myPreStepModifiers.end(); ++i)
    if ((*i)->isInternal()) {
      myPreStepModifiers.erase(i);
      deleteModifier(*i);
    }

  for (iterator i = myPreDriftOrNextModifiers.begin();
       i != myPreDriftOrNextModifiers.end(); ++i)
    if ((*i)->isInternal()) {
      myPreDriftOrNextModifiers.erase(i);
      deleteModifier(*i);
    }

  for (iterator i = myPostDriftOrNextModifiers.begin();
       i != myPostDriftOrNextModifiers.end(); ++i)
    if ((*i)->isInternal()) {
      myPostDriftOrNextModifiers.erase(i);
      deleteModifier(*i);
    }

  for (iterator i = myPreForceModifiers.begin();
       i != myPreForceModifiers.end(); ++i)
    if ((*i)->isInternal()) {
      myPreForceModifiers.erase(i);
      deleteModifier(*i);
    }

  for (iterator i = myMediForceModifiers.begin();
       i != myMediForceModifiers.end(); ++i)
    if ((*i)->isInternal()) {
      myMediForceModifiers.erase(i);
      deleteModifier(*i);
    }

  for (iterator i = myPostForceModifiers.begin();
       i != myPostForceModifiers.end(); ++i)
    if ((*i)->isInternal()) {
      myPostForceModifiers.erase(i);
      deleteModifier(*i);
    }

  for (iterator i = myPostStepModifiers.begin();
       i != myPostStepModifiers.end(); ++i)
    if ((*i)->isInternal()) {
      myPostStepModifiers.erase(i);
      deleteModifier(*i);
    }
}

void Integrator::deleteExternalModifiers() {
  report << debug(10) << "[Integrator::deleteExternalModifiers] size="
         << myListModifiers.size() << endr;

  iterator i;

  for (i = myPreStepModifiers.begin(); i != myPreStepModifiers.end(); ++i)
    if (!((*i)->isInternal())) {
      myPreStepModifiers.erase(i);
      deleteModifier(*i);
    }

  for (i = myPreDriftOrNextModifiers.begin();
       i != myPreDriftOrNextModifiers.end(); ++i)
    if (!((*i)->isInternal())) {
      myPreDriftOrNextModifiers.erase(i);
      deleteModifier(*i);
    }

  for (i = myPostDriftOrNextModifiers.begin();
       i != myPostDriftOrNextModifiers.end(); ++i)
    if (!((*i)->isInternal())) {
      myPostDriftOrNextModifiers.erase(i);
      deleteModifier(*i);
    }

  for (i = myPreForceModifiers.begin(); i != myPreForceModifiers.end(); ++i)
    if (!((*i)->isInternal())) {
      myPreForceModifiers.erase(i);
      deleteModifier(*i);
    }

  for (i = myMediForceModifiers.begin(); i != myMediForceModifiers.end(); ++i)
    if (!((*i)->isInternal())) {
      myMediForceModifiers.erase(i);
      deleteModifier(*i);
    }

  for (i = myPostForceModifiers.begin(); i != myPostForceModifiers.end(); ++i)
    if (!((*i)->isInternal())) {
      myPostForceModifiers.erase(i);
      deleteModifier(*i);
    }

  for (i = myPostStepModifiers.begin(); i != myPostStepModifiers.end(); ++i)
    if (!((*i)->isInternal())) {
      myPostStepModifiers.erase(i);
      deleteModifier(*i);
    }

  report << debug(10) << "[Integrator::deleteExternalModifiers] end size="
         << myListModifiers.size() << endr;
}

bool Integrator::removeModifier(const Modifier *modifier) {
  report << debug(10) << "[Integrator::removeModifier]" << endr;

  iterator i;
  bool ok = false;
  for (i = myPreStepModifiers.begin(); i != myPreStepModifiers.end(); ++i)
    if (modifier == (*i)) {
      myPreStepModifiers.erase(i);
      ok = true;
    }

  for (i = myPreDriftOrNextModifiers.begin();
       i != myPreDriftOrNextModifiers.end(); ++i)
    if (modifier == (*i)) {
      myPreDriftOrNextModifiers.erase(i);
      ok = true;
    }

  for (i = myPostDriftOrNextModifiers.begin();
       i != myPostDriftOrNextModifiers.end(); ++i)
    if (modifier == (*i)) {
      myPostDriftOrNextModifiers.erase(i);
      ok = true;
    }

  for (i = myPreForceModifiers.begin(); i != myPreForceModifiers.end(); ++i)
    if (modifier == (*i)) {
      myPreForceModifiers.erase(i);
      ok = true;
    }

  for (i = myMediForceModifiers.begin(); i != myMediForceModifiers.end(); ++i)
    if (modifier == (*i)) {
      myMediForceModifiers.erase(i);
      ok = true;
    }

  for (i = myPostForceModifiers.begin(); i != myPostForceModifiers.end(); ++i)
    if (modifier == (*i)) {
      myPostForceModifiers.erase(i);
      ok = true;
    }

  for (i = myPostStepModifiers.begin(); i != myPostStepModifiers.end(); ++i)
    if (modifier == (*i)) {
      myPostStepModifiers.erase(i);
      ok = true;
    }

  if (ok) {
    i = myListModifiers.find(const_cast<Modifier *>(modifier));
    if (i != myListModifiers.end()) myListModifiers.erase(i);
  }

  return ok;
}

void Integrator::initializeModifiers() {
  report << debug(10) << "[Integrator::initializeModifiers] "
         << (app != 0 ? app->topology->time : 0.0) << endr;

  for (iterator i = myListModifiers.begin(); i != myListModifiers.end(); ++i)
    (*i)->initialize(app, myForces);
}

void Integrator::addModifier(Modifier *modifier) {
  report << debug(10) << "[Integrator::addModifier] size=";

  myListModifiers.insert(modifier);

  report << myListModifiers.size() << endr;
}

void Integrator::deleteModifier(Modifier *modifier) {
  report << debug(10) << "[Integrator::deleteModifier] size="
         << myListModifiers.size() << "," << modifier->isInternal() << endr;

  iterator i = myListModifiers.find(modifier);
  if (i != myListModifiers.end()) {
    report << debug(10) << "[Integrator::deleteModifier] delete"
           << (long)(modifier) << endr;

    delete modifier;
    myListModifiers.erase(i);
  }

  report << debug(10) << "[Integrator::deleteModifier] end size="
         << myListModifiers.size() << endr;
}

//____  --------------------------------------------------------------------  //
//____  The last modifier found with modifierName, is removed.                //
//____  --------------------------------------------------------------------  //

bool Integrator::removeModifier(const string modifierName) {
  report << debug(10) << "[Integrator::removeModifier]" << endr;

  bool found = false;
  Modifier *foundMod = 0;
  iterator i;

  for (iterator i = myPostStepModifiers.begin();
       i != myPostStepModifiers.end(); ++i)

    if ((*i)->getId() == modifierName) {
      myPostStepModifiers.erase(i);
      foundMod = *i;
      found = true;
    }

  for (i = myPreStepModifiers.begin(); i != myPreStepModifiers.end(); ++i)
    if ((*i)->getId() == modifierName) {
      myPreStepModifiers.erase(i);
      foundMod = *i;
      found = true;
    }

  for (i = myPreDriftOrNextModifiers.begin();
       i != myPreDriftOrNextModifiers.end(); ++i)
    if ((*i)->getId() == modifierName) {
      myPreDriftOrNextModifiers.erase(i);
      foundMod = *i;
      found = true;
    }

  for (i = myPostDriftOrNextModifiers.begin();
       i != myPostDriftOrNextModifiers.end(); ++i)
    if ((*i)->getId() == modifierName) {
      myPostDriftOrNextModifiers.erase(i);
      foundMod = *i;
      found = true;
    }

  for (i = myPreForceModifiers.begin(); i != myPreForceModifiers.end(); ++i)
    if ((*i)->getId() == modifierName) {
      myPreForceModifiers.erase(i);
      foundMod = *i;
      found = true;
    }

  for (i = myMediForceModifiers.begin(); i != myMediForceModifiers.end(); ++i)
    if ((*i)->getId() == modifierName) {
      myMediForceModifiers.erase(i);
      foundMod = *i;
      found = true;
    }

  for (i = myPostForceModifiers.begin(); i != myPostForceModifiers.end(); ++i)
    if ((*i)->getId() == modifierName) {
      myPostForceModifiers.erase(i);
      foundMod = *i;
      found = true;
    }

  if (found) {
    i = myListModifiers.find(const_cast<Modifier *>(foundMod));
    if (i != myListModifiers.end())
      myListModifiers.erase(i);
  }

  return found;
}

//____  --------------------------------------------------------------------  //
//____  These methods will save/restore the forces.  This is probably only    //
//____  going to be used by *MCIntegrator methods, so it may one day be       //
//____  moved to a more appropriate position.                                 //
//____  --------------------------------------------------------------------  //

void Integrator::saveForces() {
  *myOldForces = *myForces;
}

void Integrator::restoreForces() {
  *myForces = *myOldForces;
}

