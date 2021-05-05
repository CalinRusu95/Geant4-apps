#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "Randomize.hh"
#include "globals.hh"

class PhysicsList : public G4VModularPhysicsList {
public:
  PhysicsList();
  virtual ~PhysicsList();
  virtual void ConstructParticle();   // construction of particles
  virtual void ConstructProcess();    // construct processes and register them to particles
  virtual void SetCuts();
};

#endif

