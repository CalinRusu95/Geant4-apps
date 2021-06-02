#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "RunInput.hh"
#include "Randomize.hh"
#include "G4Box.hh"
#include "G4IonTable.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"

/*
   The primary generator action class with a general particle source.
   The default kinematic is a 3 MeV proton.
*/

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction(RunInput*, CLHEP::HepRandomEngine*);
    virtual ~PrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);
    G4ParticleGun* GetParticleGun() { return fParticleGun; }

private:
    G4ParticleGun* fParticleGun;
    G4Box* fEnvelopeBox;
    RunInput* runInput;
    CLHEP::HepRandomEngine* randGen;
    G4RandGauss* enerDist;
    G4RandFlat* flatDist;
};

#endif
