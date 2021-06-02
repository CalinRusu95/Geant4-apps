#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunInput.hh"

// Action initialization class: the user defines here all the user action classes

class ActionInitialization : public G4VUserActionInitialization {
public:
    ActionInitialization(RunInput*, CLHEP::HepRandomEngine*);
    virtual ~ActionInitialization();

    //     virtual void BuildForMaster() const; // called by the master thread
    virtual void Build() const; // called by worker threads

private:
    CLHEP::HepRandomEngine* randGen;
    RunInput* runInput;
};

#endif