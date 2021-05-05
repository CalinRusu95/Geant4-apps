#ifndef SteppingAction_h
#define SteppingAction_h 1

#include <iostream>
#include <fstream>

#include "G4UserSteppingAction.hh"
#include "EventAction.hh"
#include "RunInput.hh"

// In UserSteppingAction() we collect the step parameters to be updated in EventAction.
class SteppingAction : public G4UserSteppingAction {
public:
  SteppingAction(EventAction*,RunInput*);
  virtual ~SteppingAction();

  virtual void UserSteppingAction(const G4Step* step);
    
private:
  EventAction*  fEventAction;
  RunInput* runInput;
  G4String fProcList[21];
  std::vector<G4String> fPartList;
  G4int partIn, partOut;

  G4int GetProcessID(const G4StepPoint* postStepPt);
  G4int GetProcessID(const G4Track* track);
  G4int GetParticleID(const G4Track* track);
};

#endif
