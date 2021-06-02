#include "Randomize.hh"

#include <time.h>
#include <algorithm>
using namespace std;

#include "RunInput.hh"
#include "PhysicsList.hh"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

// #ifdef G4MULTITHREADED
// #include "G4MTRunManager.hh"
// #else
#include "G4RunManager.hh"
// #endif

#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc, char** argv){
  RunInput* runInput = new RunInput();
  G4bool VisualOn = runInput->IsVisualOn();

  CLHEP::HepRandomEngine* randGen = new CLHEP::RanecuEngine;
  long int seed = static_cast<long int> (time(NULL));   // get seed from time
  G4Random::setTheEngine(randGen);   // choose the Random engine
  G4Random::setTheSeed(seed);

  //DetectorConstruction* detector = new DetectorConstruction(runInput);

  // get number of events
  G4int numberOfEvent = runInput->GetNumEvents();

  // construct the default run manager (use G4MTRunManager for multi-threading)
// #ifdef G4MULTITHREADED
//   G4MTRunManager* runManager = new G4MTRunManager;
//   runManager->SetNumberOfThreads(8);
//   G4cout<<"Running in multi-threading mode!"<<G4endl;
// #else
  G4RunManager* runManager = new G4RunManager;
  G4cout<<"Running in sequential mode!"<<G4endl;
// #endif

  // Set mandatory initialization classes
  runManager->SetUserInitialization(new DetectorConstruction(runInput));
  runManager->SetUserInitialization(new PhysicsList());
  runManager->SetUserInitialization(new ActionInitialization(runInput, randGen));
  runManager->Initialize();                                    // initialize G4 kernel
  
  G4VisManager* visManager = new G4VisExecutive;
  if(VisualOn) {
    visManager->Initialize();
    // interactive mode : define UI session
    G4UIExecutive* UIexec = new G4UIExecutive(argc, argv);  // then "/control/execute vis.mac"
    UIexec->SessionStart();
    delete UIexec;
  } else {
    // get the pointer to the UI manager and set verbosities
//     G4UImanager* UI = G4UImanager::GetUIpointer();
//     UI->ApplyCommand("/run/verbose 1");
//     UI->ApplyCommand("/event/verbose 1");
//     UI->ApplyCommand("/tracking/verbose 1");
//     UI->ApplyCommand("/physics/verbose 1");

    // start a run
    G4cout<<G4endl<<"***Running "<<numberOfEvent<<" events with random seed "<<seed<<G4endl<<G4endl;
    runManager->BeamOn(numberOfEvent);
  }
  
  // Free the store: user actions, physics list and detector description are owned and 
  // deleted by the run manager, so they should not be deleted in the main() program!
  delete runInput;
  delete visManager;
  delete runManager;  
  
  return 0;
}
