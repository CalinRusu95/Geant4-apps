#ifndef RunAction_h
#define RunAction_h 1

#include "RunInput.hh"

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Analysis.hh"

// Run action class
// Histograms and ntuples are created in BeginOfRunAction().
// The histograms and ntuples are saved in the output file in a format
// accoring to a selected technology in Analysis.hh.
class RunAction : public G4UserRunAction {
  public:
    RunAction(RunInput*);
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

    G4int GetStpNtuple() {return stpNtu;}

  private:
    RunInput* runInput;
    G4AnalysisManager* analysisManager;
    G4int stpNtu;

    // void PrintXSections();
};

#endif

