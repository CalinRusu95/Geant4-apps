#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "RunAction.hh"

#include <vector>

// Event action class
// It defines data members to hold the XYZ parameters to be stored.
// They are collected step by step via the functions AddXYZ().

class EventAction : public G4UserEventAction {
public:
  EventAction(RunAction*);
  virtual ~EventAction();
  
  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);

  // track functions for step ntuple:
  void AddTrkKin(G4double v1, G4double v2) { fTrkEne.push_back(v1); fTrkThe.push_back(v2); }
  void AddTrkLabel(G4int v1, G4int v2) { fParent.push_back(v1); fStep.push_back(v2); }
  void AddTrkProc(G4int v1) { fProcT.push_back(v1); }
  void AddParticle(G4int v1, G4int v2, G4int v3) { fPart.push_back(v1); fQ.push_back(v2); fAtoM.push_back(v3); }
  // step functions for step ntuple:
  void AddHitKin(G4double v1, G4double v2) { fHitEne.push_back(v1); fHitThe.push_back(v2); }
  void AddHitProc(G4int v1) { fProcP.push_back(v1); }
  void AddStepLen(G4double v1) { fSLen.push_back(v1); }
  void AddStepQ(G4double v1) { fsq.push_back(v1); }
  void AddDetector(G4int v1) { fDet.push_back(v1); }
  void AddPos(G4double v1, G4double v2, G4double v3) { fHitX.push_back(v1); fHitY.push_back(v2); fHitZ.push_back(v3); }
  void AddStepEne(G4double v1) { fEDep.push_back(v1); }
  
private:
  RunAction* fRunAction;
  // track variables in step ntuple:
  std::vector<G4double> fTrkEne;
  std::vector<G4double> fTrkThe;
  std::vector<G4int>    fParent;
  std::vector<G4int>    fStep;
  std::vector<G4int>    fProcT;
  std::vector<G4int>    fPart;
  std::vector<G4int>    fQ;
  std::vector<G4int>    fAtoM;
  // step variables in step ntuple:
  std::vector<G4double> fHitEne;
  std::vector<G4double> fHitThe;
  std::vector<G4int>    fProcP;
  std::vector<G4double> fSLen;
  std::vector<G4double> fsq;
  std::vector<G4int>    fDet;
  std::vector<G4double> fHitX;
  std::vector<G4double> fHitY;
  std::vector<G4double> fHitZ;
  std::vector<G4double> fEDep;
  
  void Print();
};

#endif

    
