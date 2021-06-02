#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include <fstream>
using namespace std;

ofstream revent("Results_Event.txt", ofstream::out);


EventAction::EventAction(RunAction* run) : G4UserEventAction(), fRunAction(run),
						 fTrkEne(), fTrkThe(), fParent(), fStep(), fProcT(), fPart(), fQ(), fAtoM(), 
						 fHitEne(), fHitThe(), fProcP(), fSLen(), fsq(), fDet(),
						 fHitX(), fHitY(), fHitZ(), fEDep() {}

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event* /*event*/) {
  fTrkEne.clear();
  fTrkThe.clear();
  fParent.clear();
  fStep.clear();
  fProcT.clear();
  fPart.clear();
  fQ.clear();
  fAtoM.clear();
  fHitEne.clear();
  fHitThe.clear();
  fProcP.clear();
  fSLen.clear();
  fsq.clear();
  fDet.clear();
  fHitX.clear();
  fHitY.clear();
  fHitZ.clear();
  fEDep.clear();
}

// Accumulate statistics
void EventAction::EndOfEventAction(const G4Event* event) {
  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4int stpNtu = fRunAction->GetStpNtuple();
  if(stpNtu<0) { G4cout<<"EventAction::EndOfEventAction - ERROR: invalid step ntuple!"<<G4endl; exit(-1); }

  //ofstream res2("Results.txt");
  //res2.open ("Results.txt");

  // fill the step NTuple from vectors:
  unsigned int NumEntries = fEDep.size();   // relies on same number of entries in each hit vector!
  for(unsigned int i=0; i<NumEntries; i++) {
    analysisManager->FillNtupleDColumn(stpNtu, 0,  fTrkEne[i]);   // track energy (at creation vertex)
    analysisManager->FillNtupleDColumn(stpNtu, 1,  fTrkThe[i]);   // track theta (at creation vertex)
    analysisManager->FillNtupleIColumn(stpNtu, 2,  fParent[i]);   // parent ID (track ID of the parent particle)
    analysisManager->FillNtupleIColumn(stpNtu, 3,  fStep[i]);     // step ID
    analysisManager->FillNtupleIColumn(stpNtu, 4,  fProcT[i]);    // track process
    analysisManager->FillNtupleIColumn(stpNtu, 5,  fPart[i]);     // particle number 
    analysisManager->FillNtupleIColumn(stpNtu, 6,  fQ[i]);        // particle PDG charge
    analysisManager->FillNtupleIColumn(stpNtu, 7,  fAtoM[i]);     // particle atomic mass
    analysisManager->FillNtupleDColumn(stpNtu, 8,  fHitEne[i]);   // hit energy (at step prePoint)
    analysisManager->FillNtupleDColumn(stpNtu, 9,  fHitThe[i]);   // hit energy (at step prePoint)
    analysisManager->FillNtupleIColumn(stpNtu, 10, fProcP[i]);    // postPoint process
    analysisManager->FillNtupleDColumn(stpNtu, 11, fHitX[i]);     // prePoint x-coordinate
    analysisManager->FillNtupleDColumn(stpNtu, 12, fHitY[i]);     // prePoint y-coordinate
    analysisManager->FillNtupleDColumn(stpNtu, 13, fHitZ[i]);     // prePoint z-coordinate
    analysisManager->FillNtupleDColumn(stpNtu, 14, fEDep[i]);     // step deposited energy
    analysisManager->FillNtupleIColumn(stpNtu, 15, fDet[i]);      // detector number
    analysisManager->FillNtupleDColumn(stpNtu, 16, fSLen[i]);     // step length
    analysisManager->FillNtupleDColumn(stpNtu, 17, fsq[i]);       // prepoint charge
    /*
    analysisManager->FillNtupleDColumn(stpNtu, 0, fHitX[i]);     // prePoint x-coordinate
    analysisManager->FillNtupleDColumn(stpNtu, 1, fHitY[i]);     // prePoint y-coordinate
    analysisManager->FillNtupleDColumn(stpNtu, 2, fHitZ[i]);     // prePoint z-coordinate
    analysisManager->FillNtupleDColumn(stpNtu, 3, fEDep[i]);     // step deposited energy
    */

    revent << fEDep[i] << G4endl;

    analysisManager->AddNtupleRow(stpNtu);
  }

  //   if(NumEntries>0) Print();
  //res2.close();
}

void EventAction::Print() {
  G4int fPrec = 5;
  for(unsigned int i=0; i<fTrkEne.size(); i++)
    G4cout<<"HIT: Det="<<fDet[i]<<", Part="<<fPart[i]<<", Parent="<<fParent[i]
	  <<", Step="<<fStep[i]<<", Proc(P/T)="<<fProcP[i]<<"/"<<fProcT[i]
	  <<", TE="<<std::setprecision(fPrec)<<G4BestUnit(fTrkEne[i],"Energy")
	  <<", HE="<<std::setprecision(fPrec)<<G4BestUnit(fHitEne[i],"Energy")
	  <<", DE="<<std::setprecision(fPrec)<<G4BestUnit(fEDep[i],"Energy")
	  <<", SL="<<std::setprecision(fPrec)<<G4BestUnit(fSLen[i],"Length")
	  <<", Q="<<fQ[i]<<", sq="<<fsq[i]<<", A="<<fAtoM[i]
	  <<", X/Y/Z="<<std::setprecision(fPrec)<<G4BestUnit(fHitX[i],"Length")
	  <<"/"<<std::setprecision(fPrec)<<G4BestUnit(fHitY[i],"Length")
	  <<"/"<<std::setprecision(fPrec)<<G4BestUnit(fHitZ[i],"Length")<<G4endl;
}
