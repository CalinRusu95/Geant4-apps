#include "TrackingAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det, RunAction* run)
:G4UserTrackingAction(),fDetector(det), fRunAction(run)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  // extract Projected Range of primary particle
  if (track->GetTrackID() == 1) {
    G4double x = track->GetPosition().z() + 0.5*fDetector->GetAbsorSizeZ();
    if(x > 0.0) fRunAction->AddProjRange(x);
    //G4AnalysisManager::Instance()->FillH1(3, x);
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

