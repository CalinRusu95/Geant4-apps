#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include <fstream>
  using namespace std;
  ofstream ofz("Distributie.txt", ofstream::out);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* RuAct)
:G4UserSteppingAction(),fDetector(det), fRunAction(RuAct), fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4double edep = step->GetTotalEnergyDeposit();
  //G4double charge = step->
  if (edep <= 0.) return;
  

  G4StepPoint* prePoint  = step->GetPreStepPoint();
  //G4StepPoint* postPoint = step->GetPostStepPoint();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  G4int copyNb = prePoint->GetTouchableHandle()->GetCopyNumber();
  if (copyNb > 0) { fRunAction->FillTallyEdep(copyNb-1, edep); }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  G4double niel = step->GetNonIonizingEnergyDeposit();
  fRunAction->FillEdep(edep, niel);
  //fRunAction->FillEdep(edep);

  
  if (step->GetTrack()->GetTrackID() == 1) {
    fRunAction->AddPrimaryStep();
    
    /*G4cout << step->GetTrack()->GetMaterial()->GetName()
           << "  E2= " << step->GetPostStepPoint()->GetKineticEnergy()
           << " Edep= " << edep 
           << " Q= " << step->GetTrack()->GetDynamicParticle()->GetCharge()
           << " Qp= " << step->GetPostStepPoint()->GetCharge()
           << G4endl; */
    
  } 
  
  
  
  
  
  if (!fScoringVolume) { 
    const DetectorConstruction* detectorConstruction
      = static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }
  
  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4double edepz = step->GetTotalEnergyDeposit();
  if (edepz <= 0.) return;
   
  if (volume != fScoringVolume) return;
  
  //G4double kinEnergyPostStep = postPoint->GetKineticEnergy();
  G4double kinEnergyPreStep = prePoint->GetKineticEnergy();
  G4double kinEnergyPosition = prePoint->GetPosition().x();
  //G4double EneDep = prePoint->GetTotalEnergyDeposit();
  
  ofz << kinEnergyPreStep << "   \t   " << kinEnergyPosition << "   \t   " << edepz << G4endl;

  //Bragg curve
  //        
  //G4double xmax = fDetector->GetAbsorSizeZ();
   
  //G4double x1 = prePoint->GetPosition().x() + xmax*0.5;
  //G4double x2 = postPoint->GetPosition().x() + xmax*0.5;
  
  //if(x1 >= 0.0 && x2 <= xmax) {  
    //G4double x  = x1 + G4UniformRand()*(x2-x1);  
	//if (step->GetTrack()->GetDefinition()->GetPDGCharge() == 0.) x = x2;
    //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    //analysisManager->FillH1(1, x, edep);  
    //analysisManager->FillH1(2, x, edep);
  //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


