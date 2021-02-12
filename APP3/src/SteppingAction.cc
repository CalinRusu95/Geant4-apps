#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"

#include <fstream>
  using namespace std;
  ofstream ofz("Distributie.txt", ofstream::out);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* RuAct)//, PrimaryGeneratorAction* gps)
:G4UserSteppingAction(),fDetector(det), fRunAction(RuAct), /*fPrimary(gps),*/ fScoringVolume(0) /*, fGPS(0)*/
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.) return;
  

  G4StepPoint* prePoint  = step->GetPreStepPoint();
  //G4StepPoint* postPoint = step->GetPostStepPoint();
  
  G4double niel = step->GetNonIonizingEnergyDeposit();
  G4double trackL = step->GetTrack()->GetTrackLength();
  fRunAction->FillEdep(edep, niel);
  fRunAction->FillTrack(trackL);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  G4int copyNb = prePoint->GetTouchableHandle()->GetCopyNumber();
  if (copyNb > 0) { fRunAction->FillTallyEdep(copyNb-1, edep);
					fRunAction->FillTallyTrack(copyNb-1,trackL);}
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////


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
   



  //const G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  //const G4ParticleDefinition* currparticle = step->GetTrack()->GetDynamicParticle()->GetParticleDefinition();

  //const G4String protname = particle->GetParticleName();
  //const G4String curpname = currparticle->GetParticleName();
	
	  
	  
  if (volume != fScoringVolume) return;
  //if (curpname != protname) return;

  
     
  //G4double kinEnergyPostStep = postPoint->GetKineticEnergy();
  G4double kinEnergyPreStep = prePoint->GetKineticEnergy();
  G4double kinEnergyPosition = prePoint->GetPosition().x();
  G4double EneDep = step->GetTotalEnergyDeposit();
  //G4double SFS = prePoint->GetControlFlag();

  ofz << kinEnergyPreStep << G4endl;// << "   \t   " << kinEnergyPosition << "   \t   " << edepz << "   \t   " << EneDep << G4endl;
  //ofz << curpname << G4endl;

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


