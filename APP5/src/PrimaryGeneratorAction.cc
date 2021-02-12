#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4GeneralParticleSource.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
:G4VUserPrimaryGeneratorAction(),                                              
 fGPS(0),
 fDetector(det),
 fEbeamCumul(0),       
 fGunMessenger(0)
{
  fGPS  = new G4GeneralParticleSource();
  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  
  fGPS->SetParticleDefinition(particle);
  //fGPS->GetCurrentSource()->SetParticleEnergy(3*MeV);  
  fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  
  fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0,0,-705.0127*mm));           //GPS position  
  //fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0,0,681.679*mm));

    
  //create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGPS;
  delete fGunMessenger;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4double envSizeXY = 0;
  G4double envSizeZ = 0;

  G4LogicalVolume* envLV
    = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope1");

  G4Box* envBox = NULL;
  if ( envLV ) envBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  if ( envBox ) {
    envSizeXY = envBox->GetXHalfLength()*2.;
    envSizeZ = envBox->GetZHalfLength()*2.;
  }
  else  {
    G4cerr << "Envelope volume of box shape not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The GPS will be place in the center." << G4endl;
  }

  G4double size = 0.8;
  G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
  G4double y0 = size * envSizeXY * (G4UniformRand()-0.5);
  G4double z0 = -0.5 * envSizeZ;

  fGPS->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  fGPS->GeneratePrimaryVertex(anEvent);
  
  fEbeamCumul += fGPS->GetParticleEnergy(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

