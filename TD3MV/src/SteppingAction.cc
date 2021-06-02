#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4PVParameterised.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include <fstream>
#include <iomanip>
using namespace std;

ofstream rstepScoring("Results_Step_Scoring.csv", ofstream::out);
//ofstream rstepScoring_U("Results_Step_Scoring_unrestricted.csv", ofstream::out);
ofstream rstepTest("Results_Step_Test.csv", ofstream::out);
ofstream rstepSetup("Results_Step_Setup.csv", ofstream::out);

SteppingAction::SteppingAction(EventAction* eventAction, RunInput* rinp) : G4UserSteppingAction(), fEventAction(eventAction), runInput(rinp) {
  // list of processes names (to be increased by hand)
  fProcList[0] = "none"; 
  fProcList[1] = "Transportation"; 
  fProcList[2] = "conv"; 
  fProcList[3] = "compt"; 
  fProcList[4] = "eIoni";
  fProcList[5] = "eBrem"; 
  fProcList[6] = "msc"; 
  fProcList[7] = "phot"; 
  fProcList[8] = "photonNuclear"; 
  fProcList[9] = "ionIoni";
  fProcList[10] = "annihil"; 
  fProcList[11] = "neutronInelastic"; 
  fProcList[12] = "hadElastic"; 
  fProcList[13] = "NeutronArgon";
  fProcList[14] = "PhotondiNeutron"; 
  fProcList[15] = "nCapture"; 
  fProcList[16] = "nFission"; 
  fProcList[17] = "GammaFission";
  fProcList[18] = "ionIonZ>3"; 
  fProcList[19] = "hIoni";  
  fProcList[20] = "CoulombScat";

  // initial list of particle names (to be increased dynamically)
  const G4String iPartList[6] = {"proton", "e-", "e+", "neutron", "gamma", "alpha"};
  for(G4int i=0; i<6; i++) fPartList.push_back(iPartList[i]);
  partIn = 0;
  partOut = 0;
}

SteppingAction::~SteppingAction() {
  for(G4int i=0; i<21; i++) G4cout<<"Process #"<<i<<" is "<<fProcList[i]<<G4endl; G4cout<<G4endl;
  for(G4int i=0; i<fPartList.size(); i++) G4cout<<"Particle #"<<i<<" is "<<fPartList[i]<<G4endl; G4cout<<G4endl;
  G4cout<<"PARTICLE ENTERING DETECTOR:"<<partIn<<G4endl;
  G4cout<<"PARTICLE EXITING DETECTOR:"<<partOut<<G4endl;
}

void SteppingAction::UserSteppingAction(const G4Step* step) {

  //res << "ala bala portocla";

  G4StepPoint* prePoint = dynamic_cast<G4StepPoint*>(step->GetPreStepPoint());
  if(!prePoint) { G4cout<<"No point associated to step - Exiting!"<<G4endl; exit(0); }
  G4StepPoint* postPoint = dynamic_cast<G4StepPoint*>(step->GetPostStepPoint());
  if(!postPoint) { G4cout<<"No point associated to step - Exiting!"<<G4endl; exit(0); }
  G4Track* track = dynamic_cast<G4Track*>(step->GetTrack());
  if(!track) { G4cout<<"No track associated to step - Exiting!"<<G4endl; exit(0); }
  G4VPhysicalVolume* Step_Volume = dynamic_cast<G4VPhysicalVolume*>(prePoint->GetTouchableHandle()->GetVolume());
  if(!Step_Volume) { G4cout<<"No volume associated to step - Exiting!"<<G4endl; exit(0); }
  G4String volume_Name = Step_Volume->GetName();

  if(volume_Name == "World") return;    // do not store hits in the world volume!

  G4int id_Detector = -1; G4String tmpStr;
  //G4cout<<"Step in "<<volume_Name<<G4endl;
  if(volume_Name == "Scattering_foil")   id_Detector = 0;
  else if(volume_Name == "Colimator")         id_Detector = 1;
  else if(volume_Name == "Si3N4_window")      id_Detector = 2;
  else if(volume_Name == "MPE")               id_Detector = 3;
  else if(volume_Name == "World")               id_Detector = 4;
  else if(volume_Name == "Envelope_1")               id_Detector = 5;
  else if(volume_Name == "Envelope_2")               id_Detector = 6;
  else if(volume_Name.contains("Scoring_Volume_V_")) {         
    tmpStr = volume_Name; tmpStr.remove(0,17);
    id_Detector = 1000 + atoi(tmpStr);
    //G4cout << id_Detector << " \t ";
  }  
 else if(volume_Name.contains("Test_Volume_V_")) {           
    tmpStr = volume_Name; tmpStr.remove(0,14);
    id_Detector = 2000 + atoi(tmpStr);
  } else { G4cout<<"SteppingAction: ERROR - "<<volume_Name<<" is not matched! Exiting..."<<G4endl; exit(0); }

  G4ThreeVector posPt = prePoint->GetPosition();                // prePoint position
  G4ParticleDefinition* particle = track->GetDefinition();      // particle PDG data
  G4int parentID = track->GetParentID();
  //  G4int trackID = track->GetTrackID();
  G4int stepID = track->GetCurrentStepNumber();
  G4double TrkE = track->GetVertexKineticEnergy()/MeV;
  G4double TrkT = track->GetVertexMomentumDirection().theta()/rad;
  G4double HitE = prePoint->GetKineticEnergy()/MeV;
  G4double HitT = prePoint->GetMomentumDirection().theta()/rad;
  G4double HitP = prePoint->GetMomentumDirection().phi()/rad;
  G4double HitX = posPt.getX()/mm; G4double HitY = posPt.getY()/mm; G4double HitZ = posPt.getZ()/mm;
  G4double preQ = prePoint->GetCharge();
  G4int    partQ = (G4int)particle->GetPDGCharge();
  G4int    partA = particle->GetAtomicMass();
  G4int    preProc = GetProcessID(prePoint);
  //  G4int    postProc = GetProcessID(postPoint);
  G4int    trackProc = GetProcessID(track);
  G4int    partID = GetParticleID(track);
  G4double eDep = step->GetTotalEnergyDeposit(); // Energy deposit
  G4double stepL = step->GetStepLength();        // Step length
  
  G4bool Setup_Volumes = (id_Detector >= 0 && id_Detector < 4);       // Step in a geometry element
  G4bool Scoring_Volumes = (id_Detector > 999 && id_Detector < 2000);  // Step in a scoring volume
  G4bool Test_Volumes = (id_Detector > 1999 && id_Detector < 3000);   // Step in a test volume
  
    // Stats for proton in cell
  // if(partID == 0) {   // proton
  //   if(Test_Volumes) {
  //     partIn++;
  //     //G4cout << "Hit in volume:" << volume_Name << ", Detector = " << id_Detector << G4endl;
  //   }  
  // }
   
    if(Scoring_Volumes){
        //if (partID == 0) {
            G4StepPoint* ScoringprePoint = step->GetPostStepPoint();
            G4double kinEnergyPreScoring = postPoint->GetKineticEnergy();
            rstepScoring << volume_Name
                << "\t" << id_Detector
                << "\t" << eDep / MeV
                << "\t" << TrkE / MeV
                << "\t" << stepL / um
                << "\t" << kinEnergyPreScoring / MeV
                << G4endl;
       // }
        //rstepScoring_U << volume_Name
        //               //<< "\t" << GetParticleID(step)
        //               << "\t" << partID
        //               << "\t" << id_Detector
        //               << "\t" << eDep / MeV
        //               << "\t" << stepL / um
        //               << G4endl;
    }
    
    if(Test_Volumes) {
      G4StepPoint* TestprePoint = step->GetPostStepPoint();
      G4double kinEnergyPreTest = postPoint->GetKineticEnergy();
      rstepTest << volume_Name
            << "\t" << id_Detector
            << "\t" << kinEnergyPreTest / MeV
            << G4endl;  }

    if(volume_Name == "MPE") { 
      G4StepPoint* MPEprePoint = step->GetPreStepPoint();
      G4double kinEnergyPreStep = prePoint->GetKineticEnergy();
      rstepSetup << volume_Name 
                   << "\t" << id_Detector 
                   << "\t" << kinEnergyPreStep / MeV
                   << G4endl;    }

  
  //if(partID == 0){ //proton
    if(Scoring_Volumes) { 
      //G4cout<<"DING DING DING"<<G4endl;
      // store track parameters:
      fEventAction->AddTrkKin(TrkE, TrkT);
      fEventAction->AddTrkLabel(parentID,stepID);
      fEventAction->AddTrkProc(trackProc);
      fEventAction->AddParticle(partID, partQ, partA);
      // store step parameters:
      fEventAction->AddHitProc(preProc);
      fEventAction->AddHitKin(HitE, HitT);
      fEventAction->AddStepQ(preQ);
      fEventAction->AddStepLen(stepL);
      fEventAction->AddDetector(id_Detector);
      fEventAction->AddStepEne(eDep);
      fEventAction->AddPos(HitX, HitY, HitZ);
      //G4cout << volume_Name << G4endl;
    }
  //}
}

// uniquely converts process string to process integer and stores the association to a list
G4int SteppingAction::GetProcessID(const G4StepPoint* StepPt) {
  G4int procID = -1;

  G4VProcess* hitProc = const_cast<G4VProcess*>(StepPt->GetProcessDefinedStep());
  G4String ProcName = "none";
  if(hitProc) ProcName = hitProc->GetProcessName();

  for(unsigned int i=0; i<21; i++)
    if(ProcName == fProcList[i]) procID = i;
  if(procID == -1) {
    //G4cout<<"SteppingAction::GetProcessID - Process:"<<ProcName<<" is not on the list! Please add it."<<G4endl;
    exit(0);
  }

  return procID;
}

// uniquely converts process string to process integer and stores the association to a list
G4int SteppingAction::GetProcessID(const G4Track* track) {
  G4int procID = -1;

  G4VProcess* hitProc = const_cast<G4VProcess*>(track->GetCreatorProcess());
  G4String ProcName = "none";
  if(hitProc) ProcName = hitProc->GetProcessName();

  for(unsigned int i=0; i<21; i++)
    if(ProcName == fProcList[i]) procID = i;
  if(procID == -1) {
    G4cout<<"SteppingAction::GetProcessID - Process:"<<ProcName<<" is not on the list! Please add it."<<G4endl;
    exit(0);
  }

  return procID;
}

// uniquely converts process string to process integer and stores the association to a list
G4int SteppingAction::GetParticleID(const G4Track* track) {
  G4int partID = -1;

  G4ParticleDefinition* hitPart = const_cast<G4ParticleDefinition*>(track->GetDefinition());
  G4String PartName = "none";
  if(hitPart) PartName = hitPart->GetParticleName();
  unsigned int NumParts = fPartList.size();

  // find a match; if no match was found, store new process string
  for(unsigned int i=0; i<NumParts; i++)
    if(PartName == fPartList[i]) partID = i;
  if(partID == -1) {
    fPartList.push_back(PartName);
    partID = fPartList.size()-1;
  }

  return partID;
}


