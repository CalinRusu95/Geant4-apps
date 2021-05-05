#include "RunAction.hh"

RunAction::RunAction(RunInput* rinp) : G4UserRunAction(),  runInput(rinp) {
  G4RunManager::GetRunManager()->SetPrintProgress(10); // set printing event number per XYZ events  

  // Create analysis manager
  // The choice of analysis technology is done via selecting a namespace in Analysis.hh
  analysisManager = G4AnalysisManager::Instance();
  G4cout<<"Using for analysis manager: "<<analysisManager->GetType()<<G4endl;
  stpNtu = -1;
}

RunAction::~RunAction() {
  delete G4AnalysisManager::Instance();  
}

void RunAction::BeginOfRunAction(const G4Run* /*run*/) {
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true); //inform the runManager to save random number seed
  
  /* Ntuples get automatically attributed an integer identifier which value is returned from the "Create" function.
     The default start value is 0 and it is incremented by 1 for each next created ntuple. The start ntuple identifier
     value can be changed with the SetFirstNtupleId(G4int) function.
  */
  stpNtu = analysisManager->CreateNtuple("StpNTuple", "Geant4 NTuple of G4Steps in CSC"); // see EventAction.cc
  analysisManager->CreateNtupleDColumn(stpNtu, "TrkE");
  analysisManager->CreateNtupleDColumn(stpNtu, "TrkT");
  analysisManager->CreateNtupleIColumn(stpNtu, "par");
  analysisManager->CreateNtupleIColumn(stpNtu, "stp");
  analysisManager->CreateNtupleIColumn(stpNtu, "procT");
  analysisManager->CreateNtupleIColumn(stpNtu, "part");
  analysisManager->CreateNtupleIColumn(stpNtu, "Q");
  analysisManager->CreateNtupleIColumn(stpNtu, "A");
  analysisManager->CreateNtupleDColumn(stpNtu, "HitE");
  analysisManager->CreateNtupleDColumn(stpNtu, "HitT");
  analysisManager->CreateNtupleIColumn(stpNtu, "procP");
  analysisManager->CreateNtupleDColumn(stpNtu, "x");
  analysisManager->CreateNtupleDColumn(stpNtu, "y");
  analysisManager->CreateNtupleDColumn(stpNtu, "z");
  analysisManager->CreateNtupleDColumn(stpNtu, "edep");
  analysisManager->CreateNtupleIColumn(stpNtu, "det");
  analysisManager->CreateNtupleDColumn(stpNtu, "slen");
  analysisManager->CreateNtupleDColumn(stpNtu, "sq");

  analysisManager->FinishNtuple(stpNtu);

  // Open an output file
  G4String fileName = "NTuple";
  analysisManager->OpenFile(fileName);
}

void RunAction::EndOfRunAction(const G4Run* /*run*/) {
  // Print info on xsections:
  // PrintXSections();
  
  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();
}

// void RunAction::PrintXSections() {
//   const G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
//   G4VProcess* process = G4ProcessTable::GetProcessTable()->FindProcess("GammaFission", particle);
//   const G4Element* elem  = G4NistManager::Instance()->FindOrBuildElement(runInput->TargName());
//   const G4Material* mate = G4NistManager::Instance()->FindOrBuildMaterial(runInput->TargNistName());
//   if(!particle || !process || !elem || !mate) {
//     if(!particle) G4cout<<"RunAction Error: Particle undefined. No XSections!."<<G4endl;
//     if(!process) G4cout<<"RunAction Error: Process undefined. No XSections!."<<G4endl;
//     if(!elem)     G4cout<<"RunAction Error: Element undefined. No XSections!."<<G4endl;
//     if(!mate)     G4cout<<"RunAction Error: Material undefined. No XSections!."<<G4endl;
//   } else {
//     G4double KE = 12.*238.*MeV;   // KarpovCrossSections parameterizes xsection vs beam KA per nucleon
//     G4double xs = G4HadronicProcessStore::Instance()->GetCrossSectionPerAtom(particle, KE,  process, elem, mate);
//     xs /= millibarn;
//     G4cout<<"Karpov XSection at "<<KE<<"MeV for beam "<<runInput->BeamName()<<" on target "<<runInput->TargName()<<" is "<<xs<<"mb"<<G4endl;
//   }
// }
