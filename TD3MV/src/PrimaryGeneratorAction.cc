#include "PrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(RunInput* rinp, CLHEP::HepRandomEngine* rgen) : G4VUserPrimaryGeneratorAction(), fParticleGun(0), fEnvelopeBox(0), runInput(rinp), randGen(rgen) {

    //RunInput* rinp = new RunInput();
    G4cout << "PrimaryGeneratorAction: beam of ";
    G4int pName_index = runInput->ParticleName();
    //G4int bZ = 0, bA = 0; G4double bQ = 0. * eplus;
    //G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    //G4String particle = "geantino";
    if (pName_index == 0) {   // proton beam
        G4cout << "protons" << G4endl;
        //G4String particle = "proton";
    }
    else if (pName_index == 1) {   // alpha beam
        G4cout << "alphas" << G4endl;
        //G4String particle = "alpha";
    }
    else if (pName_index == 2) {   // alpha beam
        G4cout << "C12 ions" << G4endl;
        //G4String particle = "alpha";
    }
    else if (pName_index == 3) {   // alpha beam
        //bZ = 3; bA = 7; bQ = 3. * eplus;
        G4cout << "Li7 ions" << G4endl;
        //G4String particle = "alpha";
    }
    else if (pName_index == 4) {   // alpha beam
        //bZ = 5; bA = 11; bQ = 5. * eplus;
        G4cout << "B11 ions" << G4endl;
        //G4String particle = "alpha";
    }
    else {
        G4cout << G4endl;
        G4cout << "PrimaryGeneratorAction - ERROR: Unknown beam!" << G4endl;
        exit(-1);
    }

    fParticleGun = new G4ParticleGun(1);    // with one particle
    //G4ParticleDefinition* part = particleTable->FindParticle(particle);
    //G4ParticleDefinition* part = G4IonTable::GetIonTable()->GetIon(bZ, bA, 0.*keV);

    //fParticleGun->SetParticleDefinition(part);
    //fParticleGun->SetParticleCharge(bQ);
    enerDist = new G4RandGauss(randGen, runInput->BeamEnergy(), runInput->BeamESpread());
    flatDist = new G4RandFlat(randGen, 0., 1.);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

    //RunInput* rinp = new RunInput();
    //G4cout << "PrimaryGeneratorAction: beam of ";
    G4int pName_index = runInput->ParticleName();
    G4int bZ = 0, bA = 0; G4double bQ = 0. * eplus;
    //G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    //G4String particle = "geantino";
    if (pName_index == 0) {   // proton beam
        bZ = 1; bA = 1; bQ = 1. * eplus;
        //G4cout << "protons" << G4endl;
        //G4String particle = "proton";
    }
    else if (pName_index == 1) {   // alpha beam
        bZ = 2; bA = 4; bQ = 2. * eplus;
        //G4cout << "alphas" << G4endl;
        //G4String particle = "alpha";
    }
    else if (pName_index == 2) {   // alpha beam
        bZ = 6; bA = 12; bQ = 0. * eplus;
        //G4cout << "C12 ions" << G4endl;
        //G4String particle = "alpha";
    }
    else if (pName_index == 3) {   // alpha beam
        bZ = 3; bA = 7; bQ = 0. * eplus;
        //G4cout << "Li7 ions" << G4endl;
        //G4String particle = "alpha";
    }
    else if (pName_index == 4) {   // alpha beam
        bZ = 5; bA = 11; bQ = 0. * eplus;
        //G4cout << "B11 ions" << G4endl;
        //G4String particle = "alpha";
    }
    else {
        G4cout << G4endl;
        G4cout << "PrimaryGeneratorAction - ERROR: Unknown beam!" << G4endl;
        exit(-1);
    }

    G4ParticleDefinition* part = G4IonTable::GetIonTable()->GetIon(bZ, bA, 0. * keV);
    fParticleGun->SetParticleDefinition(part);
    fParticleGun->SetParticleCharge(bQ);

    G4double bRad = 1.5 * mm;
    G4double radius = sqrt(bRad * flatDist->fire()); radius *= mm;
    G4double phi = 2 * (CLHEP::pi)*flatDist->fire(); phi *= rad;
    G4double px = radius * cos(phi); px *= mm;
    G4double py = radius * sin(phi); py *= mm;
    G4double pz = -2.5 * m;
    G4double bEne = enerDist->fire(); bEne *= MeV;
    fParticleGun->SetParticlePosition(G4ThreeVector(px, py, pz));
    fParticleGun->SetParticleEnergy(bEne);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));  // direction: positive z
    G4String pname = fParticleGun->GetParticleDefinition()->GetParticleName();
    //G4cout<<"Ion: " << pname << "-> (px,py,pz) [mm] = ("<<px/mm<<","<<py/mm<<","<<pz/mm<<"), E = "<< bEne/MeV << " MeV" << G4endl;

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
