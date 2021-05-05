#include "PhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"

// PARTICLES:
#include "G4Gamma.hh"
#include "G4LeptonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// PHYSICS PROCESSES:
#include "G4Decay.hh"
// GAMMA
#include "G4GammaConversion.hh"
#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"
// E+E-
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
// NEUTRONS
#include "G4HadronElasticProcess.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
// PROTONS
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4NuclearStopping.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4CascadeInterface.hh"
#include "G4ProtonInelasticCrossSection.hh"
// ALPHAS
#include "G4NuclearStopping.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiLightCrossSection.hh"
// IONS
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4IonInelasticProcess.hh"
#include "G4AtimaEnergyLossModel.hh"
#include "G4AtimaFluctuations.hh"
#include "G4hMultipleScattering.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList() {
  G4LossTableManager::Instance();
  // range in millimeters for energy cuts: EM_option1=0.7, EM_option4=0.02
  SetDefaultCutValue(0.02*mm);
}

PhysicsList::~PhysicsList() {}

void PhysicsList::ConstructParticle() {
  G4Gamma::GammaDefinition();                    // gamma
  G4LeptonConstructor::ConstructParticle();      // leptons
  G4IonConstructor::ConstructParticle();         // ions
  G4MesonConstructor::ConstructParticle();       // hadrons
  G4BaryonConstructor::ConstructParticle();
  G4ShortLivedConstructor::ConstructParticle();  // short lived
}

void PhysicsList::ConstructProcess(){
  AddTransportation(); // register the G4Transportation class with all particle classes

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4String particleName = particle->GetParticleName();
    
    G4Decay* theDecay = new G4Decay();
    if(theDecay->IsApplicable(*particle)) ph->RegisterProcess(theDecay,particle);

    if(particleName=="gamma") {    // PROCESSES FOR GAMMAS
      ph->RegisterProcess(new G4ComptonScattering(),particle);
      ph->RegisterProcess(new G4PhotoElectricEffect(),particle);
      ph->RegisterProcess(new G4GammaConversion(),particle);
    } else if(particleName=="e-" || particleName=="e+") {    // PROCESSES FOR ELECTRONS AND POSITRONS
      ph->RegisterProcess(new G4eMultipleScattering(),particle);
      ph->RegisterProcess(new G4eBremsstrahlung(),particle);
      ph->RegisterProcess(new G4eIonisation(),particle);
      if(particleName == "e+") ph->RegisterProcess(new G4eplusAnnihilation(),particle);
    }  else if(particleName=="neutron") {    // PROCESSES FOR NEUTRONS
      G4HadronElasticProcess* nElasticProc = new G4HadronElasticProcess("neutronElastic");
      G4NeutronHPElastic* nElasticModel = new G4NeutronHPElastic();
      nElasticModel->SetMinEnergy(0.001*MeV);
      nElasticProc->RegisterMe(nElasticModel);
      G4NeutronHPElasticData* nElasticData = new G4NeutronHPElasticData();
      nElasticProc->AddDataSet(nElasticData);
      ph->RegisterProcess(nElasticProc,particle);
    } else if(particleName=="proton") {
      ph->RegisterProcess(new G4hMultipleScattering(),particle);
      ph->RegisterProcess(new G4hIonisation(),particle);
      ph->RegisterProcess(new G4NuclearStopping(), particle);
      //G4ProtonInelasticProcess* theProtonInelastic = new G4ProtonInelasticProcess();
      //theProtonInelastic->RegisterMe(new G4CascadeInterface());
      //theProtonInelastic->AddDataSet(new G4ProtonInelasticCrossSection());
      //ph->RegisterProcess(theProtonInelastic,particle);
    } else if(particleName=="alpha" || particleName=="He3") {
      ph->RegisterProcess(new G4hMultipleScattering(),particle);
      ph->RegisterProcess(new G4ionIonisation(),particle);
      ph->RegisterProcess(new G4NuclearStopping(), particle);
      // Inelastic
      G4AlphaInelasticProcess* theAlphaInelastic = new G4AlphaInelasticProcess();
      theAlphaInelastic->RegisterMe(new G4BinaryLightIonReaction());
      theAlphaInelastic->AddDataSet(new G4TripathiLightCrossSection());
      ph->RegisterProcess(theAlphaInelastic,particle);
    } else if (particleName == "GenericIon") {
      G4ionIonisation* ionIoni = new G4ionIonisation("ionIonZ>3");
      ionIoni->SetEmModel(new G4AtimaEnergyLossModel());
      ionIoni->SetFluctModel(new G4AtimaFluctuations());
      ph->RegisterProcess(ionIoni,particle);   // see twiki.cern.ch/twiki/bin/view/Geant4/LoweIonParameterized
      //ph->RegisterProcess(new G4hMultipleScattering("ionmsc"), particle);
      //ph->RegisterProcess(new G4NuclearStopping(), particle);
    }
  } // iteration over all defined particles
}

void PhysicsList::SetCuts(){
  G4VModularPhysicsList::SetCuts();
  DumpCutValuesTable();
}
