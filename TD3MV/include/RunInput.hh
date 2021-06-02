#ifndef RunInput_h
#define RunInput_h 1

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "globals.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

// Class that loads input parameters
class RunInput {
public:
  RunInput();
  ~RunInput();

  G4bool IsVisualOn() { return visualOn; }
  G4int GetScFoilMat() { return scFoil_MaterialIndex; }
  G4int GetMPEMat() { return MPE_MaterialIndex; }
  G4double GetFoilThickness() { return foil_thickness; }
  G4double GetMPEThickness() { return MPE_thickness; }
  G4int GetNumberOfLayers() { return number_Layers; }
  G4int GetNumEvents() { return numberOfEvent; }
  G4double BeamEnergy() { return bEnergy; }
  G4double BeamESpread() { return bESpread; }
  G4int ParticleName() { return pName; }
  void SetNumEvents(G4int v) { numberOfEvent = v; }

private:
  std::ifstream dataFile;

  G4bool visualOn;
  G4int scFoil_MaterialIndex;
  G4int MPE_MaterialIndex;
  G4double foil_thickness;
  G4double MPE_thickness;
  G4int number_Layers;
  G4double bEnergy;
  G4double bESpread;
  G4int pName;
  G4int numberOfEvent;
};

#endif
