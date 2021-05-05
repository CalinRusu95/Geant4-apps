#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4ExtrudedSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "RunInput.hh"

#include <sstream>

// Detector construction class to define:
// Construct(): geometry, materials, and volumes.
// ConstructSDandField(): sensitive deterctors and the field.
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction(RunInput*);
  virtual ~DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();    
  
private:
  RunInput* runInput;
};

#endif
