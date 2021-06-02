#include "RunInput.hh"

using namespace std;

RunInput::RunInput() {
  dataFile.open("/home/calin/geant4-apps/TD3MV/RunInput.txt", ios::in | ios::binary);
  if(!dataFile) {G4cout<<"RunInput - Error: Can't open the input file."<<G4endl; exit(1);}
  
  G4int visualInt = 0;
  dataFile>>visualInt>>scFoil_MaterialIndex>>foil_thickness>>MPE_MaterialIndex>>MPE_thickness>>pName>>bEnergy>>bESpread>>number_Layers>>numberOfEvent;

  // NOW FIX THE UNITS FOR ALL INPUT GEOMETRY PARAMETERS:
  foil_thickness *= um; MPE_thickness *= um; bEnergy *= MeV; bESpread *= keV;
  visualOn = (visualInt!=0);
 
  // PRINT OUT VARIOUS INPUT INFO:
  G4cout<<"\nRunInput settings..."<<G4endl;

  if(visualOn) G4cout << "Visualisation: ON" << G4endl;
  else G4cout << "Visualisation: OFF" << G4endl;

  if (scFoil_MaterialIndex == 0) G4cout << "Scattering foil: OFF" << G4endl;
  else if (scFoil_MaterialIndex == 1) G4cout << "Scattering foil: Au, " << foil_thickness << " um - low flux" << G4endl;
  else if (scFoil_MaterialIndex == 2) G4cout << "Scattering foil: Al, " << foil_thickness << " um - high flux" << G4endl;

  G4String pName_s;
  if (pName == 0) pName_s = "proton";
  else if (pName == 1) pName_s = "alpha";
  else if (pName == 2) pName_s = "C12 ions";
  else if (pName == 3) pName_s = "Li7 ions";
  else if (pName == 4) pName_s = "B11 ions";
  
  G4cout << "Beam consists of " << pName_s << " of " << bEnergy << " MeV and " << bESpread << " MeV energy spread" << G4endl;

  G4cout << "Number of layers in detector: " << number_Layers << G4endl;
  G4cout << "Number of Events: " << numberOfEvent << G4endl;

}

RunInput::~RunInput() {
  dataFile.close();
}
