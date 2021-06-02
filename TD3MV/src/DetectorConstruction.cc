#include "DetectorConstruction.hh"

using namespace std;

DetectorConstruction::DetectorConstruction(RunInput* rinp): G4VUserDetectorConstruction(), runInput(rinp) {}

DetectorConstruction::~DetectorConstruction() { }

G4VPhysicalVolume* DetectorConstruction::Construct() {  
  G4NistManager* nist = G4NistManager::Instance(); // Get nist material manager
  G4bool checkOverlaps = true; // Option to switch on/off checking of volumes overlaps
  // G4double maxStep = 0.5*nm;     // limit on step size in solids (target, dump)
  // G4UserLimits* fStepLimit = new G4UserLimits(maxStep);


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ALL GEOMTRY PARAMETERS
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
 /*
     Defined materials (these will be used for the entire geometry).
     We will define define the World and Envelope_1 materials as vacuum,
     similarly with the one used in experimental conditions. The Envelope_2
     material will be a custom defined air material and the detector will
     be a custom set of defined water volumes with the mean excitation energy of
     78 eV. Additionally, the intertwined test volumes, among the water volumes
     will have the same vacuum configuration as the World/Envelope_1. All the
     other volumes will have materials similar to the experimental setup.
 */

  // Helpful elements
  G4double z, a;
  G4Element* H = new G4Element("Hydrogen", "H", z = 1, a = 1.008 * g / mole);
  G4Element* N = new G4Element("Nitrogen", "N", z = 7, a = 14.01 * g / mole);
  G4Element* O = new G4Element("Oxygen", "O", z = 8, a = 16.00 * g / mole);
  G4Element* C = new G4Element("Carbon", "C", z = 6., a = 12.01 * g / mole);
  G4Element* Si = new G4Element("Silicon", "Si", z = 14., a = 28.0855 * g / mole);
  // Vacuum (World + Envelope_1 + test volumes)
  //G4double v_density = 1.55e-5*g/cm3, v_pressure = 0.0001*pascal;
  //G4double v_density = 1.e-23*g/cm3, v_pressure = 1.e-19*pascal;
  //G4Material* vacuum = new G4Material("Vacuum", 1., 1.008*g/mole, v_density, kStateGas,  295.*kelvin,  v_pressure);
  G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
  // Air (Envelope_2)
  /* The definition below compiles but then crashes. Replaced with NIST default. Debug if really need custom!
  G4double air_temp = 298.15*kelvin, air_pres = 101325*pascal, air_density = 1.290*mg/cm3;
  G4int air_ncomp = 2, air_fractionmass;
  G4Material* Air = new G4Material("Air",  air_density,  air_ncomp,  kStateGas, air_temp,  air_pres);
  Air->AddElement(N, air_fractionmass = 0.7);
  Air->AddElement(O, air_fractionmass = 0.3);
  */
  G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");
  // Water (Scoring volumes)
  //G4Material* Water = nist->FindOrBuildMaterial("G4_Water");
  G4double w_density = 1.0*g/cm3;
  G4int w_ncomponents = 2, w_natoms;
  G4Material* Water = new G4Material("Water",  w_density = 1.0*g/cm3, w_ncomponents = 2);
  Water->AddElement(H, w_natoms = 2);
  Water->AddElement(O, w_natoms = 1);
  Water->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  // Gold - Au or Aluminium - Al (Thin scattering foil)
  G4int ScFoil_MaterialIndex = runInput->GetScFoilMat();    //1 for Au or 2 for Al.For removal, set to 0.
  G4Material* ScFoil_Material = vacuum; 
  if(ScFoil_MaterialIndex == 1) { ScFoil_Material = nist->FindOrBuildMaterial("G4_Au"); }
  else if(ScFoil_MaterialIndex == 2) { ScFoil_Material = nist->FindOrBuildMaterial("G4_Al"); }  
  // Graphite (colimator)
  G4double gr_z, gr_a, gr_density = 2.2*g/cm3;
  G4Material* Graphite = new G4Material("Graphite",  gr_z = 6., gr_a = 12.0107*g/mole, gr_density);
  // Silicon nitride - Si3N4 (vacuum - air tube interface window)
  G4double g_density = 3.44*g/cm3;
  G4int g_nel = 2, g_natoms;
  G4Material* Si3N4 = new G4Material("Si3N4",  g_density, g_nel);
  Si3N4->AddElement(Si, g_natoms = 3);                           
  Si3N4->AddElement(N, g_natoms = 4);
  // Mylar / Polyethylene (air - water interface before entrance to the scoring volume)
  G4double m_density = 1.39*g/cm3;
  G4int m_nel = 3, m_natoms;
  G4Material* Mylar = new G4Material("Mylar", m_density, m_nel);
  Mylar->AddElement(O, m_natoms = 2);                                    
  Mylar->AddElement(C, m_natoms = 5);                                    
  Mylar->AddElement(H, m_natoms = 4);
  G4Material* PE = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4int MPE_MaterialIndex = runInput->GetMPEMat();         // 1 for Mylar or 2 for PE. For removal, set to 0.
  G4Material* MPE_Material = vacuum; 
  if(MPE_MaterialIndex == 1) { MPE_Material = Mylar; }
  else if(MPE_MaterialIndex == 2) { MPE_Material = PE; }

 /*
   Dimensional parameters.
   We will define the size of all geometry components of the setup. 
   The detector will consist of an array of water + vacuum (1 um + 0.1 um) pairs, 
   built long the Z axis, therefore we need to define the number of elements.
*/ 
  // World and envelopes configuration
  G4double worldSize_Z = 6 * m;                             // World length
  G4double worldSize_XY = 0.6 * m;                          // World height and width
  G4double envelopeSize_XY = .5 * m;                        // Envelope height and width (avaliable for both)
  G4double envelope1Size_Z = 3. * m;                        // First envelope (vacuum tube) length
  G4double envelope2Size_Z = 2. * m;                        // Second envelope length
  // Scattering foil
  G4double ScFoil_X = 4. * cm;                              // Scattering foil width
  G4double ScFoil_Y = 4. * cm;                              // Scattering foil height
  G4double ScFoil_Z = runInput->GetFoilThickness();         // Scattering foil thickness (length)
  // Graphite colimator and air its air gap
  G4double GrColSize_XY = .5 * m;                           // Colimator body width and height
  G4double GrColSize_Z = 12. * mm;                          // Colimator thickness (length)
  // Air gap is cylindrical
  G4double gap_startingPhi = 0. * deg;                      // Starting angle for volume
  G4double gap_segmentAngle = 360.0 * deg;                  // Segment angle for volume
  G4double gap_innerRadius = 0 * mm;                        // Gap inner radius
  G4double gap_outerRadius = 6 * mm;                        // Gap diameter (particles passing through here)
  //G4double gap_outerRadius = .4 * m; //test modification!!!!!
  G4double gap_Z = 12.001 * mm;                             // Gap length
  //G4double gap_Z = 12.01 * mm; // test modification !!!!!
  // Silicon nitride window
  G4double Si3N4Size_XY = .5 * m;                           // Window height and width
  G4double Si3N4Size_Z = 1. * um;                           // Window thickness (length)
  // Mylar / PE (identified with MPE)
  G4double MPE_startingPhi = 0. * deg;                      // Starting angle for the volume
  G4double MPE_segmentAngle = 360.0 * deg;                  // Segment angle for the volume
  G4double MPE_innerRadius = 0. * cm;                       // MPE inner radius
  G4double MPE_outerRadius = 7. * mm;                       // Diameter of the MPE
  //G4double MPE_outerRadius = .4 * m; //test modification!!!!!
  G4double MPESize_Z = runInput -> GetMPEThickness();       // Mylar of PE foil thickness (length)

/*
   Placement parameters.
   We will define the G4ThreeVector placement parameter for every volume. 
   For the scoring and test volumes we will perform this step in a separate section
   along with the definition of their solid, logical and physical volumes.
   This section will be ordered acording to the hierarcy of the volumes, 
   each subsection containing its daughter volumes placement. The structure
   of the volumes hierarchy is:

   World:
         Envelope 1:
                    Scattering foil
                    Graphite colimator
         Envelope 2:
                    Silicon nitride window
                    MPE
                    for every NumberOfLayers -> Scoring volume
                                                Test volume
*/
  G4ThreeVector world_P = G4ThreeVector(0., 0., 0.);                // World volume placement at 0, 0, 0
  G4ThreeVector envelope1_P = G4ThreeVector(0, 0, -1.5 * m);      // First envelope - Mother: World
  G4ThreeVector ScFoil_P = G4ThreeVector(0, 0, 1. * m);      // Scattering foil - Mother: Envelope_1
  G4ThreeVector GrCol_P = G4ThreeVector(0, 0, 1.494 * m);     // Colimator - Mother: Envelope_1  
  G4ThreeVector envelope2_P = G4ThreeVector(0, 0, 1. * m);        // Second envelope - Mother: World
  G4ThreeVector Si3N4_P = G4ThreeVector(0, 0, -0.9999995 * m);// Silicon nitride window - Mother: Envelope_2
  //G4ThreeVector MPE_P = G4ThreeVector(0, 0, -0.9635 * m);   // Mylar or PE - Mother: Envelope_2 (proton and alpha configuration)
  G4ThreeVector MPE_P = G4ThreeVector(0, 0, -0.995 * m);   // Mylar or PE - Mother: Envelope_2 (heavy ions configuration)

  G4cout << " -------------------------------------------------------------------------------------- " << G4endl << G4endl;
  // World
  G4Box* World_Solid = new G4Box("World", 0.5 * worldSize_XY, 0.5 * worldSize_XY, 0.5 * worldSize_Z);
  G4LogicalVolume* World_Logical = new G4LogicalVolume(World_Solid, vacuum, "World_Vol");
  G4VPhysicalVolume* World_Physical = new G4PVPlacement(0, world_P, World_Logical, "World", 0, false, 0);
  G4cout<<"WORLD SIZE:"<<worldSize_XY/mm<<"mm X "<<worldSize_XY/mm<<"mm X "<<worldSize_Z/mm<<"mm; MATERIAL:"<<vacuum->GetName()<<"; PLACEMENT: x="<<world_P.getX()/m<<"m, y="<<world_P.getY()/m<<"m, z="<<world_P.getZ()/m<<"m"<<G4endl;
  // Envelope_1
  G4Box* Envelope1_Solid = new G4Box("Envelope_1", 0.5 * envelopeSize_XY, 0.5 * envelopeSize_XY, 0.5 * envelope1Size_Z);
  G4LogicalVolume* Envelope1_Logical = new G4LogicalVolume(Envelope1_Solid, vacuum, "Envelope_1_Vol");
  new G4PVPlacement(0, envelope1_P, Envelope1_Logical, "Envelope_1", World_Logical, false, 0, checkOverlaps);
  G4cout<<"Envelope_1 SIZE:"<<envelopeSize_XY/mm<<"mm X "<<envelopeSize_XY/mm<<"mm X "<<envelope1Size_Z/mm<<"mm; MATERIAL:"<<vacuum->GetName()<<"; PLACEMENT: x="<<envelope1_P.getX()/m<<"m, y="<<envelope1_P.getY()/m<<"m, z="<<envelope1_P.getZ()/m<<"m"<< G4endl;
  // Envelope_2
  G4Box* Envelope2_Solid = new G4Box("Envelope_2", 0.5 * envelopeSize_XY, 0.5 * envelopeSize_XY, 0.5 * envelope2Size_Z);
  G4LogicalVolume* Envelope2_Logical = new G4LogicalVolume(Envelope2_Solid, Air, "Envelope_2_Vol");
  new G4PVPlacement(0, envelope2_P, Envelope2_Logical, "Envelope_2", World_Logical, false, 0, checkOverlaps);
  G4cout<<"Envelope_2 SIZE: "<<envelopeSize_XY/mm<<"mm X "<<envelopeSize_XY/mm<<"mm X "<<envelope2Size_Z/mm<<"mm; MATERIAL: "<<Air->GetName()<<"; PLACEMENT: x="<<envelope2_P.getX()/m<<"m, y="<<envelope2_P.getY()/m<<"m, z="<<envelope2_P.getZ()/m<<"m"<< G4endl;
  //Scattering foil
  G4Box* ScFoil_Solid = new G4Box("ScFoil", 0.5 * ScFoil_X, 0.5 * ScFoil_Y, 0.5 * ScFoil_Z);
  G4LogicalVolume* ScFoil_Logical = new G4LogicalVolume(ScFoil_Solid, ScFoil_Material, "Scattering_foil_Vol");
  new G4PVPlacement(0, ScFoil_P, ScFoil_Logical, "Scattering_foil", Envelope1_Logical, false, 0, checkOverlaps);
  G4cout<<"Scattering_foil SIZE: "<<ScFoil_X/mm<<"mm X "<<ScFoil_Y/mm<<"mm X "<<ScFoil_Z/mm<<"mm; MATERIAL: "<<ScFoil_Material->GetName()<<"; PLACEMENT: x="<<ScFoil_P.getX()/m<<"m, y="<<ScFoil_P.getY()/m<<"m, z="<<ScFoil_P.getZ()/m<<"m"<< G4endl;
  // Graphite colimator and its air gap
  G4Box* GrCol_Solid = new G4Box("Graphite_colimator", 0.5 * GrColSize_XY, 0.5 * GrColSize_XY, 0.5 * GrColSize_Z);
  G4Tubs* gap_Solid = new G4Tubs("Gap", 0.5 * gap_innerRadius, 0.5 * gap_outerRadius, 0.5 * gap_Z, gap_startingPhi, gap_segmentAngle);
  G4ThreeVector shift = G4ThreeVector(0, 0, 0);
  G4VSolid* Colimator = new G4SubtractionSolid("Colimator", GrCol_Solid, gap_Solid, NULL, shift);
  G4LogicalVolume* Colimator_Logical = new G4LogicalVolume(Colimator, Graphite, "Colimator_Vol");
  new G4PVPlacement(0, GrCol_P, Colimator_Logical, "Colimator", Envelope1_Logical, false, 0, checkOverlaps);
  G4cout<<"Colimator SIZE: "<<GrColSize_XY/mm<<"mm X "<<GrColSize_XY/mm<<"mm X "<<GrColSize_Z/mm<<"mm with gap of "<<gap_outerRadius/mm<<"mm in diameter; MATERIAL: "<<Graphite->GetName()<< "; PLACEMENT: x="<<GrCol_P.getX()/m<<"m, y="<<GrCol_P.getY()/m<<"m, z="<<GrCol_P.getZ()/m<<"m"<< G4endl;
  // Silicon nitride window
  G4Box* Si3N4_Solid = new G4Box("Si3N4_window", 0.5 * Si3N4Size_XY, 0.5 * Si3N4Size_XY, 0.5 * Si3N4Size_Z);
  G4LogicalVolume* Si3N4_Logical = new G4LogicalVolume(Si3N4_Solid, Si3N4, "Si3N4_window_Vol");
  new G4PVPlacement(0, Si3N4_P, Si3N4_Logical, "Si3N4_window", Envelope2_Logical, false, 0, checkOverlaps);
  G4cout<<"Si3N4_window SIZE: "<<Si3N4Size_XY/mm<<"mm X "<<Si3N4Size_XY/mm<<"mm X "<<Si3N4Size_Z/mm<<"mm; MATERIAL: "<<Si3N4->GetName()<<"; PLACEMENT: x="<<Si3N4_P.getX()/m<<"m, y="<<Si3N4_P.getY()/m<<"m, z="<<Si3N4_P.getZ()/m<<"m"<< G4endl;
  // Mylar / PE (identified with MPE)
  G4Tubs* MPE_Solid = new G4Tubs("MPE", 0.5 * MPE_innerRadius, 0.5 * MPE_outerRadius, 0.5 * MPESize_Z, MPE_startingPhi, MPE_segmentAngle);
  G4LogicalVolume* MPE_Logical = new G4LogicalVolume(MPE_Solid, MPE_Material, "MPE_Vol");
  new G4PVPlacement(0, MPE_P, MPE_Logical, "MPE", Envelope2_Logical, false, 0, checkOverlaps);
  G4cout<<"MPE (air-water interface) SIZE: "<<MPE_outerRadius/mm<<"mm in diameter and thickness of "<<MPESize_Z/mm<<"mm; MATERIAL: "<<MPE_Material->GetName()<<"; PLACEMENT: x="<<MPE_P.getX()/m<<"m, y="<<MPE_P.getY()/m<<"m, z="<<MPE_P.getZ()/m<<"m"<< G4endl;
  G4cout << G4endl << G4endl << " -------------------------------------------------------------------------------------- " << G4endl << G4endl;

  /*
   Scoring and test volumes definition.
   Here we will define the scoring and test volumes. We will define the 
   solid, logical and physical volumes in a loop. The mother volume 
   MUST be Envelope_2.
*/

 G4int NumberOfLayers = runInput->GetNumberOfLayers();     // Number of water layer + test layer pairs
 G4double waterLayerThickness = 1.0 * um;                  // Water layer thickness
 G4double testLayerThickness = 0.1 * um;                   // Test layer thickness (for scoring fluence)
 G4double detector_innerRadius = 0. * mm;                  // Detector inner radius
 G4double detector_outerRadius = 6. * mm;                  // Diameter of the detector
 //G4double detector_outerRadius = 0.4 * m; // test modification!!!
 G4double detector_startingPhi = 0. * deg;                 // Starting angle for volume
 G4double detector_segmentAngle = 360.0 * deg;             // Segment angle for volume

 std::ostringstream solid_name(std::ostringstream::ate);
 G4LogicalVolume* Scoring_Logical[NumberOfLayers];
 G4LogicalVolume* Test_Logical[NumberOfLayers];

/*
   Useful variables
*/
 
 G4double interface_Z = MPE_P.getZ() / m;
 G4double MPESize_Z_update = (MPESize_Z / 2) / m;
 //G4double buffer = 1e-7 * m;
  
 G4double Test_ZStart = (interface_Z + MPESize_Z_update + ((testLayerThickness / 2) / m)) * 1000.0;
 G4double Scoring_ZStart = ((interface_Z + MPESize_Z_update + (testLayerThickness / m) + ((waterLayerThickness / 2) / m)) * 1000.0);
 //G4double Scoring_ZStart = (interface_Z + MPESize_Z_update + ((waterLayerThickness / 2) / m)) * 1000.0;

 G4ThreeVector LayerCenter = G4ThreeVector(0, 0, 0);

 G4Tubs* Scoring_Solid = new G4Tubs("Scoring_Volume", 0.5 * detector_innerRadius, 0.5 * detector_outerRadius, 0.5 * waterLayerThickness, detector_startingPhi, detector_segmentAngle);
 G4Tubs* Test_Solid = new G4Tubs("Test_Volume", 0.5 * detector_innerRadius, 0.5 * detector_outerRadius, 0.5 * testLayerThickness, detector_startingPhi, detector_segmentAngle);

 G4cout << "Adding " << NumberOfLayers << " scoring layer volumes made of " << Water->GetName() << ", starting from " << Scoring_ZStart / m << " m"
     << G4endl << "Scoring layer thickness: " << waterLayerThickness / m << " um"
     << G4endl << G4endl;

 G4cout << "Adding " << NumberOfLayers << " test/scoring layer volumes made of " << vacuum->GetName() << "/" << Water->GetName() << ", starting from " << Test_ZStart / m << "/" << Scoring_ZStart / m << " m"
	<< G4endl << "Test layer thickness: " << testLayerThickness / m << " um"
        << G4endl << "Scoring layer thickness: " << waterLayerThickness / m << " um"
        << G4endl << G4endl;

 for(G4int iDetector = 0; iDetector < NumberOfLayers; iDetector++) {
     solid_name.str("Test_Volume_V_"); solid_name << (iDetector);
     Test_Logical[iDetector] = new G4LogicalVolume(Test_Solid, vacuum, solid_name.str());
     LayerCenter.set(0, 0, Test_ZStart + (iDetector * waterLayerThickness) + (iDetector * testLayerThickness));
     new G4PVPlacement(0, LayerCenter, Test_Logical[iDetector], solid_name.str(), Envelope2_Logical, false, 0, checkOverlaps);
     G4cout << " Test layer volume placed at " << LayerCenter.getZ() / m << G4endl;

     solid_name.str("Scoring_Volume_V_"); solid_name << (iDetector);
     Scoring_Logical[iDetector] = new G4LogicalVolume(Scoring_Solid, Water, solid_name.str());
     LayerCenter.set(0, 0, Scoring_ZStart + (iDetector * waterLayerThickness) + (iDetector * testLayerThickness));
     new G4PVPlacement(0, LayerCenter, Scoring_Logical[iDetector], solid_name.str(), Envelope2_Logical, false, 0, checkOverlaps);
     G4cout << " Scoring layer volume placed at " << LayerCenter.getZ() / m << G4endl;
 }
  
  return World_Physical;   // Always return the physical World
}
