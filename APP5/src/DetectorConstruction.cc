//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4SubtractionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4NistManager.hh"
#include "G4UnitsTable.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh" 

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(), 
    fMagField(nullptr), 
    fLAbsor(nullptr),
    fLWorld(nullptr),
	fScoringVolume()
{
  // default parameter values
 //fAbsorSizeZ = 30*mm;
  fAbsorSizeXY = 5.*mm;
  //fAbsorSizeXY = 8*mm; //lateral measurement
	

  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.0034815*m); //1 um
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.003457*m); //50 um
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.003402*m); //100 um
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.002952*m); //1 mm
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.003302*m); //300 um
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.003447*m); //10 
  fTallyPosition[0] = G4ThreeVector(0, 0, 0.0220005*m); //1 um Si config.

   
  //fTallyPosition[0] = G4ThreeVector(0, 0, 0.0280315*m);
  fTallySize = 1.*um;
  //fTallySize = 50.*um;
  //fTallySize = 100.*um;
  //fTallySize = 1.*mm;
  //fTallySize = 300.*um;
  //fTallySize = 10.*um;
  //fTallyNumber = 2000;
  //fTallyNumber = 100;
  //fTallyNumber = 40;
  fTallyNumber = 1;
  
  
  for (G4int j=1; j<kMaxTally; j++) {
	fTallyPosition[j] = fTallyPosition[j-1]+G4ThreeVector(0.,0.,fTallySize);
    fTallyMass[j]     = 0.;
    fLTally[j]        = nullptr; 
  }
    
  DefineMaterials();

  // create commands for interactive definition of the detector  
  fDetectorMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

void DetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nist = G4NistManager::Instance();
  //
  // define Elements
  //
  G4double z, a;

  G4Element* H = new G4Element("Hydrogen", "H", z= 1, a= 1.008*g/mole);
  G4Element* N = new G4Element("Nitrogen", "N", z= 7, a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z= 8, a= 16.00*g/mole);
  G4Element* C = new G4Element("Carbon", "C", z= 6., a= 12.01*g/mole);   

  //
  // define Materials.
  //
  G4double density, temperature, pressure, air_temp, air_pres, density_vac, pressure_vac;
  G4int    ncomponents, natoms;
  G4double fractionmass;
 
  G4Material* H2O = 
    new G4Material("Water", density= 1.0*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

  // In this line both G4_WATER and Water_1.05 will be constructed
  G4NistManager::Instance()->
    BuildMaterialWithNewDensity("Water_1.05","G4_WATER",1.05*g/cm3);

  air_temp = 298.15*kelvin;
  air_pres = 101325*pascal;
  G4Material* Air = 
    new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2, kStateGas, air_temp, air_pres);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  density_vac    = 1.55e-5*g/cm3;
  pressure_vac   = 0.13332236842*pascal;
  //temperature = STP_Temperature;  // From PhysicalConstants.h .
  G4Material* vac = new G4Material( "TechVacuum", density_vac, 1, kStateGas, air_temp, pressure_vac );
  vac->AddMaterial( Air, 1. );

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material* vacuum = 
    new G4Material("Galactic",z= 1,a= 1.008*g/mole,density,
                   kStateGas,temperature,pressure);

  G4Material* matplexiglassABS = new G4Material("Plexiglass",1.19*g/cm3,3);
  	matplexiglassABS->AddElement(H,0.08);
  	matplexiglassABS->AddElement(C,0.60);
  	matplexiglassABS->AddElement(O,0.32);

  //default materials
  //G4Material* tallymat = nist->FindOrBuildMaterial("G4_WATER");
  //G4Material* videanu = nist->FindOrBuildMaterial("G4_Galactic");
  
  //G4NistManager* nist = G4NistManager::Instance(); // Get NIST material manager
  G4Material * Si = nist->FindOrBuildMaterial("G4_Si");

  
  //fAbsorMaterial = H2O;
  //fAbsorMaterial = Air;
  //fAbsorMaterial = matplexiglassABS;
  fAbsorMaterial = Si;
  
  fWorldMaterial = vacuum;
}
  
G4VPhysicalVolume* DetectorConstruction::Construct()
{
	
	G4NistManager* nist = G4NistManager::Instance(); // Get NIST material manager
	
/////////////////////////////////////////////////////////////////////////////////////////////////
//World

  G4double fWorldSizeZ = 1500*mm, fWorldSizeXY = 0.06*m;
  
  G4Box*
  sWorld = new G4Box("World",            //name
                     0.5*fWorldSizeXY,
		     0.5*fWorldSizeXY,
		     0.5*fWorldSizeZ);   //dimensions

  fLWorld = new G4LogicalVolume(sWorld,                        //shape
                                fWorldMaterial,                //material
                                "World");                      //name

  G4VPhysicalVolume*                                   
  pWorld = new G4PVPlacement(0,                           //no rotation
                             G4ThreeVector(0.,0.,0.),     //at (0,0,0)
                             fLWorld,                     //logical volume
                             "World",                     //name
                             0,                           //mother  volume
                             false,                       //no boolean operation
                             0);                          //copy number
                          
 ///////////////////////////////////////////////////////////////////////////////////////////////
//Envelope1
  
  G4double air_temp = 298.15*kelvin;
  G4double air_pres = 101325*pascal;
  G4double density_vac = 1.55e-5*g/cm3;
  G4double pressure_vac = 1.3332236842*pascal;
  
  G4double z, a;

  G4Element* H = new G4Element("Hydrogen", "H", z= 1, a= 1.008*g/mole);
  G4Element* N = new G4Element("Nitrogen", "N", z= 7, a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z= 8, a= 16.00*g/mole);
  G4Element* C = new G4Element("Carbon", "C", z= 6., a= 12.01*g/mole); 
  
  G4Material* Air = 
    new G4Material("Air"  , 1.290*mg/cm3, 2, kStateGas, air_temp, air_pres);
  Air->AddElement(N, 0.7);
  Air->AddElement(O, 0.3);
  
  //temperature = STP_Temperature;  // From PhysicalConstants.h .
  G4Material* vac = new G4Material( "TechVacuum", density_vac, 1, kStateGas, air_temp, pressure_vac );
  vac->AddMaterial( Air, 1. );
 
  G4Material* env_mat1 = nist->FindOrBuildMaterial("G4_AIR");
  
  G4double env_sizeXY = 50*mm, env_sizeZA = 120*mm;
 
  G4Box* solidEnv1 =                                     
    new G4Box("Envelope1",                               
              0.5*env_sizeXY,
              0.5*env_sizeXY,                           
              0.5*env_sizeZA);                         
      
  G4LogicalVolume* logicEnv1 =                
    new G4LogicalVolume(solidEnv1,  
						//env_mat1,
                        Air,             
                        "Envelope1");         
               
  new G4PVPlacement(0,                       
                    G4ThreeVector(0, 0, -0.71*m),
                    logicEnv1,                
                    "Envelope1",             
                    fLWorld,              
                    false,                 
                    0);                     

/////////////////////////////////////////////////////////////////////////////////////////////////
//W_LAYER_1
 
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_W");
  G4ThreeVector pos1 = G4ThreeVector(0, 0, 0.005*m);
               
  G4double shape1_pX =  50.*mm, shape1_pY = 50.*mm;
  G4double shape1_pZ = 0.025*mm;

  
  G4VSolid* solidShape1 =          
    new G4Box ("W_Layer_1",       
	       0.5*shape1_pX,          
	       0.5*shape1_pY,          
	       0.5*shape1_pZ);         
                      
  G4LogicalVolume* logicShape1 =                        
    new G4LogicalVolume(solidShape1,         
                        shape1_mat,         
                        "W_Layer_1");      
              
  new G4PVPlacement(0,                       
                    pos1,                   
                    logicShape1,           
                    "W_Layer_1",           
                    logicEnv1,               
                    false,                  
                    0); 

/////////////////////////////////////////////////////////////////////////////////////////////////
//Au_DEP
 
//  G4Material* shape1_1_mat = nist->FindOrBuildMaterial("G4_Au");
//  G4ThreeVector pos1_1 = G4ThreeVector(0, 0, 0.005015*m);
               
//  G4double shape1_1_pZ = 0.005*mm;

//  G4VSolid* solidShape1_1 =          
//    new G4Box ("Au_DEP",       
//	       0.5*shape1_pX,          
//	       0.5*shape1_pY,          
//	       0.5*shape1_1_pZ);         
                      
//  G4LogicalVolume* logicShape1_1 =                        
//    new G4LogicalVolume(solidShape1_1,         
//                        shape1_1_mat,         
//                        "Au_DEP");      
              
//  new G4PVPlacement(0,                       
//                    pos1_1,                   
//                    logicShape1_1,           
//                    "Au_DEP",           
//                    logicEnv1,               
//                    false,                  
//                    0); 
                    
///////////////////////////////////////////////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//Envelope2
  
  G4Material* env_mat2 = nist->FindOrBuildMaterial("G4_Galactic");
  
  
  G4double env_sizeZB = 1300.*mm;

   G4Box* solidEnv2 =                                     
     new G4Box("Envelope2",                               
               0.5*env_sizeXY,
               0.5*env_sizeXY,                           
               0.5*env_sizeZB);                           

    G4LogicalVolume* logicEnv2 =               
      new G4LogicalVolume(solidEnv2,
						  //env_mat2,
                          vac,             
                          "Envelope2");         

    new G4PVPlacement(0,                      
                     G4ThreeVector(0, 0, 0.*m),
                     logicEnv2,               
                     "Envelope2",        
                     fLWorld,             
                      false,                  
                     0); 

/////////////////////////////////////////////////////////////////////////////////////////////////
//Al_LAYER_1
 
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_Al");
  G4ThreeVector pos2 = G4ThreeVector(0, 0, -0.64998*m);
               
  G4double shape2_pX =  50.*mm;
  G4double shape2_pY = 50.*mm;
  G4double shape2_pZ = 0.04*mm; //Shape1 size

  
  G4VSolid* solidShape2 =          
    new G4Box ("Al_Layer_1",       
	       0.5*shape2_pX,          
	       0.5*shape2_pY,          
	       0.5*shape2_pZ);         
                      
  G4LogicalVolume* logicShape2 =                        
    new G4LogicalVolume(solidShape2,         
                        shape2_mat,         
                        "Al_Layer_1");      
              
  new G4PVPlacement(0,                       
                    pos2,                   
                    logicShape2,           
                    "Al_Layer_1",           
                    logicEnv2,               
                    false,                  
                    0);                     
/////////////////////////////////////////////////////////////////////////////////////////////////
//Al_Col_1 + Al_Col_window_1
    
  G4ThreeVector pos3 = G4ThreeVector(0,0,-0.63996*m);
  G4ThreeVector pos31 = G4ThreeVector(0,0,-0.63996*m);

  G4double shape3_pX = 50*mm;
  G4double shape3_pY = 50*mm;
  G4double shape3_pZ = 20.*mm;
  
  G4double shape31_sa = 0.*deg, shape31_ea = 360.0*deg;
  G4double shape31_ir = 0*mm, shape31_or = 5*mm, shape31_h = 20.001*mm;
  

  G4VSolid* solidShape3 =                    
    new G4Box("Al_Col_1",                    
              0.5*shape3_pX,                
              0.5*shape3_pY,                
              0.5*shape3_pZ);                
 
  G4VSolid* solidShape31 =                   
    new G4Tubs("Al_Col_window_1",                      
              0.5*shape31_ir,               
              0.5*shape31_or,               
              0.5*shape31_h,
		  shape31_sa,
	          shape31_ea);              

  G4ThreeVector shift1;
  G4VSolid* solidShape331  = new G4SubtractionSolid("Al_Col_1",solidShape3,solidShape31,NULL,shift1);

  G4LogicalVolume* logicShape331 =                   
    new G4LogicalVolume(solidShape331,                    
                        shape2_mat,               
                        "Al_Col_1");  

  new G4PVPlacement(0,                       
                    pos3,                   
                    logicShape331,               
                    "Al_Col_1",   
                    logicEnv2,             
                    false,                   
                    0); 
/////////////////////////////////////////////////////////////////////////////////////////////////
//W_LAYER_2

  G4ThreeVector pos4 = G4ThreeVector(0, 0, -0.6299475*m);
               
  G4double shape4_pX =  50.*mm;
  G4double shape4_pY = 50.*mm;
  G4double shape4_pZ = 0.025*mm;

  
  G4VSolid* solidShape4 =          
    new G4Box ("W_Layer_2",       
	       0.5*shape4_pX,          
	       0.5*shape4_pY,          
	       0.5*shape4_pZ);         
                      
  G4LogicalVolume* logicShape4 =                        
    new G4LogicalVolume(solidShape4, 
                        env_mat2,
						//vac,
						//shape1_mat,         
                        "W_Layer_2");      
              
  new G4PVPlacement(0,                       
                    pos4,                   
                    logicShape4,           
                    "W_Layer_2",           
                    logicEnv2,               
                    false,                  
                    0);  
/////////////////////////////////////////////////////////////////////////////////////////////////
//Al_Col_2 + Al_Col_window_2
    
  G4ThreeVector pos5 = G4ThreeVector(0,0,0.64496*m);
  G4ThreeVector pos51 = G4ThreeVector(0,0,0.64496*m);

  G4double shape5_pX = 50*mm;
  G4double shape5_pY = 50*mm;
  G4double shape5_pZ = 10.*mm;
  
  G4double shape51_sa = 0.*deg, shape51_ea = 360.0*deg;
  G4double shape51_ir = 0*mm, shape51_or = 30*mm, shape51_h = 10.001*mm;
  

  G4VSolid* solidShape5 =                    
    new G4Box("Al_Col_2",                    
              0.5*shape5_pX,                
              0.5*shape5_pY,                
              0.5*shape5_pZ);                
 
  G4VSolid* solidShape51 =                   
    new G4Tubs("Al_Col_window_2",                      
              0.5*shape51_ir,               
              0.5*shape51_or,               
              0.5*shape51_h,
		  shape51_sa,
	          shape51_ea);              

  G4ThreeVector shift2;
  G4VSolid* solidShape551  = new G4SubtractionSolid("Al_Col_2",solidShape5,solidShape51,NULL,shift2);

  G4LogicalVolume* logicShape551 =                   
    new G4LogicalVolume(solidShape551,                    
                        shape2_mat,               
                        "Al_Col_2");  

  new G4PVPlacement(0,                       
                    pos5,                   
                    logicShape551,               
                    "Al_Col_1",   
                    logicEnv2,             
                    false,                   
                    0); 
/////////////////////////////////////////////////////////////////////////////////////////////////
//Al_LAYER_2
 
  G4ThreeVector pos6 = G4ThreeVector(0, 0, 0.64998*m);
               
  G4double shape6_pX =  50.*mm;
  G4double shape6_pY = 50.*mm;
  G4double shape6_pZ = 0.04*mm; //Shape1 size

  
  G4VSolid* solidShape6 =          
    new G4Box ("Al_Layer_2",       
	       0.5*shape6_pX,          
	       0.5*shape6_pY,          
	       0.5*shape6_pZ);         
                      
  G4LogicalVolume* logicShape6 =                        
    new G4LogicalVolume(solidShape6,         
                        shape2_mat,         
                        "Al_Layer_2");      
              
  new G4PVPlacement(0,                       
                    pos6,                   
                    logicShape6,           
                    "Al_Layer_2",           
                    logicEnv2,               
                    false,                  
                    0); 
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//Envelope3
  
  G4double env_sizeZC = 70.*mm;

   G4Box* solidEnv3 =                                     
     new G4Box("Envelope3",                               
               0.5*env_sizeXY,
               0.5*env_sizeXY,                           
               0.5*env_sizeZC);                           

    G4LogicalVolume* logicEnv3 =               
      new G4LogicalVolume(solidEnv3, 
						  //env_mat1,
                          Air,             
                          "Envelope3");         

    new G4PVPlacement(0,                      
                     G4ThreeVector(0, 0, 0.685*m),		//Placement coordinates (at [0,0,1.*m])
                     logicEnv3,               
                     "Envelope3",        
                     fLWorld,             
                      false,                  
                     0);                      
//////////////////////////////////////////////////////////////////////////////////////////////////////// 
//PLEX_Col + PLEX_Col_window

  G4Material* matplexiglass = new G4Material("Plexiglass",1.19*g/cm3,3);
  	matplexiglass->AddElement(H,0.08);
  	matplexiglass->AddElement(C,0.60);
  	matplexiglass->AddElement(O,0.32);
    
  G4ThreeVector pos7 = G4ThreeVector(0,0,0.-0.026*m);
  G4ThreeVector pos71 = G4ThreeVector(0,0,-0.026*m);

  G4double shape7_pX = 50*mm;
  G4double shape7_pY = 50*mm;
  G4double shape7_pZ = 10.*mm;
  
  G4double shape71_sa = 0.*deg, shape71_ea = 360.0*deg;
  G4double shape71_ir = 0*mm, shape71_or = 21*mm, shape71_h = 10.001*mm;
  

  G4VSolid* solidShape7 =                    
    new G4Box("PLEX_Col",                    
              0.5*shape7_pX,                
              0.5*shape7_pY,                
              0.5*shape7_pZ);                
 
  G4VSolid* solidShape71 =                   
    new G4Tubs("PLEX_Col_window",                      
              0.5*shape71_ir,               
              0.5*shape71_or,               
              0.5*shape71_h,
		  shape71_sa,
	          shape71_ea);              

  G4ThreeVector shift3;
  G4VSolid* solidShape771  = new G4SubtractionSolid("PLEX_Col",solidShape7,solidShape71,NULL,shift3);

  G4LogicalVolume* logicShape771 =                   
    new G4LogicalVolume(solidShape771,                    
                        matplexiglass,               
                        "PLEX_Col");  

  new G4PVPlacement(0,                       
                    pos7,                   
                    logicShape771,               
                    "PLEX_Col",   
                    logicEnv3,             
                    false,                   
                    0);
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Al_Col_3 + Al_Col_window_3
    
  G4ThreeVector pos8 = G4ThreeVector(0,0,-0.0195*m);
  G4ThreeVector pos81 = G4ThreeVector(0,0,-0.0195*m);

  G4double shape8_pX = 50*mm;
  G4double shape8_pY = 50*mm;
  G4double shape8_pZ = 3.*mm;
  
  G4double shape81_sa = 0.*deg, shape81_ea = 360.0*deg;
  G4double shape81_ir = 0*mm, shape81_or = 21*mm, shape81_h = 3.001*mm;
  

  G4VSolid* solidShape8 =                    
    new G4Box("Al_Col_3",                    
              0.5*shape8_pX,                
              0.5*shape8_pY,                
              0.5*shape8_pZ);                
 
  G4VSolid* solidShape81 =                   
    new G4Tubs("Al_Col_window_3",                      
              0.5*shape81_ir,               
              0.5*shape81_or,               
              0.5*shape81_h,
		  shape81_sa,
	          shape81_ea);              

  G4ThreeVector shift4;
  G4VSolid* solidShape881  = new G4SubtractionSolid("Al_Col_3",solidShape8,solidShape81,NULL,shift4);

  G4LogicalVolume* logicShape881 =                   
    new G4LogicalVolume(solidShape881,                    
                        shape2_mat,               
                        "Al_Col_3");  

  new G4PVPlacement(0,                       
                    pos8,                   
                    logicShape881,               
                    "Al_Col_3",   
                    logicEnv3,             
                    false,                   
                    0);                 
////////////////////////////////////////////////////////////////////////////////////////////////////////
//PLEX
  
  //G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.0175*m); //1 mm
  //G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.013*m); // 10 mm
  //G4ThreeVector pos9 = G4ThreeVector (0, 0, -0.017*m); // 2 mm
  //G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.01725*m); //1.5 mm
  
  G4ThreeVector pos9 = G4ThreeVector(0, 0, 0.0005*m); //1 mm SI CONFIG.
  
  //G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.017475*m); // 1 mm + 50 um
  //G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.01745*m); // 1 mm + 100 um
  //G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.017425*m); // 1 mm + 150 um
  //G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.0174*m); // 1 mm + 200 um
  //G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.017375*m); // 1 mm + 250 um
  //G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.01735*m); // 1 mm + 300 um
  //G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.017325*m); // 1 mm + 350 um
  //G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.0173*m); // 1 mm + 400 um
  //G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.017275*m); // 1 mm + 450 um
  
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005232*m); G4double shape9_pZ = 500.*um; 	// 500 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005357*m); G4double shape9_pZ = 750.*um; 	// 750 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005382*m); G4double shape9_pZ = 800.*um;	// 800 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005407*m); G4double shape9_pZ = 850.*um;	// 850 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005432*m); G4double shape9_pZ = 900.*um;	// 900 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005457*m); G4double shape9_pZ = 950.*um;	// 950 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005482*m); G4double shape9_pZ = 1000.*um;	//1000 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005507*m); G4double shape9_pZ = 1050.*um;	//1050 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005532*m); G4double shape9_pZ = 1100.*um;	//1100 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005557*m); G4double shape9_pZ = 1150.*um;	//1150 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005582*m); G4double shape9_pZ = 1200.*um;	//1200 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005607*m); G4double shape9_pZ = 1250.*um;	//1250 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005632*m); G4double shape9_pZ = 1300.*um;	//1300 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005657*m); G4double shape9_pZ = 1350.*um;	//1350 um
//  G4ThreeVector pos9 = G4ThreeVector(0, 0, -0.005682*m); G4double shape9_pZ = 1400.*um;	//1400 um
               
			   
  G4double shape9_pX =  50.*mm;
  G4double shape9_pY = 50.*mm;
  G4double shape9_pZ = 1*mm;
  
  //G4double shape9_pZ = 10*mm;
  //G4double shape9_pZ = 2*mm;
  //G4double shape9_pZ = 1.5*mm;
  
  //G4double shape9_pZ = 1.05*mm;
  //G4double shape9_pZ = 1.1*mm;
  //G4double shape9_pZ = 1.15*mm;
  //G4double shape9_pZ = 1.2*mm;
  //G4double shape9_pZ = 1.25*mm;
  //G4double shape9_pZ = 1.3*mm;
  //G4double shape9_pZ = 1.35*mm;
  //G4double shape9_pZ = 1.4*mm;
  //G4double shape9_pZ = 1.45*mm;

  G4VSolid* solidShape9 =          
    new G4Box ("PLEX",       
	       0.5*shape9_pX,          
	       0.5*shape9_pY,          
	       0.5*shape9_pZ);         
                      
  G4LogicalVolume* logicShape9 =                        
    new G4LogicalVolume(solidShape9, 
						//env_mat1,
						//Air,	
                        matplexiglass,         
                        "PLEX");      
              
  new G4PVPlacement(0,                       
                    pos9,                   
                    logicShape9,           
                    "PLEX",           
                    logicEnv3,               
                    false,                  
                    0);  
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Mylar  
/*
  G4Material* Mylar = new G4Material("Mylar", 1.39*g/cm3, 3);  
   Mylar->AddElement(O, 2);                                     
   Mylar->AddElement(C, 5);                                   
   Mylar->AddElement(H, 4);                                    
   
//  G4ThreeVector pos10 = G4ThreeVector(0, 0, -0.006138*m); G4double shape10_pZ = 12.*um;	//1150+12 um
//  G4ThreeVector pos10 = G4ThreeVector(0, 0, -0.006144*m); G4double shape10_pZ = 24.*um;	//1150+24 um
//  G4ThreeVector pos10 = G4ThreeVector(0, 0, -0.006150*m); G4double shape10_pZ = 36.*um;	//1150+36 um

//  G4ThreeVector pos10 = G4ThreeVector(0, 0, -0.006188*m); G4double shape10_pZ = 12.*um;	//1200+12 um
  G4ThreeVector pos10 = G4ThreeVector(0, 0, -0.006194*m); G4double shape10_pZ = 24.*um;	//1200+24 um

  G4double shape10_ir =  0.*cm;
  G4double shape10_or = 5.*mm;
  G4double shape10_sa = 0.*deg;
  G4double shape10_ea = 360.0*deg;


  G4VSolid* solidShape10 =                       
      new G4Tubs("Mylar",                     
                 0.5*shape10_ir,                
      	         0.5*shape10_or,                
                 0.5*shape10_pZ,                 
		     	 shape10_sa,                    
		     	 shape10_ea);                    

  G4LogicalVolume* logicShape10 =                     
      new G4LogicalVolume(solidShape10,
                          Mylar,
						  //Air,
                          "Mylar");            

  new G4PVPlacement(0,                      
                    pos10,                 
                    logicShape10,             
                    "Mylar",                
                    logicEnv3,                
                    false,                   
                    0);  		
*/			
/////////////////////////////////////////////////////////////////////////////////////////////////
//Caseta

  G4Material* Polystyrene = new G4Material("Polystyrene", 1.032*g/cm3, 2);
  	Polystyrene->AddElement(C, 19);
  	Polystyrene->AddElement(H, 21);

  //G4ThreeVector pos11 = G4ThreeVector(0, 0, -0.004232*m);
  G4ThreeVector pos11 = G4ThreeVector(0, 0, 0.0115*m);
   
  G4double shape11_ir =  0.*cm;
  G4double shape11_or = 21.*mm;
  G4double shape11_h = 1.5*mm;  
  G4double shape11_sa = 0.*deg;
  G4double shape11_ea = 360.0*deg;

  G4VSolid* solidShape11 =                       
      new G4Tubs("Caseta",                     
                 0.5*shape11_ir,                
      	         0.5*shape11_or,                
                 0.5*shape11_h,                 
		     shape11_sa,                    
		     shape11_ea);                    

  G4LogicalVolume* logicShape11 =                     
      new G4LogicalVolume(solidShape11,
						  //env_mat1,
                          Polystyrene,
                          "Caseta");            

  new G4PVPlacement(0,                      
                    pos11,                 
                    logicShape11,             
                    "Caseta",                
                    logicEnv3,                
                    false,                   
                    0); 
					
/////////////////////////////////////////////////////////////////////////////////////////////////
//PTE
/*
  G4Material* PTE = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4ThreeVector pos12 = G4ThreeVector(0, 0, -0.003467*m);
   
  G4double shape12_ir =  0.*cm;
  G4double shape12_or = 21.*mm;
  G4double shape12_h = 0.03*mm;  
  G4double shape12_sa = 0.*deg;
  G4double shape12_ea = 360.0*deg;

  G4VSolid* solidShape12 =                       
      new G4Tubs("PTE",                     
                 0.5*shape12_ir,                
      	         0.5*shape12_or,                
                 0.5*shape12_h,                 
				 shape12_sa,                    
				 shape12_ea);                    

  G4LogicalVolume* logicShape12 =                     
      new G4LogicalVolume(solidShape12,
						  PTE,
                          "PTE");            

  new G4PVPlacement(0,                      
                    pos12,                 
                    logicShape12,             
                    "PTE",                
                    logicEnv3,                
                    false,                   
                    0); 
*/
/////////////////////////////////////////////////////////////////////////////////////////////////  
// Absorber 
   
  G4double absorSizesh_sa = 0.*deg, absorSizesh_ea = 360*deg;   //Shape 7 angles
  
  if (fTallyNumber > 0) {
    for (G4int j=0; j<fTallyNumber; ++j) {
            
       G4VSolid* sTally = new G4Tubs("Tally",
				     0.5*shape11_ir,
				     0.5*fAbsorSizeXY,
				     0.5*fTallySize,
				     absorSizesh_sa, 
				     absorSizesh_ea);

  fLTally[j] = new G4LogicalVolume(sTally,
								   fAbsorMaterial,
								   "Tally");
           
       new G4PVPlacement(0,                        
                         fTallyPosition[j],       
                         fLTally[j],               
                         "Tally",                  
						 logicEnv3,
                         false,                    
                         j+1);                     
       
      fTallyMass[j] = pi*fTallySize*fAbsorSizeXY*fAbsorSizeXY*(fAbsorMaterial->GetDensity());
    }               
  } 
  
  fScoringVolume = fLTally[0];
  //fScoringVolume = logicShape11;
  PrintParameters();

  return pWorld;
}

void DetectorConstruction::PrintParameters() const
{
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Absorber is " << G4BestUnit(fAbsorSizeZ,"Length")
         << " of " << fAbsorMaterial->GetName() << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
  
  if (fTallyNumber > 0) {
    G4cout << "---> There are " << fTallyNumber << " tallies : " << G4endl;    
    for (G4int j=0; j<fTallyNumber; ++j) {
      G4cout << "fTally " << j << ": "
             << fAbsorMaterial->GetName()
             << ",  mass = " << G4BestUnit(fTallyMass[j],"Mass")
             << " size = "   << G4BestUnit(fTallySize,"Length")           
             << " position = " << G4BestUnit(fTallyPosition[j],"Length")
             << G4endl;
    }                 
    G4cout << "\n---------------------------------------------------------\n";
  }  
}


void DetectorConstruction::SetSizeZ(G4double value)
{
  fAbsorSizeZ = value; 
}
  

void DetectorConstruction::SetSizeXY(G4double value)
{
  fAbsorSizeXY = value; 
}  


void DetectorConstruction::SetMaterial(const G4String& materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  if (pttoMaterial && pttoMaterial != fAbsorMaterial) {
    // change target material everywhere
    fAbsorMaterial = pttoMaterial;
    for (G4int j=0; j<fTallyNumber; ++j) {
      if(fLTally[j]) { 
        fLTally[j]->SetMaterial(pttoMaterial); 
        fTallyMass[j] = pi*fTallySize*fAbsorSizeXY*fAbsorSizeXY*(pttoMaterial->GetDensity());
      }
    } 
    if(fLAbsor) {
      fLAbsor->SetMaterial(fAbsorMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}


void DetectorConstruction::SetWorldMaterial(const G4String& materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  if (pttoMaterial && pttoMaterial != fWorldMaterial) {
    fWorldMaterial = pttoMaterial;
    if(fLWorld) {
      fLWorld->SetMaterial(fAbsorMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}


void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr 
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
  if (fMagField) delete fMagField;        //delete the existing magn field  

  if (fieldValue!=0.)                        // create a new one if non nul
    {
      fMagField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
      fieldMgr->SetDetectorField(fMagField);
      fieldMgr->CreateChordFinder(fMagField);
    }
   else
    {
      fMagField = nullptr;
      fieldMgr->SetDetectorField(fMagField);
    }
}

void DetectorConstruction::SetTallyNumber(G4int value)
{
  if(value >= 0 && value <kMaxTally) {
    fTallyNumber = value;
  } else {
    G4cout << "### DetectorConstruction::SetTallyNumber WARNING: wrong tally "
           << "number " << value << " is ignored" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTallySize(G4int j, G4double value)
{
  if(j >= 0 && j < kMaxTally) {
    fTallySize = value;
  } else {
    G4cout << "### DetectorConstruction::SetTallyNumber WARNING: wrong tally size. Pick a different tally size.";
  } 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTallyPosition(G4int j, const G4ThreeVector& value)
{
  if(j >= 0 && j < kMaxTally) {
    fTallyPosition[j] = value; 
  } else {
    G4cout << "### DetectorConstruction::SetTallyPosition WARNING: wrong tally "
           << "number " << j << " is ignored" << G4endl;
  } 
}  

G4double DetectorConstruction::GetTallyMass(G4int j) const
{
  if(j >= 0 && j < kMaxTally) {
    return fTallyMass[j];
  } else {
    G4cout << "### DetectorConstruction::GetTallyMass WARNING: wrong tally "
           << "number " << j << " is ignored" << G4endl;
    return 0.0;
  } 
}

const G4LogicalVolume* DetectorConstruction::GetLogicalTally(G4int j) const 
{
  if(j >= 0 && j < kMaxTally) {
    return fLTally[j];
  } else {
    G4cout << "### DetectorConstruction::GetLOgicalTally WARNING: wrong tally "
           << "number " << j << " is ignored" << G4endl;
    return nullptr;
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
