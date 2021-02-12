////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
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

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSNofStep.hh"
#include "G4PSCellFlux.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(), 
    fMagField(nullptr), 
    fLAbsor(nullptr),
    fLWorld(nullptr),
	fScoringVolume(),
	fMfdName("Target_MFD")
{
  // default parameter values
  fAbsorSizeZ = 150*um;
  fAbsorSizeXY = 6*mm;
  //fAbsorSizeXY = 8*mm; //lateral measurement

	//fTallyPosition[0] = G4ThreeVector(0, 0, -0.9899694*m);
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.9949685*m);
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.9899685*m);

 // fTallyPosition[0] = G4ThreeVector(0, 0, -0.971969*m);
  fTallyPosition[0] = G4ThreeVector(0, 0, -0.9719685*m);     //30, 6 or 4 um Mylar 
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.971966*m);        // 30 - 6 um
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.9719745*m);     //24 um Mylar
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.9799925*m);     //6 um Mylar He CONFIG
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.9949685*m);     //30 um Mylar 12C CONFIG
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.9949925*m);     //Heavy ions config. 2
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.971484*m);     //He Air
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.971469*m);     //He Air 2
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.994469*m);     //He Air 3
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.984469*m);     //He Air 4
  //fTallyPosition[0] = G4ThreeVector(0, 0, -0.971469*m); //MARKUS 1 mm
  
     
  //fTallyPosition[0] = G4ThreeVector(0, 0, 0.0280315*m);
  fTallySize = 1.*um;
  //fTallySize = 1.*mm;
  //fTallySize = 6.*um;
  fTallyNumber = 200;
  //G4double gapair = 500*um;
  
  
  for (G4int j=1; j<kMaxTally; j++) {
	fTallyPosition[j] = fTallyPosition[j-1]+G4ThreeVector(0.,0.,fTallySize);//+G4ThreeVector(0.,0.,gapair);
    //fTallyPosition[j] = fTallyPosition[j-1]-G4ThreeVector(0.,0.,fTallySize);//+G4ThreeVector(0.,0.,gapair);
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
  G4Element* C = new G4Element("Carbon", "C", z= 6., 12.01*g/mole);  

  //
  // define Materials.
  //
  G4double density, temperature, pressure;
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

  G4Material* Air = 
    new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  //density     = 1.e-5*g/cm3;
  density	  = universe_mean_density;
  //pressure    = 2.e-2*bar;
  pressure = 0.0001*pascal;
  //temperature = STP_Temperature;  // From PhysicalConstants.h .
  temperature = 2.73*kelvin;
  G4Material* vac = new G4Material( "TechVacuum", density, 1,
                           kStateGas, temperature, pressure );
  vac->AddMaterial( Air, 1. );

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material* vacuum = 
    new G4Material("Galactic",z= 1,a= 1.008*g/mole,density,
                   kStateGas,temperature,pressure);

  //default materials
  //G4Material* tallymat = nist->FindOrBuildMaterial("G4_WATER");
  //G4Material* videanu = nist->FindOrBuildMaterial("G4_Galactic");

    //Material definition
  G4Material* MylarMylar = new G4Material("Mylar", 1.39*g/cm3, 3);  
   MylarMylar->AddElement(O, natoms=2);                                    //Add 2 atoms of O  
   MylarMylar->AddElement(C, natoms=5);                                    //Add 5 atoms of C
   MylarMylar->AddElement(H, natoms=4);                                    //Add 4 atoms of H
  
  fAbsorMaterial = H2O;
  //fAbsorMaterial = MylarMylar;
  //fAbsorMaterial = Air;
  fWorldMaterial = vacuum;
}
  
G4VPhysicalVolume* DetectorConstruction::Construct()
{
	
	G4NistManager* nist = G4NistManager::Instance(); // Get NIST material manager
	
/////////////////////////////////////////////////////////////////////////////////////////////////
//World

  G4double fWorldSizeZ = 6.*m, fWorldSizeXY = .6*m;
  
  G4Box*
  sWorld = new G4Box("World",                                      //name
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
 
  G4Material* vac = nist->FindOrBuildMaterial("G4_Galactic");
  
  G4double air_temp = 298.15*kelvin;
  G4double air_pres = 101325*pascal;
  G4double density_vac = 1.55e-5*g/cm3;
  G4double pressure_vac = 0.0001*pascal;
  
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
  //G4Material* vac = new G4Material( "TechVacuum", density_vac, 1, kStateGas, air_temp, pressure_vac );
  //vac->AddMaterial( Air, 1. );


 
  //Envelope1 material
  G4double env_sizeXY = .5*m, env_sizeZA = 3.*m;
 
  G4Box* solidEnv1 =                                     
    new G4Box("Envelope1",                               
              0.5*env_sizeXY,
              0.5*env_sizeXY,                           
              0.5*env_sizeZA);                         
      
  G4LogicalVolume* logicEnv1 =                
    new G4LogicalVolume(solidEnv1,          
                        //env_mat1,
					    vac,
                        "Envelope1");         
               
  new G4PVPlacement(0,                       
                    G4ThreeVector(0, 0, -1.5*m),	//Placement coordinates (at [0,0,-1.5*m])
                    logicEnv1,                
                    "Envelope1",             
                    fLWorld,              
                    false,                 
                    0);                     

/////////////////////////////////////////////////////////////////////////////////////////////////
//Envelope2
  
  G4Material* env_mat2 = nist->FindOrBuildMaterial("G4_AIR"); //Envelope2 material
  G4double env_sizeZB = 2.*m;

   G4Box* solidEnv2 =                                     
     new G4Box("Envelope2",                               
               0.5*env_sizeXY,
               0.5*env_sizeXY,                           
               0.5*env_sizeZB);                           

    G4LogicalVolume* logicEnv2 =               
      new G4LogicalVolume(solidEnv2,           
                          env_mat2,
						  //Air,
                          "Envelope2");         

    new G4PVPlacement(0,                      
                     G4ThreeVector(0, 0, 1.*m),		//Placement coordinates (at [0,0,1.*m])
                     logicEnv2,               
                     "Envelope2",        
                     fLWorld,             
                      false,                  
                     0); 
                     
//////////////////////////////////////////////////////////////////////////////////////////////////////// 
//Gold layer
 
  //G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Au");	//Shape1 material
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Al");	//Shape1 material
  G4ThreeVector pos1 = G4ThreeVector(0, 0, -.5*m);              //Shape1 position
  
  //G4ThreeVector pos1 = G4ThreeVector(0, 0, -25000*um);                
  //G4ThreeVector pos1 = G4ThreeVector(0, 0, -1.3*m);                 


  G4double shape1_pX =  4.*cm, shape1_pY = 4.*cm;
  //G4double shape1_pZ = 6.7*um; //Shape1 size
  G4double shape1_pZ = 2.*um; //Shape1 size
  //G4double shape1_pZ = 1.*um; //Shape1 size

  
  G4VSolid* solidShape1 =          
    new G4Box ("Gold_Layer",       
	       0.5*shape1_pX,          
	       0.5*shape1_pY,          
	       0.5*shape1_pZ);         
                      
  G4LogicalVolume* logicShape1 =                        
    new G4LogicalVolume(solidShape1,         
                        shape1_mat,         
                        "Gold_Layer");      
              
  new G4PVPlacement(0,                       
                    pos1,                   
                    logicShape1,           
                    "Gold_Layer",           
                    logicEnv1,               
                    false,                  
                    0);                     

///////////////////////////////////////////////////////////////////////////////////////////////// 
//Graphite colimator + Gap

  G4Material* Graphite = new G4Material("Graphite", 6., 12.0107*g/mole, 2.2*g/cm3);	//Shape2 material
    
  G4ThreeVector pos2 = G4ThreeVector(0,0,1.494*m);	//Shape2 position
  G4ThreeVector pos3 = G4ThreeVector(0,0,1.494*m);  //Shape3 position (air gap)

  G4double shape2_pX = env_sizeXY, shape2_pY = env_sizeXY, shape2_pZ = 12.*mm;	//Shape 2 size
  
  G4double shape3_sa = 0.*deg, shape3_ea = 360.0*deg;
  G4double shape3_ir = 0*mm, shape3_or = 6*mm, shape3_h = 12.001*mm;	//Shape 3 size
  

  G4VSolid* solidShape2 =                    
    new G4Box("Graphite",                    
              0.5*shape2_pX,                
              0.5*shape2_pY,                
              0.5*shape2_pZ);                
 
  G4VSolid* solidShape3 =                   
    new G4Tubs("Gap",                      
              0.5*shape3_ir,               
              0.5*shape3_or,               
              0.5*shape3_h,
			  shape3_sa,
			  shape3_ea);              

  G4ThreeVector shift; //G4ThreeVector(0.,0.,0.)

  G4VSolid* Col  = new G4SubtractionSolid("Graphite_Colimator",solidShape2,solidShape3,NULL,shift); //Build gap

  G4LogicalVolume* logicCol =                   
    new G4LogicalVolume(Col,                    
                        Graphite,               
                        "Graphite_Colimator");  

  new G4PVPlacement(0,                       
                    pos2,                   
                    logicCol,               
                    "Graphite_Colimator",   
                    logicEnv1,             
                    false,                   
                    0);                    

/////////////////////////////////////////////////////////////////////////////////////////////////
//Si3N4 window
  
  //Material definition
  G4String name;            
  G4double density1, density2;             
  G4int natoms,nel;
  
  density1 = 3.44*g/cm3;                                                
  G4Element *  Si = new G4Element ("Silicon","Si",14., 28.0855*g/mole); 
  
  G4Material* Si3N4 = new G4Material("Si3N4", density1, nel=2);   
   Si3N4->AddElement (Si, natoms=3);                              //Add 3 atoms of Si
   Si3N4->AddElement (N , natoms=4);                              //Add 4 atoms of N
    
   G4ThreeVector pos4 = G4ThreeVector(0,0,-0.9999995*m); //Shape4 position


  G4double shape4_pX = env_sizeXY, shape4_pY = env_sizeXY, shape4_pZ = 1.*um; //Shape4 size  

  G4VSolid* solidShape4 =                   
    new G4Box("Si3N4_Window",                
              0.5*shape4_pX,                 
  	          0.5*shape4_pY,                 
              0.5*shape4_pZ);                
 
  G4LogicalVolume* logicWindow =                    
    new G4LogicalVolume(solidShape4,         
                        Si3N4,               
                        "Si3N4_Window");     

  new G4PVPlacement(0,                       
                    pos4,                   
                    logicWindow,            
                    "Si3N4_Window",          
                    logicEnv2,                
                    false,                  
                    0);                    

/////////////////////////////////////////////////////////////////////////////////////////////////
//Shape 6 - Mylar

  //Material definition
  density2 = 1.39*g/cm3;                                                
  G4Material* Mylar = new G4Material(name="Mylar", density2, nel=3);  
   Mylar->AddElement(O, natoms=2);                                    //Add 2 atoms of O  
   Mylar->AddElement(C, natoms=5);                                    //Add 5 atoms of C
   Mylar->AddElement(H, natoms=4);                                    //Add 4 atoms of H
   
   G4Material* poly = nist->FindOrBuildMaterial("G4_POLYETHYLENE");

   G4ThreeVector pos6 = G4ThreeVector(0, 0, -0.971984*m);  //pos 30 um
   G4double shape6_h = 6*um;
  
//   G4ThreeVector pos6 = G4ThreeVector(0, 0, -0.972002*m); G4double shape6_h = 6.*um;
//   G4ThreeVector pos6 = G4ThreeVector(0, 0, -0.972008*m); G4double shape6_h = 12.*um;
//   G4ThreeVector pos6 = G4ThreeVector(0, 0, -0.972014*m); G4double shape6_h = 18.*um;
//   G4ThreeVector pos6 = G4ThreeVector(0, 0, -0.972020*m); G4double shape6_h = 24.*um;
//   G4ThreeVector pos6 = G4ThreeVector(0, 0, -0.972026*m); G4double shape6_h = 30.*um;
//   G4ThreeVector pos6 = G4ThreeVector(0, 0, -0.972032*m); G4double shape6_h = 36.*um;
//   G4ThreeVector pos6 = G4ThreeVector(0, 0, -0.972038*m); G4double shape6_h = 42.*um;
//   G4ThreeVector pos6 = G4ThreeVector(0, 0, -0.972044*m); G4double shape6_h = 48.*um;
//   G4ThreeVector pos6 = G4ThreeVector(0, 0, -0.972050*m); G4double shape6_h = 54.*um;
//   G4ThreeVector pos6 = G4ThreeVector(0, 0, -0.972056*m); G4double shape6_h = 60.*um;

  
//  G4double shape6_h = 1*um;
  G4double shape6_ir =  0.*cm;
  G4double shape6_or = 7.*mm;

  
  G4double shape6_sa = 0.*deg, shape6_ea = 360.0*deg;	//Shape6 angles

  //G4double shape6_ir =  0.*cm, shape6_or = 1.*cm, shape6_h = 24.*um;  


  G4VSolid* solidShape6 =                       
      new G4Tubs("Mylar",                     
                 0.5*shape6_ir,                
      	         0.5*shape6_or,                
                 0.5*shape6_h,                 
				shape6_sa,                    
				shape6_ea);                    

  G4LogicalVolume* logicShape6 =                     
      new G4LogicalVolume(solidShape6,
						  //env_mat2,
                          Mylar,
						  //Air,
						  //poly,
						  //Si3N4,
                          "Mylar");            

  new G4PVPlacement(0,                      
                    pos6,                 
                    logicShape6,             
                    "Mylar",                
                    logicEnv2,                
                    false,                   
                    0);  					

/////////////////////////////////////////////////////////////////////////////////////////////////
//Shape 6.1 - PE

  // //G4ThreeVector pos6PE = G4ThreeVector(0, 0, -0.971972*m);  //pos 6 um
  // //G4ThreeVector pos6PE = G4ThreeVector(0, 0, -0.979996*m);  //pos 6 um He CONFIG
  // //G4ThreeVector pos6PE = G4ThreeVector(0, 0, -0.994984*m);  //pos 30 um 12C CONFIG
  // //G4ThreeVector pos6PE = G4ThreeVector(0, 0, -0.971984*m);  //pos 30 um
  // //G4ThreeVector pos6PE = G4ThreeVector(0, 0, -0.971987*m);  //pos 24 um
  // //G4ThreeVector pos6PE = G4ThreeVector(0, 0, -0.971971*m);  //pos 4 um
  // //G4ThreeVector pos6PE = G4ThreeVector(0, 0, -0.994996*m);  //pos 6 um HP config. 2
  // //G4ThreeVector pos6PE = G4ThreeVector(0, 0, -0.972014*m);  //pos 30 um - He Air  
  // //G4ThreeVector pos6PE = G4ThreeVector(0, 0, -0.994984*m);  //pos 30 um - He Air 2 
  // G4ThreeVector pos6PE = G4ThreeVector(0, 0, -0.971984*m); //pos 1 mm Markus

  //G4double shape6PE_ir =  0.*cm;
  //G4double shape6PE_or = 7.*mm;
  ////G4double shape6PE_or = 12.*mm;
  //
  //G4double shape6PE_h = 30.*um;  
  ////G4double shape6PE_h = 24.*um;
  ////G4double shape6PE_h = 6.*um;
  ////G4double shape6PE_h = 4.*um;
  //
  //G4double shape6PE_sa = 0.*deg, shape6PE_ea = 360.0*deg;	//Shape6 angles

  ////G4double shape6PE_ir =  0.*cm, shape6PE_or = 1.*cm, shape6PE_h = 24.*um;  


  //G4VSolid* solidShape6PE =                       
  //   new G4Tubs("PE",                     
  //              0.5*shape6PE_ir,                
  //   	        0.5*shape6PE_or,                
  //              0.5*shape6PE_h,                 
  //				shape6PE_sa,                    
		//		shape6PE_ea);                    

  //G4LogicalVolume* logicShape6PE =                     
	 //new G4LogicalVolume(solidShape6PE,
	 // 				     //env_mat2,
		//				 poly,
  //                       "PE");            

  //new G4PVPlacement(0,                      
  //                  pos6PE,                 
  //                  logicShape6PE,             
  //                  "PE",                
  //                  logicEnv2,                
  //                  false,                   
  //                  0);  

/////////////////////////////////////////////////////////////////////////////////////////////////  
// Absorber 
   
  G4double absorSizesh_sa = 0.*deg, absorSizesh_ea = 360*deg;   //Shape 7 angles
   
  //G4VSolid*
  
  //sAbsor = new G4Tubs("Absorber",                                 //name
   //                   0.5*shape6_ir,
	//	      0.5*fAbsorSizeXY,
	//	      0.5*fAbsorSizeZ,
	//	      absorSizesh_sa,
	//	      absorSizesh_ea); //dimensions
                                                                 
  //fLAbsor = new G4LogicalVolume(sAbsor,                   //shape
     //                           fAbsorMaterial,           //material
    //                            "Absorber");              //name
  
                              
 // new G4PVPlacement(0,                           //no rotation
   //                 G4ThreeVector(0, 0, -0.971469*m),  //30 um Mylar
	//				//G4ThreeVector(0, 0, -0.971474*m),
   //                 fLAbsor,                     //logical volume
   //                 "Absorber",                  //name
   //                 logicEnv2,                     //mother  volume
    //                false,                       //no boolean operation
    //                0);                          //copy number
  
  // Tallies (optional)
  
  
  if (fTallyNumber > 0) {
    for (G4int j=0; j<fTallyNumber; ++j) {
            
       G4VSolid* sTally = new G4Tubs("Tally",
				     0.5*shape6_ir,
				     0.5*fAbsorSizeXY,
				     0.5*fTallySize,
				     absorSizesh_sa, 
				     absorSizesh_ea);

	   fLTally[j] = new G4LogicalVolume(sTally,fAbsorMaterial,"Tally");
           
       new G4PVPlacement(0,                        //no rotation
                         fTallyPosition[j],        //position
                         fLTally[j],               //logical volume
                         "Tally",                  //name
                         //fLAbsor,                  //mother  volume
						 logicEnv2,
                         false,                    //no boolean operation
                         j+1);                     //copy number
       
		fTallyMass[j] = pi*fTallySize*fAbsorSizeXY*fAbsorSizeXY*(fAbsorMaterial->GetDensity());
	 
		  //G4MultiFunctionalDetector* fluxscore = new G4MultifunctionalDetector("fluxscore");

		  //fLTally[j]->SetSensitiveDetector(fluxscore);

		  //G4VPrimitiveScorer* TSF = new G4PSFlatSurfaceFlux("TotalSurfFlux");
		  //fluxscore->RegisterPrimitive(TSF);	  
    }               
  } 



  
  fScoringVolume = fLTally[178];
  //fScoringVolume = logicShape6;
  PrintParameters();
    
  //
  //always return the World volume
  //  
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

void DetectorConstruction::ConstructSDandField(G4int j)
{
	//------------------------------------------------//
	//            MultiFunctionalDetector             //
	//------------------------------------------------//
	// Define MultiFunctionalDetector with name.
	G4MultiFunctionalDetector* MFDet = new G4MultiFunctionalDetector(fMfdName);
	G4SDManager::GetSDMpointer()->AddNewDetector(MFDet);
	G4VPrimitiveScorer* pcur = new G4PSCellFlux("PassageCurrent");
	MFDet->RegisterPrimitive(pcur);

	// add scoring volumes
	if (j >= 0 && j < kMaxTally) {
		SetSensitiveDetector(fLTally[j], MFDet);
	}
	else {
		G4cout << "### DetectorConstruction::GetLOgicalTally WARNING: wrong tally "
			<< "number " << j << " is ignored" << G4endl;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
