//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fTargetMater(0), fLogicTarget(0),
 fDetectorMater(0), fLogicDetector(0),
 fGapMater(0), fLogicGap(0),
 fShieldMater(0), fLogicShield(0),
 fEndCapNorthDetectorMater(0), fLogicEndCapNorthDetector(0),
 fEndCapNorthGapMater(0), fLogicEndCapNorthGap(0),
 fEndCapNorthShieldMater(0), fLogicEndCapNorthShield(0),
 fEndCapSouthDetectorMater(0), fLogicEndCapSouthDetector(0),
 fEndCapSouthGapMater(0), fLogicEndCapSouthGap(0),
 fEndCapSouthShieldMater(0), fLogicEndCapSouthShield(0),
 fWorldMater(0), fPhysiWorld(0),
 fDetectorMessenger(0)
{
  // Figures given by Chris:
  // Target: radius 12m | length 24m
  // Detector: thickness 1mm-1cm | length 24m
  // Gap: thickness 10cm-5m | length 24m
  // Shield: Barlow's formula + safety factor (1.5 - 10) | length 24m
  fTargetLength      = 24*m;
  fTargetRadius      = 12*m;
  fDetectorLength    = 24*m;
  fDetectorThickness = 5*mm;
  fGapLength = 24*m;
  fGapThickness = 255*cm;
  fShieldLength = 24*m;
  fShieldThickness = 25*cm;

  fEndCapNorthDetectorLength = fDetectorThickness;
  fEndCapNorthDetectorRadius = fTargetRadius + fDetectorThickness + fGapThickness + fShieldThickness;
  fEndCapNorthGapLength = fGapThickness;
  fEndCapNorthGapRadius = fTargetRadius + fDetectorThickness + fGapThickness + fShieldThickness;
  fEndCapNorthShieldLength = fShieldThickness;
  fEndCapNorthShieldRadius = fTargetRadius + fDetectorThickness + fGapThickness + fShieldThickness;

  fEndCapSouthDetectorLength = fDetectorThickness;
  fEndCapSouthDetectorRadius = fTargetRadius + fDetectorThickness + fGapThickness + fShieldThickness;
  fEndCapSouthGapLength = fGapThickness;
  fEndCapSouthGapRadius = fTargetRadius + fDetectorThickness + fGapThickness + fShieldThickness;
  fEndCapSouthShieldLength = fShieldThickness;
  fEndCapSouthShieldRadius = fTargetRadius + fDetectorThickness + fGapThickness + fShieldThickness;

  fWorldLength = (std::max(std::max(fTargetLength,fDetectorLength), std::max(fGapLength,fShieldLength))) + (2*(fEndCapNorthDetectorLength+fEndCapNorthGapLength+fEndCapNorthShieldLength));
  fWorldRadius = fTargetRadius + fDetectorThickness + fGapThickness + fShieldThickness + (1.0*m);

  DefineMaterials();

  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // target density
  targetPressure=15;
  densityXe = ((targetPressure*molarMassXe)/(0.0821*275))/1E3;
  densityXe136 = 80.3; //kg/m3

  // build materials
  //
  mXenon = new G4Material("Xenon", 54, molarMassXe*g/mole, densityXe*g/cm3, kStateGas, 275.*kelvin, targetPressure=15*atmosphere);
  mXenon136 = new G4Material("Xenon136", 54, molarMassXe136*g/mole, densityXe136*kg/m3, kStateGas, 275.*kelvin, targetPressure=15*atmosphere);

  fDetectorMater =
  new G4Material("Copper", 29, 63.546*g/mole, 8.96*g/cm3);


  G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",   "O", 8, 15.9994*g/mole);
  //
  G4int ncomponents; G4double fractionmass;
  G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2,
                      kStateGas, 293.*kelvin, 1.*atmosphere);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);
  //
  fWorldMater = Air20;
  // fGapMater
  G4Element* H = new G4Element("Hydrogen", "H", 1., 1.00794*g/mole);
  G4double density;
  fGapMater = new G4Material("Water", density = 15*1.000*g/cm3, 2, kStateGas, 275.*kelvin, gapPressure=15*atmosphere);
  fGapMater -> AddElement(H, 2);
  fGapMater -> AddElement(O, 1);

  // fShieldMater Stainless steel
  G4Element* Cr = new G4Element("Chromium", "Cr", 24, 51.996*g/mole);
  G4Element* C = new G4Element("Carbon", "C", 6, 12.011*g/mole);
  G4Element* Mn = new G4Element("Manganese", "Mn", 25, 54.938*g/mole);
  G4Element* Si = new G4Element("Silicon", "Si", 14, 28.085*g/mole);
  G4Element* P = new G4Element("Phosphorous", "P", 15, 30.974*g/mole);
  G4Element* S = new G4Element("Sulphur", "S", 16, 32.06*g/mole);
  G4Element* Ni = new G4Element("Nickel", "Ni", 28, 58.693*g/mole);
  G4Element* Fe = new G4Element("Iron", "Fe", 26, 55.845*g/mole);
  //Stainless Steel 304. MIGHT BE ADJUSTED.
  fShieldMater = new G4Material("Stainless Steel", 8000*kg/m3, 9, kStateGas, 275.*kelvin, 55.584*atmosphere);
  fShieldMater -> AddElement(C, 0.0007);
  fShieldMater -> AddElement(Cr, 0.18);
  fShieldMater -> AddElement(Mn, 0.02);
  fShieldMater -> AddElement(Si, 0.01);
  fShieldMater -> AddElement(P, 0.00045);
  fShieldMater -> AddElement(S, 0.00015);
  fShieldMater -> AddElement(Ni, 0.09);
  fShieldMater -> AddElement(N, 0.001);
  fShieldMater -> AddElement(Fe, 0.6977);

  fTargetMater = mXenon136;

  //EndCap-North -------- Start -------- EndCap-North//

  fEndCapNorthDetectorMater = fDetectorMater;
  fEndCapNorthGapMater = fGapMater;
  fEndCapNorthShieldMater = fShieldMater;

  //EndCap-North -------- End -------- EndCap-North//

  //EndCap-South -------- Start -------- EndCap-South//

  fEndCapSouthDetectorMater = fDetectorMater;
  fEndCapSouthGapMater = fGapMater;
  fEndCapSouthShieldMater = fShieldMater;

  //EndCap-South -------- End -------- EndCap-South//

 ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // World
  //
  // (re) compute World dimensions if necessary
  fWorldLength = (std::max(std::max(fTargetLength,fDetectorLength), std::max(fGapLength,fShieldLength))) + (2*(fEndCapNorthDetectorLength+fEndCapNorthGapLength+fEndCapNorthShieldLength));
  fWorldRadius = fTargetRadius + fDetectorThickness + fGapThickness + fShieldThickness + (1.0*m);

  G4Tubs*
  sWorld = new G4Tubs("World",                                 //name
                 0.,fWorldRadius, 0.5*fWorldLength, 0.,twopi); //dimensions

  G4LogicalVolume*
  lWorld = new G4LogicalVolume(sWorld,                  //shape
                             fWorldMater,               //material
                             "World");                  //name

  fPhysiWorld = new G4PVPlacement(0,                    //no rotation
                            G4ThreeVector(),            //at (0,0,0)
                            lWorld,                     //logical volume
                            "World",                    //name
                            0,                          //mother volume
                            false,                      //no boolean operation
                            0);                         //copy number

  // Target
  //
  G4Tubs*
  sTarget = new G4Tubs("Target",                                   //name
                  0., fTargetRadius, 0.5*fTargetLength, 0.,twopi); //dimensions


  fLogicTarget = new G4LogicalVolume(sTarget,           //shape
                             fTargetMater,              //material
                             "Target");                 //name

           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(),             //at (0,0,0)
                           fLogicTarget,                //logical volume
                           "Target",                    //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number
  // Detector
  //
  G4Tubs* 
  sDetector = new G4Tubs("Detector",  
                fTargetRadius, (fTargetRadius+fDetectorThickness), (0.5*fDetectorLength) + fDetectorThickness, 0.,twopi);


  fLogicDetector = new G4LogicalVolume(sDetector,       //shape
                             fDetectorMater,            //material
                             "Detector");               //name
                               
           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(),             //at (0,0,0)
                           fLogicDetector,              //logical volume
                           "Detector",                  //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  // Gap
  //
  G4Tubs*
  sGap = new G4Tubs("Gap",
                (fDetectorThickness+fTargetRadius), (fDetectorThickness+fTargetRadius+fGapThickness), (0.5*fGapLength) + fDetectorThickness + fGapThickness, 0.,twopi);


  fLogicGap = new G4LogicalVolume(sGap,       //shape
                             fGapMater,            //material
                             "Gap");               //name

           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(),             //at (0,0,0)
                           fLogicGap,              //logical volume
                           "Gap",                  //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  // Shield
  //
  G4Tubs*
  sShield = new G4Tubs("Shield", (fDetectorThickness+fTargetRadius+fGapThickness), (fDetectorThickness+fTargetRadius+fGapThickness+fShieldThickness), (0.5*fShieldLength) + fDetectorThickness + fGapThickness + fShieldThickness, 0.,twopi);


  fLogicShield = new G4LogicalVolume(sShield,       //shape
                             fShieldMater,            //material
                             "Shield");               //name

           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(),             //at (0,0,0)
                           fLogicShield,              //logical volume
                           "Shield",                  //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  // EndCap-North-Detector ----------- Start ----------- EndCap-North-Detector //

  G4Tubs*
  sEndCapNorthDetector = new G4Tubs("EndCapNorthDetector", 0., fTargetRadius, 0.5*fEndCapNorthDetectorLength, 0.,twopi);


  fLogicEndCapNorthDetector = new G4LogicalVolume(sEndCapNorthDetector,       //shape
                             fEndCapNorthDetectorMater,            //material
                             "EndCapNorthDetector");               //name

           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(0, 0, (0.5*fEndCapNorthDetectorLength) + (0.5*(std::max(std::max(fTargetLength,fDetectorLength), std::max(fGapLength,fShieldLength))))),             //at (0,0,Z)
                           fLogicEndCapNorthDetector,              //logical volume
                           "EndCapNorthDetector",                  //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  // EndCap-North-Detector ----------- End ----------- EndCap-North-Detector //

  // EndCap-North-Gap ----------- Start ----------- EndCap-North-Gap //

  G4Tubs*
  sEndCapNorthGap = new G4Tubs("EndCapNorthGap", 0., fTargetRadius+fDetectorThickness, 0.5*fEndCapNorthGapLength, 0.,twopi);


  fLogicEndCapNorthGap = new G4LogicalVolume(sEndCapNorthGap,       //shape
                             fEndCapNorthGapMater,            //material
                             "EndCapNorthGap");               //name

           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(0, 0, ((0.5* fEndCapNorthGapLength) + fEndCapNorthDetectorLength) + (0.5*(std::max(std::max(fTargetLength,fDetectorLength), std::max(fGapLength,fShieldLength))))),             //at (0,0,Z)
                           fLogicEndCapNorthGap,              //logical volume
                           "EndCapNorthGap",                  //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  // EndCap-North-Gap ----------- End ----------- EndCap-North-Gap //

  // EndCap-North-Shield ----------- Start ----------- EndCap-North-Shield //

  G4Tubs*
  sEndCapNorthShield = new G4Tubs("EndCapNorthShield", 0., fTargetRadius+fDetectorThickness+fGapThickness, 0.5*fEndCapNorthShieldLength, 0.,twopi);


  fLogicEndCapNorthShield = new G4LogicalVolume(sEndCapNorthShield,       //shape
                             fEndCapNorthShieldMater,            //material
                             "EndCapNorthShield");               //name

           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(0, 0, ((0.5* fEndCapNorthShieldLength) + fEndCapNorthGapLength + fEndCapNorthDetectorLength) + (0.5*(std::max(std::max(fTargetLength,fDetectorLength), std::max(fGapLength,fShieldLength))))),             //at (0,0,Z)
                           fLogicEndCapNorthShield,              //logical volume
                           "EndCapNorthShield",                  //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  // EndCap-North-Shield ----------- End ----------- EndCap-North-Shield //

  // EndCap-South-Detector ----------- Start ----------- EndCap-South-Detector //

  G4Tubs*
  sEndCapSouthDetector = new G4Tubs("EndCapSouthDetector", 0., fTargetRadius, 0.5*fEndCapSouthDetectorLength, 0.,twopi);


  fLogicEndCapSouthDetector = new G4LogicalVolume(sEndCapSouthDetector,       //shape
                             fEndCapSouthDetectorMater,            //material
                             "EndCapSouthDetector");               //name

           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(0, 0, -0.5 * (std::max(std::max(fTargetLength,fDetectorLength), std::max(fGapLength,fShieldLength))) - (0.5*fEndCapSouthDetectorLength)),             //at (0,0,Z)
                           fLogicEndCapSouthDetector,              //logical volume
                           "EndCapSouthDetector",                  //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  // EndCap-South-Detector ----------- End ----------- EndCap-South-Detector //

  // EndCap-South-Gap ----------- Start ----------- EndCap-South-Gap //

  G4Tubs*
  sEndCapSouthGap = new G4Tubs("EndCapSouthGap", 0., fTargetRadius + fDetectorThickness, 0.5*fEndCapSouthGapLength, 0.,twopi);


  fLogicEndCapSouthGap = new G4LogicalVolume(sEndCapSouthGap,       //shape
                             fEndCapSouthGapMater,            //material
                             "EndCapSouthGap");               //name

           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(0, 0, -0.5 * (std::max(std::max(fTargetLength,fDetectorLength), std::max(fGapLength,fShieldLength))) - ((0.5*fEndCapSouthGapLength) + fEndCapSouthDetectorLength)),             //at (0,0,Z)
                           fLogicEndCapSouthGap,              //logical volume
                           "EndCapSouthGap",                  //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  // EndCap-South-Gap ----------- End ----------- EndCap-South-Gap //

  // EndCap-South-Shield ----------- Start ----------- EndCap-South-Shield //

  G4Tubs*
  sEndCapSouthShield = new G4Tubs("EndCapSouthShield", 0., fTargetRadius + fDetectorThickness + fGapThickness, 0.5*fEndCapSouthShieldLength, 0.,twopi);


  fLogicEndCapSouthShield = new G4LogicalVolume(sEndCapSouthShield,       //shape
                             fEndCapSouthShieldMater,            //material
                             "EndCapSouthShield");               //name

           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(0, 0, -0.5 * (std::max(std::max(fTargetLength,fDetectorLength), std::max(fGapLength,fShieldLength))) - ((0.5*fEndCapSouthShieldLength) + fEndCapSouthDetectorLength + fEndCapSouthGapLength)),             //at (0,0,Z)
                           fLogicEndCapSouthShield,              //logical volume
                           "EndCapSouthShield",                  //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  // EndCap-South-Shield ----------- End ----------- EndCap-South-Shield //


  PrintParameters();

  //always return the root volume
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4double xenonMass = (fTargetRadius/m) * (fTargetRadius/m) * pi * (fTargetLength/m) * densityXe136 * kg;

  G4double copperVolume = (((fTargetRadius/m) + (fDetectorThickness/m)) * ((fTargetRadius/m) + (fDetectorThickness/m)) * pi * ((fTargetLength/m) + (fDetectorThickness/m) * 2))-((fTargetRadius/m) * (fTargetRadius/m) * pi * (fTargetLength/m));
  G4double copperMass = copperVolume * 8960. * kg;

  G4double waterVolume = (((fTargetRadius/m) + (fDetectorThickness/m) + (fGapThickness/m)) * ((fTargetRadius/m) + (fDetectorThickness/m) + (fGapThickness/m)) * pi * ((fTargetLength/m) + 2*((fDetectorThickness/m) + (fGapThickness/m)))) - (((fTargetRadius/m) + (fDetectorThickness/m)) * ((fTargetRadius/m) + (fDetectorThickness/m)) * pi * ((fTargetLength/m) + (fDetectorThickness/m) * 2));
  G4double waterMass = waterVolume * 15000. * kg;

  G4double steelVolume = (((fTargetRadius/m) + (fDetectorThickness/m) + (fGapThickness/m) + (fShieldThickness/m)) * ((fTargetRadius/m) + (fDetectorThickness/m) + (fGapThickness/m) + (fShieldThickness/m)) * pi * ((fTargetLength/m) + 2*((fDetectorThickness/m) + (fGapThickness/m) + (fShieldThickness/m)))) - (((fTargetRadius/m) + (fDetectorThickness/m) + (fGapThickness/m)) * ((fTargetRadius/m) + (fDetectorThickness/m) + (fGapThickness/m)) * pi * ((fTargetLength/m) + 2*((fDetectorThickness/m) + (fGapThickness/m))));
  G4double steelMass = steelVolume * 8000. * kg;
  
  G4cout << "\n Target : Length = " << G4BestUnit(fTargetLength,"Length")
         << " Radius = " << G4BestUnit(fTargetRadius,"Length")
         << " Material = " << fTargetMater->GetName();
  G4cout << "\n Detector : Length = " << G4BestUnit(fDetectorLength,"Length")
         << " Thickness = " << G4BestUnit(fDetectorThickness,"Length")
         << " Material = " << fDetectorMater->GetName();
  G4cout << "\n Gap : Length = " << G4BestUnit(fGapLength,"Length")
         << " Thickness = " << G4BestUnit(fGapThickness,"Length")
         << " Material = " << fGapMater->GetName();
  G4cout << "\n Shield : Length = " << G4BestUnit(fShieldLength,"Length")
         << " Thickness = " << G4BestUnit(fShieldThickness,"Length")
         << " Material = " << fShieldMater->GetName();
  G4cout << "\n Total length = " << G4BestUnit((fTargetRadius + fDetectorThickness + fGapThickness + fShieldThickness),"Length") << G4endl;
  G4cout << "\n North End Cap: "
         << "\n Detector Thickness = " << G4BestUnit(fEndCapNorthDetectorLength,"Length")
         << "\n Gap Thickness = " << G4BestUnit(fEndCapNorthGapLength,"Length")
         << "\n Shield Thickness = " << G4BestUnit(fEndCapNorthShieldLength,"Length") << G4endl;
  G4cout << "\n South End Cap: "
         << "\n Detector Thickness = " << G4BestUnit(fEndCapSouthDetectorLength,"Length")
         << "\n Gap Thickness = " << G4BestUnit(fEndCapSouthGapLength,"Length")
         << "\n Shield Thickness = " << G4BestUnit(fEndCapSouthShieldLength,"Length")
         << G4endl;
  G4cout << "\n World : Length = " << G4BestUnit(fWorldLength,"Length")
         << " Radius = " << G4BestUnit(fWorldRadius,"Length") << G4endl;
  G4cout << "\n Target mass = " << G4BestUnit(xenonMass,"Mass")
         << "\n Detector mass = " << G4BestUnit(copperMass,"Mass")
         << "\n Gap mass = " << G4BestUnit(waterMass,"Mass")
         << "\n Shield mass = " << G4BestUnit(steelMass,"Mass") << G4endl;

  G4cout << "\n" << fTargetMater << "\n" << fDetectorMater << "\n" << fGapMater << "\n" << fShieldMater << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fTargetMater = pttoMaterial;
    if(fLogicTarget) { fLogicTarget->SetMaterial(fTargetMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetTargetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fDetectorMater = pttoMaterial;
    if(fLogicDetector) { fLogicDetector->SetMaterial(fDetectorMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
           << materialChoice << " not found" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGapMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fGapMater = pttoMaterial;
    if(fLogicGap) { fLogicGap->SetMaterial(fGapMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetGapMaterial : "
           << materialChoice << " not found" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetShieldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fShieldMater = pttoMaterial;
    if(fLogicGap) { fLogicShield->SetMaterial(fShieldMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetShieldMaterial : "
           << materialChoice << " not found" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetRadius(G4double value)
{
  fTargetRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetLength(G4double value)
{
  fTargetLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorThickness(G4double value)
{
  fDetectorThickness = value;
  fEndCapNorthDetectorLength = value;
  fEndCapSouthDetectorLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorLength(G4double value)
{
  fDetectorLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGapThickness(G4double value)
{
  fGapThickness = value;
  fEndCapNorthGapLength = value;
  fEndCapSouthGapLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGapLength(G4double value)
{
  fGapLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetShieldThickness(G4double value)
{
  fShieldThickness = value;
  fEndCapNorthShieldLength = value;
  fEndCapSouthShieldLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetShieldLength(G4double value)
{
  fShieldLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetPressure(G4double value)
{
  targetPressure = value;
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGapPressure(G4double value)
{
  gapPressure = value;
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetLength()
{
  return fTargetLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetRadius()
{
  return fTargetRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetTargetMaterial()
{
  return fTargetMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicTarget()
{
  return fLogicTarget;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorLength()
{
  return fDetectorLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorThickness()
{
  return fDetectorThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetDetectorMaterial()
{
  return fDetectorMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicDetector()
{
  return fLogicDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetGapLength()
{
  return fGapLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetGapThickness()
{
  return fGapThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetGapMaterial()
{
  return fGapMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicGap()
{
  return fLogicGap;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetShieldLength()
{
  return fShieldLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetShieldThickness()
{
  return fShieldThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetShieldMaterial()
{
  return fShieldMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicShield()
{
  return fLogicShield;
}
