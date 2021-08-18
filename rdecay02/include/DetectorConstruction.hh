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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
    virtual G4VPhysicalVolume* Construct();
    
    void SetTargetLength (G4double value);
    void SetTargetRadius (G4double value);
    void SetTargetMaterial (G4String);
    
    void SetDetectorLength(G4double value);           
    void SetDetectorThickness(G4double value);  
    void SetDetectorMaterial(G4String);               
    
    void SetGapLength(G4double value);           
    void SetGapThickness(G4double value);  
    void SetGapMaterial(G4String);

    void SetShieldLength(G4double value);           
    void SetShieldThickness(G4double value);  
    void SetShieldMaterial(G4String);

    void SetEndCapNorthDetectorLength(G4double value);
    void SetEndCapNorthGapLength(G4double value);
    void SetEndCapNorthShieldLength(G4double value);

    void SetEndCapSouthDetectorLength(G4double value);
    void SetEndCapSouthGapLength(G4double value);
    void SetEndCapSouthShieldLength(G4double value);


    void SetTargetPressure(G4double value);
    void SetGapPressure(G4double value);
                   
    void PrintParameters();
    
  public:
      
    G4double GetTargetLength();
    G4double GetTargetRadius();
    G4Material* GetTargetMaterial();       
    G4LogicalVolume* GetLogicTarget();
    
    G4double GetDetectorLength();
    G4double GetDetectorThickness();
    G4Material* GetDetectorMaterial();                 
    G4LogicalVolume* GetLogicDetector();

    G4double GetGapLength();
    G4double GetGapThickness();
    G4Material* GetGapMaterial();                 
    G4LogicalVolume* GetLogicGap();

    G4double GetShieldLength();
    G4double GetShieldThickness();
    G4Material* GetShieldMaterial();                 
    G4LogicalVolume* GetLogicShield(); 

    //EndCap-North -------- Start -------- EndCap-North//

    G4double GetEndCapNorthDetectorLength();
    G4double GetEndCapNorthDetectorRadius();
    G4Material* GetEndCapNorthDetectorMaterial();                 
    G4LogicalVolume* GetLogicEndCapNorthDetector();   

    G4double GetEndCapNorthGapLength();
    G4double GetEndCapNorthGapRadius();
    G4Material* GetEndCapNorthGapMaterial();                 
    G4LogicalVolume* GetLogicEndCapNorthGap();   

    G4double GetEndCapNorthShieldLength();
    G4double GetEndCapNorthShieldRadius();
    G4Material* GetEndCapNorthShieldMaterial();                 
    G4LogicalVolume* GetLogicEndCapNorthShield();  

    //EndCap-North -------- End -------- EndCap-North//

    //EndCap-South -------- Start -------- EndCap-South//

    G4double GetEndCapSouthDetectorLength();
    G4double GetEndCapSouthDetectorRadius();
    G4Material* GetEndCapSouthDetectorMaterial();                 
    G4LogicalVolume* GetLogicEndCapSouthDetector();   

    G4double GetEndCapSouthGapLength();
    G4double GetEndCapSouthGapRadius();
    G4Material* GetEndCapSouthGapMaterial();                 
    G4LogicalVolume* GetLogicEndCapSouthGap();   

    G4double GetEndCapSouthShieldLength();
    G4double GetEndCapSouthShieldRadius();
    G4Material* GetEndCapSouthShieldMaterial();                 
    G4LogicalVolume* GetLogicEndCapSouthShield();  

    //EndCap-South -------- End -------- EndCap-South//     


                       
  private:
  
    G4double           fTargetLength;
    G4double           fTargetRadius;
    G4Material*        fTargetMater;
    G4LogicalVolume*   fLogicTarget;
                 
    G4double           fDetectorLength;
    G4double           fDetectorThickness;
    G4Material*        fDetectorMater;
    G4LogicalVolume*   fLogicDetector;

    G4double           fGapLength;
    G4double           fGapThickness;
    G4Material*        fGapMater;
    G4LogicalVolume*   fLogicGap;

    G4double           fShieldLength;
    G4double           fShieldThickness;
    G4Material*        fShieldMater;
    G4LogicalVolume*   fLogicShield;
    
    //EndCap-North -------- Start -------- EndCap-North//
    
    G4double           fEndCapNorthDetectorLength;
    G4double           fEndCapNorthDetectorRadius;
    G4Material*        fEndCapNorthDetectorMater;
    G4LogicalVolume*   fLogicEndCapNorthDetector;

    G4double           fEndCapNorthGapLength;
    G4double           fEndCapNorthGapRadius;
    G4Material*        fEndCapNorthGapMater;
    G4LogicalVolume*   fLogicEndCapNorthGap;

    G4double           fEndCapNorthShieldLength;
    G4double           fEndCapNorthShieldRadius;
    G4Material*        fEndCapNorthShieldMater;
    G4LogicalVolume*   fLogicEndCapNorthShield;

    //EndCap-North -------- End -------- EndCap-North//

    //EndCap-South -------- Start -------- EndCap-South//
    
    G4double           fEndCapSouthDetectorLength;
    G4double           fEndCapSouthDetectorRadius;
    G4Material*        fEndCapSouthDetectorMater;
    G4LogicalVolume*   fLogicEndCapSouthDetector;

    G4double           fEndCapSouthGapLength;
    G4double           fEndCapSouthGapRadius;
    G4Material*        fEndCapSouthGapMater;
    G4LogicalVolume*   fLogicEndCapSouthGap;

    G4double           fEndCapSouthShieldLength;
    G4double           fEndCapSouthShieldRadius;
    G4Material*        fEndCapSouthShieldMater;
    G4LogicalVolume*   fLogicEndCapSouthShield;

    //EndCap-South -------- End -------- EndCap-South//
               
    G4double           fWorldLength;
    G4double           fWorldRadius;
    G4Material*        fWorldMater;
    G4VPhysicalVolume* fPhysiWorld;

    G4Material* mXenon;
    G4Material* mXenon136;
    
    G4double targetPressure;
    G4double gapPressure;
    const G4double molarMassXe = 131.292;
    const G4double molarMassXe136 = 135.91;
    G4double densityXe;
    G4double densityXe136;
                
    DetectorMessenger* fDetectorMessenger;

  private:
    
    void               DefineMaterials();
    G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

