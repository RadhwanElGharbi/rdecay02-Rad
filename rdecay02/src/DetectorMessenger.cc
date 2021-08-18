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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
:G4UImessenger(), 
 fDetector(Det), fRdecayDir(0), fDetDir(0),
 fTargMatCmd(0), fDetectMatCmd(0), fGapMatCmd(0), fShieldMatCmd(0),
 fTargRadiusCmd(0), fDetectThicknessCmd(0), fGapThicknessCmd(0), fShieldThicknessCmd(0),
 fTargLengthCmd(0), fDetectLengthCmd(0), fGapLengthCmd(0), fShieldLengthCmd(0),
 fTargPressureCmd(0), fGapPressureCmd(0)
{ 
  fRdecayDir = new G4UIdirectory("/rdecay02/");
  fRdecayDir->SetGuidance("commands specific to this example");
  
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/rdecay02/det/",broadcast);
  fDetDir->SetGuidance("detector construction commands");
        
  fTargMatCmd = new G4UIcmdWithAString("/rdecay02/det/setTargetMate",this);
  fTargMatCmd->SetGuidance("Select material of the target");
  fTargMatCmd->SetParameterName("choice",false);
  fTargMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fTargRadiusCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setTargetRadius", this);
  fTargRadiusCmd->SetGuidance("Set the Target Radius.");
  fTargRadiusCmd->SetUnitCategory("Length");
  fTargRadiusCmd->SetParameterName("choice",false);
  fTargRadiusCmd->AvailableForStates(G4State_PreInit);  

  
  fTargLengthCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setTargetLength", this);
  fTargLengthCmd->SetGuidance("Set the Target Length.");
  fTargLengthCmd->SetUnitCategory("Length");
  fTargLengthCmd->SetParameterName("choice",false);
  fTargLengthCmd->AvailableForStates(G4State_PreInit);
  

  fDetectMatCmd = new G4UIcmdWithAString("/rdecay02/det/setDetectorMate",this);
  fDetectMatCmd->SetGuidance("Select Material of the Detector.");
  fDetectMatCmd->SetParameterName("choice",false);
  fDetectMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fDetectThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setDetectorThickness",this);
  fDetectThicknessCmd->SetGuidance("Set the Detector Thickness.");
  fDetectThicknessCmd->SetUnitCategory("Length");
  fDetectThicknessCmd->SetParameterName("choice",false);
  fDetectThicknessCmd->AvailableForStates(G4State_PreInit);

  fDetectLengthCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setDetectorLength",this);
  fDetectLengthCmd->SetGuidance("Set the Detector Length.");
  fDetectLengthCmd->SetUnitCategory("Length");
  fDetectLengthCmd->SetParameterName("choice",false);
  fDetectLengthCmd->AvailableForStates(G4State_PreInit);

  //=========================================================================//

  fGapMatCmd = new G4UIcmdWithAString("/rdecay02/det/setGapMate",this);
  fGapMatCmd->SetGuidance("Select Material of the Gap.");
  fGapMatCmd->SetParameterName("choice",false);
  fGapMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fGapThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setGapThickness",this);
  fGapThicknessCmd->SetGuidance("Set the Gap Thickness.");
  fGapThicknessCmd->SetUnitCategory("Length");
  fGapThicknessCmd->SetParameterName("choice",false);
  fGapThicknessCmd->AvailableForStates(G4State_PreInit);

  fGapLengthCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setGapLength",this);
  fGapLengthCmd->SetGuidance("Set the Gap Length.");
  fGapLengthCmd->SetUnitCategory("Length");
  fGapLengthCmd->SetParameterName("choice",false);
  fGapLengthCmd->AvailableForStates(G4State_PreInit);

  //=========================================================================//

  fShieldMatCmd = new G4UIcmdWithAString("/rdecay02/det/setShieldMate",this);
  fShieldMatCmd->SetGuidance("Select Material of the Shield.");
  fShieldMatCmd->SetParameterName("choice",false);
  fShieldMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fShieldThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setShieldThickness",this);
  fShieldThicknessCmd->SetGuidance("Set the Shield Thickness.");
  fShieldThicknessCmd->SetUnitCategory("Length");
  fShieldThicknessCmd->SetParameterName("choice",false);
  fShieldThicknessCmd->AvailableForStates(G4State_PreInit);

  fShieldLengthCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setShieldLength",this);
  fShieldLengthCmd->SetGuidance("Set the Shield Length.");
  fShieldLengthCmd->SetUnitCategory("Length");
  fShieldLengthCmd->SetParameterName("choice",false);
  fShieldLengthCmd->AvailableForStates(G4State_PreInit);

  //=========================================================================//

  fTargPressureCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setTargetPressure",this);
  fTargPressureCmd->SetGuidance("Set the Target Pressure.");
  fTargPressureCmd->SetUnitCategory("Pressure");
  fTargPressureCmd->SetParameterName("choice",false);
  fTargPressureCmd->AvailableForStates(G4State_PreInit);

  fGapPressureCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setGapPressure",this);
  fGapPressureCmd->SetGuidance("Set the Gap Pressure.");
  fGapPressureCmd->SetUnitCategory("Pressure");
  fGapPressureCmd->SetParameterName("choice",false);
  fGapPressureCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fTargMatCmd;
  delete fDetectMatCmd;
  delete fGapMatCmd;
  delete fShieldMatCmd;
  delete fTargRadiusCmd;
  delete fDetectThicknessCmd;
  delete fGapThicknessCmd;
  delete fShieldThicknessCmd;
  delete fTargLengthCmd;
  delete fDetectLengthCmd;
  delete fGapLengthCmd;
  delete fShieldLengthCmd;
  delete fTargPressureCmd;
  delete fGapPressureCmd;
  delete fDetDir;
  delete fRdecayDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == fTargMatCmd )
   { fDetector->SetTargetMaterial(newValue);}
   
  if (command == fTargLengthCmd ) 
    { fDetector->SetTargetLength(fTargLengthCmd->GetNewDoubleValue(newValue));}
    
  if (command == fTargRadiusCmd ) 
    {fDetector->SetTargetRadius(fTargLengthCmd->GetNewDoubleValue(newValue));}
    
  if (command == fDetectMatCmd )
    { fDetector->SetDetectorMaterial(newValue);}
    
  if (command == fDetectLengthCmd ) 
    {fDetector->SetDetectorLength(
                     fDetectLengthCmd->GetNewDoubleValue(newValue));}

  if (command == fDetectThicknessCmd ) 
    {fDetector->SetDetectorThickness(
                     fDetectThicknessCmd->GetNewDoubleValue(newValue));}
     
  //==========================================//
    
  if (command == fGapMatCmd )
    { fDetector->SetGapMaterial(newValue);}
    
  if (command == fGapLengthCmd ) 
    {fDetector->SetGapLength(
                     fGapLengthCmd->GetNewDoubleValue(newValue));}

  if (command == fGapThicknessCmd ) 
    {fDetector->SetGapThickness(
                     fGapThicknessCmd->GetNewDoubleValue(newValue));}  

  //==========================================//
    
  if (command == fShieldMatCmd )
    { fDetector->SetShieldMaterial(newValue);}
    
  if (command == fShieldLengthCmd ) 
    {fDetector->SetShieldLength(
                     fShieldLengthCmd->GetNewDoubleValue(newValue));}

  if (command == fShieldThicknessCmd ) 
    {fDetector->SetShieldThickness(
                     fShieldThicknessCmd->GetNewDoubleValue(newValue));}  

  //==========================================//

  if (command == fTargPressureCmd ) 
    {fDetector->SetTargetPressure(
                     fTargPressureCmd->GetNewDoubleValue(newValue));}

  if (command == fGapPressureCmd ) 
    {fDetector->SetGapPressure(
                     fGapPressureCmd->GetNewDoubleValue(newValue));} 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
