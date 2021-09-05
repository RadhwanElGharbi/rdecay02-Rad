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
/// \file HitsAction.cc
/// \brief Implementation of the HitsAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HitsAction.hh"
#include "DetectorConstruction.hh"

SensitiveDetectorVol::SensitiveDetectorVol(G4String SDname) : G4VSensitiveDetector(SDname)
{
  G4cout<<"Creating SD with name: "<<SDname<<G4endl;
  collectionName.insert("SensitiveDetectorHitCollection");
}

SensitiveDetectorVol::~SensitiveDetectorVol()
{}

void SensitiveDetectorVol::Initialize(G4HCofThisEvent* HCE)
{
  hitCollection = new SensitiveDetectorHitCollection(GetName(), collectionName[0]);

  static G4int HCID = -1;
  if (HCID<0) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, hitCollection);
}

G4bool SensitiveDetectorVol::ProcessHits(G4Step *step, G4TouchableHistory *)
{
  G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();
  G4int copyNo = touchable->GetVolume(0)->GetCopyNo();
  G4double edep = step->GetTotalEnergyDeposit();

  SensitiveDetectorHit* aHit = new SensitiveDetectorHit(layerNumber);
  hitCollection->insert(aHit);
  aHit->AddEdep( edep );
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
