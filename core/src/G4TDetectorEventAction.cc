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
/// \file optical/G4TDetector/src/G4TDetectorEventAction.cc
/// \brief Implementation of the G4TDetectorEventAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4TDetectorEventAction.hh"
//#include "G4TDetectorRun.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"

extern G4String DetectFor;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4TDetectorEventAction::G4TDetectorEventAction(G4TRunAction* runaction)
  : G4UserEventAction()
{
  fRayleigh   = 0;
  fAbsorption = 0;
  fMie        = 0;
  fBoundary   = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4TDetectorEventAction::~G4TDetectorEventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4TDetectorEventAction::BeginOfEventAction(const G4Event*)
{
  fRayleigh   = 0;
  fAbsorption = 0;
  fMie        = 0;
  fBoundary   = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4TDetectorEventAction::EndOfEventAction(const G4Event* )
{

    G4cout << " EndOfEventAction 1" << G4endl;

    //runAction->AddEventsDataInDetectors(G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(),decetorEventEnerDepoMap,decetorEventEnerDepoMap2,decetorEventHitMap);
    G4cout << " EndOfEventAction 2" << G4endl;

    //if(DetectFor == "ForCNN"){
    //}else{
    //}

  //G4TDetectorRun* run = static_cast<G4TDetectorRun*>( G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  //run->AddRayleigh(fRayleigh);
  //run->AddAbsorption(fAbsorption);
  //run->AddMie(fMie);
  //run->AddBoundary(fBoundary);
}

//void setEnerDepo(G4double EE) { EnerDepo = EE; }
void G4TDetectorEventAction::setDetectorNameMap(G4String SS, double EnerDepo) {

    G4cout << " SS " << SS  << " EnerDepo " << EnerDepo << G4endl;


    decetorEventEnerDepoMap["World"] = 0.1;

    //if(decetorEventEnerDepoMap[SS] == NULL){
    //    decetorEventEnerDepoMap[SS] = 0;
    //}

    G4cout << decetorEventEnerDepoMap[SS]  << G4endl;

    G4cout << " \n\n\n\n\n\n\n\n\n\nSS " << SS  << " EnerDepo " << EnerDepo << G4endl;

    G4cout << decetorEventEnerDepoMap2[SS]  << G4endl;
    G4cout << decetorEventHitMap[SS]  << G4endl;


    decetorEventEnerDepoMap[SS] = decetorEventEnerDepoMap[SS] + EnerDepo;
    decetorEventEnerDepoMap2[SS] = decetorEventEnerDepoMap2[SS] + EnerDepo*EnerDepo;
    decetorEventHitMap[SS] = decetorEventHitMap[SS] + 1;

    G4cout << " End of setDetectorNameMap " << G4endl;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
