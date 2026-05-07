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
/// \file optical/G4TDetector/include/G4TDetectorEventAction.hh
/// \brief Definition of the G4TDetectorEventAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4TDetectorEventAction_h
#define G4TDetectorEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"
#include "G4TRunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4TDetectorEventAction : public G4UserEventAction
{
public:
    G4TDetectorEventAction(G4TRunAction*);
    virtual ~G4TDetectorEventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    std::map<G4String,unsigned long long int> decetorEventHitMap;
    std::map<G4String,G4double> decetorEventEnerDepoMap;
    std::map<G4String,G4double> decetorEventEnerDepoMap2;
    void AddRayleigh() { ++fRayleigh; }
    void AddAbsorption() { ++fAbsorption; }
    void AddMie() { ++fMie; }
    void AddBoundary() { ++fBoundary; }

    //void setEnerDepo(G4double EE) { EnerDepo = EE; }
    void setDetectorNameMap(G4String SS, double EnerDepo);

    //  G4double getEnetDepo() { ++fBoundary; }

private:
    G4double EnerDepo;
    G4int fRayleigh;
    G4int fAbsorption;
    G4int fMie;
    G4int fBoundary;
    G4TRunAction* runAction;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
