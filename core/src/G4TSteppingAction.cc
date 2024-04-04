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
// 
/// \file G4TSteppingAction.cc
/// \brief Implementation of the G4TSteppingAction class

#include "G4TSteppingAction.hh"
#include "G4TRunAction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"

//G4ThreadLocal G4double G4TSteppingAction::edep;
//G4ThreadLocal G4String G4TSteppingAction::reg;
//G4ThreadLocal unsigned int  G4TSteppingAction::CN;

G4TSteppingAction::G4TSteppingAction(G4TRunAction* ra):G4UserSteppingAction(),RunAction(ra){}

G4TSteppingAction::~G4TSteppingAction(){}

void G4TSteppingAction::UserSteppingAction(const G4Step* step)
{

    //G4String Particlename = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    //if(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "nFission"){
    //    std::cout << " EventID " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
    //          << " TrackID " << step->GetTrack()->GetTrackID()
    //          << " Particlename " << Particlename
    //          << " ProcessName " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
    //          << std::endl;
    //}

    auto edep = step->GetTotalEnergyDeposit();

    //std::cout <<" edep:" << edep << " @@@@@@@@@@@@@@@@@@@@@@@" << std::endl;

    if (edep == 0.) return;

    auto reg = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName();
    RunAction->FillRegionStepHits(reg, edep);

    //std::cout << "---- reg:" << reg << " edep:" << edep << " @@@@@@@@@@@@@@@@@@@@@@@" << std::endl;

}
