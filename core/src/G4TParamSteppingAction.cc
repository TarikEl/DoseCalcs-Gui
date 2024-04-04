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
/// \file G4TParamSteppingAction.cc
/// \brief Implementation of the G4TParamSteppingAction class

#include "G4TParamSteppingAction.hh"
#include "G4TRunAction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4RunManager.hh"

//extern G4String GenerateVoxelsResuls;
extern G4String* CopyNumberRegionNameMap;

//G4ThreadLocal G4double G4TParamSteppingAction::edep;
//G4ThreadLocal G4String G4TParamSteppingAction::reg;
//G4ThreadLocal unsigned int  G4TParamSteppingAction::CN;

G4TParamSteppingAction::G4TParamSteppingAction(G4TRunAction* ra):G4UserSteppingAction(),RunAction(ra){}

G4TParamSteppingAction::~G4TParamSteppingAction(){}


void G4TParamSteppingAction::UserSteppingAction(const G4Step* step)
{

    auto edep = step->GetTotalEnergyDeposit();
    if (edep == 0.) return;
    G4int CN = step->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
    RunAction->FillRegionStepHits(CopyNumberRegionNameMap[CN], edep);

    //RunAction->FillVoxelStepHits(CN, edep);






    //else stepCounter=0;

    //std::cout << " Energy(MeV): "<< edep << " CN: "<< CN << " CopyNumberRegionNameMap[CN]: "<< CopyNumberRegionNameMap[CN] << "  << std::endl;

    //std::cout << "---- ParticleName: " << step->GetTrack()->GetParticleDefinition()->GetParticleName() << " Energy(MeV): "<< edep << " ProcessName: "<< step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;

    //if(step->GetPreStepPoint()->GetTouchable()!= nullptr){

    //const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();

    //std::cout << "---- CN:" << touchable->GetCopyNumber() << " reg:"<< CopyNumberRegionNameMap[touchable->GetCopyNumber()] << " edep:" << edep << " @@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
    //}

    //auto reg = CopyNumberRegionNameMap[step->GetPreStepPoint()->GetTouchable()->GetCopyNumber()];
    //auto reg = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    //G4cout << "---- CN:" << step->GetPreStepPoint()->GetTouchable()->GetCopyNumber() << " reg:"<< reg << " edep:" << edep << " @@@@@@@@@@@@@@@@@@@@@@@" << G4endl;

    /*
    if(GenerateVoxelsResuls == "yes"){
        //auto CN = step->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
        RunAction->FillVoxelStepHits(step->GetPreStepPoint()->GetTouchable()->GetCopyNumber(), edep);
    }else{

    }
    */

}
