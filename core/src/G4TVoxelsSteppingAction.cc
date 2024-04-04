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
/// \file G4TVoxelsSteppingAction.cc
/// \brief Implementation of the G4TVoxelsSteppingAction class

#include "G4TVoxelsSteppingAction.hh"
#include "G4TRunAction.hh"
#include "G4Step.hh"

extern G4String GenerateVoxelsResuls;
extern G4String* CopyNumberRegionNameMap;

extern G4int VoxXNumber;
extern G4int VoxYNumber;
extern G4int VoxZNumber;
extern G4int ParamType;

//G4ThreadLocal G4double G4TVoxelsSteppingAction::edep;
//G4ThreadLocal G4String G4TVoxelsSteppingAction::reg;
//G4ThreadLocal unsigned int  G4TVoxelsSteppingAction::CN;

G4TVoxelsSteppingAction::G4TVoxelsSteppingAction(G4TRunAction* ra):G4UserSteppingAction(),RunAction(ra){}


G4TVoxelsSteppingAction::~G4TVoxelsSteppingAction(){}


void G4TVoxelsSteppingAction::UserSteppingAction(const G4Step* step)
{

    auto edep = step->GetTotalEnergyDeposit();

    if (edep == 0.) return;

    const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();

    if(ParamType == 0){
        //G4int CN = step->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
        RunAction->FillVoxelStepHits(touchable->GetCopyNumber(), edep);
    }else{
        if(step->GetPreStepPoint()->GetTouchable()->GetHistoryDepth()!=0){
            RunAction->FillVoxelStepHits(touchable->GetReplicaNumber(1) + VoxXNumber*touchable->GetReplicaNumber(2) + VoxXNumber*VoxYNumber*touchable->GetReplicaNumber(0), edep);
        }
    }

    //std::cout << "TrackID: "<< step->GetTrack()->GetTrackID() << " ParticleName: " << step->GetTrack()->GetParticleDefinition()->GetParticleName() << " Length "<< step->GetTrack()->GetTrackLength() << " Energy(MeV): "<< edep << std::endl;

    /*
                //std::cout << " Touchable " << std::endl;
                //std::cout << " Begin " << touchable->GetHistoryDepth() << " " << touchable->GetVolume()->GetName() << std::endl;
                G4int ix = touchable->GetReplicaNumber(1);
                //std::cout << " ix " << touchable->GetVolume(1)->GetName() << std::endl;
                G4int iy = touchable->GetReplicaNumber(2);
                //std::cout << " iy " << touchable->GetVolume(2)->GetName() << std::endl;
                G4int iz = touchable->GetReplicaNumber(0);
                //std::cout << " iz " << touchable->GetVolume(0)->GetName() << std::endl;
                G4int cp = ix + VoxXNumber*iy + VoxXNumber*VoxYNumber*iz;
                auto reg = CopyNumberRegionNameMap[cp];
                //std::cout << " ix " << ix << " iy " << iy << " iz " << iz << " cp " << cp <<  " reg " << reg << std::endl;
                //RunAction->FillRegionStepHits(reg, edep);
                */
    /*
        if(GenerateVoxelsResuls == "yes"){
            //auto CN = step->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
            RunAction->FillVoxelStepHits(touchable->GetReplicaNumber(1) + VoxXNumber*touchable->GetReplicaNumber(2) + VoxXNumber*VoxYNumber*touchable->GetReplicaNumber(0), edep);
        }else{

        }
        */

}
