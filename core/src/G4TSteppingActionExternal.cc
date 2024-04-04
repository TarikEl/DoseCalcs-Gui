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
/// \file G4TSteppingActionExternal.cc
/// \brief Implementation of the G4TSteppingActionExternal class

#include "G4TSteppingActionExternal.hh"
#include "G4TRunAction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

extern G4String GenerateVoxelsResuls;
extern G4String* CopyNumberRegionNameMap;

extern G4int VoxXNumber;
extern G4int VoxYNumber;
extern G4int VoxZNumber;
extern G4int ParamType;

//G4ThreadLocal G4double G4TSteppingActionExternal::edep;
//G4ThreadLocal G4String G4TSteppingActionExternal::reg;
//G4ThreadLocal unsigned int  G4TSteppingActionExternal::CN;

G4TSteppingActionExternal::G4TSteppingActionExternal(G4TRunAction* ra):G4UserSteppingAction(),RunAction(ra){}

G4TSteppingActionExternal::~G4TSteppingActionExternal(){}

void G4TSteppingActionExternal::UserSteppingAction(const G4Step* step)
{

    //G4String Particlename = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    //if(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "nFission"){
    //    std::cout << " EventID " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
    //          << " TrackID " << step->GetTrack()->GetTrackID()
    //          << " Particlename " << Particlename
    //          << " ProcessName " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
    //          << std::endl;
    //}

    const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
    auto edep = step->GetTotalEnergyDeposit();


    if(GenerateVoxelsResuls == "yes"){
        G4int CN ;

        if(ParamType == 0){
            CN = touchable->GetCopyNumber();
            if(step->GetTrack()->GetTrackID() == 1){
                //std::cout <<" edep:" << edep << " SL " << step->GetStepLength() << std::endl;
                RunAction->FillVoxelLenghts(CN, 0.1*step->GetStepLength()); // in cm
                //std::cout << " Particle:" << step->GetTrack()->GetParticleDefinition()->GetParticleName() << " ID:" << step->GetTrack()->GetTrackID() << " reg:" << reg << " Edep:" << step->GetTotalEnergyDeposit() << " Lenght:" << step->GetStepLength() << " @@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
            }
            RunAction->FillVoxelStepHits(touchable->GetCopyNumber(), edep);
        }else{
            if(touchable->GetHistoryDepth()!=0){
                CN = touchable->GetReplicaNumber(1) + VoxXNumber*touchable->GetReplicaNumber(2) + VoxXNumber*VoxYNumber*touchable->GetReplicaNumber(0);
                if(step->GetTrack()->GetTrackID() == 1){
                    //std::cout <<" edep:" << edep << " SL " << step->GetStepLength() << std::endl;
                    RunAction->FillVoxelLenghts(CN, 0.1*step->GetStepLength()); // in cm
                    //std::cout << " Particle:" << step->GetTrack()->GetParticleDefinition()->GetParticleName() << " ID:" << step->GetTrack()->GetTrackID() << " reg:" << reg << " Edep:" << step->GetTotalEnergyDeposit() << " Lenght:" << step->GetStepLength() << " @@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
                }
                RunAction->FillVoxelStepHits(CN, edep);
            }
        }
    }else{
        auto reg = touchable->GetVolume()->GetLogicalVolume()->GetName();
        if(reg != "World"){
            if(step->GetTrack()->GetTrackID() == 1){
                //std::cout <<" edep:" << edep << " SL " << step->GetStepLength() << std::endl;
                RunAction->FillRegionLenghts(reg, 0.1*step->GetStepLength()); // in cm
                //std::cout << " Particle:" << step->GetTrack()->GetParticleDefinition()->GetParticleName() << " ID:" << step->GetTrack()->GetTrackID() << " reg:" << reg << " Edep:" << step->GetTotalEnergyDeposit() << " Lenght:" << step->GetStepLength() << " @@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
            }
            RunAction->FillRegionStepHits(reg, edep);
        }
    }

    //if (edep == 0.) return;

}
