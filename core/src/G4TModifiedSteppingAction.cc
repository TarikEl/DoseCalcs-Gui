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
/// \file G4TModifiedSteppingAction.cc
/// \brief Implementation of the G4TModifiedSteppingAction class

#include "G4TModifiedSteppingAction.hh"
#include "G4TRunAction.hh"
#include "G4Step.hh"



#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"

G4TModifiedSteppingAction::G4TModifiedSteppingAction(G4TRunAction* ra):G4UserSteppingAction(),RunAction(ra){}

G4TModifiedSteppingAction::~G4TModifiedSteppingAction(){}


void G4TModifiedSteppingAction::UserSteppingAction(const G4Step* step)
{
    // For Protontherapy voxelized phantom irradiation

    auto edep = step->GetTotalEnergyDeposit();
    if (edep == 0.) return;
    const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
    if(touchable->GetVolume()->GetName() == "Voxel"){
        RunAction->FillVoxelStepHits(touchable->GetCopyNumber(), edep);
        //G4cout << "---- CN:" << touchable->GetCopyNumber() << " VolumeName" << touchable->GetVolume()->GetName() <<  " edep:" << edep << " @@@@@@@@@@@@@@@@@@@@@@@" << G4endl;
    }
}

/////////////////////////////////////  For Protontherapy voxelized phantom irradiation ///////////////////////////////////////////////////////////////////////////
/*
auto edep = step->GetTotalEnergyDeposit();
if (edep == 0.) return;
const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
if(touchable->GetVolume()->GetName() == "Voxel"){
    RunAction->FillVoxelStepHits(touchable->GetCopyNumber(), edep);
    //G4cout << "---- CN:" << touchable->GetCopyNumber() << " VolumeName" << touchable->GetVolume()->GetName() <<  " edep:" << edep << " @@@@@@@@@@@@@@@@@@@@@@@" << G4endl;
}
*/
///////////////////////////////////// For Detector, Scintillation and ////////////////////////////////////////////////////////////////////////////
/*
    #include "G4OpBoundaryProcess.hh"
    #include "G4OpticalPhoton.hh"

    static G4ParticleDefinition* opticalphoton = G4OpticalPhoton::OpticalPhotonDefinition();
    const G4ParticleDefinition* particleDef = step->GetTrack()->GetDynamicParticle()->GetParticleDefinition();
    RunAction->CountOpticalInteractions(step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName()+"_"+particleDef->GetParticleName(), step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
*/
/*
    extern G4String* CopyNumberRegionNameMap;
    extern G4int VoxXNumber;
    extern G4int VoxYNumber;
    extern G4int VoxZNumber;

    auto edep = step->GetTotalEnergyDeposit();
    if (edep == 0.) return;
    auto reg = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName();
    if(reg != "RO"){return;}

    const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
    RunAction->FillVoxelStepHits(touchable->GetReplicaNumber(2) + VoxXNumber*touchable->GetReplicaNumber(1) + VoxXNumber*VoxYNumber*touchable->GetReplicaNumber(0), edep);
*/

///////////////////////////////////// For Stylized Phantoms, Internal Dosimetry  ////////////////////////////////////////////////////////////////////////////

/*
    auto edep = step->GetTotalEnergyDeposit();
    if (edep == 0.) return;
    auto reg = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName();
    RunAction->FillRegionStepHits(reg, edep);
*/

///////////////////////////////////// For Voxelized Phantoms, Internal Dosimetry, Region Level  ////////////////////////////////////////////////////////////////////////////

/*   For Parammetrization 0
extern G4String* CopyNumberRegionNameMap;

auto edep = step->GetTotalEnergyDeposit();
if (edep == 0.) return;
G4int CN = step->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
RunAction->FillRegionStepHits(CopyNumberRegionNameMap[CN], edep);
*/


/*   For Nested Parammetrization 1
extern G4String* CopyNumberRegionNameMap;
extern G4int VoxXNumber;
extern G4int VoxYNumber;
extern G4int VoxZNumber;

auto edep = step->GetTotalEnergyDeposit();
if (edep == 0.) return;
if(step->GetPreStepPoint()->GetTouchable()->GetHistoryDepth()!=0){
    const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
    RunAction->FillRegionStepHits(CopyNumberRegionNameMap[touchable->GetReplicaNumber(1) + VoxXNumber*touchable->GetReplicaNumber(2) + VoxXNumber*VoxYNumber*touchable->GetReplicaNumber(0)], edep);
}
*/

///////////////////////////////////// For Voxelized Phantoms, Internal Dosimetry, Voxel Level  ////////////////////////////////////////////////////////////////////////////

/*   For Parammetrization 0
extern G4String* CopyNumberRegionNameMap;

auto edep = step->GetTotalEnergyDeposit();
if (edep == 0.) return;
const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
RunAction->FillVoxelStepHits(touchable->GetCopyNumber(), edep);
*/


/*   For Nested Parammetrization 1
extern G4String* CopyNumberRegionNameMap;
extern G4int VoxXNumber;
extern G4int VoxYNumber;
extern G4int VoxZNumber;

auto edep = step->GetTotalEnergyDeposit();
if (edep == 0.) return;
const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
RunAction->FillVoxelStepHits(touchable->GetCopyNumber(), edep);
if(step->GetPreStepPoint()->GetTouchable()->GetHistoryDepth()!=0){
        RunAction->FillVoxelStepHits(touchable->GetReplicaNumber(1) + VoxXNumber*touchable->GetReplicaNumber(2) + VoxXNumber*VoxYNumber*touchable->GetReplicaNumber(0), edep);
}
*/


///////////////////////////////////// For Voxelized Phantoms, External Dosimetry, Region Level  ////////////////////////////////////////////////////////////////////////////


/*   For Parammetrization 0
extern G4String* CopyNumberRegionNameMap;

const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
auto edep = step->GetTotalEnergyDeposit();
G4int CN = touchable->GetCopyNumber();
if(step->GetTrack()->GetTrackID() == 1){
    RunAction->FillRegionLenghts(CopyNumberRegionNameMap[CN], 0.1*step->GetStepLength()); // in cm
}
RunAction->FillRegionStepHits(CopyNumberRegionNameMap[CN], edep);
*/


/*   For Nested Parammetrization 1
extern G4String* CopyNumberRegionNameMap;
extern G4int VoxXNumber;
extern G4int VoxYNumber;
extern G4int VoxZNumber;

const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
auto edep = step->GetTotalEnergyDeposit();
if(touchable->GetHistoryDepth()!=0){
    G4int CN = touchable->GetReplicaNumber(1) + VoxXNumber*touchable->GetReplicaNumber(2) + VoxXNumber*VoxYNumber*touchable->GetReplicaNumber(0);
    if(step->GetTrack()->GetTrackID() == 1){
        RunAction->FillRegionLenghts(CopyNumberRegionNameMap[CN], 0.1*step->GetStepLength()); // in cm
    }
    RunAction->FillRegionStepHits(CopyNumberRegionNameMap[CN], edep);
}
*/

///////////////////////////////////// For Voxelized Phantoms, External Dosimetry, Voxel Level  ////////////////////////////////////////////////////////////////////////////

/*   For Parammetrization 0
extern G4String* CopyNumberRegionNameMap;

const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
auto edep = step->GetTotalEnergyDeposit();
G4int CN = touchable->GetCopyNumber();
if(step->GetTrack()->GetTrackID() == 1){
    RunAction->FillVoxelLenghts(CN, 0.1*step->GetStepLength()); // in cm
}
RunAction->FillVoxelStepHits(touchable->GetCopyNumber(), edep);
*/


/*   For Nested Parammetrization 1
extern G4String* CopyNumberRegionNameMap;
extern G4int VoxXNumber;
extern G4int VoxYNumber;
extern G4int VoxZNumber;

const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
auto edep = step->GetTotalEnergyDeposit();
if(touchable->GetHistoryDepth()!=0){
    G4int CN = touchable->GetReplicaNumber(1) + VoxXNumber*touchable->GetReplicaNumber(2) + VoxXNumber*VoxYNumber*touchable->GetReplicaNumber(0);
    if(step->GetTrack()->GetTrackID() == 1){
        RunAction->FillVoxelLenghts(CN, 0.1*step->GetStepLength()); // in cm
    }
    RunAction->FillVoxelStepHits(CN, edep);
}
*/


//////////////////////// For Mesh-type Phantoms, Internal Dosimetry and External Dosimetry are already separated ////////////////////////////////////////////////////////////////////////////

