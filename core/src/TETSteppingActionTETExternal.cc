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
// TETSteppingActionTETExternal.cc
// \file   MRCP_GEANT4/Internal/src/TETSteppingActionTETExternal.cc
// \author Haegin Han
//

#include "TETSteppingActionTETExternal.hh"
#include "TETRunAction.hh"

extern G4String* CopyNumberRegionNameMap;

TETSteppingActionTETExternal::TETSteppingActionTETExternal(TETRunAction* ra) : G4UserSteppingAction(), RunAction(ra), kCarTolerance(1.0000000000000002e-07), stepCounter(0), checkFlag(0)
{}

TETSteppingActionTETExternal::~TETSteppingActionTETExternal()
{}

void TETSteppingActionTETExternal::UserSteppingAction(const G4Step* step)
{
    // Slightly move the particle when the step length of five continuous steps is
    // shorter than the tolerance (0.1 nm)
    //
    G4Track* theTrack = step->GetTrack();
    G4bool CheckingLength = (step->GetStepLength() < kCarTolerance);
    if(CheckingLength)
    {
        ++stepCounter;
        if( checkFlag && stepCounter>=5 )
        {
            // kill the track if the particle is stuck even after the slight move
            // (this hardly occurs)
            theTrack->SetTrackStatus(fStopAndKill);
            stepCounter=0;
            checkFlag=0;
        }
        else if( stepCounter>=5 )
        {
            // if a particle is at the same position (step length < 0.1 nm) for five consecutive steps,
            // slightly move (0.1 nm) the stuck particle in the direction of momentum
            theTrack->SetPosition(theTrack->GetPosition() + theTrack->GetMomentumDirection()*kCarTolerance);
            checkFlag=1;
        }
    }
    else stepCounter=0;

    G4int CN = step->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
    if(step->GetTrack()->GetTrackID() == 1){
        //std::cout <<" edep:" << edep << " SL " << step->GetStepLength() << std::endl;
        RunAction->FillRegionLenghts(CopyNumberRegionNameMap[CN], 0.1*step->GetStepLength()); // in cm
        //std::cout << " Particle:" << step->GetTrack()->GetParticleDefinition()->GetParticleName() << " ID:" << step->GetTrack()->GetTrackID() << " reg:" << reg << " Edep:" << step->GetTotalEnergyDeposit() << " Lenght:" << step->GetStepLength() << " @@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
    }
}
