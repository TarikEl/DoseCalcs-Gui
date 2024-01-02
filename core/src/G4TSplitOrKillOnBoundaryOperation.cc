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
/// \file G4TSplitOrKillOnBoundaryOperation.cc
/// \brief Implementation of the G4TSplitOrKillOnBoundaryOperation class

#include "G4TSplitOrKillOnBoundaryOperation.hh"
#include "G4BiasingProcessInterface.hh"

#include "G4ParticleChangeForLoss.hh"
#include "G4ParticleChangeForGamma.hh"


G4TSplitOrKillOnBoundaryOperation::G4TSplitOrKillOnBoundaryOperation(G4String name):G4VBiasingOperation(name), fSplittingFactor(1), fParticleChange(){}
G4TSplitOrKillOnBoundaryOperation::~G4TSplitOrKillOnBoundaryOperation(){}


G4VParticleChange* G4TSplitOrKillOnBoundaryOperation::ApplyFinalStateBiasing( const G4BiasingProcessInterface* callingProcess,const G4Track* track,const G4Step* step,G4bool&)
{
    //G4cout << __FUNCTION__<< G4endl;

    // -- Collect brem. process (wrapped process) final state:
    G4VParticleChange* processFinalState = callingProcess->GetWrappedProcess()->PostStepDoIt(*track, *step);

    if ( fSplittingFactor == 1 ) return processFinalState;

    G4double initialWeight = track->GetWeight();
    G4String VolFrom = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName();
    G4String VolTo = step->GetPostStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName();
    G4String BiasVol = "Adrenal";

    //G4cout << "Track weight " << initialWeight << " - " << VolFrom << " ---To---> " << VolTo << G4endl;

    if ( VolFrom != BiasVol && VolTo != BiasVol ) return processFinalState;

    fParticleChange.Initialize(*track);


    if ( VolFrom == BiasVol && VolTo != BiasVol ) {

        G4cout << " Killing ##################### " << G4endl;

        // -- We apply Russian roulette if the track is moving backward:
        G4double random = G4UniformRand();
        G4double killingProbability = 1.0 - 1.0/fSplittingFactor;
        if ( random < killingProbability )
        {
            // We ask for the the track to be killed:
            fParticleChange.ProposeTrackStatus(fStopAndKill);
        }
        else
        {
            // In this case, the track survives. We change its weight to conserve weight among killed and survival tracks:
            fParticleChange.ProposeParentWeight( initialWeight*fSplittingFactor );
        }

    }
    else if ( VolTo == BiasVol ) {

        G4cout << " Splitting ##################### " << G4endl;

        G4cout << "Track weight " << initialWeight << " - " << VolFrom << " ---To---> " << VolTo << G4endl;


        if ( processFinalState->GetNumberOfSecondaries() == 0 )  return processFinalState;

        //G4ParticleChangeForLoss* actualParticleChange = ( G4ParticleChangeForLoss* ) processFinalState ;
        G4ParticleChangeForGamma * actualParticleChange = ( G4ParticleChangeForGamma* ) processFinalState;

        //fParticleChange.Initialize(*track);

        // -- Store electron(which is the primary gamma in interaction in compt) final state:
        //fParticleChange.ProposeTrackStatus      ( actualParticleChange->GetTrackStatus() );
        //fParticleChange.ProposeEnergy           ( actualParticleChange->GetProposedKineticEnergy() );
        //fParticleChange.ProposeMomentumDirection( actualParticleChange->GetProposedMomentumDirection() );
        //fParticleChange.ProposeParentWeight( gammaWeight );

        //fParticleChange.ProposeTrackStatus      ( track->GetTrackStatus() );
        //fParticleChange.ProposeEnergy           ( track->GetTotalEnergy() );
        //fParticleChange.ProposeMomentumDirection( track->GetMomentumDirection() );

        fParticleChange.SetSecondaryWeightByProcess(true);
        G4double gammaWeight = initialWeight/(fSplittingFactor);

        fParticleChange.SetNumberOfSecondaries(fSplittingFactor+1);// electron of compt secondary and gammas

        G4Track* electronTrack = actualParticleChange->GetSecondary(0);
        electronTrack->SetWeight( initialWeight );
        fParticleChange.AddSecondary( electronTrack ); // Store first electron:
        actualParticleChange->Clear(); // -- and clean-up the brem. process particle change:

        G4cout << "Primary   " << track->GetParticleDefinition()->GetParticleName() << " E=" << actualParticleChange->GetCurrentTrack()->GetKineticEnergy() << " Pos=" << track->GetPosition() << " MomDir=" << track->GetMomentumDirection() <<" Weight=" << track->GetWeight()<< G4endl;
        G4cout << "secondary    " << electronTrack->GetParticleDefinition()->GetParticleName() << " E=" << electronTrack->GetKineticEnergy() << " Pos=" << electronTrack->GetPosition() << " MomDir=" << electronTrack->GetMomentumDirection() <<" Weight=" << electronTrack->GetWeight()<< G4endl;

        G4int ii = 1;
        G4int SecNum = fSplittingFactor+1;
        while ( ii < SecNum )
        {

            G4Track* clone = new G4Track( *actualParticleChange->GetCurrentTrack());
            clone->SetKineticEnergy(actualParticleChange->GetProposedKineticEnergy());
            clone->SetMomentumDirection(actualParticleChange->GetProposedMomentumDirection());
            clone->SetWeight( gammaWeight );
            G4cout << "secondary " << clone->GetParticleDefinition()->GetParticleName() << " E=" << clone->GetKineticEnergy() << " Pos=" << clone->GetPosition() << " MomDir=" << clone->GetMomentumDirection() <<" Weight=" << clone->GetWeight()<< G4endl;
            fParticleChange.AddSecondary( clone );
            ii++;
        }

        G4cout << "New Number Of Secondaries " <<fParticleChange.GetNumberOfSecondaries() << G4endl;
    }

    return &fParticleChange;
}

