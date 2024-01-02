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
/// \file G4TGeometryBasedBiasingOperator.cc
/// \brief Implementation of the G4TGeometryBasedBiasingOperator class

#include "G4TGeometryBasedBiasingOperator.hh"
#include "G4TSplitOrKillOnBoundaryOperation.hh"

#include "G4BiasingProcessInterface.hh"
#include "G4GenericMessenger.hh"
//#include "G4ParticleTable.hh"


G4TGeometryBasedBiasingOperator::G4TGeometryBasedBiasingOperator():G4VBiasingOperator("ComptonGammaSplitting"), fSplittingFactor(50),fBiasPrimaryOnly(true), fBiasOnlyOnce(true){

    G4cout << __FUNCTION__<< G4endl;

    fBremSplittingOperation = new G4TSplitOrKillOnBoundaryOperation("ComptonGammaSplittingOperation");

    // -- Define messengers:
    // -- Splitting factor:
    fSplittingFactorMessenger = new G4GenericMessenger(this, "/GB04/biasing/","Biasing control" );
    G4GenericMessenger::Command& splittingFactorCmd = fSplittingFactorMessenger->DeclareProperty("setSplittingFactor", fSplittingFactor,
                                                                                                 "Define the brem. splitting factor." );
    splittingFactorCmd.SetStates(G4State_Idle);
    // -- Bias ony primary particle:
    fBiasPrimaryOnlyMessenger = new G4GenericMessenger(this, "/GB04/biasing/","Biasing control" );
    G4GenericMessenger::Command& biasPrimaryCmd = fBiasPrimaryOnlyMessenger->DeclareProperty("biasPrimaryOnly", fBiasPrimaryOnly,
                                                                                             "Chose if brem. splitting applies to primary particles only." );
    biasPrimaryCmd.SetStates(G4State_Idle);
    // -- Bias ony primary particle:
    fBiasOnlyOnceMessenger = new G4GenericMessenger(this, "/GB04/biasing/","Biasing control" );
    G4GenericMessenger::Command& biasOnlyOnceCmd = fBiasPrimaryOnlyMessenger->DeclareProperty("biasOnlyOnce", fBiasOnlyOnce,
                                                                                              "Chose if apply the brem. splitting only once for the track." );
    biasOnlyOnceCmd.SetStates(G4State_Idle);

}


void G4TGeometryBasedBiasingOperator::StartRun()
{
    G4cout << __FUNCTION__<< G4endl;

    fBremSplittingOperation->SetSplittingFactor ( fSplittingFactor );
    G4cout << GetName() << " : starting run with compt splitting factor = " << fSplittingFactor << G4endl;

}


void G4TGeometryBasedBiasingOperator::StartTracking( const G4Track* /* track */ )
{
    //G4cout << __FUNCTION__<< G4endl;

    fNInteractions = 0;

}


G4VBiasingOperation* G4TGeometryBasedBiasingOperator::ProposeFinalStateBiasingOperation(const G4Track* track , const G4BiasingProcessInterface* /* callingProcess */)
{
    //G4cout << __FUNCTION__<< G4endl;

    // -- Check if biasing of primary particle only is requested. If so, and if particle is not a primary one, don't ask for biasing:
    if ( fBiasPrimaryOnly && ( track->GetParentID() !=0 ) ) return 0;

    // -- Check if brem. splitting should be applied only once to the track, and if so, and if brem. splitting already occured, don't ask for biasing:
    //if ( fBiasOnlyOnce && ( fNInteractions > 0 )) return 0;

    //fNInteractions++;

    return fBremSplittingOperation;
}

