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
/// \file G4TSplitOrKillOnBoundaryOperation.hh
/// \brief Definition of the G4TSplitOrKillOnBoundaryOperation class
//
//---------------------------------------------------------------
//
// G4TSplitOrKillOnBoundaryOperation
//
// Class Description:
//    A G4VBiasingOperation to change a process cross-section.
//
//
//---------------------------------------------------------------
//   Initial version                         Nov. 2014 M. Verderi

#ifndef G4TSplitOrKillOnBoundaryOperation_hh
#define G4TSplitOrKillOnBoundaryOperation_hh 1

#include "G4VBiasingOperation.hh"
#include "G4ParticleChange.hh"

class G4TSplitOrKillOnBoundaryOperation : public G4VBiasingOperation {
public:
    // -- Constructor :
    G4TSplitOrKillOnBoundaryOperation(G4String name);
    // -- destructor:
    virtual ~G4TSplitOrKillOnBoundaryOperation();

public:
    // ----------------------------------------------
    // -- Methods from G4VBiasingOperation interface:
    // ----------------------------------------------

    // --Used:
    virtual G4VParticleChange*   ApplyFinalStateBiasing( const G4BiasingProcessInterface*, const G4Track*, const G4Step*, G4bool&);

    // -- Unsued:
    virtual G4double           DistanceToApplyOperation( const G4Track*, G4double, G4ForceCondition*) {return DBL_MAX;}
    virtual G4VParticleChange* GenerateBiasingFinalState( const G4Track*, const G4Step*   ) {return 0;}
    virtual const G4VBiasingInteractionLaw* ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*, G4ForceCondition& ){ return 0; }

public:
    // ----------------------------------------------
    // -- Additional methods, specific to this class:
    // ----------------------------------------------
    // -- Splitting factor:
    void     SetSplittingFactor(G4int splittingFactor) { fSplittingFactor = splittingFactor; }
    G4int GetSplittingFactor() const { return fSplittingFactor; }

private:
    G4int            fSplittingFactor;
    G4ParticleChange fParticleChange;
};

#endif
