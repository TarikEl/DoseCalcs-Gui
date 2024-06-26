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
// TETSteppingActionTETExternal.hh
// \file   MRCP_GEANT4/External/include/TETSteppingActionTETExternal.hh
// \author Haegin Han
//

#ifndef TETSteppingActionTETExternal_h
#define TETSteppingActionTETExternal_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

class G4LogicalVolume;

// *********************************************************************
// With very low probability, because of the internal bug of G4TET, the
// particles can be stuck in the vertices of tetrahedrons. This
// UserSteppingAction class was written to slightly move these stuck
// particles.
// -- UserSteppingAction: Slightly move the stuck particles.
// *********************************************************************

class TETRunAction;

class TETSteppingActionTETExternal : public G4UserSteppingAction
{
  public:
    TETSteppingActionTETExternal(TETRunAction* runAction);
    virtual ~TETSteppingActionTETExternal();

    virtual void UserSteppingAction(const G4Step*);

  private:
    G4double kCarTolerance;
    G4int    stepCounter;
    G4bool   checkFlag;
    TETRunAction* RunAction;

};

#endif
