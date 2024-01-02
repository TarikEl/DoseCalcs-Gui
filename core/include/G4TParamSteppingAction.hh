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
/// \file G4TParamSteppingAction.hh
/// \brief Definition of the G4TParamSteppingAction class

#ifndef G4TParamSteppingAction_h
#define G4TParamSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

//class G4TVolumeConstruction;
class G4TRunAction;

class G4TParamSteppingAction : public G4UserSteppingAction
{
public:
    G4TParamSteppingAction(G4TRunAction* runAction);
    virtual ~G4TParamSteppingAction();

    virtual void UserSteppingAction(const G4Step* step);
    
private:
    //const G4TVolumeConstruction* fDetConstruction;
    G4TRunAction*  RunAction;

    //G4ThreadLocal static G4double edep;
    //G4ThreadLocal static G4String reg;
    //G4ThreadLocal static unsigned int CN;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
