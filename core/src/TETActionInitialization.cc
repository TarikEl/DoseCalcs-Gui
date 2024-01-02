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
// TETActionInitialization.cc
// \file   MRCP_GEANT4/Internal/src/TETActionInitialization.cc
// \author Haegin Han
//

#include "G4TVolumeConstruction.hh"
#include "TETActionInitialization.hh"

TETActionInitialization::TETActionInitialization()
 : G4VUserActionInitialization()
{}

TETActionInitialization::~TETActionInitialization()
{}

void TETActionInitialization::BuildForMaster() const
{
    SetUserAction(new TETRunAction());
}

void TETActionInitialization::Build() const
{
	// initialise UserAction classes
    SetUserAction(new TETRunAction());
	SetUserAction(new TETSteppingAction);

    const G4TVolumeConstruction* TConstruction2 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    if( TConstruction2->getUseGeneratedData() == "read"){
        SetUserAction(new G4TReadPrimaryGeneratorAction);
    }else if(TConstruction2->getUseGeneratedData() == "save"){
        SetUserAction(new G4TDirectToFilesPrimaryGeneratorAction);
    }else{
        SetUserAction(new G4TDirectPrimaryGeneratorAction);
    }

}  
