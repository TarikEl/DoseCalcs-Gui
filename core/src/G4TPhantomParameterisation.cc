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
// Author: Tarik Elghalbzouri,  Abdelmalek Essaâdi University,
// faculty of sciences Tetouane, morocco. email : telghalbzouri@uae.ac.ma
//
// This application is based on code developed by :
// G. Guerrieri, University of Genova, Italy .
// S. Guatelli. University of Wollongong, Australia.
//


#include "G4TPhantomParameterisation.hh"

#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "Randomize.hh"
#include "CLHEP/Random/RandFlat.h"

#include <stdio.h>

// it call ReadColourData()
G4TPhantomParameterisation::G4TPhantomParameterisation() : G4PhantomParameterisation()
{
#if VERBOSE_USE
    G4cout  << "\n" << __FUNCTION__ << " with SetSkipEqualMaterials option " << G4endl;
    G4cout  << "UseLogVolColour" << " -- " << UseLogVolColour << G4endl;
#endif
    SetSkipEqualMaterials(false);
}


G4TPhantomParameterisation::~G4TPhantomParameterisation(){}


// called automatically when the creation of voxel of copyNo when parametrisation it's use just to SetVisAttributes materials are already filled the voxels
G4Material* G4TPhantomParameterisation::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume * physVol, const G4VTouchable *){

    //if(parentTouch == nullptr) return fMaterials[0];
    //std::size_t matIndex = GetMaterialIndex(copyNo);
    //static G4Material* mate = nullptr;
    //mate = fMaterials[matIndex];

    G4Material* mate = G4PhantomParameterisation::ComputeMaterial( copyNo , physVol , 0 );


    //G4cout  << "ComputeMaterial : copyNo " << " -- " << copyNo << " " << mate->GetName() << G4endl;

    if(UseLogVolColour == true) {

        //G4cout  << UseLogVolColour << " -- " << RegionCopyNumberColour[copyNo] << G4endl;

        if(physVol){
            G4VisAttributes* voxvis = new G4VisAttributes(RegionNameColour[CopyNumberRegionNameMap[copyNo]]);
            voxvis->SetVisibility(true); voxvis->SetForceSolid(ForceSolid);
            physVol->GetLogicalVolume()->SetVisAttributes(voxvis);
        }
        return mate;
    }
    else {
        return mate;
    }
}
