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
//


#include "G4TMyGeometry.hh"
#include "G4VPhysicalVolume.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialTable.hh"
#include "Randomize.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Ellipsoid.hh"
#include "G4UnionSolid.hh"
#include "G4EllipticalTube.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"

#include "G4TNuclearReactorGeometry.hh"
#include "G4TDetectorGeometry.hh"
#include "G4TPassiveProtonBeamLineGeometry.hh"
#include "G4TFlashRadioTherapyGeometry.hh"

G4TMyGeometry::G4TMyGeometry(){}
G4TMyGeometry::~G4TMyGeometry(){}

G4VPhysicalVolume* G4TMyGeometry::ConstructPhysicalVolume(){

    // //////////////////////////////// Required Material creation ////////////////////////


    //G4TDetectorGeometry* NewGeom = new G4TDetectorGeometry();
    //WorldPhysicalVolume = NewGeom->ConstructLXe();
    //WorldPhysicalVolume = NewGeom->ConstructNeutronDetector();

    //G4TPassiveProtonBeamLineGeometry* NewGeom = new G4TPassiveProtonBeamLineGeometry();
    //WorldPhysicalVolume = NewGeom->Construct();


    G4TFlashRadioTherapyGeometry* NewGeom = new G4TFlashRadioTherapyGeometry();
    WorldPhysicalVolume = NewGeom->ConstructGeometry();

    return WorldPhysicalVolume;

}

void G4TMyGeometry::ConstructLogicalVolumes(){

}
