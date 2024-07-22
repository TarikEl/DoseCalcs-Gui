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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#ifndef G4TPhantomCreationAdding_H
#define G4TPhantomCreationAdding_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

class G4VPhysicalVolume;
class G4THadrontherapyModulator;
//class G4TPhantomCreationAddingMessenger;
class G4THadrontherapyDetectorROGeometry;

class G4TPhantomCreationAdding
        //: public G4VUserDetectorConstruction
{
public:

    G4TPhantomCreationAdding();
    ~G4TPhantomCreationAdding();
    // static G4bool doCalculation;
    
    void CreateAndAddPhantom(G4VPhysicalVolume*, G4ThreeVector, G4ThreeVector);
    
private:

    G4VisAttributes* blue;
    G4VisAttributes* gray;
    G4VisAttributes* white;
    G4VisAttributes* red;
    G4VisAttributes* yellow;
    G4VisAttributes* green;
    G4VisAttributes* darkGreen;
    G4VisAttributes* darkOrange3;
    G4VisAttributes* skyBlue;

    
    G4Box *phantom , *detector;
    G4LogicalVolume *phantomLogicalVolume, *detectorLogicalVolume;
    G4VPhysicalVolume *phantomPhysicalVolume, *detectorPhysicalVolume;

    G4Box* solidVirtualLayer;
    G4LogicalVolume* logicVirtualLayer;
    G4VPhysicalVolume*  physVirtualLayer;

    G4double phantomSizeX;
    G4double phantomSizeY;
    G4double phantomSizeZ;

    G4double ROSizeX;
    G4double ROSizeY;
    G4double ROSizeZ;

    G4ThreeVector ROToWorldPosition, phantomPosition, ROPosition,
        ROToPhantomPosition; //  phantom center, detector center, detector
                                   //  to phantom relative position

    G4double sizeOfVoxelAlongX;
    G4double sizeOfVoxelAlongY;
    G4double sizeOfVoxelAlongZ;

    G4int numberOfVoxelsAlongX;
    G4int numberOfVoxelsAlongY;
    G4int numberOfVoxelsAlongZ;

    G4double volumeOfVoxel, massOfVoxel;

    G4Material *phantomMaterial, *ROMaterial;



    // RO detector

    G4Box* RODetector;
    G4Box* RODetectorXDivision;
    G4Box* RODetectorYDivision;
    G4Box* RODetectorZDivision;

    //Logical volumes used for the re-build on-the-fly
    G4LogicalVolume* RODetectorLog;
    G4LogicalVolume* RODetectorXDivisionLog;
    G4LogicalVolume* RODetectorYDivisionLog;
    G4LogicalVolume* RODetectorZDivisionLog;



    
};
#endif
