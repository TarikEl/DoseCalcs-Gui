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

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"



#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4PhysicalVolumeStore.hh"

#include "G4ThreeVector.hh"

#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"

#include "G4THadrontherapyDetectorROGeometry.hh"
//#include "HadrontherapyDummySD.hh"
//#include "HadrontherapyDetectorSD.hh"

/////////////////////////////////////////////////////////////////////////////
G4THadrontherapyDetectorROGeometry::G4THadrontherapyDetectorROGeometry(G4VPhysicalVolume* PV)
//: G4VUserParallelWorld(aString),
    : WorldPhysicalVolume(PV), RODetector(0),RODetectorXDivision(0),
      RODetectorYDivision(0),RODetectorZDivision(0),worldLogical(0),RODetectorLog(0),
      RODetectorXDivisionLog(0),RODetectorYDivisionLog(0),RODetectorZDivisionLog(0),
      sensitiveLogicalVolume(0)
{
    isBuilt = false;
    isInitialized = false;
}

/////////////////////////////////////////////////////////////////////////////

void G4THadrontherapyDetectorROGeometry::Initialize(G4ThreeVector pos,
                                                    G4double detectorDimX,
                                                    G4double detectorDimY,
                                                    G4double detectorDimZ,
                                                    G4int numberOfVoxelsX,
                                                    G4int numberOfVoxelsY,
                                                    G4int numberOfVoxelsZ)
{
    detectorToWorldPosition = pos;
    detectorSizeX = detectorDimX;
    detectorSizeY= detectorDimY;
    detectorSizeZ=detectorDimZ;
    numberOfVoxelsAlongX=numberOfVoxelsX;
    numberOfVoxelsAlongY=numberOfVoxelsY;
    numberOfVoxelsAlongZ=numberOfVoxelsZ;

    isInitialized = true;
}

void G4THadrontherapyDetectorROGeometry::UpdateROGeometry()
{
    //Nothing happens if the RO geometry is not built. But parameters are properly set.
    if (!isBuilt)
    {
        //G4Exception("G4THadrontherapyDetectorROGeometry::UpdateROGeometry","had001",
        //		  JustWarning,"Cannot update geometry before it is built");
        return;
    }

    //1) Update the dimensions of the G4Boxes
    G4double halfDetectorSizeX = detectorSizeX;
    G4double halfDetectorSizeY = detectorSizeY;
    G4double halfDetectorSizeZ = detectorSizeZ;

    RODetector->SetXHalfLength(halfDetectorSizeX);
    RODetector->SetYHalfLength(halfDetectorSizeY);
    RODetector->SetZHalfLength(halfDetectorSizeZ);

    G4double halfXVoxelSizeX = halfDetectorSizeX/numberOfVoxelsAlongX;
    G4double halfXVoxelSizeY = halfDetectorSizeY;
    G4double halfXVoxelSizeZ = halfDetectorSizeZ;
    G4double voxelXThickness = 2*halfXVoxelSizeX;

    RODetectorXDivision->SetXHalfLength(halfXVoxelSizeX);
    RODetectorXDivision->SetYHalfLength(halfXVoxelSizeY);
    RODetectorXDivision->SetZHalfLength(halfXVoxelSizeZ);

    G4double halfYVoxelSizeX = halfXVoxelSizeX;
    G4double halfYVoxelSizeY = halfDetectorSizeY/numberOfVoxelsAlongY;
    G4double halfYVoxelSizeZ = halfDetectorSizeZ;
    G4double voxelYThickness = 2*halfYVoxelSizeY;

    RODetectorYDivision->SetXHalfLength(halfYVoxelSizeX);
    RODetectorYDivision->SetYHalfLength(halfYVoxelSizeY);
    RODetectorYDivision->SetZHalfLength(halfYVoxelSizeZ);

    G4double halfZVoxelSizeX = halfXVoxelSizeX;
    G4double halfZVoxelSizeY = halfYVoxelSizeY;
    G4double halfZVoxelSizeZ = halfDetectorSizeZ/numberOfVoxelsAlongZ;
    G4double voxelZThickness = 2*halfZVoxelSizeZ;

    RODetectorZDivision->SetXHalfLength(halfZVoxelSizeX);
    RODetectorZDivision->SetYHalfLength(halfZVoxelSizeY);
    RODetectorZDivision->SetZHalfLength(halfZVoxelSizeZ);

    //Delete and re-build the relevant physical volumes
    G4PhysicalVolumeStore* store =
            G4PhysicalVolumeStore::GetInstance();

    //Delete...
    G4VPhysicalVolume* myVol = store->GetVolume("RODetectorPhys");
    store->DeRegister(myVol);
    //..and rebuild
    G4VPhysicalVolume *RODetectorPhys = new G4PVPlacement(0,
                                                          detectorToWorldPosition,
                                                          RODetectorLog,
                                                          "RODetectorPhys",
                                                          worldLogical,
                                                          false,0);

    myVol = store->GetVolume("RODetectorXDivisionPhys");
    store->DeRegister(myVol);
    G4VPhysicalVolume *RODetectorXDivisionPhys = new G4PVReplica("RODetectorXDivisionPhys",
                                                                 RODetectorXDivisionLog,
                                                                 RODetectorPhys,
                                                                 kXAxis,
                                                                 numberOfVoxelsAlongX,
                                                                 voxelXThickness);
    myVol = store->GetVolume("RODetectorYDivisionPhys");
    store->DeRegister(myVol);
    G4VPhysicalVolume *RODetectorYDivisionPhys = new G4PVReplica("RODetectorYDivisionPhys",
                                                                 RODetectorYDivisionLog,
                                                                 RODetectorXDivisionPhys,
                                                                 kYAxis,
                                                                 numberOfVoxelsAlongY,
                                                                 voxelYThickness);

    myVol = store->GetVolume("RODetectorZDivisionPhys");
    store->DeRegister(myVol);
    new G4PVReplica("RODetectorZDivisionPhys",
                    RODetectorZDivisionLog,
                    RODetectorYDivisionPhys,
                    kZAxis,
                    numberOfVoxelsAlongZ,
                    voxelZThickness);

    return;

}

/////////////////////////////////////////////////////////////////////////////
G4THadrontherapyDetectorROGeometry::~G4THadrontherapyDetectorROGeometry()
{;}

/////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* G4THadrontherapyDetectorROGeometry::Construct()
{
    // A dummy material is used to fill the volumes of the readout geometry.
    // (It will be allowed to set a NULL pointer in volumes of such virtual
    // division in future, since this material is irrelevant for tracking.)


    //
    // World
    //
    //G4VPhysicalVolume* ghostWorld = GetWorld();
    worldLogical = WorldPhysicalVolume->GetLogicalVolume();

    if (!isInitialized)
    {
        G4Exception("G4THadrontherapyDetectorROGeometry::Construct","had001",
                    FatalException,"Parameters of the RO geometry are not initialized");
        return WorldPhysicalVolume;
    }




    G4double halfDetectorSizeX = detectorSizeX;
    G4double halfDetectorSizeY = detectorSizeY;
    G4double halfDetectorSizeZ = detectorSizeZ;
    
    // World volume of ROGeometry ... SERVE SOLO PER LA ROG

    // Detector ROGeometry
    RODetector = new G4Box("RODetector",
                           halfDetectorSizeX,
                           halfDetectorSizeY,
                           halfDetectorSizeZ);

    G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    RODetectorLog = new G4LogicalVolume(RODetector,
                                        airNist,
                                        "RODetectorLog",
                                        0,0,0);



    G4VPhysicalVolume *RODetectorPhys = new G4PVPlacement(0,
                                                          detectorToWorldPosition,RODetectorLog,
                                                          "RODetectorPhys",
                                                          worldLogical,
                                                          false,0);




    // Division along X axis: the detector is divided in slices along the X axis

    G4double halfXVoxelSizeX = halfDetectorSizeX/numberOfVoxelsAlongX;
    G4double halfXVoxelSizeY = halfDetectorSizeY;
    G4double halfXVoxelSizeZ = halfDetectorSizeZ;
    G4double voxelXThickness = 2*halfXVoxelSizeX;

    std::cout << numberOfVoxelsAlongX << " " << halfXVoxelSizeX << std::endl;

    RODetectorXDivision = new G4Box("RODetectorXDivision",
                                    halfXVoxelSizeX,
                                    halfXVoxelSizeY,
                                    halfXVoxelSizeZ);

    RODetectorXDivisionLog = new G4LogicalVolume(RODetectorXDivision,
                                                 airNist,
                                                 "RODetectorXDivisionLog",
                                                 0,0,0);

    G4VPhysicalVolume *RODetectorXDivisionPhys = new G4PVReplica("RODetectorXDivisionPhys",
                                                                 RODetectorXDivisionLog,
                                                                 RODetectorPhys,
                                                                 kXAxis,
                                                                 numberOfVoxelsAlongX,
                                                                 voxelXThickness);

    // Division along Y axis: the slices along the X axis are divided along the Y axis

    G4double halfYVoxelSizeX = halfXVoxelSizeX;
    G4double halfYVoxelSizeY = halfDetectorSizeY/numberOfVoxelsAlongY;
    G4double halfYVoxelSizeZ = halfDetectorSizeZ;
    G4double voxelYThickness = 2*halfYVoxelSizeY;

    std::cout << numberOfVoxelsAlongY << " " << halfXVoxelSizeY << std::endl;

    RODetectorYDivision = new G4Box("RODetectorYDivision",
                                    halfYVoxelSizeX,
                                    halfYVoxelSizeY,
                                    halfYVoxelSizeZ);

    RODetectorYDivisionLog = new G4LogicalVolume(RODetectorYDivision,
                                                 airNist,
                                                 "RODetectorYDivisionLog",
                                                 0,0,0);

    G4VPhysicalVolume *RODetectorYDivisionPhys = new G4PVReplica("RODetectorYDivisionPhys",
                                                                 RODetectorYDivisionLog,
                                                                 RODetectorXDivisionPhys,
                                                                 kYAxis,
                                                                 numberOfVoxelsAlongY,
                                                                 voxelYThickness);

    // Division along Z axis: the slices along the Y axis are divided along the Z axis

    G4double halfZVoxelSizeX = halfXVoxelSizeX;
    G4double halfZVoxelSizeY = halfYVoxelSizeY;
    G4double halfZVoxelSizeZ = halfDetectorSizeZ/numberOfVoxelsAlongZ;
    G4double voxelZThickness = 2*halfZVoxelSizeZ;

    std::cout << numberOfVoxelsAlongZ << " " << halfXVoxelSizeZ << std::endl;

    RODetectorZDivision = new G4Box("RODetectorZDivision",
                                    halfZVoxelSizeX,
                                    halfZVoxelSizeY,
                                    halfZVoxelSizeZ);

    RODetectorZDivisionLog = new G4LogicalVolume(RODetectorZDivision,
                                                 airNist,
                                                 "RODetectorZDivisionLog",
                                                 0,0,0);

    new G4PVReplica("RODetectorZDivisionPhys",
                    RODetectorZDivisionLog,
                    RODetectorYDivisionPhys,
                    kZAxis,
                    numberOfVoxelsAlongZ,
                    voxelZThickness);

    //sensitiveLogicalVolume = RODetectorZDivisionLog;
    isBuilt = true;

    return WorldPhysicalVolume;
}

