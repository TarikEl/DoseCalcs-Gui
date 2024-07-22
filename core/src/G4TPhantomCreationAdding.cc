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

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "G4NistElementBuilder.hh"
#include "G4THadrontherapyModulator.hh"
#include "G4TPhantomCreationAdding.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalVolumeStore.hh"

extern G4String* CopyNumberRegionNameMap;
extern G4double VoxXHalfSize;
extern G4double VoxYHalfSize;
extern G4double VoxZHalfSize;
extern G4int VoxXNumber;
extern G4int VoxYNumber;
extern G4int VoxZNumber;
extern  size_t* MateIDs; // index of material of each voxel unsigned int* fMateIDs; // index of material of each voxel
extern  G4String* CopyNumberRegionNameMap;
extern  G4float* CopyNumberXPos;
extern  G4float* CopyNumberYPos;
extern  G4float* CopyNumberZPos;
extern  G4float* CopyNumberMassSize;

G4TPhantomCreationAdding::G4TPhantomCreationAdding(){}
G4TPhantomCreationAdding::~G4TPhantomCreationAdding(){}

void G4TPhantomCreationAdding::CreateAndAddPhantom(G4VPhysicalVolume* motherPhyVol, G4ThreeVector PhantomPos, G4ThreeVector PhantomSize){

    // ***********************************  Phantom Geometry

    // Set of coulors that can be used
    white = new G4VisAttributes( G4Colour());
    white -> SetVisibility(true);
    white -> SetForceSolid(true);

    blue = new G4VisAttributes(G4Colour(0. ,0. ,1.));
    blue -> SetVisibility(true);
    blue -> SetForceSolid(true);

    gray = new G4VisAttributes( G4Colour(0.5, 0.5, 0.5 ));
    gray-> SetVisibility(true);
    gray-> SetForceSolid(true);

    red = new G4VisAttributes(G4Colour(1. ,0. ,0.));
    red-> SetVisibility(true);
    red-> SetForceSolid(true);

    yellow = new G4VisAttributes(G4Colour(1., 1., 0. ));
    yellow-> SetVisibility(true);
    yellow-> SetForceSolid(true);

    green = new G4VisAttributes( G4Colour(25/255. , 255/255. ,  25/255. ));
    green -> SetVisibility(true);
    green -> SetForceSolid(true);

    darkGreen = new G4VisAttributes( G4Colour(0/255. , 100/255. ,  0/255. ));
    darkGreen -> SetVisibility(true);
    darkGreen -> SetForceSolid(true);

    darkOrange3 = new G4VisAttributes( G4Colour(205/255. , 102/255. ,  000/255. ));
    darkOrange3 -> SetVisibility(true);
    darkOrange3 -> SetForceSolid(true);

    skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
    skyBlue -> SetVisibility(true);
    skyBlue -> SetForceSolid(true);

    // ELEMENTS
    G4bool isotopes = false;
    G4Material* aluminumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al", isotopes);
    G4Material* tantalumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ta", isotopes);
    G4Material* copperNistAsMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu", isotopes);
    G4Element* zincNist = G4NistManager::Instance()->FindOrBuildElement("Zn");
    G4Element* copperNist = G4NistManager::Instance()->FindOrBuildElement("Cu");

    // COMPOUND
    G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
    G4Material* kaptonNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON", isotopes);
    G4Material* galacticNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic", isotopes);
    G4Material* PMMANist = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);
    G4Material* mylarNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR", isotopes);

    G4double d; // Density
    G4int nComponents;// Number of components
    G4double fractionmass; // Fraction in mass of an element in a material

    d = 8.40*g/cm3;
    nComponents = 2;
    G4Material* brass = new G4Material("Brass", d, nComponents);
    brass -> AddElement(zincNist, fractionmass = 30 *perCent);
    brass -> AddElement(copperNist, fractionmass = 70 *perCent);


    //***************************** PW ***************************************

    // DetectorROGeometry Material
    new G4Material("dummyMat", 1., 1.*g/mole, 1.*g/cm3);

    // Define here the material of the water phantom and of the detector

    G4Material* pMat = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", false);
    phantomMaterial  = pMat;

    //phantomSizeX = 40. *cm;
    //phantomSizeY = 40. *cm;
    //phantomSizeZ = 40. *cm;

    phantomSizeX = PhantomSize.getX();
    phantomSizeY = PhantomSize.getY();
    phantomSizeZ = PhantomSize.getZ();

    //phantomPosition = G4ThreeVector(20. *cm, 0. *cm, 0. *cm);
    phantomPosition = PhantomPos;




    // Definition of the solid volume of the Phantom
    phantom = new G4Box("Phantom",
                        phantomSizeX/2,
                        phantomSizeY/2,
                        phantomSizeZ/2);

    // Definition of the logical volume of the Phantom
    phantomLogicalVolume = new G4LogicalVolume(phantom,
                                               phantomMaterial,
                                               "phantomLog", 0, 0, 0);

    // Definition of the physics volume of the Phantom
    phantomPhysicalVolume = new G4PVPlacement(0,
                                              phantomPosition,
                                              "phantomPhys",
                                              phantomLogicalVolume,
                                              motherPhyVol,
                                              false,
                                              0);

    // Visualisation attributes of the phantom
    red = new G4VisAttributes(G4Colour(255/255., 0/255. ,0/255.));
    red -> SetVisibility(true);
    red -> SetForceSolid(true);
    red -> SetForceWireframe(true);
    phantomLogicalVolume -> SetVisAttributes(red);


    // ************************************ Voxelized tumour Region


    // World volume of ROGeometry ... SERVE SOLO PER LA ROG

    ROSizeX = 4. *cm;
    ROSizeY = 4. *cm;
    ROSizeZ = 4. *cm;


    // Adjust RO position
    ROToPhantomPosition = G4ThreeVector(-5. *cm, 0. *cm, 0. *cm);

    ROPosition.setX(ROToPhantomPosition.getX() - phantomSizeX/2 + ROSizeX/2);
    ROPosition.setY(ROToPhantomPosition.getY() - phantomSizeY/2 + ROSizeY/2);
    ROPosition.setZ(ROToPhantomPosition.getZ() - phantomSizeZ/2 + ROSizeZ/2);

    ROToWorldPosition = phantomPosition + ROPosition;

    ROMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", false);

    // Default detector RO voxels size
    // 200 slabs along the beam direction (X)

    // // by setting the half sizes it calculate the number of voxels
    //VoxXHalfSize = 200 *um;
    //VoxYHalfSize = 100 *um;
    //VoxZHalfSize = 80  *um;
    //VoxXNumber = G4lrint(ROSizeX / VoxXHalfSize);
    //VoxYNumber = G4lrint(ROSizeY / VoxYHalfSize);
    //VoxZNumber = G4lrint(ROSizeZ / VoxZHalfSize);

    // // by setting the number of voxels it calculate the half sizes
    VoxXNumber = 20;
    VoxYNumber = 1;
    VoxZNumber = 1;
    VoxXHalfSize = ROSizeX/VoxXNumber;
    VoxYHalfSize = ROSizeY/VoxYNumber;
    VoxZHalfSize = ROSizeZ/VoxZNumber;


    // Detector ROGeometry
    RODetector = new G4Box("RODetector",
                           ROSizeX,
                           ROSizeY,
                           ROSizeZ);

    RODetectorLog = new G4LogicalVolume(RODetector,
                                        ROMaterial,
                                        "RODetectorLog",
                                        0,0,0);

    G4VPhysicalVolume *RODetectorPhys = new G4PVPlacement(0,
                                                          ROToPhantomPosition,
                                                          //ROToWorldPosition,
                                                          RODetectorLog,
                                                          "RODetectorPhys",
                                                          phantomLogicalVolume,
                                                          false,0);

    // Division along X axis: the detector is divided in slices along the X axis

    G4double halfXVoxelSizeX = ROSizeX/VoxXNumber;
    G4double halfXVoxelSizeY = ROSizeY;
    G4double halfXVoxelSizeZ = ROSizeZ;
    G4double voxelXThickness = 2*halfXVoxelSizeX;

    std::cout << VoxXNumber << " " << halfXVoxelSizeX << std::endl;

    RODetectorXDivision = new G4Box("RODetectorXDivision",
                                    halfXVoxelSizeX,
                                    halfXVoxelSizeY,
                                    halfXVoxelSizeZ);

    RODetectorXDivisionLog = new G4LogicalVolume(RODetectorXDivision,
                                                 ROMaterial,
                                                 "RO",
                                                 0,0,0);

    G4VPhysicalVolume *RODetectorXDivisionPhys = new G4PVReplica("RODetectorXDivisionPhys",
                                                                 RODetectorXDivisionLog,
                                                                 RODetectorLog,
                                                                 kXAxis,
                                                                 VoxXNumber,
                                                                 voxelXThickness);

    // Division along Y axis: the slices along the X axis are divided along the Y axis

    G4double halfYVoxelSizeX = halfXVoxelSizeX;
    G4double halfYVoxelSizeY = ROSizeY/VoxYNumber;
    G4double halfYVoxelSizeZ = ROSizeZ;
    G4double voxelYThickness = 2*halfYVoxelSizeY;

    std::cout << VoxYNumber << " " << halfXVoxelSizeY << std::endl;

    RODetectorYDivision = new G4Box("RODetectorYDivision",
                                    halfYVoxelSizeX,
                                    halfYVoxelSizeY,
                                    halfYVoxelSizeZ);

    RODetectorYDivisionLog = new G4LogicalVolume(RODetectorYDivision,
                                                 ROMaterial,
                                                 "RO",
                                                 0,0,0);

    G4VPhysicalVolume *RODetectorYDivisionPhys = new G4PVReplica("RODetectorYDivisionPhys",
                                                                 RODetectorYDivisionLog,
                                                                 RODetectorXDivisionLog,
                                                                 kYAxis,
                                                                 VoxYNumber,
                                                                 voxelYThickness);

    // Division along Z axis: the slices along the Y axis are divided along the Z axis

    G4double halfZVoxelSizeX = halfXVoxelSizeX;
    G4double halfZVoxelSizeY = halfYVoxelSizeY;
    G4double halfZVoxelSizeZ = ROSizeZ/VoxZNumber;
    G4double voxelZThickness = 2*halfZVoxelSizeZ;

    std::cout << VoxZNumber << " " << halfXVoxelSizeZ << std::endl;

    RODetectorZDivision = new G4Box("RODetectorZDivision",
                                    halfZVoxelSizeX,
                                    halfZVoxelSizeY,
                                    halfZVoxelSizeZ);

    RODetectorZDivisionLog = new G4LogicalVolume(RODetectorZDivision,
                                                 ROMaterial,
                                                 "RO",
                                                 0,0,0);

    new G4PVReplica("RODetectorZDivisionPhys",
                    RODetectorZDivisionLog,
                    RODetectorYDivisionLog,
                    kZAxis,
                    VoxZNumber,
                    voxelZThickness);

    size_t VoxelNumberInPhantom = VoxXNumber*VoxYNumber*VoxZNumber;

    CopyNumberXPos = new G4float[VoxelNumberInPhantom];
    CopyNumberYPos = new G4float[VoxelNumberInPhantom];
    CopyNumberZPos = new G4float[VoxelNumberInPhantom];

    CopyNumberMassSize = new G4float[VoxelNumberInPhantom];
    CopyNumberRegionNameMap = new G4String[VoxelNumberInPhantom];

    G4double VoxelVolume = 8*(VoxXHalfSize*VoxYHalfSize*VoxZHalfSize)/mm;
#if VERBOSE_USE
    G4cout << " VoxelVolume (mm3) of RO: " << VoxelVolume << G4endl;
#endif
    // CN coordinate and position
    G4float Voxel0PosX = -((VoxXNumber*VoxXHalfSize) - VoxXHalfSize) + ROToWorldPosition.getX();
    G4float Voxel0PosY = -((VoxYNumber*VoxYHalfSize) - VoxYHalfSize) + ROToWorldPosition.getY();
    G4float Voxel0PosZ = -((VoxZNumber*VoxZHalfSize) - VoxZHalfSize) + ROToWorldPosition.getZ();

    size_t Cn_inc = 0;
    for(size_t f = 0; f < VoxZNumber ;f++ ){
        for(size_t g = 0; g < VoxYNumber ;g++ ){
            for(size_t d = 0; d < VoxXNumber ;d++ ){

                CopyNumberXPos[Cn_inc] = Voxel0PosX + d * 2 * VoxXHalfSize;
                CopyNumberYPos[Cn_inc] = Voxel0PosY + g * 2 * VoxYHalfSize;
                CopyNumberZPos[Cn_inc] = Voxel0PosZ + f * 2 * VoxZHalfSize;
                //G4cout << Cn_inc << " " << f << " " << g << " " << d << G4endl;

                CopyNumberMassSize[Cn_inc]= ROMaterial->GetDensity() * (mm3/kg) * VoxelVolume; /*density ; e+21 in g/mm3 , and e+18 g/cm3 )*/
                CopyNumberRegionNameMap[Cn_inc] = "RO";

                //G4float ff = VoxelsMaterials[MateIDs[Cn_inc]]->GetDensity() * G4Density_to_kgPerMm3 * VoxelVolume ;
                //G4double nn = VoxelsMaterials[MateIDs[Cn_inc]]->GetDensity() * G4Density_to_kgPerMm3 * VoxelVolume ;
                //G4cout << " Float " << ff << "  Vs Double " << nn  << G4endl;

                Cn_inc++;
            }
        }
    }

}


