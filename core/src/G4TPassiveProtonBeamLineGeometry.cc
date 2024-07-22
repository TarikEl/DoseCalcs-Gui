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
#include "G4TPassiveProtonBeamLineGeometry.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4TPhantomCreationAdding.hh"

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
void G4TPassiveProtonBeamLineGeometry::ConstructVoxelizedRO(){


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

//G4bool G4TPassiveProtonBeamLineGeometry::doCalculation = false;
/////////////////////////////////////////////////////////////////////////////
G4TPassiveProtonBeamLineGeometry::G4TPassiveProtonBeamLineGeometry():
    physicalTreatmentRoom(0),
    physiBeamLineSupport(0), physiBeamLineCover(0), physiBeamLineCover2(0),
    firstScatteringFoil(0), physiFirstScatteringFoil(0), physiKaptonWindow(0),
    solidStopper(0), physiStopper(0), secondScatteringFoil(0), physiSecondScatteringFoil(0),
    physiFirstCollimator(0), solidRangeShifterBox(0), logicRangeShifterBox(0),
    physiRangeShifterBox(0), physiSecondCollimator(0), physiFirstCollimatorModulatorBox(0),
    physiHoleFirstCollimatorModulatorBox(0), physiSecondCollimatorModulatorBox(0),
    physiHoleSecondCollimatorModulatorBox(0), physiMOPIMotherVolume(0),
    physiFirstMonitorLayer1(0), physiFirstMonitorLayer2(0), physiFirstMonitorLayer3(0),
    physiFirstMonitorLayer4(0), physiSecondMonitorLayer1(0), physiSecondMonitorLayer2(0),
    physiSecondMonitorLayer3(0), physiSecondMonitorLayer4(0), physiNozzleSupport(0),
    physiBrassTube(0), solidFinalCollimator(0), physiFinalCollimator(0),
    solidMod1(0),         logicMod1(0),          physiMod1(0),
    solidMod2(0),         logicMod2(0),          physiMod2(0),
    solidMod3(0),         logicMod3(0),          physiMod3(0),
    solidMod4(0),         logicMod4(0),          physiMod4(0)
{
    // Messenger to change parameters of the G4TPassiveProtonBeamLineGeometry geometry
    //passiveMessenger = new G4TPassiveProtonBeamLineGeometryMessenger(this);
    
    //***************************** PW ***************************************
    //static G4String ROGeometryName = "DetectorROGeometry";
    //RO = new G4THadrontherapyDetectorROGeometry(ROGeometryName);
    //
    //G4cout << "Going to register Parallel world...";
    //RegisterParallelWorld(RO);
    //G4cout << "... done" << G4endl;

}
/////////////////////////////////////////////////////////////////////////////
G4TPassiveProtonBeamLineGeometry::~G4TPassiveProtonBeamLineGeometry()
{
    //delete passiveMessenger;
    //delete hadrontherapyDetectorConstruction;
}

/////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* G4TPassiveProtonBeamLineGeometry::Construct()
{
    // Sets default geometry and materials
    SetDefaultDimensions();
    
    // Construct the whole Passive Beam Line
    ConstructG4TPassiveProtonBeamLineGeometry();
    G4TPhantomCreationAdding* pca = new G4TPhantomCreationAdding();

    pca->CreateAndAddPhantom(physicalTreatmentRoom,
                             G4ThreeVector(20. *cm, 0. *cm, 0. *cm),
                             G4ThreeVector(40. *cm, 40. *cm, 40. *cm));

    //ConstructPhantom();
    //ConstructVoxelizedRO();

    return physicalTreatmentRoom;
}

// In the following method the DEFAULTS used in the geometry of
// passive beam line are provided
// HERE THE USER CAN CHANGE THE GEOMETRY CHARACTERISTICS OF BEAM
// LINE ELEMENTS, ALTERNATIVELY HE/SHE CAN USE THE MACRO FILE (IF A
// MESSENGER IS PROVIDED)
//
// DEFAULT MATERIAL ARE ALSO PROVIDED
// and COLOURS ARE ALSO DEFINED
// ----------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::SetDefaultDimensions()
{
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
    
    
    // VACUUM PIPE: first track of the beam line is inside vacuum;
    // The PIPE contains the FIRST SCATTERING FOIL and the KAPTON WINDOW
    G4double defaultVacuumZoneXSize = 100.0 *mm;
    vacuumZoneXSize = defaultVacuumZoneXSize;
    
    G4double defaultVacuumZoneYSize = 52.5 *mm;
    vacuumZoneYSize = defaultVacuumZoneYSize;
    
    G4double defaultVacuumZoneZSize = 52.5 *mm;
    vacuumZoneZSize = defaultVacuumZoneZSize;
    
    G4double defaultVacuumZoneXPosition = -3010.0 *mm;
    vacuumZoneXPosition = defaultVacuumZoneXPosition;
    
    // FIRST SCATTERING FOIL: a thin foil performing a first scattering
    // of the original beam
    G4double defaultFirstScatteringFoilXSize = 0.0075 *mm;
    firstScatteringFoilXSize = defaultFirstScatteringFoilXSize;
    
    G4double defaultFirstScatteringFoilYSize = 52.5   *mm;
    firstScatteringFoilYSize = defaultFirstScatteringFoilYSize;
    
    G4double defaultFirstScatteringFoilZSize = 52.5   *mm;
    firstScatteringFoilZSize = defaultFirstScatteringFoilZSize;
    
    G4double defaultFirstScatteringFoilXPosition = 0.0 *mm;
    firstScatteringFoilXPosition = defaultFirstScatteringFoilXPosition;
    
    // KAPTON WINDOW: it prmits the passage of the beam from vacuum to air
    G4double defaultKaptonWindowXSize = 0.010*mm;
    kaptonWindowXSize = defaultKaptonWindowXSize;
    
    G4double defaultKaptonWindowYSize = 5.25*cm;
    kaptonWindowYSize = defaultKaptonWindowYSize;
    
    G4double defaultKaptonWindowZSize = 5.25*cm;
    kaptonWindowZSize = defaultKaptonWindowZSize;
    
    G4double defaultKaptonWindowXPosition = 100.0*mm - defaultKaptonWindowXSize;
    kaptonWindowXPosition = defaultKaptonWindowXPosition;
    
    // STOPPER: is a small cylinder able to stop the central component
    // of the beam (having a gaussian shape). It is connected to the SECON SCATTERING FOIL
    // and represent the second element of the scattering system
    G4double defaultInnerRadiusStopper = 0.*cm;
    innerRadiusStopper = defaultInnerRadiusStopper;
    
    G4double defaultHeightStopper = 3.5*mm;
    heightStopper = defaultHeightStopper;
    
    G4double defaultStartAngleStopper = 0.*deg;
    startAngleStopper = defaultStartAngleStopper;
    
    G4double defaultSpanningAngleStopper = 360.*deg;
    spanningAngleStopper = defaultSpanningAngleStopper;
    
    G4double defaultStopperXPosition = -2705.0 *mm;
    stopperXPosition = defaultStopperXPosition;
    
    G4double defaultStopperYPosition = 0.*m;
    stopperYPosition = defaultStopperYPosition;
    
    G4double defaultStopperZPosition = 0.*m;
    stopperZPosition = defaultStopperZPosition;
    
    G4double defaultOuterRadiusStopper = 2 *mm;
    outerRadiusStopper = defaultOuterRadiusStopper;
    
    // SECOND SCATTERING FOIL: it is another thin foil and provides the
    // final diffusion of the beam. It represents the third element of the scattering
    // system;
    G4double defaultSecondScatteringFoilXSize = 0.0125 *mm;
    secondScatteringFoilXSize = defaultSecondScatteringFoilXSize;
    
    G4double defaultSecondScatteringFoilYSize = 52.5   *mm;
    secondScatteringFoilYSize = defaultSecondScatteringFoilYSize;
    
    G4double defaultSecondScatteringFoilZSize = 52.5   *mm;
    secondScatteringFoilZSize = defaultSecondScatteringFoilZSize;
    
    G4double defaultSecondScatteringFoilXPosition = defaultStopperXPosition + defaultHeightStopper + defaultSecondScatteringFoilXSize;
    secondScatteringFoilXPosition = defaultSecondScatteringFoilXPosition;
    
    G4double defaultSecondScatteringFoilYPosition =  0 *mm;
    secondScatteringFoilYPosition = defaultSecondScatteringFoilYPosition;
    
    G4double defaultSecondScatteringFoilZPosition =  0 *mm;
    secondScatteringFoilZPosition = defaultSecondScatteringFoilZPosition;
    
    // RANGE SHIFTER: is a slab of PMMA acting as energy degreader of
    // primary beam
    
    //Default material of the range shifter
    
    G4double defaultRangeShifterXSize = 5. *mm;
    rangeShifterXSize = defaultRangeShifterXSize;
    
    G4double defaultRangeShifterYSize = 176. *mm;
    rangeShifterYSize = defaultRangeShifterYSize;
    
    G4double defaultRangeShifterZSize = 176. *mm;
    rangeShifterZSize = defaultRangeShifterZSize;
    
    G4double defaultRangeShifterXPosition = -2393.0 *mm;
    rangeShifterXPosition = defaultRangeShifterXPosition;
    
    G4double defaultRangeShifterYPosition = 0. *mm;
    rangeShifterYPosition = defaultRangeShifterYPosition;
    
    G4double defaultRangeShifterZPosition = 0. *mm;
    rangeShifterZPosition = defaultRangeShifterZPosition;
    
    // MOPI DETECTOR: two orthogonal microstrip gas detectors developed
    // by the INFN Section of Turin in collaboration with some
    // of the author of this example. It permits the
    // on-line check of the beam simmetry via the signal
    // integration of the collected charge for each strip.
    
    // Mother volume of MOPI
    
    G4double defaultMOPIMotherVolumeXSize = 12127.0 *um;
    MOPIMotherVolumeXSize = defaultMOPIMotherVolumeXSize;
    
    G4double defaultMOPIMotherVolumeYSize = 40.0 *cm;
    MOPIMotherVolumeYSize = defaultMOPIMotherVolumeYSize;
    
    G4double defaultMOPIMotherVolumeZSize = 40.0 *cm;
    MOPIMotherVolumeZSize = defaultMOPIMotherVolumeZSize;
    
    G4double defaultMOPIMotherVolumeXPosition = -1000.0 *mm;
    MOPIMotherVolumeXPosition = defaultMOPIMotherVolumeXPosition;
    
    G4double defaultMOPIMotherVolumeYPosition = 0.0 *mm;
    MOPIMotherVolumeYPosition = defaultMOPIMotherVolumeYPosition;
    
    G4double defaultMOPIMotherVolumeZPosition = 0.0 *mm;
    MOPIMotherVolumeZPosition = defaultMOPIMotherVolumeZPosition;
    
    // First Kapton Layer of MOPI
    G4double defaultMOPIFirstKaptonLayerXSize = 35 *um;
    MOPIFirstKaptonLayerXSize = defaultMOPIFirstKaptonLayerXSize;
    
    G4double defaultMOPIFirstKaptonLayerYSize = 30 *cm;
    MOPIFirstKaptonLayerYSize = defaultMOPIFirstKaptonLayerYSize;
    
    G4double defaultMOPIFirstKaptonLayerZSize = 30 *cm;
    MOPIFirstKaptonLayerZSize = defaultMOPIFirstKaptonLayerZSize;
    
    G4double defaultMOPIFirstKaptonLayerXPosition = -(MOPIMotherVolumeXSize/2 - (MOPIFirstKaptonLayerXSize/2));
    MOPIFirstKaptonLayerXPosition = defaultMOPIFirstKaptonLayerXPosition;
    
    G4double defaultMOPIFirstKaptonLayerYPosition = 0.0 *mm;
    MOPIFirstKaptonLayerYPosition = defaultMOPIFirstKaptonLayerYPosition;
    
    G4double defaultMOPIFirstKaptonLayerZPosition = 0.0 *mm;
    MOPIFirstKaptonLayerZPosition = defaultMOPIFirstKaptonLayerZPosition;
    
    //First Aluminum  Layer of MOPI
    G4double defaultMOPIFirstAluminumLayerXSize = 15 *um;
    MOPIFirstAluminumLayerXSize = defaultMOPIFirstAluminumLayerXSize;
    
    G4double defaultMOPIFirstAluminumLayerYSize = 30 *cm;
    MOPIFirstAluminumLayerYSize = defaultMOPIFirstAluminumLayerYSize;
    
    G4double defaultMOPIFirstAluminumLayerZSize = 30 *cm;
    MOPIFirstAluminumLayerZSize = defaultMOPIFirstAluminumLayerZSize;
    
    G4double defaultMOPIFirstAluminumLayerXPosition =
            MOPIFirstKaptonLayerXPosition + MOPIFirstKaptonLayerXSize/2 + MOPIFirstAluminumLayerXSize/2;
    MOPIFirstAluminumLayerXPosition = defaultMOPIFirstAluminumLayerXPosition;
    
    G4double defaultMOPIFirstAluminumLayerYPosition = 0.0 *mm;
    MOPIFirstAluminumLayerYPosition = defaultMOPIFirstAluminumLayerYPosition;
    
    G4double defaultMOPIFirstAluminumLayerZPosition = 0.0 *mm;
    MOPIFirstAluminumLayerZPosition = defaultMOPIFirstAluminumLayerZPosition;
    
    // First Air gap of MOPI
    G4double defaultMOPIFirstAirGapXSize = 6000 *um;
    MOPIFirstAirGapXSize = defaultMOPIFirstAirGapXSize;
    
    G4double defaultMOPIFirstAirGapYSize = 30 *cm;
    MOPIFirstAirGapYSize = defaultMOPIFirstAirGapYSize;
    
    G4double defaultMOPIFirstAirGapZSize = 30 *cm;
    MOPIFirstAirGapZSize = defaultMOPIFirstAirGapZSize;
    
    G4double defaultMOPIFirstAirGapXPosition =
            MOPIFirstAluminumLayerXPosition + MOPIFirstAluminumLayerXSize/2 + MOPIFirstAirGapXSize/2;
    MOPIFirstAirGapXPosition = defaultMOPIFirstAirGapXPosition;
    
    G4double defaultMOPIFirstAirGapYPosition = 0.0 *mm;
    MOPIFirstAirGapYPosition = defaultMOPIFirstAirGapYPosition;
    
    G4double defaultMOPIFirstAirGapZPosition = 0.0 *mm;
    MOPIFirstAirGapZPosition = defaultMOPIFirstAirGapZPosition;
    
    // Cathode of MOPI
    G4double defaultMOPICathodeXSize = 25.0 *um;
    MOPICathodeXSize = defaultMOPICathodeXSize;
    
    G4double defaultMOPICathodeYSize = 30.0 *cm;
    MOPICathodeYSize = defaultMOPICathodeYSize;
    
    G4double defaultMOPICathodeZSize = 30.0 *cm;
    MOPICathodeZSize = defaultMOPICathodeZSize;
    
    G4double defaultMOPICathodeXPosition =
            MOPIFirstAirGapXPosition + MOPIFirstAirGapXSize/2 + MOPICathodeXSize/2;
    MOPICathodeXPosition = defaultMOPICathodeXPosition;
    
    G4double defaultMOPICathodeYPosition = 0.0 *mm;
    MOPICathodeYPosition = defaultMOPICathodeYPosition;
    
    G4double defaultMOPICathodeZPosition = 0.0 *mm;
    MOPICathodeZPosition = defaultMOPICathodeZPosition;
    
    // Second Air gap of MOPI
    G4double defaultMOPISecondAirGapXSize = 6000 *um;
    MOPISecondAirGapXSize = defaultMOPISecondAirGapXSize;
    
    G4double defaultMOPISecondAirGapYSize = 30 *cm;
    MOPISecondAirGapYSize = defaultMOPISecondAirGapYSize;
    
    G4double defaultMOPISecondAirGapZSize = 30 *cm;
    MOPISecondAirGapZSize = defaultMOPISecondAirGapZSize;
    
    G4double defaultMOPISecondAirGapXPosition =
            MOPICathodeXPosition + MOPICathodeXSize/2 + MOPISecondAirGapXSize/2;
    MOPISecondAirGapXPosition = defaultMOPISecondAirGapXPosition;
    
    G4double defaultMOPISecondAirGapYPosition = 0.0 *mm;
    MOPISecondAirGapYPosition = defaultMOPISecondAirGapYPosition;
    
    G4double defaultMOPISecondAirGapZPosition = 0.0 *mm;
    MOPISecondAirGapZPosition = defaultMOPISecondAirGapZPosition;
    
    //Second Aluminum  Layer of MOPI
    G4double defaultMOPISecondAluminumLayerXSize = 15 *um;
    MOPISecondAluminumLayerXSize = defaultMOPISecondAluminumLayerXSize;
    
    G4double defaultMOPISecondAluminumLayerYSize = 30 *cm;
    MOPISecondAluminumLayerYSize = defaultMOPISecondAluminumLayerYSize;
    
    G4double defaultMOPISecondAluminumLayerZSize = 30 *cm;
    MOPISecondAluminumLayerZSize = defaultMOPISecondAluminumLayerZSize;
    
    G4double defaultMOPISecondAluminumLayerXPosition =
            MOPISecondAirGapXPosition + MOPISecondAirGapXSize/2 + MOPISecondAluminumLayerXSize/2;
    MOPISecondAluminumLayerXPosition = defaultMOPISecondAluminumLayerXPosition;
    
    G4double defaultMOPISecondAluminumLayerYPosition = 0.0 *mm;
    MOPISecondAluminumLayerYPosition = defaultMOPISecondAluminumLayerYPosition;
    
    G4double defaultMOPISecondAluminumLayerZPosition = 0.0 *mm;
    MOPISecondAluminumLayerZPosition = defaultMOPISecondAluminumLayerZPosition;
    
    // Second Kapton Layer of MOPI
    G4double defaultMOPISecondKaptonLayerXSize = 35 *um;
    MOPISecondKaptonLayerXSize = defaultMOPISecondKaptonLayerXSize;
    
    G4double defaultMOPISecondKaptonLayerYSize = 30 *cm;
    MOPISecondKaptonLayerYSize = defaultMOPISecondKaptonLayerYSize;
    
    G4double defaultMOPISecondKaptonLayerZSize = 30 *cm;
    MOPISecondKaptonLayerZSize = defaultMOPISecondKaptonLayerZSize;
    
    G4double defaultMOPISecondKaptonLayerXPosition =
            MOPISecondAluminumLayerXPosition + MOPISecondAluminumLayerXSize/2 + MOPISecondKaptonLayerXSize/2;
    MOPISecondKaptonLayerXPosition = defaultMOPISecondKaptonLayerXPosition;
    
    G4double defaultMOPISecondKaptonLayerYPosition = 0.0 *mm;
    MOPISecondKaptonLayerYPosition = defaultMOPISecondKaptonLayerYPosition;
    
    G4double defaultMOPISecondKaptonLayerZPosition = 0.0 *mm;
    MOPISecondKaptonLayerZPosition = defaultMOPISecondKaptonLayerZPosition;
    
    
    // FINAL COLLIMATOR: is the collimator giving the final transversal shape
    // of the beam
    G4double defaultinnerRadiusFinalCollimator = 7.5 *mm;
    innerRadiusFinalCollimator = defaultinnerRadiusFinalCollimator;
    
    // DEFAULT DEFINITION OF THE MATERIALS
    // All elements and compound definition follows the NIST database
    
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
    
    //***************************** PW ***************************************
    
    
    
    // MATERIAL ASSIGNMENT
    // Range shifter
    rangeShifterMaterial = airNist;
    
    // Support of the beam line
    beamLineSupportMaterial = aluminumNist;
    
    // Vacuum pipe
    vacuumZoneMaterial = galacticNist;
    
    // Material of the fisrt scattering foil
    firstScatteringFoilMaterial = tantalumNist;
    
    // Material of kapton window
    kaptonWindowMaterial = kaptonNist;
    
    // Material of the stopper
    stopperMaterial = brass;
    
    // Material of the second scattering foil
    secondScatteringFoilMaterial = tantalumNist;
    
    // Materials of the collimators
    firstCollimatorMaterial = PMMANist;
    holeFirstCollimatorMaterial = airNist;
    
    // Box containing the modulator wheel
    modulatorBoxMaterial = aluminumNist;
    holeModulatorBoxMaterial = airNist;
    
    // Materials of the monitor chamber
    layer1MonitorChamberMaterial = kaptonNist;
    layer2MonitorChamberMaterial = copperNistAsMaterial;
    layer3MonitorChamberMaterial = airNist;
    layer4MonitorChamberMaterial = copperNistAsMaterial;
    
    // Mother volume of the MOPI detector
    MOPIMotherVolumeMaterial = airNist;
    MOPIFirstKaptonLayerMaterial = kaptonNist;
    MOPIFirstAluminumLayerMaterial = aluminumNist;
    MOPIFirstAirGapMaterial = airNist;
    MOPICathodeMaterial = mylarNist;
    MOPISecondAirGapMaterial = airNist;
    MOPISecondAluminumLayerMaterial = aluminumNist;
    MOPISecondKaptonLayerMaterial = kaptonNist;
    
    // material of the final nozzle
    nozzleSupportMaterial = PMMANist;
    brassTubeMaterial = brassTube2Material = brassTube3Material = brass;
    holeNozzleSupportMaterial = airNist;
    
    // Material of the final collimator
    finalCollimatorMaterial = brass;




    // Define here the material of the water phantom and of the detector

    G4Material* pMat = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", false);
    phantomMaterial  = pMat;

    phantomSizeX = 40. *cm;
    phantomSizeY = 40. *cm;
    phantomSizeZ = 40. *cm;
    phantomPosition = G4ThreeVector(20. *cm, 0. *cm, 0. *cm);


}

/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::ConstructG4TPassiveProtonBeamLineGeometry()
{
    // -----------------------------
    // Treatment room - World volume
    //------------------------------
    // Treatment room sizes
    const G4double worldX = 400.0 *cm;
    const G4double worldY = 400.0 *cm;
    const G4double worldZ = 400.0 *cm;
    G4bool isotopes = false;
    
    G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
    G4Box* treatmentRoom = new G4Box("TreatmentRoom",worldX,worldY,worldZ);
    G4LogicalVolume* logicTreatmentRoom = new G4LogicalVolume(treatmentRoom,
                                                              airNist,
                                                              "logicTreatmentRoom",
                                                              0,0,0);
    physicalTreatmentRoom = new G4PVPlacement(0,
                                              G4ThreeVector(),
                                              "physicalTreatmentRoom",
                                              logicTreatmentRoom,
                                              0,false,0);
    
    
    // The treatment room is invisible in the Visualisation
    //logicTreatmentRoom -> SetVisAttributes(G4VisAttributes::GetInvisible());
    
    // Components of the Passive Proton Beam Line
    HadrontherapyBeamLineSupport();
    HadrontherapyBeamScatteringFoils();
    HadrontherapyRangeShifter();
    HadrontherapyBeamCollimators();
    HadrontherapyBeamMonitoring();
    HadrontherapyMOPIDetector();
    HadrontherapyBeamNozzle();
    HadrontherapyBeamFinalCollimator();
    
    // The following lines construc a typical modulator wheel inside the Passive Beam line.
    // Please remember to set the nodulator material (default is air, i.e. no modulator!)
    BuildModulator();
}

/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::HadrontherapyBeamLineSupport()
{
    // ------------------//
    // BEAM LINE SUPPORT //
    //-------------------//
    const G4double beamLineSupportXSize = 1.5*m;
    const G4double beamLineSupportYSize = 20.*mm;
    const G4double beamLineSupportZSize = 600.*mm;
    
    const G4double beamLineSupportXPosition = -1745.09 *mm;
    const G4double beamLineSupportYPosition = -230. *mm;
    const G4double beamLineSupportZPosition = 0.*mm;
    
    G4Box* beamLineSupport = new G4Box("BeamLineSupport",
                                       beamLineSupportXSize,
                                       beamLineSupportYSize,
                                       beamLineSupportZSize);
    
    G4LogicalVolume* logicBeamLineSupport = new G4LogicalVolume(beamLineSupport,
                                                                beamLineSupportMaterial,
                                                                "BeamLineSupport");
    physiBeamLineSupport = new G4PVPlacement(0, G4ThreeVector(beamLineSupportXPosition,
                                                              beamLineSupportYPosition,
                                                              beamLineSupportZPosition),
                                             "BeamLineSupport",
                                             logicBeamLineSupport,
                                             physicalTreatmentRoom, false, 0);
    
    // Visualisation attributes of the beam line support
    
    logicBeamLineSupport -> SetVisAttributes(gray);

    //---------------------------------//
    //  Beam line cover 1 (left panel) //
    //---------------------------------//
    const G4double beamLineCoverXSize = 1.5*m;
    const G4double beamLineCoverYSize = 750.*mm;
    const G4double beamLineCoverZSize = 10.*mm;
    
    const G4double beamLineCoverXPosition = -1745.09 *mm;
    const G4double beamLineCoverYPosition = -1000.*mm;
    const G4double beamLineCoverZPosition = 600.*mm;
    
    G4Box* beamLineCover = new G4Box("BeamLineCover",
                                     beamLineCoverXSize,
                                     beamLineCoverYSize,
                                     beamLineCoverZSize);
    
    G4LogicalVolume* logicBeamLineCover = new G4LogicalVolume(beamLineCover,
                                                              beamLineSupportMaterial,
                                                              "BeamLineCover");
    
    physiBeamLineCover = new G4PVPlacement(0, G4ThreeVector(beamLineCoverXPosition,
                                                            beamLineCoverYPosition,
                                                            beamLineCoverZPosition),
                                           "BeamLineCover",
                                           logicBeamLineCover,
                                           physicalTreatmentRoom,
                                           false,
                                           0);
    
    // ---------------------------------//
    //  Beam line cover 2 (rigth panel) //
    // ---------------------------------//
    // It has the same characteristic of beam line cover 1 but set in a different position
    physiBeamLineCover2 = new G4PVPlacement(0, G4ThreeVector(beamLineCoverXPosition,
                                                             beamLineCoverYPosition,
                                                             - beamLineCoverZPosition),
                                            "BeamLineCover2",
                                            logicBeamLineCover,
                                            physicalTreatmentRoom,
                                            false,
                                            0);
    
    logicBeamLineCover -> SetVisAttributes(blue);

}
/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::HadrontherapyBeamScatteringFoils()
{
    // ------------//
    // VACUUM PIPE //
    //-------------//
    //
    // First track of the beam line is inside vacuum;
    // The PIPE contains the FIRST SCATTERING FOIL and the KAPTON WINDOW
    G4Box* vacuumZone = new G4Box("VacuumZone", vacuumZoneXSize, vacuumZoneYSize, vacuumZoneZSize);
    G4LogicalVolume* logicVacuumZone = new G4LogicalVolume(vacuumZone, vacuumZoneMaterial, "VacuumZone");
    G4VPhysicalVolume* physiVacuumZone = new G4PVPlacement(0, G4ThreeVector(vacuumZoneXPosition, 0., 0.),
                                                           "VacuumZone", logicVacuumZone, physicalTreatmentRoom, false, 0);
    // --------------------------//
    // THE FIRST SCATTERING FOIL //
    // --------------------------//
    // A thin foil performing a first scattering
    // of the original beam
    firstScatteringFoil = new G4Box("FirstScatteringFoil",
                                    firstScatteringFoilXSize,
                                    firstScatteringFoilYSize,
                                    firstScatteringFoilZSize);
    
    G4LogicalVolume* logicFirstScatteringFoil = new G4LogicalVolume(firstScatteringFoil,
                                                                    firstScatteringFoilMaterial,
                                                                    "FirstScatteringFoil");
    
    physiFirstScatteringFoil = new G4PVPlacement(0, G4ThreeVector(firstScatteringFoilXPosition, 0.,0.),
                                                 "FirstScatteringFoil", logicFirstScatteringFoil, physiVacuumZone,
                                                 false, 0);
    
    logicFirstScatteringFoil -> SetVisAttributes(skyBlue);
    // -------------------//
    // THE KAPTON WINDOWS //
    //--------------------//
    //It prmits the passage of the beam from vacuum to air

    G4Box* solidKaptonWindow = new G4Box("KaptonWindow",
                                         kaptonWindowXSize,
                                         kaptonWindowYSize,
                                         kaptonWindowZSize);
    
    G4LogicalVolume* logicKaptonWindow = new G4LogicalVolume(solidKaptonWindow,
                                                             kaptonWindowMaterial,
                                                             "KaptonWindow");
    
    physiKaptonWindow = new G4PVPlacement(0, G4ThreeVector(kaptonWindowXPosition, 0., 0.),
                                          "KaptonWindow", logicKaptonWindow,
                                          physiVacuumZone, false,	0);
    
    logicKaptonWindow -> SetVisAttributes(darkOrange3);
    
    // ------------//
    // THE STOPPER //
    //-------------//
    // Is a small cylinder able to stop the central component
    // of the beam (having a gaussian shape). It is connected to the SECON SCATTERING FOIL
    // and represent the second element of the scattering system
    G4double phi = 90. *deg;
    // Matrix definition for a 90 deg rotation with respect to Y axis
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    solidStopper = new G4Tubs("Stopper",
                              innerRadiusStopper,
                              outerRadiusStopper,
                              heightStopper,
                              startAngleStopper,
                              spanningAngleStopper);
    
    logicStopper = new G4LogicalVolume(solidStopper,
                                       stopperMaterial,
                                       "Stopper",
                                       0, 0, 0);
    
    physiStopper = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(stopperXPosition,
                                                                     stopperYPosition,
                                                                     stopperZPosition)),
                                     "Stopper",
                                     logicStopper,
                                     physicalTreatmentRoom,
                                     false,
                                     0);
    
    logicStopper -> SetVisAttributes(red);
    
    // ---------------------------//
    // THE SECOND SCATTERING FOIL //
    // ---------------------------//
    // It is another thin foil and provides the
    // final diffusion of the beam. It represents the third element of the scattering
    // system;
    
    secondScatteringFoil = new G4Box("SecondScatteringFoil",
                                     secondScatteringFoilXSize,
                                     secondScatteringFoilYSize,
                                     secondScatteringFoilZSize);
    
    G4LogicalVolume* logicSecondScatteringFoil = new G4LogicalVolume(secondScatteringFoil,
                                                                     secondScatteringFoilMaterial,
                                                                     "SecondScatteringFoil");
    
    physiSecondScatteringFoil = new G4PVPlacement(0, G4ThreeVector(secondScatteringFoilXPosition,
                                                                   secondScatteringFoilYPosition,
                                                                   secondScatteringFoilZPosition),
                                                  "SeconScatteringFoil",
                                                  logicSecondScatteringFoil,
                                                  physicalTreatmentRoom,
                                                  false,
                                                  0);
    
    logicSecondScatteringFoil -> SetVisAttributes(skyBlue);


}
/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::HadrontherapyRangeShifter()
{
    // ---------------------------- //
    //         THE RANGE SHIFTER    //
    // -----------------------------//
    // It is a slab of PMMA acting as energy degreader of
    // primary beam

    
    solidRangeShifterBox = new G4Box("RangeShifterBox",
                                     rangeShifterXSize,
                                     rangeShifterYSize,
                                     rangeShifterZSize);
    
    logicRangeShifterBox = new G4LogicalVolume(solidRangeShifterBox,
                                               rangeShifterMaterial,
                                               "RangeShifterBox");
    physiRangeShifterBox = new G4PVPlacement(0,
                                             G4ThreeVector(rangeShifterXPosition, 0., 0.),
                                             "RangeShifterBox",
                                             logicRangeShifterBox,
                                             physicalTreatmentRoom,
                                             false,
                                             0);
    
    
    logicRangeShifterBox -> SetVisAttributes(yellow);

    
}
/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::HadrontherapyBeamCollimators()
{
    
    
    // -----------------//
    // FIRST COLLIMATOR //
    // -----------------//
    // It is a slab of PMMA with an hole in its center
    const G4double firstCollimatorXSize = 20.*mm;
    const G4double firstCollimatorYSize = 100.*mm;
    const G4double firstCollimatorZSize = 100.*mm;
    
    const G4double firstCollimatorXPosition = -2673.00*mm;
    const G4double firstCollimatorYPosition = 0.*mm;
    const G4double firstCollimatorZPosition = 0.*mm;
    
    
    G4Box* solidFirstCollimator = new G4Box("FirstCollimator",
                                            firstCollimatorXSize,
                                            firstCollimatorYSize,
                                            firstCollimatorZSize);
    
    G4LogicalVolume* logicFirstCollimator = new G4LogicalVolume(solidFirstCollimator,
                                                                firstCollimatorMaterial,
                                                                "FirstCollimator");
    
    physiFirstCollimator = new G4PVPlacement(0, G4ThreeVector(firstCollimatorXPosition,
                                                              firstCollimatorYPosition,
                                                              firstCollimatorZPosition),
                                             "FirstCollimator",
                                             logicFirstCollimator,
                                             physicalTreatmentRoom,
                                             false,
                                             0);
    // ----------------------------//
    // Hole of the first collimator//
    //-----------------------------//
    G4double innerRadiusHoleFirstCollimator   = 0.*mm;
    G4double outerRadiusHoleFirstCollimator   = 15.*mm;
    G4double hightHoleFirstCollimator         = 20.*mm;
    G4double startAngleHoleFirstCollimator    = 0.*deg;
    G4double spanningAngleHoleFirstCollimator = 360.*deg;
    
    G4Tubs* solidHoleFirstCollimator = new G4Tubs("HoleFirstCollimator",
                                                  innerRadiusHoleFirstCollimator,
                                                  outerRadiusHoleFirstCollimator,
                                                  hightHoleFirstCollimator,
                                                  startAngleHoleFirstCollimator,
                                                  spanningAngleHoleFirstCollimator);
    
    G4LogicalVolume* logicHoleFirstCollimator = new G4LogicalVolume(solidHoleFirstCollimator,
                                                                    holeFirstCollimatorMaterial,
                                                                    "HoleFirstCollimator",
                                                                    0, 0, 0);
    G4double phi = 90. *deg;
    // Matrix definition for a 90 deg rotation. Also used for other volumes
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    physiHoleFirstCollimator = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                                 "HoleFirstCollimator",
                                                 logicHoleFirstCollimator,
                                                 physiFirstCollimator,
                                                 false,
                                                 0);
    // ------------------//
    // SECOND COLLIMATOR //
    //-------------------//
    // It is a slab of PMMA with an hole in its center
    const G4double secondCollimatorXPosition = -1900.00*mm;
    const G4double secondCollimatorYPosition =  0*mm;
    const G4double secondCollimatorZPosition =  0*mm;
    
    physiSecondCollimator = new G4PVPlacement(0, G4ThreeVector(secondCollimatorXPosition,
                                                               secondCollimatorYPosition,
                                                               secondCollimatorZPosition),
                                              "SecondCollimator",
                                              logicFirstCollimator,
                                              physicalTreatmentRoom,
                                              false,
                                              0);
    
    // ------------------------------//
    // Hole of the second collimator //
    // ------------------------------//
    physiHoleSecondCollimator = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                                  "HoleSecondCollimator",
                                                  logicHoleFirstCollimator,
                                                  physiSecondCollimator,
                                                  false,
                                                  0);
    
    // --------------------------------------//
    // FIRST SIDE OF THE MODULATOR BOX      //
    // --------------------------------------//
    // The modulator box is an aluminum box in which
    // the range shifter and the energy modulator are located
    // In this example only the entrance and exit
    // faces of the box are simulated.
    // Each face is an aluminum slab with an hole in its center
    
    const G4double firstCollimatorModulatorXSize = 10.*mm;
    const G4double firstCollimatorModulatorYSize = 200.*mm;
    const G4double firstCollimatorModulatorZSize = 200.*mm;
    
    const G4double firstCollimatorModulatorXPosition = -2523.00*mm;
    const G4double firstCollimatorModulatorYPosition = 0.*mm;
    const G4double firstCollimatorModulatorZPosition = 0.*mm;
    
    G4Box* solidFirstCollimatorModulatorBox = new G4Box("FirstCollimatorModulatorBox",
                                                        firstCollimatorModulatorXSize,
                                                        firstCollimatorModulatorYSize,
                                                        firstCollimatorModulatorZSize);
    
    G4LogicalVolume* logicFirstCollimatorModulatorBox = new G4LogicalVolume(solidFirstCollimatorModulatorBox,
                                                                            modulatorBoxMaterial,
                                                                            "FirstCollimatorModulatorBox");
    
    physiFirstCollimatorModulatorBox = new G4PVPlacement(0, G4ThreeVector(firstCollimatorModulatorXPosition,
                                                                          firstCollimatorModulatorYPosition,
                                                                          firstCollimatorModulatorZPosition),
                                                         "FirstCollimatorModulatorBox",
                                                         logicFirstCollimatorModulatorBox,
                                                         physicalTreatmentRoom, false, 0);
    
    // ----------------------------------------------------//
    //   Hole of the first collimator of the modulator box //
    // ----------------------------------------------------//
    const G4double innerRadiusHoleFirstCollimatorModulatorBox = 0.*mm;
    const G4double outerRadiusHoleFirstCollimatorModulatorBox = 31.*mm;
    const G4double hightHoleFirstCollimatorModulatorBox = 10.*mm;
    const G4double startAngleHoleFirstCollimatorModulatorBox = 0.*deg;
    const G4double spanningAngleHoleFirstCollimatorModulatorBox = 360.*deg;
    
    G4Tubs* solidHoleFirstCollimatorModulatorBox  = new G4Tubs("HoleFirstCollimatorModulatorBox",
                                                               innerRadiusHoleFirstCollimatorModulatorBox,
                                                               outerRadiusHoleFirstCollimatorModulatorBox,
                                                               hightHoleFirstCollimatorModulatorBox ,
                                                               startAngleHoleFirstCollimatorModulatorBox,
                                                               spanningAngleHoleFirstCollimatorModulatorBox);
    
    G4LogicalVolume* logicHoleFirstCollimatorModulatorBox = new G4LogicalVolume(solidHoleFirstCollimatorModulatorBox,
                                                                                holeModulatorBoxMaterial,
                                                                                "HoleFirstCollimatorModulatorBox",
                                                                                0, 0, 0);
    
    physiHoleFirstCollimatorModulatorBox = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                                             "HoleFirstCollimatorModulatorBox",
                                                             logicHoleFirstCollimatorModulatorBox,
                                                             physiFirstCollimatorModulatorBox, false, 0);
    
    // --------------------------------------------------//
    //       SECOND SIDE OF THE MODULATOR BOX            //
    // --------------------------------------------------//
    const G4double secondCollimatorModulatorXSize = 10.*mm;
    const G4double secondCollimatorModulatorYSize = 200.*mm;
    const G4double secondCollimatorModulatorZSize = 200.*mm;
    
    const G4double secondCollimatorModulatorXPosition = -1953.00 *mm;
    
    const G4double secondCollimatorModulatorYPosition = 0.*mm;
    const G4double secondCollimatorModulatorZPosition = 0.*mm;
    
    G4Box* solidSecondCollimatorModulatorBox = new G4Box("SecondCollimatorModulatorBox",
                                                         secondCollimatorModulatorXSize,
                                                         secondCollimatorModulatorYSize,
                                                         secondCollimatorModulatorZSize);
    
    G4LogicalVolume* logicSecondCollimatorModulatorBox = new G4LogicalVolume(solidSecondCollimatorModulatorBox,
                                                                             modulatorBoxMaterial,
                                                                             "SecondCollimatorModulatorBox");
    
    physiSecondCollimatorModulatorBox = new G4PVPlacement(0, G4ThreeVector(secondCollimatorModulatorXPosition,
                                                                           secondCollimatorModulatorYPosition,
                                                                           secondCollimatorModulatorZPosition),
                                                          "SecondCollimatorModulatorBox",
                                                          logicSecondCollimatorModulatorBox,
                                                          physicalTreatmentRoom, false, 0);
    
    // ----------------------------------------------//
    //   Hole of the second collimator modulator box //
    // ----------------------------------------------//
    const G4double innerRadiusHoleSecondCollimatorModulatorBox = 0.*mm;
    const G4double outerRadiusHoleSecondCollimatorModulatorBox = 31.*mm;
    const G4double hightHoleSecondCollimatorModulatorBox = 10.*mm;
    const G4double startAngleHoleSecondCollimatorModulatorBox = 0.*deg;
    const G4double spanningAngleHoleSecondCollimatorModulatorBox = 360.*deg;
    
    G4Tubs* solidHoleSecondCollimatorModulatorBox  = new G4Tubs("HoleSecondCollimatorModulatorBox",
                                                                innerRadiusHoleSecondCollimatorModulatorBox,
                                                                outerRadiusHoleSecondCollimatorModulatorBox,
                                                                hightHoleSecondCollimatorModulatorBox ,
                                                                startAngleHoleSecondCollimatorModulatorBox,
                                                                spanningAngleHoleSecondCollimatorModulatorBox);
    
    G4LogicalVolume* logicHoleSecondCollimatorModulatorBox = new G4LogicalVolume(solidHoleSecondCollimatorModulatorBox,
                                                                                 holeModulatorBoxMaterial,
                                                                                 "HoleSecondCollimatorModulatorBox",
                                                                                 0, 0, 0);
    
    physiHoleSecondCollimatorModulatorBox = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                                              "HoleSecondCollimatorModulatorBox",
                                                              logicHoleSecondCollimatorModulatorBox,
                                                              physiSecondCollimatorModulatorBox, false, 0);
    
    logicFirstCollimator -> SetVisAttributes(yellow);
    logicFirstCollimatorModulatorBox -> SetVisAttributes(blue);
    logicSecondCollimatorModulatorBox -> SetVisAttributes(blue);



}
/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::HadrontherapyBeamMonitoring()
{
    // ----------------------------
    //   THE FIRST MONITOR CHAMBER
    // ----------------------------
    // A monitor chamber is a free-air  ionisation chamber
    // able to measure do proton fluence during the treatment.
    // Here its responce is not simulated in terms of produced
    // charge but only the energy losses are taked into account.
    // Each chamber consist of 9 mm of air in a box
    // that has two layers one of kapton and one
    // of copper
    const G4double monitor1XSize = 4.525022*mm;
    const G4double monitor2XSize = 0.000011*mm;
    const G4double monitor3XSize = 4.5*mm;
    const G4double monitorYSize = 10.*cm;
    const G4double monitorZSize = 10.*cm;
    const G4double monitor1XPosition = -1262.47498 *mm;
    const G4double monitor2XPosition = -4.500011*mm;
    const G4double monitor4XPosition = 4.500011*mm;

    

    G4Box* solidFirstMonitorLayer1 = new G4Box("FirstMonitorLayer1",
                                               monitor1XSize,
                                               monitorYSize,
                                               monitorZSize);
    
    G4LogicalVolume* logicFirstMonitorLayer1 = new G4LogicalVolume(solidFirstMonitorLayer1,
                                                                   layer1MonitorChamberMaterial,
                                                                   "FirstMonitorLayer1");
    
    physiFirstMonitorLayer1 = new G4PVPlacement(0,
                                                G4ThreeVector(monitor1XPosition,0.*cm,0.*cm),
                                                "FirstMonitorLayer1",
                                                logicFirstMonitorLayer1,
                                                physicalTreatmentRoom,
                                                false,
                                                0);
    
    G4Box* solidFirstMonitorLayer2 = new G4Box("FirstMonitorLayer2",
                                               monitor2XSize,
                                               monitorYSize,
                                               monitorZSize);
    
    G4LogicalVolume* logicFirstMonitorLayer2 = new G4LogicalVolume(solidFirstMonitorLayer2,
                                                                   layer2MonitorChamberMaterial,
                                                                   "FirstMonitorLayer2");
    
    physiFirstMonitorLayer2 = new G4PVPlacement(0, G4ThreeVector(monitor2XPosition,0.*cm,0.*cm),
                                                "FirstMonitorLayer2",
                                                logicFirstMonitorLayer2,
                                                physiFirstMonitorLayer1,
                                                false,
                                                0);
    
    G4Box* solidFirstMonitorLayer3 = new G4Box("FirstMonitorLayer3",
                                               monitor3XSize,
                                               monitorYSize,
                                               monitorZSize);
    
    G4LogicalVolume* logicFirstMonitorLayer3 = new G4LogicalVolume(solidFirstMonitorLayer3,
                                                                   layer3MonitorChamberMaterial,
                                                                   "FirstMonitorLayer3");
    
    physiFirstMonitorLayer3 = new G4PVPlacement(0,
                                                G4ThreeVector(0.*mm,0.*cm,0.*cm),
                                                "MonitorLayer3",
                                                logicFirstMonitorLayer3,
                                                physiFirstMonitorLayer1,
                                                false,
                                                0);
    
    G4Box* solidFirstMonitorLayer4 = new G4Box("FirstMonitorLayer4",
                                               monitor2XSize,
                                               monitorYSize,
                                               monitorZSize);
    
    G4LogicalVolume* logicFirstMonitorLayer4 = new G4LogicalVolume(solidFirstMonitorLayer4,
                                                                   layer4MonitorChamberMaterial,
                                                                   "FirstMonitorLayer4");
    
    physiFirstMonitorLayer4 = new G4PVPlacement(0, G4ThreeVector(monitor4XPosition,0.*cm,0.*cm),
                                                "FirstMonitorLayer4",
                                                logicFirstMonitorLayer4,
                                                physiFirstMonitorLayer1, false, 0);
    // ----------------------------//
    // THE SECOND MONITOR CHAMBER  //
    // ----------------------------//
    physiSecondMonitorLayer1 = new G4PVPlacement(0, G4ThreeVector(-1131.42493 *mm,0.*cm,0.*cm),
                                                 "SecondMonitorLayer1", logicFirstMonitorLayer1,physicalTreatmentRoom, false, 0);
    
    physiSecondMonitorLayer2 = new G4PVPlacement(0, G4ThreeVector( monitor2XPosition,0.*cm,0.*cm), "SecondMonitorLayer2",
                                                 logicFirstMonitorLayer2, physiSecondMonitorLayer1, false, 0);
    
    physiSecondMonitorLayer3 = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*cm,0.*cm), "MonitorLayer3",
                                                 logicFirstMonitorLayer3, physiSecondMonitorLayer1, false, 0);
    
    physiSecondMonitorLayer4 = new G4PVPlacement(0, G4ThreeVector(monitor4XPosition,0.*cm,0.*cm), "SecondMonitorLayer4",
                                                 logicFirstMonitorLayer4, physiSecondMonitorLayer1, false, 0);
    
    logicFirstMonitorLayer3 -> SetVisAttributes(white);
    

    
}
/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::HadrontherapyMOPIDetector()
{

    // --------------------------------//
    //        THE MOPI DETECTOR        //
    // --------------------------------//
    // MOPI DETECTOR: two orthogonal microstrip gas detectors developed
    // by the INFN Section of Turin in collaboration with some
    // of the author of this example. It permits the
    // on-line check of the beam simmetry via the signal
    // integration of the collected charge for each strip.
    //
    // In this example it is simulated as:
    // 1. First anode: 35 mu of kapton + 15 mu of aluminum,
    // 2. First air gap: 6 mm of air,
    // 3. The cathode: 1 mu Al + 25 mu mylar + 1 mu Al
    //    (in common with the two air gap),
    // 4. Second air gap: 6 mm of air,
    // 5  Second anode: 15 mu Al + 35 mu kapton
    // Color used in the graphical output
    
    
    // Mother volume
    solidMOPIMotherVolume = new G4Box("MOPIMotherVolume",
                                      MOPIMotherVolumeXSize/2,
                                      MOPIMotherVolumeYSize/2,
                                      MOPIMotherVolumeYSize/2);
    
    logicMOPIMotherVolume = new G4LogicalVolume(solidMOPIMotherVolume,
                                                MOPIMotherVolumeMaterial,
                                                "MOPIMotherVolume");
    physiMOPIMotherVolume = new G4PVPlacement(0,
                                              G4ThreeVector(MOPIMotherVolumeXPosition,
                                                            MOPIMotherVolumeYPosition,
                                                            MOPIMotherVolumeZPosition),
                                              "MOPIMotherVolume",
                                              logicMOPIMotherVolume,
                                              physicalTreatmentRoom,
                                              false,
                                              0);
    
    // First Kapton layer
    solidMOPIFirstKaptonLayer = new G4Box("MOPIFirstKaptonLayer",
                                          MOPIFirstKaptonLayerXSize/2,
                                          MOPIFirstKaptonLayerYSize/2 ,
                                          MOPIFirstKaptonLayerZSize/2);
    
    logicMOPIFirstKaptonLayer = new G4LogicalVolume(solidMOPIFirstKaptonLayer,
                                                    MOPIFirstKaptonLayerMaterial,
                                                    "MOPIFirstKaptonLayer");
    
    physiMOPIFirstKaptonLayer = new G4PVPlacement(0,
                                                  G4ThreeVector(MOPIFirstKaptonLayerXPosition,
                                                                MOPIFirstKaptonLayerYPosition ,
                                                                MOPIFirstKaptonLayerZPosition),
                                                  "MOPIFirstKaptonLayer",
                                                  logicMOPIFirstKaptonLayer,
                                                  physiMOPIMotherVolume,
                                                  false,
                                                  0);
    
    // First Aluminum layer
    solidMOPIFirstAluminumLayer = new G4Box("MOPIFirstAluminumLayer",
                                            MOPIFirstAluminumLayerXSize/2,
                                            MOPIFirstAluminumLayerYSize/2 ,
                                            MOPIFirstAluminumLayerZSize/2);
    
    logicMOPIFirstAluminumLayer = new G4LogicalVolume(solidMOPIFirstAluminumLayer,
                                                      MOPIFirstAluminumLayerMaterial,
                                                      "MOPIFirstAluminumLayer");
    
    physiMOPIFirstAluminumLayer = new G4PVPlacement(0,
                                                    G4ThreeVector(MOPIFirstAluminumLayerXPosition,
                                                                  MOPIFirstAluminumLayerYPosition ,
                                                                  MOPIFirstAluminumLayerZPosition),
                                                    "MOPIFirstAluminumLayer",
                                                    logicMOPIFirstAluminumLayer, physiMOPIMotherVolume, false, 0);
    
    // First Air GAP
    solidMOPIFirstAirGap = new G4Box("MOPIFirstAirGap",
                                     MOPIFirstAirGapXSize/2,
                                     MOPIFirstAirGapYSize/2,
                                     MOPIFirstAirGapZSize/2);
    
    logicMOPIFirstAirGap = new G4LogicalVolume(solidMOPIFirstAirGap,
                                               MOPIFirstAirGapMaterial,
                                               "MOPIFirstAirgap");
    
    physiMOPIFirstAirGap = new G4PVPlacement(0,
                                             G4ThreeVector(MOPIFirstAirGapXPosition,
                                                           MOPIFirstAirGapYPosition ,
                                                           MOPIFirstAirGapZPosition),
                                             "MOPIFirstAirGap",
                                             logicMOPIFirstAirGap, physiMOPIMotherVolume, false, 0);
    
    
    // The Cathode
    solidMOPICathode = new G4Box("MOPICathode",
                                 MOPICathodeXSize/2,
                                 MOPICathodeYSize/2,
                                 MOPICathodeZSize/2);
    
    logicMOPICathode = new G4LogicalVolume(solidMOPICathode,
                                           MOPICathodeMaterial,
                                           "MOPICathode");
    
    physiMOPICathode = new G4PVPlacement(0,
                                         G4ThreeVector(MOPICathodeXPosition,
                                                       MOPICathodeYPosition ,
                                                       MOPICathodeZPosition),
                                         "MOPICathode",
                                         logicMOPICathode,
                                         physiMOPIMotherVolume, false, 0);
    
    // Second Air GAP
    solidMOPISecondAirGap = new G4Box("MOPISecondAirGap",
                                      MOPISecondAirGapXSize/2,
                                      MOPISecondAirGapYSize/2,
                                      MOPISecondAirGapZSize/2);
    
    logicMOPISecondAirGap = new G4LogicalVolume(solidMOPISecondAirGap,
                                                MOPISecondAirGapMaterial,
                                                "MOPISecondAirgap");
    
    physiMOPISecondAirGap = new G4PVPlacement(0,
                                              G4ThreeVector(MOPISecondAirGapXPosition,
                                                            MOPISecondAirGapYPosition ,
                                                            MOPISecondAirGapZPosition),
                                              "MOPISecondAirGap",
                                              logicMOPISecondAirGap, physiMOPIMotherVolume, false, 0);
    
    // Second Aluminum layer
    solidMOPISecondAluminumLayer = new G4Box("MOPISecondAluminumLayer",
                                             MOPISecondAluminumLayerXSize/2,
                                             MOPISecondAluminumLayerYSize/2 ,
                                             MOPISecondAluminumLayerZSize/2);
    
    logicMOPISecondAluminumLayer = new G4LogicalVolume(solidMOPISecondAluminumLayer,
                                                       MOPISecondAluminumLayerMaterial,
                                                       "MOPISecondAluminumLayer");
    
    physiMOPISecondAluminumLayer = new G4PVPlacement(0,
                                                     G4ThreeVector(MOPISecondAluminumLayerXPosition,
                                                                   MOPISecondAluminumLayerYPosition ,
                                                                   MOPISecondAluminumLayerZPosition),
                                                     "MOPISecondAluminumLayer",
                                                     logicMOPISecondAluminumLayer,
                                                     physiMOPIMotherVolume,
                                                     false,
                                                     0);
    
    // Second Kapton layer
    solidMOPISecondKaptonLayer = new G4Box("MOPISecondKaptonLayer",
                                           MOPISecondKaptonLayerXSize/2,
                                           MOPISecondKaptonLayerYSize/2 ,
                                           MOPISecondKaptonLayerZSize/2);
    
    logicMOPISecondKaptonLayer = new G4LogicalVolume(solidMOPISecondKaptonLayer,
                                                     MOPIFirstKaptonLayerMaterial,
                                                     "MOPISecondKaptonLayer");
    
    physiMOPISecondKaptonLayer = new G4PVPlacement(0,
                                                   G4ThreeVector(MOPISecondKaptonLayerXPosition,
                                                                 MOPISecondKaptonLayerYPosition ,
                                                                 MOPISecondKaptonLayerZPosition),
                                                   "MOPISecondKaptonLayer",
                                                   logicMOPISecondKaptonLayer,
                                                   physiMOPIMotherVolume,
                                                   false,
                                                   0);
    
    logicMOPIFirstAirGap -> SetVisAttributes(darkGreen);
    logicMOPISecondAirGap -> SetVisAttributes(darkGreen);
    
    
}
/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::HadrontherapyBeamNozzle()
{
    // ------------------------------//
    // THE FINAL TUBE AND COLLIMATOR //
    //-------------------------------//
    // The last part of the transport beam line consists of
    // a 59 mm thick PMMA slab (to stop all the diffused radiation), a 370 mm brass tube
    // (to well collimate the proton beam) and a final collimator with 25 mm diameter
    // aperture (that provide the final trasversal shape of the beam)
    
    // -------------------//
    //     PMMA SUPPORT   //
    // -------------------//
    const G4double nozzleSupportXSize = 29.5 *mm;
    const G4double nozzleSupportYSize = 180. *mm;
    const G4double nozzleSupportZSize = 180. *mm;
    
    const G4double nozzleSupportXPosition = -397.50 *mm;
    
    G4double phi = 90. *deg;
    // Matrix definition for a 90 deg rotation. Also used for other volumes
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    G4Box* solidNozzleSupport = new G4Box("NozzleSupport",
                                          nozzleSupportXSize,
                                          nozzleSupportYSize,
                                          nozzleSupportZSize);
    
    G4LogicalVolume* logicNozzleSupport = new G4LogicalVolume(solidNozzleSupport,
                                                              nozzleSupportMaterial,
                                                              "NozzleSupport");
    
    physiNozzleSupport = new G4PVPlacement(0, G4ThreeVector(nozzleSupportXPosition,0., 0.),
                                           "NozzleSupport",
                                           logicNozzleSupport,
                                           physicalTreatmentRoom,
                                           false,
                                           0);
    
    logicNozzleSupport -> SetVisAttributes(yellow);
    
    
    
    //------------------------------------//
    // HOLE IN THE SUPPORT                //
    //------------------------------------//
    const G4double innerRadiusHoleNozzleSupport = 0.*mm;
    const G4double outerRadiusHoleNozzleSupport = 21.5*mm;
    const G4double hightHoleNozzleSupport = 29.5 *mm;
    const G4double startAngleHoleNozzleSupport = 0.*deg;
    const G4double spanningAngleHoleNozzleSupport = 360.*deg;
    
    G4Tubs* solidHoleNozzleSupport = new G4Tubs("HoleNozzleSupport",
                                                innerRadiusHoleNozzleSupport,
                                                outerRadiusHoleNozzleSupport,
                                                hightHoleNozzleSupport,
                                                startAngleHoleNozzleSupport,
                                                spanningAngleHoleNozzleSupport);
    
    G4LogicalVolume* logicHoleNozzleSupport = new G4LogicalVolume(solidHoleNozzleSupport,
                                                                  holeNozzleSupportMaterial,
                                                                  "HoleNozzleSupport",
                                                                  0,
                                                                  0,
                                                                  0);
    
    
    physiHoleNozzleSupport = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                               "HoleNozzleSupport",
                                               logicHoleNozzleSupport,
                                               physiNozzleSupport,
                                               false, 0);
    
    logicHoleNozzleSupport -> SetVisAttributes(darkOrange3);
    
    // ---------------------------------//
    //     BRASS TUBE 1 (phantom side)    //
    // ---------------------------------//
    const G4double innerRadiusBrassTube= 18.*mm;
    const G4double outerRadiusBrassTube = 21.5 *mm;
    const G4double hightBrassTube = 140.5*mm;
    const G4double startAngleBrassTube = 0.*deg;
    const G4double spanningAngleBrassTube = 360.*deg;
    
    const G4double brassTubeXPosition = -227.5 *mm;
    
    G4Tubs* solidBrassTube = new G4Tubs("BrassTube",
                                        innerRadiusBrassTube,
                                        outerRadiusBrassTube,
                                        hightBrassTube,
                                        startAngleBrassTube,
                                        spanningAngleBrassTube);
    
    G4LogicalVolume* logicBrassTube = new G4LogicalVolume(solidBrassTube,
                                                          brassTubeMaterial,
                                                          "BrassTube",
                                                          0, 0, 0);
    
    physiBrassTube = new G4PVPlacement(G4Transform3D(rm,
                                                     G4ThreeVector(brassTubeXPosition,
                                                                   0.,
                                                                   0.)),
                                       "BrassTube",
                                       logicBrassTube,
                                       physicalTreatmentRoom,
                                       false,
                                       0);
    
    logicBrassTube -> SetVisAttributes(darkOrange3);
    
    // ----------------------------------------------//
    //     BRASS TUBE 2 (inside the PMMA support)    //
    // ----------------------------------------------//
    const G4double innerRadiusBrassTube2= 18.*mm;
    const G4double outerRadiusBrassTube2 = 21.5 *mm;
    const G4double hightBrassTube2 = 29.5*mm;
    const G4double startAngleBrassTube2 = 0.*deg;
    const G4double spanningAngleBrassTube2 = 360.*deg;
    
    
    G4Tubs* solidBrassTube2 = new G4Tubs("BrassTube2",
                                         innerRadiusBrassTube2,
                                         outerRadiusBrassTube2,
                                         hightBrassTube2,
                                         startAngleBrassTube2,
                                         spanningAngleBrassTube2);
    
    G4LogicalVolume* logicBrassTube2 = new G4LogicalVolume(solidBrassTube2,
                                                           brassTube2Material,
                                                           "BrassTube2",
                                                           0, 0, 0);
    
    physiBrassTube2 = new G4PVPlacement(0,
                                        G4ThreeVector(0,0.,0.),
                                        logicBrassTube2,
                                        "BrassTube2",
                                        logicHoleNozzleSupport,
                                        false,
                                        0);
    
    logicBrassTube2 -> SetVisAttributes(darkOrange3);
    
    
    // --------------------------------------//
    //     BRASS TUBE 3 (beam line side)    //
    // -------------------------------------//
    const G4double innerRadiusBrassTube3= 18.*mm;
    const G4double outerRadiusBrassTube3 = 21.5 *mm;
    const G4double hightBrassTube3 = 10.0 *mm;
    const G4double startAngleBrassTube3 = 0.*deg;
    const G4double spanningAngleBrassTube3 = 360.*deg;
    
    const G4double brassTube3XPosition = -437 *mm;
    
    G4Tubs* solidBrassTube3 = new G4Tubs("BrassTube3",
                                         innerRadiusBrassTube3,
                                         outerRadiusBrassTube3,
                                         hightBrassTube3,
                                         startAngleBrassTube3,
                                         spanningAngleBrassTube3);
    
    G4LogicalVolume* logicBrassTube3 = new G4LogicalVolume(solidBrassTube3,
                                                           brassTube3Material,
                                                           "BrassTube3",
                                                           0, 0, 0);
    
    physiBrassTube3 = new G4PVPlacement(G4Transform3D(rm,
                                                      G4ThreeVector(brassTube3XPosition,
                                                                    0.,
                                                                    0.)),
                                        "BrassTube3",
                                        logicBrassTube3,
                                        physicalTreatmentRoom,
                                        false,
                                        0);
    
    logicBrassTube3 -> SetVisAttributes(darkOrange3);
}
/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::HadrontherapyBeamFinalCollimator()
{
    // -----------------------//
    //     FINAL COLLIMATOR   //
    //------------------------//
    const G4double outerRadiusFinalCollimator = 21.5*mm;
    const G4double hightFinalCollimator = 3.5*mm;
    const G4double startAngleFinalCollimator = 0.*deg;
    const G4double spanningAngleFinalCollimator = 360.*deg;
    const G4double finalCollimatorXPosition = -83.5 *mm;
    
    G4double phi = 90. *deg;
    
    // Matrix definition for a 90 deg rotation. Also used for other volumes
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    solidFinalCollimator = new G4Tubs("FinalCollimator",
                                      innerRadiusFinalCollimator,
                                      outerRadiusFinalCollimator,
                                      hightFinalCollimator,
                                      startAngleFinalCollimator,
                                      spanningAngleFinalCollimator);
    
    G4LogicalVolume* logicFinalCollimator = new G4LogicalVolume(solidFinalCollimator,
                                                                finalCollimatorMaterial,
                                                                "FinalCollimator",
                                                                0,
                                                                0,
                                                                0);
    
    physiFinalCollimator = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(finalCollimatorXPosition,0.,0.)),
                                             "FinalCollimator", logicFinalCollimator, physicalTreatmentRoom, false, 0);
    
    logicFinalCollimator -> SetVisAttributes(yellow);
}

void G4TPassiveProtonBeamLineGeometry::ConstructPhantom()
{

    // ***********************************  Phantom Geometry

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
                                              physicalTreatmentRoom,
                                              false,
                                              0);

    // Visualisation attributes of the phantom
    red = new G4VisAttributes(G4Colour(255/255., 0/255. ,0/255.));
    red -> SetVisibility(true);
    red -> SetForceSolid(true);
    red -> SetForceWireframe(true);
    phantomLogicalVolume -> SetVisAttributes(red);





    // ***********************************  Detector Geometry


    /*
    // Definition of the solid volume of the Detector
    detector = new G4Box("Detector",

                         ROSizeX/2,

                         ROSizeY/2,

                         ROSizeZ/2);

    // Definition of the logic volume of the Phantom
    detectorLogicalVolume = new G4LogicalVolume(detector,
                                                ROMaterial,
                                                "DetectorLog",
                                                0,0,0);
    // Definition of the physical volume of the Phantom
    detectorPhysicalVolume = new G4PVPlacement(0,
                                               G4ThreeVector(),//ROPosition, // Setted by displacement
                                               "DetectorPhys",
                                               detectorLogicalVolume,
                                               phantomPhysicalVolume,
                                               false,0);

    // Visualisation attributes of the detector
    skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
    skyBlue -> SetVisibility(true);
    skyBlue -> SetForceSolid(true);
    //skyBlue -> SetForceWireframe(true);
    detectorLogicalVolume -> SetVisAttributes(skyBlue);

    */

}


void G4TPassiveProtonBeamLineGeometry::BuildModulator(){


    // ********************** initialization

    pi=4*std::atan(1.);
    StepNumbers=22;
    Weight=new G4double[StepNumbers];
    StepThickness=new G4double[StepNumbers];
    StartingAngle=new G4double[StepNumbers];
    SpanningAngle=new G4double[StepNumbers];
    PositionMod=new G4ThreeVector[StepNumbers];


    solidMod=new G4Tubs *[StepNumbers];
    logicMod=new G4LogicalVolume *[StepNumbers];
    physiMod=new G4VPhysicalVolume *[(4*(StepNumbers-1)+1)];

    for (G4int i=0;i<StepNumbers;i++)
    {
        Weight[i]=0;
        StepThickness[i]=0;
        StartingAngle[i]=0;
        SpanningAngle[i]=0;
        PositionMod[i]=G4ThreeVector(0,0,0);
        solidMod[i]=0;
        logicMod[i]=0;

    }
    for (G4int i=0;i<4*(StepNumbers-1)+1;i++)
    {
        physiMod[i]=0;
    }



    // ModulatorDefaultProperties()

    /*
       Here we initialize the step properties of Modulator wheel, you can create your
       specific modulator by changing the values in this class or writing them in an external
       file and activate reading from file via a macrofile.
    */

    StepThickness[0]=0; Weight[0]=.14445;
    StepThickness[1]=.8; Weight[1]=.05665;
    StepThickness[2]=1.6; Weight[2]=.05049;
    StepThickness[3]=2.4; Weight[3]=.04239;
    StepThickness[4]=3.2; Weight[4]=.04313;
    StepThickness[5]=4.0; Weight[5]=.03879;
    StepThickness[6]=4.8; Weight[6]=.04182;
    StepThickness[7]=5.6; Weight[7]=.03422;
    StepThickness[8]=6.4; Weight[8]=.03469;
    StepThickness[9]=7.2; Weight[9]=.03589;
    StepThickness[10]=8.0; Weight[10]=.03633;
    StepThickness[11]=8.8; Weight[11]=.03842;
    StepThickness[12]=9.6; Weight[12]=.03688;
    StepThickness[13]=10.4; Weight[13]=.03705;
    StepThickness[14]=11.2; Weight[14]=.03773;
    StepThickness[15]=12.0; Weight[15]=.03968;
    StepThickness[16]=12.8; Weight[16]=.04058;
    StepThickness[17]=13.6; Weight[17]=.03903;
    StepThickness[18]=14.4; Weight[18]=.04370;
    StepThickness[19]=15.2; Weight[19]=.03981;
    StepThickness[20]=16.0; Weight[20]=.05226;
    StepThickness[21]=16.8; Weight[21]=.03603;

    // GetStepInformation();

    G4double TotalWeight=0;
    // convert the absolute weight values to relative ones
    for(G4int i=0;i<StepNumbers;i++)
    {
        TotalWeight+=Weight[i];
    }

    for(G4int i=0;i<StepNumbers;i++)
    {
        Weight[i]=Weight[i]/TotalWeight;
    }

    // To build the RMW step layers will be put one after another

    StartingAngle[0]=0 *deg;
    SpanningAngle[0]=90 *deg;
    G4double PositionModx;
    G4double WholeStartingAngle=0 *deg;
    G4double WholeThickness=0;
    for(G4int i=1;i<StepNumbers;i++)
    {
        StartingAngle[i]=WholeStartingAngle+(Weight[i-1]*(2*pi))/8;
        SpanningAngle[i]=90* deg -2*StartingAngle[i];
        StepThickness[i]=StepThickness[i]-WholeThickness;
        PositionModx=WholeThickness+StepThickness[i]/2.;
        PositionMod[i]=G4ThreeVector(0,0,PositionModx);
        WholeThickness+=StepThickness[i];
        WholeStartingAngle=StartingAngle[i];
    }

    rm = new G4RotationMatrix();
    G4double phi = 270. *deg;
    rm -> rotateY(phi);


    // BuildModulato

    G4bool isotopes = false;
    G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);


    Mod0Mater = airNist;
    ModMater = airNist; // You have to change modulator material via a macrofile (default is air)

    innerRadiusOfTheTube = 2.5 *cm;
    outerRadiusOfTheTube = 9.5 *cm;

    // Mother of the modulator wheel
    G4ThreeVector positionMotherMod = G4ThreeVector(-2160.50 *mm, 30 *mm, 50 *mm);

    G4Box* solidMotherMod = new G4Box("MotherMod", 12 *cm, 12 *cm, 12 *cm);

    logicMotherMod = new G4LogicalVolume(solidMotherMod, Mod0Mater,"MotherMod",0,0,0);

    physiMotherMod = new G4PVPlacement(rm,positionMotherMod,  "MotherMod",
                                       logicMotherMod,
                                       physicalTreatmentRoom,
                                       false,
                                       0);


    // BuildSteps()

    //----------------------------------------------------------
    // Mother volume of first quarter of the modulator
    //----------------------------------------------------------

    G4double hightOfTheTube0 = 10.0 *cm;
    G4double startAngleOfTheTube0 = 0 *deg;
    G4double spanningAngleOfTheTube0 = 90 *deg;

    G4RotationMatrix rm1;
    rm1.rotateZ(0 *deg);

    G4ThreeVector positionMod1 = G4ThreeVector(0*cm,0*cm,0*cm);

    solidMod1 = new G4Tubs("Mod1",
                           innerRadiusOfTheTube,
                           outerRadiusOfTheTube,
                           hightOfTheTube0/2.,
                           startAngleOfTheTube0,
                           spanningAngleOfTheTube0);

    logicMod1 = new G4LogicalVolume(solidMod1, Mod0Mater, "Mod1",0,0,0);

    physiMod1 = new G4PVPlacement(G4Transform3D(rm1, positionMod1),
                                  logicMod1,
                                  "Mod1",
                                  logicMotherMod,
                                  false,
                                  0);

    //----------------------------------------------------------
    //  modulator steps
    //----------------------------------------------------------
    for (G4int i=1;i<StepNumbers;i++)
    {

        solidMod[i] = new G4Tubs("Modstep",
                                 innerRadiusOfTheTube,
                                 outerRadiusOfTheTube,
                                 StepThickness[i]/2.,
                                 StartingAngle[i],
                                 SpanningAngle[i]);

        logicMod[i] = new G4LogicalVolume(solidMod[i],
                                          ModMater, "Modstep",0,0,0);

        physiMod[i] = new G4PVPlacement(0,
                                        PositionMod[i],
                                        logicMod[i],
                                        "Modstep",
                                        logicMod1,
                                        false,
                                        0);


    }

    //----------------------------------------------------------
    // Mother volume of the second modulator quarter
    //----------------------------------------------------------

    G4RotationMatrix rm2;
    rm2.rotateZ(90 *deg);

    G4ThreeVector positionMod2 = G4ThreeVector(0*cm,0*cm,0*cm);

    solidMod2 = new G4Tubs("Mod2",
                           innerRadiusOfTheTube,
                           outerRadiusOfTheTube,
                           hightOfTheTube0/2.,
                           startAngleOfTheTube0,
                           spanningAngleOfTheTube0);

    logicMod2 = new G4LogicalVolume(solidMod2,
                                    Mod0Mater, "Mod2",0,0,0);


    physiMod2 = new G4PVPlacement(G4Transform3D(rm2, positionMod2),
                                  logicMod2,
                                  "Mod2",
                                  logicMotherMod,
                                  false,
                                  0);


    for (G4int i=1;i<StepNumbers;i++)
    {

        physiMod[StepNumbers+i-1] = new G4PVPlacement(0,
                                                      PositionMod[i],
                                                      logicMod[i],
                                                      "Modstep",
                                                      logicMod2,
                                                      false,
                                                      0);

    }

    //----------------------------------------------------------
    // Mother volume of the third modulator quarter
    //----------------------------------------------------------

    G4RotationMatrix rm3;
    rm3.rotateZ(180 *deg);

    G4ThreeVector positionMod3 = G4ThreeVector(0*cm,0*cm,0*cm);

    solidMod3 = new G4Tubs("Mod3",
                           innerRadiusOfTheTube,
                           outerRadiusOfTheTube,
                           hightOfTheTube0,
                           startAngleOfTheTube0/2.,
                           spanningAngleOfTheTube0);

    logicMod3 = new G4LogicalVolume(solidMod3,
                                    Mod0Mater, "Mod3",0,0,0);


    physiMod3 = new G4PVPlacement(G4Transform3D(rm3, positionMod3),
                                  logicMod3,    // its logical volume
                                  "Mod3",        // its name
                                  logicMotherMod,  // its mother  volume
                                  false,         // no boolean operations
                                  0);            // no particular field




    for (G4int i=1;i<StepNumbers;i++)
    {

        physiMod[2*(StepNumbers-1)+i] = new G4PVPlacement(0,
                                                          PositionMod[i],
                                                          logicMod[i],
                                                          "Modstep",
                                                          logicMod3,
                                                          false,
                                                          0);

    }

    //----------------------------------------------------------
    // Mother volume of the fourth modulator quarter
    //----------------------------------------------------------


    G4RotationMatrix rm4;
    rm4.rotateZ(270 *deg);

    G4ThreeVector positionMod4 = G4ThreeVector(0*cm,0*cm,0*cm);

    solidMod4 = new G4Tubs("Mod4",
                           innerRadiusOfTheTube,
                           outerRadiusOfTheTube,
                           hightOfTheTube0,
                           startAngleOfTheTube0/2.,
                           spanningAngleOfTheTube0);

    logicMod4 = new G4LogicalVolume(solidMod4,
                                    Mod0Mater, "Mod4",0,0,0);


    physiMod4 = new G4PVPlacement(G4Transform3D(rm4, positionMod4),
                                  logicMod4,
                                  "Mod4",
                                  logicMotherMod,
                                  false,
                                  0);


    for (G4int i=1;i<StepNumbers;i++)
    {
        physiMod[3*(StepNumbers-1)+i] = new G4PVPlacement(0,
                                                          PositionMod[i],
                                                          logicMod[i],
                                                          "Modstep",
                                                          logicMod4,
                                                          false,
                                                          0);
    }
    // Inform the kernel about the new geometry
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
    G4VisAttributes * red = new G4VisAttributes( G4Colour(1. ,0. ,0.));
    red-> SetVisibility(true);
    red-> SetForceSolid(true);
    logicMotherMod -> SetVisAttributes(G4VisAttributes::GetInvisible());

    logicMod1 ->SetVisAttributes(G4VisAttributes::GetInvisible());
    logicMod2 ->SetVisAttributes(G4VisAttributes::GetInvisible());
    logicMod3 ->SetVisAttributes(G4VisAttributes::GetInvisible());
    logicMod4 ->SetVisAttributes(G4VisAttributes::GetInvisible());

    for (G4int i=1;i<StepNumbers;i++)
    {
        logicMod[i] -> SetVisAttributes(red);
    }



}
/////////////////////////////////////////////////////////////////////////////
// Messenger values
/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::SetModulatorAngle(G4double angle)
{
    G4double rotationAngle = angle;
    rm -> rotateZ(rotationAngle);
    physiMotherMod -> SetRotation(rm);
    G4cout << "MODULATOR HAS BEEN ROTATED OF " << rotationAngle/deg
           << " deg" << G4endl;
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
/////////////////////////////////////////////////////////////////////////
// Change modulator material
void G4TPassiveProtonBeamLineGeometry::SetModulatorMaterial(G4String Material)
{
    if (G4Material* NewMaterial = G4NistManager::Instance()->FindOrBuildMaterial(Material, false) )
    {
        if (NewMaterial)
        {
            for(G4int i=1;i<StepNumbers;i++)
            {
                logicMod[i] -> SetMaterial(NewMaterial);
                //  G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
                G4RunManager::GetRunManager() -> GeometryHasBeenModified();

                //  G4cout<<(logicMod[i]->GetMaterial()->GetName())<<G4endl;
            }
            G4cout << "The material of the Modulator wheel has been changed to " << Material << G4endl;
        }
    }
    else
    {
        G4cout << "WARNING: material \"" << Material << "\" doesn't exist in NIST elements/materials"
                                                        " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl;
        G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl;


    }
}

////////////////////////////////////////////////////////////////////////////////
// Change modulator position in the beam line
void G4TPassiveProtonBeamLineGeometry::SetModulatorPosition(G4ThreeVector Pos)
{
    G4ThreeVector NewModulatorPos=Pos;
    physiMotherMod -> SetTranslation( NewModulatorPos);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout << "The modulator wheel is translated to"<<  NewModulatorPos/mm <<"mm " <<G4endl;

}
/////////////////////////////////////////////////////////////////////////////////
//change modulator inner raduis
void G4TPassiveProtonBeamLineGeometry::SetModulatorInnerRadius(G4double newvalue)
{
    solidMod1 -> SetInnerRadius(newvalue);
    solidMod2 -> SetInnerRadius(newvalue);
    solidMod3 -> SetInnerRadius(newvalue);
    solidMod4 -> SetInnerRadius(newvalue);
    for(G4int i=1;i<StepNumbers;i++)
    {
        solidMod[i] -> SetInnerRadius(newvalue);}
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout << "InnerRadius of the Modulator Wheel has been changed to :"
           << newvalue/mm<<" mm"<< G4endl;
}
/////////////////////////////////////////////////////////////////////////////////
//change modulator outer raduis
void G4TPassiveProtonBeamLineGeometry::SetModulatorOuterRadius(G4double newvalue)
{
    solidMod1 -> SetOuterRadius(newvalue);
    solidMod2 -> SetOuterRadius(newvalue);
    solidMod3 -> SetOuterRadius(newvalue);
    solidMod4 -> SetOuterRadius(newvalue);
    for(G4int i=1;i<StepNumbers;i++)
    {
        solidMod[i] -> SetOuterRadius(newvalue);}
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout << "OuterRadius of the Modulator Wheel has been changed to :"
           << newvalue/mm<<" mm"<<G4endl;
}



/////////////////////////// MESSENGER ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::SetRangeShifterXPosition(G4double value)
{
    physiRangeShifterBox -> SetTranslation(G4ThreeVector(value, 0., 0.));
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout << "The Range Shifter is translated to"<< value/mm <<"mm along the X axis" <<G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::SetRangeShifterXSize(G4double value)
{
    solidRangeShifterBox -> SetXHalfLength(value) ;
    G4cout << "RangeShifter size X (mm): "<< ((solidRangeShifterBox -> GetXHalfLength())*2.)/mm
           << G4endl;
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}

/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::SetFirstScatteringFoilXSize(G4double value)
{
    firstScatteringFoil -> SetXHalfLength(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout <<"The X size of the first scattering foil is (mm):"<<
             ((firstScatteringFoil -> GetXHalfLength())*2.)/mm
          << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::SetSecondScatteringFoilXSize(G4double value)
{
    secondScatteringFoil -> SetXHalfLength(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout <<"The X size of the second scattering foil is (mm):"<<
             ((secondScatteringFoil -> GetXHalfLength())*2.)/mm
          << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::SetOuterRadiusStopper(G4double value)
{
    solidStopper -> SetOuterRadius(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout << "OuterRadius od the Stopper is (mm):"
           << solidStopper -> GetOuterRadius()/mm
           << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::SetInnerRadiusFinalCollimator(G4double value)
{
    solidFinalCollimator -> SetInnerRadius(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout<<"Inner Radius of the final collimator is (mm):"
         << solidFinalCollimator -> GetInnerRadius()/mm
         << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void G4TPassiveProtonBeamLineGeometry::SetRSMaterial(G4String materialChoice)
{
    if (G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice, false) )
    {
        if (pttoMaterial)
        {
            rangeShifterMaterial  = pttoMaterial;
            logicRangeShifterBox -> SetMaterial(pttoMaterial);
            G4cout << "The material of the Range Shifter has been changed to " << materialChoice << G4endl;
        }
    }
    else
    {
        G4cout << "WARNING: material \"" << materialChoice << "\" doesn't exist in NIST elements/materials"
                                                              " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl;
        G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl;
    }
}
