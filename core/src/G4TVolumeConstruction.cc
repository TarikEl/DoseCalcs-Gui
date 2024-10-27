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

#include "globals.hh"

#include "TETPSEnergyDeposit.hh"
#include "G4TVolumeConstruction.hh"
//#include "G4TSD.hh"

#include "G4DrawVoxels.hh"

#include "G4RunManager.hh"
#include "G4tgbVolumeMgr.hh"
#include "G4tgbVolume.hh"
#include "G4tgbGeometryDumper.hh"
#if GDML_USE
#include "G4GDMLParser.hh"
#endif

#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVParameterised.hh"
#include "G4TVolumeBuilderUsingDICOM.hh"
#include "G4TCPPGeometryFormat.hh"
#include "G4TMyGeometry.hh"

#include "G4TVolumeBuilderUsingVoxel.hh"

#include "G4TGeometryBasedBiasingOperator.hh"

#include "G4UImanager.hh"
#include "G4VVisManager.hh"

#include "G4IStore.hh"
#include "G4UserLimits.hh"

#include <iostream>
#include <fstream>
#include <ios>

#include "G4TNuclearReactorGeometry.hh"

//#ifdef G4MPI_USE
//#include "G4AutoLock.hh"
//namespace {G4Mutex	mutex = G4MUTEX_INITIALIZER;}
//#endif


#include "G4Sphere.hh"

#include "G4TDCMHandler.hh"
#ifdef DCMTK_USE
#include "DicomFileMgr.hh"
#endif

extern std::string getFileNameFromPath(std::string const & path, std::string const & delims = "/\\");
extern std::string getFileExt(const std::string& s);

G4TVolumeConstruction::G4TVolumeConstruction()
{
#if VERBOSE_USE
    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  G4TVolumeConstruction  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;
#endif

    //CreateICRPSAFsReferenceDataFile();


    InternalSourcePosition = G4ThreeVector();

    DicomCTName = "CTDATA.txt";
    DicomPETName = "PETDATA.txt";

    GeometryFileType = "Construct";
    GeometrySymbol = "phantom0";
    ThreadsNumber = 0;
    MPISimulationNum = "o";

    MaterialNameAsRegionName = true;

    GenerateCrossSectionTableFlag = false;

    BoxLVName = "BoxVolume";
    DataDirectoryPath = "EventsData";
    //ScriptsDirectoryPath  = "Scripts/";

    G4Density_to_gPerCm3 = cm3/g; // or 1/(6.24151e+18); // used for setting the values to DICOMHANDLER
    G4Density_to_kgPerMm3 = mm3/kg; // or 1/(6.24151e+24); // used to fill Volumes needed for calculation as OrganNameMassMap; CopyNumberMassSize
    G4Density_to_kgPerMm3ToGPercm3 = (kg/mm3)*(cm3/g);

    G4Mass_to_Kg = 1/kg; // or 1/(6.24151e+24); // used by construct Type for volumes and world because we can get Mass of a logical Volume

    DicomPixelsCompression = 1;
    PlanesToVisualize = "all";

    RegionFromInputInc = 0; // used for Dcm region to add data segmentation to the regions

    OtherDataMessenger = new G4TMessenger(this);
    geometryMessenger = new G4TGeometryMessenger(this);
    PhysicsMessenger = new G4TPhysicsMessenger(this);
    PGAMessenger = new G4TPrimaryGeneratorMessenger(this);

    UseDicomCumAct = false;
    GenerateVoxelsResuls = "no";
    //VOXTET_USE = false;

    CPPLogVolAreBuilt = false;
    dumpGeom = false;

    //Vacuum = new G4Material("Vacuum", 1., 1.01*g/mole, universe_mean_density,kStateGas,2.73*kelvin,3.e-18*pascal);

    //CreateICRPSAFsReferenceDataFile();
}
G4TVolumeConstruction::~G4TVolumeConstruction()
{
    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;

    delete OtherDataMessenger;
    delete geometryMessenger;
    delete PhysicsMessenger;
    delete PGAMessenger;

    if(tetData != 0){delete tetData;}

}

G4VPhysicalVolume* G4TVolumeConstruction::Construct()
{
    TestAndShowUserInputs();

    // Create data files Name

    setSourcePositionFileName();
    setSourceEnergyFileName();
    setSourceMomDirFileName();

    //ConvertVoxelPosToXYZVoxelIDs("Z", "Max", 2.6);
    //ConvertVoxelPosToXYZVoxelIDs("Z", "Min", 2.6);

    //CalcutateSizeOfVoxelizedArrays();

    //ConstructICRPPhantomVolume(); return 0;



    if(GeometryFileType == "GDML"){
#if GDML_USE
#else
        G4String msg =  "DoseCalcs should be built with \"-DWITH_GDML_USE=ON\" to construct geometry from .gdml files "; G4Exception("Geometry Data", "1", FatalErrorInArgument, msg.c_str());
#endif
    }

    if(GeometryFileType == "ADD"){
        CreateVolumesData();
    }
    else if(GeometryFileType == "GDML" || GeometryFileType == "TEXT" || GeometryFileType == "CPP" || GeometryFileType == "Construct" || GeometryFileType == "STL" || GeometryFileType == "MyGeometry" ){
        ConstructVolumes();

        CreateVolumesData();

        ShowBoxVolume();
        makeVolumeVisualization();
    }
    else if(GeometryFileType == "VOXEL"){
        ConstructVOXELVolume();
    }
    else if(GeometryFileType == "DICOM"){
        ConstructDICOMVolume();
    }
    else if(GeometryFileType == "VoxIDs") {
        ConstructPhantomFromVoxIds();
    }
    else if(GeometryFileType == "Generate") {

        //For ICRP voxelized phantoms materials Generation
        GenerateICRPMaterialsCommands();

        //For TET phantoms materials Generation

        return 0;
    }
    else if(GeometryFileType == "TET"){
        ConstructTETPhantomVolume();
    }
    else {
        return 0;
    }

    setRankDataForMPIMode();

    return WorldPhysicalVolume;
}

void G4TVolumeConstruction::ConstructVolumes(){

    if(GeometryFileType == "GDML"){
#if GDML_USE
        G4GDMLParser parser;
        parser.Read( GeometryPath , false );  //false to eliminate the xchema validation because it print a lot of lines of error validation
        WorldPhysicalVolume = parser.GetWorldVolume();
#else
        G4String msg =  "DoseCalcs should be built with \"-DWITH_GDML_USE=ON\" to construct geometry from .gdml files "; G4Exception("Geometry Data", "1", FatalErrorInArgument, msg.c_str());
#endif
    }
    else if(GeometryFileType == "TEXT"){
        G4tgbVolumeMgr* G4tgbVolumeMgrObj;
        G4tgbVolumeMgrObj = G4tgbVolumeMgr::GetInstance();
        G4tgbVolumeMgrObj->AddTextFile(GeometryPath);
        WorldPhysicalVolume = G4tgbVolumeMgrObj->ReadAndConstructDetector();
    }
    else if(GeometryFileType == "CPP"){

        G4TCPPGeometryFormat * CGF = new G4TCPPGeometryFormat();
        WorldPhysicalVolume = CGF->ConstructPhysicalVolume();

    }
    else if(GeometryFileType == "MyGeometry"){

        G4TMyGeometry * CGF = new G4TMyGeometry();
        WorldPhysicalVolume = CGF->ConstructPhysicalVolume();

    }
    else if(GeometryFileType == "Construct"){

    }
    //else if(GeometryFileType == "NR"){
    //    G4TNuclearReactorGeometry* NucRea = new G4TNuclearReactorGeometry();
    //    WorldPhysicalVolume = NucRea->ConstructNuclearReactor();
    //}

    // CreateVolumesData();

    // //boxCenterPos = CreatedPositionOrgans[SourceRegionName];
    // //G4cout << "boxCenterPos " << boxCenterPos <<G4endl;

    // ShowBoxVolume();
    // makeVolumeVisualization();
    // setRankDataForMPIMode();
    // //CreateMaterialsGDMLTags();

    //return WorldPhysicalVolume ;
}
void G4TVolumeConstruction::ConstructDICOMVolume(){
#ifdef DCMTK_USE

    // to set materials ordered in MaterialIndices[] map
    std::map<G4double,G4String> DensityMaterial;
    std::vector<G4double> densities;
    //std::vector<G4String> names;
#if VERBOSE_USE
    G4cout << "\n\nUser materials : \n" <<G4endl;
#endif
    for ( auto it = CreatedMaterials.begin(); it != CreatedMaterials.end(); ++it  ){
        densities.push_back(it->second->GetDensity());
        //names.push_back(it->first);
        DensityMaterial[it->second->GetDensity()] = it->first ;
        RegionNameColour[it->second->GetName()] = G4Colour( (G4double) G4UniformRand(),(G4double)G4UniformRand(),(G4double)G4UniformRand(), 1.);

        G4String jj =it->first; jj.resize(32);
#if VERBOSE_USE
        G4cout << std::setw(33) << std::left << jj.c_str() << std::setw(15) << std::left << it->second->GetDensity() * G4Density_to_gPerCm3  <<G4endl;
#endif
    }
    sort(densities.begin(), densities.end());
#if VERBOSE_USE
    G4cout << "\n\nUser Density ordered materials to be used by DICOM Reader: \n" <<G4endl;
    G4cout << std::setw(33) << std::left << "Material Name" << std::setw(15) << std::left << "Density(g/cm3) " <<G4endl;
#endif
    for ( G4int h = 0; h < densities.size(); h++  ){
        //DcmMaterials.push_back(CreatedMaterials[DensityMaterial[densities[h]]]);
        G4String jj = CreatedMaterials[DensityMaterial[densities[h]]]->GetName(); jj.resize(32);
        MaterialIndices[CreatedMaterials[DensityMaterial[densities[h]]]->GetDensity()*G4Density_to_gPerCm3] = CreatedMaterials[DensityMaterial[densities[h]]]->GetName();// /(6.24151e+15)
#if VERBOSE_USE
        G4cout << std::setw(33) << std::left << jj.c_str() << std::setw(15) << std::left << CreatedMaterials[DensityMaterial[densities[h]]]->GetDensity() * G4Density_to_gPerCm3 <<G4endl;
#endif
    }


    G4bool fpf = true;

    for(G4int R = 0 ; R < DicomFileTypesVec.size() ; R++){

        DicomFileType = DicomFileTypesVec[R];
        DicomFilesDirPath = DicomFilePathsVec[R];

        G4cout << " From function ConstructDICOMVolume : " << DicomFileType << " : " << DicomFilesDirPath << G4endl;

        if(DicomFileType == "CT"){
            DicomOutTextName = DicomFilesDirPath+"/CTDATA.txt";
#if VERBOSE_USE
            G4cout << "\n\n========= Start Reading, processing, converting Dicom files to Text files for CT modality ======================= \n"<<G4endl;
#endif
        }
        else if(!UseDicomCumAct){
            continue;
        }

        else if(DicomFileType == "PET"){

            DicomOutTextName = DicomFilesDirPath+"/"+DicomPETName;

#if VERBOSE_USE
            G4cout << "\n\n========= Start Reading, processing, converting Dicom files to Text files for PET modality ======================= \n"<<G4endl;
#endif
        }

        DicomFileMgr* theFileMgr = 0;
        theFileMgr = DicomFileMgr::GetInstance();
        //theFileMgr->setDicomType(DicomFileType);
        theFileMgr->setDicomFilesDirPath(DicomFilesDirPath);
        theFileMgr->setDicomOutTextName(DicomOutTextName);
        DicomPixelsCompression = SeriesCompressionXY[R];

        theFileMgr->SetCompression(DicomPixelsCompression+"");
        //theFileMgr->Convert("Data.dat.new_dens");
        theFileMgr->GenerateDICOMTextGeometryFile();

        if(DicomFileType == "PET"){
            Serie_Order_DataMap[R-1] = theFileMgr->getSerieAcquisitionsDataVec();
            if(fpf){

                PETGeneralData = theFileMgr->getPETGeneralData();
                PETGeneralDateData = theFileMgr->getPETGeneralDateData();

                fpf = false;
            }
        }

        theFileMgr = nullptr;
    }

#else
    G4TDCMHandler* G4TDCMHandlerObject = new G4TDCMHandler();
    G4TDCMHandlerObject->CheckFileFormat();
    delete G4TDCMHandlerObject;
#endif

    getDicomDataAndConstruct();

    // to fill the minmaxVectors by the new values because the number of Voxels on x y z not known by the user
    for(G4int RNL = 0 ; RNL < DcmRegionsNames.size() ; RNL++){
        if(PhantomLimits[RNL] == "all"){
            DcmRegionsMinX[RNL] = 0 ;
            DcmRegionsMaxX[RNL] = VoxXNumber-1 ;
            DcmRegionsMinY[RNL] = 0 ;
            DcmRegionsMaxY[RNL] = VoxYNumber-1 ;
            DcmRegionsMinZ[RNL] = 0 ;
            DcmRegionsMaxZ[RNL] = VoxZNumber-1 ;
        }
    }
    for(G4int RNL = 0 ; RNL < DcmRegionsNames.size() ; RNL++){

        //G4cout << DcmRegionsNames[RNL] << " " << DcmRegionsMinX[RNL] << " " << DcmRegionsMaxX[RNL] << " " << DcmRegionsMinY[RNL] << " " << DcmRegionsMaxY[RNL] << " " << DcmRegionsMinZ[RNL] << " " << DcmRegionsMaxZ[RNL] << " \n"<<G4endl;

        // segmentation data
        VoxRegionName = DcmRegionsNames[RNL];
        VoxRegionMinX = DcmRegionsMinX[RNL];
        VoxRegionMaxX = DcmRegionsMaxX[RNL];
        VoxRegionMinY = DcmRegionsMinY[RNL];
        VoxRegionMaxY = DcmRegionsMaxY[RNL];
        VoxRegionMinZ = DcmRegionsMinZ[RNL];
        VoxRegionMaxZ = DcmRegionsMaxZ[RNL];
        UseDcmRegionMinDensity = UseDcmRegionsMinDensityMap[RNL];
        UseDcmRegionMaxDensity = UseDcmRegionsMaxDensityMap[RNL];
        DcmRegionMinDensity = DcmRegionsMinDensityMap[RNL];
        DcmRegionMaxDensity = DcmRegionsMaxDensityMap[RNL];

        CreateRegionVoxelsData();
    }

    CreateDefaultVOXELRegionData();

    ShowMaterialRegionVoxelsData();

    //setRankDataForMPIMode();
    //return ConstructVoxeDcmGeometry();
    ConstructVoxeDcmGeometry();
}
void G4TVolumeConstruction::ConstructVOXELVolume(){

    CreateDefaultVOXELRegionData();
    CreateRegionsDataFile();
    //setRankDataForMPIMode();

    //return ConstructVoxeDcmGeometry();
    ConstructVoxeDcmGeometry();
}
void G4TVolumeConstruction::ConstructPhantomFromVoxIds(){
    ReadVoxelsIDsAndFillCNMatIDsMassColour();
    if(MaterialNameAsRegionName == false){
        GenerateDataForVoxelsIdsFilePhantom();
        CreateDefaultVOXELRegionData();
        CreateRegionsDataFile();
    }
    ShowMaterialRegionVoxelsData();
    //setRankDataForMPIMode();
    //return ConstructVoxeDcmGeometry();
    ConstructVoxeDcmGeometry();
}
void G4TVolumeConstruction::ConstructTETPhantomVolume(){

    //Original
    //GenerateDataFromTETPhantomFiles();
    //setRankDataForMPIMode();
    // //return ConstructTETGeometry();
    //ConstructTETGeometry();


    // modified
    GenerateDataFromTETPhantomFiles();
    //return ConstructTETGeometry();
    ConstructTETGeometry();
}


void G4TVolumeConstruction::ConstructTETGeometry(){

    // Define the phantom container (10 cm-margins from the bounding box of phantom)
    //

    G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    ContSolidVoll = new G4Box("Container", phantomSize.x() * 0.5 + 10.*cm, phantomSize.y() * 0.5 + 10.*cm,phantomSize.z() * 0.5 + 10.*cm);

    // added
    if(UseInternalSourceVolume == true){
        G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance() ;
        for (size_t iLV = 0; iLV < pvs->size(); iLV++ ) {
            G4String nmm = (*pvs)[iLV]->GetLogicalVolume()->GetName();
            //G4cout << (*pvs)[iLV]->GetLogicalVolume()->GetName() << " " << nmm << G4endl;
            if(strstr(nmm.c_str(),InternalSourceName)){
                InternalSourceLogicalSolid = (*pvs)[iLV]->GetLogicalVolume()->GetSolid();
                InternalSourceLogicalVolume = (*pvs)[iLV]->GetLogicalVolume();
                InternalSourcePosition = getRegionAbsolutePosition(InternalSourceName);
                WorldPhysicalVolume->GetLogicalVolume()->RemoveDaughter((*pvs)[iLV]);
            }
        }

        G4RotationMatrix* rm = new G4RotationMatrix();
        NewContSolidVoll = new G4SubtractionSolid( "Container" , ContSolidVoll, InternalSourceLogicalSolid, rm, InternalSourcePosition);
        ContLogicalVoll = new G4LogicalVolume(NewContSolidVoll, vacuum, "Container");
        placeInternalSourceVolume();
    }else{
        ContLogicalVoll = new G4LogicalVolume(ContSolidVoll, vacuum, "Container");
    }

    //ContLogicalVoll = new G4LogicalVolume(ContSolidVoll, vacuum, "Container");
    ContPhysicalVoll = new G4PVPlacement(0, G4ThreeVector(), ContLogicalVoll, "Container", WorldPhysicalVolume->GetLogicalVolume(), false, 0);
    //G4RotationMatrix* rm = new G4RotationMatrix(); rm->rotateX(G4ThreeVector().getX()); rm->rotateY(G4ThreeVector().y()); rm->rotateZ(G4ThreeVector().z());
    //ContPhysicalVoll = new G4PVPlacement(rm, G4ThreeVector(), ContLogicalVoll, "Container", WorldPhysicalVolume->GetLogicalVolume(), false, 0);
    ContLogicalVoll->SetOptimisation(TRUE);
    ContLogicalVoll->SetSmartless( 0.5 ); // for optimization (default=2)

    // Define the tetrahedral mesh phantom as a parameterised geometry
    //
    // solid and logical volume to be used for parameterised geometry

    G4VSolid* tetraSolid = new G4Tet("TetSolid",
                                     G4ThreeVector(),
                                     G4ThreeVector(1.*cm,0,0),
                                     G4ThreeVector(0,1.*cm,0),
                                     G4ThreeVector(0,0,1.*cm));

    tetLogic = new G4LogicalVolume(tetraSolid, vacuum, "TetLogic");
/*
    G4TTETParameterisation* param2 =  new G4TTETParameterisation(tetData);
    if(UseVoxelsColour == true){ param2->setUseLogVolColour(true);}

    physical volume (phantom) constructed as parameterised geometry
    new G4PVParameterised("wholePhantom",tetLogic,ContLogicalVoll,
                          kUndefined,
                          tetData->GetNumTetrahedron(),
                          param2);
*/
    // physical volume (phantom) constructed as parameterised geometry

    new G4PVParameterised("wholePhantom",tetLogic,ContLogicalVoll,
                          kUndefined,tetData->GetNumTetrahedron(),
                          new TETParameterisation(tetData));

    // added
    //setSourceVolumeData(G4ThreeVector(12*cm,32*cm,62*cm));
    //setSourceVolumeData(G4ThreeVector(0*cm,0*cm,0*cm));
    // // we make it here because we have to add a volume before creting replicate then
    // // getting its position

    //setRankDataForMPIMode();

    //return WorldPhysicalVolume;

}
void G4TVolumeConstruction::GenerateDataFromTETPhantomFiles()
{

    // for tetrahedral always = true because we use materials ID in Stepping and HitMap for estimating edep
    MaterialNameAsRegionName = true;


    tetData = new G4TTETModelImport();

    tetData->setIdNameMap(MaterialIDName); // readed from macros file and not from materials files as it is given in the mesh-type phantom example data files
    tetData->setNameMatMap(CreatedMaterials);

    tetData->setRegionSegmentationDataVectors(DcmRegionsNames,
                                              //UseVoxelMatForSegMap,
                                              //VoxelMatForSegMap,
                                              UseDcmRegionsMinDensityMap,
                                              UseDcmRegionsMaxDensityMap,
                                              DcmRegionsMinDensityMap,
                                              DcmRegionsMaxDensityMap);

    tetData->GenerateDataforTET();

    OrganNamesVector = tetData->GetOrganNamesVector();
    OrganNameMassMap = tetData->GetOrganNameMassMap();
    OrganNameDensityMap = tetData->GetOrganNameDensityMap();
    OrganNameVolumeMap = tetData->GetOrganNameVolumeMap();

    phantomBoxMin   = tetData -> GetPhantomBoxMin();
    phantomBoxMax   = tetData -> GetPhantomBoxMax();
    nOfTetrahedrons = tetData -> GetNumTetrahedron();
    phantomSize     = tetData -> GetPhantomSize();

#if VERBOSE_USE
    tetData->PrintMaterialInfomation();
#endif
    ShowMaterialRegionVoxelsData();
}

void G4TVolumeConstruction::placeInternalSourceVolume(){
    new G4PVPlacement(0, InternalSourcePosition, InternalSourceLogicalVolume, "Tumour", WorldPhysicalVolume->GetLogicalVolume(), false, 0, true);
}
void G4TVolumeConstruction::setInternalSourceData(G4String n ){
    UseInternalSourceVolume = true;
    InternalSourceName = n;
}

// called in the begining of construct method to see the values seten by commands
void G4TVolumeConstruction::TestAndShowUserInputs(){

    G4String msg = "";

#if VERBOSE_USE
    G4cout<<"\n\n========= Geometry, Source, Score Inputs ====================" <<G4endl;

    G4cout<<"\n\n >> Geometry Type: " << GeometryFileType <<G4endl;
    G4cout<<" >> Geometry File: " << GeometryPath <<G4endl;
    G4cout<<" >> SourcePosition: " << SourcePosition <<G4endl;
    G4cout<<" >> WorldHalfSize: " << WorldHalfSize <<G4endl;
    G4cout<<" >> WorldMaterialName: " << WorldMaterialName << G4endl;
    G4cout<<" >> Use Material Name As Region Name: " << MaterialNameAsRegionName <<"\n\n"<<G4endl;

    G4cout<<"\n\n >> Source Particle Name: " << ParticleName <<G4endl;
    G4cout<<" >> Source Position Dist: " << SourceType <<G4endl;
    G4cout<<" >> Source Energy Dist: " << EnergyDistribution <<G4endl;
    G4cout<<" >> Source Momentum Direction Dist: " << MomDirDistribution <<G4endl;
    //G4cout<<" >> Source Events Data Number: " << NumberOfGenPointsToSave <<G4endl;

    //G4cout<<"\n\n >> Quantities To Score: " << variable_To_Score <<G4endl;
    //G4cout<<" >> Regions To Score: " << organs_to_score << G4endl;
    G4cout<<"\n\n >> Simulation Number On Ranks: " << MPISimulationNum <<G4endl;
    G4cout<<" >> Number Of Threads: " << ThreadsNumber <<"\n\n"<<G4endl;

    //if(WorldPhysicalVolume == nullptr){msg = "construct the world geometry"; G4Exception("Source Data", "1", FatalException, msg.c_str());}
    //if(SourceRegionsNamesValues.size() == 0){msg = "set a source region name"; G4Exception("Source Data", "3", FatalException, msg.c_str());}
    //if(SourceEnergiesValues.size() == 0){msg = "set a source primary energy"; G4Exception("Source Data", "3", FatalException, msg.c_str());}
    //if(SourceMomDirsValues.size() == 0){msg = "set a source primary momentum direction "; G4Exception("Source Data", "3", FatalException, msg.c_str());}

    /*
    G4cout<<" >> Dicom pixels compression number " << DicomPixelsCompression <<G4endl;
    G4cout<<" >> SourceType " << SourceType <<G4endl;
    G4cout<<" >> SourceRegionName " << SourceRegionName <<G4endl;
    G4cout<<" >> ShowBox " << ShowBox <<G4endl;
    G4cout<<" >> BoxDimGene " << BoxDimGene <<G4endl;
    G4cout<<" >> EnergyDistribution " << EnergyDistribution <<G4endl;
    G4cout<<" >> MomDirDistribution " << MomDirDistribution <<G4endl;
    G4cout<<" >> NumberOfPointsToSave " << NumberOfGenPointsToSave <<G4endl;
    G4cout<<" >> GaussMean " << GaussMean <<G4endl;
    G4cout<<" >> GaussSDev " << GaussSDev <<G4endl;
    G4cout<<" >> UniformEmin " << UniformEmin <<G4endl;
    G4cout<<" >> UniformEmax " << UniformEmax <<G4endl;
    G4cout<<" >> RayleighEmax " << RayleighEmax <<G4endl;
    G4cout<<" >> ParticleName " << ParticleName <<G4endl;
    G4cout<<" >> MonoEnergy " << MonoEnergy <<G4endl;
    G4cout<<" >> Theta " << Theta <<G4endl;
    G4cout<<" >> Phi " << Phi <<G4endl;

    G4cout<<"\n >> organs_to_score " << organs_to_score <<G4endl;
    G4cout<<" >> variable_To_Score " << variable_To_Score << G4endl;
    G4cout<<" >> ThreadsNumber " << ThreadsNumber <<G4endl;
    G4cout<<" >> MPISimulationNum " << MPISimulationNum <<G4endl;

    //G4cout<<" >> EventNumberInOneThreads " << EventNumberInOneThreads <<G4endl;
    //G4cout<<" >> BatchsNumber " << BatchsNumber <<G4endl;
    G4cout<<" >> AccuracyCalculationLevel " << AccuracyCalculationLevel <<G4endl;

    G4cout<<"\n >> Graphs_Data " << Graphs_Data <<G4endl;
    G4cout<<" >> compare_type " << compare_type <<G4endl;
    G4cout<<" >> graphs_Ext " << graphs_Ext <<G4endl;
    G4cout<<" >> ref_File_Path " << ref_File_Path <<G4endl;
    G4cout<<" >> ref_Name " << ref_Name <<G4endl;
*/
#endif

}
void G4TVolumeConstruction::CreateVolumesData(){

    //G4cout << " From function : " << __FUNCTION__ << G4endl;

    G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance() ;

    // to make the name pf physical volume the same as logical volume and remove world Volume from name as in human phantoms gdml files
    for (size_t iLV = 0; iLV < pvs->size(); iLV++ ) {
        G4String nmm = (*pvs)[iLV]->GetLogicalVolume()->GetName();
        //G4cout << (*pvs)[iLV]->GetLogicalVolume()->GetName() << " " << nmm << G4endl;
        if(strstr(nmm.c_str(),"Volume")){
            (*pvs)[iLV]->GetLogicalVolume()->SetName(nmm.replace(nmm.find("Volume"), 6, ""));
            (*pvs)[iLV]->SetName((*pvs)[iLV]->GetLogicalVolume()->GetName());
        }
    }

    G4String WorldVolName = "World";
    WorldVolName = WorldPhysicalVolume->GetLogicalVolume()->GetName();

    G4String LVName= "";
    G4String PVName= "";
#if VERBOSE_USE
    G4cout<<"\n\n========= Created Volumes in Geometry ====================\n\n" <<G4endl;
#endif
    //std::fprintf(stdout,"------------------------------------------------------------------------------------------------------------------");
    //std::fprintf(stdout,"\n%-3s%-25s%-20s%-20s%-20s%-25s%-1s\n", "|", "Volume Name" , "| Mass(Kg)", "| Volume(cm3)", "| Density(g/cm3)" , "| X Y Z (mm)" ,"|" );
    //std::fprintf(stdout,"------------------------------------------------------------------------------------------------------------------\n");

#if VERBOSE_USE
    G4cout << "------------------------------------------------------------------------------------------------------------------"<< G4endl;
    G4cout << std::setw(25) << std::left << "| Volume Name"
           << std::setw(20) << std::left << "| Mass(Kg)"
           << std::setw(20) << std::left << "| Volume(cm3)"
           << std::setw(20) << std::left << "| Density(g/cm3)"
           << std::setw(20) << std::left << "| Position(mm)"
              //<< std::setw(20) << std::left << "| Rotation(rad)" << " |"
           << G4endl;
    G4cout << "------------------------------------------------------------------------------------------------------------------"<< G4endl;
#endif

    for (size_t iLV = 0; iLV < pvs->size(); iLV++ ) {

        G4VPhysicalVolume* PV = (*pvs)[iLV];
        //G4LogicalVolume* LV = (*nn)[iLV];
        PVName = PV->GetName();
        LVName = PV->GetLogicalVolume()->GetName();

        // this will be considered after setting the pre calculated results

        bool isIn = false; for ( int df = 0 ; df < OrganNamesVector.size(); df++  ){ if(LVName == OrganNamesVector[df] ){ isIn = true; }}
        if(isIn == true ){/*G4cout << " yes is in " << LVName << G4endl; */ continue; } else{OrganNamesVector.push_back(LVName);}
        //OrganNamesVector.push_back(LVName);


        // Visualization Attributes
        G4Colour colour = G4Colour((G4double)G4UniformRand(), (G4double)G4UniformRand(), (G4double)G4UniformRand());
        G4VisAttributes* organVisAtt = new G4VisAttributes(colour);
        organVisAtt->SetForceSolid(false);
        PV->GetLogicalVolume()->SetVisAttributes(organVisAtt);

        G4double Density1 = PV->GetLogicalVolume()->GetMaterial()->GetDensity() * G4Density_to_gPerCm3; // density in g/cm3;
        //G4cout << PV->GetLogicalVolume()->GetMaterial()->GetName() << " " << PV->GetLogicalVolume()->GetMaterial()->GetDensity() << " " << LVName.c_str()  << G4endl;

        if(OrganNameDensityMap[LVName] == 0. || __isinf(OrganNameDensityMap[LVName]) || __isnan(OrganNameDensityMap[LVName])){
        //if(OrganNameDensityMap[LVName] == 0.){
            OrganNameDensityMap[LVName] = Density1;
        }
        else{
            OrganNameDensityMap[LVName] = (Density1 + OrganNameDensityMap[LVName])/2;
        }
        G4double Volume1 = PV->GetLogicalVolume()->GetSolid()->GetCubicVolume()/cm3;
        OrganNameVolumeMap[LVName] += Volume1; // in cm3
        G4double Mass1 = Volume1*Density1;
        //if(createdOrganMass[LVName] == 0.){}

        OrganNameMassMap[LVName] += Mass1/1e3; // from g to Kg;
        OrganNamePositionMap[LVName] = getRegionAbsolutePosition(LVName);

        if(LVName != "World"){
            AllGeometryMass += OrganNameMassMap[LVName];
            AllGeometryVolume += OrganNameVolumeMap[LVName];
#if VERBOSE_USE
            //std::fprintf(stdout,"%-3s%-25s%-20e%-20e%-20e%-7.2f%-1s%-7.2f%-1s%-7.2f%-1s\n", "|", LVName.c_str(), (double)OrganNameMassMap[LVName] , (double)OrganNameVolumeMap[LVName], (double)OrganNameDensityMap[LVName], (double)OrganNamePositionMap[LVName].getX(), "  " , (double)OrganNamePositionMap[LVName].getY(),"  ", (double)OrganNamePositionMap[LVName].getZ(),"|");
            G4cout << std::setw(25) << std::left << LVName.c_str()
                   << std::setw(20) << std::left << OrganNameMassMap[LVName]
                      << std::setw(20) << std::left << OrganNameVolumeMap[LVName]
                         << std::setw(20) << std::left << OrganNameDensityMap[LVName]
                            << std::setw(1) << std::left << (OrganNamePositionMap[LVName])
                               //<< std::setw(20) << std::left << OrganNameRotMatrixMap[LVName] << " |"
                            << G4endl;
#endif
        }
    }

#if VERBOSE_USE
    //std::fprintf(stdout,"------------------------------------------------------------------------------------------------------------------\n");
    G4cout << "------------------------------------------------------------------------------------------------------------------"<< G4endl;
    G4cout << " Phantom Mass " << AllGeometryMass << " kg" << G4endl;
    G4cout << " Phantom Volume " << AllGeometryVolume << " cm3" << G4endl;
#endif

    for(G4int RNL = 0 ; RNL < OrganNamesVector.size() ; RNL++){
        //OrganNamePositionMap[OrganNamesVector[RNL]] = getRegionAbsolutePosition(OrganNamesVector[RNL]);
        //G4cout << "tttttttttttttttttttt " << OrganNamesVector[RNL] << " " << OrganNamePositionMap[OrganNamesVector[RNL]] << G4endl;
    }


}
void G4TVolumeConstruction::makeVolumeVisualization()
{
    if(ShowBox != "yes" && TestPointsPositions != "yes"){
        //return;
    }

    // to use correctly the "allregions" region in calculation
    for (size_t iLV = 0; iLV < SourceRegionsNamesValues.size(); iLV++ ) {
        if(SourceRegionsNamesValues[iLV]=="allregions"){
            SourceRegionsNamesToBeIgnoredValues.push_back("World");
            if(SourceRegionsNamesValues.size()>0){
                SourceRegionName = SourceRegionsNamesValues[0];
                BoxDimGene = SourceRegionsBoxDimValues[0];
            }
            if(BoxDimGene.getX()>WorldHalfSize.getX()){
                BoxDimGene.setX(WorldHalfSize.getX());
            }
            if(BoxDimGene.getY()>WorldHalfSize.getY()){
                BoxDimGene.setY(WorldHalfSize.getY());
            }
            if(BoxDimGene.getZ()>WorldHalfSize.getZ()){
                BoxDimGene.setY(WorldHalfSize.getZ());
            }
        }
    }

    G4LogicalVolumeStore* nn = G4LogicalVolumeStore::GetInstance() ;
    for (size_t iLV = 0; iLV < nn->size(); iLV++ ) {
        G4LogicalVolume* LV = (*nn)[iLV];
        G4String LVN = LV->GetName();
        if(LVN == SourceRegionName || LVN == BoxLVName){
            G4Colour colour = G4Colour((G4double)G4UniformRand(), (G4double)G4UniformRand(), (G4double)G4UniformRand());
            G4VisAttributes* VolumeVisAtt = new G4VisAttributes(colour);
            // Visualization Attributes
            VolumeVisAtt->SetVisibility(true);

            if(LVN == SourceRegionName && TestPointsPositions != "yes"){
                VolumeVisAtt->SetForceSolid(true);

            }
            LV->SetVisAttributes(VolumeVisAtt);

        }else{
            G4Colour colour = G4Colour((G4double)G4UniformRand(), (G4double)G4UniformRand(), (G4double)G4UniformRand());
            G4VisAttributes* VolumeVisAtt = new G4VisAttributes(colour);
            // Visualization Attributes
            VolumeVisAtt->SetVisibility(true);
            VolumeVisAtt->SetForceSolid(true);

            for (int ss = 0; ss < VolumesNotVisualized.size(); ss++ ) {
                if(LVN == VolumesNotVisualized[ss]){
                    VolumeVisAtt->SetForceSolid(false);
                    break;
                }
            }
            LV->SetVisAttributes(VolumeVisAtt);
        }
    }
}
G4ThreeVector G4TVolumeConstruction::getRegionAbsolutePosition(G4String RegionName){

    G4ThreeVector vectPos1(0,0,0), vectPos2(0,0,0), vectPos3(0,0,0), vectPos4(0,0,0);
    G4RotationMatrix rm, rm1, rm2, rm3, rm4;
    //G4ThreeVector rot(0,0,0);

    //double abstheta, absphi, abspsi ;

    G4VPhysicalVolume* Phy1 = 0;
    G4VPhysicalVolume* Phy2 = 0;
    G4VPhysicalVolume* Phy3 = 0;
    G4VPhysicalVolume* Phy4 = 0;

    G4ThreeVector vectPos(0,0,0);

    if(RegionName == "AllRegion"){
        return vectPos;
    }

    for (G4int aa = 0 ; aa < WorldPhysicalVolume->GetLogicalVolume()->GetNoDaughters() ; aa++) {

        WorldPhysicalVolume->GetLogicalVolume()->GetDaughter(aa);
        Phy1 = WorldPhysicalVolume->GetLogicalVolume()->GetDaughter(aa);
        vectPos1 = Phy1->GetObjectTranslation();
        rm1 = Phy1->GetObjectRotationValue();
        //theta = Phy1->GetObjectRotationValue().getTheta(), phi = Phy1->GetObjectRotationValue().getPhi() , psi = Phy1->GetObjectRotationValue().getPsi();

        if(Phy1->GetLogicalVolume()->GetName() == RegionName){

            vectPos = vectPos1;
            rm.transform(rm1);
            //rm.setTheta(Phy1->GetObjectRotationValue().getTheta()); rm.setPhi(Phy1->GetObjectRotationValue().getPhi()); rm.setPsi(Phy1->GetObjectRotationValue().getPsi());

            OrganNameRotMatrixMap[RegionName] = rm;
            //G4cout << " ---- " << RegionName << " Pos:" << vectPos  << " ---- RotationMatrix: " << rm <<  " ---- Angles:" << rm.getTheta() << " " << rm.getPhi() << " " << rm.getPsi() << G4endl;
            return vectPos;
        }

        for (G4int bb =0 ; bb < Phy1->GetLogicalVolume()->GetNoDaughters() ; bb++) {

            Phy2 = Phy1->GetLogicalVolume()->GetDaughter(bb);
            vectPos2 = Phy2->GetObjectTranslation();
            rm2 = Phy2->GetObjectRotationValue();
            if(Phy2->GetLogicalVolume()->GetName() == RegionName){

                vectPos = vectPos2.transform(Phy1->GetObjectRotationValue()) + vectPos1;
                rm.transform(rm1); rm.transform(rm2);

                //rm.setTheta(Phy1->GetObjectRotationValue().getTheta()); rm.setPhi(Phy1->GetObjectRotationValue().getPhi()); rm.setPsi(Phy1->GetObjectRotationValue().getPsi());
                //rm.setTheta(Phy2->GetObjectRotationValue().getTheta()); rm.setPhi(Phy2->GetObjectRotationValue().getPhi()); rm.setPsi(Phy2->GetObjectRotationValue().getPsi());

                OrganNameRotMatrixMap[RegionName] = rm;
                //G4cout << " ---- " << RegionName << " Pos:" << vectPos  << " ---- RotationMatrix: " << rm <<  " ---- Angles:" << rm.getTheta() << " " << rm.getPhi() << " " << rm.getPsi() << G4endl;
                return vectPos;
            }

            for (G4int cc =0 ; cc < Phy2->GetLogicalVolume()->GetNoDaughters() ; cc++) {

                Phy3 = Phy2->GetLogicalVolume()->GetDaughter(cc);
                vectPos3 = Phy3->GetObjectTranslation();
                rm3 = Phy3->GetObjectRotationValue();

                if(Phy3->GetLogicalVolume()->GetName() == RegionName){

                    vectPos = vectPos3.transform(Phy2->GetObjectRotationValue()) + vectPos2.transform(Phy1->GetObjectRotationValue()) + vectPos1;
                    rm.transform(rm1); rm.transform(rm2); rm.transform(rm3);

                    OrganNameRotMatrixMap[RegionName] = rm;
                    //G4cout << " ---- " << RegionName << " Pos:" << vectPos  << " ---- RotationMatrix: " << rm <<  " ---- Angles:" << rm.getTheta() << " " << rm.getPhi() << " " << rm.getPsi() << G4endl;
                    return vectPos;
                }

                for (G4int cc =0 ; cc < Phy3->GetLogicalVolume()->GetNoDaughters() ; cc++) {

                    Phy4 = Phy3->GetLogicalVolume()->GetDaughter(cc);
                    vectPos4 = Phy4->GetObjectTranslation();
                    rm4 = Phy4->GetObjectRotationValue();

                    if(Phy4->GetLogicalVolume()->GetName() == RegionName){

                        vectPos = vectPos4.transform(Phy3->GetObjectRotationValue()) + vectPos3.transform(Phy2->GetObjectRotationValue()) + vectPos2.transform(Phy1->GetObjectRotationValue()) + vectPos1;
                        rm.transform(rm1); rm.transform(rm2); rm.transform(rm3);  rm.transform(rm4);

                        OrganNameRotMatrixMap[RegionName] = rm;
                        //G4cout << " ---- " << RegionName << " Pos:" << vectPos  << " ---- RotationMatrix: " << rm <<  " ---- Angles:" << rm.getTheta() << " " << rm.getPhi() << " " << rm.getPsi() << G4endl;
                        return vectPos;
                    }
                }
            }
        }
    }
    return vectPos;
}
void G4TVolumeConstruction::ShowBoxVolume(){

    if(SourceType == "Voxels"){
        return;
    }

    boxCenterPos = getRegionAbsolutePosition(SourceRegionName);

    //G4cout << "\n" << boxCenterPos << " ---- " << BoxDimGene.getX()  << " " << BoxDimGene.getY() << " " << BoxDimGene.getZ() <<G4endl;

    if(ShowBox == "yes"){

        G4LogicalVolume*    LogicalBoxVolume;

        //  G4cout << "\n" << boxCenterPos << " ---- " << getRegionAbsolutePosition(SourceRegionName) << " ---- " << BoxDimGene.getX()  << " " << BoxDimGene.getY() << " " << BoxDimGene.getZ() <<G4endl;

        G4Box* box = new G4Box("Box", BoxDimGene.getX(), BoxDimGene.getY(), BoxDimGene.getZ());
        G4double Z, A;
        A = 14.01*g/mole; G4Element* elN = new G4Element("Nitrogen","N",Z = 7.,A);
        A = 16.00*g/mole; G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A);

        G4double d = 1.290*mg/cm3;
        G4Material* AirMaterial = new G4Material("Air",d,2);
        AirMaterial->AddElement(elN,0.7); AirMaterial->AddElement(elO,0.3);

        LogicalBoxVolume = new G4LogicalVolume(box , AirMaterial , BoxLVName , 0 , 0 , 0);

        G4VisAttributes* boxVisAtt = new G4VisAttributes(G4Colour::Yellow());
        boxVisAtt->SetForceSolid(true);
        LogicalBoxVolume->SetVisAttributes(boxVisAtt);

        G4RotationMatrix* BoxMatrixRot = new G4RotationMatrix();
        BoxMatrixRot->rotateX((G4ThreeVector(0. , 0. , 0.)).x()); BoxMatrixRot->rotateY((G4ThreeVector(0. , 0. , 0.)).y()); BoxMatrixRot->rotateZ((G4ThreeVector(0. , 0. , 0.)).z());
        //PhysicalBoxVolume = new G4PVPlacement( &OrganNameRotMatrixMap[SourceRegionName] , boxCenterPos , BoxLVName , LogicalBoxVolume, WorldPhysicalVolume , false , 0 , true );

        // ///////////////////////////////////////////////////////////////


        G4RotationMatrix* RegMatrixRot = new G4RotationMatrix();
        RegMatrixRot->rotateX((G4ThreeVector(0. , 0. , 0.)).x()); BoxMatrixRot->rotateY((G4ThreeVector(0. , 0. , 0.)).y()); BoxMatrixRot->rotateZ((G4ThreeVector(0. , 0. , 0.)).z());

        G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance() ;

        G4bool IsFound = false;
        G4String RegName;
        G4String PVName= "";
        G4VPhysicalVolume* PV;
        G4LogicalVolume* RegLV;

        for (size_t iLV = 0; iLV < pvs->size(); iLV++ ) {

            PV = (*pvs)[iLV];
            if(PV->GetLogicalVolume()->GetName() == SourceRegionName){
                PVName = PV->GetName();
                RegLV = PV->GetLogicalVolume();
                RegName = PV->GetLogicalVolume()->GetName();
                RegMatrixRot = PV->GetObjectRotation();
                IsFound = true;
                //G4cout << " Relative data ---- " << SourceRegionName << " Pos:" << PV->GetObjectTranslation() << " ---- RotationMatrix: " << PV->GetObjectRotationValue() << G4endl;
                //G4cout << " Absolute data ---- " << SourceRegionName << " Pos:" << OrganNamePositionMap[SourceRegionName] << " ---- RotationMatrix: " << OrganNameRotMatrixMap[SourceRegionName] << G4endl;

                delete PV;
                break;
            }
        }

        if(IsFound){
            G4cout << " \n\nRegion found to be surrended by BoxVolume\n\n" << G4endl;

            //WorldPhysicalVolume->GetLogicalVolume()->RemoveDaughter(PV);
            PhysicalBoxVolume = new G4PVPlacement( new G4RotationMatrix() , OrganNamePositionMap[SourceRegionName] , BoxLVName , LogicalBoxVolume, WorldPhysicalVolume , false , 0 , false );
            new G4PVPlacement( &OrganNameRotMatrixMap[SourceRegionName] , G4ThreeVector() , RegName , RegLV, PhysicalBoxVolume , false , 0 , false );
        }else{
            G4cout << " \n\nCanno't find the region to be surrended by BoxVolume\n\n" << G4endl;
        }
    }
}

void G4TVolumeConstruction::CreateMaterialsGDMLTags(){

    G4MaterialTable Mates = *G4Material::GetMaterialTable();
    G4ElementTable Elems = *G4Element::GetElementTable();
    const G4IsotopeTable Isots= *G4Isotope::GetIsotopeTable();

    G4String IN, EN, MN ;
    G4bool IsAbun;
    G4double IsoRelAbun, MassFra;
    G4int nAtom;

    /*
    G4Isotope* ISO = new G4Isotope( IN , MaterialDensity, MaterialCompNumber); // if we dont use g/cm3, we will find that density equal to the same as it added by user

    G4Element* ELE = new G4Element( EN , MaterialDensity, MaterialCompNumber); // if we dont use g/cm3, we will find that density equal to the same as it added by user
    for (G4int n = 0; n < Isots.size(); n++ ){
        if(ISO[n].GetName() != IN){
            ELE->AddIsotope(&ISO[n], IsoRelAbun);
        }
    }

    G4Material* MAT = new G4Material( MN , MaterialDensity, MaterialCompNumber); // if we dont use g/cm3, we will find that density equal to the same as it added by user
    for (G4int n = 0; n < Elems.size(); n++ ){
        if(Elems[n]->GetName() != EN){
            if(IsAbun){
                MAT->AddElementByMassFraction(Elems[n], MassFra);
            }
            else{
                MAT->AddElementByNumberOfAtoms(Elems[n], nAtom);
            }
        }
    }
    */
    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::ostringstream texttt;

    for (G4int n = 0; n < Elems.size(); n++ ){
        texttt << "  <element Z=\""
               << Elems[n]->GetZ() << "\" formula=\""
               << Elems[n]->GetName() << "\" name=\""
               << Elems[n]->GetSymbol() << "\" > <atom value=\"" << Elems[n]->GetA() << "\" unit=\"g/mole\" /> </element>\n";
    }

    for (G4int n = 0; n < Mates.size(); n++ ){
        texttt << "\n  <material  name=\"" << Mates[n]->GetName()
               << "\"  formula=\"" << Mates[n]->GetName()+std::to_string(ElemInc)
               << "\">  <D value=\"" << Mates[n]->GetDensity()*G4Density_to_gPerCm3 << "\"  unit=\"g/cm3\" /> \n" ;
#if VERBOSE_USE
        G4cout << "\nMates[n]->GetNumberOfElements() " << Mates[n]->GetNumberOfElements() << G4endl ;
#endif
        for (G4int m = 0; m < Mates[n]->GetNumberOfElements(); m++ ){
            if( Mates[n]->GetFractionVector()){
                texttt << "    <fraction  n=\"" << Mates[n]->GetElement(m)->GetAtomicMassAmu()*perCent << "\" ref=\"" << Mates[n]->GetElement(m)->GetName() << "\" /> \n" ;
            }else {
                texttt << "    <composite n=\"" << Mates[n]->GetElement(m)->GetAtomicMassAmu()*perCent << "\" ref=\"" << Mates[n]->GetElement(m)->GetName() << "\" /> \n" ;
            }
        }
        texttt << "  </material> \n" ;
    }

#if VERBOSE_USE
    G4cout << "\nCreating file " << texttt.str().c_str() << G4endl ;
#endif

}

// called after construct method, after the geometry is created, to set the SD to region of interest
void G4TVolumeConstruction::ConstructSDandField()
{

    if(GeometryFileType == "TET"){

        //std::cout  << "\n\n========= Set The SD to the Geometry Volume..." << std::endl;

        // Define detector (Phantom SD) and scorer (eDep)
        //
        G4String phantomSDname = "PhantomSD";

        // MultiFunctional detector
        G4MultiFunctionalDetector* MFDet = new G4MultiFunctionalDetector(phantomSDname);

        G4SDManager::GetSDMpointer()->AddNewDetector( MFDet );

        // scorer for energy depositon in each organ
        MFDet->RegisterPrimitive(new TETPSEnergyDeposit("eDep", tetData));

        // attach the detector to logical volume for parameterised geometry (phantom geometry)
        SetSensitiveDetector(tetLogic, MFDet);
    }

    //G4cout  << "************** Set The SD to the Geometry Volume..." << G4endl;
    /*
    G4TSD* SD = new G4TSD("SD", "TCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(SD);

    if(GeometryFileType == "ICRPFiles" || GeometryFileType == "VOXEL" || GeometryFileType == "DICOM" || GeometryFileType == "VoxIDs"){ //
        SetSensitiveDetector("VoxelLogical", SD);
        if(WorldAsSD){
            //G4cout  << OrganNamesVector[lk] <<" yes " << G4endl ;
            SetSensitiveDetector(WorldVolumeName,SD);
        }
    }

    else { // for TEXT GDML Construct

        for(G4int lk = 0 ; lk < (G4int)OrganNamesVector.size() ; lk++ ){ // the WorldVolume
            if(OrganNamesVector[lk] == "World"){
                if(WorldAsSD){
                    //G4cout  << OrganNamesVector[lk] <<" yes " << G4endl ;
                    SetSensitiveDetector(OrganNamesVector[lk],SD);
                }
                else {
                    //G4cout  << OrganNamesVector[lk] <<" no " << G4endl ;
                    continue;
                }
            }else{
                SetSensitiveDetector(OrganNamesVector[lk],SD);
                //G4cout  << OrganNamesVector[lk] <<" yes " << G4endl ;
            }
        }

        BiasFlag = "no";
        if(BiasFlag=="yes"){

            RegionWhereToBias = "World";
            // -- Fetch volume for biasing:
            G4TGeometryBasedBiasingOperator* GeometryBasedBiasingOperator =  new G4TGeometryBasedBiasingOperator();
            G4LogicalVolume* logicTest;
            for ( G4int nj = 0; nj < OrganNamesVector.size(); nj++  )
            {
                if(OrganNamesVector[nj] == "World"){
                    continue;
                }
                if(OrganNamesVector[nj] == SourceRegionName){
                    continue;
                }
                logicTest = G4LogicalVolumeStore::GetInstance()->GetVolume(OrganNamesVector[nj]);
                GeometryBasedBiasingOperator->AttachTo(logicTest);
                G4cout << " Attaching biasing operator " << GeometryBasedBiasingOperator->GetName() << " to logical volume " << logicTest->GetName() << G4endl;
                //}
            }
        }

    }
   */
}

// this to add Special tracking cuts to a volume (not like production cut(special is advanced))
// it's activated as a process G4StepLimiter to a particle in physics list class or as a physics in main using class new G4StepLimiterPhysics() .
void G4TVolumeConstruction::setUserLimits(){

    VolumeOfLimits = "Liver";

    max_allowed_step_size = 0.5*mm; // Limitation to step, to activate it, G4StepLimiter process must be defined in physics to affected particle types.
    max_total_track_length = DBL_MAX; //5*mm; //Limitations to track, to activate it, G4UserSpecialCuts process must be defined in physics to affected particle types.
    max_total_time_of_flight = 10*ms; //Limitations to track, to activate it, G4UserSpecialCuts process must be defined in physics to affected particle types.
    min_kinetic_energy = 1*keV; //Limitations to track, to activate it, G4UserSpecialCuts process must be defined in physics to affected particle types.
    min_remaining_range = 0.1*mm; //Limitations to track, to activate it, G4UserSpecialCuts process must be defined in physics to affected particle types.

    G4UserLimits* stepLimits = new G4UserLimits(max_allowed_step_size);

    G4LogicalVolume* logicTest;
    for ( G4int nj = 0; nj < OrganNamesVector.size(); nj++  )
    {
        if(OrganNamesVector[nj] == VolumeOfLimits){
            logicTest = G4LogicalVolumeStore::GetInstance()->GetVolume(OrganNamesVector[nj]);
            logicTest->SetUserLimits(stepLimits);
        }
    }

}
G4VIStore* G4TVolumeConstruction::CreateImportanceStore()
{
#if VERBOSE_USE
    G4cout << " ################################################# " << __FUNCTION__ <<G4endl;

    G4cout << " B01DetectorConstruction:: Creating Importance Store " << G4endl;
#endif
    //if (!fPhysicalVolumeVector.size())
    //{
    // G4Exception("B01DetectorConstruction::CreateImportanceStore" ,"exampleB01_0001",RunMustBeAborted ,"no physical volumes created yet!");
    //}

    //WorldPhysicalVolume = fPhysicalVolumeVector[0];

    // creating and filling the importance store

    G4IStore *istore = G4IStore::GetInstance();

    std::vector<G4String> VolumesWithImportance;
    VolumesWithImportance.push_back("Thyroid");
    VolumesWithImportance.push_back("Scapulae");
    VolumesWithImportance.push_back("Pnacreas");

    G4int n = 1;
    G4double imp = 2;
    istore->AddImportanceGeometryCell(1,  *WorldPhysicalVolume);
    G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance() ;
    for (size_t iLV = 0; iLV < pvs->size(); iLV++ ) {
        G4VPhysicalVolume* PV = (*pvs)[iLV];

        if (PV != WorldPhysicalVolume)
        {
        }
        else{

            for ( G4int nj = 0; nj < VolumesWithImportance.size(); nj++  )
            {
                if(PV->GetLogicalVolume()->GetName() == VolumesWithImportance[nj]){
                    imp = std::pow(2., n);
#if VERBOSE_USE
                    G4cout << "Going to assign importance: " << imp << ", to volume: " << PV->GetName() << G4endl;
#endif
                    istore->AddImportanceGeometryCell(imp, *(PV),n); // /////////////////////////////// n !!!!!!!!!!!!!!
                }else {
                    istore->AddImportanceGeometryCell(1, *(PV));
                }
            }
        }
    }

    // the remaining part pf the geometry (rest) gets the same
    // importance as the last conrete cell
    //
    //istore->AddImportanceGeometryCell(imp, *(fPhysicalVolumeVector[fPhysicalVolumeVector.size()-1]),++n);

    return istore;
}

// called after setting commands to create the data files needed by data generation class
G4String G4TVolumeConstruction::setSourcePositionFileName(){

    //G4cout << " ################################################# NumberOfGenPointsToSave : " << NumberOfGenPointsToSave <<G4endl;

    std::ostringstream a ;
    if(SourceType == "Solid"){
        a << DataDirectoryPath  << "/Pos_" << SourceRegionName << "_" << SourceSolid << "_" << SourceType << "_" << NumberOfGenPointsToSave << "_0"<< DataFilesExtension;
    }
    else if(SourceType == "Surface"){
        a << DataDirectoryPath  << "/Pos_" << SourceRegionName << "_" << SourceSurface << "_" << SourceType << "_" << NumberOfGenPointsToSave << "_0"<< DataFilesExtension;
    }
    else if(SourceType == "Plane"){
        a << DataDirectoryPath  << "/Pos_" << SourceRegionName << "_" << SourcePlane << "_" << SourceType << "_" << NumberOfGenPointsToSave << "_0"<< DataFilesExtension;
    }
    else {
        a << DataDirectoryPath  << "/Pos_" << SourceRegionName << "_" << SourceType << "_" << NumberOfGenPointsToSave << "_0"<< DataFilesExtension;
    }

    PositionDataFile = a.str().c_str();
    return PositionDataFile;
}
G4String G4TVolumeConstruction::setSourceEnergyFileName(){
    std::ostringstream a ;
    a << DataDirectoryPath  << "/Ene_" << EnergyDistribution <<"_";
    if(EnergyDistribution == "Mono"){
        a << MonoEnergy;
    }
    else if(EnergyDistribution == "Rayleigh"){
        a << RayleighEmax;
    }
    else if(EnergyDistribution == "Gauss"){
        a << GaussMean;
    }
    else if(EnergyDistribution == "Uniform"){
        a << UniformEmax;
    }
    else if(EnergyDistribution == "Spectrum"){
        a << SpectrumMaxEnergy;
    }
    else if(EnergyDistribution == "File"){
        a << FileEnergyCharacterizer;
    }
    else if(EnergyDistribution == "RadioNuclide"){
        a << RadioNuclideMaxEnergy;
    }
    a << "_" << NumberOfGenPointsToSave << "_0"<< DataFilesExtension;
    EnergyDataFile = a.str().c_str();
    return EnergyDataFile;
}
G4String G4TVolumeConstruction::setSourceMomDirFileName(){
    std::ostringstream a ;
    a << DataDirectoryPath  << "/MomDir_" << MomDirDistribution ;
    if(MomDirDistribution == "Isotropic"){
        //a << "";
    }
    else if(MomDirDistribution == "Uniform"){
        //a << "";
    }
    else if(MomDirDistribution == "Directed"){
        a << Theta << "_" << Phi ;
    }
    a <<  "_" << NumberOfGenPointsToSave << "_0"<< DataFilesExtension;
    MomDirDataFile = a.str().c_str();
    return MomDirDataFile;

}

//############################################################# Creating Materials

void G4TVolumeConstruction::createElement(G4String n, G4double Z, G4double A){
    //G4cout << __FUNCTION__<<G4endl;
    CreatedElements[n] = new G4Element( n , "ELE", Z , A );

    ElementsGDMLText << "  <element Z=\""
                     << Z << "\" formula=\""
                     << n << "\" name=\""
                     << n << "\" > <atom value=\"" << A << "\" unit=\"g/mole\" /> </element>\n";
}
void G4TVolumeConstruction::createNistMaterial(G4String n, G4int i){
    //G4cout << __FUNCTION__<<G4endl;
    G4NistManager* nist = G4NistManager::Instance();
    CreatedMaterials[n]=nist->FindOrBuildMaterial(n,true);
    MaterialIDName[i] = n ;

    //MaterialsGDMLtext << "\n  <material  name=\"" << n << "\"> </material> \n" ;

}
void G4TVolumeConstruction::createMaterialFromComponents(G4String n, G4int i, G4int nu, G4double d, G4String f){
    //G4cout << __FUNCTION__<<G4endl;

    MaterialName = n;
    MaterialIDName[i] = n ;
    FracOrNum = f;
    MaterialCompNumber = nu;

    // the density should be setted with with unit to recognized internally by geant4
    CurrentMaterial = new G4Material( n , d, nu); // if we dont use g/cm3, we will find that density equal to the same as it added by user
    //G4cout << " Material->GetDensity "<< CurrentMaterial->GetDensity()*G4Density_to_gPerCm3 << " Elem Num " << nu <<G4endl;

    ElemInc = 0;

    MaterialsGDMLtext << "\n  <material  name=\"" << n
                      << "\"  formula=\"" << n+std::to_string(ElemInc)
                      << "\">  <D value=\"" << CurrentMaterial->GetDensity()*G4Density_to_gPerCm3 << "\"  unit=\"g/cm3\" /> \n" ;
}
void G4TVolumeConstruction::AddMatElement(G4String s, G4double d, G4int n){

    //G4cout << __FUNCTION__<<G4endl;

    for ( auto it = CreatedElements.begin(); it != CreatedElements.end(); ++it  )
    {
        if(s == it->first){
            if( FracOrNum == "frac" ){
                CurrentMaterial->AddElement(CreatedElements[s], d*perCent);
                //G4cout << s << " " << d << " - percent: " << d*perCent << G4endl;

                MaterialsGDMLtext << "    <fraction n=\"" << d*perCent << "\" ref=\"" << s << "\" /> \n" ;

            }else {
                CurrentMaterial->AddElement(CreatedElements[s], n);
                MaterialsGDMLtext << "    <composite n=\"" << n<< "\" ref=\"" << s << "\" /> \n" ;
            }
        }
    }
    for ( auto it = CreatedMaterials.begin(); it != CreatedMaterials.end(); ++it  ){
        if(s == it->first){
            if( FracOrNum == "frac" ){
                CurrentMaterial->AddMaterial(CreatedMaterials[s], d*perCent);
                MaterialsGDMLtext << "    <fraction n=\"" << d*perCent << "\" ref=\"" << s << "\" /> \n" ;

            }else {
                CurrentMaterial->AddMaterial(CreatedMaterials[s],n);
                MaterialsGDMLtext << "    <composite n=\"" <<n<< "\" ref=\"" << s << "\" /> \n" ;
            }
        }
    }

    ElemInc++;
    //std::cout << ElemInc << " " << MaterialCompNumber << std::endl;

    if( MaterialCompNumber == ElemInc){ MaterialsGDMLtext << "  </material> \n" ; }

    CreatedMaterials[MaterialName]=CurrentMaterial;
}

//############################################################ World Constructing

void G4TVolumeConstruction::constructWorldVolume(G4String mn, G4ThreeVector dim){

    WorldMaterialName = mn;
    WorldHalfSize = dim;
    //WorldVolumeName = "World";
    //G4cout<<"\nConstruct World Volume "<< WorldVolumeName <<" \n" <<G4endl;

    G4Box* world = new G4Box("World", WorldHalfSize.getX()/2., WorldHalfSize.getY()/2., WorldHalfSize.getZ()/2.);
    G4LogicalVolume* logicMotherWorld = new G4LogicalVolume(world, CreatedMaterials[WorldMaterialName], "World", 0, 0,0);
    WorldPhysicalVolume = new G4PVPlacement(0,G4ThreeVector(), "World", logicMotherWorld, 0 , false, 0 , false );

/*
    G4UserLimits(G4double uStepMax = DBL_MAX,
                 G4double uTrakMax = DBL_MAX,
                 G4double uTimeMax = DBL_MAX,
                 G4double uEkinMin = 0.,
                 G4double uRangMin = 0. );
*/
    //logicMotherWorld->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX,1*MeV,DBL_MAX));

    //OrganNamesVector.push_back(WorldVolumeName);

    //OrganNameMassMap[WorldVolumeName] = WorldPhysicalVolume->GetLogicalVolume()->GetMass()*G4Mass_to_Kg; // density in g/cm3;
    //OrganNameDensityMap[WorldVolumeName] = WorldPhysicalVolume->GetLogicalVolume()->GetMaterial()->GetDensity()*G4Density_to_gPerCm3; // density in g/cm3;
    //CreatedPositionOrgans[WorldVolumeName] = G4ThreeVector();
}

//############################################################ For Construct Geometry

void G4TVolumeConstruction::setSTLVolumeMaterial(G4String n){STLVolumeMaterial = n ;}
// place a volume (create a Physical Volume)
void G4TVolumeConstruction::placeVolume(G4String n, G4String MotherVolumeName, G4ThreeVector Position, G4ThreeVector Rotation){

    G4String PhysicalVolumeName = n ;
    G4String FileName = getFileNameFromPath(n) ;
    G4String FileExtension = getFileExt(n);;

//    G4cout << "\n\n --------------- " << __FUNCTION__ << G4endl;

    //G4cout << __FUNCTION__ << G4endl;
    //G4cout << PhysicalVolumeName<<G4endl;
    //G4cout << __FUNCTION__ << " " << PhysicalVolumeName << " " << MotherVolumeName << " " << Position << " " << Rotation << G4endl;

    if(WorldPhysicalVolume == NULL){
        if(FileExtension == "gdml" || FileExtension == "geom" || FileExtension == "c++" || FileExtension == "cpp" || n == "MyGeometry"){
            FileName = "World";
        }else{
            // create world volume first
        }
    }

    //G4cout << "\n\n --------------- " << FileName << G4endl;

    if(FileName == "World"){

        if(n == "MyGeometry"){

            G4TMyGeometry * CGF = new G4TMyGeometry();
            WorldPhysicalVolume = CGF->ConstructPhysicalVolume();
        }else{
            if(FileExtension == "gdml"){
#if GDML_USE
                G4GDMLParser parser;
                parser.Read( PhysicalVolumeName , false );  //false to eliminate the xchema validation because it print a lot of lines of error validation
                WorldPhysicalVolume = parser.GetWorldVolume();
#else
            G4String msg =  "DoseCalcs should be built with \"-DWITH_GDML_USE=ON\" to construct geometry from .gdml files "; G4Exception("Geometry Data", "1", FatalErrorInArgument, msg.c_str());
#endif
            }
            else if(FileExtension == "geom"){
                G4tgbVolumeMgr* G4tgbVolumeMgrObj;
                G4tgbVolumeMgrObj = G4tgbVolumeMgr::GetInstance();
                G4tgbVolumeMgrObj->AddTextFile(PhysicalVolumeName);
                WorldPhysicalVolume = G4tgbVolumeMgrObj->ReadAndConstructDetector();
            }
            else if(FileExtension == "cpp" || FileExtension == "c++"){

                G4TCPPGeometryFormat * CGF = new G4TCPPGeometryFormat();
                WorldPhysicalVolume = CGF->ConstructPhysicalVolume();
            }
        }
        WorldPhysicalVolume->SetName("World");
        WorldPhysicalVolume->GetLogicalVolume()->SetName("World");
#if VERBOSE_USE
        //G4cout << "\n\n - Volume " << WorldPhysicalVolume->GetLogicalVolume()->GetName() << " - mother : " << WorldPhysicalVolume->GetName() << G4endl;
#endif
        //G4cout << __FUNCTION__ << " --------------- "  << G4endl;

        //G4cout << "World Logical Name " << WorldPhysicalVolume->GetLogicalVolume()->GetName() << " World Physical Name " << WorldPhysicalVolume->GetName()<< G4endl;

        return;
    }

    if(FileExtension == "gdml"){

        PhysicalVolumeName = FileName;
        std::ifstream ifile(n.c_str());
        if (ifile) {
#if GDML_USE
            G4GDMLParser parser;
            parser.Read(n, false );  //false to eliminate the xchema validation because it print a lot of lines of error validation
            CurrentLogicalVol = parser.GetVolume(PhysicalVolumeName);
#else
            G4String msg =  "DoseCalcs should be built with \"-DWITH_GDML_USE=ON\" to construct geometry from .gdml files "; G4Exception("Geometry Data", "1", FatalErrorInArgument, msg.c_str());
#endif
        }
        else{
            //G4cout << " Unable to get logical volume : " << CurrentLogicalVol->GetName() << G4endl;
        }
    }
    else if(FileExtension == "geom"){

        PhysicalVolumeName = FileName;

        std::ifstream ifile(n.c_str());
        if (ifile) {

            G4tgbVolumeMgr* G4tgbVolumeMgrObj;
            G4tgbVolumeMgrObj = G4tgbVolumeMgr::GetInstance();

            G4tgbVolumeMgrObj->AddTextFile(n);
            //G4tgbVolumeMgrObj->RegisterMe(G4tgbVolumeMgrObj->FindVolume(PhysicalVolumeName));
            //G4LogicalVolume * logicChamber = G4tgbVolumeMgr::GetInstance()->FindG4LogVol(PhysicalVolumeName, 0);

            G4tgbVolumeMgrObj->ReadAndConstructDetector()->GetLogicalVolume();
            delete G4PhysicalVolumeStore::GetInstance()->GetVolume(PhysicalVolumeName); // because when reading gdml by parser an physical volume will be registered
        }
        else{
            //G4cout << " Unable to get logical volume : " << CurrentLogicalVol->GetName() << G4endl;
        }
    }
    else if(FileExtension == "cpp" || FileExtension == "c++"){
        if(CPPLogVolAreBuilt == false){
            G4TCPPGeometryFormat * CGF = new G4TCPPGeometryFormat();
            CGF->ConstructLogicalVolumes();
            CPPLogVolAreBuilt == true;
        }
        PhysicalVolumeName = FileName;
    }
    else if(FileExtension == "stl" || FileExtension == "ast"){

#if GDML_USE
        PhysicalVolumeName = FileName;

        createStlVol = new G4TStlToGdml();
        createStlVol->CreateMaterialGDMLFile();

        createStlVol->setStlMaterialName(STLVolumeMaterial);
        createStlVol->setStructureVolumeName(PhysicalVolumeName);
        createStlVol->stl_to_gdml(PhysicalVolumeName, n);
        //CurrentLogicalVol = createStlVol->GetLogicalVolFromGDMLFile();
        createStlVol->GetLogicalVolFromGDMLFile();

        //G4VSolid::BoundingLimits()
        // because when reading gdml by parser an physical volume will be registered
        // G4PhysicalVolumeStore::GetInstance()->GetVolume(PhysicalVolumeName)->Clean();
#else
            G4String msg =  "DoseCalcs should be built with \"-DWITH_GDML_USE=ON\" to construct geometry from .gdml files "; G4Exception("Geometry Data", "1", FatalErrorInArgument, msg.c_str());
#endif
    }
    else{
        PhysicalVolumeName = FileName;
    }

    // to make the name pf physical volume the same as logical volume
    G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance() ;
    for (size_t iLV = 0; iLV < pvs->size(); iLV++ ) {if((*pvs)[iLV]->GetLogicalVolume()->GetName() != (*pvs)[iLV]->GetName()){(*pvs)[iLV]->SetName((*pvs)[iLV]->GetLogicalVolume()->GetName());}}

    G4RotationMatrix* rm = new G4RotationMatrix(); rm->rotateX(Rotation.x()); rm->rotateY(Rotation.y()); rm->rotateZ(Rotation.z());
    //CurrentPhysicalVol = new G4PVPlacement( rm , Position , PhysicalVolumeName , ls->GetVolume(PhysicalVolumeName) , ps->GetVolume(MotherVolumeName) , false , 0 , false );

    G4VPhysicalVolume* MotherVol = G4PhysicalVolumeStore::GetInstance()->GetVolume(MotherVolumeName);
    new G4PVPlacement( rm , Position , PhysicalVolumeName , G4LogicalVolumeStore::GetInstance()->GetVolume(PhysicalVolumeName) , MotherVol , false , 0 , false );

    //G4cout << "\n\n - Volume " << PhysicalVolumeName << " - mother : " << MotherVolumeName << "=" << G4PhysicalVolumeStore::GetInstance()->GetVolume(PhysicalVolumeName)->GetMotherLogical()->GetName() << " - Pos "<< Position << " - Rot " << Rotation << " - rm "<< *G4PhysicalVolumeStore::GetInstance()->GetVolume()->GetRotation() <<G4endl;

}
void G4TVolumeConstruction::PlaceVolumeFromRegionfile(G4String n, G4String MotherVolumeName){

    G4ThreeVector Position, Rotation;

    G4String FileName = getFileNameFromPath(n) ;
    G4String FileExtension = getFileExt(n);;

    //G4cout << __FUNCTION__ << " " << FileName << " " << MotherVolumeName << " " << Position << " " << Rotation << G4endl;

    if(WorldPhysicalVolume == NULL){

        if(FileExtension == "gdml" || FileExtension == "geom" || FileExtension == "c++" || FileExtension == "cpp" ){
            FileName = "World";
            if(FileExtension == "gdml" ){
#if GDML_USE
#else
                G4String msg =  "DoseCalcs should be built with \"-DWITH_GDML_USE=ON\" to construct geometry from .gdml files "; G4Exception("Geometry Data", "1", FatalErrorInArgument, msg.c_str());
#endif
            }
        }else{
            // create world volume first
        }
    }

    std::ifstream ifile(n.c_str());
    if (ifile) {

        if(FileExtension == "gdml"){

#if GDML_USE
#else
            G4String msg =  "DoseCalcs should be built with \"-DWITH_GDML_USE=ON\" to construct geometry from .gdml files "; G4Exception("Geometry Data", "1", FatalErrorInArgument, msg.c_str());
#endif

#if GDML_USE
            G4GDMLParser parser;
            parser.Read(n, false );  //false to eliminate the xchema validation because it print a lot of lines of error validation

            CurrentLogicalVol = parser.GetVolume(FileName);

            G4String sub_str = "Pos";
            G4String nm1 = FileName;
            //G4cout << "\n\n --------------- " << nm1 <<  " " << Poss << G4endl;
            sub_str = "Rot";
            G4String nm2 = FileName;
            //G4cout << "\n\n --------------- " << nm2 <<  " " << Rott << G4endl;

            //G4String Poss = nm1.replace(nm1.find("Volume"), 6, sub_str); //remove and replace from that position
            //G4String Rott = nm2.replace(nm2.find("Volume"), 6, sub_str); //remove and replace from that position
            G4String Poss = "OrganPos"; // for ORNL
            G4String Rott = "OrganRot"; // for ORNL
            Position = parser.GetPosition(Poss);
            Rotation = parser.GetRotation(Rott);

            G4PhysicalVolumeStore::DeRegister(parser.GetWorldVolume());
#endif
        }
        else if(FileExtension == "geom"){

            G4tgbVolumeMgr* G4tgbVolumeMgrObj;
            G4tgbVolumeMgrObj = G4tgbVolumeMgr::GetInstance();

            G4tgbVolumeMgrObj->AddTextFile(n);
            //G4tgbVolumeMgrObj->RegisterMe(G4tgbVolumeMgrObj->FindVolume(PhysicalVolumeName));
            //G4LogicalVolume * logicChamber = G4tgbVolumeMgr::GetInstance()->FindG4LogVol(PhysicalVolumeName, 0);

            CurrentLogicalVol = G4tgbVolumeMgrObj->ReadAndConstructDetector()->GetLogicalVolume();
            delete G4PhysicalVolumeStore::GetInstance()->GetVolume(FileName); // because when reading gdml by parser an physical volume will be registered

            G4String sub_str = "Pos";
            G4String nm1 = FileName;
            //G4cout << "\n\n --------------- " << nm1 <<  " " << Poss << G4endl;
            sub_str = "Rot";
            G4String nm2 = FileName;
            //G4cout << "\n\n --------------- " << nm2 <<  " " << Rott << G4endl;

            //G4String Poss = nm1.replace(nm1.find("Volume"), 6, sub_str); //remove and replace from that position
            //G4String Rott = nm2.replace(nm2.find("Volume"), 6, sub_str); //remove and replace from that position
            G4String Poss = "OrganPos"; // for ORNL
            G4String Rott = "OrganRot"; // for ORNL
            //Position = G4tgbVolumeMgrObj->GetInstance()-> GetPosition(Poss);
            //Rotation = parser.GetRotation(Rott);

        }
    }
    else{
        std::cout << "\nCanno't read the file " << n <<  std::endl;
        return;
    }

    G4RotationMatrix* rm = new G4RotationMatrix(); rm->rotateX(Rotation.x()); rm->rotateY(Rotation.y()); rm->rotateZ(Rotation.z());
    //CurrentPhysicalVol = new G4PVPlacement( rm , Position , FileName , ls->GetVolume(FileName) , ps->GetVolume(MotherVolumeName) , false , 0 , false );

    G4VPhysicalVolume* MotherVol = G4PhysicalVolumeStore::GetInstance()->GetVolume(MotherVolumeName);
    new G4PVPlacement( rm , Position , FileName , G4LogicalVolumeStore::GetInstance()->GetVolume(FileName) , MotherVol , false , 0 , false );



    //G4cout << "\n\n - Volume " << FileName << " - mother : " << MotherVolumeName << "=" << G4PhysicalVolumeStore::GetInstance()->GetVolume(FileName)->GetMotherLogical()->GetName() << " - Pos "<< Position << " - Rot " << Rotation << " - rm "<< *G4PhysicalVolumeStore::GetInstance()->GetVolume()->GetRotation() <<G4endl;

}

//########################################################### For Dicom Geometry

void G4TVolumeConstruction::setDcmMaterialName(G4String n){
    /*
    DcmMaterialName = n ;
    //G4cout << "DcmMaterialName : "<< DcmMaterialName <<G4endl;

    MaterialIndices[CreatedMaterials[DcmMaterialName]->GetDensity()*G4Density_to_gPerCm3] = DcmMaterialName;// /(6.24151e+15)
    //G4cout << " density in g/cm3 : " << CreatedMaterials[DcmMaterialName]->GetDensity()*G4Density_to_gPerCm3 <<    " - DcmMaterialName  " << DcmMaterialName <<G4endl;
    DcmMaterials.push_back(CreatedMaterials[DcmMaterialName]);
    */
}
void G4TVolumeConstruction::setDcmCtNumber(G4double n){
    DcmCtNumber = n ;
}
void G4TVolumeConstruction::setDcmCtDensity(G4double n){
    DcmCtDensity = n ;
    ValueCT.push_back(DcmCtNumber);
    ValueDensity.push_back(DcmCtDensity*G4Density_to_gPerCm3);
    //G4cout << "DcmCtNumber " << DcmCtNumber << " - density As it setted : " << DcmCtDensity << G4endl;
    //G4cout << "DcmCtNumber " << DcmCtNumber << " - density in g/cm3     : " << DcmCtDensity*G4Density_to_gPerCm3 << G4endl;
}

void G4TVolumeConstruction::getDicomDataAndConstruct(){

#ifdef DCMTK_USE

    for(G4int V = 0 ; V < DicomFileTypesVec.size() ; V++){

        DicomFileType = DicomFileTypesVec[V];
        G4TVolumeBuilderUsingVoxel* G4TVolumeBuilderUsingDICOMObject;

        if(DicomFileType == "CT"){ // should be readed one time

            DicomOutTextName = DicomFilesDirPath+"/"+DicomCTName;;

            //G4cout << " From function getDicomDataAndConstruct : " << DicomFileType << " : " << DicomOutTextName << G4endl;
            G4TVolumeBuilderUsingDICOMObject = new G4TVolumeBuilderUsingVoxel();
            G4TVolumeBuilderUsingDICOMObject->setDicomType(DicomFileType);
            G4TVolumeBuilderUsingDICOMObject->setDicomOutTextName(DicomOutTextName);
            G4TVolumeBuilderUsingDICOMObject->StartGettingDicomData();

            // geting Voxels data
            VoxXNumber = G4TVolumeBuilderUsingDICOMObject->GetVoxelsDataXYZ()[0];
            VoxYNumber = G4TVolumeBuilderUsingDICOMObject->GetVoxelsDataXYZ()[1];
            VoxZNumber = G4TVolumeBuilderUsingDICOMObject->GetVoxelsDataXYZ()[2];
            VoxXHalfSize = G4TVolumeBuilderUsingDICOMObject->GetVoxelsDataXYZ()[3];
            VoxYHalfSize = G4TVolumeBuilderUsingDICOMObject->GetVoxelsDataXYZ()[4];
            VoxZHalfSize = G4TVolumeBuilderUsingDICOMObject->GetVoxelsDataXYZ()[5];

            MateIDs = G4TVolumeBuilderUsingDICOMObject->GetMateIds();
            VoxelsMaterials = G4TVolumeBuilderUsingDICOMObject->GetMaterials();

            defaultMat = VoxelsMaterials[0]; // used in container Volume, there is no default material name to be set by the user for DICOM gepmetry
            VoxDefaultMaterialName = VoxelsMaterials[0]->GetName();
            MaterialColour.clear();

            for ( G4int h = 0; h < VoxelsMaterials.size(); h++  ){
                RegionNameColour[VoxelsMaterials[h]->GetName()] = G4Colour( (G4double) G4UniformRand(),(G4double)G4UniformRand(),(G4double)G4UniformRand(), 1.);
            }

            //G4TVolumeBuilderUsingDICOM* G4TVolumeBuilderUsingDICOMObject = new G4TVolumeBuilderUsingDICOM();
#if VERBOSE_USE
            G4cout << "\n\n========= Start getting Dicom data from Text file " << DicomOutTextName << " for CT modality======================= \n"<<G4endl;
            G4cout << "\n\n========= Generated DICOM Data ======================= \n"<<G4endl;

            G4cout << std::setw(15) << std::left << " VoxXNumber " << std::setw(15) << std::left << VoxXNumber << std::setw(15) << std::left << "VoxXHalfSize " << std::setw(15) << std::left << VoxXHalfSize <<G4endl;
            G4cout << std::setw(15) << std::left << " VoxYNumber " << std::setw(15) << std::left << VoxYNumber << std::setw(15) << std::left << "VoxYHalfSize " << std::setw(15) << std::left << VoxYHalfSize <<G4endl;
            G4cout << std::setw(15) << std::left << " VoxZNumber " << std::setw(15) << std::left << VoxZNumber << std::setw(15) << std::left << "VoxZHalfSize " << std::setw(15) << std::left << VoxZHalfSize <<G4endl;
            G4cout << std::setw(15) << std::left << " TotalVoxelNum " << std::setw(15) << std::left << VoxXNumber*VoxYNumber*VoxZNumber << "\n" <<  G4endl;

            G4cout << std::setw(15) << std::left << " MatIDs array \n" <<G4endl;
            G4cout << std::setw(15) << std::left << " The new voxels materials data " << G4endl;
#endif

            if(MaterialNameAsRegionName == false){
                CreateRegionsDataFile();
            }
        }
        else  if(DicomFileType == "PET"){

            DicomOutTextName = DicomFilesDirPath+"/"+DicomPETName;

            //G4cout << " From function getDicomDataAndConstruct : " << DicomFileType << " : " << DicomOutTextName << G4endl;

#if VERBOSE_USE
            G4cout << "\n\n========= Start getting " << DicomFileType << " Dicom data from Text file " << DicomOutTextName << " for PET modality ======================= \n"<<G4endl;
#endif
            G4TVolumeBuilderUsingDICOMObject = new G4TVolumeBuilderUsingVoxel();
            G4TVolumeBuilderUsingDICOMObject->setDicomType(DicomFileType);
            G4TVolumeBuilderUsingDICOMObject->setDicomOutTextName(DicomOutTextName);
            G4TVolumeBuilderUsingDICOMObject->StartGettingDicomData();

            VoxPETXNumber = G4TVolumeBuilderUsingDICOMObject->GetPETVoxelsDataXYZ()[0];
            VoxPETYNumber = G4TVolumeBuilderUsingDICOMObject->GetPETVoxelsDataXYZ()[1];
            VoxPETZNumber = G4TVolumeBuilderUsingDICOMObject->GetPETVoxelsDataXYZ()[2];
            VoxPETXHalfSize = G4TVolumeBuilderUsingDICOMObject->GetPETVoxelsDataXYZ()[3];
            VoxPETYHalfSize = G4TVolumeBuilderUsingDICOMObject->GetPETVoxelsDataXYZ()[4];
            VoxPETZHalfSize = G4TVolumeBuilderUsingDICOMObject->GetPETVoxelsDataXYZ()[5];

            //            G4cout << std::setw(15) << std::left << " VoxPETXNumber " << std::setw(15) << std::left << VoxPETXNumber << std::setw(15) << std::left << "VoxPETXHalfSize " << std::setw(15) << std::left << VoxPETXHalfSize <<G4endl;
            //            G4cout << std::setw(15) << std::left << " VoxPETYNumber " << std::setw(15) << std::left << VoxPETYNumber << std::setw(15) << std::left << "VoxPETYHalfSize " << std::setw(15) << std::left << VoxPETYHalfSize <<G4endl;
            //            G4cout << std::setw(15) << std::left << " VoxPETZNumber " << std::setw(15) << std::left << VoxPETZNumber << std::setw(15) << std::left << "VoxPETZHalfSize " << std::setw(15) << std::left << VoxPETZHalfSize <<G4endl;
            //            G4cout << std::setw(15) << std::left << " TotalPETVoxelNum " << std::setw(15) << std::left << VoxPETXNumber*VoxPETYNumber*VoxPETZNumber << "\n" <<  G4endl;

            //            G4cout << std::setw(15) << std::left << " Activity  array \n" <<G4endl;

            Activities = G4TVolumeBuilderUsingDICOMObject->GetActivities();
            ActivitiesVec.push_back(Activities);
            //delete[] Activities;

        }
        delete G4TVolumeBuilderUsingDICOMObject;
    }

    if(DicomFileType == "PET" && UseDicomCumAct){

        // calculate time for each measurement
#if VERBOSE_USE
        G4cout << "\n\n========= Estimating Voxel Cummulated Activity from " << ActivitiesVec.size() << " PET series ======================= \n"<<G4endl;
#endif
        //PETGeneralData["DCM_RadionuclideTotalDose"];
        //PETGeneralData["DCM_RadionuclideHalfLife"];
        //PETGeneralData["DCM_RadionuclidePositronFraction"];
        //PETGeneralData["DCM_RadiopharmaceuticalSpecificActivity"];
        //PETGeneralData["DCM_RadiopharmaceuticalVolume"];

        G4double y1, m1, d1, h1, min1, sec1, firstT, deltaT, TimeInHour, AdministrationT;

        for ( auto it = Serie_Order_DataMap.begin(); it != Serie_Order_DataMap.end(); ++it  ){
#if VERBOSE_USE
            G4cout << " Serie:" << it->first << " - Acquisition Date:" << it->second[0]
                   << " " << it->second[1]
                   << " " << it->second[2]
                   << " " << it->second[3]
                   << " " << it->second[4]
                   << " " << it->second[5]
                   << G4endl;
#endif
        }

        // all in hour unit
        y1   = PETGeneralDateData["DCM_AcquisitionDate"][0]*12*30*24;
        m1   = PETGeneralDateData["DCM_AcquisitionDate"][1]*30*24;
        d1   = PETGeneralDateData["DCM_AcquisitionDate"][2]*24;
        h1   = PETGeneralDateData["DCM_AcquisitionDate"][4];
        min1 = PETGeneralDateData["DCM_AcquisitionDate"][5]/60;
        sec1 = PETGeneralDateData["DCM_AcquisitionDate"][6]/(60*60);
        AdministrationT = y1+m1+d1+h1+min1+sec1; // in hour
        //AdministrationT = AdministrationT/(24);// in day


        //G4cout << " for all series , RadiopharmaceuticalStartTime - RadiopharmaceuticalStopTime = Administration time = " << AdministrationT << " hour "<< G4endl;

#if VERBOSE_USE
        G4cout << "\n\n========= Voxels Cummulated activity Calculation ======================= \n"<<G4endl;
#endif
        //std::vector<G4double> Time = PETGeneralDateData[1] ;// Radiopharmaceutique intake Date
        std::vector<G4double> Time ; Time.push_back(0); Time.push_back(0);Time.push_back(0);Time.push_back(0);Time.push_back(0);Time.push_back(0);

        for ( G4int h = 0; h < ActivitiesVec.size(); h++  ){

            y1   = (Serie_Order_DataMap[h][0]-Time[0])*12*30*24;
            m1   = (Serie_Order_DataMap[h][1]-Time[1])*30*24;
            d1   = (Serie_Order_DataMap[h][2]-Time[2])*24;
            h1   = (Serie_Order_DataMap[h][3]-Time[3]);
            min1 = (Serie_Order_DataMap[h][4]-Time[4])/60;
            sec1 = (Serie_Order_DataMap[h][5]-Time[5])/(60*60);
            TimeInHour = y1+m1+d1+h1+m1+d1;

            deltaT = TimeInHour;

            //<< " For serie "<< Serie_Order_DataMap[h][6] << " of duration(min) "<< Serie_Order_DataMap[h][7]

            Time = Serie_Order_DataMap[h];
            //G4cout << " Serie "<< h << " ActivityTimeVec: " << ActivityTimeVec[h] << G4endl;

            ActivityTimeVec.push_back(SeriesResidenceTime[h]); // just for experience
            //ActivityTimeVec.push_back(deltaT); // just for experience

        }

        // calculate cummulative activity for each voxel

        //ActivityTimeVec.push_back(0);
        //ActivityTimeVec.push_back(200);
        
        size_t VoxelNumberInPhantom = VoxPETXNumber*VoxPETYNumber*VoxPETZNumber;
        CumulativeActivities = new double[VoxelNumberInPhantom];

        //for(size_t hg = 0; hg < VoxelNumberInPhantom ; hg++ ){
        //G4cout << hg << " " << ActivitiesVec[0][hg] << " " << ActivitiesVec[1][hg] << G4endl;
        //}

        G4double TotCummVoxAct = 0;
        for ( G4int h = 0; h < ActivitiesVec.size(); h++  ){
            G4double totalVoxAct = 0;
            G4cout << " Serie "<< h ;

            for(size_t hg = 0; hg < VoxelNumberInPhantom ; hg++ ){
                //G4cout << " for serie "<< h << " voxelID : " << hg << " last CumuAct : " << CumulativeActivities[hg] ;
                totalVoxAct += ActivitiesVec[h][hg];
                CumulativeActivities[hg] += ActivitiesVec[h][hg]*ActivityTimeVec[h];
            }

            TotCummVoxAct+= totalVoxAct;
#if VERBOSE_USE
            G4cout << " - TotalVoxelsAct ("<< totalVoxAct << ") * Residence Time (" <<ActivityTimeVec[h] << " h) = "<< totalVoxAct*ActivityTimeVec[h] << " -TotCummVoxelsAct = " << TotCummVoxAct << G4endl;
#endif
        }
#if VERBOSE_USE
        G4cout << "\n\n========= CT and PET Modality Voxels Data Configuration  ======================= \n" << G4endl;

        G4cout << std::setw(15) << std::left << "| Modality" << std::setw(15) << std::left << "| Compression" << std::setw(15) << std::left << "| VoxXNumber " << std::setw(15) << std::left << "| VoxYNumber"  << std::setw(15) << std::left << "| VoxZNumber"  << std::setw(15) << std::left << "| VoxXHalfSize " << std::setw(15) << std::left << "| VoxYHalfSize " << std::setw(15) << std::left << "| VoxZHalfSize "  <<G4endl;
        G4cout << std::setw(15) << std::left << "| CT" << std::setw(15) << std::left << SeriesCompressionXY[0] << std::setw(15) << std::left << VoxXNumber << std::setw(15) << std::left << VoxYNumber << std::setw(15) << std::left << VoxZNumber  << std::setw(15) << std::left << VoxXHalfSize << std::setw(15) << std::left << VoxYHalfSize << std::setw(15) << std::left << VoxZHalfSize <<G4endl;
        G4cout << std::setw(15) << std::left << "| PET" << std::setw(15) << std::left << SeriesCompressionXY[1] << std::setw(15) << std::left << VoxPETXNumber << std::setw(15) << std::left << VoxPETYNumber << std::setw(15) << std::left << VoxPETZNumber  << std::setw(15) << std::left << VoxPETXHalfSize << std::setw(15) << std::left << VoxPETYHalfSize << std::setw(15) << std::left << VoxPETZHalfSize <<G4endl;
#endif
    }

#else

    G4TVolumeBuilderUsingDICOM* G4TVolumeBuilderUsingDICOMObject = new G4TVolumeBuilderUsingDICOM();

    // geting Voxels data
    VoxXNumber = G4TVolumeBuilderUsingDICOMObject->GetVoxelsDataXYZ()[0];
    VoxYNumber = G4TVolumeBuilderUsingDICOMObject->GetVoxelsDataXYZ()[1];
    VoxZNumber = G4TVolumeBuilderUsingDICOMObject->GetVoxelsDataXYZ()[2];
    VoxXHalfSize = G4TVolumeBuilderUsingDICOMObject->GetVoxelsDataXYZ()[3];
    VoxYHalfSize = G4TVolumeBuilderUsingDICOMObject->GetVoxelsDataXYZ()[4];
    VoxZHalfSize = G4TVolumeBuilderUsingDICOMObject->GetVoxelsDataXYZ()[5];
#if VERBOSE_USE
    G4cout << "\n Generated DICOM Data: \n"<<G4endl;
    G4cout << " VoxXNumber        " << VoxXNumber <<    "  VoxXHalfSize    " << VoxXHalfSize <<G4endl;
    G4cout << " VoxYNumber        " << VoxYNumber <<    "  VoxYHalfSize    " << VoxYHalfSize <<G4endl;
    G4cout << " VoxZNumber    	  " << VoxZNumber <<    "  VoxZHalfSize    " << VoxZHalfSize <<G4endl;
    G4cout << " TotalPixels       " << VoxXNumber*VoxYNumber*VoxZNumber << "\n" <<  G4endl;
#endif
    //MateIDs = G4TVolumeBuilderUsingDICOMObject->GetMateIds();
    VoxelsMaterials = G4TVolumeBuilderUsingDICOMObject->GetMaterials();

#endif

    // generate CN Mass, colour and coordinate data
#if VERBOSE_USE
    G4cout << "\n\n========= CopyNumbers arrays initialization ======================= \n"<<G4endl;
#endif

    size_t VoxelNumberInPhantom = VoxXNumber*VoxYNumber*VoxZNumber;

    CopyNumberXPos = new G4float[VoxelNumberInPhantom];
    CopyNumberYPos = new G4float[VoxelNumberInPhantom];
    CopyNumberZPos = new G4float[VoxelNumberInPhantom];

    CopyNumberMassSize = new G4float[VoxelNumberInPhantom];
    CopyNumberRegionNameMap = new G4String[VoxelNumberInPhantom];

    InitializeCNMassColourRegNamePosAxisIDs();
    ConstructVoxDcmContainerGeometry();
}

void G4TVolumeConstruction::setDcmRegionMinDensity(G4double n ){
    DcmRegionsMinDensityMap[RegionFromInputInc-1] = n;
    UseDcmRegionsMinDensityMap[RegionFromInputInc-1] = true;
}
void G4TVolumeConstruction::setDcmRegionMaxDensity(G4double n ){
    DcmRegionsMaxDensityMap[RegionFromInputInc-1] = n;
    UseDcmRegionsMaxDensityMap[RegionFromInputInc-1] = true;
}
void G4TVolumeConstruction::setVoxelsegmentedMaterial(G4int n ){

    //G4cout << "\n" << __FUNCTION__<<G4endl;

    if(GeometryFileType == "VOXEL"){ // in case of voxelized , the creation of region is directly here without any density segmentation
        //G4cout << "\n VOXEL TYPE"<<G4endl;

        VoxRegionMaterialName = MaterialIDName[n];
        CreateRegionVoxelsData();
    }else { // we save it to use with each region in a loop

    }

    VoxelMatForSegMap[RegionFromInputInc-1] = n;
    UseVoxelMatForSegMap[RegionFromInputInc-1] = true;

    VoxelMatForSegMapMap[RegionFromInputInc-1][n] = n;
}

//############################################################ For VOXEL and VoxIDs Geometries

void G4TVolumeConstruction::setVoxDefaultMaterialName(G4String n ){// get default material and create contaner Volume
    //G4cout << __FUNCTION__<<G4endl;

#if VERBOSE_USE
    G4cout << "\n\n========= Voxels materials filling and CopyNumbers arrays initialization ======================= \n"<<G4endl;
#endif
    VoxDefaultMaterialName = n;
    defaultMat = CreatedMaterials[VoxDefaultMaterialName]; //its used by container;

    G4int maxID=0; for ( auto it = MaterialIDName.begin(); it != MaterialIDName.end(); ++it  ){if(maxID < it->first){maxID = it->first;}}
    G4bool isDefined = false;
    VoxelsMaterials.clear(); // fill the VoxelsMaterials for VOXEL and VoxIDs
    for(int d = 0 ; d < maxID ;d++ ){
        isDefined = false;
        for ( auto it = MaterialIDName.begin(); it != MaterialIDName.end(); ++it  ){
            if( d == it->first){
                isDefined = true;
                VoxelsMaterials.push_back(CreatedMaterials[MaterialIDName[d]]);
                MaterialColour.push_back(G4Colour( (G4double)G4UniformRand(), (G4double)G4UniformRand()  , (G4double)G4UniformRand()  , 1. ));
                //RegionNameColour[MaterialIDName[d]] = G4Colour( (G4double) G4UniformRand(),(G4double)G4UniformRand(),(G4double)G4UniformRand(), 1.);
                //std::cout << "ID:"<< d <<  " - Material Name: " << MaterialIDName[d] << " Density(g/cm3): "<< VoxelsMaterials[d]->GetDensity() * G4Density_to_gPerCm3 << std::endl;
                continue;
            }
        }
        if( isDefined == false){ // in case the user jump or ignors an id we fill it by default
            VoxelsMaterials.push_back(defaultMat); //
        }
    }


    //if(GeometryFileType == "VoxIDs" || GeometryFileType == "VOXEL"){ //we dont need it fo ICRPFiles type because it filled from file

    size_t VoxelNumberInPhantom = VoxXNumber*VoxYNumber*VoxZNumber;
    MateIDs = new size_t[VoxelNumberInPhantom];

    //std::cout << "VoxelNumberInPhantom ======================= " << VoxelNumberInPhantom << std::endl;

    CopyNumberXPos = new G4float[VoxelNumberInPhantom];
    CopyNumberYPos = new G4float[VoxelNumberInPhantom];
    CopyNumberZPos = new G4float[VoxelNumberInPhantom];

    CopyNumberMassSize = new G4float[VoxelNumberInPhantom];

    CopyNumberRegionNameMap = new G4String[VoxelNumberInPhantom];

    for(size_t hg = 0; hg < VoxelNumberInPhantom ;hg++ ){
        MateIDs[hg] = 0 ;
    }
    //}

    InitializeCNMassColourRegNamePosAxisIDs();
    ConstructVoxDcmContainerGeometry();

}

//############################################################ For all voxelized Geometries

// we use vector to save the Region name despite if we have the same name of region to update ; their data will be save not like in case of
// map using, that the seconde same name will modify the data for first befor begin of regions specification generation, which gives an
// unexpected region data
void G4TVolumeConstruction::setVoxRegionName(G4String n){

    //G4cout << "\nData Initialization for "<< RegionFromInputInc << "-ieme region " << n <<G4endl;

    VoxRegionName = n ;
    DcmRegionsNames.push_back(VoxRegionName);
    UseDcmRegionsMinDensityMap.push_back(false);
    UseDcmRegionsMaxDensityMap.push_back(false);
    DcmRegionsMinDensityMap.push_back(0);
    DcmRegionsMaxDensityMap.push_back(0);

    UseVoxelMatForSegMap.push_back(false);
    VoxelMatForSegMap.push_back(0);

    std::map<G4int,G4int> b;
    VoxelMatForSegMapMap.push_back(b);

    //G4cout << "RegionFromInputInc " << RegionFromInputInc <<G4endl;

    RegionFromInputInc++;
}
void G4TVolumeConstruction::setVoxRegionMinX(G4int n){    VoxRegionMinX = n ; DcmRegionsMinX.push_back(VoxRegionMinX) ; PhantomLimits.push_back("") ;} // despite if we have the same name of region to updete ; their data
void G4TVolumeConstruction::setVoxRegionMaxX(G4int n){    VoxRegionMaxX = n ; DcmRegionsMaxX.push_back(VoxRegionMaxX) ;} // will be save not like in case of map using, that the seconde
void G4TVolumeConstruction::setVoxRegionMinY(G4int n){    VoxRegionMinY = n ; DcmRegionsMinY.push_back(VoxRegionMinY) ;} // same name will modify the data for first befor begin of
void G4TVolumeConstruction::setVoxRegionMaxY(G4int n){    VoxRegionMaxY = n ; DcmRegionsMaxY.push_back(VoxRegionMaxY) ;} // regions specification generation, which gives an unexpected
void G4TVolumeConstruction::setVoxRegionMinZ(G4int n){    VoxRegionMinZ = n ; DcmRegionsMinZ.push_back(VoxRegionMinZ) ;} // region data
void G4TVolumeConstruction::setVoxRegionMaxZ(G4int n){    VoxRegionMaxZ = n ; DcmRegionsMaxZ.push_back(VoxRegionMaxZ) ;}
void G4TVolumeConstruction::setAllGeomAsMinMaxVoxRegionLimits(){

    //G4cout << "\n\n\n\n\n\n\n\n\n\n" << __FUNCTION__<<G4endl;

    VoxRegionMinX = 0 ;            DcmRegionsMinX.push_back(VoxRegionMinX) ;
    VoxRegionMaxX = VoxXNumber-1 ; DcmRegionsMaxX.push_back(VoxRegionMaxX) ;
    VoxRegionMinY = 0 ;            DcmRegionsMinY.push_back(VoxRegionMinY) ;
    VoxRegionMaxY = VoxYNumber-1 ; DcmRegionsMaxY.push_back(VoxRegionMaxY) ;
    VoxRegionMinZ = 0 ;            DcmRegionsMinZ.push_back(VoxRegionMinZ) ;
    VoxRegionMaxZ = VoxZNumber-1 ; DcmRegionsMaxZ.push_back(VoxRegionMaxZ) ;
    PhantomLimits.push_back("all") ;
}
void G4TVolumeConstruction::setVoxRegionsFractions(G4String n, double l){
    RegionRegionFraction[VoxRegionName][n] = l;

}
void G4TVolumeConstruction::InitializeCNMassColourRegNamePosAxisIDs(){

    VoxelVolume = 8*(VoxXHalfSize*VoxYHalfSize*VoxZHalfSize);
#if VERBOSE_USE
    G4cout << " VoxelVolume (mm3) : " << VoxelVolume << G4endl;
#endif
    // CN coordinate and position
    G4float Voxel0PosX = -((VoxXNumber*VoxXHalfSize) - VoxXHalfSize) + VoxContainerPos.getX();
    G4float Voxel0PosY = -((VoxYNumber*VoxYHalfSize) - VoxYHalfSize) + VoxContainerPos.getY();
    G4float Voxel0PosZ = -((VoxZNumber*VoxZHalfSize) - VoxZHalfSize) + VoxContainerPos.getZ();

    //size_t* CnCor = 0;
    size_t Cn_inc = 0;
    for(size_t f = 0; f < VoxZNumber ;f++ ){
        for(size_t g = 0; g < VoxYNumber ;g++ ){
            for(size_t d = 0; d < VoxXNumber ;d++ ){

                CopyNumberXPos[Cn_inc] = Voxel0PosX + d * 2 * VoxXHalfSize;
                CopyNumberYPos[Cn_inc] = Voxel0PosY + g * 2 * VoxYHalfSize;
                CopyNumberZPos[Cn_inc] = Voxel0PosZ + f * 2 * VoxZHalfSize;
                //G4cout << Cn_inc << " " << f << " " << g << " " << d << G4endl;

                CopyNumberMassSize[Cn_inc]= VoxelsMaterials[MateIDs[Cn_inc]]->GetDensity() * G4Density_to_kgPerMm3 * VoxelVolume; /*density ; e+21 in g/mm3 , and e+18 g/cm3 )*/
                if(GeometryFileType == "DICOM"){
                    CopyNumberRegionNameMap[Cn_inc] = VoxelsMaterials[MateIDs[Cn_inc]]->GetName();
                }
                else{
                    CopyNumberRegionNameMap[Cn_inc] = "VOXEL";
                }
                //G4float ff = VoxelsMaterials[MateIDs[Cn_inc]]->GetDensity() * G4Density_to_kgPerMm3 * VoxelVolume ;
                //G4double nn = VoxelsMaterials[MateIDs[Cn_inc]]->GetDensity() * G4Density_to_kgPerMm3 * VoxelVolume ;
                //G4cout << " Float " << ff << "  Vs Double " << nn  << G4endl;

                Cn_inc++;
            }
        }
    }
}
// called from commands for VOXEL at each cmd and from Construct for DICOM in a loop of all registered regions and from GenerateDataForVoxelsIdsFilePhantom in a loop of all registered regions
void G4TVolumeConstruction::CreateRegionVoxelsData(){

    //std::cout << "------MateIDs.size()------- "<<sizeof(MateIDs)<< std::endl;

    std::ostringstream RegionsDataText;

    G4bool IsRegDefined = false;
    for(G4int nn = 0 ; nn < OrganNamesVector.size() ; nn++){

        if(OrganNamesVector[nn] == VoxRegionName){
            IsRegDefined = true;
            break;
        }
    }

    if(IsRegDefined){
        RegionsDataText << "\n ***************** Updating the Region : "<< VoxRegionName << " *********** "<< G4endl;
    }else {
        RegionsDataText << "\n ***************** New Region : "<< VoxRegionName << " ***************** "<< G4endl;
        RegionNumberOfCN[VoxRegionName] = 0;
    }

    std::vector<G4int> RelativenumbersInXForOrgan; //
    std::vector<G4int> RelativenumbersInYForOrgan; // relative to Z
    std::vector<G4int> RelativenumbersInZForOrgan; // relative to Y

    if(VoxRegionMaxX > VoxXNumber){VoxRegionMaxX = VoxXNumber - 1 ; }
    if(VoxRegionMinX > VoxXNumber){VoxRegionMinX = VoxXNumber - 1 ; }
    if(VoxRegionMaxY > VoxYNumber){VoxRegionMaxY = VoxYNumber - 1 ; }
    if(VoxRegionMinY > VoxYNumber){VoxRegionMinY = VoxYNumber - 1 ; }
    if(VoxRegionMaxZ > VoxZNumber){VoxRegionMaxZ = VoxZNumber - 1 ; }
    if(VoxRegionMinZ > VoxZNumber){VoxRegionMinZ = VoxZNumber - 1 ; }

    if(VoxRegionMinX == VoxRegionMaxX){ // MinCol lk=1 and maxCol lk=2
        RegionsDataText << " MinCol                     : " << VoxRegionMinX << G4endl;
        RegionsDataText << " MaxCol                     : " << VoxRegionMaxX << G4endl;

        RelativenumbersInXForOrgan.push_back(VoxRegionMinX);
        RegionsDataText << " Total X voxels             : " << RelativenumbersInXForOrgan.size() << G4endl;

    }else{

        G4int nn = VoxRegionMinX;

        RegionsDataText << " MinCol                     : " << VoxRegionMinX << G4endl;

        for(nn ; nn < VoxRegionMaxX+1 ; nn++){
            RelativenumbersInXForOrgan.push_back(nn);
        }
        RegionsDataText << " MaxCol                     : " << VoxRegionMaxX << G4endl;
        RegionsDataText << " Total X voxels             : " << RelativenumbersInXForOrgan.size() << G4endl;
    }

    //RegionsDataText << "For Y" << G4endl;
    if(VoxRegionMinY == VoxRegionMaxY){ // MinCol lk=1 and maxCol lk=2
        RegionsDataText << " MinRow                     : " << VoxRegionMinY << G4endl;
        RegionsDataText << " MaxRow                     : " << VoxRegionMaxY << G4endl;

        RelativenumbersInYForOrgan.push_back(VoxRegionMinY);
        RegionsDataText << " Total Y voxels             : " << RelativenumbersInYForOrgan.size() << G4endl;

    }else{

        G4int nn = VoxRegionMinY;

        RegionsDataText << " MinRow                     : " << VoxRegionMinY << G4endl;

        for(nn ; nn < VoxRegionMaxY+1 ; nn++){
            RelativenumbersInYForOrgan.push_back(nn);
        }
        RegionsDataText << " MaxRow                     : " << VoxRegionMaxY << G4endl;
        RegionsDataText << " Total Y voxels             : " << RelativenumbersInYForOrgan.size() << G4endl;

    }

    //G4cout << "For Z" << G4endl;
    if(VoxRegionMinZ == VoxRegionMaxZ){ // MinCol lk=1 and maxCol lk=2
        RegionsDataText << " MinSlice                   : " << VoxRegionMinZ << G4endl;
        RegionsDataText << " MaxSlice                   : " << VoxRegionMaxZ << G4endl;

        RelativenumbersInZForOrgan.push_back(VoxRegionMinZ);
        RegionsDataText << " Total Z voxels             : " << RelativenumbersInZForOrgan.size() << G4endl;

    }else{

        G4int nn = VoxRegionMinZ;

        RegionsDataText << " MinSlice                   : " << VoxRegionMinZ << G4endl;

        for(nn ; nn < VoxRegionMaxZ+1 ; nn++){
            RelativenumbersInZForOrgan.push_back(nn);
        }
        RegionsDataText << " MaxSlice                   : " << VoxRegionMaxZ << G4endl;
        RegionsDataText << " Total Z voxels             : " << RelativenumbersInZForOrgan.size() << G4endl;

    }

    //G4ThreeVector CopyNumPos;
    G4int Xcordinate = 0 ;

    RegionsDataText << " Last CN Number             : " << RegionNumberOfCN[VoxRegionName] << G4endl;
    RegionsDataText << " Last Region Volume(cm3)    : " << OrganNameVolumeMap[VoxRegionName]/1000. << G4endl;
    RegionsDataText << " Last Region Mass(Kg)       : " << OrganNameMassMap[VoxRegionName] << G4endl;
    RegionsDataText << " Last Region Density(g/cm3) : " << OrganNameDensityMap[VoxRegionName]*G4Density_to_kgPerMm3ToGPercm3 << G4endl;

    G4bool IsMatDefined = false;
    G4int MatIndex = 0;
    G4int Seg= 0;

    if(GeometryFileType == "VOXEL"){

        if(IsRegDefined){
            for(G4int nn = 0 ; nn < VoxelsMaterials.size() ; nn++){
                if(VoxelsMaterials[nn]->GetName() == VoxRegionMaterialName){
                    IsMatDefined = true;
                    MatIndex = nn ;
                }
            }
            if(!IsMatDefined){
                VoxelsMaterials.push_back(CreatedMaterials[VoxRegionMaterialName]);
                MatIndex = VoxelsMaterials.size() - 1 ;
            }
        }else {
            VoxelsMaterials.push_back(CreatedMaterials[VoxRegionMaterialName]);
            MatIndex = VoxelsMaterials.size() - 1 ;
        }

        VoxelMass = CreatedMaterials[VoxRegionMaterialName]->GetDensity()*G4Density_to_kgPerMm3 * VoxelVolume ; /*density ; e+21 in g/mm3 , and e+18 g/cm3 )*/
        //G4cout << " VoxelMass (g) : " << VoxelMass << G4endl;

        if(IsRegDefined){
        }else {
            RegionNameColour[VoxRegionName] = G4Colour( G4UniformRand(), G4UniformRand() , G4UniformRand() , 1. );
        }

    }else if (GeometryFileType == "DICOM" || GeometryFileType == "VoxIDs") {

        if(UseDcmRegionMinDensity && UseDcmRegionMaxDensity){
            Seg = 1;
            RegionsDataText << " Density Segmentation(g/cm3): " << Seg << " - Min & Max " << G4endl;
            RegionsDataText << " Density Min Max (g/cm3)    : " << DcmRegionMinDensity*G4Density_to_kgPerMm3 * G4Density_to_kgPerMm3ToGPercm3 <<  " " << DcmRegionMaxDensity*G4Density_to_kgPerMm3 * G4Density_to_kgPerMm3ToGPercm3 << G4endl;
        }
        else if(UseDcmRegionMinDensity && !UseDcmRegionMaxDensity){
            Seg = 2;
            RegionsDataText << " Density Segmentation(g/cm3): " << Seg << " - Min & !Max " << G4endl;
            RegionsDataText << " Density Min Max (g/cm3)    : " << DcmRegionMinDensity*G4Density_to_kgPerMm3 * G4Density_to_kgPerMm3ToGPercm3 <<  " ---" << G4endl;
        }
        else if(!UseDcmRegionMinDensity && UseDcmRegionMaxDensity){
            Seg = 3;
            RegionsDataText << " Density Segmentation(g/cm3): " << Seg << " - !Min & Max " << G4endl;
            RegionsDataText << " Density Min Max (g/cm3)    : " << "----- " << DcmRegionMaxDensity*G4Density_to_kgPerMm3 * G4Density_to_kgPerMm3ToGPercm3<< G4endl;
        }
        else if(!UseDcmRegionMinDensity && !UseDcmRegionMaxDensity){
            Seg = 0;
            RegionsDataText << " Density Min Max (g/cm3)    : " << Seg << " - !Min & !Max " << G4endl;
        }

        if(UseVoxelMat){
            RegionsDataText << " Material ID for Seg        : " << VoxelMatForSeg << " is used " << G4endl;
        }else {
            RegionsDataText << " Material Segmentation      : not used " << G4endl;
        }

        //G4cout << RegionsDataText.str().c_str() << G4endl;

    }

    RegionMass = 0.;
    G4int RegionGoodCN = 0;

    //std::cout << "0 " << std::endl;

    for(G4int bb = 0 ; bb < RelativenumbersInZForOrgan.size() ; bb++){
        //G4cout << " For relative Z index : " << RelativenumbersInZForOrgan[OrganChosenName][bb] << G4endl;

        for(G4int cc = 0 ; cc < RelativenumbersInYForOrgan.size() ; cc++){
            //G4cout << " For relative Y index : " << RelativenumbersInYForOrgan[OrganChosenName][cc] << G4endl;

            for(G4int dd = 0 ; dd < RelativenumbersInXForOrgan.size() ; dd++){

                // from if Z=5 "RelativenumbersInZForOrgan[orgId].size()=5"; then From 0 to 5 ; we have 0,1,2,3,4 are filled which means 5*NVoxelY*NVoxelX are filled
                G4int NumInFilledZ ; // in the filled Z .
                G4int NumInFilledY ; // in the filled Y .
                G4int NumRestOfX ; // in the not filled Y
                /*
                if(GeometryFileType == "VoxIDs"){ // because in the dta description of max an min limits, the cosidere that 1 is the min for X Y Z, and max is total number of voxels and not total number of voxels mines 1,
                    NumInFilledZ = (RelativenumbersInZForOrgan[bb]-1)*(VoxYNumber)*VoxXNumber; // in the filled Z .
                    NumInFilledY = (RelativenumbersInYForOrgan[cc]-1)*VoxXNumber; // in the filled Y .
                    NumRestOfX = RelativenumbersInXForOrgan[dd]; // in the not filled Y
                }else {
                    NumInFilledZ = (RelativenumbersInZForOrgan[bb])*(VoxYNumber)*VoxXNumber; // in the filled Z .
                    NumInFilledY = (RelativenumbersInYForOrgan[cc])*VoxXNumber; // in the filled Y .
                    NumRestOfX = RelativenumbersInXForOrgan[dd]; // in the not filled Y
                }
*/
                NumInFilledZ = (RelativenumbersInZForOrgan[bb])*(VoxYNumber)*VoxXNumber; // in the filled Z .
                NumInFilledY = (RelativenumbersInYForOrgan[cc])*VoxXNumber; // in the filled Y .
                NumRestOfX = RelativenumbersInXForOrgan[dd]; // in the not filled Y

                Xcordinate = NumInFilledZ + NumInFilledY + NumRestOfX;

                if(Xcordinate > VoxZNumber*VoxYNumber*VoxXNumber){
                    G4cout << "\n\n\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!! CN > VoxZNumber*VoxYNumber*VoxXNumbe !!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;
                    G4cout << NumInFilledZ << "+" << NumInFilledY  << "+" << NumInFilledY << " = " << Xcordinate << G4endl;
                    G4cout << VoxZNumber*VoxYNumber*VoxXNumber << " " << (RelativenumbersInZForOrgan[bb])*VoxYNumber*VoxXNumber << " " << (RelativenumbersInYForOrgan[cc])*VoxXNumber << " " << RelativenumbersInXForOrgan[dd] << G4endl;
                    G4cout << " The voxel ID " << Xcordinate << " is out of the MatIDs size " << VoxZNumber*VoxYNumber*VoxXNumber << G4endl;
                    continue;
                }else if(Xcordinate < 0){
                    G4cout << "\n\n\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!! CN < 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;
                    G4cout << NumInFilledZ << "+" << NumInFilledY  << "+" << NumInFilledY << " = " << Xcordinate << G4endl;
                    G4cout << " The voxel ID " << Xcordinate << " can't be negative" << G4endl;
                    continue;
                }

                //CopyNumPos = G4ThreeVector( CopyNumberXPos[Xcordinate], CopyNumberYPos[Xcordinate], CopyNumberZPos[Xcordinate]);
                //CopyNumPos.rotateX(VoxContainerRot.x()); CopyNumPos.rotateY(VoxContainerRot.y()); CopyNumPos.rotateZ(VoxContainerRot.z());
                //CopyNumPos.transform(VoxContainerRot);

                if(GeometryFileType == "VOXEL"){
                    MateIDs[Xcordinate] = MatIndex;
                    CopyNumberMassSize[Xcordinate] = CreatedMaterials[VoxRegionMaterialName]->GetDensity()*G4Density_to_kgPerMm3*VoxelVolume; // all region Mass must be edited for new Region
                    CopyNumberRegionNameMap[Xcordinate] = VoxRegionName;
                    RegionMass += CopyNumberMassSize[Xcordinate];
                    RegionNumberOfCN[VoxRegionName]++;
                    //G4cout <<  " Region: "<< VoxRegionName <<  " CN: " << Xcordinate << " Mat: " << VoxelsMaterials[VoxelsMaterials.size()-1]->GetName() << G4endl;
                }
                else if (GeometryFileType == "VoxIDs") { // no density segmentation available because all the mat ids are already defined by the user in material constructing

                    if(UseVoxelMat){
                        // for multiple materials
                        //for ( auto it = MaterialsIDsForSegmentationMap.begin(); it != MaterialsIDsForSegmentationMap.end(); ++it  ){

                        //std::cout << Xcordinate << " " <<MateIDs[Xcordinate] << " " << VoxelMatForSeg << G4endl;

                        if(MateIDs[Xcordinate] == VoxelMatForSeg){

                            //std::cout << "--------------------------------- "<<MateIDs[Xcordinate] << " " << VoxelMatForSeg << std::endl;

                            CopyNumberRegionNameMap[Xcordinate] = VoxRegionName;
                            RegionMass += CopyNumberMassSize[Xcordinate];
                            RegionNumberOfCN[VoxRegionName]++;
                            RegionGoodCN = Xcordinate;
                            continue;
                        }
                    }else {
                        //std::cout << "6 " << std::endl;

                        CopyNumberRegionNameMap[Xcordinate] = VoxRegionName;
                        RegionMass += CopyNumberMassSize[Xcordinate];
                        RegionNumberOfCN[VoxRegionName]++;
                    }
                }
                else if (GeometryFileType == "DICOM") {

                    G4double Densityy = CopyNumberMassSize[Xcordinate]/VoxelVolume;
                    //G4cout << DcmRegionMinDensity*G4Density_to_kgPerMm3 << " =< " <<  Densityy  << " =< " << DcmRegionMaxDensity*G4Density_to_kgPerMm3 << G4endl;

                    if(Seg == 0){ //
                        CopyNumberRegionNameMap[Xcordinate] = VoxRegionName;
                        RegionMass += CopyNumberMassSize[Xcordinate];
                        RegionNumberOfCN[VoxRegionName]++;
                    }
                    else if(Seg == 1){
                        if(DcmRegionMinDensity*G4Density_to_kgPerMm3 <= Densityy  &&  Densityy <= DcmRegionMaxDensity*G4Density_to_kgPerMm3){
                            CopyNumberRegionNameMap[Xcordinate] = VoxRegionName;
                            RegionMass += CopyNumberMassSize[Xcordinate];
                            RegionNumberOfCN[VoxRegionName]++;
                        }
                    }
                    else if(Seg == 2){
                        if(DcmRegionMinDensity*G4Density_to_kgPerMm3 <= Densityy){
                            CopyNumberRegionNameMap[Xcordinate] = VoxRegionName;
                            RegionMass += CopyNumberMassSize[Xcordinate];
                            RegionNumberOfCN[VoxRegionName]++;
                        }
                    }
                    else if(Seg == 3){
                        if(Densityy <= DcmRegionMaxDensity*G4Density_to_kgPerMm3){
                            CopyNumberRegionNameMap[Xcordinate] = VoxRegionName;
                            RegionMass += CopyNumberMassSize[Xcordinate];
                            RegionNumberOfCN[VoxRegionName]++;
                        }
                    }
                    //G4coutOrganName <<  " Region: "<< VoxRegionName <<  " CN: " << Xcordinate << " Mat: " << VoxelsMaterials[VoxelsMaterials.size()-1]->GetName() << G4endl;
                }
            }
        }
    }


    RegionVolume = VoxelVolume*RegionNumberOfCN[VoxRegionName]; // this because the voxel volume is cte the variable is the
    // number of Voxels belonging to the region( variable in updating), and the mass too is variable from voxel to another.

    if(IsRegDefined){
        OrganNameMassMap[VoxRegionName] = OrganNameMassMap[VoxRegionName] + RegionMass; // should be in Kg
        OrganNameVolumeMap[VoxRegionName] = RegionVolume/1000.; // should be in cm3
        OrganNameDensityMap[VoxRegionName] = (OrganNameMassMap[VoxRegionName]/RegionVolume)*G4Density_to_kgPerMm3ToGPercm3; // from mm3 to cm3 // should be in g/cm3

    }else {

        RegionNameColour[VoxRegionName] = G4Colour( (G4double) G4UniformRand(),(G4double)G4UniformRand(),(G4double)G4UniformRand(), 1.);
        OrganNamesVector.push_back(VoxRegionName);
        OrganNameMassMap[VoxRegionName] = RegionMass; // should be in Kg
        OrganNameVolumeMap[VoxRegionName] = RegionVolume/1000.; // should be in cm3
        OrganNameDensityMap[VoxRegionName] = (OrganNameMassMap[VoxRegionName]/RegionVolume)*G4Density_to_kgPerMm3ToGPercm3; // from mm3 to cm3 // should be in g/cm3
    }

    //G4cout << " A2 "<< G4endl;

    if(IsRegDefined){
        RegionsDataText << " Is First building          : " << "no " << G4endl;
    }else {
        RegionsDataText << " Is First building          : " << "yes" << G4endl;
    }
    if(GeometryFileType == "VOXEL"){
        RegionsDataText << " Material Index             : " << MatIndex << G4endl;
        RegionsDataText << " Material Name              : " << VoxelsMaterials[MatIndex]->GetName() << G4endl;
        if(IsRegDefined){
            if(IsMatDefined){
                RegionsDataText << " Different Material !       : " << "no" << G4endl;
            }else {
                RegionsDataText << " Different Material !       : " << "yes" << G4endl;
            }
        }
    }

    RegionsDataText << " hx, hy, hz (mm)            : " << 2*VoxXHalfSize<< " " << 2*VoxYHalfSize << " " << 2*VoxZHalfSize << G4endl;
    RegionsDataText << " Voxel Volume(cm3)          : " << VoxelVolume/1000. << G4endl;
    if(GeometryFileType == "VOXEL" ){
        RegionsDataText << " Voxel Mass(Kg)             : " << VoxelMass << G4endl;
    }
    if(GeometryFileType == "VoxIDs"){
        if(RegionNumberOfCN[VoxRegionName] != 0 && UseVoxelMat ){
            RegionsDataText << " Voxel Mass(Kg)             : " << RegionMass/RegionNumberOfCN[VoxRegionName] << G4endl;
            RegionsDataText << " MatID, Name, dens(g/cm3)   : " << MateIDs[RegionGoodCN] << " , " << VoxelsMaterials[MateIDs[RegionGoodCN]]->GetName() << " , " << VoxelsMaterials[MateIDs[RegionGoodCN]]->GetDensity()*G4Density_to_kgPerMm3*G4Density_to_kgPerMm3ToGPercm3 << G4endl;
        }
    }
    RegionsDataText << " New CN Number              : " << RegionNumberOfCN[VoxRegionName] << G4endl;
    if(IsRegDefined){
        RegionsDataText << " Region Volume(cm3)         : " << RegionVolume/1000. << G4endl;
        RegionsDataText << " Region Mass(Kg)            : " << RegionMass << G4endl;
    }
    RegionsDataText << " New Region Volume(cm3)     : " << OrganNameVolumeMap[VoxRegionName] << G4endl;
    RegionsDataText << " New Region Mass(Kg)        : " << OrganNameMassMap[VoxRegionName] << G4endl;
    RegionsDataText << " New Density(g/cm3)         : " << OrganNameDensityMap[VoxRegionName] << G4endl;
    RegionsDataText << " ****************************************************\n" << G4endl;

#if VERBOSE_USE
    if(RegionNumberOfCN[VoxRegionName] != 0. ){
        G4cout << RegionsDataText.str().c_str() << G4endl;
    }else {
        if(IsRegDefined){
            G4cout << "\n ***************** Updating the Region : "<< VoxRegionName << " *********** "<< G4endl;
        }else {
            G4cout << "\n ***************** New Region : "<< VoxRegionName << " ***************** "<< G4endl;
        }
        G4cout << " The number of voxels attached to this region is " << RegionNumberOfCN[VoxRegionName] << G4endl;
        G4cout << " ****************************************************\n" << G4endl;
    }
#endif

}
void G4TVolumeConstruction::CreateRegionsDataFile(){

    size_t VoxelNumberInPhantom = VoxXNumber*VoxYNumber*VoxZNumber;
    std::vector<G4double> MaterialsPercentage;
    for ( G4int h = 0; h < VoxelsMaterials.size(); h++  ){ MaterialsPercentage.push_back(0.);}
    for ( G4int h = 0; h < VoxelNumberInPhantom; h++  ){ MaterialsPercentage[MateIDs[h]]++;}
#if VERBOSE_USE
    G4cout << std::setw(20) << std::left << "Material ID " << std::setw(33) << std::left << "Material Name" << std::setw(15) << std::left << "Density(g/cm3) "  << std::setw(15) << std::left << "Voxels Num" << std::setw(15) << std::left << "Phantom %" <<G4endl;
#endif
    for ( G4int h = 0; h < VoxelsMaterials.size(); h++  ){
        if(GeometryFileType == "DICOM"){
            MaterialColour.push_back(G4Colour( (G4double)G4UniformRand(), (G4double)G4UniformRand()  , (G4double)G4UniformRand()  , 1. ));
        }
        G4String jj = VoxelsMaterials[h]->GetName();
        if(VoxDefaultMaterialName == VoxelsMaterials[h]->GetName() ){
            jj.resize(23); jj = jj+"(Default)";
#if VERBOSE_USE
            G4cout << std::setw(20) << std::left << h << std::setw(33) << std::left << jj.c_str() << std::setw(15) << std::left << VoxelsMaterials[h]->GetDensity()* G4Density_to_gPerCm3  << std::setw(15) << std::left << MaterialsPercentage[h] << std::setw(15) << std::left << MaterialsPercentage[h]*100/VoxelNumberInPhantom<<G4endl;
#endif
        } else{
            jj.resize(32);
#if VERBOSE_USE
            G4cout << std::setw(20) << std::left << h << std::setw(33) << std::left << jj.c_str() << std::setw(15) << std::left << VoxelsMaterials[h]->GetDensity()* G4Density_to_gPerCm3 << std::setw(15) << std::left << MaterialsPercentage[h]  << std::setw(15) << std::left << MaterialsPercentage[h]*100/VoxelNumberInPhantom<<G4endl;
#endif
        }
    }

    /*std::ostringstream filname ;
    filname << ResultDirectoryPath << "CreatedRegionsData.txt";
    G4String fileNameString = filname.str();

    if(FILE* file = fopen(fileNameString.c_str(),"a")){

        //G4cout << "Creating file  " << fileNameString.c_str() << " : \n" << G4endl;

        //std::fprintf(file,"%s", RegionsDataText.str().c_str());
        //std::fprintf(stdout,"%s",RegionsDataText.str().c_str());

    }*/


}
void G4TVolumeConstruction::ShowMaterialRegionVoxelsData(){

#if VERBOSE_USE

    if(MaterialNameAsRegionName == true){
        G4cout<<"\n\n========= Generated phantom data by material as a region ====================\n" <<G4endl;
    }else{
        G4cout<<"\n\n========= Generated phantom data by constructing new regions ====================\n" <<G4endl;
    }
    G4cout << "--------------------------------------------------------------------------------------------"<<G4endl;

    G4cout << std::setw(20) << std::left << "Region Name " << " "
           << std::setw(33) << std::left << "Mass(kg)"
           << std::setw(15) << std::left << "Volume(cm3) "
           << std::setw(15) << std::left << "Density(g/cm3) "<<G4endl;

    G4cout << "--------------------------------------------------------------------------------------------"<<G4endl;

    G4cout<<std::setiosflags(std::ios::fixed);
    G4cout.precision(3);
    for(int a = 0; a < OrganNamesVector.size() ;a++ ){
        G4cout << std::setw(20) << std::left << OrganNamesVector[a] << " "                        // organ ID
               << std::setw(33) << std::left << OrganNameMassMap[OrganNamesVector[a]]
               << std::setw(15) << std::left << OrganNameVolumeMap[OrganNamesVector[a]]
               << std::setw(15) << std::left << OrganNameDensityMap[OrganNamesVector[a]]<<G4endl;
    }

    G4cout << "--------------------------------------------------------------------------------------------"<<G4endl;

    G4cout << "Number of created materials by commands : " << MaterialIDName.size() <<G4endl;
    G4cout << "Number of created regions               : " << OrganNamesVector.size() <<G4endl;

    if(GeometryFileType == "VoxIDs" || GeometryFileType == "DICOM" ){
        G4cout << std::setw(15) << std::left << "VoxXNumber "     << std::setw(15) << std::left << VoxXNumber << std::setw(15) << std::left << "VoxXHalfSize " << std::setw(15) << std::left << VoxXHalfSize <<G4endl;
        G4cout << std::setw(15) << std::left << "VoxYNumber "     << std::setw(15) << std::left << VoxYNumber << std::setw(15) << std::left << "VoxYHalfSize " << std::setw(15) << std::left << VoxYHalfSize <<G4endl;
        G4cout << std::setw(15) << std::left << "VoxZNumber "     << std::setw(15) << std::left << VoxZNumber << std::setw(15) << std::left << "VoxZHalfSize " << std::setw(15) << std::left << VoxZHalfSize <<G4endl;
        G4cout << std::setw(15) << std::left << "TotalVoxelNum "  << std::setw(15) << std::left << VoxXNumber*VoxYNumber*VoxZNumber << "\n" <<  G4endl;
        G4cout << std::setw(15) << std::left << "MinX | MaxX : "  << std::setw(15) << std::left << 0/*-VoxXNumber*2*VoxXHalfSize*/ << " | " << VoxXNumber*2*VoxXHalfSize <<G4endl;
        G4cout << std::setw(15) << std::left << "MinY | MaxY : "  << std::setw(15) << std::left << 0/*-VoxYNumber*2*VoxYHalfSize*/ << " | " << VoxYNumber*2*VoxYHalfSize <<G4endl;
        G4cout << std::setw(15) << std::left << "MinZ | MaxZ : "  << std::setw(15) << std::left << 0/*-VoxZNumber*2*VoxZHalfSize*/ << " | " << VoxZNumber*2*VoxZHalfSize <<G4endl;
    }else if(GeometryFileType == "TET"){
        G4cout << std::setw(15) << std::left << "Number of readed tetrahedrons : " << tetData->GetNumTetrahedron() << "/" << tetData->GetTotalNumTetrahedronInFile() << " of entire phantom "<<G4endl;
        G4cout << std::setw(15) << std::left << "MinX | MaxX : "  << std::setw(15) << std::left << tetData->GetPhantomBoxMin().getX() << " | " << tetData->GetPhantomBoxMax().getX() <<G4endl;
        G4cout << std::setw(15) << std::left << "MinY | MaxY : "  << std::setw(15) << std::left << tetData->GetPhantomBoxMin().getY() << " | " << tetData->GetPhantomBoxMax().getY() <<G4endl;
        G4cout << std::setw(15) << std::left << "MinZ | MaxZ : "  << std::setw(15) << std::left << tetData->GetPhantomBoxMin().getZ() << " | " << tetData->GetPhantomBoxMax().getZ() <<G4endl;
        G4cout << std::setw(15) << std::left << "Phantom size : " << std::setw(15) << std::left << phantomSize <<G4endl;
    }
#endif
}

G4int G4TVolumeConstruction::ConvertVoxelPosToXYZVoxelIDs(G4String Axis, G4String MaxMin ,G4double val){

    G4double VHS; G4int VN, VID;
    if( Axis == "X")     { VHS = VoxXHalfSize; VN = VoxXNumber;}
    else if( Axis == "Y"){ VHS = VoxYHalfSize; VN = VoxYNumber;}
    else if( Axis == "Z"){ VHS = VoxZHalfSize; VN = VoxZNumber;}

    G4double VS = VHS*2;
    G4double diff = (VN*VS)/2. - (G4int((VN*VS)/2.));

    G4double PosX0 = -(VN*VS)/2;
    G4int VoxInc =0 ;
#if VERBOSE_USE
    G4cout << "\n\n\nAxis " << Axis << " MaxMin " << MaxMin << " val " << val << G4endl;
    G4cout << "VN " << VN << " VS " << VS << " PosX0 " << PosX0 << " diff " << diff << G4endl;
#endif
    //if( diff == 0.){

    G4double posInc = PosX0; // first VoxInc = 0, first posInc = PosX0,
    while(posInc < val ){

        posInc += VS;
        VoxInc++;
        //G4cout << "posInc " << posInc << " VoxInc " << VoxInc << G4endl;

    }
    //G4cout << "final posInc " << posInc << " final VoxInc " << VoxInc << G4endl;
    if(MaxMin == "Max"){
        if(posInc-VHS <= val && val <= posInc ){

            //G4cout << "Pour Max : " << posInc-VHS << " <= val=" << val << " <= " << posInc << G4endl;

            VID = VoxInc-1;
        }else {

            //G4cout << "Pour Max : " << val << " < " << posInc-VHS << G4endl;

            VID = VoxInc-2;
        }
    }
    if(MaxMin == "Min"){
        if(posInc-VHS >= val && val <= posInc ){

            //G4cout << "Pour Min : " << posInc-VHS << " >= val=" << val << " <= " << posInc << G4endl;

            VID = VoxInc-1;
        }else {

            //G4cout << "Pour Min : " << val << " < " << posInc-VHS << G4endl;

            VID = VoxInc;
        }
    }

    //G4cout << "VID " << VID << " posInc " << posInc << G4endl;

    //}
    return VID;
}
// not finished and not tested yet
G4int G4TVolumeConstruction::ConvertVoxelPosToXYZToID(G4double X, G4double Y, G4double Z){

    G4double XVS = VoxXHalfSize*2;
    G4double YVS = VoxYHalfSize*2;
    G4double ZVS = VoxZHalfSize*2;

    G4double Xdiff = (VoxXNumber*XVS)/2. - (G4int((VoxXNumber*XVS)/2.));
    G4double Ydiff = (VoxYNumber*YVS)/2. - (G4int((VoxYNumber*YVS)/2.));
    G4double Zdiff = (VoxZNumber*ZVS)/2. - (G4int((VoxZNumber*ZVS)/2.));

    G4double PosX0 = -(VoxXNumber*XVS)/2;
    G4double PosY0 = -(VoxYNumber*YVS)/2;
    G4double PosZ0 = -(VoxYNumber*ZVS)/2;

    G4int VoxInc =0 ;

#if VERBOSE_USE
    //G4cout << "\n\n\nAxis " << Axis << " MaxMin " << MaxMin << " val " << val << G4endl;
    G4cout << " XVS " << XVS << " PosX0 " << PosX0 << " Xdiff " << Xdiff << G4endl;
    G4cout << " YVS " << YVS << " PosY0 " << PosY0 << " Ydiff " << Ydiff << G4endl;
    G4cout << " ZVS " << ZVS << " PosZ0 " << PosZ0 << " Zdiff " << Zdiff << G4endl;
#endif

    //if( diff == 0.){

    G4double PosX = PosX0;
    G4double PosY = PosY0;
    G4double PosZ = PosZ0;

    G4double posInc = PosX0; // first VoxInc = 0, first posInc = PosX0,

    G4int ID = 0;
    G4bool Found = false;

    for(G4int z = 0; z <VoxZNumber ; z++){
        PosZ += ZVS;
        for(G4int y = 0; y <VoxYNumber ; y++){
            PosY += YVS;
            for(G4int x = 0; x <VoxXNumber ; x++){
                PosX += XVS;
                if(PosX < PosX+XVS && PosY < PosY+YVS && PosZ < PosZ+ZVS){
                    Found = true;
                    break;
                }
                //G4cout << "posInc " << posInc << " ID " << ID << G4endl;
                ID++;
            }
            if(Found == true){
                break;
            }
        }
        if(Found == true){
            break;
        }
    }
#if VERBOSE_USE
    G4cout << "Final ID " << ID << " PosX " << PosX << " PosY " << PosY << " PosZ " << PosZ << G4endl;
#endif
    return ID;
}
void G4TVolumeConstruction::CreateDefaultVOXELRegionData(){

    RegionMass = 0. ; RegionVolume = 0. ;
    G4int inc = 0;
    //G4cout << " " <<CopyNumberRegionNameMap.size() <<G4endl;

    /*
    for ( auto it = CopyNumberRegionNameMap.begin(); it != CopyNumberRegionNameMap.end(); ++it  ){
        if(it->second =="VOXEL"){
            //MateIDs[Xcordinate] = MatIndex;

            //G4cout << " " << it->first <<G4endl;

            CopyNumberMassSize[it->first] = VoxelsMaterials[MateIDs[it->first]]->GetDensity()*G4Density_to_kgPerMm3*VoxelVolume; // all region Mass must be edited for new Region
            RegionMass += CopyNumberMassSize[it->first];
            RegionCopyNumberPositionMap["VOXEL"][it->first] = G4ThreeVector(CopyNumberXPos[it->first],CopyNumberYPos[it->first],CopyNumberZPos[it->first]); // is needed in points generation with name and CN and Pos
            inc++;
        }
    }
*/
    G4int VoxNum = VoxXNumber * VoxYNumber * VoxZNumber ;

    for(G4int a = 0 ; a < VoxNum ; a++){

        if(CopyNumberRegionNameMap[a] =="VOXEL"){
            //MateIDs[Xcordinate] = MatIndex;

            //G4cout << " " << it->first <<G4endl;

            CopyNumberMassSize[a] = VoxelsMaterials[MateIDs[a]]->GetDensity()*G4Density_to_kgPerMm3*VoxelVolume; // all region Mass must be edited for new Region
            RegionMass += CopyNumberMassSize[a];
            //RegionCopyNumberPositionMap["VOXEL"][a] = G4ThreeVector(CopyNumberXPos[a],CopyNumberYPos[a],CopyNumberZPos[a]); // is needed in points generation with name and CN and Pos
            inc++;
        }
    }

    RegionVolume = VoxelVolume*inc;

    //G4bool IsRegDefined = false;
    //for(G4int nn = 0 ; nn < OrganNamesVector.size() ; nn++){if(OrganNamesVector[nn] == VOXEL){ IsRegDefined = true;break;}}
    //if(!IsRegDefined){ OrganNamesVector.push_back("VOXEL");}
    OrganNamesVector.push_back("VOXEL");
    OrganNameMassMap["VOXEL"] = RegionMass; //should be in Kg
    OrganNameVolumeMap["VOXEL"] = RegionVolume/1000.; //should be in cm3
    OrganNameDensityMap["VOXEL"] = (OrganNameMassMap["VOXEL"]/RegionVolume)*G4Density_to_kgPerMm3ToGPercm3; // from mm3 to cm3 // should be in g/cm3;

}
void G4TVolumeConstruction::ConstructVoxDcmContainerGeometry(){

#if VERBOSE_USE
    G4cout << "\n\n========= Voxels container construction and setting regions data ======================= \n"<<G4endl;
#endif

    ContSolidVoll = new G4Box("Container",VoxXNumber*VoxXHalfSize, VoxYNumber*VoxYHalfSize, VoxZNumber*VoxZHalfSize);

    G4RotationMatrix* rm = new G4RotationMatrix(); rm->rotateX(G4ThreeVector().getX()); rm->rotateY(G4ThreeVector().y()); rm->rotateZ(G4ThreeVector().z());

    if(UseInternalSourceVolume == true){
        G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance() ;
        for (size_t iLV = 0; iLV < pvs->size(); iLV++ ) {
            G4String nmm = (*pvs)[iLV]->GetLogicalVolume()->GetName();
            //G4cout << (*pvs)[iLV]->GetLogicalVolume()->GetName() << " " << nmm << G4endl;
            if(strstr(nmm.c_str(),InternalSourceName)){
                InternalSourceLogicalSolid = (*pvs)[iLV]->GetLogicalVolume()->GetSolid();
                InternalSourceLogicalVolume = (*pvs)[iLV]->GetLogicalVolume();
                InternalSourcePosition = getRegionAbsolutePosition(InternalSourceName);
                WorldPhysicalVolume->GetLogicalVolume()->RemoveDaughter((*pvs)[iLV]);
            }
        }

        G4RotationMatrix* rm = new G4RotationMatrix();
        NewContSolidVoll = new G4SubtractionSolid( "Container" , ContSolidVoll, InternalSourceLogicalSolid, rm, InternalSourcePosition);
        ContLogicalVoll = new G4LogicalVolume(NewContSolidVoll, defaultMat, "Container");
        placeInternalSourceVolume();
    }else{
        ContLogicalVoll = new G4LogicalVolume(ContSolidVoll, defaultMat, "Container");
    }

    //ContLogicalVoll = new G4LogicalVolume( ContSolidVoll, defaultMat , "Container", 0, 0, 0 ); //the material is not important, it will be fully filled by the voxels
    ContPhysicalVoll = new G4PVPlacement(rm,VoxContainerPos, ContLogicalVoll, "Container", WorldPhysicalVolume->GetLogicalVolume(), 0, false, 0);              // Copy number
}
G4VPhysicalVolume* G4TVolumeConstruction::ConstructVoxeDcmGeometry(){

    //G4cout << " - plane " << PlanesToVisualize << " from " << MinPlaneID << " to " << MaxPlaneID << G4endl;

    if(PlanesToVisualize != "all"){
#if VERBOSE_USE
        G4cout << "\n\n========= Construct voxelized geometry for planes visualization ======================= \n"<<G4endl;
#endif
        VisualizationVoxelizedGeometry();

    }else{

        //----- Define voxel logical volume
        voxel_solid = new G4Box( "Voxel", VoxXHalfSize, VoxYHalfSize, VoxZHalfSize);
        voxel_logic = new G4LogicalVolume(voxel_solid,defaultMat,"VoxelLogical", 0,0,0);

        if( ParamType == 0 ) {

#if VERBOSE_USE
            G4cout << "\n\n========= Construct voxelized geometry with G4TPhantomParameterisation (" << ParamType << ") ======================= \n"<<G4endl;
#endif

            //----- Create parameterisation
            G4TPhantomParameterisation* param = new G4TPhantomParameterisation();
            //G4cout << " 1 " << G4endl;
            if(UseVoxelsColour == true){ param->setUseLogVolColour(true);}
            param->SetVoxelXYZData(VoxXNumber, VoxYNumber, VoxZNumber);
            param->SetVoxelDimensions( VoxXHalfSize, VoxYHalfSize, VoxZHalfSize );
            //----- Set number of voxels
            param->SetNoVoxels( VoxXNumber, VoxYNumber, VoxZNumber );
            //----- Set list of materials , the materials have to be organized as the indexs for map in fMaterials are the indices that fills fMateIDs
            param->SetMaterials( VoxelsMaterials );
            //G4cout << " 3 " << G4endl;
            //----- Set list of material indices: for each voxel it is a number that correspond to the index of its material in the vector of materials defined above
            param->SetMaterialIndices( MateIDs );
            param->BuildContainerSolid(ContPhysicalVoll);
            //--- Assure yourself that the voxels are completely filling the fContainer volume
            param->CheckVoxelsFillContainer( ContSolidVoll->GetXHalfLength(),ContSolidVoll->GetYHalfLength(),ContSolidVoll->GetZHalfLength() );
            //----- The G4PVParameterised object that uses the created parameterisation should be placed in the fContainer logical volume
            G4PVParameterised* phantom_phys = new G4PVParameterised("Container",voxel_logic,ContLogicalVoll, kXAxis, VoxXNumber*VoxYNumber*VoxZNumber, param);
            //G4cout << " 4 " << G4endl;
            // if axis is set as kUndefined instead of kXAxis, GEANT4 will do an smart voxel optimisation (not needed if G4RegularNavigation is used)
            //----- Set this physical volume as having a regular structure of type 1, so that G4RegularNavigation is used
            phantom_phys->SetRegularStructureId(1); // if not set, G4VoxelNavigation will be used instead

        } else {

#if VERBOSE_USE
            G4cout << "\n\n========= Construct voxelized geometry with G4TPhantomNestedParameterisation (" << ParamType << ") ======================= \n"<<G4endl;
#endif
            voxel_logic->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));
            // Parameterisation to define the material of each voxel
            G4ThreeVector halfVoxelSize(VoxXHalfSize,VoxYHalfSize,VoxZHalfSize);

            G4TPhantomNestedParameterisation* param1 =  new G4TPhantomNestedParameterisation(halfVoxelSize, VoxelsMaterials);
            if(UseVoxelsColour == true){ param1->setUseLogVolColour(true);}

            param1->SetNoVoxel( VoxXNumber, VoxYNumber, VoxZNumber );
            //----- Set list of material indices: for each voxel it is a number that correspond to the index of its material in the vector of materials defined above
            param1->SetMaterialIndices( MateIDs );

            //--- Slice the phantom along Y axis
            G4String yRepName("RepY");
            G4VSolid* solYRep = new G4Box(yRepName,VoxXNumber*VoxXHalfSize, VoxYHalfSize, VoxZNumber*VoxZHalfSize);
            G4LogicalVolume* logYRep = new G4LogicalVolume(solYRep,defaultMat,yRepName);
            new G4PVReplica(yRepName,logYRep,ContLogicalVoll,kYAxis, VoxYNumber,VoxYHalfSize*2.);
            logYRep -> SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));

            //--- Slice the phantom along X axis
            G4String xRepName("RepX");
            G4VSolid* solXRep = new G4Box(xRepName,VoxXHalfSize,VoxYHalfSize, VoxZNumber*VoxZHalfSize);
            G4LogicalVolume* logXRep = new G4LogicalVolume(solXRep,defaultMat,xRepName);
            new G4PVReplica(xRepName,logXRep,logYRep,kXAxis,VoxXNumber,VoxXHalfSize*2.);
            logXRep -> SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));

            new G4PVParameterised("phantom", voxel_logic,logXRep,kZAxis,VoxZNumber,param1);

            //G4PVParameterised* phantom_phys = new G4PVParameterised("Container",voxel_logic,ContLogicalVoll, kXAxis, VoxXNumber*VoxYNumber*VoxZNumber, param1);
            //phantom_phys->SetRegularStructureId(1); // if not set, G4VoxelNavigation will be used instead
        }
    }
    return WorldPhysicalVolume;
}
void G4TVolumeConstruction::VisualizationVoxelizedGeometry(){

    //G4cout << " - Visualization for plane " << PlanesToVisualize << " from " << MinPlaneID << " to " << MaxPlaneID << G4endl;

    G4String axis = "z";
    voxel_solid = new G4Box( "Voxel", VoxXHalfSize, VoxYHalfSize, VoxZHalfSize);
    int NumberOfVisualizedCN;
    int* VisualizedCN;
    G4int TotalVoxelsNumber = VoxXNumber*VoxYNumber*VoxZNumber;




    if(RegionsToVisualize.size() != 0){

        std::vector<int> CPNumbs;

        for (int i = 0; i < RegionsToVisualize.size(); i++) {

            for (int n = 0; n < TotalVoxelsNumber; n++) {

                if(CopyNumberRegionNameMap[n] == RegionsToVisualize[i]){
                    //G4cout << RegionsToVisualize[i] << G4endl;
                    CPNumbs.push_back(n);
                }
            }
        }
        NumberOfVisualizedCN = CPNumbs.size();
        VisualizedCN = new int[NumberOfVisualizedCN];

        for (int n = 0; n < NumberOfVisualizedCN; n++) {
            VisualizedCN[n] = CPNumbs[n];
        }

#if VERBOSE_USE
        G4cout << " - Number of regions to visualize: " << RegionsToVisualize.size() << ": ";
        for (int i = 0; i < RegionsToVisualize.size(); i++) {
            G4cout << RegionsToVisualize[i] << " - ";
        }
        G4cout << G4endl;
        G4cout << " - Total number of visualized voxels is "<< NumberOfVisualizedCN << G4endl;
#endif

    }else {

        G4int VoxNumbAxis1 = VoxYNumber; G4int VoxNumbAxis2 = VoxXNumber;

        if(PlanesToVisualize == "xy" || PlanesToVisualize == "yx"){
            if(MaxPlaneID > VoxZNumber){MaxPlaneID = VoxZNumber - 1 ;
                //G4cout << " The set ID number " << MaxPlaneID << " exceeds the plane maximum index across Z axis " << VoxZNumber-1  << ". The new Z-Plane ID is " << VoxZNumber-1 << G4endl;
            }
            if(MinPlaneID > VoxZNumber){MinPlaneID = VoxZNumber - 1 ;
                //G4cout << " The set ID number " << MinPlaneID << " exceeds the plane maximum index across Z axis " << VoxZNumber-1  << ". The new Z-Plane ID is " << VoxZNumber-1 << G4endl;
            }

            axis = "z";

            VoxNumbAxis1 = VoxYNumber;
            VoxNumbAxis2 = VoxXNumber;
            //PlanePosition = -((VoxZNumber*VoxZHalfSize) - VoxZHalfSize) + VoxContainerPos.getZ() + PlaneID * 2 * VoxZHalfSize;
        }
        else if(PlanesToVisualize == "xz" || PlanesToVisualize == "zx"){
            if(MaxPlaneID > VoxYNumber){MaxPlaneID = VoxYNumber - 1 ;
                //G4cout << " The set ID number " << MaxPlaneID << " exceeds the plane maximum index across Y axis " << VoxYNumber-1  << ". The new Y-Plane ID is " << VoxYNumber-1 << G4endl;
            }
            if(MinPlaneID > VoxYNumber){MinPlaneID = VoxYNumber - 1 ;
                //G4cout << " The set ID number " << MinPlaneID << " exceeds the plane maximum index across Y axis " << VoxYNumber-1  << ". The new Y-Plane ID is " << VoxYNumber-1 << G4endl;
            }

            axis = "y";

            VoxNumbAxis1 = VoxZNumber;
            VoxNumbAxis2 = VoxXNumber;
            //PlanePosition = -((VoxYNumber*VoxYHalfSize) - VoxYHalfSize) + VoxContainerPos.getY() + PlaneID * 2 * VoxYHalfSize;
        }
        else if(PlanesToVisualize == "yz" || PlanesToVisualize == "zy"){
            if(MaxPlaneID > VoxXNumber){MaxPlaneID = VoxXNumber - 1 ;
                //G4cout << " The set ID number " << MaxPlaneID << " exceeds the plane maximum index across X axis " << VoxXNumber-1  << ". The new X-Plane ID is " << VoxXNumber-1 << G4endl;
            }
            if(MinPlaneID > VoxXNumber){MinPlaneID = VoxXNumber - 1 ;
                //G4cout << " The set ID number " << MinPlaneID << " exceeds the plane maximum index across X axis " << VoxXNumber-1  << ". The new X-Plane ID is " << VoxXNumber-1 << G4endl;
            }

            axis = "x";

            VoxNumbAxis1 = VoxZNumber;
            VoxNumbAxis2 = VoxYNumber;
            //PlanePosition = -((VoxXNumber*VoxXHalfSize) - VoxXHalfSize) + VoxContainerPos.getX() + PlaneID * 2 * VoxXHalfSize;
        }

        std::vector<G4int> RelativeCNNumbersForAxisPlanes; //

        if(MinPlaneID == MaxPlaneID){ // MinCol lk=1 and maxCol lk=2
            RelativeCNNumbersForAxisPlanes.push_back(MinPlaneID);
        }else{
            G4int nn = MinPlaneID;
            for(nn ; nn < MaxPlaneID+1 ; nn++){RelativeCNNumbersForAxisPlanes.push_back(nn);}
        }

        NumberOfVisualizedCN = RelativeCNNumbersForAxisPlanes.size()*VoxNumbAxis1*VoxNumbAxis2;

#if VERBOSE_USE
        G4cout << " - Number of " << PlanesToVisualize << " plans " << RelativeCNNumbersForAxisPlanes.size() << " from " << axis<< "=" << MinPlaneID << " to " << axis<< "=" << MaxPlaneID << G4endl;
        G4cout << " - Total number of visualized voxels " << NumberOfVisualizedCN << G4endl;
#endif

       VisualizedCN = new int[NumberOfVisualizedCN];
        int cn=0;
        int inc=0;

        for(size_t f = 0; f < VoxZNumber ;f++ ){
            for(size_t g = 0; g < VoxYNumber ;g++ ){
                for(size_t d = 0; d < VoxXNumber ;d++ ){

                    for(size_t m = 0; m < RelativeCNNumbersForAxisPlanes.size() ;m++ ){

                        if(PlanesToVisualize == "xy" && RelativeCNNumbersForAxisPlanes[m] == f){
                            VisualizedCN[inc] = cn;
                            inc++;
                        }
                        else if(PlanesToVisualize == "xz" && RelativeCNNumbersForAxisPlanes[m] == g){
                            VisualizedCN[inc] = cn;
                            inc++;
                        }
                        else if(PlanesToVisualize == "yz" && RelativeCNNumbersForAxisPlanes[m] == d){
                            VisualizedCN[inc] = cn;
                            inc++;
                        }
                    }
                    cn++;
                }
            }
        }
    }


    for(int inc = 0; inc < NumberOfVisualizedCN ;inc++ ){

        voxel_logic = new G4LogicalVolume(voxel_solid,VoxelsMaterials[MateIDs[VisualizedCN[inc]]],"VoxelLogical", 0,0,0);
        G4ThreeVector VoxPosition;

        //if(PlanesToVisualize == "xy" || PlanesToVisualize == "yx"){VoxPosition =      G4ThreeVector(CopyNumberXPos[VisualizedCN[inc]],CopyNumberYPos[VisualizedCN[inc]],CopyNumberZPos[VisualizedCN[inc]]);}
        //else if(PlanesToVisualize == "xz" || PlanesToVisualize == "zx"){VoxPosition = G4ThreeVector(CopyNumberXPos[VisualizedCN[inc]],CopyNumberYPos[VisualizedCN[inc]],CopyNumberZPos[VisualizedCN[inc]]);}
        //else if(PlanesToVisualize == "yz" || PlanesToVisualize == "zy"){VoxPosition = G4ThreeVector(CopyNumberXPos[VisualizedCN[inc]],CopyNumberYPos[VisualizedCN[inc]],CopyNumberZPos[VisualizedCN[inc]]);}

        VoxPosition = G4ThreeVector(CopyNumberXPos[VisualizedCN[inc]],CopyNumberYPos[VisualizedCN[inc]],CopyNumberZPos[VisualizedCN[inc]]);
        new G4PVPlacement( new G4RotationMatrix() ,VoxPosition, "PhysVox" , voxel_logic , ContPhysicalVoll , false , 0 , false );

        G4VisAttributes* voxvis = new G4VisAttributes(RegionNameColour[CopyNumberRegionNameMap[VisualizedCN[inc]]]);
        voxvis->SetVisibility(true); voxvis->SetForceSolid(true);
        voxel_logic->SetVisAttributes(voxvis);

        //G4cout << " CN ID " << VisualizedCN[inc] << " Position " << VoxPosition << "  Colour " << RegionNameColour[CopyNumberRegionNameMap[VisualizedCN[inc]]] << " Material name " << VoxelsMaterials[MateIDs[VisualizedCN[inc]]]->GetName() << G4endl;
    }
#if VERBOSE_USE
    G4cout << G4endl;
#endif

    /*
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance() ;
    if(pVVisManager)
    {
        G4UImanager::GetUIpointer()->ApplyCommand("/vis/open OGLSX");

        // Camera setting
        //G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/viewpointVector 1 0 0");
        //G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 20 20");

        G4UImanager::GetUIpointer()->ApplyCommand("/vis/drawLogicalVolume ContLog");
        G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/flush");

        // set print mode to vectored
        G4UImanager::GetUIpointer()->ApplyCommand("/vis/ogl/set/printMode vectored");

        // set print size larger than screen
        G4UImanager::GetUIpointer()->ApplyCommand("/vis/ogl/set/printSize 2000 2000");

        // print
        G4UImanager::GetUIpointer()->ApplyCommand("/vis/ogl/printEPS");
    }
    */

}
//############################################################ For VoxIDs voxelized phantom

/*
to generate A voxelized volume you should get and/or fill the :
VoxIDs voxelized
VoxXNumber VoxYNumber VoxZNumber VoxXHalfSize VoxYHalfSize  VoxZHalfSize
VoxelsMaterials MaterialColour MateIDs
VoxelVolume RegionNameColour[] CopyNumberMassSize[i] CopyNumberRegionNameMap[i] egionCopyNumberColour[i] RegionCopyNumberPosMap[Reg][i] CNXPos[i] CNYPos[i] CNZPos[i] defaultMat

ConstructVoxDcmContainerGeometry()
ConstructVoxeDcmGeometry()
*/

void G4TVolumeConstruction::ReadVoxelsIDsAndFillCNMatIDsMassColour(){

#if VERBOSE_USE
    G4cout<<"\n\n========= Read VOXELS Organs/Tissus/Regions IDs Values ====================\n" <<G4endl;
#endif
    std::ifstream fileR;
    fileR.open(GeometryPath.c_str(), std::ios::binary | std::ios::in);

    G4String line, indicator, word;
    if(fileR.is_open()){

#if VERBOSE_USE
        std::cout << "\nThe file " << GeometryPath << " is opened "<<  std::endl;
#endif
        //VoxelVolume = 8*(VoxXHalfSize*VoxYHalfSize*VoxZHalfSize); // in mm3

        unsigned int Cn_inc = 0;

        size_t VoxelNumberInPhantom = VoxXNumber*VoxYNumber*VoxZNumber;
#if VERBOSE_USE
        std::cout << "\nTotal Number of Voxels to Read: " << VoxelNumberInPhantom <<  std::endl;
#endif
        int ID;

        for(int f = 0; f < VoxZNumber ;f++ ){
            for(int g = 0; g < VoxYNumber ;g++ ){
                for(int d = 0; d < VoxXNumber ;d++ ){

                    fileR >> ID ;

                    //std::cout << Cn_inc << " MateID : " << ID << std::endl;

                    MateIDs[Cn_inc] = ID;
                    //G4cout << Cn_inc << " MateID : " << MateIDs[Cn_inc]  << G4endl;
                    //std::cout << Cn_inc << " MateID : " << ID << std::endl;

                    CopyNumberMassSize[Cn_inc] = VoxelsMaterials[ID]->GetDensity() * G4Density_to_kgPerMm3 * VoxelVolume; /*density ; e+21 in g/mm3 , and e+18 g/cm3 )*/
                    //std::cout << Cn_inc << " MateID : " << ID << std::endl;

                    if(MaterialNameAsRegionName == true){
                        //std::cout << Cn_inc << " VoxelsMaterials.size : " << VoxelsMaterials.size()
                        //          << " OrganNamesVector.size : " << OrganNamesVector.size() << std::endl;

                        bool isIn = false;
                        for ( int df = 0 ; df < OrganNamesVector.size(); df++  ){
                            //std::cout << Cn_inc << " OrganNamesVector[df] : " << OrganNamesVector[df]<< std::endl;
                            //std::cout << Cn_inc << " VoxelsMaterials[ID] : " << VoxelsMaterials[ID]<< std::endl;
                            if(VoxelsMaterials[ID]->GetName() == OrganNamesVector[df] ){ isIn = true;}
                        }
                        if(isIn == false ){
                            //std::cout << Cn_inc << " VoxelsMaterials[ID]->GetName() : " << std::endl;

                            OrganNamesVector.push_back(VoxelsMaterials[ID]->GetName());
                            //std::cout << VoxelsMaterials[ID]->GetName() << std::endl;
                        }
                        //std::cout << Cn_inc << " MateID : " << ID << std::endl;
                        //std::cout << Cn_inc << " MaterialColour[ID].size() : " << MaterialColour.size() << std::endl;

                        RegionNameColour[VoxelsMaterials[ID]->GetName()] = MaterialColour[ID];
                        //std::cout << Cn_inc << " MateID : " << ID << std::endl;
                        OrganNameMassMap[VoxelsMaterials[ID]->GetName()] += CopyNumberMassSize[Cn_inc]; // should be in Kg
                        OrganNameVolumeMap[VoxelsMaterials[ID]->GetName()] += VoxelVolume/1000; // should be in cm3
                        OrganNameDensityMap[VoxelsMaterials[ID]->GetName()] = VoxelsMaterials[ID]->GetDensity() * G4Density_to_kgPerMm3*G4Density_to_kgPerMm3ToGPercm3; // from mm3 to cm3 // should be in g/cm3
                        CopyNumberRegionNameMap[Cn_inc] = VoxelsMaterials[ID]->GetName();
                    }else{
                        CopyNumberRegionNameMap[Cn_inc] = "VOXEL";
                    }

                    //if(ID == 95){ std::cout << Cn_inc << " " << f << " " << g << " " << d << " MateID : " << ID  << " Region name " << VoxelsMaterials[ID]->GetName() << " Material name " << VoxRegionName << " Mass(kg) : " << CopyNumberMassSize[Cn_inc] << " RegionName : " << CopyNumberRegionNameMap[Cn_inc] <<  std::endl;}
                    //std::ostringstream text;
                    //text <<"\r"<< " " << Cn_inc <<"/"<<VoxelNumberInPhantom << " ";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
                    //printf(text.str().c_str());

                    Cn_inc++;
                }
            }
        }
#if VERBOSE_USE
        std::cout << "\nTotal Number of Voxels to Read: " << VoxelNumberInPhantom  << " ===> Number Of readed Voxels: " << Cn_inc <<"\n" <<  std::endl;
#endif
        fileR.close();
    }else{
#if VERBOSE_USE
        std::cout << "\nCanno't read the file " << GeometryPath <<  std::endl;
#endif
        G4String msg = "Canno't read the file \"" + GeometryPath + "\""; G4Exception("Geometry Data", "1", FatalErrorInArgument, msg.c_str());
    }

}
void G4TVolumeConstruction::GenerateDataForVoxelsIdsFilePhantom(){

    //defaultMat = VoxelsMaterials[0]; // used in container Volume
    //ConstructVoxDcmContainerGeometry();

#if VERBOSE_USE
    std::cout << "DcmRegionsNames.size() : " << DcmRegionsNames.size() <<  std::endl;
    std::cout << "DcmRegionsMinX.size() : " << DcmRegionsMinX.size() <<  std::endl;
    std::cout << "UseVoxelMatForSegMap.size() : " << UseVoxelMatForSegMap.size() <<  std::endl;
    std::cout << "VoxelMatForSegMap.size() : " << VoxelMatForSegMap.size() <<  std::endl;
#endif

    OrganNamesVector.clear();
    OrganNameMassMap.clear(); // should be in Kg
    OrganNameVolumeMap.clear(); // should be in cm3
    OrganNameDensityMap.clear();
    for(G4int RNL = 0 ; RNL < DcmRegionsNames.size() ; RNL++){

        //std::cout << "DcmRegionsNames : " << DcmRegionsNames[RNL] <<  std::endl;

        // segmentation data
        VoxRegionName = DcmRegionsNames[RNL];
        VoxRegionMinX = DcmRegionsMinX[RNL];
        VoxRegionMaxX = DcmRegionsMaxX[RNL];
        VoxRegionMinY = DcmRegionsMinY[RNL];
        VoxRegionMaxY = DcmRegionsMaxY[RNL];
        VoxRegionMinZ = DcmRegionsMinZ[RNL];
        VoxRegionMaxZ = DcmRegionsMaxZ[RNL];

        UseDcmRegionMinDensity = UseDcmRegionsMinDensityMap[RNL];
        UseDcmRegionMaxDensity = UseDcmRegionsMaxDensityMap[RNL];
        DcmRegionMinDensity = DcmRegionsMinDensityMap[RNL];
        DcmRegionMaxDensity = DcmRegionsMaxDensityMap[RNL];

        UseVoxelMat = UseVoxelMatForSegMap[RNL];
        VoxelMatForSeg = VoxelMatForSegMap[RNL];
        MaterialsIDsForSegmentationMap = VoxelMatForSegMapMap[RNL] ;

        //std::cout << "VoxelMatForSeg : " << VoxelMatForSeg <<  std::endl;

        CreateRegionVoxelsData();
    }

}

void G4TVolumeConstruction::CalcutateSizeOfVoxelizedArrays(){

    G4int numberOfVoxels = VoxXNumber*VoxYNumber*VoxZNumber;

#if VERBOSE_USE
    G4cout<<"\n\n========= Memory Size needed to hold "<< numberOfVoxels <<" voxels data elements in different arrays, maps and vectors ====================\n" <<G4endl;
#endif

    G4double Gega = 1./(1024.*1024.*1024.);
    G4int nn = 50 ;

    G4double  dataSizes = 0 ;
    G4double v = 0 ;

#if VERBOSE_USE
    G4cout << "\nGeometry construction class " << G4endl;
#endif
    v = G4double(numberOfVoxels)*sizeof(G4float ); G4cout << std::setw(nn) << std::left << " Array Size Of activities in Gega" << v * Gega << G4endl; dataSizes += v;
    v = G4double(numberOfVoxels)*sizeof(G4float ); G4cout << std::setw(nn) << std::left << " Array Size Of cumulated activities in Gega" << v * Gega << G4endl; dataSizes += v;
    v = G4double(numberOfVoxels)*sizeof(G4int); G4cout << std::setw(nn) << std::left << " Array Size Of MatIds in Gega" << v * Gega << G4endl; dataSizes += v;
    v = G4double(numberOfVoxels)*sizeof(G4float )*3; G4cout << std::setw(nn) << std::left << " Array Size Of CN Pos X Y Z in Gega" << v * Gega << G4endl; dataSizes += v;
    v = G4double(numberOfVoxels)*sizeof(G4float ); G4cout << std::setw(nn) << std::left << " Array Size Of CN Masses in Gega" << v * Gega << G4endl; dataSizes += v;
    v = G4double(numberOfVoxels)*sizeof(G4String); G4cout << std::setw(nn) << std::left << " Array Size Of CN Reg Names in Gega" << v * Gega << G4endl; dataSizes += v ;
    v = G4double(numberOfVoxels)*sizeof(G4Colour); G4cout << std::setw(nn) << std::left << " Map Size Of CN Colour in Gega" << v * Gega << G4endl; dataSizes += v;
    v = sizeof(std::map<G4int, G4ThreeVector>) +G4double(numberOfVoxels)*sizeof(G4int)+ G4double(numberOfVoxels)*sizeof(G4ThreeVector); G4cout << std::setw(nn) << std::left << " Map Map Size Of CN Pos XYZ in Gega" << v * Gega << G4endl; dataSizes += v;

    //G4cout << "\nEventsDataGeneration class " << G4endl;
    v = sizeof(std::map<G4int, G4ThreeVector>) +G4double(numberOfVoxels)*sizeof(G4int)+ G4double(numberOfVoxels)*sizeof(G4ThreeVector); G4cout << std::setw(nn) << std::left << " Map Size Of CN Pos XYZ in Gega" << v * Gega << G4endl; dataSizes += v;

    //G4cout << "\nPrimaryGeneratorAction class " << G4endl;
    v = G4double(10000000)*(sizeof(G4double)+sizeof(G4ThreeVector)+sizeof(G4ThreeVector)); G4cout << std::setw(nn) << std::left << " Arrays Size Of Ene, Pos and MomDir in Primary" << v * Gega << G4endl; dataSizes += v;

    //G4cout << "\nRun class " << G4endl;
    v = G4double(numberOfVoxels)*sizeof(G4float )*3; G4cout << std::setw(nn) << std::left << " Array Size Of CN Pos X Y Z in Gega" << v * Gega << G4endl; dataSizes += v;
    v = G4double(numberOfVoxels)*sizeof(G4String); G4cout << std::setw(nn) << std::left << " Array Size Of CN Reg Names in Gega" << v * Gega << G4endl; dataSizes += v;
    v = G4double(numberOfVoxels)*sizeof(G4float ) ; G4cout << std::setw(nn) << std::left << " Array Size Of CN Masses in Gega" << v * Gega << G4endl; dataSizes +=  v;

    v = 2*(sizeof(std::map<unsigned int,G4float >) +G4double(numberOfVoxels)*sizeof(unsigned int)+ G4double(numberOfVoxels)*sizeof(G4float )); G4cout << std::setw(nn) << std::left << " Map ED and ED2 in Gega" << v * Gega << G4endl; dataSizes += v;
    v = (sizeof(std::map<unsigned int,unsigned int>) +G4double(numberOfVoxels)*sizeof(unsigned int)+ G4double(numberOfVoxels)*sizeof(unsigned int)) ; G4cout << std::setw(nn) << std::left << " Map NOFSteps in Gega" << v * Gega << G4endl; dataSizes += v;

    G4cout << "\nResultsCalculation class " << G4endl;
    v = sizeof(std::vector<unsigned int>)+G4double(numberOfVoxels)*sizeof(unsigned int) ; G4cout << std::setw(nn) << std::left << " Vector Size Of CN ID in Gega" << v * Gega << G4endl; dataSizes += v;
    v = 3*(sizeof(std::map<G4int, G4float>)+G4double(numberOfVoxels)*sizeof(G4int)+G4double(numberOfVoxels)*sizeof(G4float)) ; G4cout << std::setw(nn) << std::left << " Map Size Of CN Pos X Y Z in Gega" << v * Gega << G4endl; dataSizes += v;
    v = sizeof(std::map<G4int, G4String>)+G4double(numberOfVoxels)*sizeof(G4int)+G4double(numberOfVoxels)*sizeof(G4String); G4cout << std::setw(nn) << std::left << " Map Size Of CN Reg Names in Gega" << v * Gega << G4endl; dataSizes += v ;
    v = sizeof(std::map<G4int, G4float>)+G4double(numberOfVoxels)*sizeof(G4int)+G4double(numberOfVoxels)*sizeof(G4float); G4cout << std::setw(nn) << std::left << " Map Size Of CN Masses in Gega" << v * Gega << G4endl; dataSizes += v;

    v = (sizeof(std::map<unsigned int,unsigned int>) +G4double(numberOfVoxels)*sizeof(unsigned int) + G4double(numberOfVoxels)*sizeof(unsigned int)) ; G4cout << std::setw(nn) << std::left << " Map NOFSteps in Gega" << v * Gega << G4endl; dataSizes += v;
    v = 13.*(sizeof(std::map<unsigned int,G4float>) +G4double(numberOfVoxels)*sizeof(unsigned int) + G4double(numberOfVoxels)*sizeof(G4float)) ; G4cout << std::setw(nn) << std::left << "Map Variables & mean, var, SDev, rel in Gega" << v * Gega << G4endl; dataSizes += v;

#if VERBOSE_USE
    G4cout << "\n"<< std::setw(nn) << std::left << "Total Size of arrays in Gega" << dataSizes*Gega << G4endl;
    G4cout << std::setw(nn) << std::left << "Size of G4float" << sizeof(G4float) << G4endl;
    G4cout << std::setw(nn) << std::left << "Size of unsigned int" << sizeof(unsigned int) << G4endl;
    G4cout << std::setw(nn) << std::left << "Size of G4double" << sizeof(G4double) << G4endl;
    G4cout << std::setw(nn) << std::left << "Size of G4int" << sizeof(G4int) << G4endl;
#endif

    G4double ICRPPhantomMass = 70;
    G4double PatientPhantomMass = 70;

    G4double ICRPXPhantomsize = 40;
    G4double ICRPYPhantomsize = 40;
    G4double ICRPZPhantomsize = 170;

    G4double PatientXPhantomsize = 40;
    G4double PatientYPhantomsize = 40;
    G4double PatientZPhantomsize = 170;

    G4double ICRPXRegionCenter = 22;
    G4double ICRPYRegionCenter = 19;
    G4double ICRPZRegionCenter = 100;

    G4double PatientXRegionCenter = (ICRPXRegionCenter/ICRPXPhantomsize)*PatientXPhantomsize;
    G4double PatientYRegionCenter = (ICRPYRegionCenter/ICRPYPhantomsize)*PatientYPhantomsize;
    G4double PatientZRegionCenter = (ICRPZRegionCenter/ICRPZPhantomsize)*PatientZPhantomsize;

    G4double ICRPRegionDensity = 1.05;
    G4double ICRPRegionDensityMin = 0.999;
    G4double ICRPRegionDensityMax = 1.05;

    G4double ICRPRegionMass = 1.;
    G4double NewRegionMass = ICRPRegionMass*(PatientPhantomMass/ICRPPhantomMass);

    G4double VoxelDensity = 0;
    G4double VoxelVolume = VoxXHalfSize*VoxYHalfSize*VoxZHalfSize*8;
    G4double RegionMass = 0;

    G4double NewXDistance=0;
    G4double NewYDistance=0;
    G4double NewZDistance=0;
    while(VoxelDensity <= ICRPRegionDensityMax && ICRPRegionDensityMin < VoxelDensity && RegionMass <= NewRegionMass ){

        NewXDistance += VoxXHalfSize;
        NewYDistance += VoxYHalfSize;
        NewZDistance += VoxZHalfSize;

        G4double VoxelXPosition = PatientXRegionCenter + NewXDistance;
        G4double VoxelYPosition = PatientYRegionCenter + NewYDistance;
        G4double VoxelZPosition = PatientZRegionCenter + NewZDistance;

        VoxelDensity = CreatedMaterials[MaterialIDName[MateIDs[ConvertVoxelPosToXYZToID(VoxelXPosition, VoxelYPosition,VoxelZPosition )]]]->GetDensity()*G4Density_to_gPerCm3;
        RegionMass += (VoxelDensity*VoxelVolume);
    }

    // ther rest voxels without segmentation will be saved with material according to the density value

}

//############################################################ For ICRP voxelized phantom

void G4TVolumeConstruction::GenerateICRPMaterialsCommands(){

    std::cout << "\n\n" << __FUNCTION__ <<  std::endl;


    std::string OrganName, MediaName ;
    int IdOrgan , IdMat;
    double HFrac, CFrac, NFrac, OFrac, NaFrac, MgFrac,PFrac,SFrac,ClFrac,KFrac, CaFrac, FeFrac, IFrac;

    G4String line, indicator, word;

    std::ifstream fileR(GeometryPath.c_str());
    G4String TextToWrite = "";
    if(fileR.is_open()){

        //G4cout  << " \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << G4endl;
        std::cout << "Reading "<< GeometryPath << " ... " << std::endl ;

        while (getline(fileR, line)) {

            //G4cout << " the line " << line << G4endl ;

            std::istringstream LineString(line);

            if(LineString.str().empty()){ continue; }

            LineString >> word;  // the white spaces more than 1 are considered as 1

            //G4cout << " the word " << word << G4endl ;
            if(word.empty() || word == "" || word == "#" ){//|| word.isNull()
                //G4cout <<  "word " << word << G4endl;
                continue; }
            if (word == ">>"){
                LineString >> indicator;
                G4cout <<  " --- data indicator : " << indicator << G4endl;
                continue; }
            else{

                if(indicator == "voxels_data"){

                    VoxYNumber = atoi(word.c_str());
                    LineString >> VoxYNumber;
                    LineString >> VoxZNumber;
                    LineString >> VoxXHalfSize;
                    LineString >> VoxYHalfSize;
                    LineString >> VoxZHalfSize;
                    LineString >> ICRPOrgansNumber;
                    LineString >> ICRPMaterialsNumber;

                    std::cout  << VoxXNumber << " "  << VoxYNumber << " " << VoxZNumber << " " << VoxXHalfSize << " " << VoxYHalfSize << " " << VoxZHalfSize << " "  << ICRPOrgansNumber << " " << ICRPMaterialsNumber << std::endl;;
                }
                else if (indicator == "materials_comp_data") {

                    IdMat = atoi(word.c_str());
                    //LineString >> MaterialIDName[IdMat];
                    LineString >> HFrac; ICRPMaterialsCompFra[IdMat].push_back(HFrac);
                    LineString >> CFrac; ICRPMaterialsCompFra[IdMat].push_back(CFrac);
                    LineString >> NFrac; ICRPMaterialsCompFra[IdMat].push_back(NFrac);
                    LineString >> OFrac; ICRPMaterialsCompFra[IdMat].push_back(OFrac);
                    LineString >> NaFrac; ICRPMaterialsCompFra[IdMat].push_back(NaFrac);
                    LineString >> MgFrac; ICRPMaterialsCompFra[IdMat].push_back(MgFrac);
                    LineString >> PFrac; ICRPMaterialsCompFra[IdMat].push_back(PFrac);
                    LineString >> SFrac; ICRPMaterialsCompFra[IdMat].push_back(SFrac);
                    LineString >> ClFrac; ICRPMaterialsCompFra[IdMat].push_back(ClFrac);
                    LineString >> KFrac; ICRPMaterialsCompFra[IdMat].push_back(KFrac);
                    LineString >> CaFrac; ICRPMaterialsCompFra[IdMat].push_back(CaFrac);
                    LineString >> FeFrac; ICRPMaterialsCompFra[IdMat].push_back(FeFrac);
                    LineString >> IFrac; ICRPMaterialsCompFra[IdMat].push_back(IFrac);
                    std::cout << IdMat <<" "<<HFrac<<" "<<CFrac<<" "<<NFrac<<" "<<OFrac<<" "<<NaFrac<<" "<<MgFrac<<" "<<PFrac<<" "<<SFrac<<" "<<ClFrac<<" "<<KFrac<<" "<<CaFrac<<" "<<FeFrac<<" "<<IFrac << std::endl;

                }
                else if (indicator == "organId_materialId_density") {

                    IdOrgan = atoi(word.c_str());
                    LineString >> OrganIDMatIDMap[IdOrgan];
                    LineString >> OrganIDDensityMap[IdOrgan];
                    std::cout  << IdOrgan << " "  << OrganIDMatIDMap[IdOrgan]  << " "<< OrganIDDensityMap[IdOrgan] << std::endl;

                }
                else if (indicator == "organId_name") {

                    IdOrgan = atoi(word.c_str());
                    LineString >> OrganIDNameMap[IdOrgan];
                    OrganNameIDMap[OrganIDNameMap[IdMat]] = IdOrgan;
                    std::cout  << IdOrgan << " "  << OrganIDNameMap[IdOrgan] << std::endl;

                }
                else if (indicator == "materialId_name") {

                    IdMat = atoi(word.c_str());
                    LineString >> MaterialIDName[IdMat];
                    MatNameIDMap[MaterialIDName[IdMat]] = IdMat;

                    std::cout  << IdMat << " "  << MaterialIDName[IdMat] << std::endl;

                }
                else if (indicator == "segmentation_data") {

                    IdOrgan = atoi(word.c_str());
                    TextToWrite += "/GeometryData/setVoxelizedRegionData " + OrganIDNameMap[IdOrgan] + " ";
                    LineString >> word; TextToWrite += std::to_string(atoi(word.c_str())-1) + " ";
                    LineString >> word; TextToWrite += std::to_string(atoi(word.c_str())-1) + " ";
                    LineString >> word; TextToWrite += std::to_string(atoi(word.c_str())-1) + " ";
                    LineString >> word; TextToWrite += std::to_string(atoi(word.c_str())-1) + " ";
                    LineString >> word; TextToWrite += std::to_string(atoi(word.c_str())-1) + " ";
                    LineString >> word; TextToWrite += std::to_string(atoi(word.c_str())-1) + " ";

                    TextToWrite += std::to_string(IdOrgan) + " ";
                    LineString >> word; TextToWrite += "null null null \n";

                    //std::cout  << TextToWrite << std::endl;

                }
                else if(indicator == "blood_mass_ratios"){

                    IdMat = atoi(word.c_str());
                    LineString >> MediaBloodFractionMap[IdMat];
                    std::cout  << IdMat << " "  << MediaBloodFractionMap[IdMat] << std::endl;
                }
                else if(indicator == "bone_mass_ratios"){

                    IdMat = atoi(word.c_str());
                    LineString >> MediaRBMFractionMap[IdMat];
                    LineString >> MediaYBMFractionMap[IdMat];
                    LineString >> MediaBoneFractionMap[IdMat];
                    std::cout  << IdMat << " "  << MediaRBMFractionMap[IdMat] << " "  << MediaYBMFractionMap[IdMat] << " "  << MediaBoneFractionMap[IdMat] << std::endl;

                }
                else if (indicator == "organName_colour") {

                    G4double a, b, c, d ;
                    IdOrgan = atoi(word.c_str());
                    LineString >> a ;
                    LineString >> b ;
                    LineString >> c ;
                    LineString >> d ;
                    OrganIDColourMap[d] = G4Colour( (G4double)G4UniformRand(), (G4double)G4UniformRand()  , (G4double)G4UniformRand()  , 1.);

                }
                else if (indicator == "voxels_matId") {
                    break;
                }
            }
        }

        fileR.close();
    }

    std::ofstream file1("MaterialsCommands.txt", std::ios_base::binary);

    if(file1.is_open()){

        G4Material* Media;

        G4double z, a, density;
        G4String name, symbol;
        G4int numberofElements;

        G4Element* elH = new G4Element( name = "Hydrogen", symbol = "H", z = 1.0, a = 1.008  * g/mole );
        file1 << "/MaterialData/createElement " << 1.0 << " " << 1.008 << " " << "Hydrogen" << "\n" ;
        G4Element* elC = new G4Element( name = "Carbon", symbol = "C", z = 6.0, a = 12.011 * g/mole );
        file1 << "/MaterialData/createElement " << 6.0 <<  " " << 12.011 << " " << "Carbon" << "\n" ;
        G4Element* elN = new G4Element( name = "Nitrogen", symbol = "N", z = 7.0, a = 14.007 * g/mole );
        file1 << "/MaterialData/createElement " << 7 <<  " " << 14.007 << " " << "Nitrogen" << "\n" ;
        G4Element* elO = new G4Element( name = "Oxygen", symbol = "O", z = 8.0, a = 16.00  * g/mole );
        file1 << "/MaterialData/createElement " << 8 <<  " " << 16 << " " << "Oxygen" << "\n" ;
        G4Element* elNa = new G4Element( name = "Sodium", symbol = "Na", z= 11.0, a = 22.98977* g/mole );
        file1 << "/MaterialData/createElement " << 11 <<  " " << 12.98977 << " " << "Sodium" << "\n" ;
        G4Element* elMg = new G4Element( name = "Magnesium", symbol = "Mg", z = 12.0, a = 24.3050* g/mole );
        file1 << "/MaterialData/createElement " << 12 <<  " " << 24.3050<< " " << "Magnesium" << "\n" ;
        G4Element* elP = new G4Element( name = "Phosphorus", symbol = "P", z = 15.0, a = 30.973976* g/mole );
        file1 << "/MaterialData/createElement " << 12 <<  " " << 30.973976 << " " << "Phosphorus" << "\n" ;
        G4Element* elS = new G4Element( name = "Sulfur", symbol = "S", z = 16.0,a = 32.065* g/mole );
        file1 << "/MaterialData/createElement " << 16 <<  " " << 32.065 << " " << "Sulfur" << "\n" ;
        G4Element* elCl = new G4Element( name = "Chlorine", symbol = "P", z = 17.0, a = 35.453* g/mole );
        file1 << "/MaterialData/createElement " << 17 <<  " " << 35.453 << " " << "Chlorine" << "\n" ;
        G4Element* elK = new G4Element( name = "Potassium", symbol = "P", z = 19.0, a = 39.0983* g/mole );
        file1 << "/MaterialData/createElement " << 19 <<  " " << 39.0983 << " " << "Potassium" << "\n" ;
        G4Element* elCa = new G4Element( name="Calcium", symbol = "Ca", z = 20.0, a = 40.078* g/mole );
        file1 << "/MaterialData/createElement " << 20 <<  " " << 40.078 << " " << "Calcium" << "\n" ;
        G4Element* elFe = new G4Element( name = "Iron", symbol = "Fe", z = 26, a = 56.845* g/mole );
        file1 << "/MaterialData/createElement " << 26 <<  " " << 56.845 << " " << "Iron" << "\n" ;
        G4Element* elI = new G4Element( name = "Iode", symbol = "I", z = 53, a = 126.90447* g/mole );
        file1 << "/MaterialData/createElement " << 53 <<  " " << 126.90447 << " " << "Iode" << "\n" ;

        file1 << "\n" << "\n" ;

        //ICRPOrgansNumber = 141;
        std::cout << "\n\n" << "ICRPOrgansNumber : " << OrganIDMatIDMap.size() <<  std::endl;

        for(int d = 0; d < OrganIDMatIDMap.size() ;d++ ){

            if(ICRPMaterialsCompFra[OrganIDMatIDMap[d]].size() < 13){
                continue;
            }
            //Media = new G4Material( OrganIDNameMap[d], density = OrganIDDensityMap[d]*g/cm3, numberofElements = 13);
            file1 << "\n/MaterialData/createMaterial " << OrganIDNameMap[d] << " " << d << " " << 13 << " " << OrganIDDensityMap[d] << " g/cm3 frac" << " \n" ;
            std::cout << "\n/MaterialData/createMaterial " << OrganIDNameMap[d] << " " << d << " " << 13 << " " << OrganIDDensityMap[d] << " g/cm3 frac" << " \n" ;

            //Media->AddElement(elH,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][0]/100);
            file1 << "/MaterialData/addElements " << elH->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][0] << " " ;
            //Media->AddElement(elC,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][1]/100);
            file1 << elC->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][1] << " " ;
            //Media->AddElement(elN,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][2]/100);
            file1 << elN->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][2] << " " ;
            //Media->AddElement(elO,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][3]/100);
            file1 << elO->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][3] << " " ;
            //Media->AddElement(elNa,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][4]/100);
            file1 << elNa->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][4] << " " ;
            //Media->AddElement(elMg,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][5]/100);
            file1 << elMg->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][5] << " " ;
            //Media->AddElement(elP,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][6]/100);
            file1 << elP->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][6] << " " ;
            //Media->AddElement(elS,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][7]/100);
            file1 << elS->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][7] << " " ;
            //Media->AddElement(elCl,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][8]/100);
            file1 << elCl->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][8] << " " ;
            //Media->AddElement(elK,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][9]/100);
            file1 << elK->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][9] << " " ;
            //Media->AddElement(elCa,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][10]/100);
            file1 << elCa->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][10] << " " ;
            //Media->AddElement(elFe,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][11]/100);
            file1 << elFe->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][11] << " " ;
            //Media->AddElement(elI,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][12]/100);
            file1 << elI->GetName() << " " << ICRPMaterialsCompFra[OrganIDMatIDMap[d]][12] << "\n" ;


            std::cout << "ID: " << d << " -Name: " << OrganIDNameMap[d] << " -Density: " << OrganIDDensityMap[d] << " g/cm3 "<<  std::endl;

            //VoxelsMaterials.push_back(Media);
            //MaterialColour.push_back(G4Colour( (G4double)G4UniformRand(), (G4double)G4UniformRand()  , (G4double)G4UniformRand()  , 1. ));
            //OrganNamesVector.push_back(OrganIDNameMap[d]);
        }
        file1.close();
    }
}

void G4TVolumeConstruction::BuildICRPMaterials(){

    std::cout << "\n ----------------------- Build organs material ----------------------- \n" <<  std::endl;

    G4Material* Media;
    G4Material* Blood;

    G4double z, a, density;
    G4String name, symbol;
    G4int numberofElements;

    G4Element* elH = new G4Element( name = "Hydrogen", symbol = "H", z = 1.0, a = 1.008  * g/mole );
    G4Element* elC = new G4Element( name = "Carbon", symbol = "C", z = 6.0, a = 12.011 * g/mole );
    G4Element* elN = new G4Element( name = "Nitrogen", symbol = "N", z = 7.0, a = 14.007 * g/mole );
    G4Element* elO = new G4Element( name = "Oxygen", symbol = "O", z = 8.0, a = 16.00  * g/mole );
    G4Element* elNa = new G4Element( name = "Sodium", symbol = "Na", z= 11.0, a = 22.98977* g/mole );
    G4Element* elMg = new G4Element( name = "Magnesium", symbol = "Mg", z = 12.0, a = 24.3050* g/mole );
    G4Element* elP = new G4Element( name = "Phosphorus", symbol = "P", z = 15.0, a = 30.973976* g/mole );
    G4Element* elS = new G4Element( name = "Sulfur", symbol = "S", z = 16.0,a = 32.065* g/mole );
    G4Element* elCl = new G4Element( name = "Chlorine", symbol = "P", z = 17.0, a = 35.453* g/mole );
    G4Element* elK = new G4Element( name = "Potassium", symbol = "P", z = 19.0, a = 39.0983* g/mole );
    G4Element* elCa = new G4Element( name="Calcium", symbol = "Ca", z = 20.0, a = 40.078* g/mole );
    G4Element* elFe = new G4Element( name = "Iron", symbol = "Fe", z = 26, a = 56.845* g/mole );
    G4Element* elI = new G4Element( name = "Iode", symbol = "I", z = 53, a = 126.90447* g/mole );

    for(int d = 0; d < ICRPOrgansNumber ;d++ ){

        Media = new G4Material( OrganIDNameMap[d], density = OrganIDDensityMap[d]*g/cm3, numberofElements = 13);
        Media->AddElement(elH,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][0]/100);
        Media->AddElement(elC,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][1]/100);
        Media->AddElement(elN,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][2]/100);
        Media->AddElement(elO,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][3]/100);
        Media->AddElement(elNa,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][4]/100);
        Media->AddElement(elMg,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][5]/100);
        Media->AddElement(elP,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][6]/100);
        Media->AddElement(elS,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][7]/100);
        Media->AddElement(elCl,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][8]/100);
        Media->AddElement(elK,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][9]/100);
        Media->AddElement(elCa,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][10]/100);
        Media->AddElement(elFe,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][11]/100);
        Media->AddElement(elI,ICRPMaterialsCompFra[OrganIDMatIDMap[d]][12]/100);

        std::cout << "ID: " << d << " -Name: " << OrganIDNameMap[d] << " -Density: " << Media->GetDensity()*cm3/g << " g/cm3 "<<  std::endl;

        VoxelsMaterials.push_back(Media);
        MaterialColour.push_back(G4Colour( (G4double)G4UniformRand(), (G4double)G4UniformRand()  , (G4double)G4UniformRand()  , 1. ));
        OrganNamesVector.push_back(OrganIDNameMap[d]);
    }


}
void G4TVolumeConstruction::CreateICRPPhantomData(){

    std::cout << "\n ----------------------- Read VOXELS data file ----------------------- \n" <<  std::endl;

    for(int d = 0; d < ICRPOrgansNumber ;d++ ){

        if(OrganIDColourMap[d] == NULL || OrganIDColourMap[d] == 0){
            MaterialColour.push_back(G4Colour( (G4double)G4UniformRand(), (G4double)G4UniformRand()  , (G4double)G4UniformRand()  , 1. ));
        }
        else {
            MaterialColour.push_back(OrganIDColourMap[d]);
        }
        OrganNamesVector.push_back(OrganIDNameMap[d]);
    }

    std::ifstream fileR(GeometryPath.c_str());
    G4String line, indicator, word;
    if(fileR.is_open()){

        //G4cout  << " \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << G4endl;
        G4cout << "Reading "<< GeometryPath << " ... " << G4endl ;

        while (getline(fileR, line)) {

            std::istringstream LineString(line);

            if(LineString.str().empty()){ continue; }

            LineString >> word;  // the white spaces more than 1 are considered as 1

            //G4cout << " the word " << word << G4endl ;
            if(word.empty() || word == "" || word == "#"){ // || word.isNull()
                //G4cout <<  "word " << word << G4endl;
                continue; }
            if (word == ">>"){
                LineString >> indicator;
                G4cout <<  " --- data indicator : " << indicator << G4endl;
                continue; }
            else{
                if (indicator == "voxels_matId") {
                    break;
                }
            }
        }

        G4cout <<  " Reading voxels MatIDs..."<< G4endl;

        //VoxelVolume = 8*(VoxXHalfSize*VoxYHalfSize*VoxZHalfSize); // in mm3

        G4float Voxel0PosX = -((VoxXNumber*VoxXHalfSize) - VoxXHalfSize);
        G4float Voxel0PosY = -((VoxYNumber*VoxYHalfSize) - VoxYHalfSize);
        G4float Voxel0PosZ = -((VoxZNumber*VoxZHalfSize) - VoxZHalfSize);

        G4ThreeVector CopyNum0pos = G4ThreeVector( Voxel0PosX , Voxel0PosY , Voxel0PosZ );
        //G4cout << "CopyNum0pos "<< CopyNum0pos << G4endl;

        G4ThreeVector CopyNumPos;

        unsigned int Cn_inc = 0;
        size_t VoxelNumberInPhantom = VoxXNumber*VoxYNumber*VoxZNumber;
        //G4cout << "Total Number of Voxels to Read : " << VoxelNumberInPhantom  <<  G4endl;

        MateIDs = new size_t[VoxelNumberInPhantom];
        int ID; //, X,Y,Z, XInc = 0,YInc = 0 ,ZInc = 0 ;

        for(int f = 0; f < VoxZNumber ;f++ ){
            for(int g = 0; g < VoxYNumber ;g++ ){
                for(int d = 0; d < VoxXNumber ;d++ ){

                    fileR >> ID ;

                    VoxRegionName = OrganIDNameMap[ID];
                    G4double CNposX = CopyNum0pos.getX()+d*2*VoxXHalfSize;
                    G4double CNposY = CopyNum0pos.getY()+g*2*VoxYHalfSize;
                    G4double CNposZ = CopyNum0pos.getZ()+f*2*VoxZHalfSize;
                    //G4cout << "CNposX "<< CNposX << G4endl;

                    MateIDs[Cn_inc] = ID;
                    CopyNumberMassSize[Cn_inc] = VoxelsMaterials[ID]->GetDensity() * G4Density_to_kgPerMm3 * VoxelVolume; /*density ; e+21 in g/mm3 , and e+18 g/cm3 )*/
                    CopyNumberRegionNameMap[Cn_inc] = VoxRegionName;
                    OrganNameMassMap[VoxRegionName] += CopyNumberMassSize[Cn_inc];

                    //CNZSize[Cn_inc]=f;
                    //CNYSize[Cn_inc]=g;
                    //CNXSize[Cn_inc]=d;
                    //std::cout << Cn_inc << " " << f << " " << g << " " << d << " MateID : " << ID  << " Organ name " << VoxRegionName << " Mass(kg) : " << CopyNumberMassSize[Cn_inc] << " RegionName : " << CopyNumberRegionNameMap[Cn_inc] <<  std::endl;

                    //std::ostringstream text;
                    //text <<"\r"<< " " << Cn_inc <<"/"<<VoxelNumberInPhantom << " ";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
                    //printf(text.str().c_str());

                    Cn_inc++;
                }
            }
        }

        fileR.close();
    }

}
void G4TVolumeConstruction::CreateICRPSAFsReferenceDataFile(){

    std::vector<G4double> EneVal;
    std::map<G4String, std::map<G4String, std::map<G4double, G4double>>> RefData;

    G4String FileName = "/home/tarik/DoseCalcs_build/core_build/rcp-am_alpha_2016-08-12.SAF";
    G4String ParticleName = "alpha";
    G4String GeomName = "ICRPAdultMale";
    G4int numE = 24;

    for(int na = 0; na < 8 ;na++ ){

        EneVal.clear();
        if(na == 0){
            FileName = "/home/tarik/DoseCalcs_build/core_build/rcp-af_photon_2016-08-12.SAF";
            ParticleName = "gamma";
            GeomName = "ICRPAdultMale";
            numE = 28;
        }
        else if(na == 1){
            FileName = "/home/tarik/DoseCalcs_build/core_build/rcp-am_electron_2016-08-12.SAF";
            ParticleName = "e-";
            GeomName = "ICRPAdultMale";
            numE = 28;
        }
        else if(na == 2){
            FileName = "/home/tarik/DoseCalcs_build/core_build/rcp-am_alpha_2016-08-12.SAF";
            ParticleName = "alpha";
            GeomName = "ICRPAdultMale";
            numE = 24;
        }
        else if(na == 3){
            FileName = "/home/tarik/DoseCalcs_build/core_build/rcp-af_photon_2016-08-12.SAF";
            ParticleName = "gamma";
            GeomName = "ICRPAdultFemale";
            numE = 28;
        }
        else if(na == 4){
            FileName = "/home/tarik/DoseCalcs_build/core_build/rcp-af_electron_2016-08-12.SAF";
            ParticleName = "e-";
            GeomName = "ICRPAdultFemale";
            numE = 28;
        }
        else if(na == 5){
            FileName = "/home/tarik/DoseCalcs_build/core_build/rcp-af_alpha_2016-08-12.SAF";
            ParticleName = "alpha";
            GeomName = "ICRPAdultFemale";
            numE = 24;
        }
        else if(na == 6){
            FileName = "/home/tarik/DoseCalcs_build/core_build/rcp-am_neutron_2016-08-12.SAF";
            ParticleName = "neutron";
            GeomName = "ICRPAdultMale";
            numE = 28;
        }
        else if(na == 7){
            FileName = "/home/tarik/DoseCalcs_build/core_build/rcp-af_neutron_2016-08-12.SAF";
            ParticleName = "neutron";
            GeomName = "ICRPAdultFemale";
            numE = 28;
        }
        std::cout << "\n ----------------------- Read ICRP Reference data file ----------------------- \n" <<  std::endl;

        std::ifstream fileR(FileName.c_str());
        G4String line, word, Src, Trg; G4double val, massfactor;

        if(fileR.is_open()){

            G4cout << "Reading "<< FileName << " ... " << G4endl ;

            G4int kk = 0;
            while (getline(fileR, line)) {

                std::istringstream LineString(line);

                //G4cout << "line : "<< LineString.str() << G4endl ;

                if(LineString.str().empty()){
                    continue;
                }
                if(kk == 0){

                    for(int f = 0; f < numE ;f++ ){
                        LineString >> val;
                        EneVal.push_back(val);
                    }
                    //G4cout << "EneVal.size() : "<< EneVal.size() << G4endl ;

                    kk++;
                    continue;
                }

                line.replace(10, 2, " ");
                //G4cout << "line without <- : "<< line << G4endl ;

                std::istringstream LineString11(line);
                LineString11 >> Trg >> Src;

                for(int f = 0; f < numE ;f++ ){

                    LineString11 >> val;
                    double valuee = val;
/*
                    if(GeomName == "ICRPAdultFemale"){
                        if(Trg == "Liver"){valuee = val*(1400+410)/1400;}
                        else if(Trg == "Brain"){valuee = val*(1300+49.2)/1300;}
                        else if(Trg == "Thymus"){valuee = val*(20+0.000615)/20;}
                        else if(Trg == "Thyroid"){valuee = val*(17+2.46)/17;}
                        else if(Trg == "Ht-wall"){valuee = val*(250+41)/250;}
                        else if(Trg == "Pancreas"){valuee = val*(120+24.6)/120;}
                        else if(Trg == "Spleen"){valuee = val*(130+57.4)/130;}
                    }else if(GeomName == "ICRPAdultMale"){
                        if(Trg == "Liver"){valuee = val*(1800+560)/1800;}
                        else if(Trg == "Brain"){valuee = val*(1450+67.2)/1450;}
                        else if(Trg == "Thymus"){valuee = val*(25+1.17)/25;}
                        else if(Trg == "Thyroid"){valuee = val*(20+3.36)/20;}
                        else if(Trg == "Ht-wall"){valuee = val*(330+56)/330;}
                        else if(Trg == "Pancreas"){valuee = val*(140+33.6)/140;}
                        else if(Trg == "Spleen"){valuee = val*(150+78.4)/150;}
                        else if(Trg == "Prostate"){valuee = val*(17+0.000797)/17;}
                    }
*/
                    // for organs name defined in our results
                    //if     (Trg == "UB-cont"){Trg = "UBCs";}
                    //else if(Src == "UB-cont"){Src = "UBCs";}
                    //else if(Trg == "Ht-wall"){Trg = "HeW";}
                    //else if(Src == "Ht-wall"){Src = "HeW";}

                    RefData[Src][Trg][EneVal[f]] = valuee;

                    if     (Src == "Liver" && Trg == "Liver"){
                           G4cout << " " << Src << " " << Trg << " " << EneVal[f] << " " << valuee << G4endl ;
                    }
                }
            }

            fileR.close();
        }

        std::ofstream fileO("ICRP133SAFs", std::ios_base::app);
        G4cout << "Writing to FemaleICRPRefData.txt ..." << G4endl ;
        if(fileO.is_open()){

            for ( auto Abeg = RefData.begin(); Abeg != RefData.end(); ++Abeg  )
            {
                Src = Abeg->first;
                fileO << "****** SAF kg-1 " << Src << " " << ParticleName << " " << GeomName << " MeV ";

                for(int f = 0; f < numE ;f++ ){
                    //fileO << std::setw(11) << std::left << EneVal[f] << " ";
                    fileO << EneVal[f] << " ";
                }
                fileO << "\n" ;

                for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  )
                {
                    Trg = Bbeg->first;

                    fileO << std::setw(30) << std::left << Trg << " " ;

                    massfactor = 1. ;
                    //for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg ){
                    for(int f = 0; f < numE ;f++ ){
                        fileO << RefData[Src][Trg][EneVal[f]] << " " ;
                    }
                    fileO << "\n" ;

                }
                fileO << "* ---------------------------------------------------------------------\n" ;
            }
            fileO.close();
        }

    }
}

//############################################################ setting or generating source data
void G4TVolumeConstruction::setRankDataForMPIMode(){

    G4int SourceParticlesNumber = SourceParticlesNamesValues.size();
    G4int SourceRegionsNumber   = SourceRegionsNamesValues.size();
    G4int SourceEnergiesNumber  = SourceEnergiesValues.size();
    G4int SourceMomDirsNumber   = SourceMomDirsValues.size();

    G4int NumberOfSimulation = SourceParticlesNumber*SourceRegionsNumber*SourceEnergiesNumber*SourceMomDirsNumber;

    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
    //          " SourceParticlesNumber " << SourceParticlesNumber << " SourceRegionsNumber " << SourceRegionsNumber
    //       << " SourceEnergiesNumber " << SourceEnergiesNumber << " SourceMomDirsNumber " << SourceMomDirsNumber << G4endl ;


    if(NumberOfSimulation==0){
        G4String msg = "Particle, region, energy or momentum direction source data are not set"; G4Exception("Source Data", "1", FatalErrorInArgument, msg.c_str());
    }

    G4int hh = 0;
    G4int a1 = 0;
    G4int b1 = 0;
    G4int c1 = 0;
    G4int d1 = 0;

    for( G4int a = 0; a < SourceParticlesNumber ; a++ ){
        for( G4int b = 0; b < SourceRegionsNumber ; b++ ){
            for( G4int c = 0; c < SourceEnergiesNumber ; c++ ){
                for( G4int d = 0; d < SourceMomDirsNumber ; d++ ){

                    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
                    //          " SourceParticlesNumber " << SourceParticlesNumber << " SourceRegionsNumber " << SourceRegionsNumber
                    //       << " SourceEnergiesNumber " << SourceEnergiesNumber << " SourceMomDirsNumber " << SourceMomDirsNumber << G4endl ;


                    if(MPISimulationNum == "m"){

                        a1 = a;
                        b1 = b;
                        c1 = c;
                        d1 = d;
                    }

                    std::vector <unsigned int> VoxelsIDsOfSourceRegion ;
                    NewRankSourceParticlesNamesValues.push_back(SourceParticlesNamesValues[a1]);
                    NewRankSourceRegionsNamesValues.push_back(SourceRegionsNamesValues[b1]);
                    if(SourceType == "Volume"){
                        NewRankSourceRegionsBoxDimValues.push_back(SourceRegionsBoxDimValues[b1]);
                        NewRankSourceRegionsPosValues.push_back(getRegionAbsolutePosition(SourceRegionsNamesValues[b1]));
                        G4Navigator* aNavigator = new G4Navigator();
                        aNavigator->SetWorldVolume(WorldPhysicalVolume);
                        NavigatorForVolumesInitialPosition.push_back(aNavigator);
                        //G4cout <<" \n\n\n\n\n\n\n\n\n\n BoxDimValues " << SourceRegionsBoxDimValues[b1]
                        //        << " NewRankSourceRegionsPosValues " << getRegionAbsolutePosition(SourceRegionsNamesValues[b1])
                        //        << " SourceRegionsNamesValues " << SourceRegionsNamesValues[b1] << G4endl;

                    }else if(SourceType == "Voxels"){
                        G4int TotalVoxelsNumber = VoxXNumber*VoxYNumber*VoxZNumber;
                        for (int n = 0; n < TotalVoxelsNumber; n++) {
                            //G4cout << n <<"/" << TotalVoxelsNumber << " " << CopyNumberRegionNameMap[n] << "***" << SourceRegionName << CopyNumberXPos[n] << " " << CopyNumberYPos[n]  <<  " " << CopyNumberZPos[n] << G4endl;

                            if(SourceRegionsNamesValues[b1] == "allregions"){

                                bool SaveCN = true;
                                for (int gg = 0 ; gg < SourceRegionsNamesToBeIgnoredValues.size() ; gg++) {
                                    if(SourceRegionsNamesToBeIgnoredValues[gg] == CopyNumberRegionNameMap[n]){
                                        SaveCN=false;
                                    }
                                }
                                if(SaveCN){
                                    //std::cout << n << " " << TotalVoxelsNumber << " " << CopyNumberRegionNameMap[n] << " X=" << CopyNumberXPos[n] << " Y=" << CopyNumberYPos[n]  <<  " Z=" << CopyNumberZPos[n] << std::endl;
                                    VoxelsIDsOfSourceRegion.push_back(n);
                                }

                            }else{

                                if(CopyNumberRegionNameMap[n] == SourceRegionsNamesValues[b1]){
                                    //G4cout << "\n\n\n\n-------------------------\n\n\n " << SourceRegionName << G4endl;
                                    VoxelsIDsOfSourceRegion.push_back(n);
                                }
                            }
                        }
                        NewRankVoxelsIDsOfSourceRegion.push_back(VoxelsIDsOfSourceRegion);
                    }else if(SourceType == "TET"){

                        std::vector<G4Tet*>      internalTetVec;
                        G4double xMin(DBL_MAX), yMin(DBL_MAX), zMin(DBL_MAX);
                        G4double xMax(DBL_MIN), yMax(DBL_MIN), zMax(DBL_MIN);

                        for(G4int n=0;n<tetData->GetNumTetrahedron();n++){

                            if(SourceRegionsNamesValues[b1] == "allregions"){

                                bool SaveCN = true;
                                for (int gg = 0 ; gg < SourceRegionsNamesToBeIgnoredValues.size() ; gg++) {
                                    if(SourceRegionsNamesToBeIgnoredValues[gg] == CopyNumberRegionNameMap[n]){
                                        SaveCN=false;
                                    }
                                }
                                if(SaveCN){
                                    G4Tet* tetSolid = tetData->GetTetrahedron(n);
                                    for(auto vertex:tetSolid->GetVertices()){
                                        if      (vertex.getX() < xMin) xMin = vertex.getX();
                                        else if (vertex.getX() > xMax) xMax = vertex.getX();
                                        if      (vertex.getY() < yMin) yMin = vertex.getY();
                                        else if (vertex.getY() > yMax) yMax = vertex.getY();
                                        if      (vertex.getZ() < zMin) zMin = vertex.getZ();
                                        else if (vertex.getZ() > zMax) zMax = vertex.getZ();
                                    }
                                    internalTetVec.push_back(tetData->GetTetrahedron(n));
                                }
                            }else{

                                if(CopyNumberRegionNameMap[n] != SourceRegionsNamesValues[b1]) continue;
                                //VoxelsIDsOfSourceRegion.push_back(n);
                                G4Tet* tetSolid = tetData->GetTetrahedron(n);
                                for(auto vertex:tetSolid->GetVertices()){
                                    if      (vertex.getX() < xMin) xMin = vertex.getX();
                                    else if (vertex.getX() > xMax) xMax = vertex.getX();
                                    if      (vertex.getY() < yMin) yMin = vertex.getY();
                                    else if (vertex.getY() > yMax) yMax = vertex.getY();
                                    if      (vertex.getZ() < zMin) zMin = vertex.getZ();
                                    else if (vertex.getZ() > zMax) zMax = vertex.getZ();
                                }
                                internalTetVec.push_back(tetData->GetTetrahedron(n));
                            }
                        }

                        //NewRankVoxelsIDsOfSourceRegion.push_back(VoxelsIDsOfSourceRegion);

                        NewRankTETBoxMinOfSourceRegion.push_back(G4ThreeVector(xMin, yMin, zMin));
                        NewRankTETBoxDimOfSourceRegion.push_back(G4ThreeVector(xMax-xMin, yMax-yMin, zMax-zMin));
                        NewRankTETOfSourceRegion.push_back(internalTetVec);
                    }
                    NewRankSourceEnergiesValues.push_back(SourceEnergiesValues[c1]);
                    NewRankSourceMomDirsValues.push_back(SourceMomDirsValues[d1]);
                    NewRankSourceMomDirsDirectedThetaValues.push_back(SourceMomDirsDirectedThetaValues[d1]);
                    NewRankSourceMomDirsDirectedPhiValues.push_back(SourceMomDirsDirectedPhiValues[d1]);

                    SourceRegionName = SourceRegionsNamesValues[b1];
                    if(EnergyDistribution == "Mono"){ MonoEnergy = SourceEnergiesValues[c1];}
                    else if(EnergyDistribution == "Rayleigh"){ RayleighEmax = SourceEnergiesValues[c1];}
                    else if(EnergyDistribution == "Gauss"){  GaussMean = SourceEnergiesValues[c1];}
                    else if(EnergyDistribution == "Uniform"){ UniformEmax = SourceEnergiesValues[c1];}
                    else if(EnergyDistribution == "Spectrum"){ SpectrumMaxEnergy = SourceEnergiesValues[c1];}
                    else if(EnergyDistribution == "File"){ FileEnergyCharacterizer = SourceEnergiesValues[c1];}
                    else if(EnergyDistribution == "RadioNuclide"){ RadioNuclideMaxEnergy = SourceEnergiesValues[c1];}

                    MomDirDistribution = SourceMomDirsValues[d1];
                    Theta = SourceMomDirsDirectedThetaValues[d1];
                    Phi = SourceMomDirsDirectedPhiValues[d1];

                    NewRankSourcePositionDataFiles.push_back(setSourcePositionFileName());
                    NewRankSourceEnergyDataFiles.push_back(setSourceEnergyFileName());
                    NewRankSourceMomDirDataFiles.push_back(setSourceMomDirFileName());

                    //G4cout << hh << " - NewRankSourceRegionsNamesValues[a]: " << NewRankSourceRegionsNamesValues[b] <<  G4endl;
                    //G4cout << hh << " - NewRankSourceRegionsBoxDimValues[a]: " << NewRankSourceRegionsBoxDimValues[b] <<  G4endl;
                    //G4cout << hh << " - NewRankSourceEnergiesValues[a]: " << NewRankSourceEnergiesValues[c] <<  G4endl;
                    //G4cout << hh << " - NewRankSourceMomDirsValues[a]: " << NewRankSourceMomDirsValues[d] <<  G4endl;
                    hh++;
                }
            }
        }
    }

    G4int SourceDataInc = 0 , MaxSourceDataSim = 1;
    MaxSourceDataSim = NewRankSourceParticlesNamesValues.size(); // for multi simulation the source data inc should be equal to the rank ID

#ifdef G4MPI_USE

    int NumberOfRanksThreads; // replace

    MPI_Comm_size(MPI_COMM_WORLD, &NumberOfRanksThreads);
    //NumberOfRanksThreads = G4MPImanager::GetManager()->GetTotalSize();

    MaxSourceDataSim = NumberOfRanksThreads;
    if(NumberOfRanksThreads > NewRankSourceParticlesNamesValues.size()){
        G4int nn = NumberOfRanksThreads - NewRankSourceParticlesNamesValues.size();
#if VERBOSE_USE
        G4cout << "\n\n !!!!!!!!!!!!!!! The number of ranks "<< NumberOfRanksThreads <<", exceeds the source configuration by " << nn << " rank \n ===> The first source data (of rank 0) will be resimulated by the additional ranks"<< G4endl;
#endif
        for( G4int d = 0; d < nn ; d++ ){
            NewRankSourceParticlesNamesValues.push_back(NewRankSourceParticlesNamesValues[0]);
            NewRankSourceRegionsNamesValues.push_back(NewRankSourceRegionsNamesValues[0]);
            if(SourceType == "Volume"){
                NewRankSourceRegionsBoxDimValues.push_back(NewRankSourceRegionsBoxDimValues[0]);
                NewRankSourceRegionsPosValues.push_back(NewRankSourceRegionsPosValues[0]);
                G4Navigator* aNavigator = new G4Navigator();
                aNavigator->SetWorldVolume(WorldPhysicalVolume);
                NavigatorForVolumesInitialPosition.push_back(aNavigator);
            }else if(SourceType == "Voxels"){
                NewRankVoxelsIDsOfSourceRegion.push_back(NewRankVoxelsIDsOfSourceRegion[0]);
            }
            NewRankSourceEnergiesValues.push_back(NewRankSourceEnergiesValues[0]);
            NewRankSourceMomDirsValues.push_back(NewRankSourceMomDirsValues[0]);
            NewRankSourceMomDirsDirectedThetaValues.push_back(NewRankSourceMomDirsDirectedThetaValues[0]);
            NewRankSourceMomDirsDirectedPhiValues.push_back(NewRankSourceMomDirsDirectedPhiValues[0]);

            NewRankSourcePositionDataFiles.push_back(NewRankSourcePositionDataFiles[0]);
            NewRankSourceEnergyDataFiles.push_back(NewRankSourceEnergyDataFiles[0]);
            NewRankSourceMomDirDataFiles.push_back(NewRankSourceMomDirDataFiles[0]);

        }

    }else if(NumberOfRanksThreads < NewRankSourceParticlesNamesValues.size()){
        G4int nn = NewRankSourceParticlesNamesValues.size() - NumberOfRanksThreads;
#if VERBOSE_USE
        G4cout << "\n\n !!!!!!!!!!!!!!! " << nn << " source configuration will not be simulated all because of the not sufficient number of ranks "<< G4endl;
#endif
    }
    /*
    std::cout << " ------- NewRankSourceParticlesNamesValues[0] " << NewRankSourceParticlesNamesValues[0] << std::endl;
    std::cout << " ------- NewRankSourceRegionsNamesValues[0] " << NewRankSourceRegionsNamesValues[0] << std::endl;
    std::cout << " ------- NewRankSourceRegionsBoxDimValues[0] " << NewRankSourceRegionsBoxDimValues[0] << std::endl;
    std::cout << " ------- NewRankSourceRegionsPosValues[0] " << NewRankSourceRegionsPosValues[0] << std::endl;
    std::cout << " ------- NewRankSourceEnergiesValues[0] " << NewRankSourceEnergiesValues[0] << std::endl;
    std::cout << " ------- NewRankSourceMomDirsValues[0] " << NewRankSourceMomDirsValues[0] << std::endl;
*/
#else
    MaxSourceDataSim = ThreadsNumber;
    if(ThreadsNumber > NewRankSourceParticlesNamesValues.size()){
        G4int nn = ThreadsNumber - NewRankSourceParticlesNamesValues.size();
#if VERBOSE_USE
        G4cout << "\n\n !!!!!!!!!!!!!!! The number of Threads "<< ThreadsNumber <<", exceeds the source configuration number " << NewRankSourceParticlesNamesValues.size() << " by " << nn << " thread \n ===> The first source data (of thread 0) will be resimulated by the additional threads"<< G4endl;
#endif
        for( G4int d = 0; d < nn ; d++ ){
            NewRankSourceParticlesNamesValues.push_back(NewRankSourceParticlesNamesValues[0]);
            NewRankSourceRegionsNamesValues.push_back(NewRankSourceRegionsNamesValues[0]);
            if(SourceType == "Volume"){
                NewRankSourceRegionsBoxDimValues.push_back(NewRankSourceRegionsBoxDimValues[0]);
                NewRankSourceRegionsPosValues.push_back(NewRankSourceRegionsPosValues[0]);
                G4Navigator* aNavigator = new G4Navigator();
                aNavigator->SetWorldVolume(WorldPhysicalVolume);
                NavigatorForVolumesInitialPosition.push_back(aNavigator);
            }else if(SourceType == "Voxels"){
                NewRankVoxelsIDsOfSourceRegion.push_back(NewRankVoxelsIDsOfSourceRegion[0]);
            }else if(SourceType == "TET"){  // is the same as we fill the array first, look for the statement
                                            // "else if(SourceType == "TET"){" before this

                std::vector<G4Tet*>      internalTetVec;
                G4double xMin(DBL_MAX), yMin(DBL_MAX), zMin(DBL_MAX);
                G4double xMax(DBL_MIN), yMax(DBL_MIN), zMax(DBL_MIN);

                for(G4int n=0;n<tetData->GetNumTetrahedron();n++){

                    if(SourceRegionsNamesValues[0] == "allregions"){

                        bool SaveCN = true;
                        for (int gg = 0 ; gg < SourceRegionsNamesToBeIgnoredValues.size() ; gg++) {
                            if(SourceRegionsNamesToBeIgnoredValues[gg] == CopyNumberRegionNameMap[n]){
                                SaveCN=false;
                            }
                        }
                        if(SaveCN){
                            G4Tet* tetSolid = tetData->GetTetrahedron(n);
                            for(auto vertex:tetSolid->GetVertices()){
                                if      (vertex.getX() < xMin) xMin = vertex.getX();
                                else if (vertex.getX() > xMax) xMax = vertex.getX();
                                if      (vertex.getY() < yMin) yMin = vertex.getY();
                                else if (vertex.getY() > yMax) yMax = vertex.getY();
                                if      (vertex.getZ() < zMin) zMin = vertex.getZ();
                                else if (vertex.getZ() > zMax) zMax = vertex.getZ();
                            }
                            internalTetVec.push_back(tetData->GetTetrahedron(n));
                        }
                    }else{

                        if(CopyNumberRegionNameMap[n] != SourceRegionsNamesValues[0]) continue;
                        //VoxelsIDsOfSourceRegion.push_back(n);
                        G4Tet* tetSolid = tetData->GetTetrahedron(n);
                        for(auto vertex:tetSolid->GetVertices()){
                            if      (vertex.getX() < xMin) xMin = vertex.getX();
                            else if (vertex.getX() > xMax) xMax = vertex.getX();
                            if      (vertex.getY() < yMin) yMin = vertex.getY();
                            else if (vertex.getY() > yMax) yMax = vertex.getY();
                            if      (vertex.getZ() < zMin) zMin = vertex.getZ();
                            else if (vertex.getZ() > zMax) zMax = vertex.getZ();
                        }
                        internalTetVec.push_back(tetData->GetTetrahedron(n));
                    }
                }

                //NewRankVoxelsIDsOfSourceRegion.push_back(VoxelsIDsOfSourceRegion);

                NewRankTETBoxMinOfSourceRegion.push_back(G4ThreeVector(xMin, yMin, zMin));
                NewRankTETBoxDimOfSourceRegion.push_back(G4ThreeVector(xMax-xMin, yMax-yMin, zMax-zMin));
                NewRankTETOfSourceRegion.push_back(internalTetVec);
            }

            NewRankSourceEnergiesValues.push_back(NewRankSourceEnergiesValues[0]);
            NewRankSourceMomDirsValues.push_back(NewRankSourceMomDirsValues[0]);
            NewRankSourceMomDirsDirectedThetaValues.push_back(NewRankSourceMomDirsDirectedThetaValues[0]);
            NewRankSourceMomDirsDirectedPhiValues.push_back(NewRankSourceMomDirsDirectedPhiValues[0]);

            NewRankSourcePositionDataFiles.push_back(NewRankSourcePositionDataFiles[0]);
            NewRankSourceEnergyDataFiles.push_back(NewRankSourceEnergyDataFiles[0]);
            NewRankSourceMomDirDataFiles.push_back(NewRankSourceMomDirDataFiles[0]);
        }
    }else if(ThreadsNumber < NewRankSourceParticlesNamesValues.size()){
        G4int nn = NewRankSourceParticlesNamesValues.size() - ThreadsNumber;
#if VERBOSE_USE
        G4cout << "\n\n !!!!!!!!!!!!!!! " << nn << " source configuration will not be simulated because of the not sufficient number of threads "<< G4endl;
#endif
    }
#endif

    for( G4int d = 0; d < MaxSourceDataSim ; d++ ){
        SourceDataInc = d;
        if(MPISimulationNum == "o"){
            SourceDataInc = 0;
        }
#ifdef G4MPI_USE
        RanksSourceData << "Rank                             " << d << " " << NewRankSourceParticlesNamesValues[SourceDataInc]<< " " << NewRankSourceRegionsNamesValues[SourceDataInc]<< " " << NewRankSourceEnergiesValues[SourceDataInc] << "\n" ;
#else
        RanksSourceData << "Thread                           " << d << " " << NewRankSourceParticlesNamesValues[SourceDataInc]<< " " << NewRankSourceRegionsNamesValues[SourceDataInc]<< " " << NewRankSourceEnergiesValues[SourceDataInc] << "\n" ;
#endif
    }


    InitializeRadiationSource();

}
void G4TVolumeConstruction::GenerateEventsDataInMPIMode(){

    int NumberOfRanksThreads;
    int DataID ;

    std::vector<G4String> Rank_GenerateType;

    G4int SourceRegionsNumber = SourceRegionsNamesValues.size();
    G4int SourceEnergiesNumber = SourceEnergiesValues.size();
    G4int SourceMomDirsNumber = SourceMomDirsValues.size();

    G4int NumberOfRanks = SourceRegionsNumber+SourceEnergiesNumber+SourceMomDirsNumber;

    G4cout << "\n\n\n\n\n\n\n\nSourceRegionsNumber: " << SourceRegionsNumber << " - SourceEnergiesNumber: " << SourceEnergiesNumber << " -SourceMomDirsNumber: " << SourceMomDirsNumber << " -NumberOfRanks: " << NumberOfRanks <<  G4endl;

    G4int RankInc;
    for(RankInc = 0 ; RankInc < SourceRegionsNumber ; RankInc++ ){
        Rank_GenerateType.push_back("Position");
    }

    G4int d = RankInc;
    for( RankInc = d; RankInc < SourceRegionsNumber+SourceEnergiesNumber ; RankInc++ ){
        Rank_GenerateType.push_back("Energy");
    }

    d = RankInc;
    for( RankInc = d; RankInc < SourceRegionsNumber+SourceEnergiesNumber+SourceMomDirsNumber ; RankInc++ ){
        Rank_GenerateType.push_back("MomDir");
    }

    for(G4int ee = 0 ; ee < SourceRegionsNamesValues.size() ; ee++ ){
        //G4cout << "RegionsNames: " << SourceRegionsNamesValues[ee] << " " << SourceRegionsBoxDimValues[ee] <<   G4endl;
    }
    for(G4int ee = 0 ; ee < SourceEnergiesValues.size() ; ee++ ){
        //G4cout << "Energies    : " << SourceEnergiesValues[ee] <<  G4endl;
    }
    for(G4int ee = 0 ; ee < SourceMomDirsValues.size() ; ee++ ){
        //G4cout << "Mom Dirs    : " << SourceMomDirsValues[ee] <<  G4endl;
    }

#ifdef G4MPI_USE

    MPI_Comm_rank(MPI_COMM_WORLD, &DataID);
    MPI_Comm_size(MPI_COMM_WORLD, &NumberOfRanksThreads);

    //NumberOfRanksThreads = G4MPImanager::GetManager()->GetTotalSize();
    //DataID = G4MPImanager::GetManager()->GetRank();


    std::vector<G4String> NN = Rank_GenerateType; Rank_GenerateType.clear();
    for(G4int ee = 0 ; ee < NN.size() ; ee++ ){
        if( ee < NumberOfRanksThreads){
            Rank_GenerateType.push_back(NN[ee]);
        }else {
            if(DataID == 0){
                G4cout << "The Not Generated Data According to Rank number is : " << NumberOfRanksThreads << " Data: " << ee << " " << Rank_GenerateType[ee] <<  G4endl;
            }
        }
    }
#endif

#ifdef G4MPI_USE
    if(DataID == 0){
        for(G4int ee = 0 ; ee < Rank_GenerateType.size() ; ee++ ){
            G4cout << "Rank ID: " << ee << " -Data To Generate: " << Rank_GenerateType[ee] <<  G4endl;
        }
    }
#else
    for(G4int ee = 0 ; ee < Rank_GenerateType.size() ; ee++ ){
        G4cout << "Rank ID: " << ee << " -Data To Generate: " << Rank_GenerateType[ee] <<  G4endl;
    }
#endif

    G4cout <<  G4endl;

    G4TPointDataGeneration* PDG = new G4TPointDataGeneration();
    PDG->SetSourceInputs();

    G4int SourceRegionInc = 0 , MomDirInc = 0 ,EnergyInc = 0 ;
    for(RankInc = 0 ; RankInc < Rank_GenerateType.size() ; RankInc++ ){

#ifdef G4MPI_USE
        if(RankInc == DataID){
            G4cout << "===== The generation data Rank: " << DataID << " which is in list index "<< RankInc << " =====" <<  G4endl;
#endif

            if(Rank_GenerateType[RankInc] == "Position" && GeneratePosFlag == "yes"){

                SourceRegionInc = RankInc;

                SourceRegionName = SourceRegionsNamesValues[SourceRegionInc];

                if(SourceType == "Volume"){
                    boxCenterPos = getRegionAbsolutePosition(SourceRegionName);
                    BoxDimGene = SourceRegionsBoxDimValues[SourceRegionInc];
                    G4cout << "SourceType "<< SourceType << " -SourceRegionName " << SourceRegionName << " -boxCenterPos " << boxCenterPos << " -BoxDimGene " << BoxDimGene <<  G4endl;
                }else {
                    G4cout << "SourceType "<< SourceType << " -SourceRegionName " << SourceRegionName <<  G4endl;
                }

                setSourcePositionFileName();
                setSourceEnergyFileName();
                setSourceMomDirFileName();
                PDG->SetSourceInputs();
                PDG->GenerateEventsPosition();
            }
            else if(Rank_GenerateType[RankInc] == "Energy" && GenerateEneFlag == "yes"){

                EnergyInc = RankInc - SourceRegionsNamesValues.size();

                if(EnergyDistribution == "Mono"){ MonoEnergy = SourceEnergiesValues[EnergyInc];
                    G4cout << "MonoEnergy " << MonoEnergy <<  G4endl;
                }
                else if(EnergyDistribution == "Rayleigh"){ RayleighEmax = SourceEnergiesValues[EnergyInc];
                    G4cout << "RayleighEmax " << RayleighEmax <<  G4endl;
                }
                else if(EnergyDistribution == "Gauss"){  GaussMean = SourceEnergiesValues[EnergyInc];
                    G4cout << "GaussMean " << GaussMean <<  G4endl;
                }
                else if(EnergyDistribution == "Uniform"){ UniformEmax = SourceEnergiesValues[EnergyInc];
                    G4cout << "UniformEmax " << UniformEmax <<  G4endl;
                }
                else if(EnergyDistribution == "Spectrum"){ SpectrumMaxEnergy = SourceEnergiesValues[EnergyInc];
                    G4cout << "SpectrumMaxEnergy " << SpectrumMaxEnergy <<  G4endl;
                }
                else if(EnergyDistribution == "File"){ FileEnergyCharacterizer = SourceEnergiesValues[EnergyInc];
                    G4cout << "FileEnergyCharacterizer " << FileEnergyCharacterizer <<  G4endl;
                }
                else if(EnergyDistribution == "RadioNuclide"){ RadioNuclideMaxEnergy = SourceEnergiesValues[EnergyInc];
                    G4cout << "RadioNuclideMaxEnergy " << RadioNuclideMaxEnergy <<  G4endl;
                }

                setSourcePositionFileName();
                setSourceEnergyFileName();
                setSourceMomDirFileName();
                PDG->SetSourceInputs();
                PDG->GenerateEventsEnergy();

            }
            else if(Rank_GenerateType[RankInc] == "MomDir" && GenerateMomDirFlag == "yes"){

                MomDirInc = RankInc - SourceRegionsNamesValues.size() - SourceEnergiesValues.size();

                MomDirDistribution = SourceMomDirsValues[MomDirInc];
                G4cout << "MomDirDistribution " << MomDirDistribution <<  G4endl;
                if(MomDirDistribution == "Directed"){
                    Theta = SourceMomDirsDirectedThetaValues[MomDirInc];
                    Phi = SourceMomDirsDirectedPhiValues[MomDirInc];
                    G4cout << "Theta " << Theta << " - Phi " << Phi <<  G4endl;
                }
                setSourcePositionFileName();
                setSourceEnergyFileName();
                setSourceMomDirFileName();
                PDG->SetSourceInputs();
                PDG->GenerateEventsMomentumDirection();
            }
#ifdef G4MPI_USE
        }
#endif

    }
}

//############################################################
void G4TVolumeConstruction::InitializeRadiationSource(){


    if(EnergyDistribution == "Spectrum"){//EnergyTypeNum = 4;

        EnergyList        = new double[EventsNumPerThreadRank];

        G4double TotalProbability = 0;

        std::map<G4double, G4double> EnergyProbabilityInterval;
        std::map<G4int, G4double> CumulatedProbabilityIndexEnergy;

        G4int index=0;

        for ( auto Abeg = EnergyValueProbability.begin(); Abeg != EnergyValueProbability.end(); ++Abeg  ){

            EnergyProbabilityInterval[Abeg->first] = Abeg->second;
            TotalProbability += EnergyProbabilityInterval[Abeg->first];

#if VERBOSE_USE
            //std::cout  << index << " TotalProbability "<< TotalProbability << " - this EnergyProbabilityInterval " << EnergyProbabilityInterval[Abeg->first] << " Energy " << std::endl;
#endif
            CumulatedProbabilityIndexEnergy[index] = Abeg->first;
            index++;
        }

        G4double* EnergyList = new G4double[EventsNumPerThreadRank];

        for(int f = 0; f < EventsNumPerThreadRank ;f++ ){

            G4double RV = TotalProbability*(G4double)G4UniformRand() ;
            G4double EnergyCumulatedProbability = 0;

            index=0;
            for ( auto Abeg = EnergyValueProbability.begin(); Abeg != EnergyValueProbability.end(); ++Abeg  ){

                EnergyCumulatedProbability += Abeg->second;

                if(RV < EnergyCumulatedProbability){
                    EnergyList[f] = CumulatedProbabilityIndexEnergy[index];
                    break;
                }
                index++;
            }

            //std::cout << f << " Energy " << EnergyValues[f] << " EnergyCumulatedProbability " << EnergyCumulatedProbability << " - this EnergyProbabilityInterval " << EnergyProbabilityInterval[EnergyValues[f]] << std::endl;
        }
    }
    else if(EnergyDistribution == "RadioNuclide" || EnergyDistribution == "File"){//EnergyTypeNum = 5;

        EnergyList        = new double[EventsNumPerThreadRank];
        ParNameList = new unsigned int[EventsNumPerThreadRank];
        //MomDirXList       = new double[EventsNumPerThreadRank];
        //MomDirYList       = new double[EventsNumPerThreadRank];
        //MomDirZList       = new double[EventsNumPerThreadRank];
        //PosXList          = new double[EventsNumPerThreadRank];
        //PosYList          = new double[EventsNumPerThreadRank];
        //PosZList          = new double[EventsNumPerThreadRank];

        G4double TotalProbability;
        for(int m = 0; m < RadioNuclideProbVec.size() ;m++ ){
            TotalProbability += RadioNuclideProbVec[m];
            //std::cout << m << " TotalProbability " << TotalProbability << " SpectrumOrDiscreteVec " << RadioNuclideSpectrumOrDiscreteVec[m] << " Proba " << RadioNuclideProbVec[m] << " EneMin " << RadioNuclideEneVec[m][0] << std::endl;
            //std::cout << m <<  " yield " << RadioNuclideProbVec[m] <<  " accumulated yield " << TotalProbability << std::endl;
        }

        for(int f = 0; f < EventsNumPerThreadRank ;f++ ){

            G4double RV = TotalProbability*(G4double)G4UniformRand() ;
            G4double EnergyCumulatedProbability = 0;

            for(int m = 0; m < RadioNuclideProbVec.size() ;m++ ){

                EnergyCumulatedProbability += RadioNuclideProbVec[m];

                if(RV < EnergyCumulatedProbability){
                    if(RadioNuclideSpectrumOrDiscreteVec[m] == 0){

                        double RandomEne = RadioNuclideEneVec[m][0] + (RadioNuclideEneVec[m][1]-RadioNuclideEneVec[m][0])*(double)G4UniformRand();
                        EnergyList[f] = RandomEne;
                        //std::cout << f << " Particle "<< RadioNuclidePartNameVec[m] << " TotalProbability "<< TotalProbability << " RandomProb " << RV << " Prob "<< RadioNuclideProbVec[m] << " EnergyCumulatedProbability " << EnergyCumulatedProbability  << " SpectrumOrDiscreteVec " << RadioNuclideSpectrumOrDiscreteVec[m] << " Min " << RadioNuclideEneVec[m][0]   << " RandomEne " << EnergyList[f] << " Max " << RadioNuclideEneVec[m][1] << std::endl;

                    }else{

                        EnergyList[f] = RadioNuclideEneVec[m][0];
                        //std::cout << f << " Particle "<< RadioNuclidePartNameVec[m] << " TotalProbability "<< TotalProbability << " RandomProb " << RV << " Prob "<< RadioNuclideProbVec[m] << " EnergyCumulatedProbability " << EnergyCumulatedProbability  << " SpectrumOrDiscreteVec " << RadioNuclideSpectrumOrDiscreteVec[m] << " selected Energy " << EnergyList[f] << std::endl;
                    }
                    ParNameList[f]=RadioNuclidePartNameVec[m];
                    break;
                }
            }
        }
    }

}


void G4TVolumeConstruction::setNumberOfThreads(G4int newNumber){
    ThreadsNumber = newNumber;
#if G4MPI_USE
    ThreadsNumber = 1;
#else
    if(ThreadsNumber == 0){
        if(G4Threading::IsMultithreadedApplication()){
            ThreadsNumber = G4Threading::G4GetNumberOfCores() ;
        }else {
            ThreadsNumber = 1;
        }
    }
#endif
    //std::cout << "\n\n\n\n ThreadsNumber "<< ThreadsNumber << std::endl;

}

// for hx, hy, hz, minx, miny, minz, maxx, maxy, maxz, we have our values and we have ICRP values
// we will find the correspondant ICRP regions minx maxx.
// 0. we have to find the coordinates of each voxel in a region
// 1. we have to find the correspondent coordinates of this voxels in our phantom,
// 2. check the density of this region voxels
