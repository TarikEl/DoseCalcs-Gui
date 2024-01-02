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

#ifndef G4TVolumeConstruction_H
#define G4TVolumeTConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4TMessenger.hh"
#include "G4TGeometryMessenger.hh"
#include "G4TPhysicsMessenger.hh"
#include "G4TPrimaryGeneratorMessenger.hh"

#include "globals.hh"
#include "G4Material.hh"
#include <map>
//#include "G4TOutputText.hh"
#include "G4TPhantomParameterisation.hh"
#include "G4TPhantomNestedParameterisation.hh"
//#include "G4TVolumeBuilderUsingGDML.hh"
//#include "G4TVolumeBuilderUsingTEXTGeomtry.hh"
#include "G4TPointDataGeneration.hh"
#include "G4TStlToGdml.hh"
#include "G4TTETModelImport.hh"
//#include "G4TTETParameterisation.hh"
#include "TETParameterisation.hh"

//#include "G4Navigator.hh"

#include "G4VIStore.hh"
#include "G4Colour.hh"

#ifdef G4MPI_USE
//#include "G4MPImanager.hh"
#include "mpi.h"
#endif

// Voxelized arrays

extern  size_t* MateIDs; // index of material of each voxel unsigned int* fMateIDs; // index of material of each voxel
extern  G4String* CopyNumberRegionNameMap;
extern  G4float* CopyNumberXPos;
extern  G4float* CopyNumberYPos;
extern  G4float* CopyNumberZPos;
extern  G4float* CopyNumberMassSize;
extern std::map<G4String, G4Colour> RegionNameColour;

extern std::vector<G4String> SourceRegionsNamesToBeIgnoredValues;

extern  G4double AllGeometryVolume;
extern  G4double AllGeometryMass;

extern double* CumulativeActivities;

// Source Data
extern G4String GeometryFileType;

extern G4String DataFilesExtension;

extern G4String ParticleName;

extern G4String SourceType;
extern G4String SourceRegionName;
extern G4ThreeVector boxCenterPos;
extern G4int NumberOfGenPointsToSave;
extern G4ThreeVector BoxDimGene;
extern G4ThreeVector newOrgPos3Vector;
extern G4ThreeVector SourcePosition;
extern G4ThreeVector SourceRotation;
extern G4String SourceSolid;
extern G4String SourceSurface;
extern G4String SourcePlane;
extern G4String SourceAxis;
extern G4double Radius;
extern G4double RadiusIn;
extern G4double BeamSDev;
extern G4double HalfX, HalfY, HalfZ;
extern G4double ThetaMin;
extern G4double ThetaMax;
extern G4double PhiMin;
extern G4double PhiMax;
extern G4ThreeVector SourceRotVector1;
extern G4ThreeVector SourceRotVector2;

extern G4String EnergyDistribution;
extern G4double GaussSDev;
extern G4double GaussMean;
extern G4double UniformEmin;
extern G4double UniformEmax;
extern G4double RayleighEmax;
extern G4double MonoEnergy;
extern G4double SpectrumMaxEnergy;
extern G4double FileEnergyCharacterizer;
extern G4double RadioNuclideMaxEnergy;

extern std::map<G4double, G4double> EnergyValueProbability;
extern std::vector<double> RadioNuclideProbVec;
extern std::vector<unsigned int> RadioNuclidePartNameVec;
extern std::vector<unsigned int> RadioNuclideSpectrumOrDiscreteVec;
extern std::vector<double*> RadioNuclideEneVec;

extern G4String MomDirDistribution;
extern G4double Theta;
extern G4double Phi;

extern G4String GeneratePosFlag ;
extern G4String GenerateEneFlag ;
extern G4String GenerateMomDirFlag ;
extern G4String UseGeneratedData ;
extern G4String ShowBox ;
extern G4String TestPointsPositions ;
extern G4String GeometrySymbol ;

extern G4String PositionDataFile, EnergyDataFile, MomDirDataFile;

extern G4String ScriptsDirectoryPath, DataDirectoryPath;

extern G4String ResultDirectoryPath;

// Run And Score
extern G4String MPISimulationNum;
extern G4int BatchsNumber;
extern G4int ThreadsNumber;
extern G4int EventNumberInOneThreads;
extern G4String AccuracyCalculationLevel;
extern G4String organs_to_score ;
extern G4String variable_To_Score ;
extern G4String PlanesToVisualize;
extern bool ForceSolid;
extern bool UseVoxelsColour;
extern G4int MinPlaneID;
extern G4int MaxPlaneID;

// Voxelized Geometry
extern G4bool VOX_USE;
extern G4double VoxXHalfSize;
extern G4double VoxYHalfSize;
extern G4double VoxZHalfSize;
extern G4int VoxXNumber;
extern G4int VoxYNumber;
extern G4int VoxZNumber;
extern G4bool UseDicomCumAct;
extern G4String GenerateVoxelsResuls;
extern G4int VoxPETXNumber;
extern G4int VoxPETYNumber;
extern G4int VoxPETZNumber;
extern G4double VoxPETXHalfSize;
extern G4double VoxPETYHalfSize;
extern G4double VoxPETZHalfSize;

//extern G4Navigator* aNavigator;

extern G4VPhysicalVolume*  WorldPhysicalVolume;

extern std::vector<G4String> NewRankSourceParticlesNamesValues;
extern std::vector<G4String> NewRankSourceRegionsNamesValues;
extern std::vector<G4ThreeVector> NewRankSourceRegionsBoxDimValues;
extern std::vector<G4ThreeVector> NewRankSourceRegionsPosValues;
extern std::vector<std::vector<unsigned int>> NewRankVoxelsIDsOfSourceRegion;
extern std::vector<std::vector<G4Tet*>> NewRankTETOfSourceRegion;
extern std::vector<G4ThreeVector> NewRankTETBoxMinOfSourceRegion;
extern std::vector<G4ThreeVector> NewRankTETBoxDimOfSourceRegion;
extern std::vector<G4double> NewRankSourceEnergiesValues;
extern std::vector<G4String> NewRankSourceMomDirsValues;
extern std::vector<G4double> NewRankSourceMomDirsDirectedThetaValues;
extern std::vector<G4double> NewRankSourceMomDirsDirectedPhiValues;
extern std::vector<G4String> NewRankSourcePositionDataFiles;
extern std::vector<G4String> NewRankSourceEnergyDataFiles;
extern std::vector<G4String> NewRankSourceMomDirDataFiles;

extern G4int ParamType;

// TET Geometry
extern G4String TETNodeDataFile, TETEleDataFile, TETMatDataFile;
extern G4double MinTETPhantom;
extern G4double MaxTETPhantom;
extern bool MaterialNameAsRegionName;


// Physics
extern G4String CutInRangeData;
extern G4String EnergyThresholdsData;
extern G4bool IsEcutsSet;
extern G4bool IsDcutsSet;
extern G4double CutsEnergy;
extern G4double CutsDistance;
extern G4String ParticlePysics;
extern G4String PhotoElectricEffectModel;
extern G4String PolarizedPhotoElectricEffectModel;
extern G4String ComptonScatteringModel;
extern G4String PolarizedComptonModel;
extern G4String GammaConversionModel;
extern G4String PolarizedGammaConversionModel;
extern G4String RayleighScatteringModel;
extern G4String GammaConversionToMuonModel;
extern G4String ElectronIonisationModel;
extern G4String ElectronBremModel;
extern G4String HadronIonisationModel;
extern G4String HadronBremModel;
extern G4String IonIonisationModel;
extern G4bool GenerateCrossSectionTableFlag;
extern G4String ParticleForCrossSection;
extern std::vector<G4double> EnergiesForCrossSectionValues;

class G4VPhysicalVolume;
class G4TVolumeConstruction : public G4VUserDetectorConstruction
{
public:

    G4TVolumeConstruction();
    ~G4TVolumeConstruction();
    G4VPhysicalVolume* Construct();

    G4VPhysicalVolume* ConstructVolumes();
    G4VPhysicalVolume* ConstructDICOMVolume();
    G4VPhysicalVolume* ConstructVOXELVolume();
    G4VPhysicalVolume* ConstructTETPhantomVolume();
    G4VPhysicalVolume* ConstructPhantomFromVoxIds();

    void ShowBoxVolume();
    void makeVolumeVisualization();
    void CreateMaterialsGDMLTags();

    virtual void ConstructSDandField();

    void CreateVolumesData();

    void setUserLimits();
    G4VIStore* CreateImportanceStore();

    void GenerateEventsDataInMPIMode();
    void setRankDataForMPIMode();

private:

    std::ostringstream RanksSourceData ;

    G4String BoxLVName;

    void GenerateICRPMaterialsCommands();
    void GenerateDataForVoxelsIdsFilePhantom();
    void ReadVoxelsIDsAndFillCNMatIDsMassColour();

    //World and Geometry types commands

    G4double G4Density_to_gPerCm3;
    G4double G4Density_to_kgPerMm3;
    G4double G4Density_to_kgPerMm3ToGPercm3;
    G4double G4Mass_to_Kg;

    G4String SaveGeoToTextFile;

    G4ThreeVector WorldHalfSize ;
    G4String WorldMaterialName ;

    G4String GeometryPath;
    G4String GeometryFileOrDir;
    G4String RegionsMassDataPath;
    G4bool UseRegionsMassDataPath;

    std::map<G4String , G4bool> organBoolToGen;
    G4VPhysicalVolume* WorldPhysicalVolume;
    std::vector<G4String> OrganNamesVector;
    std::map<G4String, G4double> OrganNameDensityMap;
    std::map<G4String, G4double> OrganNameMassMap;
    std::map<G4String, G4double> OrganNameVolumeMap;
    std::map<G4String, G4ThreeVector> OrganNamePositionMap;
    std::map<G4String, G4RotationMatrix> OrganNameRotMatrixMap;

    G4int ElemInc = 0;
    std::ostringstream ElementsGDMLText;
    std::ostringstream MaterialsGDMLtext;


public:

    G4String getGeometrySymbol() const {return GeometrySymbol;}

    G4String getRanksSourceData() const {return RanksSourceData.str();}

    G4String getAllMaterialsGDMLText() const {
        return (ElementsGDMLText.str()+MaterialsGDMLtext.str()).c_str();
    }

    void setMPISimulationNum(G4String nn ) { MPISimulationNum = nn;}
    G4String getMPISimulationNum() const { return MPISimulationNum;}

    G4double getG4Density_to_gPerCm3() const { return G4Density_to_gPerCm3;}
    G4double getG4Mass_to_Kg() const { return G4Mass_to_Kg;}
    G4double getG4Density_to_kgPerMm3() const { return G4Density_to_kgPerMm3;}

    void setGeometryTETDataFiles(G4String a, G4String b, G4String c){TETNodeDataFile = a; TETEleDataFile = b; TETMatDataFile = c;}

    void setGeometryFileType(G4String str){    GeometryFileType = str;}
    void setGeometryFileOrDir(G4String str){    GeometryFileOrDir = str;}
    void setRegionsMassDataPath(G4String File, G4bool j){    RegionsMassDataPath = File; UseRegionsMassDataPath = j; }
    void setGeometryPath(G4String GDMLDir){    GeometryPath = GDMLDir;}

    void constructWorldVolume(G4String , G4ThreeVector);
    G4String getGeometryPath() const { return GeometryPath;}
    G4String getGeometryFileType() const {return GeometryFileType;}

    void TestAndShowUserInputs();

    G4VPhysicalVolume* getWorldPhyVolume()const {return WorldPhysicalVolume;}
    std::vector<G4String> GetOrganNamesVector() const { return OrganNamesVector;}
    std::map<G4String, G4double> GetOrganNameDensityMap() const {return OrganNameDensityMap;}
    std::map<G4String, G4double> GetOrganNameMassMap() const {return OrganNameMassMap;}
    std::map<G4String, G4double> GetOrganNameVolumeMap() const {return OrganNameVolumeMap;}
    std::map<G4String, G4ThreeVector> GetOrganNamePositionMap() const {return OrganNamePositionMap;}

    std::vector<G4String> SimulatedCouplesData;

    std::vector<G4String> SourceParticlesNamesValues;
    std::vector<G4String> SourceRegionsNamesValues;
    std::vector<G4ThreeVector> SourceRegionsBoxDimValues;
    std::vector<G4double> SourceEnergiesValues;
    std::vector<G4String> SourceMomDirsValues;
    std::vector<G4double> SourceMomDirsDirectedThetaValues;
    std::vector<G4double> SourceMomDirsDirectedPhiValues;


    void setInitializationSourceForCommandSourceInputs(){
        SourceParticlesNamesValues.clear();SourceParticlesNamesValues.empty();
        SourceRegionsNamesValues.clear();SourceRegionsNamesValues.empty();
        SourceRegionsBoxDimValues.clear();SourceRegionsBoxDimValues.empty();
        SourceEnergiesValues.clear();SourceEnergiesValues.empty();
        SourceMomDirsValues.clear();SourceMomDirsValues.empty();
    }

    void setInitializationSourceParticlesInputs(){
        SourceParticlesNamesValues.clear();SourceParticlesNamesValues.empty();
    }
    void setInitializationSourceEnergiesInputs(){
        SourceEnergiesValues.clear();SourceEnergiesValues.empty();
    }
    void setInitializationSourcePositionsInputs(){
        SourceRegionsNamesValues.clear();SourceRegionsNamesValues.empty();
        SourceRegionsBoxDimValues.clear();SourceRegionsBoxDimValues.empty();
    }
    void setInitializationSourceMomDirsInputs(){
        SourceMomDirsValues.clear();SourceMomDirsValues.empty();
    }

    std::vector<G4String> getSourceParticlesNamesValues() const { return SourceParticlesNamesValues;}
    std::vector<G4String> getSourceRegionsNamesValues() const { return SourceRegionsNamesValues;}
    std::vector<G4ThreeVector> getSourceRegionsBoxDimValues() const { return SourceRegionsBoxDimValues;}
    std::vector<G4double> getSourceEnergiesValues() const { return SourceEnergiesValues;}
    std::vector<G4String> getSourceMomDirsValues(){ return SourceMomDirsValues; }

    void setSourceParticlesNamesValues(G4String n ){ SourceParticlesNamesValues.push_back(n); }
    void setSourceRegionsNamesValues(G4String n ){ SourceRegionsNamesValues.push_back(n); }
    void setSourceRegionsBoxDimValues(G4ThreeVector n ){ SourceRegionsBoxDimValues.push_back(n); }
    void setSourceEnergiesValues(G4double n ){ SourceEnergiesValues.push_back(n); }
    void setSourceMomDirsValues(G4String n ){ SourceMomDirsValues.push_back(n); }
    void setSourceMomDirsDirectedThetaValues(G4double n ){ SourceMomDirsDirectedThetaValues.push_back(n); }
    void setSourceMomDirsDirectedPhiValues(G4double n ){ SourceMomDirsDirectedPhiValues.push_back(n); }

    void setSourceRegionsNamesToBeIgnoredValues(G4String n ){ SourceRegionsNamesToBeIgnoredValues.push_back(n); }

    std::vector<G4String> VolumesNotVisualized;
    void setVolumesNotVisualized(G4String n ){ VolumesNotVisualized.push_back(n); }
    void setGeometrySymbol(G4String n ){ GeometrySymbol = n; }

    struct SourceDataStruct
    {
        G4int ID;
        G4String SrcPartName;
        G4String srcRegName;
        G4ThreeVector SrcRegPos;
        G4ThreeVector SrcBoxDim;
        std::vector <unsigned int> SrcRegVoxelsIDs;
        G4String SrnEneDistName;
        G4double SrnEneVal;
        G4String SrnMomDirDistName;
        G4double SrnMomDirDistTheta;
        G4double SrnMomDirDistPhi;
        G4String SrcPosDataFilePath;
        G4String SrcEneDataFilePath;
        G4String SrcMomDirDataFilePath;
    };

    std::vector<SourceDataStruct> RankThreadSourceData;
    std::vector<SourceDataStruct> GetRankThreadSourceData() const { return RankThreadSourceData;}

private:

    // for biasing

    G4String RegionWhereToBias;
    G4String BiasFlag;

    // for UserLimits

    G4String VolumeOfLimits;
    G4double max_allowed_step_size;
    G4double max_total_track_length;
    G4double max_total_time_of_flight;
    G4double min_kinetic_energy;
    G4double min_remaining_range;

    // for Source geometry

public:

    G4String getBiasFlag() const { return BiasFlag;}

    G4String getPositionDataFile() const { return PositionDataFile;}
    G4String getEnergyDataFile() const { return EnergyDataFile;}
    G4String getMomDirDataFile() const { return MomDirDataFile;}

    G4String getGeneratePostions() const { return GeneratePosFlag;}
    G4String getGenerateEnergies() const { return GenerateEneFlag;}
    G4String getGenerateMomDirs() const { return GenerateMomDirFlag;}
    G4String getUseGeneratedData() const { return UseGeneratedData;}


    G4String getSourceType() const { return SourceType;}
    G4ThreeVector getSourcePosition() const { return SourcePosition;}
    G4ThreeVector getSourceRotation() const { return SourceRotVector1;}
    G4ThreeVector getSourceRotVector1() const { return SourceRotVector1;}
    G4ThreeVector getSourceRotVector2() const { return SourceRotVector2;}

    G4String getSourceSolid() const { return SourceSolid ;}
    G4String getSourcePlane() const { return SourcePlane;}
    G4String getSourceSurface() const { return SourceSurface;}
    G4String getSourceAxis() const { return SourceAxis;}

    G4double getRadiusIn() const { return RadiusIn;}
    G4double getBeamSDev() const { return BeamSDev;}

    G4double getRadius() const { return Radius;}
    G4double getHalfX() const { return HalfX;}
    G4double getHalfY() const { return HalfY;}
    G4double getHalfZ() const { return HalfZ;}

    G4double getThetaMin() const { return ThetaMin;}
    G4double getThetaMax() const { return ThetaMax;}
    G4double getPhiMin() const { return PhiMin;}
    G4double getPhiMax() const { return PhiMax;}

    G4String getShowBox() const { return ShowBox;}
    G4String getTestPointsPositions() const { return TestPointsPositions;}
    G4int getPointNumberToSave() const { return NumberOfGenPointsToSave;}
    G4ThreeVector getOrgPosThreeVector() const { return newOrgPos3Vector;}
    G4ThreeVector getBoxDimGene() const { return BoxDimGene;}
    G4String getOrganSource() const { return SourceRegionName;}
    G4String getParticleName() const { return ParticleName;}
    G4double getMonoEnergy() const { return MonoEnergy ;}
    G4String getMomDirDistribution() const { return MomDirDistribution;}
    G4String getEnergyDistribution() const { return EnergyDistribution;}
    G4double getTheta() const { return Theta ;}
    G4double getPhi() const { return Phi ;}
    G4double getGaussSDev() const { return GaussSDev ;}
    G4double getUniformEmin() const { return UniformEmin ;}
    G4double getUniformEmax() const { return UniformEmax ;}
    G4double getRayleighEmax() const { return RayleighEmax ;}
    G4double getGaussMean() const { return GaussMean ;}
    std::map<G4double, G4double> getEnergyValueProbability() const { return EnergyValueProbability ;}
    G4double getSpectrumMaxEnergy() const { return SpectrumMaxEnergy ;}
    G4double getFileEnergyCharacterizer() const { return FileEnergyCharacterizer ;}
    G4double getRadioNuclideMaxEnergy() const { return RadioNuclideMaxEnergy ;}


    void setSourceType(G4String nn){    SourceType = nn ;}
    void setSourcePosition(G4ThreeVector newPos){   SourcePosition = newPos;}
    void setSourceRotVector1(G4ThreeVector newPos){   SourceRotVector1 = newPos;}
    void setSourceRotVector2(G4ThreeVector newPos){   SourceRotVector2 = newPos;}

    void setSourceSolid(G4String newPos){   SourceSolid = newPos;}
    void setSourcePlane(G4String nn ){    SourcePlane = nn;}

    void setSourceSurface(G4String nn ){    SourceSurface = nn;}
    void setSourceAxis(G4String nn ){    SourceAxis = nn;}

    void setSourceRotation(G4ThreeVector nn ){    SourceRotation = nn;}
    void setRadius(G4double nn ){ Radius = nn;}
    void setHalfX(G4double nn ){ HalfX = nn;}
    void setHalfY(G4double nn ){ HalfY = nn;}
    void setHalfZ(G4double nn ){ HalfZ = nn;}

    void setThetaMin(G4double nn ){ ThetaMin = nn;}
    void setThetaMax(G4double nn ){ ThetaMax = nn;}
    void setPhiMin(G4double nn ){ PhiMin = nn;}
    void setPhiMax(G4double nn ){ PhiMax = nn;}

    void setRadiusIn(G4double nn ){ RadiusIn = nn;}
    void setBeamSDev(G4double nn ){ BeamSDev = nn;}

    void setSourceRegionName(G4String nn){   SourceRegionName = nn;}
    // PositionDataFile = SourceRegionName + "PositionData.bin";

    G4String  setSourcePositionFileName();
    G4String  setSourceEnergyFileName();
    G4String  setSourceMomDirFileName();

    void setPointNumberToSave(G4int nn){    NumberOfGenPointsToSave = nn ;}

    void setOrgPosThreeVector(G4ThreeVector newPos ){    newOrgPos3Vector = newPos; }
    void setXYZOfBox(G4ThreeVector nn , G4bool newBoxWidthUse){    BoxDimGene = nn; }
    void setParticleName(G4String nn){   ParticleName = nn;}
    void setMonoEnergy(G4double nn){   MonoEnergy = nn;}

    void setEnergyValueProbability(G4double nn, G4double ll){   EnergyValueProbability[nn] = ll;}
    void setProbabilityParticleNameEnergy(double aa, unsigned int bb, unsigned int cc, double* dd){
        RadioNuclideProbVec.push_back(aa);
        RadioNuclidePartNameVec.push_back(bb);
        RadioNuclideSpectrumOrDiscreteVec.push_back(cc);
        RadioNuclideEneVec.push_back(dd);
    }

    void setResultDirectoryPath(G4String nn){
        ResultDirectoryPath = nn ;
    }

    void setInitialDirectionModel(G4String nn){   MomDirDistribution = nn ;}
    void setPhi(G4double nn){   Phi = nn;}
    void setTheta(G4double nn){   Theta = nn;}
    void setEnergyDistribution(G4String nn){   EnergyDistribution = nn ;}
    void setGaussMean(G4double nn){   GaussMean = nn ;}
    void setGaussSDev(G4double nn){    GaussSDev = nn ;}
    void setSpectrumMaxEnergy(G4double nn){    SpectrumMaxEnergy = nn ;}
    void setFileEnergyCharacterizer(G4double nn){    FileEnergyCharacterizer = nn ;}
    void setRadioNuclideMaxEnergy(G4double nn){    RadioNuclideMaxEnergy = nn ;}

    void setUniformEmin(G4double nn){    UniformEmin = nn ;}
    void setUniformEmax(G4double nn){    UniformEmax = nn ;}
    void setRayleighEmax(G4double nn){   RayleighEmax = nn ;}
    void setShowBox(G4String nn){    ShowBox = nn ;}
    void setTestPointsPositions(G4String nn){    TestPointsPositions = nn ;}

    void setGeneratePositions(G4String nn ){    GeneratePosFlag = nn;}
    void setGenerateEnergies(G4String nn ){    GenerateEneFlag = nn;}
    void setGenerateMomDirs(G4String nn ){    GenerateMomDirFlag = nn;}
    void setUseGeneratedData(G4String nn ){    UseGeneratedData = nn;}

private:


public:

    void setOrgans_to_score(G4String lll){    organs_to_score = lll;}
    void setVariable_To_Score(G4String lll){    variable_To_Score = lll;}
    void setAccuracyCalculationLevel(G4String sss){    AccuracyCalculationLevel = sss;}
    void setNumberOfThreads(G4int);
    void setNumberOfEventsPerThread(G4int newNumber){    EventNumberInOneThreads = newNumber;}

    G4int    getNumberOfThreads() const { return ThreadsNumber;}
    G4int    getNumberOfEventsPerThread() const { return EventNumberInOneThreads;}
    G4String setAccuracyCalculationLevel() const { return AccuracyCalculationLevel;}

    G4String getorgans_to_score() const { return organs_to_score;}
    G4String getvariable_To_Score()const {return variable_To_Score;}

    // for analysis

private:

    G4String Graphs_Data ;
    G4String compare_type ;
    G4String ref_File_Path ;
    G4String ref_Name ;

    G4String RegionVariableName ;
    //G4double EnergyGraphValue;
    G4String GenerateRelativeErrGraph;
    G4String DifferenceMethod;
    G4String GenerateRelativeSDevGraph;
    G4String GenerateRegionsVariableGraph;

    G4String TissueFactors;
    G4String RadiationFactors;
    G4String QuantitiesUnits;

    std::vector<G4String> RadioTracerDataVector;
    std::vector<G4String> RadioTracerBiokineticVector;

    //G4double InjectedActivity;
    //G4String SAFDataReference;

    G4String UseLogE;
    G4String UseLogVariable;
    G4String UseGridXY;
    G4String PrintTitle;
    G4String LegendPos;
    G4double LegendXWidth ;
    G4double LegendYHeight ;
    G4String AddErrorBarInGraphs ;
    G4String graphs_Ext ;

    G4int SliceID;
    G4String SliceFor2DGraph;
    G4String BeamAxis;
    G4String DoseProfilQuantity;
    G4String EventsDataHistograms;
    G4String EventsPositionHistogram;
    G4String EventsEnergyHistogram;
    G4String EventsMomDirHistogram;

public:

    // for analysis data arguments

    void setGraphs_Data(G4String lll){   Graphs_Data = lll;}
    void setCompare_type(G4String lll){    compare_type = lll;}
    void setRef_File_Path(G4String lll){    ref_File_Path = lll;}
    void setRef_Name(G4String lll){    ref_Name = lll;}

    void setRegionVariableName(G4String l){    RegionVariableName = l;}
    //void setEnergyGraphValue(G4double l){    EnergyGraphValue = l;}
    void setGenerateRelativeErrGraph(G4String lll){    GenerateRelativeErrGraph = lll;}
    void setDifferenceMethod(G4String lll){    DifferenceMethod = lll;}
    void setGenerateRelativeSDevGraph(G4String lll){    GenerateRelativeSDevGraph = lll;}
    void setGenerateRegionsVariableGraph(G4String lll){    GenerateRegionsVariableGraph = lll;}

    void setRadiationFactors(G4String lll){    RadiationFactors = lll;}
    void setQuantitiesUnits(G4String lll){    QuantitiesUnits = lll;}
    void setTissueFactors(G4String lll){    TissueFactors = lll;}

    void setRadioTracerData(G4String lll){    RadioTracerDataVector.push_back(lll);}
    void setRadioTracerBiokinetic(G4String lll){    RadioTracerBiokineticVector.push_back(lll);}

    //void setInjectedActivity(G4double lll){    InjectedActivity = lll;}
    //void setSAFDataReference(G4String lll){    SAFDataReference = lll;}

    void setSliceID(G4int n){SliceID = n;}
    void setSliceFor2DGraph(G4String n){SliceFor2DGraph = n;}
    void setBeamAxis(G4String n){BeamAxis = n;}
    void setDoseProfilQuantity(G4String n){DoseProfilQuantity = n;}

    void setEventsDataHistograms(G4String n){EventsDataHistograms = n;}
    void setEventsPositionHistogram(G4String n){EventsPositionHistogram = n;}
    void setEventsEnergyHistogram(G4String n){EventsEnergyHistogram = n;}
    void setEventsMomDirHistogram(G4String n){EventsMomDirHistogram = n;}

    G4String getGraphs_Data() const { return Graphs_Data;}
    G4String getcompare_type() const { return compare_type;}
    G4String getref_File_Path() const { return ref_File_Path ;}
    G4String getref_Name() const { return ref_Name;}

    G4String getUseLogE() const { return UseLogE;}
    G4String getUseLogVariable() const { return UseLogVariable;}
    G4String getUseGridXY() const { return UseGridXY;}
    G4String getPrintTitle() const { return PrintTitle;}
    G4String getLegendPos() const { return LegendPos;}
    G4double getLegendXWidth() const { return LegendXWidth;}
    G4double getLegendYHeight() const { return LegendYHeight;}
    G4String getAddErrorBarInGraphs() const { return AddErrorBarInGraphs;}
    G4String getgraphs_Ext() const { return graphs_Ext;}

    void setUseLogE(G4String lll){    UseLogE = lll;}
    void setUseLogVariable(G4String lll){    UseLogVariable = lll;}
    void setUseGridXY(G4String lll){    UseGridXY = lll;}
    void setPrintTitle(G4String lll){    PrintTitle = lll;}
    void setLegendXWidth(G4double lll){    LegendXWidth = lll;}
    void setLegendYHeight(G4double lll){    LegendYHeight = lll;}
    void setLegendPos(G4String lll){    LegendPos = lll;}
    void setAddErrorBarInGraphs(G4String lll){    AddErrorBarInGraphs = lll;}
    void setGraphs_Ext(G4String lll){    graphs_Ext = lll;}

    G4String getRegionVariableName() const { return RegionVariableName;}
    //G4double getEnergyGraphValue() const { return EnergyGraphValue;}
    G4String getGenerateRegionsVariableGraph() const { return GenerateRegionsVariableGraph;}
    G4String getGenerateRelativeSDevGraph() const { return GenerateRelativeSDevGraph;}
    G4String getGenerateRelativeErrGraph() const { return GenerateRelativeErrGraph;}
    G4String getDifferenceMethod() const { return DifferenceMethod;}

    G4String getRadiationFactors() const { return RadiationFactors;}
    G4String getQuantitiesUnits() const { return QuantitiesUnits;}
    G4String getTissueFactors() const { return TissueFactors;}

    std::vector<G4String> getRadioTracerData() const { return RadioTracerDataVector;}
    std::vector<G4String> getRadioTracerBiokinetic() const { return RadioTracerBiokineticVector;}
    //G4double getInjectedActivity() const { return InjectedActivity;}
    //G4String getSAFDataReference() const { return SAFDataReference;}

    G4int    getSliceID()const {return SliceID;}
    G4String getSliceFor2DGraph()const {return SliceFor2DGraph;}
    G4String getBeamAxis()const {return BeamAxis;}
    G4String getDoseProfilQuantity()const {return DoseProfilQuantity;}

    G4String getEventsDataHistograms   ()const {return EventsDataHistograms;}
    G4String getEventsPositionHistogram()const {return EventsPositionHistogram;}
    G4String getEventsEnergyHistogram  ()const {return EventsEnergyHistogram;}
    G4String getEventsMomDirHistogram  ()const {return EventsMomDirHistogram;}


    // for Material constructing
private:

    G4String MaterialName, FracOrNum; G4int MaterialCompNumber;
    std::map<G4int,G4String> MaterialIDName;

public:
    std::map<G4int,G4String> getMaterialIDName  ()const {return MaterialIDName;}

    void createElement(G4String, G4double, G4double);
    void createNistMaterial(G4String, G4int);
    void createMaterialFromComponents(G4String, G4int, G4int, G4double, G4String);
    void AddMatElement(G4String, G4double, G4int);


    // for Construct geometry
private:

    G4String STLVolumeMaterial;
    G4bool CPPLogVolAreBuilt;

    G4String DumpGeomTextFile;
    G4bool dumpGeom;

public:

    void dumpGeometryToTextFile(G4String n ){ dumpGeom = true; DumpGeomTextFile= n;}
    void setSTLVolumeMaterial(G4String n);
    void placeVolume(G4String, G4String, G4ThreeVector, G4ThreeVector);
    void PlaceVolumeFromRegionfile(G4String, G4String);

private:

    G4Material* CurrentMaterial;
    G4VSolid* CurrentSolid;
    G4LogicalVolume* CurrentLogicalVol;
    G4VPhysicalVolume* CurrentPhysicalVol;
    G4VPhysicalVolume* CurrentMotherPhysicalVol;

    std::map<G4String,G4Material*> CreatedMaterials;
    std::map<G4String,G4Element*> CreatedElements;

    // For Dicom and Voxelized geometry

private:

    G4int VoxRegionMinX;
    G4int VoxRegionMaxX;
    G4int VoxRegionMinY;
    G4int VoxRegionMaxY;
    G4int VoxRegionMinZ;
    G4int VoxRegionMaxZ;
    G4String  VoxRegionName ;
    G4ThreeVector VoxContainerPos;
    G4ThreeVector VoxContainerRot;

    G4double VoxelVolume;
    G4double VoxelMass;
    G4double RegionVolume ;
    G4double RegionMass ;

public:

    G4int getVoxXNumber() const { return VoxXNumber ;}
    G4int getVoxYNumber() const { return VoxYNumber ;}
    G4int getVoxZNumber() const { return VoxZNumber ;}

    G4ThreeVector getVoxContainerPos() const { return VoxContainerPos;}

    G4String getGenerateVoxelsResuls() const { return GenerateVoxelsResuls ;}
    void setGenerateVoxelsResuls(G4String n ){GenerateVoxelsResuls = n;}
    void setPlanesToVisualize(G4String n ){PlanesToVisualize = n;}
    void setForcedSolid(G4bool n ){ForceSolid = n;}

    void setMinPlaneID(G4int n){    MinPlaneID = n ;}
    void setMaxPlaneID(G4int n){    MaxPlaneID = n ;}

    void setVoxContainerPos(G4ThreeVector n){    VoxContainerPos = n ;}
    void setVoxContainerRot(G4ThreeVector n){    VoxContainerRot = n ;}

    void setVoxXNumber(G4int n){    VoxXNumber = n ;}
    void setVoxYNumber(G4int n){    VoxYNumber = n ;}
    void setVoxZNumber(G4int n){    VoxZNumber = n ;}
    void setVoxRegionName(G4String);
    void setVoxRegionMinX(G4int );
    void setVoxRegionMaxX(G4int );
    void setVoxRegionMinY(G4int );
    void setVoxRegionMaxY(G4int );
    void setVoxRegionMinZ(G4int );
    void setVoxRegionMaxZ(G4int );
    void setAllGeomAsMinMaxVoxRegionLimits();

    void setVoxRegionsFractions(G4String, double);

    void setTETRegionName(G4String);

    G4ThreeVector getVoxelsSize() const {
        return G4ThreeVector(VoxXHalfSize,VoxYHalfSize,VoxZHalfSize) ;
    }


    std::map<G4String,std::map<G4String, G4double>>RegionRegionFraction ;

    std::vector<G4Colour> MaterialColour;
    std::map<G4String,unsigned int> RegionNumberOfCN;
    std::vector<G4Material*> VoxelsMaterials;

    std::map<G4String,unsigned int> getRegionNumberOfCNMap() const { return RegionNumberOfCN ;}
    /*
    size_t* MateIDs; // index of material of each voxel unsigned int* fMateIDs; // index of material of each voxel
    G4String* CopyNumberRegionNameMap; G4String* GetCopyNumberRegionNameMap() const { return CopyNumberRegionNameMap;}
    G4float* CopyNumberXPos;  G4float* getCopyNumberXPos() const { return CopyNumberXPos;}
    G4float* CopyNumberYPos;  G4float* getCopyNumberYPos() const { return CopyNumberYPos;}
    G4float* CopyNumberZPos;  G4float* getCopyNumberZPos() const { return CopyNumberZPos;}
    G4float* CopyNumberMassSize; G4float* getCopyNumberMassSize() const { return CopyNumberMassSize;}
    */

    void ConstructVoxDcmContainerGeometry();
    void InitializeCNMassColourRegNamePosAxisIDs();
    G4VPhysicalVolume* ConstructVoxeDcmGeometry();

    G4VPhysicalVolume* ConstructTETGeometry();

    void VisualizationVoxelizedGeometry();

    G4Material* defaultMat;
    G4Box* voxel_solid;
    G4LogicalVolume* voxel_logic;

    G4Box* ContSolidVoll;
    G4LogicalVolume* ContLogicalVoll;
    G4VPhysicalVolume* ContPhysicalVoll;

    // for Voxelized geometry

    G4String  VoxRegionMaterialName ;
    G4String  VoxDefaultMaterialName ;

    G4double getVoxXHalfSize() const { return VoxXHalfSize ;}
    G4double getVoxYHalfSize() const { return VoxYHalfSize ;}
    G4double getVoxZHalfSize() const { return VoxZHalfSize ;}

    void setVoxXHalfSize(G4double n ){VoxXHalfSize = n;}
    void setVoxYHalfSize(G4double n ){VoxYHalfSize = n;}
    void setVoxZHalfSize(G4double n ){VoxZHalfSize = n;}
    void CreateRegionVoxelsData();
    void CreateDefaultVOXELRegionData();

    //std::ostringstream RegionsDataText;

    void ShowMaterialRegionVoxelsData();
    void CreateRegionsDataFile();
    G4int ConvertVoxelPosToXYZVoxelIDs(G4String, G4String, G4double);
    G4int ConvertVoxelPosToXYZToID(G4double, G4double, G4double);

    void setVoxDefaultMaterialName(G4String);
    void CalcutateSizeOfVoxelizedArrays();

    // for DICOM geometry

    void getDicomDataAndConstruct();

    G4int RegionFromInputInc;

    G4String DcmMaterialName;
    G4double DcmCtNumber;
    G4double DcmCtDensity;
    G4double DcmRegionMinDensity;
    G4bool UseDcmRegionMinDensity;
    G4bool UseDcmRegionMaxDensity;
    G4bool UseVoxelMat;
    G4int VoxelMatForSeg;

    std::vector<G4int> VoxelMatForSegVec ;

    G4double DcmRegionMaxDensity;
    int DicomPixelsCompression;

    void setDcmMaterialName(G4String);
    void setDcmCtNumber(G4double);
    void setDcmCtDensity(G4double);

    void setDcmPixelsCompression(G4int num){ DicomPixelsCompression = num;}

    void setDcmPETSerieResidenceTime(G4double nn){ SeriesResidenceTime.push_back(nn);}
    void setDcmPixelsCompressionX(G4int num){ SeriesCompressionX.push_back(num);}
    void setDcmPixelsCompressionY(G4int num){ SeriesCompressionY.push_back(num);}
    void setDcmPixelsCompressionXY(G4int num){ SeriesCompressionXY.push_back(num);}

    std::vector<G4double> SeriesResidenceTime ;
    std::vector<G4int> SeriesCompressionX ;
    std::vector<G4int> SeriesCompressionY ;
    std::vector<G4int> SeriesCompressionXY ;

    G4String getDcmMaterialName() const { return DcmMaterialName ;}
    G4double getDcmCtNumber() const { return DcmCtNumber ;}
    G4double getDcmCtDensity() const { return DcmCtDensity ;}
    G4double getDcmRegionMinDensity() const { return DcmRegionMinDensity ;}
    G4double getDcmRegionMaxDensity() const { return DcmRegionMaxDensity ;}

    G4int getDicomPixelsCompression() const { return DicomPixelsCompression;}

    std::vector<G4String> DicomFileTypesVec ;
    G4String DicomFileType;
    void setDicomFileType(G4String n ){DicomFileType = n ; DicomFileTypesVec.push_back(n) ; }
    G4String getDicomFileType() const { return DicomFileType ; }

    G4String DicomOutTextName, DicomCTName, DicomPETName;
    G4String getDicomOutTextName() const { return DicomOutTextName;}
    G4String getDicomCTName() const { return DicomCTName;}
    G4String getDicomPETName() const { return DicomPETName;}

    std::vector<G4String> DicomFilePathsVec ;
    G4String DicomFilesDirPath;
    void setDicomDataDirPath(G4String Dir){ DicomFilesDirPath = Dir ; DicomFilePathsVec.push_back(Dir) ; }
    G4String getDicomDataDirPath() const { return DicomFilesDirPath;}

    double* Activities;
    std::vector<double*> ActivitiesVec ;
    std::vector<G4double> ActivityTimeVec ;
    std::map<G4String,G4double> PETGeneralData ;
    std::map<G4String,std::vector<double>> PETGeneralDateData ;
    std::map<G4int ,std::vector<double>> Serie_Order_DataMap ;

    //double* getCumulativeActivities() const { return CumulativeActivities;}
    void setUsePETCumulativeAct(G4String n){ if(n == "yes"){UseDicomCumAct = true;} else {UseDicomCumAct = false;}}
    G4bool getUseDicomCumAct() const { return UseDicomCumAct;}

    void setParamType(G4int n){ParamType = n;
                               //G4cout << " ParamType " << ParamType << " fffffffffffffff \n"<<G4endl;
                              }


    G4double Time;
    G4double BiologicPeriode;
    G4double PhysicalPeriode;
    G4double EffectifPeriode;
    G4double Lambda;

    std::vector<G4Material*> DcmMaterials;
    std::vector<G4Material*> getDcmMaterialsVector() const {
        return DcmMaterials;
    }
    std::vector<G4double> ValueCT; std::vector<G4double> getValueCTVector() const {return ValueCT;}
    std::vector<G4double> ValueDensity; std::vector<G4double> getValueDensityVector() const {return ValueDensity;}
    std::map<G4float,G4String> MaterialIndices; std::map<G4float,G4String> getMaterialIndicesMap() const { return MaterialIndices;}

    std::vector<G4String> DcmRegionsNames;
    std::vector<G4int> DcmRegionsMinX;
    std::vector<G4int> DcmRegionsMaxX;
    std::vector<G4int> DcmRegionsMinY;
    std::vector<G4int> DcmRegionsMaxY;
    std::vector<G4int> DcmRegionsMinZ;
    std::vector<G4int> DcmRegionsMaxZ;

    std::vector<G4String> PhantomLimits;

    void setDcmRegionMinDensity(G4double);
    void setDcmRegionMaxDensity(G4double);
    void setVoxelsegmentedMaterial(G4int);

    std::vector<G4double> DcmRegionsMinDensityMap;
    std::vector<G4double> DcmRegionsMaxDensityMap;
    std::vector<G4bool> UseDcmRegionsMinDensityMap;
    std::vector<G4bool> UseDcmRegionsMaxDensityMap;

    std::vector<std::map<G4int,G4int>> VoxelMatForSegMapMap;
    std::map<G4int,G4int> MaterialsIDsForSegmentationMap;

    std::vector<G4int> VoxelMatForSegMap;
    std::vector<G4bool> UseVoxelMatForSegMap;

    // ICRP geometry

    std::map<int,double> MediaBloodFractionMap;
    std::map<int,double> MediaRBMFractionMap;
    std::map<int,double> MediaYBMFractionMap;
    std::map<int,double> MediaBoneFractionMap;

    std::map<int,std::vector<G4double>> ICRPMaterialsCompFra;
    std::map<int,G4String> OrganIDNameMap;
    std::map<int,G4Colour> OrganIDColourMap;

    std::map<G4String,int> MatNameIDMap;
    std::map<G4String,int> OrganNameIDMap;

    std::map<int,int> OrganIDMatIDMap;
    std::map<int,G4double> OrganIDDensityMap;
    G4int ICRPOrgansNumber;
    G4int ICRPMaterialsNumber;

    void BuildICRPMaterials();
    void CreateICRPPhantomData();
    void CreateICRPSAFsReferenceDataFile();

    // for GDML geometry

    G4ThreeVector getRegionAbsolutePosition(G4String);

    std::map<G4String,G4ThreeVector> CreatedPositionOrgans;
    std::map<G4String,G4ThreeVector> getOrgansPositionCreated() const { return CreatedPositionOrgans;}
    std::map<G4String,G4ThreeVector> CreatedRotationOrgans;
    std::map<G4String,G4ThreeVector> getOrgansRotationCreated() const { return CreatedRotationOrgans;}

    G4VPhysicalVolume*  PhysicalBoxVolume;


    bool getMaterialNameAsRegionName() const { return MaterialNameAsRegionName;}
    void setMaterialNameAsRegionName(bool n){MaterialNameAsRegionName = n;}

    // for TET geometry

    G4TTETModelImport*    tetData;
    G4ThreeVector      phantomBoxMin, phantomBoxMax;
    G4int              nOfTetrahedrons;

    G4LogicalVolume*   tetLogic;


    //G4String phantomDataPath;
    //G4String phantomName;

    G4ThreeVector phantomSize;

    void GenerateDataFromTETPhantomFiles();

public:

    G4TTETModelImport* gettetData() const {return tetData;}

    void setTETPhantomLimits(G4double a, G4double b){
        MinTETPhantom = a ;
        MaxTETPhantom = b ;

        //MinXTET = a ;
        //MaxXTET = b ;
        //MinYTET = c ;
        //MaxYTET = d ;
        //MinZTET = e ;
        //MaxZTET = f ;
    }

    /*
    G4ThreeVector boundingBox_Min;
    G4ThreeVector boundingBox_Max;

    std::vector<G4ThreeVector> vertexVector;
    std::vector<G4Tet*>        tetVector;
    std::vector<G4int*>        eleVector;
    std::vector<G4int>         materialVector;
    std::map<G4int, G4int>     numTetMap;
    std::map<G4int, G4double>  volumeMap;
    std::map<G4int, G4double>  massMap;
    std::map<G4int, G4Colour>  colourMap;
*/
    //std::map<G4int, std::vector<std::pair<G4int, G4double>>> materialIndexMap;
    //std::vector<G4int>                                       materialIndex;
    //std::map<G4int, G4Material*>                             materialMap;
    //std::map<G4int, G4double>                                densityMap;
    //std::map<G4int, G4String>                                organNameMap;


    //////////////////////////////////////////// begin Physics Commands

    void setCutInRangeData(G4String newVal){
        CutInRangeData= newVal;
        IsDcutsSet = true ;
        //G4cout<<" >> New CutsEnergy " << CutsEnergy <<G4endl;
    }

    void setEnergyThresholdsData(G4String newVal){
        EnergyThresholdsData= newVal;
        IsEcutsSet = true ;
        //G4cout<<" >> New CutsEnergy " << CutsEnergy <<G4endl;
    }

    void setCutsEnergy(G4double newVal){
        CutsEnergy= newVal*MeV;
        IsEcutsSet = true ;
        //G4cout<<" >> New CutsEnergy " << CutsEnergy <<G4endl;
    }
    void setCutsDistance(G4double newVal){
        CutsDistance = newVal*mm;
        IsDcutsSet = true ;
        //G4cout<<" >> New CutsDistance " << CutsDistance <<G4endl;
    }
    void setPhysics(G4String newVal){
        ParticlePysics = newVal;
        //G4cout<<" >> New ParticlePysics " << ParticlePysics <<G4endl;
    }
    void setPhotoElectricEffectModel(G4String newVal){
        PhotoElectricEffectModel = newVal;
        //G4cout<<" >> New PhotoElectricEffectModel " << PhotoElectricEffectModel <<G4endl;
    }
    void setPolarizedPhotoElectricEffectModel(G4String newVal){
        PolarizedPhotoElectricEffectModel = newVal;
        //G4cout<<" >> New PolarizedPhotoElectricEffectModel " << PolarizedPhotoElectricEffectModel <<G4endl;
    }
    void setComptonScatteringModel(G4String newVal){
        ComptonScatteringModel = newVal;
        //G4cout<<" >> New ComptonScatteringModel " << ComptonScatteringModel <<G4endl;
    }
    void setPolarizedComptonModel(G4String newVal){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;
        PolarizedComptonModel = newVal;
        //G4cout<<" >> New PolarizedComptonModel " << PolarizedComptonModel <<G4endl;
    }
    void setGammaConversionModel(G4String newVal){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;
        GammaConversionModel = newVal;
        //G4cout<<" >> New GammaConversionModel " << GammaConversionModel <<G4endl;
    }
    void setPolarizedGammaConversionModel(G4String newVal){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;
        PolarizedGammaConversionModel = newVal;
        //G4cout<<" >> New PolarizedGammaConversionModel " << PolarizedGammaConversionModel <<G4endl;
    }
    void setRayleighScatteringModel(G4String newVal){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;
        RayleighScatteringModel = newVal;
        //G4cout<<" >> New ParticlePysics " << RayleighScatteringModel <<G4endl;
    }
    void setGammaConversionToMuonModel(G4String newVal){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;
        GammaConversionToMuonModel = newVal;
        //G4cout<<" >> New GammaConversionToMuonModel " << GammaConversionToMuonModel <<G4endl;
    }
    void setElectronIonisationModel(G4String newVal){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;
        ElectronIonisationModel = newVal;
        //G4cout<<" >> New ElectronIonisationModel " << ElectronIonisationModel <<G4endl;
    }
    void setElectronBremModel(G4String newVal){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;
        ElectronBremModel = newVal;
        //G4cout<<" >> New ElectronBremModel " << ElectronBremModel <<G4endl;
    }

    void setHadronIonisationModel(G4String newVal){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;
        HadronIonisationModel = newVal;
        //G4cout<<" >> New HadronIonisationModel " << HadronIonisationModel <<G4endl;
    }

    void setIonIonisationModel(G4String newVal){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;
        IonIonisationModel = newVal;
        //G4cout<<" >> New IonIonisationModel " << IonIonisationModel <<G4endl;
    }

    void setGenerateCrossSectionTableFlag(G4bool n ){ GenerateCrossSectionTableFlag = n; }
    G4bool getGenerateCrossSectionTableFlag() const {return GenerateCrossSectionTableFlag;}

    void setParticleForCrossSection(G4String n ){ ParticleForCrossSection = n; }
    G4String getParticleForCrossSectionValues() const {return ParticleForCrossSection;}

    void setEnergiesForCrossSectionValues(G4double n ){ EnergiesForCrossSectionValues.push_back(n); }
    std::vector<G4double> getEnergiesForCrossSectionValues() const {return EnergiesForCrossSectionValues;}


    //////////////////////////////////////////// end Physics Commands

private:

#ifdef G4MPI_USE
    //G4MPImanager* g4MPI1 ;
#endif

    G4TStlToGdml* createStlVol ;

    //G4TPhantomParameterisation* param;
    //G4TPhantomNestedParameterisation* param1;

    //G4TVolumeBuilderUsingGDML* G4TVolumeBuilderUsingGDMLObject;
    //G4TVolumeBuilderUsingTEXTGeomtry * G4TVolumeBuilderUsingTEXTGeomtryObject;

    G4TPointDataGeneration* PDG;

    G4TMessenger* OtherDataMessenger;
    G4TGeometryMessenger* geometryMessenger;
    G4TPhysicsMessenger* PhysicsMessenger;

    G4TPrimaryGeneratorMessenger* PGAMessenger;



};

#endif
