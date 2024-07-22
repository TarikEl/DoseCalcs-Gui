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

//#include <stdexcept>

//#include "globals.hh"
//#include "G4SystemOfUnits.hh"

#include "G4TVolumeConstruction.hh"
#include "G4TUserPhysicsList.hh"
#include "G4TNeutronPhysicsList.hh"

#include "G4TActionInitialization.hh"

#include "G4PhysListFactory.hh"

#include "G4RunManagerFactory.hh"

//#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
//#endif

//#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#include "G4UIterminal.hh"
//#endif

#include "TETPhysicsList.hh"
#include "TETActionInitialization.hh"

#ifdef G4MPI_USE
#include "mpi.h"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
//#include "G4RunManager.hh"
#else
#ifdef G4MULTITHREADED
//#include "G4MTRunManager.hh"
#else
//#include "G4RunManager.hh"
#endif
#endif

#include "G4UImanager.hh"

#include <dirent.h>
#include <libgen.h>         // dirname
#include <unistd.h>         // readlink
#include <linux/limits.h>   // PATH_MAX


#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"

size_t* MateIDs; // index of material of each voxel unsigned int* fMateIDs; // index of material of each voxel
G4String* CopyNumberRegionNameMap;
std::map<G4String, G4Colour> RegionNameColour;

G4float* CopyNumberXPos;
G4float* CopyNumberYPos;
G4float* CopyNumberZPos;
G4float* CopyNumberMassSize;
double* CumulativeActivities;

std::map<unsigned int,G4double* >      EnergyListForCriticality;
std::map<unsigned int,G4ThreeVector* > PositionsListForCriticality;
std::map<unsigned int,G4ParticleMomentum* > MomDirecsListForCriticality;
std::map<G4int,std::map<G4int,std::vector<G4double>>> FissionCapturesOfThreadsRanks;
std::map<G4int,std::map<G4String,std::map<G4String,G4int>>> OpticalPhotonInteractionRate;
std::map<G4int,std::map<G4int,G4bool>> TerminatedThreadBatch;
std::vector<G4double> KeffectiveInEachBatch;

G4int NumberOfEventInBatch;
G4int NumberOfBatch;

G4double AllGeometryVolume;
G4double AllGeometryMass;

G4String GeometryFileType;

G4String ParticleName;
G4String SourceType;
G4String SourceRegionName;
G4ThreeVector boxCenterPos;
G4int NumberOfGenPointsToSave;
G4ThreeVector BoxDimGene;
G4ThreeVector newOrgPos3Vector;
G4ThreeVector SourcePosition;
G4ThreeVector SourceRotation;
G4String SourceSolid;
G4String SourceSurface;
G4String SourcePlane;
G4String SourceAxis;
G4double Radius;
G4double RadiusIn;
G4double BeamSDev;
G4double RotTheta;
G4double RotPhi;
G4String RotPosAxis;
G4double HalfX, HalfY, HalfZ;
G4double ThetaMin;
G4double ThetaMax;
G4double PhiMin;
G4double PhiMax;
G4String DirectedParallelAxis;
G4String SimulationIntExtNeutDet;

G4double ToVolumeX;
G4double ToVolumeY;
G4double ToVolumeZ;
G4double DirectedToX;
G4double DirectedToY;
G4double DirectedToZ;
G4String MomDirDirectedHow; // Volume, Point, ParallelToAxis
G4ThreeVector SourceRotVector1;
G4ThreeVector SourceRotVector2;

G4String EnergyDistribution;
G4double GaussSDev;
G4double GaussMean;
G4double UniformEmin;
G4double UniformEmax;
G4double RayleighEmax;
G4double MonoEnergy;
G4double SpectrumMaxEnergy;
G4double FileEnergyCharacterizer;
G4double RadioNuclideMaxEnergy;

std::map<G4double, G4double> EnergyValueProbability;
std::vector<double> RadioNuclideProbVec;
std::vector<unsigned int> RadioNuclidePartNameVec;
std::vector<double*> RadioNuclideEneVec;
std::vector<unsigned int> RadioNuclideSpectrumOrDiscreteVec;

bool UseVoxelsColour;
bool ForceSolid;


G4String MomDirDistribution;
G4double Theta;
G4double Phi;

G4String GeneratePosFlag ;
G4String GenerateEneFlag ;
G4String GenerateMomDirFlag ;
G4String UseGeneratedData ;
G4String ShowBox ;
G4String TestPointsPositions ;

G4String PositionDataFile, EnergyDataFile, MomDirDataFile;

G4String ScriptsDirectoryPath, DataDirectoryPath;

// Run And Score
G4String MPISimulationNum;
G4int BatchsNumber;
G4int ThreadsNumber;
G4String AccuracyCalculationLevel;
G4String organs_to_score ;
G4String variable_To_Score ;

// Voxelized Geometry
G4bool VOXTET_USE;
G4int ParamType;

G4String DataFilesExtension = ".bin";

G4double VoxXHalfSize;
G4double VoxYHalfSize;
G4double VoxZHalfSize;
G4int VoxXNumber;
G4int VoxYNumber;
G4int VoxZNumber;
G4int VoxPETXNumber;
G4int VoxPETYNumber;
G4int VoxPETZNumber;
G4double VoxPETXHalfSize;
G4double VoxPETYHalfSize;
G4double VoxPETZHalfSize;
G4bool UseDicomCumAct;
G4String GenerateVoxelsResuls;

G4String GeometrySymbol ;

std::vector<G4String> RegionsToVisualize;

G4String PlanesToVisualize;
G4int MinPlaneID;
G4int MaxPlaneID;
G4int TotalNumberOfSimulations;

G4String TETNodeDataFile, TETEleDataFile, TETMatDataFile;
std::map<G4int,G4String> MaterialIDName;

G4double MinTETPhantom;
G4double MaxTETPhantom;

bool MaterialNameAsRegionName;

double* SourceEnergies;
double* SourcePositions;
double* SourceMomDirs;


G4String CutInRangeData;
G4String EnergyThresholdsData;
G4bool IsEcutsSet;
G4bool IsDcutsSet;
G4double CutsEnergy;
G4double CutsDistance;
G4String ParticlePysics;
G4String PhotoElectricEffectModel;
G4String PolarizedPhotoElectricEffectModel;
G4String ComptonScatteringModel;
G4String PolarizedComptonModel;
G4String GammaConversionModel;
G4String PolarizedGammaConversionModel;
G4String RayleighScatteringModel;
G4String GammaConversionToMuonModel;
G4String ElectronIonisationModel;
G4String ElectronBremModel;
G4String HadronIonisationModel;
G4String HadronBremModel;
G4String IonIonisationModel;
G4bool GenerateCrossSectionTableFlag;
G4String ParticleForCrossSection;
std::vector<G4double> EnergiesForCrossSectionValues;

//G4Navigator* aNavigator;

G4VPhysicalVolume*  WorldPhysicalVolume;

std::vector<G4String> NewRankSourceParticlesNamesValues;
std::vector<G4String> NewRankSourceRegionsNamesValues;
std::vector<G4ThreeVector> NewRankSourceRegionsBoxDimValues;
std::vector<G4ThreeVector> NewRankSourceRegionsPosValues;
std::vector<std::vector<unsigned int>> NewRankVoxelsIDsOfSourceRegion;
std::vector<std::vector<G4Tet*>> NewRankTETOfSourceRegion;
std::vector<G4ThreeVector> NewRankTETBoxMinOfSourceRegion;
std::vector<G4ThreeVector> NewRankTETBoxDimOfSourceRegion;

std::vector<G4double> NewRankSourceEnergiesValues;
std::vector<G4String> NewRankSourceMomDirsValues;
std::vector<G4double> NewRankSourceMomDirsDirectedThetaValues;
std::vector<G4double> NewRankSourceMomDirsDirectedPhiValues;
std::vector<G4String> NewRankSourcePositionDataFiles;
std::vector<G4String> NewRankSourceEnergyDataFiles;
std::vector<G4String> NewRankSourceMomDirDataFiles;
std::vector<G4String> SourceRegionsNamesToBeIgnoredValues;
std::vector<G4Navigator*> NavigatorForVolumesInitialPosition;

std::vector<G4String> RadioNuclideParticleNames;

G4int EventsNumPerThreadRank;

G4String ResultDirectoryPath ;

unsigned int* ParNameList ;
double* EnergyList        ;
double* MomDirXList       ;
double* MomDirYList       ;
double* MomDirZList       ;
double* PosXList          ;
double* PosYList          ;
double* PosZList          ;


std::string getProgpath()
{
    char result[ PATH_MAX ];
    ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
    return std::string( result, (count > 0) ? count : 0 );
}

G4String MacrosStartingFile;

char result[PATH_MAX];
ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
std::string getProgDirpath()
{
    std::string path;
    if (count != -1) { path = dirname(result); }
    return path;
}

std::string getFileNameFromPath(std::string const & path, std::string const & delims = "/\\") {
    std::string const & filename = path.substr(path.find_last_of(delims) + 1);
    typename std::string::size_type const p(filename.find_last_of('.'));
    return p > 0 && p != std::string::npos ? filename.substr(0, p) : filename;
}

std::string getFileExt(const std::string& s) {
    size_t i = s.rfind('.', s.length());
    if (i != std::string::npos) {
        return(s.substr(i+1, s.length() - i));
    }
    return("");
}

bool DirectoryExists( const char* pzPath )
{
    if ( pzPath == NULL) return false;

    DIR *pDir;
    bool bExists = false;

    pDir = opendir (pzPath);

    if (pDir != NULL)
    {
        bExists = true;
        (void) closedir (pDir);
    }

    return bExists;
}

G4String appBuildDir;

int main(int argc,char** argv){

    // ./core B inputFile.mac 10000

    appBuildDir = getProgDirpath().c_str(); //(getenv("PWD"));
    ResultDirectoryPath = appBuildDir+"/Results";

    G4cout << "\n========= DoseCalcs Build Directory "<< appBuildDir << " ========= "<< "\n" << G4endl;

    G4String InteractiveMode = "B";
    MacrosStartingFile = "inputFile.mac" ;
    G4String MacrosVisFile = "v" ; //"openGLVis.mac" ;
    EventsNumPerThreadRank = 0;

    for (int fd = 0 ;fd < argc;fd++) {
        G4cout << "################# : " << fd+1 << "/" << argc << " "<< argv[fd] << G4endl;
    }

    if(argc > 4){

        InteractiveMode = argv[1];
        MacrosStartingFile = argv[2];
        MacrosVisFile = argv[3];
        EventsNumPerThreadRank = atoi(argv[3]);
        if( strcmp(InteractiveMode.c_str(),"B") != 0 && strcmp(InteractiveMode.c_str(),"G") != 0 && strcmp(InteractiveMode.c_str(),"T") != 0 && strcmp(InteractiveMode.c_str(),"Gen") != 0){   // another value is entered for mode1
            InteractiveMode = "B" ;
        }
    }


    if(argc == 4){

        InteractiveMode = argv[1];
        MacrosStartingFile = argv[2];
        if(InteractiveMode == "G"){ MacrosVisFile = argv[3];G4cout << "################# MacrosVisFile         : " << MacrosVisFile << G4endl;}
        else { EventsNumPerThreadRank = atoi(argv[3]);}

        if( strcmp(InteractiveMode.c_str(),"B") != 0 && strcmp(InteractiveMode.c_str(),"G") != 0 && strcmp(InteractiveMode.c_str(),"T") != 0 && strcmp(InteractiveMode.c_str(),"Gen") != 0){   // another value is entered for mode1
            InteractiveMode = "B" ;
        }
    }
    if(argc == 3){

        InteractiveMode = argv[1];
        MacrosStartingFile = argv[2];
        if(InteractiveMode == "G"){ MacrosVisFile = "v";}
        else { EventsNumPerThreadRank = 0;}

        if( strcmp(InteractiveMode.c_str(),"B") != 0 && strcmp(InteractiveMode.c_str(),"G") != 0 && strcmp(InteractiveMode.c_str(),"T") != 0 && strcmp(InteractiveMode.c_str(),"Gen") != 0){   // another value is entered for mode1
            InteractiveMode = "B" ;
        }
    }
    else if(argc == 2){

        InteractiveMode = argv[1] ;
        MacrosStartingFile = "inputFile.mac";
        EventsNumPerThreadRank = 0;
        MacrosVisFile = "v";
        if( strcmp(InteractiveMode.c_str(),"B") != 0 && strcmp(InteractiveMode.c_str(),"G") != 0 && strcmp(InteractiveMode.c_str(),"T") != 0 && strcmp(InteractiveMode.c_str(),"Gen") != 0){
            InteractiveMode = "B" ;
        }
    }
    else if(argc == 1){

        InteractiveMode = "B" ;
        MacrosStartingFile = "inputFile.mac";
        EventsNumPerThreadRank = 0;
        MacrosVisFile = "v";
    }

    G4cout << "################# InteractiveMode       : " << InteractiveMode << G4endl;
    G4cout << "################# MacrosStartingFile    : " << MacrosStartingFile << G4endl;
    G4cout << "################# EventsNumPerThreadRank: " << EventsNumPerThreadRank << G4endl;

    // Choose the Random engine
    //G4long seed = 100001;
    //G4Random::setTheEngine(new CLHEP::MixMaxRng());
    //G4Random::setTheEngine(new CLHEP::RanecuEngine());
    //G4long seed = time(NULL); G4Random::setTheSeed(seed);

    G4Random::setTheEngine(new CLHEP::MixMaxRng);
    unsigned ins = unsigned(clock());
    G4Random::setTheSeed(ins);
    G4cout << "Seed " << ins << G4endl;

    if(InteractiveMode == "G"){

#ifdef G4MPI_USE
        G4cout << "################# No Graphical run in MPI mode ..." << G4endl;
#else

        auto* runManager = G4RunManagerFactory::CreateRunManager();

        G4TVolumeConstruction* VolCon = new G4TVolumeConstruction;
        UseVoxelsColour = true;

        runManager->SetUserInitialization(VolCon);

        G4UImanager* UImanager = G4UImanager::GetUIpointer();
        UImanager->SetIgnoreCmdNotFound(true);
        UImanager->SetVerboseLevel(0);
        UImanager->ExecuteMacroFile(MacrosStartingFile.c_str());

        if(ParticlePysics.contains("FACTORY")){
        //if(ParticleName == "neutron" && ParticlePysics.contains("FACTORY")){
            std::vector<std::string> RefPhyLists;
            RefPhyLists.push_back("FACTORY_FTFP_BERT"); RefPhyLists.push_back("FACTORY_FTFP_BERT_ATL"); RefPhyLists.push_back("FACTORY_FTFP_BERT_TRV"); RefPhyLists.push_back("FACTORY_QGSP_FTFP_BERT"); RefPhyLists.push_back("FACTORY_QGSP_BERT"); RefPhyLists.push_back("FACTORY_QGSP_BERT_HP"); RefPhyLists.push_back("FACTORY_QGSP_BIC"); RefPhyLists.push_back("FACTORY_QGSP_BIC_AllHP"); RefPhyLists.push_back("FACTORY_Shielding"); RefPhyLists.push_back("FACTORY_ShieldingLEND");
            //RefPhyLists.push_back("FACTORY_INCLXX");
            bool isIn = false; for ( int df = 0 ; df < RefPhyLists.size(); df++  ){ if(ParticlePysics == RefPhyLists[df] ){ isIn = true; }}

            if(isIn == false ){
                ParticlePysics = "FACTORY_QGSP_BERT_HP";
                G4Exception("Reference Physics List", "InvalidSetup", JustWarning, "Invalid reference physics choice, the QGSP_BERT_HP will be used.");
            }

            G4String PhyyName = ParticlePysics;
            std::string s = "FACTORY_";
            std::string::size_type i = PhyyName.find(s);
            if (i != std::string::npos){PhyyName.erase(i, s.length());}

            //std::cout << " ParticlePysics " << ParticlePysics << " PhyyName " << PhyyName << std::endl;

            G4PhysListFactory* physFactory = new G4PhysListFactory();
            G4VModularPhysicsList* physicsList = physFactory->GetReferencePhysList(PhyyName);

            if(CutInRangeData == ""){
                physicsList->SetCutsWithDefault();
            }else{
                std::istringstream LineString(CutInRangeData);
                //G4cout << CutInRangeData << G4endl;

                G4double Ene; G4String pn, Unit;
                while(LineString >> pn ){
                    LineString >> Ene;
                    LineString >> Unit;
                    //G4cout << " pn " << pn << " Ene " << Ene << " Unit " << Unit << G4endl;
                    physicsList->SetCutValue(Ene*G4UnitDefinition::GetValueOf(Unit), pn);
                }
            }
            if(EnergyThresholdsData != ""){
                std::istringstream LineString(EnergyThresholdsData);
                //G4cout << EnergyThresholdsData << G4endl;
                std::string ww; std::vector<std::string> sl;

                while(LineString >> ww ){
                    sl.push_back(ww);
                }

                G4double lowLimit = 1*keV;
                G4double highLimit = 100. * GeV;
                if(sl.size()>1){
                    G4double lowLimit = std::stod(sl[0])*G4UnitDefinition::GetValueOf(sl[1]);
                    if(sl.size()>3){
                        highLimit = std::stod(sl[2])*G4UnitDefinition::GetValueOf(sl[3]);
                    }
                    //G4cout << CutInRangeData << " - lowLimit " << lowLimit << " highLimit " << highLimit << G4endl;
                    if(lowLimit != 0. && highLimit > lowLimit ){
                        G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);
                    }
                }
            }
            runManager->SetUserInitialization(physicsList);
        }
        else if(ParticlePysics.contains("OpticalPhysics")){

            G4PhysListFactory* physFactory = new G4PhysListFactory();
            G4VModularPhysicsList* physicsList = physFactory->GetReferencePhysList("FTFP_BERT");
            physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());

            G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
            auto opticalParams               = G4OpticalParameters::Instance();

            opticalParams->SetWLSTimeProfile("delta");

            opticalParams->SetScintTrackSecondariesFirst(true);

            opticalParams->SetCerenkovMaxPhotonsPerStep(100);
            opticalParams->SetCerenkovMaxBetaChange(10.0);
            opticalParams->SetCerenkovTrackSecondariesFirst(true);

            physicsList->RegisterPhysics(opticalPhysics);
            runManager->SetUserInitialization(physicsList);
        }
        else{
            runManager->SetUserInitialization(new G4TUserPhysicsList());
        }

        // depend on the input to the VolumeConstructor then it called after ExecuteMacroFile(MacrosStartingFile.c_str());
        runManager->SetUserInitialization(new G4TActionInitialization());

#ifdef G4MULTITHREADED
        runManager->SetNumberOfThreads(1);
        // For B and T modes we need a normal simulation without Box or Points Testing
#endif

        //#ifdef G4VIS_USE
        //#endif

        runManager->Initialize();
        //runManager->SetEventModulo(100);
        //runManager->BeamOn(400);

        G4cout << "\n\n========= Graphical session starts for visualization ======================= \n"<<G4endl;

        //G4cout << "\n\n\nTest Points " << VolCon->getTestPointsPositions()<< G4endl;

        //#ifdef G4UI_USE

        if(MacrosVisFile == "v"){

            G4VisExecutive* visManager = new G4VisExecutive;
            visManager->Initialize();

            G4UIExecutive* ui = new G4UIExecutive(argc, argv);

            G4String Plane;
            G4String ThetaPhi = "";

            if(argc > 4){
                Plane = argv[4];
            }else{
                Plane = PlanesToVisualize;
            }

            Plane.toLower();

            G4cout << " Plane " << Plane << G4endl;

            if(Plane == "xy" || Plane == "yx"){ //
                ThetaPhi = " 0 0 deg";
                G4cout << " ThetaPhi " << ThetaPhi << G4endl;
            }
            else if(Plane == "xz" || Plane == "zx"){
                ThetaPhi = " 90 90 deg";
                G4cout << " ThetaPhi " << ThetaPhi << G4endl;
            }
            else if(Plane == "zy" || Plane == "yz"){//
                ThetaPhi = " 90 0 deg";
                G4cout << " ThetaPhi " << ThetaPhi << G4endl;
            }
            else{
                ThetaPhi = " 90 0 deg";
                G4cout << " Using default theta and phi " << ThetaPhi << G4endl;
            }

            UImanager->ApplyCommand("/vis/open OGL");
            UImanager->ApplyCommand("/vis/viewer/set/autoRefresh false");
            UImanager->ApplyCommand("/vis/verbose errors");
            UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi "+ ThetaPhi);
            UImanager->ApplyCommand("/vis/drawVolume");
            UImanager->ApplyCommand("/vis/scene/add/axes 0 0 0 10 cm");
            UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
            UImanager->ApplyCommand("/vis/modeling/trajectories/create/drawByCharge");
            UImanager->ApplyCommand("/vis/modeling/trajectories/drawByCharge-0/default/setDrawStepPts true");
            UImanager->ApplyCommand("/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2");
            UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
            UImanager->ApplyCommand("/vis/viewer/set/autoRefresh true");
            UImanager->ApplyCommand("/vis/verbose warnings");
            UImanager->ApplyCommand("/vis/viewer/flush");

            ui->SessionStart();
            delete ui;
        }
        else if(MacrosVisFile == "d"){

            G4String Plane;
            G4String ThetaPhi = "";

            if(argc > 4){
                Plane = argv[4];
            }else{
                Plane = PlanesToVisualize;
            }

            Plane.toLower();

            G4cout << " Plane " << Plane << G4endl;

            if(Plane == "xy" || Plane == "yx"){ //
                ThetaPhi = " 0 0 deg";
                G4cout << " ThetaPhi " << ThetaPhi << G4endl;
            }
            // reset the command "/vis/viewer/set/viewpointThetaPhi 90 91 deg" in Qt viewer
            // to view the right angle
            else if(Plane == "xz" || Plane == "zx"){
                ThetaPhi = " 90 91 deg ";
                G4cout << " ThetaPhi " << ThetaPhi << G4endl;
            }
            else if(Plane == "zy" || Plane == "yz"){//
                ThetaPhi = " 90 0 deg";
                G4cout << " ThetaPhi " << ThetaPhi << G4endl;
            }
            else{
                ThetaPhi = " 90 0 deg";
                G4cout << " Using default theta and phi " << ThetaPhi << G4endl;
            }

            G4cout << "\n" << G4endl;

            G4VisExecutive* visManager = new G4VisExecutive;
            visManager->Initialize();

            // necessite installation of dawn, install dawn // all steps are described here : https://twiki.cern.ch/twiki/bin/view/CLIC/DawnVisualization
            // in terminal export G4DAWNFILE_VIEWER="dawn -d"

            //UImanager->ApplyCommand("/vis/scene/add/axes 1 1 1 10 cm");

            UImanager->ApplyCommand("/vis/open DAWNFILE");
            //UImanager->ApplyCommand("/vis/open OGLSX");

            UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi "+ ThetaPhi);

            if(PlanesToVisualize != "all"){UImanager->ApplyCommand("/vis/drawVolume TETVoxelContainer");}
            else{UImanager->ApplyCommand("/vis/drawVolume World");}

            //UImanager->ApplyCommand("/vis/viewer/set/upVector "+ ThetaPhi);

            UImanager->ApplyCommand("/vis/viewer/flush");
            UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi "+ ThetaPhi);
        }
        else {
            UImanager->ApplyCommand("/control/execute "+MacrosVisFile); // this is were added after ExecuteMacroFile(MacrosStartingFile.c_str()) to execute according to the user choice graphical or basic without changing the input.mac by adding or removing the command /control/execute openGLVis.mac
            //delete visManager;
        }

        //#endif
        delete runManager;
#endif
        return 0;
    }
    if(InteractiveMode == "Gen"){

        // or make it run Of seq for all modes
        G4cout << "################# Gen-in-SEQ mode ..." << G4endl;
        //G4RunManager* runManager = new G4RunManager;

        auto* runManager = G4RunManagerFactory::CreateRunManager();

        G4TVolumeConstruction* VolCon= new G4TVolumeConstruction();
        VolCon->setMPISimulationNum("m");
        VolCon->setShowBox("no");

        if(VolCon->getGeometryFileType() != "TET"){ // because we TET require a high computational materials, then in all cases we will simulate just
                                                    // a part by setting the command "/GeometryData/setTETPhantomLimits xy -15 15"
            PlanesToVisualize = "all";
        }

        runManager->SetUserInitialization(VolCon);
        G4UImanager* UImanager = G4UImanager::GetUIpointer();

        UImanager->ExecuteMacroFile(MacrosStartingFile.c_str());
        if(ParticlePysics.contains("FACTORY")){
        //if(ParticleName == "neutron" && ParticlePysics.contains("FACTORY")){
            std::vector<std::string> RefPhyLists;
            RefPhyLists.push_back("FACTORY_FTFP_BERT"); RefPhyLists.push_back("FACTORY_FTFP_BERT_ATL"); RefPhyLists.push_back("FACTORY_FTFP_BERT_TRV"); RefPhyLists.push_back("FACTORY_QGSP_FTFP_BERT"); RefPhyLists.push_back("FACTORY_QGSP_BERT"); RefPhyLists.push_back("FACTORY_QGSP_BERT_HP"); RefPhyLists.push_back("FACTORY_QGSP_BIC"); RefPhyLists.push_back("FACTORY_QGSP_BIC_AllHP"); RefPhyLists.push_back("FACTORY_Shielding"); RefPhyLists.push_back("FACTORY_ShieldingLEND");
            //RefPhyLists.push_back("FACTORY_INCLXX");
            bool isIn = false; for ( int df = 0 ; df < RefPhyLists.size(); df++  ){ if(ParticlePysics == RefPhyLists[df] ){ isIn = true; }}

            if(isIn == false ){
                ParticlePysics = "FACTORY_QGSP_BERT_HP";
                G4Exception("Reference Physics List", "InvalidSetup", JustWarning, "Invalid reference physics choice, the QGSP_BERT_HP will be used.");
            }

            G4String PhyyName = ParticlePysics;
            std::string s = "FACTORY_";
            std::string::size_type i = PhyyName.find(s);
            if (i != std::string::npos){PhyyName.erase(i, s.length());}

            //std::cout << " ParticlePysics " << ParticlePysics << " PhyyName " << PhyyName << std::endl;

            G4PhysListFactory* physFactory = new G4PhysListFactory();
            G4VModularPhysicsList* physicsList = physFactory->GetReferencePhysList(PhyyName);

            if(CutInRangeData == ""){
                physicsList->SetCutsWithDefault();
            }else{
                std::istringstream LineString(CutInRangeData);
                //G4cout << CutInRangeData << G4endl;

                G4double Ene; G4String pn, Unit;
                while(LineString >> pn ){
                    LineString >> Ene;
                    LineString >> Unit;
                    //G4cout << " pn " << pn << " Ene " << Ene << " Unit " << Unit << G4endl;
                    physicsList->SetCutValue(Ene*G4UnitDefinition::GetValueOf(Unit), pn);
                }
            }
            if(EnergyThresholdsData != ""){
                std::istringstream LineString(EnergyThresholdsData);
                //G4cout << EnergyThresholdsData << G4endl;
                std::string ww; std::vector<std::string> sl;

                while(LineString >> ww ){
                    sl.push_back(ww);
                }

                G4double lowLimit = 1*keV;
                G4double highLimit = 100. * GeV;
                if(sl.size()>1){
                    G4double lowLimit = std::stod(sl[0])*G4UnitDefinition::GetValueOf(sl[1]);
                    if(sl.size()>3){
                        highLimit = std::stod(sl[2])*G4UnitDefinition::GetValueOf(sl[3]);
                    }
                    //G4cout << CutInRangeData << " - lowLimit " << lowLimit << " highLimit " << highLimit << G4endl;
                    if(lowLimit != 0. && highLimit > lowLimit ){
                        G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);
                    }
                }
            }
            runManager->SetUserInitialization(physicsList);
        }
        else if(ParticlePysics.contains("OpticalPhysics")){

            G4PhysListFactory* physFactory = new G4PhysListFactory();
            G4VModularPhysicsList* physicsList = physFactory->GetReferencePhysList("FTFP_BERT");
            physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());

            G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
            auto opticalParams               = G4OpticalParameters::Instance();

            opticalParams->SetWLSTimeProfile("delta");

            opticalParams->SetScintTrackSecondariesFirst(true);

            opticalParams->SetCerenkovMaxPhotonsPerStep(100);
            opticalParams->SetCerenkovMaxBetaChange(10.0);
            opticalParams->SetCerenkovTrackSecondariesFirst(true);

            physicsList->RegisterPhysics(opticalPhysics);
            runManager->SetUserInitialization(physicsList);
        }
        else{
            runManager->SetUserInitialization(new G4TUserPhysicsList());
        }

        //G4VModularPhysicsList* physicsList = new G4TNeutronPhysicsList("NeutronHP");
        //runManager->SetUserInitialization(physicsList);

        // depend on the input to the VolumeConstructor then it called after ExecuteMacroFile(MacrosStartingFile.c_str());
        runManager->SetUserInitialization(new G4TActionInitialization());
        runManager->Initialize();

        VolCon->setUseGeneratedData("generate");
        VolCon->setGeneratePositions("yes");
        VolCon->setGenerateEnergies("yes");
        VolCon->setGenerateMomDirs("yes");
        VolCon->setPointNumberToSave(EventsNumPerThreadRank);

        VolCon->GenerateEventsDataInMPIMode();

        delete runManager;
        return 0;
    }


#ifdef G4MPI_USE
    //MPI_Init(&argc,&argv);
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if(provided < MPI_THREAD_MULTIPLE){MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);/*printf("The threading support level is lesser than that demanded.\n");*/}
    else{/*printf("The threading support level corresponds to that demanded.\n");   */ }
    int DataID;
    MPI_Comm_rank(MPI_COMM_WORLD, &DataID);
    //G4MPImanager* g4MPI = new G4MPImanager(/*argc,argv*/); // if we use the parameter passed from terminal (argv) then the constructor must be empty
    //G4MPIsession* session = g4MPI-> GetMPIsession();
#endif

    UseVoxelsColour = false;

    //G4RunManager* runManager = new G4RunManager;
    //auto* runManager = G4RunManagerFactory::CreateRunManager();
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

    G4TVolumeConstruction* VolCon= new G4TVolumeConstruction();
    runManager->SetUserInitialization(VolCon);

    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    UImanager->SetIgnoreCmdNotFound(true);
    UImanager->SetVerboseLevel(0);
    UImanager->ExecuteMacroFile(MacrosStartingFile.c_str());

    if(ParticlePysics.contains("FACTORY")){
    //if(ParticleName == "neutron" && ParticlePysics.contains("FACTORY")){
        std::vector<std::string> RefPhyLists;
        RefPhyLists.push_back("FACTORY_FTFP_BERT"); RefPhyLists.push_back("FACTORY_FTFP_BERT_ATL"); RefPhyLists.push_back("FACTORY_FTFP_BERT_TRV"); RefPhyLists.push_back("FACTORY_QGSP_FTFP_BERT"); RefPhyLists.push_back("FACTORY_QGSP_BERT"); RefPhyLists.push_back("FACTORY_QGSP_BERT_HP"); RefPhyLists.push_back("FACTORY_QGSP_BIC"); RefPhyLists.push_back("FACTORY_QGSP_BIC_AllHP"); RefPhyLists.push_back("FACTORY_Shielding"); RefPhyLists.push_back("FACTORY_ShieldingLEND");
        //RefPhyLists.push_back("FACTORY_INCLXX");
        bool isIn = false; for ( int df = 0 ; df < RefPhyLists.size(); df++  ){ if(ParticlePysics == RefPhyLists[df] ){ isIn = true; }}

        if(isIn == false ){
            ParticlePysics = "FACTORY_QGSP_BERT_HP";
            G4Exception("Reference Physics List", "InvalidSetup", JustWarning, "Invalid reference physics choice, the QGSP_BERT_HP will be used.");
        }

        G4String PhyyName = ParticlePysics;
        std::string s = "FACTORY_";
        std::string::size_type i = PhyyName.find(s);
        if (i != std::string::npos){PhyyName.erase(i, s.length());}

        //std::cout << " ParticlePysics " << ParticlePysics << " PhyyName " << PhyyName << std::endl;

        G4PhysListFactory* physFactory = new G4PhysListFactory();
        G4VModularPhysicsList* physicsList = physFactory->GetReferencePhysList(PhyyName);

        if(CutInRangeData == ""){
            physicsList->SetCutsWithDefault();
        }else{
            std::istringstream LineString(CutInRangeData);
            G4cout << " CutInRangeData " << CutInRangeData << G4endl;

            G4double Ene; G4String pn, Unit;
            while(LineString >> pn ){
                LineString >> Ene;
                LineString >> Unit;
                //G4cout << " pn " << pn << " Ene " << Ene << " Unit " << Unit << G4endl;
                physicsList->SetCutValue(Ene*G4UnitDefinition::GetValueOf(Unit), pn);
            }
        }
        if(EnergyThresholdsData != ""){
            std::istringstream LineString(EnergyThresholdsData);
            //G4cout << EnergyThresholdsData << G4endl;
            std::string ww; std::vector<std::string> sl;

            while(LineString >> ww ){
                sl.push_back(ww);
            }

            G4double lowLimit = 1*keV;
            G4double highLimit = 100. * GeV;
            if(sl.size()>1){
                G4double lowLimit = std::stod(sl[0])*G4UnitDefinition::GetValueOf(sl[1]);
                if(sl.size()>3){
                    highLimit = std::stod(sl[2])*G4UnitDefinition::GetValueOf(sl[3]);
                }
                //G4cout << CutInRangeData << " - lowLimit " << lowLimit << " highLimit " << highLimit << G4endl;
                if(lowLimit != 0. && highLimit > lowLimit ){
                    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);
                }
            }
        }
        runManager->SetUserInitialization(physicsList);
    }
    else if(ParticlePysics.contains("OpticalPhysics")){

        G4PhysListFactory* physFactory = new G4PhysListFactory();
        G4VModularPhysicsList* physicsList = physFactory->GetReferencePhysList("FTFP_BERT");
        physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());

        G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
        auto opticalParams               = G4OpticalParameters::Instance();

        opticalParams->SetWLSTimeProfile("delta");

        opticalParams->SetScintTrackSecondariesFirst(true);

        opticalParams->SetCerenkovMaxPhotonsPerStep(100);
        opticalParams->SetCerenkovMaxBetaChange(10.0);
        opticalParams->SetCerenkovTrackSecondariesFirst(true);

        physicsList->RegisterPhysics(opticalPhysics);
        runManager->SetUserInitialization(physicsList);
    }
    else {

        //G4PhysListFactory* physFactory = new G4PhysListFactory();
        //G4VModularPhysicsList* physicsList = physFactory->GetReferencePhysList("FTFP_BERT");
        //physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
        ////physicsList->RemovePhysics(new G4EmStandardPhysics_option4());
        //G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
        //physicsList->RegisterPhysics(opticalPhysics);
        //runManager->SetUserInitialization(physicsList);
        runManager->SetUserInitialization(new G4TUserPhysicsList());
    }

/*
    if(ParticlePysics == "EMS" || ParticlePysics == "EMS1" || ParticlePysics == "EMS2" || ParticlePysics == "EMS3"|| ParticlePysics == "EMS4"|| ParticlePysics == "Livermore"|| ParticlePysics == "Penelope"){
        if(ParticleName == "neutron"){
            G4VModularPhysicsList* physicsList = new G4TNeutronPhysicsList("NeutronHP");
            runManager->SetUserInitialization(physicsList);
            if(CutInRangeData == ""){
                physicsList->SetCutsWithDefault();
            }else{
                std::istringstream LineString(CutInRangeData);
                G4cout << CutInRangeData << G4endl;

                G4double Ene; G4String pn, Unit;
                while(LineString >> pn ){
                    LineString >> Ene;
                    LineString >> Unit;
                    //G4cout << " pn " << pn << " Ene " << Ene << " Unit " << Unit << G4endl;
                    physicsList->SetCutValue(Ene*G4UnitDefinition::GetValueOf(Unit), pn);
                }
            }
        }else{
            runManager->SetUserInitialization(new G4TUserPhysicsList());
        }
    }else{
        if(ParticleName == "neutron"){
            std::vector<std::string> RefPhyLists;
            RefPhyLists.push_back("FTFP_BERT"); RefPhyLists.push_back("FTFP_BERT_ATL"); RefPhyLists.push_back("FTFP_BERT_TRV"); RefPhyLists.push_back("QGSP_FTFP_BERT"); RefPhyLists.push_back("QGSP_BERT"); RefPhyLists.push_back("QGSP_BERT_HP"); RefPhyLists.push_back("QGSP_BIC"); RefPhyLists.push_back("QGSP_BIC_AllHP"); RefPhyLists.push_back("INCLXX"); RefPhyLists.push_back("Shielding"); RefPhyLists.push_back("ShieldingLEND");
            bool isIn = false; for ( int df = 0 ; df < RefPhyLists.size(); df++  ){ if(ParticlePysics == RefPhyLists[df] ){ isIn = true; }}

            if(isIn == false ){
                ParticlePysics = "QGSP_BERT_HP";
                G4Exception("Reference Physics List", "InvalidSetup", JustWarning, "Invalid reference physics choice, the QGSP_BERT_HP will be used.");
            }

            G4PhysListFactory* physFactory = new G4PhysListFactory();
            G4VModularPhysicsList* physicsList = physFactory->GetReferencePhysList(ParticlePysics);

            if(CutInRangeData == ""){
                physicsList->SetCutsWithDefault();
            }else{
                std::istringstream LineString(CutInRangeData);
                G4cout << CutInRangeData << G4endl;

                G4double Ene; G4String pn, Unit;
                while(LineString >> pn ){
                    LineString >> Ene;
                    LineString >> Unit;
                    //G4cout << " pn " << pn << " Ene " << Ene << " Unit " << Unit << G4endl;
                    physicsList->SetCutValue(Ene*G4UnitDefinition::GetValueOf(Unit), pn);
                }
            }
    if(EnergyThresholdsData != ""){
        std::istringstream LineString(EnergyThresholdsData);
        //G4cout << EnergyThresholdsData << G4endl;
        std::string ww; std::vector<std::string> sl;

        while(LineString >> ww ){
            sl.push_back(ww);
        }

        G4double lowLimit = 1*keV;
        G4double highLimit = 100. * GeV;
        if(sl.size()>1){
            G4double lowLimit = std::stod(sl[0])*G4UnitDefinition::GetValueOf(sl[1]);
            if(sl.size()>3){
                highLimit = std::stod(sl[2])*G4UnitDefinition::GetValueOf(sl[3]);
            }
            //G4cout << CutInRangeData << " - lowLimit " << lowLimit << " highLimit " << highLimit << G4endl;
            if(lowLimit != 0. && highLimit > lowLimit ){
                G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);
            }
        }
    }

            runManager->SetUserInitialization(physicsList);
        }else{
            ParticlePysics = "EMS3";
            runManager->SetUserInitialization(new G4TUserPhysicsList());
        }
    }
*/

    G4int par = ParamType ; // because ParamType is unset accidentally !!
    VolCon->setParamType(par);  // because ParamType is unset accidentally !!

    if( GeometryFileType == "TET"){
        runManager->SetUserInitialization(new TETActionInitialization());
    }else {
        runManager->SetUserInitialization(new G4TActionInitialization());
    }

    runManager->SetNumberOfThreads(VolCon->getNumberOfThreads());
    runManager->SetNumberOfEventsToBeProcessed(VolCon->getNumberOfThreads()*EventsNumPerThreadRank);
    static_cast<G4MTRunManager*>(runManager)->SetEventModulo(EventsNumPerThreadRank);

    VolCon->setShowBox("no");
    VolCon->setTestPointsPositions("no");
    if(VolCon->getGeometryFileType() != "TET"){ // because we TET require a high computational materials, then in all cases we will simulate just
        // a part by setting the command "/GeometryData/setTETPhantomLimits xy -15 15"
        PlanesToVisualize = "all";
    }

    VolCon->setPointNumberToSave(EventsNumPerThreadRank); // because we eliminate the NumberOfGenPointsToSave

    std::cout << "\n\n========= Geometry, physics and radiation source data initialization ======================= \n"<<std::endl;

    runManager->Initialize();

    if(InteractiveMode == "B" || InteractiveMode == ""){

#ifdef G4MPI_USE
        std::cout << "\n\n========= MPI Calculation session starts. " << VolCon->getNumberOfThreads() << " ranks, " << EventsNumPerThreadRank << " in each rank ======================= \n"<<std::endl;
#else
        std::cout << "\n\n========= MT Calculation session starts. " << VolCon->getNumberOfThreads() << " threads, " << EventsNumPerThreadRank << " in each thread ======================= \n"<<std::endl;
#endif
        runManager->BeamOn(EventsNumPerThreadRank*VolCon->getNumberOfThreads());
    }
    else if(InteractiveMode == "T" ){

        G4cout << "\n Terminal session starts ..." << G4endl;

        G4UIsession * session = new G4UIterminal;
        session->SessionStart();
        delete session;
    }

    //delete g4MPI;
    delete runManager;

#ifdef G4MPI_USE
    MPI_Finalize();
#endif

    std::cout << "\n Terminated_Session" << std::endl;

    return 0;
}
