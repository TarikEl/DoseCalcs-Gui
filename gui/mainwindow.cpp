#include "gui/mainwindow.h"
#include "gui/ui_mainwindow.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>
#include <cstdlib>
#include <unistd.h>
#include <QProcess>
#include <QDir>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <QFuture>
#include <QtConcurrent/QtConcurrent>
//#include <QThread>


extern QStringList CompleterWords;

extern QString  TestPackagesPathsBeforToRun();

extern int NumberOfCPUCores ;
extern qint64 pidOfPreviousNotKilled;

extern QString OSNameVersion;

extern QString DoseCalcsDownloadURL;
extern QString DoseCalcsCore_build_dir_path;
extern QString UserCurrentResultsDirPath;
extern QString DoseCalcs_build_dir_name;

extern QString ResultDirectoryName;
extern QString ScriptDirectoryName;
extern QString DoseCalcsExecutableName;
extern QString MergeExecutableName;
extern QString GraphExecutableName;
extern QString MacroFileName;
extern QString ResultFileName;
extern QString ReferenceFileName;
extern QString GraphsOutDirName;
extern QString DoseCalcsExecutingFileName;
extern QString MacroFilePath;
extern QString ICRPDATAPath;

extern QString geant4_Lib_dir_path ;
extern QString DCMTK_Lib_dir_path ;
extern QString MPI_Lib_dir_path ;
extern QString CMAKE_Lib_dir_path ;
extern QString Root_Lib_dir_path ;
extern QString DoseCalcsCore_source_dir_path;
extern QString DoseCalcsGui_source_dir_path;

extern QString GUIConfigFileName;
extern QString GUIPackagesAndFilesDirPath;
extern QString GUIPackagesAndFilesDirName;

extern QString Execution_setEventNumber;

extern QString CENTOS_ROCKS_CLUSTER ;
extern QString MPI_USE ;
extern QString ROOT_USE ;
extern QString DCMTK_USE ;
extern QString VERBOSE_USE ;
extern QString GDML_USE;
extern QString ConfigDataText;

extern QString cmakeTruePath;

extern QString QueueName;
extern QString ExeDataText;
extern QString ExeFileName;

extern QVector <QString> MaterialCommands ;
extern QVector <QString> GeometryCommands ;
extern QVector <QString> VOXELCommands ;
extern QVector <QString> DICOMCommands ;
extern QVector <QString> CONSCommands ;
extern QVector <QString> PhysicsCommands ;
extern QVector <QString> SourceCommands ;
extern QVector <QString> RunAndScoreCommands ;
extern QVector <QString> AnalysisCommands ;

extern QStringList DoseCalcsQuantities ;
extern QMap <QString,QStringList> QuantitiesUnitsLists ;
extern QMap <QString,QMap <QString,double>> QuantitiesConversionFromDefault ;

extern std::string getFileNameFromPath(std::string const & path, std::string const & delims = "/\\");
extern std::string getFileExt(const std::string& s);

extern void PrintToConsole(QString text, int level);

// it call initializeVariable(), ReadConfigurationFile(ConfigFile),
MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)  // it look to a ui file form MainWindow to initialize Qwidget created by Qt creator and that will be used by MainWindow class
{

    ui->setupUi(this);

    //ui->groupBoxExternal->setVisible(false);

    QStringList ICRPPhantoms=(QStringList()<<"ICRPAdultMale"<<"ICRPAdultFemale");
    ui->comboBoxPhantom->addItems(ICRPPhantoms);

    QStringList PreDefinedGeom=(QStringList()<<"ICRP110 Voxel-Type Adult Male"
                                <<"ICRP110 Voxel-Type Adult Female"
                                <<"ICRP143 Voxel-Type Male 15-years-old"
                                <<"ICRP143 Voxel-Type Female 15-years-old"
                                <<"ICRP143 Voxel-Type Male 10-years-old"
                                <<"ICRP143 Voxel-Type Female 10-years-old"
                                <<"ICRP143 Voxel-Type Male 05-years-old"
                                <<"ICRP143 Voxel-Type Female 05-years-old"
                                <<"ICRP143 Voxel-Type Male 01-years-old"
                                <<"ICRP143 Voxel-Type Female 01-years-old"
                                <<"ICRP143 Voxel-Type Male Newborn"
                                <<"ICRP143 Voxel-Type Female Newborn"
                                <<"ICRP145 Mesh-Type Adult Male"
                                <<"ICRP145 Mesh-Type Adult Female"
                                <<"Phantom Constructed using DICOM/CT Files"
                                <<"Simple Voxelized Geometry Using DoseCalcs Commands"
                                <<"Stylized MIRD Adult Female Phantom using CPP (../DoseCalcs/core/src/G4TCPPGeometryFormat.cc) Geometry Method"
                                <<"Stylized MIRD Adult Female Phantom using GDML (.gdml) Geometry Method"
                                <<"Stylized MIRD Adult Female Phantom using TEXT (.geom) Geometry Method"
                                <<"Phantom Constructed using a Combination of GDML, CPP, TEXT, STL, and STANDARD methods"
                                <<"Simple Geometry using STL Solids and DoseCalcs Commands"
                                <<"Simple Geometry using Geant4 Standard Solids and DoseCalcs Commands"
                                <<"Simple Voxelized Water Phantom"
                                <<"MyGeometry"
                                );

    PreDefinedGeomMap[PreDefinedGeom[0]] ="VoxICRPAdultMale";
    PreDefinedGeomMap[PreDefinedGeom[1]] ="VoxICRPAdultFemale";
    PreDefinedGeomMap[PreDefinedGeom[2]] ="VoxICRPMale15";
    PreDefinedGeomMap[PreDefinedGeom[3]] ="VoxICRPFemale15";
    PreDefinedGeomMap[PreDefinedGeom[4]] ="VoxICRPMale10";
    PreDefinedGeomMap[PreDefinedGeom[5]] ="VoxICRPFemale10";
    PreDefinedGeomMap[PreDefinedGeom[6]] ="VoxICRPMale05";
    PreDefinedGeomMap[PreDefinedGeom[7]] ="VoxICRPFemale05";
    PreDefinedGeomMap[PreDefinedGeom[8]] ="VoxICRPMale01";
    PreDefinedGeomMap[PreDefinedGeom[9]] ="VoxICRPFemale01";
    PreDefinedGeomMap[PreDefinedGeom[10]]="VoxICRPMale00";
    PreDefinedGeomMap[PreDefinedGeom[11]]="VoxICRPFemale00";
    PreDefinedGeomMap[PreDefinedGeom[12]]="TetICRPAdultMale";
    PreDefinedGeomMap[PreDefinedGeom[13]]="TetICRPAdultFemale";
    PreDefinedGeomMap[PreDefinedGeom[14]]="DICOM_CT_Phantom";
    PreDefinedGeomMap[PreDefinedGeom[15]]="DoseCalcs_Voxelized_Geometry";
    PreDefinedGeomMap[PreDefinedGeom[16]]="Stylized_CPP_MIRD_AdultFemale";
    PreDefinedGeomMap[PreDefinedGeom[17]]="Stylized_GDML_MIRD_AdultFemale";
    PreDefinedGeomMap[PreDefinedGeom[18]]="Stylized_TEXT_MIRD_AdultFemale";
    PreDefinedGeomMap[PreDefinedGeom[19]]="Constructed_GDML_CPP_TEXT_STL_STAN_Example";
    PreDefinedGeomMap[PreDefinedGeom[20]]="STL_Solids_Example";
    PreDefinedGeomMap[PreDefinedGeom[21]]="Geant4_STANDARDs_Example";
    PreDefinedGeomMap[PreDefinedGeom[22]]="SimpleVoxelizedWaterPhantom";
    PreDefinedGeomMap[PreDefinedGeom[23]]="MyGeometry";

    ui->comboBoxPreDefinedGeom->addItems(PreDefinedGeom);

    //QStringList SourceOrTargetOrCombination=(QStringList()<<"Region as a source and target"<<"In a specific target from all sources"<<"For a target<--source combination"<<"In all targets (phantom) from all sources" ); //<<"In a specific targets from a specific sources");
    //ui->comboBoxTotalorSourceOrCombNucl->addItems(SourceOrTargetOrCombination);

    QStringList SourceOrTargetOrCombination1=(QStringList()<<"Region as a source and target"<<"In each target from a specific source"<<"In each target from all sources" ); //<<"In each target from a number of sources");
    ui->comboBoxSourceOrTargetsForBiokinetics->addItems(SourceOrTargetOrCombination1);

    ui->progressBarReadingCalcData->setWindowIconText("Reading");
    ui->progressBarReadingCalcData->setValue(0);

    QStringList DDD=(QStringList()<<"File"<<"P-E-Y");
    ui->comboBoxRadionuclidedataType->addItems(DDD);
    RadDataType["File"]= "File";
    RadDataType["P-E-Y"]= "P-E-Y";

    ui->comboBoxPeriodUnit->addItems(QuantitiesUnitsLists["T"]);
    ui->comboBoxQuantityNucl->addItems(DoseCalcsQuantities);
    ui->comboBoxActivityAdministered->addItems(QuantitiesUnitsLists["A"]);
    ui->comboBoxEffDoseUnit->clear();ui->comboBoxEffDoseUnit->addItems(QuantitiesUnitsLists["AE"]);

    ui->comboBoxRedionuclidesRadEmmiterType->addItems((QStringList()<<"All"<<"Gamma or X-ray emitters"<<"Electron emitters"<<"Positron emitters"<<"Alpha emitters"<<"Neutron emitters"));

    ui->RadPhotonCheckBox->setChecked(true);
    ui->RadelectronCheckBox->setChecked(true);
    ui->RadPositronCheckBox->setChecked(true);
    ui->RadAlphaCheckBox->setChecked(true);
    ui->RadNeutronCheckBox->setChecked(true);
    ui->RadFFCheckBox->setChecked(true);


    // ///////////////////////////////////////

    MaterialsDataCommandsString= "";
    GeometryCommandsString= "";
    DICOMCommandsString= "";
    VoxIDsCommandsString= "";
    VOXELCommandsString= "";
    CONSCommandsString= "";
    PhysicsCommandsString= "";
    AnalysisCommandsString= "";
    RunAndScoreCommandsString= "";
    SourceCommandsString = "";
    MacrosCommandsString = "";

    macrosfileinc=0;

    initializeVariable();
    CommandsInitialization();

    MacroFilePath = GUIPackagesAndFilesDirPath+"/"+GUIConfigFileName;
    ReadConfigurationFile(MacroFilePath);

    //ui->Tab->setVisible(true);
    //ui->tabNJOY->setVisible(false);
    //ui->tabReactor->setVisible(false);
    //ui->tabDetector->setVisible(false);

    ui->checkBoxFixGeometryCommands->setChecked(true);
    ui->checkBoxFixPhysicsCommands->setChecked(true);
    ui->checkBoxFixSourceCommands->setChecked(true);
    ui->checkBoxFixScoreCommands->setChecked(true);
    ui->checkBoxFixRunCommands->setChecked(true);
    ui->checkBoxFixGraphsCommands->setChecked(true);

    QStringList GraphsDatalist=(QStringList()<<"none"<<"Result"<<"Reference"<<"Reference_Result");
    ui->AnalysisComboBoxGraphData->addItems(GraphsDatalist);

    QStringList GraphsTypelist=(QStringList()<<"Self"<<"Cross"<<"Self_Cross");
    ui->AnalysisComboBoxGraphsType->addItems(GraphsTypelist);

    QStringList GraphsExtlist=(QStringList()<<".root"<<".png"<<".jpg"<<".pdf"<<".ps"<<".tex");
    ui->AnalysisComboBoxGraphsExt->addItems(GraphsExtlist);

    //QStringList VarToScorellist=(QStringList()<<"AE"<<"AD"<<"AF"<<"SAF"<<"S"<<"ED"<<"EDE");
    //ui->AnalysisComboBoxVarToScore->addItems(VarToScorellist);


    // source Tab

    QStringList EnergyDistlist=(QStringList()<<"Mono"<<"Uniform"<<"Gauss"<<"Rayleigh"<<"Spectrum"<<"File");
    ui->SourceComboBoxEnergyDist->addItems(EnergyDistlist);

    QStringList AngleDistlist=(QStringList()<<"Isotropic"<<"Uniform"<<"Directed");
    ui->SourceComboBoxAngleDist->addItems(AngleDistlist);

    QStringList Physicslist=(QStringList()<<"EMS"<<"EMS1"<<"EMS2"<<"EMS3"<<"EMS4"<<"Penelope"<<"Livermore"<<"Construct"<<"OpticalPhysics"
                             <<"HADRON_FTFP_BERT"<<"HADRON_FTFP_BERT_ATL"<<"HADRON_FTFP_BERT_TRV"<<"HADRON_QGSP_FTFP_BERT"<<"HADRON_QGSP_BERT"
                             <<"HADRON_QGSP_BERT_HP"<<"HADRON_QGSP_BIC"<<"HADRON_QGSP_BIC_AllHP"<<"HADRON_INCLXX"<<"HADRON_Shielding"<<"HADRON_ShieldingLEND"
                             <<"FACTORY_FTFP_BERT"<<"FACTORY_FTFP_BERT_ATL"<<"FACTORY_FTFP_BERT_TRV"<<"FACTORY_QGSP_FTFP_BERT"<<"FACTORY_QGSP_BERT"
                             <<"FACTORY_QGSP_BERT_HP"<<"FACTORY_QGSP_BIC"<<"FACTORY_QGSP_BIC_AllHP"<<"FACTORY_Shielding"<<"FACTORY_ShieldingLEND");

    ui->SourceComboBoxPhysUsed->addItems(Physicslist);

    AngleUnits=(QStringList()<<"degree"<<"rad");
    AngleSizeUnits=(QStringList()<<"degree"<<"rad"<<"mm"<<"cm"<<"m");
    SizeUnits=(QStringList()<<"mm"<<"cm"<<"m");
    EnergyUnits=(QStringList()<<"eV"<<"keV"<<"MeV");

    ui->comboBoxSizeUnit->addItems(SizeUnits);
    ui->comboBoxAngleUnit->addItems(AngleSizeUnits);
    ui->comboBoxEnergyUnit->addItems(EnergyUnits);
    ui->comboBoxWorldSizeUnit->addItems(SizeUnits);
    ui->comboBoxEnergyUnitsForCrossSection->addItems(EnergyUnits);

    //ui->comboBoxCutsSizeUnit->addItems(SizeUnits);
    //ui->comboBoxCutsEnergyUnit->addItems(EnergyUnits);

    ui->comboBoxSizeUnit->setCurrentIndex(1);
    ui->comboBoxAngleUnit->setCurrentIndex(0);
    ui->comboBoxEnergyUnit->setCurrentIndex(2);

    // PEE
    ValuesOfInputs["G4PEEffectFluoModel"] = "1";
    ValuesOfInputs["G4LivermorePhotoElectricModel"] = "2";
    ValuesOfInputs["G4LivermorePolarizedPhotoElectricModel"] = "3";
    ValuesOfInputs["G4PenelopePhotoElectricModel"] = "4";
    QStringList EPEModelslist=(QStringList()
                               <<"G4PEEffectFluoModel"
                               <<"G4LivermorePhotoElectricModel"
                               <<"G4LivermorePolarizedPhotoElectricModel"
                               <<"G4PenelopePhotoElectricModel");
    ui->comboBoxPEEModels->addItems(EPEModelslist);

    // CE
    ValuesOfInputs["G4KleinNishinaCompton"] = "1";
    ValuesOfInputs["G4KleinNishinaModel"] = "2";
    ValuesOfInputs["G4LowEPComptonModel"] = "3";
    ValuesOfInputs["G4LivermoreComptonModel"] = "4";
    ValuesOfInputs["G4LivermoreComptonModifiedModel"] = "5";
    ValuesOfInputs["G4LivermorePolarizedComptonModel"] = "6";
    ValuesOfInputs["G4PenelopeComptonModel"] = "7";
    ValuesOfInputs["G4TKleinNishinaCompton"] = "8";
    QStringList ComptonModelslist=(QStringList()
                                   <<"G4KleinNishinaCompton"
                                   <<"G4KleinNishinaModel"
                                   <<"G4LowEPComptonModel"
                                   <<"G4LivermoreComptonModel"
                                   <<"G4LivermoreComptonModifiedModel"
                                   <<"G4LivermorePolarizedComptonModel"
                                   <<"G4PenelopeComptonModel"
                                   <<"G4TKleinNishinaCompton");
    ui->comboBoxComptonModels->addItems(ComptonModelslist);

    // GC
    ValuesOfInputs["G4BetheHeitlerModel"] = "1";
    ValuesOfInputs["G4BetheHeitler5DModel"] = "2";
    ValuesOfInputs["G4PairProductionRelModel"] = "3";
    ValuesOfInputs["G4LivermoreGammaConversionModel"] = "4";
    ValuesOfInputs["G4BoldyshevTripletModel"] = "5";
    ValuesOfInputs["G4LivermoreNuclearGammaConversionModel"] = "6";
    ValuesOfInputs["G4LivermorePolarizedGammaConversionModel"] = "7";
    ValuesOfInputs["G4PenelopeGammaConversionModel"] = "8";
    QStringList GammaConversionModelslist=(QStringList()
                                           <<"G4BetheHeitlerModel"
                                           <<"G4BetheHeitler5DModel"
                                           <<"G4PairProductionRelModel"
                                           <<"G4LivermoreGammaConversionModel"
                                           <<"G4BoldyshevTripletModel"
                                           <<"G4LivermoreNuclearGammaConversionModel"
                                           <<"G4LivermorePolarizedGammaConversionModel"
                                           <<"G4PenelopeGammaConversionModel");
    ui->comboBoxGammaConversionModels->addItems(GammaConversionModelslist);

    // RE
    ValuesOfInputs["G4LivermoreRayleighModel"] = "1";
    ValuesOfInputs["G4LivermorePolarizedRayleighModel"] = "2";
    ValuesOfInputs["G4PenelopeRayleighModel"] = "3";
    QStringList RayleighScatteringModelslist=(QStringList()
                                              <<"G4LivermoreRayleighModel"
                                              <<"G4LivermorePolarizedRayleighModel"
                                              <<"G4PenelopeRayleighModel");
    ui->comboBoxRayleighScatteringModels->addItems(RayleighScatteringModelslist);

    // just e- ionization
    ValuesOfInputs["G4MollerBhabhaModel"] = "1";
    ValuesOfInputs["G4LivermoreIonisationModel"] = "2";
    QStringList electrIonilist=(QStringList()
                                <<"G4MollerBhabhaModel"
                                <<"G4LivermoreIonisationModel");
    ui->comboBoxElectronIonisationModels->addItems(electrIonilist);

    // e- and e+ brem
    ValuesOfInputs["G4SeltzerBergerModel"] = "7";
    ValuesOfInputs["G4eBremsstrahlungRelModel"] = "7";
    ValuesOfInputs["G4LivermoreBremsstrahlungModel"] = "7";
    ValuesOfInputs["G4PenelopeBremsstrahlungModel"] = "7";
    QStringList electrBremlist=(QStringList()
                                <<"G4SeltzerBergerModel"
                                <<"G4eBremsstrahlungRelModel"
                                <<"G4LivermoreBremsstrahlungModel"
                                <<"G4PenelopeBremsstrahlungModel");
    ui->comboBoxElectronBremModels->addItems(electrBremlist);


    // ion ionization
    ValuesOfInputs["G4BetheBlochModel"] = "1";
    ValuesOfInputs["G4BetheBlochIonGasModel"] = "2";
    ValuesOfInputs["G4BraggIonModel"] = "3";
    ValuesOfInputs["G4BraggIonGasModel"] = "4";
    ValuesOfInputs["G4IonParametrisedLossModel"] = "5";
    ValuesOfInputs["G4AtimaEnergyLossModel"] = "6";
    ValuesOfInputs["G4LindhardSorensenIonModel"] = "7";
    QStringList IonIonilist = (QStringList()
                               <<"G4BetheBlochModel"
                               <<"G4BetheBlochIonGasModel"
                               <<"G4BraggIonModel"
                               <<"G4BraggIonGasModel"
                               <<"G4IonParametrisedLossModel"
                               <<"G4AtimaEnergyLossModel"
                               <<"G4LindhardSorensenIonModel");

    ValuesOfInputs["Phantom parametrization"] = "0";
    ValuesOfInputs["Phantom nested parametrization"] = "1";

    // hadron ionization
    ValuesOfInputs["G4BetheBlochModel"] = "1";
    ValuesOfInputs["G4BraggModel"] = "2";
    ValuesOfInputs["G4ICRU73QOModel"] = "3";
    QStringList hadronIonilist=(QStringList()
                                <<"G4BetheBlochModel"
                                <<"G4BraggModel"
                                <<"G4ICRU73QOModel");
    ui->comboBoxHadronIonisationModels->addItems(hadronIonilist);

    QStringList StandartSourceslist=(QStringList()<<"Volume"<<"Voxels"<<"TET"<<"Point"<<"Beam"<<"Plane"<<"Surface"<<"Solid");
    ui->comboBoxTypeOfSources->addItems(StandartSourceslist);

    // score

    //QStringList AccuracyLevellist=(QStringList()<<"StepLevel"<<"EventLevel"<<"BatchLevel");
    //ui->ScoreComboBoxAcuuracyLevel->addItems(AccuracyLevellist);

    ValuesOfInputs["One-source Configuration (o)"] = "o";
    ValuesOfInputs["Multi-source Configurations (m)"] = "m";

    ValuesOfInputs["0"] = "Phantom Parameterisation";
    ValuesOfInputs["1"] = "Phantom Nested Parameterisation";

    QStringList SimNumOnRanksList=(QStringList()
                                   <<"One-source Configuration (o)"
                                   <<"Multi-source Configurations (m)");
    ui->ScoreCombobowSimNumOnRanksLineEdit->addItems(SimNumOnRanksList);

    //to MeV, to mm
    UnitsConversionFactors["mm"] = 1;
    UnitsConversionFactors["cm"] = 10;
    UnitsConversionFactors["m"] = 1000;
    UnitsConversionFactors["eV"] = 1e-6;
    UnitsConversionFactors["keV"] = 1e-3;
    UnitsConversionFactors["MeV"] = 1;


    //QStringList RadNuc=(QStringList()<<"F18"<<"I131"<<"I123");
    //ui->lineEditRadioNucleid->addItems(RadNuc);

    QStringList ActUnits=(QStringList()<<"MBq"<<"Bq");
    ui->comboBoxActivityUnits->addItems(ActUnits);

    QStringList TimeUnits=(QStringList()<<"y"<<"h"<<"min"<<"s");
    ui->comboBoxTimeUnit->addItems(TimeUnits);

    QStringList AEUnits=(QStringList()<<"eV"<<"keV"<<"MeV");
    ui->comboBoxAEUnits->addItems(AEUnits);
    QStringList SAFUnits=(QStringList()<<"g-1"<<"kg-1");
    ui->comboBoxSAFUnits->addItems(SAFUnits);
    QStringList ADUnits=(QStringList()<<"MeV/kg"<<"nGy"<<"miGy"<<"mGy"<<"mGy/MBq"<<"Gy"<<"kGy"<<"MGy");
    ui->comboBoxADUnits->addItems(ADUnits);
    ui->comboBoxSUnits->addItems(ADUnits);
    QStringList HUnits=(QStringList()<<"MeV/kg"<<"nGy"<<"miGy"<<"mGy"<<"mGy/MBq"<<"Gy"<<"kGy"<<"MGy"<<"mSv"<<"mSv/MBq"<<"Sv");
    ui->comboBoxHUnits->addItems(HUnits);
    ui->comboBoxEUnits->addItems(HUnits);


    // ROOT analysis

    QStringList Analy1=(QStringList()<<"Mass"<<"Volume"<<"Density"<<"Distance");
    ui->AnalysisComboBoxRegionVariable->addItems(Analy1);

    QStringList Analy2=(QStringList()<<"X"<<"Y"<<"Z");
    ui->AnalysisComboBoxBeamAxis->addItems(Analy2);

    QStringList Analy3=(QStringList()<<"XY"<<"YZ"<<"XZ");
    ui->AnalysisComboBoxSliceFor2DGraph->addItems(Analy3);


    ValuesOfInputs["RightTop"] = 1;
    ValuesOfInputs["RightBottom"] = 2;
    ValuesOfInputs["LeftTop"] = 3;
    ValuesOfInputs["RightBottom"] = 4;
    ValuesOfInputs["MiddleTop"] = 5;
    ValuesOfInputs["MiddleBottom"] = 6;

    QStringList LegPosList=(QStringList()
                            <<"RightTop"
                            <<"RightBottom"
                            <<"LeftTop"
                            <<"RightBottom"
                            <<"MiddleTop"
                            <<"MiddleBottom");

    ui->AnalysisLegendPoscomboBox->addItems(LegPosList);

    LegPosList=(QStringList()
                <<"MPI"
                <<"MT");

    ui->MPIOrMTOnRockscomboBox->addItems(LegPosList);

    LegPosList=(QStringList()
                <<""
                <<"RA"
                <<"LRD"
                <<"RD");

    ui->comboBoxRelDiff->addItems(LegPosList);


    QStringList DataFilesUses=(QStringList()
                               <<""
                               <<"read"
                               <<"save"
                               <<"generate"
                               <<"MyCPPGenerator");

    ui->UseDataFilesFor->addItems(DataFilesUses);

    QStringList RunForList=(QStringList()
                               <<"InternalDosimetry"
                               <<"ExternalDosimetry"
                               <<"MyCPPSteppingAction"
                               //<<"Neutronic"
                               //<<"Detector"
                               );

    ui->comboBoxSimulationRunFor->addItems(RunForList);

    QStringList PredefinedSources=(QStringList()
                                   <<""
                                   <<"External RLAT(Ang=-90) for Stylized Phantoms"
                                   <<"External LLAT(Ang=90) for Stylized Phantoms"
                                   <<"External PA(Ang=+-180) for Stylized Phantoms"
                                   <<"External AP(Ang=0) for Stylized Phantoms"
                                   <<"External ROT for Stylized Phantoms"
                                   <<"External ISO for Stylized Phantoms"
                                   <<"External RLAT(Ang=-90) for Voxelized and Mesh-type Phantoms"
                                   <<"External LLAT(Ang=90) for Voxelized and Mesh-type Phantoms"
                                   <<"External PA(Ang=+-180) for Voxelized and Mesh-type Phantoms"
                                   <<"External AP(Ang=0) for Voxelized and Mesh-type Phantoms"
                                   <<"External ROT for Voxelized and Mesh-type Phantoms"
                                   <<"External ISO for Voxelized and Mesh-type Phantoms"
                                   <<"External Ang=-75 for Stylized Phantoms"
                                   <<"External Ang=-60 for Stylized Phantoms"
                                   <<"External Ang=-45 for Stylized Phantoms"
                                   <<"External Ang=-30 for Stylized Phantoms"
                                   <<"External Ang=-15 for Stylized Phantoms"
                                   <<"External Ang=+15 for Stylized Phantoms"
                                   <<"External Ang=+30 for Stylized Phantoms"
                                   <<"External Ang=+45 for Stylized Phantoms"
                                   <<"External Ang=+60 for Stylized Phantoms"
                                   <<"External Ang=+75 for Stylized Phantoms"
                               );

    ui->comboBoxPredefinedSources->addItems(PredefinedSources);

    LegPosList=(QStringList()
                <<"DCC"
                <<"Ambient Dose"
                <<"Personal Dose"
                <<"Directional or Personal Absorbed Dose in the Lens of the Eye"
                <<"Directional or Personal Absorbed Dose in Local Skin"
                );

    ui->comboBoxExternalDosQuantities->addItems(LegPosList);

    LegPosList=(QStringList()
                <<"e-"
                <<"gamma"
                <<"e+"
                <<"alpha"
                <<"proton"
                <<"neutron"
                );

    ui->comboBoxExternalDoseParticles->addItems(LegPosList);

    LegPosList=(QStringList()
                <<"For Each Source"
                <<"For Each Target"
                <<"For Each Energy"
                );

    ui->comboBoxGenerateExternalFor->addItems(LegPosList);

    ElementsSymbolSym["1-H-1.008"]="H";     ElementsSymbolZ["1-H-1.008"]="1";ElementsSymbolA["1-H-1.008"]="1.008";ElementsSymbolName["1-H-1.008"]="Hydrogen";
    ElementsSymbolSym["2-He-4.003"]="He";   ElementsSymbolZ["2-He-4.003"]="2";ElementsSymbolA["2-He-4.003"]="4.003";ElementsSymbolName["2-He-4.003"]="Helium";
    ElementsSymbolSym["3-Li-6.941"]="Li";   ElementsSymbolZ["3-Li-6.941"]="3";ElementsSymbolA["3-Li-6.941"]="6.941";ElementsSymbolName["3-Li-6.941"]="Lithium";
    ElementsSymbolSym["4-Be-9.012"]="Be";   ElementsSymbolZ["4-Be-9.012"]="4";ElementsSymbolA["4-Be-9.012"]="9.012";ElementsSymbolName["4-Be-9.012"]="Beryllium";
    ElementsSymbolSym["5-B-10.811"]="B";    ElementsSymbolZ["5-B-10.811"]="5";ElementsSymbolA["5-B-10.811"]="10.811";ElementsSymbolName["5-B-10.811"]="Boron";
    ElementsSymbolSym["6-C-12.011"]="C";    ElementsSymbolZ["6-C-12.011"]="6";ElementsSymbolA["6-C-12.011"]="12.011";ElementsSymbolName["6-C-12.011"]="Carbon";
    ElementsSymbolSym["7-N-14.007"]="N";    ElementsSymbolZ["7-N-14.007"]="7";ElementsSymbolA["7-N-14.007"]="14.007";ElementsSymbolName["7-N-14.007"]="Nitrogen";
    ElementsSymbolSym["8-O-15.999"]="O";    ElementsSymbolZ["8-O-15.999"]="8";ElementsSymbolA["8-O-15.999"]="15.999";ElementsSymbolName["8-O-15.999"]="Oxygen";
    ElementsSymbolSym["9-F-18.998"]="F";    ElementsSymbolZ["9-F-18.998"]="9";ElementsSymbolA["9-F-18.998"]="18.998";ElementsSymbolName["9-F-18.998"]="Fluorine";
    ElementsSymbolSym["10-Ne-20.18"]="Ne";  ElementsSymbolZ["10-Ne-20.18"]="1";ElementsSymbolA["10-Ne-20.18"]="20.18";ElementsSymbolName["10-Ne-20.18"]="Neon";
    ElementsSymbolSym["11-Na-22.99"]="Na";  ElementsSymbolZ["11-Na-22.99"]="11";ElementsSymbolA["11-Na-22.99"]="22.99";ElementsSymbolName["11-Na-22.99"]="Sodium";
    ElementsSymbolSym["12-Mg-24.305"]="Mg"; ElementsSymbolZ["12-Mg-24.305"]="12";ElementsSymbolA["12-Mg-24.305"]="24.305";ElementsSymbolName["12-Mg-24.305"]="Magnesium";
    ElementsSymbolSym["13-Al-26.982"]="Al"; ElementsSymbolZ["13-Al-26.982"]="13";ElementsSymbolA["13-Al-26.982"]="26.982";ElementsSymbolName["13-Al-26.982"]="Aluminum";
    ElementsSymbolSym["14-Si-28.086"]="Si"; ElementsSymbolZ["14-Si-28.086"]="14";ElementsSymbolA["14-Si-28.086"]="28.086";ElementsSymbolName["14-Si-28.086"]="Silicon";
    ElementsSymbolSym["15-P-30.974"]="P";   ElementsSymbolZ["15-P-30.974"]="15";ElementsSymbolA["15-P-30.974"]="30.974";ElementsSymbolName["15-P-30.974"]="Phosphorus";
    ElementsSymbolSym["16-S-32.065"]="S";   ElementsSymbolZ["16-S-32.065"]="16";ElementsSymbolA["16-S-32.065"]="32.065";ElementsSymbolName["16-S-32.065"]="Sulfur";
    ElementsSymbolSym["17-Cl-35.453"]="Cl"; ElementsSymbolZ["17-Cl-35.453"]="17";ElementsSymbolA["17-Cl-35.453"]="35.453";ElementsSymbolName["17-Cl-35.453"]="Chlorine";
    ElementsSymbolSym["18-Ar-39.948"]="Ar"; ElementsSymbolZ["18-Ar-39.948"]="18";ElementsSymbolA["18-Ar-39.948"]="39.948";ElementsSymbolName["18-Ar-39.948"]="Argon";
    ElementsSymbolSym["19-K-39.098"]="K";   ElementsSymbolZ["19-K-39.098"]="19";ElementsSymbolA["19-K-39.098"]="39.098";ElementsSymbolName["19-K-39.098"]="Potassium";
    ElementsSymbolSym["20-Ca-40.078"]="Ca"; ElementsSymbolZ["20-Ca-40.078"]="2";ElementsSymbolA["20-Ca-40.078"]="40.078";ElementsSymbolName["20-Ca-40.078"]="Calcium";
    ElementsSymbolSym["21-Sc-44.956"]="Sc"; ElementsSymbolZ["21-Sc-44.956"]="21";ElementsSymbolA["21-Sc-44.956"]="44.956";ElementsSymbolName["21-Sc-44.956"]="Scandium";
    ElementsSymbolSym["22-Ti-47.867"]="Ti"; ElementsSymbolZ["22-Ti-47.867"]="22";ElementsSymbolA["22-Ti-47.867"]="47.867";ElementsSymbolName["22-Ti-47.867"]="Titanium";
    ElementsSymbolSym["23-V-50.942"]="V";   ElementsSymbolZ["23-V-50.942"]="23";ElementsSymbolA["23-V-50.942"]="50.942";ElementsSymbolName["23-V-50.942"]="Vanadium";
    ElementsSymbolSym["24-Cr-51.996"]="Cr"; ElementsSymbolZ["24-Cr-51.996"]="24";ElementsSymbolA["24-Cr-51.996"]="51.996";ElementsSymbolName["24-Cr-51.996"]="Chromium";
    ElementsSymbolSym["25-Mn-54.938"]="Mn"; ElementsSymbolZ["25-Mn-54.938"]="25";ElementsSymbolA["25-Mn-54.938"]="54.938";ElementsSymbolName["25-Mn-54.938"]="Manganese";
    ElementsSymbolSym["26-Fe-55.845"]="Fe"; ElementsSymbolZ["26-Fe-55.845"]="26";ElementsSymbolA["26-Fe-55.845"]="55.845";ElementsSymbolName["26-Fe-55.845"]="Iron";
    ElementsSymbolSym["27-Co-58.933"]="Co"; ElementsSymbolZ["27-Co-58.933"]="27";ElementsSymbolA["27-Co-58.933"]="58.933";ElementsSymbolName["27-Co-58.933"]="Cobalt";
    ElementsSymbolSym["28-Ni-58.693"]="Ni"; ElementsSymbolZ["28-Ni-58.693"]="28";ElementsSymbolA["28-Ni-58.693"]="58.693";ElementsSymbolName["28-Ni-58.693"]="Nickel";
    ElementsSymbolSym["29-Cu-63.546"]="Cu"; ElementsSymbolZ["29-Cu-63.546"]="29";ElementsSymbolA["29-Cu-63.546"]="63.546";ElementsSymbolName["29-Cu-63.546"]="Copper";
    ElementsSymbolSym["30-Zn-65.39"]="Zn";  ElementsSymbolZ["30-Zn-65.39"]="3";ElementsSymbolA["30-Zn-65.39"]="65.39";ElementsSymbolName["30-Zn-65.39"]="Zinc";
    ElementsSymbolSym["31-Ga-69.723"]="Ga"; ElementsSymbolZ["31-Ga-69.723"]="31";ElementsSymbolA["31-Ga-69.723"]="69.723";ElementsSymbolName["31-Ga-69.723"]="Gallium";
    ElementsSymbolSym["32-Ge-72.64"]="Ge";  ElementsSymbolZ["32-Ge-72.64"]="32";ElementsSymbolA["32-Ge-72.64"]="72.64";ElementsSymbolName["32-Ge-72.64"]="Germanium";
    ElementsSymbolSym["33-As-74.922"]="As"; ElementsSymbolZ["33-As-74.922"]="33";ElementsSymbolA["33-As-74.922"]="74.922";ElementsSymbolName["33-As-74.922"]="Arsenic";
    ElementsSymbolSym["34-Se-78.96"]="Se";  ElementsSymbolZ["34-Se-78.96"]="34";ElementsSymbolA["34-Se-78.96"]="78.96";ElementsSymbolName["34-Se-78.96"]="Selenium";
    ElementsSymbolSym["35-Br-79.904"]="Br"; ElementsSymbolZ["35-Br-79.904"]="35";ElementsSymbolA["35-Br-79.904"]="79.904";ElementsSymbolName["35-Br-79.904"]="Bromine";
    ElementsSymbolSym["36-Kr-83.8"]="Kr";   ElementsSymbolZ["36-Kr-83.8"]="36";ElementsSymbolA["36-Kr-83.8"]="83.8";ElementsSymbolName["36-Kr-83.8"]="Krypton";
    ElementsSymbolSym["37-Rb-85.468"]="Rb"; ElementsSymbolZ["37-Rb-85.468"]="37";ElementsSymbolA["37-Rb-85.468"]="85.468";ElementsSymbolName["37-Rb-85.468"]="Rubidium";
    ElementsSymbolSym["38-Sr-87.62"]="Sr";  ElementsSymbolZ["38-Sr-87.62"]="38";ElementsSymbolA["38-Sr-87.62"]="87.62";ElementsSymbolName["38-Sr-87.62"]="Strontium";
    ElementsSymbolSym["39-Y-88.906"]="Y";   ElementsSymbolZ["39-Y-88.906"]="39";ElementsSymbolA["39-Y-88.906"]="88.906";ElementsSymbolName["39-Y-88.906"]="Yttrium";
    ElementsSymbolSym["40-Zr-91.224"]="Zr"; ElementsSymbolZ["40-Zr-91.224"]="4";ElementsSymbolA["40-Zr-91.224"]="91.224";ElementsSymbolName["40-Zr-91.224"]="Zirconium";
    ElementsSymbolSym["41-Nb-92.906"]="Nb"; ElementsSymbolZ["41-Nb-92.906"]="41";ElementsSymbolA["41-Nb-92.906"]="92.906";ElementsSymbolName["41-Nb-92.906"]="Niobium";
    ElementsSymbolSym["42-Mo-95.94"]="Mo";  ElementsSymbolZ["42-Mo-95.94"]="42";ElementsSymbolA["42-Mo-95.94"]="95.94";ElementsSymbolName["42-Mo-95.94"]="Molybdenum";
    ElementsSymbolSym["43-Tc-98"]="Tc";     ElementsSymbolZ["43-Tc-98"]="43";ElementsSymbolA["43-Tc-98"]="98";ElementsSymbolName["43-Tc-98"]="Technetium";
    ElementsSymbolSym["44-Ru-101.07"]="Ru"; ElementsSymbolZ["44-Ru-101.07"]="44";ElementsSymbolA["44-Ru-101.07"]="101.07";ElementsSymbolName["44-Ru-101.07"]="Ruthenium";
    ElementsSymbolSym["45-Rh-102.906"]="Rh";ElementsSymbolZ["45-Rh-102.906"]="45";ElementsSymbolA["45-Rh-102.906"]="102.906";ElementsSymbolName["45-Rh-102.906"]="Rhodium";
    ElementsSymbolSym["46-Pd-106.42"]="Pd"; ElementsSymbolZ["46-Pd-106.42"]="46";ElementsSymbolA["46-Pd-106.42"]="106.42";ElementsSymbolName["46-Pd-106.42"]="Palladium";
    ElementsSymbolSym["47-Ag-107.868"]="Ag";ElementsSymbolZ["47-Ag-107.868"]="47";ElementsSymbolA["47-Ag-107.868"]="107.868";ElementsSymbolName["47-Ag-107.868"]="Silver";
    ElementsSymbolSym["48-Cd-112.411"]="Cd";ElementsSymbolZ["48-Cd-112.411"]="48";ElementsSymbolA["48-Cd-112.411"]="112.411";ElementsSymbolName["48-Cd-112.411"]="Cadmium";
    ElementsSymbolSym["49-In-114.818"]="In";ElementsSymbolZ["49-In-114.818"]="49";ElementsSymbolA["49-In-114.818"]="114.818";ElementsSymbolName["49-In-114.818"]="Indium";
    ElementsSymbolSym["50-Sn-118.71"]="Sn"; ElementsSymbolZ["50-Sn-118.71"]="5";ElementsSymbolA["50-Sn-118.71"]="118.71";ElementsSymbolName["50-Sn-118.71"]="Tin";
    ElementsSymbolSym["51-Sb-121.76"]="Sb"; ElementsSymbolZ["51-Sb-121.76"]="51";ElementsSymbolA["51-Sb-121.76"]="121.76";ElementsSymbolName["51-Sb-121.76"]="Antimony";
    ElementsSymbolSym["52-Te-127.6"]="Te";  ElementsSymbolZ["52-Te-127.6"]="52";ElementsSymbolA["52-Te-127.6"]="127.6";ElementsSymbolName["52-Te-127.6"]="Tellurium";
    ElementsSymbolSym["53-I-126.905"]="I";  ElementsSymbolZ["53-I-126.905"]="53";ElementsSymbolA["53-I-126.905"]="126.905";ElementsSymbolName["53-I-126.905"]="Iodine";
    ElementsSymbolSym["54-Xe-131.293"]="Xe";ElementsSymbolZ["54-Xe-131.293"]="54";ElementsSymbolA["54-Xe-131.293"]="131.293";ElementsSymbolName["54-Xe-131.293"]="Xenon";
    ElementsSymbolSym["55-Cs-132.906"]="Cs";ElementsSymbolZ["55-Cs-132.906"]="55";ElementsSymbolA["55-Cs-132.906"]="132.906";ElementsSymbolName["55-Cs-132.906"]="Cesium";
    ElementsSymbolSym["56-Ba-137.327"]="Ba";ElementsSymbolZ["56-Ba-137.327"]="56";ElementsSymbolA["56-Ba-137.327"]="137.327";ElementsSymbolName["56-Ba-137.327"]="Barium";
    ElementsSymbolSym["57-La-138.906"]="La";ElementsSymbolZ["57-La-138.906"]="57";ElementsSymbolA["57-La-138.906"]="138.906";ElementsSymbolName["57-La-138.906"]="Lanthanum";
    ElementsSymbolSym["58-Ce-140.116"]="Ce";ElementsSymbolZ["58-Ce-140.116"]="58";ElementsSymbolA["58-Ce-140.116"]="140.116";ElementsSymbolName["58-Ce-140.116"]="Cerium";
    ElementsSymbolSym["59-Pr-140.908"]="Pr";ElementsSymbolZ["59-Pr-140.908"]="59";ElementsSymbolA["59-Pr-140.908"]="140.908";ElementsSymbolName["59-Pr-140.908"]="Praseodymium";
    ElementsSymbolSym["60-Nd-144.24"]="Nd"; ElementsSymbolZ["60-Nd-144.24"]="6";ElementsSymbolA["60-Nd-144.24"]="144.24";ElementsSymbolName["60-Nd-144.24"]="Neodymium";
    ElementsSymbolSym["61-Pm-145"]="Pm";    ElementsSymbolZ["61-Pm-145"]="61";ElementsSymbolA["61-Pm-145"]="145";ElementsSymbolName["61-Pm-145"]="Promethium";
    ElementsSymbolSym["62-Sm-150.36"]="Sm"; ElementsSymbolZ["62-Sm-150.36"]="62";ElementsSymbolA["62-Sm-150.36"]="150.36";ElementsSymbolName["62-Sm-150.36"]="Samarium";
    ElementsSymbolSym["63-Eu-151.964"]="Eu";ElementsSymbolZ["63-Eu-151.964"]="63";ElementsSymbolA["63-Eu-151.964"]="151.964";ElementsSymbolName["63-Eu-151.964"]="Europium";
    ElementsSymbolSym["64-Gd-157.25"]="Gd"; ElementsSymbolZ["64-Gd-157.25"]="64";ElementsSymbolA["64-Gd-157.25"]="157.25";ElementsSymbolName["64-Gd-157.25"]="Gadolinium";
    ElementsSymbolSym["65-Tb-158.925"]="Tb";ElementsSymbolZ["65-Tb-158.925"]="65";ElementsSymbolA["65-Tb-158.925"]="158.925";ElementsSymbolName["65-Tb-158.925"]="Terbium";
    ElementsSymbolSym["66-Dy-162.5"]="Dy";  ElementsSymbolZ["66-Dy-162.5"]="66";ElementsSymbolA["66-Dy-162.5"]="162.5";ElementsSymbolName["66-Dy-162.5"]="Dysprosium";
    ElementsSymbolSym["67-Ho-164.93"]="Ho"; ElementsSymbolZ["67-Ho-164.93"]="67";ElementsSymbolA["67-Ho-164.93"]="164.93";ElementsSymbolName["67-Ho-164.93"]="Holmium";
    ElementsSymbolSym["68-Er-167.259"]="Er";ElementsSymbolZ["68-Er-167.259"]="68";ElementsSymbolA["68-Er-167.259"]="167.259";ElementsSymbolName["68-Er-167.259"]="Erbium";
    ElementsSymbolSym["69-Tm-168.934"]="Tm";ElementsSymbolZ["69-Tm-168.934"]="69";ElementsSymbolA["69-Tm-168.934"]="168.934";ElementsSymbolName["69-Tm-168.934"]="Thulium";
    ElementsSymbolSym["70-Yb-173.04"]="Yb"; ElementsSymbolZ["70-Yb-173.04"]="7";ElementsSymbolA["70-Yb-173.04"]="173.04";ElementsSymbolName["70-Yb-173.04"]="Ytterbium";
    ElementsSymbolSym["71-Lu-174.967"]="Lu";ElementsSymbolZ["71-Lu-174.967"]="71";ElementsSymbolA["71-Lu-174.967"]="174.967";ElementsSymbolName["71-Lu-174.967"]="Lutetium";
    ElementsSymbolSym["72-Hf-178.49"]="Hf"; ElementsSymbolZ["72-Hf-178.49"]="72";ElementsSymbolA["72-Hf-178.49"]="178.49";ElementsSymbolName["72-Hf-178.49"]="Hafnium";
    ElementsSymbolSym["73-Ta-180.948"]="Ta";ElementsSymbolZ["73-Ta-180.948"]="73";ElementsSymbolA["73-Ta-180.948"]="180.948";ElementsSymbolName["73-Ta-180.948"]="Tantalum";
    ElementsSymbolSym["74-W-183.84"]="W";   ElementsSymbolZ["74-W-183.84"]="74";ElementsSymbolA["74-W-183.84"]="183.84";ElementsSymbolName["74-W-183.84"]="Tungsten";
    ElementsSymbolSym["75-Re-186.207"]="Re";ElementsSymbolZ["75-Re-186.207"]="75";ElementsSymbolA["75-Re-186.207"]="186.207";ElementsSymbolName["75-Re-186.207"]="Rhenium";
    ElementsSymbolSym["76-Os-190.23"]="Os"; ElementsSymbolZ["76-Os-190.23"]="76";ElementsSymbolA["76-Os-190.23"]="190.23";ElementsSymbolName["76-Os-190.23"]="Osmium";
    ElementsSymbolSym["77-Ir-192.217"]="Ir";ElementsSymbolZ["77-Ir-192.217"]="77";ElementsSymbolA["77-Ir-192.217"]="192.217";ElementsSymbolName["77-Ir-192.217"]="Iridium";
    ElementsSymbolSym["78-Pt-195.078"]="Pt";ElementsSymbolZ["78-Pt-195.078"]="78";ElementsSymbolA["78-Pt-195.078"]="195.078";ElementsSymbolName["78-Pt-195.078"]="Platinum";
    ElementsSymbolSym["79-Au-196.967"]="Au";ElementsSymbolZ["79-Au-196.967"]="79";ElementsSymbolA["79-Au-196.967"]="196.967";ElementsSymbolName["79-Au-196.967"]="Gold";
    ElementsSymbolSym["80-Hg-200.59"]="Hg"; ElementsSymbolZ["80-Hg-200.59"]="8";ElementsSymbolA["80-Hg-200.59"]="200.59";ElementsSymbolName["80-Hg-200.59"]="Mercury";
    ElementsSymbolSym["81-Tl-204.383"]="Tl";ElementsSymbolZ["81-Tl-204.383"]="81";ElementsSymbolA["81-Tl-204.383"]="204.383";ElementsSymbolName["81-Tl-204.383"]="Thallium";
    ElementsSymbolSym["82-Pb-207.2"]="Pb";  ElementsSymbolZ["82-Pb-207.2"]="82";ElementsSymbolA["82-Pb-207.2"]="207.2";ElementsSymbolName["82-Pb-207.2"]="Lead";
    ElementsSymbolSym["83-Bi-208.98"]="Bi"; ElementsSymbolZ["83-Bi-208.98"]="83";ElementsSymbolA["83-Bi-208.98"]="208.98";ElementsSymbolName["83-Bi-208.98"]="Bismuth";
    ElementsSymbolSym["84-Po-209"]="Po";    ElementsSymbolZ["84-Po-209"]="84";ElementsSymbolA["84-Po-209"]="209";ElementsSymbolName["84-Po-209"]="Polonium";
    ElementsSymbolSym["85-At-21"]="At";     ElementsSymbolZ["85-At-21"]="85";ElementsSymbolA["85-At-21"]="21";ElementsSymbolName["85-At-21"]="Astatine";
    ElementsSymbolSym["86-Rn-222"]="Rn";    ElementsSymbolZ["86-Rn-222"]="86";ElementsSymbolA["86-Rn-222"]="222";ElementsSymbolName["86-Rn-222"]="Radon";
    ElementsSymbolSym["87-Fr-223"]="Fr";    ElementsSymbolZ["87-Fr-223"]="87";ElementsSymbolA["87-Fr-223"]="223";ElementsSymbolName["87-Fr-223"]="Francium";
    ElementsSymbolSym["88-Ra-226"]="Ra";    ElementsSymbolZ["88-Ra-226"]="88";ElementsSymbolA["88-Ra-226"]="226";ElementsSymbolName["88-Ra-226"]="Radium";
    ElementsSymbolSym["89-Ac-227"]="Ac";    ElementsSymbolZ["89-Ac-227"]="89";ElementsSymbolA["89-Ac-227"]="227";ElementsSymbolName["89-Ac-227"]="Actinium";
    ElementsSymbolSym["90-Th-232.038"]="Th";ElementsSymbolZ["90-Th-232.038"]="9";ElementsSymbolA["90-Th-232.038"]="232.038";ElementsSymbolName["90-Th-232.038"]="Thorium";
    ElementsSymbolSym["91-Pa-231.036"]="Pa";ElementsSymbolZ["91-Pa-231.036"]="91";ElementsSymbolA["91-Pa-231.036"]="231.036";ElementsSymbolName["91-Pa-231.036"]="Protactinium";
    ElementsSymbolSym["92-U-238.029"]="U";  ElementsSymbolZ["92-U-238.029"]="92";ElementsSymbolA["92-U-238.029"]="238.029";ElementsSymbolName["92-U-238.029"]="Uranium";
    ElementsSymbolSym["93-Np-237"]="Np";    ElementsSymbolZ["93-Np-237"]="93";ElementsSymbolA["93-Np-237"]="237";ElementsSymbolName["93-Np-237"]="Neptunium";
    ElementsSymbolSym["94-Pu-244"]="Pu";    ElementsSymbolZ["94-Pu-244"]="94";ElementsSymbolA["94-Pu-244"]="244";ElementsSymbolName["94-Pu-244"]="Plutonium";
    ElementsSymbolSym["95-Am-243"]="Am";    ElementsSymbolZ["95-Am-243"]="95";ElementsSymbolA["95-Am-243"]="243";ElementsSymbolName["95-Am-243"]="Americium";
    ElementsSymbolSym["96-Cm-247"]="Cm";    ElementsSymbolZ["96-Cm-247"]="96";ElementsSymbolA["96-Cm-247"]="247";ElementsSymbolName["96-Cm-247"]="Curium";
    ElementsSymbolSym["97-Bk-247"]="Bk";    ElementsSymbolZ["97-Bk-247"]="97";ElementsSymbolA["97-Bk-247"]="247";ElementsSymbolName["97-Bk-247"]="Berkelium";
    ElementsSymbolSym["98-Cf-251"]="Cf";    ElementsSymbolZ["98-Cf-251"]="98";ElementsSymbolA["98-Cf-251"]="251";ElementsSymbolName["98-Cf-251"]="Californium";
    ElementsSymbolSym["99-Es-252"]="Es";    ElementsSymbolZ["99-Es-252"]="99";ElementsSymbolA["99-Es-252"]="252";ElementsSymbolName["99-Es-252"]="Einsteinium";
    ElementsSymbolSym["100-Fm-257"]="Fm";   ElementsSymbolZ["100-Fm-257"]="1";ElementsSymbolA["100-Fm-257"]="257";ElementsSymbolName["100-Fm-257"]="Fermium";
    ElementsSymbolSym["101-Md-258"]="Md";   ElementsSymbolZ["101-Md-258"]="101";ElementsSymbolA["101-Md-258"]="258";ElementsSymbolName["101-Md-258"]="Mendelevium";
    ElementsSymbolSym["102-No-259"]="No";   ElementsSymbolZ["102-No-259"]="102";ElementsSymbolA["102-No-259"]="259";ElementsSymbolName["102-No-259"]="Nobelium";
    ElementsSymbolSym["103-Lr-262"]="Lr";   ElementsSymbolZ["103-Lr-262"]="103";ElementsSymbolA["103-Lr-262"]="262";ElementsSymbolName["103-Lr-262"]="Lawrencium";
    ElementsSymbolSym["104-Rf-261"]="Rf";   ElementsSymbolZ["104-Rf-261"]="104";ElementsSymbolA["104-Rf-261"]="261";ElementsSymbolName["104-Rf-261"]="Rutherfordium";
    ElementsSymbolSym["105-Db-262"]="Db";   ElementsSymbolZ["105-Db-262"]="105";ElementsSymbolA["105-Db-262"]="262";ElementsSymbolName["105-Db-262"]="Dubnium";
    ElementsSymbolSym["106-Sg-266"]="Sg";   ElementsSymbolZ["106-Sg-266"]="106";ElementsSymbolA["106-Sg-266"]="266";ElementsSymbolName["106-Sg-266"]="Seaborgium";
    ElementsSymbolSym["107-Bh-264"]="Bh";   ElementsSymbolZ["107-Bh-264"]="107";ElementsSymbolA["107-Bh-264"]="264";ElementsSymbolName["107-Bh-264"]="Bohrium";
    ElementsSymbolSym["108-Hs-277"]="Hs";   ElementsSymbolZ["108-Hs-277"]="108";ElementsSymbolA["108-Hs-277"]="277";ElementsSymbolName["108-Hs-277"]="Hassium";
    ElementsSymbolSym["109-Mt-268"]="Mt";   ElementsSymbolZ["109-Mt-268"]="109";ElementsSymbolA["109-Mt-268"]="268";ElementsSymbolName["109-Mt-268"]="Meitnerium";

    PeriodicTableElementsSymbol=(QStringList()
                                 <<"1-H-1.008"
                                 <<"2-He-4.003"
                                 <<"3-Li-6.941"
                                 <<"4-Be-9.012"
                                 <<"5-B-10.811"
                                 <<"6-C-12.011"
                                 <<"7-N-14.007"
                                 <<"8-O-15.999"
                                 <<"9-F-18.998"
                                 <<"10-Ne-20.18"
                                 <<"11-Na-22.99"
                                 <<"12-Mg-24.305"
                                 <<"13-Al-26.982"
                                 <<"14-Si-28.086"
                                 <<"15-P-30.974"
                                 <<"16-S-32.065"
                                 <<"17-Cl-35.453"
                                 <<"18-Ar-39.948"
                                 <<"19-K-39.098"
                                 <<"20-Ca-40.078"
                                 <<"21-Sc-44.956"
                                 <<"22-Ti-47.867"
                                 <<"23-V-50.942"
                                 <<"24-Cr-51.996"
                                 <<"25-Mn-54.938"
                                 <<"26-Fe-55.845"
                                 <<"27-Co-58.933"
                                 <<"28-Ni-58.693"
                                 <<"29-Cu-63.546"
                                 <<"30-Zn-65.39"
                                 <<"31-Ga-69.723"
                                 <<"32-Ge-72.64"
                                 <<"33-As-74.922"
                                 <<"34-Se-78.96"
                                 <<"35-Br-79.904"
                                 <<"36-Kr-83.8"
                                 <<"37-Rb-85.468"
                                 <<"38-Sr-87.62"
                                 <<"39-Y-88.906"
                                 <<"40-Zr-91.224"
                                 <<"41-Nb-92.906"
                                 <<"42-Mo-95.94"
                                 <<"43-Tc-98"
                                 <<"44-Ru-101.07"
                                 <<"45-Rh-102.906"
                                 <<"46-Pd-106.42"
                                 <<"47-Ag-107.868"
                                 <<"48-Cd-112.411"
                                 <<"49-In-114.818"
                                 <<"50-Sn-118.71"
                                 <<"51-Sb-121.76"
                                 <<"52-Te-127.6"
                                 <<"53-I-126.905"
                                 <<"54-Xe-131.293"
                                 <<"55-Cs-132.906"
                                 <<"56-Ba-137.327"
                                 <<"57-La-138.906"
                                 <<"58-Ce-140.116"
                                 <<"59-Pr-140.908"
                                 <<"60-Nd-144.24"
                                 <<"61-Pm-145"
                                 <<"62-Sm-150.36"
                                 <<"63-Eu-151.964"
                                 <<"64-Gd-157.25"
                                 <<"65-Tb-158.925"
                                 <<"66-Dy-162.5"
                                 <<"67-Ho-164.93"
                                 <<"68-Er-167.259"
                                 <<"69-Tm-168.934"
                                 <<"70-Yb-173.04"
                                 <<"71-Lu-174.967"
                                 <<"72-Hf-178.49"
                                 <<"73-Ta-180.948"
                                 <<"74-W-183.84"
                                 <<"75-Re-186.207"
                                 <<"76-Os-190.23"
                                 <<"77-Ir-192.217"
                                 <<"78-Pt-195.078"
                                 <<"79-Au-196.967"
                                 <<"80-Hg-200.59"
                                 <<"81-Tl-204.383"
                                 <<"82-Pb-207.2"
                                 <<"83-Bi-208.98"
                                 <<"84-Po-209"
                                 <<"85-At-21"
                                 <<"86-Rn-222"
                                 <<"87-Fr-223"
                                 <<"88-Ra-226"
                                 <<"89-Ac-227"
                                 <<"90-Th-232.038"
                                 <<"91-Pa-231.036"
                                 <<"92-U-238.029"
                                 <<"93-Np-237"
                                 <<"94-Pu-244"
                                 <<"95-Am-243"
                                 <<"96-Cm-247"
                                 <<"97-Bk-247"
                                 <<"98-Cf-251"
                                 <<"99-Es-252"
                                 <<"100-Fm-257"
                                 <<"101-Md-258"
                                 <<"102-No-259"
                                 <<"103-Lr-262"
                                 <<"104-Rf-261"
                                 <<"105-Db-262"
                                 <<"106-Sg-266"
                                 <<"107-Bh-264"
                                 <<"108-Hs-277"
                                 <<"109-Mt-268");

    NistMaterialsDataBase=(QStringList()
                           <<"G4_H"
                           <<"G4_He"
                           <<"G4_Li"
                           <<"G4_Be"
                           <<"G4_B"
                           <<"G4_C"
                           <<"G4_N"
                           <<"G4_O"
                           <<"G4_F"
                           <<"G4_Ne"
                           <<"G4_Na"
                           <<"G4_Mg"
                           <<"G4_Al"
                           <<"G4_Si"
                           <<"G4_P"
                           <<"G4_S"
                           <<"G4_Cl"
                           <<"G4_Ar"
                           <<"G4_K"
                           <<"G4_Ca"
                           <<"G4_Sc"
                           <<"G4_Ti"
                           <<"G4_V"
                           <<"G4_Cr"
                           <<"G4_Mn"
                           <<"G4_Fe"
                           <<"G4_Co"
                           <<"G4_Ni"
                           <<"G4_Cu"
                           <<"G4_Zn"
                           <<"G4_Ga"
                           <<"G4_Ge"
                           <<"G4_As"
                           <<"G4_Se"
                           <<"G4_Br"
                           <<"G4_Kr"
                           <<"G4_Rb"
                           <<"G4_Sr"
                           <<"G4_Y"
                           <<"G4_Zr"
                           <<"G4_Nb"
                           <<"G4_Mo"
                           <<"G4_Tc"
                           <<"G4_Ru"
                           <<"G4_Rh"
                           <<"G4_Pd"
                           <<"G4_Ag"
                           <<"G4_Cd"
                           <<"G4_In"
                           <<"G4_Sn"
                           <<"G4_Sb"
                           <<"G4_Te"
                           <<"G4_I"
                           <<"G4_Xe"
                           <<"G4_Cs"
                           <<"G4_Ba"
                           <<"G4_La"
                           <<"G4_Ce"
                           <<"G4_Pr"
                           <<"G4_Nd"
                           <<"G4_Pm"
                           <<"G4_Sm"
                           <<"G4_Eu"
                           <<"G4_Gd"
                           <<"G4_Tb"
                           <<"G4_Dy"
                           <<"G4_Ho"
                           <<"G4_Er"
                           <<"G4_Tm"
                           <<"G4_Yb"
                           <<"G4_Lu"
                           <<"G4_Hf"
                           <<"G4_Ta"
                           <<"G4_W"
                           <<"G4_Re"
                           <<"G4_Os"
                           <<"G4_Ir"
                           <<"G4_Pt"
                           <<"G4_Au"
                           <<"G4_Hg"
                           <<"G4_Tl"
                           <<"G4_Pb"
                           <<"G4_Bi"
                           <<"G4_Po"
                           <<"G4_At"
                           <<"G4_Rn"
                           <<"G4_Fr"
                           <<"G4_Ra"
                           <<"G4_Ac"
                           <<"G4_Th"
                           <<"G4_Pa"
                           <<"G4_U"
                           <<"G4_Np"
                           <<"G4_Pu"
                           <<"G4_Am"
                           <<"G4_Cm"
                           <<"G4_Bk"
                           <<"G4_Cf"
                           <<"====================================="
                           <<"G4_A-150_TISSUE"
                           <<"G4_ACETONE"
                           <<"G4_ACETYLENE"
                           <<"G4_ADENINE"
                           <<"G4_ADIPOSE_TISSUE_ICRP"
                           <<"G4_AIR"
                           <<"G4_ALANINE"
                           <<"G4_ALUMINUM_OXIDE"
                           <<"G4_AMBER"
                           <<"G4_AMMONIA"
                           <<"G4_ANILINE"
                           <<"G4_ANTHRACENE"
                           <<"G4_B-100_BONE"
                           <<"G4_BAKELITE"
                           <<"G4_BARIUM_FLUORIDE"
                           <<"G4_BARIUM_SULFATE"
                           <<"G4_BENZENE"
                           <<"G4_BERYLLIUM_OXIDE"
                           <<"G4_BGO"
                           <<"G4_BLOOD_ICRP"
                           <<"G4_BONE_COMPACT_ICRU"
                           <<"G4_BONE_CORTICAL_ICRP"
                           <<"G4_BORON_CARBIDE"
                           <<"G4_BORON_OXIDE"
                           <<"G4_BRAIN_ICRP"
                           <<"G4_BUTANE"
                           <<"G4_N-BUTYL_ALCOHOL"
                           <<"G4_C-552"
                           <<"G4_CADMIUM_TELLURIDE"
                           <<"G4_CADMIUM_TUNGSTATE"
                           <<"G4_CALCIUM_CARBONATE"
                           <<"G4_CALCIUM_FLUORIDE"
                           <<"G4_CALCIUM_OXIDE"
                           <<"G4_CALCIUM_SULFATE"
                           <<"G4_CALCIUM_TUNGSTATE"
                           <<"G4_CARBON_DIOXIDE"
                           <<"G4_CARBON_TETRACHLORIDE"
                           <<"G4_CELLULOSE_CELLOPHANE"
                           <<"G4_CELLULOSE_BUTYRATE"
                           <<"G4_CELLULOSE_NITRATE"
                           <<"G4_CERIC_SULFATE"
                           <<"G4_CESIUM_FLUORIDE"
                           <<"G4_CESIUM_IODIDE"
                           <<"G4_CHLOROBENZENE"
                           <<"G4_CHLOROFORM"
                           <<"G4_CONCRETE"
                           <<"G4_CYCLOHEXANE"
                           <<"G4_1,2-DICHLOROBENZENE"
                           <<"G4_DICHLORODIETHYL_ETHER"
                           <<"G4_1,2-DICHLOROETHANE"
                           <<"G4_DIETHYL_ETHER"
                           <<"G4_N,N-DIMETHYL_FORMAMIDE"
                           <<"G4_DIMETHYL_SULFOXIDE"
                           <<"G4_ETHANE"
                           <<"G4_ETHYL_ALCOHOL"
                           <<"G4_ETHYL_CELLULOSE"
                           <<"G4_ETHYLENE"
                           <<"G4_EYE_LENS_ICRP"
                           <<"G4_FERRIC_OXIDE"
                           <<"G4_FERROBORIDE"
                           <<"G4_FERROUS_OXIDE"
                           <<"G4_FERROUS_SULFATE"
                           <<"G4_FREON-12"
                           <<"G4_FREON-12B2"
                           <<"G4_FREON-13"
                           <<"G4_FREON-13B1"
                           <<"G4_FREON-13I1"
                           <<"G4_GADOLINIUM_OXYSULFIDE"
                           <<"G4_GALLIUM_ARSENIDE"
                           <<"G4_GEL_PHOTO_EMULSION"
                           <<"G4_Pyrex_Glass"
                           <<"G4_GLASS_LEAD"
                           <<"G4_GLASS_PLATE"
                           <<"G4_GLUTAMINE"
                           <<"G4_GLYCEROL"
                           <<"G4_GUANINE"
                           <<"G4_GYPSUM"
                           <<"G4_N-HEPTANE"
                           <<"G4_N-HEXANE"
                           <<"G4_KAPTON"
                           <<"G4_LANTHANUM_OXYBROMIDE"
                           <<"G4_LANTHANUM_OXYSULFIDE"
                           <<"G4_LEAD_OXIDE"
                           <<"G4_LITHIUM_AMIDE"
                           <<"G4_LITHIUM_CARBONATE"
                           <<"G4_LITHIUM_FLUORIDE"
                           <<"G4_LITHIUM_HYDRIDE"
                           <<"G4_LITHIUM_IODIDE"
                           <<"G4_LITHIUM_OXIDE"
                           <<"G4_LITHIUM_TETRABORATE"
                           <<"G4_LUNG_ICRP"
                           <<"G4_M3_WAX"
                           <<"G4_MAGNESIUM_CARBONATE"
                           <<"G4_MAGNESIUM_FLUORIDE"
                           <<"G4_MAGNESIUM_OXIDE"
                           <<"G4_MAGNESIUM_TETRABORATE"
                           <<"G4_MERCURIC_IODIDE"
                           <<"G4_METHANE"
                           <<"G4_METHANOL"
                           <<"G4_MIX_D_WAX"
                           <<"G4_MS20_TISSUE"
                           <<"G4_MUSCLE_SKELETAL_ICRP"
                           <<"G4_MUSCLE_STRIATED_ICRU"
                           <<"G4_MUSCLE_WITH_SUCROSE"
                           <<"G4_MUSCLE_WITHOUT_SUCROSE"
                           <<"G4_NAPHTHALENE"
                           <<"G4_NITROBENZENE"
                           <<"G4_NITROUS_OXIDE"
                           <<"G4_NYLON-8062"
                           <<"G4_NYLON-6-6"
                           <<"G4_NYLON-6-10"
                           <<"G4_NYLON-11_RILSAN"
                           <<"G4_OCTANE"
                           <<"G4_PARAFFIN"
                           <<"G4_N-PENTANE"
                           <<"G4_PHOTO_EMULSION"
                           <<"G4_PLASTIC_SC_VINYLTOLUENE"
                           <<"G4_PLUTONIUM_DIOXIDE"
                           <<"G4_POLYACRYLONITRILE"
                           <<"G4_POLYCARBONATE"
                           <<"G4_POLYCHLOROSTYRENE"
                           <<"G4_POLYETHYLENE"
                           <<"G4_MYLAR"
                           <<"G4_PLEXIGLASS"
                           <<"G4_POLYOXYMETHYLENE"
                           <<"G4_POLYPROPYLENE"
                           <<"G4_POLYSTYRENE"
                           <<"G4_TEFLON"
                           <<"G4_POLYTRIFLUOROCHLOROETHYLENE"
                           <<"G4_POLYVINYL_ACETATE"
                           <<"G4_POLYVINYL_ALCOHOL"
                           <<"G4_POLYVINYL_BUTYRAL"
                           <<"G4_POLYVINYL_CHLORIDE"
                           <<"G4_POLYVINYLIDENE_CHLORIDE"
                           <<"G4_POLYVINYLIDENE_FLUORIDE"
                           <<"G4_POLYVINYL_PYRROLIDONE"
                           <<"G4_POTASSIUM_IODIDE"
                           <<"G4_POTASSIUM_OXIDE"
                           <<"G4_PROPANE"
                           <<"G4_lPROPANE"
                           <<"G4_N-PROPYL_ALCOHOL"
                           <<"G4_PYRIDINE"
                           <<"G4_RUBBER_BUTYL"
                           <<"G4_RUBBER_NATURAL"
                           <<"G4_RUBBER_NEOPRENE"
                           <<"G4_SILICON_DIOXIDE"
                           <<"G4_SILVER_BROMIDE"
                           <<"G4_SILVER_CHLORIDE"
                           <<"G4_SILVER_HALIDES"
                           <<"G4_SILVER_IODIDE"
                           <<"G4_SKIN_ICRP"
                           <<"G4_SODIUM_CARBONATE"
                           <<"G4_SODIUM_IODIDE"
                           <<"G4_SODIUM_MONOXIDE"
                           <<"G4_SODIUM_NITRATE"
                           <<"G4_STILBENE"
                           <<"G4_SUCROSE"
                           <<"G4_TERPHENYL"
                           <<"G4_TESTIS_ICRP"
                           <<"G4_TETRACHLOROETHYLENE"
                           <<"G4_THALLIUM_CHLORIDE"
                           <<"G4_TISSUE_SOFT_ICRP"
                           <<"G4_TISSUE_SOFT_ICRU-4"
                           <<"G4_TISSUE-METHANE"
                           <<"G4_TISSUE-PROPANE"
                           <<"G4_TITANIUM_DIOXIDE"
                           <<"G4_TOLUENE"
                           <<"G4_TRICHLOROETHYLENE"
                           <<"G4_TRIETHYL_PHOSPHATE"
                           <<"G4_TUNGSTEN_HEXAFLUORIDE"
                           <<"G4_URANIUM_DICARBIDE"
                           <<"G4_URANIUM_MONOCARBIDE"
                           <<"G4_URANIUM_OXIDE"
                           <<"G4_UREA"
                           <<"G4_VALINE"
                           <<"G4_VITON"
                           <<"G4_WATER_VAPOR"
                           <<"G4_XYLENE"
                           <<"G4_GRAPHITE"
                           <<"G4_WATER"
                           <<"G4_lH2"
                           <<"G4_lN2"
                           <<"G4_lO2"
                           <<"G4_lAr"
                           <<"G4_lBr"
                           <<"G4_lKr"
                           <<"G4_lXe"
                           <<"G4_PbWO4"
                           <<"G4_Galactic"
                           <<"G4_GRAPHITE_POROUS"
                           <<"G4_LUCITE"
                           <<"G4_BRASS"
                           <<"G4_BRONZE"
                           <<"G4_STAINLESS-STEEL"
                           <<"G4_CR39"
                           <<"G4_OCTADECANOL"
                           <<"====================================="
                           <<"G4_KEVLAR"
                           <<"G4_DACRON"
                           <<"G4_NEOPRENE"
                           <<"====================================="
                           <<"G4_CYTOSINE"
                           <<"G4_THYMINE"
                           <<"G4_URACIL"
                           <<"G4_DEOXYRIBOSE"
                           <<"G4_DNA_DEOXYRIBOSE"
                           <<"G4_DNA_PHOSPHATE"
                           <<"G4_DNA_ADENINE"
                           <<"G4_DNA_GUANINE"
                           <<"G4_DNA_CYTOSINE"
                           <<"G4_DNA_THYMINE"
                           <<"G4_DNA_URACIL"
                           );

    on_checkBoxUseMacroCommandFile_stateChanged(0);
    on_checkBoxWorldConst_clicked(true);
    //ui->lineEditNumberOfEvent->setText("1000");
    //ui->tabWidgetTerm->removeTab(1);
    //ui->tabWidgetRootAnalysis->removeTab(1);
    on_checkBoxRocks_clicked(false);
    ui->pushButtonLoadExe->setVisible(false);
    ui->pushButtonGenerateExe->setVisible(false);
    ui->RESUBButton_2->setVisible(false);
    ui->pushButton_2LoadICRPSpectrum->setVisible(false);

    highlighter1 = new Highlighter(ui->outputTextConsole->document());
    highlighter2 = new Highlighter(ui->GeometryFileTextEdit->document());
    ui->frame_19->setVisible(false);

    FillCoomponentByDefaultData();

    ui->tabWidget->setCurrentIndex(1);
    FillComponentsFromInputsFile(MacroFilePath); // this is called here after making and filling combobox lists

    on_checkBoxHalflivesNucl_clicked(false);
    on_checkBoxEffDoseLimNucl_clicked(false);

    ui->comboBoxNohupFiles->setVisible(false);
    ui->pushButtonShowOutputs->setVisible(false);

    if(ui->radioButtonDICOM->isChecked() || ui->radioButtonVoxel->isChecked() || ui->radioButtonVoxIDs->isChecked()){
        ui->checkBoxVoxelOrRegionLevel->setVisible(true);
        ui->checkBoxVoxelOrRegionLevel->setChecked(true);
    }else{
        if(ui->comboBoxPreDefinedGeom->currentText() == "MyGeometry"){
            ui->checkBoxVoxelOrRegionLevel->setVisible(true);
            ui->checkBoxVoxelOrRegionLevel->setChecked(true);
        }
    }

    ui->checkBoxInterpolationType->setChecked(true);
    setCompleters();

    ui->scrollArea_2->setVisible(false);
    ui->Tab->setCurrentIndex(1);
    BashCommandsForExecuting = "#! /bin/bash \n "
                               "cd "+DoseCalcsCore_build_dir_path+"\n"
                               "bash \n ";
    fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);
    ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);

    updateApplicationTabs();

}
MainWindow::~MainWindow()
{

    QProcess::startDetached("kill", {QString::number(pidOfPreviousNotKilled)}); // because the last xterm process not killed , then it will be killed at program terminating

    CoreProcess.kill();

    delete ui;
}

void MainWindow::UsePackagesMethods(){


    QMap<QString,QVector<QString>> TableWords = fileManagerObject->ReadLinesFromFileWithFirstWordIndicatorAndVector("/home/tarik/Bureau/DoseCalcs_Files/Scripts/TextChangeFromToWords.txt");
    QString Text = fileManagerObject->ReadTextFromFileInOneString("/home/tarik/Bureau/DoseCalcs_Files/Scripts/TextToChange.txt");

    for ( auto Abeg = TableWords.begin(); Abeg != TableWords.end(); ++Abeg  ){

        for (int ff = 0 ; ff < Abeg.value().size(); ff++) {
            //showResultsOutput(Abeg.value()[ff] + " -  key:"+ Abeg.key());

            Text.replace(Abeg.value()[ff],Abeg.key());
        }
    }

    showResultsOutput(Text , 4);
    fileManagerObject->WriteTextToFile("/home/tarik/Bureau/DoseCalcs_Files/Scripts/TextResultedWithChanges.txt",Text);

}

//called from MainWindow() and ...
void MainWindow::initializeVariable(){

    showResultsOutput("Initializing variables...", 0);


    //process = new QProcess(this);
    //process->setReadChannel(QProcess::StandardOutput);

    //connect(Process, &QProcess::readyReadStandardOutput, [process, this](){ processOutput()});
    //connect(process, &QProcess::readyReadStandardError, [=](){ ui->outputTextConsole->appendPlainText(process->readAllStandardError()); });

    connect (&CoreProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(processOutput()));  // connect process signals with your code
    connect (&CoreProcess, SIGNAL(readyReadStandardError()), this, SLOT(processOutput()));

    //connect(ui->customPlot, SIGNAL(plottableClick()), this, SLOT(graphClicked()));

    //connect (process, SIGNAL(readyReadStandardOutput()), this , SLOT(processOutput()));  // connect process signals with your code
    //QObject::connect (process, SIGNAL(readyReadStandardError()), this, SLOT(processOutput()));  // same here

    NumberOfExecutions = 1;

    Path1 = "user1@10.13.1.17:/home/user1/DoseCalcs/Results";
    Path2 = "user2@10.13.1.17:/home/user2/DoseCalcs";

    isInexec = true;

    IsICRPFilesAreRead = false;

    fileManagerObject = new filesManager;

    BashCommandsForExecuting = "";

    download_Directory_Path="";

    autoSimuFileName = "autoSimu.sh";
    autoMacrosFileName = "autoMacros.mac";
    autoMacrosFilePath="";

    DEFAULT_INPUTS = "";

    CENTOS_ROCKS_CLUSTER = "";
    MPI_USE = "";
    ROOT_USE = "";
    DCMTK_USE = "";
    VERBOSE_USE = "";
    GDML_USE="";

    // World
    Geometry_setWorldHalfSize = "";
    Geometry_setWorldMaterialName = "";
    Geometry_setWorldHalfSizeunit = "";

    // Geometry

    Geometry_CreateVolume_GeometryFileType = "";
    //Geometry_CreateVolume_GeometryPath = "";

    // Physics
    SourceData_setParticleName = "";
    PhysicsData_setPhysicsName = "";

    PhysicsData_setPhotoElectricEffectModel = "";
    PhysicsData_setComptonScatteringModel = "";
    PhysicsData_setGammaConversionModel = "";
    PhysicsData_setRayleighScatteringModel = "";
    PhysicsData_setElectronIonisationModel = "";
    PhysicsData_setElectronBremModel = "";
    PhysicsData_setHadronIonisationModel = "";

    PhysicsData_setCutsEnergy = "";
    PhysicsData_setCutsDistance = "";

    //PhysicsData_setIonIonisationModel = "";

    // source
    SourceData_setSourceType = "";
    SourceData_setSourceData = "";
    SourceData_setEnergyDistribution = "";
    SourceData_setEnergyData = "";
    SourceData_setAngleDistribution = "";
    SourceData_setMomDirData = "";
    SourceData_setEventsNumForDataGen = "";
    SourceData_GeneratePositions = "";
    SourceData_GenerateEnergies = "";
    SourceData_GenerateMomDirs = "";
    SourceData_setSourceSizeUnit = "";
    SourceData_setSourceEnergyUnit = "";
    SourceData_setSourceAngleUnit = "";

    // Run and score
    Score_setVolumesToScore = "";
    Score_setVariableToScore = "";
    //Score_setAccuracyCalculationLevel = "";
    Score_setSimNumOnRanksLineEdit = "";
    Score_setRadioNucleidDataLineEdit = "";
    //Score_setRadiationFactors = "";
    Score_setTissueFactors = "";

    Score_setQuantitiesUnits = "";

    Score_setRadioNucleidBiokineticsLineEdit = "";
    Execution_setNumberOfRanksOrThreads = "";
    Execution_setEventNumber = "1000";

    // Root analysis
    Analysis_setGraphsData = "";
    Analysis_setCompareType = "";
    Analysis_setGraphsExt = "";
    Analysis_setRefFilePath = "";
    Analysis_setRefName = "";
    Analysis_setRegionVariableName = "";
    Analysis_GenerateRelativeSDevGraph = "" ;
    Analysis_GenerateRelativeErrGraph = "" ;
    Analysis_GenerateRegionsVariableGraph = "" ;
    Analysis_GenerateEventsDataHisto = "" ;
    Analysis_setSliceFor2DGraph = "" ;
    Analysis_setBeamAxis = "" ;
    Analysis_setSliceID = "" ;

    Analysis_UseLogE = "" ;
    Analysis_UseLogVariable = "" ;
    Analysis_UseGridXY = "" ;
    Analysis_PrintTitle = "" ;
    Analysis_LegendPos = "" ;
    Analysis_LegendWidth ="";
    Analysis_AddErrorBar ="";

    PhysicsData_ParticleForCrossSection = "";
    PhysicsData_EUnitForCrossSection = "";
    PhysicsData_EnergiesForCrossSection = "";


    PhysicsCutEnergyUnit = "KeV";
    PhysicsCutDistanceUnit = "mm";
    SourceEnergyUnit = "MeV";
    SourceMomDirAngleUnit = "degree";
    SourceBoxHalfSizeUnit = "cm";
    AnalysisRegionEnergyUnit = "MeV";

    ExecutionCommand = "" ;

    FirstTimeAdd = true;

    if(OSNameVersion.toLower().contains("centos")){
        ui->checkBoxRocks->setChecked(true);
    }
}
void MainWindow::CommandsInitialization(){

    // !!!!!!!!!!!!!!!!!!!!! if you add a new command to a vector category, you have to add it at the last index.

    // Material commands
    MaterialCommands.push_back("/MaterialData/createElement");
    MaterialCommands.push_back("/MaterialData/createMaterial");
    MaterialCommands.push_back("/MaterialData/addElements");
    MaterialCommands.push_back("/MaterialData/addMaterial");
    MaterialCommands.push_back("/MaterialData/setNistMaterialNameAndID");

    // Geometry commands
    GeometryCommands.push_back("/GeometryData/createWorld");
    GeometryCommands.push_back("/GeometryData/createSolid");
    GeometryCommands.push_back("/GeometryData/createVolume");
    GeometryCommands.push_back("/GeometryData/setVolumesVisToNotForced");
    GeometryCommands.push_back("/GeometryData/setGeometrySymbol");

    CONSCommands.push_back("/GeometryData/createSolid");
    CONSCommands.push_back("/GeometryData/createVolume");
    CONSCommands.push_back("/GeometryData/setVolumesVisToNotForced");

    //GeometryCommands.push_back("/GeometryData/setGeometryFileType");
    //GeometryCommands.push_back("/GeometryData/geometryFileOrDir");
    //GeometryCommands.push_back("/GeometryData/setGeometryPath");
    //GeometryCommands.push_back("/GeometryData/setRegionsMassDataPath");

    // VOXEL commands
    VOXELCommands.push_back("/GeometryData/setVoxelsData");
    VOXELCommands.push_back("/GeometryData/setVoxContainerPos");
    VOXELCommands.push_back("/GeometryData/setVoxContainerRot");
    VOXELCommands.push_back("/GeometryData/setVoxDefaultMaterialName");

    VOXELCommands.push_back("/GeometryData/setVoxelizedRegionData");
    VOXELCommands.push_back("/GeometryData/setCtDensityValues");
    VOXELCommands.push_back("/GeometryData/setMatNumberAndMaterials");
    VOXELCommands.push_back("/GeometryData/setDcmTypeAndDataPath");
    VOXELCommands.push_back("/GeometryData/setDcmMaterialName");
    VOXELCommands.push_back("/GeometryData/visualizeVoxelizedPlanes");
    VOXELCommands.push_back("/GeometryData/setTETPhantomLimits");
    VOXELCommands.push_back("/GeometryData/setMaterialNameAsRegionName");
    VOXELCommands.push_back("/GeometryData/setTETRegionData");

    // Physics commands

    PhysicsCommands.push_back("/PhysicsData/setPhysicsData");
    PhysicsCommands.push_back("/PhysicsData/setCutsInRange");
    PhysicsCommands.push_back("/PhysicsData/generateCrossSectionFor");
    PhysicsCommands.push_back("/PhysicsData/setEnergyRange");

    // Source commands

    SourceCommands.push_back("/SourceData/setEventsParticleNameData");
    SourceCommands.push_back("/SourceData/setEventsInitialPosData");
    SourceCommands.push_back("/SourceData/setEventsInitialEneData");
    SourceCommands.push_back("/SourceData/setEventsInitialMomDirData");
    SourceCommands.push_back("/SourceData/useDataGenerationFiles");
    SourceCommands.push_back("/SourceData/testEventsInitialPositions");
    SourceCommands.push_back("/SourceData/showSourceBox");

    // Run and Score

    RunAndScoreCommands.push_back("/RunAndScoreData/setVolumesToScore");
    RunAndScoreCommands.push_back("/RunAndScoreData/setQuantitiesToScore");
    RunAndScoreCommands.push_back("/RunAndScoreData/setNumberOfThreads");
    RunAndScoreCommands.push_back("/RunAndScoreData/setAccuracyCalculationLevel");
    RunAndScoreCommands.push_back("/RunAndScoreData/setEventNumberPerThread");
    RunAndScoreCommands.push_back("/RunAndScoreData/setSimNumOnRanks");
    RunAndScoreCommands.push_back("/RunAndScoreData/setRadioTracerData");
    RunAndScoreCommands.push_back("/RunAndScoreData/setRadioTracerBiokinetic");
    RunAndScoreCommands.push_back("/RunAndScoreData/setQuantitiesUnits");
    RunAndScoreCommands.push_back("/RunAndScoreData/setRadiationFactors");
    RunAndScoreCommands.push_back("/RunAndScoreData/setResultDirectoryPath");
    RunAndScoreCommands.push_back("/RunAndScoreData/setTissueFactors");
    RunAndScoreCommands.push_back("/RunAndScoreData/generateVoxelsResults");
    RunAndScoreCommands.push_back("/RunAndScoreData/RunFor");

    // Analysis commands

    AnalysisCommands.push_back("/AnalysisData/generateSelfCrossGraphs");
    AnalysisCommands.push_back("/AnalysisData/generateRelativeErrGraph");
    AnalysisCommands.push_back("/AnalysisData/generateRelativeSDevGraph");
    AnalysisCommands.push_back("/AnalysisData/generateVariableRegionGraph");
    AnalysisCommands.push_back("/AnalysisData/generateEventsDataHisto");
    AnalysisCommands.push_back("/AnalysisData/setSliceFor2DGraph");
    AnalysisCommands.push_back("/AnalysisData/setBeamAxis");
    AnalysisCommands.push_back("/AnalysisData/setSliceID");
    AnalysisCommands.push_back("/AnalysisData/setGraphsParameters");
    AnalysisCommands.push_back("/AnalysisData/generateVoxelizedHistograms");

    /*
    for(int cc=0; cc < MaterialCommands.size();cc++){
        QTextStream(stdout) << " QStringLiteral(\"\\\\b" << MaterialCommands[cc] << "\\\\b\"), ";
    }
    for(int cc=0; cc < GeometryCommands.size();cc++){
        QTextStream(stdout) << " QStringLiteral(\"\\\\b" << GeometryCommands[cc] << "\\\\b\"), ";
    }
    for(int cc=0; cc < VOXELCommands.size();cc++){
        QTextStream(stdout) << " QStringLiteral(\"\\\\b" << VOXELCommands[cc] << "\\\\b\"), ";
    }
    for(int cc=0; cc < CONSCommands.size();cc++){
        QTextStream(stdout) << " QStringLiteral(\"\\\\b" << CONSCommands[cc] << "\\\\b\"), ";
    }
    for(int cc=0; cc < PhysicsCommands.size();cc++){
        QTextStream(stdout) << " QStringLiteral(\"\\\\b" << PhysicsCommands[cc] << "\\\\b\"), ";
    }
    for(int cc=0; cc < SourceCommands.size();cc++){
        QTextStream(stdout) << " QStringLiteral(\"\\\\b" << SourceCommands[cc] << "\\\\b\"), ";
    }
    for(int cc=0; cc < RunAndScoreCommands.size();cc++){
        QTextStream(stdout) << " QStringLiteral(\"\\\\b" << RunAndScoreCommands[cc] << "\\\\b\"), ";
    }
    for(int cc=0; cc < AnalysisCommands.size();cc++){
        QTextStream(stdout) << " QStringLiteral(\"\\\\b" << AnalysisCommands[cc] << "\\\\b\"), ";
    }
*/

}
void MainWindow::ReadConfigurationFile(QString ConfigFilePath){

    showResultsOutput("Reading configuration file...", 1);

    QMap <QString,QString> lines;

    QFile* filee = new QFile(ConfigFilePath);
    if(!filee->exists()){
        showResultsOutput("Cannot find the default configuration file : " + ConfigFilePath , 3);
        fileManagerObject->WriteTextToFile(ConfigFilePath,ConfigDataText);
        showResultsOutput("Check the new configuration file : " + ConfigFilePath + ", and fill it with needed paths and parameters", 1);

        return;
    }

    lines = fileManagerObject->ReadLinesFromFileWithFirstWordIndicator(ConfigFilePath);

    // supposed that are installed in /usr/local
    if(QFile::exists(lines["CMAKE_INSTALL_DIR"])){}else{lines["CMAKE_INSTALL_DIR"] = "/usr/bin";}
    if(QFile::exists(lines["GEANT4_INSTALL_DIR"])){}else{lines["GEANT4_INSTALL_DIR"] = "/usr/local/bin";}
    if(QFile::exists(lines["MPI_INSTALL_DIR"])){}else{lines["MPI_INSTALL_DIR"] = "/usr/local/bin";}
    if(QFile::exists(lines["ROOT_INSTALL_DIR"])){}else{lines["ROOT_INSTALL_DIR"] = "/usr/local/bin";}
    if(QFile::exists(lines["DCMTK_INSTALL_DIR"])){}else{lines["DCMTK_INSTALL_DIR"] = "/usr/local/lib/cmake/dcmtk";}
    if(QFile::exists(lines["DoseCalcs_SOURCE_DIR"])){}else{lines["DoseCalcs_SOURCE_DIR"] = QDir(QCoreApplication::applicationDirPath()).absolutePath()+"/core";}
    if(QFile::exists(lines["DEFAULT_DoseCalcs_INPUTS"])){
        QString s1 = "../"+GUIPackagesAndFilesDirName;
        if(lines["DEFAULT_DoseCalcs_INPUTS"].contains(s1)){ // /GeometryData/createVolume
            lines["DEFAULT_DoseCalcs_INPUTS"] = lines["DEFAULT_DoseCalcs_INPUTS"].replace(s1,GUIPackagesAndFilesDirPath);
        }
    }else{lines["DEFAULT_DoseCalcs_INPUTS"] = QDir(QCoreApplication::applicationDirPath()).absolutePath()+"/"+DoseCalcs_build_dir_name+"/macros.mac";}

    CMAKE_Lib_dir_path=lines["CMAKE_INSTALL_DIR"];
    MPI_Lib_dir_path=lines["MPI_INSTALL_DIR"];
    Root_Lib_dir_path=lines["ROOT_INSTALL_DIR"];
    DCMTK_Lib_dir_path=lines["DCMTK_INSTALL_DIR"];
    geant4_Lib_dir_path=lines["GEANT4_INSTALL_DIR"];
    DEFAULT_INPUTS=lines["DEFAULT_DoseCalcs_INPUTS"];
    MacroFilePath = DEFAULT_INPUTS ;

    OpenMacrosFilePath = DEFAULT_INPUTS;
    SaveMacrosFilePath = DEFAULT_INPUTS;

    if(QFile::exists(CMAKE_Lib_dir_path)){ // create install dir for DCMTK
        cmakeTruePath = "\n " + CMAKE_Lib_dir_path +"/cmake ";
    }
    else{
        cmakeTruePath = "\n cmake " ;
    }

    DoseCalcsCore_source_dir_path=lines["DoseCalcs_SOURCE_DIR"];
    if(!QFile::exists(DoseCalcsCore_source_dir_path)){ // create install dir for DCMTK
        DoseCalcsGui_source_dir_path = QDir::currentPath()+"/../DoseCalcs-Gui-main";
        DoseCalcsCore_source_dir_path = DoseCalcsGui_source_dir_path+"/core";
        showResultsOutput("Default DoseCalcs source directory " + DoseCalcsCore_source_dir_path + " is used", 4);
    }else{
        QDir dir = QDir(DoseCalcsCore_source_dir_path); dir.cdUp();
        DoseCalcsGui_source_dir_path = dir.path();
    }

    QDir dir = QDir(QCoreApplication::applicationDirPath());
    DoseCalcsCore_build_dir_path = dir.absolutePath()+"/"+DoseCalcs_build_dir_name;
    if(!QFile::exists(DoseCalcsCore_build_dir_path)){ // create install dir for DCMTK
        dir.mkdir(DoseCalcs_build_dir_name);
        showResultsOutput("DoseCalcs build directory " + DoseCalcsCore_build_dir_path + " is created, but you can't run until you build DoseCalcs code", 4);
    }

    UserCurrentResultsDirPath = DoseCalcsCore_build_dir_path+"/"+ResultDirectoryName;

    ui->openResultsDirButton->setToolTip("Click to choose result directory for simulation. The current directory is " + UserCurrentResultsDirPath);
    ui->pushButtonChooseResultsDir->setToolTip("Click to choose result directory for simulation. The current directory is " + UserCurrentResultsDirPath);
    ui->pushButton->setToolTip("Open " + UserCurrentResultsDirPath +" directory");

    MPI_USE = lines["MPI_USE"];
    ROOT_USE = lines["ROOT_USE"];
    DCMTK_USE = lines["DCMTK_USE"];
    VERBOSE_USE = lines["VERBOSE_USE"];
    GDML_USE = lines["GDML_USE"];
    CENTOS_ROCKS_CLUSTER = lines["CENTOS_ROCKS_CLUSTER"];

    if(MPI_USE == "YES"){
        ui->lineEditNumberOfRanksOrThreads->setPlaceholderText("Ranks Number");
        ui->lineEditNumberOfEvent->setPlaceholderText("Events Number Per Rank");
        ui->MPIOrMTOnRockscomboBox->setCurrentText("MPI");
        //ui->SimPerThreadOrRankLabel->setText("Events Number Per Rank");
    }else{
        ui->lineEditNumberOfRanksOrThreads->setPlaceholderText("Threads Number");
        ui->lineEditNumberOfEvent->setPlaceholderText("Events Number Per Thread");
        ui->MPIOrMTOnRockscomboBox->setCurrentText("MT");
        //ui->SimPerThreadOrRankLabel->setText("Events Number Per Thread");
    }

    if(CENTOS_ROCKS_CLUSTER == "YES"){
        ui->checkBoxRocks->setChecked(true);
    }

    ui->lineEditNumberOfEvent->setText(Execution_setEventNumber);

}
void MainWindow::updateApplicationTabs(){
    ui->Geometry->repaint(); ui->Geometry->update();
    ui->Physics->repaint();ui->Physics->update();
    ui->Root->repaint(); ui->Root->update();
}
void MainWindow::setCompleters(){

    QCompleter *completer;

    DefinedParticlesNames=(QStringList()
                           <<"gamma"
                           <<"alpha"
                           <<"e-"
                           <<"e+"
                           <<"proton"
                           <<"neutron");
    completer = new QCompleter(DefinedParticlesNames, this);
    completer->setCaseSensitivity(Qt::CaseInsensitive);
    ui->lineEditParticleNamesForCrossSection->setCompleter(completer);

    completer = new QCompleter(DefinedParticlesNames, this);
    completer->setCaseSensitivity(Qt::CaseInsensitive);
    ui->SourcelineEditParName->setCompleter(completer);

    completer = new QCompleter(DefinedParticlesNames, this);
    completer->setCaseSensitivity(Qt::CaseInsensitive);
    //ui->radiationEnergyFactor->setCompleter(completer);

    //completer = new QCompleter(EnergyUnits, this);
    //completer->setCaseSensitivity(Qt::CaseInsensitive);
    ui->SourceLineEditEnergyCut->setCompleter(completer);

    completer = new QCompleter(SizeUnits, this);
    completer->setCaseSensitivity(Qt::CaseInsensitive);
    ui->AnalysisLineEditVarToScore->setCompleter(completer);

    VarToScorellist=(QStringList()<<"AE"<<"AD"<<"AF"<<"SAF"<<"S"<<"H"<<"E"<<"DR"<<"DCC");
    completer = new QCompleter(VarToScorellist, this);
    completer->setCaseSensitivity(Qt::CaseInsensitive);
    ui->AnalysisLineEditVarToScore->setCompleter(completer);

    for(int aa=0; aa < MaterialCommands.size();aa++){CompleterWords << MaterialCommands[aa];}
    for(int aa=0; aa < GeometryCommands.size();aa++){CompleterWords << GeometryCommands[aa];}
    for(int aa=0; aa < VOXELCommands.size();aa++){CompleterWords << VOXELCommands[aa];}
    for(int aa=0; aa < DICOMCommands.size();aa++){CompleterWords << DICOMCommands[aa];}
    for(int aa=0; aa < CONSCommands.size();aa++){CompleterWords << CONSCommands[aa];}
    for(int aa=0; aa < PhysicsCommands.size();aa++){CompleterWords << PhysicsCommands[aa];}
    for(int aa=0; aa < SourceCommands.size();aa++){CompleterWords << SourceCommands[aa];}
    for(int aa=0; aa < RunAndScoreCommands.size();aa++){CompleterWords << RunAndScoreCommands[aa];}
    for(int aa=0; aa < AnalysisCommands.size();aa++){CompleterWords << AnalysisCommands[aa];}

    CompleterWords << "GeometryData" << "MaterialData" <<
                      "SourceData" << "RunAndScoreData" <<
                      "AnalysisData" << "PhysicsData"
                      ;

    CompleterWords << "createElement" <<
                      "createMaterial" <<
                      "addElements" <<
                      "addMaterial" <<
                      "setNistMaterialNameAndID" <<
                      "createWorld" <<
                      "createSolid" <<
                      "createVolume" <<
                      "setTETPhantomLimits" <<
                      "setMaterialNameAsRegionName" <<
                      "setTETRegionData" <<
                      "setGeometrySymbol" <<
                      "setVolumesVisToNotForced" <<
                      "setVoxelsData" <<
                      "setVoxContainerPos" <<
                      "setVoxContainerRot" <<
                      "setVoxDefaultMaterialName" <<
                      "setVoxelizedRegionData" <<
                      "setCtDensityValues" <<
                      "setMatNumberAndMaterials" <<
                      "setDcmTypeAndDataPath" <<
                      "setDcmMaterialName" <<
                      "visualizeVoxelizedPlanes" <<
                      "createSolid" <<
                      "createVolume" <<
                      "setPhysicsData" <<
                      "setCutsInRange" <<
                      "setEnergyRange" <<
                      "generateCrossSectionFor" <<
                      "setEventsParticleNameData" <<
                      "setEventsInitialPosData" <<
                      "setEventsInitialEneData" <<
                      "setEventsInitialMomDirData" <<
                      "useDataGenerationFiles" <<
                      "testEventsInitialPositions" <<
                      "showSourceBox" <<
                      "setVolumesToScore" <<
                      "setQuantitiesToScore" <<
                      "setNumberOfThreads" <<
                      "setAccuracyCalculationLevel" <<
                      "setQuantitiesUnits" <<
                      "setRadiationFactors" <<
                      "setEventNumberPerThread" <<
                      "setSimNumOnRanks" <<
                      "setRadioTracerData" <<
                      "setTissueFactors" <<
                      "setResultDirectoryPath" <<
                      "setRadioTracerBiokinetic" <<
                      "generateSelfCrossGraphs" <<
                      "generateRelativeErrGraph" <<
                      "generateRelativeSDevGraph" <<
                      "generateVariableRegionGraph" <<
                      "generateEventsDataHisto" <<
                      "setSliceFor2DGraph" <<
                      "setBeamAxis" <<
                      "setSliceID" <<
                      "setGraphsParameters" ;

    CompleterWords << "VOXEL" << "VoxIDs" <<
                      "GDML" << "STL" <<
                      "C++" << "TEXT" <<
                      "DICOM" <<

                      "gdml" << "geom" <<
                      "ast" << "stl" <<
                      "c++" <<

                      "xy" << "yz" <<
                      "yx" << "zy" <<
                      "xz" << "zx" <<

                      "Box" << "Tubs" <<
                      "CutTubs" << "Cons" <<
                      "Para" << "Trd" <<
                      "Sphere" << "Orb" <<
                      "Torus" << "Ellipsoid" <<
                      "Union" << "Intersection" <<
                      "Subtraction" <<

                      "Livermore" << "Penelope" <<
                      "EMS" << "EMS1" <<
                      "EMS2" << "EMS3" <<
                      "EMS4" << "Construct" <<

                      "numb" << "frac" <<

                      "Voxels" << "Volume" <<
                      "Point" << "Beam" <<
                      "Surface" << "Plane" << "Solid" <<
                      "AllRegions" <<"allregions" <<"null" <<

                      "Mass" << "Density" << "Volume" <<

                      "Mono" << "Rayleigh" <<
                      "File" << "Uniform" <<
                      "Gauss" << "Spectrum" <<

                      "Isotropic" << "Uniform" <<
                      "Directed" <<

                      "gamma" << "e-" <<
                      "e+" << "proton" <<
                      "alpha" << "neutron" <<

                      "SAF" << "AF" <<
                      "AE" << "AD" <<
                      "S" << "AD" <<
                      "DR" << "E" <<
                      "H" <<"DCC"<<

                      "all" << "All" <<
                      "m" << "o" <<

                      "yes" << "no" <<

                      "Result" << "Reference_Result" <<
                      "Self" << "Self_Cross" <<
                      "Cross" << "RA" <<
                      "RD" << "LRD" <<

                      "RightBottom" << "LeftBottom" <<
                      "RightTop" << "LeftTop" <<
                      "MiddleBottom" << "MiddleTop" <<

                      "XY" << "YZ" <<
                      "YX" << "ZY" <<
                      "ZX" << "XZ" <<

                      "root" << "pdf" <<
                      "ps" << "png" <<
                      "jpeg" ;

    CompleterWords << "mm" <<
                      "cm" <<
                      "m" <<
                      "g/cm3" <<
                      "mg/cm3" <<
                      "mg/mm3" <<
                      "kg/m3" <<
                      "mg/mL" <<
                      "g/mL" <<
                      "degree" <<
                      "radian" <<
                      "eV" <<
                      "keV" <<
                      "MeV" <<
                      "Bq" <<
                      "kBq" <<
                      "MBq" <<
                      "GBq" <<
                      "J" << "kg-1" << "g-1" <<
                      "MeV/kg" << "Gy" << "mGy" <<
                      "MGy" << "miGy" << "nGy" <<
                      "MeV" << "MeV" << "MeV" <<
                      "Sv" << "mSv"
                      ;

    CompleterWords << "DWITH_GDML_USE" << "DWITH_MPI_USE" <<
                      "DWITH_ANALYSIS_USE" << "DROOT_DIR" <<
                      "DWITH_DCMTK_USE" << "DDCMTK_DIR" <<
                      "DWITH_VERBOSE_USE" << "DCMAKE_CXX_COMPILER" <<
                      "DCMAKE_C_COMPILER" <<
                      "ON" <<"OFF" <<
                      "CMAKE_INSTALL_DIR" <<
                      "GEANT4_INSTALL_DIR" <<
                      "MPI_INSTALL_DIR" <<
                      "ROOT_INSTALL_DIR" <<
                      "DCMTK_INSTALL_DIR" <<
                      "DoseCalcs_SOURCE_DIR" <<
                      "DEFAULT_DoseCalcs_INPUTS" <<
                      "CMAKE_DOWNLOAD_URL" <<
                      "GEANT4_DOWNLOAD_URL" <<
                      "XERCES_DOWNLOAD_URL" <<
                      "MPI_DOWNLOAD_URL" <<
                      "ROOT_DOWNLOAD_URL" <<
                      "DCMTK_DOWNLOAD_URL" <<
                      "VERBOSE_USE" <<
                      "GDML_USE" <<
                      "CENTOS_ROCKS_CLUSTER" <<
                      "MPI_USE" <<
                      "ROOT_USE" <<
                      "DCMTK_USE"
                      ;

    completer1 = new QCompleter(this);
    completer1->setModel(new QStringListModel(CompleterWords, completer1));
    completer1->setModelSorting(QCompleter::CaseInsensitivelySortedModel);
    completer1->setCaseSensitivity(Qt::CaseInsensitive);
    completer1->setWrapAround(false);
    ui->GeometryFileTextEdit->setCompleter(completer1);

}
void MainWindow::FillCoomponentByDefaultData(){

    // World data

    ui->PhantomWorldHalfSizeslineEdit->setText("3 3 3");
    ui->comboBoxWorldSizeUnit->setCurrentText("m");

    // Geometry data

    ui->radioButtonConstruct->setChecked(true);
    ui->lineEditGeometrySymbole->setText("phantom0");

    // Physics

    ui->SourceLineEditEnergyCut->setText("1 keV 10 GeV");
    ui->SourceLineEditDistanceCut->setText("e- 1 mm e+ 1 mm gamma 1 mm proton 1 mm");

    ui->lineEditParticleNamesForCrossSection->setText("gamma");
    ui->comboBoxEnergyUnitsForCrossSection->setCurrentText("MeV");
    ui->lineEditEnergiesForCrossSection->setText("0.01 0.02 0.05 0.1 0.2 0.5 1 2");

    // radiation source

    ui->SourcelineEditParName->setText("e-");
    ui->comboBoxTypeOfSources->setCurrentText("Volume");
    ui->comboBoxSizeUnit->setCurrentText("cm");
    ui->lineEditChosenSourceTypeData->setText("Liver 15 15 15");

    ui->SourceComboBoxEnergyDist->setCurrentText("mono");
    ui->comboBoxEnergyUnit->setCurrentText("MeV");
    ui->lineEditSpecialEnergyDistributionParameter->setText("1");

    ui->SourceComboBoxAngleDist->setCurrentText("Isotropic");

    ui->UseDataFilesFor->setCurrentText("");

    // run and score

    ui->AnalysisLineEditVarToScore->setText("AE AF SAF AD S H E DR");

    ui->comboBoxAEUnits->setCurrentText("MeV");
    ui->comboBoxSAFUnits->setCurrentText("kg-1");
    ui->comboBoxADUnits->setCurrentText("mGy");
    ui->comboBoxSUnits->setCurrentText("mGy");
    ui->comboBoxHUnits->setCurrentText("mSv");
    ui->comboBoxEUnits->setCurrentText("mSv");

    ui->TissueFactorLineEdit->setText("Other 0.12 Gonads 0.08 Bladder 0.04 Esophagus 0.04 Liver 0.04 Thyroid 0.04 BoneSurface 0.1 Brain 0.1 Skin 0.1 SalivaryGlands 0.1 Kidneys 0.1");
    //ui->radiationEnergyFactor->setText("gamma 0.511 1 e+ 0.6335 1 alpha 2 20");

    ui->lineEditRadioNucleid->setText("F-18");
    ui->comboBoxRadionuclidedataType->setCurrentIndex(1);
    ui->RadioNucleidlineEmissionDataEdit->setText("e- 0.249776 96.73 gamma 0.511 193.46");

    ui->RadioNucleidAdmActivitylineEdit->setText("1");
    ui->comboBoxActivityUnits->setCurrentText("MBq");
    ui->RadioNucleidSourceResTimelineEdit->setText("Brain 0.21 HeartWall 0.11 Lungs 0.079 Liver 0.13 UBCs 0.26 allregions 1.7");
    ui->comboBoxTimeUnit->setCurrentText("h");

    ui->lineEditNumberOfRanksOrThreads->setText("2");
    ui->lineEditNumberOfEvent->setText("1000");

}

// to build the application core
void MainWindow::on_BuildButton_clicked()
{

    installationDialogObj = new InstallationDialog(this);
    installationDialogObj->show();

    if(MPI_USE == "YES"){
        ui->lineEditNumberOfRanksOrThreads->setPlaceholderText("Ranks Number");
        ui->lineEditNumberOfEvent->setPlaceholderText("Events Number Per Rank");
    }else{
        ui->lineEditNumberOfRanksOrThreads->setPlaceholderText("Threads Number");
        ui->lineEditNumberOfEvent->setPlaceholderText("Events Number Per Thread");
        //ui->SimPerThreadOrRankLabel->setText("Events Number Per Thread");
    }
}
void MainWindow::on_pushButtonCalculateNumSimulation_clicked()
{
    if(ValuesOfInputs[ui->ScoreCombobowSimNumOnRanksLineEdit->currentText()] == "o"){
        return;
    }

    int num = 1;
    QString tx = "(";

    QStringList InputsVals;

    InputsVals = ui->SourcelineEditParName->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
    num = num * InputsVals.size();

    tx = tx + "" + QString::number(num) + "*";

    if(ui->comboBoxTypeOfSources->currentText() == "Voxels" || ui->comboBoxTypeOfSources->currentText() == "TET"){
        InputsVals = ui->lineEditChosenSourceTypeData->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"

        int mines = 0;
        for(int b=0; b < InputsVals.size();b++){
            if(InputsVals[b].toLower() == "allregions"){
                if(InputsVals.size() >= b ){
                    mines = InputsVals[b+1].toInt() + 1;
                }
                break;
            }
        }

        num = num * (InputsVals.size()-mines);
        tx = tx + "" + QString::number(InputsVals.size()-mines) + "*";
    }
    else if(ui->comboBoxTypeOfSources->currentText() == "Volume"){
        InputsVals = ui->lineEditChosenSourceTypeData->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        num = num * (InputsVals.size()/4);
        tx = tx + "" + QString::number(InputsVals.size()/4) + "*";
    }

    if(ui->SourceComboBoxEnergyDist->currentText() != "File" && ui->SourceComboBoxEnergyDist->currentText() != "Spectrum"  ){
        InputsVals = ui->lineEditSpecialEnergyDistributionParameter->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        num = num * InputsVals.size();
    }else{
        num = num * 1;
    }
    tx = tx + "" + QString::number(InputsVals.size()) + "="+QString::number(num)+")";

    ui->pushButtonCalculateNumSimulation->setToolTip("Calculate the number of simulations according to the entered radiation source inputs. " + tx);
    ui->lineEditNumberOfRanksOrThreads->setText(QString::number(num));

}
void MainWindow::on_RunButton_clicked()
{
    if(ui->checkBoxSimulateGeometriesWitOneSource->isChecked()){
        RunForMultiGeomeries();
        return;
    }

    if(!TestSimulateExecutableInputsToRun()){
        return;
    }

    if(!QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName)){
        showResultsOutput("Cannot find DoseCalcs executable "+DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName + " . Please build DoseCalcs before run", 3);
        on_BuildButton_clicked();
    }

    bool pvr = false;
    QDir dir(UserCurrentResultsDirPath);
    foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "AE@*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {{pvr = true; break;}}
    if(pvr){if (QMessageBox::Yes == QMessageBox::question(this, tr("Simulation result directory"), tr("The current directory is the result of previous simulation, you can choose another directory for results, or create a new one (yes) "))){on_openResultsDirButton_clicked();}}


    initializeVariable();
    SaveDataFromInputComponents(); // get the same componenet data but the four under variable values are related to the SaveEnePharOrgLists() data that is called one time when the the run in pushed
    CreateUserCommands(); //fill the variables with all values to save it to user file and not inputFile executed by the geant4 application core

    if(ui->checkBoxUseMacroCommandFile->isChecked()){
        if(EditFlag == 8){
            on_pushButtonEditGeomFile_clicked(); // save data to MacroFileName
        }
        //fill the components by the text in the saved command file
        FillComponentsFromInputsFile(DoseCalcsCore_build_dir_path+"/"+MacroFileName); // fill data from MacroFileName

        if(!TestSimulateExecutableInputsToRun()){
            return;
        }
    }

    initializeVariable();
    SaveDataFromInputComponents(); // get the same componenet data but the four under variable values are related to the SaveEnePharOrgLists() data that is called one time when the the run in pushed
    CreateUserCommands(); //fill the variables with all values to save it to user file and not inputFile executed by the geant4 application core

    QsubSeparatedMacroFilePath = DoseCalcsCore_build_dir_path+"/Macros"+ConstructDoseCalcsJobName().remove("DoseCalcs")+".mac";
    fileManagerObject->WriteTextToFile(QsubSeparatedMacroFilePath , generateInputUserTextForinputFile());

    if(ui->checkBoxRocks->isChecked()){

        QString SpecificNodes = SetToASpecificMPIRank();

        BashCommandsForExecuting = "#! /bin/bash \ncd " + DoseCalcsCore_build_dir_path
                + "\n qsub "+ SpecificNodes +" "+ DoseCalcsCore_build_dir_path+"/"+ExeFileName;

        if(false){
            BashCommandsForExecuting = "#! /bin/bash \ncd " + DoseCalcsCore_build_dir_path
                    + "\n sbatch "+ SpecificNodes +" "+ DoseCalcsCore_build_dir_path+"/"+ExeFileName;
        }

        BashCommandsForExecuting += "\n bash \n";

        showResultsOutput("Writing Run Commands : \n", 0);
        showResultsOutput(BashCommandsForExecuting , 0);
        showResultsOutput("to --> " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , 4);
        fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

        //if(EditFlag == 1){ // if you late it open, is not saved, can be for the last execution
        on_pushButtonGenerateExe_clicked();
        on_pushButtonEditGeomFile_clicked();
        //}
        ui->outputTextConsole->setPlainText(fileManagerObject->ReadTextFromFileInOneString(DoseCalcsCore_build_dir_path+"/"+ExeFileName));
        ui->tabWidget->setCurrentIndex(1);

        if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName) && QFile::exists(DoseCalcsCore_build_dir_path+"/"+ExeFileName)){

            if(!ShowImportantSimulationData()){return;}
            showResultsOutput("Computation Run" , 1);
            ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
        }
    }else{

        if(QFile::exists(QsubSeparatedMacroFilePath) && QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName)){

            if(MPI_USE=="YES"){

                if(!QFile::exists(MPI_Lib_dir_path+"/mpirun")){
                    MPI_Lib_dir_path = ShowMessageBoxAndGetPath("Directory containing mpirun or mpiexec Not Found, Click OK to Choose the Directory");
                }

                QString nohup = "";
                QString andd = "";
                if(ui->checkBoxnohup->isChecked()){
                    LastRunOutputFile = "nohup_"+ConstructDoseCalcsJobName();

                    if(!QFile::exists("/usr/bin/nohup")){
                        QMessageBox::information(this, tr(""), "/usr/bin/nohup Not Found, Please install nohup and check the /usr/bin/nohup path.");
                        return;
                    }
                    nohup = "/usr/bin/nohup ";
                    andd = "  > " + LastRunOutputFile + " & ";
                }

                BashCommandsForExecuting = "#! /bin/bash \ncd " + DoseCalcsCore_build_dir_path + "\n"
                        + nohup + " " + MPI_Lib_dir_path+"/mpirun " + Execution_setNumberOfRanksOrThreads + " "
                        + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName + " B " + QsubSeparatedMacroFilePath + " " + Execution_setEventNumber + " " +andd
                        ;
            }else{

                QString f = "";

                if(!QFile::exists(geant4_Lib_dir_path+"/geant4.sh")){
                    geant4_Lib_dir_path = ShowMessageBoxAndGetPath("Directory containing geant4.sh Not Found, Click OK to Choose the Directory");
                }

                QString nohup = "";
                QString andd = "";
                if(ui->checkBoxnohup->isChecked()){
                    ui->RunButton->setText("Run("+QString::number(macrosfileinc)+")");macrosfileinc++;

                    LastRunOutputFile = "nohup_"+ConstructDoseCalcsJobName();

                    if(!QFile::exists("/usr/bin/nohup")){
                        QMessageBox::information(this, tr(""), "/usr/bin/nohup Not Found, Please install nohup and check the /usr/bin/nohup path.");
                        return;
                    }
                    nohup = "/usr/bin/nohup ";
                    andd = "  > " + LastRunOutputFile + " & ";
                }

                BashCommandsForExecuting = "#! /bin/bash \n" +f
                        + "cd " +geant4_Lib_dir_path +"\n"+
                        + ". ./geant4.sh\n" +f +
                        + "cd " + DoseCalcsCore_build_dir_path + "\n"
                        + nohup + " "+DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName + " B " + QsubSeparatedMacroFilePath + " " + Execution_setEventNumber + " " +andd
                        ;
            }

            BashCommandsForExecuting += "\n bash \n";

            showResultsOutput("Writing Run Commands : \n", 0);
            showResultsOutput(BashCommandsForExecuting , 0);
            showResultsOutput("to --> " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , 4);

            fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

            if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName)){

                if(!ShowImportantSimulationData()){return;}
                showResultsOutput("Computation Run" , 1);
                ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
            }
            else{
                showResultsOutput("Cannot find file containing execution commands "+ DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName + " , you should build DoseCalcs with ROOT Analysis option" , 3);
            }

            isInexec = true;

        }else{
            showResultsOutput("Verify that the simulate binary and the inputs.mac are existed in " + QsubSeparatedMacroFilePath, 0);
        }
    }

    ui->RunButton->setText("Run("+QString::number(macrosfileinc)+")");macrosfileinc++;
    if (ui->checkBoxnohup->isChecked()){
        ui->comboBoxNohupFiles->addItem(LastRunOutputFile);
        ui->comboBoxNohupFiles->setCurrentText(LastRunOutputFile);
    }
}
void MainWindow::on_checkBoxnohup_clicked()
{
    if (ui->checkBoxnohup->isChecked()){
        ui->comboBoxNohupFiles->setVisible(true);
        ui->pushButtonShowOutputs->setVisible(true);
        QDir dir(DoseCalcsCore_build_dir_path);
        ui->comboBoxNohupFiles->clear();
        foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "nohup_DoseCalcs*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {
            ui->comboBoxNohupFiles->addItem(entry.fileName());
        }
    }else{
        ui->comboBoxNohupFiles->setVisible(false);
        ui->pushButtonShowOutputs->setVisible(false);
        QDir dir(DoseCalcsCore_build_dir_path);
        ui->comboBoxNohupFiles->clear();
    }
}
void MainWindow::on_pushButtonShowOutputs_clicked()
{

    BashCommandsForExecuting = "#! /bin/bash \n cd " + DoseCalcsCore_build_dir_path+ "\n tail -f "+ ui->comboBoxNohupFiles->currentText() +" -n +0";
    BashCommandsForExecuting += "\n bash \n";

    showResultsOutput("Writing Check Commands : \n", 0);
    showResultsOutput(BashCommandsForExecuting , 0);
    showResultsOutput("to --> " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , 4);
    fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);
    showResultsOutput("Checking DoseCalcs simulation output file" , 1);
    ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
}
void MainWindow::on_pushButtonShowOutputsAndMacros_clicked()
{
    QDir dir(DoseCalcsCore_build_dir_path);
    QString filesnames = "\n";
    foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "DoseCalcs_*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {
        filesnames += entry.fileName() + "\n";
    }
    foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "nohup_*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {
        filesnames += entry.fileName() + "\n";
    }
    foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "Macros_*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {
        filesnames += entry.fileName() + "\n";
    }

    if (QMessageBox::Yes == QMessageBox::question(this, "Warning!", "Do you want to remove these files: \n"+filesnames , QMessageBox::Yes | QMessageBox::No))
    {

        foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "DoseCalcs_*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {
            dir.remove(entry.fileName());
        }
        foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "nohup_*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {
            dir.remove(entry.fileName());
        }
        foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "Macros_*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {
            dir.remove(entry.fileName());
        }
    }

}
void MainWindow::on_pushButtonStopCurrentProcess_clicked()
{
    if (QMessageBox::Yes == QMessageBox::question(this, "Warning!", "Do you want to stop all DoseCalcs simulations in progress (\"simulate\" processes)? in case there is more than one simulation, you should verify the current simulations in progress \n1) Type \"top\" in terminal;\n"
                                                  "2) Then take the process_ID of the \"simulate\" process to be killed;\n3) execute \"kill process_ID\" on terminal.\n\n"
                                                  "If you click on \"yes\", all \"simulate\" processes will be killed" , QMessageBox::Yes | QMessageBox::No))
    {
        QProcess p;
        p.start("pkill simulate");
        p.waitForFinished();
    }
}
void MainWindow::on_openResultsDirButton_clicked()
{
    QFileDialog dialog;
    dialog.setFileMode(QFileDialog::Directory);
    dialog.setOption(QFileDialog::ShowDirsOnly, false);

    QString pp = dialog.getExistingDirectory(0, ("Choose the results directory where ResultsData, Ref1 and Ref2 can be found, or create one to be used in simulations"), UserCurrentResultsDirPath);
    if(pp != "" || !pp.isEmpty()){
        UserCurrentResultsDirPath = pp;
    }

    ui->openResultsDirButton->setToolTip("Click to choose result directory for simulation. The current directory is " + UserCurrentResultsDirPath);
    ui->pushButtonChooseResultsDir->setToolTip("Click to choose result directory for simulation. The current directory is " + UserCurrentResultsDirPath);
    ui->pushButton->setToolTip("Open " + UserCurrentResultsDirPath +" directory");

}
void MainWindow::on_pushButtonChooseResultsDir_clicked()
{
    QFileDialog dialog;
    dialog.setFileMode(QFileDialog::Directory);
    dialog.setOption(QFileDialog::ShowDirsOnly, false);

    QString pp = dialog.getExistingDirectory(0, ("Choose the results directory where ResultsData, Ref1 and Ref2 can be found, or create one to be used in simulations"), UserCurrentResultsDirPath);
    if(pp != "" || !pp.isEmpty()){
        UserCurrentResultsDirPath = pp;
    }

    ui->openResultsDirButton->setToolTip("Click to choose result directory for simulation. The current directory is " + UserCurrentResultsDirPath);
    ui->pushButtonChooseResultsDir->setToolTip("Click to choose result directory for simulation. The current directory is " + UserCurrentResultsDirPath);
    ui->pushButton->setToolTip("Open " + UserCurrentResultsDirPath +" directory");

}
void MainWindow::on_pushButtonMerge_clicked()
{

    if(!TestMergeExecutableInputsToRun()){
        return;
    }

    //UserCurrentResultsDirPath = QFileDialog::getExistingDirectory(0, ("Choose Results Directory"), DoseCalcsCore_build_dir_path + "/" + ResultDirectoryName);
    if(UserCurrentResultsDirPath == "" || UserCurrentResultsDirPath.isEmpty()){UserCurrentResultsDirPath = DoseCalcsCore_build_dir_path + "/" + ResultDirectoryName;}

    QString MacrosFileForAnalysis;

    if(QFile::exists(MacroFilePath)){
        MacrosFileForAnalysis = MacroFilePath;
    }else{
        MacrosFileForAnalysis = DoseCalcsCore_build_dir_path+"/"+MacroFileName;
    }


    if(QFile(UserCurrentResultsDirPath+"/"+ResultFileName).size() != 0){
        QDir dir(UserCurrentResultsDirPath);
        QString nn = "The "+UserCurrentResultsDirPath+"/"+ResultFileName+" not empty. Click on \"yes\" if you want erasing its data before generating the new results data";
        if(QMessageBox::Yes == QMessageBox::question(this, tr(""),nn)){
            dir.remove("ResultsData");
        }
    }

    initializeVariable();
    SaveDataFromInputComponents(); // get the same componenet data but the four under variable values are related to the SaveEnePharOrgLists() data that is called one time when the the run in pushed
    CreateUserCommands(); //fill the variables with all values to save it to user file and not inputFile executed by the geant4 application core
    fileManagerObject->WriteTextToFile( MacrosFileForAnalysis , generateInputUserTextForinputFile());

    if(!QFile::exists(MacrosFileForAnalysis)){
        QMessageBox::information(this, tr(""), "Cannot find macros file for data files merging in this directory." );
        MacrosFileForAnalysis = QFileDialog::getOpenFileName( this, tr("Open a macros file for data files merging"), DoseCalcsCore_build_dir_path, "All files (*.*)" );
    }

    if(!QFile::exists(geant4_Lib_dir_path+"/geant4.sh")){
        geant4_Lib_dir_path = ShowMessageBoxAndGetPath("Directory containing geant4.sh Not Found, Click OK to Choose the Directory");
    }

    //CoreProcess.execute("sh " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
    BashCommandsForExecuting = "\ncd " +geant4_Lib_dir_path +"\n" + ". ./geant4.sh \n"+
            "cd " + DoseCalcsCore_build_dir_path + "\n" +
            DoseCalcsCore_build_dir_path+"/"+MergeExecutableName + " " + MacrosFileForAnalysis ;
    ;
    BashCommandsForExecuting += "\n bash \n";

    fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

    showResultsOutput("Writing merging Commands : \n", 0);
    showResultsOutput(BashCommandsForExecuting , 0);
    showResultsOutput("to --> " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , 0);

    if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName)){

        showResultsOutput("Merging results files " , 1);

        ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);

    }
    else{
        showResultsOutput("Cannot find file containing execution commands "+ DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName + " , you should build DoseCalcs with ROOT Analysis option" , 3);
    }

}
void MainWindow::on_pushButtonOpenResultsFile_clicked()
{

    if(UserCurrentResultsDirPath == "" || UserCurrentResultsDirPath.isEmpty()){
        UserCurrentResultsDirPath = DoseCalcsCore_build_dir_path + "/" + ResultDirectoryName;
    }

    on_pushButtonEditGeomFile_clicked();

    ui->tabWidget->setTabText(0,ResultFileName);
    ui->GeometryFileTextEdit->clear();
    showResultsOutput("getting "+ResultFileName+" data", 4);
    ui->GeometryFileTextEdit->setPlainText(fileManagerObject->ReadTextFromFileInOneString(UserCurrentResultsDirPath+"/"+ResultFileName));
    ui->tabWidget->setCurrentIndex(0);

    EditFlag = 9;

}
void MainWindow::on_pushButton_clicked()
{

    if(QFile::exists(UserCurrentResultsDirPath)){
        QString command = UserCurrentResultsDirPath;
        QProcess process;
        QStringList qsl = {command};
        process.startDetached("nautilus", qsl);
    }else{
        QMessageBox::information(this, tr(""), "Cannot find "+ UserCurrentResultsDirPath +" directory");
    }

    /*

    QDialog * d = new QDialog(); d->setWindowTitle("Create a new results directory for DoseCalcs simulations");
    QGridLayout* GraphLayout = new QGridLayout;

    QLabel* lb = new QLabel("New Results Directory Name");
    QLineEdit * IDs = new QLineEdit();
    IDs->setToolTip("Add the name of directory");

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    int ii = 0, jj=0;
    GraphLayout->addWidget(lb, jj,ii,1,1);
    GraphLayout->addWidget(IDs, jj,++ii,1,1);
    GraphLayout->addWidget(buttonBox);

    d->setLayout(GraphLayout);

    int result = d->exec();

    if(result == QDialog::Accepted)
    {

        if(IDs->text()!=""){

            if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+ResultDirectoryName+"/"+IDs->text())){
                QMessageBox::information(this, tr(""), "The directory name ("+IDs->text()+") is already exists. Choose another name for results directory !" );
            }else{

                showResultsOutput("Creating DoseCalcs Results Directory: "+ IDs->text() , 1);

                QProcess process;
                QStringList qsl = {DoseCalcsCore_build_dir_path+"/"+ResultDirectoryName+"/"+IDs->text()};
                process.startDetached("mkdir", qsl);

                UserCurrentResultsDirPath = DoseCalcsCore_build_dir_path+"/"+IDs->text();
    ui->openResultsDirButton->setToolTip("Click to choose result directory for simulation. The current directory is " + UserCurrentResultsDirPath);
        ui->pushButton->setToolTip("Open " + UserCurrentResultsDirPath +" directory");

            }
        }
    }

    */
}
void MainWindow::on_checkBoxRocks_clicked(bool checked)
{
    if (ui->checkBoxRocks->isChecked()){
        ui->checkBoxnohup->setEnabled(false);
        ui->checkBoxnohup->setChecked(false);
        ui->comboBoxNohupFiles->setVisible(true);
        ui->pushButtonShowOutputs->setVisible(true);
        QDir dir(DoseCalcsCore_build_dir_path);
        ui->comboBoxNohupFiles->clear();
        foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "DoseCalcs_*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {
            ui->comboBoxNohupFiles->addItem(entry.fileName());
        }
        foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "nohup_*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {
            ui->comboBoxNohupFiles->addItem(entry.fileName());
        }
    }else{
        ui->comboBoxNohupFiles->setVisible(false);
        ui->pushButtonShowOutputs->setVisible(false);
        QDir dir(DoseCalcsCore_build_dir_path);
        ui->comboBoxNohupFiles->clear();
    }


    if(ui->checkBoxRocks->isChecked()){
        ui->pushButtonLoadExe->setEnabled(true);
        ui->MPIOrMTOnRockscomboBox->setEnabled(true);
        ui->pushButtonGenerateExe->setEnabled(true);
        ui->pushButtonQstat->setEnabled(true);
        ui->pushButtonChecQsubMPIsim->setEnabled(true);
        ui->pushButtonStopJob->setEnabled(true);

        if(MPI_USE == "YES"){
            ui->MPIOrMTOnRockscomboBox->setCurrentText("MPI");
        }else{
            ui->MPIOrMTOnRockscomboBox->setCurrentText("MT");
        }
        timer.singleShot(60000, this,SLOT(removeHugFiles_slot())); // check in 3 minutes
    }else{
        ui->checkBoxnohup->setEnabled(true);
        ui->pushButtonLoadExe->setEnabled(false);
        ui->MPIOrMTOnRockscomboBox->setEnabled(false);
        ui->pushButtonGenerateExe->setEnabled(false);
        ui->pushButtonQstat->setEnabled(false);
        ui->pushButtonChecQsubMPIsim->setEnabled(false);
        ui->pushButtonStopJob->setEnabled(false);
        timer.singleShot(0, this,SLOT(removeHugFiles_slot()));
    }
}
void MainWindow::on_pushButtonQstat_clicked()
{

    BashCommandsForExecuting = "#! /bin/bash \ncd " + DoseCalcsCore_build_dir_path
            + "\n qstat -r | grep DoseCalcs";
    BashCommandsForExecuting += "\n bash \n";

    fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

    ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);

}
void MainWindow::on_pushButtonChecQsubMPIsim_clicked()
{

    StopOrCheck = false;

    getRockcsDoseCalcsJobs();

}
void MainWindow::on_pushButtonStopJob_clicked()
{

    StopOrCheck = true;

    getRockcsDoseCalcsJobs();

}
void MainWindow::on_pushButtonLoadExe_clicked()
{
    on_pushButtonEditGeomFile_clicked();

    ui->tabWidget->setTabText(0,ExeFileName);
    ui->GeometryFileTextEdit->clear();
    showResultsOutput("getting exe.sh data", 4);
    ui->GeometryFileTextEdit->setPlainText(fileManagerObject->ReadTextFromFileInOneString(DoseCalcsCore_build_dir_path+"/"+ExeFileName));
    ui->tabWidget->setCurrentIndex(0);
    EditFlag = 1;
}
void MainWindow::on_pushButtonGenerateExe_clicked()
{

    on_pushButtonEditGeomFile_clicked();

    QString DoseCalcsJobName = ConstructDoseCalcsJobName();

    QString numbOfRanks;
    if(ui->MPIOrMTOnRockscomboBox->currentText() == "MPI"){
        numbOfRanks = Execution_setNumberOfRanksOrThreads;
    }else{
        numbOfRanks = "1";
    }

    if(!QFile::exists(geant4_Lib_dir_path+"/geant4.sh")){
        geant4_Lib_dir_path = ShowMessageBoxAndGetPath("Directory containing geant4.sh Not Found, Click OK to Choose the Directory");
    }

    if(!QFile::exists(MPI_Lib_dir_path+"/mpirun")){
        MPI_Lib_dir_path = ShowMessageBoxAndGetPath("Directory containing mpirun or mpiexec Not Found, Click OK to Choose the Directory");
    }

    ExeDataText = "#!/bin/bash \n"
                  "#$ -S /bin/bash \n"
                  "#$ -cwd \n"
                  "#$ -j y \n"
                  "#$ -q all.q \n"
                  "#$ -N " + DoseCalcsJobName + " \n"
                                                "#$ -pe mpi "+ Execution_setNumberOfRanksOrThreads + " \n"
                                                                                                     "#$ -M imttarikk@gmail.com \n"
                                                                                                     "#$ -m e \n"
                                                                                                     ". "+geant4_Lib_dir_path+"/geant4.sh \n"+
            MPI_Lib_dir_path + "/mpirun -np "+ numbOfRanks + " " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName + " B " + QsubSeparatedMacroFilePath + " " + Execution_setEventNumber +" \n";

    if(false){
        ExeDataText = "#!/bin/bash \n"
                  "#SBATCH --job-name=" + DoseCalcsJobName + " \n"
                  "#SBATCH --ntasks=1 --nodes=1 \n"
                  "# #SBATCH --mem-per-cpu=5G \n"
                  "# #SBATCH --partition=defq \n"
                  "# #SBATCH --time=12:00:00 \n"
                  //"#$ -pe mpi "+ Execution_setNumberOfRanksOrThreads + " \n"
                  //"#$ -M imttarikk@gmail.com \n"
                  ". "+geant4_Lib_dir_path+"/geant4.sh \n"+
            //MPI_Lib_dir_path + "/mpirun -np "+ numbOfRanks + " " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName + " B " + QsubSeparatedMacroFilePath + " " + Execution_setEventNumber +" \n";
            DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName + " -p "+ numbOfRanks + " B " + QsubSeparatedMacroFilePath + " " + Execution_setEventNumber +" \n";
    }

    ui->tabWidget->setTabText(0,ExeFileName);
    ui->GeometryFileTextEdit->clear();
    showResultsOutput("getting exe.sh data", 4);
    ui->GeometryFileTextEdit->setPlainText(ExeDataText);
    ui->tabWidget->setCurrentIndex(0);
    EditFlag = 1;
}

void MainWindow::on_RESUBButton_2_clicked()
{
    QDialog * d = new QDialog(); d->setWindowTitle("Add macros file and number of simulations");
    QVBoxLayout * vbox = new QVBoxLayout();

    QPushButton * openqsubMacros = new QPushButton(); vbox->addWidget(openqsubMacros);
    openqsubMacros->setText("Choose macros file"); connect(openqsubMacros, SIGNAL(clicked()), this, SLOT(openqsubMacros_slot()));
    openqsubMacros->setToolTip("/../../macros.mac");

    QLineEdit * numbOfSim= new QLineEdit(); vbox->addWidget(numbOfSim);
    numbOfSim->setToolTip("Add number of process (each process for a simulation run)");

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));
    vbox->addWidget(buttonBox);

    d->setLayout(vbox);

    int result = d->exec();

    if(result == QDialog::Accepted)
    {

    }
}
void MainWindow::openqsubMacros_slot(){
    OpenMacrosFilePath = QFileDialog::getOpenFileName(
                this,
                tr("Open a macros file"),
                DoseCalcsCore_build_dir_path, //pathBuildApp,
                "All files (*.*);;Text files (*.txt)"
                );
}

void MainWindow::getRockcsDoseCalcsJobs()
{
    // all the list of job files
    QDir dir(DoseCalcsCore_build_dir_path);
    QStringList listOfFile ;
    listOfFilesName.clear(); listOfFilesName.empty();

    //foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "DoseCalcs.o*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {
    //QString ID = entry.fileName().remove(0,11);
    //listOfFile << ID;
    foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "DoseCalcs_*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {
        QStringList InputsVals = entry.fileName().split(".o", QString::SkipEmptyParts);
        if(InputsVals.size() > 0 ){listOfFilesName[InputsVals[0]] = InputsVals[1]; listOfFile << InputsVals[0];}

        showResultsOutput("DoseCalcs Job File ID "+ InputsVals[1] , 1);
        showResultsOutput("DoseCalcs Job File Name "+ InputsVals[0] , 1);
    }


    QStringList InputsVals;
    QVector<QString> Commlines;


    // read runing job names by file interface
    terminal *term = new terminal;
    term->scroll(20, 20);
    BashCommandsForExecuting = "#! /bin/bash \ncd " + DoseCalcsCore_build_dir_path + "\n";
    BashCommandsForExecuting += "qstat -r | grep DoseCalcs > " + DoseCalcsCore_build_dir_path+"/F" + "\n";
    BashCommandsForExecuting += "\n bash \n";
    fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);
    term->executeLocalCommand(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
    Commlines = fileManagerObject->ReadTextFromFile(DoseCalcsCore_build_dir_path+"/F");
    // remove the file interface
    BashCommandsForExecuting = "#! /bin/bash \ncd " + DoseCalcsCore_build_dir_path + "\n";
    BashCommandsForExecuting += "rm -r " + DoseCalcsCore_build_dir_path+"/F" + "\n";
    BashCommandsForExecuting += "\n bash \n";
    fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);
    term->executeLocalCommand(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);

    //QTextStream(stdout) << " \n\n\n\n\n -------------------------------------- " << "\n";
    //// read runing job names by process output read
    //QProcess process;
    //process.start("ls");
    ////if (!process.waitForFinished()) {
    ////    //qDebug() << "Error: Unable to execute Rocks command.";
    ////    //return nodes;
    ////}
    //QString output = process.readAllStandardOutput();
    //QTextStream(stdout) << " output " << output << "\n";

    //QStringList lines = output.split("\n", Qt::SkipEmptyParts);
    //for(int dd=0; dd < lines.size();dd++){
    //    Commlines.push_back(lines[dd]);
    //    QTextStream(stdout) << " lines[dd] " << lines[dd] << "\n";

    //}




    listOfRun.empty();listOfRun.clear();
    listOfRunIDs.empty();listOfRunIDs.clear();

    for(int aa=0; aa < Commlines.size();aa++){
        // InputsVals = Commlines[aa].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
        // listOfRun << InputsVals[2];
        // listOfRunIDs << InputsVals[0];

        InputsVals = Commlines[aa].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
        listOfRunIDs << InputsVals[0];
        aa++;
        InputsVals = Commlines[aa].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
        listOfRun << InputsVals[2];

        //showResultsOutput("listOfRun "+ InputsVals[0] , 1);
    }

    // list of terminated job
    QStringList listOfTerminated ;
    for(int aa=0; aa < listOfFile.size();aa++){
        bool isin = false;
        for(int b=0; b < listOfRun.size();b++){
            if(listOfRun[b] == listOfFile[aa]){
                isin = true; break;
            }
        }
        if(isin == false){
            listOfTerminated << listOfFile[aa];
            //showResultsOutput("listOfTerminated "+ listOfFile[aa] , 1);
        }
    }

    // list of shown jobs
    QStringList ComboboxListOfJobs ; ComboboxListOfJobs << "" << "all" << "run" ;
    for(int aa=0; aa < listOfRun.size();aa++){
        ComboboxListOfJobs << listOfRun[aa];
    }
    ComboboxListOfJobs << "terminated";
    for(int aa=0; aa < listOfTerminated.size();aa++){
        ComboboxListOfJobs << listOfTerminated[aa];
    }

    QDialog * d = new QDialog();
    QGridLayout* GraphLayout = new QGridLayout;

    QLabel* lb = new QLabel("DoseCalcs Job IDs");
    DoseCalcsJobIDsCombobox = new QComboBox(); DoseCalcsJobIDsCombobox->addItems(ComboboxListOfJobs);

    if(StopOrCheck){
        DoseCalcsJobIDsCombobox->setToolTip("Choose the ID of job to be deleted");
        d->setWindowTitle("Choose Job ID to be deleted");
    }else{
        DoseCalcsJobIDsCombobox->setToolTip("Choose the ID of job to be checked");
        d->setWindowTitle("Choose Job ID to be checked");
    }

    DoseCalcsJobIDsCombobox->setCurrentText("");
    if(QdelID != "" || !QdelID.isEmpty()){
        DoseCalcsJobIDsCombobox->setCurrentText(QdelID);
    }

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    int ii = 0, jj=0;
    GraphLayout->addWidget(lb, jj,ii,1,1);
    GraphLayout->addWidget(DoseCalcsJobIDsCombobox, jj,++ii,1,1);
    QCheckBox* cb = new QCheckBox("Remove");;
    QPushButton* pb = new QPushButton("Check");;

    if(StopOrCheck){
        cb->setToolTip("Remove the job output file corresponding to the selected run job ID. If terminated is chosen or a specific finished job, the job file will be deleted besides without checking ");
        connect(pb, SIGNAL(clicked()), this, SLOT(CheckJobOutputsFromStopSlot()));
        pb->setToolTip("Check the DoseCalcs job before deleting it.");
        GraphLayout->addWidget(pb, jj,++ii,1,1);
        GraphLayout->addWidget(cb, jj,++ii,1,1);
    }
    GraphLayout->addWidget(buttonBox);
    d->setLayout(GraphLayout);

    int result = d->exec();

    if(result == QDialog::Accepted){

        if(StopOrCheck){

            if(DoseCalcsJobIDsCombobox->currentText()!=""){

                showResultsOutput("Deleting DoseCalcs Job "+ DoseCalcsJobIDsCombobox->currentText() , 1);

                BashCommandsForExecuting = "#! /bin/bash \ncd " + DoseCalcsCore_build_dir_path + "\n";

                if(DoseCalcsJobIDsCombobox->currentText() == "all"){

                    for(int dd = 0 ; dd < listOfRun.size() ; dd++){
                        BashCommandsForExecuting += "qdel "+listOfFilesName[listOfRun[dd]] + "\n";
                        if(cb->isChecked()){
                            QString nnn = listOfRun[dd] ; nnn.remove("DoseCalcs");
                            BashCommandsForExecuting += "rm "+listOfRun[dd]+".o"+listOfFilesName[listOfRun[dd]] + "\n" +
                                    "rm Macros"+nnn+".mac\n";
                            listOfRun[dd] = "DoseCalcs" + listOfRun[dd] ;

                        }
                    }
                    if(cb->isChecked()){
                        for(int aa = 0; aa < listOfTerminated.size() ; aa++){
                            QString nnn = listOfTerminated[aa] ; nnn.remove("DoseCalcs");
                            BashCommandsForExecuting += "rm "+ listOfTerminated[aa]+".o"+listOfFilesName[listOfTerminated[aa]] + "\n" +
                                    "rm Macros"+nnn+".mac\n";
                            listOfTerminated[aa] = "DoseCalcs" + listOfTerminated[aa] ;
                        }
                    }
                }
                else if(DoseCalcsJobIDsCombobox->currentText() == "run"){
                    for(int dd = 0 ; dd < listOfRun.size() ; dd++){
                        BashCommandsForExecuting += "qdel "+listOfFilesName[listOfRun[dd]] + "\n";
                        if(cb->isChecked()){
                            QString nnn = listOfRun[dd] ; nnn.remove("DoseCalcs");
                            BashCommandsForExecuting += "rm " + listOfRun[dd] +".o"+listOfFilesName[listOfRun[dd]] + "\n" +
                                    "rm Macros"+nnn+".mac\n";
                            listOfRun[dd] = "DoseCalcs" + listOfRun[dd] ;
                        }
                    }
                }
                else if(DoseCalcsJobIDsCombobox->currentText() == "terminated"){
                    if(cb->isChecked()){
                        for(int aa = 0; aa < listOfTerminated.size() ; aa++){
                            QString nnn = listOfTerminated[aa] ; nnn.remove("DoseCalcs");
                            BashCommandsForExecuting += "rm " +listOfTerminated[aa]+ ".o"+ listOfFilesName[listOfTerminated[aa]] + "\n" +
                                    "rm Macros"+nnn+".mac\n";
                            listOfTerminated[aa] = "DoseCalcs" + listOfTerminated[aa] ;
                        }
                    }
                }
                else{

                    bool isin = false;

                    // look in the run
                    for(int b=0; b < listOfRun.size();b++){
                        if(listOfRun[b] == DoseCalcsJobIDsCombobox->currentText()){
                            isin = true; break;
                        }
                    }
                    if(isin == true){
                        BashCommandsForExecuting += "qdel "+listOfFilesName[DoseCalcsJobIDsCombobox->currentText() ]+ "\n";
                        if(cb->isChecked()){
                            BashCommandsForExecuting += "rm "+DoseCalcsJobIDsCombobox->currentText()+".o"+listOfFilesName[DoseCalcsJobIDsCombobox->currentText()] + "\n" +
                                    "rm Macros"+DoseCalcsJobIDsCombobox->currentText().remove("DoseCalcs")+".mac\n";
                        }
                    }

                    isin = false;

                    // look in the terminated
                    for(int b=0; b < listOfTerminated.size();b++){
                        if(listOfTerminated[b] == DoseCalcsJobIDsCombobox->currentText()){
                            isin = true; break;
                        }
                    }
                    if(isin == true){
                        if(cb->isChecked()){
                            BashCommandsForExecuting += "rm "+DoseCalcsJobIDsCombobox->currentText()+".o"+listOfFilesName[DoseCalcsJobIDsCombobox->currentText()] + "\n" +
                                    "rm Macros"+DoseCalcsJobIDsCombobox->currentText().remove("DoseCalcs")+".mac\n";
                        }
                    }
                }

                BashCommandsForExecuting += "\n bash \n";
                fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

                ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
            }
        }
        else{
            CheckJobOutputsFromStopSlot();
        }
    }
}
void MainWindow::CheckJobOutputsFromStopSlot(){

    if(DoseCalcsJobIDsCombobox->currentText() == "all" || DoseCalcsJobIDsCombobox->currentText() == "run" || DoseCalcsJobIDsCombobox->currentText() == "terminated" || DoseCalcsJobIDsCombobox->currentText() == ""){
        return;
    }

    QString fileToCheck = DoseCalcsCore_build_dir_path+"/"+DoseCalcsJobIDsCombobox->currentText()+".o" + listOfFilesName[DoseCalcsJobIDsCombobox->currentText()];

    showResultsOutput("Checking the : " + fileToCheck , 0);

    QdelID = DoseCalcsJobIDsCombobox->currentText() ;

    BashCommandsForExecuting = "#! /bin/bash \ncd " + DoseCalcsCore_build_dir_path +"\n";
    //+ "\n qstat -j "+QdelID

    bool isin = false;

    // look in the run
    for(int b=0; b < listOfRun.size();b++){
        if(listOfRun[b] == DoseCalcsJobIDsCombobox->currentText()){
            isin = true; break;
        }
    }
    if(isin == true){
        BashCommandsForExecuting += "\n qstat -j "+listOfFilesName[QdelID]+" | grep job_name"
                + "\n qstat -j "+listOfFilesName[QdelID]+" | grep script_file"
                + "\n qstat -j "+listOfFilesName[QdelID]+" | grep submission_time"
                + "\n qstat -j "+listOfFilesName[QdelID]+" | grep range"
                + "\n qstat -j "+listOfFilesName[QdelID]+" | grep usage";
    }

    BashCommandsForExecuting += "\n tail -f "+fileToCheck;
    BashCommandsForExecuting += "\n bash \n";

    showResultsOutput("Writing Check Commands : \n", 0);
    showResultsOutput(BashCommandsForExecuting , 0);
    showResultsOutput("to --> " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , 4);
    fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

    if(QFile::exists(fileToCheck)){

        showResultsOutput("Checking DoseCalcs output file" , 1);

        //ShowTerminal("tail -f " + fileToCheck);
        ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
    }

}
QString MainWindow::ConstructDoseCalcsJobName(){

    NumberOfSourceParticle = 1;
    NumberOfSourceEnergies = 1;
    NumberOfSourceRegions = 1;
    NumberOfSourceMomentumDirection = 1;

    QString DoseCalcsJobName = "";

    QStringList InputsVals0,InputsVals1,InputsVals2,InputsVals3;

    InputsVals0 = ui->SourcelineEditParName->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
    NumberOfSourceParticle = InputsVals0.size();

    if(ui->comboBoxTypeOfSources->currentText() == "Voxels"){
        InputsVals1 = ui->lineEditChosenSourceTypeData->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        NumberOfSourceRegions = InputsVals1.size();
    }
    else if(ui->comboBoxTypeOfSources->currentText() == "Volume"){

        QStringList InputsVals11 = ui->lineEditChosenSourceTypeData->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        if(InputsVals11.size()>0){
            int nnn = 0;
            while( nnn < InputsVals11.size() ){
                InputsVals1 << InputsVals11[nnn];
                QTextStream(stdout) << " nnn " << nnn << " -- " << InputsVals11[nnn] << "\n";

                nnn = nnn + 4;
            }
        }
        NumberOfSourceRegions = InputsVals1.size();
    }

    if(ui->SourceComboBoxEnergyDist->currentText() != "File"){
        InputsVals2 = ui->lineEditSpecialEnergyDistributionParameter->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        NumberOfSourceEnergies = InputsVals2.size();
    }

    if(ValuesOfInputs[ui->ScoreCombobowSimNumOnRanksLineEdit->currentText()] == "o"){
        if(InputsVals0.size() > 0 && InputsVals1.size() > 0 && InputsVals2.size() > 0){
            DoseCalcsJobName += InputsVals0[0] +"_"+InputsVals1[0] +"_"+InputsVals2[0];
        }
    }else{

        for(int a = 0 ; a < InputsVals0.size() ; a++){
            DoseCalcsJobName += InputsVals0[a][0]+"";
        }

        //DoseCalcsJobName += "_";
        for(int c = 0 ; c < InputsVals2.size() ; c++){
            //DoseCalcsJobName += InputsVals2[c][0];
        }

        DoseCalcsJobName += "_";
        for(int b = 0 ; b < InputsVals1.size() ; b++){
            DoseCalcsJobName += InputsVals1[b][0];
        }
    }

    Geometry_setGeometrySymbole = ui->lineEditGeometrySymbole->text();
    Execution_setEventNumber = ui->lineEditNumberOfEvent->text();
    Execution_setNumberOfRanksOrThreads = ui->lineEditNumberOfRanksOrThreads->text();

    if(ui->lineEditGeometrySymbole->text()==""){ ui->lineEditGeometrySymbole->setText("phantom0") ; Geometry_setGeometrySymbole = "phantom0";}
    if(ui->lineEditNumberOfRanksOrThreads->text()==""){ ui->lineEditNumberOfRanksOrThreads->setText("3") ; Execution_setNumberOfRanksOrThreads = "3";}
    if(ui->lineEditNumberOfEvent->text()==""){ ui->lineEditNumberOfEvent->setText("10000") ; Execution_setEventNumber = "10000";}

    DoseCalcsJobName = "DoseCalcs_"+QString::number(macrosfileinc)+"_"+ ui->lineEditGeometrySymbole->text()+"_"+DoseCalcsJobName+ "_"+QString::number(Execution_setNumberOfRanksOrThreads.toInt());
    //DoseCalcsJobName = "DoseCalcs_"+ ui->lineEditGeometrySymbole->text()+"_"+Execution_setNumberOfRanksOrThreads;

    return DoseCalcsJobName;
}

void MainWindow::GenerateMacrosFilesForBatch(){

    MacrosNamesForBatch->clear();
    ChosenEnergiesVector->clear();
    ChosenRegionVector->clear();
    ChosenParticlesVector->clear();

    QStringList InputsVals0,InputsVals1,InputsVals2,InputsVals3;

    InputsVals0 = ui->SourcelineEditParName->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"

    if(ui->comboBoxTypeOfSources->currentText() == "Voxels" || ui->comboBoxTypeOfSources->currentText() == "TET"){
        for(int qq=0; qq < InputsVals1.size();qq++){
            ChosenRegionVector->push_back(InputsVals1[qq]);
        }
    }
    else if(ui->comboBoxTypeOfSources->currentText() == "Volume"){

        QStringList InputsVals11 = ui->lineEditChosenSourceTypeData->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        if(InputsVals11.size()>0){
            int nnn = 0;
            while( nnn < InputsVals11.size() ){
                ChosenRegionVector->push_back(InputsVals11[nnn]);
                //QTextStream(stdout) << " nnn " << nnn << " -- " << InputsVals11[nnn] << "\n";

                nnn = nnn + 4;
            }
        }
    }

    if(ui->SourceComboBoxEnergyDist->currentText() != "File"){
        for(int qq=0; qq < InputsVals2.size();qq++){
            ChosenEnergiesVector->push_back(InputsVals2[qq]);
        }
    }

    for(int qq=0; qq < InputsVals0.size();qq++){
        ChosenParticlesVector->push_back(InputsVals0[qq]);
    }

    Geometry_setGeometrySymbole = ui->lineEditGeometrySymbole->text();
    Execution_setEventNumber = ui->lineEditNumberOfEvent->text();
    Execution_setNumberOfRanksOrThreads = ui->lineEditNumberOfRanksOrThreads->text();

    if(ui->lineEditGeometrySymbole->text()==""){ ui->lineEditGeometrySymbole->setText("phantom0") ; Geometry_setGeometrySymbole = "phantom0";}
    if(ui->lineEditNumberOfRanksOrThreads->text()==""){ ui->lineEditNumberOfRanksOrThreads->setText("3") ; Execution_setNumberOfRanksOrThreads = "3";}
    if(ui->lineEditNumberOfEvent->text()==""){ ui->lineEditNumberOfEvent->setText("10000") ; Execution_setEventNumber = "10000";}


    if(ValuesOfInputs[ui->ScoreCombobowSimNumOnRanksLineEdit->currentText()] == "o"){

    }else{

        for(int a=0; a < ChosenParticlesVector->size(); a++){
            for(int b=0; b < ChosenEnergiesVector->size();b++){
                for(int c=0; c < ChosenRegionVector->size();c++){

                    initializeVariable();
                    SaveDataFromInputComponents(); // get the same componenet data but the four under variable values are related to the SaveEnePharOrgLists() data that is called one time when the the run in pushed
                    CreateUserCommands(); //fill the variables with all values to save it to user file and not inputFile executed by the geant4 application core

                    QsubSeparatedMacroFilePath = DoseCalcsCore_build_dir_path+"/Macros"+ConstructDoseCalcsJobName().remove("DoseCalcs")+".mac";
                    fileManagerObject->WriteTextToFile(QsubSeparatedMacroFilePath , generateInputUserTextForinputFile());

                }
            }
        }
    }

}

void MainWindow::ShowTerminal(QString Command){

    fileManagerObject->WriteTextToFile(Command, fileManagerObject->ReadTextFromFileInOneString(Command)+"\n bash");

    terminal *term = new terminal;
    term->scroll(20, 20);
    term->executeCommand(Command);
    ui->tabWidget->removeTab(2);
    ui->tabWidget->insertTab(2, term, "Terminal");
    ui->tabWidget->setCurrentIndex(2);
}

// Return Buttons
void MainWindow::on_BtnPhantomReturn_clicked()
{
    RefillGeometryInputs();
}
void MainWindow::on_BtnsourceReturn_clicked()
{
    RefillSourcePhysicsInputs();
}
void MainWindow::on_btnanalyseReturn_clicked()
{
}
void MainWindow::on_RootAnalyseBtnReturn_clicked()
{
    RefillROOTAnalysisInputs();
}
// Reset buttons
void MainWindow::on_PhantomBtnReset_clicked()
{
    SaveGeometryInputs();
    InitialiseGeometryInputs();
}
void MainWindow::on_SourceBtnReset_clicked()
{

    SaveSourcePhysicsInputs();
    InitialiseSourcePhysicsInputs();

}
void MainWindow::on_AnalyseBtnReset_clicked()
{


}
void MainWindow::on_RootAnalyseBtnReset_clicked()
{
    SaveROOTAnalysisInputs();
    InitialiseROOTAnalysisInputs();

}

void MainWindow::SaveDataFromInputComponents(){

    showResultsOutput("Saving data from Components to the variables...", 0);

    SaveGeometryInputs();
    SaveSourcePhysicsInputs();
    SaveROOTAnalysisInputs();

}
// called from ResetButtons and SaveInput SaveDataFromInputComponents()
void MainWindow::SaveGeometryInputs(){

    Geometry_setWorldHalfSize = ui->PhantomWorldHalfSizeslineEdit->text();
    Geometry_setWorldHalfSizeunit = ui->comboBoxWorldSizeUnit->currentText();

    Geometry_setWorldMaterialName = ui->PhantomWorldMaterialLineEdit->currentText();

    if(ui->radioButtonGDML->isChecked()){
        ui->radioButtonGDML->setChecked(true);
        Geometry_CreateVolume_GeometryFileType = "GDML";
        //Geometry_CreateVolume_GeometryPath = ui->GeometryFilePathLineEdit->text();
    }
    if(ui->radioButtonTEXT->isChecked()){
        ui->radioButtonTEXT->setChecked(true);
        Geometry_CreateVolume_GeometryFileType = "TEXT";
        //Geometry_CreateVolume_GeometryPath = ui->GeometryFilePathLineEdit->text();;
    }
    else if(ui->radioButtonDICOM->isChecked()){
        ui->radioButtonDICOM->setChecked(true);
        Geometry_CreateVolume_GeometryFileType = "DICOM";
    }
    else if(ui->radioButtonVoxel->isChecked()){
        ui->radioButtonVoxel->setChecked(true);
        Geometry_CreateVolume_GeometryFileType = "VOXEL";
    }
    else if(ui->radioButtonVoxIDs->isChecked()){
        ui->radioButtonVoxIDs->setChecked(true);
        Geometry_CreateVolume_GeometryFileType = "VoxIDs";
        Geometry_CreateVolume_GeometryPath = ui->lineEditVoxIDsFilePath->text();
    }
    else if(ui->radioButtonTET->isChecked()){
        ui->radioButtonTET->setChecked(true);
        Geometry_CreateVolume_GeometryFileType = "TET";
        Geometry_CreateVolume_GeometryPath = ui->lineEditVoxIDsFilePath->text();
    }
    else if(ui->radioButtonConstruct->isChecked()){
        ui->radioButtonConstruct->setChecked(true);
        Geometry_CreateVolume_GeometryFileType = "Construct";
    }

    Geometry_setGeometrySymbole = ui->lineEditGeometrySymbole->text();

}
void MainWindow::SaveSourcePhysicsInputs(){

    SourceData_setParticleName = ui->SourcelineEditParName->text();
    PhysicsData_setPhysicsName = ui->SourceComboBoxPhysUsed->currentText();
    PhysicsData_setPhotoElectricEffectModel = ui->comboBoxPEEModels->currentText();
    PhysicsData_setComptonScatteringModel = ui->comboBoxComptonModels->currentText();
    PhysicsData_setGammaConversionModel = ui->comboBoxGammaConversionModels->currentText();
    PhysicsData_setRayleighScatteringModel = ui->comboBoxRayleighScatteringModels->currentText();
    PhysicsData_setElectronIonisationModel = ui->comboBoxElectronIonisationModels->currentText();
    PhysicsData_setElectronBremModel = ui->comboBoxElectronBremModels->currentText();
    PhysicsData_setHadronIonisationModel = ui->comboBoxHadronIonisationModels->currentText();

    PhysicsData_setCutsEnergy = ui->SourceLineEditEnergyCut->text();
    PhysicsData_setCutsDistance = ui->SourceLineEditDistanceCut->text();

    PhysicsData_ParticleForCrossSection = ui->lineEditParticleNamesForCrossSection->text() ;
    PhysicsData_EUnitForCrossSection = ui->comboBoxEnergyUnitsForCrossSection->currentText() ;
    PhysicsData_EnergiesForCrossSection = ui->lineEditEnergiesForCrossSection->text() ;

    //SourceData_setBoxHalfSize = ui->PointsLineEditCubDim->text();

    SourceData_setSourceType = ui->comboBoxTypeOfSources->currentText();
    SourceData_setSourceData = ui->lineEditChosenSourceTypeData->text();
    SourceData_setSourceSizeUnit = ui->comboBoxSizeUnit->currentText();

    SourceData_setEnergyDistribution = ui->SourceComboBoxEnergyDist->currentText();
    SourceData_setEnergyData = ui->lineEditSpecialEnergyDistributionParameter->text();
    SourceData_setSourceEnergyUnit = ui->comboBoxEnergyUnit->currentText();

    SourceData_setAngleDistribution = ui->SourceComboBoxAngleDist->currentText();
    SourceData_setMomDirData = ui->lineEditSpecialAngulatDistributionParameter->text();
    SourceData_setSourceAngleUnit = ui->comboBoxAngleUnit->currentText();

    SourceData_setEventsNumForDataGen = ui->UseDataFilesFor->currentText();

    // Run
    Execution_setEventNumber = ui->lineEditNumberOfEvent->text();
    Execution_setNumberOfRanksOrThreads = ui->lineEditNumberOfRanksOrThreads->text();
    Score_setSimNumOnRanksLineEdit = ui->ScoreCombobowSimNumOnRanksLineEdit->currentText();

    if(MPI_USE != "YES"){
        // the MPI and multithreading computation modes should not depasse the use of all physical cores - 1
        if(Execution_setNumberOfRanksOrThreads.toInt() > NumberOfCPUCores){
            //Execution_setNumberOfRanksOrThreads = QString::number(NumberOfCPUCores);
            //ui->lineEditNumberOfRanksOrThreads->setText(QString::number(NumberOfCPUCores));
            showResultsOutput("Thread Number should not exceed the number of physical cores, lefting 1 thread to other applications " , 3);
        }
    }


}
void MainWindow::SaveROOTAnalysisInputs(){

    //Score
    Score_setVolumesToScore = ui->AnalysisLineEdit_ScorOrg->text();
    Score_setVariableToScore = ui->AnalysisLineEditVarToScore->text();
    //Score_setAccuracyCalculationLevel = ui->ScoreComboBoxAcuuracyLevel->currentText();

    Score_setRadioNucleidDataLineEdit = ui->lineEditRadioNucleid->text() + " " + RadDataType[ui->comboBoxRadionuclidedataType->currentText()] + " " + ui->RadioNucleidlineEmissionDataEdit->text();
    Score_setRadioNucleidBiokineticsLineEdit = ui->lineEditRadioNucleid->text() + " " + ui->RadioNucleidAdmActivitylineEdit->text() + " " + ui->comboBoxActivityUnits->currentText() + " " + ui->comboBoxTimeUnit->currentText() + " " + ui->RadioNucleidSourceResTimelineEdit->text();

    //Score_setRadiationFactors = ui->radiationEnergyFactor->text();
    Score_setTissueFactors = ui->TissueFactorLineEdit->text();
    Score_setQuantitiesUnits = " AE " +ui->comboBoxAEUnits->currentText() + " SAF " +ui->comboBoxSAFUnits->currentText() +" AD " +ui->comboBoxADUnits->currentText() +" S " +ui->comboBoxSUnits->currentText() +" H " +ui->comboBoxHUnits->currentText() + " E " +ui->comboBoxEUnits->currentText();


    // Analysis
    Analysis_setGraphsData = ui->AnalysisComboBoxGraphData->currentText();
    Analysis_setCompareType =  ui->AnalysisComboBoxGraphsType->currentText();
    Analysis_setRefFilePath = ui->AnalysisLineEditRefFile->text();
    Analysis_setRefName = ui->AnalysisLineEdit_RefName->text();
    Analysis_setRegionVariableName = ui->AnalysisComboBoxRegionVariable->currentText();

    if(ui->checkBoxRelSDev->isChecked()){
        Analysis_GenerateRelativeSDevGraph = "yes" ;
    }
    else{
        Analysis_GenerateRelativeSDevGraph = "no" ;
    }

    //if(ui->comboBoxRelDiff->currentText() != "" ){
    //if(ui->checkBoxRelErr->isChecked()){
    //}
    Analysis_GenerateRelativeErrGraph = ui->comboBoxRelDiff->currentText() ;

    if(ui->checkBoxRegionParameter->isChecked()){
        Analysis_GenerateRegionsVariableGraph = "yes" ;
    }
    else{
        Analysis_GenerateRegionsVariableGraph = "no" ;
    }

    if(ui->checkBoxEventsDataHisto->isChecked()){
        Analysis_GenerateEventsDataHisto = "yes" ;
    }
    else{
        Analysis_GenerateEventsDataHisto = "no" ;
    }


    PositionsDataFilePath = ui->AnalysisLineEditPositionsFile->text();
    EnergiesDataFilePath = ui->AnalysisLineEditEnergiesFile->text();
    MomDirsDataFilePath = ui->AnalysisLineEditMomDirsFile->text();

    Analysis_setSliceFor2DGraph = ui->AnalysisComboBoxSliceFor2DGraph->currentText();
    Analysis_setBeamAxis = ui->AnalysisComboBoxBeamAxis->currentText();
    Analysis_setSliceID = ui->AnalysisLineEdit_SliceID->text();

    if(ui->checkBoxUseLogE->isChecked()){
        Analysis_UseLogE = "yes" ;
    }
    else{
        Analysis_UseLogE = "no" ;
    }

    if(ui->checkBoxUseLogVar->isChecked()){
        Analysis_UseLogVariable = "yes" ;
    }
    else{
        Analysis_UseLogVariable = "no" ;
    }


    if(ui->checkBoxUseGrid->isChecked()){
        Analysis_UseGridXY = "yes" ;
    }
    else{
        Analysis_UseGridXY = "no" ;
    }

    if(ui->checkBoxPrintTitle->isChecked()){
        Analysis_PrintTitle = "yes" ;
    }else{
        Analysis_PrintTitle = "no" ;
    }

    Analysis_LegendPos = ui->AnalysisLegendPoscomboBox->currentText() ;

    if(ui->checkBoxAddErrorBar->isChecked()){
        Analysis_AddErrorBar = "yes" ;
    }
    else{
        Analysis_AddErrorBar = "no" ;
    }


    Analysis_LegendWidth = ui->lineEditLegendXYWidth->text();

    Analysis_setGraphsExt = ui->AnalysisComboBoxGraphsExt->currentText();
}

void MainWindow::InitialiseGeometryInputs(){

    ui->PhantomWorldHalfSizeslineEdit->setText("");
    ui->PhantomWorldMaterialLineEdit->setCurrentText("");
    ui->comboBoxWorldSizeUnit->setCurrentIndex(1);

    ui->radioButtonGDML->setChecked(false);
    ui->radioButtonTEXT->setChecked(false);
    ui->radioButtonVoxIDs->setChecked(false);
    ui->radioButtonTET->setChecked(false);
    ui->radioButtonDICOM->setChecked(false);
    ui->radioButtonVoxel->setChecked(false);
    ui->radioButtonConstruct->setChecked(false);

    ui->lineEditVoxIDsFilePath->setText("");

    ui->lineEditGeometrySymbole->setText("");

}
void MainWindow::InitialiseSourcePhysicsInputs(){

    ui->SourcelineEditParName->setText("");
    ui->SourceComboBoxPhysUsed->setCurrentText("");
    ui->comboBoxPEEModels->setCurrentIndex(0);
    ui->comboBoxComptonModels->setCurrentIndex(0);
    ui->comboBoxGammaConversionModels->setCurrentIndex(0);
    ui->comboBoxRayleighScatteringModels->setCurrentIndex(0);
    ui->comboBoxElectronIonisationModels->setCurrentIndex(0);
    ui->comboBoxElectronBremModels->setCurrentIndex(0);
    ui->comboBoxHadronIonisationModels->setCurrentIndex(0);

    ui->SourceLineEditEnergyCut->setText("");
    ui->SourceLineEditDistanceCut->setText("");

    ui->lineEditParticleNamesForCrossSection->setText("");
    ui->comboBoxEnergyUnitsForCrossSection->setCurrentText("");
    ui->lineEditEnergiesForCrossSection->setText(0);

    ui->comboBoxTypeOfSources->setCurrentIndex(0);
    ui->lineEditChosenSourceTypeData->setText("");
    ui->comboBoxSizeUnit->setCurrentIndex(1);

    ui->SourceComboBoxEnergyDist->setCurrentIndex(0);
    ui->lineEditSpecialEnergyDistributionParameter->setText("");
    ui->comboBoxEnergyUnit->setCurrentIndex(2);

    ui->SourceComboBoxAngleDist->setCurrentIndex(0);
    ui->lineEditSpecialAngulatDistributionParameter->setText("");
    ui->comboBoxAngleUnit->setCurrentIndex(0);

    ui->UseDataFilesFor->setCurrentIndex(0);

    // Run
    //ui->LineEditDoseCalcsCore_build_dir_pathPath->setText("");
    ui->lineEditNumberOfRanksOrThreads->setText("");
    ui->lineEditNumberOfEvent->setText("");
    ui->ScoreCombobowSimNumOnRanksLineEdit->setCurrentIndex(0);

    if(MPI_USE == "YES"){
        ui->lineEditNumberOfRanksOrThreads->setPlaceholderText("Ranks Number");
        ui->lineEditNumberOfEvent->setPlaceholderText("Events Number Per Rank");
    }else{
        ui->lineEditNumberOfRanksOrThreads->setPlaceholderText("Threads Number");
        ui->lineEditNumberOfEvent->setPlaceholderText("Events Number Per Thread");
    }
}
void MainWindow::InitialiseROOTAnalysisInputs(){

    // Score
    ui->AnalysisLineEdit_ScorOrg->setText("");
    ui->AnalysisLineEditVarToScore->setText("");
    //ui->ScoreComboBoxAcuuracyLevel->setCurrentIndex(0);

    ui->lineEditRadioNucleid->setText("");
    ui->RadioNucleidlineEmissionDataEdit->setText("");
    ui->comboBoxActivityUnits->setCurrentIndex(0);
    ui->RadioNucleidSourceResTimelineEdit->setText("");
    ui->comboBoxTimeUnit->setCurrentIndex(0);

    //ui->radiationEnergyFactor->setText("");
    ui->TissueFactorLineEdit->setText("");

    // analysis
    ui->AnalysisComboBoxGraphData->setCurrentIndex(0);
    ui->AnalysisComboBoxGraphsType->setCurrentIndex(0);
    ui->AnalysisComboBoxGraphsExt->setCurrentIndex(0);
    ui->AnalysisLineEditRefFile->setText("");
    ui->AnalysisLineEdit_RefName->setText("");
    ui->AnalysisComboBoxRegionVariable->setCurrentIndex(0);

    ui->checkBoxRelSDev->setChecked(true);
    //ui->checkBoxRelErr->setChecked(true);
    ui->comboBoxRelDiff->setCurrentText("");

    ui->checkBoxRegionParameter->setChecked(true);
    ui->checkBoxEventsDataHisto->setChecked(true);

    ui->AnalysisLineEditPositionsFile->setText("");
    ui->AnalysisLineEditEnergiesFile->setText("");
    ui->AnalysisLineEditMomDirsFile->setText("");

    ui->AnalysisComboBoxSliceFor2DGraph->setCurrentText("");
    ui->AnalysisComboBoxBeamAxis->setCurrentText("");
    ui->AnalysisLineEdit_SliceID->setText("");

    ui->checkBoxUseLogE->setChecked(true);
    ui->checkBoxUseLogVar->setChecked(true);
    ui->checkBoxUseGrid->setChecked(true);
    ui->checkBoxPrintTitle->setChecked(true);
    ui->AnalysisLegendPoscomboBox->setCurrentText("") ;

    ui->lineEditLegendXYWidth->setText("");

    ui->checkBoxAddErrorBar->setChecked(true);

}

void MainWindow::RefillGeometryInputs(){

    ui->PhantomWorldHalfSizeslineEdit->setText(Geometry_setWorldHalfSize);
    ui->comboBoxWorldSizeUnit->setCurrentText(Geometry_setWorldHalfSizeunit);

    ui->PhantomWorldMaterialLineEdit->setCurrentText(Geometry_setWorldMaterialName);

    if(Geometry_CreateVolume_GeometryFileType == "GDML"){
        ui->radioButtonGDML->setChecked(true);
        //ui->GeometryFilePathLineEdit->setText(Geometry_CreateVolume_GeometryPath);
    }
    else if(Geometry_CreateVolume_GeometryFileType == "TEXT"){
        ui->radioButtonTEXT->setChecked(true);
        //ui->GeometryFilePathLineEdit->setText(Geometry_CreateVolume_GeometryPath);
    }
    else if(Geometry_CreateVolume_GeometryFileType == "VoxIDs"){
        ui->radioButtonVoxIDs->setChecked(true);
        ui->lineEditVoxIDsFilePath->setText(Geometry_CreateVolume_GeometryPath);
    }
    else if(Geometry_CreateVolume_GeometryFileType == "TET"){
        ui->radioButtonTET->setChecked(true);
        ui->lineEditVoxIDsFilePath->setText(Geometry_CreateVolume_GeometryPath);
    }
    else if(Geometry_CreateVolume_GeometryFileType == "DICOM"){
        ui->radioButtonDICOM->setChecked(true);
    }
    else if(Geometry_CreateVolume_GeometryFileType == "VOXEL"){
        ui->radioButtonVoxel->setChecked(true);
    }
    else if(Geometry_CreateVolume_GeometryFileType ==  "Construct"){
        ui->radioButtonConstruct->setChecked(true);
    }

    ui->lineEditGeometrySymbole->setText(Geometry_setGeometrySymbole);

    updateApplicationTabs();
}
void MainWindow::RefillSourcePhysicsInputs(){

    ui->SourceComboBoxPhysUsed->setCurrentText(PhysicsData_setPhysicsName);
    ui->comboBoxPEEModels->setCurrentText(PhysicsData_setPhotoElectricEffectModel);
    ui->comboBoxComptonModels->setCurrentText(PhysicsData_setComptonScatteringModel);
    ui->comboBoxGammaConversionModels->setCurrentText(PhysicsData_setGammaConversionModel);
    ui->comboBoxRayleighScatteringModels->setCurrentText(PhysicsData_setRayleighScatteringModel);
    ui->comboBoxElectronIonisationModels->setCurrentText(PhysicsData_setElectronIonisationModel);
    ui->comboBoxElectronBremModels->setCurrentText(PhysicsData_setElectronBremModel);
    ui->comboBoxHadronIonisationModels->setCurrentText(PhysicsData_setHadronIonisationModel);

    ui->SourceLineEditEnergyCut->setText(PhysicsData_setCutsEnergy);
    ui->SourceLineEditDistanceCut->setText(PhysicsData_setCutsDistance);

    ui->lineEditParticleNamesForCrossSection->setText(PhysicsData_ParticleForCrossSection);
    ui->comboBoxEnergyUnitsForCrossSection->setCurrentText(PhysicsData_EUnitForCrossSection);
    ui->lineEditEnergiesForCrossSection->setText(PhysicsData_EnergiesForCrossSection);

    ui->comboBoxTypeOfSources->setCurrentText(SourceData_setSourceType);
    ui->lineEditChosenSourceTypeData->setText(SourceData_setSourceData);
    ui->comboBoxSizeUnit->setCurrentText(SourceData_setSourceSizeUnit);

    ui->SourceComboBoxAngleDist->setCurrentText(SourceData_setAngleDistribution);
    ui->lineEditSpecialAngulatDistributionParameter->setText(SourceData_setMomDirData) ;
    ui->comboBoxAngleUnit->setCurrentText(SourceData_setSourceAngleUnit);

    ui->SourceComboBoxEnergyDist->setCurrentText(SourceData_setEnergyDistribution);
    ui->lineEditSpecialEnergyDistributionParameter->setText(SourceData_setEnergyData);
    ui->comboBoxEnergyUnit->setCurrentText(SourceData_setSourceEnergyUnit);

    ui->UseDataFilesFor->setCurrentText(SourceData_setEventsNumForDataGen);

    /*

    ui->PointsLineEditCubDim->setText(SourceData_setBoxHalfSize);

    ui->PointsOrganSourceLineEdit->setText(SourceData_setSourceInPoint);

    ui->comboBoxTypeOfSources->setCurrentText(SourceData_setSourceType);
    ui->lineEditChosenSourceTypeData->setText(SourceData_setSourcePosition);


    if(ui->SourceComboBoxAngleDist->currentText()=="Directed"){
        ui->lineEditSpecialAngulatDistributionParameter->setText(SourceData_setDirectionTheta + " " + SourceData_setDirectionPhi) ;
    }

    if(ui->SourceComboBoxEnergyDist->currentText()=="Gauss"){
        ui->lineEditSpecialEnergyDistributionParameter->setText(SourceData_setGaussEData);
    }
    else if(ui->SourceComboBoxEnergyDist->currentText()=="Rayleigh"){
        ui->lineEditSpecialEnergyDistributionParameter->setText(SourceData_setRayleighEmax);
    }
    else if(ui->SourceComboBoxEnergyDist->currentText()=="Mono"){
        ui->lineEditSpecialEnergyDistributionParameter->setText(SourceData_setMonoEnergy);
    }
    else if(ui->SourceComboBoxEnergyDist->currentText()=="Uniform"){
        ui->lineEditSpecialEnergyDistributionParameter->setText(SourceData_setUniformEData);
    }
*/

    // Run
    ui->lineEditNumberOfEvent->setText(Execution_setEventNumber);
    ui->lineEditNumberOfRanksOrThreads->setText(Execution_setNumberOfRanksOrThreads);
    ui->ScoreCombobowSimNumOnRanksLineEdit->setCurrentText(Score_setSimNumOnRanksLineEdit);

    //ui->LineEditDoseCalcsCore_build_dir_pathPath->setText(DoseCalcsCore_build_dir_path);

    updateApplicationTabs();
}
void MainWindow::RefillROOTAnalysisInputs(){

    // Score

    ui->AnalysisLineEdit_ScorOrg->setText(Score_setVolumesToScore);
    ui->AnalysisLineEditVarToScore->setText(Score_setVariableToScore);
    //ui->ScoreComboBoxAcuuracyLevel->setCurrentText(Score_setAccuracyCalculationLevel);

    //ui->radiationEnergyFactor->setText(Score_setRadiationFactors);
    ui->TissueFactorLineEdit->setText(Score_setTissueFactors);

    //ui->lineEditRadioNucleid->setText(Score_setRadioNucleidDataLineEdit);
    //ui->RadioNucleidlineEmissionDataEdit->setText(Score_setRadioNucleidDataLineEdit);

    // analysis
    ui->AnalysisComboBoxGraphData->setCurrentText(Analysis_setCompareType);
    ui->AnalysisComboBoxGraphsType->setCurrentText(Analysis_setGraphsData);
    ui->AnalysisLineEditRefFile->setText(Analysis_setRefFilePath);
    ui->AnalysisLineEdit_RefName->setText(Analysis_setRefName);
    ui->AnalysisComboBoxRegionVariable->setCurrentText(Analysis_setRegionVariableName);

    if(Analysis_GenerateRelativeSDevGraph == "yes" ){
        ui->checkBoxRelSDev->setChecked(true);
    }else {
        ui->checkBoxRelSDev->setChecked(false);
    }

    ui->comboBoxRelDiff->setCurrentText(Analysis_GenerateRelativeErrGraph);
    //if(ui->comboBoxRelDiff->currentText() != "" ){
    //if(Analysis_GenerateRelativeErrGraph == "yes" ){
    //ui->comboBoxRelDiff->setCurrentText(Analysis_GenerateRelativeErrGraph);
    //}else {
    //ui->checkBoxRelErr->setChecked(false);
    //}
    if(Analysis_GenerateRegionsVariableGraph == "yes" ){
        ui->checkBoxRegionParameter->setChecked(true);
    }else {
        ui->checkBoxRegionParameter->setChecked(false);
    }
    if(Analysis_GenerateEventsDataHisto == "yes" ){
        ui->checkBoxEventsDataHisto->setChecked(true);
    }else {
        ui->checkBoxEventsDataHisto->setChecked(false);
    }

    ui->AnalysisLineEditPositionsFile->setText(PositionsDataFilePath);
    ui->AnalysisLineEditEnergiesFile->setText(EnergiesDataFilePath);
    ui->AnalysisLineEditMomDirsFile->setText(MomDirsDataFilePath);

    ui->AnalysisLineEditPositionsFile->setText(PositionsDataFilePath);
    ui->AnalysisLineEditMomDirsFile->setText(PositionsDataFilePath);
    ui->AnalysisLineEditEnergiesFile->setText(PositionsDataFilePath);

    ui->AnalysisComboBoxSliceFor2DGraph->setCurrentText(Analysis_setSliceFor2DGraph);
    ui->AnalysisComboBoxBeamAxis->setCurrentText(Analysis_setBeamAxis);
    ui->AnalysisLineEdit_SliceID->setText(Analysis_setSliceID);

    if(Analysis_UseLogE == "yes" ){
        ui->checkBoxUseLogE->setChecked(true);
    }else {
        ui->checkBoxUseLogE->setChecked(false);
    }
    if(Analysis_UseLogVariable == "yes" ){
        ui->checkBoxUseLogVar->setChecked(true);
    }else {
        ui->checkBoxUseLogVar->setChecked(false);
    }
    if(Analysis_UseGridXY == "yes" ){
        ui->checkBoxUseGrid->setChecked(true);
    }else {
        ui->checkBoxUseGrid->setChecked(false);
    }
    if(Analysis_PrintTitle == "yes" ){
        ui->checkBoxPrintTitle->setChecked(true);
    }else {
        ui->checkBoxPrintTitle->setChecked(false);
    }

    ui->AnalysisLegendPoscomboBox->setCurrentText(Analysis_LegendPos) ;
    ui->AnalysisComboBoxGraphsExt->setCurrentText(Analysis_setGraphsExt);

    if(Analysis_AddErrorBar == "yes" ){
        ui->checkBoxAddErrorBar->setChecked(true);
    }else {
        ui->checkBoxAddErrorBar->setChecked(false);
    }

    ui->lineEditLegendXYWidth->setText(Analysis_LegendWidth) ;

    updateApplicationTabs();
}

void MainWindow::on_checkBoxUseMacroCommandFile_stateChanged(int arg1)
{
    if(arg1 == 0){
        ui->checkBoxUseMacroCommandFile->setCheckState(Qt::Unchecked);
        ui->pushButtonEditMacros->setText("View");
        ui->actionEditFile->setEnabled(false);
        ui->actionUseTextInput->setChecked(false);
    }else{
        ui->checkBoxUseMacroCommandFile->setCheckState(Qt::Checked);
        ui->pushButtonEditMacros->setText("Edit");
        ui->actionEditFile->setEnabled(true);
        ui->actionUseTextInput->setChecked(true);
        ui->Tab->setCurrentIndex(1);
    }
}
void MainWindow::on_pushButtonEditMacros_clicked()
{
    initializeVariable();
    SaveDataFromInputComponents();
    CreateUserCommands();

    if(ui->actionUseTextInput->isChecked()){

        fileManagerObject->WriteTextToFile(DoseCalcsCore_build_dir_path+"/"+MacroFileName , generateInputUserTextForinputFile());

        on_pushButtonEditGeomFile_clicked();

        EditFlag = 8;
        ui->tabWidget->setTabText(0,MacroFileName);
        ui->GeometryFileTextEdit->clear();
        showResultsOutput("", 4);
        ui->GeometryFileTextEdit->setPlainText(fileManagerObject->ReadTextFromFileInOneString(DoseCalcsCore_build_dir_path+"/"+MacroFileName));
        ui->tabWidget->setCurrentIndex(0);
    }else{
        ui->outputTextConsole->setPlainText("");
        showResultsOutput("Macros file data",1);
        ui->outputTextConsole->appendPlainText(generateInputUserTextForinputFile());
        ui->tabWidget->setCurrentIndex(1);
    }
}
void MainWindow::on_actionUseTextInput_triggered()
{
    if(!ui->actionUseTextInput->isChecked()){
        ui->checkBoxUseMacroCommandFile->setCheckState(Qt::Unchecked);
        ui->pushButtonEditMacros->setText("View");
        ui->actionEditFile->setEnabled(false);
    }else{
        ui->checkBoxUseMacroCommandFile->setCheckState(Qt::Checked);
        ui->pushButtonEditMacros->setText("Edit");
        ui->actionEditFile->setEnabled(true);
    }
}
void MainWindow::on_actionEditFile_triggered()
{
    initializeVariable();
    SaveDataFromInputComponents();
    CreateUserCommands();

    if(ui->actionUseTextInput->isChecked()){

        fileManagerObject->WriteTextToFile(DoseCalcsCore_build_dir_path+"/"+MacroFileName , generateInputUserTextForinputFile());

        on_pushButtonEditGeomFile_clicked();

        EditFlag = 8;
        ui->tabWidget->setTabText(0,MacroFileName);
        ui->GeometryFileTextEdit->clear();
        showResultsOutput("", 4);
        ui->GeometryFileTextEdit->setPlainText(fileManagerObject->ReadTextFromFileInOneString(DoseCalcsCore_build_dir_path+"/"+MacroFileName));
        ui->tabWidget->setCurrentIndex(0);
    }else{
        ui->outputTextConsole->setPlainText("");
        showResultsOutput("Macros file data",1);
        ui->outputTextConsole->appendPlainText(generateInputUserTextForinputFile());
        ui->tabWidget->setCurrentIndex(1);
    }
}

void MainWindow::on_actionStop_triggered()
{

    QTextStream(stdout) << CoreProcess.processId() << "\n";

    //QProcess p; p.start("pkill simulate"); p.waitForFinished(); // you cant kill a detached process from your application directly, it will be decoupled from your application. You could get the the PID of it and then send a SIGSTOP on Linux,

    if(CoreProcess.isOpen() /*isInexec == true*/){

        QTextStream(stdout) << isInexec << "\n";

        CoreProcess.kill();
        //disconnect(&CoreProcess, 0, 0, 0); CoreProcess.kill();

        isInexec = false;
    }

}
void MainWindow::on_actionAnalysis_triggered()
{

    plotDialogObject = new PlotDialog(this);
    /*plotDialogObject->setDataAnalysis(ui->AnalysisComboBoxGraphData->currentText() , ui->AnalysisComboBoxGraphsType->currentText(),
                                      ui->AnalysisComboBoxGraphsExt->currentText() , ui->AnalysisLineEditRefFile->text(),
                                      ui->AnalysisLineEdit_RefName->text() , ui->AnalysisLineEdit_ScorOrg->text() ,
                                      ui->AnalysisLineEditVarToScore->text(), 9, 9);
    */

    plotDialogObject->show();
}
void MainWindow::on_actionInstallations_triggered()
{
    installationDialogObj = new InstallationDialog(this);
    installationDialogObj->show();

    UserCurrentResultsDirPath = DoseCalcsCore_build_dir_path;

}
void MainWindow::on_actionVisualization_triggered()
{

    if(!TestSimulateExecutableInputsToRun()){
        return;
    }

    QDialog * d = new QDialog(); d->setWindowTitle("Visualization");
    QGridLayout* GraphLayout = new QGridLayout;

    QLabel* lb = new QLabel("Choose visualization editor");

    QStringList ComboboxListOfJobs ;
    if(ui->radioButtonCpp->isChecked() ||
            ui->radioButtonConstruct->isChecked() ||
            ui->radioButtonTEXT->isChecked() ||
            ui->radioButtonGDML->isChecked() ||
            ui->radioButtonSTL->isChecked()){
        ComboboxListOfJobs << "QT Editor";
    }
    else{
        ComboboxListOfJobs << "QT Editor" << "Voxels visualization" ;
    }

    QComboBox * VisType = new QComboBox(); VisType->addItems(ComboboxListOfJobs);
    VisType->setToolTip("Choose visualization editor");

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    int ii = 0, jj=0;
    GraphLayout->addWidget(lb, jj,ii,1,1);
    GraphLayout->addWidget(VisType, jj,++ii,1,1);
    GraphLayout->addWidget(buttonBox);

    d->setLayout(GraphLayout);
    int result = d->exec();

    if(result == QDialog::Accepted)
    {
        if(VisType->currentText()=="QT Editor"){

            if(ui->radioButtonTET->isChecked() || ui->radioButtonVoxIDs->isChecked() || ui->radioButtonDICOM->isChecked()){
                if(!TestVisualizingInputsToRun()){
                    return;
                }
            }

            //initializeVariable();
            SaveDataFromInputComponents();
            CreateUserCommands();

            if(!QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName)){
                on_BuildButton_clicked();
            }

            QString MacrosFilePath = DoseCalcsCore_build_dir_path+"/VisTestMacros.mac";
            QString commandsText="";


            // Geometry

            commandsText += MaterialsDataCommandsString + "\n\n" +
                    GeometryCommands[0] + " " + GeometryData_CreateWorld + "\n";

            if(ui->radioButtonVoxIDs->isChecked()){
                commandsText += GeometryCommands[2] + " " + Geometry_CreateVolume_GeometryFileType + " " + Geometry_CreateVolume_GeometryPath + "\n" +
                        VoxIDsCommandsString;
            }
            else if(ui->radioButtonTET->isChecked()){
                commandsText += GeometryCommands[2] + " " + Geometry_CreateVolume_GeometryFileType + " " + Geometry_CreateVolume_GeometryPath + "\n" +
                        VoxIDsCommandsString;
            }
            else if(ui->radioButtonDICOM->isChecked()){
                commandsText += GeometryCommands[2] + " " + Geometry_CreateVolume_GeometryFileType + "\n" +
                        DICOMCommandsString;
            }
            else if(ui->radioButtonVoxel->isChecked()){
                commandsText += GeometryCommands[2] + " " + Geometry_CreateVolume_GeometryFileType + "\n" +
                        VOXELCommandsString;
            }
            else {
                commandsText += CONSCommandsString;
            }

            commandsText += "\n" + NotForcedVolumesCommand + "\n";
            commandsText += "\n" +GeometryCommands[4]+ " "+ Geometry_setGeometrySymbole + "\n";

            commandsText += "\n"+

                    PhysicsCommands[0] + " " + PhysicsData_setPhysicsData + "\n" +
                    PhysicsCommands[1] + " " + PhysicsData_setCutsDistance+ "\n" +
                    PhysicsCommands[3] + " " + PhysicsData_setCutsEnergy+ "\n\n" +

                    SourceCommands[0] + " " + SourceData_setParticleName + "\n" +
                    SourceCommands[1] + " " + SourceData_setSourcePosData + "\n" +
                    SourceCommands[2] + " " + SourceData_setSourceEneData + "\n" +
                    SourceCommands[3] + " " + SourceData_setSourceMomDirData + "\n" ;

            if(ui->UseDataFilesFor->currentIndex() != 0){
                commandsText += SourceCommands[4] + " " + SourceData_UseDataFiles + "\n\n" ;
            }
            else{
                commandsText += "# " + SourceCommands[4] + " " + SourceData_UseDataFiles + "\n\n" ;
            }

            commandsText += "\n" + RunAndScoreCommands[0] + " " + Score_setVolumesToScore + "\n" +
                    RunAndScoreCommands[1] + " " + Score_setVariableToScore + "\n" +
                    RunAndScoreCommands[2] + " " + Execution_setNumberOfRanksOrThreads + "\n" +
                    //RunAndScoreCommands[3] + " " + Score_setAccuracyCalculationLevel + "\n" +
                    //RunAndScoreCommands[4] + " " + Execution_setEventNumber + "\n" +
                    RunAndScoreCommands[5] + " " + ValuesOfInputs[Score_setSimNumOnRanksLineEdit] + "\n\n" ;

            if(ui->radioButtonDICOM->isChecked() || ui->radioButtonVoxel->isChecked() || ui->radioButtonVoxIDs->isChecked()){
                if(!ui->checkBoxVoxelOrRegionLevel->isChecked()){
                    commandsText += "\n" + RunAndScoreCommands[12] + "\n";
                }
            }else{
                if(ui->comboBoxPreDefinedGeom->currentText() == "MyGeometry"){
                    commandsText += "\n" + RunAndScoreCommands[12] + "\n";
                }
            }

            commandsText += "\n" + RunAndScoreCommands[13] + " " + ui->comboBoxSimulationRunFor->currentText()+"\n";

            fileManagerObject->WriteTextToFile(MacrosFilePath , commandsText);

            ui->outputTextConsole->clear();
            showResultsOutput("Visualization commands", 1);
            ui->outputTextConsole->setPlainText(commandsText);
            ui->tabWidget->setCurrentIndex(1);

            if(MPI_USE=="YES"){
                return;
            }else{

                QString f = "";

                if(!QFile::exists(geant4_Lib_dir_path+"/geant4.sh")){
                    geant4_Lib_dir_path = ShowMessageBoxAndGetPath("Directory containing geant4.sh Not Found, Click OK to Choose the Directory");
                }

                BashCommandsForExecuting = "#! /bin/bash \n" +f
                        + "cd " +geant4_Lib_dir_path +"\n"+
                        + ". ./geant4.sh\n" +f +
                        + "cd " + DoseCalcsCore_build_dir_path + "\n"
                        + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName + " G " + MacrosFilePath
                        ;
            }
            BashCommandsForExecuting += "\n bash \n";

            showResultsOutput("Writing Visualisation Commands : \n", 0);
            showResultsOutput(BashCommandsForExecuting , 0);
            showResultsOutput("to --> " + MacrosFilePath , 0);
            fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

            if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName)){

                showResultsOutput("Visualization Run" , 1);

                ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
            }
        }
        else{

            geometriesvisualizationObj = new geometriesvisualization(this);
            geometriesvisualizationObj->show();

        }
    }

}
void MainWindow::on_actionClose_triggered()
{
    QApplication::quit();
}
void MainWindow::on_actionClear_Output_triggered()
{
    if( ui->tabWidget->currentIndex() == 2){
        BashCommandsForExecuting = "#! /bin/bash \n "
                                   "cd "+DoseCalcsCore_build_dir_path+"\n"
                                   "bash \n ";
        fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);
        ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
    }
    else if( ui->tabWidget->currentIndex() == 0){
        ui->outputTextConsole->clear();
    }
    else if( ui->tabWidget->currentIndex() == 1){
        ui->outputTextConsole->clear();
    }

    //ui->outputTextConsole->undo();
}
void MainWindow::on_actionClear_Inputs_triggered()
{
    SaveGeometryInputs();
    InitialiseGeometryInputs();
    SaveSourcePhysicsInputs();
    InitialiseSourcePhysicsInputs();
    SaveROOTAnalysisInputs();
    InitialiseROOTAnalysisInputs();
}
void MainWindow::on_actionReturn_triggered()
{
    RefillGeometryInputs();
    RefillSourcePhysicsInputs();
    RefillROOTAnalysisInputs();
}
void MainWindow::on_actionOpen_triggered()
{
    QString WorkDirectory = QFileInfo(OpenMacrosFilePath).absoluteDir().absolutePath();
    if(QFile(WorkDirectory).exists()){

    }else{
        WorkDirectory = DoseCalcsCore_build_dir_path;
    }
    OpenMacrosFilePath = QFileDialog::getOpenFileName(
                this,
                tr("Open a macros file to fill GUI components"),
                WorkDirectory, //pathBuildApp,
                "All files (*.*);;Text files (*.txt)"
                );

    if(OpenMacrosFilePath != NULL && OpenMacrosFilePath != "" ){

        showResultsOutput("Chosen file is " + OpenMacrosFilePath , 0);
        FillComponentsFromInputsFile(OpenMacrosFilePath);
        MacroFilePath = OpenMacrosFilePath;
    }else{
        //QMessageBox::information(this, tr("File..!!"),"No file is choosen");
        showResultsOutput("No file is chosen, please choose an input file to fill the components automatically ! ", 3);
    }
}
void MainWindow::on_actionSave_triggered()
{
    if(!TestSimulateExecutableInputsToRun()){
        return;
    }

    initializeVariable();
    SaveDataFromInputComponents();
    CreateUserCommands();


    QString WorkDirectory = QFileInfo(OpenMacrosFilePath).absoluteDir().absolutePath();
    if(QFile(WorkDirectory).exists()){

    }else{
        WorkDirectory = DoseCalcsCore_build_dir_path;
    }
    SaveMacrosFilePath = QFileDialog::getOpenFileName(
                this,
                tr("Choose a file to save macros from GUI components"),
                WorkDirectory, //pathBuildApp,
                "All files (*.*);;Text files (*.txt)"
                );

    if(SaveMacrosFilePath == NULL && SaveMacrosFilePath == "" ){
        SaveMacrosFilePath = OpenMacrosFilePath;
    }
    if(QFile(SaveMacrosFilePath).exists()){
        SaveMacrosFilePath = DoseCalcsCore_build_dir_path+"/"+MacroFileName;
        MacroFilePath = DoseCalcsCore_build_dir_path+"/"+MacroFileName;
    }

    if(SaveMacrosFilePath != NULL && SaveMacrosFilePath != "" ){

        showResultsOutput("Chosen file is " + SaveMacrosFilePath , 0);
        QMessageBox::information(this, tr(""), "You can edit the macro commands in File editor \" " + SaveMacrosFilePath + " \" and then click Save/Close button.");

        on_pushButtonEditGeomFile_clicked();
        ui->tabWidget->setTabText(0,QFileInfo(SaveMacrosFilePath).fileName());
        ui->GeometryFileTextEdit->clear();
        ui->GeometryFileTextEdit->setPlainText(generateInputUserTextForinputFile());
        ui->tabWidget->setCurrentIndex(0);

        EditFlag = 10;
    }

}
void MainWindow::on_actionRead_File_triggered()
{

    QString WorkDirectory = QFileInfo(AnyOpenedFilePath).absoluteDir().absolutePath();
    if(QFile(WorkDirectory).exists()){

    }else{
        WorkDirectory = DoseCalcsCore_build_dir_path;
    }

    AnyOpenedFilePath = QFileDialog::getOpenFileName( this, tr("Choose a file to edit in File Editor"), WorkDirectory, "All files (*.*);;Text files (*.txt)" );

    if(AnyOpenedFilePath != NULL){

        showResultsOutput("Chosen file is " + AnyOpenedFilePath , 0);

        QMessageBox::information(this, tr(""), "You can edit the file in File editor \" " + AnyOpenedFilePath + " \" and then click Save/Close button.");

        on_pushButtonEditGeomFile_clicked();
        ui->tabWidget->setTabText(0,QFileInfo(AnyOpenedFilePath).fileName());
        ui->GeometryFileTextEdit->clear();
        ui->GeometryFileTextEdit->setPlainText(fileManagerObject->ReadTextFromFileInOneString(AnyOpenedFilePath));
        ui->tabWidget->setCurrentIndex(0);

        EditFlag = 11;
    }

}
void MainWindow::on_actionCheck_Inputs_triggered()
{

    if(ui->checkBoxSimulateGeometriesWitOneSource->isChecked()){
        bool matgeomedatagood = true;
        for(int qq=0; qq < MacrosFiles.size();qq++){
            MacrosFilePathReadForMultiGeometryAndOneSource = MacrosFiles[qq];
            MacrosCommandsString = CreateMaterialAndGeometryDataFromMacrosFile(MacrosFiles[qq]);
            if(ui->checkBoxSimulateGeometriesWitOneSource->isChecked()){
                matgeomedatagood = CheckInputsForGeometryFiles();
            }
        }
        if(matgeomedatagood){
            QMessageBox::information(this, tr(""), "The first inputs check \"GOOD\"");
        }
    }else{
        if(TestSimulateExecutableInputsToRun()){
            QMessageBox::information(this, tr(""), "The first inputs check \"GOOD\"");
        }
    }
}
void MainWindow::on_actionRun_2_triggered()
{
    on_RunButton_clicked();
}
void MainWindow::on_actionGeometryModelling_triggered()
{
    GeometryModellingEditorObject = new GeometryModellingEditor(this);
    GeometryModellingEditorObject->show();
}
void MainWindow::on_actionAbout_triggered()
{


    //UsePackagesMethods();
    QMessageBox msgBox;
    msgBox.setWindowTitle("About");
    msgBox.setWindowIcon(QIcon(QDir(QCoreApplication::applicationDirPath()).filePath(GUIPackagesAndFilesDirName+"/AppIcon.png")));
    msgBox.setText("*************************************************************\n"
                   "*** DoseCalcs v1.0 is a Geant4-based code with a GUI interface, for dosimetry Monte Carlo simulations and calculations.\n"
                   "*** DoseCalcs for internal dosimetry calculations is developed by Tarik El Ghalbzouri in a thesis directed by Prof Tarek El Bardouni, and co-directed by Prof Jaafar El Bakkali ,at ERSN, faculty of sicences, university Abdelmalek Essaadi, Tetouan, Morocco.  \n"
                   "*** Application source code can be Downloaded from https://github.com/TarikEl/DoseCalcs-gui \n"
                   "*** Documentation can be found on https://dosecalcs-gui.readthedocs.io/en/latest/ \n"
                   "*** Developer contact email: imttarikk@gmail.com\n"
                   "*************************************************************\n");
    if(msgBox.exec() == QDialog::Accepted){}
}
void MainWindow::on_actionHow_To_Use_triggered()
{
    /*
    QMessageBox msgBox;
    msgBox.setWindowTitle("Documentation");
    msgBox.setWindowIcon(QIcon(QDir(QCoreApplication::applicationDirPath()).filePath(GUIPackagesAndFilesDirName+"/AppIcon.png")));
    msgBox.setText("*************************************************************\n"
"*** Documentation can be found on https://dosecalcs-gui.readthedocs.io/en/latest/ \n"
"*************************************************************\n");
    msgBox.exec();
*/
    QMessageBox::about(0, "How to use !", " For documentation, visite the link <a href='https://dosecalcs-gui.readthedocs.io/en/latest/'>DoseCalcs</a>");

}
void MainWindow::on_actionUpdate_triggered()
{

    if(QMessageBox::Yes == QMessageBox::question(this, tr("Update"), "Automatic download and installation of the latest version of DoseCalcs on the terminal. Modify the download and installation commands if necessary and click Save/Close button. Restart DoseCalcs-Gui after finishing.")){

        BashCommandsForExecuting = "cd "+QDir::currentPath()+"/..\n"+
                "rm -r main.tar.gz \n"+
                "wget "+ DoseCalcsDownloadURL +"\n"+
                "tar xvf main.tar.gz -C " + QDir::currentPath()+"/..\n"+
                "cd "+ QDir::currentPath() + "\n"+
                "rm -r Makefile gui DoseCalcs \n"+
                "qmake "+QDir::currentPath()+"/../DoseCalcs-Gui-main/DoseCalcs.pro \n"+
                "make -j" + QString::number(NumberOfCPUCores) + "\n"
                "bash \n"
                ;
        EditFlag = 12;
        ui->tabWidget->setTabText(0,"DoseCalcs Updating commands");
        ui->GeometryFileTextEdit->clear();
        showResultsOutput("", 4);
        ui->GeometryFileTextEdit->setPlainText(BashCommandsForExecuting);
        ui->tabWidget->setCurrentIndex(0);
    }
}
void MainWindow::on_actionRestart_triggered()
{
    // Start a new process to launch the updated application
    QProcess::startDetached(QApplication::applicationFilePath());

    // Close the current instance of the application
    QApplication::quit();
}
void MainWindow::on_actionPackages_Statut_triggered()
{
    QMessageBox::information(this, tr(""), TestPackagesPathsBeforToRun());
}
void MainWindow::on_actionSave_Inputs_To_Default_File_triggered()
{
    if(!TestSimulateExecutableInputsToRun()){
        return;
    }

    initializeVariable();
    SaveDataFromInputComponents();
    CreateUserCommands();

    fileManagerObject->WriteTextToFile( DEFAULT_INPUTS, generateInputUserTextForinputFile());

}
void MainWindow::on_actionSend_Results_triggered()
{
    if(!QFile::exists("/usr/bin/scp")){

        if(QMessageBox::Yes == QMessageBox::question(this, tr("scp not found")," Install \"scp\" by typing the command below in your terminal: \n\"openssh-client openssh-server\" on ubuntu, \n or \n\"sudo yum install -y openssh-clients openssh\". Install it on xterm ?")){

            if(OSNameVersion.toLower().contains("centos")){
                BashCommandsForExecuting = "#! /bin/bash \n"
                                           "sudo yum install -y openssh-clients openssh\n"
                                           "\n bash \n"
                        ;
            }
            else if(OSNameVersion.toLower().contains("ubuntu")){
                BashCommandsForExecuting = "#! /bin/bash \n"
                                           "sudo apt-get install -y openssh-client openssh-server\n"
                                           "\n bash \n"
                        ;
            }

            showResultsOutput("Writing Run Commands : \n", 0);
            showResultsOutput(BashCommandsForExecuting , 0);
            showResultsOutput("to --> " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , 4);

            fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

            if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName)){
                ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
            }
        }else{
            return;
        }
    }
    QDialog * d = new QDialog(); d->setWindowTitle("Send results using scp command");

    QGridLayout* GraphLayout = new QGridLayout;

    QPushButton* chooseFromDirectory = new QPushButton("."); connect(chooseFromDirectory, SIGNAL(clicked()), this, SLOT(btnChooseADir1_slot()));
    FromPath = new QLineEdit(Path1); FromPath->setToolTip("user@10.13.1.17:/home/user/DoseCalcs or select directory");

    QLabel* fromtolbl= new QLabel("from ---> to");

    QPushButton* chooseToDirectory = new QPushButton("."); connect(chooseToDirectory, SIGNAL(clicked()), this, SLOT(btnChooseADir2_slot()));
    ToPath = new QLineEdit(Path2); FromPath->setToolTip("user@10.13.1.17:/home/user/DoseCalcs or select directory");

    int ii = 0, jj=0;

    GraphLayout->addWidget(chooseFromDirectory, jj,ii,1,1);
    GraphLayout->addWidget(FromPath, jj,++ii,1,1);
    GraphLayout->addWidget(fromtolbl, jj,++ii,1,1);
    GraphLayout->addWidget(ToPath, jj,++ii,1,1);
    GraphLayout->addWidget(chooseToDirectory, jj,++ii,1,1);

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    GraphLayout->addWidget(buttonBox);

    d->setLayout(GraphLayout);

    int result = d->exec();

    if(result == QDialog::Accepted)
    {

        Path1 = FromPath->text();
        Path2 = ToPath->text();

        BashCommandsForExecuting = "#! /bin/bash \n"
                                   "scp -r " +FromPath->text() +" " + ToPath->text() + "\n"
                ;

        BashCommandsForExecuting += "\n bash \n";

        showResultsOutput("Writing Run Commands : \n", 0);
        showResultsOutput(BashCommandsForExecuting , 0);
        showResultsOutput("to --> " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , 4);

        fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

        if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName)){
            ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
        }
        else{
            showResultsOutput("Cannot find file containing execution commands "+ DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName + " , you should build DoseCalcs with ROOT Analysis option" , 3);
        }
    }

}
void MainWindow::btnChooseADir1_slot(){
    QString chosen_DirPath = QFileDialog::getExistingDirectory(0, ("Choose Directory"), UserCurrentResultsDirPath) ;
    if(chosen_DirPath.isEmpty() || chosen_DirPath.isNull()){

    }else{
        Path1 = chosen_DirPath;
        FromPath->setText(Path1);
    }
}
void MainWindow::btnChooseADir2_slot(){
    QString chosen_DirPath = QFileDialog::getExistingDirectory(0, ("Choose Directory"), UserCurrentResultsDirPath) ;
    if(chosen_DirPath.isEmpty() || chosen_DirPath.isNull()){

    }else{
        Path2 = chosen_DirPath;
        ToPath->setText(Path2);

    }
}

int MainWindow::FillComponentsFromInputsFile(QString FilePathString){

    QFile* filee = new QFile(FilePathString);
    if(!filee->exists()){
        showResultsOutput("Cannot find the default inputs file : " + FilePathString , 3);
        return 0;
    }

    QStringList InputsVals;
    QString sss ;

    if(ui->checkBoxFixGeometryCommands->isChecked()){
        MaterialsDataCommandsString = "";
        GeometryCommandsString = "";
        VOXELCommandsString = "";
        CONSCommandsString = "";

        // fill the vectors of dynamique components
        ElementsNames.clear();
        MaterialsNames.clear(); MaterialsNameIDs.clear();
        SolidsNames.clear();
        VolsNames.clear();
        DefinedRegionsNames.clear();
        MaterialRegionsNames.clear();
    }
    if(ui->checkBoxFixPhysicsCommands->isChecked()){
        PhysicsCommandsString = "";
    }
    if(ui->checkBoxFixSourceCommands->isChecked()){
        SourceCommandsString = "";
    }
    if(ui->checkBoxFixRunCommands->isChecked() && ui->checkBoxFixScoreCommands->isChecked()){
        RunAndScoreCommandsString = "";
    }
    if(ui->checkBoxFixGraphsCommands->isChecked()){
        AnalysisCommandsString = "";
    }

    MacrosCommandsString = "";

    QVector< QPair<QString,QString>> Commlines = fileManagerObject->ReadTextFromFileInQStringList(FilePathString);
    QMap <QString,QString> lines = fileManagerObject->ReadLinesFromFileWithFirstWordIndicator(FilePathString);

    // to fixe the paths when using predefined geometries from PackagesAndFiles directory
    if(ui->checkBoxUsePreDefinedGeom->isChecked() && ui->checkBoxFixGeometryCommands->isChecked()){
        QString s1 = "../"+GUIPackagesAndFilesDirName;
        for(int dd=0; dd < Commlines.size();dd++){
            if(Commlines[dd].second.contains(s1)){ // /GeometryData/createVolume
                Commlines[dd].second = Commlines[dd].second.replace(s1,GUIPackagesAndFilesDirPath);
            }
        }
        for ( auto Abeg = lines.begin(); Abeg != lines.end(); ++Abeg  ){

            if(Abeg.value().contains(s1)){ // "/AnalysisData/generateRelativeErrGraph"
                lines[Abeg.key()] = lines[Abeg.key()].replace(s1,GUIPackagesAndFilesDirPath);
            }
        }
    }

    for(int dd=0; dd < Commlines.size();dd++){
        QStringList InputsVals;
        if(Commlines[dd].first == MaterialCommands[0]){ // /MaterialData/createElement
            InputsVals = Commlines[dd].second.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
            if(InputsVals.size() > 1){
                ElementsNames.push_back(InputsVals[2]);
            }
            continue;
        }
        if(Commlines[dd].first == MaterialCommands[1] || Commlines[dd].first == MaterialCommands[4] ){ // /MaterialData/createMaterial or /MaterialData/createNISTMaterial
            InputsVals = Commlines[dd].second.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
            if(InputsVals.size() > 0){
                MaterialsNames.push_back(InputsVals[0]);
                ElementsNames.push_back(InputsVals[0]);
                MaterialRegionsNames.push_back(InputsVals[0]);

            }
            continue;
        }
        if(Commlines[dd].first == VOXELCommands[4] || Commlines[dd].first == VOXELCommands[12] ){ // /MaterialData/createMaterial or /MaterialData/createNISTMaterial
            InputsVals = Commlines[dd].second.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
            if(InputsVals.size() > 0){
                DefinedRegionsNames.push_back(InputsVals[0]);
            }
            continue;
        }
        if(Commlines[dd].first == GeometryCommands[0] || Commlines[dd].first == GeometryCommands[2]){ // /GeometryData/createWorld or /GeometryData/createVolume
            //showResultsOutput(GeometryCommands[0] + " "+ Commlines[dd].second, 4)
            InputsVals = Commlines[dd].second.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);

            //QMessageBox::information(this, tr(""),"Commlines[dd].second: "+ Commlines[dd].second+ " InputsVals[0]: " +InputsVals[0]);

            if(InputsVals.size() > 0){
                if(InputsVals[0] == "GDML" || InputsVals[0] == "TEXT" || InputsVals[0] == "VoxIDs" || InputsVals[0] == "VOXEL" || InputsVals[0] == "DICOM"){
                    continue;
                }
                QString fe =  QString::fromLocal8Bit(getFileExt(InputsVals[0].toStdString()).c_str());
                QString fn =  QString::fromLocal8Bit(getFileNameFromPath(InputsVals[0].toStdString()).c_str());

                if(fe == "gdml" || fe == "geom" || fe == "c++" || fe == "cpp" || fe == "cc" || fe == "stl" || fe == "ast"){
                    if(fn == "World" || Commlines[dd].first == GeometryCommands[0]){ // for World components data
                        ui->checkBoxWorldConst->setChecked(false);
                        ui->checkBoxWorldConst->setCheckState(Qt::Unchecked);
                        on_checkBoxWorldConst_clicked(false);
                        ui->PhantomWorldHalfSizeslineEdit->setText(InputsVals[0]);
                        VolsNames.push_back("World");
                    }else{
                        VolsNames.push_back(QString::fromLocal8Bit(getFileNameFromPath(InputsVals[0].toStdString()).c_str()));
                    }
                }
                else{
                    if(Commlines[dd].first == GeometryCommands[0]){
                        if(InputsVals[0] == "MyGeometry"){
                            ui->checkBoxWorldConst->setChecked(false);
                            ui->checkBoxWorldConst->setCheckState(Qt::Unchecked);
                            on_checkBoxWorldConst_clicked(false);
                            ui->PhantomWorldHalfSizeslineEdit->setText(InputsVals[0]);
                            ui->radioButtonCpp->setChecked(true);
                            Geometry_CreateVolume_GeometryFileType = "CPP";
                        }
                        VolsNames.push_back("World");
                    }
                }
            }


            continue;
        }
        if(Commlines[dd].first == GeometryCommands[1]){ // /GeometryData/createSolid
            InputsVals = Commlines[dd].second.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
            SolidsNames.push_back(InputsVals[1]);
            continue;
        }
    }

    showResultsOutput("---ElementsNames : ", 4);for(int dd=0; dd < ElementsNames.size();dd++){ showResultsOutput(ElementsNames[dd], 4);}
    showResultsOutput("---MaterialsNames : ", 4);for(int dd=0; dd < MaterialsNames.size();dd++){ showResultsOutput(MaterialsNames[dd], 4);}
    showResultsOutput("---SolidsNames : ", 4);for(int dd=0; dd < SolidsNames.size();dd++){ showResultsOutput(SolidsNames[dd], 4);}
    showResultsOutput("---VolumesNames : ", 4);for(int dd=0; dd < VolsNames.size();dd++){ showResultsOutput(VolsNames[dd], 4);}
    showResultsOutput("---DefinedRegionsNames : ", 4);for(int dd=0; dd < DefinedRegionsNames.size();dd++){ showResultsOutput(DefinedRegionsNames[dd], 4);}
    showResultsOutput("---MaterialRegionsNames : ", 4);for(int dd=0; dd < MaterialRegionsNames.size();dd++){ showResultsOutput(MaterialRegionsNames[dd], 4);}

    for(int dd=0; dd < Commlines.size();dd++){

        if(ui->checkBoxFixGeometryCommands->isChecked()){

            // remove world volume creation(in gdml, text, c++ ) and VoxIDs, VOXEL, DICOM volume indicator(/createVolume) from geometry DataCommandString
            if(Commlines[dd].first == GeometryCommands[2]){ // /GeometryData/createVolume
                InputsVals = Commlines[dd].second.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
                if(InputsVals.size() < 1){continue;}
                QString fn =  QString::fromLocal8Bit(getFileNameFromPath(InputsVals[0].toStdString()).c_str());
                if(fn == "World" || InputsVals[0] == "VoxIDs" || InputsVals[0] == "VOXEL" || InputsVals[0] == "DICOM"){ // for World components data
                    continue;
                }
            }

            if(Commlines[dd].first == GeometryCommands[0]){ // /GeometryData/createWorld or /GeometryData/createWorld
                continue;
            }

            for(int cc=0; cc < MaterialCommands.size();cc++){
                if(Commlines[dd].first == MaterialCommands[cc]){
                    MaterialsDataCommandsString += Commlines[dd].first + " " + Commlines[dd].second +"\n" ;
                    continue;
                }
            }
            for(int cc=0; cc < GeometryCommands.size();cc++){
                if(Commlines[dd].first == GeometryCommands[cc]){
                    GeometryCommandsString += Commlines[dd].first + " " + Commlines[dd].second +"\n" ;
                    continue;
                }
            }
            for(int cc=0; cc < VOXELCommands.size();cc++){
                if(Commlines[dd].first == VOXELCommands[cc]){
                    VOXELCommandsString += Commlines[dd].first + " " + Commlines[dd].second +"\n" ;
                    continue;
                }
            }
            for(int cc=0; cc < CONSCommands.size();cc++){
                if(Commlines[dd].first == CONSCommands[cc]){
                    CONSCommandsString += Commlines[dd].first + " " + Commlines[dd].second +"\n" ;
                    continue;
                }
            }
        }
        if(ui->checkBoxFixPhysicsCommands->isChecked()){
            for(int cc=0; cc < PhysicsCommands.size();cc++){
                if(Commlines[dd].first == PhysicsCommands[cc]){
                    PhysicsCommandsString += Commlines[dd].first + " " + Commlines[dd].second +"\n" ;
                    continue;
                }
            }
        }
        if(ui->checkBoxFixSourceCommands->isChecked()){
            for(int cc=0; cc < SourceCommands.size();cc++){
                if(Commlines[dd].first == SourceCommands[cc]){
                    SourceCommandsString += Commlines[dd].first + " " + Commlines[dd].second +"\n" ;
                    continue;
                }
            }
        }
        if(ui->checkBoxFixRunCommands->isChecked() || ui->checkBoxFixScoreCommands->isChecked()){
            for(int cc=0; cc < RunAndScoreCommands.size();cc++){
                if(Commlines[dd].first == RunAndScoreCommands[cc]){
                    RunAndScoreCommandsString += Commlines[dd].first + " " + Commlines[dd].second +"\n" ;
                    continue;
                }
            }
        }
        if(ui->checkBoxFixGraphsCommands->isChecked()){
            for(int cc=0; cc < AnalysisCommands.size();cc++){
                if(Commlines[dd].first == AnalysisCommands[cc]){
                    AnalysisCommandsString += Commlines[dd].first + " " + Commlines[dd].second +"\n" ;
                    continue;
                }
            }
        }
    }

    showResultsOutput("================= Materials Data Commands ================= \n\n"+MaterialsDataCommandsString, 4);
    showResultsOutput("================= Geometry Data Commands ================= \n\n"+GeometryCommandsString, 4);
    showResultsOutput("================= Construct Geometry Data Commands ================= \n\n"+CONSCommandsString, 4);
    showResultsOutput("================= Voxelized Data Commands ================= \n\n"+VOXELCommandsString, 4);
    showResultsOutput("================= Physics Data Commands ================= \n\n"+PhysicsCommandsString, 4);
    showResultsOutput("================= Source Data Commands ================= \n\n"+SourceCommandsString, 4);
    showResultsOutput("================= Run And Score Data Commands ================= \n\n"+RunAndScoreCommandsString, 4);
    showResultsOutput("================= Analysis Data Commands =================\n\n"+AnalysisCommandsString, 4);

    if(ui->checkBoxFixGeometryCommands->isChecked()){

        ui->PhantomWorldMaterialLineEdit->clear();for(int dd=0; dd < MaterialsNames.size();dd++){ ui->PhantomWorldMaterialLineEdit->addItem(MaterialsNames[dd]);}
        ui->lineEditVoxIDsFilePath->setText("");

        InputsVals = lines[GeometryCommands[2]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
        if(InputsVals.size() > 0){
            Geometry_CreateVolume_GeometryFileType = InputsVals[0]; // "/GeometryData/createVolume first parameter"
        }

        if(Geometry_CreateVolume_GeometryFileType == "GDML"){

            //showResultsOutput("\n\n\n\n\n\n\n\n\n\n 333333333333333 \n\n\n\n\n\n\n\n\n\n", 4);

            ui->radioButtonGDML->setChecked(true);
            on_radioButtonGDML_clicked(true);
            //ui->GeometryFilePathLineEdit->setText(InputsVals[1]); // "/GeometryData/createVolume second parameter"
        }
        else if(Geometry_CreateVolume_GeometryFileType == "TEXT"){
            ui->radioButtonTEXT->setChecked(true);
            on_radioButtonTEXT_clicked(true);
            //ui->GeometryFilePathLineEdit->setText(InputsVals[1]); // "/GeometryData/createVolume second parameter"
        }
        else if(Geometry_CreateVolume_GeometryFileType == "CPP" || Geometry_CreateVolume_GeometryFileType == "C++"){

            //showResultsOutput("\n\n\n\n\n\n\n\n\n\n 333333333333333 \n\n\n\n\n\n\n\n\n\n", 4);

            ui->radioButtonCpp->setChecked(true);
            on_radioButtonCpp_clicked(true);
            //ui->GeometryFilePathLineEdit->setText(InputsVals[1]); // "/GeometryData/createVolume second parameter"
        }
        else {

            if(Geometry_CreateVolume_GeometryFileType == "VOXEL"){
                ui->radioButtonVoxel->setChecked(true);
                on_radioButtonVoxel_clicked(true);
            }
            else if(Geometry_CreateVolume_GeometryFileType == "DICOM"){
                ui->radioButtonDICOM->setChecked(true);
                on_radioButtonDICOM_clicked(true);
                InputsVals = lines[VOXELCommands[7]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/setDcmDataDir"

                //ui->GeometryFilePathLineEdit->setText(lines[InputsVals[1]]); // second parameter of "/GeometryData/setDcmDataDir"
                DICOMCommandsString = VOXELCommandsString;
            }
            else if(Geometry_CreateVolume_GeometryFileType == "VoxIDs"){
                InputsVals = lines[GeometryCommands[2]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
                ui->radioButtonVoxIDs->setChecked(true);
                on_radioButtonVoxIDs_clicked(true);
                if(InputsVals.size()>1){
                    Geometry_CreateVolume_GeometryPath = InputsVals[1];
                    ui->lineEditVoxIDsFilePath->setText(Geometry_CreateVolume_GeometryPath);
                }
                //ui->GeometryFilePathLineEdit->setText(InputsVals[1]); // "/GeometryData/createVolume second parameter"
                VoxIDsCommandsString = VOXELCommandsString;
            }
            else if(Geometry_CreateVolume_GeometryFileType == "TET"){
                InputsVals = lines[GeometryCommands[2]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
                ui->radioButtonTET->setChecked(true);
                on_radioButtonVoxIDs_clicked(true);
                if(InputsVals.size()>1){
                    Geometry_CreateVolume_GeometryPath = InputsVals[1] + " " + InputsVals[2];
                    ui->lineEditVoxIDsFilePath->setText(Geometry_CreateVolume_GeometryPath);
                }
                //ui->GeometryFilePathLineEdit->setText(InputsVals[1]); // "/GeometryData/createVolume second parameter"
                VoxIDsCommandsString = VOXELCommandsString;
            }
            else {
                ui->radioButtonConstruct->setChecked(true);
                on_radioButtonConstruct_clicked(true);
            }
        }

        InputsVals = lines[GeometryCommands[0]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
        if(InputsVals.size() > 0){
            if(InputsVals.size() == 5){
                on_checkBoxWorldConst_clicked(true);
                ui->PhantomWorldMaterialLineEdit->setCurrentText(InputsVals[0]);
                ui->PhantomWorldHalfSizeslineEdit->setText(InputsVals[1] + " " + InputsVals[2] + " " + InputsVals[3] );
                ui->comboBoxWorldSizeUnit->setCurrentText(InputsVals[4]);
            }else{

                QString fe =  QString::fromLocal8Bit(getFileExt(InputsVals[0].toStdString()).c_str());

                if(fe == "c++" || fe == "cpp" || fe == "gdml" || fe == "geom"){
                    ui->radioButtonConstruct->setChecked(true);
                }
                on_checkBoxWorldConst_clicked(false);
                ui->PhantomWorldHalfSizeslineEdit->setText(InputsVals[0]);
            }
        }

        ui->lineEditGeometrySymbole->setText("Phantom");
        if(!lines[GeometryCommands[4]].isEmpty() || lines[GeometryCommands[4]] != ""){ ui->lineEditGeometrySymbole->setText(lines[GeometryCommands[4]]);        }
    }
    if(ui->checkBoxFixPhysicsCommands->isChecked()){

        showResultsOutput("Filling Physics data...", 0);
        // Physics

        InputsVals = lines[PhysicsCommands[0]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setPhysicsData"
        if(InputsVals.size() == 1){ui->SourceComboBoxPhysUsed->setCurrentText(InputsVals[0]);}
        if(InputsVals.size() == 8){
            ui->comboBoxPEEModels->setCurrentIndex(InputsVals[1].toInt());
            ui->comboBoxComptonModels->setCurrentIndex(InputsVals[2].toInt());
            ui->comboBoxGammaConversionModels->setCurrentIndex(InputsVals[3].toInt());
            ui->comboBoxRayleighScatteringModels->setCurrentIndex(InputsVals[4].toInt());
            ui->comboBoxElectronIonisationModels->setCurrentIndex(InputsVals[5].toInt());
            ui->comboBoxElectronBremModels->setCurrentIndex(InputsVals[6].toInt());
            ui->comboBoxHadronIonisationModels->setCurrentIndex(InputsVals[7].toInt());
        }

        ui->SourceLineEditDistanceCut->setText(lines[PhysicsCommands[1]]);
        ui->SourceLineEditEnergyCut->setText(lines[PhysicsCommands[3]]);

        InputsVals = lines[PhysicsCommands[2]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setEventsInitialEneData"
        if(InputsVals.size() >= 2){
            ui->lineEditParticleNamesForCrossSection->setText(InputsVals[0]);
            ui->comboBoxEnergyUnitsForCrossSection->setCurrentText(InputsVals[1]);
            sss = ""; for(int cc = 2; cc < InputsVals.size();cc++){ sss += InputsVals[cc] + " " ; }
            ui->lineEditEnergiesForCrossSection->setText(sss);
        }
    }
    if(ui->checkBoxFixSourceCommands->isChecked()){
        // source
        showResultsOutput("Filling Source data...", 0);

        InputsVals = lines[SourceCommands[0]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setEventsInitialPosData"
        if(InputsVals.size() > 0){
            ui->SourcelineEditParName->setText(lines[SourceCommands[0]]); // "/SourceData/setEventsParticleNameData"
        }
        InputsVals = lines[SourceCommands[1]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setEventsInitialPosData"
        if(InputsVals.size() > 1){
            ui->comboBoxTypeOfSources->setCurrentText(InputsVals[1]);
            ui->comboBoxSizeUnit->setCurrentText(InputsVals[0]);
            sss = ""; for(int cc = 2; cc < InputsVals.size();cc++){ sss += InputsVals[cc] + " " ; }
            ui->lineEditChosenSourceTypeData->setText(sss);
        }

        QStringList args = getArgumentOfACommandFromText(VOXELCommands[11], 5);
        AddNewMaterialOrRegionToListOfSourceRegionNames();

        InputsVals = lines[SourceCommands[2]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setEventsInitialEneData"
        if(InputsVals.size() > 1){
            ui->SourceComboBoxEnergyDist->setCurrentText(InputsVals[1]);
            ui->comboBoxEnergyUnit->setCurrentText(InputsVals[0]);
            sss = ""; for(int cc = 2; cc < InputsVals.size();cc++){ sss += InputsVals[cc] + " " ; }
            ui->lineEditSpecialEnergyDistributionParameter->setText(sss);
        }
        InputsVals = lines[SourceCommands[3]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setEventsInitialMomDirData"
        if(InputsVals.size() > 1){
            ui->SourceComboBoxAngleDist->setCurrentText(InputsVals[1]);
            ui->comboBoxAngleUnit->setCurrentText(InputsVals[0]);
            sss = ""; for(int cc = 2; cc < InputsVals.size();cc++){ sss += InputsVals[cc] + " " ; }
            ui->lineEditSpecialAngulatDistributionParameter->setText(sss);
        }
        InputsVals = lines[SourceCommands[4]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        if(InputsVals.size() > 0){
            ui->UseDataFilesFor->setCurrentText(InputsVals[0]);
        }
    }else{
        AddNewMaterialOrRegionToListOfSourceRegionNames();
    }
    if(ui->checkBoxFixScoreCommands->isChecked()){

        showResultsOutput("Filling Score data...", 0);

        ui->AnalysisLineEdit_ScorOrg->setText(lines[RunAndScoreCommands[0]]); // "/RunAndScoreData/setRegionsToScore"
        ui->AnalysisLineEditVarToScore->setText(lines[RunAndScoreCommands[1]]); //"/RunAndScoreData/setQuantitiesToScore"
        InputsVals = lines[RunAndScoreCommands[6]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        if(InputsVals.size() > 0){
            ui->lineEditRadioNucleid->setText(InputsVals[0]);
            if(InputsVals.size() > 1){
                ui->comboBoxRadionuclidedataType->setCurrentText(InputsVals[1]);
                QString aa = "";
                for(int cc=2; cc < InputsVals.size();cc++){
                    aa += InputsVals[cc] + " ";
                }
                ui->RadioNucleidlineEmissionDataEdit->setText(aa);
            }
        }
        InputsVals = lines[RunAndScoreCommands[7]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        if(InputsVals.size() > 2){
            ui->RadioNucleidAdmActivitylineEdit->setText(InputsVals[1]);
            ui->comboBoxActivityUnits->setCurrentText(InputsVals[2]);
            if(InputsVals.size() > 4){
                ui->comboBoxTimeUnit->setCurrentText(InputsVals[3]);
                QString aa = "";
                for(int cc=4; cc < InputsVals.size();cc++){
                    aa += InputsVals[cc] + " ";
                }
                ui->RadioNucleidSourceResTimelineEdit->setText(aa);
            }
        }
        //ui->ScoreComboBoxAcuuracyLevel->setCurrentText(lines[RunAndScoreCommands[3]]); // "/RunAndScoreData/setAccuracyCalculationLevel"
        //ui->lineEditNumberOfEvent->setText(lines[RunAndScoreCommands[4]]); //"/RunAndScoreData/setEventNumberPerThread"

        InputsVals = lines[RunAndScoreCommands[8]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        for(int cc=0; cc < InputsVals.size();cc++){
            if(InputsVals[cc] == "AE"){ cc++; ui->comboBoxAEUnits->setCurrentText(InputsVals[cc]);}
            if(InputsVals[cc] == "SAF"){ cc++; ui->comboBoxSAFUnits->setCurrentText(InputsVals[cc]);}
            if(InputsVals[cc] == "AD"){ cc++; ui->comboBoxADUnits->setCurrentText(InputsVals[cc]);}
            if(InputsVals[cc] == "S"){ cc++; ui->comboBoxSUnits->setCurrentText(InputsVals[cc]);}
            if(InputsVals[cc] == "H"){ cc++; ui->comboBoxHUnits->setCurrentText(InputsVals[cc]);}
            if(InputsVals[cc] == "E"){ cc++; ui->comboBoxEUnits->setCurrentText(InputsVals[cc]);}
            if(InputsVals[cc] == "DCC"){ cc++; ui->comboBoxEUnits->setCurrentText(InputsVals[cc]);}
        }

        //ui->radiationEnergyFactor->setText(lines[RunAndScoreCommands[9]]);
        ui->TissueFactorLineEdit->setText(lines[RunAndScoreCommands[11]]);

    }
    if(ui->checkBoxFixRunCommands->isChecked()){

        showResultsOutput("Filling Run data...", 0);

        InputsVals = lines[RunAndScoreCommands[2]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        if(InputsVals.size() > 0){
            ui->lineEditNumberOfRanksOrThreads->setText(lines[RunAndScoreCommands[2]]); // "/RunAndScoreData/setNumberOfThreads"
        }

        InputsVals = lines[RunAndScoreCommands[5]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        if(InputsVals.size() > 0){
            if(lines[RunAndScoreCommands[5]]== "m"){
                ui->ScoreCombobowSimNumOnRanksLineEdit->setCurrentIndex(1); // "/RunAndScoreData/setSimNumOnRanks"
            }else{
                ui->ScoreCombobowSimNumOnRanksLineEdit->setCurrentIndex(0);
            }
        }
        if(QFile::exists(lines[RunAndScoreCommands[10]])){
            UserCurrentResultsDirPath = lines[RunAndScoreCommands[10]];
            ui->openResultsDirButton->setToolTip("Click to choose result directory for simulation. The current directory is " + UserCurrentResultsDirPath);
            ui->pushButtonChooseResultsDir->setToolTip("Click to choose result directory for simulation. The current directory is " + UserCurrentResultsDirPath);
            ui->pushButton->setToolTip("Open " + UserCurrentResultsDirPath +" directory");
        }
    }
    if(ui->checkBoxFixGraphsCommands->isChecked()){

        // Analyse
        showResultsOutput("Filling Analysis data...", 0);

        InputsVals = lines[AnalysisCommands[0]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/Analysis/GenerateSelfCrossGraphs"
        if(InputsVals.size() > 1){
            ui->AnalysisComboBoxGraphData->setCurrentText(InputsVals[0]);
            ui->AnalysisComboBoxGraphsType->setCurrentText(InputsVals[1]);
            if(InputsVals.size() > 3){
                ui->AnalysisLineEdit_RefName->setText(InputsVals[2]);
                ui->AnalysisLineEditRefFile->setText(InputsVals[3]);
                if(InputsVals.size() > 5){
                    ui->AnalysisLineEdit_RefName->setText(ui->AnalysisLineEdit_RefName->text()+" "+InputsVals[4]);
                    ui->AnalysisLineEditRefFile->setText(ui->AnalysisLineEditRefFile->text()+" "+InputsVals[5]);
                }
            }
        }


        for ( auto Abeg = lines.begin(); Abeg != lines.end(); ++Abeg  ){

            if(AnalysisCommands[1] == Abeg.key()){ // "/AnalysisData/generateRelativeErrGraph"
                //ui->checkBoxRelErr->setChecked(true) ;
                ui->comboBoxRelDiff->setCurrentText(lines[AnalysisCommands[1]]);
            }

            if(AnalysisCommands[2] == Abeg.key()){ // /AnalysisData/generateRelativeSDevGraph
                ui->checkBoxRelSDev->setChecked(true) ;
            }

            if(AnalysisCommands[3] == Abeg.key()){ // /AnalysisData/generateVariableRegionGraph
                ui->checkBoxRegionParameter->setChecked(true) ;
                ui->AnalysisComboBoxRegionVariable->setCurrentText(lines[AnalysisCommands[3]]);
            }
        }

        InputsVals = lines[AnalysisCommands[4]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/Analysis/GenerateSelfCrossGraphs"
        if(InputsVals.size() > 2){
            ui->AnalysisLineEditPositionsFile->setText(InputsVals[0]);
            ui->AnalysisLineEditEnergiesFile->setText(InputsVals[1]);
            ui->AnalysisLineEditMomDirsFile->setText(InputsVals[2]);
            ui->checkBoxEventsDataHisto->setChecked(true);
        }

        ui->AnalysisComboBoxSliceFor2DGraph->setCurrentText(lines["/Analysis/setSliceFor2DGraph"]);
        ui->AnalysisComboBoxBeamAxis->setCurrentText(lines["/Analysis/setBeamAxis"]);
        ui->AnalysisLineEdit_SliceID->setText(lines["/Analysis/setSliceID"]);

        InputsVals = lines[AnalysisCommands[8]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"
        if(InputsVals.size() > 8){

            if(InputsVals[0].toLower() == "yes"){
                ui->checkBoxPrintTitle->setChecked(true);
            }else{
                ui->checkBoxPrintTitle->setChecked(false);
            }
            if(InputsVals[1].toLower() == "yes"){
                ui->checkBoxUseLogE->setChecked(true);
            }else{
                ui->checkBoxUseLogE->setChecked(false);
            }
            if(InputsVals[2].toLower() == "yes"){
                ui->checkBoxUseLogVar->setChecked(true);
            }else{
                ui->checkBoxUseLogVar->setChecked(false);
            }
            if(InputsVals[3].toLower() == "yes"){
                ui->checkBoxUseGrid->setChecked(true);
            }else{
                ui->checkBoxUseGrid->setChecked(false);
            }

            ui->AnalysisLegendPoscomboBox->setCurrentText(InputsVals[4]);
            ui->lineEditLegendXYWidth->setText(InputsVals[5] +" "+ InputsVals[6]);

            if(InputsVals[7].toLower() == "yes"){
                ui->checkBoxAddErrorBar->setChecked(true);
            }else{
                ui->checkBoxAddErrorBar->setChecked(false);
            }

            ui->AnalysisComboBoxGraphsExt->setCurrentText(InputsVals[8]);
        }
    }

    updateApplicationTabs();

    return 1;
}

void MainWindow::CreateUserCommands(){

    showResultsOutput("Creating user Inputs Commands...", 1);

    // World and Materials

    if(ui->checkBoxWorldConst->isChecked()){
        GeometryData_CreateWorld = " " +
                Geometry_setWorldMaterialName + " " +
                Geometry_setWorldHalfSize + " " +
                Geometry_setWorldHalfSizeunit;
    }else {
        GeometryData_CreateWorld = " " + Geometry_setWorldHalfSize + " " ;
    }

    if(ui->radioButtonDICOM->isChecked() ){
        DICOMCommandsString =  VOXELCommandsString;
    }
    if(ui->radioButtonVoxIDs->isChecked() || ui->radioButtonTET->isChecked()){
        VoxIDsCommandsString =  VOXELCommandsString;
    }


    // Physics

    if(ui->SourceComboBoxPhysUsed->currentText() == "Construct"){
        PhysicsData_setPhysicsData = PhysicsData_setPhysicsName + " " +
                ValuesOfInputs[PhysicsData_setPhotoElectricEffectModel] + " " +
                ValuesOfInputs[PhysicsData_setComptonScatteringModel] + " " +
                ValuesOfInputs[PhysicsData_setGammaConversionModel] + " " +
                ValuesOfInputs[PhysicsData_setRayleighScatteringModel] + " " +
                ValuesOfInputs[PhysicsData_setElectronIonisationModel] + " " +
                ValuesOfInputs[PhysicsData_setElectronBremModel] + " " +
                ValuesOfInputs[PhysicsData_setHadronIonisationModel] ;
    }else {
        PhysicsData_setPhysicsData = PhysicsData_setPhysicsName ;
    }

    PhysicsData_setCutsDistance = PhysicsData_setCutsDistance;
    PhysicsData_setCutsEnergy = PhysicsData_setCutsEnergy;

    // Source
    SourceData_setSourcePosData = SourceData_setSourceSizeUnit + " " +
            SourceData_setSourceType + " " +
            SourceData_setSourceData ;
    SourceData_setSourceEneData = SourceData_setSourceEnergyUnit + " " +
            SourceData_setEnergyDistribution + " " +
            SourceData_setEnergyData ;
    SourceData_setSourceMomDirData = SourceData_setSourceAngleUnit + " " +
            SourceData_setAngleDistribution + " " +
            SourceData_setMomDirData ;
    SourceData_UseDataFiles = SourceData_setEventsNumForDataGen +" " +
            SourceData_GeneratePositions +" " +
            SourceData_GenerateEnergies +" " +
            SourceData_GenerateMomDirs
            // + " " + GeometryData_TestPointsPos +" " +
            //GeometryData_ShowBox
            ;

    // other are saved;

}
QString MainWindow::generateInputUserTextForinputFile(){

    showResultsOutput("Generating Geometry,Physics, and Radiation Source Commands text...", 0);

    QString MaterialsAndWorldData = "";
    QString GeometryData = "";
    QString PhysicsData = "";
    QString SourceData = "";
    QString ScoreData = "";
    QString AnalysisData = "";

    // World and Materials

    MaterialsAndWorldData = MaterialsDataCommandsString + "\n\n" +
            GeometryCommands[0] + " " + GeometryData_CreateWorld ;


    // Geometry

    if(ui->radioButtonDICOM->isChecked() || ui->radioButtonVoxel->isChecked() || ui->radioButtonVoxIDs->isChecked() || ui->radioButtonTET->isChecked()){
        if(ui->radioButtonVoxIDs->isChecked() || ui->radioButtonTET->isChecked()){
            GeometryData += GeometryCommands[2] + " " + Geometry_CreateVolume_GeometryFileType + " " + Geometry_CreateVolume_GeometryPath + "\n" +
                    VoxIDsCommandsString;
        }else{
            GeometryData += GeometryCommands[2] + " " + Geometry_CreateVolume_GeometryFileType;
            if(ui->radioButtonDICOM->isChecked()){
                GeometryData += "\n" + DICOMCommandsString;
            }else if(ui->radioButtonVoxel->isChecked()){
                GeometryData += "\n" + VOXELCommandsString;
            }
        }
    }
    else{
        //showResultsOutput("Geometry...---------------------------\n" + CONSCommandsString , 4);
        GeometryData += CONSCommandsString;
    }

    GeometryData += "\n" +GeometryCommands[4]+ " "+ Geometry_setGeometrySymbole + "\n";

    // Physics
    PhysicsData = PhysicsCommands[0] + " " + PhysicsData_setPhysicsData + "\n" +
            PhysicsCommands[1]+ " " + PhysicsData_setCutsDistance + "\n" ;
    if(ui->checkBoxGenerateCrossSection->isChecked()){
        if(ui->lineEditParticleNamesForCrossSection->text()!="" && ui->lineEditEnergiesForCrossSection->text()!=""){
            PhysicsData += PhysicsCommands[2]+ " " + PhysicsData_ParticleForCrossSection + " " + PhysicsData_EUnitForCrossSection + " " + PhysicsData_EnergiesForCrossSection ;
        }else{
            PhysicsData += "# " + PhysicsCommands[2]+ " " + PhysicsData_ParticleForCrossSection + " " + PhysicsData_EUnitForCrossSection + " " + PhysicsData_EnergiesForCrossSection ;
        }
    }
    PhysicsData += PhysicsCommands[3]+ " " + PhysicsData_setCutsEnergy + "\n" ;


    // Source
    SourceData = SourceCommands[0] + " " + SourceData_setParticleName + "\n" +
            SourceCommands[1] + " " + SourceData_setSourcePosData + "\n" +
            SourceCommands[2] + " " + SourceData_setSourceEneData + "\n" +
            SourceCommands[3] + " " + SourceData_setSourceMomDirData + "\n";

    if(ui->UseDataFilesFor->currentIndex() != 0 ){
        SourceData += SourceCommands[4] + " " + SourceData_UseDataFiles;
    }else{
        SourceData += "# " +SourceCommands[4] + " " + SourceData_UseDataFiles;
    }

    // Score
    ScoreData = RunAndScoreCommands[0] + " " + Score_setVolumesToScore + "\n" +
            RunAndScoreCommands[1] + " " + Score_setVariableToScore + "\n" +
            RunAndScoreCommands[2] + " " + Execution_setNumberOfRanksOrThreads + "\n" +
            //RunAndScoreCommands[3] + " " + Score_setAccuracyCalculationLevel + "\n" +
            //RunAndScoreCommands[4] + " " + Execution_setEventNumber + "\n"+
            RunAndScoreCommands[5] + " " + ValuesOfInputs[Score_setSimNumOnRanksLineEdit]+ "\n"+
            RunAndScoreCommands[10] + " " + UserCurrentResultsDirPath + "\n" +
            RunAndScoreCommands[8] + " " + Score_setQuantitiesUnits + "\n"
            //RunAndScoreCommands[9] + " " + Score_setRadiationFactors + "\n" +
            ;

    if(ui->radioButtonDICOM->isChecked() || ui->radioButtonVoxel->isChecked() || ui->radioButtonVoxIDs->isChecked()){
        if(!ui->checkBoxVoxelOrRegionLevel->isChecked()){
            ScoreData += "\n" + RunAndScoreCommands[12] + "\n";
        }
    }else{
        if(ui->comboBoxPreDefinedGeom->currentText() == "MyGeometry"){
            ScoreData += "\n" + RunAndScoreCommands[12] + "\n";
        }
    }

    ScoreData += "\n" + RunAndScoreCommands[13] + " " + ui->comboBoxSimulationRunFor->currentText()+"\n";

    if(Score_setTissueFactors != "" ){
        ScoreData += RunAndScoreCommands[11] + " " + Score_setTissueFactors + "\n" ;
    }
    if(ui->lineEditRadioNucleid->text() != "" && ui->RadioNucleidAdmActivitylineEdit->text().toDouble() != 0.){
        ScoreData +=
                RunAndScoreCommands[6] + " " + Score_setRadioNucleidDataLineEdit + "\n" +
                RunAndScoreCommands[7] + " " + Score_setRadioNucleidBiokineticsLineEdit + "\n" ;
    }

    // Analysis


    Analysis_GenerateSelfCrossGraphs = Analysis_setGraphsData +" " +
            Analysis_setCompareType +" " +
            Analysis_setRefName +" " +
            Analysis_setRefFilePath
            ;
    AnalysisData += AnalysisCommands[0] + " " + Analysis_GenerateSelfCrossGraphs + "\n" ;

    if(Analysis_GenerateRelativeErrGraph != ""){
        AnalysisData += AnalysisCommands[1] + " " + Analysis_GenerateRelativeErrGraph + "\n" ;
    }
    if(Analysis_GenerateRelativeSDevGraph == "yes"){
        AnalysisData += AnalysisCommands[2] + "\n" ;
    }
    if(Analysis_GenerateRegionsVariableGraph == "yes"){
        AnalysisData += AnalysisCommands[3] + " " + Analysis_setRegionVariableName + "\n" ;
    }
    if(Analysis_GenerateEventsDataHisto == "yes"){
        AnalysisData += AnalysisCommands[4] + " " + ui->AnalysisLineEditPositionsFile->text() + " " +
                " " + ui->AnalysisLineEditEnergiesFile->text() + " " +
                ui->AnalysisLineEditMomDirsFile->text() +" \n" ;
    }

    AnalysisData += AnalysisCommands[8] + " " + Analysis_UseLogE + " " +
            Analysis_UseLogVariable + " " +
            Analysis_UseGridXY + " " +
            Analysis_PrintTitle + " " +
            Analysis_LegendPos + " " +
            Analysis_LegendWidth + " " +
            Analysis_AddErrorBar + " " +
            Analysis_setGraphsExt + "\n"
            ;

    AnalysisData += AnalysisCommands[9]
            + " " + Score_setVariableToScore
            + " " + Analysis_setBeamAxis
            + " " + Analysis_setSliceFor2DGraph
            + " " + Analysis_setSliceID;

    //}

    QString commandsText = "\n"+
            MaterialsAndWorldData + "\n\n" +
            GeometryData + "\n\n" +
            PhysicsData + "\n\n" +
            SourceData + "\n\n" +
            ScoreData + "\n\n"+
            AnalysisData ;


    return commandsText ;

}

bool MainWindow::TestSimulateExecutableInputsToRun(){

    QStringList InputsVals;
    QString sizes = "";

    // check Computation mode and number of events
    if(ui->checkBoxRocks->isChecked()){
        if(ui->MPIOrMTOnRockscomboBox->currentText()=="MT"){

            double TotNumOfEvt = ui->lineEditNumberOfRanksOrThreads->text().toLong()*ui->lineEditNumberOfEvent->text().toLong();
            if(TotNumOfEvt > INT32_MAX ){
                QMessageBox::information(this, tr(""), "In Geant4, the total number of events \""+QString::number(TotNumOfEvt)+"\" in MT simulations should not exceed \"" + QString::number(INT32_MAX) + "\". Reduce the number of simulations \""+ ui->lineEditNumberOfRanksOrThreads->text() +"\" or number of events per simulation \""+ui->lineEditNumberOfEvent->text()+"\"");
                return false;
            }
        }
        else if(ui->MPIOrMTOnRockscomboBox->currentText()=="MPI"){

            long TotNumOfEvt = ui->lineEditNumberOfEvent->text().toLong();
            if(TotNumOfEvt > INT32_MAX ){
                QMessageBox::information(this, tr(""), "In Geant4, the total number of events \""+QString::number(TotNumOfEvt)+"\" in an MPI simulation should not exceed \"" + QString::number(INT32_MAX) + "\". Please reduve the number of events per simulation \""+ui->lineEditNumberOfEvent->text()+"\"");
                return false;
            }
        }
    }
    else{
        if(MPI_USE =="YES"){

            long TotNumOfEvt = ui->lineEditNumberOfEvent->text().toLong();
            if(TotNumOfEvt > INT32_MAX ){
                QMessageBox::information(this, tr(""), "In Geant4, the total number of events \""+QString::number(TotNumOfEvt)+"\" in an MPI simulation should not exceed \"" + QString::number(INT32_MAX) + "\". Please reduve the number of events per simulation \""+ui->lineEditNumberOfEvent->text()+"\"");
                return false;
            }
        }
        else{

            double TotNumOfEvt = ui->lineEditNumberOfRanksOrThreads->text().toLong()*ui->lineEditNumberOfEvent->text().toLong();
            if(TotNumOfEvt > INT32_MAX ){
                QMessageBox::information(this, tr(""), "In Geant4, the total number of events \""+QString::number(TotNumOfEvt)+"\" in MT simulations should not exceed \"" + QString::number(INT32_MAX) + "\". Reduce the number of simulations \""+ ui->lineEditNumberOfRanksOrThreads->text() +"\" or number of events per simulation \""+ui->lineEditNumberOfEvent->text()+"\"");
                return false;
            }
        }

    }

    // check World data
    if(ui->radioButtonVoxIDs->isChecked() || ui->radioButtonVoxel->isChecked()){
        InputsVals = ui->PhantomWorldHalfSizeslineEdit->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
        double xl,yl,zl,vxl,vyl,vzl;
        if(InputsVals.size() == 3){

            xl = QString(InputsVals[0]).toDouble() * UnitsConversionFactors[ui->comboBoxWorldSizeUnit->currentText()];
            yl = QString(InputsVals[1]).toDouble() * UnitsConversionFactors[ui->comboBoxWorldSizeUnit->currentText()];
            zl = QString(InputsVals[2]).toDouble() * UnitsConversionFactors[ui->comboBoxWorldSizeUnit->currentText()];

            //showResultsOutput(InputsVals[0] + " " + InputsVals[1] +InputsVals[2] + " " + ui->comboBoxWorldSizeUnit->currentText() , 4);

            InputsVals = getArgumentOfACommandFromText(VOXELCommands[0], 5);
            if(InputsVals.size() == 9){
                vxl = QString(InputsVals[0]).toDouble() *2* QString(InputsVals[5]).toDouble() * UnitsConversionFactors[InputsVals[8]];
                vyl = QString(InputsVals[1]).toDouble() *2* QString(InputsVals[6]).toDouble() * UnitsConversionFactors[InputsVals[8]];
                vzl = QString(InputsVals[2]).toDouble() *2* QString(InputsVals[7]).toDouble() * UnitsConversionFactors[InputsVals[8]];
                //showResultsOutput(InputsVals[0] + "*"+InputsVals[5]+ " " + InputsVals[1] + "*" +InputsVals[6]+ " " + InputsVals[2] + "*" +InputsVals[7] + " " + InputsVals[8], 4);
            }

            if(xl < vxl || yl < vyl || zl < vzl){
                sizes = "("+QString::number(xl)+";"+QString::number(vxl)+"), ("+QString::number(yl)+";"+QString::number(vyl)+"), ("+QString::number(zl)+";"+QString::number(vzl)+")";
                ui->Tab->setCurrentIndex(0);
                QMessageBox::information(this, tr(""), "The Voxelized phantom geometrical limits exceed the simulation world dimensions " + sizes +" . Edit the world dimension (i.e. 100 100 100) or voxel data (half x, y, and z).");
                return false;
            }
        }
    }
    if(ui->checkBoxWorldConst->isChecked()){
        InputsVals = ui->PhantomWorldHalfSizeslineEdit->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
        QString mat = "";mat = ui->PhantomWorldMaterialLineEdit->currentText();
        if(mat == "" || InputsVals.size() < 3){
            ui->Tab->setCurrentIndex(0);
            QMessageBox::information(this, tr(""), "Add materials, specify the world material, and world box dimension (i.e. 100 100 100)");
            return false;
        }
    }else{
        InputsVals = ui->PhantomWorldHalfSizeslineEdit->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
        if(InputsVals.size() == 0){
            ui->Tab->setCurrentIndex(0);
            QMessageBox::information(this, tr(""), "Add world geometry file (i.e. World.c++, World.cpp, World.gdml, or World.geom )");
            return false;
        }else{

            QString dd = QString::fromLocal8Bit(getFileExt(ui->PhantomWorldHalfSizeslineEdit->text().toStdString()).c_str());
            if(dd == "gdml" || dd == "geom" || dd == "c++" || dd == "cpp" || dd == "cc"){}
            else{
                if(ui->PhantomWorldHalfSizeslineEdit->text() == "MyGeometry"){

                }else{
                    ui->Tab->setCurrentIndex(0);
                    QMessageBox::information(this, tr(""), "The name of world geometry file should be one of (i.e.World.c++, World.cpp, World.gdml, World.geom, MyGeometry )");
                    return false;
                }
            }
        }
    }

    // check geometry data
    if(ui->radioButtonVoxIDs->isChecked()){
        if(QFile::exists(ui->lineEditVoxIDsFilePath->text()) ||
                QFile::exists(DoseCalcsCore_build_dir_path+"/"+ui->lineEditVoxIDsFilePath->text()) ||
                QFile::exists(DoseCalcsCore_source_dir_path+"/"+ui->lineEditVoxIDsFilePath->text())){
        }else{
            ui->Tab->setCurrentIndex(0);
            QMessageBox::information(this, tr(""), "Canno't find the voxels IDs data file. Add the file path to the LineEdit widget");
            on_pushButtonChooseVoxIDsFile_clicked();
            return false;
        }
    }
    if(ui->radioButtonTET->isChecked()){

        InputsVals = ui->lineEditVoxIDsFilePath->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
        if(InputsVals.size() < 2){
            ui->Tab->setCurrentIndex(0);
            QMessageBox::information(this, tr(""), "Add the two file paths to the LineEdit widget, the first fro tetrahedrons nodes data and the seconde for tetrahedrons elements data");
            on_pushButtonChooseVoxIDsFile_clicked();
            return false;
        }else{
            if((QFile::exists(InputsVals[0]) || QFile::exists(DoseCalcsCore_build_dir_path+"/"+InputsVals[0])|| QFile::exists(DoseCalcsCore_source_dir_path+"/"+InputsVals[0])) && (QFile::exists(InputsVals[1]) || QFile::exists(DoseCalcsCore_build_dir_path+"/"+InputsVals[1]) || QFile::exists(DoseCalcsCore_source_dir_path+"/"+InputsVals[1]))){
            }else{ui->Tab->setCurrentIndex(0);
                QMessageBox::information(this, tr(""), "Add the two file paths to the LineEdit widget, the first fro tetrahedrons nodes data and the seconde for tetrahedrons elements data");
                on_pushButtonChooseVoxIDsFile_clicked();
                return false;
            }
        }
    }

    // check source regions data
    InputsVals = ui->lineEditChosenSourceTypeData->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    if(ui->comboBoxTypeOfSources->currentText()=="Volume"){
        if(InputsVals.size() < 4){
            ui->Tab->setCurrentIndex(1);
            QMessageBox::information(this, tr(""), "Add region name with box dimension (i.e. 10 25 15)");
            return false;
        }
    }else if(ui->comboBoxTypeOfSources->currentText()=="Voxels" || ui->comboBoxTypeOfSources->currentText()=="TET"){
        if(InputsVals.size() < 1){
            ui->Tab->setCurrentIndex(1);
            QMessageBox::information(this, tr(""), "Add region/s name (i.e. Liver)");
            return false;
        }else{

            for(int aa=0; aa < InputsVals.size();aa++){

                QString ss = InputsVals[aa];
                if("allregions" == ss.toLower()){

                    aa++; for(int bb=0; bb < InputsVals[aa].toInt();bb++){
                        aa++;
                        if(aa == InputsVals.size()){
                            break;
                        }
                    }
                }
                else{
                    QStringList args = getArgumentOfACommandFromText(VOXELCommands[11], 5);
                    bool IsIn = false;

                    if(args.size() > 0){
                        if(args[0] == "yes"){for(int dd=0; dd < MaterialRegionsNames.size();dd++){ if(MaterialRegionsNames[dd] == InputsVals[aa]){ IsIn = true;}}}
                        else{ for(int dd=0; dd < DefinedRegionsNames.size();dd++){ if(DefinedRegionsNames[dd] == InputsVals[aa]){ IsIn = true;}}}
                    }else{
                        for(int dd=0; dd < MaterialRegionsNames.size();dd++){ if(MaterialRegionsNames[dd] == InputsVals[aa]){ IsIn = true;}}
                    }

                    if(IsIn == false){
                        ui->Tab->setCurrentIndex(1);
                        QMessageBox::information(this, tr(""), "The region\""+ InputsVals[aa] +"\" is not defined. Remove it from the source region LineEdit widget or add the \""+ InputsVals[aa] +"\" region data");
                        return false;
                    }
                }
            }
        }
    }
    if(ui->radioButtonGDML->isChecked() || ui->radioButtonTEXT->isChecked() || ui->radioButtonCpp->isChecked() || ui->radioButtonSTL->isChecked() || ui->radioButtonConstruct->isChecked()){
        if(ui->comboBoxTypeOfSources->currentText()=="Voxels" || ui->comboBoxTypeOfSources->currentText()=="TET"){
            ui->Tab->setCurrentIndex(1);
            QMessageBox::information(this, tr(""), "You can't use Voxels or TET source with standart volumes geometry, use Volume or other source types with standart volumes geometry");
            return false;
        }
    }else if(ui->radioButtonDICOM->isChecked() || ui->radioButtonVoxIDs->isChecked() || ui->radioButtonVoxel->isChecked() ){
        if(ui->comboBoxTypeOfSources->currentText() == "Volume" || ui->comboBoxTypeOfSources->currentText() == "TET"){
            ui->Tab->setCurrentIndex(1);
            QMessageBox::information(this, tr(""), "You can't use Volume or TET source, use Voxels or other source types with Voxelized geometry");
            return false;
        }
    }
    else if(ui->radioButtonTET->isChecked() ){
        if(ui->comboBoxTypeOfSources->currentText()=="Volume" || ui->comboBoxTypeOfSources->currentText()=="Voxels" ){
            ui->Tab->setCurrentIndex(1);
            QMessageBox::information(this, tr(""), "You can't use Volume or Voxels source, use TET or other source types with Voxelized geometry");
            return false;
        }
    }

    // check particle data
    InputsVals = ui->SourcelineEditParName->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    if(InputsVals.size() < 1){
        ui->Tab->setCurrentIndex(1);
        QMessageBox::information(this, tr(""), "Add particle name (i.e. gamma, e-, e+, alpha and proton)");
        return false;
    }else{
        for(int aa=0; aa < InputsVals.size();aa++){
            bool IsIn = false;
            for(int dd=0; dd < DefinedParticlesNames.size();dd++){ if(DefinedParticlesNames[dd] == InputsVals[aa]){ IsIn = true;}}
            if(IsIn == false){
                ui->Tab->setCurrentIndex(1);
                QMessageBox::information(this, tr(""), "The particle\""+ InputsVals[aa] +"\" is not known by DoseCalcs. Remove it from the source particle LineEdit widget");
                return false;
            }
        }

    }

    // check energy data
    if(ui->SourceComboBoxEnergyDist->currentText()=="Mono"){
        InputsVals = ui->lineEditSpecialEnergyDistributionParameter->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
        if(InputsVals.size() < 1){
            ui->Tab->setCurrentIndex(1);
            QMessageBox::information(this, tr(""), "Add an energy value");
            return false;
        }
        for(int aa=0; aa < InputsVals.size();aa++){
            QString val = InputsVals[aa];
            if(val.toDouble() == 0. || val.isEmpty() || val.isNull()){
                ui->Tab->setCurrentIndex(1);
                QMessageBox::information(this, tr(""), "The value: \"" + val + "\" not accepted.");
                return false;
            }
        }
    }

    // check momentum direction data
    InputsVals = ui->lineEditSpecialEnergyDistributionParameter->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    if(InputsVals.size() == 0){
        ui->Tab->setCurrentIndex(1);
        QMessageBox::information(this, tr(""), "Add energy value (i.e. 0.5)");
        return false;
    }

    if(ui->SourceComboBoxAngleDist->currentText()=="Directed"){
        InputsVals = ui->lineEditSpecialAngulatDistributionParameter->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
        if(InputsVals.size() > 2){
            if(InputsVals[0] != "ToPoint" && InputsVals[0] != "ThetaPhi" && InputsVals[0] != "ParallelTo" && InputsVals[0] != "ToVolume"){
                ui->Tab->setCurrentIndex(1);
                QMessageBox::information(this, tr(""), "Add how directed momentum distribution is (i.e. ThetaPhi 90 270)");
                return false;
            }
        }else{
            ui->Tab->setCurrentIndex(1);
            QMessageBox::information(this, tr(""), "Add Directed momentum distribution Data");
            return false;
        }
    }

    // check run end score data
    QString val = ui->lineEditNumberOfEvent->text();
    if(val.toDouble() == 0. || val.isEmpty() || val.isNull() || val.toInt() >= INT32_MAX){
        ui->Tab->setCurrentIndex(1);
        QMessageBox::information(this, tr(""), "The number of events: \"" + val + "\" not accepted.");
        return false;
    }
    val = ui->lineEditNumberOfRanksOrThreads->text();
    if(val.toInt() == 0. || val.isEmpty() || val.isNull()){
        ui->Tab->setCurrentIndex(1);
        QMessageBox::information(this, tr(""), "The number of sub-simulations: \"" + val + "\" not accepted.");
        return false;
    }

    return true;
}
bool MainWindow::TestMergeExecutableInputsToRun(){

    QStringList InputsVals;

    InputsVals = ui->AnalysisLineEdit_ScorOrg->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    if(InputsVals.size() == 1){
        if (ui->AnalysisLineEdit_ScorOrg->text() == "all" || ui->AnalysisLineEdit_ScorOrg->text() == "All"){}
        else{
            ui->Tab->setCurrentIndex(2);
            QMessageBox::information(this, tr(""), "set \"all\" value in volume to score");
            return false;
        }
    }
    if (InputsVals.size() > 1 ){
        bool IsIn = false;
        bool IsIn2 = false;
        for(int aa=0; aa < InputsVals.size();aa++){

            if(InputsVals[aa] == "source"){
                IsIn = true;
            }
            if(InputsVals[aa] == "target"){
                IsIn2 = true;
            }
        }
        if(!IsIn || !IsIn2){
            ui->Tab->setCurrentIndex(2);
            QMessageBox::information(this, tr(""), "Please add \"source\" word followed by source volumes, and \"target\" word followed by target volumes");
            return false;
        }
    }else if(InputsVals.size() == 0){
        ui->Tab->setCurrentIndex(2);
        QMessageBox::information(this, tr(""), "Please add \"source\" word followed by source volumes, and \"target\" word followed by target volumes");
        return false;
    }

    InputsVals = ui->AnalysisLineEditVarToScore->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    if(InputsVals.size() < 1){
        ui->Tab->setCurrentIndex(2);
        QMessageBox::information(this, tr(""), "Add a quantity to score (i.e. SAF)");
        return false;
    }else{
        for(int aa=0; aa < InputsVals.size();aa++){
            bool IsIn = false;
            for(int dd=0; dd < VarToScorellist.size();dd++){ if(VarToScorellist[dd] == InputsVals[aa]){ IsIn = true;}}
            if(IsIn == false){
                ui->Tab->setCurrentIndex(2);
                QMessageBox::information(this, tr(""), "The quantity\""+ InputsVals[aa] +"\" is not known by DoseCalcs. Remove it from the quantities to score LineEdit widget");
                return false;
            }
        }
    }

    InputsVals = ui->lineEditRadioNucleid->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    QStringList InputsVals2 = ui->RadioNucleidlineEmissionDataEdit->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);

    if(InputsVals.size() < 1){}else{
        if(ui->comboBoxRadionuclidedataType->currentText() == "File"){
            if(InputsVals2.size() > 1){
                if(!QFile::exists(InputsVals2[0])){
                    ui->Tab->setCurrentIndex(2);
                    QMessageBox::information(this, tr(""), "There is no file contains the "+InputsVals[0]+ " radiation data");
                    //return false;
                }
            }
        }else{
            for(int aa=0; aa < InputsVals2.size();aa++){
                bool IsIn = false;
                for(int dd=0; dd < DefinedParticlesNames.size();dd++){ if(DefinedParticlesNames[dd] == InputsVals2[aa]){ IsIn = true;}}
                if(IsIn == false){
                    ui->Tab->setCurrentIndex(2);
                    QMessageBox::information(this, tr(""), "Warning: The particle\""+ InputsVals2[aa] +"\" is not known by DoseCalcs. Remove it from the radiation factors LineEdit widget");
                    break;
                    //return false;
                }
                if(IsIn == false){
                    break;
                }
                aa+=2;
            }
        }
    }

    if(QFile::exists(UserCurrentResultsDirPath)){}else{
        while(!QFile::exists(UserCurrentResultsDirPath)){
            ui->Tab->setCurrentIndex(1);
            QMessageBox::information(this, tr(""), "Canno't find the results directory \""+ UserCurrentResultsDirPath +"\". Please choose a new directory path.");
            on_openResultsDirButton_clicked();
        }
    }

    return true;
}
bool MainWindow::TestVisualizingInputsToRun(){

    QStringList InputsVals;

    if(ui->radioButtonDICOM->isChecked() || ui->radioButtonVoxIDs->isChecked() ){
        InputsVals = getArgumentOfACommandFromText(VOXELCommands[9], 5);
    }else if(ui->radioButtonTET->isChecked() ){
        InputsVals = getArgumentOfACommandFromText(VOXELCommands[10], 5);
    }

    if(InputsVals.size() < 1){
        if(QMessageBox::Yes == QMessageBox::question(this, tr("Voxelized or tetrahedral phantom visualization"),
                                                     tr("For voxelized or tetrahedral phantoms, it is preferred to use the limits of planes to visualize in Qt Editor. To set the limits, click yes, then edit geometry, add the limits data. Add the limits ?")))
        {
            ui->Tab->setCurrentIndex(0);
            on_EditGeometyDataBtn_clicked();
            return false;
        }
        else{return true;}
    }

    return true;
}
bool MainWindow::ShowImportantSimulationData(){

    //showResultsOutput("31 : \n", 0);

    QString ImportantSimulationInputs = "";

    if(ui->checkBoxWorldConst->isChecked()){
        ImportantSimulationInputs += "*** World Data : " + ui->PhantomWorldMaterialLineEdit->currentText() + " " + ui->PhantomWorldHalfSizeslineEdit->text() + " (" +ui->comboBoxWorldSizeUnit->currentText()+")\n\n";
    }else{
        ImportantSimulationInputs += "*** World Data are imported from : " + ui->PhantomWorldHalfSizeslineEdit->text() + "\n\n";
    }

    if(ui->radioButtonVoxIDs->isChecked() || ui->radioButtonVoxel->isChecked()){
        ImportantSimulationInputs += "*** Parammetrization : " + ValuesOfInputs[getArgumentOfACommandFromText(VOXELCommands[0], 5)[3]]+"("+getArgumentOfACommandFromText(VOXELCommands[0],5)[3]+")\n\n";
    }
    else if(ui->radioButtonTET->isChecked()){

        ImportantSimulationInputs += "*** Tetrahedral phantom \n";
    }
    else if(ui->radioButtonGDML->isChecked()|| ui->radioButtonTEXT->isChecked()|| ui->radioButtonCpp->isChecked() || ui->radioButtonConstruct->isChecked() ){

        ImportantSimulationInputs += "\n";
        ImportantSimulationInputs += "Volumes : ";
        for(int cc=0; cc < VolsNames.size() ;cc++){
            ImportantSimulationInputs += " " + VolsNames[cc] ;
        }

        ImportantSimulationInputs += "\n";
        ImportantSimulationInputs += "Solids : ";
        for(int cc=0; cc < SolidsNames.size() ;cc++){
            ImportantSimulationInputs += " " + SolidsNames[cc] ;
        }
    }

    ImportantSimulationInputs +=
            "*** Geometry Symbole : " + ui->lineEditGeometrySymbole->text() + "\n\n";

    ImportantSimulationInputs +=
            "*** Physics : " + ui->SourceComboBoxPhysUsed->currentText() + "\n"+
            "*** Cuts In Range: " + ui->SourceLineEditDistanceCut->text()+ "\n"+
            "*** Energy Range: " + ui->SourceLineEditEnergyCut->text()+ "\n\n";
    ImportantSimulationInputs +=
            "*** Initial Particles : "+ui->SourcelineEditParName->text() + "\n"+
            "*** Initial Position : "+ui->comboBoxTypeOfSources->currentText() +" : " + ui->lineEditChosenSourceTypeData->text() +"\n" +
            "*** Initial Energy : "+ui->SourceComboBoxEnergyDist->currentText() + " : " + QString::number(ui->lineEditSpecialEnergyDistributionParameter->text().toDouble()) + "\n" +
            "*** Initial Momentum Direction : "+ui->SourceComboBoxAngleDist->currentText() + " : " + ui->lineEditSpecialAngulatDistributionParameter->text() + "\n\n";
    ImportantSimulationInputs +=
            "*** Simulation Per Thread or Rank: " + ui->ScoreCombobowSimNumOnRanksLineEdit->currentText() + "\n"+
            "*** Number Of Simulations : " + QString::number(ui->lineEditNumberOfRanksOrThreads->text().toInt()) + "\n"+
            "*** Number Of Events : " + QString::number(ui->lineEditNumberOfEvent->text().toInt()) + "\n"+
            "*** Result directory : " + UserCurrentResultsDirPath + "\n\n";
    ImportantSimulationInputs +=
            "*** DoseCalcs Job Name: " + ConstructDoseCalcsJobName()+ "\n\n";

    //ImportantSimulationInputs += "*** simulation number: " + SimulationsNumberDataString+ "\n\n";

    if(QMessageBox::Yes == QMessageBox::question(this, tr("The simulation will be executed for:"), ImportantSimulationInputs)){
        return true;
    }else{
        return false;
    }
}
QString MainWindow::ShowMessageBoxAndGetPath(QString message ){

    QString DirPath = "";

    QDialog * d = new QDialog(); d->setWindowTitle("A Needed Package Not Found");
    QVBoxLayout * vbox = new QVBoxLayout();
    //QStringList RegionVars=;
    //QPushButton * DirPathButton = new QPushButton("."); vbox->addWidget(DirPathButton);
    QLabel * Path = new QLabel(message); vbox->addWidget(Path);

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));
    vbox->addWidget(buttonBox);
    d->setLayout(vbox);

    if(d->exec() == QDialog::Accepted)
    {
        DirPath = QFileDialog::getExistingDirectory(0, (message), DoseCalcsCore_build_dir_path);
        //DirPath = QFileDialog::getOpenFileName( this, message, DoseCalcsCore_build_dir_path, tr("All files (*.sh*)") );

        if(DirPath.isEmpty() || DirPath.isNull()){

        }else{

        }

    }else{
        DirPath = "";
    }

    return DirPath;
}
// this Slot gets called whenever the process has something to say... then in header declare it in the SLOT part, and the process is initialized and connected to system calls in the constructor of this class
void MainWindow::processOutput()
{
    //showResultsOutput(" .............................. ", 4);

    QString outputText = CoreProcess.readAllStandardOutput();
    ui->outputTextConsole->appendPlainText(outputText); ui->outputTextConsole->update();
    /*
    bytes = CoreProcess.readAllStandardOutput();

    QStringList lines = QString(bytes).split(" ");
    foreach (QString line, lines) {
        ui->outputTextConsole->appendPlainText(line); ui->outputTextConsole->update();

        qDebug() << line;
    }
    //showResultsOutput( process->readAllStandardOutput() , 4 );
    */
}
void MainWindow::showResultsOutput(QString text, int level){

    if(level == 0){
        QTextStream(stdout) << text << "\n";
        return;
    }
    if(level == 1){
        QTextStream(stdout) << " ---------------------> " << text << "\n";
        ui->outputTextConsole->appendPlainText(" ---------------------> "+text);
    }
    if(level == 4){
        //QTextStream(stdout) << text << "\n";
        ui->outputTextConsole->appendPlainText("\n"+text);
    }
    if(level == 3){
        QTextStream(stdout) << "!!!!!!!!!!!!!!!!! " << text << "\n";
        ui->outputTextConsole->appendPlainText("!!!!!!!!!!!!!!!!! "+text);
    }

}
QString MainWindow::SetToASpecificMPIRank(){



    //QProcess process;

    //// Run the Rocks command to list hosts
    ////process.start("rocks list host");
    //process.start("ls");

    //if (!process.waitForFinished()) {
    //    //qDebug() << "Error: Unable to execute Rocks command.";
    //    //return nodes;
    //}

    //// Read the output from the process
    //QString output = process.readAllStandardOutput();

    //// Split the output into lines
    //QStringList lines = output.split("\n", Qt::SkipEmptyParts);
    //QStringList nodesstringlist;

    ////QTextStream(stdout) << " lines size " << lines.size() << "\n";

    //// Display the hostnames
    ////qDebug() << "List of Hostnames:";
    //for(int dd=0; dd < lines.size();dd++){
    //    QStringList fields = lines[dd].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
    //    if (fields.length() > 0) {
    //        //qDebug() << " fields.at(1) " << fields.at(0) << "\n";
    //        nodesstringlist.push_back(fields.at(0));
    //    }
    //}
    //ListNodes = new QComboBox(); ListNodes->addItems(nodesstringlist);
    //ListNodes->setToolTip("Choose a node from the listed nodes");
    //ListNodes->setCurrentText("");
    //connect(ListNodes, SIGNAL(textActivated(QString)), this, SLOT(on_ListNodes_textActivated(QString)));
    //GraphLayout->addWidget(ListNodes, jj,++ii,1,1);

    QDialog * d = new QDialog(); d->setWindowTitle("Specify hostnames \"i.e., node-0-2\"");

    QGridLayout* GraphLayout = new QGridLayout;

    Textnodes = new QLineEdit(); Textnodes->setToolTip("Add specific hostnames");
    QPushButton* qstatfbutton = new QPushButton(); qstatfbutton->setText("List hostnames");
    qstatfbutton->setToolTip("List all nodes (hostnames) in Terminal ");
    connect(qstatfbutton, SIGNAL(clicked()), this, SLOT(runTerminalCommandSlot()));


    int ii = 0, jj=0;
    GraphLayout->addWidget(Textnodes, jj,ii,1,1);
    GraphLayout->addWidget(qstatfbutton, jj,++ii,1,1);

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    GraphLayout->addWidget(buttonBox);

    d->setLayout(GraphLayout);

    int result = d->exec();
    QString nodes = "";

    if(result == QDialog::Accepted)
    {
        nodes = Textnodes->text();

        if(nodes.isEmpty() || nodes == " "|| nodes == "  " || nodes == "   "){
            return "";
        }else{
            return " -l hostname="+nodes;
        }
    }else{
        return "";
    }


}

// world
void MainWindow::on_checkBoxWorldConst_clicked(bool checked)
{
    if(checked == true){
        ui->checkBoxWorldConst->setCheckState(Qt::Checked);
        ui->checkBoxWorldConst->setText("Construct");
        ui->pushButtonChoosWorldFile->setVisible(false);
        ui->PhantomWorldMaterialLineEdit->setVisible(true);
        ui->comboBoxWorldSizeUnit->setVisible(true);
        ui->PhantomWorldHalfSizeslineEdit->setPlaceholderText("100 100 100");
        ui->PhantomWorldHalfSizeslineEdit->setText("100 100 100");
        ui->pushButtonEditGeometryFile->setVisible(false);

    }else if(checked == false){

        ui->checkBoxWorldConst->setCheckState(Qt::Unchecked);
        ui->checkBoxWorldConst->setText("Import");

        ui->PhantomWorldMaterialLineEdit->setVisible(false);
        ui->comboBoxWorldSizeUnit->setVisible(false);

        ui->pushButtonChoosWorldFile->setVisible(true);
        //ui->PhantomWorldHalfSizeslineEdit->setText("");
        ui->PhantomWorldHalfSizeslineEdit->setPlaceholderText("World.gdml/World.geom/World.c++/World.cpp");
        ui->PhantomWorldHalfSizeslineEdit->setText(DoseCalcsCore_source_dir_path+"/src/G4TCPPGeometryFormat.cc");
        ui->pushButtonEditGeometryFile->setVisible(true);

    }
}
void MainWindow::on_pushButtonEditGeometryFile_clicked()
{
    on_pushButtonEditGeomFile_clicked(); // if we click edit again, save before.

    QString flnm = QString::fromLocal8Bit(getFileNameFromPath(ui->PhantomWorldHalfSizeslineEdit->text().toStdString()).c_str())+"."+QString::fromLocal8Bit(getFileExt(ui->PhantomWorldHalfSizeslineEdit->text().toStdString()).c_str());
    QString ext = QString::fromLocal8Bit(getFileExt(ui->PhantomWorldHalfSizeslineEdit->text().toStdString()).c_str());

    if(ext == "gdml" || ext == "geom" || ext == "text" || ext == "c++" || ext == "cpp"  || ext == "cc"){

        if(!QFile::exists(ui->PhantomWorldHalfSizeslineEdit->text())){
            on_pushButtonChoosWorldFile_clicked();
        }

        EditFlag = 5;
        ui->tabWidget->setTabText(0,flnm);
        ui->GeometryFileTextEdit->clear();
        showResultsOutput("", 4);
        ui->GeometryFileTextEdit->setPlainText(fileManagerObject->ReadTextFromFileInOneString(ui->PhantomWorldHalfSizeslineEdit->text()));
        ui->tabWidget->setCurrentIndex(0);
    }
}
void MainWindow::on_WorldDataShowButtonpushButton_clicked()
{
    QString WorldData = "";
    if(ui->checkBoxWorldConst->isChecked()){
        WorldData = GeometryCommands[0] + " " +
                ui->PhantomWorldMaterialLineEdit->currentText() + " " +
                ui->PhantomWorldHalfSizeslineEdit->text() + " " +
                ui->comboBoxWorldSizeUnit->currentText();
    }else {
        WorldData = GeometryCommands[0] + " " + ui->PhantomWorldHalfSizeslineEdit->text() + " " ;
    }

    ui->outputTextConsole->appendPlainText("\n============ World Volume Command ============\n");
    ui->outputTextConsole->appendPlainText(WorldData);
    ui->tabWidget->setCurrentIndex(1);
}
void MainWindow::on_pushButtonChoosWorldFile_clicked()
{
    QString chosen_DirPath = QFileDialog::getOpenFileName( this, tr("Choose VoxIDs file (region or material IDs file) \".dat\""), DoseCalcsCore_build_dir_path+"/"+ScriptDirectoryName, tr("All files (*.*)") );

    if(chosen_DirPath.isEmpty() || chosen_DirPath.isNull()){

    }else{
        ui->PhantomWorldHalfSizeslineEdit->setText(chosen_DirPath);
    }
}

// edit and save from and to files buttons
void MainWindow::on_EditGeometyDataBtn_clicked()
{
    on_pushButtonEditGeomFile_clicked();

    if(ui->radioButtonGDML->isChecked() || ui->radioButtonTEXT->isChecked() || ui->radioButtonCpp->isChecked()){

        ui->tabWidget->setTabText(0,"Volumes reading macros");
        ui->GeometryFileTextEdit->clear();
        showResultsOutput("", 4);
        ui->GeometryFileTextEdit->setPlainText(CONSCommandsString);

        RemoveDynamiqueGeomAndMatFrame();
        showConstructVolumeGDMLTEXTCPPCommandsFrame();
    }
    if(ui->radioButtonSTL->isChecked()){

        ui->tabWidget->setTabText(0,"Volumes creation by .stl solid macros");
        ui->GeometryFileTextEdit->clear();
        showResultsOutput("", 4);
        ui->GeometryFileTextEdit->setPlainText(CONSCommandsString);

        RemoveDynamiqueGeomAndMatFrame();
        showConstructSTLCommandsFrame();
    }
    if(ui->radioButtonVoxIDs->isChecked() || ui->radioButtonDICOM->isChecked() || ui->radioButtonVoxel->isChecked()){

        ui->tabWidget->setTabText(0,"Voxelized data macros");
        ui->GeometryFileTextEdit->clear();
        showResultsOutput("getting the saved Voxelized Geometry Commands", 4);
        ui->GeometryFileTextEdit->setPlainText(VOXELCommandsString);

        RemoveDynamiqueGeomAndMatFrame();
        showConstructVoxelizedCommandsFrame();
    }
    if(ui->radioButtonTET->isChecked()){

        ui->tabWidget->setTabText(0,"TET data macros");
        ui->GeometryFileTextEdit->clear();
        showResultsOutput("getting the saved TET Geometry Commands", 4);
        ui->GeometryFileTextEdit->setPlainText(VOXELCommandsString);

        RemoveDynamiqueGeomAndMatFrame();
        showConstructVoxelizedCommandsFrame();
    }
    if(ui->radioButtonConstruct->isChecked()){

        ui->tabWidget->setTabText(0,"Volumes creation macros");
        ui->GeometryFileTextEdit->clear();
        showResultsOutput("getting the saved Volumes Constructing Commands", 4);
        ui->GeometryFileTextEdit->setPlainText(CONSCommandsString);

        RemoveDynamiqueGeomAndMatFrame();
        showConstructSolidAndVolumeCommandsFrame();
    }
    ui->tabWidget->setCurrentIndex(0);

    EditFlag = 2;
}
void MainWindow::on_ViewGeometyDataBtn_clicked()
{
    on_actionVisualization_triggered();
}
void MainWindow::on_GeometryDataShowButton_clicked()
{
    ui->outputTextConsole->appendPlainText("\n============ Geometry Commands ============\n");
    if(ui->radioButtonTET->isChecked() || ui->radioButtonVoxIDs->isChecked() || ui->radioButtonDICOM->isChecked() || ui->radioButtonVoxel->isChecked()){
        ui->outputTextConsole->appendPlainText(VOXELCommandsString);
    }else{
        ui->outputTextConsole->appendPlainText(CONSCommandsString);
    }
    ui->tabWidget->setCurrentIndex(1);
}
void MainWindow::on_comboBoxPreDefinedGeom_currentTextChanged(const QString &arg1)
{
    if(ui->checkBoxUsePreDefinedGeom->isChecked()){
        MacroFilePath = GUIPackagesAndFilesDirPath+"/PreDefinedGeometry/macros"+PreDefinedGeomMap[arg1]+".mac";
        FillComponentsFromInputsFile(MacroFilePath);
    }

    if(ui->radioButtonDICOM->isChecked() || ui->radioButtonVoxel->isChecked() || ui->radioButtonVoxIDs->isChecked()){
        ui->checkBoxVoxelOrRegionLevel->setVisible(true);
        ui->checkBoxVoxelOrRegionLevel->setChecked(true);
    }else{
        if(ui->comboBoxPreDefinedGeom->currentText() == "MyGeometry"){
            ui->checkBoxVoxelOrRegionLevel->setVisible(true);
            ui->checkBoxVoxelOrRegionLevel->setChecked(true);
        }else{
            ui->checkBoxVoxelOrRegionLevel->setVisible(false);
            ui->checkBoxVoxelOrRegionLevel->setChecked(true);
        }
    }

    //ReadWTFactor(ICRPDATAPath+"/ICRP110RegionsData");
}
void MainWindow::on_checkBoxVoxelOrRegionLevel_stateChanged(int arg1)
{
    if(ui->checkBoxVoxelOrRegionLevel->isChecked()){
        //ui->checkBoxVoxelOrRegionLevel->setCheckState(Qt::Checked);
        ui->checkBoxVoxelOrRegionLevel->setText("Region Level");
    }else{
        //ui->checkBoxVoxelOrRegionLevel->setCheckState(Qt::Unchecked);
        ui->checkBoxVoxelOrRegionLevel->setText("Voxel Level");
    }
}

void MainWindow::on_comboBoxPredefinedSources_currentTextChanged(const QString &arg1)
{
    if(arg1 == "External RLAT(Ang=-90) for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("RLAT 25 0 0 Rectangle X 14 90");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 180 ");
    }
    else if(arg1=="External LLAT(Ang=90) for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("LLAT 25 0 0 Rectangle X 14 90 Z 180");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 360");
    }
    else if(arg1=="External PA(Ang=+-180) for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("PA 0 -25 0 Rectangle Y 90 20");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 90");
    }
    else if(arg1=="External AP(Ang=0) for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("AP 25 0 0 Rectangle X 20 90 Z 90");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 270");
    }
    else if(arg1=="External ROT for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Surface");
        ui->lineEditChosenSourceTypeData->setText("ROT 0 0 0 Cylinder 40 90");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ParallelTo Z 0 0");
    }
    else if(arg1=="External ISO for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Surface");
        ui->lineEditChosenSourceTypeData->setText("ISO 0 0 0 Sphere 110");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ToVolume 0 0 0 20 12 90");
    }
    else if(arg1=="External RLAT(Ang=-90) for Voxelized and Mesh-type Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("RLAT -40 0 0 Rectangle X 26 100");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 360");
    }
    else if(arg1=="External LLAT(Ang=90) for Voxelized and Mesh-type Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("LLAT 40 0 0 Rectangle X 26 100");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 180 ");
    }
    else if(arg1=="External PA(Ang=+-180) for Voxelized and Mesh-type Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("PA 30 0 0 Rectangle X 41 100 Z 90");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 270");
    }
    else if(arg1=="External AP(Ang=0) for Voxelized and Mesh-type Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("AP 0 -30 0 Rectangle Y 100 41");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 90");
    }
    else if(arg1=="External ROT for Voxelized and Mesh-type Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Surface");
        ui->lineEditChosenSourceTypeData->setText("ROT 0 0 0 Cylinder 50 100");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ParallelTo Z 0 0");
    }
    else if(arg1=="External ISO for Voxelized and Mesh-type Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Surface");
        ui->lineEditChosenSourceTypeData->setText("ISO 0 0 0 Sphere 110");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ToVolume 0 0 0 40 20 100");
    }
    else if(arg1=="External Ang=-75 for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("Ang=-75 25 0 0 Rectangle X 14 90 Z 15");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 195");
    }
    else if(arg1=="External Ang=-60 for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("Ang=-60 25 0 0 Rectangle X 16 90 Z 30");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 210 ");
    }
    else if(arg1=="External Ang=-45 for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("Ang=-45 25 0 0 Rectangle X 19 90 Z 45");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 225");
    }
    else if(arg1=="External Ang=-30 for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("Ang=-30 25 0 0 Rectangle X 20 90 Z 60");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 240 ");
    }
    else if(arg1=="External Ang=-15 for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("Ang=-15 25 0 0 Rectangle X 20 90 Z 75");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 255");
    }
    else if(arg1=="External Ang=+15 for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("Ang=+15 25 0 0 Rectangle X 20 90 Z 105");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 285");
    }
    else if(arg1=="External Ang=+30 for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("Ang=+30 25 0 0 Rectangle X 20 90 Z 120");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 300");
    }
    else if(arg1=="External Ang=+45 for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("Ang=+45 25 0 0 Rectangle X 20 90 Z 135");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 315");
    }
    else if(arg1=="External Ang=+60 for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("Ang=+60 25 0 0 Rectangle X 16 90 Z 150");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 330");
    }
    else if(arg1=="External Ang=+75 for Stylized Phantoms"){
        ui->comboBoxTypeOfSources->setCurrentText("Plane");
        ui->lineEditChosenSourceTypeData->setText("Ang=+75 25 0 0 Rectangle X 14 90 Z 165");
        ui->SourceComboBoxAngleDist->setCurrentText("Directed");
        ui->lineEditSpecialAngulatDistributionParameter->setText("ThetaPhi 90 345");
    }
}
void MainWindow::on_checkBoxUsePreDefinedGeom_stateChanged(int arg1)
{
    if(ui->checkBoxUsePreDefinedGeom->isChecked()){
        if(!QFile::exists(GUIPackagesAndFilesDirPath+"/PreDefinedGeometry")){
            QString n = "To use the predefined geometries macros, you should download the files "
                        "unzip them under the directory \"" + GUIPackagesAndFilesDirPath+
                    "\". This can be automatically from \"Installations\" window by clicking on "
                    "\"Download DoseCalcs Supplementary Files \"button.\n"
                    "You can check this files by clicking on \"Open Files and Installation Packages Directory\" button"
                    "\n\n You read the message, Download PreDefinedGeometry Data ?"
                    ;
            if (QMessageBox::Yes == QMessageBox::question(this, tr("PreDefined Geometry Files not found"), n)){
                on_actionInstallations_triggered();
                if(!QFile::exists(GUIPackagesAndFilesDirPath+"/PreDefinedGeometry")){
                    ui->checkBoxUsePreDefinedGeom->setCheckState(Qt::Unchecked);
                }
                else{
                    on_comboBoxPreDefinedGeom_currentTextChanged(ui->comboBoxPreDefinedGeom->currentText());
                }
            }else{
                ui->checkBoxUsePreDefinedGeom->setCheckState(Qt::Unchecked);
            }
        }else{
            on_comboBoxPreDefinedGeom_currentTextChanged(ui->comboBoxPreDefinedGeom->currentText());
        }
    }
}

// edit and save from and to String buttons
void MainWindow::on_MaterialsDataEditButton_clicked()
{
    on_pushButtonEditGeomFile_clicked();

    ui->tabWidget->setTabText(0,"Materials creation macros");
    ui->GeometryFileTextEdit->clear();
    showResultsOutput("getting the saved Materials Data Commands", 4);
    ui->GeometryFileTextEdit->setPlainText(MaterialsDataCommandsString);

    ui->tabWidget->setCurrentIndex(0);

    RemoveDynamiqueGeomAndMatFrame();
    showConstructMaterialsFrame();

    EditFlag = 3;
}
void MainWindow::on_MaterialsDataShowButton_clicked()
{
    ui->outputTextConsole->appendPlainText("\n============ Materials Commands ============\n");
    ui->outputTextConsole->appendPlainText(MaterialsDataCommandsString);
    ui->tabWidget->setCurrentIndex(1);
}

void MainWindow::RemoveDynamiqueGeomAndMatFrame(){

    if ( ui->fraHelpGUIInput->layout() != NULL )
    {
        QLayoutItem* item;
        while ( ( item = ui->fraHelpGUIInput->layout()->takeAt( 0 ) ) != NULL )
        {
            delete item->widget();
            delete item;
        }

        delete ui->fraHelpGUIInput->layout();
        ui->fraHelpGUIInput->setVisible(false);
    }

}

void MainWindow::showConstructMaterialsFrame(){

    //QLayout* framLay = new QGridLayout(ui->fraHelpGUIInput);
    framLay = new QGridLayout;

    // first row

    MatElemComb = new QComboBox(); MatElemComb->addItems(PeriodicTableElementsSymbol);
    btnSaveElem= new QPushButton(); btnSaveElem->setText("Save Elem"); connect(btnSaveElem, SIGNAL(clicked()), this, SLOT(btnSaveElem_slot()));
    MatElemComb->setToolTip("Choose the elements to be used after in material creation"); btnSaveElem->setToolTip("Add element to the created list of elements");

    MaterialName = new QLineEdit(); MaterialName->setPlaceholderText("Mat Name");
    MaterialID = new QLineEdit(); MaterialID->setPlaceholderText("Mat ID");
    MaterialEleNumber = new QLineEdit(); MaterialEleNumber->setPlaceholderText("Numb");
    MaterialDensity = new QLineEdit(); MaterialDensity->setPlaceholderText("Density");
    NumFraCheckBox = new QCheckBox(); connect(NumFraCheckBox, SIGNAL(clicked()), this, SLOT(CheckBoxFracNum_slot()));
    btnAddMat = new QPushButton(); btnAddMat->setText("Add Mat"); connect(btnAddMat, SIGNAL(clicked()), this, SLOT(btnAddMat_slot()));

    MaterialName->setToolTip("Add material name"); MaterialID->setToolTip("Add material specific ID"); MaterialEleNumber->setToolTip("Add number of elements that will constitue the material"); MaterialDensity->setToolTip("Add material density value and units");
    NumFraCheckBox->setToolTip("Choose by by which elements method will be added to the material, by mass fraction or number"); btnAddMat->setToolTip("Add material to the created list of materials");

    ElemComb = new QComboBox(); for (int kk =0 ;kk < ElementsNames.size() ; kk++) {ElemComb->addItem(ElementsNames[kk]);}
    EleFraOrNumInMat = new QLineEdit(); EleFraOrNumInMat->setPlaceholderText("Numb");
    btnAddElem = new QPushButton(); btnAddElem->setText("Add Elem"); connect(btnAddElem, SIGNAL(clicked()), this, SLOT(btnAddElem_slot()));
    ElemComb->setToolTip("Choose the element constituant to the material"); EleFraOrNumInMat->setToolTip("material mass fraction(ie. 50) or number (ie. 2) "); btnAddElem->setToolTip("Add element to the material");

    NistMatComb = new QComboBox(); NistMatComb->addItems(NistMaterialsDataBase);
    NistMaterialID = new QLineEdit(); NistMaterialID->setPlaceholderText("Mat ID");
    btnAddNistMat = new QPushButton(); btnAddNistMat->setText("add NIST Mat");
    connect(btnAddNistMat, SIGNAL(clicked()), this, SLOT(btnAddNistMat_slot()));
    NistMatComb->setToolTip("Add NIST material name"); btnAddNistMat->setToolTip("Add NIST material name to the created list of materials");

    int ii = 0, jj = 0;
    framLay->addWidget(MatElemComb, jj,ii,1,1);
    framLay->addWidget(btnSaveElem, jj,++ii,1,1);

    ii = 0; jj++;
    framLay->addWidget(MaterialName, jj,ii,1,1);
    framLay->addWidget(MaterialID, jj,++ii,1,1);
    framLay->addWidget(MaterialEleNumber, jj,++ii,1,1);
    framLay->addWidget(MaterialDensity, jj,++ii,1,1);
    framLay->addWidget(NumFraCheckBox , jj,++ii,1,1);
    framLay->addWidget(btnAddMat, jj,++ii,1,1);

    framLay->addWidget(ElemComb, jj,++ii,1,1);
    framLay->addWidget(EleFraOrNumInMat, jj,++ii,1,1);
    framLay->addWidget(btnAddElem, jj,++ii,1,1);

    ii = 0; jj++;
    framLay->addWidget(NistMatComb, jj,ii,1,1);
    framLay->addWidget(NistMaterialID, jj,++ii,1,1);
    framLay->addWidget(btnAddNistMat, jj,++ii,1,1);

    ui->fraHelpGUIInput->setLayout(framLay);
    ui->fraHelpGUIInput->setVisible(true);

    MatCommandToAddElem = "";

}
void MainWindow::showConstructVoxelizedCommandsFrame(){

    framLay = new QGridLayout;
    int ii=0, jj=0;
    /*
    if(ui->radioButtonVoxIDs->isChecked()){

        BtnAddVolPath = new QPushButton(); BtnAddVolPath->setText("Add IDs file"); connect(BtnAddVolPath, SIGNAL(clicked()), this, SLOT(BtnAddVolPath_slot()));
        PhyVolName = new QLineEdit();
        if(QFile::exists(Geometry_CreateVolume_GeometryPath)){
            PhyVolName->setText(Geometry_CreateVolume_GeometryPath);
        }else{
            PhyVolName->setText(DoseCalcsCore_build_dir_path+"/Scripts/VoxIDs.dat");
        }
        BtnEditVolFile = new QPushButton(); BtnEditVolFile->setText("set IDs File"); connect(BtnEditVolFile, SIGNAL(clicked()), this, SLOT(BtnEditVolFile_slot()));
        framLay->addWidget(BtnAddVolPath, jj,ii,1,1);
        framLay->addWidget(PhyVolName, jj,++ii,1,1);
        framLay->addWidget(BtnEditVolFile, jj,++ii,1,1);
        ii = 0; jj++;
    }
    else{
        ii=0, jj=0;
    }
*/
    if(!ui->radioButtonDICOM->isChecked() && !ui->radioButtonTET->isChecked()){

        XYZVoxelsNumb = new QLineEdit(); XYZVoxelsNumb->setPlaceholderText("nX nY nZ");
        ParamType = new QComboBox(); QStringList NIST_Mat=(QStringList()<<"Phantom parametrization"<<"Phantom nested parametrization"); ParamType->addItems(NIST_Mat);
        XYZVoxelsHalfSize = new QLineEdit();  XYZVoxelsHalfSize->setPlaceholderText("Hx Hy Hz unit");
        LogVolMatName = new QComboBox(); for (int kk =0 ;kk < MaterialsNames.size() ; kk++) {LogVolMatName->addItem(MaterialsNames[kk]);}
        VoxContainerPos = new QLineEdit(); VoxContainerPos->setPlaceholderText("ContPos unit");
        VoxContainerRot = new QLineEdit(); VoxContainerRot->setPlaceholderText("ContRot unit");
        btnAddVoxelsData = new QPushButton(); btnAddVoxelsData->setText("Add Voxels Data"); connect(btnAddVoxelsData, SIGNAL(clicked()), this, SLOT(btnAddVoxelsData_slot()));

        QStringList InputsVals = getArgumentOfACommandFromText(VOXELCommands[0], 5);
        if(InputsVals.size() == 9){
            XYZVoxelsNumb->setText(InputsVals[0]+" "+InputsVals[1]+" "+InputsVals[2]);
            if(InputsVals[3].toInt() == 0){ParamType->setCurrentIndex(0);}else{ParamType->setCurrentIndex(1);}
            LogVolMatName->setCurrentText(InputsVals[4]);
            XYZVoxelsHalfSize->setText(InputsVals[5]+" "+InputsVals[6]+" "+InputsVals[7]+" "+InputsVals[8]);
        }else{
            XYZVoxelsNumb->setText("3 4 5");
            ParamType->setCurrentIndex(0);
            LogVolMatName->setCurrentIndex(0);
            XYZVoxelsHalfSize->setText("10 12 15 cm");
        }

        InputsVals = getArgumentOfACommandFromText(VOXELCommands[1], 5);
        if(InputsVals.size() == 4){
            VoxContainerPos->setText(InputsVals[0]+" "+InputsVals[1]+" "+InputsVals[2]+" "+InputsVals[3]);
        }else{
            VoxContainerPos->setText("0 0 0 cm");
        }
        InputsVals = getArgumentOfACommandFromText(VOXELCommands[2], 5);
        if(InputsVals.size() == 4){
            VoxContainerRot->setText(InputsVals[0]+" "+InputsVals[1]+" "+InputsVals[2]+" "+InputsVals[3]);
        }else{
            VoxContainerRot->setText("0 0 0 degree");
        }

        framLay->addWidget(XYZVoxelsNumb, jj,ii,1,1);
        framLay->addWidget(ParamType, jj,++ii,1,1);
        framLay->addWidget(XYZVoxelsHalfSize, jj,++ii,1,1);
        framLay->addWidget(LogVolMatName, jj,++ii,1,1);
        framLay->addWidget(VoxContainerPos, jj,++ii,1,1);
        framLay->addWidget(VoxContainerRot, jj,++ii,1,1);
        framLay->addWidget(btnAddVoxelsData, jj,++ii,1,1);

    }

    if(ui->radioButtonDICOM->isChecked()){

        AllOrSpecificQComboBox = new QComboBox(); QStringList Elem=(QStringList()<<""<<"All"); AllOrSpecificQComboBox->addItems(Elem);
        NumbOfDCMMat = new QLineEdit(); NumbOfDCMMat->setText("2"); NumbOfDCMMat->setPlaceholderText("DCM Mats numb");
        LogVolMatName = new QComboBox(); for (int kk =0 ;kk < MaterialsNames.size() ; kk++) {LogVolMatName->addItem(MaterialsNames[kk]);}
        btnAddDcmMatOrdered = new QPushButton(); btnAddDcmMatOrdered->setText("Add Mat"); connect(btnAddDcmMatOrdered, SIGNAL(clicked()), this, SLOT(btnAddDcmMatOrdered_slot()));

        DcmType = new QComboBox(); QStringList DcmElem=(QStringList()<<"CT"<<"PET"); DcmType->addItems(DcmElem);
        btnChooseDcmDir = new QPushButton(); btnChooseDcmDir->setText("..."); connect(btnChooseDcmDir, SIGNAL(clicked()), this, SLOT(btnChooseDcmDir_slot()));
        DcmDataDir = new QLineEdit(); DcmDataDir->setPlaceholderText("DCM Dir");
        btnAddDcmTypeDir = new QPushButton(); btnAddDcmTypeDir->setText("Add Moda-Dir"); connect(btnAddDcmTypeDir, SIGNAL(clicked()), this, SLOT(btnAddDcmTypeDir_slot()));

        DcmCtNumber = new QLineEdit(); DcmCtNumber->setText("-4000"); DcmCtNumber->setPlaceholderText("CT Numb");
        DcmCtDensity = new QLineEdit(); DcmCtDensity->setText("0 0.5 g/cm3"); DcmCtDensity->setPlaceholderText("CT Dens unit");
        btnAddDcmCTNumDen = new QPushButton(); btnAddDcmCTNumDen->setText("Add CTNum-Dens"); connect(btnAddDcmCTNumDen, SIGNAL(clicked()), this, SLOT(btnAddDcmCTNumDen_slot()));

        ii = 0; jj++;
        framLay->addWidget(AllOrSpecificQComboBox, jj,ii,1,1);
        framLay->addWidget(NumbOfDCMMat, jj,++ii,1,1);
        framLay->addWidget(LogVolMatName, jj,++ii,1,1);
        framLay->addWidget(btnAddDcmMatOrdered, jj,++ii,1,1);

        ii = 0; jj++;;
        framLay->addWidget(DcmType, jj,ii,1,1);
        framLay->addWidget(btnChooseDcmDir, jj,++ii,1,1);
        framLay->addWidget(DcmDataDir, jj,++ii,1,1);
        framLay->addWidget(btnAddDcmTypeDir, jj,++ii,1,1);

        ii = 0; jj++;
        framLay->addWidget(DcmCtNumber, jj,ii,1,1);
        framLay->addWidget(DcmCtDensity, jj,++ii,1,1);
        framLay->addWidget(btnAddDcmCTNumDen, jj,++ii,1,1);

        FirstTimeAdd = true;
    }


    UseMaterialsAsRegionNames = new QCheckBox(); UseMaterialsAsRegionNames->setText(""); UseMaterialsAsRegionNames->setToolTip("Check to define new regions with names, or uncheck to use materials as regions names"); connect(UseMaterialsAsRegionNames, SIGNAL(clicked()), this, SLOT(UseMaterialsAsRegionNames_slot()));
    VoxRegionName = new QLineEdit(); VoxRegionName->setText("Region1"); VoxRegionName->setPlaceholderText("Region Name");
    if(!ui->radioButtonTET->isChecked()){
        VoxRegionMinMaxX = new QLineEdit(); VoxRegionMinMaxX->setText("1 2"); VoxRegionMinMaxX->setPlaceholderText("MinX MaxX");
        VoxRegionMinMaxY = new QLineEdit(); VoxRegionMinMaxY->setText("0 2"); VoxRegionMinMaxY->setPlaceholderText("MinY MinZ");
        VoxRegionMinMaxZ = new QLineEdit(); VoxRegionMinMaxZ->setText("3 3"); VoxRegionMinMaxZ->setPlaceholderText("MinZ MaxZ");
        MaterialsIDForVox = new QComboBox(); for (int kk =0 ;kk < MaterialsNames.size() ; kk++) {MaterialsIDForVox->addItem(MaterialsNames[kk]);}
    }
    VoxRegionMinMaxDensity = new QLineEdit(); VoxRegionMinMaxDensity->setText("null null g/cm3"); VoxRegionMinMaxDensity->setPlaceholderText("MinD MaxD");

    btnAddVoxelsRegionData = new QPushButton(); btnAddVoxelsRegionData->setText("Add Region Data"); connect(btnAddVoxelsRegionData, SIGNAL(clicked()), this, SLOT(btnAddVoxelsRegionData_slot()));

    ii = 0; jj++;

    framLay->addWidget(UseMaterialsAsRegionNames, jj,ii,1,1);
    framLay->addWidget(VoxRegionName, jj,++ii,1,1);
    if(!ui->radioButtonTET->isChecked()){
        framLay->addWidget(VoxRegionMinMaxX, jj,++ii,1,1);
        framLay->addWidget(VoxRegionMinMaxY, jj,++ii,1,1);
        framLay->addWidget(VoxRegionMinMaxZ, jj,++ii,1,1);
        framLay->addWidget(MaterialsIDForVox, jj,++ii,1,1);
    }
    framLay->addWidget(VoxRegionMinMaxDensity, jj,++ii,1,1);
    framLay->addWidget(btnAddVoxelsRegionData, jj,++ii,1,1);

    // fourth row

    LimitsPlan = new QComboBox(); QStringList DcmElem=(QStringList()<<"all"<<"xy"<<"xz"<<"yz"<<"regions"); LimitsPlan->addItems(DcmElem);
    connect(LimitsPlan, SIGNAL(textActivated(QString)), this, SLOT(on_LimitsPlan_textActivated(QString)));

    LimitsValues = new QLineEdit(); LimitsValues->setPlaceholderText("add min and max plan IDs (i.e. 134 146)");

    RegionToVisualizeComb = new QComboBox(); for (int kk =0 ;kk < ui->comboBoxSourceRegionProposed->count() ; kk++) {RegionToVisualizeComb->addItem(ui->comboBoxSourceRegionProposed->itemText(kk));}
    connect(RegionToVisualizeComb, SIGNAL(textActivated(QString)), this, SLOT(on_comboBoxRegionToVisualize_textActivated(QString)));

    btnAddVisualizationLimits = new QPushButton(); btnAddVisualizationLimits->setText("Add Visualization Limits"); connect(btnAddVisualizationLimits, SIGNAL(clicked()), this, SLOT(btnAddVisualizationLimits_slot()));btnAddVisualizationLimits->setToolTip("Add min slice ID and max slice ID over the selected axis for voxelized geometry.Instead, add min and max position in mm over the selected axis");

    QStringList InputsVals;
    if(ui->radioButtonTET->isChecked()){
        InputsVals = getArgumentOfACommandFromText(VOXELCommands[10], 5);
    }else{
        InputsVals = getArgumentOfACommandFromText(VOXELCommands[9], 5);
    }

    if(InputsVals.size() > 0){
        InputsVals[0].toLower();
        if(InputsVals[0] == "regions"){
            if(InputsVals.size() > 0){
            }
        }
    }

    if(InputsVals.size() > 2){
        InputsVals[0].toLower();
        if(InputsVals[0] == "regions"){
            RegionToVisualizeComb->setVisible(true);
            LimitsPlan->setCurrentIndex(4);
            QString ttt = "";
            for (int kk =0 ;kk < InputsVals.size() ; kk++){
                ttt += InputsVals[kk] + " ";
            }
            LimitsValues->setText(ttt);
        }else{
            RegionToVisualizeComb->setVisible(false);
            if(InputsVals[0] == "xy" || InputsVals[0] == "yx"){LimitsPlan->setCurrentIndex(1);}
            if(InputsVals[0] == "xz" || InputsVals[0] == "zx"){LimitsPlan->setCurrentIndex(2);}
            if(InputsVals[0] == "yz" || InputsVals[0] == "zy"){LimitsPlan->setCurrentIndex(3);}
            LimitsValues->setText(InputsVals[1] + " " + InputsVals[2]);
        }
    }else{
        RegionToVisualizeComb->setVisible(false);
        LimitsPlan->setCurrentIndex(0);
    }

    InputsVals = getArgumentOfACommandFromText(VOXELCommands[11], 5);
    if(InputsVals.size() > 0){
        if(InputsVals[0] == "yes"){
            UseMaterialsAsRegionNames->setChecked(false);
            UseMaterialsAsRegionNames_slot();
        }else{
            UseMaterialsAsRegionNames->setChecked(true);
            UseMaterialsAsRegionNames_slot();
        }
    }else{
        UseMaterialsAsRegionNames_slot();
    }

    ii = 0; jj++;
    framLay->addWidget(LimitsPlan, jj,ii,1,1);
    framLay->addWidget(LimitsValues, jj,++ii,1,1);
    framLay->addWidget(RegionToVisualizeComb, jj,++ii,1,1);
    framLay->addWidget(btnAddVisualizationLimits, jj,++ii,1,1);

    ui->fraHelpGUIInput->setLayout(framLay);
    ui->fraHelpGUIInput->setVisible(true);

}
void MainWindow::showConstructSolidAndVolumeCommandsFrame(){

    //ui->GeometryFileTextEdit->appendPlainText("ComboxLogVolSolid_slot " );

    framLay = new QGridLayout;

    ComboxSolidType = new QComboBox(); QStringList SolType=(QStringList()<<"Box"<<"Tubs"<<"CutTubs"<<"Cons"<<"Para"<<"Trd"<<"Sphere"<<"Orb"<<"Torus"<<"Ellipsoid"<<"Union"<<"Intersection"<<"Subtraction"); ComboxSolidType->addItems(SolType);
    connect(ComboxSolidType, SIGNAL(currentTextChanged(QString)), this, SLOT(ComboxSolidType_slot()));

    SolidSpecications = new QLineEdit(); SolidSpecications->setPlaceholderText("Solid Data");
    SolidName = new QLineEdit(); SolidName->setText("SolName1");
    BtnAddSol = new QPushButton(); BtnAddSol->setText("Add Solid"); connect(BtnAddSol, SIGNAL(clicked()), this, SLOT(btnAddSol_slot()));
    ComboxLogVolSolid = new QComboBox();for (int kk =0 ;kk < SolidsNames.size() ; kk++) {ComboxLogVolSolid->addItem(SolidsNames[kk]);}
    LogVolMatName = new QComboBox(); for (int kk =0 ;kk < MaterialsNames.size() ; kk++) {LogVolMatName->addItem(MaterialsNames[kk]);}
    //LogVolName = new QLineEdit(); LogVolName->setPlaceholderText("Log Vol Name");
    ComboxMotherVol = new QComboBox(); for (int kk =0 ;kk < VolsNames.size() ; kk++) {ComboxMotherVol->addItem(VolsNames[kk]);}
    SolidName->setToolTip("Add solid name to be used after"); BtnAddSol->setToolTip("Click to add the solid name to the created list of solids"); ComboxLogVolSolid->setToolTip("Add solid name to construct volume"); LogVolMatName->setToolTip("Add material that will fill the solid to construct volume"); ComboxMotherVol->setToolTip("Add mother volume where this new volume will be placed");

    PhyVolPos = new QLineEdit(); PhyVolPos->setText("0 0 0");
    PhyVolRot = new QLineEdit(); PhyVolRot->setText("0 0 0");
    PosRotUnits = new QLineEdit(); PosRotUnits->setText("cm degree");
    PhyVolName = new QLineEdit(); PhyVolName->setText("VolName1");
    ViewVolCheckBox = new QCheckBox(); ViewVolCheckBox->setText("View"); ViewVolCheckBox->setChecked(false);
    BtnAddVol = new QPushButton(); BtnAddVol->setText("Add Volume"); connect(BtnAddVol, SIGNAL(clicked()), this, SLOT(btnAddVol_slot()));
    PhyVolPos->setToolTip("Add volume position in its mother volume"); PhyVolRot->setToolTip("Add volume rotation in its mother volume"); PosRotUnits->setToolTip("Add position and rotation units"); PhyVolName->setToolTip("Click to add the volume name to the list of Volumes"); BtnAddVol->setToolTip("Add volume name to the created list of volumes");
    int ii = 0; int jj = 0;
    framLay->addWidget(ComboxSolidType, jj,ii,1,1);
    framLay->addWidget(SolidSpecications, jj,++ii,1,1);
    framLay->addWidget(SolidName, jj,++ii,1,1);
    framLay->addWidget(BtnAddSol, jj,++ii,1,1);

    ii = 0; jj++;
    framLay->addWidget(ComboxLogVolSolid, jj,ii,1,1);
    framLay->addWidget(LogVolMatName, jj,++ii,1,1);
    //framLay->addWidget(LogVolName, jj,++ii,1,1);
    framLay->addWidget(ComboxMotherVol, jj,++ii,1,1);

    ii = 0; jj++;
    framLay->addWidget(PhyVolPos, jj,ii,1,1);
    framLay->addWidget(PhyVolRot, jj,++ii,1,1);
    framLay->addWidget(PosRotUnits, jj,++ii,1,1);
    framLay->addWidget(PhyVolName, jj,++ii,1,1);
    framLay->addWidget(ViewVolCheckBox, jj,++ii,1,1);
    framLay->addWidget(BtnAddVol, jj,++ii,1,1);

    //ii = 0, jj++;
    //framLay->addWidget(BtnAddVol, jj,ii,1,1);

    ui->fraHelpGUIInput->setLayout(framLay);
    ui->fraHelpGUIInput->setVisible(true);

}
void MainWindow::showConstructSTLCommandsFrame(){

    //ui->GeometryFileTextEdit->appendPlainText("ComboxLogVolSolid_slot " );

    framLay = new QGridLayout;

    BtnAddSolPath = new QPushButton(); BtnAddSolPath->setText("Add Solid"); connect(BtnAddSolPath, SIGNAL(clicked()), this, SLOT(BtnAddSolPath_slot()));
    PhyVolName = new QLineEdit(); PhyVolName->setText(DoseCalcsCore_build_dir_path+"/Scripts/solid.stl");
    LogVolMatName = new QComboBox(); for (int kk =0 ;kk < MaterialsNames.size() ; kk++) {LogVolMatName->addItem(MaterialsNames[kk]);}
    ComboxMotherVol = new QComboBox(); for (int kk = 0 ;kk < VolsNames.size() ; kk++) {ComboxMotherVol->addItem(VolsNames[kk]);}
    PhyVolPos = new QLineEdit(); PhyVolPos->setText("0 0 0");
    PhyVolRot = new QLineEdit(); PhyVolRot->setText("0 0 0");
    PosRotUnits = new QLineEdit(); PosRotUnits->setText("cm degree");
    BtnAddVol = new QPushButton(); BtnAddVol->setText("Add Volume"); connect(BtnAddVol, SIGNAL(clicked()), this, SLOT(btnAddVol_slot()));
    ViewVolCheckBox = new QCheckBox(); ViewVolCheckBox->setText("View"); ViewVolCheckBox->setChecked(false);
    BtnAddSolPath->setToolTip("Add solid file path that contains solid data (.stl or .ast)"); PhyVolPos->setToolTip("Add volume position in its mother volume"); PhyVolRot->setToolTip("Add volume rotation in its mother volume"); PosRotUnits->setToolTip("Add position and rotation units"); PhyVolName->setToolTip("This solid name will be the the volume name"); BtnAddVol->setToolTip("Add volume name to the created list of volumes"); LogVolMatName->setToolTip("Add material that will fill the solid to construct volume");

    int ii = 0; int jj = 0;
    framLay->addWidget(BtnAddSolPath, jj,ii,1,1);
    framLay->addWidget(PhyVolName, jj,++ii,1,1);
    framLay->addWidget(LogVolMatName, jj,++ii,1,1);

    ii = 0; jj++;
    framLay->addWidget(ComboxMotherVol, jj,ii,1,1);
    framLay->addWidget(PhyVolPos, jj,++ii,1,1);
    framLay->addWidget(PhyVolRot, jj,++ii,1,1);
    framLay->addWidget(PosRotUnits, jj,++ii,1,1);
    framLay->addWidget(ViewVolCheckBox, jj,++ii,1,1);
    framLay->addWidget(BtnAddVol, jj,++ii,1,1);

    ui->fraHelpGUIInput->setLayout(framLay);
    ui->fraHelpGUIInput->setVisible(true);

}
void MainWindow::showConstructVolumeGDMLTEXTCPPCommandsFrame(){

    framLay = new QGridLayout;

    BtnAddVolPath = new QPushButton(); BtnAddVolPath->setText("Add Vol file"); connect(BtnAddVolPath, SIGNAL(clicked()), this, SLOT(BtnAddVolPath_slot()));
    BtnEditVolFile = new QPushButton(); BtnEditVolFile->setText("Edit Vol file"); connect(BtnEditVolFile, SIGNAL(clicked()), this, SLOT(BtnEditVolFile_slot()));

    PhyVolName = new QLineEdit();

    if(ui->radioButtonTEXT->isChecked()){
        PhyVolName->setText(DoseCalcsCore_build_dir_path+"/Scripts/LogVol.geom");
    }else if(ui->radioButtonCpp->isChecked()){
        PhyVolName->setText(DoseCalcsCore_build_dir_path+"/Scripts/LogVol.c++");
    }else{
        PhyVolName->setText(DoseCalcsCore_build_dir_path+"/Scripts/LogVol.gdml");
    }


    ComboxMotherVol = new QComboBox(); for (int kk =0 ;kk < VolsNames.size() ; kk++) {ComboxMotherVol->addItem(VolsNames[kk]);}
    PhyVolPos = new QLineEdit(); PhyVolPos->setText("0 0 0");
    PhyVolRot = new QLineEdit(); PhyVolRot->setText("0 0 0");
    PosRotUnits = new QLineEdit(); PosRotUnits->setText("cm degree");
    BtnAddVol = new QPushButton(); BtnAddVol->setText("Add Volume"); connect(BtnAddVol, SIGNAL(clicked()), this, SLOT(btnAddVol_slot()));
    ViewVolCheckBox = new QCheckBox(); ViewVolCheckBox->setText("View"); ViewVolCheckBox->setChecked(false);

    BtnEditVolFile->setToolTip("Edit volume geometry file that contains volume data (.gdml, .geom, .c++)");BtnAddVolPath->setToolTip("Add volume path that contains volume data (.gdml, .geom, .c++)"); PhyVolPos->setToolTip("Add volume position in its mother volume"); PhyVolRot->setToolTip("Add volume rotation in its mother volume"); PosRotUnits->setToolTip("Add position and rotation units"); PhyVolName->setToolTip("Click to add the volume name to the list of Volumes"); BtnAddVol->setToolTip("Add volume name to the created list of volumes");

    int ii = 0; int jj = 0;
    framLay->addWidget(BtnAddVolPath, jj,ii,1,1);
    framLay->addWidget(PhyVolName, jj,++ii,1,1);
    framLay->addWidget(BtnEditVolFile, jj,++ii,1,1);

    ii = 0; jj = jj+1;
    framLay->addWidget(ComboxMotherVol, jj,ii,1,1);
    framLay->addWidget(PhyVolPos, jj,++ii,1,1);
    framLay->addWidget(PhyVolRot, jj,++ii,1,1);
    framLay->addWidget(PosRotUnits, jj,++ii,1,1);
    framLay->addWidget(ViewVolCheckBox, jj,++ii,1,1);
    framLay->addWidget(BtnAddVol, jj,++ii,1,1);

    ui->fraHelpGUIInput->setLayout(framLay);
    ui->fraHelpGUIInput->setVisible(true);
}

// slots of dynamic widgets
void MainWindow::btnSaveElem_slot(){

    QString ElemCommandToAdd = "";

    bool isIn = false;

    for ( int df = 0 ; df < ElementsNames.size(); df++  ){
        if(ElementsNames[df] == MatElemComb->currentText() ||
                ElementsNames[df] == ElementsSymbolName[MatElemComb->currentText()] ||
                ElementsNames[df] == ElementsSymbolSym[MatElemComb->currentText()]){
            isIn = true;
        }
    }
    if(isIn == true ){
    }
    else{

        ElementsNames.push_back(MatElemComb->currentText());
        ElemComb->addItem(MatElemComb->currentText());

        ElemCommandToAdd = MaterialCommands[0]
                + " " + ElementsSymbolZ[MatElemComb->currentText()]
                + " " + ElementsSymbolA[MatElemComb->currentText()]
                + " " + ElementsSymbolName[MatElemComb->currentText()];

        ui->GeometryFileTextEdit->appendPlainText(ElemCommandToAdd );
    }


}
void MainWindow::btnAddElem_slot(){

    if(MatCommandToAddElem == ""){
        MatCommandToAddElem = "C";
        ui->GeometryFileTextEdit->appendPlainText( "\n" + MaterialCommands[2] + " " + ElemComb->currentText() + " " + EleFraOrNumInMat->text() + " " );
    }
    else{
        QString tx = ui->GeometryFileTextEdit->toPlainText();
        ui->GeometryFileTextEdit->setPlainText(tx + " " + ElemComb->currentText() + " " + EleFraOrNumInMat->text() + " ");
    }
    ui->tabWidget->setCurrentIndex(0);
}
void MainWindow::btnAddMat_slot(){

    QString MatCommandToAdd = "";
    MatCommandToAddElem = "";
    bool isIn = false;
    for ( int df = 0 ; df < MaterialsNames.size(); df++  ){
        if(MaterialsNames[df] == MaterialName->text() ){
            isIn = true;
        }
    }

    if(isIn == true ){ }
    else{
        MaterialsNames.push_back(MaterialName->text());
        MaterialRegionsNames.push_back(MaterialName->text());
        MaterialsNameIDs[MaterialName->text()] = MaterialID->text();

        QString frOrNum = "Numb";
        if( NumFraCheckBox->isChecked() ){
            frOrNum = "frac";
        }
        MatCommandToAdd += "\n" + MaterialCommands[1]
                + " " + MaterialName->text()
                + " " + MaterialID->text()
                + " " + MaterialEleNumber->text()
                + " " + MaterialDensity->text()
                + " " + frOrNum;

        ui->GeometryFileTextEdit->appendPlainText(MatCommandToAdd);
    }
    ui->tabWidget->setCurrentIndex(0);
}
void MainWindow::btnAddNistMat_slot(){

    bool isIn = false;
    QString MatCommandToAdd = "";

    for ( int df = 0 ; df < MaterialsNames.size(); df++  ){
        if(MaterialsNames[df] == NistMatComb->currentText() ){
            isIn = true;
        }
    }

    if(isIn == true ){}
    else{
        MaterialsNames.push_back(NistMatComb->currentText());
        MaterialRegionsNames.push_back(NistMatComb->currentText());
        MaterialsNameIDs[MaterialName->text()] = NistMaterialID->text();
        ui->GeometryFileTextEdit->appendPlainText("\n"+MaterialCommands[4] + " " + NistMatComb->currentText()+ " " + NistMaterialID->text());
    }
    ui->tabWidget->setCurrentIndex(0);
}
void MainWindow::CheckBoxFracNum_slot(){

    if(NumFraCheckBox->isChecked()){
        EleFraOrNumInMat->setPlaceholderText("Fraction");
    }else{
        EleFraOrNumInMat->setPlaceholderText("Number");
    }
}
void MainWindow::on_LimitsPlan_textActivated(QString arg)
{
    if(arg == "regions"){
        RegionToVisualizeComb->setVisible(true);
    }else{
        RegionToVisualizeComb->setVisible(false);
    }
}
void MainWindow::runTerminalCommandSlot()
{
    BashCommandsForExecuting = "#! /bin/bash \n "
                               "qstat -f\n"
                               //"ls \n"
                               "bash \n ";
    fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);
    ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);

}
void MainWindow::on_ListNodes_textActivated(QString arg)
{
    bool IsIn = false;
    QStringList args = Textnodes->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
    for(int dd=0; dd < args.size();dd++){ if(args[dd] == ListNodes->currentText()){ IsIn = true;}}
    if(IsIn == false){
        if(ListNodes->currentText() != ""){
            QString nn = Textnodes->text();
            nn.replace(" ","");
            if(nn == ""){
                Textnodes->setText(ListNodes->currentText());
            }
            else{
                Textnodes->setText(nn + "+" + ListNodes->currentText());
            }
        }
    }
}
void MainWindow::on_comboBoxRegionToVisualize_textActivated(QString arg)
{
    bool IsIn = false;
    QStringList args = LimitsValues->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
    for(int dd=0; dd < args.size();dd++){ if(args[dd] == RegionToVisualizeComb->currentText()){ IsIn = true;}}
    if(IsIn == false){

        if(RegionToVisualizeComb->currentText() != ""){
            QString nn = LimitsValues->text();
            LimitsValues->setText(nn + " " + RegionToVisualizeComb->currentText());
        }
    }
}
void MainWindow::btnAddVoxelsRegionData_slot(){

    if(ui->radioButtonTET->isChecked()){
        RegionDataCommands = VOXELCommands[4]
                + " " + VoxRegionName->text()
                + " " + VoxRegionMinMaxDensity->text() +"\n" ;
    }else{
        RegionDataCommands = VOXELCommands[4]
                + " " + VoxRegionName->text()
                + " " + VoxRegionMinMaxX->text()
                + " " + VoxRegionMinMaxY->text()
                + " " + VoxRegionMinMaxZ->text()
                + " " + MaterialsNameIDs[MaterialsIDForVox->currentText()]
                + " " + VoxRegionMinMaxDensity->text() +"\n" ;
    }

    DefinedRegionsNames.push_back(VoxRegionName->text());

    ui->tabWidget->setCurrentIndex(0);
    ui->GeometryFileTextEdit->appendPlainText(RegionDataCommands);

}
void MainWindow::UseMaterialsAsRegionNames_slot(){

    ui->tabWidget->setCurrentIndex(0);
    if(!UseMaterialsAsRegionNames->isChecked()){
        ui->GeometryFileTextEdit->setPlainText(changeArgumentOfACommandFromText(ui->GeometryFileTextEdit->toPlainText(),VOXELCommands[11], VOXELCommands[11] + " yes "));
        VoxRegionName->setEnabled(false);
        if(!ui->radioButtonTET->isChecked()){
            VoxRegionMinMaxX->setEnabled(false);
            VoxRegionMinMaxY->setEnabled(false);
            VoxRegionMinMaxZ->setEnabled(false);
            MaterialsIDForVox->setEnabled(false);
        }
        VoxRegionMinMaxDensity->setEnabled(false);
        btnAddVoxelsRegionData->setEnabled(false);
    }
    else{
        VoxRegionName->setEnabled(true);
        if(!ui->radioButtonTET->isChecked()){
            VoxRegionMinMaxX->setEnabled(true);
            VoxRegionMinMaxY->setEnabled(true);
            VoxRegionMinMaxZ->setEnabled(true);
            MaterialsIDForVox->setEnabled(true);
        }
        VoxRegionMinMaxDensity->setEnabled(true);
        btnAddVoxelsRegionData->setEnabled(true);

        ui->GeometryFileTextEdit->setPlainText(changeArgumentOfACommandFromText(ui->GeometryFileTextEdit->toPlainText(),VOXELCommands[11], VOXELCommands[11] + " no "));
    }
}

void MainWindow::btnAddVoxelsData_slot(){

    if(LogVolMatName->count() != 0){
        ui->GeometryFileTextEdit->setPlainText(changeArgumentOfACommandFromText(ui->GeometryFileTextEdit->toPlainText(),VOXELCommands[0], VOXELCommands[0] + " " + XYZVoxelsNumb->text() + " " + ValuesOfInputs[ParamType->currentText()] + " " + LogVolMatName->currentText() + " " + XYZVoxelsHalfSize->text()));
    }else{
        QMessageBox::information(this, tr("Voxels default material"),"Cannot create voxels data command. Build a material to be added as a default material to voxels");
    }
    ui->GeometryFileTextEdit->setPlainText(changeArgumentOfACommandFromText(ui->GeometryFileTextEdit->toPlainText(),VOXELCommands[1], VOXELCommands[1] + " " + VoxContainerPos->text()));
    ui->GeometryFileTextEdit->setPlainText(changeArgumentOfACommandFromText(ui->GeometryFileTextEdit->toPlainText(),VOXELCommands[2], VOXELCommands[2] + " " + VoxContainerRot->text()));

    ui->tabWidget->setCurrentIndex(0);
    //ui->GeometryFileTextEdit->appendPlainText(VoxString);
}
void MainWindow::btnAddDcmMatOrdered_slot(){

    if(AllOrSpecificQComboBox->currentText() == "All"){

        DcmCtDensCommands = VOXELCommands[6] + " " + QString::number(MaterialsNames.size());

        for (int kk =0 ;kk < MaterialsNames.size() ; kk++)
        {
            DcmCtDensCommands += " " + MaterialsNames[kk];
        }
        ui->GeometryFileTextEdit->appendPlainText(DcmCtDensCommands + " \n");

    }else{
        if(FirstTimeAdd){
            DcmCtDensCommands = VOXELCommands[6] + " " + NumbOfDCMMat->text() + " " + LogVolMatName->currentText() ;
            ui->GeometryFileTextEdit->appendPlainText(DcmCtDensCommands);
            FirstTimeAdd = false;
        }else{
            QString tx = ui->GeometryFileTextEdit->toPlainText();
            ui->GeometryFileTextEdit->setPlainText( tx + " " + LogVolMatName->currentText());
        }
    }
    ui->tabWidget->setCurrentIndex(0);
}
void MainWindow::btnAddDcmCTNumDen_slot(){
    DcmCtDensCommands  = VOXELCommands[5]
            + " " + DcmCtNumber->text()
            + " " + DcmCtDensity->text();
    ui->tabWidget->setCurrentIndex(0);
    ui->GeometryFileTextEdit->appendPlainText(DcmCtDensCommands);
}
void MainWindow::btnAddDcmTypeDir_slot(){
    DcmTypeDirCommands = VOXELCommands[7]
            + " " + DcmType->currentText()
            + " " + DcmDataDir->text();
    ui->tabWidget->setCurrentIndex(0);
    ui->GeometryFileTextEdit->appendPlainText(DcmTypeDirCommands);
}
void MainWindow::btnChooseDcmDir_slot(){
    QString chosen_DirPath = QFileDialog::getExistingDirectory(0, ("Choose modality files"), DoseCalcsCore_build_dir_path+"/"+ScriptDirectoryName) ;
    if(chosen_DirPath.isEmpty() || chosen_DirPath.isNull()){

    }else{
        DcmDataDir->setText(chosen_DirPath);
    }
}
void MainWindow::btnAddVisualizationLimits_slot(){
    //  ui->GeometryFileTextEdit->appendPlainText(VOXELCommands[9] + " " + LimitsPlan->currentText() + " " + LimitsValues->text() + "\n" );
    ui->tabWidget->setCurrentIndex(0);

    if(ui->radioButtonTET->isChecked()){
        ui->GeometryFileTextEdit->setPlainText(changeArgumentOfACommandFromText(ui->GeometryFileTextEdit->toPlainText(),VOXELCommands[10], VOXELCommands[10] + " " + LimitsPlan->currentText() + " " + LimitsValues->text()));
    }else{
        ui->GeometryFileTextEdit->setPlainText(changeArgumentOfACommandFromText(ui->GeometryFileTextEdit->toPlainText(),VOXELCommands[9], VOXELCommands[9] + " " + LimitsPlan->currentText() + " " + LimitsValues->text()));
    }

}
void MainWindow::BtnAddSolPath_slot(){

    QString chosen_DirPath = QFileDialog::getOpenFileName( this, tr("Choose Solid .stl file"), DoseCalcsCore_build_dir_path+"/"+ScriptDirectoryName, tr("All files (*.*)") );

    if(chosen_DirPath.isEmpty() || chosen_DirPath.isNull()){

    }else{
        SolidName->setText(chosen_DirPath);
    }
}
void MainWindow::BtnAddVolPath_slot(){

    QString chosen_DirPath = QFileDialog::getOpenFileName( this, tr("Choose Volume file"), DoseCalcsCore_build_dir_path+"/"+ScriptDirectoryName, tr("All files (*.*)") );

    if(chosen_DirPath.isEmpty() || chosen_DirPath.isNull()){

    }else{
        PhyVolName->setText(chosen_DirPath);
    }
}
void MainWindow::BtnEditVolFile_slot(){

    on_pushButtonEditGeomFile_clicked(); // if we click edit again, save before.

    QString ext = QString::fromLocal8Bit(getFileExt(PhyVolName->text().toStdString()).c_str());
    QString flnm = QString::fromLocal8Bit(getFileNameFromPath(PhyVolName->text().toStdString()).c_str())+"."+QString::fromLocal8Bit(getFileExt(ui->PhantomWorldHalfSizeslineEdit->text().toStdString()).c_str());

    if(ext == "gdml" || ext == "geom" || ext == "stl" || ext == "ast"){
        EditFlag = 7;
        ui->tabWidget->setTabText(0,flnm);
        ui->GeometryFileTextEdit->clear();
        showResultsOutput("", 4);
        ui->GeometryFileTextEdit->setPlainText(fileManagerObject->ReadTextFromFileInOneString(PhyVolName->text()));
        ui->tabWidget->setCurrentIndex(0);
    }
    else if(ext == "c++" || ext == "cpp"  || ext == "cc"){
        EditFlag = 7;
        ui->tabWidget->setTabText(0,"G4TCPPGeometryFormat.cc");
        ui->GeometryFileTextEdit->clear();
        showResultsOutput("", 4);
        ui->GeometryFileTextEdit->setPlainText(fileManagerObject->ReadTextFromFileInOneString(DoseCalcsCore_source_dir_path+"/src/G4TCPPGeometryFormat.cc"));
        ui->tabWidget->setCurrentIndex(0);
    }
}
void MainWindow::btnAddSol_slot(){

    bool isIn = false;

    for ( int df = 0 ; df < SolidsNames.size(); df++  ){
        if(SolidsNames[df] == SolidName->text() ){
            isIn = true;
        }
    }
    if(isIn == true ){
    }
    else{

        SolidsNames.push_back(SolidName->text());
        ComboxLogVolSolid->addItem(SolidName->text());

        QString SolidCommandToAdd = GeometryCommands[1]
                + " " + ComboxSolidType->currentText()
                + " " + SolidName->text()
                + " " + SolidSpecications->text();

        ui->tabWidget->setCurrentIndex(0);

        ui->GeometryFileTextEdit->appendPlainText(SolidCommandToAdd );
        //CONSCommandsString += SolidCommandToAdd + "\n";

    }
}
void MainWindow::btnAddVol_slot(){

    QString PhysVoNm = PhyVolName->text();
    QString fe =  QString::fromLocal8Bit(getFileExt(PhyVolName->text().toStdString()).c_str());

    if(fe == "gdml" || fe == "geom" || fe == "c++" || fe == "cpp" || fe == "stl" || fe == "ast"){
        PhysVoNm = QString::fromLocal8Bit(getFileNameFromPath(PhysVoNm.toStdString()).c_str());
    }

    bool isIn = false;
    for ( int df = 0 ; df < VolsNames.size(); df++  ){
        if(VolsNames[df] == PhysVoNm ){
            isIn = true;
        }
    }

    if(isIn == true ){

        if(ViewVolCheckBox->isChecked()){
            QString TransparentVol = "";
            for ( int df = 0 ; df < VolsNames.size(); df++  ){
                if(VolsNames[df] == PhysVoNm ){}
                else{ TransparentVol += VolsNames[df] + " "; }
            }

            NotForcedVolumesCommand = GeometryCommands[3] + " " + TransparentVol;
            on_actionVisualization_triggered();
        }else{
            NotForcedVolumesCommand = GeometryCommands[3] + " World ";
        }

    }
    else{
        ComboxMotherVol->addItem(PhysVoNm);

        QString VolumeCommandToAdd;
        if(fe == "stl" || fe == "ast"){
            if(LogVolMatName->count() != 0){

                VolumeCommandToAdd = GeometryCommands[2]
                        + " " + PhyVolName->text()
                        + " " + LogVolMatName->currentText()
                        + " " + ComboxMotherVol->currentText()
                        + " " + PhyVolPos->text()
                        + " " + PhyVolRot->text()
                        + " " + PosRotUnits->text();
            }else{
                QMessageBox::information(this, tr("No material found"),"Cannot create volume command. Build a material first to fill the STL solid");
            }

        }
        else if(fe == "gdml" || fe == "geom" || fe == "c++" || fe == "cpp"){
            VolumeCommandToAdd = GeometryCommands[2]
                    + " " + PhyVolName->text()
                    + " " + ComboxMotherVol->currentText()
                    + " " + PhyVolPos->text()
                    + " " + PhyVolRot->text()
                    + " " + PosRotUnits->text();
        }
        else {
            VolumeCommandToAdd = GeometryCommands[2]
                    + " " + PhyVolName->text()
                    + " " + ComboxLogVolSolid->currentText()
                    + " " + LogVolMatName->currentText()
                    + " " + ComboxMotherVol->currentText()
                    + " " + PhyVolPos->text()
                    + " " + PhyVolRot->text()
                    + " " + PosRotUnits->text();
        }

        if(ViewVolCheckBox->isChecked()){
            QString TransparentVol = "";
            for ( int df = 0 ; df < VolsNames.size(); df++  ){
                if(VolsNames[df] == PhysVoNm ){}
                else{ TransparentVol += VolsNames[df] + " "; }
            }

            NotForcedVolumesCommand = GeometryCommands[3] + " " + TransparentVol;
            on_actionVisualization_triggered();
        }else{
            NotForcedVolumesCommand = GeometryCommands[3] + " World ";
        }

        ui->GeometryFileTextEdit->appendPlainText(VolumeCommandToAdd );

        ui->tabWidget->setCurrentIndex(0);

        //CONSCommandsString += VolumeCommandToAdd + "\n";

    }
}
void MainWindow::ComboxSolidType_slot(){

    if(ComboxSolidType->currentText() == "Box"){
        SolidSpecications->setText("100 100 100 cm");
        SolidSpecications->setToolTip("The " + ComboxSolidType->currentText() +" parameters are: hx hy hz LengthUnit AngleUnit");
        SolidSpecications->setPlaceholderText("hx hy hz cm");
    }
    else if(ComboxSolidType->currentText() == "CutTubs"){
        SolidSpecications->setText("  cm degree");
        SolidSpecications->setToolTip("The " + ComboxSolidType->currentText() +" parameters are: Rmin Rmax Dz SPhi DPhi LowNorm HighNorm LengthUnit AngleUnit");
        SolidSpecications->setPlaceholderText(" Rmin Rmax Dz SPhi DPhi LowNorm HighNorm cm degree ");
    }
    else if(ComboxSolidType->currentText() == "Tubs"){
        SolidSpecications->setText("  cm degree");
        SolidSpecications->setToolTip("The " + ComboxSolidType->currentText() +" parameters are: Rmin Rmax Dz SPhi DPhi LengthUnit AngleUnit");
        SolidSpecications->setPlaceholderText(" Rmin Rmax Dz SPhi DPhi cm degree ");
    }
    else if(ComboxSolidType->currentText() == "Cons"){
        SolidSpecications->setText("  cm degree");
        SolidSpecications->setToolTip("The " + ComboxSolidType->currentText() +" parameters are: Rmin1 Rmax1 Rmin2 Rmax2 Dz SPhi DPhi LengthUnit AngleUnit");
        SolidSpecications->setPlaceholderText("  ");
    }
    else if(ComboxSolidType->currentText() == "Para"){
        SolidSpecications->setText("  cm degree");
        SolidSpecications->setToolTip("The " + ComboxSolidType->currentText() +" parameters are: Dx Dy Dz Alpha Theta0 Phi0 LengthUnit AngleUnit");
        SolidSpecications->setPlaceholderText(" Dx Dy Dz Alpha Theta0 Phi0 cm degree ");
    }
    else if(ComboxSolidType->currentText() == "Trd"){
        SolidSpecications->setText("  cm degree");
        SolidSpecications->setToolTip("The " + ComboxSolidType->currentText() +" parameters are: Dx1 Dx2 Dy1 Dy2 Dz LengthUnit AngleUnit");
        SolidSpecications->setPlaceholderText(" Dx1 Dx2 Dy1 Dy2 Dz cm degree ");
    }
    else if(ComboxSolidType->currentText() == "Sphere"){
        SolidSpecications->setText("  cm degree");
        SolidSpecications->setToolTip("The " + ComboxSolidType->currentText() +" parameters are: Rmin Rmax SPhi DPhi STheta DTheta LengthUnit AngleUnit");
        SolidSpecications->setPlaceholderText(" Rmin Rmax SPhi DPhi STheta DTheta cm degree ");
    }
    else if(ComboxSolidType->currentText() == "Orb"){
        SolidSpecications->setText("  cm degree");
        SolidSpecications->setToolTip("The " + ComboxSolidType->currentText() +" parameters are: Rmax LengthUnit");
        SolidSpecications->setPlaceholderText(" Rmax cm ");
    }
    else if(ComboxSolidType->currentText() == "Torus"){
        SolidSpecications->setText("  cm degree");
        SolidSpecications->setToolTip("The " + ComboxSolidType->currentText() +" parameters are: Rmin Rmax Rtor SPhi DPhi LengthUnit AngleUnit");
        SolidSpecications->setPlaceholderText(" Rmin Rmax Rtor SPhi DPhi cm degree ");
    }
    else if(ComboxSolidType->currentText() == "Ellipsoid"){
        SolidSpecications->setText("  cm degree");
        SolidSpecications->setToolTip("The " + ComboxSolidType->currentText() +" parameters are: xSemiAxis ySemiAxis zSemiAxis zBottomCut zBottomCut LengthUnit AngleUnit");
        SolidSpecications->setPlaceholderText(" xSemiAxis ySemiAxis zSemiAxis zBottomCut zBottomCut cm degree ");
    }
    else if(ComboxSolidType->currentText() == "Union" || ComboxSolidType->currentText() == "Intersection" || ComboxSolidType->currentText() == "Subtraction"){
        SolidSpecications->setText("  cm degree");
        SolidSpecications->setToolTip("The " + ComboxSolidType->currentText() +" parameters are: Solid1 Solid2 XPos YPos ZPos XRot YRot ZRot LengthUnit AngleUnit");
        SolidSpecications->setPlaceholderText(" Solid1 Solid2 XPos YPos ZPos XRot YRot ZRot cm degree ");
    }

}
void MainWindow::on_groupBoxRootGraphs_clicked()
{
    ui->groupBoxRootGraphs->hide();
}

// Analysis box inputs handle
void MainWindow::on_btnPositionsFile_clicked()
{
    PositionsDataFilePath = QFileDialog::getOpenFileName(
                this,
                tr("Open file"),
                DoseCalcsCore_build_dir_path, //pathBuildApp,
                "All files (*.*);;Text files (*.txt)"
                );

    //QMessageBox::information(this, tr("File Path is : "), PositionsDataFilePath);

    if(PositionsDataFilePath != NULL && PositionsDataFilePath != "" ){

        showResultsOutput("Chosen file is " + PositionsDataFilePath , 4);
        ui->AnalysisLineEditPositionsFile->setText(PositionsDataFilePath);

    }else{
        //QMessageBox::information(this, tr("File..!!"),"No file is choosen");
        showResultsOutput("No file is chosen, please choose an input file to fill the components automatically ! ", 0);
    }
}
void MainWindow::on_btnEnergiessFile_clicked()
{
    EnergiesDataFilePath = QFileDialog::getOpenFileName(
                this,
                tr("Open file"),
                DoseCalcsCore_build_dir_path, //pathBuildApp,
                "All files (*.*);;Text files (*.txt)"
                );

    //QMessageBox::information(this, tr("File Path is : "), EnergiesDataFilePath);

    if(EnergiesDataFilePath != NULL && EnergiesDataFilePath != "" ){

        showResultsOutput("Chosen file is " + EnergiesDataFilePath , 4);
        ui->AnalysisLineEditEnergiesFile->setText(EnergiesDataFilePath);

    }else{
        //QMessageBox::information(this, tr("File..!!"),"No file is choosen");
        showResultsOutput("No file is chosen, please choose an input file to fill the components automatically ! ", 0);
    }
}
void MainWindow::on_btnMomDirsFile_clicked()
{
    MomDirsDataFilePath = QFileDialog::getOpenFileName(
                this,
                tr("Open file"),
                DoseCalcsCore_build_dir_path, //pathBuildApp,
                "All files (*.*);;Text files (*.txt)"
                );

    //QMessageBox::information(this, tr("File Path is : "), MomDirsDataFilePath);

    if(MomDirsDataFilePath != NULL && MomDirsDataFilePath != "" ){

        showResultsOutput("Chosen file is " + MomDirsDataFilePath , 4);
        ui->AnalysisLineEditMomDirsFile->setText(MomDirsDataFilePath);

    }else{
        //QMessageBox::information(this, tr("File..!!"),"No file is choosen");
        showResultsOutput("No file is chosen, please choose an input file to fill the components automatically ! ", 0);
    }
}
void MainWindow::on_btnAnalysiInputFileFile_clicked()
{

    showResultsOutput("---------- Filling analysis inputs by default values ", 4);

    ui->AnalysisComboBoxGraphData->setCurrentText("Result");
    ui->AnalysisComboBoxGraphsType->setCurrentText("Self_Cross");
    ui->AnalysisLineEdit_RefName->setText("ICRP");
    ui->AnalysisLineEditRefFile->setText(UserCurrentResultsDirPath+"/ICRPValues");

    ui->comboBoxRelDiff->setCurrentText("RA");
    ui->checkBoxRelSDev->setChecked(true) ;
    ui->checkBoxRegionParameter->setChecked(true) ;
    ui->AnalysisComboBoxRegionVariable->setCurrentText("Mass");

    ui->AnalysisLineEditPositionsFile->setText(DoseCalcsCore_build_dir_path+"/EventsData/Pos_Liver_Voxels_100000000_0.bin");
    ui->AnalysisLineEditEnergiesFile->setText(DoseCalcsCore_build_dir_path+"/EventsData/Ene_Mono_0.01_100000000_0.bin");
    ui->AnalysisLineEditMomDirsFile->setText(DoseCalcsCore_build_dir_path+"/EventsData/MomDir_Isotropic_100000000_0.bin");
    ui->checkBoxEventsDataHisto->setChecked(true);

    ui->AnalysisComboBoxSliceFor2DGraph->setCurrentText("XY");
    ui->AnalysisComboBoxBeamAxis->setCurrentText("Z");
    ui->AnalysisLineEdit_SliceID->setText("5");

    ui->checkBoxPrintTitle->setChecked(true);
    ui->checkBoxUseLogE->setChecked(true);
    ui->checkBoxUseLogVar->setChecked(true);
    ui->checkBoxUseGrid->setChecked(true);

    ui->AnalysisLegendPoscomboBox->setCurrentText("RightTop");
    ui->lineEditLegendXYWidth->setText("0.15 0.23");

    ui->checkBoxAddErrorBar->setChecked(false);

    ui->AnalysisComboBoxGraphsExt->setCurrentText(".root");

    updateApplicationTabs();
    /*
        QVector< QPair<QString,QString>> Commlines1 = fileManagerObject->ReadTextFromFileInQStringList(DoseCalcsCore_build_dir_path+"/"+MacroFileName);

        for ( int ds = 0 ; ds < Commlines1.size(); ds++  ){

            if(Commlines1[ds].first == "GraphsData" ){ui->AnalysisComboBoxGraphData->setCurrentText(Commlines1[ds].second) ; }
            if(Commlines1[ds].first == "CompareType" ){ui->AnalysisComboBoxGraphsType->setCurrentText(Commlines1[ds].second) ; }
            if(Commlines1[ds].first == "RefFilePath" ){ui->AnalysisLineEditRefFile->setText(Commlines1[ds].second) ; }
            if(Commlines1[ds].first == "RefName" ){ui->AnalysisLineEdit_RefName->setText(Commlines1[ds].second) ; }


            if(Commlines1[ds].first == "RegionVariableName" ){ui->AnalysisComboBoxRegionVariable->setCurrentText(Commlines1[ds].second) ; }

            if(Commlines1[ds].first == "GenerateRegionsVariableGraph" ){
                if(Commlines1[ds].second == "yes" ){
                    ui->checkBoxRegionParameter->setChecked(true);
                }
            }
            if(Commlines1[ds].first == "GenerateRelativeSDevGraph"){
                if(Commlines1[ds].second == "yes" ){
                    ui->checkBoxRelSDev->setChecked(true);
                }
            }
            if(Commlines1[ds].first == "GenerateRelativeErrGraph" ){
                if(Commlines1[ds].second == "yes" ){
                    ui->checkBoxRelErr->setChecked(true);
                }
            }

            if(Commlines1[ds].first == "SliceFor2DGraph" ){ui->AnalysisComboBoxSliceFor2DGraph->setCurrentText(Commlines1[ds].second) ; }
            if(Commlines1[ds].first == "BeamAxis" ){ui->AnalysisComboBoxBeamAxis->setCurrentText(Commlines1[ds].second) ; }
            if(Commlines1[ds].first == "SliceID" ){ui->AnalysisLineEdit_SliceID->setText(Commlines1[ds].second) ; }

            if(Commlines1[ds].first == "PositionDataFile" ){ui->AnalysisLineEditPositionsFile->setText(Commlines1[ds].second) ; }
            if(Commlines1[ds].first == "EnergyDataFile" ){ui->AnalysisLineEditEnergiesFile->setText(Commlines1[ds].second) ; }
            if(Commlines1[ds].first == "MomDirDataFile" ){ui->AnalysisLineEditMomDirsFile->setText(Commlines1[ds].second) ; }

            if(Commlines1[ds].first == "EventsDataHistograms" ){
                if(Commlines1[ds].second == "yes" ){
                    ui->checkBoxEventsDataHisto->setChecked(true);
                }
            }

            if(Commlines1[ds].first == "UseLogE" ){
                if(Commlines1[ds].second == "yes" ){
                    ui->checkBoxUseLogE->setChecked(true);
                }
            }
            if(Commlines1[ds].first == "UseLogVariable" ){
                if(Commlines1[ds].second == "yes" ){
                    ui->checkBoxUseLogVar->setChecked(true);
                }
            }
            if(Commlines1[ds].first == "UseGridXY" ){
                if(Commlines1[ds].second == "yes" ){
                    ui->checkBoxUseGrid->setChecked(true);
                }
            }
            if(Commlines1[ds].first == "PrintTitle" ){
                if(Commlines1[ds].second == "yes" ){
                    ui->checkBoxPrintTitle->setChecked(true);
                }
            }

            if(Commlines1[ds].first == "LegendPos" ){ui->AnalysisLegendPoscomboBox->setCurrentText(Commlines1[ds].second) ; }

            if(Commlines1[ds].first == "LegendXWidth" ){LegWidth += Commlines1[ds].second ; ui->lineEditLegendXYWidth->setText(LegWidth) ; }
            if(Commlines1[ds].first == "LegendYHeight" ){LegWidth += Commlines1[ds].second ; ui->lineEditLegendXYWidth->setText(LegWidth) ; }

            if(Commlines1[ds].first == "AddErrorBarInGraphs" ){
                if(Commlines1[ds].second == "yes" ){
                    ui->checkBoxAddErrorBar->setChecked(true);
                }
            }

            if(Commlines1[ds].first == "GraphsExt" ){ui->AnalysisComboBoxGraphsExt->setCurrentText(Commlines1[ds].second) ; }

            showResultsOutput(Commlines1[ds].first + " " + Commlines1[ds].second , 4);
            //AnalysisData.push_back()
        }
    */

}
void MainWindow::on_AnalysisbtnGraphsDirPath_clicked()
{
    if(QFile::exists(UserCurrentResultsDirPath+"/"+GraphsOutDirName)){
        QString command = UserCurrentResultsDirPath+"/"+GraphsOutDirName;
        QProcess process;
        QStringList qsl = {command};
        process.startDetached("nautilus", qsl);
    }else{
        showResultsOutput("cannot find analysis graphs directory", 4);
    }
}
void MainWindow::on_AnalysisbtnGenerate_clicked()
{

    if(UserCurrentResultsDirPath == "" || UserCurrentResultsDirPath.isEmpty()){UserCurrentResultsDirPath = DoseCalcsCore_build_dir_path + "/" + ResultDirectoryName;}

    if(!QFile::exists(UserCurrentResultsDirPath+"/"+ResultFileName)){
        //QMessageBox::information(this, tr(""), "Cannot find ResultsData file in the given result directory. Choose a result directory that contains ResultsData file to perform root analysis tasks." );
        UserCurrentResultsDirPath = QFileDialog::getExistingDirectory(0, ("Cannot find ResultsData file in the given result directory. Choose a result directory that contains ResultsData file to perform root analysis tasks."), UserCurrentResultsDirPath);
        if(UserCurrentResultsDirPath.isEmpty()){
            return;
        }
    }

    if(!QFile::exists(UserCurrentResultsDirPath+"/"+GraphsOutDirName)){
        QDir* dir = new QDir(UserCurrentResultsDirPath+"/"+GraphsOutDirName);
        dir->mkdir(UserCurrentResultsDirPath+"/"+GraphsOutDirName);
    }

    if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+GraphExecutableName)){

        QDir dir(UserCurrentResultsDirPath+"/"+GraphsOutDirName);
        if(dir.count() != 0){
            QString nn = "The directory "+UserCurrentResultsDirPath+"/"+GraphsOutDirName+" not empty. Click on \"yes\" if you want removing these existed files before generating the new files";

            if(QMessageBox::Yes == QMessageBox::question(this, tr(""),nn)){
                foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {dir.remove(entry.fileName());}
            }
        }

        initializeVariable();
        SaveDataFromInputComponents(); // get the same componenet data but the four under variable values are related to the SaveEnePharOrgLists() data that is called one time when the the run in pushed
        CreateUserCommands(); //fill the variables with all values to save it to user file and not inputFile executed by the geant4 application core
        fileManagerObject->WriteTextToFile(DoseCalcsCore_build_dir_path+"/"+MacroFileName, generateInputUserTextForinputFile());

        if(!QFile::exists(DoseCalcsCore_build_dir_path+"/"+MacroFileName)){
            QMessageBox::information(this, tr(""), "Cannot create macros file "+DoseCalcsCore_build_dir_path+"/"+MacroFileName+" for analysis" );
        }


        if(!QFile::exists(Root_Lib_dir_path+"/thisroot.sh")){
            Root_Lib_dir_path = ShowMessageBoxAndGetPath("Directory Containing thisroot.sh Not Found, Click OK to Choose the Directory, else, you will just generate .tex and .csv tables without graphs");
        }

        //CoreProcess.execute("sh " + UserCurrentResultsDirPath+"/"+DoseCalcsExecutingFileName);
        BashCommandsForExecuting = "#! /bin/bash \n . " +Root_Lib_dir_path +"/thisroot.sh \n"+
                "cd " + DoseCalcsCore_build_dir_path + "\n" +
                DoseCalcsCore_build_dir_path+"/"+GraphExecutableName + " " + DoseCalcsCore_build_dir_path+"/"+MacroFileName
                ;
        BashCommandsForExecuting += "\n bash \n";
        fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

        showResultsOutput("Writing Analysis Commands : \n", 0);
        showResultsOutput(BashCommandsForExecuting , 0);
        showResultsOutput("to --> " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , 0);

        if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName)){

            showResultsOutput("Generating DoseCalcs Results in tables .tex and .csv and ROOT graphs and histograms" , 1);

            ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);

        }
        else{
            showResultsOutput("Cannot find file containing execution commands "+ DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName + " , you should build DoseCalcs with ROOT Analysis option" , 3);
        }

    }else{
        QMessageBox::information(this, tr("Performing root analysis tasks"), "Cannot find analysis executable. Reinstall DoseCalcs-Core with  -DWITH_ANALYSIS_USE=ON option in cmake command. Also if you will generate graphs, set ROOT analysis system bin path in the DoseCalcs-Core installation cmake command -DROOT_DIR=/home/../root_install/bin" );
    }
}
void MainWindow::on_RootAnalysispushButtonOpenROOTGUI_clicked()
{


    if(!QFile::exists(UserCurrentResultsDirPath+"/"+GraphsOutDirName)){
        QString nn = "Generate graphs first. "+UserCurrentResultsDirPath+"/"+GraphsOutDirName+" directory not Found. Open ROOT browser GUI";
        if (QMessageBox::Yes == QMessageBox::question(this, tr(""),nn)){}
        else{return;}
    }

    if(!QFile::exists(Root_Lib_dir_path+"/thisroot.sh")){
        Root_Lib_dir_path = ShowMessageBoxAndGetPath("Directory containing thisroot.sh Not Found, Click OK to Choose the Directory");
    }

    //CoreProcess.execute("sh " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
    if(QFile::exists(Root_Lib_dir_path) && QFile::exists(Root_Lib_dir_path+"/thisroot.sh") ){
        BashCommandsForExecuting = "\ncd " +Root_Lib_dir_path +"\n" + ". ./thisroot.sh \n"
                +"cd "+UserCurrentResultsDirPath+"/"+GraphsOutDirName+" \n"
                +"root --web=off -e \"TBrowser x\" "
                ;
    }else{
        BashCommandsForExecuting = "cd "+UserCurrentResultsDirPath+"/"+GraphsOutDirName+" \n root --web=off -e \"TBrowser x\" "
                ;
    }

    BashCommandsForExecuting += "\n bash \n";
    fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

    if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName)){

        showResultsOutput("Opening ROOT graphical user interface" , 1);

        ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);

    }
    else{
        showResultsOutput("Cannot find file containing execution commands "+ DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName + " , you should build DoseCalcs with ROOT Analysis option" , 3);
    }
}
void MainWindow::AddNewMaterialOrRegionToListOfSourceRegionNames(){


    QStringList args = getArgumentOfACommandFromText(VOXELCommands[11], 5);
    QStringList ll;
    if(args.size() > 0){
        ll.push_back("");
        if(args[0] == "yes"){

            for(int dd=0; dd < MaterialRegionsNames.size();dd++){
                ll.push_back(MaterialRegionsNames[dd]);
            }
        }
        else{
            for(int dd=0; dd < MaterialRegionsNames.size();dd++){
                ll.push_back(MaterialRegionsNames[dd]);
            }
        }

    }else{
        for(int dd=0; dd < MaterialRegionsNames.size();dd++){ ll.push_back(MaterialRegionsNames[dd]);}
    }
    ui->comboBoxSourceRegionProposed->clear();
    ui->comboBoxSourceRegionProposed->addItems(ll);
}
void MainWindow::on_comboBoxSourceRegionProposed_textActivated(const QString &arg1)
{
    bool IsIn = false;
    QStringList args = ui->lineEditChosenSourceTypeData->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
    for(int dd=0; dd < args.size();dd++){ if(args[dd] == arg1){ IsIn = true;}}
    if(IsIn == false){

        if(arg1 != ""){
            QString nn = ui->lineEditChosenSourceTypeData->text();
            ui->lineEditChosenSourceTypeData->setText(nn + " " + arg1);
        }
    }
}

void MainWindow::on_pushButton_2LoadICRPSpectrum_clicked()
{
    ReadLoadICRPSpectrumData();
    //QFuture<void> future = QtConcurrent::run(this, &MainWindow::ReadLoadICRPSpectrumData);
}
void MainWindow::ReadLoadICRPSpectrumData(){

    if(!IsICRPFilesAreRead){

        Read_ICRP107SpectrumRadiationFiles(ICRPDATAPath);

        if(ICRPRadioNuclideDataDiscSpec.size() != 0){
            IsICRPFilesAreRead = true;
        }
    }

    QDialog * d = new QDialog(); d->setWindowTitle("ICRP Spectrums");

    QGridLayout* GraphLayout = new QGridLayout;

    QStringList radlist ;
    for ( auto it = ICRPRadioNuclideDataDiscSpec.begin(); it != ICRPRadioNuclideDataDiscSpec.end(); ++it  ){

        QString RadioTracer_NAME = it.key();
        //        QTextStream(stdout) << "--Geometry_NAME " << Geometry_NAME <<"\n";

        radlist.push_back(RadioTracer_NAME);

    }

    QComboBox * RadioList = new QComboBox(); RadioList->addItems(radlist);
    RadioList->setToolTip("Choose a radionuclide");
    RadioList->setCurrentText("");

    QComboBox * RadioPart = new QComboBox(); RadioPart->addItems(DefinedParticlesNames);
    RadioPart->setToolTip("Choose a particle");
    RadioPart->setCurrentText("");

    int ii = 0, jj=0;
    GraphLayout->addWidget(RadioList, jj,ii,1,1);
    GraphLayout->addWidget(RadioPart, jj,++ii,1,1);

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    GraphLayout->addWidget(buttonBox);

    d->setLayout(GraphLayout);

    int result = d->exec();

    if(result == QDialog::Accepted)
    {

        if(ICRPRadioNuclideDataDiscSpec[RadioList->currentText()][RadioPart->currentText()]["Spectrum"].size()==0)
        {
            QMessageBox::information(this, tr(""), "Warning! There is no spectrum data for the particle "+RadioPart->currentText()+" emitted from " +RadioList->currentText());
            return;
        }

        double mean = 0.22222;
        QString textt = " 0 0";
        for ( auto it = ICRPRadioNuclideDataDiscSpec[RadioList->currentText()][RadioPart->currentText()]["Spectrum"].begin(); it != ICRPRadioNuclideDataDiscSpec[RadioList->currentText()][RadioPart->currentText()]["Spectrum"].end(); ++it  ){

            double Energy  = it.key(); mean = mean + Energy;
            double SAFValue  = it.value();
            textt += " " + QString::number(Energy) + " " + QString::number(SAFValue) ;
            //QTextStream(stdout) << " Particle_NAME " << RadioPart->currentText() << " Energy " << Energy << " Prob " << SAFValue <<"\n";
        }

        if(ui->SourcelineEditParName->text() != RadioPart->currentText() ){
            QMessageBox::information(this, tr(""), "Warning! The particle of spectrum not the same as the particle given in the source radiation particle inputs, the source particle will be changed to "+RadioPart->currentText()+".");
            ui->SourcelineEditParName->setText(RadioPart->currentText());
        }
        ui->lineEditSpecialEnergyDistributionParameter->setText(QString::number(mean/ICRPRadioNuclideDataDiscSpec[RadioList->currentText()][RadioPart->currentText()]["Spectrum"].size())+" " +textt);
    }
}
void MainWindow::on_pushButton_2RootEnergyView_clicked()
{
    CreateROOTFile();
}
void MainWindow::OpenROOTCanvas()
{

    if(!QFile::exists(Root_Lib_dir_path+"/thisroot.sh")){
        Root_Lib_dir_path = ShowMessageBoxAndGetPath("Directory containing thisroot.sh Not Found, Click OK to Choose the Directory");
    }

    if( QFile::exists(Root_Lib_dir_path) && QFile::exists(Root_Lib_dir_path+"/thisroot.sh") && QFile::exists(DoseCalcsCore_build_dir_path+"/"+GraphExecutableName) ){

        BashCommandsForExecuting = "#! /bin/bash \ncd " +Root_Lib_dir_path +"\n" + ". ./thisroot.sh \n"
                +"cd "+DoseCalcsCore_build_dir_path+" \n"
                +"root -l -q Canvas.C"
                ;
        BashCommandsForExecuting += "\n bash \n";
        fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

        if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName)){

            showResultsOutput("Opening ROOT graphical Canvas" , 1);

            ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
        }
    }else{
        QMessageBox::information(this, tr(""), "Install ROOT Analysis System, and rebuild DoseCalcs-Core with analysis option before run the ROOT-View");
    }
}
void MainWindow::CreateROOTFile()
{

    QDialog * d = new QDialog(); d->setWindowTitle("Samples Number in ROOT Histogram");

    QGridLayout* GraphLayout = new QGridLayout;

    QStringList radlist ;
    for ( auto it = ICRPRadioNuclideDataDiscSpec.begin(); it != ICRPRadioNuclideDataDiscSpec.end(); ++it  ){

        QString RadioTracer_NAME = it.key();
        //        QTextStream(stdout) << "--Geometry_NAME " << Geometry_NAME <<"\n";

        radlist.push_back(RadioTracer_NAME);

    }


    QLabel * Labell = new QLabel("Sample Number");

    QDoubleSpinBox * doubleSpinBoxNum = new QDoubleSpinBox();
    doubleSpinBoxNum->setMaximum(1000000000);
    doubleSpinBoxNum->setToolTip("Add number (default is 1000) of samples to generate an histogram with the chosen distribution (i.e. 10000)");

    int ii = 0, jj=0;
    GraphLayout->addWidget(Labell, jj,ii,1,1);
    GraphLayout->addWidget(doubleSpinBoxNum, jj,++ii,1,1);

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    GraphLayout->addWidget(buttonBox);

    d->setLayout(GraphLayout);

    int result = d->exec();
    int numberOfData = 1000;

    if(result == QDialog::Accepted)
    {
        numberOfData = doubleSpinBoxNum->value();

    }


    QString text = "";
    bool open = true;
    QStringList InputsVals;

    //QString nnn = "0.2222222 0.000 0.00E+00 4.14e-07 1.51834e-14 1e-06 9.53797e-13 1e-05 1.05004e-11 5e-05 2.12085e-11 0.0001 6.00885e-11 0.0002 1.70091e-10 0.0004 3.46901e-10 0.0007 4.30984e-10 0.001 4.3634e-09 0.003 9.8687e-09 0.006 1.75521e-08 0.01 5.9753e-08 0.02 1.67755e-07 0.04 2.15148e-07 0.06 2.52307e-07 0.08 2.83391e-07 0.1 8.19421e-07 0.15 9.46597e-07 0.2 1.04714e-06 0.25 1.12894e-06 0.3 1.19649e-06 0.35 1.25268e-06 0.4 1.29951e-06 0.45 1.33846e-06 0.5 1.37064e-06 0.55 1.39695e-06 0.6 2.85286e-06 0.7 2.90359e-06 0.8 2.92682e-06 0.9 2.92785e-06 1 5.78989e-06 1.2 5.61695e-06 1.4 5.3727e-06 1.6 5.08312e-06 1.8 4.76705e-06 2 6.53296e-06 2.3 5.79237e-06 2.6 6.62417e-06 3 6.68499e-06 3.5 5.1738e-06 4 3.94373e-06 4.5 2.96873e-06 5 3.84379e-06 6 2.06529e-06 7 1.08069e-06 8 5.53594e-07 9 2.78648e-07 10 1.38192e-07 11 6.76674e-08 12 3.27686e-08 13 1.5714e-08 14 7.47027e-09 15 0 ";
    //InputsVals = nnn.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
    InputsVals = ui->lineEditSpecialEnergyDistributionParameter->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);

    if(ui->SourceComboBoxEnergyDist->currentText() == "Spectrum" ){

        text += "void Canvas(){\n"
                "TH1F* h = new TH1F(\"Data\", \"Energy Spectrum\","+QString::number(numberOfData)+",0.0,"+InputsVals[InputsVals.size()-2]+");\n";

        //QString nnn = "0.2222222 0.000 0.00E+00 4.14e-07 1.51834e-14 1e-06 9.53797e-13 1e-05 1.05004e-11 5e-05 2.12085e-11 0.0001 6.00885e-11 0.0002 1.70091e-10 0.0004 3.46901e-10 0.0007 4.30984e-10 0.001 4.3634e-09 0.003 9.8687e-09 0.006 1.75521e-08 0.01 5.9753e-08 0.02 1.67755e-07 0.04 2.15148e-07 0.06 2.52307e-07 0.08 2.83391e-07 0.1 8.19421e-07 0.15 9.46597e-07 0.2 1.04714e-06 0.25 1.12894e-06 0.3 1.19649e-06 0.35 1.25268e-06 0.4 1.29951e-06 0.45 1.33846e-06 0.5 1.37064e-06 0.55 1.39695e-06 0.6 2.85286e-06 0.7 2.90359e-06 0.8 2.92682e-06 0.9 2.92785e-06 1 5.78989e-06 1.2 5.61695e-06 1.4 5.3727e-06 1.6 5.08312e-06 1.8 4.76705e-06 2 6.53296e-06 2.3 5.79237e-06 2.6 6.62417e-06 3 6.68499e-06 3.5 5.1738e-06 4 3.94373e-06 4.5 2.96873e-06 5 3.84379e-06 6 2.06529e-06 7 1.08069e-06 8 5.53594e-07 9 2.78648e-07 10 1.38192e-07 11 6.76674e-08 12 3.27686e-08 13 1.5714e-08 14 7.47027e-09 15 0 ";
        //InputsVals = nnn.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
        QVector<double> RadioNuclideProbVec;
        QVector<QVector<double>> RadioNuclideEneVec;

        int m = 0;
        for (int i = 1; i < InputsVals.size(); i++) {

            if(i == InputsVals.size()-2){
                break;
            }

            double V = InputsVals[i].toDouble();
            QVector<double> V12;
            V12.push_back(V);
            i++;
            V12.push_back(InputsVals[i+1].toDouble());
            RadioNuclideEneVec.push_back(V12);

            double W = InputsVals[i].toDouble();

            RadioNuclideProbVec.push_back(W);
            //QTextStream(stdout) << " m= " << m  << " E1= " << V12[0]<< " E2= " << V12[1]  << " W= " << W << "\n";

            m++;
        }

        double TotalProbability;
        for(int m = 0; m < RadioNuclideProbVec.size() ;m++ ){
            TotalProbability += RadioNuclideProbVec[m];
            //QTextStream(stdout) << " m= " << m  << " TotalProbability= " << TotalProbability << "\n";

        }

        for(int f = 0; f < numberOfData ;f++ ){

            double RV = TotalProbability*((float)rand() / RAND_MAX) ;
            double EnergyCumulatedProbability = 0;

            //QTextStream(stdout)  << " ((float)rand() / RAND_MAX)= " << ((float)rand() / RAND_MAX) << " TotalProbability= " << TotalProbability  << " EnergyCumulatedProbability= " << EnergyCumulatedProbability  << " RV= " << RV  << "\n";

            for(int m = 0; m < RadioNuclideProbVec.size() ;m++ ){

                if(RV < EnergyCumulatedProbability){

                    double RandomEne = RadioNuclideEneVec[m][0] + (RadioNuclideEneVec[m][1]-RadioNuclideEneVec[m][0])*((float)rand() / RAND_MAX);
                    text += "h->Fill("+QString::number(RandomEne)+");\n";

                    //QTextStream(stdout) << " RandomEne= " << RandomEne << " Min-Max " << RadioNuclideEneVec[m][0] << " " << RadioNuclideEneVec[m][1] << "\n";

                    //std::cout << f << " Particle "<< RadioNuclidePartNameVec[m] << " TotalProbability "<< TotalProbability << " RandomProb " << RV << " Prob "<< RadioNuclideProbVec[m] << " EnergyCumulatedProbability " << EnergyCumulatedProbability  << " SpectrumOrDiscreteVec " << RadioNuclideSpectrumOrDiscreteVec[m] << " Min " << RadioNuclideEneVec[m][0]   << " RandomEne " << EnergyList[f] << " Max " << RadioNuclideEneVec[m][1] << std::endl;

                    break;
                }
                //QTextStream(stdout) << " EnergyCumulatedProbability= " << EnergyCumulatedProbability  << " RV= " << RV  << "\n";

                EnergyCumulatedProbability += RadioNuclideProbVec[m];
            }
        }

    } else if(ui->SourceComboBoxEnergyDist->currentText() == "Gauss" ){

        if(InputsVals.size() > 1){

            text += "void Canvas(){\n"
                    "TH1F* h = new TH1F(\"Data\", \"Energy Spectrum\","+QString::number(numberOfData)+",0.0,"+QString::number(InputsVals[1].toDouble()*10)+");\n";


            //QTextStream(stdout) << " InputsVals[0].toDouble()= " << InputsVals[0].toDouble()  << " InputsVals[1].toDouble()= " << InputsVals[1].toDouble() << "\n";

            text +="for (int i = 0; i < "+QString::number(numberOfData)+"; i++) {\n";
            text +="h->Fill(gRandom->Gaus("+QString::number(InputsVals[1].toDouble())+", "+QString::number(InputsVals[0].toDouble())+"));\n";
            text +="}\n";
        }
    } else if(ui->SourceComboBoxEnergyDist->currentText() == "Rayleigh" ){

        if(InputsVals.size() > 0){

            text += "void Canvas(){\n"
                    "TH1F* h = new TH1F(\"Data\", \"Energy Spectrum\","+QString::number(numberOfData)+",0.0,"+QString::number((InputsVals[0].toDouble()/3)*10)+");\n";

            //QTextStream(stdout) << " InputsVals[0].toDouble()= " << InputsVals[0].toDouble()  << " InputsVals[1].toDouble()= " << InputsVals[1].toDouble() << "\n";

            text +="double mean = "+QString::number(InputsVals[0].toDouble()/3.)+" ;\n";
            text +="double sigma = mean * std::sqrt(1/M_PI) ;\n";

            text +="for (int i = 0; i < "+QString::number(numberOfData)+"; i++) {\n";

            text +="double ENERGY = sigma*std::sqrt(-2*std::log(1-gRandom->Rndm()));\n";

            //text +="double sigma = ("+QString::number(InputsVals[0].toDouble()/3)+") / sqrt(2*TMath::Pi());";
            //text +="double ENERGY = sqrt(-2.0 * log(gRandom->Rndm())) * sigma;";

            text +="h->Fill(ENERGY);\n";
            text +="}\n";
        }
    } else if(ui->SourceComboBoxEnergyDist->currentText() == "Uniform" ){

        if(InputsVals.size() > 0){

            text += "void Canvas(){\n"
                    "TH1F* h = new TH1F(\"Data\", \"Energy Spectrum\","+QString::number(numberOfData)+",0.0,"+QString::number((InputsVals[0].toDouble()/3)*10)+");\n";

            //QTextStream(stdout) << " InputsVals[0].toDouble()= " << InputsVals[0].toDouble()  << " InputsVals[1].toDouble()= " << InputsVals[1].toDouble() << "\n";

            text +="double min = "+QString::number(InputsVals[0].toDouble())+" ;\n";
            text +="double max = "+QString::number(InputsVals[1].toDouble())+" ;\n";

            text +="for (int i = 0; i < "+QString::number(numberOfData)+"; i++) {\n";

            text +="double ENERGY = min+(max-min)*gRandom->Rndm();\n";

            text +="h->Fill(ENERGY);\n";
            text +="}\n";
        }
    }else{
        open = false;
        QMessageBox::information(this, tr(""), "Please add spectrum, Rayleigh or Gauss distributions data in energy line input");
    }

    //QTextStream(stdout) << " InputsVals[0].toDouble()= " << InputsVals[0].toDouble()  << " InputsVals[1].toDouble()= " << InputsVals[1].toDouble() << "\n";

    text += "h->GetXaxis()->SetTitle(\"Energy ("+ui->comboBoxEnergyUnit->currentText()+")\");\n"
            "h->GetYaxis()->SetTitle(\"Sample Number\");\n"
            "TCanvas* c = new TCanvas(\"c\", \"Energy Distribution\");\n"
            "//h->Draw();\n"
            "h->Draw(\"HIST SAME\");\n"
            // Set the line color and width\n"
            "h->SetLineColor(kBlue);\n"
            "h->SetLineWidth(1);\n"
            "//c->WaitPrimitive();\n"
            "while (1) {gSystem->ProcessEvents();}\n"
            "//gApplication.Run();\n"
            "}";

    fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/Canvas.C" , text);
    if(open){
        OpenROOTCanvas();
    }
}


// RadioButton value changed
void MainWindow::on_pushButtonChooseVoxIDsFile_clicked()
{
    if(ui->radioButtonVoxIDs->isChecked()){
        QString chosen_DirPath = QFileDialog::getOpenFileName( this, tr("Choose VoxIDs file (region or material IDs file) \".dat\""), DoseCalcsCore_build_dir_path+"/"+ScriptDirectoryName, tr("All files (*.*)") );

        if(chosen_DirPath.isEmpty() || chosen_DirPath.isNull()){}else{
            ui->lineEditVoxIDsFilePath->setText(chosen_DirPath);
        }
    }else if(ui->radioButtonTET->isChecked()){

        QString chosen_DirPath = QFileDialog::getOpenFileName( this, tr("Choose Tetrahedrons nodes, then elements file"), DoseCalcsCore_build_dir_path+"/"+ScriptDirectoryName, tr("All files (*.*)") );

        if(chosen_DirPath.isEmpty() || chosen_DirPath.isNull()){

        }else{
            if(QFile::exists(ui->lineEditVoxIDsFilePath->text())){
                ui->lineEditVoxIDsFilePath->setText(ui->lineEditVoxIDsFilePath->text() + " " + chosen_DirPath);
            }else{
                ui->lineEditVoxIDsFilePath->setText(chosen_DirPath);
            }
        }
    }
}
void MainWindow::on_radioButtonTET_clicked()
{
    ui->comboBoxTypeOfSources->setCurrentText("TET");
    ui->checkBoxVoxelOrRegionLevel->setVisible(false);
}
void MainWindow::on_radioButtonVoxel_clicked(bool checked)
{
    if(checked == true){

        ui->frame_19->setVisible(false);
        ui->checkBoxVoxelOrRegionLevel->setVisible(true);
        ui->checkBoxVoxelOrRegionLevel->setChecked(true);

        ui->comboBoxTypeOfSources->setCurrentText("TET");
        //ui->PhantomWorldMaterialLineEdit->setEnabled(true);
        //ui->PhantomWorldHalfSizeslineEdit->setEnabled(true);

        //ui->frame_World->setEnabled(true);
        //ui->frame_materials->setEnabled(true);
        //ui->GeometryFilePathLineEdit->setVisible(false);
        //ui->btnOpenGeometryFilePath->setEnabled(false);

        RemoveDynamiqueGeomAndMatFrame();
        //showConstructVoxelizedCommandsFrame();
    }

}
void MainWindow::on_radioButtonDICOM_clicked(bool checked)
{

    FirstTimeAdd = true;
    if(checked == true){

        ui->frame_19->setVisible(false);
        ui->checkBoxVoxelOrRegionLevel->setVisible(true);
        ui->checkBoxVoxelOrRegionLevel->setChecked(true);

        ui->comboBoxTypeOfSources->setCurrentText("Voxels");
        //ui->frame_World->setEnabled(true);
        //ui->frame_materials->setEnabled(true);
        //ui->GeometryFilePathLineEdit->setVisible(false);
        //ui->btnOpenGeometryFilePath->setEnabled(false);

        RemoveDynamiqueGeomAndMatFrame();
        //showConstructVoxelizedCommandsFrame();

    }
}
void MainWindow::on_radioButtonVoxIDs_clicked(bool checked)
{
    if(checked == true){

        ui->frame_19->setVisible(true);
        ui->checkBoxVoxelOrRegionLevel->setVisible(true);
        ui->checkBoxVoxelOrRegionLevel->setChecked(true);

        ui->comboBoxTypeOfSources->setCurrentText("Voxels");
        //ui->frame_World->setEnabled(true);
        //ui->frame_materials->setEnabled(true);
        //ui->btnOpenGeometryFilePath->setEnabled(true);
        //ui->GeometryFilePathLineEdit->setVisible(true);

        RemoveDynamiqueGeomAndMatFrame();
        //showConstructVoxelizedCommandsFrame();

    }
}
void MainWindow::on_radioButtonConstruct_clicked(bool checked)
{
    if(checked == true){

        ui->frame_19->setVisible(false);

        ui->comboBoxTypeOfSources->setCurrentText("Volume");
        ui->checkBoxVoxelOrRegionLevel->setVisible(false);
        //ui->frame_World->setEnabled(true);
        //ui->frame_materials->setEnabled(true);
        //ui->GeometryFilePathLineEdit->setVisible(false);
        //ui->btnOpenGeometryFilePath->setEnabled(false);

        RemoveDynamiqueGeomAndMatFrame();
        //showConstructSolidAndVolumeCommandsFrame();

    }

}
void MainWindow::on_radioButtonGDML_clicked(bool checked)
{
    if(checked == true){

        ui->comboBoxTypeOfSources->setCurrentText("Volume");
        ui->frame_19->setVisible(false);
        ui->checkBoxVoxelOrRegionLevel->setVisible(false);

        //ui->frame_World->setEnabled(false);
        //ui->frame_materials->setEnabled(false);
        //ui->GeometryFilePathLineEdit->setVisible(true);
        //ui->btnOpenGeometryFilePath->setEnabled(true);

        RemoveDynamiqueGeomAndMatFrame();
    }
}
void MainWindow::on_radioButtonTEXT_clicked(bool checked)
{
    if(checked == true){

        ui->frame_19->setVisible(false);
        ui->checkBoxVoxelOrRegionLevel->setVisible(false);

        ui->comboBoxTypeOfSources->setCurrentText("Volume");
        //ui->frame_World->setEnabled(false);
        //ui->frame_materials->setEnabled(false);
        //ui->GeometryFilePathLineEdit->setVisible(true);
        //ui->btnOpenGeometryFilePath->setEnabled(true);

        RemoveDynamiqueGeomAndMatFrame();

    }
}
void MainWindow::on_radioButtonCpp_clicked(bool checked)
{
    if(checked == true){

        ui->frame_19->setVisible(false);
        ui->checkBoxVoxelOrRegionLevel->setVisible(false);

        ui->comboBoxTypeOfSources->setCurrentText("Volume");
        //ui->frame_World->setEnabled(true);
        //ui->frame_materials->setEnabled(true);
        //ui->GeometryFilePathLineEdit->setVisible(false);
        //ui->btnOpenGeometryFilePath->setEnabled(false);

        RemoveDynamiqueGeomAndMatFrame();
        //showConstructSolidAndVolumeCommandsFrame();

    }
}
void MainWindow::on_radioButtonSTL_clicked(bool checked)
{
    if(checked == true){

        ui->frame_19->setVisible(false);
        ui->checkBoxVoxelOrRegionLevel->setVisible(false);

        ui->comboBoxTypeOfSources->setCurrentText("Volume");
        //ui->frame_World->setEnabled(true);
        //ui->frame_materials->setEnabled(true);
        //ui->GeometryFilePathLineEdit->setVisible(false);
        //ui->btnOpenGeometryFilePath->setEnabled(false);

        RemoveDynamiqueGeomAndMatFrame();
        //showConstructSolidAndVolumeCommandsFrame();

    }
}

// combobox Values changed
void MainWindow::on_SourceComboBoxPhysUsed_currentIndexChanged(const QString &arg1)
{
    if(arg1 == "Construct"){

        ui->comboBoxComptonModels->setVisible(true);
        ui->comboBoxPEEModels->setVisible(true);
        ui->comboBoxGammaConversionModels->setVisible(true);
        ui->comboBoxRayleighScatteringModels->setVisible(true);
        ui->comboBoxElectronBremModels->setVisible(true);
        ui->comboBoxElectronIonisationModels->setVisible(true);
        ui->comboBoxHadronIonisationModels->setVisible(true);

        ui->labelPEE->setVisible(true);
        ui->label_3->setVisible(true);
        ui->label_4->setVisible(true);
        ui->label_5->setVisible(true);
        ui->label_6->setVisible(true);
        ui->label_7->setVisible(true);
        ui->label_8->setVisible(true);

    }else{
        ui->comboBoxComptonModels->setVisible(false);
        ui->comboBoxPEEModels->setVisible(false);
        ui->comboBoxGammaConversionModels->setVisible(false);
        ui->comboBoxRayleighScatteringModels->setVisible(false);
        ui->comboBoxElectronBremModels->setVisible(false);
        ui->comboBoxElectronIonisationModels->setVisible(false);
        ui->comboBoxHadronIonisationModels->setVisible(false);

        ui->labelPEE->setVisible(false);
        ui->label_3->setVisible(false);
        ui->label_4->setVisible(false);
        ui->label_5->setVisible(false);
        ui->label_6->setVisible(false);
        ui->label_7->setVisible(false);
        ui->label_8->setVisible(false);
    }
}
void MainWindow::on_comboBoxTypeOfSources_currentIndexChanged(const QString &arg1)
{
    if(ui->comboBoxTypeOfSources->currentText() == "Voxels" || ui->comboBoxTypeOfSources->currentText() == "TET" ){
        ui->comboBoxSourceRegionProposed->setVisible(true);
    }else{
        ui->comboBoxSourceRegionProposed->setVisible(false);
    }


    if(ui->comboBoxTypeOfSources->currentText() == "Volume"){
        ui->lineEditChosenSourceTypeData->setPlaceholderText("VolumeName1 hx hy hz VolumeName2 hx2 hy2 hz2 ... \n");
        ui->lineEditChosenSourceTypeData->setToolTip("VolumeName1 hx hy hz VolumeName2 hx2 hy2 hz2 VolumeName2 hx3 hy3 hz3... :\n"
                                                     "Liver 15. 8. 8. Thyroid 5 4 3.5 Brain 12 14 13.5)"
                                                     );
    }
    else if(ui->comboBoxTypeOfSources->currentText() == "Voxels"){
        ui->lineEditChosenSourceTypeData->setPlaceholderText("RegionName1 RegionName2 RegionName3 ...\n");
        ui->lineEditChosenSourceTypeData->setToolTip("RegionName1 RegionName2 RegionName3 ... :\n"
                                                     "Liver Thyroid Brain \n"
                                                     );
    }
    else if(ui->comboBoxTypeOfSources->currentText() == "TET"){
        ui->lineEditChosenSourceTypeData->setPlaceholderText("RegionName1 RegionName2 RegionName3 ...\n");
        ui->lineEditChosenSourceTypeData->setToolTip("RegionName1 RegionName2 RegionName3 ... :\n"
                                                     "Liver Thyroid Brain \n"
                                                     );
    }
    else if(ui->comboBoxTypeOfSources->currentText() == "Point"){
        ui->lineEditChosenSourceTypeData->setPlaceholderText("SourceName SourcePosition(X Y Z) RotationAxis RotationAngle");
        ui->lineEditChosenSourceTypeData->setToolTip("SourceName SourcePosition(X Y Z) RotationAxis RotationAngle: "
                                                     "P1 2 15 3.2 Z 45"
                                                     );
    }
    else if(ui->comboBoxTypeOfSources->currentText() == "Beam"){
        ui->lineEditChosenSourceTypeData->setPlaceholderText("SourceName SourcePosition(X Y Z) BeamSDev RotationAxis RotationAngle");
        ui->lineEditChosenSourceTypeData->setToolTip("SourceName SourcePosition(X Y Z) BeamAxis BeamSDev RotationAxis RotationAngle:\n"
                                                     "B1 10 5 8 Z 2 Y 45\n"
                                                     );
    }
    else if(ui->comboBoxTypeOfSources->currentText() == "Plane"){
        ui->lineEditChosenSourceTypeData->setPlaceholderText("SourceName SourcePosition(X Y Z) PlaneShape(Axis) ShapeData RotationAxis RotationAngle");
        ui->lineEditChosenSourceTypeData->setToolTip("SourceName SourcePosition(X Y Z) PlaneShape(Axis) ShapeData RotationAxis RotationAngle: \n"
                                                     "Source1 10 5 8 Square X 16 Z 45\n"
                                                     "Source1 10 5 8 Rectangle X 11 12 Z 45\n"
                                                     "Source1 10 5 8 Circle X 16 Z 45\n"
                                                     "Source1 10 5 8 Ellipse X 11 12 Z 45\n"
                                                     "Source1 10 5 8 Annulus X 11 15 Z 45\n"
                                                     );
    }
    else if(ui->comboBoxTypeOfSources->currentText() == "Surface"){
        ui->lineEditChosenSourceTypeData->setPlaceholderText("SourceName SourcePosition(X Y Z) SurfaceShape ShapeData RotationAxis RotationAngle");
        ui->lineEditChosenSourceTypeData->setToolTip("SourceName SourcePosition(X Y Z) SurfaceShape ShapeData RotationAxis RotationAngle: \n"
                                                     "Source1 10 5 8 Sphere 11 Z 45\n"
                                                     "Source1 10 5 8 Ellipsoid 11 13 5 Z 45\n"
                                                     "Source1 10 5 8 Cylinder 11 13 Z 45\n"
                                                     );
    }
    else if(ui->comboBoxTypeOfSources->currentText() == "Solid"){
        ui->lineEditChosenSourceTypeData->setPlaceholderText("SourceName SourcePosition(X Y Z) SolidShape ShapeData RotationAxis RotationAngle");
        ui->lineEditChosenSourceTypeData->setToolTip("SourceName SourcePosition(X Y Z) SolidShape ShapeData RotationAxis RotationAngle: \n"
                                                     "Source1 10 5 8 Sphere 11 Z 45\n"
                                                     "Source1 10 5 8 Ellipsoid 11 13 5 Z 45\n"
                                                     "Source1 10 5 8 Para 11 13 5 Z 45\n"
                                                     "Source1 10 5 8 Cylinder 11 13 Z 45\n"
                                                     "Source1 10 5 8 EllipticCylinder 11 13 5 Z 45\n"
                                                     );
    }
}
void MainWindow::on_SourceComboBoxEnergyDist_currentIndexChanged(const QString &arg1)
{
    if(arg1 == "Mono" || arg1 == "File" ){
        ui->pushButton_2RootEnergyView->setVisible(false);
    }else{
        ui->pushButton_2RootEnergyView->setVisible(true);
    }


    if(arg1 == "Spectrum"){
        ui->pushButton_2LoadICRPSpectrum->setVisible(true);
    }
    else{
        ui->pushButton_2LoadICRPSpectrum->setVisible(false);
    }


    if(arg1 == "Mono"){
        ui->lineEditSpecialEnergyDistributionParameter->setPlaceholderText("");
        ui->lineEditSpecialEnergyDistributionParameter->setEnabled(true);
    }
    else if(arg1 == "Uniform"){
        ui->lineEditSpecialEnergyDistributionParameter->setPlaceholderText("minE maxE ");
        ui->lineEditSpecialEnergyDistributionParameter->setEnabled(true);
    }
    else if(arg1 == "Gauss"){
        ui->lineEditSpecialEnergyDistributionParameter->setPlaceholderText("Emean SDev");
        ui->lineEditSpecialEnergyDistributionParameter->setEnabled(true);
    }
    else if(arg1 == "Rayleigh"){
        ui->lineEditSpecialEnergyDistributionParameter->setPlaceholderText("Emax");
        ui->lineEditSpecialEnergyDistributionParameter->setEnabled(true);
    }
    else if(arg1 == "Spectrum"){
        ui->lineEditSpecialEnergyDistributionParameter->setPlaceholderText("Energy-probability list");
        ui->lineEditSpecialEnergyDistributionParameter->setEnabled(true);
    }
    else if(arg1 == "File"){
        ui->lineEditSpecialEnergyDistributionParameter->setPlaceholderText("File path");
        ui->lineEditSpecialEnergyDistributionParameter->setEnabled(true);
    }
}
void MainWindow::on_SourceComboBoxAngleDist_currentIndexChanged(const QString &arg1)
{
    if(arg1 == "Directed"){
        ui->lineEditSpecialAngulatDistributionParameter->setPlaceholderText("ThetaPhi 90 270");
        ui->lineEditSpecialAngulatDistributionParameter->setToolTip("ThetaPhi Theta Phi (ThetaPhi 90 270) (unit in degree or rad)\n"
                                                     "ToPoint x y z (ToPoint 12 58 90) (unit in cm or mm)\n"
                                                     "ParallelTo Plane Axis2Val Axis3Val (ParallelTo XY 12 58 or ParallelTo Y 90 12 or ParallelTo X 58 90) (unit in cm or mm)\n"
                                                     "ToVolume x y z hx hy hz (ToVolume 0 0 0 20 12 90) (unit in cm or mm)");
        //ui->lineEditSpecialAngulatDistributionParameter->setEnabled(true);
    }else{
        ui->lineEditSpecialAngulatDistributionParameter->setPlaceholderText("");
        //ui->lineEditSpecialAngulatDistributionParameter->setEnabled(false);
    }
}
void MainWindow::on_AnalysisComboBoxGraphData_currentIndexChanged(const QString &arg1)
{
    if(arg1 == "Reference_Result"){
            ui->AnalysisLineEditRefFile->setEnabled(true);
            ui->AnalysisLineEdit_RefName->setEnabled(true);
            ui->AnalysisComboBoxGraphsType->setEnabled(true);
            ui->comboBoxRelDiff->setEnabled(true);

    }else{
        ui->AnalysisLineEdit_RefName->setEnabled(false);
        ui->AnalysisLineEditRefFile->setEnabled(false);
        ui->btnReferenceFile->setEnabled(false);
        ui->comboBoxRelDiff->setEnabled(false);
    }

    if(arg1 == "none"){
        //ui->AnalysisLineEditRefFile->setEnabled(false);
        //ui->AnalysisLineEdit_RefName->setEnabled(false);
        //ui->AnalysisComboBoxGraphsExt->setEnabled(false);
        //ui->AnalysisComboBoxGraphsType->setEnabled(false);
        //ui->btnReferenceFile->setEnabled(false);
    }
    else if(arg1 == "Result"){
        //ui->AnalysisLineEditRefFile->setEnabled(false);
        //ui->AnalysisLineEdit_RefName->setEnabled(false);
        //ui->btnReferenceFile->setEnabled(false);
        //ui->AnalysisComboBoxGraphsExt->setEnabled(true);
        //ui->AnalysisComboBoxGraphsType->setEnabled(true);
    }
    else if(arg1 == "Reference_Result"){
        //ui->AnalysisLineEditRefFile->setEnabled(true);
        //ui->AnalysisLineEdit_RefName->setEnabled(true);
        //ui->AnalysisComboBoxGraphsExt->setEnabled(true);
        //ui->AnalysisComboBoxGraphsType->setEnabled(true);
        //ui->btnReferenceFile->setEnabled(true);

    }
}


void MainWindow::on_pushButtonEditGeomFile_clicked()
{
    if(EditFlag == 0){

    }
    else if (EditFlag == 7){ // .gdml, .text, c++ file text saving

        if(QFile::exists(PhyVolName->text())){
            fileManagerObject->WriteTextToFile( PhyVolName->text() , ui->GeometryFileTextEdit->toPlainText());
        }
        else{
            QString ext = QString::fromLocal8Bit(getFileExt(PhyVolName->text().toStdString()).c_str());
            QString fn = QString::fromLocal8Bit(getFileNameFromPath(PhyVolName->text().toStdString()).c_str());

            fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+ScriptDirectoryName+"/"+fn+"."+ext , ui->GeometryFileTextEdit->toPlainText());
        }
    }
    else if (EditFlag == 1){ // ExeFileName file text saving
        fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+ExeFileName , ui->GeometryFileTextEdit->toPlainText());
    }
    else if (EditFlag == 2){ // Geometry data saving

        if(ui->radioButtonGDML->isChecked() || ui->radioButtonTEXT->isChecked() || ui->radioButtonCpp->isChecked() || ui->radioButtonSTL->isChecked() || ui->radioButtonConstruct->isChecked()){

            if(ui->GeometryFileTextEdit->toPlainText() != ""){

                showResultsOutput("Saving Volumes Constructing commands ", 4);
                CONSCommandsString = ui->GeometryFileTextEdit->toPlainText();
                UpdateMatVolSolList(2);
            }
        }
        if(ui->radioButtonTET->isChecked() || ui->radioButtonVoxIDs->isChecked() || ui->radioButtonDICOM->isChecked() || ui->radioButtonVoxel->isChecked()){

            if(ui->GeometryFileTextEdit->toPlainText() != ""){

                showResultsOutput("Saving Voxelized/TET  Geometry commands ", 4);
                VOXELCommandsString = ui->GeometryFileTextEdit->toPlainText();
            }
        }
        if(ui->radioButtonConstruct->isChecked()){

            if(ui->GeometryFileTextEdit->toPlainText() != ""){

                showResultsOutput("Saving Volumes Constructing commands ", 4);
                CONSCommandsString = ui->GeometryFileTextEdit->toPlainText();
            }
        }

        ui->GeometryFileTextEdit->clear();
        RemoveDynamiqueGeomAndMatFrame();
    }
    else if (EditFlag == 3){ // Materials data saving

        if(ui->GeometryFileTextEdit->toPlainText() != ""){
            showResultsOutput("Saving Materials Data commands ", 4);
            MaterialsDataCommandsString = ui->GeometryFileTextEdit->toPlainText();
            UpdateMatVolSolList(3);
        }
        ui->GeometryFileTextEdit->clear();
        ui->PhantomWorldMaterialLineEdit->clear();
        for(int dd=0; dd < MaterialsNames.size();dd++){
            ui->PhantomWorldMaterialLineEdit->addItem(MaterialsNames[dd]);
        }
        RemoveDynamiqueGeomAndMatFrame();
    }
    else if (EditFlag == 4){ // macros file for analysis file text saving
        fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+MacroFileName , ui->GeometryFileTextEdit->toPlainText());
    }
    else if (EditFlag == 8){ // MacroFileName file text saving
        fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+MacroFileName , ui->GeometryFileTextEdit->toPlainText());
        if(ui->checkBoxUseMacroCommandFile->isChecked()){
            FillComponentsFromInputsFile(DoseCalcsCore_build_dir_path+"/"+MacroFileName);
        }
    }
    else if (EditFlag == 9){ // ResultFileName file text saving
        fileManagerObject->WriteTextToFile( UserCurrentResultsDirPath+"/"+ResultFileName , ui->GeometryFileTextEdit->toPlainText());
    }
    else if (EditFlag == 10){ // SaveMacrosFilePath file text saving
        fileManagerObject->WriteTextToFile( SaveMacrosFilePath , ui->GeometryFileTextEdit->toPlainText());
    }
    else if (EditFlag == 5){ // World(.gdml, .text, .c++) file text saving

        //CheckIfIsInputsUsed();

        QString ext = QString::fromLocal8Bit(getFileExt(ui->PhantomWorldHalfSizeslineEdit->text().toStdString()).c_str());

        if(ext == "gdml" || ext == "geom"){
            fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/Scripts/"+ui->PhantomWorldHalfSizeslineEdit->text() , ui->GeometryFileTextEdit->toPlainText());
        }
        else if(ext == "c++" || ext == "cpp" || ext == "cc"){

            QString oldd = fileManagerObject->ReadTextFromFileInOneString(DoseCalcsCore_source_dir_path+"/src/G4TCPPGeometryFormat.cc");
            QString neww = ui->GeometryFileTextEdit->toPlainText();

            if(oldd != neww){

                fileManagerObject->WriteTextToFile( DoseCalcsCore_source_dir_path+"/src/G4TCPPGeometryFormat.cc" , ui->GeometryFileTextEdit->toPlainText());

                on_BuildButton_clicked();
            /*
                if(!QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName)){
                    showResultsOutput("Cannot find DoseCalcs executable "+DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName + " . Please build DoseCalcs to run it", 3);
                    on_BuildButton_clicked();
                }

                QString f = "";

                BashCommandsForExecuting = "#! /bin/bash \n" +f
                        + "cd " + DoseCalcsCore_build_dir_path + "\n"
                        + cmakeTruePath + " --build . --target " + DoseCalcsCore_source_dir_path+"/src/G4TCPPGeometryFormat.cc.o"
                        ;
                BashCommandsForExecuting += "\n bash \n";

                showResultsOutput("Writing building G4TCPPGeometryFormat class commands: \n", 0);
                showResultsOutput(BashCommandsForExecuting , 0);
                fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

                ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
            */
            }
        }
    }
    else if (EditFlag == 11){ // SaveMacrosFilePath file text saving
        fileManagerObject->WriteTextToFile( AnyOpenedFilePath , ui->GeometryFileTextEdit->toPlainText());
    }
    else if (EditFlag == 12){ // DoseCalcs update (download and install) commands file text saving
        fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , ui->GeometryFileTextEdit->toPlainText());
        ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
    }

    ui->tabWidget->setTabText(0,"Inputs Edit");
    ui->GeometryFileTextEdit->clear();
    EditFlag = 0;
}
void MainWindow::on_pushButtonCloseFile_clicked()
{
    ui->tabWidget->setTabText(0,"Inputs Edit");
    ui->GeometryFileTextEdit->clear();
    EditFlag = 0;
}
void MainWindow::UpdateMatVolSolList(int inn){
    QStringList InputsVals;
    QVector< QPair<QString,QString>> Commlines ;

    if(inn==3){
        ElementsNames.clear();
        MaterialsNames.clear(); MaterialsNameIDs.clear();

        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(MaterialsDataCommandsString);
        for(int dd=0; dd < Commlines.size();dd++){
            if(Commlines[dd].first == MaterialCommands[0]){ // /MaterialData/createElement
                InputsVals = Commlines[dd].second.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
                ElementsNames.push_back(InputsVals[2]);
                continue;
            }
            if(Commlines[dd].first == MaterialCommands[1] || Commlines[dd].first == MaterialCommands[4] ){ // /MaterialData/createMaterial or /MaterialData/createNISTMaterial
                InputsVals = Commlines[dd].second.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
                MaterialsNames.push_back(InputsVals[0]);
                ElementsNames.push_back(InputsVals[0]);
                continue;
            }
        }

    }else if(inn == 2){
        SolidsNames.clear();
        VolsNames.clear();

        VolsNames.push_back("World");

        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(CONSCommandsString);
        for(int dd=0; dd < Commlines.size();dd++){
            if(Commlines[dd].first == GeometryCommands[2]){ // /GeometryData/createWorld or /GeometryData/createVolume

                InputsVals = Commlines[dd].second.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);

                if(InputsVals.size() < 1){continue;}

                QString fe =  QString::fromLocal8Bit(getFileExt(InputsVals[0].toStdString()).c_str());
                QString fn =  QString::fromLocal8Bit(getFileNameFromPath(InputsVals[0].toStdString()).c_str());

                if(fe == "gdml" || fe == "geom" || fe == "c++" || fe == "cpp" || fe == "cc" || fe == "stl" || fe == "ast"){
                    VolsNames.push_back(fn);
                }
                else if( InputsVals[0] == "GDML" || InputsVals[0] == "TEXT" || InputsVals[0] == "VoxIDs" || InputsVals[0] == "VOXEL" || InputsVals[0] == "DICOM"){
                    continue;
                }
                else{
                    VolsNames.push_back(fn);
                }

                continue;
            }
            if(Commlines[dd].first == GeometryCommands[1]){ // /GeometryData/createSolid
                InputsVals = Commlines[dd].second.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
                SolidsNames.push_back(InputsVals[1]);
                continue;
            }
        }
    }
}
QStringList MainWindow::getArgumentOfACommandFromText(QString command, int whichData){

    QStringList InputsVals;
    QVector< QPair<QString,QString>> Commlines ;

    if(whichData == 1){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(MaterialsDataCommandsString);
    }
    else if(whichData == 2){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(GeometryCommandsString);
    }
    else if(whichData == 3){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(VoxIDsCommandsString);
    }
    else if(whichData == 4){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(DICOMCommandsString);
    }
    else if(whichData == 5){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(VOXELCommandsString);
    }
    else if(whichData == 6){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(PhysicsCommandsString);
    }
    else if(whichData == 7){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(SourceCommandsString);
    }
    else if(whichData == 8){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(RunAndScoreCommandsString);
    }
    else if(whichData == 9){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(AnalysisCommandsString);
    }
    else if(whichData == 10){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(MacrosCommandsString);
    }

    for(int dd=0; dd < Commlines.size();dd++){
        if(Commlines[dd].first == command){
            InputsVals = Commlines[dd].second.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
            return InputsVals;
            break;
        }
    }

    return InputsVals;
}
QString MainWindow::changeArgumentOfACommandFromText(QString Text,QString command, QString NewCommandArg){

    QString NewText = "";
    QVector< QPair<QString,QString>> Commlines ;
    /*
    if(whichData == 1){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(MaterialsDataCommandsString);
    }
    else if(whichData == 2){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(GeometryCommandsString);
    }
    else if(whichData == 3){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(VoxIDsCommandsString);
    }
    else if(whichData == 4){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(DICOMCommandsString);
    }
    else if(whichData == 5){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(VOXELCommandsString);
    }
    else if(whichData == 6){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(PhysicsCommandsString);
    }
    else if(whichData == 7){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(SourceCommandsString);
    }
    else if(whichData == 8){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(RunAndScoreCommandsString);
    }
    else if(whichData == 9){
        Commlines = fileManagerObject->ReadTextFromQStringInQStringList(AnalysisCommandsString);
    }
*/
    Commlines = fileManagerObject->ReadTextFromQStringInQStringList(Text);

    bool ComFound = false;
    for(int dd=0; dd < Commlines.size();dd++){

        if(Commlines[dd].first == command){
            ComFound = true;
            NewText += NewCommandArg + "\n";
        }else{
            NewText += Commlines[dd].first + " " + Commlines[dd].second + "\n";
        }
    }
    if(ComFound == false){
        NewText += NewCommandArg + "\n";
    }

    return NewText;
}




// for multi Geometry files and one physicssourcce data
void MainWindow::on_pushButtonAddGeometryMacros_clicked()
{

    if(ui->checkBoxSimulateGeometriesWitOneSource->isChecked())
    {
        MacrosFiles.clear();

        QFileDialog dialog(this);
        dialog.setDirectory(QDir::homePath());
        dialog.setFileMode(QFileDialog::ExistingFiles);
        dialog.setNameFilter(trUtf8("Splits (*.mac)"));
        if (dialog.exec()){
            MacrosFiles = dialog.selectedFiles();

            QString GoodFilesNames;
            int GoodFilesNum = 0;
            for(int dd=0; dd < MacrosFiles.size();dd++){

                MacrosCommandsString = CreateMaterialAndGeometryDataFromMacrosFile(MacrosFiles[dd]);

                if(MacrosCommandsString != "" || MacrosCommandsString.isEmpty()){

                    Geometry_setGeometrySymbole = "Phantom"+QString::number(dd);
                    if(getArgumentOfACommandFromText(GeometryCommands[4], 10).size() > 0){Geometry_setGeometrySymbole = getArgumentOfACommandFromText(GeometryCommands[4], 10)[0];
                    }else{changeArgumentOfACommandFromText(MacrosCommandsString,GeometryCommands[4],GeometryCommands[4]+" "+Geometry_setGeometrySymbole);}

                    GoodFilesNames = GoodFilesNames + QString::number(dd) + ": " + QString::fromLocal8Bit(getFileNameFromPath(MacrosFiles[dd].toStdString()).c_str()) +" -Symbol: "+ Geometry_setGeometrySymbole + "\n";
                    GoodFilesNum++;
                }
            }

            GoodFilesNames = GoodFilesNames + "\nRead files:" + QString::number(GoodFilesNum);
            ui->LabelRededMacrosFiles->setText(GoodFilesNames);
        }
    }

}
QString MainWindow::CreateMaterialAndGeometryDataFromMacrosFile(QString FilePathString){

    initializeVariable();

    QFile* filee = new QFile(FilePathString);
    if(!filee->exists()){
        showResultsOutput("Cannot find the default inputs file : " + FilePathString , 3);
        return "";
    }

    QStringList InputsVals;

    QString MaterialsAndGeometriesCommands = "";

    QVector< QPair<QString,QString>> Commlines = fileManagerObject->ReadTextFromFileInQStringList(FilePathString);

    for(int dd=0; dd < Commlines.size();dd++){

        for(int cc=0; cc < MaterialCommands.size();cc++){
            if(Commlines[dd].first == MaterialCommands[cc]){
                MaterialsAndGeometriesCommands += Commlines[dd].first + " " + Commlines[dd].second +"\n" ;
                continue;
            }
        }
        for(int cc=0; cc < GeometryCommands.size();cc++){
            if(Commlines[dd].first == GeometryCommands[cc]){
                MaterialsAndGeometriesCommands += Commlines[dd].first + " " + Commlines[dd].second +"\n" ;
                continue;
            }
        }
        for(int cc=0; cc < VOXELCommands.size();cc++){
            if(Commlines[dd].first == VOXELCommands[cc]){
                MaterialsAndGeometriesCommands += Commlines[dd].first + " " + Commlines[dd].second +"\n" ;
                continue;
            }
        }
        for(int cc=0; cc < CONSCommands.size();cc++){
            if(Commlines[dd].first == CONSCommands[cc]){
                MaterialsAndGeometriesCommands += Commlines[dd].first + " " + Commlines[dd].second +"\n" ;
                continue;
            }
        }
    }

    showResultsOutput("================= Materials and Geometries Commands ================= \n\n"+MaterialsAndGeometriesCommands, 4);
    return MaterialsAndGeometriesCommands;

}
void MainWindow::RunForMultiGeomeries()
{

    if(!QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName)){
        showResultsOutput("Cannot find DoseCalcs executable "+DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName + " . Please build DoseCalcs before run", 3);
        on_BuildButton_clicked();
    }

    bool pvr = false;
    QDir dir(UserCurrentResultsDirPath);
    foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "AE@*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {{pvr = true; break;}}
    if(pvr){if (QMessageBox::Yes == QMessageBox::question(this, tr("Simulation result directory"), tr("The current directory is the result of previous simulation, you can choose another directory for results, or create a new one (yes) "))){on_openResultsDirButton_clicked();}}


    for(int qq=0; qq < MacrosFiles.size();qq++){

        MacrosFilePathReadForMultiGeometryAndOneSource = MacrosFiles[qq];
        MacrosCommandsString = CreateMaterialAndGeometryDataFromMacrosFile(MacrosFiles[qq]);

        // to set
        Geometry_setGeometrySymbole = "Phantom"+QString::number(macrosfileinc);
        if(getArgumentOfACommandFromText(GeometryCommands[4], 10).size() > 0){
            Geometry_setGeometrySymbole = getArgumentOfACommandFromText(GeometryCommands[4], 10)[0];
        }

        QString MacrosText = MacrosCommandsString;

        if(MacrosText == ""){
            continue;
        }

        if(!CheckInputsForGeometryFiles()){
            continue;
        }

        SaveSourcePhysicsInputs();
        CreateUserCommands();

        // Physics
        MacrosText += PhysicsCommands[0] + " " + PhysicsData_setPhysicsData + "\n" +
                PhysicsCommands[1]+ " " + PhysicsData_setCutsDistance + "\n" +
                PhysicsCommands[3]+ " " + PhysicsData_setCutsEnergy + "\n" ;

        if(ui->checkBoxGenerateCrossSection->isChecked()){
            if(PhysicsData_ParticleForCrossSection !="" && PhysicsData_EnergiesForCrossSection !=""){
                MacrosText += PhysicsCommands[2]+ " " + PhysicsData_ParticleForCrossSection + " " + PhysicsData_EUnitForCrossSection + " " + PhysicsData_EnergiesForCrossSection ;
            }else{
                MacrosText += "# " + PhysicsCommands[2]+ " " + PhysicsData_ParticleForCrossSection + " " + PhysicsData_EUnitForCrossSection + " " + PhysicsData_EnergiesForCrossSection ;
            }
        }
        MacrosText += "\n\n";

        // Source
        MacrosText += "\n\n" + SourceCommands[0] + " " + SourceData_setParticleName + "\n" +
                SourceCommands[1] + " " + SourceData_setSourcePosData + "\n" +
                SourceCommands[2] + " " + SourceData_setSourceEneData + "\n" +
                SourceCommands[3] + " " + SourceData_setSourceMomDirData + "\n";

        if(ui->UseDataFilesFor->currentIndex() != 0 ){
            MacrosText += SourceCommands[4] + " " + SourceData_UseDataFiles;
        }else{
            MacrosText += "# " +SourceCommands[4] + " " + SourceData_UseDataFiles;
        }
        MacrosText += "\n\n";

        // Score
        MacrosText += RunAndScoreCommands[2] + " " + Execution_setNumberOfRanksOrThreads + "\n" +
                RunAndScoreCommands[5] + " " + ValuesOfInputs[Score_setSimNumOnRanksLineEdit]+ "\n"+
                RunAndScoreCommands[10] + " " + UserCurrentResultsDirPath + "\n\n";

        QsubSeparatedMacroFilePath = DoseCalcsCore_build_dir_path+"/Macros"+ConstructDoseCalcsJobName().remove("DoseCalcs")+".mac";
        ui->RunButton->setText("Run("+QString::number(macrosfileinc)+")");macrosfileinc++;
        if(ui->checkBoxnohup->isChecked()){
            ui->comboBoxNohupFiles->setVisible(true);
            ui->pushButtonShowOutputs->setVisible(true);
            ui->comboBoxNohupFiles->addItem(LastRunOutputFile);
            ui->comboBoxNohupFiles->setCurrentText(LastRunOutputFile);
        }else{
            ui->comboBoxNohupFiles->setVisible(false);
            ui->pushButtonShowOutputs->setVisible(false);
        }

        fileManagerObject->WriteTextToFile(QsubSeparatedMacroFilePath , MacrosText);
        FillComponentsFromInputsFile(QsubSeparatedMacroFilePath);

        if(ui->checkBoxRocks->isChecked()){

            // this file name is used just in exe.sh file in rocks cluster simulations
            QString SpecificNodes = SetToASpecificMPIRank();
            BashCommandsForExecuting = "#! /bin/bash \ncd " + DoseCalcsCore_build_dir_path
                    + "\n qsub "+ SpecificNodes +" "+ DoseCalcsCore_build_dir_path+"/"+ExeFileName;

            if(false){
                BashCommandsForExecuting = "#! /bin/bash \ncd " + DoseCalcsCore_build_dir_path
                        + "\n sbatch "+ SpecificNodes +" "+ DoseCalcsCore_build_dir_path+"/"+ExeFileName;
            }
            BashCommandsForExecuting += "\n bash \n";

            showResultsOutput("Writing Run Commands : \n", 0);
            showResultsOutput(BashCommandsForExecuting , 0);
            showResultsOutput("to --> " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , 4);
            fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

            //if(EditFlag == 1){ // if you late it open, is not saved, can be for the last execution
            on_pushButtonGenerateExe_clicked();
            on_pushButtonEditGeomFile_clicked();
            //}
            ui->outputTextConsole->setPlainText(fileManagerObject->ReadTextFromFileInOneString(DoseCalcsCore_build_dir_path+"/"+ExeFileName));
            ui->tabWidget->setCurrentIndex(0);

            if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName) && QFile::exists(DoseCalcsCore_build_dir_path+"/"+ExeFileName)){

                if(!ShowImportantSimulationData()){return;}
                showResultsOutput("Computation Run" , 1);

                //ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
            }
        }else{

            if(QFile::exists(QsubSeparatedMacroFilePath) && QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName)){

                if(MPI_USE=="YES"){

                    if(!QFile::exists(MPI_Lib_dir_path+"/mpirun")){
                        MPI_Lib_dir_path = ShowMessageBoxAndGetPath("Directory containing mpirun or mpiexec Not Found, Click OK to Choose the Directory");
                    }

                    QString nohup = "";
                    QString andd = "";
                    if(ui->checkBoxnohup->isChecked()){
                        LastRunOutputFile = "nohup_"+ConstructDoseCalcsJobName();

                        if(!QFile::exists("/usr/bin/nohup")){
                            QMessageBox::information(this, tr(""), "/usr/bin/nohup Not Found, Please install nohup and check the /usr/bin/nohup path.");
                            return;
                        }
                        nohup = "/usr/bin/nohup ";
                        andd = "  > " + LastRunOutputFile + " & ";
                    }

                    BashCommandsForExecuting = "#! /bin/bash \ncd " + DoseCalcsCore_build_dir_path + "\n"
                            + nohup + " "+MPI_Lib_dir_path+"/mpirun " + Execution_setNumberOfRanksOrThreads + " "
                            + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName + " B " + QsubSeparatedMacroFilePath + " " + Execution_setEventNumber + " " +andd
                            ;
                }else{

                    QString f = "";

                    if(!QFile::exists(geant4_Lib_dir_path+"/geant4.sh")){
                        geant4_Lib_dir_path = ShowMessageBoxAndGetPath("Directory containing geant4.sh Not Found, Click OK to Choose the Directory");
                    }

                    QString nohup = "";
                    QString andd = "";
                    if(ui->checkBoxnohup->isChecked()){
                        LastRunOutputFile = "nohup_"+ConstructDoseCalcsJobName();

                        if(!QFile::exists("/usr/bin/nohup")){
                            QMessageBox::information(this, tr(""), "/usr/bin/nohup Not Found, Please install nohup and check the /usr/bin/nohup path.");
                            return;
                        }
                        nohup = "/usr/bin/nohup ";
                        andd = "  > " + LastRunOutputFile + " & ";
                    }

                    BashCommandsForExecuting = "#! /bin/bash \n" +f
                            + "cd " +geant4_Lib_dir_path +"\n"+
                            + ". ./geant4.sh\n" +f +
                            + "cd " + DoseCalcsCore_build_dir_path + "\n"
                            + nohup +" "+DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName + " B " + QsubSeparatedMacroFilePath + " " + Execution_setEventNumber + " " +andd
                            ;
                }
                BashCommandsForExecuting += "\n bash \n";

                showResultsOutput("Writing Run Commands : \n", 0);
                showResultsOutput(BashCommandsForExecuting , 0);
                showResultsOutput("to --> " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , 4);
                fileManagerObject->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

                if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName)){

                    if(!ShowImportantSimulationData()){return;}
                    showResultsOutput("Computation Run" , 1);
                    ShowTerminal(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
                }
                else{
                    showResultsOutput("Cannot find file containing execution commands "+ DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName + " , you should build DoseCalcs with ROOT Analysis option" , 3);
                }

                isInexec = true;

                //process->close();

            }else{
                showResultsOutput("Verify that the simulate binary and the inputs.mac are existed in " + QsubSeparatedMacroFilePath, 0);
            }
        }
    }
}
bool MainWindow::CheckInputsForGeometryFiles(){

    QStringList InputsVals;
    QString sizes = "";

    // check World data

    //QMessageBox::information(this, tr(""), getArgumentOfACommandFromText(GeometryCommands[2], 10)[0]);
    //QMessageBox::information(this, tr(""), getArgumentOfACommandFromText(GeometryCommands[2], 10)[1]);
    //QMessageBox::information(this, tr(""), getArgumentOfACommandFromText(GeometryCommands[0], 10)[0]);
    //QMessageBox::information(this, tr(""), getArgumentOfACommandFromText(GeometryCommands[0], 10)[1]);

    if(getArgumentOfACommandFromText(GeometryCommands[2], 10).size() > 0){

        if(getArgumentOfACommandFromText(GeometryCommands[2], 10)[0] == "VoxIDs"  ||
                getArgumentOfACommandFromText(GeometryCommands[2], 10)[0] == "VOXEL" ||
                getArgumentOfACommandFromText(GeometryCommands[2], 10)[0] == "DICOM" ||
                getArgumentOfACommandFromText(GeometryCommands[2], 10)[0] == "TET" ){
            InputsVals = getArgumentOfACommandFromText(GeometryCommands[0], 10); // "/GeometryData/createWorld"
            double xl,yl,zl,vxl,vyl,vzl;
            if(InputsVals.size() > 3){

                xl = QString(InputsVals[1]).toDouble() * UnitsConversionFactors[ui->comboBoxWorldSizeUnit->currentText()];
                yl = QString(InputsVals[2]).toDouble() * UnitsConversionFactors[ui->comboBoxWorldSizeUnit->currentText()];
                zl = QString(InputsVals[3]).toDouble() * UnitsConversionFactors[ui->comboBoxWorldSizeUnit->currentText()];

                //showResultsOutput(InputsVals[0] + " " + InputsVals[1] +InputsVals[2] + " " + ui->comboBoxWorldSizeUnit->currentText() , 4);

                InputsVals = getArgumentOfACommandFromText(VOXELCommands[0], 5);
                if(InputsVals.size() == 9){
                    vxl = QString(InputsVals[0]).toDouble() *2* QString(InputsVals[5]).toDouble() * UnitsConversionFactors[InputsVals[8]];
                    vyl = QString(InputsVals[1]).toDouble() *2* QString(InputsVals[6]).toDouble() * UnitsConversionFactors[InputsVals[8]];
                    vzl = QString(InputsVals[2]).toDouble() *2* QString(InputsVals[7]).toDouble() * UnitsConversionFactors[InputsVals[8]];
                    //showResultsOutput(InputsVals[0] + "*"+InputsVals[5]+ " " + InputsVals[1] + "*" +InputsVals[6]+ " " + InputsVals[2] + "*" +InputsVals[7] + " " + InputsVals[8], 4);
                }

                if(xl < vxl || yl < vyl || zl < vzl){
                    sizes = "("+QString::number(xl)+";"+QString::number(vxl)+"), ("+QString::number(yl)+";"+QString::number(vyl)+"), ("+QString::number(zl)+";"+QString::number(vzl)+")";
                    ui->Tab->setCurrentIndex(0);
                    QMessageBox::information(this, tr(""), "The Voxelized phantom geometrical limits in command("+ VOXELCommands[0] +") exceed the simulation world dimensions " + sizes +". Edit the world dimension (i.e. 100 100 100) using (" + GeometryCommands[0] +") or voxel data (half x, y, and z) in " + MacrosFilePathReadForMultiGeometryAndOneSource);
                    return false;
                }
            }else{
                ui->Tab->setCurrentIndex(0);
                QMessageBox::information(this, tr(""), "Canno't find the World volume data (" + GeometryCommands[0] +") in " + MacrosFilePathReadForMultiGeometryAndOneSource);
                return false;
            }
        }
    }


    // check geometry data
    if(getArgumentOfACommandFromText(GeometryCommands[2], 10).size() > 0){
        if(getArgumentOfACommandFromText(GeometryCommands[2], 10)[0] == "VoxIDs"){

            if(getArgumentOfACommandFromText(GeometryCommands[2], 10).size() > 1){
                if(QFile::exists(getArgumentOfACommandFromText(GeometryCommands[2], 10)[1])){

                }else{
                    ui->Tab->setCurrentIndex(0);
                    QMessageBox::information(this, tr(""), "Canno't find the voxels IDs data file in macros file command " + MacrosFilePathReadForMultiGeometryAndOneSource);
                    return false;
                }
            }else{
                QMessageBox::information(this, tr(""), "Canno't find the voxels IDs data file in macros file command " + MacrosFilePathReadForMultiGeometryAndOneSource);
                return false;
            }
        }
    }
    if(getArgumentOfACommandFromText(GeometryCommands[2], 10).size() > 0){
        if(getArgumentOfACommandFromText(GeometryCommands[0], 10)[0] == "TET"){
            if(getArgumentOfACommandFromText(GeometryCommands[0], 10).size() < 2){
                ui->Tab->setCurrentIndex(0);
                QMessageBox::information(this, tr(""), "Add the two geometry file paths in macros file "  + MacrosFilePathReadForMultiGeometryAndOneSource +", the first fro tetrahedrons nodes data and the seconde for tetrahedrons elements data");
                return false;
            }else{
                if(QFile::exists(getArgumentOfACommandFromText(GeometryCommands[0], 10)[1]) && QFile::exists(InputsVals[1])){
                }else{ui->Tab->setCurrentIndex(0);
                    QMessageBox::information(this, tr(""), "Add the two geometry file paths in macros file "  + MacrosFilePathReadForMultiGeometryAndOneSource +", the first fro tetrahedrons nodes data and the seconde for tetrahedrons elements data");
                    return false;
                }
            }
        }
    }

    // check particle data
    InputsVals = ui->SourcelineEditParName->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    if(InputsVals.size() < 1){
        ui->Tab->setCurrentIndex(1);
        QMessageBox::information(this, tr(""), "Add particle name (i.e. gamma, e-, e+, alpha and proton)");
        return false;
    }else{
        for(int aa=0; aa < InputsVals.size();aa++){
            bool IsIn = false;
            for(int dd=0; dd < DefinedParticlesNames.size();dd++){ if(DefinedParticlesNames[dd] == InputsVals[aa]){ IsIn = true;}}
            if(IsIn == false){
                ui->Tab->setCurrentIndex(1);
                QMessageBox::information(this, tr(""), "The particle\""+ InputsVals[aa] +"\" is not known by DoseCalcs. Remove it from the source particle LineEdit widget");
                return false;
            }
        }

    }

    // check source region data
    InputsVals = ui->lineEditChosenSourceTypeData->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    if(ui->comboBoxTypeOfSources->currentText()=="Volume"){
        if(InputsVals.size() < 4){
            ui->Tab->setCurrentIndex(1);
            QMessageBox::information(this, tr(""), "Add region name with box dimension (i.e. 10 25 15)");
            return false;
        }
    }else if(ui->comboBoxTypeOfSources->currentText()=="Voxels" || ui->comboBoxTypeOfSources->currentText()=="TET"){
        if(InputsVals.size() < 1){
            ui->Tab->setCurrentIndex(1);
            QMessageBox::information(this, tr(""), "Add region/s name (i.e. Liver)");
            return false;
        }else{

            for(int aa=0; aa < InputsVals.size();aa++){

                QString ss = InputsVals[aa];
                if("allregions" == ss.toLower()){

                    aa++; for(int bb=0; bb < InputsVals[aa].toInt();bb++){
                        aa++;
                        if(aa == InputsVals.size()){
                            break;
                        }
                    }
                }
                else{
                    QStringList args = getArgumentOfACommandFromText(VOXELCommands[11], 5);
                    bool IsIn = false;

                    if(args.size() > 0){
                        if(args[0] == "yes"){for(int dd=0; dd < MaterialRegionsNames.size();dd++){ if(MaterialRegionsNames[dd] == InputsVals[aa]){ IsIn = true;}}}
                        else{ for(int dd=0; dd < DefinedRegionsNames.size();dd++){ if(DefinedRegionsNames[dd] == InputsVals[aa]){ IsIn = true;}}}
                    }else{
                        for(int dd=0; dd < MaterialRegionsNames.size();dd++){ if(MaterialRegionsNames[dd] == InputsVals[aa]){ IsIn = true;}}
                    }

                    if(IsIn == false){
                        ui->Tab->setCurrentIndex(1);
                        QMessageBox::information(this, tr(""), "The region\""+ InputsVals[aa] +"\" is not defined. Remove it from the source region LineEdit widget or add the \""+ InputsVals[aa] +"\" region data");
                        return false;
                    }
                }
            }
        }
    }


    // check energy data
    if(ui->SourceComboBoxEnergyDist->currentText()=="Mono"){
        InputsVals = ui->lineEditSpecialEnergyDistributionParameter->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
        if(InputsVals.size() < 1){
            ui->Tab->setCurrentIndex(1);
            QMessageBox::information(this, tr(""), "Add an energy value");
            return false;
        }
        for(int aa=0; aa < InputsVals.size();aa++){
            QString val = InputsVals[aa];
            if(val.toDouble() == 0. || val.isEmpty() || val.isNull()){
                ui->Tab->setCurrentIndex(1);
                QMessageBox::information(this, tr(""), "The value: \"" + val + "\" not accepted.");
                return false;
            }
        }
    }


    // check momentum direction data
    InputsVals = ui->lineEditSpecialEnergyDistributionParameter->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    if(InputsVals.size() == 0){
        ui->Tab->setCurrentIndex(1);
        QMessageBox::information(this, tr(""), "Add energy value (i.e. 0.5)");
        return false;
    }

    if(ui->SourceComboBoxAngleDist->currentText()=="Directed"){
        InputsVals = ui->lineEditSpecialAngulatDistributionParameter->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
        if(InputsVals.size() < 2){
            ui->Tab->setCurrentIndex(1);
            QMessageBox::information(this, tr(""), "Add theta and phi values for directed momentum distribution (i.e. Liver)");
            return false;
        }
    }

    // check run end score data
    QString val = ui->lineEditNumberOfEvent->text();
    if(val.toDouble() == 0. || val.isEmpty() || val.isNull() || val.toInt() >= INT32_MAX){
        ui->Tab->setCurrentIndex(1);
        QMessageBox::information(this, tr(""), "The number of events: \"" + val + "\" not accepted.");
        return false;
    }
    val = ui->lineEditNumberOfRanksOrThreads->text();
    if(val.toInt() == 0. || val.isEmpty() || val.isNull()){
        ui->Tab->setCurrentIndex(1);
        QMessageBox::information(this, tr(""), "The number of sub-simulations: \"" + val + "\" not accepted.");
        return false;
    }

    return true;

}
void MainWindow::on_checkBoxSimulateGeometriesWitOneSource_clicked(bool checked)
{
    if(checked){
        ui->checkBoxSimulateGeometriesWitOneSource->setText("ON: Choose macros files (.mac)");
    }else{
        ui->checkBoxSimulateGeometriesWitOneSource->setText("OFF");
    }
}
void MainWindow::on_pushButtonCheckMatGeoCommandsOfFiles_clicked()
{
    bool matgeomedatagood = true;
    for(int qq=0; qq < MacrosFiles.size();qq++){
        MacrosFilePathReadForMultiGeometryAndOneSource = MacrosFiles[qq];
        MacrosCommandsString = CreateMaterialAndGeometryDataFromMacrosFile(MacrosFiles[qq]);
        if(ui->checkBoxSimulateGeometriesWitOneSource->isChecked()){
            matgeomedatagood = CheckInputsForGeometryFiles();
        }
    }
    if(matgeomedatagood){
        QMessageBox::information(this, tr(""), "The first inputs check \"GOOD\"");
    }
}

void MainWindow::removeHugFiles_slot(){
//    QMessageBox::information(this, tr(""), "Remove File");

    if(ui->checkBoxRocks->isChecked()){

        bool pvr = false;
        QString files = "";
        QDir dir(DoseCalcsCore_build_dir_path);
        foreach( const QFileInfo& entry, dir.entryInfoList( QStringList() << "core*", QDir::Files | QDir::Hidden | QDir::NoSymLinks ) ) {
            {
                pvr = true;
                files += "* " + entry.filePath() + " ;\n";
                dir.remove(entry.filePath());
            }
        }
        if(pvr){
            QMessageBox::information(this, tr(""), "The files: \n" + files + " are removed. Please check the simulations output files to resubmite the failed simulations");
        }

        timer.singleShot(60000, this, SLOT(removeHugFiles_slot()));
    }else{

    }
}




void MainWindow::on_pushButtonReadICRP107128_clicked()
{
    Read_ICRP107_108Files(ICRPDATAPath);
    if(ui->checkBoxReadSpectrum->isChecked()){
        ReadNeutronSpectrum = false; // the neutron spectrum are not read for SAFs calculations, beacause the SAFs of this neutrons
                                     // are estimated by simulation of spectrum and one SAF value is given for each spectrum
                                     // then, the spectrum is read just to simulate the spectrum in 'MC Simulation' window.
                                     // Each neutron spectrum (cad for each radionuclide) is identified here by a number (111 222 333 ....)
        Read_ICRP107SpectrumRadiationFiles(ICRPDATAPath);
    }
    //CombineSpectrumWithMonoData();

    ui->comboBoxRadioPharmaceutiques->clear();
    for ( auto it = RadiotracerradionucleidMap.begin(); it != RadiotracerradionucleidMap.end(); ++it  ){
            ui->comboBoxRadioPharmaceutiques->addItem(it->first);
    }

    ui->comboBoxSourceOrTargetsForBiokinetics->setCurrentIndex(2);
    ui->doubleSpinBoxAdministeredActivity->setValue(1);

}
void MainWindow::on_pushButtonReadICRPData_clicked()
{

    Read_Reference_file(ICRPDATAPath+"/ICRP133SAFsData");
    ui->pushButtonReadICRPData->setToolTip("Current readed file path: "+ICRPDATAPath+"/ICRP133SAFsData");

    setBiokineticsDefaulsInputs();


    /*
    QStringList args=(QStringList()<<"Pancreas"<<"Spleen"<<"Thyroid");

    QMap<QString,QMap<QString,QMap<double,QMap<QString,QMap<QString,double>>>>> SAFMap = ICRPSAFs["SAF"];


    for ( auto it = SAFMap.begin(); it != SAFMap.end(); ++it  ){

        QString Geometry_NAME = it.key();
        //        QTextStream(stdout) << "--Geometry_NAME " << Geometry_NAME <<"\n";

        for ( auto it2 = it.value().begin(); it2 != it.value().end(); ++it2  ){
            QString Particle_NAME = it2.key();
            //            QTextStream(stdout) << "---Particle_NAME " << Particle_NAME <<"\n";

            for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
                double Energy_Val = it3.key();
                //                QTextStream(stdout) << "----Energy_Val " << Energy_Val <<"\n";

                for ( auto DD = it3.value().begin(); DD != it3.value().end(); ++DD  ){
                    QString Source_NAME  = DD.key();
                    //                    QTextStream(stdout) << "-----Source_NAME " << Source_NAME <<"\n";

                    bool IsIn = false;
                    for(int dd=0; dd < args.size();dd++){ if(args[dd] == Source_NAME){ IsIn = true;}}
                    if(IsIn == false){

                        for ( auto CC = DD.value().begin(); CC != DD.value().end(); ++CC  ){
                            QString  Target_NAME  = CC.key();
                            double SAFValue  = CC.value();

                            ICRPSAFs["SAF"][Geometry_NAME][Particle_NAME][Energy_Val]["OtherTissues"][Target_NAME] += SAFValue;

                            //                        QTextStream(stdout) << "------Target_NAME " << Target_NAME <<"\n";
                            //                        QTextStream(stdout) << "-------SAFValue " << SAFValue <<"\n";
                            //if(Source_NAME == "Liver" && Target_NAME == "Brain"){
                            //QTextStream(stdout) << " Geometry_NAME " << Geometry_NAME << " Particle_NAME " << Particle_NAME << " Energy_Val " << Energy_Val << " Source_NAME " << Source_NAME << " Target_NAME " << Target_NAME << " SAFValue " << SAFValue <<"\n";
                            //}
                        }
                    }
                }
            }
        }
    }
*/


    // to run the reading in another thread you should add this code to a new method (ReadICRPilesAndGetData()) and this method will be run on new thread for example:
    //QFuture<void> future = QtConcurrent::run(this, &MainWindow::ReadICRPilesAndGetData);
}
void MainWindow::on_pushButton_ReadDoseCalcsSAFs_clicked()
{
    Read_DoseCalcs_file(ICRPDATAPath+"/DoseCalcsSAFsData");
    ReadMassesAndWTFactor(ICRPDATAPath+"/DoseCalcsRegionsData");

    if(ICRPSAFs.size() == 0){
        ui->pushButtonReadICRPData->setText("Read ICRP133 SAFs");
        ui->pushButton_ReadDoseCalcsSAFs->setText("Read DoseCalcs SAFs");
        ui->pushButtonReadUSERData->setText("Read User SAFs");

    }else{
        ui->pushButtonReadICRPData->setText("Read ICRP133 SAFs");
        ui->pushButton_ReadDoseCalcsSAFs->setText("Read DoseCalcs SAFs *");
        ui->pushButtonReadUSERData->setText("Read User SAFs");
    }


    ui->pushButton_ReadDoseCalcsSAFs->setToolTip("Current readed file path: "+ICRPDATAPath+"/DoseCalcsSAFsData");
    setBiokineticsDefaulsInputs();

}
void MainWindow::on_pushButtonReadUSERData_clicked()
{
    QString DataDirName = QFileDialog::getOpenFileName( this, tr("Choose DoseCalcs ResultsData file"), UserCurrentResultsDirPath, tr("All files (*)") );

    if(DataDirName.isEmpty() || DataDirName.isNull()){
    }else{
        if(QFile::exists(DataDirName)){
            Read_DoseCalcs_file(DataDirName);
            ReadMassesAndWTFactor(ICRPDATAPath+"/DoseCalcsRegionsData");

            if(ICRPSAFs.size() == 0){
                ui->pushButtonReadICRPData->setText("Read ICRP133 SAFs");
                ui->pushButton_ReadDoseCalcsSAFs->setText("Read DoseCalcs SAFs");
                ui->pushButtonReadUSERData->setText("Read User SAFs");

            }else{
                ui->pushButtonReadICRPData->setText("Read ICRP133 SAFs");
                ui->pushButton_ReadDoseCalcsSAFs->setText("Read DoseCalcs SAFs");
                ui->pushButtonReadUSERData->setText("Read User SAFs *");
            }
            ui->pushButtonReadUSERData->setToolTip("Current readed file path: "+DataDirName);
            setBiokineticsDefaulsInputs();
        }
    }
    //QFuture<void> future = QtConcurrent::run(this, &MainWindow::ReadUserilesAndGetData);
}

void MainWindow::setBiokineticsDefaulsInputs(){

    QTextStream(stdout) <<"Data reading terminated \n";

    if( ICRPSAFs.size() == 0){
        QMessageBox::information(this, tr(""), "No SAF data were registered");
    }

    QStringList Geometries;
    for ( auto it = ICRPSAFs.begin(); it != ICRPSAFs.end(); ++it  ){
        for ( auto it2 = it.value().begin(); it2 != it.value().end(); ++it2  ){
            bool isIn = false;
            for(int gg = 0 ; gg < Geometries.size() ; gg++){
                if(it2.key() == Geometries[gg]){
                    isIn = true; break;
                }
            }
            if(isIn == false){
                Geometries.push_back(it2.key());
            }
        }
    }
    ui->comboBoxPhantom->clear();
    ui->comboBoxPhantom->addItems(Geometries);

    //QStringList Organs;
    //for ( auto it = PhantomSourcesMap[ui->comboBoxPhantom->currentText()].begin(); it != PhantomSourcesMap[ui->comboBoxPhantom->currentText()].end(); ++it  ){Organs.push_back(it->first);}
    //ui->comboBoxSources->clear();
    ui->comboBoxSources->clear();
    ui->comboBoxSources->addItems(PhantomSourcesMap[ui->comboBoxPhantom->currentText()].toList());


    QStringList Organs1;
    for ( auto it = RegionParameterValueMap[ui->comboBoxPhantom->currentText()]["Mass"].begin(); it != RegionParameterValueMap[ui->comboBoxPhantom->currentText()]["Mass"].end(); ++it  ){Organs1.push_back(it.key());}
    ui->comboBoxTargets->clear();
    //ui->comboBoxTargets->addItem("all");
    ui->comboBoxTargets->addItems(Organs1);


    // these energies should be sorted one time for all source particle to be used in the calculations of E1 E2 for E value finding
    for ( auto Abeg = SourceParticleEnergyValues.begin(); Abeg != SourceParticleEnergyValues.end(); ++Abeg  )
    {
        for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  ){
            std::sort(Bbeg->second.begin(), Bbeg->second.end());
        }
    }

    //ui->comboBoxQuantityNucl->setCurrentText("SAF");
    ui->comboBoxSources->setCurrentText("Liver");
    //ui->comboBoxTargets->setCurrentText("Brain");

    ConstructOtherTissuesSource();
}
void MainWindow::CalculateQuantitiesBasedOnICRPData()
{
    // quantity, geometry, radionuclide, source, target, value

    RadionuclideValuePairVectorAlphabetical.clear(); // Radionuclide, QuantityValue
    RadionuclideValuePairVectorValuesAscending.clear(); // Radionuclide, QuantityValue
    RadionuclideValuePairVectorValuesDescending.clear(); // Radionuclide, QuantityValue

    QString Quantity_NAME = ui->comboBoxQuantityNucl->currentText();
    QString Geometry_NAME = ui->comboBoxPhantom->currentText();

    QString RadioTracer_NAME;
    QString Particle_NAME;
    QString Source_NAME;
    QString Target_NAME;
    double Energy_Val;

    //QTextStream(stdout) <<"The quantities values will be calculated for "<< ICRPRadioNuclideData.size()<<" radionuclide\n";

    double QuantityMinValue;
    double QuantityMaxValue;
    if(ui->checkBoxEffDoseLimNucl->isChecked()){
        QuantityMinValue = ui->doubleSpinBoxMinEffDose->value();
        if(ui->doubleSpinBoxMaxEffDose->value() == 0.){
            QuantityMaxValue = 1.7e+308; // maxdouble value
        }
        else{
            QuantityMaxValue = ui->doubleSpinBoxMaxEffDose->value();
        }
    }

    double TMinValue;
    double TMaxValue;
    if(ui->checkBoxHalflivesNucl->isChecked()){

        //QTextStream(stdout) << ui->doubleSpinBoxMaxHalflives->value() << " " << ui->doubleSpinBoxMinHalflives->value() << " " << QuantitiesConversionFromDefault["T"][ui->comboBoxPeriodUnit->currentText()] << "\n";

        if(ui->doubleSpinBoxMaxHalflives->value() == 0.){
            TMaxValue = 1.7e+308;// maxdouble value
        }
        else{
            TMaxValue = ui->doubleSpinBoxMaxHalflives->value()/QuantitiesConversionFromDefault["T"][ui->comboBoxPeriodUnit->currentText()];
        }
        TMinValue = ui->doubleSpinBoxMinHalflives->value()/QuantitiesConversionFromDefault["T"][ui->comboBoxPeriodUnit->currentText()];
    }

    //QTextStream(stdout) <<"Number of readed radionuclides data "<< ICRPRadioNuclideData.size() << " " << ICRPRadioNuclideHalfLives.size()<<"\n";

    //QTextStream(stdout) <<"Filter on Emitters and particles and half-lives min-max\n";

    ui->progressBarReadingCalcData->setRange(0, 100);
    ui->progressBarReadingCalcData->setValue(0);
    ui->progressBarReadingCalcData->show();
    double radincc = 0;
    double percent = 0;


    for ( auto it = ICRPRadioNuclideData.begin(); it != ICRPRadioNuclideData.end(); ++it  ){

        percent = (radincc/ICRPRadioNuclideData.size())*100;
        ui->progressBarReadingCalcData->setValue(percent);
        radincc++;

        RadioTracer_NAME = it.key();
        double Value = 0;

        bool SaveOrNotRadionuclide = true;
        if(ui->checkBoxHalflivesNucl->isChecked()){
            if(ICRPRadioNuclideHalfLives[RadioTracer_NAME] >= TMinValue && ICRPRadioNuclideHalfLives[RadioTracer_NAME] < TMaxValue ){
            }else{
                SaveOrNotRadionuclide = false;
            }
        }
        if(SaveOrNotRadionuclide == false){
            continue;
            //QTextStream(stdout) << "1 Remove RadioTracer_NAME " << RadioTracer_NAME << " because its not emitter of " << pn <<"\n";
        }

        if(ui->comboBoxRedionuclidesRadEmmiterType->currentIndex() != 0){

            QString pn="";

            if(ui->comboBoxRedionuclidesRadEmmiterType->currentIndex() == 1){
                pn="gamma";
            }
            else if(ui->comboBoxRedionuclidesRadEmmiterType->currentIndex() == 2){
                pn="e-";
            }
            else if(ui->comboBoxRedionuclidesRadEmmiterType->currentIndex() == 3){
                pn="e+";
            }
            else if(ui->comboBoxRedionuclidesRadEmmiterType->currentIndex() == 4){
                pn="alpha";
            }
            else if(ui->comboBoxRedionuclidesRadEmmiterType->currentIndex() == 5){
                pn="neutron";
            }

            //QTextStream(stdout) << RadioTracer_NAME << " " << it.value().size() <<"\n";
            SaveOrNotRadionuclide = false;
            for ( auto it2 = it.value().begin(); it2 != it.value().end(); ++it2  ){
                //QTextStream(stdout) << it2.key() <<" ";
                if(pn == it2.key()){
                    SaveOrNotRadionuclide = true;
                }
            }
        }

        if(SaveOrNotRadionuclide == false){
            continue;
            //QTextStream(stdout) << "2 Remove RadioTracer_NAME " << RadioTracer_NAME << " because its not emitter of " << pn <<"\n";
        }

        for ( auto it2 = it.value().begin(); it2 != it.value().end(); ++it2  ){
            Particle_NAME = it2.key();

            if(!ui->RadPhotonCheckBox->isChecked() && Particle_NAME == "gamma"){
                continue;
            }
            if(!ui->RadelectronCheckBox->isChecked() && Particle_NAME == "e-"){
                continue;
            }
            if(!ui->RadPositronCheckBox->isChecked() && Particle_NAME == "e+"){
                continue;
            }
            if(!ui->RadAlphaCheckBox->isChecked() && Particle_NAME == "alpha"){
                continue;
            }
            if(!ui->RadNeutronCheckBox->isChecked() && Particle_NAME == "neutron"){
                continue;
            }
            if(!ui->RadFFCheckBox->isChecked() && Particle_NAME == "FF"){
                continue;
            }

            //QTextStream(stdout) << "Particle_NAME " << Particle_NAME << " is checked " << ui->RadNeutronCheckBox->isChecked() <<"\n";

            QString Partial_Particle_NAME = Particle_NAME;
            if(Particle_NAME == "e+"){
                Partial_Particle_NAME = "e-";
            }

            //QTextStream(stdout) << "Particle_NAME " << Particle_NAME <<"\n";
            for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
                Energy_Val = it3.key();
                if(Particle_NAME == "neutron"){
                   Energy_Val = EnergyIDRadionuclideForNeutronSAF[RadioTracer_NAME];
                }

                // here we considere that energy emitted from FF is absorbed locally in the source region
                if(Particle_NAME == "FF" && Target_NAME == Source_NAME){
                    ICRPSAFs["SAF"][Geometry_NAME][Particle_NAME][Energy_Val][Source_NAME][Source_NAME] = 1/RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME];
                    //QTextStream(stdout) << "Particle_NAME " << Particle_NAME << " ICRPSAFs " << ICRPSAFs["AE"][Geometry_NAME][Particle_NAME][Energy_Val][Source_NAME][Source_NAME] <<"\n";
                }

                //QTextStream(stdout) << " RadioTracer_NAME " << RadioTracer_NAME <<  " Particle_NAME " << Particle_NAME <<  " Partial_Particle_NAME " << Partial_Particle_NAME  <<  " Spectrum Energy_Val " << Energy_Val  << " yield " << ICRPRadioNuclideData[RadioTracer_NAME][Particle_NAME][Energy_Val] << " SAF "<<  ICRPSAFs[Quantity_NAME][Geometry_NAME][Partial_Particle_NAME][Energy_Val][ui->comboBoxSources->currentText()][ui->comboBoxTargets->currentText()] << " ValueInc " << Value  << "\n" ;

                GenerateRadiotracerQuantitiesByInterpolationInDefaultUnit(Partial_Particle_NAME,Energy_Val);

                Value += it3.value()
                        *ICRPSAFs[Quantity_NAME][Geometry_NAME][Partial_Particle_NAME][Energy_Val][ui->comboBoxSources->currentText()][ui->comboBoxTargets->currentText()];

                //if(RadioTracer_NAME == "F-18"){
                //    QTextStream(stdout) << " RadioTracer_NAME " << RadioTracer_NAME <<  " Particle_NAME " << Particle_NAME <<  " Partial_Particle_NAME " << Partial_Particle_NAME  <<  " Spectrum Energy_Val " << Energy_Val << " percentt " << percentt << " yield " << ICRPRadioNuclideData[RadioTracer_NAME][Particle_NAME][energyval] << " ValueInc " << Value  << "\n" ;
                //}

                /*
                double yieldd = it3.value();

                bool isin = false;
                for ( auto Oterat = RadioParticleEnergyYieldForSpectrum[RadioTracer_NAME][Particle_NAME].begin(); Oterat != ICRPRadioNuclideData[RadioTracer_NAME][Particle_NAME].end(); ++Oterat  ){
                    if(Energy_Val == Oterat.key() && yieldd == Oterat.value()){
                        isin = true;
                        break;
                    }
                }

                if (isin == true){

                    if(ICRPRadioNuclideDataDiscSpec[RadioTracer_NAME][Partial_Particle_NAME]["Spectrum"].size() == 0){
                        for ( auto it3 = ICRPRadioNuclideDataDiscSpec[RadioTracer_NAME][Partial_Particle_NAME]["Spectrum"].begin(); it3 != ICRPRadioNuclideDataDiscSpec[RadioTracer_NAME][Partial_Particle_NAME]["Spectrum"].end(); ++it3  ){
                            double energyval = it3.key();

                            //QTextStream(stdout) << "energyval " << energyval << " RadiationPerCent " << RadiationPerCent <<"\n";

                            //if(ICRPSAFs[Quantity_NAME][Geometry_NAME][Partial_Particle_NAME][energyval][ui->comboBoxSources->currentText()][ui->comboBoxTargets->currentText()] == 0){
                            //QTextStream(stdout) << "\n-----> For Radiotracer " << RadioTracer_NAME << ", the data for particle "<< Particle_NAME << " with energy "<< energyval << " are not found (this source configuration not simulated), the related quantities values will be generated by interpolating the already existed data of " << Particle_NAME << " with energies surround the " << energyval << " from the existing source regions." << "\n";
                            //  GenerateRadiotracerQuantitiesByInterpolationInDefaultUnit(Partial_Particle_NAME,energyval);
                            //}

                            GenerateRadiotracerQuantitiesByInterpolationInDefaultUnit(Partial_Particle_NAME,energyval);
                            double percentt = it3.value()/ICRPRadioNuclideDataDiscSpec[RadioTracer_NAME][Partial_Particle_NAME]["Total"][1];

                            Value += percentt
                                    *ICRPRadioNuclideData[RadioTracer_NAME][Particle_NAME][energyval]
                                    *ICRPSAFs[Quantity_NAME][Geometry_NAME][Partial_Particle_NAME][energyval][ui->comboBoxSources->currentText()][ui->comboBoxTargets->currentText()];

                            if(RadioTracer_NAME == "F-18"){
                                QTextStream(stdout) << " RadioTracer_NAME " << RadioTracer_NAME <<  " Particle_NAME " << Particle_NAME <<  " Partial_Particle_NAME " << Partial_Particle_NAME  <<  " Spectrum Energy_Val " << energyval << " percentt " << percentt << " yield " << ICRPRadioNuclideData[RadioTracer_NAME][Particle_NAME][energyval] << " ValueInc " << Value  << "\n" ;
                            }
                        }
                    }
                }
                */
                //QTextStream(stdout) << "Energy_Val " << Energy_Val << " RadiationPerCent " << RadiationPerCent <<"\n";

                //if(ICRPSAFs[Quantity_NAME][Geometry_NAME][Partial_Particle_NAME][Energy_Val][ui->comboBoxSources->currentText()][ui->comboBoxTargets->currentText()] == 0){
                //QTextStream(stdout) << "\n-----> For Radiotracer " << RadioTracer_NAME << ", the data for particle "<< Particle_NAME << " with energy "<< Energy_Val << " are not found (this source configuration not simulated), the related quantities values will be generated by interpolating the already existed data of " << Particle_NAME << " with energies surround the " << Energy_Val << " from the existing source regions." << "\n";
                //  GenerateRadiotracerQuantitiesByInterpolationInDefaultUnit(Partial_Particle_NAME,Energy_Val);
                //}
            }
        }

        bool save0 = false;
        bool save1 = false;

        if( !__isnan(Value) && !__isinf(Value) && Value != 0 && Value != NULL){
            save0 = true;
        }
        if(ui->checkBoxEffDoseLimNucl->isChecked()){
            if(Value >= QuantityMinValue && Value < QuantityMaxValue){
                save1 = true;
            }
        }else{
            save1 = true;
        }

        QPair<double,QString> RadVal;

        if(save0 == true && save1 == true){
            RadVal.first = Value;
            RadVal.second = RadioTracer_NAME;
            RadionuclideValuePairVectorAlphabetical.push_back(RadVal);
        }
    }

    ui->progressBarReadingCalcData->setValue(100);

    RadionuclideValuePairVectorValuesAscending = RadionuclideValuePairVectorAlphabetical;
    std::sort(RadionuclideValuePairVectorValuesAscending.begin(), RadionuclideValuePairVectorValuesAscending.end());

    RadionuclideValuePairVectorValuesDescending.clear();
    for(int gg = (int)RadionuclideValuePairVectorValuesAscending.size()-1 ; gg > -1 ; gg--){
        RadionuclideValuePairVectorValuesDescending.push_back(RadionuclideValuePairVectorValuesAscending[gg]);
    }
    ui->pushButtonReverseData->setText("Descending Table, show Alphabetical");

    on_pushButtonReverseData_clicked();

    //QTextStream(stdout) << "The calculation terminated, showing data in table" << "\n";
}
void MainWindow::GenerateSAFFromNewSource(QString SSOURCE){

    QString Quantity_NAME = "SAF";

    QMap<QString,QMap<QString,QMap<double,QMap<QString,QMap<QString,double>>>>> SAFMap = ICRPSAFs[Quantity_NAME];

    ui->progressBarReadingCalcData->setRange(0, 100);
    ui->progressBarReadingCalcData->setValue(0);
    ui->progressBarReadingCalcData->show();
    double percent = 0;
    int qq = 0;

    for ( auto it = SAFMap.begin(); it != SAFMap.end(); ++it  ){

        QString Geometry_NAME = it.key();
        //        QTextStream(stdout) << "--Geometry_NAME " << Geometry_NAME <<"\n";


        percent = ((qq+1)/it->size())*100;
        ui->progressBarReadingCalcData->setValue(percent);
        qq++;

        for ( auto it2 = it.value().begin(); it2 != it.value().end(); ++it2  ){
            QString Particle_NAME = it2.key();
            //            QTextStream(stdout) << "---Particle_NAME " << Particle_NAME <<"\n";

            for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
                double Energy_Val = it3.key();
                //                QTextStream(stdout) << "----Energy_Val " << Energy_Val <<"\n";

                for ( auto DD = it3.value().begin(); DD != it3.value().end(); ++DD  ){
                    QString Source_NAME  = DD.key();

                    bool isin = false;
                    for (int dd = 0 ; dd < CurrentSources.size(); dd++) {
                        if(Source_NAME == CurrentSources[dd]){
                            isin = true;break;
                        }
                    }
                    if(isin == false){
                        continue;
                    }
                    //QTextStream(stdout) << "-----Source_NAME " << Source_NAME << "-----fraction " << RegionFraction[SSOURCE] <<"\n";

                    for ( auto CC = DD.value().begin(); CC != DD.value().end(); ++CC  ){
                        QString  Target_NAME  = CC.key();
                        double SAFValue  = CC.value();

                        ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][SSOURCE][Target_NAME] += SAFValue*RegionFraction[SSOURCE][Source_NAME];
                        ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][SSOURCE][SSOURCE] += SAFValue*RegionFraction[SSOURCE][Source_NAME];
                        bool isin = false;
                        for (int dd = 0 ; dd < SourceParticleEnergyValues[SSOURCE][Particle_NAME].size(); dd++) {
                            if(Energy_Val == SourceParticleEnergyValues[SSOURCE][Particle_NAME][dd]){isin = true;break;}}
                        if(isin == false){SourceParticleEnergyValues[SSOURCE][Particle_NAME].push_back(Energy_Val);}
                    }
                }
            }
        }
    }

    ui->progressBarReadingCalcData->setValue(100);

    QMap<QString,QMap<QString,QMap<QString,double>>> RegionParameterValueMap1 = RegionParameterValueMap;
    for ( auto it = RegionParameterValueMap.begin(); it != RegionParameterValueMap.end(); ++it  ){
        for (int dd = 0 ; dd < CurrentSources.size(); dd++) {
            QString Geometry_NAME = it.key();
            RegionParameterValueMap[Geometry_NAME]["Mass"][SSOURCE] += RegionParameterValueMap1[Geometry_NAME]["Mass"][CurrentSources[dd]]*RegionFraction[SSOURCE][CurrentSources[dd]];
        }

        // The WT are defined in a file for DoseCalcs ICRP Phantoms
        //TissueFactorMap[SSOURCE] += TissueFactorMap[CurrentSources[dd]];
        //TissueFactorMap[CurrentSources[dd]]=0;
    }

}
void MainWindow::GenerateSAFFromNewTarget(QString TTARGET){

    QString Quantity_NAME = "SAF";

    QMap<QString,QMap<QString,QMap<double,QMap<QString,QMap<QString,double>>>>> SAFMap = ICRPSAFs[Quantity_NAME];

    ui->progressBarReadingCalcData->setRange(0, 100);
    ui->progressBarReadingCalcData->setValue(0);
    ui->progressBarReadingCalcData->show();
    double percent = 0;
    int qq = 0;

    for ( auto it = SAFMap.begin(); it != SAFMap.end(); ++it  ){

        QString Geometry_NAME = it.key();
        //        QTextStream(stdout) << "--Geometry_NAME " << Geometry_NAME <<"\n";

        percent = ((qq+1)/it->size())*100;
        ui->progressBarReadingCalcData->setValue(percent);
        qq++;

        for ( auto it2 = it.value().begin(); it2 != it.value().end(); ++it2  ){
            QString Particle_NAME = it2.key();
            //            QTextStream(stdout) << "---Particle_NAME " << Particle_NAME <<"\n";

            for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
                double Energy_Val = it3.key();
                //                QTextStream(stdout) << "----Energy_Val " << Energy_Val <<"\n";

                for ( auto DD = it3.value().begin(); DD != it3.value().end(); ++DD  ){
                    QString Source_NAME  = DD.key();

                    //QTextStream(stdout) << "-----Source_NAME " << Source_NAME <<"\n";
                    for ( auto CC = DD.value().begin(); CC != DD.value().end(); ++CC  ){
                        QString  Target_NAME  = CC.key();


                        bool isin = false;
                        for (int dd = 0 ; dd < CurrentTargets.size(); dd++) {
                            if(Target_NAME == CurrentTargets[dd]){
                                isin = true;break;
                            }
                        }
                        if(isin == false){
                            continue;
                        }

                        //QTextStream(stdout) << "-----Target_NAME " << Target_NAME << "-----fraction " << RegionFraction[TTARGET] <<"\n";

                        double SAFValue  = CC.value();

                        ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][Source_NAME][TTARGET] += SAFValue*RegionFraction[TTARGET][Target_NAME];
                    }
                }
            }
        }
    }

    ui->progressBarReadingCalcData->setValue(100);

    QMap<QString,QMap<QString,QMap<QString,double>>> RegionParameterValueMap1 = RegionParameterValueMap;
    for ( auto it = RegionParameterValueMap.begin(); it != RegionParameterValueMap.end(); ++it  ){
        QString Geometry_NAME = it.key();
        for (int dd = 0 ; dd < CurrentTargets.size(); dd++) {
            RegionParameterValueMap[Geometry_NAME]["Mass"][TTARGET] += RegionParameterValueMap1[Geometry_NAME]["Mass"][CurrentTargets[dd]]*RegionFraction[TTARGET][CurrentTargets[dd]];
        }

        // The WT are defined in a file for DoseCalcs ICRP Phantoms
        //TissueFactorMap[TTARGET] += TissueFactorMap[CurrentTargets[dd]];
        //TissueFactorMap[CurrentTargets[dd]]=0;
    }
}

void MainWindow::GenerateRadiotracerQuantitiesByInterpolationInDefaultUnit(QString ParticleName, double Energy){

    QString Geometry_NAME = ui->comboBoxPhantom->currentText();
    QString Quantity_Name = ui->comboBoxQuantityNucl->currentText();
    QString Source_NAME = ui->comboBoxSources->currentText();
    QString Target_NAME = ui->comboBoxTargets->currentText();

    double Val1 = 0;
    double Val2 = 0;
    double Val = 0;
    double Energy1;
    double Energy2;

    if(ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME] != 0){

        if(Quantity_Name == "AE"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]*Energy);
        }
        else if(Quantity_Name == "AF"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]);
        }
        else if(Quantity_Name == "SAF"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME];
        }
        else if(Quantity_Name == "S"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(Energy);
        }
        else if(Quantity_Name == "H"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(Energy)*GenerateRadiationFactor(ParticleName,Energy);
        }
        else if(Quantity_Name == "E"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(Energy)*GenerateRadiationFactor(ParticleName,Energy)*TissueFactorMap[Target_NAME];
        }

        if(Source_NAME == "Liver" && Target_NAME == "Liver" && Energy == 0.249776){
            QTextStream(stdout) << " Val " << Val << "\n" ;
        }

    }else{

        Energy1 = 0;
        Energy2 = 0;

        for(int ss = 0 ; ss < SourceParticleEnergyValues[Source_NAME][ParticleName].size() ; ss++){

            int ff = ss+1;

            int da = SourceParticleEnergyValues[Source_NAME][ParticleName].size()-1; if(ff == da){break;}

            double E1 = SourceParticleEnergyValues[Source_NAME][ParticleName][ss];
            double E2 = SourceParticleEnergyValues[Source_NAME][ParticleName][ff];
            //QTextStream(stdout) << ss << " " << E1 << " " << ff << " " << E2 << "\n" ;

            if(E1 < Energy && Energy < E2){
                Energy1 = E1;
                Energy2 = E2;

                //QTextStream(stdout) << " ParticleName " << ParticleName << " Energy1 " << Energy1 << " Energy " << Energy << " Energy2 " << Energy2 << "\n" ;
                break;
            }
        }

        if(Quantity_Name == "AE"){
            Val1 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy1][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]*Energy1);
            Val2 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy2][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]*Energy2);
        }
        else if(Quantity_Name == "AF"){
            Val1 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy1][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]);
            Val2 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy2][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]);
        }
        else if(Quantity_Name == "SAF"){
            Val1 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy1][Source_NAME][Target_NAME];
            Val2 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy2][Source_NAME][Target_NAME];
        }
        else if(Quantity_Name == "S"){
            Val1 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy1][Source_NAME][Target_NAME]*(Energy1);
            Val2 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy2][Source_NAME][Target_NAME]*(Energy2);
        }
        else if(Quantity_Name == "H"){
            Val1 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy1][Source_NAME][Target_NAME]*(Energy1)*GenerateRadiationFactor(ParticleName,Energy1);
            Val2 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy2][Source_NAME][Target_NAME]*(Energy2)*GenerateRadiationFactor(ParticleName,Energy2);
        }
        else if(Quantity_Name == "E"){
            Val1 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy1][Source_NAME][Target_NAME]*(Energy1)*GenerateRadiationFactor(ParticleName,Energy1)*TissueFactorMap[Target_NAME];
            Val2 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy2][Source_NAME][Target_NAME]*(Energy2)*GenerateRadiationFactor(ParticleName,Energy2)*TissueFactorMap[Target_NAME];
        }

        Val  = Val1 + (Energy-Energy1)*((Val2-Val1)/(Energy2-Energy1));

        if(Source_NAME == "Liver" && Target_NAME == "Liver" && Energy == 0.249776){
            QTextStream(stdout) << " Val " << Val << "\n" ;
        }

    }

    if( !__isnan(Val) && !__isinf(Val) && Val != 0 && Val != NULL){
        ICRPSAFs[Quantity_Name][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME] = Val;
    }

}
double MainWindow::GenerateRadiationFactor(QString ParticleName, double Energy){ // energy in MeV

    double factor = 1;

    if(ParticleName == "gamma"){
        factor = 1; // all energies
    }
    else if (ParticleName == "e-" || ParticleName == "e+"){
        factor = 1; // all energies
    }
    else if (ParticleName == "alpha"){
        factor = 20; // all energies for alpha and heavy nuclei
    }
    else if (ParticleName == "proton"){
        if(Energy <= 2.){
            factor = 1;
        }else{
            factor = 5;
        }
    }
    else if (ParticleName == "neutron"){

        // this for ICRP radiation factors
        if(Energy == 111   ){factor=17.49;}
        else if(Energy == 222   ){factor=16.67;}
        else if(Energy == 333   ){factor=16.99;}
        else if(Energy == 444   ){factor=17.12;}
        else if(Energy == 555   ){factor=17.07;}
        else if(Energy == 666   ){factor=17.37;}
        else if(Energy == 777   ){factor=16.47;}
        else if(Energy == 888   ){factor=16.88;}
        else if(Energy == 999   ){factor=16.85;}
        else if(Energy == 101010){factor=16.84;}
        else if(Energy == 111111){factor=16.92;}
        else if(Energy == 121212){factor=17.09;}
        else if(Energy == 131313){factor=17.27;}
        else if(Energy == 141414){factor=16.57;}
        else if(Energy == 151515){factor=16.56;}
        else if(Energy == 161616){factor=16.57;}
        else if(Energy == 171717){factor=16.57;}
        else if(Energy == 181818){factor=16.57;}
        else if(Energy == 191919){factor=16.57;}
        else if(Energy == 202020){factor=17.02;}
        else if(Energy == 212121){factor=17.02;}
        else if(Energy == 222222){factor=17.02;}
        else if(Energy == 232323){factor=17.02;}
        else if(Energy == 242424){factor=17.02;}
        else if(Energy == 252525){factor=17.02;}
        else if(Energy == 262626){factor=17.02;}
        else if(Energy == 272727){factor=17.02;}
        else if(Energy == 282828){factor=17.02;}
        // for general case
        else if(Energy < 0.01){
            factor = 5;
        }
        else if( 0.01 <= Energy && Energy <= 0.1){
            factor = 10;
        }
        else if( 0.1 < Energy && Energy <= 2){
            factor = 20;
        }
        else if( 2 < Energy && Energy <= 10){
            factor = 20;
        }
        else if( 20 < Energy){
            factor = 5;
        }
    }
    return factor;
}
double MainWindow::GenerateTissueFactor(QString OrganName){


    double factor = 0.12;
    if(TissueFactorMap["Others"] == 0.){
        TissueFactorMap["Others"] = 0;
        factor = TissueFactorMap["Others"];
    }

    bool isIn = false;
    for ( auto it = TissueFactorMap.begin(); it != TissueFactorMap.end(); ++it  ){

        if(it->first == OrganName){
            isIn = true;
            factor = TissueFactorMap[OrganName];
            break;
        }
    }

    if(isIn == false){
        TissueFactorMap[OrganName] = TissueFactorMap["Others"];
        factor = TissueFactorMap[OrganName];
    }

    return factor;
}
void MainWindow::on_pushButtonGetResultsNucl_clicked()
{
    //QFuture<void> future = QtConcurrent::run(this, &MainWindow::CalculateQuantitiesBasedOnICRPData);
    CalculateQuantitiesBasedOnICRPData();
}
void MainWindow::on_pushButtonReverseData_clicked()
{
    GenerateDataInTableView();
}
void MainWindow::GenerateDataInTableView(){ // energy in MeV

    QVector<QPair<double, QString>> Data;
    //QVector<QString> Data1;


    if(ui->pushButtonReverseData->text() == "Ascending Table, show Descending"){
        ui->pushButtonReverseData->setText("Descending Table, show Alphabetical");
        Data = RadionuclideValuePairVectorValuesDescending;

    }else if(ui->pushButtonReverseData->text() == "Descending Table, show Alphabetical"){
        ui->pushButtonReverseData->setText("Alphabetical Table, show Ascending");
        Data = RadionuclideValuePairVectorAlphabetical;
    }
    else if(ui->pushButtonReverseData->text() == "Alphabetical Table, show Ascending"){
        ui->pushButtonReverseData->setText("Ascending Table, show Descending");
        Data = RadionuclideValuePairVectorValuesAscending;
    }

    ui->tableWidgetForOneGraph->clear();
    ui->tableWidgetForOneGraph->setRowCount(0);
    ui->tableWidgetForOneGraph->setColumnCount(0);

    QStringList headers;

    headers.append(tr("RadioNuclide(T(s))"));
    headers.append(tr((ui->comboBoxQuantityNucl->currentText()+ "(" + ui->comboBoxEffDoseUnit->currentText()+")").toStdString().c_str()));

    ui->tableWidgetForOneGraph->setColumnCount(2);
    ui->tableWidgetForOneGraph->setShowGrid(true);
    ui->tableWidgetForOneGraph->setSelectionMode(QAbstractItemView::SingleSelection);
    ui->tableWidgetForOneGraph->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->tableWidgetForOneGraph->setHorizontalHeaderLabels(headers);
    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);

    //QTextStream(stdout) << "Table to show rows " << Data.size() << "\n";

    for(int row = 0 ; row < Data.size(); row++){
        // Insert row

        double valuee = ICRPRadioNuclideHalfLives[Data[row].second];
        QString Tval = QString::number(valuee)+"s";

        //QTextStream(stdout) << row << "  "<< Data[row] << " " << Data1[row] << "\n";

        ui->tableWidgetForOneGraph->insertRow(row);
        ui->tableWidgetForOneGraph->setItem(row,0, new QTableWidgetItem(Data[row].second+" ("+Tval+")"));
        ui->tableWidgetForOneGraph->setItem(row,1, new QTableWidgetItem(QString::number(Data[row].first/QuantitiesConversionFromDefault[ui->comboBoxQuantityNucl->currentText()][ui->comboBoxEffDoseUnit->currentText()])));
    }

    ui->tableWidgetForOneGraph->resizeColumnsToContents();

    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->viewport()->installEventFilter(this);

}
void MainWindow::on_checkBoxInterpolationType_stateChanged(int arg1)
{
    if(ui->checkBoxInterpolationType->isChecked()){
        ui->checkBoxInterpolationType->setText("Log");
    }else{
        ui->checkBoxInterpolationType->setText("Linear");
    }
}
void MainWindow::on_comboBoxQuantityNucl_currentTextChanged(const QString &arg1)
{
    ui->comboBoxEffDoseUnit->clear();
    ui->comboBoxEffDoseUnit->addItems(QuantitiesUnitsLists[arg1]);
    ui->checkBoxEffDoseLimNucl->setText(ui->comboBoxQuantityNucl->currentText()+" limits");
}
void MainWindow::on_comboBoxPhantom_currentTextChanged(const QString &arg1)
{
    ui->comboBoxSources->clear();
    ui->comboBoxSources->addItems(PhantomSourcesMap[ui->comboBoxPhantom->currentText()].toList());

    QStringList Organs1;
    for ( auto it = RegionParameterValueMap[ui->comboBoxPhantom->currentText()]["Mass"].begin(); it != RegionParameterValueMap[ui->comboBoxPhantom->currentText()]["Mass"].end(); ++it  ){Organs1.push_back(it.key());}
    ui->comboBoxTargets->clear();
    //ui->comboBoxTargets->addItem("all");
    ui->comboBoxTargets->addItems(Organs1);
}
// used with activated it not work with currentTextChanged
void MainWindow::on_comboBoxPhantom_textActivated(const QString &arg1)
{
    ConstructOtherTissuesSource();
}
void MainWindow::on_comboBoxSources_currentTextChanged(const QString &arg1)
{
    /*
    if(ui->comboBoxTotalorSourceOrCombNucl->currentIndex() == 4 || ui->comboBoxSourceOrTargetsForBiokinetics->currentIndex() == 3){
        bool IsIn = false;
        QStringList args = ui->lineEditSourcesToUse->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
        for(int dd=0; dd < args.size();dd++){ if(args[dd] == arg1){ IsIn = true;}}
        if(IsIn == false){

            if(arg1 != ""){
                QString nn = ui->lineEditSourcesToUse->text();
                ui->lineEditSourcesToUse->setText(nn + " " + arg1);
            }
        }
    }
    */
}
void MainWindow::on_comboBoxTargets_currentTextChanged(const QString &arg1)
{
    /*
    if(ui->comboBoxTotalorSourceOrCombNucl->currentIndex() == 4){
        bool IsIn = false;
        QStringList args = ui->lineEdittargetsToUse->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
        for(int dd=0; dd < args.size();dd++){ if(args[dd] == arg1){ IsIn = true;}}
        if(IsIn == false){

            if(arg1 != ""){
                QString nn = ui->lineEdittargetsToUse->text();
                ui->lineEdittargetsToUse->setText(nn + " " + arg1);
            }
        }
    }
    */
}
void MainWindow::on_comboBoxRadioPharmaceutiques_textActivated(const QString &arg1)
{
    ConstructOtherTissuesSource();
}
void MainWindow::on_comboBoxSourceOrTargetsForBiokinetics_currentTextChanged(const QString &arg1)
{
    if(ui->comboBoxSourceOrTargetsForBiokinetics->currentIndex() == 3){
        ui->lineEditSourcesToUse->setVisible(true);
    }else{
        ui->lineEditSourcesToUse->setVisible(false);
    }
}
void MainWindow::on_pushButton_selectSource_clicked()
{

    if(ui->comboBoxSources->count() == 0){
        QMessageBox::information(this, tr(""), "Canno't find source regions, please read SAFs from ICRP or DoseCalcs");
        return;
    }

    QDialog * d = new QDialog();
    QGridLayout* GraphLayout = new QGridLayout;
    QGridLayout* ElementLayou = new QGridLayout();

    int maxx = 14;
    int h = 0;
    int jj=0;

    QPushButton* CheckAllBtn = new QPushButton("Check-All"); QObject::connect(CheckAllBtn, SIGNAL(clicked()), this, SLOT(btnCheckALL()));
    QPushButton* UNCheckAllBtn = new QPushButton("Uncheck-All"); QObject::connect(UNCheckAllBtn, SIGNAL(clicked()), this, SLOT(btnUNCheckALL()));

    GraphLayout->addWidget(CheckAllBtn, jj,0,1,1);
    GraphLayout->addWidget(UNCheckAllBtn, jj,1,1,1);

    jj++;

    //QTextStream(stdout) << "ZZZ" << "\n" ;

    QLineEdit* regioname ;
    if (OpenSourceRegionsDialogFor == "UnknownBiokineticSource"){
        d->setWindowTitle("The source region \"" + SourceRegionNameInNewRegiondialog + "\" not found in source region list, construct it! ");
        regioname = new QLineEdit(); regioname->setText(SourceRegionNameInNewRegiondialog);
        regioname->setEnabled(false);
        regioname->setToolTip("The name is disabled and will be used as in biokinetics model of radiopharmaceutical, This source should be as apre-defined source region or will be considered as a combination of already simulated region shown below");
    }else{
        d->setWindowTitle("Construct new source region");
        regioname = new QLineEdit(); regioname->setPlaceholderText("RegionName");
        regioname->setToolTip("Add a name of new source that will be considered as a combination of already simulated region shown below");

    }

    GraphLayout->addWidget(regioname, jj,0,1,2);
    jj++;

    //QTextStream(stdout) << "ZZZ " << CurrentSources.size() << "\n" ;

    checkboxes.clear();
    SpinBoxes.clear();
    for(int ii = 0; ii < ui->comboBoxSources->count(); ii++){

        QCheckBox* CB = new QCheckBox(ui->comboBoxSources->itemText(ii));
        QDoubleSpinBox* SB = new QDoubleSpinBox();SB->setMaximum(1);SB->setMinimum(0);SB->setValue(1); SB->setSingleStep(0.0001);SB->setDecimals(4);
        SB->setToolTip("Add mass fraction of "+ui->comboBoxSources->itemText(ii)+" in the new source region");

        //QTextStream(stdout) << "AAA " << ui->comboBoxSources->itemText(ii) << "\n" ;

        bool isin = false;
        for (int dd = 0 ; dd < CurrentSources.size(); dd++) {
            //QTextStream(stdout) << "CurrentSources[dd] " << CurrentSources[dd] << "CurrentSourcesFractions[dd] " << CurrentSourcesFractions[dd] << "\n" ;

            if(CurrentSources[dd] == ui->comboBoxSources->itemText(ii)){
                CB->setChecked(true);
                if(CurrentSourcesFractions.size()-1 >= dd){
                    SB->setValue(CurrentSourcesFractions[dd]);
                    //QTextStream(stdout) << "CurrentSources[dd] " << CurrentSources[dd] << "CurrentSourcesFractions[dd] " << CurrentSourcesFractions[dd] << "\n" ;
                }

                break;
            }
        }

        //QTextStream(stdout) << "BBB" << "\n" ;

        GraphLayout->addWidget(CB, jj,h,1,1);
        h++;
        GraphLayout->addWidget(SB, jj,h,1,1);
        h++;

        //QTextStream(stdout) << "CCC" << "\n" ;

        checkboxes.push_back(CB);
        SpinBoxes.push_back(SB);

        //QTextStream(stdout) << "DDD " << jj << " " << "\n" ;

        if (h >= maxx){
            jj++;
            h = 0;
            continue;
        }
    }
    //QTextStream(stdout) << "ZZZ" << "\n" ;

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    jj++;
    GraphLayout->addWidget(buttonBox,jj,0,2,1);
    d->setLayout(GraphLayout);

    int result = d->exec();

    QString RegionNamee = "";

    if(result == QDialog::Accepted){

        CurrentSources.clear();
        CurrentSourcesFractions.clear();
        RegionFraction.clear();


        RegionNamee = regioname->text();
        for(int ii = 0; ii < checkboxes.size(); ii++){
            if (checkboxes[ii]->isChecked()){
                CurrentSources.push_back(checkboxes[ii]->text());
                CurrentSourcesFractions.push_back(SpinBoxes[ii]->value());
                RegionFraction[RegionNamee][checkboxes[ii]->text()] = SpinBoxes[ii]->value();
            }
        }

        //QTextStream(stdout) << "RegionNamee "<< RegionNamee << "== CurrentSources[0] " << CurrentSources[0] << " " << CurrentSourcesFractions[0] << "\n" ;

        if(CurrentSources.size() == 0 ){
            //CurrentSources.push_back(ui->comboBoxSources->currentText());
            return;
        }
        else if(CurrentSources.size() == 1 && (CurrentSources[0] == RegionNamee || CurrentSourcesFractions[0] == 1)){

            ui->comboBoxSources->setCurrentText(CurrentSources[0]);
            return;
        }
        else{

            if(RegionNamee != ""){

                bool isin = false;
                for (int dd = 0 ; dd < ui->comboBoxSources->count(); dd++) {
                    if(RegionNamee == ui->comboBoxSources->itemText(dd)){
                        isin = true;break;
                    }
                }
                if(isin == false){
                    ui->comboBoxSources->addItem(RegionNamee);
                    ui->comboBoxTargets->addItem(RegionNamee);
                    if (OpenSourceRegionsDialogFor == "UnknownBiokineticSource"){
                        ui->comboBoxSources->setCurrentText(RegionNamee);
                    }
                    GenerateSAFFromNewSource(RegionNamee);
                }
            }
        }
    }
}
void MainWindow::on_pushButton_selectTarget_clicked()
{
    if(ui->comboBoxTargets->count() == 0){
        QMessageBox::information(this, tr(""), "Canno't find target regions, please read SAFs from ICRP or DoseCalcs");
        return;
    }
    QDialog * d = new QDialog();
    QGridLayout* GraphLayout = new QGridLayout;

    d->setWindowTitle("Construct new target region");

    int maxx = 14;
    int h = 0;
    int jj=0;

    QPushButton* CheckAllBtn = new QPushButton("Check-All"); QObject::connect(CheckAllBtn, SIGNAL(clicked()), this, SLOT(btnCheckALL()));
    QPushButton* UNCheckAllBtn = new QPushButton("Uncheck-All"); QObject::connect(UNCheckAllBtn, SIGNAL(clicked()), this, SLOT(btnUNCheckALL()));

    GraphLayout->addWidget(CheckAllBtn, jj,0,1,1);
    GraphLayout->addWidget(UNCheckAllBtn, jj,1,1,1);

    jj++;

    QLineEdit* regioname = new QLineEdit(); regioname->setPlaceholderText("RegionName");
    regioname->setToolTip("Add a name of new target region that will be considered as a combination of already simulated region shown below");

    GraphLayout->addWidget(regioname, jj,0,1,2);

    jj++;

    checkboxes.clear();
    SpinBoxes.clear();
    for(int ii = 0; ii < ui->comboBoxTargets->count(); ii++){

        QCheckBox* CB = new QCheckBox(ui->comboBoxTargets->itemText(ii));
        QDoubleSpinBox* SB = new QDoubleSpinBox();SB->setMaximum(1);SB->setMinimum(0);SB->setValue(1);SB->setSingleStep(0.0001);SB->setDecimals(4);
        SB->setToolTip("Add mass fraction of "+ui->comboBoxSources->itemText(ii)+" in the new target region");

        bool isin = false;
        for (int dd = 0 ; dd < CurrentTargets.size(); dd++) {
            if(CurrentTargets[dd] == ui->comboBoxTargets->itemText(ii)){
                CB->setChecked(true);
                SB->setValue(CurrentTargetsFractions[dd]);

                break;
            }
        }

        GraphLayout->addWidget(CB, jj,h,1,1);
        h++;
        GraphLayout->addWidget(SB, jj,h,1,1);
        h++;

        checkboxes.push_back(CB);
        SpinBoxes.push_back(SB);

        if (h >= maxx){
            jj++;
            h = 0;
            continue;
        }
    }

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    jj++;
    GraphLayout->addWidget(buttonBox,jj,0,2,1);
    d->setLayout(GraphLayout);

    int result = d->exec();

    QString RegionNamee = "";
    if(result == QDialog::Accepted){

        CurrentTargets.clear();
        CurrentTargetsFractions.clear();
        RegionFraction.clear();

        RegionNamee = regioname->text();
        for(int ii = 0; ii < checkboxes.size(); ii++){
            if (checkboxes[ii]->isChecked()){
                CurrentTargets.push_back(checkboxes[ii]->text());
                CurrentTargetsFractions.push_back(SpinBoxes[ii]->value());
                RegionFraction[RegionNamee][checkboxes[ii]->text()] = SpinBoxes[ii]->value();
            }
        }
        if(CurrentTargets.size() == 0 ){
            //CurrentTargets.push_back(ui->comboBoxTargets->currentText());
            return;
        }
        else if(CurrentTargets.size() == 1 && (CurrentTargets[0] == RegionNamee || CurrentTargetsFractions[0] == 1)){
            ui->comboBoxTargets->setCurrentText(CurrentTargets[0]);
            return;
        }
        else{

            if(RegionNamee != ""){

                bool isin = false;
                for (int dd = 0 ; dd < ui->comboBoxTargets->count(); dd++) {
                    if(RegionNamee == ui->comboBoxTargets->itemText(dd)){
                        isin = true;break;
                    }
                }
                if(isin == false){
                    ui->comboBoxTargets->addItem(RegionNamee);
                    ui->comboBoxTargets->setCurrentText(RegionNamee);

                    GenerateSAFFromNewTarget(RegionNamee);
                }
            }
        }
    }
}
void MainWindow::btnCheckALL(){
    for(int ii = 0; ii < checkboxes.size(); ii++){
        checkboxes[ii]->setChecked(true);
    }
}
void MainWindow::btnUNCheckALL(){
    for(int ii = 0; ii < checkboxes.size(); ii++){
        checkboxes[ii]->setChecked(false);
    }
}
void MainWindow::on_checkBoxHalflivesNucl_clicked(bool checked)
{
    if(checked == true){
        ui->doubleSpinBoxMinHalflives->setEnabled(true);
        ui->doubleSpinBoxMaxHalflives->setEnabled(true);
        ui->comboBoxPeriodUnit->setEnabled(true);
    }else{
        ui->doubleSpinBoxMinHalflives->setEnabled(false);
        ui->doubleSpinBoxMaxHalflives->setEnabled(false);
        ui->comboBoxPeriodUnit->setEnabled(false);
    }
}
void MainWindow::on_checkBoxEffDoseLimNucl_clicked(bool checked)
{
    if(checked == true){
        ui->doubleSpinBoxMaxEffDose->setEnabled(true);
        ui->doubleSpinBoxMinEffDose->setEnabled(true);
    }else{
        ui->doubleSpinBoxMaxEffDose->setEnabled(false);
        ui->doubleSpinBoxMinEffDose->setEnabled(false);
    }
}

void MainWindow::on_pushButtonAddRowToTable_clicked()
{


    if(ui->spinBoxNumbeOfRows->value() == 0){return;}
    int inirow = ui->tableWidgetForOneGraph->rowCount();
    for(int ii = 0; ii < ui->spinBoxNumbeOfRows->value(); ii++){
        ui->tableWidgetForOneGraph->insertRow(inirow+ii);
        for(int jj = 0; jj < ui->tableWidgetForOneGraph->columnCount(); jj++){
            ui->tableWidgetForOneGraph->setItem(inirow+ii,jj, new QTableWidgetItem(""));
        }
    }


}

void MainWindow::on_pushButton_SaveADEDResultsInResultsData_clicked()
{
    QString label = ui->tableWidgetForOneGraph->horizontalHeaderItem(1)->text();
    if(ui->tableWidgetForOneGraph->horizontalHeaderItem(0)->text() == "RadioNuclide(T(s))" || label == "Mass (kg)" || label == "WT" ){
        return;
    }

    int SZ = 17;
    int VarSZ = 8;

    QString RadiTracerBiokineticsData = " - ";
    QString Quantity_NAME = ui->comboBoxQuantityNucl->currentText();
    QString Quantity_NAME_UNIT = ui->comboBoxQuantityNucl->currentText()+"("+ ui->comboBoxEffDoseUnit->currentText()+")";
    //double convfac = QuantitiesConversionFromDefault[ui->comboBoxQuantityNucl->currentText()][ui->comboBoxEffDoseUnit->currentText()];
    QString RadioTracer_NAME = RadiotracerradionucleidMap[ui->comboBoxRadioPharmaceutiques->currentText()];
    QString Geometry_NAME = ui->comboBoxPhantom->currentText();
    QString Particle_NAME;
    QString Source_NAME;
    QString Target_NAME;

    std::ostringstream filname ;
    filname << UserCurrentResultsDirPath.toStdString() << "/" << ResultFileName.toStdString();
    QString fileNameString = filname.str().c_str();
    std::ofstream file(fileNameString.toStdString(), std::ios_base::app);

    if(file.is_open()){

        std::ostringstream OutS;
        OutS << "****** "
             << ui->comboBoxQuantityNucl->currentText().toStdString() << " "
             << "IntakeIntoBody "
             << "RadioTracer "
             << RadioTracer_NAME.toStdString() << " "
             << Geometry_NAME.toStdString() <<  " "
             << "Physics "
             << " TotalEmittedEnergy! "
             << RadiTracerBiokineticsData.toStdString() << " "
             << "InjectedActivity(Bq)=" << ui->doubleSpinBoxAdministeredActivity->value()/QuantitiesConversionFromDefault["A"][ui->comboBoxActivityAdministered->currentText()] << " " << Quantity_NAME_UNIT.toStdString() << "_Total"<< "=!" <<" " ;

        OutS << "\n";

        QString headerText= OutS.str().c_str();

        file << headerText.toStdString();

        file << std::setw(24) << std::left << "# Volume" << " "
             << std::setw(SZ) << std::left << Quantity_NAME_UNIT.toStdString() << " "
             << std::setw(SZ) << std::left << "SDev"  << " "
             << std::setw(SZ) << std::left << "Rel_SDev(%)" << " "
             << std::setw(SZ) << std::left << "Values Num" << " "
             << std::setw(SZ) << std::left << "Mass[kg]" << " "
             << std::setw(SZ) << std::left << "Volume[cm3]" << " "
             << std::setw(SZ) << std::left << "Density[g/cm3] " << " \n";

        int rowNum = ui->tableWidgetForOneGraph->rowCount();
        QMap<QString,double> RegionValueMap;

        for(int ii = 0; ii < rowNum; ii++){

            Target_NAME = ui->tableWidgetForOneGraph->item( ii, 0 )->text();Target_NAME.remove(" ");
            double value = ui->tableWidgetForOneGraph->item( ii, 1 )->text().toDouble();
            double sdv = 1;
            unsigned long long int sm = 1;
            double rsdv = ((sdv*sm)/100)*100;

            file << std::setw(24) << std::left << Target_NAME.toStdString() << " "
                 << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << value << " "
                 << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << sdv << " "
                 << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << rsdv << " "
                 << std::resetiosflags(file.basefield) << std::resetiosflags( file.floatfield) << std::resetiosflags( file.flags())
                 << std::setw(SZ) << std::left << sm << " "
                 << std::setw(SZ) << std::left << RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME] << " "
                 << std::setw(SZ) << std::left << RegionParameterValueMap[Geometry_NAME]["Volume"][Target_NAME] << " "
                 << std::setw(SZ) << std::left << RegionParameterValueMap[Geometry_NAME]["Density"][Target_NAME] << " \n";
        }

        file << "* ----------------------------------------------------------------------------------------------------------------------------------------------\n";

        file.close();
    }
}
void MainWindow::on_pushButtonSaveTableInPDF_clicked()
{
/*
    QTableWidget* tbl ;
    tbl = ui->tableWidgetForOneGraph;

    const int columns = tbl->columnCount() ;
    const int rows = tbl->rowCount();
    QTextDocument doc;
    QTextCursor cursor(&doc);
    QTextTableFormat tableFormat;
    tableFormat.setHeaderRowCount(1);
    tableFormat.setAlignment(Qt::AlignHCenter);
    tableFormat.setCellPadding(3);
    tableFormat.setCellSpacing(3);
    tableFormat.setBorder(1);
    tableFormat.setBorderBrush(QBrush(Qt::SolidPattern));
    tableFormat.clearColumnWidthConstraints();

    QTextTable *textTable = cursor.insertTable(rows + 1, columns, tableFormat);
    QTextCharFormat tableHeaderFormat;
    tableHeaderFormat.setBackground(QColor("#DADADA"));
    //QTextStream(stdout) << "coloumn" << columns << "\n";

    for (int i = 0; i < columns; i++) {
        //QTextStream(stdout) << "i" << i << "\n";
        QTextTableCell cell = textTable->cellAt(0, i);
        cell.setFormat(tableHeaderFormat);
        QTextCursor cellCursor = cell.firstCursorPosition();
        cellCursor.insertText(tbl->horizontalHeaderItem(i)->data(Qt::DisplayRole).toString());
    }
    //QTextStream(stdout) << "---------------- 3 ------------------" << "\n";

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {

            QTextTableCell cell = textTable->cellAt(i+1, j);
            QTextCursor cellCursor = cell.firstCursorPosition();
            cellCursor.insertText(tbl->item(i , j)->text());
        }
    }

    //QTextStream(stdout) << "---------------- 4 ------------------" << "\n";

    cursor.movePosition(QTextCursor::End);
#if QT_PRINTSUPPORT_LIB
    QPrinter printer(QPrinter::PrinterResolution);
    printer.setOutputFormat(QPrinter::PdfFormat);
    printer.setPaperSize(QPrinter::A4);
    printer.setOrientation(QPrinter::Landscape);
#endif

    QString FileName = tr((UserCurrentResultsDirPath+"/" +GraphsOutDirName+"/" + "RadionuclidesClassification").toStdString().c_str());

    bool ok;
    QString newLabel = QInputDialog::getText(this, "Enter the File Name", "New File Name:", QLineEdit::Normal, UserCurrentResultsDirPath+"/" +GraphsOutDirName, &ok);
    if (ok)
    {
        FileName = newLabel;
    }

#if QT_PRINTSUPPORT_LIB
    printer.setOutputFileName( FileName );

    doc.setDocumentMargin(0);
    doc.setTextWidth(10);

    doc.print(&printer);
    //QTextStream(stdout) << "---------------- 5 ------------------" << "\n";
#endif
*/

    int rowNum = ui->tableWidgetForOneGraph->rowCount();
    int colNum = ui->tableWidgetForOneGraph->columnCount();

    if(!QFile::exists(UserCurrentResultsDirPath+"/"+GraphsOutDirName)){
        QDir* dir = new QDir(UserCurrentResultsDirPath+"/"+GraphsOutDirName);
        dir->mkdir(UserCurrentResultsDirPath+"/"+GraphsOutDirName);
    }

    QString fileName = tr((UserCurrentResultsDirPath+"/" +GraphsOutDirName+"/" + "Data.csv").toStdString().c_str());

    QString TextToWrite="";
    TextToWrite += "\n\n\n";
    TextToWrite += "------------------------------------------------------------ , \n";

    for(int ii = 0; ii < rowNum; ii++){
        for(int jj = 0; jj < colNum; jj++){
             TextToWrite += ui->tableWidgetForOneGraph->item( ii, jj )->text() + ", ";
        }
        TextToWrite += "\n";
    }
    TextToWrite += "------------------------------------------------------------ , \n";

    fileManagerObject->WriteTextToFile(fileName , TextToWrite);

    ui->pushButtonSaveTableInPDF->setToolTip("The data in table are saved to file: "+fileName);
}

void MainWindow::on_pushButtonChangeMasses_clicked()
{
    ui->tableWidgetForOneGraph->clear();
    ui->tableWidgetForOneGraph->setRowCount(0);
    ui->tableWidgetForOneGraph->setColumnCount(0);

    QStringList headers;

    headers.append(tr("Region Name"));
    headers.append(tr("Mass (kg)"));

    ui->tableWidgetForOneGraph->setColumnCount(2);
    ui->tableWidgetForOneGraph->setShowGrid(true);
    ui->tableWidgetForOneGraph->setSelectionMode(QAbstractItemView::SingleSelection);
    ui->tableWidgetForOneGraph->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->tableWidgetForOneGraph->setHorizontalHeaderLabels(headers);
    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);

    //QTextStream(stdout) << "Table to show rows " << Data.size() << "\n";

    int row = 0;
    for ( auto it = RegionParameterValueMap[ui->comboBoxPhantom->currentText()]["Mass"].begin(); it != RegionParameterValueMap[ui->comboBoxPhantom->currentText()]["Mass"].end(); ++it  ){

        ui->tableWidgetForOneGraph->insertRow(row);
        ui->tableWidgetForOneGraph->setItem(row,0, new QTableWidgetItem(it.key()));
        ui->tableWidgetForOneGraph->setItem(row,1, new QTableWidgetItem(QString::number(it.value())));
        row++;
    }

    ui->tableWidgetForOneGraph->resizeColumnsToContents();

    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->viewport()->installEventFilter(this);

}
void MainWindow::on_pushButtonSaveMasses_clicked()
{
    QString Quantity_NAME = "SAF";
    int rowNum = ui->tableWidgetForOneGraph->rowCount();
    QMap<QString,double> TargetMassMap = RegionParameterValueMap[ui->comboBoxPhantom->currentText()]["Mass"];

    for(int ii = 0; ii < rowNum; ii++){

    /*QTextStream(stdout) << ii <<" For Phantom with " << ui->comboBoxPhantom->currentText()
                            <<" organ " << ui->tableWidgetForOneGraph->item( ii, 0 )->text()
                           <<" mass " << ui->tableWidgetForOneGraph->item( ii, 1 )->text().toDouble()
                          << "\n";
    */

        TargetMassMap[ui->tableWidgetForOneGraph->item( ii, 0 )->text()] = ui->tableWidgetForOneGraph->item( ii, 1 )->text().toDouble();
    }

    QMap<QString,QMap<QString,QMap<double,QMap<QString,QMap<QString,double>>>>> SAFMap = ICRPSAFs[Quantity_NAME];

    ui->progressBarReadingCalcData->setRange(0, 100);
    ui->progressBarReadingCalcData->setValue(0);
    ui->progressBarReadingCalcData->show();
    double radincc = 0;
    double percent = 0;


    for ( auto it2 = SAFMap[ui->comboBoxPhantom->currentText()].begin(); it2 != SAFMap[ui->comboBoxPhantom->currentText()].end(); ++it2  ){
        QString Particle_NAME = it2.key();
        //            QTextStream(stdout) << "---Particle_NAME " << Particle_NAME <<"\n";

        for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
            double Energy_Val = it3.key();
            //                QTextStream(stdout) << "----Energy_Val " << Energy_Val <<"\n";

            for ( auto DD = it3.value().begin(); DD != it3.value().end(); ++DD  ){
                QString Source_NAME  = DD.key();
                //                    QTextStream(stdout) << "-----Source_NAME " << Source_NAME <<"\n";

                for ( auto CC = DD.value().begin(); CC != DD.value().end(); ++CC  ){
                    QString  Target_NAME  = CC.key();
                    double SAFValue  = CC.value();

                    percent = (radincc/ICRPRadioNuclideData.size())*100;
                    ui->progressBarReadingCalcData->setValue(percent);
                    radincc++;

                    ICRPSAFs[Quantity_NAME][ui->comboBoxPhantom->currentText()][Particle_NAME][Energy_Val][Source_NAME][Target_NAME] = SAFValue * RegionParameterValueMap[ui->comboBoxPhantom->currentText()]["Mass"][Target_NAME]/TargetMassMap[Target_NAME];

                    //if(Source_NAME == "Liver" && Target_NAME == "Liver"){
                    //    QTextStream(stdout) << " Geometry_NAME " << ui->comboBoxPhantom->currentText() << " Particle_NAME " << Particle_NAME << " Energy_Val " << Energy_Val << " Source_NAME " << Source_NAME << " Target_NAME " << Target_NAME << " SAFValue " << SAFValue << " " << RegionParameterValueMap[ui->comboBoxPhantom->currentText()]["Mass"][Target_NAME]/TargetMassMap[Target_NAME] << " "<< ICRPSAFs["SAF"][ui->comboBoxPhantom->currentText()][Particle_NAME][Energy_Val][Source_NAME][Target_NAME]  << "\n";
                    //}
                }
            }
        }
    }

    RegionParameterValueMap[ui->comboBoxPhantom->currentText()]["Mass"] = TargetMassMap;

    ui->progressBarReadingCalcData->setValue(100);

    ui->tableWidgetForOneGraph->clear();
    ui->tableWidgetForOneGraph->setRowCount(0);
    ui->tableWidgetForOneGraph->setColumnCount(0);
}

void MainWindow::on_pushButton_TissueWeithingFactor_clicked()
{
    ui->tableWidgetForOneGraph->clear();
    ui->tableWidgetForOneGraph->setRowCount(0);
    ui->tableWidgetForOneGraph->setColumnCount(0);

    QStringList headers;

    headers.append(tr("Region Name"));
    headers.append(tr("WT"));

    ui->tableWidgetForOneGraph->setColumnCount(2);
    ui->tableWidgetForOneGraph->setShowGrid(true);
    ui->tableWidgetForOneGraph->setSelectionMode(QAbstractItemView::SingleSelection);
    ui->tableWidgetForOneGraph->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->tableWidgetForOneGraph->setHorizontalHeaderLabels(headers);
    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);

    double totalWT = 0;
    //for ( auto it = TissueFactorMap.begin(); it != TissueFactorMap.end(); ++it  ){
    for(int ii = 0; ii < ui->comboBoxTargets->count(); ii++){

        ui->tableWidgetForOneGraph->insertRow(ii);
        ui->tableWidgetForOneGraph->setItem(ii,0, new QTableWidgetItem(ui->comboBoxTargets->itemText(ii)));
        ui->tableWidgetForOneGraph->setItem(ii,1, new QTableWidgetItem(QString::number(TissueFactorMap[ui->comboBoxTargets->itemText(ii)])));
        totalWT += TissueFactorMap[ui->comboBoxTargets->itemText(ii)];
    }

    ui->tableWidgetForOneGraph->insertRow(ui->comboBoxTargets->count());
    ui->tableWidgetForOneGraph->setItem(ui->comboBoxTargets->count(),0, new QTableWidgetItem("WT_Total"));
    ui->tableWidgetForOneGraph->setItem(ui->comboBoxTargets->count(),1, new QTableWidgetItem(QString::number(totalWT)));


    ui->tableWidgetForOneGraph->resizeColumnsToContents();

    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->viewport()->installEventFilter(this);
}
void MainWindow::on_pushButton_SaveWT_clicked()
{

    if(ui->pushButtonReadUSERData->text() == "Read User SAFs *"){

    }else if(ui->pushButton_ReadDoseCalcsSAFs->text() == "Read DoseCalcs SAFs *"){
        QMessageBox::warning(this, tr(""), "You canno't edit the WT (tissue-weigthing factors) when you use DoseCalcs or ICRP133 data, you can do it when you choose User SAFs" );
        return;
    }else if(ui->pushButtonReadICRPData->text() == "Read ICRP133 SAFs *"){
        QMessageBox::warning(this, tr(""), "You canno't edit the WT (tissue-weigthing factors) when you use DoseCalcs or ICRP133 data, you can do it when you choose User SAFs" );
        return;
    }

    int rowNum = ui->tableWidgetForOneGraph->rowCount();
    std::map<QString,double> TargetMassMap = TissueFactorMap;

    for(int ii = 0; ii < rowNum; ii++){
        if(ui->tableWidgetForOneGraph->item( ii, 0 )->text() == "WT_Total"){
            continue;
        }
        TargetMassMap[ui->tableWidgetForOneGraph->item( ii, 0 )->text()] = ui->tableWidgetForOneGraph->item( ii, 1 )->text().toDouble();
    }

    ui->progressBarReadingCalcData->setRange(0, 100);
    ui->progressBarReadingCalcData->setValue(0);
    ui->progressBarReadingCalcData->show();

    TissueFactorMap = TargetMassMap ;

    ui->progressBarReadingCalcData->setValue(100);

    ui->tableWidgetForOneGraph->clear();
    ui->tableWidgetForOneGraph->setRowCount(0);
    ui->tableWidgetForOneGraph->setColumnCount(0);
    ConstructOtherTissuesSource();
}

void MainWindow::on_pushButtonAddBiokineticModel_clicked()
{
    QString Xtitle = "Source Region"; QString Ytitle = "Residence Time (h)";

    QDialog * d = new QDialog(); d->setWindowTitle("Biokinetic Table inputs for " + ui->comboBoxRadioPharmaceutiques->currentText());
    QVBoxLayout * vbox = new QVBoxLayout();

    //QStringList RegionVars=;
    QLineEdit * RadiPharmaceutic = new QLineEdit(); RadiPharmaceutic->setText(ui->comboBoxRadioPharmaceutiques->currentText()); vbox->addWidget(RadiPharmaceutic);

    QComboBox * Radionucleids = new QComboBox();
    for ( auto it = ICRPRadioNuclideHalfLives.begin(); it != ICRPRadioNuclideHalfLives.end(); ++it  ){
        Radionucleids->addItem(it.key());
    }
    vbox->addWidget(Radionucleids);

    QLineEdit * lineEditXTitle = new QLineEdit(); lineEditXTitle->setText(Xtitle); vbox->addWidget(lineEditXTitle);
    QLineEdit * lineEditYTitle = new QLineEdit(); lineEditYTitle->setText(Ytitle); vbox->addWidget(lineEditYTitle);
    QLineEdit * lineEditRowsNumber = new QLineEdit(); lineEditRowsNumber->setText("10"); vbox->addWidget(lineEditRowsNumber);

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);

    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    vbox->addWidget(buttonBox);

    d->setLayout(vbox);

    int RowNum = 10;

    int result = d->exec();

    if(result == QDialog::Accepted)
    {
        RowNum = lineEditRowsNumber->text().toInt();

        Xtitle = lineEditXTitle->text();
        Ytitle = lineEditYTitle->text();
        RadiotracerradionucleidMap[RadiPharmaceutic->text()] = Radionucleids->currentText();

        ui->tableWidgetForOneGraph->clear();
        ui->tableWidgetForOneGraph->setRowCount(0);
        ui->tableWidgetForOneGraph->setColumnCount(0);

        QStringList headers;

        headers.append(Xtitle);
        headers.append(Ytitle);

        ui->tableWidgetForOneGraph->setColumnCount(2);
        ui->tableWidgetForOneGraph->setShowGrid(true);
        ui->tableWidgetForOneGraph->setSelectionMode(QAbstractItemView::SingleSelection);
        ui->tableWidgetForOneGraph->setSelectionBehavior(QAbstractItemView::SelectRows);
        ui->tableWidgetForOneGraph->setHorizontalHeaderLabels(headers);
        ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);

        ui->tableWidgetForOneGraph->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);

        bool isin = false;
        for (int dd = 0 ; dd < ui->comboBoxRadioPharmaceutiques->count(); dd++) {
            if(RadiPharmaceutic->text() == ui->comboBoxRadioPharmaceutiques->itemText(dd)){
                isin = true;break;
            }
        }
        if(isin == false){
            ui->comboBoxRadioPharmaceutiques->addItem(RadiPharmaceutic->text());
        }

        ui->comboBoxRadioPharmaceutiques->setCurrentText(RadiPharmaceutic->text());

        for(int row = 0 ; row < RowNum; row++){
            // Insert row

            ui->tableWidgetForOneGraph->insertRow(row);

            ui->tableWidgetForOneGraph->setItem(row,0, new QTableWidgetItem(""));
            ui->tableWidgetForOneGraph->setItem(row,1, new QTableWidgetItem(""));
        }

        ui->tableWidgetForOneGraph->resizeColumnsToContents();
        //ui->tableWidgetForOneGraph->horizontalHeader()->viewport()->installEventFilter(this);
    }

}
void MainWindow::on_pushButtonSaveBiokinetikModelData_clicked()
{
    int rowNum = ui->tableWidgetForOneGraph->rowCount();

    RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadiotracerradionucleidMap[ui->comboBoxRadioPharmaceutiques->currentText()]].clear();

    QString line = "" ;
    line += ui->comboBoxRadioPharmaceutiques->currentText() + " " + RadiotracerradionucleidMap[ui->comboBoxRadioPharmaceutiques->currentText()] + " ";
    for(int ii = 0; ii < rowNum; ii++){

        if(ui->tableWidgetForOneGraph->item( ii, 1 )->text().toDouble() != 0.){
            RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadiotracerradionucleidMap[ui->comboBoxRadioPharmaceutiques->currentText()]][ui->tableWidgetForOneGraph->item( ii, 0 )->text()] = ui->tableWidgetForOneGraph->item( ii, 1 )->text().toDouble();
            line += ui->tableWidgetForOneGraph->item( ii, 0 )->text() + " " + ui->tableWidgetForOneGraph->item( ii, 1 )->text() + " ";
        }

        /*QTextStream(stdout) << ii <<" radiotracer with " << ui->comboBoxRadioPharmaceutiques->currentText()
                            <<" source " << ui->tableWidgetForOneGraph->item( ii, 0 )->text()
                           <<" value " << ui->tableWidgetForOneGraph->item( ii, 1 )->text().toDouble()
                           << "\n";
        */
    }

    if (QMessageBox::Yes == QMessageBox::question(this, tr("Save Data ?"), tr("Save Radiopharmaceuticals Biokinetics to File"))){
        fileManagerObject->WriteTextToFile(ICRPDATAPath+"/ICRP128RadioPharmaceuticalsData",fileManagerObject->ReadTextFromFileInOneString(ICRPDATAPath+"/ICRP128RadioPharmaceuticalsData")+"\n"+line);
    }


    ui->tableWidgetForOneGraph->clear();
    ui->tableWidgetForOneGraph->setRowCount(0);
    ui->tableWidgetForOneGraph->setColumnCount(0);

}
void MainWindow::on_pushButtonShowBiokineticData_clicked()
{
    ui->tableWidgetForOneGraph->clear();
    ui->tableWidgetForOneGraph->setRowCount(0);
    ui->tableWidgetForOneGraph->setColumnCount(0);

    QStringList headers;

    QString Xtitle = "Source Region"; QString Ytitle = "Residence Time (h)";

    headers.append(Xtitle);
    headers.append(Ytitle);

    ui->tableWidgetForOneGraph->setColumnCount(2);
    ui->tableWidgetForOneGraph->setShowGrid(true);
    ui->tableWidgetForOneGraph->setSelectionMode(QAbstractItemView::SingleSelection);
    ui->tableWidgetForOneGraph->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->tableWidgetForOneGraph->setHorizontalHeaderLabels(headers);
    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);

    //QTextStream(stdout) << "Table to show rows " << "\n";

    int row = 0;
    for ( auto it = RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadiotracerradionucleidMap[ui->comboBoxRadioPharmaceutiques->currentText()]].begin(); it != RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadiotracerradionucleidMap[ui->comboBoxRadioPharmaceutiques->currentText()]].end(); ++it  ){

            ui->tableWidgetForOneGraph->insertRow(row);
            ui->tableWidgetForOneGraph->setItem(row,0, new QTableWidgetItem(it.key()));
            ui->tableWidgetForOneGraph->setItem(row,1, new QTableWidgetItem(QString::number(it.value())));
            row++;
    }

    ui->tableWidgetForOneGraph->resizeColumnsToContents();

    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->viewport()->installEventFilter(this);

}
void MainWindow::on_pushButtonGenerateQuantitiesWithBiokineticData_clicked()
{
    if(ui->doubleSpinBoxAdministeredActivity->value() == 0.){
        QMessageBox::information(this, tr(""), "Please add the administered activity and choose the corresponding unit");
        return;
    }

    QString Quantity_NAME = ui->comboBoxQuantityNucl->currentText();
    double convfac = QuantitiesConversionFromDefault[ui->comboBoxQuantityNucl->currentText()][ui->comboBoxEffDoseUnit->currentText()];
    QString RadioTracer_NAME = RadiotracerradionucleidMap[ui->comboBoxRadioPharmaceutiques->currentText()];
    QString Geometry_NAME = ui->comboBoxPhantom->currentText();
    QString Particle_NAME;
    QString Source_NAME;
    QString Target_NAME;
    double Energy_Val;

    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>>> DATAMap; // Quantity, Geometry, Radionuclide, organ, value

    QVector<QString> Targets; Targets.push_back("All Targets");
    QVector<QString> Particles; Particles.push_back("All Particles");

    // For the radionuclides we calculate other AE, AF, S values, H, E for each geometry, radionuclide, source target energy

    // check if the specific sources are in the radiotracer biokinetics
    //if(ui->comboBoxSourceOrTargetsForBiokinetics->currentIndex() == 1 || ui->comboBoxSourceOrTargetsForBiokinetics->currentIndex() == 0){ // from a specifc sources
        /*
        QStringList srcs = ui->lineEditSourcesToUse->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
        QString nain = "";
        for(int b=0; b < srcs.size();b++){
            bool isin = false;
            for ( auto DD = RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME].begin(); DD != RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME].end(); ++DD  ){
                Source_NAME  = DD.key();
                if(srcs[b] == Source_NAME){
                    isin = true;
                    break;
                }
            }
            if(isin == false){
                nain += srcs[b] + " ";
            }
        }
        if(nain != ""){
            QMessageBox::information(this, tr("Warning"), "The specific sources ("+ nain +")are not defined in the biokinetics data of "+ ui->comboBoxRadioPharmaceutiques->currentText());
        }
        */
    //}

    //QTextStream(stdout) << " 1 " <<"\n";

    // check if the source region name in the biokinetics is defined in the SAF data file with that name
    // modify the name to be used in all session until it is modified by edit biokinetics data button
    QMap<QString,QMap<QString,QMap<QString,double>>> RadioTracerSourceOrganResidenceTime1 = RadioTracerSourceOrganResidenceTime;
    for ( auto DD = RadioTracerSourceOrganResidenceTime1[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME].begin(); DD != RadioTracerSourceOrganResidenceTime1[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME].end(); ++DD  ){
        Source_NAME = DD.key();

        bool isin = false;
        for ( auto CC = SourceParticleEnergyValues.begin(); CC != SourceParticleEnergyValues.end(); ++CC  ){
            if(CC->first == Source_NAME){
                //QTextStream(stdout) << " Source_NAME " << Source_NAME << " CC->first " << CC->first <<"\n";
                isin = true; break;
            }
        }
        //QTextStream(stdout) << Source_NAME <<"\n";


        if(isin == false){

            OpenSourceRegionsDialogFor = "UnknownBiokineticSource";
            SourceRegionNameInNewRegiondialog = Source_NAME;
            //QTextStream(stdout) << "on_pushButton_selectSource_clicked" <<"\n";
            on_pushButton_selectSource_clicked();
            //QTextStream(stdout) << "on_pushButton_selectSource_clicked" <<"\n";

            if(CurrentSources.size() == 1 ){
                //QTextStream(stdout) << "Source_NAME " << Source_NAME << " Changed to " << CurrentSources[0]  <<"\n";

                RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME][Source_NAME] = NULL;
                RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME].remove(Source_NAME);
                RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME][CurrentSources[0]] = RadioTracerSourceOrganResidenceTime1[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME][Source_NAME];
            }
            OpenSourceRegionsDialogFor = "";

        }
    }

    //QTextStream(stdout) << " 1 " <<"\n";

    for ( auto it2 = ICRPRadioNuclideData[RadioTracer_NAME].begin(); it2 != ICRPRadioNuclideData[RadioTracer_NAME].end(); ++it2  ){
        Particle_NAME = it2.key();
        QString Partial_Particle_NAME = Particle_NAME;
        if(Particle_NAME == "e+"){
            Partial_Particle_NAME = "e-";
            //QTextStream(stdout) << "Particle_NAME " << Particle_NAME << " Partial_Particle_NAME " << Partial_Particle_NAME <<"\n";
        }

        bool isin = false;
        for(int b=0; b < Particles.size();b++){ if(Particles[b] == Particle_NAME){isin = true; break;}}
        if(isin == false){Particles.push_back(Particle_NAME);}

        for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
            Energy_Val = it3.key();
            // fOR NEUTRON the yield is registered for reel mono energy of neutron in .rad file
            double RadiationPerCent = it3.value();
            //ICRPRadioNuclideData[RadioTracer_NAME][Particle_NAME][Energy_Val];
            if(Particle_NAME == "neutron"){
               Energy_Val = EnergyIDRadionuclideForNeutronSAF[RadioTracer_NAME];
            }

            //QTextStream(stdout) << "Energy_Val " << Energy_Val << " RadiationPerCent " << RadiationPerCent <<"\n";

            // this is to provide results for all sources that can user simulate without taking into account the RadioTracer source Organs

            for ( auto DD = RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME].begin(); DD != RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME].end(); ++DD  ){
                Source_NAME  = DD.key();

                double CumulatedActivityInSource = (ui->doubleSpinBoxAdministeredActivity->value()/QuantitiesConversionFromDefault["A"][ui->comboBoxActivityAdministered->currentText()])*(RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME][Source_NAME]/QuantitiesConversionFromDefault["T"]["h"]);

                //QTextStream(stdout) << " A " << ui->doubleSpinBoxAdministeredActivity->value()/QuantitiesConversionFromDefault["A"][ui->comboBoxActivityAdministered->currentText()]
                  //      << " Tho of source " << Source_NAME << " " << RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME][Source_NAME]/QuantitiesConversionFromDefault["T"]["h"]
                    //    << " CumulatedActivityInSource " << CumulatedActivityInSource << " " << ui->doubleSpinBoxAdministeredActivity->value()/QuantitiesConversionFromDefault["A"][ui->comboBoxActivityAdministered->currentText()]<< " " << RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME][Source_NAME]/QuantitiesConversionFromDefault["T"]["h"] <<"\n";

                if( !__isnan(CumulatedActivityInSource) && !__isinf(CumulatedActivityInSource) && CumulatedActivityInSource != 0 && CumulatedActivityInSource != NULL){}
                else{continue;}

                for ( auto CC = RegionParameterValueMap[Geometry_NAME]["Mass"].begin(); CC != RegionParameterValueMap[Geometry_NAME]["Mass"].end(); ++CC  ){
                    Target_NAME  = CC.key();

                    // here we considere that energy emitted from FF is absorbed locally in the source region
                    if(Particle_NAME == "FF" && Target_NAME == Source_NAME){
                        ICRPSAFs["SAF"][Geometry_NAME][Particle_NAME][Energy_Val][Source_NAME][Source_NAME] = 1/RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME];
                        //QTextStream(stdout) << "Particle_NAME " << Particle_NAME << " ICRPSAFs " << ICRPSAFs["AE"][Geometry_NAME][Particle_NAME][Energy_Val][Source_NAME][Source_NAME] <<"\n";
                    }

                    double ccc = GenerateRadiotracerQuantitiesByInterpolationInDefaultUnitForBiokinetic(Partial_Particle_NAME, Energy_Val, Source_NAME, Target_NAME);;
                    //double ccc = (ICRPSAFs["AE"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]/RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME])*GenerateRadiationFactor(ParticleName,Energy);

                    if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){

                        bool isin = false;
                        for(int b=0; b < Targets.size();b++){ if(Targets[b] == Target_NAME){isin = true; break;}}
                        if(isin == false){Targets.push_back(Target_NAME);}

                        //if(Source_NAME == "Liver" && Target_NAME == "Liver"){
                            //QTextStream(stdout) << " AE=" << ICRPSAFs["AE"][Geometry_NAME][Partial_Particle_NAME][Energy_Val][Source_NAME][Target_NAME] << " CumulatedActivityInSource="  << CumulatedActivityInSource << " RadiationPerCent="  <<  RadiationPerCent << " ccc="  << ccc << " "<< Quantity_NAME << "="<< RadiationPerCent*CumulatedActivityInSource*ccc <<"\n";
                        //}

                        DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particle_NAME][Source_NAME][Target_NAME] +=
                                RadiationPerCent
                                *CumulatedActivityInSource
                                *ccc;

                        DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME]["All Particles"][Source_NAME][Target_NAME] +=
                                RadiationPerCent
                                *CumulatedActivityInSource
                                *ccc;

                        DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particle_NAME]["All Sources"][Target_NAME] +=
                                RadiationPerCent
                                *CumulatedActivityInSource
                                *ccc;

                        DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particle_NAME][Source_NAME]["All Targets"] +=
                                RadiationPerCent
                                *CumulatedActivityInSource
                                *ccc;

                        DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME]["All Particles"][Source_NAME]["All Targets"] +=
                                RadiationPerCent
                                *CumulatedActivityInSource
                                *ccc;

                        DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME]["All Particles"]["All Sources"][Target_NAME] +=
                                RadiationPerCent
                                *CumulatedActivityInSource
                                *ccc;

                        DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particle_NAME]["All Sources"]["All Targets"] +=
                                RadiationPerCent
                                *CumulatedActivityInSource
                                *ccc;

                        DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME]["All Particles"]["All Sources"]["All Targets"] +=
                                RadiationPerCent
                                *CumulatedActivityInSource
                                *ccc;

                        /*
                        if(Source_NAME == "Liver" && Target_NAME == "Brain"){
                           QTextStream(stdout) << Quantity_NAME << " " << RadioTracer_NAME<< " " << Geometry_NAME  << " " << Particle_NAME  << " " << Target_NAME  << "<--" << Source_NAME << " = " << DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particle_NAME][Source_NAME][Target_NAME]
                                                   << " for Target Organ = " << DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particle_NAME]["All Sources"][Target_NAME] << "\n";
                        }
                        */
                    }
                }
            }
        }
    }

    ui->tableWidgetForOneGraph->clear();
    ui->tableWidgetForOneGraph->setRowCount(0);
    ui->tableWidgetForOneGraph->setColumnCount(0);

    if(DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME].size() == 0){
        QMessageBox::information(this, tr(""), "No data were registered for this configuration, "+ui->comboBoxSourceOrTargetsForBiokinetics->currentText()+", "+ui->comboBoxPhantom->currentText()+", source "+ui->comboBoxSources->currentText()+", or "+ ui->comboBoxQuantityNucl->currentText() +" not calculated for this configuration, or the biokinetics data table of " + ui->comboBoxRadioPharmaceutiques->currentText() + " not defined for this configuration");
        return;
    }

    QStringList headers;

    if(ui->comboBoxSourceOrTargetsForBiokinetics->currentIndex() == 0 ){ // In all sources
        Targets.clear();
        for ( auto DD = RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME].begin(); DD != RadioTracerSourceOrganResidenceTime[ui->comboBoxRadioPharmaceutiques->currentText()][RadioTracer_NAME].end(); ++DD  ){
            Source_NAME = DD.key();
            Targets.push_back(Source_NAME);
            //QTextStream(stdout) << "Source_NAME " << Source_NAME <<"\n";
        }
        headers.append("Source-Target Region");
    }else{
        headers.append("Target Region");
    }

    for(int b=0; b < Particles.size();b++){
        headers.append(Particles[b] +" "+ui->comboBoxQuantityNucl->currentText()+"("+ui->comboBoxEffDoseUnit->currentText()+")");
    }

    ui->tableWidgetForOneGraph->setColumnCount(Particles.size()+1);
    ui->tableWidgetForOneGraph->setShowGrid(true);
    ui->tableWidgetForOneGraph->setSelectionMode(QAbstractItemView::SingleSelection);
    ui->tableWidgetForOneGraph->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->tableWidgetForOneGraph->setHorizontalHeaderLabels(headers);
    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    //ui->tableWidgetForOneGraph->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    ui->tableWidgetForOneGraph->resizeColumnsToContents();
    //ui->tableWidgetForOneGraph->horizontalHeader()->sectionResizeMode(0, QHeaderView::Stretch);

    //QTextStream(stdout) << "Table to show rows " << Data.size() << "\n";

    int row = 0;
    for(int a=0; a < Targets.size();a++){

        if(Quantity_NAME == "E" && TissueFactorMap[Targets[a]] == 0 && Targets[a] != "All Targets"){continue;}

        ui->tableWidgetForOneGraph->insertRow(row);
        ui->tableWidgetForOneGraph->setItem(row,0, new QTableWidgetItem(Targets[a]));

        if(ui->comboBoxSourceOrTargetsForBiokinetics->currentIndex() == 0 ){ // In all sources with each source as a target
            for(int b=0; b < Particles.size();b++){
                //QTextStream(stdout) << "Col "<< b+1 << " Particle " << Particles[b] << " val " << QString::number(DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particles[b]][Targets[a]][Targets[a]]/convfac) <<"\n";
                ui->tableWidgetForOneGraph->setItem(row,b+1, new QTableWidgetItem(QString::number(DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particles[b]][Targets[a]][Targets[a]]/convfac)));
            }
        }
        else if(ui->comboBoxSourceOrTargetsForBiokinetics->currentIndex() == 1){ // From a Source or to each of targets
            for(int b=0; b < Particles.size();b++){
                //QTextStream(stdout) << "Col "<< b+1 << " Particle " << Particles[b] << " val " << QString::number(DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particles[b]][ui->comboBoxSources->currentText()][Targets[a]]/convfac) <<"\n";
                ui->tableWidgetForOneGraph->setItem(row,b+1, new QTableWidgetItem(QString::number(DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particles[b]][ui->comboBoxSources->currentText()][Targets[a]]/convfac)));
            }
        }
        else if(ui->comboBoxSourceOrTargetsForBiokinetics->currentIndex() == 2 || ui->comboBoxSourceOrTargetsForBiokinetics->currentIndex() == 3){ // From all Sources to each of targets (all targets also included in list of targets)
            for(int b=0; b < Particles.size();b++){
                //QTextStream(stdout) << "Col "<< b+1 << " Particle " << Particles[b] << " val " << QString::number(DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particles[b]][ui->comboBoxSources->currentText()][Targets[a]]/convfac) <<"\n";
                ui->tableWidgetForOneGraph->setItem(row,b+1, new QTableWidgetItem(QString::number(DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particles[b]]["All Sources"][Targets[a]]/convfac)));
            }
        }

        row++;
    }

    ui->tableWidgetForOneGraph->resizeColumnsToContents();
    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->viewport()->installEventFilter(this);

}
double MainWindow::GenerateRadiotracerQuantitiesByInterpolationInDefaultUnitForBiokinetic(QString ParticleName, double Energy, QString Source_NAME, QString Target_NAME){

    QString Geometry_NAME = ui->comboBoxPhantom->currentText();
    QString Quantity_Name = ui->comboBoxQuantityNucl->currentText();

    double Val = 0;

    if(ICRPSAFs["AE"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME] != 0){

        if(Quantity_Name == "AE"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]*Energy);
        }
        else if(Quantity_Name == "AF"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]);
        }
        else if(Quantity_Name == "SAF"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME];
        }
        else if(Quantity_Name == "S"){
            Val = (ICRPSAFs["AE"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]/RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]);
        }
        else if(Quantity_Name == "H"){
            //Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(Energy)*GenerateRadiationFactor(ParticleName,Energy);
            Val = (ICRPSAFs["AE"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]/RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME])*GenerateRadiationFactor(ParticleName,Energy);
        }
        else if(Quantity_Name == "E"){
            Val = (ICRPSAFs["AE"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]/RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME])*GenerateRadiationFactor(ParticleName,Energy)*TissueFactorMap[Target_NAME];
        }

        if(Source_NAME == "Liver" && Target_NAME == "Liver"){
            QTextStream(stdout) << " Val " << Val << "\n" ;
        }

    }

    //if(Source_NAME == "Liver" && Target_NAME == "Liver"){
      //  QTextStream(stdout) << " Val befor entering " << ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME] << "\n" ;
    //}


    else if(ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME] != 0){

        if(Quantity_Name == "AE"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]*Energy);
        }
        else if(Quantity_Name == "AF"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]);
        }
        else if(Quantity_Name == "SAF"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME];
        }
        else if(Quantity_Name == "S"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(Energy);
        }
        else if(Quantity_Name == "H"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(Energy)*GenerateRadiationFactor(ParticleName,Energy);
        }
        else if(Quantity_Name == "E"){
            Val = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME]*(Energy)*GenerateRadiationFactor(ParticleName,Energy)*TissueFactorMap[Target_NAME];
        }

        //if(Source_NAME == "Liver" && Target_NAME == "Liver"){
          //  QTextStream(stdout) << " Val " << Val << "\n" ;
        //}

    }
    else {

        double Val1 = 0;
        double Val2 = 0;
        double Energy1;
        double Energy2;

        Energy1 = 0;
        Energy2 = 0;

        std::map<QString,std::map<QString, std::vector<double>>> Sourecee = SourceParticleEnergyValues ;

        for(int ss = 0 ; ss < Sourecee[Source_NAME][ParticleName].size() ; ss++){

            int ff = ss+1;

            int da = Sourecee[Source_NAME][ParticleName].size()-1; if(ff == da){break;}

            double E1 = Sourecee[Source_NAME][ParticleName][ss];
            double E2 = Sourecee[Source_NAME][ParticleName][ff];
            //QTextStream(stdout) << ss << " " << E1 << " " << ff << " " << E2 << "\n" ;

            if(E1 < Energy && Energy < E2){
                Energy1 = E1;
                Energy2 = E2;
                if(Source_NAME == "Liver" && Target_NAME == "Liver"){
                    //QTextStream(stdout) << " ParticleName " << ParticleName << " Energy1 " << Energy1 << " Energy " << Energy << " Energy2 " << Energy2 << "\n" ;
                }
                break;
            }
        }

        if(Quantity_Name == "AE"){
            Val1 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy1][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]*Energy1);
            Val2 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy2][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]*Energy2);
        }
        else if(Quantity_Name == "AF"){
            Val1 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy1][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]);
            Val2 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy2][Source_NAME][Target_NAME]*(RegionParameterValueMap[Geometry_NAME]["Mass"][Target_NAME]);
        }
        else if(Quantity_Name == "SAF"){
            Val1 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy1][Source_NAME][Target_NAME];
            Val2 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy2][Source_NAME][Target_NAME];
        }
        else if(Quantity_Name == "S"){
            Val1 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy1][Source_NAME][Target_NAME]*(Energy1);
            Val2 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy2][Source_NAME][Target_NAME]*(Energy2);
        }
        else if(Quantity_Name == "H"){
            Val1 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy1][Source_NAME][Target_NAME]*(Energy1)*GenerateRadiationFactor(ParticleName,Energy1);
            Val2 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy2][Source_NAME][Target_NAME]*(Energy2)*GenerateRadiationFactor(ParticleName,Energy2);
        }
        else if(Quantity_Name == "E"){
            Val1 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy1][Source_NAME][Target_NAME]*(Energy1)*GenerateRadiationFactor(ParticleName,Energy1)*TissueFactorMap[Target_NAME];
            Val2 = ICRPSAFs["SAF"][Geometry_NAME][ParticleName][Energy2][Source_NAME][Target_NAME]*(Energy2)*GenerateRadiationFactor(ParticleName,Energy2)*TissueFactorMap[Target_NAME];
        }

        if(ui->checkBoxInterpolationType->isChecked()){
            Val  = std::exp(std::log(Val1) + (std::log(Energy) - std::log(Energy1)) * (std::log(Val2) - std::log(Val1)) / (std::log(Energy2) - std::log(Energy1)));
        }else{
            Val  = Val1 + (Energy-Energy1)*((Val2-Val1)/(Energy2-Energy1));
        }

        //if(Source_NAME == "Liver" && Target_NAME == "Liver"){
          //  QTextStream(stdout) << " Val1 " << Val1  << " Val " << Val  << " Val2 " << Val2 << "\n" ;
        //}

    }

    if( !__isnan(Val) && !__isinf(Val) && Val != 0 && Val != NULL){
        return Val;
        //ICRPSAFs[Quantity_Name][Geometry_NAME][ParticleName][Energy][Source_NAME][Target_NAME] = Val;
    }else{
        return 0.;
    }
}


void MainWindow::Read_ICRP107_108Files(QString DataDirName ){

    // this for ICRP radionuclides, this numbers should be the indicators in SAFs files
    EnergyIDRadionuclideForNeutronSAF["U-238"  ]  =111;
    EnergyIDRadionuclideForNeutronSAF["Pu-236" ] =222;
    EnergyIDRadionuclideForNeutronSAF["Pu-238" ] =333;
    EnergyIDRadionuclideForNeutronSAF["Pu-240" ] =444;
    EnergyIDRadionuclideForNeutronSAF["Pu-242" ] =555;
    EnergyIDRadionuclideForNeutronSAF["Pu-244" ] =666;
    EnergyIDRadionuclideForNeutronSAF["Cm-240" ] =777;
    EnergyIDRadionuclideForNeutronSAF["Cm-242" ] =888;
    EnergyIDRadionuclideForNeutronSAF["Cm-244" ] =999;
    EnergyIDRadionuclideForNeutronSAF["Cm-245" ] =101010;
    EnergyIDRadionuclideForNeutronSAF["Cm-246" ] =111111;
    EnergyIDRadionuclideForNeutronSAF["Cm-248" ] =121212;
    EnergyIDRadionuclideForNeutronSAF["Cm-250" ] =131313;
    EnergyIDRadionuclideForNeutronSAF["Cf-246" ] =141414;
    EnergyIDRadionuclideForNeutronSAF["Cf-248" ] =151515;
    EnergyIDRadionuclideForNeutronSAF["Cf-249" ] =161616;
    EnergyIDRadionuclideForNeutronSAF["Cf-250" ] =171717;
    EnergyIDRadionuclideForNeutronSAF["Cf-252" ] =181818;
    EnergyIDRadionuclideForNeutronSAF["Cf-254" ] =191919;
    EnergyIDRadionuclideForNeutronSAF["Es-253" ] =202020;
    EnergyIDRadionuclideForNeutronSAF["Es-254" ] =212121;
    EnergyIDRadionuclideForNeutronSAF["Es-254m"] =222222;
    EnergyIDRadionuclideForNeutronSAF["Es-255" ] =232323;
    EnergyIDRadionuclideForNeutronSAF["Fm-252" ] =242424;
    EnergyIDRadionuclideForNeutronSAF["Fm-254" ] =252525;
    EnergyIDRadionuclideForNeutronSAF["Fm-255" ] =262626;
    EnergyIDRadionuclideForNeutronSAF["Fm-256" ] =272727;
    EnergyIDRadionuclideForNeutronSAF["Fm-257" ] =282828;

    RadioParticleEnergyYieldForSpectrum.clear();
    ICRPRadioNuclideData.clear();
    ICRPRadioNuclideHalfLives.clear();
    RadionuclidesParticles.clear();
    RadioTracerSourceOrganResidenceTime.clear();
    RadiotracerradionucleidMap.clear();

    ui->progressBarReadingCalcData->setRange(0, 100);
    ui->progressBarReadingCalcData->setValue(0);
    ui->progressBarReadingCalcData->show();
    double percent = 0;

    for (int zzz = 0 ; zzz < 2; zzz++) { // ICRP files  "< 14"

        percent = ((zzz+1)/2.)*100;

        ui->progressBarReadingCalcData->setValue(percent);

        QVector< double > particleEnergies;

        QString Geometry ,
                Quantity = "SAF",
                SrcRegionName ,
                organTargetname,
                ParticleName,
                RadioNuclideName,
                FilePathString,
                word;

        if(zzz==0){
            FilePathString = "ICRP-07.RAD";
        }
        else if(zzz==1){
            FilePathString = "ICRP128RadioPharmaceuticalsData";
        }

        QString fm = DataDirName+"/"+FilePathString;
        QFile file(fm);

        if(!file.open(QIODevice::ReadOnly)) {
            //QMessageBox::information(0, "error", file.errorString());
        }else{

            //QTextStream(stdout) <<zzz << " -Reading file: " << fm << "....\n";

            QTextStream in(&file);

            int lineInc = 0;
            int NumOfDataLines = 0;
            int DataLineInc = 0;
            int RadioNuclideDataInc = 0;
            int NumberOfEneRows = 0;
            double LastEnergy = 0;
            double Prob = 0;
            while(!in.atEnd()) {

                QRegExp space("\\s++");
                QString line = in.readLine().remove(space);

                //if(lineInc == 10000){
                //  break;
                //}

                if(!line.isEmpty()){

                    //QTextStream(stdout) << "line: " << line << "\n";

                    if(zzz==0){ // for Radiation files

                        QStringList fields = line.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);

                        if(fields.size() == 3){ // radionu data

                            RadioNuclideName = fields[0];
/*
                            if(
                                    //positron emitters for PET
                                    RadioNuclideName != "F-18" &&
                                    RadioNuclideName != "O-15" &&
                                    RadioNuclideName != "C-11" &&
                                    RadioNuclideName != "N-13" &&

                                    //Targeted therapy with alpha emitters
                                    RadioNuclideName != "At-211" &&
                                    RadioNuclideName != "Bi-213" &&
                                    RadioNuclideName != "Ac-225" &&
                                    RadioNuclideName != "Ra-223" &&
                                    RadioNuclideName != "Th-227" &&
                                    //Spontaneous neutron emitters
                                    RadioNuclideName != "U-238" &&
                                    RadioNuclideName != "Pu-240" &&
                                    RadioNuclideName != "Fm-258" &&
                                    RadioNuclideName != "Es-254" &&
                                    RadioNuclideName != "Cm-244" &&
                                    RadioNuclideName != "Cf-252" &&
                                    // For cancer imaging
                                    RadioNuclideName != "I-131" &&
                                    RadioNuclideName != "Ga-68" &&
                                    RadioNuclideName != "Tc-99m"

                                    ){
                                continue;
                            }
*/
                            int HLConversion = 1;
                            QString HLunit = "s";
                            if(fields[1].contains("s")){HLConversion = 1;HLunit = "s";}
                            else if(fields[1].contains("m")){HLConversion = 60;HLunit = "m";}
                            else if(fields[1].contains("h")){HLConversion = 60*60;HLunit = "h";}
                            else if(fields[1].contains("d")){HLConversion = 60*60*24;HLunit = "d";}
                            else if(fields[1].contains("y")){HLConversion = 60*60*24*365;HLunit = "y";}

                            ICRPRadioNuclideHalfLives[RadioNuclideName] = HLConversion*fields[1].remove(HLunit).toDouble() ;

                            //QTextStream(stdout) << "RadioNuclideName " << RadioNuclideName << " Half-life "<< ICRPRadioNuclideHalfLives[RadioNuclideName] << "\n";

                        }else if(fields.size() == 4) { //radiation data ene yield par...

/*
                            if(
                                    //positron emitters for PET
                                    RadioNuclideName != "F-18" &&
                                    RadioNuclideName != "O-15" &&
                                    RadioNuclideName != "C-11" &&
                                    RadioNuclideName != "N-13" &&

                                    //Targeted therapy with alpha emitters
                                    RadioNuclideName != "At-211" &&
                                    RadioNuclideName != "Bi-213" &&
                                    RadioNuclideName != "Ac-225" &&
                                    RadioNuclideName != "Ra-223" &&
                                    RadioNuclideName != "Th-227" &&
                                    //Spontaneous neutron emitters
                                    RadioNuclideName != "U-238" &&
                                    RadioNuclideName != "Pu-240" &&
                                    RadioNuclideName != "Fm-258" &&
                                    RadioNuclideName != "Es-254" &&
                                    RadioNuclideName != "Cm-244" &&
                                    RadioNuclideName != "Cf-252" &&
                                    // For cancer imaging
                                    RadioNuclideName != "I-131" &&
                                    RadioNuclideName != "Ga-68" &&
                                    RadioNuclideName != "Tc-99m"

                                    ){
                                continue;
                            }
*/
                            QStringList fields = line.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);

                            if(fields[3] == "G"){ParticleName = "gamma";}
                            else if(fields[3] == "PG"){ParticleName = "gamma";}
                            else if(fields[3] == "DG"){ParticleName = "gamma";}
                            else if(fields[3] == "X"){ParticleName = "gamma";}
                            else if(fields[3] == "AQ"){ParticleName = "gamma";}
                            else if(fields[3] == "B+"){ParticleName = "e+";}
                            else if(fields[3] == "B-"){ParticleName = "e-";}
                            else if(fields[3] == "BD"){ParticleName = "e-";}
                            else if(fields[3] == "IE"){ParticleName = "e-";}
                            else if(fields[3] == "AE"){ParticleName = "e-";}
                            else if(fields[3] == "A"){ParticleName = "alpha";}
                            else if(fields[3] == "AR"){ParticleName = "alpha";}
                            else if(fields[3] == "FF"){ParticleName = "FF";}
                            else if(fields[3] == "N"){ParticleName = "neutron";}

                            //QTextStream(stdout) << "RadioNuclideName " << RadioNuclideName << " fields[3] " << fields[3]  << " ParticleName " << ParticleName << " fields[2].toDouble() " << fields[2].toDouble() << " fields[1].toDouble() " << fields[1].toDouble() << "\n";

                            // because beta spectrum is readed from beta files
                            if(ui->checkBoxReadSpectrum->isChecked()){
                                if(fields[3] == "B+" || fields[3] == "B-"){
                                }else{
                                  ICRPRadioNuclideData[RadioNuclideName][ParticleName][fields[2].toDouble()] = fields[1].toDouble();
                                }
                            }else{
                                ICRPRadioNuclideData[RadioNuclideName][ParticleName][fields[2].toDouble()] = fields[1].toDouble();
                            }

                        }
                    }
                    else if(zzz==1){ // for Radiotracer data

                        std::string c1; std::string c2; std::string c3; double val;
                        std::string line1 = line.toStdString();
                        std::istringstream LineString(line1);

                        if(LineString.str().empty()){
                            continue;
                        }

                        LineString  >> c1 ;
                        LineString  >> c2 ;

                        RadiotracerradionucleidMap[c1.c_str()] = c2.c_str();

                        while(LineString >> c3 ){
                           LineString >> val;
                           RadioTracerSourceOrganResidenceTime[c1.c_str()][c2.c_str()][c3.c_str()] = val;
                        }
                    }
                }
                lineInc++;
            }
            file.close();
            //QTextStream(stdout) << "Closing file " << file.fileName() << "\n";
        }
    }

    ui->progressBarReadingCalcData->setValue(100);

    if(ICRPRadioNuclideHalfLives.size() != 0){
        ui->pushButtonReadICRP107128->setText("Read ICRP Radiation Data *");
    }else{
        ui->pushButtonReadICRP107128->setText("Read ICRP Radiation Data");
    }

}
void MainWindow::Read_ICRP107SpectrumRadiationFiles(QString DataDirName ){

    ICRPRadioNuclideDataDiscSpec.clear();

    ui->progressBarReadingCalcData->setRange(0, 100);
    ui->progressBarReadingCalcData->setValue(0);
    ui->progressBarReadingCalcData->show();
    double percent = 0;

    for (int zzz = 0 ; zzz < 3; zzz++) { // ICRP files  "< 14"

        percent = ((zzz+1)/3)*100;

        ui->progressBarReadingCalcData->setValue(percent);

        QVector< double > particleEnergies;

        QString Geometry ,
                Quantity = "SAF",
                SrcRegionName ,
                organTargetname,
                ParticleName,
                RadioNuclideName,
                FilePathString,
                word;

        if(zzz==0){
            FilePathString = "ICRP-07.BET";
        }
        else if(zzz==1){
            FilePathString = "ICRP-07.ACK";
        }
        else if(zzz==2){
            FilePathString = "ICRP-07.NSF";
        }

        QString fm = DataDirName+"/"+FilePathString;
        QFile file(fm);

        if(!file.open(QIODevice::ReadOnly)) {
            //QMessageBox::information(0, "error", file.errorString());

        }else{

            //QTextStream(stdout) <<zzz << " -Reading file: " << fm << "....\n";

            QTextStream in(&file);

            int lineInc = 0;
            int NumOfDataLines = 0;
            int DataLineInc = 0;
            int RadioNuclideDataInc = 0;
            int NumberOfEneRows = 0;
            double LastEnergy = 0;
            double Prob = 0;
            while(!in.atEnd()) {

                QRegExp space("\\s++");
                QString line = in.readLine().remove(space);

                //if(lineInc == 10000){
                //  break;
                //}

                if(!line.isEmpty()){

                    //QTextStream(stdout) << "line: " << line << "\n";

                    if(zzz==0){ // for beta (e- and e+) Radiation files
                        // in this filee the e- and e+ spectrums are combined in one spectrum, if a e+ is defined in radiation file, then the AQ is defined too

                        QStringList fields = line.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);

                        if(fields.size() == 2){ // radionu data

                            if(NumberOfEneRows == 0){ //
                                RadioNuclideName = fields[0];
                                NumberOfEneRows = fields[1].toInt();

                            }else if(RadioNuclideDataInc > NumberOfEneRows){

                                RadioNuclideName = fields[0];
                                NumberOfEneRows = fields[1].toInt();
                                RadioNuclideDataInc = 0;

                                //QTextStream(stdout) << "RadioNuclideName " << RadioNuclideName << " NumberOfEneRows " << NumberOfEneRows  << " RadioNuclideDataInc " << RadioNuclideDataInc  << "\n";
                            }
                            else{

                                if(fields[0].toDouble() == 0 && fields[1].toDouble() == 0){

                                }else{

                                    //QTextStream(stdout) << RadioNuclideDataInc << " RadioNuclideName " << RadioNuclideName << " ParticleName " << ParticleName << " RadioNuclideDataInc " << RadioNuclideDataInc << " fields[0] " << fields[0].toDouble()  << " fields[1] " << fields[1].toDouble()  << "\n";

                                    //double Energy = LastEnergy + (fields[0].toDouble()-LastEnergy)*QRandomGenerator::global()->generateDouble();
                                    double Energy = fields[0].toDouble();
                                    Prob = fields[1].toDouble();

                                    //QTextStream(stdout) << RadioNuclideDataInc << " RadioNuclideName " << RadioNuclideName << " RadioNuclideDataInc " << RadioNuclideDataInc << " Random Energy " << Energy << " Prob " << Prob << " fields[0] " << fields[0].toDouble()  << " fields[1] " << fields[1].toDouble()  << "\n";
/*
                                    if(
                                    //positron emitters for PET
                                    RadioNuclideName != "F-18" &&
                                    RadioNuclideName != "O-15" &&
                                    RadioNuclideName != "C-11" &&
                                    RadioNuclideName != "N-13" &&

                                    //Targeted therapy with alpha emitters
                                    RadioNuclideName != "At-211" &&
                                    RadioNuclideName != "Bi-213" &&
                                    RadioNuclideName != "Ac-225" &&
                                    RadioNuclideName != "Ra-223" &&
                                    RadioNuclideName != "Th-227" &&
                                    //Spontaneous neutron emitters
                                    RadioNuclideName != "U-238" &&
                                    RadioNuclideName != "Pu-240" &&
                                    RadioNuclideName != "Fm-258" &&
                                    RadioNuclideName != "Es-254" &&
                                    RadioNuclideName != "Cm-244" &&
                                    RadioNuclideName != "Cf-252" &&
                                    // For cancer imaging
                                    RadioNuclideName != "I-131" &&
                                    RadioNuclideName != "Ga-68" &&
                                    RadioNuclideName != "Tc-99m"
                                            ){
                                        continue;
                                    }
*/

                                    //if(RadioNuclideName == "F-18"){
                                    //    QTextStream(stdout) << " RadioNuclideName " << RadioNuclideName << " Energy " << Energy << " Prob " << Prob << " Energy*Prob " << Energy*Prob  << "\n";
                                    //}
                                    ICRPRadioNuclideData        [RadioNuclideName]["e-"][Energy] = (Energy-LastEnergy)*Prob;
                                    ICRPRadioNuclideDataDiscSpec[RadioNuclideName]["e-"]["Spectrum"][Energy] = (Energy-LastEnergy)*Prob;
                                    //ICRPRadioNuclideDataDiscSpec[RadioNuclideName]["e-"]["Total"][1] += Energy*Prob;

                                    //ICRPRadioNuclideDataDiscSpec[RadioNuclideName]["e+"]["Spectrum"][Energy] = Prob;
                                    //ICRPRadioNuclideDataDiscSpec[RadioNuclideName]["e+"]["Total"][1] += Prob;

                                    LastEnergy = fields[0].toDouble();
                                }
                            }
                            RadioNuclideDataInc++;
                        }
                    }
                    else if(zzz==1){ // for Auger elect...

                        QStringList fields = line.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);

                        if(NumberOfEneRows == 0){ //
                            if(fields.size() == 3){ // radionu data
                                RadioNuclideName = fields[0];
                                NumberOfEneRows = fields[2].toInt();
                            }

                        }else if(NumberOfEneRows == RadioNuclideDataInc){

                            if(fields.size() == 3){
                                RadioNuclideName = fields[0];
                                NumberOfEneRows = fields[2].toInt();
                                RadioNuclideDataInc = 0;
                            }
                            //QTextStream(stdout) << "RadioNuclideName " << RadioNuclideName << " NumberOfEneRows " << NumberOfEneRows  << " RadioNuclideDataInc " << RadioNuclideDataInc  << "\n";

                        }else{
/*
                            if(
                                    //positron emitters for PET
                                    RadioNuclideName != "F-18" &&
                                    RadioNuclideName != "O-15" &&
                                    RadioNuclideName != "C-11" &&
                                    RadioNuclideName != "N-13" &&

                                    //Targeted therapy with alpha emitters
                                    RadioNuclideName != "At-211" &&
                                    RadioNuclideName != "Bi-213" &&
                                    RadioNuclideName != "Ac-225" &&
                                    RadioNuclideName != "Ra-223" &&
                                    RadioNuclideName != "Th-227" &&
                                    //Spontaneous neutron emitters
                                    RadioNuclideName != "U-238" &&
                                    RadioNuclideName != "Pu-240" &&
                                    RadioNuclideName != "Fm-258" &&
                                    RadioNuclideName != "Es-254" &&
                                    RadioNuclideName != "Cm-244" &&
                                    RadioNuclideName != "Cf-252" &&
                                    // For cancer imaging
                                    RadioNuclideName != "I-131" &&
                                    RadioNuclideName != "Ga-68" &&
                                    RadioNuclideName != "Tc-99m"
                                    ){
                                continue;
                            }
*/
                            //QTextStream(stdout) << RadioNuclideDataInc << " RadioNuclideName " << RadioNuclideName << " RadioNuclideDataInc " << RadioNuclideDataInc << " fields[0] " << fields[0].toDouble()  << " fields[1] " << fields[1].toDouble()  << "\n";
                            //ICRPRadioNuclideDataDiscSpec[RadioNuclideName][ParticleName]["Discrete"][fields[0].toDouble()] = fields[1].toDouble();

                            RadioNuclideDataInc ++;
                        }
                    }
                    else if(zzz==2 && ReadNeutronSpectrum){ // for Neutron spectrum fission file

                        QStringList fields = line.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
                            //QTextStream(stdout) << "size "<< fields.size() << " -line: " << line << "\n";

                        if((lineInc == 0) || (lineInc != 0 && NumOfDataLines == DataLineInc)){ // radionu data

                            RadioNuclideName = fields[0];
                            NumOfDataLines = fields[2].toInt();
                            DataLineInc = 0;
                            //QTextStream(stdout) << "\n\n\n#"<< RadioNuclideName << "\n/SourceData/setEventsInitialEneData MeV Spectrum 0.2222222 0.000 0.00E+00 ";

                            //QTextStream(stdout) << "RadioNuclideName " << RadioNuclideName << " fields[1].toDouble() " << fields[1].toDouble() << " fields[2].toDouble() " << fields[2].toDouble()  << " DataLineInc " << DataLineInc<< " NumOfDataLines " << NumOfDataLines << "\n";

                        }else{

                            //QTextStream(stdout) << "RadioNuclideName " << RadioNuclideName << " fields[0].toDouble() " << fields[0].toDouble() << " fields[1].toDouble() " << fields[1].toDouble() << " fields[2].toDouble() " << fields[2].toDouble()<< "\n";
                            //QTextStream(stdout) << fields[0].toDouble() << " " << fields[2].toDouble() << " " ;
                            if(NumOfDataLines-1 == DataLineInc){ // radionu data
                                //QTextStream(stdout) << fields[1].toDouble() << " 0 ";
                            }
/*
                            if(
                                    //positron emitters for PET
                                    RadioNuclideName != "F-18" &&
                                    RadioNuclideName != "O-15" &&
                                    RadioNuclideName != "C-11" &&
                                    RadioNuclideName != "N-13" &&

                                    //Targeted therapy with alpha emitters
                                    RadioNuclideName != "At-211" &&
                                    RadioNuclideName != "Bi-213" &&
                                    RadioNuclideName != "Ac-225" &&
                                    RadioNuclideName != "Ra-223" &&
                                    RadioNuclideName != "Th-227" &&
                                    //Spontaneous neutron emitters
                                    RadioNuclideName != "U-238" &&
                                    RadioNuclideName != "Pu-240" &&
                                    RadioNuclideName != "Fm-258" &&
                                    RadioNuclideName != "Es-254" &&
                                    RadioNuclideName != "Cm-244" &&
                                    RadioNuclideName != "Cf-252" &&
                                    // For cancer imaging
                                    RadioNuclideName != "I-131" &&
                                    RadioNuclideName != "Ga-68" &&
                                    RadioNuclideName != "Tc-99m"
                                    ){
                                continue;
                            }
*/
                            Prob = fields[2].toDouble();
                            //ICRPRadioNuclideFSNData[RadioNuclideName][fields[0].toDouble()][fields[1].toDouble()] = fields[2].toDouble();
                            //ICRPRadioNuclideData[RadioNuclideName]["neutron"][fields[0].toDouble()] = fields[2].toDouble();
                            ICRPRadioNuclideDataDiscSpec[RadioNuclideName]["neutron"]["Spectrum"][fields[0].toDouble()] = fields[2].toDouble();
                            ICRPRadioNuclideData        [RadioNuclideName]["neutron"][fields[1].toDouble()] = Prob;

                            DataLineInc++;

                            if(NumOfDataLines == DataLineInc){
                                ICRPRadioNuclideDataDiscSpec[RadioNuclideName]["neutron"]["Spectrum"][fields[1].toDouble()] = 0;
                            }
                        }
                    }
                }
                lineInc++;
            }
            file.close();
            //QTextStream(stdout) << "Closing file " << file.fileName() << "\n";
        }
    }

    ReadNeutronSpectrum = true;
    ui->progressBarReadingCalcData->setValue(100);

    /*

    for ( auto it = ICRPRadioNuclideDataDiscSpec.begin(); it != ICRPRadioNuclideDataDiscSpec.end(); ++it  ){

        QString RadioTracer_NAME = it.key();
        //        QTextStream(stdout) << "--Geometry_NAME " << Geometry_NAME <<"\n";

        for ( auto it2 = it.value().begin(); it2 != it.value().end(); ++it2  ){
            QString Particle_NAME = it2.key();
            //            QTextStream(stdout) << "---Particle_NAME " << Particle_NAME <<"\n";

            for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
                QString Source_NAME  = it3.key();
                //                QTextStream(stdout) << "----Energy_Val " << Energy_Val <<"\n";

                for ( auto DD = it3.value().begin(); DD != it3.value().end(); ++DD  ){
                    double Energy  = DD.key();
                    //                  QTextStream(stdout) << "-----Source_NAME " << Source_NAME <<"\n";

                    double SAFValue  = DD.value();
                    //QTextStream(stdout) << " RadioTracer_NAME " << RadioTracer_NAME << " Period(s) " << ICRPRadioNuclideHalfLives[RadioTracer_NAME] << " Particle_NAME " << Particle_NAME << " Energy " << Energy << " Value " << SAFValue <<"\n";
                }
            }
        }
    }
    */
}
void MainWindow::CombineSpectrumWithMonoData(){

    QMap<QString,QMap<QString,QMap<double,double>>> ICRPRadioNuclideData1 = ICRPRadioNuclideData;

    for ( auto it = ICRPRadioNuclideData.begin(); it != ICRPRadioNuclideData.end(); ++it  ){
        QString RadioTracer_NAME = it.key();

        for ( auto it2 = it.value().begin(); it2 != it.value().end(); ++it2  ){
            QString Particle_NAME = it2.key();

            for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
                double Energy_Val = it3.key();


                if(Particle_NAME == "B-" || Particle_NAME == "B+"){

                    QTextStream(stdout) << " RadioTracer_NAME "<< RadioTracer_NAME << " Particle_NAME "<< Particle_NAME  << " Energy_Val "<< Energy_Val << " yield "<< it3.value() <<  "\n";

                    /*
                    if(RadioTracer_NAME == "F-18"){
                        QTextStream(stdout) << " RadioTracer_NAME "<< RadioTracer_NAME << " Particle_NAME "<< Particle_NAME  << " Energy_Val "<< Energy_Val << " it3.value() "<< it3.value() <<  "\n";
                    }
                    */
                    for ( auto AA = ICRPRadioNuclideDataDiscSpec[RadioTracer_NAME]["e-"]["Spectrum"].begin(); AA != ICRPRadioNuclideDataDiscSpec[RadioTracer_NAME]["e-"]["Spectrum"].end(); ++AA  ){
                        double energy = AA.key();
                        ICRPRadioNuclideData1[RadioTracer_NAME]["e-"][energy] = it3.value() * (AA.key()/ICRPRadioNuclideDataDiscSpec[RadioTracer_NAME]["e-"]["Total"][1]);
                    }
                }
            }
        }
    }

    for ( auto it = ICRPRadioNuclideData.begin(); it != ICRPRadioNuclideData.end(); ++it  ){
        QString RadioTracer_NAME = it.key();

        for ( auto it2 = it.value().begin(); it2 != it.value().end(); ++it2  ){
            QString Particle_NAME = it2.key();

            for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
                double Energy_Val = it3.key();

                if(Particle_NAME == "B-" || Particle_NAME == "B+"){

                    ICRPRadioNuclideData1[RadioTracer_NAME].remove("B+");
                    ICRPRadioNuclideData1[RadioTracer_NAME].remove("B-");
                }

                if(Particle_NAME == "G"){ICRPRadioNuclideData1[RadioTracer_NAME]["gamma"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "PG"){ICRPRadioNuclideData1[RadioTracer_NAME]["gamma"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "DG"){ICRPRadioNuclideData1[RadioTracer_NAME]["gamma"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "X"){ICRPRadioNuclideData1[RadioTracer_NAME]["gamma"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "AQ"){ICRPRadioNuclideData1[RadioTracer_NAME]["gamma"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "B+"){ICRPRadioNuclideData1[RadioTracer_NAME]["e+"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "B-"){ICRPRadioNuclideData1[RadioTracer_NAME]["e-"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "BD"){ICRPRadioNuclideData1[RadioTracer_NAME]["e-"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "IE"){ICRPRadioNuclideData1[RadioTracer_NAME]["e-"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "AE"){ICRPRadioNuclideData1[RadioTracer_NAME]["e-"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "A"){ICRPRadioNuclideData1[RadioTracer_NAME]["alpha"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "AR"){ICRPRadioNuclideData1[RadioTracer_NAME]["alpha"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "FF"){ICRPRadioNuclideData1[RadioTracer_NAME]["FF"][Energy_Val]=it3.value();}
                else if(Particle_NAME == "N"){ICRPRadioNuclideData1[RadioTracer_NAME]["neutron"][Energy_Val]=it3.value();}
            }
        }
    }

    ICRPRadioNuclideData.clear();
    ICRPRadioNuclideData = ICRPRadioNuclideData1;

}
void MainWindow::Read_DoseCalcs_file(QString FilePath){

    QTextStream(stdout) << "---------------- Read_DoseCalcs_file() ------------------" << "\n";
    showResultsOutput("Reading Result file... " + FilePath, 4);

    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> ResultTable;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> ErrorTables;

    QString line, Geometry, Quantity,SrcRegionName , organTargetname, RadTracerName , ParticleName , word;
    std::string Quantity1, EnergyDist, Geometry1, SrcRegionName1 , organTargetname1,RadTracerName1 ,  ParticleName1 , word1;
    double mass, volume, density, particleE, Val, StandardDeviation, err, CompuTime;
    unsigned long long int ival;
    bool isADataLines = false ;

    std::ifstream fileR;

    fileR.open(FilePath.toStdString(), std::ios::binary | std::ios::in);
    fileR.seekg(0, std::ios::end);
    long file_size = fileR.tellg();
    fileR.close();

    fileR.open(FilePath.toStdString(), std::ios::binary | std::ios::in);

    long curr, end;
    int percent;
    end = file_size;
    ui->progressBarReadingCalcData->setRange(0, 100);
    ui->progressBarReadingCalcData->setValue(0);
    ui->progressBarReadingCalcData->show();

    if(fileR.is_open()){

        TissueFactorMap.clear();
        PhantomSourcesMap.clear();
        RegionParameterValueMap.clear();
        ICRPSourceMassMap.clear();
        ICRPSAFs.clear();
        SourceParticleEnergyValues.clear();
        ui->comboBoxPhantom->clear();
        ui->comboBoxTargets->clear();
        ui->comboBoxSources->clear();

        std::string linee = line.toStdString() ;

        // to fill the map array DataTables[SrcRegionName][organTargetname][ParticleName][particleE] = Val;
        while ( getline(fileR, linee)) {

            //QTextStream(stdout) << " the line " << linee.c_str() << "\n" ;

            curr = fileR.tellg();
            if (curr != -1)
                percent = curr * 100 / end;
            else
                percent = 100;

            ui->progressBarReadingCalcData->setValue(percent);

            std::istringstream LineString(linee);
            LineString >> word1;
            word = QString::fromStdString(word1);

            //G4cout << " the word " << word << G4"\n" ;

            if(isADataLines == true){

                if (word == "*") {

                    isADataLines = false;
                    continue;
                }
                if (word == "#") {
                    continue;
                }
                organTargetname = word;

                LineString >> Val;
                LineString >> StandardDeviation;
                LineString >> err;
                LineString >> ival; // for Num Steps
                LineString >> mass;
                LineString >> volume;
                LineString >> density;

                if(ParticleName == "RadioTracer"){

                }
                else{

                    if(Quantity == "AE"){
                        ICRPSAFs[Quantity][Geometry][ParticleName][particleE][SrcRegionName][organTargetname] = Val/200000000;
                    }else{
                        ICRPSAFs[Quantity][Geometry][ParticleName][particleE][SrcRegionName][organTargetname] = Val;
                    }
                    //QTextStream(stdout) << " Quantity " << Quantity << " Geometry " << Geometry << " SrcRegionName " << SrcRegionName<< " ParticleName " << ParticleName << " particleE " << particleE << " Val " << Val << "\n" ;

                    bool isin = false;
                    for (int dd = 0 ; dd < SourceParticleEnergyValues[SrcRegionName][ParticleName].size(); dd++) {
                        if(particleE == SourceParticleEnergyValues[SrcRegionName][ParticleName][dd]){isin = true;break;}}
                    if(isin == false){SourceParticleEnergyValues[SrcRegionName][ParticleName].push_back(particleE);                                    }
                }

                RegionParameterValueMap[Geometry]["Mass"][organTargetname] = mass;
                RegionParameterValueMap[Geometry]["Volume"][organTargetname] = volume;
                RegionParameterValueMap[Geometry]["Density"][organTargetname] = density;
                ICRPSourceMassMap[Geometry]["Mass"][organTargetname] = mass;
            }

            if (word == "******" && isADataLines == false) {

                LineString >> Quantity1 >> SrcRegionName1 >> ParticleName1 ;

                //std::cout << " Quantity " << Quantity << " Geometry " << Geometry << " SrcRegionName " << SrcRegionName<< " ParticleName " << ParticleName  << std::"\n" ;

                Quantity = QString::fromStdString(Quantity1);
                ParticleName = QString::fromStdString(ParticleName1);
                SrcRegionName = QString::fromStdString(SrcRegionName1);

                if(ParticleName == "RadioTracer"){
                    LineString >> RadTracerName1 >> Geometry1 ;
                    RadTracerName = QString::fromStdString(RadTracerName1);
                    Geometry = QString::fromStdString(Geometry1);
                }
                else{
                    LineString >> particleE
                            >> Geometry1 >> word1 >> word1 >> word1 >> word1
                            >> word1 >> word1 >> word1 >> word1 >> word1
                            >> word1 >> word1 >> word1 >> word1 >> word1
                            >> word1 >> CompuTime ;

                    Geometry = QString::fromStdString(Geometry1);
                }

                //LineString >> EnergyDist >> EnergyDist;

                bool isin = false;
                for (int dd = 0 ; dd < PhantomSourcesMap[Geometry].size(); dd++) {
                    if(SrcRegionName == PhantomSourcesMap[Geometry][dd]){isin = true;break;}}
                if(isin == false){PhantomSourcesMap[Geometry].push_back(SrcRegionName);}

                isADataLines = true;
            }
        }

        fileR.close();
    }
    else{
        showResultsOutput("canno't open the file : "+FilePath , 4 );
    }

    ui->progressBarReadingCalcData->setValue(100);

}
void MainWindow::Read_Reference_file(QString FilePath){

    QTextStream(stdout) << "---------------- Read_Comparison_file() ------------------" << "\n";
    showResultsOutput("Reading Reference file... " + FilePath, 4);

    // organSource, organTarget, particle , energy, SAF
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> ReferenceTable ;

    QString line ,Geometry ,Quantity, SrcRegionName , organTargetname, ParticleName, RadTracerName, word;
    double Val;

    bool isADataLines = false ;
    int NumOfRefEne;

    std::string  ParticleName1 ,Geometry1, RadTracerName1, SrcRegionName1, organTargetname1 ,Quantity1, word1, QuantityUnit, EnergyUnit;

    QVector< double > particleEnergies;

    std::ifstream fileR;

    fileR.open(FilePath.toStdString(), std::ios::binary | std::ios::in);
    fileR.seekg(0, std::ios::end);
    long file_size = fileR.tellg();
    fileR.close();

    fileR.open(FilePath.toStdString(), std::ios::binary | std::ios::in);


    long curr, end;
    int percent;
    end = file_size;
    ui->progressBarReadingCalcData->setRange(0, 100);
    ui->progressBarReadingCalcData->setValue(0);
    ui->progressBarReadingCalcData->show();


    TissueFactorMap.clear();
    PhantomSourcesMap.clear();
    RegionParameterValueMap.clear();
    ICRPSourceMassMap.clear();
    ICRPSAFs.clear();
    SourceParticleEnergyValues.clear();
    ui->comboBoxPhantom->clear();
    ui->comboBoxTargets->clear();
    ui->comboBoxSources->clear();

    if(fileR.is_open()){

        std::string linee = line.toStdString() ;

        // to fill the map array DataTables[SrcRegionName][organTargetname][ParticleName][particleE] = Val;
        while (getline(fileR, linee)) {

            curr = fileR.tellg();
            if (curr != -1)
                percent = curr * 100 / end;
            else
                percent = 100;

            ui->progressBarReadingCalcData->setValue(percent);


            std::istringstream LineString(linee);
            LineString >> word1;
            word = QString::fromStdString(word1);

            if (word == "******" && isADataLines == false) {

                particleEnergies.clear();

                LineString >> Quantity1 >> QuantityUnit >> SrcRegionName1 >> ParticleName1;
                //G4cout << " SrcRegionName " << SrcRegionName << " ParticleName " << ParticleName  << G4"\n" ;
                SrcRegionName = QString::fromStdString(SrcRegionName1);
                ParticleName = QString::fromStdString(ParticleName1);
                Quantity = QString::fromStdString(Quantity1);

                if(ParticleName == "RadioTracer"){

                }else{
                    LineString >> Geometry1 >> EnergyUnit ;
                    Geometry = QString::fromStdString(Geometry1);

                    double Enee;
                    int numWords = 0;
                    while (LineString >> Enee) {
                        if(EnergyUnit == "keV" || EnergyUnit == "KeV" || EnergyUnit == "KEV" || EnergyUnit == "kev") {Enee = Enee * 0.001;}
                        else if(EnergyUnit == "eV" || EnergyUnit == "EV" || EnergyUnit == "Ev" || EnergyUnit == "ev") {Enee = Enee * 0.000001;}
                        particleEnergies.push_back(Enee);
                        //if(SrcRegionName == "Liver" && (ParticleName == "neutron" || ParticleName == "alpha" || ParticleName == "e-")){
                        //    QTextStream(stdout) << ParticleName << " - Enee " <<  Enee  <<  "\n" ;
                        //}
                        ++numWords;
                    }
                    NumOfRefEne = numWords;
                }


                bool isin = false;
                for (int dd = 0 ; dd < PhantomSourcesMap[Geometry].size(); dd++) {
                    if(SrcRegionName == PhantomSourcesMap[Geometry][dd]){isin = true;break;}}
                if(isin == false){PhantomSourcesMap[Geometry].push_back(SrcRegionName);}

                isADataLines = true;
            }

            // read line of data that contains the name of target organ and the SAF correspondants
            if(isADataLines == true && word != "******"){

                if (word == "*") {

                    isADataLines = false;
                    continue;
                }

                if (word == "#") {
                    continue;
                }

                organTargetname = word;
                //std::cout << " - SrcRegionName " << SrcRegionName << " - organTargetname " << organTargetname << " - ParticleName " << ParticleName << std::"\n" ;
                if(ParticleName == "RadioTracer"){

                }
                else{

                    for( int zas = 0 ; zas < NumOfRefEne ; zas++ ){

                        LineString >> Val ;
                        if(Quantity == "SAF"){ if(QuantityUnit == "g-1" || QuantityUnit == "G-1" || QuantityUnit == "gram-1" || QuantityUnit == "Gram-1" || QuantityUnit == "GRAM-1") {Val = Val * 1000;}}
                        if(Quantity == "S"){if(QuantityUnit == "rad/uCi-h" || QuantityUnit == "RAD/UCI-H" || QuantityUnit == "rad/uci-h") {Val = Val*1e-2/(37000000000*1e-6*3.6e+3);}}
                        ICRPSAFs[Quantity][Geometry][ParticleName][particleEnergies[zas]][SrcRegionName][organTargetname] = Val;

                        bool isin = false;
                        for (int dd = 0 ; dd < SourceParticleEnergyValues[SrcRegionName][ParticleName].size(); dd++) {
                            if(particleEnergies[zas] == SourceParticleEnergyValues[SrcRegionName][ParticleName][dd]){isin = true;break;}}
                        if(isin == false){SourceParticleEnergyValues[SrcRegionName][ParticleName].push_back(particleEnergies[zas]);}

                        //if(SrcRegionName == "Liver" && organTargetname == "Liver" && ParticleName == "neutron" ){
                        //  QTextStream(stdout) << " - SrcRegionName " << SrcRegionName << " - organTargetname " << organTargetname << " - ParticleName " << ParticleName << " particleEnergies[zas] " << particleEnergies[zas] << " - Val " <<  Val  <<  "\n" ;
                        //}
                    }
                }
            }
        }

        ui->progressBarReadingCalcData->setValue(100);

        if(ICRPSAFs.size() == 0){
            ui->pushButtonReadICRPData->setText("Read ICRP133 SAFs");
            ui->pushButton_ReadDoseCalcsSAFs->setText("Read DoseCalcs SAFs");
            ui->pushButtonReadUSERData->setText("Read User SAFs");
        }else{
            ui->pushButtonReadICRPData->setText("Read ICRP133 SAFs *");
            ui->pushButton_ReadDoseCalcsSAFs->setText("Read DoseCalcs SAFs");
            ui->pushButtonReadUSERData->setText("Read User SAFs");

        }
    }
    else{
        showResultsOutput("Canno't Open Reference file... " + FilePath, 4);
    }

    ReadMassesAndWTFactor(ICRPDATAPath+"/ICRP110RegionsData");

}
void MainWindow::ReadMassesAndWTFactor(QString FilePath){

    std::ifstream file(FilePath.toStdString() , std::ios::binary);

    TissueFactorMap.clear();
    SpecialOrganSoftMassFraction.clear();

    if(file.is_open()){

        double WT, frac, frac1, Mass1 , Mass2;

        std::string line, indicator, organ, word;

        while (getline(file, line)) {

            //G4cout << " the line " << line << G4"\n" ;

            std::istringstream A(line);

            if(A.str().empty()){
                continue;
            }

            //QTextStream(stdout) << word.c_str() << " "<< Mass1 <<  " " << Mass2 <<"\n";

            A >> word ;
            if(word == "#"){
                continue;
            }
            else if(word == "Source"){
                indicator = "Source";
                continue;
            } else if (word == "Target"){
                indicator = "Target";
                continue;
            } else if (word == "WT-factor"){
                indicator = "WT-factor";
                continue;
            } else if (word == "OtherTissues"){
                //QTextStream(stdout) << " the word " << word.c_str() << "\n" ;
                indicator = "OtherTissues";
                continue;
            } else if (word == "TotalBody"){
                //QTextStream(stdout) << " the word " << word.c_str() << "\n" ;
                indicator = "TotalBody";
                continue;
            } else if (word == "NewSourceRegions"){
                //QTextStream(stdout) << " the word " << word.c_str() << "\n" ;
                indicator = "NewSourceRegions";
                continue;
            } else{

                if (indicator == "WT-factor"){

                    CurrentTargets.clear();
                    CurrentTargetsFractions.clear();
                    RegionFraction.clear();

                    organ = word;
                    A >> WT ;
                    if(organ.c_str() == "#"){
                        continue;
                    }

                    //QTextStream(stdout) << organ.c_str() << " " << WT << " ";

                    while (A >> word && A >> frac){
                        CurrentTargets.push_back(word.c_str());
                        CurrentTargetsFractions.push_back(frac);
                        RegionFraction[organ.c_str()][word.c_str()] = frac;
                        //QTextStream(stdout) << word.c_str() << " "<< frac <<  " ";
                    }

                    //QTextStream(stdout) <<" ------ \n";

                    bool isin = false;

                    for ( auto it = RegionParameterValueMap.begin(); it != RegionParameterValueMap.end(); ++it  ){
                        QString Geometry_NAME = it.key();
                        if(RegionParameterValueMap[Geometry_NAME]["Mass"][organ.c_str()] != 0.){
                            isin = true;break;
                        }
                    }
                    if(isin == false){
                        GenerateSAFFromNewTarget(organ.c_str());
                    }
                    TissueFactorMap[organ.c_str()] = WT;
                }
                else if (indicator == "NewSourceRegions"){

                    organ = word;

                    A >> word ;// add or all_except
                    CurrentSources.clear();
                    while (A >> word ){
                        CurrentSources.push_back(word.c_str());

                        //QTextStream(stdout) << word.c_str() << " "<< frac <<  " ";
                    }
                    CreateNewSourceRegion(organ.c_str());

                }
                else if (indicator == "OtherTissues"){
                    A >> frac >> frac1;
                    SpecialOrganSoftMassFraction["OtherTissues"]["ICRPAdultMale"][word.c_str()] = frac;
                    SpecialOrganSoftMassFraction["OtherTissues"]["ICRPAdultFemale"][word.c_str()] = frac1;
                    //QTextStream(stdout) << " organ " << word.c_str() << " frac " << frac << " frac1 " << frac1 <<"\n";
                }
                else if (indicator == "TotalBody"){
                    A >> frac >> frac1;
                    SpecialOrganSoftMassFraction["TotalBody"]["ICRPAdultMale"][word.c_str()] = frac;
                    SpecialOrganSoftMassFraction["TotalBody"]["ICRPAdultFemale"][word.c_str()] = frac1;
                    //QTextStream(stdout) << " organ " << word.c_str() << " frac " << frac << " frac1 " << frac1 <<"\n";
                }
                else if(indicator == "Target"){
                    //QTextStream(stdout) << word.c_str() << " "<< Mass1 <<  " " << Mass2 <<"\n";
                    A >> Mass1 >> Mass2 ;

                    RegionParameterValueMap["ICRPAdultMale"]["Mass"][word.c_str()] = Mass1;
                    RegionParameterValueMap["ICRPAdultFemale"]["Mass"][word.c_str()] = Mass2;

                }else if(indicator == "Source"){
                    A >> Mass1 >> Mass2 ;

                    ICRPSourceMassMap["ICRPAdultMale"]["Mass"][word.c_str()] = Mass1;
                    ICRPSourceMassMap["ICRPAdultFemale"]["Mass"][word.c_str()] = Mass2;
                }
            }
        }

        QMap<QString,QMap<QString,QMap<QString,double>>> RegionParameterValueMap1 = RegionParameterValueMap;
        for ( auto it = RegionParameterValueMap1.begin(); it != RegionParameterValueMap1.end(); ++it  ){
            QString Geometry_NAME = it.key();
            if(Geometry_NAME.contains("Female")){
                RegionParameterValueMap[Geometry_NAME]["Mass"].remove("BoneEndosteumAdultMale") ;
                RegionParameterValueMap[Geometry_NAME]["Mass"].remove("ActiveBoneMarrowAdultMale") ;
                RegionParameterValueMap[Geometry_NAME]["Mass"].remove("Prostate") ;
            }
            else if(Geometry_NAME.contains("Male")){
                RegionParameterValueMap[Geometry_NAME]["Mass"].remove("BoneEndosteumAdultFemale") ;
                RegionParameterValueMap[Geometry_NAME]["Mass"].remove("ActiveBoneMarrowAdultFemale") ;
                RegionParameterValueMap[Geometry_NAME]["Mass"].remove("Uterus") ;
            }
        }
    }
}
void MainWindow::ReadMassesAndWTFactorForDCC(QString FilePath){

    std::ifstream file(FilePath.toStdString() , std::ios::binary);

    TissueFactorMap.clear();
    //SpecialOrganSoftMassFraction.clear();

    if(file.is_open()){

        double WT, frac, frac1, Mass1 , Mass2;

        std::string line, indicator, organ, word;

        while (getline(file, line)) {

            //G4cout << " the line " << line << G4"\n" ;

            std::istringstream A(line);

            if(A.str().empty()){
                continue;
            }

            //QTextStream(stdout) << word.c_str() << " "<< Mass1 <<  " " << Mass2 <<"\n";

            A >> word ;
            if(word == "#"){
                continue;
            }
            else if(word == "Source"){
                indicator = "Source";
                continue;
            } else if (word == "Target"){
                indicator = "Target";
                continue;
            } else if (word == "WT-factor"){
                indicator = "WT-factor";
                continue;
            } else if (word == "OtherTissues"){
                //QTextStream(stdout) << " the word " << word.c_str() << "\n" ;
                indicator = "OtherTissues";
                continue;
            } else if (word == "TotalBody"){
                //QTextStream(stdout) << " the word " << word.c_str() << "\n" ;
                indicator = "TotalBody";
                continue;
            } else if (word == "NewSourceRegions"){
                //QTextStream(stdout) << " the word " << word.c_str() << "\n" ;
                indicator = "NewSourceRegions";
                continue;
            } else{

                if (indicator == "WT-factor"){

                    CurrentTargets.clear();
                    CurrentTargetsFractions.clear();
                    RegionFraction.clear();

                    organ = word;
                    A >> WT ;
                    if(organ.c_str() == "#"){
                        continue;
                    }

                    //QTextStream(stdout) << organ.c_str() << " " << WT << " ";

                    while (A >> word && A >> frac){
                        CurrentTargets.push_back(word.c_str());
                        CurrentTargetsFractions.push_back(frac);
                        RegionFraction[organ.c_str()][word.c_str()] = frac;
                        //QTextStream(stdout) << word.c_str() << " "<< frac <<  " ";
                    }

                    //QTextStream(stdout) <<" ------ \n";

                    bool isin = false;

                    for ( auto it = RegionParameterValueMap.begin(); it != RegionParameterValueMap.end(); ++it  ){
                        QString Geometry_NAME = it.key();
                        if(RegionParameterValueMap[Geometry_NAME]["Mass"][organ.c_str()] != 0.){
                            isin = true;break;
                        }
                    }
                    if(isin == false){
                        GenerateDCCFromNewTarget(organ.c_str());
                    }
                    TissueFactorMap[organ.c_str()] = WT;
                }
                /*
                else if (indicator == "NewSourceRegions"){

                    organ = word;

                    A >> word ;// add or all_except
                    CurrentSources.clear();
                    while (A >> word ){
                        CurrentSources.push_back(word.c_str());

                        //QTextStream(stdout) << word.c_str() << " "<< frac <<  " ";
                    }
                    CreateNewSourceRegion(organ.c_str());

                }
                else if (indicator == "OtherTissues"){
                    A >> frac >> frac1;
                    SpecialOrganSoftMassFraction["OtherTissues"]["ICRPAdultMale"][word.c_str()] = frac;
                    SpecialOrganSoftMassFraction["OtherTissues"]["ICRPAdultFemale"][word.c_str()] = frac1;
                    //QTextStream(stdout) << " organ " << word.c_str() << " frac " << frac << " frac1 " << frac1 <<"\n";
                }
                else if (indicator == "TotalBody"){
                    A >> frac >> frac1;
                    SpecialOrganSoftMassFraction["TotalBody"]["ICRPAdultMale"][word.c_str()] = frac;
                    SpecialOrganSoftMassFraction["TotalBody"]["ICRPAdultFemale"][word.c_str()] = frac1;
                    //QTextStream(stdout) << " organ " << word.c_str() << " frac " << frac << " frac1 " << frac1 <<"\n";
                }
                else if(indicator == "Source"){
                    A >> Mass1 >> Mass2 ;

                    ICRPSourceMassMap["ICRPAdultMale"]["Mass"][word.c_str()] = Mass1;
                    ICRPSourceMassMap["ICRPAdultFemale"]["Mass"][word.c_str()] = Mass2;
                }
                */
                else if(indicator == "Target"){
                    //QTextStream(stdout) << word.c_str() << " "<< Mass1 <<  " " << Mass2 <<"\n";
                    A >> Mass1 >> Mass2 ;

                    RegionParameterValueMap["ICRPAdultMale"]["Mass"][word.c_str()] = Mass1;
                    RegionParameterValueMap["ICRPAdultFemale"]["Mass"][word.c_str()] = Mass2;

                }
            }
        }

        QMap<QString,QMap<QString,QMap<QString,double>>> RegionParameterValueMap1 = RegionParameterValueMap;
        for ( auto it = RegionParameterValueMap1.begin(); it != RegionParameterValueMap1.end(); ++it  ){
            QString Geometry_NAME = it.key();
            if(Geometry_NAME.contains("Female")){
                RegionParameterValueMap[Geometry_NAME]["Mass"].remove("BoneEndosteumAdultMale") ;
                RegionParameterValueMap[Geometry_NAME]["Mass"].remove("ActiveBoneMarrowAdultMale") ;
                RegionParameterValueMap[Geometry_NAME]["Mass"].remove("Prostate") ;
            }
            else if(Geometry_NAME.contains("Male")){
                RegionParameterValueMap[Geometry_NAME]["Mass"].remove("BoneEndosteumAdultFemale") ;
                RegionParameterValueMap[Geometry_NAME]["Mass"].remove("ActiveBoneMarrowAdultFemale") ;
                RegionParameterValueMap[Geometry_NAME]["Mass"].remove("Uterus") ;
            }
        }
    }
}


void MainWindow::ReadWTFactor(QString FilePath){

    std::ifstream file(FilePath.toStdString() , std::ios::binary);

    if(file.is_open()){

        double WT, frac, frac1, Mass1 , Mass2;

        std::string line, indicator, organ, word;

        if(!ui->checkBoxFixScoreCommands->isChecked()){
            return;
        }
        QString sstt= "";

        while (getline(file, line)) {

            //G4cout << " the line " << line << G4"\n" ;

            std::istringstream A(line);

            if(A.str().empty()){
                continue;
            }

            //QTextStream(stdout) << word.c_str() << " "<< Mass1 <<  " " << Mass2 <<"\n";

            A >> word ;
            if(word == "#"){
                continue;
            }
            else if (word == "WT-factor"){
                indicator = "WT-factor";
                continue;
            }
            else{
                if (indicator == "WT-factor"){

                    organ = word;  QString orr = organ.c_str();
                    A >> WT ;

                    sstt += orr + " " + QString::number(WT) + " ";

                    if(organ.c_str() == "#"){
                        continue;
                    }

                    //QTextStream(stdout) <<" ------ \n";
                }
            }
        }
        ui->TissueFactorLineEdit->setText(sstt);

    }
}
void MainWindow::ConstructOtherTissuesSource(){

    QStringList OtherTissuesSourceOrgans;

    if(ui->pushButton_ReadDoseCalcsSAFs->text() == "Read DoseCalcs SAFs *"){
        return;
    }

    QString Quantity_NAME = "SAF";
    QString Radiotracer_Name = ui->comboBoxRadioPharmaceutiques->currentText();
    QString Geometry_NAME = ui->comboBoxPhantom->currentText();

    QString OtherTissuesSource_NAME = "OtherTissues";
    for ( auto AAA = RadioTracerSourceOrganResidenceTime[Radiotracer_Name][RadiotracerradionucleidMap[Radiotracer_Name]].begin(); AAA != RadioTracerSourceOrganResidenceTime[Radiotracer_Name][RadiotracerradionucleidMap[Radiotracer_Name]].end(); ++AAA  ){
        if("OtherTissues" == AAA.key()){OtherTissuesSource_NAME = "OtherTissues"; break;}
        else if("TotalBody" == AAA.key()){OtherTissuesSource_NAME = "TotalBody"; break;}
        RegionParameterValueMap[Geometry_NAME]["Mass"][OtherTissuesSource_NAME] = 0.;
    }


    QMap<QString,QMap<double,QMap<QString,QMap<QString,double>>>> SAFMap = ICRPSAFs[Quantity_NAME][Geometry_NAME];

    // remove the older SAF of OtherTissues
    for ( auto it2 = SAFMap.begin(); it2 != SAFMap.end(); ++it2  ){
        QString Particle_NAME = it2.key();
        for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
            double Energy_Val = it3.key();
            it3.value()[OtherTissuesSource_NAME].clear();
            ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][OtherTissuesSource_NAME].clear();
        }
    }
    //QTextStream(stdout) << "---NewSourceTissue_NAME " << OtherTissuesSource_NAME <<"\n";

    // choose the sources according to radiopharmaceutical and geometry, with elimination of some sources that are not soft tissues
    for ( auto it2 = SAFMap.begin(); it2 != SAFMap.end(); ++it2  ){
        QString Particle_NAME = it2.key();
        // QTextStream(stdout) << "---Particle_NAME " << Particle_NAME <<"\n";

        for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
            double Energy_Val = it3.key();
            // QTextStream(stdout) << "----Energy_Val " << Energy_Val <<"\n";

            for ( auto DD = it3.value().begin(); DD != it3.value().end(); ++DD  ){
                QString Source_NAME  = DD.key();

                // eliminate radiopharmaceuticals source organs
                bool isin = false;
                for ( auto AAA = RadioTracerSourceOrganResidenceTime[Radiotracer_Name][RadiotracerradionucleidMap[Radiotracer_Name]].begin(); AAA != RadioTracerSourceOrganResidenceTime[Radiotracer_Name][RadiotracerradionucleidMap[Radiotracer_Name]].end(); ++AAA  ){
                    if(Source_NAME == AAA.key()){isin = true;break;}
                }
                if(isin == true){
                    continue;
                }

                // eliminate repeated source organs
                isin = false;
                for (int dd = 0 ; dd < OtherTissuesSourceOrgans.size(); dd++) {
                    if(Source_NAME == OtherTissuesSourceOrgans[dd]){isin = true;break;}}
                if(isin == false){
                    OtherTissuesSourceOrgans.push_back(Source_NAME);
                }

                //QTextStream(stdout) << "---Registered Source_NAME " << Source_NAME <<"\n";

                // eliminate the non soft tissue source organs
            }
        }
    }
    SAFMap.clear();

    // calculate the total mass of OtherTissues source organs
    for (int dd = 0 ; dd < OtherTissuesSourceOrgans.size(); dd++) {
        RegionParameterValueMap[Geometry_NAME]["Mass"][OtherTissuesSource_NAME] +=
                ICRPSourceMassMap[Geometry_NAME]["Mass"][OtherTissuesSourceOrgans[dd]]*
                SpecialOrganSoftMassFraction[OtherTissuesSource_NAME][Geometry_NAME][OtherTissuesSourceOrgans[dd]];
        //QTextStream(stdout) << " SourceRegion " << OtherTissuesSourceOrgans[dd] << " Mass " << ICRPSourceMassMap[Geometry_NAME]["Mass"][OtherTissuesSourceOrgans[dd]] << " FRACTION "<< SpecialOrganSoftMassFraction[OtherTissuesSource_NAME][Geometry_NAME][OtherTissuesSourceOrgans[dd]] <<  " mass ValTotal " << RegionParameterValueMap[Geometry_NAME]["Mass"][OtherTissuesSource_NAME] <<"\n";
    }

    QMap<QString,QMap<double,QMap<QString,QMap<QString,double>>>> SAFMap1 = ICRPSAFs[Quantity_NAME][Geometry_NAME];

    for ( auto it2 = SAFMap1.begin(); it2 != SAFMap1.end(); ++it2  ){
        QString Particle_NAME = it2.key();
        // QTextStream(stdout) << "---Particle_NAME " << Particle_NAME <<"\n";

        for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
            double Energy_Val = it3.key();
            // QTextStream(stdout) << "----Energy_Val " << Energy_Val <<"\n";

            for (int dd = 0 ; dd < OtherTissuesSourceOrgans.size(); dd++) {

                for ( auto DD = it3.value()[OtherTissuesSourceOrgans[dd]].begin(); DD != it3.value()[OtherTissuesSourceOrgans[dd]].end(); ++DD  ){
                    QString  Target_NAME  = DD.key();
                    double SAFValue  = DD.value();
                    if( !__isnan(SAFValue) && !__isinf(SAFValue) && SAFValue != 0 && SAFValue != NULL){
                    }else{
                        continue;
                    }
                    ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][OtherTissuesSource_NAME][Target_NAME] +=
                            SAFValue*
                            (
                                (SpecialOrganSoftMassFraction[OtherTissuesSource_NAME][Geometry_NAME][OtherTissuesSourceOrgans[dd]]*
                                ICRPSourceMassMap[Geometry_NAME]["Mass"][OtherTissuesSourceOrgans[dd]])/
                            RegionParameterValueMap[Geometry_NAME]["Mass"][OtherTissuesSource_NAME]
                            )
                            ;

                    bool isin = false;
                    for (int dd = 0 ; dd < SourceParticleEnergyValues[OtherTissuesSource_NAME][Particle_NAME].size(); dd++) {
                        if(Energy_Val == SourceParticleEnergyValues[OtherTissuesSource_NAME][Particle_NAME][dd]){isin = true;break;}}
                    if(isin == false){SourceParticleEnergyValues[OtherTissuesSource_NAME][Particle_NAME].push_back(Energy_Val);}

                    //if(Target_NAME == "Liver" && Particle_NAME == "e-"){
                    //    QTextStream(stdout) << "----SRC " << OtherTissuesSourceOrgans[dd]  << " Target_NAME " << Target_NAME << " SAFValue " << SAFValue << " Related_SAFValue_Othertissues " << ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][OtherTissuesSource_NAME][Target_NAME] <<"\n";
                    //}

                }

            }
        }
    }

    // to add the OtherTissues name to the list of sources
    bool isin = false;
    for (int dd = 0 ; dd < ui->comboBoxSources->count(); dd++) {
        if(OtherTissuesSource_NAME == ui->comboBoxSources->itemText(dd)){isin = true;break;}}
    if(isin == false){
        ui->comboBoxSources->addItem(OtherTissuesSource_NAME);
    }

}
void MainWindow::CreateNewSourceRegion(QString SSOURCE){

    QString Quantity_NAME = "SAF";

    QMap<QString,QMap<QString,QMap<double,QMap<QString,QMap<QString,double>>>>> SAFMap = ICRPSAFs[Quantity_NAME];

    QMap<QString,QMap<QString,QMap<QString,double>>> RegionParameterValueMap1 = RegionParameterValueMap;
    for ( auto it = RegionParameterValueMap.begin(); it != RegionParameterValueMap.end(); ++it  ){
        for (int dd = 0 ; dd < CurrentSources.size(); dd++) {
            QString Geometry_NAME = it.key();
            RegionParameterValueMap[Geometry_NAME]["Mass"][SSOURCE] += RegionParameterValueMap1[Geometry_NAME]["Mass"][CurrentSources[dd]];

            bool isin = false;
            for (int dd = 0 ; dd < PhantomSourcesMap[Geometry_NAME].size(); dd++) {
                if(SSOURCE == PhantomSourcesMap[Geometry_NAME][dd]){isin = true;break;}}
            if(isin == false){PhantomSourcesMap[Geometry_NAME].push_back(SSOURCE);}
        }
    }

    ui->progressBarReadingCalcData->setRange(0, 100);
    ui->progressBarReadingCalcData->setValue(0);
    ui->progressBarReadingCalcData->show();
    double percent = 0;
    int qq = 0;

    for ( auto it = SAFMap.begin(); it != SAFMap.end(); ++it  ){

        QString Geometry_NAME = it.key();
        //        QTextStream(stdout) << "--Geometry_NAME " << Geometry_NAME <<"\n";


        percent = ((qq+1)/it->size())*100;
        ui->progressBarReadingCalcData->setValue(percent);
        qq++;

        for ( auto it2 = it.value().begin(); it2 != it.value().end(); ++it2  ){
            QString Particle_NAME = it2.key();
            //            QTextStream(stdout) << "---Particle_NAME " << Particle_NAME <<"\n";

            for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
                double Energy_Val = it3.key();
                //                QTextStream(stdout) << "----Energy_Val " << Energy_Val <<"\n";

                for ( auto DD = it3.value().begin(); DD != it3.value().end(); ++DD  ){
                    QString Source_NAME  = DD.key();

                    bool isin = false;
                    for (int dd = 0 ; dd < CurrentSources.size(); dd++) {
                        if(Source_NAME == CurrentSources[dd]){
                            isin = true;break;
                        }
                    }
                    if(isin == false){
                        continue;
                    }
                    //QTextStream(stdout) << "-----Source_NAME " << Source_NAME << "-----fraction " << RegionFraction[SSOURCE] <<"\n";

                    for ( auto CC = DD.value().begin(); CC != DD.value().end(); ++CC  ){
                        QString  Target_NAME  = CC.key();
                        double SAFValue  = CC.value();

                        ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][SSOURCE][Target_NAME] += SAFValue*(RegionParameterValueMap[Geometry_NAME]["Mass"][Source_NAME]/RegionParameterValueMap[Geometry_NAME]["Mass"][SSOURCE]);
                        ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][SSOURCE][SSOURCE] += SAFValue*(RegionParameterValueMap[Geometry_NAME]["Mass"][Source_NAME]/RegionParameterValueMap[Geometry_NAME]["Mass"][SSOURCE]);

                        bool isin = false;
                        for (int dd = 0 ; dd < SourceParticleEnergyValues[SSOURCE][Particle_NAME].size(); dd++) {
                            if(Energy_Val == SourceParticleEnergyValues[SSOURCE][Particle_NAME][dd]){isin = true;break;}}
                        if(isin == false){SourceParticleEnergyValues[SSOURCE][Particle_NAME].push_back(Energy_Val);}
                    }
                }
            }
        }
    }

    ui->progressBarReadingCalcData->setValue(100);

}








// External Dosimetry coefficients
void MainWindow::on_pushButtonReadDCCsDoseCalcs_clicked()
{
    Read_DoseCalcs_file  (ICRPDATAPath+"/DoseCalcsDCCsData");
    ReadMassesAndWTFactor(ICRPDATAPath+"/DoseCalcsRegionsData");

    CalculateExternalDosimetryQuantities();

    //ui->pushButton_ReadDoseCalcsSAFs->setToolTip("Current readed file path: "+ICRPDATAPath+"/DoseCalcsSAFsData");

}
void MainWindow::on_pushButtonGenerateExternalDosValues_clicked()
{

    QString Quantity_NAME = ui->comboBoxExternalDosQuantities->currentText();
    //double convfac = QuantitiesConversionFromDefault[ui->comboBoxQuantityNucl->currentText()][ui->comboBoxEffDoseUnit->currentText()];
    double convfac = 1;
    QString Geometry_NAME    = ui->comboBoxPhantom->currentText();
    QString Particle_NAME    = ui->comboBoxExternalDoseParticles->currentText();
    QString Source_NAME      = ui->comboBoxSources->currentText();
    QString Target_NAME      = ui->comboBoxTargets->currentText();
    double  Energy_Val;


    ui->tableWidgetForOneGraph->clear();
    ui->tableWidgetForOneGraph->setRowCount(0);
    ui->tableWidgetForOneGraph->setColumnCount(0);


    QStringList headers;

    headers.append("Energy\Source");

    QVector<QString> Sources;
    if(ui->comboBoxExternalDosQuantities->currentText() == "Ambient Dose"){
        headers.append("Ambient Dose");
        Sources.push_back("Ambient Dose");
    }else{
        for(int b=0; b < PhantomSourcesMap[Geometry_NAME].size();b++){
            headers.append(PhantomSourcesMap[Geometry_NAME][b]);
            Sources.push_back(PhantomSourcesMap[Geometry_NAME][b]);
        }
    }

    QVector<double> energies;
    for(int b=0; b < PhantomSourcesMap[Geometry_NAME].size();b++){
        for(int a=0; a < SourceParticleEnergyValues[PhantomSourcesMap[Geometry_NAME][b]][Particle_NAME].size();a++){

            double particleE = SourceParticleEnergyValues[PhantomSourcesMap[Geometry_NAME][b]][Particle_NAME][a];
            bool isin = false;
            for (int dd = 0 ; dd < energies.size(); dd++) {
                if(particleE == energies[dd]){isin = true;break;}}
            if(isin == false){energies.push_back(particleE);}
        }
    }

    ui->tableWidgetForOneGraph->setColumnCount(PhantomSourcesMap[Geometry_NAME].size()+1);
    ui->tableWidgetForOneGraph->setShowGrid(true);
    ui->tableWidgetForOneGraph->setSelectionMode(QAbstractItemView::SingleSelection);
    ui->tableWidgetForOneGraph->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->tableWidgetForOneGraph->setHorizontalHeaderLabels(headers);
    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    //ui->tableWidgetForOneGraph->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    ui->tableWidgetForOneGraph->resizeColumnsToContents();
    //ui->tableWidgetForOneGraph->horizontalHeader()->sectionResizeMode(0, QHeaderView::Stretch);

    //QTextStream(stdout) << "Table to show rows " << Data.size() << "\n";

    int row = 0;

    for(int a=0; a < energies.size();a++){

        ui->tableWidgetForOneGraph->insertRow(row);
        ui->tableWidgetForOneGraph->setItem(row,0, new QTableWidgetItem(QString::number(energies[a])));
        for(int b=0; b < Sources.size();b++){

            double value;

            if(ui->comboBoxExternalDosQuantities->currentText() == "DCC"){
                value = ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][energies[a]][Source_NAME][Target_NAME];
            }
            else if(ui->comboBoxExternalDosQuantities->currentText() == "Ambient Dose"){
                value = ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][energies[a]]["Ambient Dose"]["Ambient Dose"];
            }
            else{
                value = ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][energies[a]][Source_NAME][ui->comboBoxExternalDosQuantities->currentText()];
            }

            //QTextStream(stdout) << "Col "<< b+1 << " Particle " << Particles[b] << " val " << QString::number(DATAMap[Quantity_NAME][Geometry_NAME][RadioTracer_NAME][Particles[b]][Targets[a]][Targets[a]]/convfac) <<"\n";
            ui->tableWidgetForOneGraph->setItem(row,b+1, new QTableWidgetItem(QString::number(value/convfac)));
        }
        row++;
    }

    ui->tableWidgetForOneGraph->resizeColumnsToContents();
    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->viewport()->installEventFilter(this);

}
void MainWindow::GenerateDCCFromNewTarget(QString TTARGET){

    QString Quantity_NAME = "DCC";

    QMap<QString,QMap<QString,QMap<double,QMap<QString,QMap<QString,double>>>>> SAFMap = ICRPSAFs[Quantity_NAME];

    ui->progressBarReadingCalcData->setRange(0, 100);
    ui->progressBarReadingCalcData->setValue(0);
    ui->progressBarReadingCalcData->show();
    double percent = 0;
    int qq = 0;

    for ( auto it = SAFMap.begin(); it != SAFMap.end(); ++it  ){

        QString Geometry_NAME = it.key();
        //        QTextStream(stdout) << "--Geometry_NAME " << Geometry_NAME <<"\n";

        percent = ((qq+1)/it->size())*100;
        ui->progressBarReadingCalcData->setValue(percent);
        qq++;

        for ( auto it2 = it.value().begin(); it2 != it.value().end(); ++it2  ){
            QString Particle_NAME = it2.key();
            //            QTextStream(stdout) << "---Particle_NAME " << Particle_NAME <<"\n";

            for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
                double Energy_Val = it3.key();
                //                QTextStream(stdout) << "----Energy_Val " << Energy_Val <<"\n";

                for ( auto DD = it3.value().begin(); DD != it3.value().end(); ++DD  ){
                    QString Source_NAME  = DD.key();

                    //QTextStream(stdout) << "-----Source_NAME " << Source_NAME <<"\n";
                    for ( auto CC = DD.value().begin(); CC != DD.value().end(); ++CC  ){
                        QString  Target_NAME  = CC.key();


                        bool isin = false;
                        for (int dd = 0 ; dd < CurrentTargets.size(); dd++) {
                            if(Target_NAME == CurrentTargets[dd]){
                                isin = true;break;
                            }
                        }
                        if(isin == false){
                            continue;
                        }

                        //QTextStream(stdout) << "-----Target_NAME " << Target_NAME << "-----fraction " << RegionFraction[TTARGET] <<"\n";

                        double SAFValue  = CC.value();

                        ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][Source_NAME][TTARGET] += SAFValue*RegionFraction[TTARGET][Target_NAME];
                    }
                }
            }
        }
    }

    ui->progressBarReadingCalcData->setValue(100);

    QMap<QString,QMap<QString,QMap<QString,double>>> RegionParameterValueMap1 = RegionParameterValueMap;
    for ( auto it = RegionParameterValueMap.begin(); it != RegionParameterValueMap.end(); ++it  ){
        QString Geometry_NAME = it.key();
        for (int dd = 0 ; dd < CurrentTargets.size(); dd++) {
            RegionParameterValueMap[Geometry_NAME]["Mass"][TTARGET] += RegionParameterValueMap1[Geometry_NAME]["Mass"][CurrentTargets[dd]]*RegionFraction[TTARGET][CurrentTargets[dd]];
        }

        // The WT are defined in a file for DoseCalcs ICRP Phantoms
        //TissueFactorMap[TTARGET] += TissueFactorMap[CurrentTargets[dd]];
        //TissueFactorMap[CurrentTargets[dd]]=0;
    }
}
void MainWindow::CalculateExternalDosimetryQuantities(){

    //double convfac = QuantitiesConversionFromDefault[ui->comboBoxQuantityNucl->currentText()][ui->comboBoxEffDoseUnit->currentText()];
    double convfac = 1;
    QString Quantity_NAME = "DCC";
    QString Geometry_NAME    = ui->comboBoxPhantom->currentText();
    QString Particle_NAME    = ui->comboBoxExternalDoseParticles->currentText();
    QString Source_NAME      = ui->comboBoxSources->currentText();
    QString Target_NAME      = ui->comboBoxTargets->currentText();
    double  Energy_Val;
    double Val;


    QMap<QString,QMap<QString,QMap<double,QMap<QString,QMap<QString,double>>>>> SAFMap = ICRPSAFs[Quantity_NAME];

    for ( auto it = SAFMap.begin(); it != SAFMap.end(); ++it  ){

        QString Geometry_NAME = it.key();
        //        QTextStream(stdout) << "--Geometry_NAME " << Geometry_NAME <<"\n";

        for ( auto it2 = it.value().begin(); it2 != it.value().end(); ++it2  ){
            QString Particle_NAME = it2.key();
            //            QTextStream(stdout) << "---Particle_NAME " << Particle_NAME <<"\n";

            for ( auto it3 = it2.value().begin(); it3 != it2.value().end(); ++it3  ){
                double Energy_Val = it3.key();
                //                QTextStream(stdout) << "----Energy_Val " << Energy_Val <<"\n";
                double ambienteffedose = 0;

                for ( auto DD = it3.value().begin(); DD != it3.value().end(); ++DD  ){
                    QString Source_NAME  = DD.key();
                    //QTextStream(stdout) << "-----Source_NAME " << Source_NAME << "-----fraction " << RegionFraction[SSOURCE] <<"\n";

                    for ( auto CC = DD.value().begin(); CC != DD.value().end(); ++CC  ){
                        QString  Target_NAME  = CC.key();
                        double SAFValue  = CC.value();

                        ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][Source_NAME]["Ambient Dose"] += SAFValue*GenerateRadiationFactor(Particle_NAME,Energy_Val)*TissueFactorMap[Target_NAME];
                        ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][Source_NAME]["Personal Dose"] += SAFValue*GenerateRadiationFactor(Particle_NAME,Energy_Val)*TissueFactorMap[Target_NAME];

                        if(        Target_NAME == "ELL"
                                || Target_NAME == "EBL"
                                || Target_NAME == "ELR"
                                || Target_NAME == "EBR"

                                || Target_NAME == "Eye_lens_sensitive_left"
                                || Target_NAME == "Eye_lens_insensitive_left"
                                || Target_NAME == "Eye_lens_sensitive_right"
                                || Target_NAME == "Eye_lens_insensitive_right"
                                ){

                            ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][Source_NAME]["Directional or Personal Absorbed Dose in the Lens of the Eye"] += SAFValue;
                        }
                        if(        Target_NAME == "SkH"
                                || Target_NAME == "SkT"
                                || Target_NAME == "SkL"
                                || Target_NAME == "SkA"
                                || Target_NAME == "STB"

                                || Target_NAME == "Skin_head_insensitive"
                                || Target_NAME == "Skin_head_sensitive(50-100)"
                                || Target_NAME == "Skin_trunk_insensitive"
                                || Target_NAME == "Skin_trunk_sensitive(50-100)"
                                || Target_NAME == "Skin_arms_insensitive"
                                || Target_NAME == "Skin_arms_sensitive(50-100)"
                                || Target_NAME == "Skin_legs_insensitive"
                                || Target_NAME == "Skin_legs_sensitive(50-100)"

                                ){

                            ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][Source_NAME]["Directional or Personal Absorbed Dose in Local Skin"] += SAFValue;
                        }
                    }

                    if(ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][Source_NAME]["Ambient Dose"] > ambienteffedose){
                        ambienteffedose = ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val][Source_NAME]["Ambient Dose"];
                    }
                }

                ICRPSAFs[Quantity_NAME][Geometry_NAME][Particle_NAME][Energy_Val]["Ambient Dose"]["Ambient Dose"] = ambienteffedose;
            }
        }
    }














    if(ui->comboBoxExternalDosQuantities->currentText() == "DCC"){


    }
    else if(ui->comboBoxExternalDosQuantities->currentText() == "Ambient Dose"){



    }
    else if(ui->comboBoxExternalDosQuantities->currentText() == "Personal Dose"){


    }
    else if(ui->comboBoxExternalDosQuantities->currentText() == "Directional or Personal Absorbed Dose in the Lens of the Eye"){


    }
    else if(ui->comboBoxExternalDosQuantities->currentText() == "Directional or Personal Absorbed Dose in Local Skin"){


    }
}

