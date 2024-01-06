#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>
#include <QString>
#include "filesManager.h"
#include <QFile>
#include <QProcess>

#include "gui/customtextedit.h"
#include "gui/organschooserdialog.h"
#include "gui/installationDialog.h"
#include "gui/geometriesvisualization.h"
#include "gui/plotDialog.h"
#include "gui/visualizationManager.h"
#include "gui/geometrymodellingeditor.h"
#include "gui/highlighter.h"
#include "gui/runmpisystem.h"
#include "gui/filesManager.h"
#include "gui/terminal.h"

//#include "mainGeant4AppCore.h"

namespace Ui {

class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:

    QTimer timer;

    bool CheckInputsForGeometryFiles();

    void UpdateMatVolSolList(int);
    QStringList getArgumentOfACommandFromText(QString, int);
    QString changeArgumentOfACommandFromText(QString, QString, QString);

    void OpenROOTCanvas();
    void CreateROOTFile();
    void AddNewMaterialOrRegionToListOfSourceRegionNames();

    void UsePackagesMethods();
    void initializeVariable();
    void FillCoomponentByDefaultData();
    void CommandsInitialization();
    void updateApplicationTabs();
    void ReadConfigurationFile(QString);
    int FillComponentsFromInputsFile(QString);

    QVector <QString> ReadFileWithText(QString) ;
    void SaveDataFromInputComponents();

    void SaveGeometryInputs();
    void SaveSourcePhysicsInputs();
    void SaveROOTAnalysisInputs();

    void InitialiseGeometryInputs();
    void InitialiseSourcePhysicsInputs();
    void InitialiseROOTAnalysisInputs();

    void RefillGeometryInputs();
    void RefillSourcePhysicsInputs();
    void RefillROOTAnalysisInputs();

    void CreateUserCommands();
    QString generateInputUserTextForinputFile();
    bool TestSimulateExecutableInputsToRun();
    bool TestMergeExecutableInputsToRun();
    bool TestVisualizingInputsToRun();
    QString ShowMessageBoxAndGetPath(QString);

    bool ShowImportantSimulationData();

    void getRockcsDoseCalcsJobs();

    //CustomTextEdit * completingTextEdit;
    QCompleter *completer1 = nullptr;
    void setCompleters();

    QString QsubSeparatedMacroFilePath;

    bool StopOrCheck = false;

    QStringList VarToScorellist;

    QString SimulationsNumberDataString;
    QString QdelID;
    int NumberOfExecutions;
    int NumberOfSourceParticle;
    int NumberOfSourceEnergies;
    int NumberOfSourceRegions;
    int NumberOfSourceMomentumDirection;

    QString OpenMacrosFilePath;
    QString SaveMacrosFilePath;
    QString AnyOpenedFilePath;

    QStringList AngleUnits;
    QStringList SizeUnits;
    QStringList EnergyUnits;

    QMap <QString , QString> RadDataType;
    QMap <QString , QString> ValuesOfInputs;
    QMap <QString , double> UnitsConversionFactors;

    QVector<QCheckBox*> checkboxes;
    QVector<QDoubleSpinBox*> SpinBoxes;

    QComboBox *MatElemComb;
    QPushButton *btnSaveElem;
    QLineEdit *MaterialName;
    QLineEdit *MaterialID;
    QLineEdit *NistMaterialID;
    QLineEdit *MaterialDensity;
    QLineEdit *MaterialEleNumber;
    QCheckBox *NumFraCheckBox;
    QComboBox *ElemComb;
    QComboBox *RegionToVisualizeComb;
    QLineEdit *EleFraOrNumInMat;
    QPushButton *btnAddElem;
    QPushButton *btnAddMat;
    QComboBox* NistMatComb;
    QPushButton *btnAddNistMat;

    QLineEdit *XYZVoxelsNumb;
    QComboBox* ParamType;
    QLineEdit *XYZVoxelsHalfSize;

    QLineEdit *VoxContainerPos;
    QLineEdit *VoxContainerRot;
    QComboBox *DcmType;
    QPushButton *btnChooseDcmDir;
    QLineEdit *DcmDataDir;
    QPushButton *btnAddDcmTypeDir;
    QLineEdit *DcmCtNumber;
    QLineEdit *DcmCtDensity;
    QPushButton *btnAddDcmCTNumDen;

    QCheckBox *UseMaterialsAsRegionNames;
    QLineEdit *VoxRegionName;
    QLineEdit *VoxRegionMinMaxX;
    QLineEdit *VoxRegionMinMaxY;
    QLineEdit *VoxRegionMinMaxZ;
    QLineEdit *VoxRegionMinMaxDensity;
    QLineEdit *VoxelsegmentedMaterial;
    QPushButton *btnAddVoxelsRegionData;

    bool FirstTimeAdd = true;

    QComboBox * LimitsPlan ;
    QLineEdit* LimitsValues;
    QPushButton * btnAddVisualizationLimits ;

    QLineEdit* NumbOfDCMMat;

    QPushButton *btnAddVoxelsData;
    QPushButton *btnAddDcmMatOrdered;

    QString NotForcedVolumesCommand;


    QString DataAlreadyInEdit;

    QString DcmTypeDirCommands;
    QString DcmCtDensCommands;
    QString DcmNumbMatsAndMatsCommands;
    QString RegionDataCommands;
    QString MatCommandToAddElem;

    QComboBox *ComboxSolidType;
    QLineEdit *SolidSpecications;
    QLineEdit *SolidName;
    QComboBox *ComboxLogVolSolid;
    QComboBox *AllOrSpecificQComboBox;
    QComboBox *LogVolMatName;
    QComboBox *MaterialsIDForVox;
    QLineEdit *LogVolName;
    QComboBox *ComboxMotherVol;
    QLineEdit *PhyVolPos;
    QLineEdit *PhyVolRot;
    QLineEdit *PosRotUnits;
    QLineEdit *PhyVolName;
    QCheckBox *ViewVolCheckBox;
    QPushButton *BtnAddVol;
    QPushButton *BtnAddSol;
    QPushButton *BtnAddVolPath;
    QPushButton *BtnEditVolFile;
    QPushButton *BtnAddSolPath;

    QComboBox* DoseCalcsJobIDsCombobox;

    QStringList listOfRun;
    QStringList listOfRunIDs;
    QMap <QString,QString> listOfFilesName ;

    void showConstructMaterialsFrame();
    void RemoveDynamiqueGeomAndMatFrame();
    void showConstructVoxelizedCommandsFrame();
    void showConstructSolidAndVolumeCommandsFrame();
    void showConstructSTLCommandsFrame();
    void showConstructVolumeGDMLTEXTCPPCommandsFrame();

    QVector <QString> ElementsNames ;
    QVector <QString> MaterialsNames ;
    QVector <QString> DefinedRegionsNames ;
    QVector <QString> MaterialRegionsNames ;
    QMap <QString,QString> MaterialsNameIDs ;
    QVector <QString> SolidsNames ;
    QVector <QString> VolsNames ;
    QVector <QString> PhyVolsNames ;
    QStringList DefinedParticlesNames;

    QMap <QString,QString> PreDefinedGeomMap ;

    QMap <QString,QString> ElementsSymbolZ ;
    QMap <QString,QString> ElementsSymbolA ;
    QMap <QString,QString> ElementsSymbolSym ;
    QMap <QString,QString> ElementsSymbolName ;
    QStringList PeriodicTableElementsSymbol;
    QStringList NistMaterialsDataBase;

    QMap <int,QString> RegionsSourceList ;
    QMap <QString , QString> RegionBoxDimMap ;
    QMap <int,QString> ParticlesList;
    QMap <int,QString> energiesList ;

    QString PhysicsCutEnergyUnit;
    QString PhysicsCutDistanceUnit;
    QString SourceEnergyUnit;
    QString SourceMomDirAngleUnit;
    QString SourceBoxHalfSizeUnit;
    QString AnalysisRegionEnergyUnit;

    QString MaterialsDataCommandsString;
    QString GeometryCommandsString;
    QString VoxIDsCommandsString;
    QString DICOMCommandsString;
    QString VOXELCommandsString;
    QString CONSCommandsString;
    QString PhysicsCommandsString;
    QString SourceCommandsString;
    QString AnalysisCommandsString;
    QString RunAndScoreCommandsString;
    QString MacrosCommandsString;


    bool isInexec;
    QProcess* process;
    QProcess CoreProcess;

    QByteArray bytes;
    void showResultsOutput(QString , int);

    //QVector <QString> organNamesToChooseArray = {"Head","Brain","Thyroid","SkullCheckBox","UpperSpine","Trunk","Arm","ArmBone","Heart","Lung","Kidney","Liver","Pancreas","Spleen","UrinaryBladder","StomachCheckBox"+"Pelvis","MiddleLowerSpine","LargeIntestine","Uterus","Ovary","Breast","MaleGenitalia","Legs","Leg","LegBone"};

    QString BashCommandsForExecuting;

    QRegExp space ;

    QString download_Directory_Path ;

    QString autoMacrosFileName;
    QString autoMacrosFilePath;
    QString autoSimuFileName ;

    QString PositionsDataFilePath ;
    QString EnergiesDataFilePath ;
    QString MomDirsDataFilePath ;

    // component data Strings.this are just to save the temporary inputs of user, for run or execution, the values are taken directly from components to the FileWriterReader object

    // World
    QString GeometryData_CreateWorld ;
    QString Geometry_setWorldHalfSize ;
    QString Geometry_setWorldHalfSizeunit;
    QString Geometry_setWorldMaterialName ;

    // Geometry
    QString Geometry_setGeometrySymbole ;
    QString Geometry_CreateVolume_GeometryFileType ;
    QString Geometry_CreateVolume_GeometryPath ;

    QString SourceData_setParticleName ;

    // Physics
    QString PhysicsData_setPhysicsData;

    QString PhysicsData_setPhysicsName ;
    QString PhysicsData_setPhotoElectricEffectModel ;
    QString PhysicsData_setComptonScatteringModel ;
    QString PhysicsData_setGammaConversionModel ;
    QString PhysicsData_setRayleighScatteringModel ;
    QString PhysicsData_setElectronIonisationModel ;
    QString PhysicsData_setElectronBremModel ;
    QString PhysicsData_setHadronIonisationModel ;

    QString PhysicsData_setCutsEnergy ;
    QString PhysicsData_setCutsDistance ;

    QString LastRunOutputFile ;

    // Source

    QString SourceData_setSourcePosData ;
    QString SourceData_setSourceEneData ;
    QString SourceData_setSourceMomDirData ;
    QString SourceData_UseDataFiles ;

    QString SourceData_setSourceType ;
    QString SourceData_setSourceData ;
    QString SourceData_setEnergyDistribution;
    QString SourceData_setEnergyData ;
    QString SourceData_setAngleDistribution ;
    QString SourceData_setMomDirData ;
    QString SourceData_setSourceSizeUnit ;
    QString SourceData_setSourceEnergyUnit ;
    QString SourceData_setSourceAngleUnit ;

    QString SourceData_setEventsNumForDataGen ;

    QString SourceData_GeneratePositions ;
    QString SourceData_GenerateEnergies ;
    QString SourceData_GenerateMomDirs ;
    /*
    QString SourceData_setSourceName ;
    QString SourceData_setBoxHalfSize ;
    QString SourceData_setSourcePosition;
    QString SourceData_setDirectionTheta ;
    QString SourceData_setDirectionPhi ;
    QString SourceData_setGaussMean ;
    QString SourceData_setRayleighEmax ;
    QString SourceData_setGaussEData = "";
    QString SourceData_setUniformEData = "";
    QString SourceData_setUniformEmax ;
    QString SourceData_setMonoEnergy ;
*/
    // Score

    int macrosfileinc;
    QString Score_setVariableToScore ;
    QString Score_setVolumesToScore ;
    //QString Score_setAccuracyCalculationLevel ;
    QString Score_setSimNumOnRanksLineEdit;
    QString Score_setRadioNucleidDataLineEdit;
    QString Score_setRadioNucleidBiokineticsLineEdit;
    //QString Score_setRadiationFactors;
    QString Score_setTissueFactors;
    QString Score_setQuantitiesUnits;

    // Analysis

    QString Analysis_GenerateSelfCrossGraphs ;
    QString Analysis_setGraphsData ;
    QString Analysis_setCompareType ;
    QString Analysis_setGraphsExt ;
    QString Analysis_setRefFilePath ;
    QString Analysis_setRefName ;

    QString Analysis_setRegionVariableName ;

    QString  Analysis_GenerateRelativeSDevGraph  ;
    QString  Analysis_GenerateRelativeErrGraph  ;
    QString  Analysis_GenerateRegionsVariableGraph  ;
    QString  Analysis_GenerateEventsDataHisto  ;
    QString  Analysis_setSliceFor2DGraph  ;
    QString  Analysis_setBeamAxis  ;
    QString  Analysis_setSliceID  ;

    // //////////////////////////////////

    QString Analysis_UseLogE ;
    QString Analysis_UseLogVariable ;
    QString Analysis_UseGridXY ;
    QString Analysis_PrintTitle ;
    QString Analysis_LegendPos ;
    QString Analysis_LegendWidth ;
    QString Analysis_AddErrorBar ;

    QString PhysicsData_ParticleForCrossSection;
    QString PhysicsData_EUnitForCrossSection;
    QString PhysicsData_EnergiesForCrossSection;

    std::vector<double> PhysicsData_EnergiesForCrossSectionValues;

    // //////////////////////////////////

    QString Execution_setNumberOfRanksOrThreads ;
    QString DEFAULT_INPUTS;

    QString ExecutionCommand;

    QGridLayout* framLay;
    QLayoutItem* LayItem;
    QWidget* HelpGUIMatVoxWDG;

    QString ConstructDoseCalcsJobName();
    void ShowTerminal(QString);
    QString CreateMaterialAndGeometryDataFromMacrosFile(QString);
    void RunForMultiGeomeries();


    void setBiokineticsDefaulsInputs();
    void CalculateQuantitiesBasedOnICRPData();
    void GenerateSAFFromNewSource(QString);
    void GenerateSAFFromNewTarget(QString);
    void GenerateRadiotracerQuantitiesByInterpolationInDefaultUnit(QString, double);
    double GenerateRadiotracerQuantitiesByInterpolationInDefaultUnitForBiokinetic(QString, double, QString, QString);
    double  GenerateRadiationFactor(QString, double);
    void GenerateDataInTableView();
    void ReadLoadICRPSpectrumData();
    void Read_ICRP107_108Files(QString);
    void Read_ICRP107SpectrumRadiationFiles(QString);
    void CombineSpectrumWithMonoData();
    void Read_DoseCalcs_file(QString);
    void Read_Reference_file(QString);
    void ReadMassesAndWTFactor(QString);
    void ConstructOtherTissuesSource();
    void CreateNewSourceRegion(QString);

    double  GenerateTissueFactor(QString);
    QMap<QString,QMap<QString,QMap<QString,QMap<double,QMap<QString,QMap<QString,double>>>>>> ICRPSAFs ;
    QMap<QString,QMap<QString,QMap<double,double>>> ICRPRadioNuclideData ; // RadioNuclide Particle MonoOrSpectrum MonoEnergyOrSpectumEnergy YieldForMonoEnergyOrSpectumEnergy
    QMap<QString,QMap<QString,QMap<double,double>>> RadioParticleEnergyYieldForSpectrum ; // RadioNuclide Particle MonoOrSpectrum MonoEnergyOrSpectumEnergy YieldForMonoEnergyOrSpectumEnergy
    QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>> ICRPRadioNuclideDataDiscSpec;
    QMap<QString,double> ICRPRadioNuclideHalfLives ;
    QMap<QString,QMap<QString,QMap<QString,double>>> RegionParameterValueMap;
    QMap<QString,QMap<QString,QMap<QString,double>>> ICRPSourceMassMap;
    QMap<QString,QVector<QString>> PhantomSourcesMap;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>>> QuatititesRadioNuclidesParticlesCalculatedData; // Quantity, Geometry, Radionuclide, source, target, value
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> QuatititesRadioNuclidesParticlesCalculatedDataInOrgan; // Quantity, Geometry, Radionuclide, organ, value
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> QuatititesRadioNuclidesCalculatedData; // Quantity, Geometry, Radionuclide, source, target, value
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>> QuatititesRadioNuclidesCalculatedDataInOrgan; // Quantity, Geometry, Radionuclide, organ, value
    QVector<QPair<double,QString>> RadionuclideValuePairVectorAlphabetical; // Radionuclide, QuantityValue
    QVector<QPair<double,QString>> RadionuclideValuePairVectorValuesAscending; // Radionuclide, QuantityValue
    QVector<QPair<double,QString>> RadionuclideValuePairVectorValuesDescending; // Radionuclide, QuantityValue
    std::map<QString,std::map<QString, std::vector<double>>> SourceParticleEnergyValues ;
    std::map<QString,double> TissueFactorMap ;
    //std::map<QString,std::map<double,double>> RadiationFactorMap ;
    QMap<QString,QVector<QString>> RadionuclidesParticles ;
    QMap<QString,QMap<QString,QMap<QString,double>>> RadioTracerSourceOrganResidenceTime;
    std::map<QString,QString> RadiotracerradionucleidMap;
    //bool IsICRPCalculationTerminated;
    QMap<QString,QMap<QString,QMap<QString,double>>> SpecialOrganSoftMassFraction;

private slots:

    //void insertCompletion();
    void removeHugFiles_slot();

    void btnCheckALL();
    void btnUNCheckALL();

    void CheckJobOutputsFromStopSlot();
    void btnSaveElem_slot();
    void btnAddElem_slot();
    void btnAddMat_slot();
    void btnAddNistMat_slot();
    void CheckBoxFracNum_slot();
    void on_comboBoxRegionToVisualize_textActivated(QString);
    void on_LimitsPlan_textActivated(QString);
    void btnAddVoxelsRegionData_slot();
    void UseMaterialsAsRegionNames_slot();
    void btnAddVoxelsData_slot();
    void btnAddDcmMatOrdered_slot();
    void btnAddDcmCTNumDen_slot();
    void ComboxSolidType_slot();
    void btnAddDcmTypeDir_slot();
    void btnChooseDcmDir_slot();
    void btnAddVol_slot();
    void openqsubMacros_slot();

    void btnAddVisualizationLimits_slot();

    void on_BtnPhantomReturn_clicked();

    void btnChooseADir1_slot();
    void btnChooseADir2_slot();
    void btnAddSol_slot();
    void BtnAddSolPath_slot();
    void BtnAddVolPath_slot();
    void BtnEditVolFile_slot();
    void on_BtnsourceReturn_clicked();
    void on_btnanalyseReturn_clicked();
    void on_AnalyseBtnReset_clicked();
    void on_SourceBtnReset_clicked();
    void on_PhantomBtnReset_clicked();
    void on_actionOpen_triggered();
    void on_actionSave_triggered();
    void on_actionRun_2_triggered();
    void processOutput();
    void on_actionInstallations_triggered();
    void on_actionClose_triggered();
    void on_actionAnalysis_triggered();
    void on_actionStop_triggered();
    void on_actionVisualization_triggered();
    void on_SourceComboBoxPhysUsed_currentIndexChanged(const QString &arg1);
    void on_comboBoxTypeOfSources_currentIndexChanged(const QString &arg1);
    void on_SourceComboBoxEnergyDist_currentIndexChanged(const QString &arg1);
    void on_radioButtonVoxel_clicked(bool checked);
    void on_radioButtonGDML_clicked(bool checked);
    void on_radioButtonConstruct_clicked(bool checked);
    void on_AnalysisComboBoxGraphData_currentIndexChanged(const QString &arg1);
    void on_radioButtonDICOM_clicked(bool checked);
    void on_SourceComboBoxAngleDist_currentIndexChanged(const QString &arg1);
    void on_actionClear_Output_triggered();
    void on_EditGeometyDataBtn_clicked();
    void on_MaterialsDataEditButton_clicked();
    void on_actionReturn_triggered();
    void on_actionClear_Inputs_triggered();
    void on_BuildButton_clicked();
    void on_groupBoxRootGraphs_clicked();
    void on_btnAnalysiInputFileFile_clicked();
    void on_btnPositionsFile_clicked();
    void on_btnEnergiessFile_clicked();
    void on_btnMomDirsFile_clicked();
    void on_AnalysisbtnGraphsDirPath_clicked();
    void on_radioButtonVoxIDs_clicked(bool checked);
    void on_radioButtonTEXT_clicked(bool checked);
    void on_AnalysisbtnGenerate_clicked();
    void on_RootAnalyseBtnReset_clicked();
    void on_RootAnalyseBtnReturn_clicked();
    void on_RunButton_clicked();
    void on_MaterialsDataShowButton_clicked();
    void on_GeometryDataShowButton_clicked();
    void on_checkBoxWorldConst_clicked(bool checked);
    void on_openResultsDirButton_clicked();
    void on_RootAnalysispushButtonOpenROOTGUI_clicked();

    void on_pushButtonMerge_clicked();

    void on_pushButtonOpenResultsFile_clicked();

    void on_actionGeometryModelling_triggered();

    void on_pushButtonQstat_clicked();

    void on_pushButtonLoadExe_clicked();

    void on_pushButtonGenerateExe_clicked();

    void on_checkBoxRocks_clicked(bool checked);

    void on_ViewGeometyDataBtn_clicked();

    void on_WorldDataShowButtonpushButton_clicked();

    void on_pushButtonEditGeometryFile_clicked();

    void on_pushButtonEditGeomFile_clicked();

    void on_pushButtonEditMacros_clicked();

    void on_checkBoxUseMacroCommandFile_stateChanged(int arg1);

    void on_actionEditFile_triggered();

    void on_actionUseTextInput_triggered();

    void on_radioButtonCpp_clicked(bool checked);

    void on_radioButtonSTL_clicked(bool checked);

    void on_pushButtonChooseVoxIDsFile_clicked();

    void on_actionAbout_triggered();

    void on_pushButtonCloseFile_clicked();

    void on_pushButton_clicked();

    void on_pushButtonChecQsubMPIsim_clicked();

    void on_pushButtonStopJob_clicked();

    void on_pushButtonCalculateNumSimulation_clicked();

    void on_actionHow_To_Use_triggered();

    void on_RESUBButton_2_clicked();

    void on_pushButtonChoosWorldFile_clicked();

    void on_actionRead_File_triggered();

    void on_actionCheck_Inputs_triggered();

    void on_pushButtonAddGeometryMacros_clicked();

    void on_checkBoxSimulateGeometriesWitOneSource_clicked(bool checked);

    void on_pushButtonCheckMatGeoCommandsOfFiles_clicked();

    void on_pushButtonReadICRPData_clicked();

    void on_pushButtonGetResultsNucl_clicked();

    void on_pushButtonReverseData_clicked();

    void on_comboBoxQuantityNucl_currentTextChanged(const QString &arg1);

    void on_pushButtonSaveTableInPDF_clicked();

    void on_pushButtonReadUSERData_clicked();

    void on_pushButtonAddBiokineticModel_clicked();

    void on_pushButtonSaveBiokinetikModelData_clicked();

    void on_pushButtonShowBiokineticData_clicked();

    void on_pushButtonGenerateQuantitiesWithBiokineticData_clicked();

    void on_actionPackages_Statut_triggered();

    void on_pushButton_2RootEnergyView_clicked();

    void on_pushButton_2LoadICRPSpectrum_clicked();

    void on_comboBoxPreDefinedGeom_currentTextChanged(const QString &arg1);

    void on_comboBoxSourceRegionProposed_textActivated(const QString &arg1);

    void on_actionSave_Inputs_To_Default_File_triggered();

    void on_actionSend_Results_triggered();

    void on_checkBoxUsePreDefinedGeom_stateChanged(int arg1);

    void on_pushButtonChangeMasses_clicked();

    void on_pushButtonSaveMasses_clicked();

    void on_actionUpdate_triggered();

    void on_actionRestart_triggered();

    void on_radioButtonTET_clicked();

    void on_comboBoxSources_currentTextChanged(const QString &arg1);

    void on_comboBoxTargets_currentTextChanged(const QString &arg1);

    void on_comboBoxSourceOrTargetsForBiokinetics_currentTextChanged(const QString &arg1);

    void on_pushButtonShowOutputs_clicked();

    void on_checkBoxnohup_clicked();

    void on_pushButtonShowOutputsAndMacros_clicked();

    void on_pushButtonStopCurrentProcess_clicked();

    void on_pushButtonChooseResultsDir_clicked();

    void on_comboBoxPhantom_currentTextChanged(const QString &arg1);

    void on_pushButtonReadICRP107128_clicked();

    void on_pushButton_selectSource_clicked();

    void on_pushButton_selectTarget_clicked();

    void on_pushButton_TissueWeithingFactor_clicked();

    void on_pushButton_SaveWT_clicked();

    void on_pushButton_ReadDoseCalcsSAFs_clicked();

    void on_checkBoxInterpolationType_stateChanged(int arg1);

    void on_pushButton_SaveADEDResultsInResultsData_clicked();

    void on_checkBoxHalflivesNucl_clicked(bool checked);

    void on_checkBoxEffDoseLimNucl_clicked(bool checked);

    void on_pushButtonAddRowToTable_clicked();

    void on_comboBoxPhantom_textActivated(const QString &arg1);

    void on_comboBoxRadioPharmaceutiques_textActivated(const QString &arg1);

private:

    QString OpenSourceRegionsDialogFor;
    QString SourceRegionNameInNewRegiondialog;
    QString OpenTargetRegionsDialogFor;
    QStringList CurrentSources;
    QVector<double> CurrentSourcesFractions;
    QMap<QString, QMap<QString, double>> RegionFraction;
    QStringList CurrentTargets;
    QVector<double> CurrentTargetsFractions;

    QMap <QString, QStringList> NewSourcesSources;
    QMap <QString, QStringList> NewTargetsTargets;

    QString Path1;
    QString Path2;

    QLineEdit* FromPath;
    QLineEdit* ToPath;

    bool IsICRPFilesAreRead ;
    QStringList MacrosFiles;
    QString MacrosFilePathReadForMultiGeometryAndOneSource ;

    Highlighter* highlighter ;
    Highlighter* highlighter1 ;
    Highlighter* highlighter2 ;

    Ui::MainWindow *ui;
    OrgansChooserDialog *OrgChoDialog;

    int EditFlag ;
    InstallationDialog *installationDialogObj;
    geometriesvisualization *geometriesvisualizationObj;

    filesManager* fileManagerObject;
    //mainGeant4AppCore* mainGeant4AppCoreObject;

    PlotDialog* plotDialogObject;

    visualizationManager* visualizationManagerObject;
    GeometryModellingEditor* GeometryModellingEditorObject;
};

#endif // MAINWINDOW_H
