#ifndef PLOTDIALOG_H
#define PLOTDIALOG_H

#include <QDialog>
#include <QMainWindow>
#include <QInputDialog>

#include "gui/qcustomplot.h"
#include "gui/graphsDataGeneration.h"
#include "gui/filesManager.h"
#include <QCloseEvent>
#include <QKeyEvent>

namespace Ui {
class PlotDialog;
}

class PlotDialog : public QDialog
{
    Q_OBJECT

public:

    void keyPressEvent(QKeyEvent*);
    explicit PlotDialog(QWidget *parent = 0);
    ~PlotDialog();
    void InitializeCustomPlotForGraphs();
    void Data_initialization();
    void ConnectCustomPlotSIGNALS_and_SLOTS();
    void create_graphs(QMap<double, double>,QString, QString, int);
    void create_RatioData_For_ParticleAndRadiotracer_In_RatioPlot();
    void create_QuantitiesData_For_Radiotracer_In_Plot_Bars();
    void create_QuantitiesData_For_Radiotracer_In_Plot_Graphs();
    void InitializeCustomPlotForStyledPlot();

    QString RadioTraceBarPlotData;

    int NumbOfScoVar;
    int NumbOfSrcReg;
    int NumbOfTrgReg;
    int NumbOfParNme;

    QString EnergyDistriburionFilePath;
    QString ReferenceFilePath;
    QString ReferenceFilePath_2;
    QString ResultFilePath;
    QString AnalysisInputFilePath;
    QString CrossSectionFilePath;

    QString Compare_Type;
    QString GraphsData;
    QString QuantitiesToScore;
    double GraphForEnergy;
    QString CompareReferenceName;
    QString CompareReferenceName_2;
    QString ParticleName;
    QString GeometrySymbol;
    QString SourceOrgan;
    QString TargetOrgan;
    QString dataOf; // MIRD or USER or...
    QString FileExtToSave;

    QString graphs_Title;
    bool UseErrorBar;

    int number_of_energies;
    int number_of_energies_in_compareFile ;

    int Reference_or_Result;
    int Self_Cross;

    QMap<double,double> ErrorMap;

    QString tblFrom;
    QString GraphName;

    //int index;

    QCPTextElement *title;

    QString GraphTypeShown;

    QVector<QString> Organs_to_score;

    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> ReferenceTable;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> ReferenceTable_2;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> ResultTable;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> ErrorTable;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> StandartDeviation;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>> ResRefErrCompTables;
    QMap<QString,QMap<QString,QMap<QString,double>>> RegionParameterValueMap;
    QMap<QString,QString> AnalysisInputMap;
    QMap<QString,QMap<QString,double>> SourceTargetDefinedEnergyMap;

    QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>> ParticleMaterialProcessEnergySigmaMap;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>> ResultParticleSourceEnergyTime;

    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> ReferenceQuantityGeometryRadioTracerSourceTargetValues ;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>> ResultQuantityRadioTracerSourceTargetValues ;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> ResultQuantityGeometryRadioTracerSourceTargetValues ;

    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,int>>>>> QuantityGeometryParticleSourceTargetsNames ;

    QString RadioTracerName;
    QString RadioTracerSourceName;

    QVector<double> xticks;
    QVector<QString> xlabels;

    QString ParticleForCrossSection;
    QString MaterialForCrossSection;
    QString ProcessForCrossSection;

    QString RegionVariableNameWithUnit;
    QString RegionVariableName;

    QString ScoreVariableUnit;

    QString DiffExp = "";
    QString DiffSym = "";

    void setComboboxInputs();

    void closeEvent(QCloseEvent *event);  // show prompt when user wants to close app

    void GetPlotInputDataAndReadFiles();
    void create_Data_for_Cross_and_Self_graphs(QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> , QString);
    void create_Data_for_Cross_and_Self_ForOneEnergy_graphs();
    void create_QuantitiesData_Cross_Self_FromRadioTracerSource_graphs();
    void create_QuantitiesData_Cross_Self_FromRadioTracerIntake_graphs();
    void create_Data_for_RelativeErr_graphs();
    void create_Data_for_MCError_graphs();
    void create_Data_for_RegionParametre_graphs();
    void create_Data_for_CrossSection_graphs();
    void create_Data_for_SourceTimeSimulation_graphs(QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>> , QString);

    double RelativeDifferenceCalculation(double, double);

    void TableViewFillingWithSelectedGraph();
    void TableViewFillingWithAllGraphs();
    void TableWidgetFillingWithAllGraphs();
    void TableWidgetFillingWithSelectedGraph();

    void SaveTableFromTableWidgetToFile();
    void SaveTableFromTableViewToFile();
    void GenerateTableInTableView();

    void getEnergyDistributionFromFile(QString);
    void OpenDialogOfCombinationChooser();


    //void GenerateLatexTableResultReferenceForOneEnergy();
    //void GenerateLatexTableResultReferenceRadioTracerSValues();

    QVector <QString> readFilesAndGetFilesNamesVector();

    void showResultsOutput(QString text, int level);

    /*
    void CalculateQuantitiesBasedOnICRPData();
    void GenerateRadiotracerQuantitiesByInterpolationInDefaultUnit(QString, double);
    double  GenerateRadiationFactor(QString, double);
    void GenerateDataInTableView();
    double  GenerateTissueFactor(QString);
    QMap<QString,double> RadioTracerInjectedActivity;
    QMap<QString,QMap<QString,QMap<QString,QMap<double,QMap<QString,QMap<QString,double>>>>>> ICRPSAFs ;
    QMap<QString,QMap<QString,QMap<double,double>>> ICRPRadioNuclideData ; // RadioNuclide Particle MonoOrSpectrum MonoEnergyOrSpectumEnergy YieldForMonoEnergyOrSpectumEnergy
    QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>> ICRPRadioNuclideDataDiscSpec;
    QMap<QString,double> ICRPRadioNuclideHalfLives ;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>>> QuatititesRadioNuclidesParticlesCalculatedData; // Quantity, Geometry, Radionuclide, source, target, value
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> QuatititesRadioNuclidesParticlesCalculatedDataInOrgan; // Quantity, Geometry, Radionuclide, organ, value
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> QuatititesRadioNuclidesCalculatedData; // Quantity, Geometry, Radionuclide, source, target, value
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>> QuatititesRadioNuclidesCalculatedDataInOrgan; // Quantity, Geometry, Radionuclide, organ, value
    QMap<QString,QVector<double>> ValuesOrderedData; // Radionuclide, QuantityValue
    QMap<double,QString> ValueRadioNuclideOrderedData; // QuantityValue, Radionuclide
*/


    std::map<QString,std::map<QString, std::vector<double>>> SourceParticleEnergyValues ;
    std::map<QString,double> TissueFactorMap ;
    std::map<QString,std::map<double,double>> RadiationFactorMap ;

    double DefaultErrorDistance ;
    double MinValForLog ;

    QVector<QCPScatterStyle::ScatterShape> shapes;
    QStringList ShapesName;
    QMap<QString,QCPScatterStyle::ScatterShape> ShapesMap;
    QStringList lineStylesName;
    QMap<QString, QCPGraph::LineStyle> lineStylesMap;

    QMap<QString,QString> QuantityUnit;

    QString WhichComponent;

    QColor TextColor;
    QColor LineColor;
    QColor ChosenBackgroundColor;
    QPushButton * TextColorBtn;
    QPushButton * LineColorBtn;
    QPushButton * BackgroundColorBtn;
    QVector<double> SelectedGraphEnergy;
    QVector<double> SelectedGraphScoreVar;


    QVector<QString> SourceOrganForCombination;
    QVector<QString> TargetOrganForCombination;

    QString XNumberFormat;
    QString YNumberFormat;

    QComboBox * TargetListForCombinations;
    QComboBox * SourceListForCombinations;
    QTextEdit * CombinationsText;
    QRadioButton* GraphsRadio;

    filesManager* fileManagerObjectPlot;
private slots:

    void LineColorMapSlot();
    void TextColorMapSlot();
    void BackgroundColorMapSlot();
    void AddCombinationBtnSlot();
    void AddAllCombinationBtnSlot();

    void titleDoubleClick(QMouseEvent *event);
    void editPlotParameters();
    void editGraphsStyles();
    void axisLabelDoubleClick(QCPAxis* axis, QCPAxis::SelectablePart part);
    void legendDoubleClick(QCPLegend* legend, QCPAbstractLegendItem* item);
    void selectionChanged();
    void mousePress();
    void mouseWheel();
    void addRandomGraph();
    void removeSelectedGraph();
    void removeAllGraphs();
    void contextMenuRequest(QPoint pos);
    void moveLegend();
    void graphClicked(QCPAbstractPlottable *plottable, int dataIndex);
    void ChangeXaxisScale();
    void ChangeYaxisScale();
    void ChangeXNumberFormat();
    void ChangeYNumberFormat();

    void on_PlotcomboBoxGraphData_currentIndexChanged(const QString &arg1);
    void on_plotcomboBoxCompareType_currentIndexChanged(const QString &arg1);
    void ActionWhenChoosingCompareType();
    void ActionWhenChoosingGraphData();

    void on_PlotButtonChoseReferenceFileDialog_clicked();
    void on_PlotButtonChoseResultFileDialog_clicked();
    void on_plotButtonSaveGraph_clicked();
    void on_plotButtonGenerateGraphs_clicked();

    void on_plotButtonSaveTableData_clicked();

    void on_PlotButtonReadResultFile_clicked();
    void on_PlotButtonReadReferenceFile_clicked();
    void on_plotButtonSaveTableOfAllGraphsData_clicked();
    void on_plotButtonGenerateTableOfAllGraphsData_clicked();
    void on_plotButtonGenerateRelErrGraphs_clicked();
    void on_plotButtonGenerateMCSimErrorGraphs_clicked();
    void on_plotButtonGenerateRegionVarGraphs_clicked();
    void on_PlotButtonReadSimulationInpFile_clicked();
    void on_PlotButtonChoseSimulationInpFileDialog_clicked();
    void on_plotButtonGenerateCrossSectionGraphGraphs_2_clicked();
    void on_PlotButtonReadCrossSectionFile_clicked();
    void on_PlotButtonChoseCrossSectionFileDialog_2_clicked();
    void on_pushButtonOpenSaveGraphsDir_clicked();
    void on_PlotButtonChoseReferenceFileDialog_2_clicked();
    void on_PlotButtonReadReferenceFile_2_clicked();
    void on_plotButtonGenerateSimulationTimeGraph_clicked();
    void on_pushButtonGenerateNewGraph_clicked();
    void on_pushButtonAddGraph_clicked();
    void on_pushButtonSaveGraphData_clicked();
    void on_pushButtonFillTable_clicked();
    //void on_pushButtonGenerateRegionDataTable_clicked();
    //void on_pushButtonGenerateSelfCrossTable_clicked();
    //void on_pushButtonSelfCrossLatex_clicked();
    //void on_pushButtonRegionLatex_clicked();
    //void on_pushButtonCrossSectionLatex_clicked();
    //void on_pushButtonGenerateCSV_clicked();
    //void on_pushButtonOpenCSV_clicked();

    void on_pushButtonGenerateGraphsForEnergy_clicked();

    void on_PlotButtonReadAllData_clicked();

    void on_PlotButtonEditSimData_clicked();

    void on_PlotButtonEditResultsData_clicked();

    void on_pushButtonRadioTracerBarsAndRatiosPlots_clicked();

    void on_plotcomboBoxScoredVariable_activated(const QString &arg1);

    void on_comboBoxGeometrySymbole_activated(const QString &arg1);

    void on_plotcomboBoxParticle_activated(const QString &arg1);

    void on_comboBoxSourceOrgan_activated(const QString &arg1);

    void on_checkBoxDiffForRadiotracerOr_stateChanged(int arg1);


    void on_pushButtonOpenInRoot_clicked();
/*
    void on_pushButtonReadICRPData_clicked();
    void on_comboBoxQuantityNucl_currentTextChanged(const QString &arg1);
    void on_pushButtonGetResultsNucl_clicked();
    void on_pushButtonReverseData_clicked();
*/
private:

    QStringList Geometrylist;
    QStringList SourceOrganlist;
    QStringList TargetOrganlist;
    QStringList Particlelist;
    QStringList ScoreVarlist;
    //QStringList RadioTracerlist;

    QMap<QString,QMap<QString,QMap<QString,QVector<double>>>> ResEnergies ;
    QMap<QString,QMap<QString,QVector<QString>>> OrganNamesToScoreForQuantityAndParticle ;
    QVector<QString> OrganNamesToScore;

    QLineEdit* header_editor;
    int editor_index;
    bool eventFilter(QObject*, QEvent*);

    GraphsDataGeneration* graphsDataGenerationObj;

    Ui::PlotDialog *ui;
    QStandardItemModel *model;
};


#endif // PLOTDIALOG_H
