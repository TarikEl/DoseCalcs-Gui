#include "gui/plotDialog.h"
#include "gui/ui_plotDialog.h"

#include "QTextStream"
#include <sstream>
#include <fstream>
#include "QTextTable"
#include "gui/terminal.h"


extern QString DoseCalcsCore_source_dir_path;
extern QString DoseCalcsCore_build_dir_path;
extern QString Root_Lib_dir_path;
extern QString DoseCalcsExecutingFileName;
extern QString UserCurrentResultsDirPath;
extern QString GUIPackagesAndFilesDirPath;

extern QString ResultDirectoryName;
extern QString ScriptDirectoryName;
extern QString DoseCalcsExecutableName;
extern QString MergeExecutableName;
extern QString GraphExecutableName;
extern QString MacroFileName;
extern QString MacroFilePath;
extern QString ResultFileName;
extern QString ReferenceFileName;
extern QString GraphsOutDirName;
extern QString BuildingFileName;

extern QVector <QString> RunAndScoreCommands ;
extern QVector <QString> AnalysisCommands ;

/*
extern QString ICRPDATAPath;
extern QStringList DoseCalcsQuantities ;
extern QMap <QString,QStringList> QuantitiesUnitsLists ;
extern QMap <QString,QMap <QString,double>> QuantitiesConversionFromDefault ;
*/

PlotDialog::~PlotDialog()
{
    delete ui;
}

PlotDialog::PlotDialog(QWidget *parent) : QDialog(parent), ui(new Ui::PlotDialog) {
    
    QTextStream(stdout) << "---------------- PlotDialog ------------------" << "\n";
    
    showResultsOutput("DoseCalcsCore_source_dir_path is " + DoseCalcsCore_source_dir_path , 4);
    showResultsOutput("DoseCalcsCore_build_dir_path is " + DoseCalcsCore_build_dir_path , 4);
    
    srand(QDateTime::currentDateTime().toTime_t());
    ui->setupUi(this);
    
    ui->PlotDataTableView->setVisible(false); //or true - later in the code
    ui->plotcomboxTableFileExtToSave_2->setVisible(false);
    ui->plotcomboxTableFileExtToSave->setVisible(false);
    ui->radioButtonLRD->setChecked(true);
    Data_initialization();

    InitializeCustomPlotForGraphs();
    ConnectCustomPlotSIGNALS_and_SLOTS();
}

// called rom PlotDialog constructor()
void PlotDialog::Data_initialization(){
    
    MinValForLog = 1e+12;
    
    fileManagerObjectPlot = new filesManager;
    


    QStringList GraphDatalist=(QStringList()<<"Result"<<"Reference_Result");
    ui->PlotcomboBoxGraphData->addItems(GraphDatalist);
    
    QStringList CompareTypelist=(QStringList()<<"Self"<<"Cross");
    ui->plotcomboBoxCompareType->addItems(CompareTypelist);
    
    QStringList SourceOrganlist=(QStringList()<<"Spleen"<<"Kidney"<<"Liver"<<"Pancreas"<<"Lung");
    //ui->comboBoxSourceOrgan->addItems(SourceOrganlist);
    
    QStringList TargetOrganlist=(QStringList()<<"Spleen"<<"Kidney"<<"Liver"<<"Pancreas"<<"Lung");
    //ui->PlotcomboBoxTargetOrgan->addItems(TargetOrganlist);
    
    QStringList Particlelist=(QStringList()<<"gamma"<<"e-");
    //ui->plotcomboBoxParticle->addItems(Particlelist);
    
    QStringList ScoredVariablelist=(QStringList()<<"SAF"<<"AF"<<"AD"<<"AE"<<"S"<<"E"<<"H"<<"DR"<<"ER");
    ui->plotcomboBoxScoredVariable->addItems(ScoredVariablelist);
    
    QStringList FileExtlist=(QStringList()<<"Pdf"<<"Jpg"<<"Png"<<"Bmp");
    ui->plotcomboxFileExtToSave->addItems(FileExtlist);
    //ui->plotcomboxFileExtToSave->setVisible(false);
    
    QStringList TableFileExtlist=(QStringList()<< "" << "Pdf"<<"Jpg"<<"Png"<<"Bmp");
    ui->plotcomboxTableFileExtToSave->addItems(TableFileExtlist);
    //ui->plotcomboxTableFileExtToSave->setVisible(false);
    
    QStringList RegionVars=(QStringList() << "Mass"<<"Density"<<"Volume"<<"Distance");
    ui->PlotComboboxRegionVariable->addItems(RegionVars);
    

    /*
    QStringList ICRPPhantoms=(QStringList()<<"ICRPAdultMale"<<"ICRPAdultFemale");
    ui->comboBoxPhantom->addItems(ICRPPhantoms);

    QStringList SourceOrTargetOrCombination=(QStringList()<<"Sources"<<"Targets"<<"Combinations"<<"Phantom");
    ui->comboBoxTotalorSourceOrCombNucl->addItems(SourceOrTargetOrCombination);

    ui->comboBoxPeriodUnit->addItems(QuantitiesUnitsLists["T"]);
    ui->comboBoxQuantityNucl->addItems(DoseCalcsQuantities);
    ui->comboBoxEffDoseUnit->clear();ui->comboBoxEffDoseUnit->addItems(QuantitiesUnitsLists["AE"]);

    ui->RadPhotonCheckBox->setChecked(true);
    ui->RadelectronCheckBox->setChecked(true);
    ui->RadPositronCheckBox->setChecked(true);
    ui->RadAlphaCheckBox->setChecked(true);
    ui->RadNeutronCheckBox->setChecked(true);
*/
    //ui->PlotLineEditResultsDir->setText(UserCurrentResultsDirPath);
    ui->PlotLineEditReferenceFile->setText(UserCurrentResultsDirPath+"/"+ReferenceFileName);
    ui->PlotLineEditReferenceFile_2->setText(UserCurrentResultsDirPath+"/"+ReferenceFileName);
    ui->PlotLineEditResultFile->setText(UserCurrentResultsDirPath+"/"+ResultFileName);
    ui->PlotLineEditSimulationInpFile->setText(MacroFilePath);
    ui->PlotLineEditCrossSectionFile->setText(UserCurrentResultsDirPath+"/CrossSectionData");

    ActionWhenChoosingCompareType();
    ActionWhenChoosingGraphData();
}

// called rom PlotDialog constructor()
void PlotDialog::InitializeCustomPlotForGraphs(){
    
    QTextStream(stdout) << "---------------- InitializeCustomPlotForGraphs() ------------------" << "\n";
    
    lineStylesName=(QStringList()
                    <<""
                    <<"lsImpulse"
                    <<"lsLine"
                    <<"lsNone"
                    <<"lsStepCenter"
                    <<"lsStepLeft"
                    <<"lsStepRight"
                    );

    lineStylesMap["lsImpulse"] = QCPGraph::lsImpulse ;
    lineStylesMap["lsLine"] = QCPGraph::lsLine ;
    lineStylesMap["lsNone"] = QCPGraph::lsNone ;
    lineStylesMap["lsStepCenter"] = QCPGraph::lsStepCenter ;
    lineStylesMap["lsStepLeft"] = QCPGraph::lsStepLeft ;
    lineStylesMap["lsStepRight"] = QCPGraph::lsStepRight ;

    ShapesName=(QStringList()
                <<""
                <<"ssCross"
                <<"ssPlus"
                <<"ssCircle"
                <<"ssDisc"
                <<"ssSquare"
                <<"ssDiamond"
                <<"ssStar"
                <<"ssTriangle"
                <<"ssTriangleInverted"
                <<"ssCrossSquare"
                <<"ssPlusSquare"
                <<"ssCrossCircle"
                <<"ssPlusCircle"
                <<"ssPeace"
                <<"ssCustom"
                );

    ShapesMap["ssCross"] = QCPScatterStyle::ssCross;
    ShapesMap["ssPlus"] = QCPScatterStyle::ssPlus;
    ShapesMap["ssCircle"] = QCPScatterStyle::ssCircle;
    ShapesMap["ssDisc"] = QCPScatterStyle::ssDisc;
    ShapesMap["ssSquare"] = QCPScatterStyle::ssSquare;
    ShapesMap["ssDiamond"] = QCPScatterStyle::ssDiamond;
    ShapesMap["ssStar"] = QCPScatterStyle::ssStar;
    ShapesMap["ssTriangle"] = QCPScatterStyle::ssTriangle;
    ShapesMap["ssTriangleInverted"] = QCPScatterStyle::ssTriangleInverted;
    ShapesMap["ssCrossSquare"] = QCPScatterStyle::ssCrossSquare;
    ShapesMap["ssPlusSquare"] = QCPScatterStyle::ssPlusSquare;
    ShapesMap["ssCrossCircle"] = QCPScatterStyle::ssCrossCircle;
    ShapesMap["ssPlusCircle"] = QCPScatterStyle::ssPlusCircle;
    ShapesMap["ssPeace"] = QCPScatterStyle::ssPeace;
    ShapesMap["ssCustom"] = QCPScatterStyle::ssCustom;

    shapes << QCPScatterStyle::ssCross;
    shapes << QCPScatterStyle::ssPlus;
    shapes << QCPScatterStyle::ssCircle;
    shapes << QCPScatterStyle::ssDisc;
    shapes << QCPScatterStyle::ssSquare;
    shapes << QCPScatterStyle::ssDiamond;
    shapes << QCPScatterStyle::ssStar;
    shapes << QCPScatterStyle::ssTriangle;
    shapes << QCPScatterStyle::ssTriangleInverted;
    shapes << QCPScatterStyle::ssCrossSquare;
    shapes << QCPScatterStyle::ssPlusSquare;
    shapes << QCPScatterStyle::ssCrossCircle;
    shapes << QCPScatterStyle::ssPlusCircle;
    shapes << QCPScatterStyle::ssPeace;
    shapes << QCPScatterStyle::ssCustom;
    
    ui->customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes | QCP::iSelectLegend | QCP::iSelectPlottables);

    //ui->customPlot->xAxis->setRange(-8, 8);
    //ui->customPlot->yAxis->setRange(-5, 5);
    ui->customPlot->axisRect()->setupFullAxesBox();
    
    ui->customPlot->plotLayout()->insertRow(0);
    title = new QCPTextElement(ui->customPlot, "Graphs Title", QFont("sans", 17, QFont::Bold));
    ui->customPlot->plotLayout()->addElement(0, 0, title);
    
    ui->customPlot->xAxis->setLabel("E(MeV)");
    ui->customPlot->yAxis->setLabel("f(E)");
    ui->customPlot->yAxis->grid()->setSubGridVisible(true);
    ui->customPlot->xAxis->grid()->setSubGridVisible(true);
    
    // log for y
    QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
    XNumberFormat = "gb"; YNumberFormat = "gb";

    ui->customPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    ui->customPlot->yAxis->setTicker(logTicker);
    ui->customPlot->yAxis->setNumberFormat(YNumberFormat); // eEfgG  e = exponential, b = beautiful decimal powers
    ui->customPlot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
    
    // log for x
    ui->customPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    ui->customPlot->xAxis->setTicker(logTicker);
    ui->customPlot->xAxis->setNumberFormat(XNumberFormat); // e = exponential, b = beautiful decimal powers
    ui->customPlot->xAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
    
    ui->customPlot->legend->setVisible(true);
    QFont legendFont = font();
    legendFont.setPointSize(10);
    ui->customPlot->legend->setFont(legendFont);
    ui->customPlot->legend->setSelectedFont(legendFont);
    ui->customPlot->legend->setSelectableParts(QCPLegend::spItems); // legend box shall not be selectable, only legend items
    
}
// called rom PlotDialog constructor()
void PlotDialog::ConnectCustomPlotSIGNALS_and_SLOTS(){

    QTextStream(stdout) << "---------------- ConnectCustomPlotSIGNALS_and_SLOTS() ------------------" << "\n";

    ui->customPlot->rescaleAxes();

    // connect slot that ties some axis selections together (especially opposite axes):
    connect(ui->customPlot, SIGNAL(selectionChangedByUser()), this, SLOT(selectionChanged()));

    // connect slots that takes care that when an axis is selected, only that direction can be dragged and zoomed:
    connect(ui->customPlot, SIGNAL(mousePress(QMouseEvent*)), this, SLOT(mousePress()));
    connect(ui->customPlot, SIGNAL(mouseWheel(QWheelEvent*)), this, SLOT(mouseWheel()));

    // make bottom and left axes transfer their ranges to top and right axes:
    connect(ui->customPlot->xAxis, SIGNAL(rangeChanged(QCPRange)), ui->customPlot->xAxis2, SLOT(setRange(QCPRange)));
    connect(ui->customPlot->yAxis, SIGNAL(rangeChanged(QCPRange)), ui->customPlot->yAxis2, SLOT(setRange(QCPRange)));

    // connect some interaction slots:
    connect(ui->customPlot, SIGNAL(axisDoubleClick(QCPAxis*,QCPAxis::SelectablePart,QMouseEvent*)), this, SLOT(axisLabelDoubleClick(QCPAxis*,QCPAxis::SelectablePart)));
    connect(ui->customPlot, SIGNAL(legendDoubleClick(QCPLegend*,QCPAbstractLegendItem*,QMouseEvent*)), this, SLOT(legendDoubleClick(QCPLegend*,QCPAbstractLegendItem*)));
    connect(title, SIGNAL(doubleClicked(QMouseEvent*)), this, SLOT(titleDoubleClick(QMouseEvent*)));

    // connect slot that shows a message in the status bar when a graph is clicked:
    connect(ui->customPlot, SIGNAL(plottableClick(QCPAbstractPlottable*,int,QMouseEvent*)), this, SLOT(graphClicked(QCPAbstractPlottable*,int)));

    // setup policy and connect slot for context menu popup:
    ui->customPlot->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(ui->customPlot, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(contextMenuRequest(QPoint)));
}

void PlotDialog::InitializeCustomPlotForStyledPlot(){
    // prepare data:
    QVector<double> x1(20), y1(20);
    QVector<double> x2(100), y2(100);
    QVector<double> x3(20), y3(20);
    QVector<double> x4(20), y4(20);
    for (int i=0; i<x1.size(); ++i)
    {
        x1[i] = i/(double)(x1.size()-1)*10;
        y1[i] = qCos(x1[i]*0.8+qSin(x1[i]*0.16+1.0))*qSin(x1[i]*0.54)+1.4;
    }
    for (int i=0; i<x2.size(); ++i)
    {
        x2[i] = i/(double)(x2.size()-1)*10;
        y2[i] = qCos(x2[i]*0.85+qSin(x2[i]*0.165+1.1))*qSin(x2[i]*0.50)+1.7;
    }
    for (int i=0; i<x3.size(); ++i)
    {
        x3[i] = i/(double)(x3.size()-1)*10;
        y3[i] = 0.05+3*(0.5+qCos(x3[i]*x3[i]*0.2+2)*0.5)/(double)(x3[i]+0.7)+qrand()/(double)RAND_MAX*0.01;
    }
    for (int i=0; i<x4.size(); ++i)
    {
        x4[i] = x3[i];
        y4[i] = (0.5-y3[i])+((x4[i]-2)*(x4[i]-2)*0.02);
    }

    // create and configure plottables:
    QCPGraph *graph1 = ui->customPlot->addGraph();
    graph1->setData(x1, y1);
    graph1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 1.5), QBrush(Qt::white), 9));
    graph1->setPen(QPen(QColor(120, 120, 120), 2));

    QCPGraph *graph2 = ui->customPlot->addGraph();
    graph2->setData(x2, y2);
    graph2->setPen(Qt::NoPen);
    graph2->setBrush(QColor(200, 200, 200, 20));
    graph2->setChannelFillGraph(graph1);

    QCPBars *bars1 = new QCPBars(ui->customPlot->xAxis, ui->customPlot->yAxis);
    bars1->setWidth(9/(double)x3.size());
    bars1->setData(x3, y3);
    bars1->setPen(Qt::NoPen);
    bars1->setBrush(QColor(10, 140, 70, 160));

    QCPBars *bars2 = new QCPBars(ui->customPlot->xAxis, ui->customPlot->yAxis);
    bars2->setWidth(9/(double)x4.size());
    bars2->setData(x4, y4);
    bars2->setPen(Qt::NoPen);
    bars2->setBrush(QColor(10, 100, 50, 70));
    bars2->moveAbove(bars1);

    // move bars above graphs and grid below bars:
    ui->customPlot->addLayer("abovemain", ui->customPlot->layer("main"), QCustomPlot::limAbove);
    ui->customPlot->addLayer("belowmain", ui->customPlot->layer("main"), QCustomPlot::limBelow);
    graph1->setLayer("abovemain");
    ui->customPlot->xAxis->grid()->setLayer("belowmain");
    ui->customPlot->yAxis->grid()->setLayer("belowmain");

    // set some pens, brushes and backgrounds:
    ui->customPlot->xAxis->setBasePen(QPen(Qt::white, 1));
    ui->customPlot->yAxis->setBasePen(QPen(Qt::white, 1));
    ui->customPlot->xAxis->setTickPen(QPen(Qt::white, 1));
    ui->customPlot->yAxis->setTickPen(QPen(Qt::white, 1));
    ui->customPlot->xAxis->setSubTickPen(QPen(Qt::white, 1));
    ui->customPlot->yAxis->setSubTickPen(QPen(Qt::white, 1));
    ui->customPlot->xAxis->setTickLabelColor(Qt::white);
    ui->customPlot->yAxis->setTickLabelColor(Qt::white);
    ui->customPlot->xAxis->grid()->setPen(QPen(QColor(140, 140, 140), 1, Qt::DotLine));
    ui->customPlot->yAxis->grid()->setPen(QPen(QColor(140, 140, 140), 1, Qt::DotLine));
    ui->customPlot->xAxis->grid()->setSubGridPen(QPen(QColor(80, 80, 80), 1, Qt::DotLine));
    ui->customPlot->yAxis->grid()->setSubGridPen(QPen(QColor(80, 80, 80), 1, Qt::DotLine));
    ui->customPlot->xAxis->grid()->setSubGridVisible(true);
    ui->customPlot->yAxis->grid()->setSubGridVisible(true);
    ui->customPlot->xAxis->grid()->setZeroLinePen(Qt::NoPen);
    ui->customPlot->yAxis->grid()->setZeroLinePen(Qt::NoPen);
    ui->customPlot->xAxis->setUpperEnding(QCPLineEnding::esSpikeArrow);
    ui->customPlot->yAxis->setUpperEnding(QCPLineEnding::esSpikeArrow);
    QLinearGradient plotGradient;
    plotGradient.setStart(0, 0);
    plotGradient.setFinalStop(0, 350);
    plotGradient.setColorAt(0, QColor(80, 80, 80));
    plotGradient.setColorAt(1, QColor(50, 50, 50));
    ui->customPlot->setBackground(plotGradient);
    QLinearGradient axisRectGradient;
    axisRectGradient.setStart(0, 0);
    axisRectGradient.setFinalStop(0, 350);
    axisRectGradient.setColorAt(0, QColor(80, 80, 80));
    axisRectGradient.setColorAt(1, QColor(30, 30, 30));
    ui->customPlot->axisRect()->setBackground(axisRectGradient);

    ui->customPlot->rescaleAxes();
    ui->customPlot->yAxis->setRange(0, 2);
}

// called from create_Data_for_Cross_and_Self_graphs()
void PlotDialog::create_graphs(QMap<double, double> E_Scor,QString GraphLeg, QString Graph_Title, int WhichGraph){
    
    if(E_Scor.isEmpty()){
        QMessageBox::information(this, tr(""), " Canno't find data of " + GraphLeg + " for graph: "+ Graph_Title);
        return;
    }

    //QTextStream(stdout) << "---------------- create_graphs() ------------------" << "\n";

    // to set the new graph Title
    int index = ui->customPlot->graphCount();
    ui->customPlot->plotLayout()->remove(title);
    title = new QCPTextElement(ui->customPlot, Graph_Title , QFont("sans", 17, QFont::Bold));
    ui->customPlot->plotLayout()->addElement(0, 0, title);
    ui->customPlot->plotLayout()->updateLayout();

    ui->customPlot->legend->setVisible(true); // because in ratio plot, we dont use it and set it to false

    connect(title, SIGNAL(doubleClicked(QMouseEvent*)), this, SLOT(titleDoubleClick(QMouseEvent*)));
    
    QMap<double,double>::iterator beginn = E_Scor.begin();
    QMap<double,double>::iterator endd = E_Scor.end();
    
    //QVector<double> x(number_of_energies), y(number_of_energies);
    QVector<double> x(E_Scor.size()), y(E_Scor.size());
    
    QTextStream(stdout) << "E_Scor.size() : " << E_Scor.size() << "\n";
    
    int k = 0;
    while(beginn != endd)
    {
        
        x[k] = beginn.key();
        y[k] = beginn.value();
        QTextStream(stdout) << x[k] << "  "<< y[k] << "\n";
        
        k++;
        beginn++;
    }
    
    QPen pen;
    ui->customPlot->setNoAntialiasingOnDrag(true); // more performance/responsiveness during dragging
    // add graphs with different scatter styles:
    
    int shapeIndex = int(((rand() % 10000) / 10000.0)*14); /* ////////////////////////////////////////////////////////////////////////// */

    ui->customPlot->addGraph();
    pen.setColor(QColor(qSin(shapeIndex*0.3)*100+100, qSin(shapeIndex*0.6+0.7)*100+100, qSin(shapeIndex*0.4+0.6)*100+100));
    pen.setWidth(2);

    ui->customPlot->graph()->setData(x, y);
    if(GraphTypeShown == "OneE" || GraphTypeShown == "RadioTracerQuantityGraph"){
        // prepare x axis with country labels:
        ui->customPlot->xAxis->setScaleType(QCPAxis::stLinear);
        QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
        textTicker->addTicks(xticks, xlabels);
        ui->customPlot->xAxis->setTicker(textTicker);
        ui->customPlot->xAxis->setTickLabelRotation(60);
        ui->customPlot->xAxis->setSubTicks(false);
    }else{

        QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);

        //ui->customPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
        ui->customPlot->xAxis->setTicker(logTicker);
        //ui->customPlot->xAxis->setNumberFormat(XNumberFormat); // e = exponential, b = beautiful decimal powers
        //ui->customPlot->xAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"

    }

    ui->customPlot->graph()->rescaleAxes(true);
    
    ui->customPlot->graph()->setPen(pen);
    //ui->customPlot->graph()->setName(QCPScatterStyle::staticMetaObject.enumerator(QCPScatterStyle::staticMetaObject.indexOfEnumerator("ScatterShape")).valueToKey(shapes.at(shapeIndex)));
    ui->customPlot->graph()->setName(GraphLeg);
    ui->customPlot->graph()->setLineStyle(QCPGraph::lsLine);
    // set scatter style:
    if (shapes.at(shapeIndex) != QCPScatterStyle::ssCustom){
        
        ui->customPlot->graph()->setScatterStyle(QCPScatterStyle(shapes.at(shapeIndex), 10));// 15 is the size of shapes vector
    }
    else{
        
        QPainterPath customScatterPath;
        for (int i=0; i<3; ++i)
            customScatterPath.cubicTo(qCos(2*M_PI*i/3.0)*9, qSin(2*M_PI*i/3.0)*9, qCos(2*M_PI*(i+0.9)/3.0)*9, qSin(2*M_PI*(i+0.9)/3.0)*9, 0, 0);
        ui->customPlot->graph()->setScatterStyle(QCPScatterStyle(customScatterPath, QPen(Qt::black, 0), QColor(40, 70, 255, 50), 10));
    }

    if(ui->checkBoxAddErrorBar->isChecked() && UseErrorBar){

        QVector<double> err(ErrorMap.size());
        QMap<double,double>::iterator beginn = ErrorMap.begin();
        QMap<double,double>::iterator endd = ErrorMap.end();

        int k = 0;
        while(beginn != endd)
        {
            err[k] = beginn.value();
            QTextStream(stdout) << err[k] << "\n";
            k++; beginn++;
        }

        QCPErrorBars *errorBars = new QCPErrorBars(ui->customPlot->xAxis, ui->customPlot->yAxis);
        errorBars->removeFromLegend();
        errorBars->setAntialiased(false);
        errorBars->setDataPlottable(ui->customPlot->graph(index));

        QPen graphPen;
        graphPen.setColor(QColor(0,0,0));
        graphPen.setWidth(2);
        errorBars->setPen(graphPen);
        errorBars->setData(err);
        ErrorMap.clear();
        UseErrorBar = false;
    }

    ui->customPlot->replot();
    
    //QTextStream(stdout) << " ----------------------------------------- " << "\n";
}

void PlotDialog::GetPlotInputDataAndReadFiles(){

    ReferenceFilePath = ui->PlotLineEditReferenceFile->text() ;
    ReferenceFilePath_2 = ui->PlotLineEditReferenceFile_2->text() ;
    ResultFilePath = ui->PlotLineEditResultFile->text();
    AnalysisInputFilePath = ui->PlotLineEditSimulationInpFile->text();
    CrossSectionFilePath = ui->PlotLineEditCrossSectionFile->text();
    
    QTextStream(stdout) << "ReferenceFilePath " << ReferenceFilePath << "\n";
    QTextStream(stdout) << "ReferenceFilePath2 " << ReferenceFilePath_2 << "\n";
    QTextStream(stdout) << "ResultFilePath " << ResultFilePath << "\n";
    QTextStream(stdout) << "AnalysisInputFilePath " << AnalysisInputFilePath << "\n";
    QTextStream(stdout) << "CrossSectionFilePath " << CrossSectionFilePath << "\n";
    
    SourceOrgan = ui->comboBoxSourceOrgan->currentText();
    TargetOrgan = ui->PlotcomboBoxTargetOrgan->currentText();
    ParticleName = ui->plotcomboBoxParticle->currentText();
    RadioTracerName = ui->plotcomboBoxParticle->currentText() ;
    Compare_Type = ui->plotcomboBoxCompareType->currentText() ;
    GraphsData = ui->PlotcomboBoxGraphData->currentText() ;
    
    QStringList InputsVals = ui->plotlineEditReferenceName->text().split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    if(InputsVals.size() > 0){
        CompareReferenceName = InputsVals[0];
        if(InputsVals.size() > 1){
            CompareReferenceName_2 = InputsVals[1];
        }else{
            CompareReferenceName_2= "";
        }
    }
    
    QuantitiesToScore = ui->plotcomboBoxScoredVariable->currentText() ;
    GraphForEnergy = ui->plotComboBoxEnergies->currentText().toDouble();
    GeometrySymbol = ui->comboBoxGeometrySymbole->currentText();

    RegionVariableName = ui->PlotComboboxRegionVariable->currentText();
    
    ParticleForCrossSection = ui->comboBoxParticleForCrossSection->currentText();
    MaterialForCrossSection = ui->comboBoxMaterialForCrossSection->currentText();
    ProcessForCrossSection = ui->comboBoxProcessForCrossSection->currentText();
    
    if(ui->radioButtonLRD->isChecked()){
        DiffSym = "LRD";
        DiffExp = "Logarithmic Relative Difference (\\%)";
    }else if(ui->radioButtonRa->isChecked()){
        DiffSym = "RA";
        DiffExp = "Ratio";
    }else{
        DiffSym = "RD";
        DiffExp = "Relative Difference (\\%)";
    }

    QTextStream(stdout) << "SourceOrgan " << SourceOrgan << "\n";
    QTextStream(stdout) << "TargetOrgan " << TargetOrgan << "\n";
    QTextStream(stdout) << "ParticleName " <<  ParticleName << "\n";
    QTextStream(stdout) << "Compare_Type " << Compare_Type << "\n";
    QTextStream(stdout) << "GraphsData " << GraphsData << "\n";
    QTextStream(stdout) << "CompareReferenceName " << CompareReferenceName << "\n";
    QTextStream(stdout) << "CompareReferenceName_2 " << CompareReferenceName_2 << "\n";
    QTextStream(stdout) << "QuantitiesToScore " << QuantitiesToScore << "\n";
    QTextStream(stdout) << "RegionVariableName " << RegionVariableName << "\n";
    QTextStream(stdout) << "ParticleForCrossSection " << ParticleForCrossSection << "\n";
    QTextStream(stdout) << "MaterialForCrossSection " << MaterialForCrossSection << "\n";
    QTextStream(stdout) << "ProcessForCrossSection " << ProcessForCrossSection << "\n";
    
    if(RegionVariableName == "Volume"){
        RegionVariableNameWithUnit="Volume(cm3)";
    }
    else if(RegionVariableName == "Mass"){
        RegionVariableNameWithUnit="Mass(kg)";
    }
    else if(RegionVariableName == "Density"){
        RegionVariableNameWithUnit="Density(mg/cm3)";
    }
    else if(RegionVariableName == "Distance"){
        RegionVariableNameWithUnit="Distance(mm)";
    }
    
    if(ui->checkBoxAppendGraph->isChecked()){ // in this case ,the graphs are appended to the graph plot
    }else{
        //index = 0;
        ui->customPlot->legend->clear();
        ui->customPlot->clearGraphs();
        ui->customPlot->clearPlottables();

        QTextStream(stdout) << "RegionVariableNameWithUnit " << RegionVariableNameWithUnit << "\n";

        if(QuantitiesToScore == "SAF"){
            ScoreVariableUnit = QuantitiesToScore + "(Kg-1)";
        }
        else if(QuantitiesToScore == "AF"){
            ScoreVariableUnit = QuantitiesToScore;
        }
        else if(QuantitiesToScore == "AD"){
            ScoreVariableUnit = QuantitiesToScore + "(MeV/Kg)";
        }
        else if(QuantitiesToScore == "AE"){
            ScoreVariableUnit = QuantitiesToScore + "(MeV)";
        }
        else if(QuantitiesToScore == "S"){
            ScoreVariableUnit = QuantitiesToScore + "(MeV/kg-emission)";
        }
        else if(QuantitiesToScore == "E"){
            ScoreVariableUnit = QuantitiesToScore + "(Sv)";
        }
        else if(QuantitiesToScore == "H"){
            ScoreVariableUnit = QuantitiesToScore + "(Sv)";
        }
        else if(QuantitiesToScore == "DR"){
            ScoreVariableUnit = QuantitiesToScore + "(%)";
        }
        else if(QuantitiesToScore == "ER"){
            ScoreVariableUnit = QuantitiesToScore + "(%)";
        }
        else {
            ScoreVariableUnit = "f(E) (unit)";
        }

        QTextStream(stdout) << "ScoreVariableUnit " << ScoreVariableUnit << "\n";

        QTextStream(stdout) << "GraphTypeShown " << GraphTypeShown << "\n";

        if(GraphTypeShown == "SC"){
            ui->customPlot->xAxis->setLabel("E(MeV)"); // to set the new axis quantity
            ui->customPlot->yAxis->setLabel(ScoreVariableUnit); // to set the new axis quantity
        }
        else if(GraphTypeShown == "RadioTracerQuantityGraph"){
            ui->customPlot->xAxis->setLabel("Regions names"); // to set the new axis quantity
            ui->customPlot->yAxis->setLabel(ScoreVariableUnit); // to set the new axis quantity
        }
        else if(GraphTypeShown == "OneE"){
            ui->customPlot->xAxis->setLabel("Regions names"); // to set the new axis quantity
            ui->customPlot->yAxis->setLabel(ScoreVariableUnit); // to set the new axis quantity
        }
        else if(GraphTypeShown == "MCS"){
            ui->customPlot->xAxis->setLabel("E(MeV)"); // to set the new axis quantity
            ui->customPlot->yAxis->setLabel("Macroscopic cross section (cm-1)"); // to set the new axis quantity
        }
        else if(GraphTypeShown == "V"){
            ui->customPlot->xAxis->setLabel(RegionVariableNameWithUnit); // to set the new axis quantity
            ui->customPlot->yAxis->setLabel(ScoreVariableUnit); // to set the new axis quantity
        }
        else if(GraphTypeShown == "RE"){
            ui->customPlot->xAxis->setLabel("E(MeV)"); // to set the new axis quantity
            ui->customPlot->yAxis->setLabel(DiffExp); // to set the new axis quantity
        }
        else if(GraphTypeShown == "PB"){
            ui->customPlot->yAxis->setLabel(ScoreVariableUnit); // to set the new axis quantity
            ui->customPlot->xAxis->setLabel(""); // to set the new axis quantity
        }
        else if(GraphTypeShown == "Rat"){
            ui->customPlot->yAxis->setLabel("Ratio"); // to set the new axis quantity
            ui->customPlot->xAxis->setLabel(""); // to set the new axis quantity
        }
        else if(GraphTypeShown == "RSD"){
            ui->customPlot->xAxis->setLabel("E(MeV)"); // to set the new axis quantity
            ui->customPlot->yAxis->setLabel("Rel SDv %"); // to set the new axis quantity
        }
        else if(GraphTypeShown == "T"){
            ui->customPlot->xAxis->setLabel("E(MeV)"); // to set the new axis quantity
            ui->customPlot->yAxis->setLabel("Computation Time (min)"); // to set the new axis quantity
        }
        else if(GraphTypeShown == "XY"){
            ui->customPlot->xAxis->setLabel("X(Unit)"); // to set the new axis quantity
            ui->customPlot->yAxis->setLabel("Y(Unit)"); // to set the new axis quantity
        }
        else{
            ui->customPlot->xAxis->setLabel("E(MeV)"); // to set the new axis quantity
            ui->customPlot->yAxis->setLabel("f(E)"); // to set the new axis quantity
        }
    }
}
void PlotDialog::setComboboxInputs()
{
    ScoreVarlist.clear();ScoreVarlist.empty();
    Geometrylist.clear();Geometrylist.empty();
    Particlelist.clear();Particlelist.empty();
    SourceOrganlist.clear();SourceOrganlist.empty();
    TargetOrganlist.clear();TargetOrganlist.empty();
    QuantityGeometryParticleSourceTargetsNames.clear();QuantityGeometryParticleSourceTargetsNames.empty();

    if(ResultTable.size() == 0){
        QMessageBox::information(this, tr(""), "Canno't find result data to fill the components for particles, check the file path or file data");
        ui->plotMessageLabel->setText("Canno't save data from file, check the file path or file data");
    }
    else if(ResultQuantityGeometryRadioTracerSourceTargetValues.size() == 0 ){

        QMessageBox::information(this, tr(""), "Canno't find result data to fill the components for radiotracers, check the file path or file data");
        ui->plotMessageLabel->setText("Canno't save data from file, check the file path or file data");
    }
    else{
        ui->plotMessageLabel->setText("The source, target and particle chooser are filled according to the content of result file data");
    }

    int ff = 0;
    for ( auto obeg = ResultTable.begin(); obeg != ResultTable.end(); ++obeg  ){

        QString VAR_NAME = obeg.key();
        QTextStream(stdout) <<"VAR_NAME : " << VAR_NAME << "\n";

        ff=0;
        for(int i=0; i<ScoreVarlist.size(); i++)
        {
            if(ScoreVarlist.indexOf(QRegExp(VAR_NAME, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                ff++;
            }
        }
        if(ff == 0){
            ScoreVarlist.push_back(VAR_NAME);
        }

        for ( auto abeg = obeg.value().begin(); abeg != obeg.value().end(); ++abeg  ){

            QString Geo_NAME = abeg.key();
            //QTextStream(stdout) <<"Geo_NAME : " << Geo_NAME << "\n";

            ff=0;
            for(int i=0; i<Geometrylist.size(); i++)
            {
                if(Geometrylist.indexOf(QRegExp(Geo_NAME, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                    ff++;
                }
            }
            if(ff == 0){
                Geometrylist.push_back(Geo_NAME);
            }
            for ( auto Abeg = abeg.value().begin(); Abeg != abeg.value().end(); ++Abeg  ){

                QString PARTICLE_NAME = Abeg.key();
                //QTextStream(stdout) <<"PARTICLE_NAME : " << PARTICLE_NAME << "\n";

                ff=0;
                for(int i=0; i<Particlelist.size(); i++)
                {
                    if(Particlelist.indexOf(QRegExp(PARTICLE_NAME, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                        ff++;
                    }
                }
                if(ff == 0){
                    Particlelist.push_back(PARTICLE_NAME);
                }

                // iterations on source name
                for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){

                    QString Source_ORG = Bbeg.key();
                    //QTextStream(stdout) <<"Source_ORG : " << Source_ORG << "\n";

                    ff=0;
                    for(int i=0; i<SourceOrganlist.size(); i++)
                    {
                        if(SourceOrganlist.indexOf(QRegExp(Source_ORG, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                            ff++;
                        }
                    }
                    if(ff == 0){
                        SourceOrganlist.push_back(Source_ORG);
                    }

                    // iterations on target name
                    for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg  ){

                        QString Target_ORG = Cbeg.key();
                        //QTextStream(stdout) <<"Target_ORG : " << Target_ORG << "\n";

                        QuantityGeometryParticleSourceTargetsNames[VAR_NAME][Geo_NAME][PARTICLE_NAME][Source_ORG][Target_ORG]=1;

                        ff=0;
                        for(int i=0; i<TargetOrganlist.size(); i++)
                        {
                            if(TargetOrganlist.indexOf(QRegExp(Target_ORG, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                                ff++;
                            }
                        }
                        if(ff == 0){
                            TargetOrganlist.push_back(Target_ORG);
                        }
                    }
                }
            }
        }
    }

    for ( auto obeg = ResultQuantityGeometryRadioTracerSourceTargetValues.begin(); obeg != ResultQuantityGeometryRadioTracerSourceTargetValues.end(); ++obeg  ){

        QString VAR_NAME = obeg.key();
        //QTextStream(stdout) <<"VAR_NAME : " << VAR_NAME << "\n";

        ff=0;
        for(int i=0; i<ScoreVarlist.size(); i++)
        {
            if(ScoreVarlist.indexOf(QRegExp(VAR_NAME, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                ff++;
            }
        }
        if(ff == 0){
            ScoreVarlist.push_back(VAR_NAME);
        }

        for ( auto abeg = obeg.value().begin(); abeg != obeg.value().end(); ++abeg  ){

            QString Geo_NAME = abeg.key();
            //QTextStream(stdout) <<"Geo_NAME : " << Geo_NAME << "\n";

            ff=0;
            for(int i=0; i<Geometrylist.size(); i++)
            {
                if(Geometrylist.indexOf(QRegExp(Geo_NAME, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                    ff++;
                }
            }
            if(ff == 0){
                Geometrylist.push_back(Geo_NAME);
            }

            for ( auto Abeg = abeg.value().begin(); Abeg != abeg.value().end(); ++Abeg  ){

                QString PARTICLE_NAME = Abeg.key();
                //QTextStream(stdout) <<"PARTICLE_NAME : " << PARTICLE_NAME << "\n";

                ff=0;
                for(int i=0; i<Particlelist.size(); i++)
                {
                    if(Particlelist.indexOf(QRegExp(PARTICLE_NAME, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                        ff++;
                    }
                }
                if(ff == 0){
                    Particlelist.push_back(PARTICLE_NAME);
                }

                // iterations on source name
                for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){

                    QString Source_ORG = Bbeg.key();
                    //QTextStream(stdout) <<"Source_ORG : " << Source_ORG << "\n";

                    ff=0;
                    for(int i=0; i<SourceOrganlist.size(); i++)
                    {
                        if(SourceOrganlist.indexOf(QRegExp(Source_ORG, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                            ff++;
                        }
                    }
                    if(ff == 0){
                        SourceOrganlist.push_back(Source_ORG);
                    }

                    // iterations on target name
                    for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg  ){

                        QString Target_ORG = Cbeg.key();
                        //QTextStream(stdout) <<"Target_ORG : " << Target_ORG << "\n";

                        QuantityGeometryParticleSourceTargetsNames[VAR_NAME][Geo_NAME][PARTICLE_NAME][Source_ORG][Target_ORG]=1;

                        ff=0;
                        for(int i=0; i<TargetOrganlist.size(); i++)
                        {
                            if(TargetOrganlist.indexOf(QRegExp(Target_ORG, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                                ff++;
                            }
                        }
                        if(ff == 0){
                            TargetOrganlist.push_back(Target_ORG);
                        }
                    }
                }
            }
        }
    }

    QTextStream(stdout) <<"----------------------------- The comboboxes list of data ----------------------------------" << "\n";

    /*
    for(int i=0; i<ScoreVarlist.size(); i++)
    {
        //QTextStream(stdout) <<"ScoreVarlist[] : " << ScoreVarlist[i] << "\n";
    }
    for(int i=0; i<Geometrylist.size(); i++)
    {
        //QTextStream(stdout) <<"Geometrylist[] : " << Geometrylist[i] << "\n";
    }
    for(int i=0; i<Particlelist.size(); i++)
    {
        //QTextStream(stdout) <<"Particlelist[] : " << Particlelist[i] << "\n";
    }
    for(int i=0; i<SourceOrganlist.size(); i++)
    {
        if(OrganNamesToScore.size() == 0){
            OrganNamesToScore.push_back(SourceOrganlist[i]);
        }
        //QTextStream(stdout) <<"SourceOrganlist[] : " << SourceOrganlist[i] << "\n";
    }
    for(int i=0; i<TargetOrganlist.size(); i++)
    {
        //QTextStream(stdout) <<"TargetOrganlist[] : " << TargetOrganlist[i] << "\n";
    }
*/

    ui->plotcomboBoxScoredVariable->clear();
    ui->comboBoxGeometrySymbole->clear();
    ui->plotcomboBoxParticle->clear();
    ui->comboBoxSourceOrgan->clear();
    ui->PlotcomboBoxTargetOrgan->clear();

    ui->plotcomboBoxScoredVariable->addItems(ScoreVarlist);
    ui->comboBoxGeometrySymbole->addItems(Geometrylist);
    ui->plotcomboBoxParticle->addItems(Particlelist);
    ui->comboBoxSourceOrgan->addItems(SourceOrganlist);
    ui->PlotcomboBoxTargetOrgan->addItems(TargetOrganlist);

    int dd = 0;

    QTextStream(stdout) << "\n\n\n\n\n QuantityGeometryParticleSourceTargetsNames size : " << QuantityGeometryParticleSourceTargetsNames.size() << "\n";

    for ( auto obeg = QuantityGeometryParticleSourceTargetsNames.begin(); obeg != QuantityGeometryParticleSourceTargetsNames.end(); ++obeg  ){
        QString VAR_NAME = obeg.key();
        if(dd==1){break;}
        //QTextStream(stdout) << " VAR_NAME " << VAR_NAME << "\n";

        for ( auto abeg = obeg.value().begin(); abeg != obeg.value().end(); ++abeg  ){
            QString Geo_NAME = abeg.key();
            //QTextStream(stdout) << " Geo_NAME " << Geo_NAME << "\n";
            if(dd==1){break;}

            for ( auto Abeg = abeg.value().begin(); Abeg != abeg.value().end(); ++Abeg  ){
                QString PARTICLE_NAME = Abeg.key();
                //QTextStream(stdout) << " PARTICLE_NAME " << PARTICLE_NAME << "\n";
                if(dd==1){break;}

                // iterations on source name
                for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                    QString Source_ORG = Bbeg.key();
                    //QTextStream(stdout) << " Source_ORG " << Source_ORG << "\n";
                    if(dd==1){break;}

                    // iterations on target name
                    for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg  ){
                        QString Target_ORG = Cbeg.key();
                        //QTextStream(stdout) << " Target_ORG " << Target_ORG << "\n";
                        if(dd==1){break;}

                        ui->plotcomboBoxScoredVariable->setCurrentText(VAR_NAME);
                        ui->comboBoxGeometrySymbole->setCurrentText(Geo_NAME);
                        ui->plotcomboBoxParticle->setCurrentText(PARTICLE_NAME);
                        ui->comboBoxSourceOrgan->setCurrentText(Source_ORG);
                        ui->PlotcomboBoxTargetOrgan->setCurrentText(Target_ORG);
                        /*
                        QTextStream(stdout) << " \n\n\n \n\n\n \n\n\n ------------------------------ \n\n\n " << "\n";

                        QTextStream(stdout) << " VAR_NAME " << VAR_NAME
                                            << " Geo_NAME " << Geo_NAME
                                             << " PARTICLE_NAME " << PARTICLE_NAME
                                              << " Source_ORG " << Source_ORG
                                               << " Target_ORG " << Target_ORG << "\n";
*/
                        dd++;
                    }
                }
            }
        }
    }



}

// get the graphs data for each source target particle for cross and self and call create_graphs for each graph
// it call create_graphs()
// called when generate graphs is clicked
void PlotDialog::create_Data_for_Cross_and_Self_graphs(QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> table, QString dataReferenceFileName){ // 1 for Self , 2 for cross
    
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << "() ------------------" << "\n";
    
    ui->customPlot->update();
    
    QString Source_ORG;
    QString Target_ORG;
    double ENERGY_VAL;
    double VAL;
    
    //double MinY = 0.00001;
    
    QTextStream(stdout) << "Table size " << table.size() << " table[QuantitiesToScore][GeometrySymbol][ParticleName].size(): " << table[QuantitiesToScore][GeometrySymbol][ParticleName].size() << "\n";

    // iterations on source name
    for ( auto Bbeg = table[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Bbeg != table[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Bbeg  ){
        
        Source_ORG = Bbeg.key();
        QTextStream(stdout) <<"Source_ORG : " << Source_ORG << "\n";
        
        // iterations on target name
        for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg  ){
            
            Target_ORG = Cbeg.key();
            QTextStream(stdout) <<"Target_ORG : " << Target_ORG << "\n";
            
            QMap<double, double> Energy_ScoreVar;
            QMap<double, double> ScoreVar_Err;

            // iterations on particle energies
            //for ( auto Dbeg = Cbeg.value().begin(); Dbeg != Cbeg.value().end(); ++Dbeg  ){
            //ENERGY_VAL = Dbeg.key();
            //VAL = Dbeg.value() ;
            
            for(int ii = 0; ii < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size(); ii++){

                ErrorMap[ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]] = ErrorTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]];
                if(table[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]] == 0 || table[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]] == MinValForLog ){
                    continue;
                }
                else{
                    Energy_ScoreVar[ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]] = table[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]] ;
                }
            }
            
            QTextStream(stdout) <<"Source_ORG " << Source_ORG << " SourceOrgan "<< SourceOrgan << " Target_ORG "<< Target_ORG << "\n";
            
            if(GraphsData == "Result"){
                if(Compare_Type == "Self"){ // we dont need the source organ because we will show all source organs, we have to just Source_ORG == Target_ORG
                    if(Source_ORG == Target_ORG){
                        
                        QTextStream(stdout) <<GraphsData << " " << Compare_Type << "\n";
                        
                        QString graph_legend = Source_ORG;
                        create_graphs( Energy_ScoreVar , graph_legend, graphs_Title , 0);
                    }
                }
                else if(Compare_Type == "Cross"){
                    if(Source_ORG == SourceOrgan){ // we need the source but all the target organs
                        
                        QTextStream(stdout) <<GraphsData << " " << Compare_Type << "\n";
                        
                        QString graph_legend = Target_ORG;
                        create_graphs( Energy_ScoreVar , graph_legend , graphs_Title , 0);
                    }
                }
            }
            else if(GraphsData == "Reference_Result"){
                if(Compare_Type == "Self"){
                    if(Source_ORG == SourceOrgan && Source_ORG == Target_ORG){ // the Source_ORG and Target_ORG have to be the same
                        
                        QTextStream(stdout) <<GraphsData << " " << Compare_Type << "\n";
                        
                        QString graph_legend = dataReferenceFileName;
                        create_graphs( Energy_ScoreVar , graph_legend, graphs_Title , 0);
                    }
                }
                else if(Compare_Type == "Cross"){
                    if(Source_ORG == SourceOrgan && Target_ORG == TargetOrgan){
                        
                        QTextStream(stdout) <<GraphsData << " " << Compare_Type << "\n";
                        
                        QString graph_legend = dataReferenceFileName;
                        create_graphs( Energy_ScoreVar , graph_legend, graphs_Title , 0);
                    }
                }
            }
        }
    }
}
void PlotDialog::create_Data_for_Cross_and_Self_ForOneEnergy_graphs(){ // 1 for Self , 2 for cross

    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << "() ------------------" << "\n";

    ui->customPlot->update();

    QMap<double, double> Energy_ScoreVar;
    QMap<double, double> RefEnergy_ScoreVar;
    QMap<double, double> Ref2Energy_ScoreVar;

    // iterations on source name
    for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  ){

        QString Source_ORG = Abeg.key();
        //QTextStream(stdout) <<"Source_ORG : " << Source_ORG << "\n";

        // for just Reference_Result we will show just the selected source organ
        if(GraphsData == "Reference_Result" && SourceOrgan != Source_ORG){continue;}

        xticks.clear();xticks.empty();
        xlabels.clear();xlabels.empty();
        double iii = 1;
        for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){

            QString Target_ORG = Bbeg.key();
            //Energy_ScoreVar[RegionParameterValueMap[GeometrySymbol][RegionVariableName][Target_ORG]] = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][GraphForEnergy];
            Energy_ScoreVar[iii] = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][GraphForEnergy];

            QTextStream(stdout) << " RadioTracerName " << RadioTracerName << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " iii " << iii << " " << RegionParameterValueMap[GeometrySymbol][RegionVariableName][Target_ORG] << "  "<< ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][GraphForEnergy] << "\n";

            if(GraphsData == "Reference_Result"){
                //RefEnergy_ScoreVar[RegionParameterValueMap[GeometrySymbol][RegionVariableName][Target_ORG]] = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][GraphForEnergy];
                RefEnergy_ScoreVar[iii] = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][GraphForEnergy];
                if(CompareReferenceName_2 != "" && QFile::exists(ui->PlotLineEditReferenceFile_2->text())){
                    //Ref2Energy_ScoreVar[RegionParameterValueMap[GeometrySymbol][RegionVariableName][Target_ORG]] = ReferenceTable_2[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][GraphForEnergy];
                    Ref2Energy_ScoreVar[iii] = ReferenceTable_2[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][GraphForEnergy];
                }

            }
            //xticks.push_back(RegionParameterValueMap[GeometrySymbol][RegionVariableName][Target_ORG]);
            xticks.push_back(iii);
            iii = iii + 1. ;
            xlabels.push_back(Target_ORG);
        }

        if(GraphsData == "Result"){
            graphs_Title = QuantitiesToScore +" from " + ParticleName + " of energy " + QString::number(GraphForEnergy) + " from sources to targets";
        }
        else{
            graphs_Title = QuantitiesToScore +" from " + ParticleName + " of energy " + QString::number(GraphForEnergy) + " from " +  Source_ORG + " source to targets ";
        }

        QString graph_legend = Source_ORG;
        if(GraphsData == "Reference_Result"){
            graph_legend = "DoseCalcs "+ Source_ORG;
        }
        QTextStream(stdout) << "Energy_ScoreVar.size() : " << Energy_ScoreVar.size() << "\n";

        create_graphs( Energy_ScoreVar , graph_legend, graphs_Title , 0);

        if(GraphsData == "Reference_Result"){
            graph_legend = CompareReferenceName +" "+ Source_ORG;
            create_graphs( RefEnergy_ScoreVar , graph_legend, graphs_Title , 0);
            if(CompareReferenceName_2 != "" && QFile::exists(ui->PlotLineEditReferenceFile_2->text())){
                graph_legend = CompareReferenceName_2 +" "+ Source_ORG;
                create_graphs( RefEnergy_ScoreVar , graph_legend, graphs_Title , 0);
            }
        }
    }
}
void PlotDialog::create_Data_for_RelativeErr_graphs(){ // if one value is empty or 0 the difference will be set to 0 and not 100 like in other difference calculations, because that 100 in graphs is not convienient
    
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << "() ------------------" << "\n";
    
    ui->customPlot->update();

    //QTextStream(stdout) << "ResultTable[QuantitiesToScore][GeometrySymbol].size " << ResultTable[QuantitiesToScore][GeometrySymbol].size() << "\n" ;
    
    // iterations on source name
    for ( auto Bbeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Bbeg != ResultTable[QuantitiesToScore][GeometrySymbol][ ParticleName].end(); ++Bbeg  ){
        
        QString Source_ORG = Bbeg.key();
        //QTextStream(stdout) << "Source_ORG " << Source_ORG << "\n" ;
        
        // iterations on target name
        for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg ){
            
            QString Target_ORG = Cbeg.key();
            //QTextStream(stdout) << "Target_ORG " << Target_ORG << "\n" ;
            
            // iterations on particle energies
            for(int ii = 0; ii < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size(); ii++){
                
                double a3 = RelativeDifferenceCalculation(ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]],ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]]);
                if(a3 != NULL){
                    ResRefErrCompTables[QuantitiesToScore][ParticleName][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]] = a3;
                }
            }
        }
    }
    
    // iterations on target name
    for ( auto Bbeg = ResRefErrCompTables[QuantitiesToScore][ParticleName].begin(); Bbeg != ResRefErrCompTables[QuantitiesToScore][ParticleName].end(); ++Bbeg  )
    {
        
        QString Source_ORG = Bbeg.key();
        
        // iterations on particle name
        for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg  )
        {
            
            QString Target_ORG = Cbeg.key();
            
            if(GraphsData == "Reference_Result"){
                if(Compare_Type == "Self"){
                    if(Source_ORG == Target_ORG){ // the Source_ORG and Target_ORG have to be the same
                        
                        //QTextStream(stdout) <<Source_ORG << " " << SourceOrgan << "\n";
                        
                        graphs_Title = "Self Relative Difference with " + CompareReferenceName + " for " + ParticleName;
                        //Energy(MeV); Comparison Difference Factor %";
                        QString graph_legend = Source_ORG;
                        create_graphs( Cbeg.value() , graph_legend, graphs_Title , 0);
                    }
                }
                else if(Compare_Type == "Cross"){
                    if(Source_ORG == SourceOrgan ){
                        
                        //QTextStream(stdout) <<Source_ORG << " " << SourceOrgan << "\n";
                        
                        graphs_Title = "Cross Relative Difference with " + CompareReferenceName + " for " + ParticleName;
                        QString graph_legend = Source_ORG+"-->"+Target_ORG;
                        create_graphs( Cbeg.value() , graph_legend, graphs_Title , 0);
                    }
                }
            }
        }
    }
}
void PlotDialog::create_Data_for_MCError_graphs(){
    
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << "() ------------------" << "\n";
    
    ui->customPlot->update();
    
    QString Source_ORG;
    QString Target_ORG;
    
    
    // iterations on source name
    for ( auto Bbeg = ErrorTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Bbeg != ErrorTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Bbeg  ){
        
        Source_ORG = Bbeg.key();
        //QTextStream(stdout) << "Source_ORG " << Source_ORG << "\n" ;
        
        // iterations on target name
        for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg ){
            
            Target_ORG = Cbeg.key();
            //QTextStream(stdout) << "Target_ORG " << Target_ORG << "\n" ;
            
            if(GraphsData == "Reference_Result" || GraphsData == "Result"){
                if(Compare_Type == "Self"){
                    if(Source_ORG == Target_ORG){ // the Source_ORG and Target_ORG have to be the same
                        
                        //QTextStream(stdout) <<Source_ORG << " " << SourceOrgan << "\n";
                        
                        graphs_Title = "Self Relative SDv for " +  ParticleName;
                        
                        QString graph_legend = Source_ORG;
                        create_graphs( Cbeg.value() , graph_legend, graphs_Title , 0);
                    }
                }
                else if(Compare_Type == "Cross"){
                    if(Source_ORG == SourceOrgan ){
                        
                        //QTextStream(stdout) <<Source_ORG << " " << SourceOrgan << "\n";
                        
                        graphs_Title = Source_ORG + " Cross Relative SDv for " + ParticleName;
                        QString graph_legend = Source_ORG+"-->"+Target_ORG;
                        create_graphs( Cbeg.value() , graph_legend, graphs_Title , 0);
                    }
                }
            }
        }
    }
}
void PlotDialog::create_Data_for_RegionParametre_graphs(){
    
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << " ------------------" << "\n";
    
    ui->customPlot->update();

    QTextStream(stdout) << "AnalysisInputFilePath    : " << AnalysisInputFilePath << "\n" ;
    QTextStream(stdout) << "RegionVariableName    : " << RegionVariableName << "\n" ;
    QTextStream(stdout) << "RegionVariableNameWithUnit    : " << RegionVariableNameWithUnit << "\n" ;
    
    QMap<double,double> MapValToGraph;
    
    if(Compare_Type == "Self"){
        
        for ( int A = 0; A < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; A++  ){
            QTextStream(stdout) << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][A] << "\n";
            
            double xi[OrganNamesToScore.size()], yi[OrganNamesToScore.size()] ;
            for ( int B = 0; B < OrganNamesToScore.size() ; B++  ){
                
                xi[B] = RegionParameterValueMap[GeometrySymbol][RegionVariableName][OrganNamesToScore[B]] ;
                yi[B] = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][OrganNamesToScore[B]][OrganNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][A]];
                MapValToGraph[xi[B]] = yi[B] ;
                QTextStream(stdout) << xi[B] << " : " << yi[B] << "\n";
            }
            
            graphs_Title =  QuantitiesToScore + " versus " + RegionVariableName + " evolution in Self Absorption for "+ ParticleName;
            std::ostringstream nm; nm << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][A] << " MeV" ;
            create_graphs( MapValToGraph , nm.str().c_str() , graphs_Title , 0);
        }
    }
    else if(Compare_Type == "Cross"){
        
        for ( int A = 0; A < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; A++  ){
            QTextStream(stdout) << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][A] << "\n";
            
            double xi[OrganNamesToScore.size()], yi[OrganNamesToScore.size()] ;
            for ( int B = 0; B < OrganNamesToScore.size() ; B++  ){
                
                if(OrganNamesToScore[B] == SourceOrgan){
                    //continue;
                }
                xi[B] = RegionParameterValueMap[GeometrySymbol][RegionVariableName][OrganNamesToScore[B]] ;
                yi[B] = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourceOrgan][OrganNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][A]];
                MapValToGraph[xi[B]] = yi[B] ;
                QTextStream(stdout) << xi[B] << " : " << yi[B] << "\n";
            }
            
            graphs_Title =  QuantitiesToScore + " versus " + RegionVariableName + " evolution in Cross Orradiation from "+SourceOrgan + ", for" + ParticleName;
            std::ostringstream nm; nm << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][A] << " MeV" ;
            create_graphs( MapValToGraph , nm.str().c_str() , graphs_Title , 0);
        }
    }
    
}
void PlotDialog::create_Data_for_CrossSection_graphs(){
    
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << " ------------------" << "\n";
    
    ui->customPlot->update();
    
    QTextStream(stdout) << "MaterialForCrossSection    : " << MaterialForCrossSection << "\n" ;
    
    QMap<double,double> Energy_Sigma;
    
    //ParticleMaterialProcessEnergySigmaMap[ParticleForCrossSection][MaterialForCrossSection][ProcessForCrossSection];
    
    if(ProcessForCrossSection == "All"){
        
        graphs_Title = "Sigma for " + ParticleForCrossSection + " in " + MaterialForCrossSection + " Material";
        
        for ( auto Bbeg = ParticleMaterialProcessEnergySigmaMap[ParticleForCrossSection][MaterialForCrossSection].begin(); Bbeg != ParticleMaterialProcessEnergySigmaMap[ParticleForCrossSection][MaterialForCrossSection].end(); ++Bbeg  ){
            
            QString SigmaProc = Bbeg.key();
            QTextStream(stdout) <<"SigmaProc : " << SigmaProc << "\n";
            
            for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg  ){
                
                double Ene = Cbeg.key();
                double Sigma = Cbeg.value();
                
                Energy_Sigma[Ene] = Sigma ;
            }
            
            QString graph_legend = SigmaProc;
            create_graphs( Energy_Sigma , graph_legend, graphs_Title , 0);
        }
    }
    else{
        
        for ( auto Bbeg = ParticleMaterialProcessEnergySigmaMap[ParticleForCrossSection][MaterialForCrossSection][ProcessForCrossSection].begin(); Bbeg != ParticleMaterialProcessEnergySigmaMap[ParticleForCrossSection][MaterialForCrossSection][ProcessForCrossSection].end(); ++Bbeg  ){
            
            double Ene = Bbeg.key();
            double Sigma = Bbeg.value();
            
            Energy_Sigma[Ene] = Sigma ;
        }
        
        graphs_Title = ProcessForCrossSection+ " Sigma " + ParticleForCrossSection + " in " + MaterialForCrossSection + " Material";
        QString graph_legend = ProcessForCrossSection;
        create_graphs( Energy_Sigma , graph_legend, graphs_Title , 0);
    }
    
}
void PlotDialog::create_Data_for_SourceTimeSimulation_graphs(QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>> table, QString dataReferenceFileName){ // 1 for Self , 2 for cross
    
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << "() ------------------" << "\n";
    
    ui->customPlot->update();
    
    QString Source_ORG;
    
    
    // iterations on source name
    for ( auto Bbeg = ResultParticleSourceEnergyTime[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Bbeg != ResultParticleSourceEnergyTime[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Bbeg  ){
        
        Source_ORG = Bbeg.key();
        QTextStream(stdout) <<"Source_ORG : " << Source_ORG << "\n";
        
        QMap<double, double> Energy_ScoreVar;
        
        for(int ii = 0; ii < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size(); ii++){
            
            if( ResultParticleSourceEnergyTime[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]] == 0){
                Energy_ScoreVar[ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]] = MinValForLog;
                
            }else{
                Energy_ScoreVar[ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]] = ResultParticleSourceEnergyTime[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]] ;
            }
        }
        
        QTextStream(stdout) <<"Source_ORG " << Source_ORG << "\n";
        
        QTextStream(stdout) <<GraphsData << " " << Compare_Type << "\n";
        
        graphs_Title = "Simulation time by " +  ParticleName + " Irradiarion From Each source ";
        QString graph_legend = Source_ORG;
        create_graphs( Energy_ScoreVar , graph_legend, graphs_Title , 0);
        
    }
}
void PlotDialog::create_QuantitiesData_Cross_Self_FromRadioTracerIntake_graphs(){ // 1 for Self , 2 for cross

    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << "() ------------------" << "\n";

    if(!ui->checkBoxAppendGraph->isChecked()){
        ui->customPlot->clearPlottables();
        ui->customPlot->legend->clear();
        ui->customPlot->clearGraphs();
    }
    ui->customPlot->update();
    QMap<double, double> Energy_ScoreVar;
    QMap<double, double> RefEnergy_ScoreVar;
    QMap<double, double> Ref2Energy_ScoreVar;
    xticks.clear();xticks.empty();
    xlabels.clear();xlabels.empty();
    double iii = 1;

    //if(GraphsData == "Reference_Result"){
    //    ReferenceTable = fileManagerObjectPlot->Read_Comparison_file(ReferenceFilePath, 12);
    //    ReferenceQuantityGeometryRadioTracerSourceTargetValues = fileManagerObjectPlot->getReferenceQuantityGeometryRadioTracerSourceTargetValues();
    //}

    // iterations on source name
    for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][SourceOrgan].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][SourceOrgan].end(); ++Abeg  ){

        // for just Reference_Result we will show just the selected source organ

        QString Target_ORG = Abeg.key();

        if(GraphsData == "Reference_Result"){

            QTextStream(stdout) << "ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][SourceOrgan].size() : " << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][SourceOrgan].size() << "\n";
            double vv = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][SourceOrgan][Target_ORG];
            double dd = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][SourceOrgan][Target_ORG];
            if(vv != 0. || dd != 0.){
                //RefEnergy_ScoreVar[RegionParameterValueMap[GeometrySymbol][RegionVariableName][Target_ORG]] = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][SourceOrgan][Target_ORG];
                RefEnergy_ScoreVar[iii] = dd ;
                Energy_ScoreVar[iii] = vv ;
            }

            QTextStream(stdout) << " Result Reference RadioTracerName " << RadioTracerName << " SourceOrgan " << SourceOrgan << " Target_ORG " << Target_ORG << " " << RegionParameterValueMap[GeometrySymbol][RegionVariableName][Target_ORG] << "  "<< vv << "  "<< dd << "\n";

            if(CompareReferenceName_2 != "" && QFile::exists(ui->PlotLineEditReferenceFile_2->text())){
                //Ref2Energy_ScoreVar[RegionParameterValueMap[GeometrySymbol][RegionVariableName][Target_ORG]] = ReferenceTable_2[QuantitiesToScore][GeometrySymbol][ParticleName][SourceOrgan][Target_ORG][GraphForEnergy];
            }
        }
        else if (GraphsData == "Result"){
            double vv = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][SourceOrgan][Target_ORG];
            if(vv != 0.){
                Energy_ScoreVar[iii] = vv ;
            }
        }

        //xticks.push_back(RegionParameterValueMap[GeometrySymbol][RegionVariableName][Target_ORG]);
        xticks.push_back(iii);
        iii = iii + 1. ;
        xlabels.push_back(Target_ORG);
    }

    if(GraphsData == "Result"){
        graphs_Title = RadioTracerName + " " + QuantitiesToScore + " from "+ SourceOrgan+" to targets in " + GeometrySymbol +" phantom";
    }
    else{
        graphs_Title = RadioTracerName + " " + QuantitiesToScore + " from "+ SourceOrgan+" to targets in " + GeometrySymbol +" phantom Compared to " + CompareReferenceName;
    }

    QString graph_legend = SourceOrgan;
    if(GraphsData == "Reference_Result"){
        graph_legend = "DoseCalcs "+ SourceOrgan;
    }

    QTextStream(stdout) << "Energy_ScoreVar.size() : " << Energy_ScoreVar.size() << "\n";
    create_graphs( Energy_ScoreVar , graph_legend, graphs_Title , 0);

    if(GraphsData == "Reference_Result"){
        graph_legend = CompareReferenceName +" "+ SourceOrgan;
        create_graphs( RefEnergy_ScoreVar , graph_legend, graphs_Title , 0);
        if(CompareReferenceName_2 != "" && QFile::exists(ui->PlotLineEditReferenceFile_2->text())){
            //graph_legend = CompareReferenceName_2 +" "+ SourceOrgan;
            //create_graphs( RefEnergy_ScoreVar , graph_legend, graphs_Title , 0);
        }
    }
}
void PlotDialog::create_RatioData_For_ParticleAndRadiotracer_In_RatioPlot(){

    QTextStream(stdout) << " ----------------- 111------------------------------- \n";

    QMap<QString,QVector<double>> ResultsDataforRatioPlot ;
    QMap<QString,QMap<QString,double>> RationInYforXData ;

    // //////////////////////////////////// For particles

    double ratio;

    QString Xaxislabel = "x";
    ui->customPlot->yAxis->setLabel(DiffExp); // to set the new axis quantity

    if(RadioTraceBarPlotData == "f(x)=Sources, x=SourceTargets, for for Particle energies and Geometry"){ // for a specific geometry, point for each source-source combination and x for radiotracer
        Xaxislabel = "Source Region (Energy)";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " " + ParticleName + " for all energies and in " + GeometrySymbol ;
        for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();
                for ( auto Dbeg = Bbeg.value().begin(); Dbeg != Bbeg.value().end(); ++Dbeg  ){
                    double Ene = Dbeg.key();

                    ratio = RelativeDifferenceCalculation(ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][Second2][Ene],ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][Second2][Ene]);
                    if(ratio == NULL){continue;}

                    ResultsDataforRatioPlot[First1].push_back(ratio);
                    QString graphname = First1+"("+QString::number(Ene)+" MeV)";
                    RationInYforXData[Second2][graphname]=ratio;
                    QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " " << ratio << " size " << ResultsDataforRatioPlot[First1].size() << "\n";
                }
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Geometries, x=Sources, for Particle energies"){ // for a specific geometry, point for each source-source combination and x for radiotracer
        Xaxislabel = "Source<-Source";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " " + ParticleName + " for all energies and in all geometries";
        for ( auto Abeg = ResultTable[QuantitiesToScore].begin(); Abeg != ResultTable[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value()[ParticleName].begin(); Bbeg != Abeg.value()[ParticleName].end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();
                for ( auto Dbeg = Bbeg.value()[Second2].begin(); Dbeg != Bbeg.value()[Second2].end(); ++Dbeg  ){
                    double Ene = Dbeg.key();
                    QString graphname = Second2+"<-"+Second2+"("+QString::number(Ene)+" MeV)";

                    ratio = RelativeDifferenceCalculation(ResultTable[QuantitiesToScore][First1][ParticleName][Second2][Second2][Ene],ReferenceTable[QuantitiesToScore][First1][ParticleName][Second2][Second2][Ene]);
                    if(ratio == NULL){continue;}

                    ResultsDataforRatioPlot[First1].push_back(ratio);
                    RationInYforXData[First1][graphname]=ratio;
                    QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " " << ratio << "\n";
                }
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Geometries, x=SourceTargets, for Particle energies"){ // for a specific geometry, point for each source-source combination and x for radiotracer
        Xaxislabel = "Target<-Source (Energy (MeV))";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " " + ParticleName + " for all energies and in all geometries";
        OpenDialogOfCombinationChooser();
        for ( auto Abeg = ResultTable[QuantitiesToScore].begin(); Abeg != ResultTable[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for(int dd=0; dd < SourceOrganForCombination.size();dd++){
                QString Second2 = TargetOrganForCombination[dd];
                SourceOrgan = SourceOrganForCombination[dd];
                for ( auto Dbeg = ResultTable[QuantitiesToScore][First1][ParticleName][SourceOrgan][Second2].begin(); Dbeg != ResultTable[QuantitiesToScore][First1][ParticleName][SourceOrgan][Second2].end(); ++Dbeg  ){
                    double Ene = Dbeg.key();
                    QString graphname = TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd]+"("+QString::number(Ene)+" MeV)";

                    ratio = RelativeDifferenceCalculation(ResultTable[QuantitiesToScore][First1][ParticleName][SourceOrgan][Second2][Ene],ReferenceTable[QuantitiesToScore][First1][ParticleName][SourceOrgan][Second2][Ene]);
                    if(ratio == NULL){continue;}

                    ResultsDataforRatioPlot[First1].push_back(ratio);
                    RationInYforXData[First1][graphname]=ratio;
                    QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " " << ratio << "\n";
                }
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Geometries, x=SourceTargets, for Particle and Energy"){ // for a specific geometry, point for each source-source combination and x for radiotracer
        Xaxislabel = "Target<-Source";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " " + ParticleName + " ("+ui->plotComboBoxEnergies->currentText()+" MeV) and in all geometries";
        OpenDialogOfCombinationChooser();
        for ( auto Abeg = ResultTable[QuantitiesToScore].begin(); Abeg != ResultTable[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for(int dd=0; dd < SourceOrganForCombination.size();dd++){
                QString Second2 = TargetOrganForCombination[dd];
                SourceOrgan = SourceOrganForCombination[dd];
                double Ene = ui->plotComboBoxEnergies->currentText().toDouble();
                QString graphname = TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd];

                ratio = RelativeDifferenceCalculation(ResultTable[QuantitiesToScore][First1][ParticleName][SourceOrgan][Second2][Ene],ReferenceTable[QuantitiesToScore][First1][ParticleName][SourceOrgan][Second2][Ene]);
                if(ratio == NULL){continue;}

                ResultsDataforRatioPlot[First1].push_back(ratio);
                RationInYforXData[First1][graphname]=ratio;
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " " << ratio << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Energies, x=SourceTargets , for Geometry and Particle"){ // for a specific geometry, point for each source-source combination and x for radiotracer
        Xaxislabel = "Target<-Source";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " for all simulated " + ParticleName + " energies in " + GeometrySymbol ;
        OpenDialogOfCombinationChooser();
        for(int dd=0; dd < SourceOrganForCombination.size();dd++){
            QString First1 = SourceOrganForCombination[dd];
            QString Second2 = TargetOrganForCombination[dd];
            QString graphname = Second2+"<-"+First1;
            for ( auto Cbeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][Second2].begin(); Cbeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][Second2].end(); ++Cbeg  ){
                double Ene = Cbeg.key();

                ratio = RelativeDifferenceCalculation(ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][Second2][Ene],ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][Second2][Ene]);
                if(ratio == NULL){continue;}

                ResultsDataforRatioPlot[QString::number(Ene)].push_back(ratio);
                RationInYforXData[QString::number(Ene)][graphname]=ratio;
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " QString::number(Ene) " << QString::number(Ene) << " Ene " << Ene << " Ratio " << ratio << " Size of vector " << ResultsDataforRatioPlot[QString::number(Ene)].size() << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=SourceTargets, x=Energies , for Geometry and Particle"){ // for a specific geometry, point for each source-source combination and x for radiotracer
        Xaxislabel = "Energy (MeV)";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " for all simulated " + ParticleName + " energies in " + GeometrySymbol ;
        OpenDialogOfCombinationChooser();
        for(int dd=0; dd < SourceOrganForCombination.size();dd++){
            QString First1 = SourceOrganForCombination[dd];
            QString Second2 = TargetOrganForCombination[dd];
            QString graphname = Second2+"<-"+First1;
            for ( auto Cbeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][Second2].begin(); Cbeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][Second2].end(); ++Cbeg  ){
                double Ene = Cbeg.key();

                ratio = RelativeDifferenceCalculation(ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][Second2][Ene],ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][Second2][Ene]);
                if(ratio == NULL){continue;}

                ResultsDataforRatioPlot[QString::number(Ene)].push_back(ratio);
                RationInYforXData[graphname][QString::number(Ene)]=ratio;
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " QString::number(Ene) " << QString::number(Ene) << " Ene " << Ene << " Ratio " << ratio << " Size of vector " << ResultsDataforRatioPlot[QString::number(Ene)].size() << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Energies, x=Sources , for Geometry and Particle"){ // for a specific geometry, point for each source-source combination and x for radiotracer
        Xaxislabel = "Source Region";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " for all simulated " + ParticleName + " energies in " + GeometrySymbol ;
        for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Cbeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][First1].begin(); Cbeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][First1].end(); ++Cbeg  ){
                double Ene = Cbeg.key();

                ratio = RelativeDifferenceCalculation(ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][First1][Ene],ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][First1][Ene]);
                if(ratio == NULL){continue;}

                ResultsDataforRatioPlot[QString::number(Ene)].push_back(ratio);
                RationInYforXData[QString::number(Ene)+ " MeV"][First1]=ratio;
                QTextStream(stdout) << " First1 " << First1 <<  " QString::number(Ene) " << QString::number(Ene) << " Ene " << Ene << " Ratio " << ratio << " Size of vector " << ResultsDataforRatioPlot[QString::number(Ene)].size() << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Sources, x=Energies , for Geometry and Particle"){ // for a specific geometry, point for each source-source combination and x for radiotracer
        Xaxislabel = "Energy (MeV)";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " for all simulated " + ParticleName + " energies in " + GeometrySymbol ;
        for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Cbeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][First1].begin(); Cbeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][First1].end(); ++Cbeg  ){
                double Ene = Cbeg.key();

                ratio = RelativeDifferenceCalculation(ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][First1][Ene],ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][First1][First1][Ene]);
                if(ratio == NULL){continue;}

                ResultsDataforRatioPlot[QString::number(Ene)].push_back(ratio);
                RationInYforXData[First1][QString::number(Ene)]=ratio;
                QTextStream(stdout) << " First1 " << First1 << " QString::number(Ene) " << QString::number(Ene) << " Ene " << Ene << " Ratio " << ratio << " Size of vector " << ResultsDataforRatioPlot[QString::number(Ene)].size() << "\n";
            }
        }
    }

    // //////////////////////////////////////// for Radiotracer

    else if(RadioTraceBarPlotData == "f(x)=Sources, x=RadioTracers, for Geometry"){ // for a specific geometry, point for each source-source combination and x for radiotracer
        Xaxislabel = "Radiotracer";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " " + GeometrySymbol ;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            //QTextStream(stdout) << " First1 " << First1 << "\n";
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();

                ratio = RelativeDifferenceCalculation(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][First1][Second2][Second2],ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][First1][Second2][Second2]);
                if(ratio == NULL){continue;}

                ResultsDataforRatioPlot[First1].push_back(ratio);
                RationInYforXData[Second2][First1]=ratio;
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " " << ratio << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=SourceTargets, x=RadioTracers, for Geometry"){ // for a specific geometry, point for each source-target combination and x for radiotracer
        Xaxislabel = "Radiotracer";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " " + GeometrySymbol ;
        OpenDialogOfCombinationChooser();
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for(int dd=0; dd < SourceOrganForCombination.size();dd++){
                QString Second2 = SourceOrganForCombination[dd];
                QString Third3 = TargetOrganForCombination[dd];
                QString graphname= TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd];

                ratio = RelativeDifferenceCalculation(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][First1][Second2][Third3],ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][First1][Second2][Third3]);
                if(ratio == NULL){continue;}

                ResultsDataforRatioPlot[First1].push_back(ratio);
                RationInYforXData[graphname][First1]=ratio;
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " Third3 " << Third3 << " " << ratio << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Sources, x=Geometries, for RadioTracer"){ // for a specific radiotracer, point for each source-source combination and x for radiotracer
        Xaxislabel = "Phantom Geometry";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " " + RadioTracerName ;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value()[RadioTracerName].begin(); Bbeg != Abeg.value()[RadioTracerName].end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();

                ratio = RelativeDifferenceCalculation(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][RadioTracerName][Second2][Second2],ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][RadioTracerName][Second2][Second2]);
                if(ratio == NULL){continue;}

                ResultsDataforRatioPlot[First1].push_back(ratio);
                RationInYforXData[Second2][First1]=ratio;
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " " << ratio << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=SourceTargets, x=Geometries, for RadioTracer"){ // for a specific radiotracer, point for each source-target combination and x for geometry
        Xaxislabel = "Phantom Geometry";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " " + RadioTracerName ;
        OpenDialogOfCombinationChooser();
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for(int dd=0; dd < SourceOrganForCombination.size();dd++){
                QString Second2 = SourceOrganForCombination[dd];
                QString Third3 = TargetOrganForCombination[dd];
                QString graphname= TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd];

                ratio = RelativeDifferenceCalculation(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][RadioTracerName][Second2][Third3],ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][RadioTracerName][Second2][Third3]);
                if(ratio == NULL){continue;}

                ResultsDataforRatioPlot[First1].push_back(ratio);
                RationInYforXData[graphname][First1]=ratio;
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " Third3 " << Third3 << " " << ratio << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Geometries, x=SourceTargets, for RadioTracer"){ // for a specific radiotracer, point for each source-target combination and x for geometry
        Xaxislabel = "Target<-Source";
        graphs_Title = QuantitiesToScore + " Ratios for " + RadioTraceBarPlotData + " " + RadioTracerName ;
        OpenDialogOfCombinationChooser();
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString graphname = Abeg.key();
            for(int dd=0; dd < SourceOrganForCombination.size();dd++){
                QString Second2 = SourceOrganForCombination[dd];
                QString Third3 = TargetOrganForCombination[dd];
                QString First1= TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd];

                ratio = RelativeDifferenceCalculation(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][graphname][RadioTracerName][Second2][Third3],ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][graphname][RadioTracerName][Second2][Third3]);
                if(ratio == NULL){continue;}

                ResultsDataforRatioPlot[First1].push_back(ratio);
                RationInYforXData[graphname][First1]=ratio;
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " Third3 " << Third3 << " " << ratio << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=RadioTracers, x=Geometries, for Source"){ // for a specific sources, point for each Radiotracer and x for geometries
        Xaxislabel = "Phantom Geometry";
        graphs_Title = QuantitiesToScore + " Ratios " + RadioTraceBarPlotData + " " + SourceOrgan +"->"+SourceOrgan;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();

                ratio = RelativeDifferenceCalculation(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][SourceOrgan],ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][SourceOrgan]);
                if(ratio == NULL){continue;}

                ResultsDataforRatioPlot[First1].push_back(ratio);
                RationInYforXData[Second2][First1]=ratio;
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " " << ratio << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=RadioTracers, x=Geometries, for SourceTarget"){ // for a specific source-targets, point for each Radiotracer and x for geometries

        Xaxislabel = "Phantom Geometry";
        graphs_Title = QuantitiesToScore + " Ratios " + RadioTraceBarPlotData + " " + SourceOrgan +"->"+ TargetOrgan;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();

                ratio = RelativeDifferenceCalculation(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][TargetOrgan],ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][TargetOrgan]);
                if(ratio == NULL){continue;}

                ResultsDataforRatioPlot[First1].push_back(ratio);
                RationInYforXData[Second2][First1]=ratio;
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << " " << ratio << "\n";
            }
        }
    }

    QTextStream(stdout) << " graphs_Title " << graphs_Title << " data size " << RationInYforXData.size() << " " << ResultsDataforRatioPlot.size() << "\n";

    if(GraphsRadio->isChecked()){

        QTextStream(stdout) << " ----------------- 222 ------------------------------- \n";

        if(RationInYforXData.size() == 0 ){
            QMessageBox::information(this, tr(""), "No ratio data for " + QuantitiesToScore +", geometry, radiotracer (or particle), source and/or target are registered for plot, check the reference file");
            return;
        }

        if(!ui->checkBoxAppendGraph->isChecked()){
            ui->customPlot->clearPlottables();
            ui->customPlot->legend->clear();
            ui->customPlot->clearGraphs();
        }
        ui->customPlot->update();

        xticks.clear();
        xlabels.clear();

        double iii = 1. ;

        ui->customPlot->xAxis->setLabel(Xaxislabel);

        for ( auto Abeg = RationInYforXData.begin(); Abeg != RationInYforXData.end(); ++Abeg  ){
            QString graphname = Abeg.key();

            QMap<double, double> XYMap;

            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString label = Bbeg.key();

                bool isIn = false;
                for ( int df = 0 ; df < xlabels.size(); df++  ){
                    if(xlabels[df] == label ){ // if already exist for older data
                        iii = xticks[df];
                        isIn = true;
                    }
                }

                XYMap[iii] = Bbeg.value();

                if(isIn == false){
                    xlabels.push_back(label);
                    xticks.push_back(iii);
                    iii = iii +1. ;
                }
            }
            create_graphs( XYMap, graphname, graphs_Title , 0);
        }
    }
    else{

        QTextStream(stdout) << " ----------------- 222 ------------------------------- \n";

        if(ResultsDataforRatioPlot.size() == 0){
            QMessageBox::information(this, tr(""), "No ratio data for " + QuantitiesToScore +", geometry, radiotracer (or particle), source and/or target are registered for plot, check the reference file");
            return;
        }

        QCPStatisticalBox *statistical = new QCPStatisticalBox(ui->customPlot->xAxis, ui->customPlot->yAxis);
        QBrush boxBrush(QColor(60, 60, 255, 100));
        boxBrush.setStyle(Qt::Dense6Pattern); // make it look oldschool
        statistical->setBrush(boxBrush);

        QString axislabel = "Ratio DoseCalcs/"+CompareReferenceName;
        ui->customPlot->yAxis->setLabel(axislabel);

        // prepare manual x axis labels:

        ui->customPlot->xAxis->setScaleType(QCPAxis::stLinear);
        ui->customPlot->xAxis->setSubTicks(false);
        ui->customPlot->xAxis->setTickLength(0, 4);
        ui->customPlot->xAxis->setTickLabelRotation(20);
        QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
        //ui->customPlot->xAxis->setScaleType(QCPAxis::stLinear);

        int inc = 1;
        for ( auto Abeg = ResultsDataforRatioPlot.begin(); Abeg != ResultsDataforRatioPlot.end(); ++Abeg  ){

            QString First1 = Abeg.key();

            QTextStream(stdout) << " First1 " << First1 <<  " vector of size " << ResultsDataforRatioPlot[First1].size() << "\n";
            double min, max, median, lowerQt, upperQt;
            int mid , lower, upper;

            if(ResultsDataforRatioPlot[First1].size() == 1){
                min = max = median = lowerQt = upperQt = ResultsDataforRatioPlot[First1][0];
            }else{

                std::sort(ResultsDataforRatioPlot[First1].begin(), ResultsDataforRatioPlot[First1].end());
                min = *std::min_element(ResultsDataforRatioPlot[First1].constBegin(), ResultsDataforRatioPlot[First1].constEnd());
                max = *std::max_element(ResultsDataforRatioPlot[First1].constBegin(), ResultsDataforRatioPlot[First1].constEnd());

                mid; if(ResultsDataforRatioPlot[First1].size()%2 == 0){ mid = ResultsDataforRatioPlot[First1].size()/2;}else{mid = (ResultsDataforRatioPlot[First1].size()/2);  if(mid >= ResultsDataforRatioPlot[First1].size()){mid = ResultsDataforRatioPlot[First1].size() - 1;}}
                lower; if(mid%2 == 0){ lower = mid/2;}else{lower = (mid/2);if(lower >= ResultsDataforRatioPlot[First1].size()){lower = ResultsDataforRatioPlot[First1].size() - 1;}}
                upper; if(mid%2 == 0){ upper = ResultsDataforRatioPlot[First1].size() - (mid/2);}else{upper = ResultsDataforRatioPlot[First1].size() - ((mid/2)); if(upper >= ResultsDataforRatioPlot[First1].size()){upper = ResultsDataforRatioPlot[First1].size() - 1;}}

                QTextStream(stdout) << ResultsDataforRatioPlot[First1].size() << " " << lower<< " " << mid<< " " << upper << " " << " \n";

                //if(mid == 0 || lower == 0 || upper == 0){continue;}

                if(mid > ResultsDataforRatioPlot[First1].size()){mid = mid - 1;} median = ResultsDataforRatioPlot[First1][mid];
                if(lower > ResultsDataforRatioPlot[First1].size()){lower = lower - 1;} lowerQt =ResultsDataforRatioPlot[First1][lower];
                if(upper > ResultsDataforRatioPlot[First1].size()){upper = upper - 1;} upperQt =ResultsDataforRatioPlot[First1][upper];
            }

            QTextStream(stdout) << min << " " << max << " " << median << " " << lowerQt << " " << upperQt << " " << " \n";

            //if(RadioTraceBarPlotData == "f(x)=Energies, x=SourceTargets , for Geometry and Particle"){ // for a specific geometry, point for each source-source combination and x for radiotracer
            //statistical->addData(First1.toDouble(), min, lowerQt, median, upperQt, max, ResultsDataforRatioPlot[First1]);
            //textTicker->addTick(First1.toDouble(), First1);
            //}

            statistical->addData(inc, min, lowerQt, median, upperQt, max, ResultsDataforRatioPlot[First1]);
            textTicker->addTick(inc, First1);

            QTextStream(stdout) << " First1 " << First1 << " First1.toDouble() " << First1.toDouble() << " " << " \n";

            inc++;

        }

        QTextStream(stdout) << " ----------------- 333 ------------------------------- \n";

        ui->customPlot->plotLayout()->remove(title);
        title = new QCPTextElement(ui->customPlot, graphs_Title , QFont("sans", 17, QFont::Bold));
        ui->customPlot->plotLayout()->addElement(0, 0, title);
        ui->customPlot->plotLayout()->updateLayout();
        connect(title, SIGNAL(doubleClicked(QMouseEvent*)), this, SLOT(titleDoubleClick(QMouseEvent*)));

        ui->customPlot->xAxis->setTicker(textTicker);

        QSharedPointer<QCPAxisTicker> FixedTicker(new QCPAxisTicker);
        ui->customPlot->yAxis->setScaleType(QCPAxis::stLinear);
        ui->customPlot->yAxis->setTicker(FixedTicker);
        ui->customPlot->yAxis->setNumberFormat(YNumberFormat); // e = exponential, b = beautiful decimal powers
        ui->customPlot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"

        // prepare axes:
        ui->customPlot->rescaleAxes();
        ui->customPlot->xAxis->scaleRange(1.7, ui->customPlot->xAxis->range().center());
        ui->customPlot->yAxis->setRange(0, 7);
        ui->customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
        //ui->customPlot->yAxis->setScaleType(QCPAxis::stLinear);

        ui->customPlot->legend->setVisible(false);

        ui->customPlot->replot();

        QTextStream(stdout) << " ----------------- 444 ------------------------------- \n";
    }
}
void PlotDialog::create_QuantitiesData_For_Radiotracer_In_Plot_Bars(){

    QVector<QColor> mycolors;
    //QStringList cl = QColor::colorNames();
    //for(int i=0; i < cl.size(); i++){
    //  mycolors.push_back(QColor(cl[i]));
    //QTextStream(stdout) << " Color Name " << cl[i] << "\n";
    //}

    mycolors << "red" << "cyan"<< "purple"<< "green"<< "yellow"<< "blue"<< "orange"
             << "brown" << "pink"<< "plum"<< "magenta"<< "olive"<< "gray"<< "orchid"
             << "indigo"  << "peru" << "khaki"<< "gold"<< "aqua"<< "chocolate"<< "black";

    QMap<QString,QMap<QString,double>> ResultsDataforBarPlot ;

    QTextStream(stdout) << " ------------------------------------------------ \n";

    QString axislabel = ScoreVariableUnit;

    if(RadioTraceBarPlotData == "f(x)=RadioTracers, x=Sources, for Geometry"){
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " " + GeometrySymbol ;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();
                ResultsDataforBarPlot[First1][Second2] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][First1][Second2][Second2];
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Sources, x=RadioTracers, for Geometry"){
        graphs_Title = QuantitiesToScore + " " +RadioTraceBarPlotData + " " + GeometrySymbol ;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();
                ResultsDataforBarPlot[Second2][First1] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][First1][Second2][Second2];
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Geometries, x=Sources, for RadioTracer"){
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " " + RadioTracerName ;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value()[RadioTracerName].begin(); Bbeg != Abeg.value()[RadioTracerName].end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();
                ResultsDataforBarPlot[First1][Second2] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][RadioTracerName][Second2][Second2];
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Sources, x=Geometries, for RadioTracer"){
        graphs_Title = QuantitiesToScore + " " +RadioTraceBarPlotData + " " + RadioTracerName ;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value()[RadioTracerName].begin(); Bbeg != Abeg.value()[RadioTracerName].end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();
                ResultsDataforBarPlot[Second2][First1] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][RadioTracerName][Second2][Second2];
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=SourceTargets, x=RadioTracers, for Geometry"){ // for a quantity, radiotracer (or particle), geometry.
        OpenDialogOfCombinationChooser();
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " " + GeometrySymbol ;

        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for(int dd=0; dd < SourceOrganForCombination.size();dd++){

                QString labelX= TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd];
                ResultsDataforBarPlot[labelX][First1] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][First1][SourceOrganForCombination[dd]][TargetOrganForCombination[dd]];
                QTextStream(stdout) << " First1 " << First1 <<  " labelX " << labelX << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=RadioTracers, x=SourceTargets, for Geometry"){ // for a quantity, radiotracer (or particle), geometry.
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " " + GeometrySymbol ;
        OpenDialogOfCombinationChooser();

        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for(int dd=0; dd < SourceOrganForCombination.size();dd++){

                QString labelX= TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd];
                ResultsDataforBarPlot[First1][labelX] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][First1][SourceOrganForCombination[dd]][TargetOrganForCombination[dd]];
                QTextStream(stdout) << " First1 " << First1 <<  " labelX " << labelX << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=SourceTargets, x=Geometries, for RadioTracer"){ // for a quantity, radiotracer (or particle), geometry.
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " " + RadioTracerName ;
        OpenDialogOfCombinationChooser();
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){

            QString First1 = Abeg.key();
            for(int dd=0; dd < SourceOrganForCombination.size();dd++){

                QString labelX= TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd];
                ResultsDataforBarPlot[labelX][First1] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][RadioTracerName][SourceOrganForCombination[dd]][TargetOrganForCombination[dd]];
                QTextStream(stdout) << " First1 " << First1 <<  " labelX " << labelX << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Geometries, x=SourceTargets, for RadioTracer"){ // for a quantity, radiotracer (or particle), geometry.
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " " + RadioTracerName ;
        OpenDialogOfCombinationChooser();
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){

            QString First1 = Abeg.key();
            for(int dd=0; dd < SourceOrganForCombination.size();dd++){

                QString labelX= TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd];
                ResultsDataforBarPlot[First1][labelX] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][RadioTracerName][SourceOrganForCombination[dd]][TargetOrganForCombination[dd]];
                QTextStream(stdout) << " First1 " << First1 <<  " labelX " << labelX << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Geometries, x=RadioTracers, for Sources"){ // for a quantity, radiotracer (or particle), geometry.
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " for (self) source " + SourceOrgan ;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();
                ResultsDataforBarPlot[First1][Second2] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][SourceOrgan];
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=RadioTracers, x=Geometries, for Sources"){ // for a quantity, radiotracer (or particle), geometry.
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " for (self) source " + SourceOrgan ;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();
                ResultsDataforBarPlot[Second2][First1] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][SourceOrgan];
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Geometries, x=RadioTracers, for SourceTarget"){ // for a quantity, radiotracer (or particle), geometry.
        graphs_Title = QuantitiesToScore + " " +RadioTraceBarPlotData + " " + SourceOrgan+"->"+TargetOrgan;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();
                ResultsDataforBarPlot[First1][Second2] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][TargetOrgan];
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=RadioTracers, x=Geometries, for SourceTarget"){ // for a quantity, radiotracer (or particle), geometry.
        graphs_Title = QuantitiesToScore + " " +RadioTraceBarPlotData + " " + SourceOrgan+"->"+TargetOrgan;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key();
                ResultsDataforBarPlot[Second2][First1] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][TargetOrgan];
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";
            }
        }
    }

    if(ResultsDataforBarPlot.size() == 0){
        QMessageBox::information(this, tr(""), "No data for " + QuantitiesToScore +", geometry, radiotracer (or particle), source and/or target are registered for plot, check the result file");
        return;
    }

    ui->customPlot->plotLayout()->remove(title);
    title = new QCPTextElement(ui->customPlot, graphs_Title , QFont("sans", 17, QFont::Bold));
    ui->customPlot->plotLayout()->addElement(0, 0, title);
    ui->customPlot->plotLayout()->updateLayout();
    connect(title, SIGNAL(doubleClicked(QMouseEvent*)), this, SLOT(titleDoubleClick(QMouseEvent*)));

    QLinearGradient gradient(0, 0, 0, 0);
    gradient.setColorAt(0, QColor(255, 255, 255));
    gradient.setColorAt(0.38, QColor(255, 255, 255));
    gradient.setColorAt(1, QColor(255, 255, 255));
    ui->customPlot->setBackground(QBrush(gradient));

    QVector<QCPBars*> Bars;
    QVector<double>  ticks;
    QVector<QString> labels;

    QTextStream(stdout) << " ------------------------------------------------ \n";

    int incr = 1;
    for ( auto Abeg = ResultsDataforBarPlot.begin(); Abeg != ResultsDataforBarPlot.end(); ++Abeg  ){

        QString First1 = Abeg.key();
        for ( auto Bbeg = ResultsDataforBarPlot[First1].begin(); Bbeg != ResultsDataforBarPlot[First1].end(); ++Bbeg  ){
            QString Second2 = Bbeg.key();
            bool isin = false;
            for(int i=0; i<labels.size(); i++){
                if(labels[i] == Second2){
                    isin = true; break;
                }
            }
            if(isin == false){
                ticks << incr; incr ++;
                labels << Second2;
            }
        }
    }

    ui->customPlot->xAxis->setScaleType(QCPAxis::stLinear);
    //ui->customPlot->yAxis->setScaleType(QCPAxis::stLinear);

    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
    textTicker->addTicks(ticks, labels);
    ui->customPlot->xAxis->setTicker(textTicker);
    ui->customPlot->xAxis->setTickLabelRotation(60);
    ui->customPlot->xAxis->setSubTicks(false);
    ui->customPlot->xAxis->setTickLength(0, 4);
    ui->customPlot->xAxis->setRange(0, 8);
    ui->customPlot->xAxis->setBasePen(QPen(Qt::black));
    ui->customPlot->xAxis->setTickPen(QPen(Qt::black));
    ui->customPlot->xAxis->grid()->setVisible(true);
    ui->customPlot->xAxis->grid()->setPen(QPen(QColor(130, 130, 130), 0, Qt::DotLine));
    ui->customPlot->xAxis->setTickLabelColor(Qt::black);
    ui->customPlot->xAxis->setLabelColor(Qt::black);

    // prepare y axis:

    QSharedPointer<QCPAxisTicker> FixedTicker(new QCPAxisTicker);
    ui->customPlot->yAxis->setScaleType(QCPAxis::stLinear);
    ui->customPlot->yAxis->setTicker(FixedTicker);
    ui->customPlot->yAxis->setNumberFormat(YNumberFormat); // e = exponential, b = beautiful decimal powers
    ui->customPlot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
    ui->customPlot->yAxis->setLabel(axislabel);

    ui->customPlot->yAxis->setRange(0, 12.1);
    ui->customPlot->yAxis->setPadding(5); // a bit more space to the left border
    ui->customPlot->yAxis->setBasePen(QPen(Qt::black));
    ui->customPlot->yAxis->setTickPen(QPen(Qt::black));
    ui->customPlot->yAxis->setSubTickPen(QPen(Qt::black));
    ui->customPlot->yAxis->grid()->setSubGridVisible(true);
    ui->customPlot->yAxis->setTickLabelColor(Qt::black);
    ui->customPlot->yAxis->setLabelColor(Qt::black);
    ui->customPlot->yAxis->grid()->setPen(QPen(QColor(130, 130, 130), 0, Qt::SolidLine));
    ui->customPlot->yAxis->grid()->setSubGridPen(QPen(QColor(130, 130, 130), 0, Qt::DotLine));

    QTextStream(stdout) << " ------------------------------------------------ \n";

    int inclr = 0;
    for ( auto Abeg = ResultsDataforBarPlot.begin(); Abeg != ResultsDataforBarPlot.end(); ++Abeg  ){

        QString First1 = Abeg.key();
        QCPBars *abar = new QCPBars(ui->customPlot->xAxis, ui->customPlot->yAxis);
        abar->setAntialiased(false); // gives more crisp, pixel aligned bar borders
        abar->setStackingGap(1);
        // set names and colors:
        abar->setName(First1);

        int clrID = qrand() % (mycolors.size()- 1); //int(rand()*(mycolors.size()-1));
        if(inclr == mycolors.size()){inclr = 0;}
        //abar->setPen(QPen(mycolors[inclr].lighter(170)));
        abar->setBrush(mycolors[inclr]);
        inclr++;

        QVector<double> data;
        for(int i=0; i<labels.size(); i++){

            data.push_back(ResultsDataforBarPlot[First1][labels[i]]);
            //QTextStream(stdout) << " -- ColorID " << clrID << " -- "<<  QuantitiesToScore << " -- "<< First1 << " -- "<< labels[i] << " -- " << ResultsDataforBarPlot[First1][labels[i]] << "\n";

            //if(RadioTraceBarPlotData == "RadioTracer_Source"){
            //              }
            //else if(RadioTraceBarPlotData == "Source_RadioTracer"){
            //  data.push_back(ResultsDataforBarPlot[First1][labels[i]]);
            //QTextStream(stdout) << " -- ColorID " << clrID << " -- "<<  QuantitiesToScore << " -- "<< First1 << " -- "<< labels[i] << " -- " << ResultsDataforBarPlot[First1][labels[i]] << "\n";
            //}
        }
        abar->setData(ticks,data);

        Bars.push_back(abar);
    }

    QCPBars * recBar;
    for(int i=0; i<Bars.size(); i++){
        if(i != 0){ Bars[i]->moveAbove(recBar);}
        recBar = Bars[i];
    }

    // setup legend:
    ui->customPlot->legend->setVisible(true);
    ui->customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignHCenter);
    ui->customPlot->legend->setBrush(QColor(255, 255, 255, 100));
    //ui->customPlot->legend->setBorderPen(Qt::NoPen);
    QFont legendFont = font();
    legendFont.setPointSize(10);
    ui->customPlot->legend->setFont(legendFont);
    ui->customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);

    ui->customPlot->replot();

}
void PlotDialog::create_QuantitiesData_For_Radiotracer_In_Plot_Graphs(){

    if(!ui->checkBoxAppendGraph->isChecked()){
        ui->customPlot->clearPlottables();
        ui->customPlot->legend->clear();
        ui->customPlot->clearGraphs();
    }
    ui->customPlot->update();
    xticks.clear();xticks.empty();
    xlabels.clear();xlabels.empty();
    double iii = 1.;

    QMap<QString,QMap<double,double>> GraphNameXYValues;
    QVector<QString> GraphsNames ;

    QString Xaxislabel = "x";
    ui->customPlot->yAxis->setLabel(DiffExp); // to set the new axis quantity

    QTextStream(stdout) << " ------------------------------------------------ \n";

    if(RadioTraceBarPlotData      == "f(x)=RadioTracers, x=Sources, for Geometry"){
        Xaxislabel = "Source Region";
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " " + GeometrySymbol ;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  ){
            QString First1 = Abeg.key();

            iii = 1.;

            for(int dd=0; dd < SourceOrganlist.size();dd++){

                QString labelX= SourceOrganlist[dd];
                GraphNameXYValues[First1][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][First1][labelX][labelX];

                bool isIn = false;
                for ( int df = 0 ; df < xticks.size(); df++  ){if(xticks[df] == iii ){ isIn = true; } }
                if(isIn == false){ xlabels.push_back(labelX); xticks.push_back(iii);}

                isIn = false;
                for ( int df = 0 ; df < GraphsNames.size(); df++  ){if(GraphsNames[df] == First1 ){ isIn = true; } }
                if(isIn == false){GraphsNames.push_back(First1);}

                iii = iii + 1.;
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Sources, x=RadioTracers, for Geometry"){
        Xaxislabel = "Radiotracer";
        graphs_Title = QuantitiesToScore + " " +RadioTraceBarPlotData + " " + GeometrySymbol ;
        for(int dd=0; dd < SourceOrganlist.size();dd++){

            QString First1 = SourceOrganlist[dd];

            iii = 1.;

            for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  ){

                QString labelX = Abeg.key();
                GraphNameXYValues[First1][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][labelX][First1][First1];

                bool isIn = false;
                for ( int df = 0 ; df < xticks.size(); df++  ){if(xticks[df] == iii ){ isIn = true; } }
                if(isIn == false){ xlabels.push_back(labelX); xticks.push_back(iii);}

                isIn = false;
                for ( int df = 0 ; df < GraphsNames.size(); df++  ){if(GraphsNames[df] == First1 ){ isIn = true; } }
                if(isIn == false){GraphsNames.push_back(First1);}

                iii = iii + 1.;
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Geometries, x=Sources, for RadioTracer"){
        Xaxislabel = "Source Region";
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " " + RadioTracerName ;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();

            iii = 1.;

            for(int dd=0; dd < SourceOrganlist.size();dd++){

                QString labelX= SourceOrganlist[dd];
                GraphNameXYValues[First1][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][RadioTracerName][labelX][labelX];

                bool isIn = false;
                for ( int df = 0 ; df < xticks.size(); df++  ){if(xticks[df] == iii ){ isIn = true; } }
                if(isIn == false){ xlabels.push_back(labelX); xticks.push_back(iii);}

                isIn = false;
                for ( int df = 0 ; df < GraphsNames.size(); df++  ){if(GraphsNames[df] == First1 ){ isIn = true; } }
                if(isIn == false){GraphsNames.push_back(First1);}

                iii = iii + 1.;
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Sources, x=Geometries, for RadioTracer"){
        Xaxislabel = "Phantom Geometry";
        graphs_Title = QuantitiesToScore + " " +RadioTraceBarPlotData + " " + RadioTracerName ;
        for(int dd=0; dd < SourceOrganlist.size();dd++){

            QString First1 = SourceOrganlist[dd];

            iii = 1.;

            for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){

                QString labelX = Abeg.key();
                GraphNameXYValues[First1][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][labelX][RadioTracerName][First1][First1];

                bool isIn = false;
                for ( int df = 0 ; df < xticks.size(); df++  ){if(xticks[df] == iii ){ isIn = true; } }
                if(isIn == false){ xlabels.push_back(labelX); xticks.push_back(iii);}

                isIn = false;
                for ( int df = 0 ; df < GraphsNames.size(); df++  ){if(GraphsNames[df] == First1 ){ isIn = true; } }
                if(isIn == false){GraphsNames.push_back(First1);}

                iii = iii + 1.;
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Geometries, x=RadioTracers, for Sources"){ // for a quantity, radiotracer (or particle), geometry.
        Xaxislabel = "Radiotracer";
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " for (self) source " + SourceOrgan ;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key(); // geometry
            iii = 1.;
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key(); // radiotracer
                GraphNameXYValues[First1][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][SourceOrgan];
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";

                QString labelX= Second2;

                bool isIn = false;
                for ( int df = 0 ; df < xticks.size(); df++  ){if(xticks[df] == iii ){ isIn = true; } }
                if(isIn == false){ xlabels.push_back(labelX); xticks.push_back(iii);}

                isIn = false;
                for ( int df = 0 ; df < GraphsNames.size(); df++  ){if(GraphsNames[df] == First1 ){ isIn = true; } }
                if(isIn == false){GraphsNames.push_back(First1);}

                iii = iii + 1.;
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=RadioTracers, x=Geometries, for Sources"){ // for a quantity, radiotracer (or particle), geometry.
        Xaxislabel = "Phantom Geometry";
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " for (self) source " + SourceOrgan ;
        iii = 1.;

        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key(); // geometry

            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key(); // radiotracer
                GraphNameXYValues[Second2][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][SourceOrgan];
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";

                QString labelX= First1;

                bool isIn = false;
                for ( int df = 0 ; df < xticks.size(); df++  ){if(xticks[df] == iii ){ isIn = true; } }
                if(isIn == false){ xlabels.push_back(labelX); xticks.push_back(iii);}

                isIn = false;
                for ( int df = 0 ; df < GraphsNames.size(); df++  ){if(GraphsNames[df] == Second2 ){ isIn = true; } }
                if(isIn == false){GraphsNames.push_back(Second2);}
            }
            iii = iii + 1.;
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Geometries, x=RadioTracers, for SourceTarget"){ // for a quantity, radiotracer (or particle), geometry.
        Xaxislabel = "Radiotracer";
        graphs_Title = QuantitiesToScore + " " +RadioTraceBarPlotData + " " + SourceOrgan+"->"+TargetOrgan;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key(); // geometry
            iii = 1.;
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key(); // radiotracer
                GraphNameXYValues[First1][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][TargetOrgan];
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";

                QString labelX= Second2;

                bool isIn = false;
                for ( int df = 0 ; df < xticks.size(); df++  ){if(xticks[df] == iii ){ isIn = true; } }
                if(isIn == false){ xlabels.push_back(labelX); xticks.push_back(iii);}

                isIn = false;
                for ( int df = 0 ; df < GraphsNames.size(); df++  ){if(GraphsNames[df] == First1 ){ isIn = true; } }
                if(isIn == false){GraphsNames.push_back(First1);}

                iii = iii + 1.;
            }

            //for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
            //  QString Second2 = Bbeg.key();
            //ResultsDataforBarPlot[First1][Second2] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][TargetOrgan];
            //QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";
            //}
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=RadioTracers, x=Geometries, for SourceTarget"){ // for a quantity, radiotracer (or particle), geometry.
        Xaxislabel = "Phantom Geometry";
        graphs_Title = QuantitiesToScore + " " +RadioTraceBarPlotData + " " + SourceOrgan+"->"+TargetOrgan;

        iii = 1.;

        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key(); // geometry

            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                QString Second2 = Bbeg.key(); // radiotracer
                GraphNameXYValues[Second2][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][Second2][SourceOrgan][TargetOrgan];
                QTextStream(stdout) << " First1 " << First1 <<  " Second2 " << Second2 << "\n";

                QString labelX= First1;

                bool isIn = false;
                for ( int df = 0 ; df < xticks.size(); df++  ){if(xticks[df] == iii ){ isIn = true; } }
                if(isIn == false){ xlabels.push_back(labelX); xticks.push_back(iii);}

                isIn = false;
                for ( int df = 0 ; df < GraphsNames.size(); df++  ){if(GraphsNames[df] == Second2 ){ isIn = true; } }
                if(isIn == false){GraphsNames.push_back(Second2);}
            }
            iii = iii + 1.;
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=Geometries, x=SourceTargets, for RadioTracer"){ // for a quantity, radiotracer (or particle), geometry.
        OpenDialogOfCombinationChooser();
        Xaxislabel = "Target<-Source";
        graphs_Title = QuantitiesToScore + " " +RadioTraceBarPlotData + " for combinations for radiotracer" + RadioTracerName;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            iii = 1.;
            for(int dd=0; dd < SourceOrganForCombination.size();dd++){

                GraphNameXYValues[First1][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][First1][RadioTracerName][SourceOrganForCombination[dd]][TargetOrganForCombination[dd]];

                QString labelX= TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd];
                bool isIn = false;
                for ( int df = 0 ; df < xticks.size(); df++  ){if(xticks[df] == iii ){ isIn = true; } }
                if(isIn == false){ xlabels.push_back(labelX); xticks.push_back(iii);}

                isIn = false;
                for ( int df = 0 ; df < GraphsNames.size(); df++  ){if(GraphsNames[df] == First1 ){ isIn = true; } }
                if(isIn == false){GraphsNames.push_back(First1);}

                iii = iii + 1.;

                //    QTextStream(stdout) << " First1 " << First1 <<  " labelX " << labelX << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=SourceTargets, x=Geometries, for RadioTracer"){ // for a quantity, radiotracer (or particle), geometry.
        OpenDialogOfCombinationChooser();
        Xaxislabel = "Phantom Geometry";
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " " + RadioTracerName ;

        for(int dd=0; dd < SourceOrganForCombination.size();dd++){

            iii = 1.;
            for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].end(); ++Abeg  ){

                QString labelX = Abeg.key();
                QString First1= TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd];

                GraphNameXYValues[First1][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][labelX][RadioTracerName][SourceOrganForCombination[dd]][TargetOrganForCombination[dd]];

                bool isIn = false;
                for ( int df = 0 ; df < xticks.size(); df++  ){if(xticks[df] == iii ){ isIn = true; } }
                if(isIn == false){ xlabels.push_back(labelX); xticks.push_back(iii);}

                isIn = false;
                for ( int df = 0 ; df < GraphsNames.size(); df++  ){if(GraphsNames[df] == First1 ){ isIn = true; } }
                if(isIn == false){GraphsNames.push_back(First1);}

                //    QTextStream(stdout) << " First1 " << First1 <<  " labelX " << labelX << "\n";
                iii = iii + 1.;
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=RadioTracers, x=SourceTargets, for Geometry"){ // for a quantity, radiotracer (or particle), geometry.
        OpenDialogOfCombinationChooser();
        Xaxislabel = "Target<-Source";
        graphs_Title = QuantitiesToScore + " " +RadioTraceBarPlotData + " for combinations in geometry" + GeometrySymbol;
        for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  ){
            QString First1 = Abeg.key();
            iii = 1.;
            for(int dd=0; dd < SourceOrganForCombination.size();dd++){

                GraphNameXYValues[First1][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][First1][SourceOrganForCombination[dd]][TargetOrganForCombination[dd]];

                QString labelX= TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd];
                bool isIn = false;
                for ( int df = 0 ; df < xticks.size(); df++  ){if(xticks[df] == iii ){ isIn = true; } }
                if(isIn == false){ xlabels.push_back(labelX); xticks.push_back(iii);}

                isIn = false;
                for ( int df = 0 ; df < GraphsNames.size(); df++  ){if(GraphsNames[df] == First1 ){ isIn = true; } }
                if(isIn == false){GraphsNames.push_back(First1);}

                iii = iii + 1.;

                //              QTextStream(stdout) << " First1 " << First1 <<  " labelX " << labelX << "\n";
            }
        }
    }
    else if(RadioTraceBarPlotData == "f(x)=SourceTargets, x=RadioTracers, for Geometry"){ // for a quantity, radiotracer (or particle), geometry.
        OpenDialogOfCombinationChooser();
        Xaxislabel = "Radiotracer";
        graphs_Title = QuantitiesToScore + " " + RadioTraceBarPlotData + " " + GeometrySymbol ;

        for(int dd=0; dd < SourceOrganForCombination.size();dd++){

            iii = 1.;
            for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  ){

                QString labelX = Abeg.key();
                QString First1= TargetOrganForCombination[dd]+"<-"+SourceOrganForCombination[dd];

                GraphNameXYValues[First1][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][labelX][SourceOrganForCombination[dd]][TargetOrganForCombination[dd]];

                bool isIn = false;
                for ( int df = 0 ; df < xticks.size(); df++  ){if(xticks[df] == iii ){ isIn = true; } }
                if(isIn == false){ xlabels.push_back(labelX); xticks.push_back(iii);}

                isIn = false;
                for ( int df = 0 ; df < GraphsNames.size(); df++  ){if(GraphsNames[df] == First1 ){ isIn = true; } }
                if(isIn == false){GraphsNames.push_back(First1);}

                //    QTextStream(stdout) << " First1 " << First1 <<  " labelX " << labelX << "\n";
                iii = iii + 1.;
            }
        }
    }

    if(GraphNameXYValues.size() == 0){
        QMessageBox::information(this, tr(""), "No data for " + QuantitiesToScore +", geometry, radiotracer (or particle), source and/or target are registered for plot, check the result file");
        return;
    }

    /*
    for ( auto Abeg = GraphNameXYValues.begin(); Abeg != GraphNameXYValues.end(); ++Abeg  ){
        QString First1 = Abeg.key();
        for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
            double Second2 = Bbeg.key();

        iii = 1.;

        for(int dd=0; dd < SourceOrganlist.size();dd++){

            QString labelX= SourceOrganlist[dd];
            GraphNameXYValues[First1][iii] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][First1][labelX][labelX];
            }
        }
    }
    */
    ui->customPlot->xAxis->setLabel(Xaxislabel);
    for ( int df = 0 ; df < GraphsNames.size(); df++  ){
        create_graphs( GraphNameXYValues[GraphsNames[df]] , GraphsNames[df], graphs_Title , 0);
    }
}

void PlotDialog::on_plotButtonGenerateGraphs_clicked()
{
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << "() ------------------" << "\n";
    
    GraphTypeShown = "SC";
    
    GetPlotInputDataAndReadFiles();
    
    graphs_Title = "Graph Title";
    
    if( GraphsData == "Result"){
        
        if( Compare_Type == "Self"){
            graphs_Title = "Self Irradiation - " +  ParticleName;
        }
        
        if( Compare_Type == "Cross"){
            graphs_Title = SourceOrgan + " - Cross Irradiation - " +  ParticleName;
        }

        UseErrorBar = true;
        create_Data_for_Cross_and_Self_graphs( ResultTable , "DoseCalcs" );
        
    }else if(GraphsData == "Reference_Result" ){
        
        if( Compare_Type == "Self"){
            graphs_Title = SourceOrgan + " - Self Irradiation - " +  ParticleName;
        }
        
        if( Compare_Type == "Cross"){

            graphs_Title = SourceOrgan + " -> "+ TargetOrgan + " - " +  ParticleName;
        }
        
        UseErrorBar = true;
        create_Data_for_Cross_and_Self_graphs( ResultTable , "DoseCalcs" );
        create_Data_for_Cross_and_Self_graphs( ReferenceTable , CompareReferenceName );
        if(CompareReferenceName_2 != "" && QFile::exists(ui->PlotLineEditReferenceFile_2->text())){
            create_Data_for_Cross_and_Self_graphs( ReferenceTable_2 , CompareReferenceName_2 );
        }
    }
}
void PlotDialog::on_pushButtonGenerateGraphsForEnergy_clicked()
{
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << " ------------------" << "\n";

    GraphTypeShown = "OneE";

    GetPlotInputDataAndReadFiles();

    create_Data_for_Cross_and_Self_ForOneEnergy_graphs();
}
void PlotDialog::on_plotButtonGenerateMCSimErrorGraphs_clicked()
{
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << "() ------------------" << "\n";
    
    GraphTypeShown = "RSD";
    
    GetPlotInputDataAndReadFiles();
    create_Data_for_MCError_graphs();
}
void PlotDialog::on_plotButtonGenerateRegionVarGraphs_clicked()
{
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << "() ------------------" << "\n";
    
    GraphTypeShown = "V";
    
    GetPlotInputDataAndReadFiles();
    create_Data_for_RegionParametre_graphs();
}
void PlotDialog::on_plotButtonGenerateCrossSectionGraphGraphs_2_clicked()
{
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << "() ------------------" << "\n";
    
    GraphTypeShown = "MCS";
    
    GetPlotInputDataAndReadFiles();
    create_Data_for_CrossSection_graphs();
}
void PlotDialog::on_plotButtonGenerateSimulationTimeGraph_clicked()
{
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << "() ------------------" << "\n";
    
    GraphTypeShown = "T";
    
    GetPlotInputDataAndReadFiles();
    create_Data_for_SourceTimeSimulation_graphs(ResultParticleSourceEnergyTime, "");
}
void PlotDialog::on_pushButtonRadioTracerBarsAndRatiosPlots_clicked()
{
    QDialog * d = new QDialog();
    QGridLayout* GraphLayout = new QGridLayout;
    QLabel* lb = new QLabel("Plot Types");
    QRadioButton* BarsRadio = new QRadioButton("Bars");
    QRadioButton* GraphsRadio = new QRadioButton("Graphs");

    GraphsRadio->setChecked(true);

    QStringList PlotData;
    PlotData  << ""
              << "f(x)=RadioTracers, x=Sources, for Geometry"
              << "f(x)=Sources, x=RadioTracers, for Geometry"
              << "f(x)=Geometries, x=Sources, for RadioTracer"
              << "f(x)=Sources, x=Geometries, for RadioTracer"
              << "f(x)=Geometries, x=RadioTracers, for Sources"
              << "f(x)=RadioTracers, x=Geometries, for Sources"
              << "f(x)=SourceTargets, x=RadioTracers, for Geometry"
              << "f(x)=RadioTracers, x=SourceTargets, for Geometry"
              << "f(x)=SourceTargets, x=Geometries, for RadioTracer"
              << "f(x)=Geometries, x=SourceTargets, for RadioTracer"
              << "f(x)=Geometries, x=RadioTracers, for SourceTarget"
              << "f(x)=RadioTracers, x=Geometries, for SourceTarget"

              << "f(x)=Quantity, x=SourceTargets, for RadioTracer and Geometry"
                 ;

    QComboBox* BarsData = new QComboBox(); BarsData->addItems(PlotData);
    BarsData->setToolTip("Choose the plot data");
    d->setWindowTitle("RadioTracer Plots");
    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));
    int ii = 0, jj=0;
    GraphLayout->addWidget(lb, jj,ii,1,1);
    GraphLayout->addWidget(BarsRadio, jj,++ii,1,1);
    GraphLayout->addWidget(GraphsRadio, jj,++ii,1,1);
    GraphLayout->addWidget(BarsData, jj,++ii,1,1);
    GraphLayout->addWidget(buttonBox);
    d->setLayout(GraphLayout);

    int result = d->exec();

    if(result == QDialog::Accepted){
        if(BarsData->currentText() != "" || BarsData->currentText().isEmpty()){
            RadioTraceBarPlotData = BarsData->currentText();
            //QString firstword = RadioTraceBarPlotData;
            //firstword.resize(5); //"Graph"



            //if(firstword == "Graph"){

            if(GraphsRadio->isChecked()){

                GraphTypeShown = "RadioTracerQuantityGraph";
                GetPlotInputDataAndReadFiles();
                if(RadioTraceBarPlotData == "f(x)=Quantity, x=SourceTargets, for RadioTracer and Geometry"){
                    create_QuantitiesData_Cross_Self_FromRadioTracerIntake_graphs();
                }else{
                    create_QuantitiesData_For_Radiotracer_In_Plot_Graphs();
                }
            }
            else if (BarsRadio->isChecked()){
                GraphTypeShown = "PB";
                GetPlotInputDataAndReadFiles();
                create_QuantitiesData_For_Radiotracer_In_Plot_Bars();
            }
        }
    }
}
void PlotDialog::on_plotButtonGenerateRelErrGraphs_clicked()
{
    QTextStream(stdout) << "---------------- "<< __FUNCTION__ << " ------------------" << "\n";

    if(ui->PlotcomboBoxGraphData->currentText() == "Result"){
        ui->PlotcomboBoxGraphData->setCurrentText("Reference_Result");
    }

    //else if ( ui->checkBoxDiffForRadiotracerOr->isChecked() || (!ui->checkBoxDiffForRadiotracerOr->isChecked() && ui->radioButtonRa->isChecked()) ){ // for radiotracer
    //else if ( ui->checkBoxDiffForRadiotracerOr->isChecked()){

    QStringList PlotData;

    if(ui->checkBoxDiffForRadiotracerOr->isChecked()){
        if(ReferenceQuantityGeometryRadioTracerSourceTargetValues.size() == 0){
            QMessageBox::information(this, tr(""), "Canno't save radiotracers reference 1 data from file, check the file path or file data");
            return;
        }

        PlotData  << ""
                  << "f(x)=Sources, x=RadioTracers, for Geometry"
                  << "f(x)=SourceTargets, x=RadioTracers, for Geometry"
                  << "f(x)=Sources, x=Geometries, for RadioTracer"
                  << "f(x)=SourceTargets, x=Geometries, for RadioTracer"
                  << "f(x)=Geometries, x=SourceTargets, for RadioTracer"
                  << "f(x)=RadioTracers, x=Geometries, for Source"
                  << "f(x)=RadioTracers, x=Geometries, for SourceTarget";
    }
    //else if ( !ui->checkBoxDiffForRadiotracerOr->isChecked() && ui->radioButtonRa->isChecked() ){ // for radiotracer
    else if ( !ui->checkBoxDiffForRadiotracerOr->isChecked() ){ // for particles

        if(ReferenceTable.size() == 0){
            QMessageBox::information(this, tr(""), "Canno't save particles reference 1 data from file, check the file path or file data");
            return;
        }

        PlotData  << ""
                  << "f(x)=Sources, x=SourceTargets, for for Particle energies and Geometry"
                  << "f(x)=Geometries, x=Sources, for Particle energies"
                  << "f(x)=Geometries, x=SourceTargets, for Particle energies"
                  << "f(x)=Geometries, x=SourceTargets, for Particle and Energy"
                  << "f(x)=Energies, x=SourceTargets , for Geometry and Particle"
                  << "f(x)=SourceTargets, x=Energies , for Geometry and Particle"
                  << "f(x)=Energies, x=Sources , for Geometry and Particle"
                  << "f(x)=Sources, x=Energies , for Geometry and Particle"

                     //<< "Ratio_EnergyForSourceAndSourceTarget"
                     ;

    }

    QDialog * d = new QDialog();
    QGridLayout* GraphLayout = new QGridLayout;
    QLabel* lb = new QLabel("Plot Types");

    QRadioButton* RatioRadio = new QRadioButton("RatioPlot");
    GraphsRadio = new QRadioButton("Graphs");

    GraphsRadio->setChecked(true);

    QComboBox* BarsData = new QComboBox(); BarsData->addItems(PlotData);
    BarsData->setToolTip("Choose the plot data");
    d->setWindowTitle("Ratio Plots");
    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));
    int ii = 0, jj=0;
    GraphLayout->addWidget(lb, jj,ii,1,1);
    GraphLayout->addWidget(RatioRadio, jj,++ii,1,1);
    GraphLayout->addWidget(GraphsRadio, jj,++ii,1,1);
    GraphLayout->addWidget(BarsData, jj,++ii,1,1);
    GraphLayout->addWidget(buttonBox);
    d->setLayout(GraphLayout);
    int result = d->exec();

    if(result == QDialog::Accepted){
        if(BarsData->currentText() != "" || BarsData->currentText().isEmpty()){

            RadioTraceBarPlotData = BarsData->currentText();
            /*
            if(RadioTraceBarPlotData == "Ratio_EnergyForSourceAndSourceTarget"){

                GraphTypeShown = "RE";
                GetPlotInputDataAndReadFiles();
                create_Data_for_RelativeErr_graphs();
            }
            else{
*/
            if(GraphsRadio->isChecked()){
                GraphTypeShown = "RadioTracerQuantityGraph";
            }else{
                GraphTypeShown = "Rat";
            }
            GetPlotInputDataAndReadFiles();
            create_RatioData_For_ParticleAndRadiotracer_In_RatioPlot();
        }
    }
    //}
}

double PlotDialog::RelativeDifferenceCalculation(double val1, double val2){

    double a1 = val1;
    double a2 = val2;
    double a3;
    if(ui->radioButtonLRD->isChecked()){
        if( a1 != 0. && a2 != 0. ){
            a3 = std::log(a1/a2)*100;
        }else{
            a3= NULL;
        }
    }
    else if(ui->radioButtonRD->isChecked()){
        if( a1 != 0. && a2 != 0. ){
            a3 = (a1-a2/a2)*100;
        }else{
            a3= NULL;
        }
    }
    else if(ui->radioButtonRa->isChecked()){
        if( a1 != 0. && a2 != 0. ){
            a3 = (a1/a2);
        }
        else{
            a3 = NULL;
        }
    }else{
        if( a1 != 0. && a2 != 0. ){
            a3 = (a1-a2/a2)*100;
        }else{
            a3 = NULL;
        }
    }

    return a3;
}

// SLOT related to the dialog components interaction
// executed when the read result, Reference and simulation inputs files buttons clicked
void PlotDialog::on_PlotButtonReadAllData_clicked()
{
    on_PlotButtonReadSimulationInpFile_clicked();
    on_PlotButtonReadResultFile_clicked();
    on_PlotButtonReadReferenceFile_clicked();
    on_PlotButtonReadReferenceFile_2_clicked();
    //on_PlotButtonReadCrossSectionFile_clicked();
}
void PlotDialog::on_PlotButtonReadResultFile_clicked()
{

    QTextStream(stdout) <<"-----------------------------on_PlotButtonReadResultFile_clicked----------------------------------" << "\n";

    SourceOrganlist.clear();SourceOrganlist.empty();
    TargetOrganlist.clear();TargetOrganlist.empty();
    Particlelist.clear();Particlelist.empty();
    ScoreVarlist.clear();ScoreVarlist.empty();
    Geometrylist.clear();Geometrylist.empty();
    //RadioTracerlist.clear();RadioTracerlist.empty();

    ResultFilePath = ui->PlotLineEditResultFile->text();
    ResultTable = fileManagerObjectPlot->Read_final_result_file(ResultFilePath, 1);// 1 for result
    ResultQuantityGeometryRadioTracerSourceTargetValues = fileManagerObjectPlot->getResultQuantityGeometryRadioTracerSourceTargetValues();
    ErrorTable = fileManagerObjectPlot->Read_final_result_file(ResultFilePath, 2);// else for Error
    StandartDeviation = fileManagerObjectPlot->getStandartDeviationMap();
    RegionParameterValueMap = fileManagerObjectPlot->getRegionParameterValueMap();
    ResultParticleSourceEnergyTime = fileManagerObjectPlot->getResultParticleSourceEnergyTime();

    //OrganNamesToScore.clear();
    ResEnergies.clear();
    ResEnergies = fileManagerObjectPlot->getResEnergies();
    if(fileManagerObjectPlot->getMinValForLog() < MinValForLog){
        MinValForLog = fileManagerObjectPlot->getMinValForLog();
    }
    
    setComboboxInputs();

}
void PlotDialog::on_PlotButtonReadReferenceFile_clicked()
{
    ReferenceFilePath = ui->PlotLineEditReferenceFile->text();
    ReferenceTable = fileManagerObjectPlot->Read_Comparison_file(ReferenceFilePath, 12);
    ReferenceQuantityGeometryRadioTracerSourceTargetValues = fileManagerObjectPlot->getReferenceQuantityGeometryRadioTracerSourceTargetValues();

    if(fileManagerObjectPlot->getMinValForLog() < MinValForLog){
        MinValForLog = fileManagerObjectPlot->getMinValForLog();
    }
    ui->plotMessageLabel->setText("Reading the reference file");

    if( CompareReferenceName == "" && CompareReferenceName.isEmpty() && ReferenceFilePath.isEmpty() && ReferenceFilePath==""){
        return;
    }

    if(ReferenceTable.size() == 0){

        QMessageBox::information(this, tr(""), "Canno't save reference 1 data from file, check the file path or file data");
        ui->plotMessageLabel->setText("Canno't save data from file, check the file path or file data of Reference");
    }else{
        ui->plotMessageLabel->setText("Now you can compare the result with reference, if reference graph not shown, please check if required data exists in reference file(The chosen organs should be listed in source and target list in file)");
    }

}
void PlotDialog::on_PlotButtonReadReferenceFile_2_clicked()
{
    
    ReferenceFilePath_2 = ui->PlotLineEditReferenceFile_2->text();
    ReferenceTable_2 = fileManagerObjectPlot->Read_Comparison_file(ReferenceFilePath_2, 12);
    if(fileManagerObjectPlot->getMinValForLog() < MinValForLog){
        MinValForLog = fileManagerObjectPlot->getMinValForLog();
    }
    if( CompareReferenceName_2 == "" && CompareReferenceName_2.isEmpty() && ReferenceFilePath_2.isEmpty() && ReferenceFilePath_2==""){
        return;
    }
    if(ReferenceFilePath_2.size() == 0){

        QMessageBox::information(this, tr(""), "Canno't save reference 2 data from file, check the file path or file data");

        ui->plotMessageLabel->setText("Canno't save data from file, check the file path or file data of Reference_2");
    }else{
        ui->plotMessageLabel->setText("Now you can compare the result and reference_1 with reference_2, if reference_2 graph not shown, please check if required data exists in reference_2 file");
    }
}
void PlotDialog::on_PlotButtonReadSimulationInpFile_clicked()
{
    AnalysisInputFilePath = ui->PlotLineEditSimulationInpFile->text();
    AnalysisInputMap = fileManagerObjectPlot->ReadLinesFromFileWithFirstWordIndicator(AnalysisInputFilePath);

    QStringList InputsVals ;

    if(QFile::exists(AnalysisInputMap[RunAndScoreCommands[10]])){
        //ui->PlotLineEditResultsDir->setText(AnalysisInputMap[RunAndScoreCommands[10]]);
        UserCurrentResultsDirPath = AnalysisInputMap[RunAndScoreCommands[10]];
        ui->PlotLineEditResultFile->setText(UserCurrentResultsDirPath+"/"+ResultFileName);
        if(!QFile::exists(UserCurrentResultsDirPath+"/"+GraphsOutDirName)){
            QDir* dir = new QDir(UserCurrentResultsDirPath);
            dir->mkdir(GraphsOutDirName);
        }
    }

    InputsVals = AnalysisInputMap[AnalysisCommands[0]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    if(InputsVals.size() > 0){
        ui->PlotcomboBoxGraphData->setCurrentText(InputsVals[0]);
    }

    InputsVals = AnalysisInputMap[AnalysisCommands[0]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"

    if(InputsVals.size() > 3){
        CompareReferenceName = InputsVals[2];
        ReferenceFilePath = InputsVals[3];
        ui->plotlineEditReferenceName->setText(CompareReferenceName);
        ui->PlotLineEditReferenceFile->setText(ReferenceFilePath);
        if(InputsVals.size() > 5){
            CompareReferenceName_2 = InputsVals[4];
            ReferenceFilePath_2 = InputsVals[5];
            ui->PlotLineEditReferenceFile_2->setText(ReferenceFilePath_2);
            ui->plotlineEditReferenceName->setText(ui->plotlineEditReferenceName->text()+" "+CompareReferenceName_2);
        }
    }

    InputsVals = AnalysisInputMap[AnalysisCommands[3]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    if(InputsVals.size() > 0){
        ui->PlotComboboxRegionVariable->setCurrentText(InputsVals[0]);
    }

    OrganNamesToScore.clear();
    InputsVals = AnalysisInputMap[RunAndScoreCommands[0]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/GeometryData/createWorld"
    for (int var = 0; var < InputsVals.size(); ++var) {OrganNamesToScore.push_back(InputsVals[var]);}

    if(AnalysisInputMap.size() == 0){

        QMessageBox::information(this, tr(""), "Canno't save simulation data from macros file, check the file path or macros data");

        ui->plotMessageLabel->setText("Canno't save Analysis input from file, check the file path or file data of Reference");
    }

    if(AnalysisInputMap.size() != 0){
        ui->plotMessageLabel->setText("Now you can use region data and analysis inputs, please check the data in macros file to ensure your data to be used");
    }

}
void PlotDialog::on_PlotButtonReadCrossSectionFile_clicked()
{
    CrossSectionFilePath = ui->PlotLineEditCrossSectionFile->text();
    ParticleMaterialProcessEnergySigmaMap = fileManagerObjectPlot->Read_CrossSection_file(CrossSectionFilePath);
    //ui->plotMessageLabel->setText("Reading the Cross Section file");
    
    QStringList ParticleForSigmalist;
    QStringList MaterialsForSigmalist;
    QStringList ProcessForSigmalist;
    
    int ff = 0;
    
    ProcessForSigmalist.push_back("All"); // this is used for printing all process cross section
    
    for ( auto abeg = ParticleMaterialProcessEnergySigmaMap.begin(); abeg != ParticleMaterialProcessEnergySigmaMap.end(); ++abeg  ){
        QString Particle = abeg.key();
        
        ff=0;
        for(int i=0; i<ParticleForSigmalist.size(); i++)
        {
            if(ParticleForSigmalist.indexOf(QRegExp(Particle, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                ff++;
            }
        }
        if(ff == 0){
            ParticleForSigmalist.push_back(Particle);
        }
        
        for ( auto Abeg = abeg.value().begin(); Abeg != abeg.value().end(); ++Abeg  ){
            
            QString MatName = Abeg.key();
            
            ff=0;
            for(int i=0; i<MaterialsForSigmalist.size(); i++)
            {
                if(MaterialsForSigmalist.indexOf(QRegExp(MatName, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                    ff++;
                }
            }
            if(ff == 0){
                MaterialsForSigmalist.push_back(MatName);
            }
            
            // iterations on source name
            for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){
                
                QString ProcName = Bbeg.key();
                
                ff=0;
                for(int i=0; i<ProcessForSigmalist.size(); i++)
                {
                    if(ProcessForSigmalist.indexOf(QRegExp(ProcName, Qt::CaseInsensitive, QRegExp::Wildcard), i)-i == 0){
                        ff++;
                    }
                }
                if(ff == 0){
                    ProcessForSigmalist.push_back(ProcName);
                }
            }
        }
    }
    
    for(int i=0; i<ParticleForSigmalist.size(); i++)
    {
        QTextStream(stdout) <<"ParticleForSigmalist[] : " << ParticleForSigmalist[i] << "\n";
    }
    for(int i=0; i<MaterialsForSigmalist.size(); i++)
    {
        QTextStream(stdout) <<"MaterialsForSigmalist[] : " << MaterialsForSigmalist[i] << "\n";
    }
    for(int i=0; i<ProcessForSigmalist.size(); i++)
    {
        QTextStream(stdout) <<"ProcessForSigmalist[] : " << ProcessForSigmalist[i] << "\n";
    }
    
    ui->comboBoxParticleForCrossSection->clear();
    ui->comboBoxMaterialForCrossSection->clear();
    ui->comboBoxProcessForCrossSection->clear();
    
    ui->comboBoxParticleForCrossSection->addItems(ParticleForSigmalist);
    ui->comboBoxMaterialForCrossSection->addItems(MaterialsForSigmalist);
    ui->comboBoxProcessForCrossSection->addItems(ProcessForSigmalist);
    
    if(ParticleMaterialProcessEnergySigmaMap.size() == 0){

        QMessageBox::information(this, tr(""), "Canno't save cross section data from file, check the file path or file data");

        ui->plotMessageLabel->setText("Canno't save cross section data from file, check the file path or file data of Reference");
    }else{
        ui->plotMessageLabel->setText("Now you can generate cross section graphs, please check if required data exists in cross section data file");
    }
}

void PlotDialog::on_PlotButtonEditSimData_clicked()
{
    QString command = DoseCalcsCore_build_dir_path+"/"+MacroFileName;
    QProcess process;
    QStringList qsl = {command};
    process.startDetached("gedit", qsl);
}
void PlotDialog::on_PlotButtonEditResultsData_clicked()
{
    QString command = UserCurrentResultsDirPath+"/"+ResultFileName;
    QProcess process;
    QStringList qsl = {command};
    process.startDetached("gedit", qsl);
}

void PlotDialog::on_PlotButtonChoseSimulationInpFileDialog_clicked()
{    
    AnalysisInputFilePath = QFileDialog::getOpenFileName( this, tr("Choose the simulation macros file (input file)"), DoseCalcsCore_build_dir_path, tr("Simulation data files (*)") );
    MacroFilePath = AnalysisInputFilePath;
    ui->PlotLineEditSimulationInpFile->setText(AnalysisInputFilePath);
    on_PlotButtonReadSimulationInpFile_clicked();
}
void PlotDialog::on_PlotButtonChoseReferenceFileDialog_clicked()
{
    ReferenceFilePath = QFileDialog::getOpenFileName( this, tr("Choose the first reference file"), UserCurrentResultsDirPath, tr("Reference files (*)"));
    ui->PlotLineEditReferenceFile->setText(ReferenceFilePath);
    on_PlotButtonReadReferenceFile_clicked();
}
void PlotDialog::on_PlotButtonChoseReferenceFileDialog_2_clicked()
{
    ReferenceFilePath_2 = QFileDialog::getOpenFileName( this, tr("Choose the second reference file"), UserCurrentResultsDirPath, tr("Reference files (*)"));
    ui->PlotLineEditReferenceFile_2->setText(ReferenceFilePath_2);
    on_PlotButtonReadReferenceFile_clicked();
}
void PlotDialog::on_PlotButtonChoseResultFileDialog_clicked()
{
    ResultFilePath = QFileDialog::getOpenFileName( this, tr("Choose the result file"), UserCurrentResultsDirPath, tr("Results data files (ResultsData*)") );
    ui->PlotLineEditResultFile->setText(ResultFilePath);
    on_PlotButtonReadResultFile_clicked();
}
void PlotDialog::on_PlotButtonChoseCrossSectionFileDialog_2_clicked()
{
    CrossSectionFilePath = QFileDialog::getOpenFileName( this, tr("Choose the Cross Section file"), UserCurrentResultsDirPath, tr("Cross section data files(CrossSection*)") );
    ui->PlotLineEditCrossSectionFile->setText(CrossSectionFilePath);
    on_PlotButtonReadCrossSectionFile_clicked();
}

// graphs and generating table, saving in a format
void PlotDialog::on_plotButtonSaveGraph_clicked()
{
    if(!QFile::exists(UserCurrentResultsDirPath+"/"+GraphsOutDirName)){
        QDir* dir = new QDir(UserCurrentResultsDirPath);
        dir->mkdir(GraphsOutDirName);
    }

    bool ok;
    QString newLabel = QInputDialog::getText(this, "Enter the File Name", "New File Name:", QLineEdit::Normal, graphs_Title, &ok);
    if (ok)
    {
        FileExtToSave = ui->plotcomboxFileExtToSave->currentText();
        if(FileExtToSave == "Bmp"){
            ui->customPlot->saveBmp(UserCurrentResultsDirPath+"/"+GraphsOutDirName +"/"+ newLabel);
        }else if(FileExtToSave == "Jpg"){
            ui->customPlot->saveJpg(UserCurrentResultsDirPath+"/" +GraphsOutDirName +"/"+ newLabel);
        }else if(FileExtToSave == "Pdf"){
            ui->customPlot->savePdf(UserCurrentResultsDirPath+"/" +GraphsOutDirName +"/"+ newLabel);
        }else if(FileExtToSave == "Png"){
            ui->customPlot->savePng(UserCurrentResultsDirPath+"/" +GraphsOutDirName +"/"+ newLabel);
        }
    }
}
void PlotDialog::on_plotButtonGenerateTableOfAllGraphsData_clicked()
{
    TableWidgetFillingWithAllGraphs();
}
void PlotDialog::on_plotButtonSaveTableOfAllGraphsData_clicked()
{
    if(!QFile::exists(UserCurrentResultsDirPath+"/"+GraphsOutDirName)){
        QDir* dir = new QDir(UserCurrentResultsDirPath);
        dir->mkdir(GraphsOutDirName);
    }

    tblFrom = "all";
    SaveTableFromTableWidgetToFile();
}
void PlotDialog::on_plotButtonSaveTableData_clicked()
{
    tblFrom = "one";
    SaveTableFromTableWidgetToFile();
}

void PlotDialog::on_PlotcomboBoxGraphData_currentIndexChanged(const QString &arg1)
{
    ActionWhenChoosingGraphData();
}
void PlotDialog::on_plotcomboBoxCompareType_currentIndexChanged(const QString &arg1)
{
    ActionWhenChoosingCompareType();
}
void PlotDialog::ActionWhenChoosingGraphData() { // called from SLOT and DataInitialization()
    
    if(ui->PlotcomboBoxGraphData->currentText()=="Result"){
        
        if(ui->plotcomboBoxCompareType->currentText()=="Self"){
            ui->comboBoxSourceOrgan->setEnabled(false);
        }else if(ui->plotcomboBoxCompareType->currentText()=="Cross"){
            ui->comboBoxSourceOrgan->setEnabled(true);
        }
        
        //ui->PlotcomboBoxTargetOrgan->setEnabled(false);
        //ui->PlotLineEditReferenceFile->setEnabled(false);
        //ui->PlotButtonChoseReferenceFileDialog->setEnabled(false);
        //ui->plotlineEditReferenceName->setEnabled(false);
        
    }else if(ui->PlotcomboBoxGraphData->currentText()=="Reference_Result" ){
        
        if(ui->plotcomboBoxCompareType->currentText()=="Self"){
            //ui->PlotcomboBoxTargetOrgan->setEnabled(false);
        }else if(ui->plotcomboBoxCompareType->currentText()=="Cross"){
            ui->PlotcomboBoxTargetOrgan->setEnabled(true);
        }
        
        ui->comboBoxSourceOrgan->setEnabled(true);

        on_PlotButtonReadReferenceFile_clicked();
        on_PlotButtonReadReferenceFile_2_clicked();

        //ui->PlotLineEditReferenceFile->setEnabled(true);
        //ui->PlotButtonChoseReferenceFileDialog->setEnabled(true);
        //ui->plotlineEditReferenceName->setEnabled(true);
    }
    
    //ui->PlotLineEditReferenceFile->update();
    
}
void PlotDialog::ActionWhenChoosingCompareType() { // called from SLOT and DataInitialization()
    
    if(ui->plotcomboBoxCompareType->currentText()=="Self"){
        
        if(ui->PlotcomboBoxGraphData->currentText()=="Reference_Result" ){
            //ui->PlotcomboBoxTargetOrgan->setEnabled(false); // i dont need source because it show for all sources
            ui->comboBoxSourceOrgan->setEnabled(true); // i need source
            //ui->PlotLineEditReferenceFile->setEnabled(true);
            //ui->PlotButtonChoseReferenceFileDialog->setEnabled(true);
            //ui->plotlineEditReferenceName->setEnabled(true);
            
        }else if (ui->PlotcomboBoxGraphData->currentText()=="Result" ){
            //ui->PlotcomboBoxTargetOrgan->setEnabled(false); // i dont need target
            ui->comboBoxSourceOrgan->setEnabled(false); // i dont need source
            //ui->PlotLineEditReferenceFile->setEnabled(false);
            //ui->PlotButtonChoseReferenceFileDialog->setEnabled(false);
            //ui->plotlineEditReferenceName->setEnabled(false);
        }
        
    }else if(ui->plotcomboBoxCompareType->currentText()=="Cross"){
        if(ui->PlotcomboBoxGraphData->currentText()=="Reference_Result" ){
            ui->PlotcomboBoxTargetOrgan->setEnabled(true); // i need target
            ui->comboBoxSourceOrgan->setEnabled(true); // i need source
            //ui->PlotLineEditReferenceFile->setEnabled(true);
            //ui->PlotButtonChoseReferenceFileDialog->setEnabled(true);
            //ui->plotlineEditReferenceName->setEnabled(true);
            
        }else if (ui->PlotcomboBoxGraphData->currentText()=="Result" ){
            //ui->PlotcomboBoxTargetOrgan->setEnabled(false); // i dont need target because it show for all for all targets
            ui->comboBoxSourceOrgan->setEnabled(true); // i need source
            //ui->PlotLineEditReferenceFile->setEnabled(false);
            //ui->PlotButtonChoseReferenceFileDialog->setEnabled(false);
            //ui->plotlineEditReferenceName->setEnabled(false);
        }
    }
}

// Table Handling : Print to table and save to file
void PlotDialog::TableWidgetFillingWithAllGraphs(){
    
    ui->tableWidget->clear();
    ui->tableWidget->setRowCount(0);
    ui->tableWidget->setColumnCount(0);
    
    QStringList headers;
    QVector<double> Energies ;
    QMap<int,QMap<double,double>> SelectedGraphIndexEnergyScoreVal;
    
    
    QTextStream(stdout) << "RegionVariableNameWithUnit " << RegionVariableNameWithUnit << "\n";
    
    if(GraphTypeShown == "V" ){
        headers.append(tr((RegionVariableNameWithUnit).toStdString().c_str()));
    }else if(GraphTypeShown == "OneE" || GraphTypeShown == "RadioTracerQuantityGraph" ){
        headers.append(tr(("Regions " +RegionVariableNameWithUnit).toStdString().c_str()));
    }else { // for all graphs E(MeV) is the x axix label except for V type
        headers.append(tr("E(MeV)"));
    }
    
    //headers.append(tr("E(MeV)"));
    for(int fg = 0; fg < ui->customPlot->graphCount() ; fg++){
        
        //headers.append(tr((ScoreVariableUnit+ui->customPlot->graph(fg)->name()).toStdString().c_str()));
        
        if(GraphTypeShown == "SC" ){
            headers.append(tr((ui->customPlot->graph(fg)->name()).toStdString().c_str()));
        }else if(GraphTypeShown == "RSD" ){
            headers.append(tr((ui->customPlot->graph(fg)->name()).toStdString().c_str()));
        }else if(GraphTypeShown == "RE" ){
            headers.append(tr((ui->customPlot->graph(fg)->name()).toStdString().c_str()));
        }else if(GraphTypeShown == "V" ){
            headers.append(tr((ScoreVariableUnit+" for "+ui->customPlot->graph(fg)->name()).toStdString().c_str()));
        }else if(GraphTypeShown == "MCS" ){
            headers.append(tr((ui->customPlot->graph(fg)->name()+"(cm-1)").toStdString().c_str()));
        }else if(GraphTypeShown == "T" ){
            headers.append(tr((ui->customPlot->graph(fg)->name()+" Time (min)").toStdString().c_str()));
        }else if(GraphTypeShown == "RadioTracerQuantityGraph" ){
            headers.append(tr((ui->customPlot->graph(fg)->name()).toStdString().c_str()));
        }else if(GraphTypeShown == "OneE" ){
            headers.append(tr((ui->customPlot->graph(fg)->name()).toStdString().c_str()));
        }else if(GraphTypeShown == "XY" ){
            headers.append(tr("Y(Unit)"));
        }else{
            headers.append(tr("f(E)"));
        }
        
        for ( auto it = ui->customPlot->graph(fg)->data()->begin(); it != ui->customPlot->graph(fg)->data()->end(); ++it  ){
            SelectedGraphIndexEnergyScoreVal[fg][it->key] = it->value;
            if(fg == 0){
                Energies.push_back(it->key);
            }
        }
    }

    if(GraphsData == "Reference_Result" && GraphTypeShown == "SC"){
        QString axilbl ;
        if(ui->radioButtonLRD->isChecked()){
            axilbl = "Logarithmic Relative Difference (%)";
        }else if(ui->radioButtonRa->isChecked()){
            axilbl = "Ratio";
        }else{
            axilbl = "Relative Difference (%)";
        }
        headers.append(axilbl);
        ui->tableWidget->setColumnCount(ui->customPlot->graphCount()+2); // we add (2) last column for relative Difference between reference and result
    }else{
        ui->tableWidget->setColumnCount(ui->customPlot->graphCount()+1);
    }
    ui->tableWidget->setShowGrid(true);
    ui->tableWidget->setSelectionMode(QAbstractItemView::SingleSelection);
    ui->tableWidget->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->tableWidget->setHorizontalHeaderLabels(headers);
    ui->tableWidget->horizontalHeader()->setStretchLastSection(true);
    
    //ui->tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    
    int EnergyIncre = 0 ;
    int graphIncre = 0 ;
    int RefRescolNum = ui->tableWidget->columnCount();
    int LastColindex = ui->tableWidget->columnCount()-1;
    
    for(int row = 0 ; row < Energies.size(); row++){
        // Insert row
        
        graphIncre = 0;
        
        ui->tableWidget->insertRow(row);
        
        for(int col = 0 ; col < RefRescolNum ; col++){
            //QTextStream(stdout) << "col" << col << "\n";
            
            if(col == 0){
                ui->tableWidget->setItem(row,col, new QTableWidgetItem(QString::number(Energies[EnergyIncre])));
                EnergyIncre++;
            }
            else if(LastColindex == col && GraphsData == "Reference_Result" && GraphTypeShown == "SC" ){                // fill the error column
                //QTextStream(stdout) << "RefRescolNum == col " << RefRescolNum  << " " << col << "\n";
                
                int RefGraIndex = 1; int ResGraIndex = 0;
                if(ui->customPlot->graph(0)->name() == CompareReferenceName){
                    RefGraIndex = 0 ; ResGraIndex = 1;
                }
                else{
                    RefGraIndex = 1 ; ResGraIndex = 0;
                }
                
                //double ReltiveError = (std::abs(SelectedGraphIndexEnergyScoreVal[RefGraIndex][Energies[EnergyIncre-1]]-SelectedGraphIndexEnergyScoreVal[ResGraIndex][Energies[EnergyIncre-1]])/SelectedGraphIndexEnergyScoreVal[RefGraIndex][Energies[EnergyIncre-1]])*100;

                double a1 = SelectedGraphIndexEnergyScoreVal[ResGraIndex][Energies[EnergyIncre-1]];
                double a2 = SelectedGraphIndexEnergyScoreVal[RefGraIndex][Energies[EnergyIncre-1]];
                double a3 ;
                if(ui->radioButtonLRD->isChecked()){
                    if( a1 != 0. && a2 != 0. ){
                        a3 = std::log(a1/a2)*100;
                    }else{
                        a3= 100;
                    }
                }
                else if(ui->radioButtonRD->isChecked()){
                    if( a1 != 0. && a2 != 0. ){
                        a3 = (a1-a2/a2)*100;
                    }else{
                        a3= 100;
                    }
                }
                else if(ui->radioButtonRa->isChecked()){
                    if( a2 != 0. ){
                        a3 = a1/a2;
                    }else{
                        a3 = 100;
                    }
                    /*if( a1 != 0. && a2 != 0. ){
                        if( a2 < a1){
                            a3 = (a2/a1)*100;
                        }
                        else{
                            a3 = -(a1/a2)*100;
                        }
                    }else{
                        a3 = 100;
                    }
                    */
                }else{
                    if( a1 != 0. && a2 != 0. ){
                        a3 = (a1-a2/a2)*100;
                    }else{
                        a3= 100;
                    }
                }


                ui->tableWidget->setItem(row,col, new QTableWidgetItem(QString::number(a3)));
            }
            else{
                ui->tableWidget->setItem(row,col, new QTableWidgetItem(QString::number(SelectedGraphIndexEnergyScoreVal[col-1][Energies[EnergyIncre-1]])));
                graphIncre ++;
            }
            
        }
    }
    
    ui->tableWidget->resizeColumnsToContents();
    //ui->tableWidget->resizeRowsToContents();
    
    ui->tableWidget->horizontalHeader()->viewport()->installEventFilter(this);
    ui->tableWidget->verticalHeader()->viewport()->installEventFilter(this);
    
}
void PlotDialog::TableWidgetFillingWithSelectedGraph(){ // called from graphClicked()
    
    ui->tableWidgetForOneGraph->clear();
    ui->tableWidgetForOneGraph->setRowCount(0);
    ui->tableWidgetForOneGraph->setColumnCount(0);
    
    QStringList headers;
    QVector<double> Energies ;
    QVector<double> ScoreVals ;
    
    GraphName = ui->customPlot->selectedGraphs().first()->name();
    
    if(GraphTypeShown == "SC" ){
        headers.append(tr("E(MeV)"));
        headers.append(tr(ScoreVariableUnit.toStdString().c_str()));
    }else if(GraphTypeShown == "RSD"){
        headers.append(tr("E(MeV)"));
        headers.append(tr(("RSD\%"+ui->customPlot->selectedGraphs().first()->name()).toStdString().c_str()));
    }else if(GraphTypeShown == "RE"){
        headers.append(tr("E(MeV)"));
        headers.append(tr(("Err\%"+ui->customPlot->selectedGraphs().first()->name()).toStdString().c_str()));
    }else if(GraphTypeShown == "V" ){
        headers.append(tr((RegionVariableNameWithUnit).toStdString().c_str()));
        headers.append(tr((ScoreVariableUnit+" for "+ui->customPlot->selectedGraphs().first()->name()).toStdString().c_str()));
    }else if(GraphTypeShown == "MCS"){
        headers.append(tr("E(MeV)"));
        headers.append(tr((ui->customPlot->selectedGraphs().first()->name()+"(cm-1)").toStdString().c_str()));
    }else if(GraphTypeShown == "T" ){
        headers.append(tr("E(MeV)"));
        headers.append(tr((ui->customPlot->selectedGraphs().first()->name()+" Time (min)").toStdString().c_str()));
    }else if(GraphTypeShown == "OneE" ){
        headers.append(tr(("Regions "+RegionVariableNameWithUnit).toStdString().c_str()));
        headers.append(tr((ScoreVariableUnit+ui->customPlot->selectedGraphs().first()->name()).toStdString().c_str()));
    }else if(GraphTypeShown == "RadioTracerQuantityGraph" ){
        headers.append(tr(("Regions "+RegionVariableNameWithUnit).toStdString().c_str()));
        headers.append(tr((ScoreVariableUnit + " from "+ui->customPlot->selectedGraphs().first()->name()).toStdString().c_str()));
    }else if(GraphTypeShown == "RDL" ){
        headers.append(tr("Radionuclide"));
        headers.append(tr((ScoreVariableUnit).toStdString().c_str()));
    }else{
        headers.append(tr("E(MeV)"));
        headers.append(tr((ScoreVariableUnit+ui->customPlot->selectedGraphs().first()->name()).toStdString().c_str()));
    }
    

    for ( auto it = ui->customPlot->selectedGraphs().first()->data()->begin(); it != ui->customPlot->selectedGraphs().first()->data()->end(); ++it  ){
        Energies.push_back(it->key);
        ScoreVals.push_back(it->value);
    }
    
    ui->tableWidgetForOneGraph->setColumnCount(2);
    ui->tableWidgetForOneGraph->setShowGrid(true);
    ui->tableWidgetForOneGraph->setSelectionMode(QAbstractItemView::SingleSelection);
    ui->tableWidgetForOneGraph->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->tableWidgetForOneGraph->setHorizontalHeaderLabels(headers);
    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    
    ui->tableWidgetForOneGraph->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    
    for(int row = 0 ; row < Energies.size(); row++){
        // Insert row
        
        ui->tableWidgetForOneGraph->insertRow(row);
        
        ui->tableWidgetForOneGraph->setItem(row,0, new QTableWidgetItem(QString::number(Energies[row])));
        ui->tableWidgetForOneGraph->setItem(row,1, new QTableWidgetItem(QString::number(ScoreVals[row])));
    }
    
    ui->tableWidgetForOneGraph->resizeColumnsToContents();

    ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidgetForOneGraph->horizontalHeader()->viewport()->installEventFilter(this);

    //ui->tableWidgetForOneGraph->horizontalHeaderItem(0)->setFlags(Qt::ItemIsEditable);

    //QTableWidgetItem* itm = ui->tableWidgetForOneGraph->item( 1, 1 );
    //if (itm) QTextStream(stdout) << " \n\n\n ------------------ " << itm->text() << "\n";

    //ui->tableWidgetForOneGraph->horizontalHeader()->viewport()->installEventFilter(this);
    //ui->tableWidgetForOneGraph->verticalHeader()->viewport()->installEventFilter(this);
    
}
void PlotDialog::SaveTableFromTableWidgetToFile() {
    
    QTableWidget* tbl ;
    if(tblFrom == "one"){
        tbl = ui->tableWidgetForOneGraph;
    }else if(tblFrom == "all"){
        tbl = ui->tableWidget;
    }
    
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

    QString FileName = tr((UserCurrentResultsDirPath+"/" +GraphsOutDirName+"/" + tblFrom +"_"+ ScoreVariableUnit +"_"+ GraphName).toStdString().c_str());
    
    if(tblFrom == "one"){
        FileName = tr((tblFrom +"_"+ ScoreVariableUnit +"_"+ GraphName).toStdString().c_str());
    }
    else if(tblFrom == "all"){
        FileName = tr((tblFrom +"_"+ ScoreVariableUnit ).toStdString().c_str());
    }
    
    bool ok;
    QString newLabel = QInputDialog::getText(this, "Enter the File Name", "New File Name:", QLineEdit::Normal, UserCurrentResultsDirPath+"/" +GraphsOutDirName, &ok);
    if (ok)
    {
        FileName = newLabel;
    }
    
    printer.setOutputFileName( FileName );
    
    doc.setDocumentMargin(0);
    doc.setTextWidth(10);
    doc.print(&printer);
    //QTextStream(stdout) << "---------------- 5 ------------------" << "\n";
#endif

}

// the tableView is setted invisible in the constructor(), we should set it visible to test it
void PlotDialog::TableViewFillingWithAllGraphs(){
    //ui->customPlot->selectedGraphs()[ui->customPlot->selectedGraphs().size()];
    
    QTextStream(stdout) << "ui->customPlot->graphCount() : " << ui->customPlot->graphCount() << "\n";
    
    ui->customPlot->graph();
    //for(int fg = 0; fg < ui->customPlot->selectedGraphs().size()  ;fg++){
    
    for(int fg = 0; fg < ui->customPlot->graphCount() ; fg++){
        
        QCPDataContainer<QCPGraphData>::iterator Abeg = ui->customPlot->graph(fg)->data()->begin();
        QCPDataContainer<QCPGraphData>::iterator Aend = ui->customPlot->graph(fg)->data()->end();
        
        //QCPDataContainer<QCPGraphData>::iterator Abeg = ui->customPlot->selectedGraphs()[fg]->data()->begin();
        //QCPDataContainer<QCPGraphData>::iterator Aend = ui->customPlot->selectedGraphs()[fg]->data()->end();
        
        QVector<double> SelectedGraphEnergyForThisGraph;
        QVector<double> SelectedGraphScoreVarForThisGraph;
        
        while(Abeg!=Aend)
        {
            SelectedGraphEnergyForThisGraph.push_back(Abeg->key);
            SelectedGraphScoreVarForThisGraph.push_back(Abeg->value);
            
            Abeg++;
        }
        
        SelectedGraphEnergy = SelectedGraphEnergyForThisGraph;
        SelectedGraphScoreVar = SelectedGraphScoreVarForThisGraph;
        
        for(int RR = 0; RR< SelectedGraphEnergy.size() ;RR++){
            QTextStream(stdout) << SelectedGraphEnergy[RR] << " - " << SelectedGraphScoreVar[RR] << "\n";
        }
        
        if(fg==0){
            // Create a new model. QStandardItemModel(int rows, int columns, QObject * parent = 0)
            model = new QStandardItemModel(SelectedGraphEnergy.size(),ui->customPlot->graphCount()+1,this);
            
            // Attach the model to the view
            ui->PlotDataTableView->setModel(model);
            
            model->setHeaderData(fg, Qt::Horizontal, tr("E(MeV)"));
        }
        //model->setHeaderData(fg+1, Qt::Horizontal, tr((ScoreVariableUnit+ui->customPlot->selectedGraphs()[fg]->name()).toStdString().c_str()));
        model->setHeaderData(fg+1, Qt::Horizontal, tr((ScoreVariableUnit+ui->customPlot->graph(fg)->name()).toStdString().c_str()));
        
        // insert data in the tableView
        for(int row = 0; row < SelectedGraphEnergy.size(); row++)
        {
            if(fg==0){
                model->setData(model->index(row,0),QString::number(SelectedGraphEnergy[row], 'f', 4));
            }
            model->setData(model->index(row,fg+1),QString::number(SelectedGraphScoreVar[row], 'f', 5));
        }
        
    }
    
}
void PlotDialog::SaveTableFromTableViewToFile(){
    
    QString filename="/home/tarik/Desktop/WorkSpace/Projects/QtCreator-build-DoseCalcsCore/AAA1.pdf";
    
#if QT_PRINTSUPPORT_LIB
    QPrinter printer(QPrinter::HighResolution);
    printer.setOutputFileName(filename);
    printer.setPaperSize(QPrinter::A4);
    printer.setOutputFormat(QPrinter::PdfFormat);

    QPainter painter(&printer);
    painter.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing | QPainter::SmoothPixmapTransform);
    ui->PlotDataTableView->render( &painter );
    painter.end();
#endif

}
void PlotDialog::TableViewFillingWithSelectedGraph(){ // called from graphClicked()
    
    QCPDataContainer<QCPGraphData>::iterator Abeg = ui->customPlot->selectedGraphs().first()->data()->begin();
    QCPDataContainer<QCPGraphData>::iterator Aend = ui->customPlot->selectedGraphs().first()->data()->end();
    
    QVector<double> SelectedGraphEnergyForThisGraph;
    QVector<double> SelectedGraphScoreVarForThisGraph;
    
    while(Abeg!=Aend)
    {
        SelectedGraphEnergyForThisGraph.push_back(Abeg->key);
        SelectedGraphScoreVarForThisGraph.push_back(Abeg->value);
        
        Abeg++;
    }
    SelectedGraphEnergy = SelectedGraphEnergyForThisGraph;
    SelectedGraphScoreVar = SelectedGraphScoreVarForThisGraph;
    
    for(int RR = 0; RR< SelectedGraphEnergy.size() ;RR++){
        QTextStream(stdout) << SelectedGraphEnergy[RR] << " - " << SelectedGraphScoreVar[RR] << "\n";
    }
    
    // Create a new model
    // QStandardItemModel(int rows, int columns, QObject * parent = 0)
    model = new QStandardItemModel(SelectedGraphEnergy.size(),2,this);
    
    // Attach the model to the view
    ui->PlotDataTableView->setModel(model);
    
    if(GraphTypeShown == "SC" ){
        model->setHeaderData(0, Qt::Horizontal, tr("E(MeV)"));
        model->setHeaderData(1, Qt::Horizontal, tr(ScoreVariableUnit.toStdString().c_str()));
    }else if(GraphTypeShown == "RSD" ){
        model->setHeaderData(0, Qt::Horizontal, tr("E(MeV)"));
        model->setHeaderData(1, Qt::Horizontal, tr(("RSD\%"+ui->customPlot->selectedGraphs().first()->name()).toStdString().c_str()));
    }else if(GraphTypeShown == "RE" ){
        model->setHeaderData(0, Qt::Horizontal, tr("E(MeV)"));
        model->setHeaderData(1, Qt::Horizontal, tr(("Err\%"+ui->customPlot->selectedGraphs().first()->name()).toStdString().c_str()));
    }else if(GraphTypeShown == "V" ){
        model->setHeaderData(0, Qt::Horizontal, tr((RegionVariableNameWithUnit).toStdString().c_str()));
        model->setHeaderData(1, Qt::Horizontal, tr((ScoreVariableUnit+" for "+ui->customPlot->selectedGraphs().first()->name()).toStdString().c_str()));
    }else if(GraphTypeShown == "MCS" ){
        model->setHeaderData(0, Qt::Horizontal, tr("E(MeV)"));
        model->setHeaderData(1, Qt::Horizontal, tr((ui->customPlot->selectedGraphs().first()->name()+"(cm-1)").toStdString().c_str()));
    }else if(GraphTypeShown == "T" ){
        model->setHeaderData(0, Qt::Horizontal, tr("E(MeV)"));
        model->setHeaderData(1, Qt::Horizontal, tr((ui->customPlot->selectedGraphs().first()->name()+" Time (min)").toStdString().c_str()));
    }else if(GraphTypeShown == "OneE" ){
        model->setHeaderData(0, Qt::Horizontal, tr("Regions Names"));
        model->setHeaderData(1, Qt::Horizontal, tr(ScoreVariableUnit.toStdString().c_str()));
    }else if(GraphTypeShown == "RadioTracerQuantityGraph" ){
        model->setHeaderData(0, Qt::Horizontal, tr("Regions Names"));
        model->setHeaderData(1, Qt::Horizontal, tr((ScoreVariableUnit + " from "+ui->customPlot->selectedGraphs().first()->name()).toStdString().c_str()));
    }else if(GraphTypeShown == "XY" ){
        model->setHeaderData(0, Qt::Horizontal, tr("X(Unit)"));
        model->setHeaderData(1, Qt::Horizontal, tr((ui->customPlot->selectedGraphs().first()->name()).toStdString().c_str()));
    }

    // insert data in the tableView
    for(int row = 0; row < SelectedGraphEnergy.size(); row++)
    {
        model->setData(model->index(row,0),QString::number(SelectedGraphEnergy[row], 'f', 4));
        model->setData(model->index(row,1),QString::number(SelectedGraphScoreVar[row], 'f', 5));
    }
    
}

// SLOT related to the graph interaction
void PlotDialog::titleDoubleClick(QMouseEvent* event)
{
    Q_UNUSED(event)
    if (QCPTextElement *title = qobject_cast<QCPTextElement*>(sender()))
    {
        QTextStream(stdout) << "titleDoubleClick " << "\n";
        
        // Set the plot title by double clicking on it
        bool ok;
        QString newTitle = QInputDialog::getText(this, "Title", "New plot title:", QLineEdit::Normal, title->text(), &ok);
        if (ok)
        {
            title->setText(newTitle);
            ui->customPlot->replot();
        }
    }
}
void PlotDialog::editPlotParameters()
{

    QDialog * d = new QDialog(); d->setWindowTitle("Background, graphs line and axis and text parameters");

    QGridLayout* GraphLayout = new QGridLayout;

    QSpinBox* TextSize = new QSpinBox(); TextSize->setToolTip("Add line size of Text In Plot (the maximum value is 100)");
    TextSize->setMinimum(12);
    TextSize->setMaximum(100);
    TextSize->setSingleStep(1);

    QSpinBox* LineSize = new QSpinBox(); LineSize->setToolTip("Add line size of graphs (the maximum value is 100)");
    LineSize->setMinimum(3);
    LineSize->setMaximum(100);
    LineSize->setSingleStep(1);// Will increment the current value with 1 (if you use up arrow key) (if you use down arrow key => -1)

    TextColorBtn = new QPushButton(); connect(TextColorBtn, SIGNAL(clicked()), this, SLOT(TextColorMapSlot()));
    TextColorBtn->setText("...");
    TextColorBtn->setToolTip("Choose a new color for text in Plot");
    QPalette Pal1(palette());
    TextColor = ui->customPlot->xAxis->tickLabelColor();
    QTextStream(stdout) << "ChosenBackgroundColor "<< ChosenBackgroundColor.name() << "\n";
    Pal1.setColor(QPalette::Button, ChosenBackgroundColor);
    TextColorBtn->setAutoFillBackground(true);
    TextColorBtn->setPalette(Pal1);

    BackgroundColorBtn = new QPushButton(); connect(BackgroundColorBtn, SIGNAL(clicked()), this, SLOT(BackgroundColorMapSlot()));
    BackgroundColorBtn->setText("...");
    BackgroundColorBtn->setToolTip("Choose a new color for Background");
    QPalette Pal2(palette());
    //ChosenBackgroundColor = ui->customPlot->background();
    //QTextStream(stdout) << "ChosenBackgroundColor "<< ChosenBackgroundColor.name() << "\n";
    Pal2.setColor(QPalette::Button, ChosenBackgroundColor);
    BackgroundColorBtn->setAutoFillBackground(true);
    BackgroundColorBtn->setPalette(Pal2);

    QComboBox * lineStyles = new QComboBox(); lineStyles->addItems(lineStylesName);
    lineStyles->setToolTip("Choose a new line style for all Graphs");
    lineStyles->setCurrentText("");

    int ii = 0, jj=0;
    GraphLayout->addWidget(lineStyles, jj,ii,1,1);
    GraphLayout->addWidget(LineSize, jj,++ii,1,1);
    GraphLayout->addWidget(TextSize, jj,++ii,1,1);
    GraphLayout->addWidget(TextColorBtn, jj,++ii,1,1);
    GraphLayout->addWidget(BackgroundColorBtn, jj,++ii,1,1);

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    GraphLayout->addWidget(buttonBox);

    d->setLayout(GraphLayout);

    int result = d->exec();

    if(result == QDialog::Accepted)
    {
        for(int fg = 0; fg < ui->customPlot->graphCount() ; fg++){

            if(lineStyles->currentText()==""){}else{ui->customPlot->graph(fg)->setLineStyle(lineStylesMap[lineStyles->currentText()]);}

            QPen graphPen;
            graphPen.setColor(ui->customPlot->graph(fg)->pen().color());
            graphPen.setWidth(LineSize->value());
            ui->customPlot->graph(fg)->setPen(graphPen);
        }
        ui->customPlot->replot();
    }

}
void PlotDialog::editGraphsStyles()
{

}
void PlotDialog::axisLabelDoubleClick(QCPAxis *axis, QCPAxis::SelectablePart part)
{
    // Set an axis label by double clicking on it
    if (part == QCPAxis::spAxisLabel) // only react when the actual axis label is clicked, not tick label or axis backbone
    {
        QTextStream(stdout) << "axisLabelDoubleClick " << "\n";

        bool ok;
        QString newLabel = QInputDialog::getText(this, "Axis", "New axis label:", QLineEdit::Normal, axis->label(), &ok);
        if (ok)
        {
            axis->setLabel(newLabel);
            ui->customPlot->replot();
        }
    }else if (part == QCPAxis::spTickLabels) // only react when the actual axis label is clicked, not tick label or axis backbone
    {
        QTextStream(stdout) << "axisTickLabelDoubleClick " << "\n";

        bool ok;
        QString newLabel = QInputDialog::getText(this, "Axis", "New Tick label:", QLineEdit::Normal, axis->label(), &ok);
        if (ok)
        {
            /*
            QVector<double> x = QCPAxis::tickVector();
            QVector<double> y = QCPAxis::tickVectorLabels();

            QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
            for(int row = 1; row < x.size(); row++){
                textTicker->addTick(row, First1);
            }
            ui->customPlot->xAxis->setTicker(textTicker);

            ui->customPlot->replot();
*/
        }
    }
}
void PlotDialog::legendDoubleClick(QCPLegend *legend, QCPAbstractLegendItem *item)
{
    // Rename a graph by double clicking on its legend item
    Q_UNUSED(legend)
    if (item) // only react if item was clicked (user could have clicked on border padding of legend where there is no item, then item is 0)
    {
        QTextStream(stdout) << "legendDoubleClick " << "\n";

        QCPPlottableLegendItem *plItem = qobject_cast<QCPPlottableLegendItem*>(item);
        bool ok;
        QString newName = QInputDialog::getText(this, "Legend", "New graph name:", QLineEdit::Normal, plItem->plottable()->name(), &ok);
        if (ok)
        {
            plItem->plottable()->setName(newName);
            ui->customPlot->replot();
        }
    }
    else{
        bool ok;
        int NewSizeVal = QInputDialog::getInt(this, tr("Titles Size"), tr("New Size:"), 10, 0, 100, 1, &ok);
        QFont legendFont = font();

        legendFont.setPointSize(NewSizeVal);
        ui->customPlot->legend->setFont(legendFont);
        ui->customPlot->legend->setSelectedFont(legendFont);
        ui->customPlot->legend->setSelectableParts(QCPLegend::spItems); // legend box shall not be selectable, only legend items
        ui->customPlot->replot();

    }
}
void PlotDialog::selectionChanged()
{
    /*
   normally, axis base line, axis tick labels and axis labels are selectable separately, but we want
   the user only to be able to select the axis as a whole, so we tie the selected states of the tick labels
   and the axis base line together. However, the axis label shall be selectable individually.
   
   The selection state of the left and right axes shall be synchronized as well as the state of the
   bottom and top axes.
   
   Further, we want to synchronize the selection of the graphs with the selection state of the respective
   legend item belonging to that graph. So the user can select a graph by either clicking on the graph itself
   or on its legend item.
  */

    // make top and bottom axes be selected synchronously, and handle axis and tick labels as one selectable object:
    if (ui->customPlot->xAxis->selectedParts().testFlag(QCPAxis::spAxis) || ui->customPlot->xAxis->selectedParts().testFlag(QCPAxis::spTickLabels) ||
            ui->customPlot->xAxis2->selectedParts().testFlag(QCPAxis::spAxis) || ui->customPlot->xAxis2->selectedParts().testFlag(QCPAxis::spTickLabels))
    {
        ui->customPlot->xAxis2->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
        ui->customPlot->xAxis->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
    }
    // make left and right axes be selected synchronously, and handle axis and tick labels as one selectable object:
    if (ui->customPlot->yAxis->selectedParts().testFlag(QCPAxis::spAxis) || ui->customPlot->yAxis->selectedParts().testFlag(QCPAxis::spTickLabels) ||
            ui->customPlot->yAxis2->selectedParts().testFlag(QCPAxis::spAxis) || ui->customPlot->yAxis2->selectedParts().testFlag(QCPAxis::spTickLabels))
    {
        ui->customPlot->yAxis2->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
        ui->customPlot->yAxis->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
    }

    // synchronize selection of graphs with selection of corresponding legend items:
    for (int i=0; i<ui->customPlot->graphCount(); ++i)
    {
        QCPGraph *graph = ui->customPlot->graph(i);
        QCPPlottableLegendItem *item = ui->customPlot->legend->itemWithPlottable(graph);
        if (item->selected() || graph->selected())
        {
            item->setSelected(true);
            graph->setSelection(QCPDataSelection(graph->data()->dataRange()));
        }
    }
}
void PlotDialog::mousePress()
{
    // if an axis is selected, only allow the direction of that axis to be dragged
    // if no axis is selected, both directions may be dragged

    if (ui->customPlot->xAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->customPlot->axisRect()->setRangeDrag(ui->customPlot->xAxis->orientation());
    else if (ui->customPlot->yAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->customPlot->axisRect()->setRangeDrag(ui->customPlot->yAxis->orientation());
    else
        ui->customPlot->axisRect()->setRangeDrag(Qt::Horizontal|Qt::Vertical);
}
void PlotDialog::mouseWheel()
{
    // if an axis is selected, only allow the direction of that axis to be zoomed
    // if no axis is selected, both directions may be zoomed

    if (ui->customPlot->xAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->customPlot->axisRect()->setRangeZoom(ui->customPlot->xAxis->orientation());
    else if (ui->customPlot->yAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->customPlot->axisRect()->setRangeZoom(ui->customPlot->yAxis->orientation());
    else
        ui->customPlot->axisRect()->setRangeZoom(Qt::Horizontal|Qt::Vertical);
}
void PlotDialog::addRandomGraph()
{
    int n = 50; // number of points in graph
    double xScale = (rand()/(double)RAND_MAX + 0.5)*2;
    double yScale = (rand()/(double)RAND_MAX + 0.5)*2;
    double xOffset = (rand()/(double)RAND_MAX - 0.5)*4;
    double yOffset = (rand()/(double)RAND_MAX - 0.5)*10;
    double r1 = (rand()/(double)RAND_MAX - 0.5)*2;
    double r2 = (rand()/(double)RAND_MAX - 0.5)*2;
    double r3 = (rand()/(double)RAND_MAX - 0.5)*2;
    double r4 = (rand()/(double)RAND_MAX - 0.5)*2;
    QVector<double> x(n), y(n);
    for (int i=0; i<n; i++)
    {
        x[i] = (i/(double)n-0.5)*10.0*xScale + xOffset;
        y[i] = (qSin(x[i]*r1*5)*qSin(qCos(x[i]*r2)*r4*3)+r3*qCos(qSin(x[i])*r4*2))*yScale + yOffset;
    }

    ui->customPlot->addGraph();
    ui->customPlot->graph()->setName(QString("New graph %1").arg(ui->customPlot->graphCount()-1));
    ui->customPlot->graph()->setData(x, y);
    ui->customPlot->graph()->setLineStyle((QCPGraph::LineStyle)(rand()%5+1));
    if (rand()%100 > 50)
        ui->customPlot->graph()->setScatterStyle(QCPScatterStyle((QCPScatterStyle::ScatterShape)(rand()%14+1)));
    QPen graphPen;
    graphPen.setColor(QColor(rand()%245+10, rand()%245+10, rand()%245+10));
    graphPen.setWidthF(rand()/(double)RAND_MAX*2+1);
    ui->customPlot->graph()->setPen(graphPen);
    ui->customPlot->replot();
}
void PlotDialog::removeSelectedGraph()
{
    if (ui->customPlot->selectedGraphs().size() > 0)
    {
        ui->customPlot->removeGraph(ui->customPlot->selectedGraphs().first());
        ui->customPlot->replot();
    }
}
void PlotDialog::removeAllGraphs()
{
    ui->customPlot->clearGraphs();
    ui->customPlot->replot();
}
void PlotDialog::contextMenuRequest(QPoint pos)
{
    QMenu *menu = new QMenu(this);
    menu->setAttribute(Qt::WA_DeleteOnClose);

    if (ui->customPlot->legend->selectTest(pos, false) >= 0) // context menu on legend requested
    {
        menu->addAction("Move to top left", this, SLOT(moveLegend()))->setData((int)(Qt::AlignTop|Qt::AlignLeft));
        menu->addAction("Move to top center", this, SLOT(moveLegend()))->setData((int)(Qt::AlignTop|Qt::AlignHCenter));
        menu->addAction("Move to top right", this, SLOT(moveLegend()))->setData((int)(Qt::AlignTop|Qt::AlignRight));
        menu->addAction("Move to bottom right", this, SLOT(moveLegend()))->setData((int)(Qt::AlignBottom|Qt::AlignRight));
        menu->addAction("Move to bottom left", this, SLOT(moveLegend()))->setData((int)(Qt::AlignBottom|Qt::AlignLeft));
    } else {
        menu->addAction("Plot Parameters", this, SLOT(editPlotParameters()));
        menu->addAction("Add random graph", this, SLOT(addRandomGraph()));
        if (ui->customPlot->selectedGraphs().size() > 0){
            menu->addAction("Remove selected graph", this, SLOT(removeSelectedGraph()));
        }
        if (ui->customPlot->graphCount() > 0){
            menu->addAction("Remove all graphs", this, SLOT(removeAllGraphs()));
        }
        menu->addAction("Change x axis scale", this, SLOT(ChangeXaxisScale()));
        menu->addAction("Change x number format", this, SLOT(ChangeXNumberFormat()));
        menu->addAction("Change y axis scale", this, SLOT(ChangeYaxisScale()));
        menu->addAction("Change y number format", this, SLOT(ChangeYNumberFormat()));

    }

    menu->popup(ui->customPlot->mapToGlobal(pos));
}
void PlotDialog::ChangeXaxisScale(){
    // log or linear for y

    QCPAxis::ScaleType ScTy = ui->customPlot->xAxis->scaleType();
    QSharedPointer<QCPAxisTicker> FixedTicker(new QCPAxisTicker);

    if(ScTy == QCPAxis::stLinear){
        QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
        ui->customPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
        ui->customPlot->xAxis->setTicker(logTicker);
        ui->customPlot->xAxis->setNumberFormat(XNumberFormat); // e = exponential, b = beautiful decimal powers
        ui->customPlot->xAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"

    } else if (ScTy == QCPAxis::stLogarithmic){
        QSharedPointer<QCPAxisTicker> FixedTicker(new QCPAxisTicker);
        ui->customPlot->xAxis->setScaleType(QCPAxis::stLinear);
        ui->customPlot->xAxis->setTicker(FixedTicker);
        ui->customPlot->xAxis->setNumberFormat(XNumberFormat); // e = exponential, b = beautiful decimal powers
        ui->customPlot->xAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
    }

    ui->customPlot->rescaleAxes();
    ui->customPlot->replot();

}
void PlotDialog::ChangeYaxisScale(){
    // log or linear for y

    QCPAxis::ScaleType ScTy = ui->customPlot->yAxis->scaleType();
    QSharedPointer<QCPAxisTicker> FixedTicker(new QCPAxisTicker);

    if(ScTy == QCPAxis::stLinear){
        QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
        ui->customPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
        ui->customPlot->yAxis->setTicker(logTicker);
        ui->customPlot->yAxis->setNumberFormat(YNumberFormat); // e = exponential, b = beautiful decimal powers
        ui->customPlot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"

    } else if (ScTy == QCPAxis::stLogarithmic){
        QSharedPointer<QCPAxisTicker> FixedTicker(new QCPAxisTicker);
        ui->customPlot->yAxis->setScaleType(QCPAxis::stLinear);
        ui->customPlot->yAxis->setTicker(FixedTicker);
        ui->customPlot->yAxis->setNumberFormat(YNumberFormat); // e = exponential, b = beautiful decimal powers
        ui->customPlot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
    }

    ui->customPlot->rescaleAxes();
    ui->customPlot->replot();
}
void PlotDialog::ChangeXNumberFormat(){
    bool ok;
    XNumberFormat = QInputDialog::getText(this, "Axis X Number Format", "New format (g gb ebc fb e E f G) :", QLineEdit::Normal, XNumberFormat, &ok);

    if (ok && XNumberFormat != "" && !XNumberFormat.isEmpty() )
    {
        ui->customPlot->xAxis->setNumberFormat(XNumberFormat);
        ui->customPlot->rescaleAxes();
        ui->customPlot->replot();
    }
}
void PlotDialog::ChangeYNumberFormat(){
    bool ok;
    YNumberFormat = QInputDialog::getText(this, "Axis Y Number Format", "New format (g gb ebc fb e E f G) :", QLineEdit::Normal, YNumberFormat, &ok);

    if (ok && YNumberFormat != "" && !YNumberFormat.isEmpty() )
    {
        ui->customPlot->yAxis->setNumberFormat(YNumberFormat);
        ui->customPlot->yAxis->setNumberPrecision(1);

        ui->customPlot->rescaleAxes();
        ui->customPlot->replot();
    }

}

void PlotDialog::moveLegend()
{
    if (QAction* contextAction = qobject_cast<QAction*>(sender())) // make sure this slot is really called by a context menu action, so it carries the data we need
    {
        bool ok;
        int dataInt = contextAction->data().toInt(&ok);
        if (ok)
        {
            ui->customPlot->axisRect()->insetLayout()->setInsetAlignment(0, (Qt::Alignment)dataInt);
            ui->customPlot->replot();
        }
    }
}
void PlotDialog::graphClicked(QCPAbstractPlottable *plottable, int dataIndex)
{
    // since we know we only have QCPGraphs in the plot, we can immediately access interface1D()
    // usually it's better to first check whether interface1D() returns non-zero, and only then use it.

    if(GraphTypeShown == "Rat" || GraphTypeShown == "PB"){
        return;
    }
    double dataValue = plottable->interface1D()->dataMainValue(dataIndex);
    //QString message = QString("Clicked on graph '%1' at data point #%2 with value %3.").arg(plottable->name()).arg(dataIndex).arg(dataValue);
    //ui->statusBar->showMessage(message, 2500);


    TableWidgetFillingWithSelectedGraph();

    QDialog * d = new QDialog(); d->setWindowTitle("Selected Graph Style");

    QGridLayout* GraphLayout = new QGridLayout;

    GraphName = ui->customPlot->selectedGraphs().first()->name();

    QLineEdit * GraphLabel = new QLineEdit(GraphName); GraphLabel->setToolTip("Choose a new label for selected graph");

    QSpinBox* LineSize = new QSpinBox(); LineSize->setToolTip("Add line size of graphs (the maximum value is 100)");
    LineSize->setMinimum(0);
    LineSize->setMaximum(100);
    LineSize->setSingleStep(1);// Will increment the current value with 1 (if you use up arrow key) (if you use down arrow key => -1)
    LineSize->setValue(ui->customPlot->selectedGraphs().first()->pen().width());// Default/begining value

    LineColorBtn = new QPushButton(); connect(LineColorBtn, SIGNAL(clicked()), this, SLOT(LineColorMapSlot()));
    LineColorBtn->setToolTip("Choose a new color for selected graph");
    QPalette Pal(palette());
    LineColor = ui->customPlot->selectedGraphs().first()->pen().color();
    QTextStream(stdout) << "LineColor "<< LineColor.name() << "\n";
    Pal.setColor(QPalette::Button, LineColor);
    LineColorBtn->setAutoFillBackground(true);
    LineColorBtn->setPalette(Pal);

    QComboBox * scatterStyle = new QComboBox(); scatterStyle->addItems(ShapesName);
    scatterStyle->setToolTip("Choose a new scatter point style for selected graph");
    scatterStyle->setCurrentText("");

    QComboBox * lineStyles = new QComboBox(); lineStyles->addItems(lineStylesName);
    lineStyles->setToolTip("Choose a new line style for selected graph");
    lineStyles->setCurrentText("");

    int ii = 0, jj=0;
    GraphLayout->addWidget(GraphLabel, jj,ii,1,1);
    GraphLayout->addWidget(scatterStyle, jj,++ii,1,1);
    GraphLayout->addWidget(lineStyles, jj,++ii,1,1);
    GraphLayout->addWidget(LineSize, jj,++ii,1,1);
    GraphLayout->addWidget(LineColorBtn, jj,++ii,1,1);

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    GraphLayout->addWidget(buttonBox);

    d->setLayout(GraphLayout);

    int result = d->exec();

    if(result == QDialog::Accepted)
    {

        if(scatterStyle->currentText()==""){ ui->customPlot->selectedGraphs().first()->setScatterStyle(ui->customPlot->selectedGraphs().first()->scatterStyle());
        }else{ui->customPlot->selectedGraphs().first()->setScatterStyle(ShapesMap[scatterStyle->currentText()]);}

        if(lineStyles->currentText()==""){ ui->customPlot->selectedGraphs().first()->setLineStyle(ui->customPlot->selectedGraphs().first()->lineStyle());
        }else{ui->customPlot->selectedGraphs().first()->setLineStyle(lineStylesMap[lineStyles->currentText()]);}

        QPen pen;
        pen.setColor(LineColor);
        pen.setWidth(LineSize->value());
        ui->customPlot->selectedGraphs().first()->setPen(pen);

        QTextStream(stdout) << "LineColor "<< LineColor.name() << "\n";

        ui->customPlot->replot();
    }


}
void PlotDialog::LineColorMapSlot()
{

    LineColor = QColorDialog::getColor(LineColor, this );
    if( LineColor.isValid() )
    {

        QTextStream(stdout) << "LineColor "<< LineColor.name() << "\n";

        QPalette Pal(palette());
        Pal.setColor(QPalette::Button, LineColor);
        LineColorBtn->setAutoFillBackground(true);
        LineColorBtn->setPalette(Pal);

        LineColorBtn->repaint();
        ui->customPlot->replot();

    }
}
void PlotDialog::BackgroundColorMapSlot()
{

    ChosenBackgroundColor = QColorDialog::getColor(ChosenBackgroundColor, this );
    if( ChosenBackgroundColor.isValid() )
    {

        QTextStream(stdout) << "ChosenBackgroundColor "<< ChosenBackgroundColor.name() << "\n";

        QLinearGradient gradient(0, 0, 0, 400);
        gradient.setColorAt(0, ChosenBackgroundColor);
        gradient.setColorAt(0.9, ChosenBackgroundColor);
        gradient.setColorAt(1, ChosenBackgroundColor);
        ui->customPlot->setBackground(QBrush(gradient));

        QPalette Pal(palette());
        Pal.setColor(QPalette::Button, ChosenBackgroundColor);
        BackgroundColorBtn->setAutoFillBackground(true);
        BackgroundColorBtn->setPalette(Pal);

        BackgroundColorBtn->repaint();
        ui->customPlot->replot();
    }
}

void PlotDialog::TextColorMapSlot()
{

    TextColor = QColorDialog::getColor(TextColor, this );
    if( TextColor.isValid() )
    {

        QTextStream(stdout) << "TextColor "<< TextColor.name() << "\n";

        ui->customPlot->xAxis->setBasePen(QPen(TextColor));
        ui->customPlot->xAxis->setTickPen(QPen(TextColor));
        //ui->customPlot->xAxis->grid()->setVisible(true);
        ui->customPlot->xAxis->grid()->setPen(QPen(TextColor, 0, Qt::DotLine));
        ui->customPlot->xAxis->setTickLabelColor(TextColor);
        ui->customPlot->xAxis->setLabelColor(TextColor);

        ui->customPlot->yAxis->setBasePen(QPen(TextColor));
        ui->customPlot->yAxis->setTickPen(QPen(TextColor));
        //ui->customPlot->yAxis->grid()->setVisible(true);
        ui->customPlot->yAxis->grid()->setPen(QPen(TextColor, 0, Qt::DotLine));
        ui->customPlot->yAxis->setTickLabelColor(TextColor);
        ui->customPlot->yAxis->setLabelColor(TextColor);

        QPalette Pal(palette());
        Pal.setColor(QPalette::Button, TextColor);
        TextColorBtn->setAutoFillBackground(true);
        TextColorBtn->setPalette(Pal);

        TextColorBtn->repaint();
        ui->customPlot->replot();

    }
}

// used as a slot for the editting the header of tableWidget
// https://stackoverflow.com/questions/29147137/how-do-i-implement-the-ability-to-edit-a-qtablewidgets-vertical-and-horizontal
bool PlotDialog::eventFilter(QObject* object, QEvent* event) {
    if ((object == ui->tableWidget->horizontalHeader()->viewport() || object == ui->tableWidget->verticalHeader()->viewport()) && event->type() == QEvent::MouseButtonDblClick) {
        if (header_editor) { //delete previous editor just in case
            QTextStream(stdout) << 1 << "\n";

            //header_editor->deleteLater();
            QTextStream(stdout) << 1 << "\n";

            header_editor = 0;
        }
        QTextStream(stdout) << 2 << "\n";

        QMouseEvent* e = static_cast<QMouseEvent*>(event);
        QHeaderView* header = static_cast<QHeaderView*>(object->parent());
        int mouse_pos = header->orientation() == Qt::Horizontal ? e->x() : e->y();
        int logical_index = header->logicalIndex(header->visualIndexAt(mouse_pos));
        if (logical_index >= 0) { // if mouse is over an item
            QTextStream(stdout) << 3 << "\n";

            QRect rect; // line edit rect in header's viewport's coordinates
            if (header->orientation() == Qt::Horizontal) {
                rect.setLeft(header->sectionPosition(logical_index));
                rect.setWidth(header->sectionSize(logical_index));
                rect.setTop(0);
                rect.setHeight(header->height());
            } else {
                rect.setTop(header->sectionPosition(logical_index));
                rect.setHeight(header->sectionSize(logical_index));
                rect.setLeft(0);
                rect.setWidth(header->width());
            }
            QTextStream(stdout) << 4 << "\n";

            rect.adjust(1, 1, -1, -1);
            header_editor = new QLineEdit(header->viewport());
            header_editor->move(rect.topLeft());
            header_editor->resize(rect.size());
            header_editor->setFrame(false);
            //get current item text
            QString text = header->model()->headerData(logical_index, header->orientation()).toString();
            header_editor->setText(text);
            header_editor->setFocus();
            editor_index = logical_index; //save for future use
            header_editor->installEventFilter(this); //catch focus out event
            //if user presses Enter it should close editor
            connect(header_editor, SIGNAL(returnPressed()), ui->tableWidget, SLOT(setFocus()));
            header_editor->show();
        }
        return true; // filter out event
    } else if (object == header_editor && event->type() == QEvent::FocusOut) {

        QHeaderView* header = static_cast<QHeaderView*>(header_editor->parentWidget()->parentWidget());
        //save item text
        header->model()->setHeaderData(editor_index, header->orientation(), header_editor->text());
        header_editor->deleteLater(); //safely delete editor
        header_editor = 0;

    }
    return false;
}

// read file of energy distribution and get the Map data to create graph
void PlotDialog::getEnergyDistributionFromFile(QString FilePath){

    QMap <double, double> E_Proba= fileManagerObjectPlot->ReadLinesFromFileWithFirstDoubleValueIndicator(FilePath);

    QString graph_Leg = "Events Energy Distribution";
    graphs_Title = "Energy Distribution";

    create_graphs(E_Proba , graphs_Title , graph_Leg , 1);

    ui->plotMessageLabel->setText("open the energy distribution file, read the energy interval and the events number correspondent generate a graph");

}

void PlotDialog::showResultsOutput(QString text, int level){

    if(level == 4){
        QTextStream(stdout) << text << "\n";
    }
    //ui->outputTextConsole->append(text);

}

void PlotDialog::on_pushButtonOpenSaveGraphsDir_clicked()
{

    if(QFile::exists(UserCurrentResultsDirPath+"/" +GraphsOutDirName)){
        QString command = UserCurrentResultsDirPath+"/" +GraphsOutDirName;
        QProcess process;
        QStringList qsl = {command};
        process.startDetached("nautilus", qsl);
    }else{
        showResultsOutput("Verify that the " + UserCurrentResultsDirPath+"/" +GraphsOutDirName + " existed", 4);
    }
}

void PlotDialog::on_pushButtonGenerateNewGraph_clicked()
{

    GraphTypeShown = "XY";
    GetPlotInputDataAndReadFiles();

    // to set the new graph Title
    graphs_Title = "Graph Title";
    ui->customPlot->plotLayout()->remove(title);
    title = new QCPTextElement(ui->customPlot, graphs_Title , QFont("sans", 17, QFont::Bold));
    ui->customPlot->plotLayout()->addElement(0, 0, title);
    ui->customPlot->plotLayout()->updateLayout();

    ui->customPlot->replot();


}
void PlotDialog::on_pushButtonAddGraph_clicked()
{
    QDialog * d = new QDialog(); d->setWindowTitle("Add Graph Legend");
    QVBoxLayout * vbox = new QVBoxLayout();

    QLineEdit * leg= new QLineEdit(); vbox->addWidget(leg);

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));
    vbox->addWidget(buttonBox);

    d->setLayout(vbox);

    int result = d->exec();

    QString graph_legend = "graph";
    if(result == QDialog::Accepted)
    {
        graph_legend = leg->text();
        int rowNum = ui->tableWidgetForOneGraph->rowCount();

        QMap<double, double> XYmap;
        QVector<double> Xvec ;

        for(int ii = 0; ii < rowNum; ii++){
            Xvec.push_back(ui->tableWidgetForOneGraph->item( ii, 0 )->text().toDouble());
        }

        for(int ii = 0; ii < Xvec.size(); ii++){
            XYmap[Xvec[ii]] = ui->tableWidgetForOneGraph->item( ii, 1 )->text().toDouble();
            QTextStream(stdout) << ii <<" Xvec[ii] " << Xvec[ii] <<" XYmap[Xvec[ii]] " << XYmap[Xvec[ii]] << "\n";
            QTextStream(stdout) << SourceOrgan << " " << TargetOrgan << " E "<< Xvec[ii] << " Total E emmitted " << Xvec[ii]*100000000 << " Target Mass "
                                << RegionParameterValueMap[GeometrySymbol]["Mass"][TargetOrgan] << " SAF " << XYmap[Xvec[ii]] << " AE " << XYmap[Xvec[ii]]*Xvec[ii]*100000000*RegionParameterValueMap[GeometrySymbol]["Mass"][TargetOrgan]
                    << "\n";
        }


        //QTextStream(stdout) <<GraphsData << " " << Compare_Type << "\n";

        create_graphs( XYmap , graph_legend , graphs_Title , 0);

        //ui->customPlot->replot();
        //ui->customPlot->plotLayout()->updateLayout();
    }
}
void PlotDialog::on_pushButtonSaveGraphData_clicked()
{

    QDialog * d = new QDialog(); d->setWindowTitle("Inputs To Save Data");
    QVBoxLayout * vbox = new QVBoxLayout();

    QComboBox * GeoCombobox = new QComboBox(); GeoCombobox->addItems(Geometrylist);
    QComboBox * VarCombobox = new QComboBox(); VarCombobox->addItems(ScoreVarlist);
    QComboBox * ParCombobox = new QComboBox(); ParCombobox->addItems(Particlelist);
    QComboBox * SrcCombobox = new QComboBox(); SrcCombobox->addItems(SourceOrganlist);
    QComboBox * TrgCombobox = new QComboBox(); TrgCombobox->addItems(TargetOrganlist);

    if(GraphTypeShown == "SC" ){

        GeoCombobox->setCurrentText(ui->comboBoxGeometrySymbole->currentText()); vbox->addWidget(GeoCombobox);
        VarCombobox->setCurrentText(ui->plotcomboBoxScoredVariable->currentText()); vbox->addWidget(VarCombobox);
        ParCombobox->setCurrentText(ui->plotcomboBoxParticle->currentText()); vbox->addWidget(ParCombobox);
        SrcCombobox->setCurrentText(ui->comboBoxSourceOrgan->currentText()); vbox->addWidget(SrcCombobox);
        TrgCombobox->setCurrentText(ui->PlotcomboBoxTargetOrgan->currentText()); vbox->addWidget(TrgCombobox);

    }
    else if(GraphTypeShown == "T" ){
        ParCombobox->setCurrentText(ui->plotcomboBoxParticle->currentText()); vbox->addWidget(ParCombobox);
        SrcCombobox->setCurrentText(ui->comboBoxSourceOrgan->currentText()); vbox->addWidget(SrcCombobox);
    }

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));
    vbox->addWidget(buttonBox);

    d->setLayout(vbox);

    int result = d->exec();

    if(result == QDialog::Accepted)
    {
        int rowNum = ui->tableWidgetForOneGraph->rowCount();

        QVector<double> Xvec ;

        for(int ii = 0; ii < rowNum; ii++){
            Xvec.push_back(ui->tableWidgetForOneGraph->item( ii, 0 )->text().toDouble());
        }

        if(GraphTypeShown == "SC" ){
            for(int ii = 0; ii < Xvec.size(); ii++){
                ResultTable[GeoCombobox->currentText()][VarCombobox->currentText()][ParCombobox->currentText()][SrcCombobox->currentText()][TrgCombobox->currentText()][Xvec[ii]] = ui->tableWidgetForOneGraph->item( ii, 1 )->text().toDouble();
                //QTextStream(stdout) << ii <<" Xvec[ii] " << Xvec[ii] <<" XYmap[Xvec[ii]] " << XYmap[Xvec[ii]] << "\n";
            }
        }
        else if(GraphTypeShown == "T" ){
            for(int ii = 0; ii < Xvec.size(); ii++){
                ResultParticleSourceEnergyTime[ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxScoredVariable->currentText()][ParCombobox->currentText()][SrcCombobox->currentText()][Xvec[ii]] = ui->tableWidgetForOneGraph->item( ii, 1 )->text().toDouble();
                //QTextStream(stdout) << ii <<" Xvec[ii] " << Xvec[ii] <<" XYmap[Xvec[ii]] " << XYmap[Xvec[ii]] << "\n";
            }
        }
    }

}
void PlotDialog::on_pushButtonFillTable_clicked()
{
    QString Xtitle = "X(Unit)"; QString Ytitle = "Y(Unit)";

    QDialog * d = new QDialog(); d->setWindowTitle("Table inputs");
    QVBoxLayout * vbox = new QVBoxLayout();

    //QStringList RegionVars=;
    QComboBox * GraphType = new QComboBox(); GraphType->addItems((QStringList() << "Quantity"<<"Time"<<"Other")); vbox->addWidget(GraphType);
    QLineEdit * lineEditXTitle = new QLineEdit(); lineEditXTitle->setText("E(MeV)"); vbox->addWidget(lineEditXTitle);
    QLineEdit * lineEditYTitle = new QLineEdit(); lineEditYTitle->setText("f(E)"); vbox->addWidget(lineEditYTitle);
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

        ui->tableWidgetForOneGraph->clear();
        ui->tableWidgetForOneGraph->setRowCount(0);
        ui->tableWidgetForOneGraph->setColumnCount(0);

        QStringList headers;

        if(GraphType->currentText() == "Quantity"){
            headers.append(tr("E(MeV)"));
            headers.append(tr(ScoreVariableUnit.toStdString().c_str()));
        }else if(GraphType->currentText() == "Time"){
            headers.append(tr("E(MeV)"));
            headers.append(tr("T(min)"));
        }else{
            headers.append(tr(Xtitle.toStdString().c_str()));
            headers.append(tr(Ytitle.toStdString().c_str()));
        }

        ui->tableWidgetForOneGraph->setColumnCount(2);
        ui->tableWidgetForOneGraph->setShowGrid(true);
        ui->tableWidgetForOneGraph->setSelectionMode(QAbstractItemView::SingleSelection);
        ui->tableWidgetForOneGraph->setSelectionBehavior(QAbstractItemView::SelectRows);
        ui->tableWidgetForOneGraph->setHorizontalHeaderLabels(headers);
        ui->tableWidgetForOneGraph->horizontalHeader()->setStretchLastSection(true);

        ui->tableWidgetForOneGraph->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);

        for(int row = 0 ; row < RowNum; row++){
            // Insert row

            ui->tableWidgetForOneGraph->insertRow(row);

            ui->tableWidgetForOneGraph->setItem(row,0, new QTableWidgetItem(""));
            ui->tableWidgetForOneGraph->setItem(row,1, new QTableWidgetItem(""));
        }

        ui->tableWidgetForOneGraph->resizeColumnsToContents();
    }
}

/*
void PlotDialog::on_pushButtonGenerateRegionDataTable_clicked()
{

    GetPlotInputDataAndReadFiles();

    if(!QFile::exists(UserCurrentResultsDirPath+"/"+GraphsOutDirName)){
        QDir* dir = new QDir(UserCurrentResultsDirPath);
        dir->mkdir(GraphsOutDirName);
    }

    QString FileName = UserCurrentResultsDirPath+"/"+ GraphsOutDirName+"/RegionsLatexTables";

    std::ostringstream LatexText;

    LatexText << "=============================== Latex Table for regions names in all the result defined geometries \n\n";


    LatexText << "\\begin{table}[H] \n"
              << "\\centering \n"
              << "\\caption{Regions names implemented in the in " << GeometrySymbol.toStdString() << " phantom geometry and the corresponding volume, mass and density.} \n"
              << "\\begin{tabular}{llll} \n";
    LatexText << "\\hline \n";
    LatexText << "\\textbf{Target region}         & \\textbf{Volume (cm3)}      & \\textbf{Mass (kg)}   & \\textbf{Density (g/cm3)}              ";

    for ( auto Bbeg = RegionParameterValueMap[GeometrySymbol].begin(); Bbeg != RegionParameterValueMap[GeometrySymbol].end(); ++Bbeg  ){

        QString reg = Bbeg.key();
        if(reg == "World"){
            continue;
        }
        LatexText << "       \\\\\\hline\n";
        LatexText << reg.toStdString() <<"     & "<< RegionParameterValueMap[GeometrySymbol]["Volume"][reg] <<"      & "<< RegionParameterValueMap[GeometrySymbol]["Mass"][reg] <<"     & " << RegionParameterValueMap[GeometrySymbol]["Density"][reg];
    }

    LatexText << "\\\\ \\hline\n\\end{tabular} \n";
    LatexText << "\\label{RegionImpleData}\n";
    LatexText << "\\end{table}";
    LatexText << "\n\n\n\n\n";

    std::ofstream outfile(FileName.toStdString() , std::ios::app);
    if(outfile.is_open()){

        //std::cout << "\nCreating file " << FileName.toStdString() << std::"\n" ;
        outfile << LatexText.str();
        outfile.close();
    }
}
void PlotDialog::on_pushButtonGenerateSelfCrossTable_clicked()
{
    GetPlotInputDataAndReadFiles();

    if(!QFile::exists(UserCurrentResultsDirPath+"/"+GraphsOutDirName)){
        QDir* dir = new QDir(UserCurrentResultsDirPath);
        dir->mkdir(GraphsOutDirName);
    }

    for (int SelfOfCross = 1 ; SelfOfCross < 4 ; SelfOfCross++) { // 0 for Self and 1 for Cross

        for (int gg = 0 ; gg < ScoreVarlist.size() ; gg++) {

            QuantitiesToScore = ScoreVarlist[gg];

            QString FileName = UserCurrentResultsDirPath+"/"+ GraphsOutDirName+"/ResRefLatexTables";

            QStringList SourcesName = SourceOrganlist;
            std::ostringstream LatexText;

            int NumberOfEne = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size();
            int NumberOfCol = NumberOfEne + 2;

            int RowNum = 2;
            //std::string ValNm = QuantitiesToScore;
            std::string ValNm = "DoseCalcs";

            if(GraphsData == "Reference_Result"){
                RowNum = 4;
                //ValNm = ReferenceFileName;
            }

            if(SelfOfCross == 1 || SelfOfCross == 3){

                //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << SourcesName[srcInc] << "\n\n";

                LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
                LatexText << "\\usepackage{threeparttable}\n";
                LatexText << "\\usepackage{multirow}\n";
                LatexText << "\n\n";

                LatexText << "\\begin{table}[H] \n"
                          << "\\centering \n";
                if(SelfOfCross == 1){
                    if(RowNum = 4){
                        LatexText << "\\caption{Self-absorption " << QuantitiesToScore.toStdString() << " values of " << ParticleName.toStdString() << " calculated in " << GeometrySymbol.toStdString() << " by DoseCalcs and compared to the " << CompareReferenceName.toStdString() << " reference. } \n";
                    }else{
                        LatexText << "\\caption{Self-absorption " << QuantitiesToScore.toStdString() << " values of " << ParticleName.toStdString() << " calculated in " << GeometrySymbol.toStdString() << " by DoseCalcs. } \n";
                    }
                }
                else if(SelfOfCross == 3){
                    if(RowNum = 4){
                        LatexText << "\\caption{Self- and cross-irradiation " << QuantitiesToScore.toStdString() << " values of " << ParticleName.toStdString() << " for each source-target combination calculated in " << GeometrySymbol.toStdString() << " by DoseCalcs and compared to the " << CompareReferenceName.toStdString() << " reference. } \n";
                    }else{
                        LatexText << "\\caption{Self- and cross-irradiation " << QuantitiesToScore.toStdString() << " values of " << ParticleName.toStdString() << " for each source-target combination calculated in " << GeometrySymbol.toStdString() << " by DoseCalcs. } \n";
                    }
                }
                LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                          << "\\begin{threeparttable}\n"
                          << "\\begin{tabular}{";

                for ( int A = 0; A < NumberOfCol ; A++  )
                {
                    LatexText << "l";
                }
                LatexText << "} \\hline \n";
                if(SelfOfCross == 1){
                    LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                }
                else if(SelfOfCross == 3){
                    LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                }

                LatexText << NumberOfEne << "}{c}{ ";
                LatexText << " \\textbf{ Energies in MeV}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";

                for ( int A = 0; A < NumberOfEne ; A++  )
                {
                    LatexText << "   & \\textbf{" << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][A]<< "}";
                }
            }

            for ( int srcInc = 0; srcInc < SourcesName.size() ; srcInc++  ){

                if(SelfOfCross == 1){

                    LatexText << "       \\\\\\hline \n";
                    LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<< SourcesName[srcInc].toStdString() << "}}       & " <<ValNm  << "                                         ";
                    for ( int B = 0; B < NumberOfEne ; B++  ){
                        LatexText << "& " << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                    }

                    if(RowNum == 4){
                        LatexText << "\\\\ \n                             & "<< CompareReferenceName.toStdString() << "                                        ";

                        for ( int B = 0; B < NumberOfEne ; B++  ){
                            LatexText << "& " << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                        }

                        LatexText << "\\\\ \n                             & " << DiffSym.toStdString() << "(\\%)\\tnote{a}" << "                                        ";

                        for ( int B = 0; B < NumberOfEne ; B++  ){
                            double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                            double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                            double a3;
                            if(ui->radioButtonLRD->isChecked()){
                                if( a1 != 0. && a2 != 0. ){
                                    a3 = std::log(a1/a2)*100;
                                }else{
                                    a3= 100;
                                }
                            }
                            else if(ui->radioButtonRD->isChecked()){
                                if( a1 != 0. && a2 != 0. ){
                                    a3 = (a1-a2/a2)*100;
                                }else{
                                    a3= 100;
                                }
                            }
                            else if(ui->radioButtonRa->isChecked()){
                                if( a2 != 0. ){
                                    a3 = a1/a2;
                                }else{
                                    a3 = 100;
                                }

                            }else{
                                if( a1 != 0. && a2 != 0. ){
                                    a3 = (a1-a2/a2)*100;
                                }else{
                                    a3= 100;
                                }
                            }
                            LatexText << "& " << a3 << "        ";
                        }
                    }
                    LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";

                    for ( int B = 0; B < NumberOfEne ; B++  ){
                        double a1 = ErrorTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                        LatexText << "& " << a1 << "        ";
                    }
                }
                if(SelfOfCross == 2){

                    //std::cout << "\n\n                                                          ========= Creation of Results Latex Table For " << SourcesName[srcInc].toStdString() << " Cross-irradiation Data ========= "<< "\n" << std::"\n";

                    //LatexText << "=============================== Latex Table for  ParticleName " <<  ParticleName.toStdString() << " and Source " << SourcesName[srcInc] << "\n\n";

                    LatexText << "\\begin{table}[H] \n"
                              << "\\centering \n";
                    if(RowNum = 4){
                        LatexText << "\\caption{Cross-irradiation " << QuantitiesToScore.toStdString() << " Values from " << SourcesName[srcInc].toStdString() << " calculated in " << GeometrySymbol.toStdString() << " by DoseCalcs and compared to the " << CompareReferenceName.toStdString() << " reference. } \n";
                    }else{
                        LatexText << "\\caption{Cross-irradiation " << QuantitiesToScore.toStdString() << " Values from " << SourcesName[srcInc].toStdString() << " calculated in " << GeometrySymbol.toStdString() << " by DoseCalcs. } \n";
                    }

                    LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                              << "\\begin{threeparttable}\n"
                              << "\\begin{tabular}{";

                    for ( int A = 0; A < NumberOfCol ; A++  )
                    {
                        LatexText << "l";
                    }
                    LatexText << "} \\hline \n";
                    LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Target region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";

                    LatexText << NumberOfEne << "}{c}{ ";
                    LatexText << " \\textbf{" << SourcesName[srcInc].toStdString() << " source with "<< ParticleName.toStdString() << " energies in MeV}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";

                    for ( int A = 0; A < NumberOfEne ; A++  )
                    {
                        LatexText << "   & \\textbf{" << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][A]<< "}";
                    }

                    for ( int A = 0; A < SourceOrganlist.size() ; A++  ){

                        if(SourceOrganlist[A] == SourcesName[srcInc]){
                            continue;
                        }

                        LatexText << "       \\\\\\hline \n";
                        LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<< SourceOrganlist[A].toStdString() << "}}       & " <<ValNm  << "                                         ";
                        for ( int B = 0; B < NumberOfEne ; B++  ){
                            LatexText << "& " << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourceOrganlist[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                        }

                        if(RowNum == 4){
                            LatexText << "\\\\ \n                             & "<< CompareReferenceName.toStdString() << "                                        ";

                            for ( int B = 0; B < NumberOfEne ; B++  ){
                                LatexText << "& " << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourceOrganlist[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                            }

                            LatexText << "\\\\ \n                             & " << DiffSym.toStdString() << "(\\%)\\tnote{a}" << "                                        ";

                            for ( int B = 0; B < NumberOfEne ; B++  ){
                                double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourceOrganlist[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourceOrganlist[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                double a3;
                                if(ui->radioButtonLRD->isChecked()){
                                    if( a1 != 0. && a2 != 0. ){
                                        a3 = std::log(a1/a2)*100;
                                    }else{
                                        a3= 100;
                                    }
                                }
                                else if(ui->radioButtonRD->isChecked()){
                                    if( a1 != 0. && a2 != 0. ){
                                        a3 = (a1-a2/a2)*100;
                                    }else{
                                        a3= 100;
                                    }
                                }
                                else if(ui->radioButtonRa->isChecked()){
                                    if( a2 != 0. ){
                                        a3 = a1/a2;
                                    }else{
                                        a3 = 100;
                                    }

                                }else{
                                    if( a1 != 0. && a2 != 0. ){
                                        a3 = (a1-a2/a2)*100;
                                    }else{
                                        a3= 100;
                                    }
                                }
                                LatexText << "& " << a3 << "        ";
                            }
                        }

                        LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";

                        for ( int B = 0; B < NumberOfEne ; B++  ){
                            double a1 = ErrorTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourceOrganlist[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                            LatexText << "& " << a1 << "        ";
                        }
                    }

                    LatexText << "\\\\ \\hline\n\\end{tabular} \n";
                    LatexText << "\\begin{tablenotes}\\footnotesize\n";
                    LatexText << "\\item[a] " << DiffExp.toStdString() << "\n";
                    LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                    LatexText << "\\end{tablenotes}\n";
                    LatexText << "\\end{threeparttable}\n";
                    LatexText << "\\end{adjustbox}\n";
                    LatexText << "\\label{tab:CrossFrom"<< SourcesName[srcInc].toStdString() <<"}\n";
                    LatexText << "\\end{table}";
                    LatexText << "\n\n";

                }
                if(SelfOfCross == 3){

                    for ( int A = 0; A < OrganNamesToScore.size() ; A++  ){

                        if(OrganNamesToScore[A] == SourcesName[srcInc]){
                            //continue;
                        }

                        LatexText << "       \\\\\\hline \n";
                        LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<SourcesName[srcInc].toStdString()<<"$\\to$" << OrganNamesToScore[A].toStdString() << "}}       & " <<ValNm  << "                                         ";
                        for ( int B = 0; B < NumberOfEne ; B++  ){
                            LatexText << "& " << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][OrganNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                        }

                        if(RowNum == 4){
                            LatexText << "\\\\ \n                             & "<< CompareReferenceName.toStdString() << "                                        ";

                            for ( int B = 0; B < NumberOfEne ; B++  ){
                                LatexText << "& " << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][OrganNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                            }

                            LatexText << "\\\\ \n                             & " << DiffSym.toStdString() << " \\tnote{a}" << "                                        ";

                            for ( int B = 0; B < NumberOfEne ; B++  ){
                                double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourceOrganlist[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourceOrganlist[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                double a3;
                                if(ui->radioButtonLRD->isChecked()){
                                    if( a1 != 0. && a2 != 0. ){
                                        a3 = std::log(a1/a2)*100;
                                    }else{
                                        a3= 100;
                                    }
                                }
                                else if(ui->radioButtonRD->isChecked()){
                                    if( a1 != 0. && a2 != 0. ){
                                        a3 = (a1-a2/a2)*100;
                                    }else{
                                        a3= 100;
                                    }
                                }
                                else if(ui->radioButtonRa->isChecked()){
                                    if( a2 != 0. ){
                                        a3 = a1/a2;
                                    }else{
                                        a3 = 100;
                                    }

                                }else{
                                    if( a1 != 0. && a2 != 0. ){
                                        a3 = (a1-a2/a2)*100;
                                    }else{
                                        a3= 100;
                                    }
                                }
                                LatexText << "& " << a3 << "        ";
                            }
                        }

                        LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";

                        for ( int B = 0; B < NumberOfEne ; B++  ){
                            double a1 = ErrorTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourceOrganlist[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                            LatexText << "& " << a1 << "        ";
                        }
                    }
                }
            }

            if(SelfOfCross == 1 || SelfOfCross == 3){
                LatexText << "\\\\ \\hline\n\\end{tabular} \n";
                LatexText << "\\begin{tablenotes}\\footnotesize\n";
                LatexText << "\\item[a] " << DiffExp.toStdString() << "\n";
                LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                LatexText << "\\end{tablenotes}\n";
                LatexText << "\\end{threeparttable}\n";
                LatexText << "\\end{adjustbox}\n";
                if(SelfOfCross == 1){LatexText << "\\label{tab:SelfIrr}\n";
                }else if(SelfOfCross == 3){LatexText << "\\label{tab:SelfCrossIrr}\n";}
                LatexText << "\\end{table}";
                LatexText << "\n\n";
            }
            std::ofstream outfile(FileName.toStdString() , std::ios::app);
            if(outfile.is_open()){

                //std::cout << "\nCreating file " << FileName << std::"\n" ;
                outfile << LatexText.str();
                outfile.close();
            }
        }
    }
    GenerateLatexTableResultReferenceForOneEnergy();

    if(ResultQuantityGeometryRadioTracerSourceTargetValues.size() != 0){
        GenerateLatexTableResultReferenceRadioTracerSValues();
    }
}
void PlotDialog::GenerateLatexTableResultReferenceForOneEnergy(){

    if(!QFile::exists(UserCurrentResultsDirPath+"/"+GraphsOutDirName)){
        QDir* dir = new QDir(UserCurrentResultsDirPath);
        dir->mkdir(GraphsOutDirName);
    }

    if(ui->radioButtonLRD->isChecked()){
        DiffSym = "LRD";
        DiffExp = "Logarithmic Relative Difference (\\%)";
    }else if(ui->radioButtonRa->isChecked()){
        DiffSym = "RA";
        DiffExp = "Ratio (\\%)";
    }else{
        DiffSym = "RD";
        DiffExp = "Relative Difference (\\%)";
    }

    for (int ss = 0 ; ss < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; ss++) {

        QString FileName = UserCurrentResultsDirPath+"/"+ GraphsOutDirName+"/ResRefLatexTables";

        QVector<QString> SourcesName = OrganNamesToScore;
        std::ostringstream LatexText;

        int NumberOfTargets = OrganNamesToScore.size();
        int NumberOfCol = NumberOfTargets + 2;

        int RowNum = 2;
        //std::string ValNm = QuantitiesToScore;
        std::string ValNm = "DoseCalcs";

        if(GraphsData == "Reference_Result"){
            RowNum = 4;
        }

        //std::cout << "\n\n                                                          ========= Creation of Latex Table For " << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss] << " Self-absorption and Cross-irradiation Data ========= "<< "\n" << std::"\n";

        //LatexText << "=============================== Latex Table for     ParticleName " <<     ParticleName << " and Source " << SourcesName[A] << "\n\n";

        LatexText << "\\begin{table}[H] \n"
                  << "\\centering \n";

        if(RowNum = 4){
            LatexText << "\\caption{Cross-irradiation " << QuantitiesToScore.toStdString() << " values for " <<     ParticleName.toStdString() << " of energy " << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss] << " MeV calculated in " << GeometrySymbol.toStdString() << " by DoseCalcs and compared to the " << CompareReferenceName.toStdString() << " reference. } \n";
        }else{
            LatexText << "\\caption{Cross-irradiation " << QuantitiesToScore.toStdString() << " values for " <<     ParticleName.toStdString() << " of energy " << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss] << " MeV calculated in " << GeometrySymbol.toStdString() << " by DoseCalcs. } \n";
        }
        LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                  << "\\begin{threeparttable}\n"
                  << "\\begin{tabular}{";

        for ( int A = 0; A < NumberOfCol ; A++  )
        {
            LatexText << "l";
        }
        LatexText << "} \\hline \n";
        LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";

        LatexText << NumberOfTargets << "}{c}{ ";
        LatexText << " \\textbf{Target region}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";

        for ( int A = 0; A < NumberOfTargets ; A++  )
        {
            LatexText << "   & \\textbf{" << OrganNamesToScore[A].toStdString()<< "}";
        }

        for ( int A = 0; A < SourcesName.size() ; A++  ){

            LatexText << "       \\\\\\hline \n";
            LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<< OrganNamesToScore[A].toStdString() << "}}       & " <<ValNm  << "                                         ";
            for ( int B = 0; B < NumberOfTargets ; B++  ){
                LatexText << "& " << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[A]][OrganNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]] << "        ";
                //QTextStream(stdout) << " " << QuantitiesToScore << " " <<     ParticleName << " " << SourcesName[A] << " " << OrganNamesToScore[B] << " " << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss] << " " << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][OrganNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]] << "\n" ;
            }

            if(RowNum == 4){
                LatexText << "\\\\ \n                             & "<< CompareReferenceName.toStdString() << "                                        ";

                for ( int B = 0; B < NumberOfTargets ; B++  ){
                    LatexText << "& " << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[A]][OrganNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]] << "        ";
                }

                LatexText << "\\\\ \n                             & " << DiffSym.toStdString() << "(\\%)\\tnote{a}" << "                                        ";

                for ( int B = 0; B < NumberOfTargets ; B++  ){
                    double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[A]][OrganNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]];
                    double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[A]][OrganNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]];
                    double a3;
                    if(ui->radioButtonLRD->isChecked()){
                        if( a1 != 0. && a2 != 0. ){
                            a3 = std::log(a1/a2)*100;
                        }else{
                            a3= 100;
                        }
                    }
                    else if(ui->radioButtonRD->isChecked()){
                        if( a1 != 0. && a2 != 0. ){
                            a3 = (a1-a2/a2)*100;
                        }else{
                            a3= 100;
                        }
                    }
                    else if(ui->radioButtonRa->isChecked()){
                        if( a2 != 0. ){
                            a3 = a1/a2;
                        }else{
                            a3 = 100;
                        }
                    }else{
                        if( a1 != 0. && a2 != 0. ){
                            a3 = (a1-a2/a2)*100;
                        }else{
                            a3= 100;
                        }
                    }
                    LatexText << "& " << a3 << "        ";
                }
            }

            LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";

            for ( int B = 0; B < NumberOfTargets ; B++  ){
                double a1 = ErrorTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[A]][OrganNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]];
                LatexText << "& " << a1 << "        ";
            }
        }

        LatexText << "\\\\ \\hline\n\\end{tabular} \n";
        LatexText << "\\begin{tablenotes}\\footnotesize\n";
        LatexText << "\\item[a] " << DiffExp.toStdString() << "\n";
        LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
        LatexText << "\\end{tablenotes}\n";
        LatexText << "\\end{threeparttable}\n";
        LatexText << "\\end{adjustbox}\n";
        LatexText << "\\label{"<< QuantitiesToScore.toStdString() << "_" <<ParticleName.toStdString() <<"_SelfCross_"<< ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss] <<"}\n";
        LatexText << "\\end{table}";
        LatexText << "\n\n";

        std::ofstream outfile(FileName.toStdString() , std::ios::app);
        if(outfile.is_open()){
            //std::cout << "\nCreating file " << FileName << std::"\n" ;
            outfile << LatexText.str();
            outfile.close();
        }
    }
}
void PlotDialog::GenerateLatexTableResultReferenceRadioTracerSValues(){

    if(!QFile::exists(UserCurrentResultsDirPath+"/"+GraphsOutDirName)){
        QDir* dir = new QDir(UserCurrentResultsDirPath);
        dir->mkdir(GraphsOutDirName);
    }

    for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  ){

        QString FileName = UserCurrentResultsDirPath+"/"+ GraphsOutDirName+"/ResRefLatexTables";

        std::ostringstream LatexText;

        int NumberOfTargets = OrganNamesToScore.size();
        int NumberOfCol = NumberOfTargets + 2;

        int RowNum = 2;

        std::string ValNm = "DoseCalcs";

        if(GraphsData == "Reference_Result"){
            RowNum = 3;
        }

        //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << SourcesName[A] << "\n\n";

        LatexText << "\\begin{table}[H] \n"
                  << "\\centering \n";
        if(RowNum = 3){
            LatexText << "\\caption{" << " S-values for " << Abeg.key().toStdString() << " calculated in " << GeometrySymbol.toStdString() << " by DoseCalcs and compared to the " << CompareReferenceName.toStdString() << " reference. } \n";
        }
        else{
            LatexText << "\\caption{" << " S-values for " << Abeg.key().toStdString() << " calculated in " << GeometrySymbol.toStdString() << " by DoseCalcs. } \n";
        }
        LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                  << "\\begin{threeparttable}\n"
                  << "\\begin{tabular}{";

        for ( int A = 0; A < NumberOfCol ; A++  )
        {
            LatexText << "l";
        }
        LatexText << "} \\hline \n";
        LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";

        LatexText << NumberOfTargets << "}{c}{ ";
        LatexText << " \\textbf{Target region}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";

        for ( int A = 0; A < NumberOfTargets ; A++  )
        {
            LatexText << "   & \\textbf{" << OrganNamesToScore[A].toStdString()<< "}";
        }


        for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  ){

            LatexText << "       \\\\\\hline \n";
            LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<< Bbeg.key().toStdString()<< "}}       & " <<ValNm  << "                                         ";

            //for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg  ){}

            for ( int B = 0; B < NumberOfTargets ; B++  ){
                LatexText << "& " << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][Abeg.key()][Bbeg.key()][OrganNamesToScore[B]] << "        ";
                //QTextStream(stdout) << " " << QuantitiesToScore << " " << ParticleName << " " << Bbeg.key() << " " << OrganNamesToScore[B] << " " << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss] << " " << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][OrganNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]] << " " << G4"\n" ;
            }

            if(RowNum == 3){
                LatexText << "\\\\ \n                             & "<< CompareReferenceName.toStdString() << "            ";

                for ( int B = 0; B < NumberOfTargets ; B++  ){
                    LatexText << "& " << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][Abeg.key()][Bbeg.key()][OrganNamesToScore[B]] << "        ";
                }

                LatexText << "\\\\ \n                             & "<< DiffSym.toStdString() <<" \\tnote{a}" << "                                        ";

                for ( int B = 0; B < NumberOfTargets ; B++  ){
                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][Abeg.key()][Bbeg.key()][OrganNamesToScore[B]];
                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][Abeg.key()][Bbeg.key()][OrganNamesToScore[B]];
                    double a3;
                    if(ui->radioButtonLRD->isChecked()){
                        if( a1 != 0. && a2 != 0. ){
                            a3 = std::log(a1/a2)*100;
                        }else{
                            a3= 100;
                        }
                    }
                    else if(ui->radioButtonRD->isChecked()){
                        if( a1 != 0. && a2 != 0. ){
                            a3 = (a1-a2/a2)*100;
                        }else{
                            a3= 100;
                        }
                    }
                    else if(ui->radioButtonRa->isChecked()){
                        if( a2 != 0. ){
                            a3 = a1/a2;
                        }else{
                            a3 = 100;
                        }
                    }else{
                        if( a1 != 0. && a2 != 0. ){
                            a3 = (a1-a2/a2)*100;
                        }else{
                            a3= 100;
                        }
                    }

                    LatexText << "& " << a3 << "        ";
                }
            }

        }

        LatexText << "\\\\ \\hline\n\\end{tabular} \n";
        LatexText << "\\begin{tablenotes}\\footnotesize\n";
        LatexText << "\\item[a] " << DiffExp.toStdString() << " \n";
        LatexText << "\\end{tablenotes}\n";
        LatexText << "\\end{threeparttable}\n";
        LatexText << "\\end{adjustbox}\n";
        LatexText << "\\label{"<< "S-values_" <<Abeg.key().toStdString() <<"_SelfCross" <<"}\n";
        LatexText << "\\end{table}";
        LatexText << "\n\n";

        std::ofstream outfile(FileName.toStdString() , std::ios::app);
        if(outfile.is_open()){

            outfile << LatexText.str();
            outfile.close();
        }
    }
}
void PlotDialog::on_pushButtonGenerateCSV_clicked()
{

    if(!QFile::exists(UserCurrentResultsDirPath+"/"+GraphsOutDirName)){
        QDir* dir = new QDir(UserCurrentResultsDirPath);
        dir->mkdir(GraphsOutDirName);
    }

    QString Quantity;
    QString Geometry;
    QString PARTICLE_NAME;
    QString Source_ORG;
    QString Target_ORG;

    // file for each quantity (a table for each geometry and particle and source and SD)

    for ( auto Obeg = ResultTable.begin(); Obeg != ResultTable.end(); ++Obeg  )
    {
        Quantity = Obeg.key();

        QString FileName = UserCurrentResultsDirPath+"/"+GraphsOutDirName +"/DoseCalcs_Particles_"+Quantity+".csv";
        std::ofstream outfile(FileName.toStdString() , std::ios::binary);
        if(outfile.is_open()){
            //std::cout << "\nCreating file " << FileName << std::endl ;
        }


        std::ostringstream LatexText;
        LatexText << Quantity.toStdString() <<"("<< QuantityUnit[Quantity].toStdString() <<")"<< "\n";

        for ( auto Mbeg = Obeg.value().begin(); Mbeg != Obeg.value().end(); ++Mbeg  )
        {
            Geometry = Mbeg.key();


            for ( auto Abeg = Mbeg.value().begin(); Abeg != Mbeg.value().end(); ++Abeg  )
            {
                PARTICLE_NAME = Abeg.key();

                for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  )
                {
                    Source_ORG = Bbeg.key();

                    LatexText << "\n" << Geometry.toStdString() << "," << PARTICLE_NAME.toStdString() << ","<< Source_ORG.toStdString() << ","<<  "\n";

                    LatexText << "Target/Energies-SDev" << ",";
                    for (int gg = 0 ; gg < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; gg++) {
                        double val = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg];
                        LatexText << val  << "," << " ,";
                    }
                    LatexText << "\n";

                    for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg  )
                    {

                        Target_ORG = Cbeg.key();
                        LatexText << Target_ORG.toStdString()  << ",";
                        for (int gg = 0 ; gg < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; gg++) {
                            double val = ResultTable[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg]];
                            //double error = RelativeStandartDeviationPerCent[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg]];
                            double error = StandartDeviation[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg]];
                            if(val == 0 || __isinf(val) || __isnan(val)){
                                LatexText << " ," << " ,";
                            }else{
                                LatexText << val << "," << error << ",";
                            }
                        }
                        LatexText << "\n";

                    }

                    LatexText << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";

                }
            }

            LatexText << "\n\n*****************************************************************************************************************************************************************************************************************************************************************\n";

        }
        if(outfile.is_open()){
            outfile << LatexText.str();
            outfile.close();
        }
    }

    // file for each quantity (a table for each geometry and particle and source and SD) with reference

    for ( auto Obeg = ResultTable.begin(); Obeg != ResultTable.end(); ++Obeg  )
    {
        Quantity = Obeg.key();

        if(ReferenceTable[Quantity].size() == 0){continue;}

        QString FileName = UserCurrentResultsDirPath+"/"+GraphsOutDirName + "/DoseCalcs-" + CompareReferenceName + "_Particles_"+Quantity+".csv";
        std::ofstream outfile(FileName.toStdString() , std::ios::binary);
        if(outfile.is_open()){
            //std::cout << "\nCreating file " << FileName << std::endl ;
        }


        std::ostringstream LatexText;
        LatexText << Quantity.toStdString() <<"("<< QuantityUnit[Quantity].toStdString() <<")"<< "\n";

        for ( auto Mbeg = Obeg.value().begin(); Mbeg != Obeg.value().end(); ++Mbeg  )
        {
            Geometry = Mbeg.key();

            if(ReferenceTable[Quantity][Geometry].size() == 0){continue;}

            for ( auto Abeg = Mbeg.value().begin(); Abeg != Mbeg.value().end(); ++Abeg  )
            {
                PARTICLE_NAME = Abeg.key();

                if(ReferenceTable[Quantity][Geometry][PARTICLE_NAME].size() == 0){continue;}

                for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  )
                {
                    Source_ORG = Bbeg.key();

                    if(ReferenceTable[Quantity][Geometry][PARTICLE_NAME][Source_ORG].size() == 0){continue;}

                    LatexText << "\n" << Geometry.toStdString() << "," << PARTICLE_NAME.toStdString() << ","<< Source_ORG.toStdString() << ","<<  "\n";

                    LatexText << "Target/Energies(DoseCalcs-"<< CompareReferenceName.toStdString() <<"-Ratio)" << ",";
                    for (int gg = 0 ; gg < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; gg++) {
                        double val = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg];
                        LatexText << val  << "," << " ,"<< " Ratio,";
                    }
                    LatexText << "\n";

                    for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg  )
                    {

                        Target_ORG = Cbeg.key();
                        LatexText << Target_ORG.toStdString()  << ",";
                        for (int gg = 0 ; gg < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; gg++) {
                            double val = ResultTable[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg]];
                            double val2 = ReferenceTable[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg]];
                            if(val == MinValForLog || val == 0 || __isinf(val) || __isnan(val)){

                                LatexText << " ,";
                                if(val2 == MinValForLog || val2 == 0 || __isinf(val2) || __isnan(val2))
                                {
                                    LatexText << " ," << " ,";

                                }else{
                                    double ratio = val/val2;
                                    LatexText << val2 << " ,"  << ratio << " ,";
                                }

                            }else{

                                LatexText << val << ",";
                                if(val2 == MinValForLog || val2 == 0 || __isinf(val2) || __isnan(val2))
                                {
                                    LatexText << " ," << " ,";

                                }else{
                                    double ratio = val/val2;
                                    LatexText << val2 << " ,"  << ratio << " ,";
                                }
                            }
                        }
                        LatexText << "\n";

                    }

                    LatexText << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";

                }
            }

            LatexText << "\n\n*****************************************************************************************************************************************************************************************************************************************************************\n";

        }
        if(outfile.is_open()){
            outfile << LatexText.str();
            outfile.close();
        }
    }

    // file for each quantity (a table for each geometry and radiotracer and source and SD)

    for ( auto Obeg = ResultQuantityGeometryRadioTracerSourceTargetValues.begin(); Obeg != ResultQuantityGeometryRadioTracerSourceTargetValues.end(); ++Obeg  )
    {
        Quantity = Obeg.key();

        QString FileName = UserCurrentResultsDirPath+"/"+GraphsOutDirName +"/DoseCalcs_RadioNucleides_"+Quantity+".csv";
        std::ofstream outfile(FileName.toStdString() , std::ios::binary);
        if(outfile.is_open()){
            //std::cout << "\nCreating file " << FileName << std::endl ;
        }


        std::ostringstream LatexText;
        LatexText << Quantity.toStdString() <<"("<< QuantityUnit[Quantity].toStdString() <<")"<< "\n";

        for ( auto Mbeg = Obeg.value().begin(); Mbeg != Obeg.value().end(); ++Mbeg  )
        {
            Geometry = Mbeg.key();


            for ( auto Abeg = Mbeg.value().begin(); Abeg != Mbeg.value().end(); ++Abeg  )
            {
                PARTICLE_NAME = Abeg.key();

                for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  )
                {
                    Source_ORG = Bbeg.key();

                    LatexText << "\n" << Geometry.toStdString() << "," << PARTICLE_NAME.toStdString() << ","<< Source_ORG.toStdString() << ","<<  "\n";

                    LatexText << "Target" << "," << Quantity.toStdString() << "\n";

                    for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg  )
                    {

                        Target_ORG = Cbeg.key();
                        LatexText << Target_ORG.toStdString()  << ",";
                        double val = ResultQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG];
                        if(val == 0 || __isinf(val) || __isnan(val)){
                            LatexText << "," ;
                        }else{
                            LatexText << val << "," ;
                        }

                        LatexText << "\n";

                    }

                    LatexText << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";

                }
            }

            LatexText << "\n\n*****************************************************************************************************************************************************************************************************************************************************************\n";

        }
        if(outfile.is_open()){
            outfile << LatexText.str();
            outfile.close();
        }
    }

    // file for each quantity (a table for each geometry and radiotracer and source and SD) with reference

    for ( auto Obeg = ResultQuantityGeometryRadioTracerSourceTargetValues.begin(); Obeg != ResultQuantityGeometryRadioTracerSourceTargetValues.end(); ++Obeg  )
    {
        Quantity = Obeg.key();

        if(ReferenceQuantityGeometryRadioTracerSourceTargetValues[Quantity].size() == 0){continue;}

        QString FileName = UserCurrentResultsDirPath+"/"+GraphsOutDirName + "/DoseCalcs-" + CompareReferenceName + "_Radionucleide_"+Quantity+".csv";
        std::ofstream outfile(FileName.toStdString() , std::ios::binary);
        if(outfile.is_open()){
            //std::cout << "\nCreating file " << FileName << std::endl ;
        }

        std::ostringstream LatexText;
        LatexText << Quantity.toStdString() <<"("<< QuantityUnit[Quantity].toStdString() <<")"<< "\n";

        for ( auto Mbeg = Obeg.value().begin(); Mbeg != Obeg.value().end(); ++Mbeg  )
        {
            Geometry = Mbeg.key();

            if(ReferenceQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry].size() == 0){continue;}

            for ( auto Abeg = Mbeg.value().begin(); Abeg != Mbeg.value().end(); ++Abeg  )
            {
                PARTICLE_NAME = Abeg.key();

                if(ReferenceQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][PARTICLE_NAME].size() == 0){continue;}

                for ( auto Bbeg = Abeg.value().begin(); Bbeg != Abeg.value().end(); ++Bbeg  )
                {
                    Source_ORG = Bbeg.key();

                    if(ReferenceQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][PARTICLE_NAME][Source_ORG].size() == 0){continue;}

                    LatexText << "\n" << Geometry.toStdString() << "," << PARTICLE_NAME.toStdString() << ","<< Source_ORG.toStdString() << ","<<  "\n";

                    LatexText << "Target" << "," << Quantity.toStdString() << "-DoseCalcs,"<< Quantity.toStdString() <<"-"<< CompareReferenceName.toStdString() <<", Ratio" << "\n";

                    for ( auto Cbeg = Bbeg.value().begin(); Cbeg != Bbeg.value().end(); ++Cbeg  )
                    {

                        Target_ORG = Cbeg.key();
                        LatexText << Target_ORG.toStdString() << ",";
                        double val = ResultQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG];
                        double val2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG];
                        if(val == MinValForLog || val == 0 || __isinf(val) || __isnan(val)){

                            LatexText << " ,";
                            if(val2 == MinValForLog || val2 == 0 || __isinf(val2) || __isnan(val2))
                            {
                                LatexText << " ," << " ,";

                            }else{
                                double ratio = val/val2;
                                LatexText << val2 << " ,"  << ratio << " ,";
                            }

                        }else{

                            LatexText << val << ",";
                            if(val2 == MinValForLog || val2 == 0 || __isinf(val2) || __isnan(val2))
                            {
                                LatexText << " ," << " ,";

                            }else{
                                double ratio = val/val2;
                                LatexText << val2 << " ,"  << ratio << " ,";
                            }
                        }
                        LatexText << "\n";
                    }

                    LatexText << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";

                }
            }

            LatexText << "\n\n*****************************************************************************************************************************************************************************************************************************************************************\n";

        }
        if(outfile.is_open()){
            outfile << LatexText.str();
            outfile.close();
        }
    }

}

void PlotDialog::on_pushButtonSelfCrossLatex_clicked()
{
    QString command = UserCurrentResultsDirPath+"/"+ GraphsOutDirName+"/ResRefLatexTables";
    QProcess process;
    QStringList qsl = {command};
    process.startDetached("gedit", qsl);
}
void PlotDialog::on_pushButtonRegionLatex_clicked()
{
    QString command = UserCurrentResultsDirPath+"/"+ GraphsOutDirName+"/RegionsLatexTables";
    QProcess process;
    QStringList qsl = {command};
    process.startDetached("gedit", qsl);
}
void PlotDialog::on_pushButtonCrossSectionLatex_clicked()
{    
    QString command = UserCurrentResultsDirPath+"/CrossSectionData";
    QProcess process;
    QStringList qsl = {command};
    process.startDetached("gedit", qsl);
}
void PlotDialog::on_pushButtonOpenCSV_clicked()
{
    QString command = UserCurrentResultsDirPath+"/"+GraphsOutDirName;
    QProcess process;
    QStringList qsl = {command};
    process.startDetached("nautilus", qsl);
}
*/

void PlotDialog::on_checkBoxDiffForRadiotracerOr_stateChanged(int arg1)
{
    if(ui->checkBoxDiffForRadiotracerOr->isChecked()){
        ui->checkBoxDiffForRadiotracerOr->setToolTip("Click on the \"Relative Diff Graph\" to calculate and plot the ratio between the radiotracer DoseCalcs results and the given reference");
    }
    else{
        ui->checkBoxDiffForRadiotracerOr->setToolTip("Click on the \"Relative Diff Graph\" to calculate and plot the ratio, relative difference or logarithmic relative difference between the particle DoseCalcs results and the given reference");
    }
}

void PlotDialog::on_plotcomboBoxScoredVariable_activated(const QString &arg1)
{
    ui->comboBoxGeometrySymbole->clear();    ui->comboBoxGeometrySymbole->update();
    for ( auto Abeg = QuantityGeometryParticleSourceTargetsNames[ui->plotcomboBoxScoredVariable->currentText()].begin(); Abeg != QuantityGeometryParticleSourceTargetsNames[ui->plotcomboBoxScoredVariable->currentText()].end(); ++Abeg  ){
        QString GeomSymb = Abeg.key();
        if(GeomSymb != "" || !GeomSymb.isEmpty()){ui->comboBoxGeometrySymbole->addItem(GeomSymb);}
        showResultsOutput(GeomSymb, 4);
    }
    ui->comboBoxGeometrySymbole->update();

    ui->plotComboBoxEnergies->clear();
    for( int ss = 0; ss < ResEnergies[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()].size(); ss++){
        ui->plotComboBoxEnergies->addItem(QString::number(ResEnergies[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()][ss]));
    }
}
void PlotDialog::on_comboBoxGeometrySymbole_activated(const QString &arg1)
{
    ui->plotcomboBoxParticle->clear();     ui->plotcomboBoxParticle->update();
    for ( auto Abeg = QuantityGeometryParticleSourceTargetsNames[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()].begin(); Abeg != QuantityGeometryParticleSourceTargetsNames[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()].end(); ++Abeg  ){
        QString PARTICLE_NAME = Abeg.key();
        if(PARTICLE_NAME != "" || !PARTICLE_NAME.isEmpty()){ui->plotcomboBoxParticle->addItem(PARTICLE_NAME);}
        showResultsOutput(PARTICLE_NAME, 4);
    }
    ui->plotcomboBoxParticle->update();

    ui->plotComboBoxEnergies->clear();
    for( int ss = 0; ss < ResEnergies[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()].size(); ss++){
        ui->plotComboBoxEnergies->addItem(QString::number(ResEnergies[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()][ss]));
    }
}
void PlotDialog::on_plotcomboBoxParticle_activated(const QString &arg1)
{
    ui->comboBoxSourceOrgan->clear(); ui->comboBoxSourceOrgan->update();
    //QTextStream(stdout)  << " \n\n\n\n Particle " << ui->plotcomboBoxParticle->currentText() << " QuantityGeometryParticleSourceTargetsNames size : " << QuantityGeometryParticleSourceTargetsNames[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()].size() << "\n";
    for ( auto Bbeg = QuantityGeometryParticleSourceTargetsNames[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()].begin(); Bbeg != QuantityGeometryParticleSourceTargetsNames[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()].end(); ++Bbeg  ){
        QString Source_ORG = Bbeg.key();
        showResultsOutput(Source_ORG, 4);
        if(Source_ORG != "" || !Source_ORG.isEmpty()){ui->comboBoxSourceOrgan->addItem(Source_ORG);}
    }
    ui->comboBoxSourceOrgan->update();

    ui->plotComboBoxEnergies->clear();
    for( int ss = 0; ss < ResEnergies[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()].size(); ss++){
        ui->plotComboBoxEnergies->addItem(QString::number(ResEnergies[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()][ss]));
    }
}
void PlotDialog::on_comboBoxSourceOrgan_activated(const QString &arg1)
{
    ui->PlotcomboBoxTargetOrgan->clear(); ui->PlotcomboBoxTargetOrgan->update();
    //QTextStream(stdout)  << " \n\n\n\n Source " << ui->PlotcomboBoxTargetOrgan->currentText() << " QuantityGeometryParticleSourceTargetsNames size : " << QuantityGeometryParticleSourceTargetsNames[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()].size() << "\n";
    //showResultsOutput(ui->plotcomboBoxParticle->currentText(), 4);
    for ( auto Bbeg = QuantityGeometryParticleSourceTargetsNames[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()][ui->comboBoxSourceOrgan->currentText()].begin(); Bbeg != QuantityGeometryParticleSourceTargetsNames[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()][ui->comboBoxSourceOrgan->currentText()].end(); ++Bbeg  ){
        QString Target_ORG = Bbeg.key();
        showResultsOutput(Target_ORG, 4);
        if(Target_ORG != "" || !Target_ORG.isEmpty()){ui->PlotcomboBoxTargetOrgan->addItem(Target_ORG);}
    }
    ui->PlotcomboBoxTargetOrgan->update();

    ui->plotComboBoxEnergies->clear();
    for( int ss = 0; ss < ResEnergies[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()].size(); ss++){
        ui->plotComboBoxEnergies->addItem(QString::number(ResEnergies[ui->plotcomboBoxScoredVariable->currentText()][ui->comboBoxGeometrySymbole->currentText()][ui->plotcomboBoxParticle->currentText()][ss]));
    }
}

void PlotDialog::on_pushButtonOpenInRoot_clicked()
{
    QMap<QString,QMap<double,double>> DataForROOTMap;
    QVector<double> EnergiesVec;
    QVector<QString> XLabels;

    int graphtype = 0;
    if(GraphTypeShown == "OneE" || GraphTypeShown == "RadioTracerQuantityGraph"){

        graphtype = 1;
    }else{
        graphtype = 0;
    }

    for(int fg = 0; fg < ui->customPlot->graphCount() ; fg++){

        for ( auto it = ui->customPlot->graph(fg)->data()->begin(); it != ui->customPlot->graph(fg)->data()->end(); ++it  ){

            QString GraphName = ui->customPlot->graph(fg)->name();
            GraphName.replace(" ","_");

            DataForROOTMap[GraphName][it->key] = it->value;
            QTextStream(stdout) << fg << " " << ui->customPlot->graph(fg)->name() << " it->key " << it->key << " it->value "<< it->value << "\n";

            bool isin = false;
            for(int i=0; i<EnergiesVec.size(); i++){
                if(EnergiesVec[i] == it->key){
                    isin = true; break;
                }
            }
            if(isin == false){
                EnergiesVec.push_back(it->key);
            }
        }
    }

    if(graphtype == 1){

        for(int a=0; a<EnergiesVec.size(); a++){
            double vv = EnergiesVec[a];
            for(int i=0; i<xticks.size(); i++){
                double nn = xticks[i];
                if(nn == vv){
                    XLabels.push_back(xlabels[i]);
                }
            }
        }
    }

    QString FileName = UserCurrentResultsDirPath+"/"+GraphsOutDirName + "/ROOTGraphData";
    std::ofstream outfile(FileName.toStdString() , std::ios::binary);
    if(outfile.is_open()){
        //std::cout << "\nCreating file " << FileName << std::endl ;
    }
    std::ostringstream LatexText;


    QString titletext = title->text(); titletext.replace(" ","_");
    QString xlbl = ui->customPlot->xAxis->label(); xlbl.replace(" ","_");
    QString ylbl = ui->customPlot->yAxis->label(); ylbl.replace(" ","_");

    LatexText  << graphtype << " " << titletext.toStdString() << " " << xlbl.toStdString() << " " << ylbl.toStdString() << " yes " << " yes " << " yes " << " RightTop " << " no \n";

    LatexText << "GraphName/Energy ";

    for(int i=0; i<EnergiesVec.size(); i++){
        if(graphtype == 0){
            LatexText << EnergiesVec[i] << " ";
        }else if(graphtype == 1){
            LatexText << EnergiesVec[i] << " " << XLabels[i].toStdString() << " ";
        }
    }
    LatexText << "\n";

    for ( auto it = DataForROOTMap.begin(); it != DataForROOTMap.end(); ++it  ){
        QString GraphName = it.key();

        LatexText << GraphName.toStdString() << " ";

        for(int i=0; i<EnergiesVec.size(); i++){
            if(__isnan(DataForROOTMap[GraphName][EnergiesVec[i]]) || __isinf(DataForROOTMap[GraphName][EnergiesVec[i]])){
                DataForROOTMap[GraphName][EnergiesVec[i]] = MinValForLog;
            }
            LatexText << DataForROOTMap[GraphName][EnergiesVec[i]] << " ";
        }
        LatexText << "\n";
    }

    outfile << LatexText.str();
    outfile.close();

    if( QFile::exists(FileName) && QFile::exists(Root_Lib_dir_path) && QFile::exists(Root_Lib_dir_path+"/thisroot.sh") && QFile::exists(DoseCalcsCore_build_dir_path+"/"+GraphExecutableName) ){

        QString BashCommandsForExecuting =
                "#! /bin/bash \n . " +Root_Lib_dir_path +"/thisroot.sh \n"+
                "cd " + DoseCalcsCore_build_dir_path + "\n" +
                DoseCalcsCore_build_dir_path+"/"+GraphExecutableName + " a b " + UserCurrentResultsDirPath+"/"+GraphsOutDirName+"\n"
                "cd " + UserCurrentResultsDirPath+"/"+GraphsOutDirName + "\n" +
                "root --web=off -e \"TBrowser x\" ";
        //"root -e \"TBrowser b(\"Graph.root\");\" ";

        fileManagerObjectPlot->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

        showResultsOutput("Writing Analysis Commands : \n", 0);
        showResultsOutput(BashCommandsForExecuting , 0);
        showResultsOutput("to --> " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , 0);

        if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName)){
            terminal *term = new terminal;
            term->executeLocalCommand(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
        }
        else{
            showResultsOutput("Cannot find file containing execution commands "+ DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName + " , you should build DoseCalcs with ROOT Analysis option" , 3);
        }

    }else{


        QMessageBox::information(this, tr(""), "Canno't find ROOT file and/or analysis executable, check the files paths.");


        showResultsOutput("Cannot find the analysis executable \""+ DoseCalcsCore_build_dir_path+"/"+GraphExecutableName + " , you should build DoseCalcs with ROOT Analysis option" , 3);
    }


}

void PlotDialog::OpenDialogOfCombinationChooser()
{
    QDialog * d = new QDialog(); d->setWindowTitle("Specify the combinations");

    QGridLayout* GraphLayout = new QGridLayout;

    CombinationsText = new QTextEdit();
    CombinationsText->setPlaceholderText("Target <- Source");

    for(int dd=0; dd < SourceOrganForCombination.size();dd++){
        CombinationsText->append(TargetOrganForCombination[dd]+" <- "+SourceOrganForCombination[dd]);
    }
    TargetOrganForCombination.clear();
    SourceOrganForCombination.clear();

    TargetListForCombinations = new QComboBox(); TargetListForCombinations->addItems(TargetOrganlist);
    TargetListForCombinations->setToolTip("Target Region");
    TargetListForCombinations->setCurrentText("");

    QLabel* ToLabel = new QLabel();
    ToLabel->setText("<-");

    SourceListForCombinations = new QComboBox(); SourceListForCombinations->addItems(SourceOrganlist);
    SourceListForCombinations->setToolTip("Source Region");
    SourceListForCombinations->setCurrentText("");

    QPushButton* AddCombinationBtn = new QPushButton(); connect(AddCombinationBtn, SIGNAL(clicked()), this, SLOT(AddCombinationBtnSlot()));
    AddCombinationBtn->setText("Add");
    AddCombinationBtn->setToolTip("Add combination to the specific list");

    QPushButton* AddAllCombinationBtn = new QPushButton(); connect(AddAllCombinationBtn, SIGNAL(clicked()), this, SLOT(AddAllCombinationBtnSlot()));
    AddAllCombinationBtn->setText("Add All");
    AddAllCombinationBtn->setToolTip("Add all combinations");

    int ii = 0, jj=0;
    GraphLayout->addWidget(TargetListForCombinations, jj,ii,1,1);
    GraphLayout->addWidget(ToLabel, jj,++ii,1,1);
    GraphLayout->addWidget(SourceListForCombinations, jj,++ii,1,1);
    GraphLayout->addWidget(AddCombinationBtn, jj,++ii,1,1);
    GraphLayout->addWidget(AddAllCombinationBtn, jj,++ii,1,1);
    GraphLayout->addWidget(CombinationsText, 1,0,3,5);

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    GraphLayout->addWidget(buttonBox);

    d->setLayout(GraphLayout);

    int result = d->exec();

    if(result == QDialog::Accepted)
    {
        QStringList InputsVals;
        QVector< QPair<QString,QString>> Commlines ;

        Commlines = fileManagerObjectPlot->ReadTextFromQStringInQStringList(CombinationsText->toPlainText());

        for(int dd=0; dd < Commlines.size();dd++){
            Commlines[dd].second.replace("<-","");
            InputsVals = Commlines[dd].second.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);

            TargetOrganForCombination.push_back(Commlines[dd].first);
            if(InputsVals.size() == 1 ){
                SourceOrganForCombination.push_back(InputsVals[0]);
                QTextStream(stdout) << Commlines[dd].first << " !!<-!! " << InputsVals[0] << "\n";
            }
        }
    }
}
void PlotDialog::AddCombinationBtnSlot(){
    CombinationsText->append(TargetListForCombinations->currentText()+" <- "+SourceListForCombinations->currentText());
}
void PlotDialog::AddAllCombinationBtnSlot(){
    for(int cc=0; cc < SourceOrganlist.size();cc++){
        for(int dd=0; dd < SourceOrganlist.size();dd++){
            CombinationsText->append(TargetOrganlist[dd]+" <- "+SourceOrganlist[cc]);
        }
    }
}

void PlotDialog::closeEvent(QCloseEvent *event)  // show prompt when user wants to close app
{
    event->ignore();
    if (QMessageBox::Yes == QMessageBox::question(this, "Close Confirmation", "Exit?", QMessageBox::Yes | QMessageBox::No))
    {
        event->accept();
    }
}
void PlotDialog::keyPressEvent(QKeyEvent *e) {
    if(e->key() != Qt::Key_Escape) {QDialog::keyPressEvent(e);}
}
