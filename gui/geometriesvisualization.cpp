#include "geometriesvisualization.h"
#include "ui_geometriesvisualization.h"

extern QString DoseCalcsCore_build_dir_path;
extern QString MacroFilePath;
extern QString GUIPackagesAndFilesDirName;
extern QString GUIPackagesAndFilesDirPath;

extern QVector <QString> GeometryCommands ;
extern QVector <QString> VOXELCommands ;

//#include <QVector>
#include <iostream>
#include <fstream>
#include <QFuture>
#include <QtConcurrent/QtConcurrent>

geometriesvisualization::geometriesvisualization(QWidget *parent) : QDialog(parent), ui(new Ui::geometriesvisualization)
{
    ui->setupUi(this);

    // configure axis rect:
    ui->customPlot->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom); // this will also allow rescaling the color scale by dragging/zooming
    ui->customPlot->axisRect()->setupFullAxesBox(true);
    ui->customPlot->xAxis->setLabel("X(mm)");
    ui->customPlot->yAxis->setLabel("Y(mm)");

    /*
    // set up the QCPColorMap:
    QCPColorMap *colorMap = new QCPColorMap(ui->customPlot->xAxis, ui->customPlot->yAxis);
    int nx = 200;
    int ny = 200;
    colorMap->data()->setSize(nx, ny); // we want the color map to have nx * ny data points
    colorMap->data()->setRange(QCPRange(-4, 4), QCPRange(-4, 4)); // and span the coordinate range -4..4 in both key (x) and value (y) dimensions
    // now we assign some data, by accessing the QCPColorMapData instance of the color map:
    double x, y, z;
    for (int xIndex=0; xIndex<nx; ++xIndex)
    {
        for (int yIndex=0; yIndex<ny; ++yIndex)
        {
            colorMap->data()->cellToCoord(xIndex, yIndex, &x, &y);
            double r = 3*qSqrt(x*x+y*y)+1e-2;
            z = 2*x*(qCos(r+2)/r-qSin(r+2)/r); // the B field strength of dipole radiation (modulo physical constants)
            colorMap->data()->setCell(xIndex, yIndex, z);
        }
    }

    // add a color scale:
    QCPColorScale *colorScale = new QCPColorScale(ui->customPlot);
    ui->customPlot->plotLayout()->addElement(0, 1, colorScale); // add it to the right of the main axis rect
    colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
    colorMap->setColorScale(colorScale); // associate the color map with the color scale
    //colorScale->axis()->setLabel("Magnetic Field Strength");

    // set the color gradient of the color map to one of the presets:
    colorMap->setGradient(QCPColorGradient::gpPolar);
    // we could have also created a QCPColorGradient instance and added own colors to
    // the gradient, see the documentation of QCPColorGradient for what's possible.

    // rescale the data dimension (color) such that all data points lie in the span visualized by the color gradient:
    colorMap->rescaleDataRange();

    // make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
    QCPMarginGroup *marginGroup = new QCPMarginGroup(ui->customPlot);
    ui->customPlot->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
    colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);

    // rescale the key (x) and value (y) axes so the whole color map is visible:
    ui->customPlot->rescaleAxes();
*/

    ui->ChangeAxiss->setVisible(false);
    ui->GenerateMap->setVisible(false);

    ui->radioButtonVoxIDs->setChecked(true);

    QStringList Axis=(QStringList()<<"Z"<<"Y"<<"X");
    ui->Axis->addItems(Axis);

    QStringList Ext=(QStringList()<<"bmp"<<"jpg"<<"pdf"<<"png");
    ui->comboBoxSaveExt->addItems(Ext);

    QStringList ColorGrad=(QStringList()<<"gpGrayscale"<<"gpHot"<<"gpCold"<<"gpNight"<<"gpCandy"<<"gpGeography"
                           <<"gpIon"<<"gpThermal"<<"gpPolar"<<"gpSpectrum"<<"gpJet"<<"gpHues");
    ui->comboBoxColorGradien->addItems(ColorGrad);

    ui->progressBarFileReading->setRange(0, 100);
    ui->progressBarFileReading->setValue(0);
    ui->progressBarFileReading->show();

    QTextStream(stdout) << MacroFilePath << "\n";

    QMap <QString,QString> lines = fileManagerObject->ReadLinesFromFileWithFirstWordIndicator(MacroFilePath);

    QTextStream(stdout) << GeometryCommands[2] << "\n";
    QTextStream(stdout) << VOXELCommands[0] << "\n";

    QStringList InputsVals = lines[GeometryCommands[2]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);

    if(InputsVals.size() > 0){
        if( InputsVals[0] == "VoxIDs" ){
            if(InputsVals.size()>1){

                QString s1 = "../"+GUIPackagesAndFilesDirName;
                if(InputsVals[1].contains(s1)){InputsVals[1] = InputsVals[1].replace(s1,GUIPackagesAndFilesDirPath);}
                ui->VoxelsDataFile->setText(InputsVals[1]);

                InputsVals = lines[VOXELCommands[0]].split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts);
                if(InputsVals.size() > 7){

                    ui->spinBoxNumberX->setValue(QString(InputsVals[0]).toInt());
                    ui->spinBoxNumberY->setValue(QString(InputsVals[1]).toInt());
                    ui->spinBoxNumberZ->setValue(QString(InputsVals[2]).toInt());
                    ui->doubleSpinBoxSizeX->setValue(QString(InputsVals[5]).toDouble()*2);
                    ui->doubleSpinBoxSizeY->setValue(QString(InputsVals[6]).toDouble()*2);
                    ui->doubleSpinBoxSizeZ->setValue(QString(InputsVals[7]).toDouble()*2);

                }
            }
        }else{

        }
    }
}

geometriesvisualization::~geometriesvisualization()
{
    delete ui;
}

void geometriesvisualization::getInputData()
{
    VoxelsDataFilePath = ui->VoxelsDataFile->text();
    Axis = ui->Axis->currentText();
    SliceID = ui->horizontalSliderSliceID->value();
    VoxXNumber = ui->spinBoxNumberX->value();
    VoxYNumber = ui->spinBoxNumberY->value();
    VoxZNumber = ui->spinBoxNumberZ->value();
    VoxXSize = ui->doubleSpinBoxSizeX->value();
    VoxYSize = ui->doubleSpinBoxSizeY->value();
    VoxZSize = ui->doubleSpinBoxSizeZ->value();
}

void geometriesvisualization::on_GenerateMap_clicked()
{

    getInputData();

    int ni;
    int nj;
    int nk;
    double smini; double smaxi;
    double sminj; double smaxj;
    //QMap <int , QMap <int , QMap <int , double>>> mapdata;
    double *** mapdata;

    ni = 200;
    nj = 200;
    nk = 200;
    smini = -4; smaxi = 4;
    sminj = -4; smaxj = 4;
    if(Axis == "Z"){
        ni = VoxYNumber;
        nj = VoxXNumber;
        nk = VoxZNumber;
        smini = -(VoxYNumber*VoxYSize)/2; smaxi = (VoxYNumber*VoxYSize)/2;
        sminj = -(VoxXNumber*VoxXSize)/2; smaxj = (VoxXNumber*VoxXSize)/2;
        mapdata = Z_YX ;
        ui->customPlot->xAxis->setLabel("Y(mm)");
        ui->customPlot->yAxis->setLabel("X(mm)");
        ui->horizontalSliderSliceID->setMaximum(VoxZNumber);
    }
    else if(Axis == "Y"){
        ni = VoxZNumber;
        nj = VoxXNumber;
        nk = VoxYNumber;
        smini = -(VoxZNumber*VoxZSize)/2; smaxi = (VoxZNumber*VoxZSize)/2;
        sminj = -(VoxXNumber*VoxXSize)/2; smaxj = (VoxXNumber*VoxXSize)/2;
        mapdata = Y_ZX ;
        ui->customPlot->xAxis->setLabel("Z(mm)");
        ui->customPlot->yAxis->setLabel("X(mm)");
        ui->horizontalSliderSliceID->setMaximum(VoxYNumber);
    }
    else if(Axis == "X"){
        ni = VoxZNumber;
        nj = VoxYNumber;
        nk = VoxXNumber;
        smini = -(VoxZNumber*VoxZSize)/2; smaxi = (VoxZNumber*VoxZSize)/2;
        sminj = -(VoxYNumber*VoxYSize)/2; smaxj = (VoxYNumber*VoxYSize)/2;
        mapdata = X_ZY ;
        ui->customPlot->xAxis->setLabel("Z(mm)");
        ui->customPlot->yAxis->setLabel("Y(mm)");
        ui->horizontalSliderSliceID->setMaximum(VoxXNumber);
    }

    // set up the QCPColorMap:
    QCPColorMap *colorMap = new QCPColorMap(ui->customPlot->xAxis, ui->customPlot->yAxis);
    colorMap->data()->setSize(ni, nj); // we want the color map to have ni * nj data points
    colorMap->data()->setRange(QCPRange(smini, smaxi), QCPRange(sminj, smaxj)); // and span the coordinate range -4..4 in both key (x) and value (y) dimensions
    // now we assign some data, by accessing the QCPColorMapData instance of the color map:
    double i, j, k;
    for (int iIndex=0; iIndex<ni; ++iIndex)
    {
        for (int jIndex=0; jIndex<nj; ++jIndex)
        {
            colorMap->data()->cellToCoord(iIndex, jIndex, &i, &j);
            k = mapdata[SliceID][iIndex][jIndex]; // the B field strength of dipole radiation (modulo physical constants)
            colorMap->data()->setCell(iIndex, jIndex, k);

            //QTextStream(stdout) << iIndex << " " << jIndex << " " << SliceID << " " << k << "\n";
        }
    }

    // add a color scale:
    //QCPColorScale *colorScale = new QCPColorScale(ui->customPlot);
    //ui->customPlot->plotLayout()->addElement(0, 1, colorScale); // add it to the right of the main axis rect
    //colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
    //colorMap->setColorScale(colorScale); // associate the color map with the color scale
    //colorScale->axis()->setLabel("Material color");

    // set the color gradient of the color map to one of the presets:

    if(ui->comboBoxColorGradien->currentText()=="gpGrayscale"){colorMap->setGradient(QCPColorGradient::gpGrayscale);}
    else if(ui->comboBoxColorGradien->currentText()=="gpHot"){colorMap->setGradient(QCPColorGradient::gpHot);}
    else if(ui->comboBoxColorGradien->currentText()=="gpCold"){colorMap->setGradient(QCPColorGradient::gpCold);}
    else if(ui->comboBoxColorGradien->currentText()=="gpNight"){colorMap->setGradient(QCPColorGradient::gpNight);}
    else if(ui->comboBoxColorGradien->currentText()=="gpGeography"){colorMap->setGradient(QCPColorGradient::gpGeography);}
    else if(ui->comboBoxColorGradien->currentText()=="gpIon"){colorMap->setGradient(QCPColorGradient::gpIon);}
    else if(ui->comboBoxColorGradien->currentText()=="gpThermal"){colorMap->setGradient(QCPColorGradient::gpThermal);}
    else if(ui->comboBoxColorGradien->currentText()=="gpPolar"){colorMap->setGradient(QCPColorGradient::gpPolar);}
    else if(ui->comboBoxColorGradien->currentText()=="gpSpectrum"){colorMap->setGradient(QCPColorGradient::gpSpectrum);}
    if(ui->comboBoxColorGradien->currentText()=="gpJet"){colorMap->setGradient(QCPColorGradient::gpJet);}
    if(ui->comboBoxColorGradien->currentText()=="gpHues"){colorMap->setGradient(QCPColorGradient::gpHues);}

    // we could have also created a QCPColorGradient instance and added own colors to
    // the gradient, see the documentation of QCPColorGradient for what's possible.

    // rescale the data dimension (color) such that all data points lie in the span visualized by the color gradient:
    colorMap->rescaleDataRange();

    // make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
    QCPMarginGroup *marginGroup = new QCPMarginGroup(ui->customPlot);
    ui->customPlot->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
    //colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);

    // rescale the key (x) and value (y) axes so the whole color map is visible:
    ui->customPlot->rescaleAxes();
    ui->customPlot->xAxis->rescale();
    ui->customPlot->yAxis->rescale();
    ui->customPlot->replot();
}
void geometriesvisualization::on_Save_clicked()
{
    bool ok;
    QString newLabel = QInputDialog::getText(this, "Enter the File Name", "New File Name:", QLineEdit::Normal, "Phantom", &ok);
    if (ok)
    {
        FileExtToSave = ui->comboBoxSaveExt->currentText();
        if(FileExtToSave == "bmp"){
            ui->customPlot->saveBmp(DoseCalcsCore_build_dir_path +"/"+ newLabel);
        }else if(FileExtToSave == "jpg"){
            ui->customPlot->saveJpg(DoseCalcsCore_build_dir_path +"/"+ newLabel);
        }else if(FileExtToSave == "bdf"){
            ui->customPlot->savePdf(DoseCalcsCore_build_dir_path +"/"+ newLabel);
        }else if(FileExtToSave == "png"){
            ui->customPlot->savePng(DoseCalcsCore_build_dir_path +"/"+ newLabel);
        }
    }


}
void geometriesvisualization::on_AddVoxelsDataFile_clicked()
{
    VoxelsDataFilePath = QFileDialog::getOpenFileName(
                this,
                tr("Choose voxels data file"),
                DoseCalcsCore_build_dir_path, //pathBuildApp,
                "All files (*.*);;Text files (*.txt)"
                );
    ui->VoxelsDataFile->setText(VoxelsDataFilePath);

}


void geometriesvisualization::on_pushButtonReadVoxelsData_clicked(){
    //QFuture<void> future = QtConcurrent::run(this, &geometriesvisualization::ReadVoxelsFilesAndFillMaps);
    ReadVoxelsFilesAndFillMaps();
}
void geometriesvisualization::ReadVoxelsFilesAndFillMaps()
{
    getInputData();

    //X_ZY.clear();
    //Y_ZX.clear();
    //Z_YX.clear();

    std::ifstream fileR;

    fileR.open(VoxelsDataFilePath.toStdString(), std::ios::binary | std::ios::in);
    fileR.seekg(0, std::ios::end);
    long file_size = fileR.tellg();
    fileR.close();

    fileR.open(VoxelsDataFilePath.toStdString(), std::ios::binary | std::ios::in);
    if(fileR.is_open()){

        if(ui->radioButtonDicom->isChecked()){

            double fMinX , fMaxX , fMinY , fMaxY , fMinZ , fMaxZ;

            int nMaterials;
            fileR >> nMaterials;
            std::string mateName;
            int nmate;
            for( int ii = 0; ii < nMaterials; ii++ ){
                fileR >> nmate;
                fileR >> mateName;
            }

            fileR >> VoxXNumber >> VoxYNumber >> VoxZNumber >> fMinX >> fMaxX >> fMinY >> fMaxY >> fMinZ >> fMaxZ;
            VoxXSize = (fMaxX-fMinX)/VoxXNumber;
            VoxYSize = (fMaxY-fMinY)/VoxYNumber;
            VoxZSize = (fMaxZ-fMinZ)/VoxZNumber;

            ui->spinBoxNumberX->setValue(VoxXNumber);
            ui->spinBoxNumberY->setValue(VoxYNumber);
            ui->spinBoxNumberZ->setValue(VoxXNumber);
            ui->doubleSpinBoxSizeX->setValue(VoxXSize);
            ui->doubleSpinBoxSizeY->setValue(VoxYSize);
            ui->doubleSpinBoxSizeZ->setValue(VoxZSize);

        }
        else if(ui->radioButtonPET->isChecked()){

            double fMinX , fMaxX , fMinY , fMaxY , fMinZ , fMaxZ;

            fileR >> VoxXNumber >> VoxYNumber >> VoxZNumber >> fMinX >> fMaxX >> fMinY >> fMaxY >> fMinZ >> fMaxZ;
            VoxXSize = (fMaxX-fMinX)/VoxXNumber;
            VoxYSize = (fMaxY-fMinY)/VoxYNumber;
            VoxZSize = (fMaxZ-fMinZ)/VoxZNumber;

            ui->spinBoxNumberX->setValue(VoxXNumber);
            ui->spinBoxNumberY->setValue(VoxYNumber);
            ui->spinBoxNumberZ->setValue(VoxXNumber);
            ui->doubleSpinBoxSizeX->setValue(VoxXSize);
            ui->doubleSpinBoxSizeY->setValue(VoxYSize);
            ui->doubleSpinBoxSizeZ->setValue(VoxZSize);
        }
        else if (ui->radioButtonVoxIDs->isChecked()){

        }


        long curr, end;
        int percent;
        end = file_size;

        ui->progressBarFileReading->setValue(0);

        std::cout << "\n\n\n\n--------- Reading file started --------- " << std::endl;
        double ID;

        // initialize the arrays

        Z_YX = new double**[VoxZNumber];
        for(int f = 0; f < VoxZNumber ;f++ ){
            double**AA = new double*[VoxYNumber];
            for(int g = 0; g < VoxYNumber ;g++ ){
                double*A = new double[VoxXNumber];
                for(int d = 0; d < VoxXNumber ;d++ ){
                    A[d] = 0.;
                }
                AA[g]=A;
            }
            Z_YX[f] = AA;
        }

        Y_ZX = new double**[VoxYNumber];
        for(int f = 0; f < VoxYNumber ;f++ ){
            double**AA = new double*[VoxZNumber];
            for(int g = 0; g < VoxZNumber ;g++ ){
                double*A = new double[VoxXNumber];
                for(int d = 0; d < VoxXNumber ;d++ ){
                    A[d] = 0.;
                }
                AA[g]=A;
            }
            Y_ZX[f] = AA;
        }

        X_ZY = new double**[VoxXNumber];
        for(int f = 0; f < VoxXNumber ;f++ ){
            double**AA = new double*[VoxZNumber];
            for(int g = 0; g < VoxZNumber ;g++ ){
                double*A = new double[VoxYNumber];
                for(int d = 0; d < VoxYNumber ;d++ ){
                    A[d] = 0.;
                }
                AA[g]=A;
            }
            X_ZY[f] = AA;
        }


        //std::cout << VoxZNumber*VoxYNumber*VoxXNumber << " " << fileR.tellg() << std::endl;

        //unsigned int voxnum = VoxZNumber*VoxYNumber*VoxXNumber;
        //unsigned int voxinc = 0;
        //colormaxvalue = 0;
        for(int f = 0; f < VoxZNumber ;f++ ){

            curr = fileR.tellg();
            if (curr != -1)
                percent = curr * 100 / end;
            else
                percent = 100;
            ui->progressBarFileReading->setValue(percent);

            for(int g = 0; g < VoxYNumber ;g++ ){
                for(int d = 0; d < VoxXNumber ;d++ ){

                    fileR >> ID ;

                    //if(colormaxvalue < ID){colormaxvalue = ID;}

                    //if(voxnum > voxinc){
                      //  QMessageBox msgBox;
                        //msgBox.setWindowTitle("Error in reading IDs file");
                        //msgBox.setText("The number of readed voxels IDs ("+QString::number(voxinc)+") and exceds the declared number ("+QString::number(voxnum)+")");
                        //msgBox.exec();
                        //return;
                    //}else{
                        Z_YX[f][g][d] = ID;
                        Y_ZX[g][f][d] = ID;
                        X_ZY[d][f][g] = ID;
                    //}

                    //voxinc++;

                    //std::cout << f << " " << g << " " << d << " ID: " << ID << std::endl;
                }
            }
        }

        ui->progressBarFileReading->setValue(100);

        getInputData();
        on_Axis_activated(0);

        std::cout << "--------- Reading file terminated --------- \n\n\n\n" << std::endl;
    }


}



void geometriesvisualization::on_radioButtonVoxIDs_clicked(bool checked)
{
    if(checked){
        ui->spinBoxNumberX->setEnabled(true);
        ui->spinBoxNumberY->setEnabled(true);
        ui->spinBoxNumberZ->setEnabled(true);
        ui->doubleSpinBoxSizeX->setEnabled(true);
        ui->doubleSpinBoxSizeY->setEnabled(true);
        ui->doubleSpinBoxSizeZ->setEnabled(true);
    }
    else{
        ui->spinBoxNumberX->setEnabled(false);
        ui->spinBoxNumberY->setEnabled(false);
        ui->spinBoxNumberZ->setEnabled(false);
        ui->doubleSpinBoxSizeX->setEnabled(false);
        ui->doubleSpinBoxSizeY->setEnabled(false);
        ui->doubleSpinBoxSizeZ->setEnabled(false);
    }
}
void geometriesvisualization::on_radioButtonDicom_toggled(bool checked)
{
    if(checked){
        ui->spinBoxNumberX->setEnabled(false);
        ui->spinBoxNumberY->setEnabled(false);
        ui->spinBoxNumberZ->setEnabled(false);
        ui->doubleSpinBoxSizeX->setEnabled(false);
        ui->doubleSpinBoxSizeY->setEnabled(false);
        ui->doubleSpinBoxSizeZ->setEnabled(false);
    }
    else{
        ui->spinBoxNumberX->setEnabled(true);
        ui->spinBoxNumberY->setEnabled(true);
        ui->spinBoxNumberZ->setEnabled(true);
        ui->doubleSpinBoxSizeX->setEnabled(true);
        ui->doubleSpinBoxSizeY->setEnabled(true);
        ui->doubleSpinBoxSizeZ->setEnabled(true);
    }

}
void geometriesvisualization::on_radioButtonPET_clicked(bool checked)
{
    if(checked){
        ui->spinBoxNumberX->setEnabled(false);
        ui->spinBoxNumberY->setEnabled(false);
        ui->spinBoxNumberZ->setEnabled(false);
        ui->doubleSpinBoxSizeX->setEnabled(false);
        ui->doubleSpinBoxSizeY->setEnabled(false);
        ui->doubleSpinBoxSizeZ->setEnabled(false);
    }
    else{
        ui->spinBoxNumberX->setEnabled(true);
        ui->spinBoxNumberY->setEnabled(true);
        ui->spinBoxNumberZ->setEnabled(true);
        ui->doubleSpinBoxSizeX->setEnabled(true);
        ui->doubleSpinBoxSizeY->setEnabled(true);
        ui->doubleSpinBoxSizeZ->setEnabled(true);
    }
}

void geometriesvisualization::on_pushButtonOpenResDir_clicked()
{
    if(QFile::exists(DoseCalcsCore_build_dir_path)){
        QString command = DoseCalcsCore_build_dir_path;
        QProcess process;
        QStringList qsl = {command};
        process.startDetached("nautilus", qsl);
    }else{
        std::cout << "cannot find directory of the saved document" << std::endl;
    }
}



void geometriesvisualization::on_horizontalSliderSliceID_sliderMoved(int position)
{
    //ui->horizontalSliderSliceID->setToolTip(QString::number(position));
    ui->labelVal->setText(QString::number(position));
    //ui->labelVal->repaint();

    on_GenerateMap_clicked();
    ui->customPlot->replot();
    //std::cout << position << std::endl;
}


void geometriesvisualization::showTETrahedrons()
{



}


void geometriesvisualization::on_Axis_activated(int index)
{
    if(index == 0){
        ui->horizontalSliderSliceID->setMaximum(VoxZNumber);
        ui->labelMax->setText(QString::number(VoxZNumber));
    }
    else if(index == 1){
        ui->horizontalSliderSliceID->setMaximum(VoxYNumber);
        ui->labelMax->setText(QString::number(VoxYNumber));
    }
    else if(index == 2){
        ui->horizontalSliderSliceID->setMaximum(VoxXNumber);
        ui->labelMax->setText(QString::number(VoxXNumber));
    }
    ui->horizontalSliderSliceID->setValue(0);
    //ui->horizontalSliderSliceID->repaint();
    on_horizontalSliderSliceID_sliderMoved(0);
    ui->customPlot->replot();
}
void geometriesvisualization::on_comboBoxColorGradien_currentTextChanged(const QString &arg1)
{
    on_GenerateMap_clicked();
}
void geometriesvisualization::on_horizontalSlider_valueChanged(int value)
{
    on_GenerateMap_clicked();
}

