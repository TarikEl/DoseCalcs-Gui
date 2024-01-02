#include "geometrymodellingeditor.h"
#include "ui_geometrymodellingeditor.h"

#include "QFileDialog"
#include "QFile"
#include "QMessageBox"
#include "QTextStream"
#include "QFont"
#include "QFontDialog"
#include "QFile"
//#include <QtXml/>

GeometryModellingEditor::GeometryModellingEditor(QWidget *parent) : QDialog(parent), ui(new Ui::GeometryModellingEditor)
{
    ui->setupUi(this);

    QStringList GeomMeth=(QStringList()<<"GDML"<<"TEXT"<<"C++"<<"Command");
    ui->comboBoxGeometrySyntax->addItems(GeomMeth);
    highlighter = new Highlighter(ui->textEdit->document());

}

GeometryModellingEditor::~GeometryModellingEditor()
{
    delete ui;
}


void GeometryModellingEditor::on_comboBoxGeometrySyntax_currentTextChanged(const QString &arg1)
{
    if(arg1 == "GDML"){

    }
    else if(arg1 == "TEXT"){

    }
    else if(arg1 == "C++"){

    }
    else if(arg1 == "Command"){

    }
    showGDMLHelperFrame();
}

void GeometryModellingEditor::showGDMLHelperFrame(){

    framLay = new QGridLayout;

    btnSaveIsotope= new QPushButton(); btnSaveIsotope->setText("Add Isotope"); connect(btnSaveIsotope, SIGNAL(clicked()), this, SLOT(btnSaveIsotope_slot()));
    btnSaveElem= new QPushButton(); btnSaveElem->setText("Add Element"); connect(btnSaveElem, SIGNAL(clicked()), this, SLOT(btnSaveElem_slot()));
    btnSaveMaterial= new QPushButton(); btnSaveMaterial->setText("Add Material"); connect(btnSaveMaterial, SIGNAL(clicked()), this, SLOT(btnSaveMaterial_slot()));
    btnAddElementForMaterial= new QPushButton(); btnAddElementForMaterial->setText("Add Element"); connect(btnAddElementForMaterial, SIGNAL(clicked()), this, SLOT(btnAddElementForMaterial_slot()));

    btnSaveConstants= new QPushButton(); btnSaveConstants->setText("Add Constant"); connect(btnSaveConstants, SIGNAL(clicked()), this, SLOT(btnSaveConstants_slot()));
    btnSaveQuantities= new QPushButton(); btnSaveQuantities->setText("Add Quantity"); connect(btnSaveQuantities, SIGNAL(clicked()), this, SLOT(btnSaveQuantities_slot()));
    btnSaveVariables= new QPushButton(); btnSaveVariables->setText("Add Variable"); connect(btnSaveVariables, SIGNAL(clicked()), this, SLOT(btnSaveVariables_slot()));
    btnSaveScales= new QPushButton(); btnSaveScales->setText("Add scales"); connect(btnSaveScales, SIGNAL(clicked()), this, SLOT(btnSaveScales_slot()));
    btnSaveMatrices= new QPushButton(); btnSaveMatrices->setText("Add Matrice"); connect(btnSaveMatrices, SIGNAL(clicked()), this, SLOT(btnSaveMatrices_slot()));
    btnSavePosition= new QPushButton(); btnSavePosition->setText("Add Position"); connect(btnSavePosition, SIGNAL(clicked()), this, SLOT(btnSavePosition_slot()));
    btnSaveRotation= new QPushButton(); btnSaveRotation->setText("Add rotation"); connect(btnSaveRotation, SIGNAL(clicked()), this, SLOT(btnSaveRotation_slot()));

    btnSaveSolid= new QPushButton(); btnSaveSolid->setText("Add Solid"); connect(btnSaveSolid, SIGNAL(clicked()), this, SLOT(btnSaveSolid_slot()));
    btnSaveLogicalVolume= new QPushButton(); btnSaveLogicalVolume->setText("Add Log Vol"); connect(btnSaveLogicalVolume, SIGNAL(clicked()), this, SLOT(btnSaveLogicalVolume_slot()));
    btnSavePhysicalVolume= new QPushButton(); btnSavePhysicalVolume->setText("Add Phy Vol"); connect(btnSavePhysicalVolume, SIGNAL(clicked()), this, SLOT(btnSavePhysicalVolume_slot()));
    btnSaveTopVolume= new QPushButton(); btnSaveTopVolume->setText("Add Volume"); connect(btnSaveTopVolume, SIGNAL(clicked()), this, SLOT(btnSaveTopVolume_slot()));

    IsotopeData = new QLineEdit(); IsotopeData->setPlaceholderText(""); IsotopeData->setToolTip("");
    ElementData = new QLineEdit(); ElementData->setPlaceholderText(""); ElementData->setToolTip("");
    MaterialData = new QLineEdit(); MaterialData->setPlaceholderText(""); MaterialData->setToolTip("");
    ElementMaterialData = new QLineEdit(); ElementMaterialData->setPlaceholderText(""); ElementMaterialData->setToolTip("");

    ConstantsData = new QLineEdit(); ConstantsData->setPlaceholderText(""); ConstantsData->setToolTip("");
    QuantitiesData = new QLineEdit(); QuantitiesData->setPlaceholderText(""); QuantitiesData->setToolTip("");
    VariablesData = new QLineEdit(); VariablesData->setPlaceholderText(""); VariablesData->setToolTip("");
    ScalesData = new QLineEdit(); ScalesData->setPlaceholderText(""); ScalesData->setToolTip("");
    MatricesData = new QLineEdit(); MatricesData->setPlaceholderText(""); MatricesData->setToolTip("");
    PositionData = new QLineEdit(); PositionData->setPlaceholderText(""); PositionData->setToolTip("");
    RotationData = new QLineEdit(); RotationData->setPlaceholderText(""); RotationData->setToolTip("");

    SolidSpecications = new QLineEdit(); SolidSpecications->setPlaceholderText(""); SolidSpecications->setToolTip("");
    LogicalVolumeData = new QLineEdit(); LogicalVolumeData->setPlaceholderText(""); LogicalVolumeData->setToolTip("");
    PhysicalVolumeData = new QLineEdit(); PhysicalVolumeData->setPlaceholderText(""); PhysicalVolumeData->setToolTip("");
    TopVolumeData = new QLineEdit(); TopVolumeData->setPlaceholderText(""); TopVolumeData->setToolTip("");

    ComboxSolidType = new QComboBox(); QStringList SolType=(QStringList()<<"Box"<<"Tubs"<<"CutTubs"<<"Cons"<<"Para"<<"Trd"<<"Sphere"<<"Orb"<<"Torus"<<"Ellipsoid"<<"Union"<<"Intersection"<<"Subtraction"); ComboxSolidType->addItems(SolType);
    connect(ComboxSolidType, SIGNAL(currentTextChanged(QString)), this, SLOT(ComboxSolidType_slot()));

    ComboxElementsNames = new QComboBox();
    ComboxSolidNames = new QComboBox();
    ComboxMaterialsNames = new QComboBox();
    ComboxPositionsNames = new QComboBox();
    ComboxRotationsNames = new QComboBox();
    ComboxLogicalVolumesNames = new QComboBox();

    int ii = 0, jj = 0;
    framLay->addWidget(IsotopeData, jj,ii,1,1);
    framLay->addWidget(btnSaveIsotope, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(ElementData, jj,ii,1,1);
    framLay->addWidget(btnSaveElem, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(MaterialData, jj,ii,1,1);
    framLay->addWidget(btnSaveMaterial, jj,++ii,1,1);
    framLay->addWidget(ComboxElementsNames, jj,++ii,1,1);
    framLay->addWidget(ElementMaterialData, jj,++ii,1,1);
    framLay->addWidget(btnAddElementForMaterial, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(ConstantsData, jj,ii,1,1);
    framLay->addWidget(btnSaveConstants, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(QuantitiesData, jj,ii,1,1);
    framLay->addWidget(btnSaveQuantities, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(VariablesData, jj,ii,1,1);
    framLay->addWidget(btnSaveVariables, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(ScalesData, jj,ii,1,1);
    framLay->addWidget(btnSaveScales, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(MatricesData, jj,ii,1,1);
    framLay->addWidget(btnSaveMatrices, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(PositionData, jj,ii,1,1);
    framLay->addWidget(btnSavePosition, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(RotationData, jj,ii,1,1);
    framLay->addWidget(btnSaveRotation, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(ComboxSolidType, jj,ii,1,1);
    framLay->addWidget(SolidSpecications, jj,++ii,1,1);
    framLay->addWidget(btnSaveSolid, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(LogicalVolumeData, jj,ii,1,1);
    framLay->addWidget(ComboxSolidNames, jj,++ii,1,1);
    framLay->addWidget(ComboxMaterialsNames, jj,++ii,1,1);
    framLay->addWidget(btnSaveLogicalVolume, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(PhysicalVolumeData, jj,ii,1,1);
    framLay->addWidget(ComboxLogicalVolumesNames, jj,++ii,1,1);
    framLay->addWidget(ComboxPositionsNames, jj,++ii,1,1);
    framLay->addWidget(ComboxRotationsNames, jj,++ii,1,1);
    framLay->addWidget(btnSavePhysicalVolume, jj,++ii,1,1);
    ii = 0, jj++;
    framLay->addWidget(TopVolumeData, jj,ii,1,1);
    framLay->addWidget(btnSaveTopVolume, jj,++ii,1,1);

    ui->frameSyntaxGenerator->setLayout(framLay);
    ui->frameSyntaxGenerator->setVisible(true);
}

void GeometryModellingEditor::ComboxSolidType_slot(){

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
void GeometryModellingEditor::btnSaveIsotope_slot(){}
void GeometryModellingEditor::btnSaveElem_slot(){


}
void GeometryModellingEditor::btnSaveMaterial_slot(){}
void GeometryModellingEditor::btnAddElementForMaterial_slot(){}
void GeometryModellingEditor::btnSaveConstants_slot(){}
void GeometryModellingEditor::btnSaveQuantities_slot(){}
void GeometryModellingEditor::btnSaveVariables_slot(){}
void GeometryModellingEditor::btnSaveScales_slot(){}
void GeometryModellingEditor::btnSaveMatrices_slot(){}
void GeometryModellingEditor::btnSavePosition_slot(){}
void GeometryModellingEditor::btnSaveRotation_slot(){}
void GeometryModellingEditor::btnSaveSolid_slot(){}
void GeometryModellingEditor::btnSaveLogicalVolume_slot(){}
void GeometryModellingEditor::btnSavePhysicalVolume_slot(){}
void GeometryModellingEditor::btnSaveTopVolume_slot(){}


void GeometryModellingEditor::on_NewBtn_clicked()
{

}
void GeometryModellingEditor::on_OpenBtn_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, "Open the file");
    QFile file(fileName);
    currentFile = fileName;
    if (!file.open(QIODevice::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning", "Cannot open file: " + file.errorString());
        return;
    }
    setWindowTitle(fileName);
    QTextStream in(&file);
    QString text = in.readAll();
    ui->textEdit->setText(text);
    file.close();

}
void GeometryModellingEditor::on_SaveBtn_clicked()
{
    QString fileName;
    // If we don't have a filename from before, get one.
    if (currentFile.isEmpty()) {
        fileName = QFileDialog::getSaveFileName(this, "Save");
        currentFile = fileName;
    } else {
        fileName = currentFile;
    }
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning", "Cannot save file: " + file.errorString());
        return;
    }
    setWindowTitle(fileName);
    QTextStream out(&file);
    QString text = ui->textEdit->toPlainText();
    out << text;
    file.close();
}
void GeometryModellingEditor::on_SaveAsBtn_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, "Save as");
    QFile file(fileName);

    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning", "Cannot save file: " + file.errorString());
        return;
    }
    currentFile = fileName;
    setWindowTitle(fileName);
    QTextStream out(&file);
    QString text = ui->textEdit->toPlainText();
    out << text;
    file.close();
}
void GeometryModellingEditor::on_PrintBtn_clicked()
{
    QString text = "<FONT  class   = \" result\"                    > tarik elghalbzouri name=\"tarik\"  <       /FONT >";
    QRegularExpression regex("<FONT class=\"result\">(.*)</FONT>");
    QRegularExpressionMatch match = regex.match(text);
    QString textYouWant = match.captured(1);
    ui->textEdit->setText(textYouWant);

/*
#if defined(QT_PRINTSUPPORT_LIB) && QT_CONFIG(printer)
    QPrinter printDev;
#if QT_CONFIG(printdialog)
    QPrintDialog dialog(&printDev, this);
    if (dialog.exec() == QDialog::Rejected)
        return;
#endif // QT_CONFIG(printdialog)
    ui->textEdit->print(&printDev);
#endif // QT_CONFIG(printer)
*/
}
void GeometryModellingEditor::on_CopyBtn_clicked()
{
#if QT_CONFIG(clipboard)
    ui->textEdit->copy();
#endif
}
void GeometryModellingEditor::on_CutBtn_clicked()
{
#if QT_CONFIG(clipboard)
    ui->textEdit->cut();
#endif
}
void GeometryModellingEditor::on_PasteBtn_clicked()
{
#if QT_CONFIG(clipboard)
    ui->textEdit->paste();
#endif
}
void GeometryModellingEditor::on_UndoBtn_clicked()
{
    ui->textEdit->undo();
}
void GeometryModellingEditor::on_RedoBtn_clicked()
{
    ui->textEdit->redo();
}
void GeometryModellingEditor::on_FontBtn_clicked()
{
    bool fontSelected;
    QFont font = QFontDialog::getFont(&fontSelected, this);
    if (fontSelected)
        ui->textEdit->setFont(font);
}
void GeometryModellingEditor::on_BoldBtn_clicked()
{

}
void GeometryModellingEditor::on_ItalicBtn_clicked()
{

}
void GeometryModellingEditor::on_UnderlineBtn_clicked()
{

}
void GeometryModellingEditor::on_AboutBtn_clicked()
{

}
void GeometryModellingEditor::on_ExitBtn_clicked()
{
    //QCoreApplication::quit();
}

