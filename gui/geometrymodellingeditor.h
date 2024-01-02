#ifndef GEOMETRYMODELLINGEDITOR_H
#define GEOMETRYMODELLINGEDITOR_H

#include <QDialog>
#include "gui/highlighter.h"
#include "QGridLayout"
#include "QLineEdit"
#include "QComboBox"

namespace Ui {
class GeometryModellingEditor;
}

class GeometryModellingEditor : public QDialog
{
    Q_OBJECT

public:
    explicit GeometryModellingEditor(QWidget *parent = nullptr);
    ~GeometryModellingEditor();

private slots:

    void on_NewBtn_clicked();
    void on_OpenBtn_clicked();
    void on_SaveBtn_clicked();
    void on_SaveAsBtn_clicked();
    void on_PrintBtn_clicked();
    void on_CopyBtn_clicked();
    void on_CutBtn_clicked();
    void on_PasteBtn_clicked();
    void on_UndoBtn_clicked();
    void on_RedoBtn_clicked();
    void on_FontBtn_clicked();
    void on_BoldBtn_clicked();
    void on_ItalicBtn_clicked();
    void on_UnderlineBtn_clicked();
    void on_AboutBtn_clicked();
    void on_ExitBtn_clicked();

    // slots for dynamique widgets

    void btnSaveIsotope_slot();
    void btnSaveElem_slot();
    void btnSaveMaterial_slot();

    void btnSaveConstants_slot();
    void btnSaveQuantities_slot();
    void btnSaveVariables_slot();
    void btnSaveScales_slot();
    void btnSaveMatrices_slot();
    void btnSavePosition_slot();
    void btnSaveRotation_slot();

    void ComboxSolidType_slot();
    void btnSaveSolid_slot();
    void btnSaveLogicalVolume_slot();
    void btnSavePhysicalVolume_slot();
    void btnAddElementForMaterial_slot();

    void btnSaveTopVolume_slot();

    void on_comboBoxGeometrySyntax_currentTextChanged(const QString &arg1);

private:

    void showGDMLHelperFrame();

    // for GDML and TEXT

    QGridLayout *framLay;

    QPushButton *btnSaveIsotope;
    QPushButton *btnSaveElem;
    QPushButton *btnSaveMaterial;
    QPushButton *btnAddElementForMaterial;

    QPushButton *btnSaveConstants;
    QPushButton *btnSaveQuantities;
    QPushButton *btnSaveVariables;
    QPushButton *btnSaveScales;
    QPushButton *btnSaveMatrices;
    QPushButton *btnSavePosition;
    QPushButton *btnSaveRotation;

    QPushButton *btnSaveSolid;
    QPushButton *btnSaveLogicalVolume;
    QPushButton *btnSavePhysicalVolume;
    QPushButton *btnSaveTopVolume;

    QLineEdit *IsotopeData;
    QLineEdit *ElementData;
    QLineEdit *ElementMaterialData;
    QLineEdit *MaterialData;

    QLineEdit *ConstantsData;
    QLineEdit *QuantitiesData;
    QLineEdit *VariablesData;
    QLineEdit *ScalesData;
    QLineEdit *MatricesData;
    QLineEdit *PositionData;
    QLineEdit *RotationData;

    QComboBox *ComboxElementsNames;

    QComboBox *ComboxSolidType;
    QLineEdit *SolidSpecications;

    QComboBox *ComboxSolidNames;
    QComboBox *ComboxMaterialsNames;
    QLineEdit *LogicalVolumeData;

    QComboBox *ComboxLogicalVolumesNames;
    QComboBox *ComboxPositionsNames;
    QComboBox *ComboxRotationsNames;
    QLineEdit *PhysicalVolumeData;

    QLineEdit *TopVolumeData;

    QVector <QString> PositionsNames ;
    QVector <QString> RotationsNames ;
    QVector <QString> IsotopesNames ;
    QVector <QString> ElementsNames ;
    QVector <QString> MaterialsNames ;
    QVector <QString> SolidsNames ;
    QVector <QString> VolsNames ;
    QVector <QString> PhyVolsNames ;

    Ui::GeometryModellingEditor *ui;
    QString currentFile;
    Highlighter *highlighter;
};

#endif // GEOMETRYMODELLINGEDITOR_H
