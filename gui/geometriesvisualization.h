#ifndef GEOMETRIESVISUALIZATION_H
#define GEOMETRIESVISUALIZATION_H

#include <QDialog>
#include "gui/qcustomplot.h"
#include "gui/filesManager.h"

namespace Ui {
class geometriesvisualization;
}

class geometriesvisualization : public QDialog
{
    Q_OBJECT

public:
    explicit geometriesvisualization(QWidget *parent = nullptr);
    ~geometriesvisualization();

    void getInputData();
    void ReadVoxelsFilesAndFillMaps();
    void showTETrahedrons();

private slots:
    void on_GenerateMap_clicked();

    void on_Save_clicked();

    void on_AddVoxelsDataFile_clicked();

    void on_pushButtonReadVoxelsData_clicked();

    void on_horizontalSlider_valueChanged(int value);

    void on_radioButtonVoxIDs_clicked(bool checked);

    void on_radioButtonDicom_toggled(bool checked);

    void on_pushButtonOpenResDir_clicked();

    void on_radioButtonPET_clicked(bool checked);

    void on_Axis_activated(int index);

    void on_horizontalSliderSliceID_sliderMoved(int position);

    void on_comboBoxColorGradien_currentTextChanged(const QString &arg1);

private:

    filesManager* fileManagerObject;
    //double colormaxvalue;


    //QMap <int , QMap <int , QMap <int , double>>> X_ZY;
    //QMap <int , QMap <int , QMap <int , double>>> Y_ZX;
    //QMap <int , QMap <int , QMap <int , double>>> Z_YX;
    double *** X_ZY;
    double *** Y_ZX;
    double *** Z_YX;


    QString Axis;
    int SliceID;
    int VoxXNumber;
    int VoxYNumber;
    int VoxZNumber;
    double VoxXSize;
    double VoxYSize;
    double VoxZSize;

    QString VoxelsDataFilePath;
    QString FileExtToSave;

    Ui::geometriesvisualization *ui;

public:

};

#endif // GEOMETRIESVISUALIZATION_H
