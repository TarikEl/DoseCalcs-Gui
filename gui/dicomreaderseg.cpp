#include "gui/dicomreaderseg.h"
#include "gui/ui_dicomreaderseg.h"

/*
#include <QtWidgets>
#include <vtkSmartPointer.h>
#include <vtkImageViewer2.h>
#include <vtkDICOMImageReader.h>
#include <vtkImageData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkInteractorStyleDrawPolygon.h>
#include <dcmtk/dcmdata/dctk.h>
#include <dcmtk/dcmimgle/dcmimage.h>
*/

dicomreaderseg::dicomreaderseg(QWidget *parent) :QDialog(parent),ui(new Ui::dicomreaderseg)
{
    ui->setupUi(this);
}

dicomreaderseg::~dicomreaderseg()
{
    delete ui;
}

