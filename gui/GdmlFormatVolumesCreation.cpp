#include "GdmlFormatVolumesCreation.h"
#include "ui_GdmlFormatVolumesCreation.h"

#include "QFileDialog"
#include "QFile"
#include "QMessageBox"
#include "QTextStream"
#include "QFont"
#include "QFontDialog"

GDMLFormatVolumesCreation::GDMLFormatVolumesCreation(QWidget *parent) : QDialog(parent), ui(new Ui::GDMLFormatVolumesCreation)
{
        ui->setupUi(this);
}

GDMLFormatVolumesCreation::~GDMLFormatVolumesCreation()
{
    delete ui;
}


void GDMLFormatVolumesCreation::newDocument()
{
    currentFile.clear();
    ui->textEdit->setText(QString());
}

void GDMLFormatVolumesCreation::open()
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

void GDMLFormatVolumesCreation::save()
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

void GDMLFormatVolumesCreation::saveAs()
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

void GDMLFormatVolumesCreation::print()
{
#if defined(QT_PRINTSUPPORT_LIB) && QT_CONFIG(printer)
    QPrinter printDev;
#if QT_CONFIG(printdialog)
    QPrintDialog dialog(&printDev, this);
    if (dialog.exec() == QDialog::Rejected)
        return;
#endif // QT_CONFIG(printdialog)
    ui->textEdit->print(&printDev);
#endif // QT_CONFIG(printer)
}

void GDMLFormatVolumesCreation::exit()
{
    QCoreApplication::quit();
}

void GDMLFormatVolumesCreation::copy()
{
#if QT_CONFIG(clipboard)
    ui->textEdit->copy();
#endif
}

void GDMLFormatVolumesCreation::cut()
{
#if QT_CONFIG(clipboard)
    ui->textEdit->cut();
#endif
}

void GDMLFormatVolumesCreation::paste()
{
#if QT_CONFIG(clipboard)
    ui->textEdit->paste();
#endif
}

void GDMLFormatVolumesCreation::undo()
{
     ui->textEdit->undo();
}

void GDMLFormatVolumesCreation::redo()
{
    ui->textEdit->redo();
}

void GDMLFormatVolumesCreation::selectFont()
{
    bool fontSelected;
    QFont font = QFontDialog::getFont(&fontSelected, this);
    if (fontSelected)
        ui->textEdit->setFont(font);
}
