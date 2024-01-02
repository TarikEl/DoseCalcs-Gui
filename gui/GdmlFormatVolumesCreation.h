#ifndef GDMLFORMATVOLUMESCREATION_H
#define GDMLFORMATVOLUMESCREATION_H

#include <QDialog>

namespace Ui {
class GDMLFormatVolumesCreation;
}

class GDMLFormatVolumesCreation : public QDialog
{
    Q_OBJECT

public:
    explicit GDMLFormatVolumesCreation(QWidget *parent = 0);
    ~GDMLFormatVolumesCreation();


private slots:
    void newDocument();

    void open();

    void save();

    void saveAs();

    void print();

    void exit();

    void copy();

    void cut();

    void paste();

    void undo();

    void redo();

    void selectFont();

    void setFontBold(bool bold);

    void setFontUnderline(bool underline);

    void setFontItalic(bool italic);

    void about();

private:
    Ui::GDMLFormatVolumesCreation *ui;
    QString currentFile;
};

#endif // GDMLFORMATVOLUMESCREATION_H
