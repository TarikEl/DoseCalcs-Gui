#ifndef DICOMREADERSEG_H
#define DICOMREADERSEG_H

#include <QDialog>

namespace Ui {
class dicomreaderseg;
}

class dicomreaderseg : public QDialog
{
    Q_OBJECT

public:
    explicit dicomreaderseg(QWidget *parent = nullptr);
    ~dicomreaderseg();

private:
    Ui::dicomreaderseg *ui;
};

#endif // DICOMREADERSEG_H
