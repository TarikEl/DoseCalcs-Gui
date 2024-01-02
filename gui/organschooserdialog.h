#ifndef ORGANSCHOOSERDIALOG_H
#define ORGANSCHOOSERDIALOG_H

#include <QDialog>

namespace Ui {
class OrgansChooserDialog;
}

class OrgansChooserDialog : public QDialog
{
    Q_OBJECT

public:
    explicit OrgansChooserDialog(QWidget *parent = 0);
    ~OrgansChooserDialog();

private:
    Ui::OrgansChooserDialog *ui;
};

#endif // ORGANSCHOOSERDIALOG_H
