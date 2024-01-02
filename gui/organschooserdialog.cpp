#include "organschooserdialog.h"
#include "ui_organschooserdialog.h"

OrgansChooserDialog::OrgansChooserDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::OrgansChooserDialog)
{
    ui->setupUi(this);
}

OrgansChooserDialog::~OrgansChooserDialog()
{
    delete ui;
}
