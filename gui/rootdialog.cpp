#include "gui/rootdialog.h"
#include "ui_rootdialog.h"

RootDialog::RootDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RootDialog)
{
    ui->setupUi(this);
}

RootDialog::~RootDialog()
{
    delete ui;
}
