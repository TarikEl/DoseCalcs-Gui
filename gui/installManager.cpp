#include "installManager.h"


installManager::installManager(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::installDialog)
{
    ui->setupUi(this);

}

installManager::~installManager()
{
    delete ui;
}

