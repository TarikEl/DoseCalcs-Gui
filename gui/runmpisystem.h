#ifndef RUNMPISYSTEM_H
#define RUNMPISYSTEM_H

#include <QDialog>
#include "gui/filesManager.h"

namespace Ui {
class runmpisystem;
}

class runmpisystem : public QDialog
{
    Q_OBJECT

public:
    explicit runmpisystem(QWidget *parent = nullptr);
    ~runmpisystem();
    void getConfigurationData();
    void generateExeFile();
private slots:
    void on_pushButtonSaveToExe_clicked();

private:
    Ui::runmpisystem *ui;
    filesManager* filesManagerObj;

};

#endif // RUNMPISYSTEM_H
