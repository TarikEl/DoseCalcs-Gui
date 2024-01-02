#ifndef ROOTDIALOG_H
#define ROOTDIALOG_H

#include <QDialog>

namespace Ui {
class RootDialog;
}

class RootDialog : public QDialog
{
    Q_OBJECT

public:
    explicit RootDialog(QWidget *parent = nullptr);
    ~RootDialog();

private:
    Ui::RootDialog *ui;
};

#endif // ROOTDIALOG_H
