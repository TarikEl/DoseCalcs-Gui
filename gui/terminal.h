#ifndef TERMINAL_H
#define TERMINAL_H

#include <QWidget>
#include <QProcess>

class terminal : public QWidget
{
    Q_OBJECT
public:
    explicit terminal(QWidget *parent = nullptr);
    void executeCommand(QString);
    void executeLocalCommand(QString);
    void executeNormalCommand(QString);
    void executeNormalCommandWithoutKillingPre(QString);
private:
    QProcess process;

signals:

};

#endif // TERMINAL_H
