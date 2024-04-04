#include "gui/terminal.h"
#include "QTextStream"
#include "QFile"
#include "QMessageBox"

extern qint64 pidOfPreviousNotKilled;

terminal::terminal(QWidget *parent) : QWidget(parent)
{
    if(!QFile::exists("/usr/bin/xterm")){
        QMessageBox::information(this, tr("xterm not found")," Install \"xterm\" by typing the command below in your terminal: \n\"sudo apt-get install xterm\" on ubuntu, \n or \n\"sudo yum install xterm\"");
    }
}
void terminal::executeCommand(QString command){

    QProcess::startDetached("kill", {QString::number(pidOfPreviousNotKilled)});

    //QStringList qsl = {"-into", QString::number(winId()), "-bg", "black", "-fg", "green", "-maximized", "-hold", "-e", "sh", command};
    QStringList qsl = {"-into", QString::number(winId()), "-bg", "black", "-fg", "green" , "-hold" ,"-geometry", "800x800+620+2", "-e", "bash", command};
    //800x800 is set as a big value, and is known that terminal will be just fill the mother, then it will be suitable to the
    process.setProgram("xterm");
    process.setArguments(qsl);

    process.startDetached(&pidOfPreviousNotKilled);

    QTextStream(stdout) << " pidOfPreviousNotKilled: " << pidOfPreviousNotKilled << "\n";
}
void terminal::executeLocalCommand(QString command){

    //QStringList qsl = {"-into", QString::number(winId()), "-bg", "black", "-fg", "green", "-maximized", "-hold", "-e", "sh", command};
    QStringList qsl = {"-into", QString::number(winId()), "-bg", "black", "-fg", "green" , "-hold" ,"-geometry", "800x800+620+2", "-e", "bash", command};
    //800x800 is set as a big value, and is known that terminal will be just fill the mother, then it will be suitable to the
    process.setProgram("xterm");
    process.setArguments(qsl);
    process.start();
    //while(process.pid()> 0){}
    process.waitForFinished(2000); // until the file is read and writed data to is
    //process.waitForFinished();
}
void terminal::executeNormalCommandWithoutKillingPre(QString command){

    //QProcess::startDetached("kill", {QString::number(pidOfPreviousNotKilled)});

    //QStringList qsl = {"-into", QString::number(winId()), "-bg", "black", "-fg", "green", "-maximized", "-hold", "-e", "sh", command};
    QStringList qsl = {"-into", QString::number(winId()), "-bg", "black", "-fg", "green" , "-hold" ,"-geometry", "800x800+620+2", "-e", command};
    //800x800 is set as a big value, and is known that terminal will be just fill the mother, then it will be suitable to the
    process.setProgram("xterm");
    process.setArguments(qsl);

    process.startDetached(&pidOfPreviousNotKilled);

    QTextStream(stdout) << " pidOfPreviousNotKilled: " << pidOfPreviousNotKilled << "\n";
}
