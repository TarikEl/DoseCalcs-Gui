#ifndef GRAPHSDATAGENERATION_H
#define GRAPHSDATAGENERATION_H

#include <QVector>
#include <QString>
#include <QMap>
#include "filesManager.h"

class GraphsDataGeneration
{
public:

    GraphsDataGeneration();
    ~GraphsDataGeneration();

    QString ResultDirectoryPath;
    QString appBuildDir;
    QString ResultDirectoryName;
    QString CompareFileName;
    QString MacroFileName;

    QString Compare_Type;
    QString Graphs_data;
    QString Graphs_Ext;
    QString Compare_file ;
    QString Score_variables;
    QString CompareReferenceName;

    int number_of_energies;
    int number_of_energies_in_compareFile ;

    //QMap <QString, QMap <QString, QMap <QString,QMap <double, QVector(double)> > > > graphs;

    QVector<QString> Organs_to_score;

    void Analyse_Main();
    void Data_initialization();

    QString ScoreVariableUnit;

    void get_data_from_Maps_simulationAndAnalyseWithCompare(QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>> table, QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>> CompareTable , int whatKindOf);
    void get_data_from_Maps_simulationAndAnalyse(QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>> table, int whatKindOf);
    QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>> Read_Comparison_file();
    QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>> Read_final_result_file();
    QVector <QString> readFilesAndGetFilesNamesVector();

    void getArguments(QString, QString, QString, QString, QString, QVector<QString> , QString);
    void showResultsOutput(QString text, int level);


    QMap <QString,QMap<QString,QMap<double,double>>> ResSelfGraphs;
    QMap <QString,QMap<QString,QMap<double,double>>> getResSelfGraphs(){return ResSelfGraphs;}
    QMap <QString,QMap<QString,QMap<QString,QMap<double,double>>>> ResCrossGraphs;
    QMap <QString,QMap<QString,QMap<QString,QMap<double,double>>>> getResCrossGraphs(){return ResCrossGraphs;}
    QMap <QString,QMap<QString,QMap<double,double>>> ResRefSelfGraphs;
    QMap <QString,QMap<QString,QMap<double,double>>> getResRefSelfGraphs(){return ResRefSelfGraphs;}
    QMap <QString,QMap<QString,QMap<QString,QMap<double,double>>>> ResRefCrossGraphs;
    QMap <QString,QMap<QString,QMap<QString,QMap<double,double>>>> getResRefCrossGraphs(){return ResRefCrossGraphs;}


    QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>> ErrorTables;

    double DefaultErrorDistance ;

    filesManager* filesManagerObject;
};

#endif // GRAPHSDATAGENERATION_H
