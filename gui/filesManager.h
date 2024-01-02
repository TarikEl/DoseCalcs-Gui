#ifndef FILESMANAGER_H
#define FILESMANAGER_H

#include <QVector>
#include <QString>
#include <QTextEdit>
#include <QDialog>

class filesManager
{
public:
    filesManager();
    ~filesManager();

    QString ReadTextFromFileInOneString(QString);
    QVector<QString> ReadTextFromFile(QString ) ;
    QMap<QString,QString> ReadLinesFromFileWithFirstWordIndicator(QString ) ;
    QMap<QString,QVector<QString>> ReadLinesFromFileWithFirstWordIndicatorAndVector(QString ) ;
    QMap<double,double> ReadLinesFromFileWithFirstDoubleValueIndicator(QString);
    QVector< QPair<QString,QString>> ReadTextFromFileInQStringList(QString) ;
    QVector< QPair<QString,QString>> ReadTextFromQStringInQStringList(QString) ;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> Read_final_result_file(QString,int);
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> Read_Comparison_file(QString,int);
    QMap<QString,QMap<QString,QMap<QString,QMap<double,double >>>> Read_CrossSection_file(QString);
    QMap<QString,QMap<QString,double>> ReadRegionsData(QString ) ;

    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> StandartDeviation;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>>> getStandartDeviationMap(){return StandartDeviation ;}

    void WriteTextToFile(QString FileName , QString TextToWrite) ;
    void showResultsOutput(QString, int);


    // called from open buttons and install buttons
    QString GetChoosenDirFromDialog(QString, int);

    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>> getResultParticleSourceEnergyTime(){return ResultParticleSourceEnergyTime ;}
    QMap<QString,QMap<QString,QMap<QString,QVector<double>>>> getResEnergies(){return ResEnergies ;}
    double getMinValForLog(){return MinValForLog ;}

    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> ResultQuantityGeometryRadioTracerSourceTargetValues ;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> getResultQuantityGeometryRadioTracerSourceTargetValues(){return ResultQuantityGeometryRadioTracerSourceTargetValues ;}

    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> ResultQuantityGeometryRadioTracerSourceTargetStandardDeviation ;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> getResultQuantityGeometryRadioTracerSourceTargetStandardDeviation(){return ResultQuantityGeometryRadioTracerSourceTargetStandardDeviation ;}

    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> ResultQuantityGeometryRadioTracerSourceTargetRelativeStandardDeviation ;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> getResultQuantityGeometryRadioTracerSourceTargetRelativeStandardDeviation(){return ResultQuantityGeometryRadioTracerSourceTargetRelativeStandardDeviation ;}

    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> ReferenceQuantityGeometryRadioTracerSourceTargetValues ;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<QString,double>>>>> getReferenceQuantityGeometryRadioTracerSourceTargetValues(){return ReferenceQuantityGeometryRadioTracerSourceTargetValues ;}

    QMap<QString,QMap<QString,QMap<QString,double>>> RegionParameterValueMap;
    QMap<QString,QMap<QString,QMap<QString,double>>> getRegionParameterValueMap(){return RegionParameterValueMap ;}

private:
    double MinValForLog;
    QMap<QString,QMap<QString,QMap<QString,QMap<QString,QMap<double,double>>>>> ResultParticleSourceEnergyTime;
    QMap<QString,QMap<QString,QMap<QString,QVector<double>>>> ResEnergies ;
    //QVector<QString> OrganNamesToScore ;

};


#endif // FILESMANAGER_H
