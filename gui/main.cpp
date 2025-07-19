#include "gui/mainwindow.h"
#include <QApplication>
#include <QCoreApplication>

int NumberOfCPUCores;

qint64 pidOfPreviousNotKilled;

QString ExampleGeant4Code_source_dir_path ;
QString ExampleGeant4Code_build_dir_path ;
QString ExampleGeant4Code_ExecutionCommand ;
QString DoseCalcsCore_source_dir_path ;
QString DoseCalcsGui_source_dir_path ;
QString DoseCalcsCore_build_dir_path;
QString UserCurrentResultsDirPath;
QString DoseCalcs_build_dir_name;


QStringList CompleterWords;

QString ResultDirectoryName;
QString ScriptDirectoryName;
QString DoseCalcsExecutableName;
QString MergeExecutableName;
QString GraphExecutableName;
QString MacroFileName;
QString ResultFileName;
QString AnalysisInputFileName;
QString ReferenceFileName;
QString GraphsOutDirName;
QString DoseCalcsExecutingFileName;
QString MacroFilePath;
QString DoseCalcsDownloadURL;
QString RocksProcessNumber ;

QString GUIPackagesAndFilesDirName;
QString GUIConfigFileName;
QString GUIPackagesAndFilesDirPath ;
QString ICRPDATAPath ;

QString geant4_Lib_dir_path ;
QString DCMTK_Lib_dir_path ;
QString MPI_Lib_dir_path ;
QString CMAKE_Lib_dir_path ;
QString Root_Lib_dir_path ;

QString cmakeTruePath;

QString CENTOS_ROCKS_CLUSTER ;
QString MPI_USE ;
QString ROOT_USE ;
QString DCMTK_USE ;
QString GDML_USE;
QString VERBOSE_USE ;

QString ExeFileName;

QString Execution_setEventNumber ;
QString QueueName;
QString ExeDataText;
QString ConfigDataText ;

QVector <QString> MaterialCommands ;
QVector <QString> GeometryCommands ;
QVector <QString> VOXELCommands ;
QVector <QString> DICOMCommands ;
QVector <QString> CONSCommands ;
QVector <QString> PhysicsCommands ;
QVector <QString> SourceCommands ;
QVector <QString> RunAndScoreCommands ;
QVector <QString> AnalysisCommands ;

QStringList DoseCalcsQuantities ;
QMap <QString,QStringList> QuantitiesUnitsLists ;
QMap <QString,QMap <QString,double>> QuantitiesConversionFromDefault ;


QString OSNameVersion ;

void getOSName(){

    FILE *fp;
    char buffer[50] = " ";
    fp = popen("lsb_release -ds", "r");
    fgets(buffer, 50, fp);
    pclose(fp);
    printf("Welcome to DoseCalcs on %s",buffer);
    OSNameVersion = buffer;

    //if(OSNameVersion.toLower().contains("ubuntu")){QTextStream(stdout) << "Welcome from Ubuntu: (" << buffer << ")"<< "\n";}
    //else if(OSNameVersion.toLower().contains("centos")){QTextStream(stdout) << "Welcome from CentOs: (" << buffer << ")"<< "\n";}
    //else if(OSNameVersion.toLower().contains("windows")){QTextStream(stdout) << "Welcome from Windows: (" << buffer << ")"<< "\n";}
}

std::string getFileNameFromPath(std::string const & path, std::string const & delims = "/\\") {
    std::string const & filename = path.substr(path.find_last_of(delims) + 1);
    typename std::string::size_type const p(filename.find_last_of('.'));
    return p > 0 && p != std::string::npos ? filename.substr(0, p) : filename;
}
std::string getFileExt(const std::string& s) {
    size_t i = s.rfind('.', s.length());
    if (i != std::string::npos) {
        return(s.substr(i+1, s.length() - i));
    }
    return("");
}

extern QString TestPackagesPathsBeforToRun(){

    QString PackagesInfo = "The needed files and directories for DoseCalcs simulation:\n\n";

    if(!QFile::exists("/usr/bin/xterm")){
        PackagesInfo += "***** xterm (Not Found) - /usr/bin/xterm\n\n";
    }else{
        PackagesInfo += "***** xterm (Found) - /usr/bin/xterm \n\n";
    }

    if(!QFile::exists(DoseCalcsCore_build_dir_path+"/simulate")){
        PackagesInfo += "***** DoseCalcs Simulation Executable (Not Found) - "+DoseCalcsCore_build_dir_path+"/simulate"+"\n\n";
    }else{
        PackagesInfo += "***** DoseCalcs Simulation Executable (Found) - "+DoseCalcsCore_build_dir_path+"/simulate"+"\n\n";
    }

    if(!QFile::exists(geant4_Lib_dir_path+"/geant4.sh")){
        PackagesInfo += "**** GEANT4 : (Not Found) - "+geant4_Lib_dir_path+"/geant4.sh \n\n";
    }else{
        PackagesInfo += "**** GEANT4 : (Found) - "+geant4_Lib_dir_path+"/geant4.sh \n\n";
    }

    if(!QFile::exists(CMAKE_Lib_dir_path+"/cmake")){
        PackagesInfo += "**** cmake (Not Found) - "+CMAKE_Lib_dir_path+"/cmake\n\n";
    }else{
        PackagesInfo += "**** cmake (Found) - "+CMAKE_Lib_dir_path+"/cmake \n\n";
    }

    if(!QFile::exists(DoseCalcsCore_build_dir_path+"/"+MergeExecutableName)){
        PackagesInfo += "**** DoseCalcs Merging Executable (Not Found) - "+DoseCalcsCore_build_dir_path+"/"+MergeExecutableName+"\n\n";
    }else{
        PackagesInfo += "**** DoseCalcs Merging Executable (Found) - "+DoseCalcsCore_build_dir_path+"/"+MergeExecutableName+"\n\n";
    }

    if(!QFile::exists(DoseCalcsCore_build_dir_path+"/"+GraphExecutableName)){
        PackagesInfo += "*** DoseCalcs Analysis Executable (Not Found) - "+DoseCalcsCore_build_dir_path+"/"+GraphExecutableName+"\n\n";
    }else{
        PackagesInfo += "*** DoseCalcs Analysis Executable (Found) - "+DoseCalcsCore_build_dir_path+"/"+GraphExecutableName+"\n\n";
    }

    if(!QFile::exists(Root_Lib_dir_path+"/thisroot.sh")){
        PackagesInfo += "*** ROOT : (Not Found) - "+Root_Lib_dir_path+"/thisroot.sh \n\n";
    }else{
        PackagesInfo += "*** ROOT : (Found) - "+Root_Lib_dir_path+"/thisroot.sh \n\n";
    }

    if(!QFile::exists(MPI_Lib_dir_path+"/mpirun")){
        PackagesInfo += "*** mpirun : (Not Found) - "+MPI_Lib_dir_path+"/mpirun \n\n";
    }else{
        PackagesInfo += "*** mpirun : (Found) - "+MPI_Lib_dir_path+"/mpirun \n\n";
    }

    if(!QFile::exists("/usr/bin/nohup")){
        PackagesInfo += "*** nohup (Not Found) - /usr/bin/nohup \n\n";
    }else{
        PackagesInfo += "*** nohup (Found) - /usr/bin/nohup \n\n";
    }

    if(!QFile::exists(DCMTK_Lib_dir_path+"/cmake")){
        PackagesInfo += "*** DCMTK (Not Found) - "+DCMTK_Lib_dir_path+"\n\n";
    }else{
        PackagesInfo += "*** DCMTK (Found) - "+DCMTK_Lib_dir_path+"\n\n";
    }

    if(!QFile::exists("/usr/bin/scp")){
        PackagesInfo += "*** scp (Not Found) - /usr/bin/scp\n\n";
    }else{
        PackagesInfo += "*** scp (Found) - /usr/bin/scp \n\n";
    }

    return PackagesInfo;
}


int main(int argc, char *argv[])
{

    getOSName();

    cmakeTruePath = "cmake";

    ConfigDataText =
    "\n CMAKE_INSTALL_DIR \n "
    "GEANT4_INSTALL_DIR \n "
    "MPI_INSTALL_DIR \n "
    "ROOT_INSTALL_DIR \n "
    "DCMTK_INSTALL_DIR \n "
    "DoseCalcs_SOURCE_DIR \n\n "
    ""
    "DEFAULT_DoseCalcs_INPUTS \n\n"
    ""
    "CMAKE_DOWNLOAD_URL \n "
    "GEANT4_DOWNLOAD_URL \n "
    "XERCES_DOWNLOAD_URL \n "
    "MPI_DOWNLOAD_URL \n "
    "ROOT_DOWNLOAD_URL \n "
    "DCMTK_DOWNLOAD_URL \n\n "
    ""
    "VERBOSE_USE \n "
    "CENTOS_ROCKS_CLUSTER \n "
    "MPI_USE \n "
    "ROOT_USE \n "
    "DCMTK_USE \n ";

    NumberOfCPUCores = std::thread::hardware_concurrency();
    if(NumberOfCPUCores > 1){
        NumberOfCPUCores = NumberOfCPUCores - 1;
    }

    ResultDirectoryName = "Results";
    ScriptDirectoryName = "Scripts";
    DoseCalcsExecutableName = "simulate";
    MergeExecutableName = "merge";
    GraphExecutableName = "analysis";
    MacroFileName = "macros.mac";
    ResultFileName="ResultsData";
    ReferenceFileName="Ref1";
    GraphsOutDirName = "PlotsAndTables";
    DoseCalcsExecutingFileName = "c.sh";
    ExeFileName = "exe.sh";

    GUIConfigFileName = "Config";
    GUIPackagesAndFilesDirName = "PackagesAndFiles";

    GUIPackagesAndFilesDirPath = QDir::currentPath()+"/"+GUIPackagesAndFilesDirName;
    if(!QFile::exists(GUIPackagesAndFilesDirPath)){
        QDir::current().mkdir(GUIPackagesAndFilesDirName);
    }

    DoseCalcsCore_build_dir_path = "";
    DoseCalcs_build_dir_name = "core_build";
    UserCurrentResultsDirPath = "";
    ICRPDATAPath = GUIPackagesAndFilesDirPath+"/ICRPDATA";

    DoseCalcsDownloadURL = "https://github.com/TarikEl/DoseCalcs-Gui/archive/main.tar.gz";

    DoseCalcsQuantities=(QStringList()<<"AE"<<"AF"<<"SAF"<<"S"<<"H"<<"E");

    QuantitiesUnitsLists["SAF"] = (QStringList()<<"kg-1"<<"g-1");
    QuantitiesUnitsLists["AF"] = (QStringList()<<"");
    QuantitiesUnitsLists["AE"] = (QStringList()<<"MeV"<<"keV"<<"eV"<<"J");
    QuantitiesUnitsLists["AD"] = (QStringList()<<"MeV/kg"<<"Gy");
    QuantitiesUnitsLists["S"] = (QStringList()<<"MeV/kg"<<"nGy"<<"miGy"<<"mGy"<<"Gy"<<"kGy"<<"MGy");
    QuantitiesUnitsLists["H"] = (QStringList()<<"MeV/kg"<<"nGy"<<"miGy"<<"mGy"<<"Gy"<<"kGy"<<"MGy"<<"mSv"<<"Sv");
    QuantitiesUnitsLists["E"] = (QStringList()<<"MeV/kg"<<"nGy"<<"miGy"<<"mGy"<<"Gy"<<"kGy"<<"MGy"<<"mSv"<<"Sv");
    QuantitiesUnitsLists["T"] = (QStringList()<<"s"<<"min"<<"h"<<"d"<<"y");
    QuantitiesUnitsLists["A"] = (QStringList()<<"MBq"<<"Bq"<<"kBq");

    double MeV_to_J = 1.60218e-13;
    double Gy_to_Sv = 1. ;
    double Bq_to_MBq = 1e-6 ;

    QuantitiesConversionFromDefault["SAF"]["kg-1"] = 1.;
    QuantitiesConversionFromDefault["SAF"]["g-1"] = 1e+3;

    QuantitiesConversionFromDefault["AF"][""] = 1.;

    QuantitiesConversionFromDefault["AE"]["MeV"] = 1.;
    QuantitiesConversionFromDefault["AE"]["keV"] = 1e-3;
    QuantitiesConversionFromDefault["AE"]["eV"] = 1e-6;
    QuantitiesConversionFromDefault["AE"]["J"] = 1/MeV_to_J;

    QuantitiesConversionFromDefault["AD"]["MeV/kg"] = 1.;
    QuantitiesConversionFromDefault["AD"]["Gy"] = 1./MeV_to_J;

    QuantitiesConversionFromDefault["S"]["MeV/kg"] = 1.;
    QuantitiesConversionFromDefault["S"]["Gy"] = 1./MeV_to_J;
    QuantitiesConversionFromDefault["S"]["miGy"] = 1./(MeV_to_J*1e+6);
    QuantitiesConversionFromDefault["S"]["nGy"] = 1./(MeV_to_J*1e+9);
    QuantitiesConversionFromDefault["S"]["mGy"] = 1./(MeV_to_J*1e+3);
    QuantitiesConversionFromDefault["S"]["MGy"] = 1./(MeV_to_J*1e-6);
    QuantitiesConversionFromDefault["S"]["kGy"] = 1./(MeV_to_J*1e-3);

    QuantitiesConversionFromDefault["H"]["MeV/kg"] = 1.;
    QuantitiesConversionFromDefault["H"]["Gy"] = 1./MeV_to_J;
    QuantitiesConversionFromDefault["H"]["miGy"] = 1./(MeV_to_J*1e+6);
    QuantitiesConversionFromDefault["H"]["nGy"] = 1./(MeV_to_J*1e+9);
    QuantitiesConversionFromDefault["H"]["mGy"] = 1./(MeV_to_J*1e+3);
    QuantitiesConversionFromDefault["H"]["MGy"] = 1./(MeV_to_J*1e-6);
    QuantitiesConversionFromDefault["H"]["kGy"] = 1./(MeV_to_J*1e-3);
    QuantitiesConversionFromDefault["H"]["Sv"] = 1./(MeV_to_J*Gy_to_Sv);
    QuantitiesConversionFromDefault["H"]["mSv"] = 1./(MeV_to_J*1e+3*Gy_to_Sv);

    QuantitiesConversionFromDefault["E"]["MeV/kg"] = 1.;
    QuantitiesConversionFromDefault["E"]["Gy"] = 1./MeV_to_J;
    QuantitiesConversionFromDefault["E"]["miGy"] = 1./(MeV_to_J*1e+6);
    QuantitiesConversionFromDefault["E"]["nGy"] = 1./(MeV_to_J*1e+9);
    QuantitiesConversionFromDefault["E"]["mGy"] = 1./(MeV_to_J*1e+3);
    QuantitiesConversionFromDefault["E"]["MGy"] = 1./(MeV_to_J*1e-6);
    QuantitiesConversionFromDefault["E"]["kGy"] = 1./(MeV_to_J*1e-3);
    QuantitiesConversionFromDefault["E"]["Sv"] = 1./(MeV_to_J*Gy_to_Sv);
    QuantitiesConversionFromDefault["E"]["mSv"] = 1./(MeV_to_J*1e+3*Gy_to_Sv);

    QuantitiesConversionFromDefault["T"]["s"]   = 1.;
    QuantitiesConversionFromDefault["T"]["min"] = 1./(60);
    QuantitiesConversionFromDefault["T"]["h"]   = 1./(60*60);
    QuantitiesConversionFromDefault["T"]["d"]   = 1./(60*60*24);
    QuantitiesConversionFromDefault["T"]["y"]   = 1./(60*60*24*365);

    QuantitiesConversionFromDefault["A"]["Bq"]   = 1.;
    QuantitiesConversionFromDefault["A"]["MBq"] = 1./1e+6;
    QuantitiesConversionFromDefault["A"]["kBq"]   = 1./1e+3;

    //QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    //env.insert("QT_QPA_PLATFORM", "xcb");
    //QProcess p ;p.setProcessEnvironment(env);

    QApplication a(argc, argv);
    MainWindow w; // this is the main GUI object
    //w.show();
    w.showMaximized();

    a.setWindowIcon(QIcon(QDir(QCoreApplication::applicationDirPath()).filePath(GUIPackagesAndFilesDirName+"/AppIcon.png")));
    return a.exec();
}

void showResultsOutput(QString text, int level){

    if(level == 4){
        QTextStream(stdout) << text << endl;
    }
    //ui->outputTextConsole->append("\n"+text);
}
