#include "gui/runmpisystem.h"
#include "ui_runmpisystem.h"

extern QString geant4_Lib_dir_path ;
extern QString DCMTK_Lib_dir_path ;
extern QString MPI_Lib_dir_path ;
extern QString CMAKE_Lib_dir_path ;
extern QString Root_Lib_dir_path ;
extern QString DoseCalcs_source_dir_path ;
extern QString DoseCalcsCore_build_dir_path;
extern QString UserCurrentResultsDirPath;
extern QString ResultDirectoryName;

extern QString GUIConfigFileName;
extern QString GUIPackagesAndFilesDirPath;

extern QString DoseCalcsExecutableName;
extern QString MacroFileName;
extern QString Execution_setEventNumber;

extern QString MPI_USE ;
extern QString ROOT_USE ;
extern QString DCMTK_USE ;
extern QString VERBOSE_USE ;
extern QString GDML_USE;

extern QString QueueName;
extern QString ExeDataText;
extern QString ExeFileName;

extern QString RocksProcessNumber ;

runmpisystem::runmpisystem(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::runmpisystem)
{
    ui->setupUi(this);
    getConfigurationData();
}

runmpisystem::~runmpisystem()
{
    delete ui;
}

void runmpisystem::getConfigurationData(){

    QMap<QString, QString> lines = filesManagerObj->ReadLinesFromFileWithFirstWordIndicator(GUIPackagesAndFilesDirPath+"/"+GUIConfigFileName);

    //showResultsOutput( "\n-------------------------------------------- Paths and Urls readed from the "+ GUIPackagesAndFilesDirPath+"/"+GUIConfigFileName , 4);

    CMAKE_Lib_dir_path = lines["CMAKE_INSTALL_DIR"];
    geant4_Lib_dir_path = lines["GEANT4_INSTALL_DIR"];
    MPI_Lib_dir_path = lines["MPI_INSTALL_DIR"];
    Root_Lib_dir_path = lines["ROOT_INSTALL_DIR"];
    DCMTK_Lib_dir_path = lines["DCMTK_INSTALL_DIR"];

}


void runmpisystem::generateExeFile()
{


}

void runmpisystem::on_pushButtonSaveToExe_clicked()
{

    QueueName = ui->lineEditQueueName->text();
    ExeDataText = "#!/bin/bash \n"
                  "#$ -S /bin/bash \n"
                  "#$ -cwd \n"
                  "#$ -j y \n"
                  "#$ -q  "+ QueueName + " \n"
                  "#$ -N Geant4.sh \n"
                  "#$ -pe mpi "+ Execution_setEventNumber + " \n"
                  "#$ -M imttarikk@gmail.com \n"
                  "#$ -m e \n"
                  ". "+geant4_Lib_dir_path+"/bin/geant4.sh \n"+
                  MPI_Lib_dir_path + "/bin/mpirun -np "+ Execution_setEventNumber + " " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutableName + " B " + MacroFileName + " " + Execution_setEventNumber +" \n";

    filesManagerObj->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+ExeFileName , ExeDataText);
}
