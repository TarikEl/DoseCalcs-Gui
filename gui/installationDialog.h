#ifndef INSTALLATIONMIALOG_H
#define INSTALLATIONMIALOG_H

#include <QDialog>
#include <QProcess>
#include <QCompleter>

#include "gui/filesManager.h"
#include "gui/httpdownload.h"
#include "gui/highlighter.h"
#include <QCloseEvent>
#include <QKeyEvent>

namespace Ui {
class InstallationDialog;
}

class InstallationDialog : public QDialog
{
    Q_OBJECT

public:
    explicit InstallationDialog(QWidget *parent = 0);
    ~InstallationDialog();

    void keyPressEvent(QKeyEvent*);
    void closeEvent(QCloseEvent *event);  // show prompt when user wants to close app
    void getConfigurationDataForInstallation();
    void create_Xerces_and_PrerequestGeant4Install_sh_file();
    void create_Geant4Install_sh_file();
    void create_MPIInstall_sh_file();
    void create_DCMTK_Install_sh_file();
    void create_RootInstall_sh_file();
    void create_RootInstall_Dir_file();
    void create_AppInstall_sh_file();

    void get_file_using_Downloading(QString Url, QString download_Directory, QString FileName);
    void getPathsFromEditLines();

    void DownloadFileToDirectory(QString,QString);

    void execShProcess(QString);
    void execBackgroundProcess(QString);

    QString GetChoosenDirFromDialog(int type);
    void showResultsOutput(QString str, int num);

    //ReadableByteChannel readableByteChannel ;  // an open url channel urlConnection to read the byte from when download
    //FileOutputStream fileOuputStreamFromDownload ;
    void download_data_UsingNIO(QString URL_string, QString file_path_string);

    QString Xerces_dir_src, Xerces_dir_install;
    QString geant4_dir_src , geant4_dir_build , geant4_dir_install;
    QString MPI_dir_src, MPI_dir_install;
    QString CMAKE_dir_src, CMAKE_dir_install;
    QString DCMTK_dir_src , DCMTK_dir_build,DCMTK_dir_install;
    QString ROOT_dir_src , ROOT_dir_build , ROOT_dir_install ;

    QCompleter *completer1 = nullptr;
    void setCompleters();

    int EditFlag;

    QString PackageName;

    QString geant4_package_name;
    QString MPI_package_name;
    QString CMAKE_package_name;
    QString xerces_package_name;
    QString Root_package_name ;
    QString DCMTK_package_name;

    QString Geant4_text_shFile;
    QString MPI_text_shFile;
    QString CMAKE_text_shFile;
    QString App_text_shFile;
    QString DCMTK_text_shFile;
    QString Root_text_shFile;
    QString Prequest_text_shFile;
    QString InstallAll_text_shFile;

    QString Geant4_sh_path ;
    QString MPI_sh_path ;
    QString CMAKE_sh_path ;
    QString DCMTK_sh_path ;
    QString Root_sh_path ;
    QString App_sh_path ;
    QString Supplementary_sh_path ;
    QString Prequests_sh_path ;
    QString AllPackagesInstall_sh_path;

    int coreNumber ;

    QString Cmake_Command ;

    QString chosen_Dir;

    QString Geant4_Url_String ;
    //QString Supplementary_Url_String ;
    QString MPI_Url_String;
    QString CMAKE_Url_String;
    QString Xerces_Url_String;
    QString Root_Url_String ;
    QString DCMTK_Url_String ;

    QString Geant4_zip_path ;
    QString MPI_zip_path;
    QString CMAKE_zip_path;
    QString Xerces_zip_path;
    QString DCMTK_zip_path ;
    QString Root_source_dir_path ;

    QString Xerces_install_dir_path;
    QString Geant4_install_dir_path ;
    QString MPI_install_dir_path;
    QString CMAKE_install_dir_path;
    QString Root_install_dir_path ;
    QString DCMTK_install_dir ;

private slots:

    void on_installXercesButton_clicked();

    void on_InstallGeant4Button_clicked();

    void on_InstallDCMTKButton_clicked();

    void on_installCodeGeant4Button_clicked();

    void on_DownloadGeant4Button_clicked();

    void on_InstallMPIButton_clicked();

    void on_DownloadMPIButton_clicked();

    void on_ReadPathFilepushButton_clicked();

    void on_EditPathsFilepushButton_clicked();

    void on_DownloadROOTButton_clicked();

    void on_InstallROOTButton_clicked();

    void on_openAndInstallPrerequesiteBtn_clicked();

    void on_openAndInstallGeant4Btn_clicked();

    void on_openAndInstallMPIBtn_clicked();

    void on_openAndInstallDCMTKBtn_clicked();

    void on_openAndInstallROOTBtn_clicked();

    void on_openAndInstallAPPBtn_clicked();

    void on_DownloadDCMTKButton_clicked();

    void on_pushButtonInstallAllPrerequisites_clicked();

    void on_pushButtonGenerateInstallerFoAll_clicked();

    void on_pushButton_7_clicked();

    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_5_clicked();

    void on_checkBoxUseCompiledSource_clicked();

    void on_pushButton_8_clicked();

    void on_pushButtonSaveInput_clicked();

    void on_pushButtonGenerate_clicked();

    void on_pushButtonGEANTGenerate_clicked();

    void on_pushButtonMPIGenerate_clicked();

    void on_pushButtonDCMTKGenerate_clicked();

    void on_pushButtonROOTGenerate_clicked();

    void on_pushButtonLoadForAll_clicked();

    void on_pushButtonGenerateDoseCalcsBuilder_clicked();

    void on_pushButton_9_clicked();

    void on_pushButton_6_clicked();

    void on_ClearTerButton_clicked();

    void on_pushButtonUpdatePaths_clicked();

    void on_pushButtonOpenPackagesFilesDir_clicked();

    void on_pushButtonDownloadSupplement_clicked();

private:
    Highlighter* highlighter ;
    QProcess *process;
    Ui::InstallationDialog *ui;
    filesManager* filesManagerObj;
    HttpDownload* HttpDownloadObj;
};

#endif // INSTALLATIONMIALOG_H
