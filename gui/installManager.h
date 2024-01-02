#ifndef INSTALLMANAGER_H
#define INSTALLMANAGER_H

#include <QString>
#include <QDialog>

namespace Ui {
class installManager;
}

class installManager : public QDialog
{
    Q_OBJECT

public:

    explicit installManager(QWidget *parent = 0);
    ~installManager();

    void get_Paths_from_PathsFile();
    void create_PrerequestGeant4Install_sh_file();
    void create_Geant4Install_sh_file();
    void create_G4mpi_Install_sh_file();
    void create_RootInstall_sh_file();
    void create_AppInstall_sh_file();
    void get_file_using_Downloading(QString Url, QString download_Directory, QString FileName);

    //ReadableByteChannel readableByteChannel ;  // an open url channel urlConnection to read the byte from when download
    //FileOutputStream fileOuputStreamFromDownload ;
    void download_data_UsingNIO(QString URL_string, QString file_path_string);

    QString Geant4_text_shFile ;
    QString App_text_shFile;
    QString g4Mpi_text_shFile;
    QString Root_text_shFile;
    QString Prequest_text_shFile;

    QString Geant4_sh_name ;
    QString Root_sh_name ;
    QString App_sh_name ;
    QString Prequests_sh_name ;
    QString G4mpi_sh_name ;

    int coreNumber ;

    QString Cmake_Command ;

    QString workSpace_dir_path ;
    QString App_source_dir_path ;
    QString Geant4_source_dir_path ;
    QString Root_source_dir_path ;
    QString G4mpi_source_dir ;

    QString geant4_Lib_dir_path ;
    QString g4Mpi_Lib_dir_path ;
    QString g4Root_Lib_dir_path ;

    QString mpi_install_dir ;
    QString XERCESC_Lib_dir_path ;

private:

    Ui::installDialog *ui;
};

#endif // INSTALLMANAGER_H
