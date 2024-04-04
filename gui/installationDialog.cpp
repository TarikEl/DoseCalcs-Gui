#include "gui/installationDialog.h"
#include "gui/ui_installationDialog.h"

//#include<QtGlobal>

#include<QFileDialog>
#include<QTextStream>
#include<QProcess>
#include<QThread>
#include<QInputDialog>
#include<QStringListModel>
#include<QDialogButtonBox>
#include<gui/terminal.h>
//#include <QMessageBox>

extern QStringList CompleterWords;

extern QString geant4_Lib_dir_path ;
extern QString DCMTK_Lib_dir_path ;
extern QString MPI_Lib_dir_path ;
extern QString CMAKE_Lib_dir_path ;
extern QString Root_Lib_dir_path ;
extern QString DoseCalcsCore_source_dir_path ;
extern QString ExampleGeant4Code_source_dir_path ;
extern QString ExampleGeant4Code_build_dir_path ;
extern QString ExampleGeant4Code_ExecutionCommand ;
extern QString DoseCalcsCore_build_dir_path;
extern QString DoseCalcs_build_dir_name;

extern QString UserCurrentResultsDirPath;
extern QString ResultDirectoryName;

extern QString DoseCalcsExecutingFileName;

extern QString GUIConfigFileName;
extern QString GUIPackagesAndFilesDirPath;
extern QString GUIPackagesAndFilesDirName;

extern QString CENTOS_ROCKS_CLUSTER ;
extern QString MPI_USE ;
extern QString ROOT_USE ;
extern QString DCMTK_USE ;
extern QString VERBOSE_USE ;
extern QString GDML_USE;

extern QString cmakeTruePath;

extern QString ConfigDataText;

extern int NumberOfCPUCores ;

extern QString OSNameVersion ;

extern QString TestPackagesPathsBeforToRun();

InstallationDialog::InstallationDialog(QWidget *parent) : QDialog(parent), ui(new Ui::InstallationDialog)
{
    ui->setupUi(this);
    
    filesManagerObj = new filesManager;
    
    PackageName = "Geant4";

    xerces_package_name = "xerces-c-3.1.2";
    geant4_package_name = "xerces-c-3.1.2";
    MPI_package_name = "xerces-c-3.1.2";
    DCMTK_package_name = "xerces-c-3.1.2";
    Root_package_name = "xerces-c-3.1.2";
    
    Geant4_text_shFile ="" ;
    App_text_shFile="" ;
    DCMTK_text_shFile="" ;
    Root_text_shFile="" ;
    Prequest_text_shFile="" ;
    
    //NumberOfCPUCores = QThread::idealThreadCount();
    Cmake_Command = "" ;
    
    //showResultsOutput("NumberOfCPUCores : " + QString::number(NumberOfCPUCores) , 4);

    //GUIPackagesAndFilesDirPath+"/"+GUIConfigFileName = GUIPackagesAndFilesDirPath+"/"+GUIConfigFileName;
    
    Geant4_zip_path = "";
    MPI_zip_path = "";
    Xerces_zip_path = "";
    DCMTK_zip_path = "";
    DoseCalcsCore_source_dir_path = GUIPackagesAndFilesDirPath+"/DoseCalcsCore";
    
    Geant4_install_dir_path = "";
    DCMTK_install_dir = "";
    Root_install_dir_path = "";
    Xerces_install_dir_path = "";
    MPI_install_dir_path = "";
    
    geant4_Lib_dir_path = "";
    DCMTK_Lib_dir_path = "";
    Root_Lib_dir_path = "";
    MPI_Lib_dir_path = "";
    
    Geant4_Url_String = "https/....";
    //Supplementary_Url_String = "https/....";
    MPI_Url_String = "https/....";
    Root_Url_String = "https/....";
    DCMTK_Url_String = "https/....";
    Xerces_Url_String = "https/....";
    
    // first read the paths and Urls
    getConfigurationDataForInstallation();

    Supplementary_sh_path = GUIPackagesAndFilesDirPath+"/DownloadSupplementary_installer.sh";
    Prequests_sh_path = GUIPackagesAndFilesDirPath+"/Prerequesite_installer.sh";
    Geant4_sh_path  = GUIPackagesAndFilesDirPath+"/geant4_installer.sh";
    MPI_sh_path  = GUIPackagesAndFilesDirPath+"/MPI_installer.sh";
    Root_sh_path = GUIPackagesAndFilesDirPath+"/Root_installer.sh";
    DCMTK_sh_path  = GUIPackagesAndFilesDirPath+"/DCMTK_installer.sh";
    App_sh_path = GUIPackagesAndFilesDirPath+"/App_installer.sh";
    Geant4Example_sh_path = GUIPackagesAndFilesDirPath+"/Geant4Example_installer.sh";
    AllPackagesInstall_sh_path  = GUIPackagesAndFilesDirPath+"/AllPackages_installer.sh";

    ui->groupBox_2->setVisible(false);
    ui->groupBox_4->setVisible(false);

    highlighter = new Highlighter(ui->textEditInput->document());

    showResultsOutput("Installation dir : " + GUIPackagesAndFilesDirPath , 0);

    QString BashCommandsForExecuting = "#! /bin/bash \n "
                                           "cd "+DoseCalcsCore_build_dir_path+"\n"
                                           "bash \n ";
    filesManagerObj->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);
    execShProcess(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);

    setCompleters();
}

InstallationDialog::~InstallationDialog()
{
    UserCurrentResultsDirPath = DoseCalcsCore_build_dir_path + "/" + ResultDirectoryName;
    delete ui;
    delete filesManagerObj;
}

// open the files of paths and urls needed to download and open and install the required package
void InstallationDialog::getConfigurationDataForInstallation(){
    
    QMap<QString, QString> lines = filesManagerObj->ReadLinesFromFileWithFirstWordIndicator(GUIPackagesAndFilesDirPath+"/"+GUIConfigFileName);
    
    //showResultsOutput( "\n-------------------------------------------- Paths and Urls readed from the "+ GUIPackagesAndFilesDirPath+"/"+GUIConfigFileName , 4);
    
    CMAKE_Lib_dir_path = lines["CMAKE_INSTALL_DIR"];
    geant4_Lib_dir_path = lines["GEANT4_INSTALL_DIR"];
    MPI_Lib_dir_path = lines["MPI_INSTALL_DIR"];
    Root_Lib_dir_path = lines["ROOT_INSTALL_DIR"];
    DCMTK_Lib_dir_path = lines["DCMTK_INSTALL_DIR"];

    DoseCalcsCore_source_dir_path = lines["DoseCalcs_SOURCE_DIR"];
    if(!QFile::exists(DoseCalcsCore_source_dir_path)){ // create install dir for DCMTK
        DoseCalcsCore_source_dir_path = GUIPackagesAndFilesDirPath+"/DoseCalcsCore";
        showResultsOutput("Default DoseCalcs source directory " + DoseCalcsCore_source_dir_path + " is used", 4);
    }

    ExampleGeant4Code_source_dir_path = lines["Geant4Example_SOURCE_DIR"];
    ExampleGeant4Code_build_dir_path  = lines["Geant4Example_BUILD_DIR"];
    ExampleGeant4Code_ExecutionCommand  = lines["Geant4Example_EXECUTION_COMMAND"];
    ui->pushButtonAppBuildDir->setToolTip(ExampleGeant4Code_build_dir_path);
    ui->pushButtonAppSourceDir->setToolTip(ExampleGeant4Code_source_dir_path);
    ui->pushButtonRunExample->setToolTip(ExampleGeant4Code_ExecutionCommand);

    QDir dir = QDir(QCoreApplication::applicationDirPath());
    DoseCalcsCore_build_dir_path = dir.absolutePath()+"/"+DoseCalcs_build_dir_name;
    if(!QFile::exists(DoseCalcsCore_build_dir_path)){ // create install dir for DCMTK
        dir.mkdir(DoseCalcs_build_dir_name);
        showResultsOutput("DoseCalcs build directory " + DoseCalcsCore_build_dir_path + " is created, but you can't run until you build DoseCalcs code", 4);
    }

    CMAKE_Url_String = lines["CMAKE_DOWNLOAD_URL"] ;
    //Xerces_Url_String = lines["XERCES_DOWNLOAD_URL"] ;
    Geant4_Url_String = lines["GEANT4_DOWNLOAD_URL"] ;
    //Supplementary_Url_String = lines["SUPPLEMENTARY_DOWNLOAD_URL"] ;
    MPI_Url_String = lines["MPI_DOWNLOAD_URL"] ;
    Root_Url_String = lines["ROOT_DOWNLOAD_URL"] ;
    DCMTK_Url_String = lines["DCMTK_DOWNLOAD_URL"] ;

    ui->lineEdit_13->setText(CMAKE_Lib_dir_path);
    ui->lineEdit->setText(geant4_Lib_dir_path);
    ui->lineEdit_2->setText(MPI_Lib_dir_path);
    ui->lineEdit_3->setText(Root_Lib_dir_path);
    ui->lineEdit_4->setText(DCMTK_Lib_dir_path);
    ui->lineEdit_5->setText(DoseCalcsCore_source_dir_path);

    if(lines["CENTOS_ROCKS_CLUSTER"] == "YES"){ ui->checkBoxRocks->setChecked(true); }else{ui->checkBoxRocks->setChecked(false);}
    if(lines["MPI_USE"] == "YES"){ ui->checkBox->setChecked(true); }else{ui->checkBox->setChecked(false);}
    if(lines["ROOT_USE"] == "YES"){ ui->checkBox_2->setChecked(true); }else{ui->checkBox_2->setChecked(false);}
    if(lines["DCMTK_USE"] == "YES"){ ui->checkBox_3->setChecked(true); }else{ui->checkBox_3->setChecked(false);}
    if(lines["VERBOSE_USE"] == "YES"){ ui->checkBox_4->setChecked(true); }else{ui->checkBox_4->setChecked(false);}
    if(lines["GDML_USE"] == "YES"){ ui->checkBox_5->setChecked(true); }else{ui->checkBox_5->setChecked(false);}

    CENTOS_ROCKS_CLUSTER  = lines["CENTOS_ROCKS_CLUSTER"];
    MPI_USE = lines["MPI_USE"];
    ROOT_USE = lines["ROOT_USE"];
    DCMTK_USE = lines["DCMTK_USE"];
    VERBOSE_USE = lines["VERBOSE_USE"];
    GDML_USE = lines["GDML_USE"];
    
    if(QFile::exists(CMAKE_Lib_dir_path)){ // create install dir for DCMTK
        cmakeTruePath = "\n " + CMAKE_Lib_dir_path +"/cmake ";
    }
    else{
        cmakeTruePath = "\n cmake " ;
    }

    QString s1 = "../"+GUIPackagesAndFilesDirName;
    if(lines["DEFAULT_DoseCalcs_INPUTS"].contains(s1)){ // /GeometryData/createVolume
        lines["DEFAULT_DoseCalcs_INPUTS"] = lines["DEFAULT_DoseCalcs_INPUTS"].replace(s1,GUIPackagesAndFilesDirPath);
        ui->lineEdit_11->setText(lines["DEFAULT_DoseCalcs_INPUTS"]);
    }
    ui->lineEdit_11->setText(lines["DEFAULT_DoseCalcs_INPUTS"]);
}

// creating shell file for each package to install
void InstallationDialog::create_Xerces_and_PrerequestGeant4Install_sh_file(){
    
    //showResultsOutput("\nWorkSpace Dir : " + GUIPackagesAndFilesDirPath , 0);

    bool ok;
    QString newLabel = QInputDialog::getText(this, "CMAKE Download URL", "CMAKE Download URL", QLineEdit::Normal, CMAKE_Url_String, &ok);
    if (ok){CMAKE_Url_String = newLabel;}

    //CMAKE_dir_install = GUIPackagesAndFilesDirPath +"/CMAKE_install";
    //QDir(GUIPackagesAndFilesDirPath).mkdir("CMAKE_install");

    CMAKE_dir_install = "/usr/local/bin";

    if(OSNameVersion.toLower().contains("ubuntu"))
    {
        Prequest_text_shFile =

                "\n# -------- Prerequisites installation commands \n\n"

                "sudo apt-get install -y build-essential  #Installation of build-essential \n"
                "sudo apt-get install -y libxaw7-dev libxaw7  #Installation of X11 Xmu library and/or headers \n"
                "//sudo apt-get install -y qt5-default   #Installation of Qt User Interface and Visualization \n"
                "sudo apt-get install -y mesa-common-dev ; sudo apt-get install -y libglu1-mesa-dev  #Installation of MesaGL headers and libraries \n"
                "sudo apt-get install -y libxerces-c-dev ; # xerces-c for GDML headers and libraries \n"
                "sudo apt-get install -y libssl-dev \n"
                "sudo apt-get install -y libexpat1-dev ; sudo apt-get install libglu1-mesa-dev -y  #Installation of MesaGL headers and libraries \n"
                "sudo apt-get install -y wget \n"

                "\n# -------- CMAKE installation commands \n\n"

                "cd "+GUIPackagesAndFilesDirPath+"\n"+
                "wget "+ CMAKE_Url_String +"\n"+
                "tar xvf "+ GUIPackagesAndFilesDirPath+"/"+QUrl(CMAKE_Url_String).fileName() + " -C " + GUIPackagesAndFilesDirPath+"\n"+
                "cd "+ QUrl(CMAKE_Url_String).fileName().replace(".tar.gz","").replace(".zip", "") + "\n"+
                "./configure \n"
                //+"-"+"-"+"­prefix=" + CMAKE_dir_install + "\n" + // it give an error when copying --prefix to text file
                "make -j" + QString::number(NumberOfCPUCores) + "# or gmake\n" +
                "sudo make install # or gmake"
                ;

    }
    else if(OSNameVersion.toLower().contains("centos"))
    {
        Prequest_text_shFile =

                "\n# -------- Prerequisites installation commands \n\n"

                "sudo yum install -y build-essential  #Installation of build-essential \n"
                "sudo yum install -y libxaw7-dev libxaw7  #Installation of X11 Xmu library and/or headers \n"
                "//sudo yum install -y qt5-default   #Installation of Qt User Interface and Visualization \n"
                "sudo yum install -y mesa-common-dev ; sudo yum install -y libglu1-mesa-dev  #Installation of MesaGL headers and libraries \n"
                "sudo yum install -y libxerces-c-dev ; \n"
                "sudo yum install -y libexpat1-dev \n"
                "sudo yum install -y wget \n"

                "\n# -------- CMAKE installation commands \n\n"

                "cd "+GUIPackagesAndFilesDirPath+"\n"+
                "wget "+ CMAKE_Url_String +"\n"+
                "tar xvf "+ GUIPackagesAndFilesDirPath+"/"+QUrl(CMAKE_Url_String).fileName() + " -C " + GUIPackagesAndFilesDirPath+"\n"+
                "cd "+ QUrl(CMAKE_Url_String).fileName().replace(".tar.gz","").replace(".zip", "") + "\n"+
                "./configure \n"
                //+"-"+"-"+"­prefix=" + CMAKE_dir_install + "\n" + // it give an error when copying --prefix to text file
                "make -j" + QString::number(NumberOfCPUCores) + "# or gmake\n" +
                "sudo make install # or gmake"
                ;
    }

    showResultsOutput("CMAKE and Prerequest installing Command : " + Prequest_text_shFile , 4); // is the same where we will install it

}
void InstallationDialog::create_Geant4Install_sh_file(){
    
    // we will work now just with the zip file

    bool ok;
    QString newLabel = QInputDialog::getText(this, "GEANT4 Download URL", "GEANT4 Download URL", QLineEdit::Normal, Geant4_Url_String, &ok);
    if (ok){Geant4_Url_String = newLabel;}

    QDir(GUIPackagesAndFilesDirPath).mkdir("Geant4_install");
    geant4_dir_install = GUIPackagesAndFilesDirPath +"/Geant4_install";

    QDir(GUIPackagesAndFilesDirPath).mkdir("Geant4_build");
    geant4_dir_build = GUIPackagesAndFilesDirPath +"/Geant4_build";

    geant4_dir_src = GUIPackagesAndFilesDirPath + "/"+QUrl(Geant4_Url_String).fileName().replace(".tar.gz","").replace(".zip", "");

    // cmake command to install the geant4 libraries
    Cmake_Command = cmakeTruePath+" -DCMAKE_INSTALL_PREFIX=" + geant4_dir_install+
            //" -DXERCESC_ROOT_DIR="+ Xerces_dir_install +
            " -DGEANT4_INSTALL_DATA=ON " +
            " ­‐DGEANT4_BUILD_MULTITHREADED=ON "+
            " -DGEANT4_USE_GDML=ON "+
            " ­-DGEANT4_USE_OPENGL_X11=ON "+
            " -DGEANT4_USE_QT=ON "+
            " ­-DGEANT4_USE_SYSTEM_EXPAT=OFF "+
            geant4_dir_src ;

    Geant4_text_shFile =
            "cd "+GUIPackagesAndFilesDirPath+"\n"+
            "wget "+ Geant4_Url_String +"\n"+
            "tar xvf "+ GUIPackagesAndFilesDirPath+"/"+QUrl(Geant4_Url_String).fileName() + " -C " + GUIPackagesAndFilesDirPath+"\n"+
            "cd "+ geant4_dir_build +"\n"+
            Cmake_Command + "\n" +
            "make -j" + QString::number(NumberOfCPUCores) + "\n" +
            "make install \n"
            ;

    showResultsOutput("Geant4 installing command : " + Geant4_text_shFile , 4); // is the same where we will install it

    /*
    PackageName = "Geant4 source zip";
    Geant4_zip_path = GetChoosenDirFromDialog(0);
    
    if(QFile::exists(Geant4_zip_path)){
        
        QFileInfo fileInfo(Geant4_zip_path);
        geant4_package_name = fileInfo.fileName();
        
        QString geant4_get_install_cmds;
        geant4_get_install_cmds = "tar xvf "+ Geant4_zip_path + " -C " + GUIPackagesAndFilesDirPath;
        showResultsOutput("Unzip command : " + geant4_get_install_cmds , 0);
        execBackgroundProcess(geant4_get_install_cmds);
        
        PackageName = "Geant4 Source "; geant4_dir_src = GetChoosenDirFromDialog(1);
        PackageName = "Xerces install "; Xerces_dir_install = GetChoosenDirFromDialog(1);
        
        showResultsOutput("Xerces_dir_install : " + Xerces_dir_install , 4);
        showResultsOutput("Geant4_dir_src : " + geant4_dir_src , 4); // is the same where we will install it

        geant4_dir_build = GUIPackagesAndFilesDirPath +"/Geant4_build";
        if(!QFile::exists(geant4_dir_build)){ // create install dir for geant4
            QDir* geant4_dir = new QDir(GUIPackagesAndFilesDirPath) ;
            geant4_dir->mkdir("Geant4_build");
        }

        geant4_dir_install = GUIPackagesAndFilesDirPath +"/Geant4_install";
        if(!QFile::exists(geant4_dir_install)){ // create install dir for geant4
            QDir* geant4_dir = new QDir(GUIPackagesAndFilesDirPath) ;
            geant4_dir->mkdir("Geant4_install");
        }
        
        showResultsOutput("geant4_dir_install : " + geant4_dir_install , 4); // is the same where we will install it

        // cmake command to install the geant4 libraries
        Cmake_Command = cmakeTruePath+" -DCMAKE_INSTALL_PREFIX=" + geant4_dir_install+
                " -DXERCESC_ROOT_DIR="+ Xerces_dir_install +
                " -DGEANT4_INSTALL_DATA=ON " +
                " ­‐DGEANT4_BUILD_MULTITHREADED=ON "+
                " -DGEANT4_USE_GDML=ON "+
                " ­-DGEANT4_USE_OPENGL_X11=ON "+
                " -DGEANT4_USE_QT=ON "+
                " ­-DGEANT4_USE_SYSTEM_EXPAT=OFF "+
                geant4_dir_src ;
        
        Geant4_text_shFile =
                "cd "+ geant4_dir_build +"\n"+
                Cmake_Command + "\n" +
                "make -j" + QString::number(NumberOfCPUCores) + "\n" +
                "make install";
        
        showResultsOutput("Geant4 installing command : " + Prequest_text_shFile , 4); // is the same where we will install it
    }
    else{
        showResultsOutput("Choose the Geant4 zip file", 3);
    }
   */
}
void InstallationDialog::create_MPIInstall_sh_file(){
    
    //showResultsOutput("\nWorkSpace Dir : " + GUIPackagesAndFilesDirPath , 0);
    
    // we will work now just with the zip file

    bool ok;
    QString newLabel = QInputDialog::getText(this, "MPI Download URL", "MPI Download URL", QLineEdit::Normal, MPI_Url_String, &ok);
    if (ok){MPI_Url_String = newLabel;}

    QDir(GUIPackagesAndFilesDirPath).mkdir("MPI_install");
    MPI_dir_install = GUIPackagesAndFilesDirPath +"/MPI_install";

    MPI_dir_src = GUIPackagesAndFilesDirPath + "/"+QUrl(MPI_Url_String).fileName().replace(".tar.gz","").replace(".zip", "");

    MPI_text_shFile = "cd "+ MPI_dir_src + "\n"+
            "cd "+GUIPackagesAndFilesDirPath+"\n"+
            "wget "+ MPI_Url_String +"\n"+
            "tar xvf "+ GUIPackagesAndFilesDirPath+"/"+QUrl(MPI_Url_String).fileName() + " -C " + GUIPackagesAndFilesDirPath+"\n"+
            "cd "+ QUrl(MPI_Url_String).fileName().replace(".tar.gz","").replace(".zip", "") + "\n"+
            "./configure --disable-fortran  --with-device=ch4:ofi --­prefix=" + MPI_dir_install + " 2>&1 | tee c.txt \n" +
            "make -j" + QString::number(NumberOfCPUCores) + "\n" +
            "sudo make install";

    showResultsOutput("MPI Installing Command : " + MPI_text_shFile , 0); // is the same where we will install it

  /*

    PackageName = "MPI source zip";
    MPI_zip_path = GetChoosenDirFromDialog(0);
    //ui->MPIDirlineEdit->setText(Xerces_zip_path);
    
    if(QFile::exists(MPI_zip_path)){
        
        QString MPI_get_install_cmds;
        MPI_get_install_cmds = "tar xvf "+ MPI_zip_path + " -C " + GUIPackagesAndFilesDirPath;
        showResultsOutput("Unzip command : " + MPI_get_install_cmds , 0);
        execBackgroundProcess(MPI_get_install_cmds);
        
        //QString MPI_dir_install ;//= xerces_package_name;
        
        PackageName = "MPI Source Dir"; MPI_dir_src = GetChoosenDirFromDialog(1);
        
        if(!QFile::exists(MPI_dir_install) || !MPI_dir_install.isNull() || !MPI_dir_install.isEmpty() ){
            //return;
        }
        
        showResultsOutput("MPI_src_dir : " + MPI_dir_src , 4); // is the same where we will install it

        MPI_dir_install = GUIPackagesAndFilesDirPath +"/MPI_install";
        if(!QFile::exists(MPI_dir_install)){ // create install dir for geant4
            QDir* geant4_dir = new QDir(GUIPackagesAndFilesDirPath) ;
            geant4_dir->mkdir("MPI_install");
        }

        showResultsOutput("MPI_dir_install : " + MPI_dir_install , 4); // is the same where we will install it

        MPI_text_shFile = "cd "+ MPI_dir_src + "\n"+
                "./configure --disable-fortran  --with-device=ch4:ofi --­prefix=" + MPI_dir_install + " 2>&1 | tee c.txt \n" +
                "make -j" + QString::number(NumberOfCPUCores) + "\n" +
                "sudo make install";
    }
    else{
        showResultsOutput("Choose the MPI zip file", 3);
    }
    */
    
}
void InstallationDialog::create_RootInstall_sh_file(){
    
    bool ok;
    QString newLabel = QInputDialog::getText(this, "ROOT Download URL", "ROOT Download URL", QLineEdit::Normal, Root_Url_String, &ok);
    if (ok){Root_Url_String = newLabel;}

    QDir(GUIPackagesAndFilesDirPath).mkdir("ROOT_install");
    ROOT_dir_install = GUIPackagesAndFilesDirPath +"/ROOT_install";

    QDir(GUIPackagesAndFilesDirPath).mkdir("ROOT_build");
    ROOT_dir_build = GUIPackagesAndFilesDirPath +"/ROOT_build";

    ROOT_dir_src = GUIPackagesAndFilesDirPath + "/"+QUrl(Root_Url_String).fileName().replace(".tar.gz","").replace(".zip", "");

    QString nl = "\n";
    Root_text_shFile = "sudo apt-get install ca-certificates curl davix-dev dcap-dev fonts-freefont-ttf g++ gcc gfortran libfcgi-dev libfftw3-dev libfreetype6-dev libftgl-dev libgfal2-dev libgif-dev libgl2ps-dev libglew-dev libglu-dev libgraphviz-dev libgsl-dev libjpeg-dev liblz4-dev liblzma-dev libmysqlclient-dev libpcre++-dev libpng-dev libpq-dev libsqlite3-dev libssl-dev libtbb-dev libtiff-dev libx11-dev libxext-dev libxft-dev libxml2-dev libxpm-dev libxxhash-dev libz-dev libzstd-dev locales python3-dev python3-numpy srm-ifce-dev unixodbc-dev " + nl+
            "sudo apt-get install dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev libxext-dev python libssl-dev gfortran libpcre3-dev xlibmesa-glu-dev libglew1.5-dev libftgl-dev libmysqlclient-dev libfftw3-dev libcfitsio-dev graphviz-dev libavahi-compat-libdnssd-dev libldap2-dev python-dev libxml2-dev libkrb5-dev libgsl0-dev"+ nl +
            "sudo apt-get install git libblas-dev libxft-dev"+ nl +
            "cd "+GUIPackagesAndFilesDirPath+"\n"+
            "wget "+ Root_Url_String +"\n"+
            "tar xvf "+ GUIPackagesAndFilesDirPath+"/"+QUrl(Root_Url_String).fileName() + " -C " + GUIPackagesAndFilesDirPath+"\n"+
            "cd "+ ROOT_dir_build + "\n"+
            cmakeTruePath +" -DCMAKE_INSTALL_PREFIX=" + ROOT_dir_install + " " + ROOT_dir_src + "\n" +
            "make -j" + QString::number(NumberOfCPUCores) + "\n" +
            "make install";

    showResultsOutput("ROOT Installing Command : " + Root_text_shFile , 0); // is the same where we will install it

    /*

    PackageName = "ROOT source zip";
    QString ROOT_zip_path = GetChoosenDirFromDialog(0);
    //ui->ROOTDirlineEdit->setText(Xerces_zip_path);
    
    if(QFile::exists(ROOT_zip_path)){
        
        QString cmds;
        cmds = "tar xvf "+ ROOT_zip_path + " -C " + GUIPackagesAndFilesDirPath;
        showResultsOutput("Unzip command : " + cmds , 0);
        execBackgroundProcess(cmds);
        
        PackageName = "ROOT Source Dir"; ROOT_dir_src = GetChoosenDirFromDialog(1);
        
        if(!QFile::exists(ROOT_dir_src) || !ROOT_dir_src.isNull() || !ROOT_dir_src.isEmpty() ){
            //return;
        }
        
        showResultsOutput("ROOT_dir_src : " + ROOT_dir_src , 4); // is the same where we will install it

        ROOT_dir_build = GUIPackagesAndFilesDirPath +"/ROOT_build";
        if(!QFile::exists(ROOT_dir_build)){ // create install dir for geant4
            QDir* geant4_dir = new QDir(GUIPackagesAndFilesDirPath) ;
            geant4_dir->mkdir("ROOT_build");
        }

        ROOT_dir_install = GUIPackagesAndFilesDirPath +"/ROOT_install";
        if(!QFile::exists(ROOT_dir_install)){ // create install dir for geant4
            QDir* geant4_dir = new QDir(GUIPackagesAndFilesDirPath) ;
            geant4_dir->mkdir("ROOT_install");
        }

        showResultsOutput("ROOT_dir_build : " + ROOT_dir_build , 4); // is the same where we will install it
        showResultsOutput("ROOT_dir_install : " + ROOT_dir_install , 4); // is the same where we will install it

        QString nl = "\n";
        Root_text_shFile = "#sudo apt-get install ca-certificates curl davix-dev dcap-dev fonts-freefont-ttf g++ gcc gfortran libfcgi-dev libfftw3-dev libfreetype6-dev libftgl-dev libgfal2-dev libgif-dev libgl2ps-dev libglew-dev libglu-dev libgraphviz-dev libgsl-dev libjpeg-dev liblz4-dev liblzma-dev libmysqlclient-dev libpcre++-dev libpng-dev libpq-dev libsqlite3-dev libssl-dev libtbb-dev libtiff-dev libx11-dev libxext-dev libxft-dev libxml2-dev libxpm-dev libxxhash-dev libz-dev libzstd-dev locales python3-dev python3-numpy srm-ifce-dev unixodbc-dev " + nl+
                "# sudo apt-get install dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev python libssl-dev gfortran libpcre3-dev xlibmesa-glu-dev libglew1.5-dev libftgl-dev libmysqlclient-dev libfftw3-dev libcfitsio-dev graphviz-dev libavahi-compat-libdnssd-dev libldap2-dev python-dev libxml2-dev libkrb5-dev libgsl0-dev "+ nl +
                "cd "+ ROOT_dir_build + "\n"+
                cmakeTruePath +" -DCMAKE_INSTALL_PREFIX=" + ROOT_dir_install + " " + ROOT_dir_src + "\n" +
                "make -j" + QString::number(NumberOfCPUCores) + "\n" +
                "make install";
        
        filesManagerObj->WriteTextToFile(Root_sh_path,Root_text_shFile);
        
        if(QFile::exists(Root_sh_path)){
            ui->ROOTDirlineEdit->setText(Root_sh_path);
            showResultsOutput("the installation file : " + Root_sh_path + " is created." , 0);
        }else{
            ui->ROOTDirlineEdit->setText("");
            showResultsOutput("the installation file not been created.", 3);
        }
        
    }
    else{
        showResultsOutput("Choose the ROOT zip file", 3);
    }
    
    */
}
void InstallationDialog::create_RootInstall_Dir_file(){
    
    PackageName = "ROOT source zip";
    QString ROOT_zip_path = GetChoosenDirFromDialog(0);
    //ui->MPIDirlineEdit->setText(Xerces_zip_path);
    
    if(QFile::exists(ROOT_zip_path)){
        
        QString cmds;
        cmds = "tar xvf "+ ROOT_zip_path + " -C " + GUIPackagesAndFilesDirPath;
        showResultsOutput("Unzip command : " + cmds , 0);
        execBackgroundProcess(cmds);
        
        PackageName = "ROOT Source Dir"; ROOT_dir_src = GetChoosenDirFromDialog(1);
        
        showResultsOutput("ROOT_dir_src : " + ROOT_dir_src , 4); // is the same where we will install it
        
    }
    else{
        showResultsOutput("Choose the ROOT zip file", 3);
    }
    
}
void InstallationDialog::create_DCMTK_Install_sh_file(){
    
    bool ok;
    QString newLabel = QInputDialog::getText(this, "DCMTK Download URL", "DCMTK Download URL", QLineEdit::Normal, DCMTK_Url_String, &ok);
    if (ok){DCMTK_Url_String = newLabel;}

    QDir(GUIPackagesAndFilesDirPath).mkdir("DCMTK_install");
    DCMTK_dir_install = GUIPackagesAndFilesDirPath +"/DCMTK_install";

    QDir(GUIPackagesAndFilesDirPath).mkdir("DCMTK_build");
    DCMTK_dir_build = GUIPackagesAndFilesDirPath +"/DCMTK_build";

    DCMTK_dir_src = GUIPackagesAndFilesDirPath + "/"+QUrl(DCMTK_Url_String).fileName().replace(".tar.gz","").replace(".zip", "");

    Cmake_Command = cmakeTruePath +" " + DCMTK_dir_src ;

    DCMTK_text_shFile =
            "cd "+GUIPackagesAndFilesDirPath+"\n"+
            "wget "+ DCMTK_Url_String +"\n"+
            "tar xvf "+ GUIPackagesAndFilesDirPath+"/"+QUrl(DCMTK_Url_String).fileName() + " -C " + GUIPackagesAndFilesDirPath+"\n"+
            "cd "+ DCMTK_dir_build +"\n"+
            Cmake_Command + "\n" +
            "make -j" + QString::number(NumberOfCPUCores) + "\n" +
            "make -DESTDIR=" + DCMTK_dir_install+ " install";

    showResultsOutput("DCMTK Installing Command : " + DCMTK_text_shFile , 0); // is the same where we will install it

    // we will work now just with the zip file

/*
    PackageName = "DCMTK source zip";
    DCMTK_zip_path = GetChoosenDirFromDialog(0);
    
    if(QFile::exists(DCMTK_zip_path)){
        
        QFileInfo fileInfo(DCMTK_zip_path);
        DCMTK_package_name = fileInfo.fileName();
        
        QString DCMTK_get_install_cmds;
        DCMTK_get_install_cmds = "tar xvf "+ DCMTK_zip_path + " -C " + GUIPackagesAndFilesDirPath;
        showResultsOutput("Unzip command : " + DCMTK_get_install_cmds , 0);
        execBackgroundProcess(DCMTK_get_install_cmds);
        
        PackageName = "DCMTK Source Dir"; DCMTK_dir_src = GetChoosenDirFromDialog(1);
        
        if(!QFile::exists(DCMTK_dir_src) || !DCMTK_dir_src.isNull() || !DCMTK_dir_src.isEmpty() ){
            //return;
        }
        

        showResultsOutput("DCMTK_dir_src : " + DCMTK_dir_src , 4); // is the same where we will install it


        DCMTK_dir_build = GUIPackagesAndFilesDirPath +"/DCMTK_build";
        if(!QFile::exists(DCMTK_dir_build)){ // create install dir for geant4
            QDir* geant4_dir = new QDir(GUIPackagesAndFilesDirPath) ;
            geant4_dir->mkdir("DCMTK_build");
        }

        DCMTK_dir_install = GUIPackagesAndFilesDirPath +"/DCMTK_install";
        if(!QFile::exists(DCMTK_dir_install)){ // create install dir for geant4
            QDir* geant4_dir = new QDir(GUIPackagesAndFilesDirPath) ;
            geant4_dir->mkdir("DCMTK_install");
        }

        showResultsOutput("DCMTK_dir_install : " + DCMTK_dir_install , 4); // is the same where we will install it

        Cmake_Command = cmakeTruePath +" " + DCMTK_dir_src ;
        
        DCMTK_text_shFile =
                "cd "+ DCMTK_dir_build +"\n"+
                Cmake_Command + "\n" +
                "make -j" + QString::number(NumberOfCPUCores) + "\n" +
                "make -DESTDIR=" + DCMTK_dir_install+ " install";
        
        showResultsOutput("DCMTK Installing Command : " + Prequest_text_shFile , 0); // is the same where we will install it
    }
    else{
        showResultsOutput("Choose the DCMTK zip file", 3);
    }
    */
}
void InstallationDialog::create_AppInstall_sh_file(){
    
    getPathsFromEditLines();
    
    if(!QFile::exists(DoseCalcsCore_source_dir_path)){ // create install dir for DCMTK
        PackageName = "DoseCalcs core source";
        DoseCalcsCore_source_dir_path = GetChoosenDirFromDialog(1);
    }

    /*QDir* dir = new QDir(DoseCalcsCore_source_dir_path) ; dir->cdUp();
    DoseCalcsCore_build_dir_path = dir->absolutePath()+"/"+DoseCalcs_build_dir_name;
    if(!QFile::exists(DoseCalcsCore_build_dir_path)){ // create install dir for DCMTK
        QDir* dir = new QDir(DoseCalcsCore_source_dir_path) ; dir->cdUp();
        dir->mkdir(DoseCalcs_build_dir_name);
        showResultsOutput("DoseCalcs build directory " + DoseCalcsCore_build_dir_path + " is created, but you can't run until you build DoseCalcs code", 4);
    }
    */

    QDir dir = QDir(QCoreApplication::applicationDirPath());
    DoseCalcsCore_build_dir_path = dir.absolutePath()+"/"+DoseCalcs_build_dir_name;
    if(!QFile::exists(DoseCalcsCore_build_dir_path)){ // create install dir for DCMTK
        dir.mkdir(DoseCalcs_build_dir_name);
        showResultsOutput("DoseCalcs build directory " + DoseCalcsCore_build_dir_path + " is created, but you can't run until you build DoseCalcs code", 4);
    }

    QString PackagesToAdd = "";
    QString PathToAdd = "";
    if(ui->checkBox->isChecked()){
        PackagesToAdd += " -DWITH_G4MPI_USE=ON -DCMAKE_CXX_COMPILER=" + MPI_Lib_dir_path+"/mpicxx -DCMAKE_C_COMPILER=" + MPI_Lib_dir_path+"/mpicc ";
        //+ " -DMPI_DIR=" + MPI_Lib_dir_path
    }
    else{
        PackagesToAdd += " -DWITH_G4MPI_USE=OFF ";
    }
    
    if(ui->checkBox_2->isChecked()){
        PackagesToAdd += " -DWITH_ANALYSIS_USE=ON -DROOT_DIR="+ Root_Lib_dir_path + " ";
        PathToAdd +=  "\ncd " +Root_Lib_dir_path +"\n" + ". ./thisroot.sh \n";
    }
    else{
        PackagesToAdd += " -DWITH_ANALYSIS_USE=OFF ";
    }
    
    if(ui->checkBox_3->isChecked()){
        PackagesToAdd += " -DWITH_DCMTK_USE=ON -DDCMTK_DIR="+DCMTK_Lib_dir_path + " ";
    }
    else{
        PackagesToAdd += " -DWITH_DCMTK_USE=OFF ";
    }
    
    if(ui->checkBox_4->isChecked()){
        PackagesToAdd += " -DWITH_VERBOSE_USE=ON ";
    }
    else{
        PackagesToAdd += " -DWITH_VERBOSE_USE=OFF ";
    }

    if(ui->checkBox_5->isChecked()){
        PackagesToAdd += " -DWITH_GDML_USE=ON ";
    }
    else{
        PackagesToAdd += " -DWITH_GDML_USE=OFF ";
    }


    Cmake_Command = cmakeTruePath + " " + PackagesToAdd
            + DoseCalcsCore_source_dir_path ;
    
    App_text_shFile =
            PathToAdd +
            "\ncd " +geant4_Lib_dir_path +
            "\n. ./geant4.sh\n" +
            "cd "+ DoseCalcsCore_build_dir_path +
            "\nmake clean \nrm -r CMake* cmake_install.cmake Makefile simulate merge analysis \n"+
            Cmake_Command + "\n"+
            "make -j" + QString::number(NumberOfCPUCores) +"\n"+
            ""+"\n";
}
void InstallationDialog::create_Geant4ExampleInstall_sh_file(){

    //getPathsFromEditLines();

    if(!QFile::exists(ExampleGeant4Code_source_dir_path)){ // create install dir for DCMTK
        PackageName = "Example Geant4 Code Source";
        ExampleGeant4Code_source_dir_path = GetChoosenDirFromDialog(1);
    }

    if(!QFile::exists(ExampleGeant4Code_build_dir_path)){ // create install dir for DCMTK

        QDir(QCoreApplication::applicationDirPath()).mkdir(ExampleGeant4Code_build_dir_path);
        showResultsOutput("Geant4 Example build directory " + ExampleGeant4Code_build_dir_path + " is created, but you can't run until you build the Geant4 example code", 4);
    }

    QString PackagesToAdd = "";
    QString PathToAdd = "";
    if(ui->checkBox->isChecked()){
        PackagesToAdd += " -DWITH_G4MPI_USE=ON -DCMAKE_CXX_COMPILER=" + MPI_Lib_dir_path+"/mpicxx -DCMAKE_C_COMPILER=" + MPI_Lib_dir_path+"/mpicc ";
        //+ " -DMPI_DIR=" + MPI_Lib_dir_path
    }
    else{
        PackagesToAdd += " -DWITH_G4MPI_USE=OFF ";
    }

    if(ui->checkBox_2->isChecked()){
        PackagesToAdd += " -DWITH_ANALYSIS_USE=ON -DROOT_DIR="+ Root_Lib_dir_path + " ";
        PathToAdd +=  "\ncd " +Root_Lib_dir_path +"\n" + ". ./thisroot.sh \n";
    }
    else{
        PackagesToAdd += " -DWITH_ANALYSIS_USE=OFF ";
    }

    if(ui->checkBox_3->isChecked()){
        PackagesToAdd += " -DWITH_DCMTK_USE=ON -DDCMTK_DIR="+DCMTK_Lib_dir_path + " ";
    }
    else{
        PackagesToAdd += " -DWITH_DCMTK_USE=OFF ";
    }

    if(ui->checkBox_4->isChecked()){
        PackagesToAdd += " -DWITH_VERBOSE_USE=ON ";
    }
    else{
        PackagesToAdd += " -DWITH_VERBOSE_USE=OFF ";
    }

    if(ui->checkBox_5->isChecked()){
        PackagesToAdd += " -DWITH_GDML_USE=ON ";
    }
    else{
        PackagesToAdd += " -DWITH_GDML_USE=OFF ";
    }


    Cmake_Command = cmakeTruePath + " "
            //+ PackagesToAdd
            + ExampleGeant4Code_source_dir_path ;

    Geant4Example_text_shFile =
            PathToAdd +
            "\ncd " +geant4_Lib_dir_path +
            "\n. ./geant4.sh\n" +
            "cd "+ ExampleGeant4Code_build_dir_path +
            "\n#make clean \n#rm -r * \n"+
            Cmake_Command + "\n"+
            "make -j" + QString::number(NumberOfCPUCores) +"\n"+
            ""+"\n";
}

// for generate/install buttons
void InstallationDialog::on_installXercesButton_clicked()
{
    
    if(!QFile::exists(Prequests_sh_path)){
        create_Xerces_and_PrerequestGeant4Install_sh_file();
        filesManagerObj->WriteTextToFile(Prequests_sh_path,Prequest_text_shFile);
    }
    
    if(QFile::exists(Prequests_sh_path)){
        ui->XercesDirlineEdit->setText(Prequests_sh_path);
        showResultsOutput("Installation file : " + Prequests_sh_path + " is created." , 0);
        
        execShProcess( Prequests_sh_path);
        
        CMAKE_Lib_dir_path = CMAKE_dir_install;
        cmakeTruePath = CMAKE_Lib_dir_path+"/cmake";
        ui->lineEdit_13->setText(CMAKE_Lib_dir_path);

        //showResultsOutput("Installation terminated", 0);
        
    }else{
        ui->XercesDirlineEdit->setText("");
        showResultsOutput("Installation file not found", 3);
    }
}
void InstallationDialog::on_InstallGeant4Button_clicked()
{
    
    if(!QFile::exists(Geant4_sh_path)){
        create_Geant4Install_sh_file();
        filesManagerObj->WriteTextToFile( Geant4_sh_path , Geant4_text_shFile );
    }
    
    if(QFile::exists(Geant4_sh_path)){
        ui->Geant4DirlineEdit->setText(Geant4_sh_path);
        showResultsOutput("Installation file : " + Geant4_sh_path + " is created." , 0);
        execShProcess( Geant4_sh_path);
        
        geant4_Lib_dir_path = geant4_dir_install+"/bin";
        ui->lineEdit->setText(geant4_Lib_dir_path);
        
    }else{
        ui->Geant4DirlineEdit->setText("");
        showResultsOutput("Installation file not found", 3);
    }
}
void InstallationDialog::on_InstallMPIButton_clicked()
{
    if(!QFile::exists(MPI_sh_path)){
        create_MPIInstall_sh_file();
        filesManagerObj->WriteTextToFile(MPI_sh_path , MPI_text_shFile);
    }
    
    if(QFile::exists(MPI_sh_path)){
        ui->MPIDirlineEdit->setText(MPI_sh_path);
        showResultsOutput("Installation file : " + MPI_sh_path + " is created." , 4);
        execShProcess( MPI_sh_path);
        
        // we dont need the library , we need just the mpicc and mpicxx binaries, and this can point to allk libraries and files needed by DoseCalcs
        MPI_Lib_dir_path = MPI_dir_install+"/bin";
        ui->lineEdit_2->setText(MPI_Lib_dir_path);
    }else{
        ui->MPIDirlineEdit->setText("");
        showResultsOutput("Installation file not found", 3);
    }
    
}
void InstallationDialog::on_InstallROOTButton_clicked()
{

    if(ui->checkBoxUseCompiledSource->isChecked()){
        create_RootInstall_Dir_file();
    }else{
        
        if(!QFile::exists(Root_sh_path)){
            create_RootInstall_sh_file();
            filesManagerObj->WriteTextToFile(Root_sh_path , Root_text_shFile);
        }

        if(QFile::exists(Root_sh_path)){
            ui->ROOTDirlineEdit->setText(Root_sh_path);
            showResultsOutput("Installation file : " + Root_sh_path + " is created." , 4);
            execShProcess( Root_sh_path);

            // we dont need the library , we need just the mpicc and mpicxx binaries, and this can point to allk libraries and files needed by DoseCalcs
            Root_Lib_dir_path = ROOT_dir_install+"/bin";
            ui->lineEdit_3->setText(Root_Lib_dir_path);

        }else{
            ui->ROOTDirlineEdit->setText("");
            showResultsOutput("Installation file not found", 3);
        }
    }    
    
}
void InstallationDialog::on_InstallDCMTKButton_clicked()
{
    
    if(!QFile::exists(DCMTK_sh_path)){
        create_DCMTK_Install_sh_file();
        filesManagerObj->WriteTextToFile( DCMTK_sh_path , DCMTK_text_shFile );
    }
    
    if(QFile::exists(DCMTK_sh_path)){
        ui->DCMTKDirlineEdit->setText(DCMTK_sh_path);
        showResultsOutput("Installation file : " + DCMTK_sh_path + " is created." , 4);
        execShProcess( DCMTK_sh_path);
        
        DCMTK_Lib_dir_path = DCMTK_dir_install+"/usr/local/lib/cmake/dcmtk";
        ui->lineEdit_4->setText(DCMTK_Lib_dir_path);

    }else{
        ui->MPIDirlineEdit->setText("");
        showResultsOutput("Installation file not found", 3);
    }
    
}
void InstallationDialog::on_installCodeGeant4Button_clicked()
{
    if(!QFile::exists(App_sh_path)){
        create_AppInstall_sh_file();
        filesManagerObj->WriteTextToFile(App_sh_path,App_text_shFile);
    }


    if(!QFile::exists(Root_Lib_dir_path)){ // create install dir for geant4
        //showResultsOutput("The library : " + Root_Lib_dir_path + " is not found. Please select Root library path in configuration variables" , 3);
    }else{
        //showResultsOutput("The library : " + Root_Lib_dir_path + " is found. DoseCalcs now can use Root library" , 4);
    }
    ui->lineEdit_3->setText(Root_Lib_dir_path);


    if(QFile::exists(App_sh_path)){
        //ui->lineEdit_5->setText(App_sh_path);
        showResultsOutput("Installation file : " + App_sh_path + " is created." , 4);
        
        execShProcess( App_sh_path);

    }else{
        ui->lineEdit_5->setText("");
        showResultsOutput("Installation file not found", 3);
    }
    
}

void InstallationDialog::on_pushButtonUpdatePaths_clicked()
{
    QMessageBox::information(this, tr(""), TestPackagesPathsBeforToRun());
}


// for open/install buttons
void InstallationDialog::on_openAndInstallPrerequesiteBtn_clicked()
{
    EditFlag = 1;
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(Prequests_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"CMAKE and Prerequisites installer");

    //PackageName = "Xerces Installer.sh file"; Prequests_sh_path = GetChoosenDirFromDialog(0);
    //GenerateShellForInstall = false;
    //on_installXercesButton_clicked(); // should not blocked here, it will make GenerateShellForInstall = false; for all timesn then you can generate/install until oyu reopen the installation dialog
    //GenerateShellForInstall = true;
}
void InstallationDialog::on_openAndInstallGeant4Btn_clicked()
{
    
    EditFlag = 3;
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(Geant4_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"GEANT4 installer");

    //PackageName = "Geant4 Installer.sh file"; Geant4_sh_path = GetChoosenDirFromDialog(0);
    //GenerateShellForInstall = false;
    //on_InstallGeant4Button_clicked();// should not blocked here, it will make GenerateShellForInstall = false; for all timesn then you can generate/install until oyu reopen the installation dialog
    //GenerateShellForInstall = true;
}
void InstallationDialog::on_openAndInstallMPIBtn_clicked()
{

    EditFlag = 4;
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(MPI_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"MPI installer");

    //PackageName = "MPI Installer.sh file"; MPI_sh_path = GetChoosenDirFromDialog(0);
    //GenerateShellForInstall = false;
    //on_InstallMPIButton_clicked();// should not blocked here, it will make GenerateShellForInstall = false; for all timesn then you can generate/install until oyu reopen the installation dialog
    //GenerateShellForInstall = true;
}
void InstallationDialog::on_openAndInstallROOTBtn_clicked()
{

    EditFlag = 5;
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(Root_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"ROOT installer");

    //PackageName = "ROOT Installer.sh file"; Root_sh_path = GetChoosenDirFromDialog(0);
    //GenerateShellForInstall = false;
    //on_InstallROOTButton_clicked(); // should not blocked here, it will make GenerateShellForInstall = false; for all timesn then you can generate/install until oyu reopen the installation dialog
    //GenerateShellForInstall = true;
}
void InstallationDialog::on_openAndInstallDCMTKBtn_clicked()
{

    EditFlag = 6;
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(DCMTK_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"DCMTK installer");

    //PackageName = "DCMTK Installer.sh file"; DCMTK_sh_path = GetChoosenDirFromDialog(0);
    //GenerateShellForInstall = false;
    //on_InstallDCMTKButton_clicked(); // should not blocked here, it will make GenerateShellForInstall = false; for all timesn then you can generate/install until oyu reopen the installation dialog
    //GenerateShellForInstall = true;
}
void InstallationDialog::on_openAndInstallAPPBtn_clicked()
{

    EditFlag = 7;
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(App_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"DoseCalcs Builder");

    //PackageName = "DoseCalcs Installer .sh file"; App_sh_path = GetChoosenDirFromDialog(0);
    //GenerateShellForInstall = false;
    //on_installCodeGeant4Button_clicked(); // should not blocked here, it will make GenerateShellForInstall = false; for all timesn then you can generate/install until oyu reopen the installation dialog
    //GenerateShellForInstall = true;
}

// for download buttons
void InstallationDialog::DownloadFileToDirectory(QString Dow_dir, QString URL)
{


#if QT_NETWORK_LIB
    HttpDownloadObj = new HttpDownload(this);
    HttpDownloadObj->setDownloadDir(Dow_dir);
    HttpDownloadObj->setUrl(URL);
    HttpDownloadObj->show();
#endif

}
void InstallationDialog::on_DownloadGeant4Button_clicked()
{
    DownloadFileToDirectory(GUIPackagesAndFilesDirPath, Geant4_Url_String);
}
void InstallationDialog::on_DownloadMPIButton_clicked()
{
    DownloadFileToDirectory(GUIPackagesAndFilesDirPath, MPI_Url_String);
}
void InstallationDialog::on_DownloadROOTButton_clicked()
{
    DownloadFileToDirectory(GUIPackagesAndFilesDirPath, Root_Url_String);
}
void InstallationDialog::on_DownloadDCMTKButton_clicked()
{
    DownloadFileToDirectory(GUIPackagesAndFilesDirPath, DCMTK_Url_String);
}

// called from open buttons and install buttons
QString InstallationDialog::GetChoosenDirFromDialog(int type){
    
    if(type == 0){
        PackageName = "Choose the "+ PackageName;
    }
    else if(type == 1){
        PackageName = "Choose the "+ PackageName + " directory";
    }
    
    if(type == 0){ // to select File
        chosen_Dir = QFileDialog::getOpenFileName(
                    this,
                    tr(PackageName.toStdString().c_str()),
                    GUIPackagesAndFilesDirPath,
                    "All files (*.*)"//;;Text files (*.txt)" // extentions to show
                    );
    }
    else if(type == 1){ // to select Directory
        chosen_Dir = QFileDialog::getExistingDirectory(0, (PackageName), GUIPackagesAndFilesDirPath);
    }
    else{
        chosen_Dir = QFileDialog::getExistingDirectory(0, (PackageName), GUIPackagesAndFilesDirPath);
    }
    
    //QMessageBox::information(this, tr("File Path is : "), defaultFileName);
    
    if(chosen_Dir != NULL && chosen_Dir != "" ){
        
        showResultsOutput("The chosen file is " + chosen_Dir , 0);
        
    }else{
        //QMessageBox::information(this, tr("File..!!"),"No file is choosen");
        showResultsOutput("No file or directory is chosen ", 3);
        
        chosen_Dir = "";
    }
    
    return chosen_Dir;
    
}

void InstallationDialog::showResultsOutput(QString text , int level){
    
    if(level == 0){
        QTextStream(stdout) << text << "\n";
        return;
    }
    if(level == 1){
        QTextStream(stdout) << " ---------------------> " << text << "\n";
        ui->installationDialogTextOutput->appendPlainText(" ---------------------> "+text);
    }
    if(level == 4){
        QTextStream(stdout) << text << "\n";
        ui->installationDialogTextOutput->appendPlainText("\n"+text);
    }
    if(level == 3){
        QTextStream(stdout) << "!!!!!!!!!!!!!!!!! " << text << "\n";
        ui->installationDialogTextOutput->appendPlainText("!!!!!!!!!!!!!!!!! "+text);
    }
}

// called from on_install... SLOT to execute installation shell file for a package
void InstallationDialog::execShProcess(QString command){
    

    filesManagerObj->WriteTextToFile(command, filesManagerObj->ReadTextFromFileInOneString(command)+"\n bash");

    //command += "\n bash \n";

    ui->tabWidget->removeTab(0);
    terminal *term = new terminal;
    term->executeCommand(command);
    ui->tabWidget->insertTab(0, term, "Terminal");
    ui->tabWidget->setCurrentIndex(0);
    
    //QProcess process;
    //QStringList ff=(QStringList() << "-hold" << "-e" << "sh" << command);
    //process.startDetached("xterm", ff);
    
    /* it works tww but it block mother process
    process = new QProcess();
    //process->start(program, params);
    process->execute("xterm -hold -e sh " + command);
    process->close();
    */
}
void InstallationDialog::execBackgroundProcess(QString command){
    
    QProcess process;
    process.execute(command);
    //QStringList ff=(QStringList() << "-hold" << "-e" << "sh" << command);
    //process.start("xterm", ff);
}

void InstallationDialog::on_pushButtonGenerateInstallerFoAll_clicked()
{
    if(ui->checkBoxInstallXerces->isChecked()){
        create_Xerces_and_PrerequestGeant4Install_sh_file();
        InstallAll_text_shFile += " \n ############################## Prerequisites ##############################\n" +
                Prequest_text_shFile;
    }
    if(ui->checkBoxInstallGeant4->isChecked()){
        create_Geant4Install_sh_file();
        InstallAll_text_shFile += " \n ############################## GEANT4 ##############################\n" +
                Geant4_text_shFile;
    }
    if(ui->checkBoxInstallMPI->isChecked()){
        create_MPIInstall_sh_file();
        InstallAll_text_shFile += " \n ############################## MPI ##############################\n" +
                MPI_text_shFile;
    }
    if(ui->checkBoxInstallROOT->isChecked()){
        create_RootInstall_sh_file();
        InstallAll_text_shFile += " \n ############################## ROOT ##############################\n" +
                Root_text_shFile;
    }
    if(ui->checkBoxInstallDCMTK->isChecked()){
        create_DCMTK_Install_sh_file();
        InstallAll_text_shFile += " \n ############################## DCMTK ##############################\n" +
                DCMTK_text_shFile;
    }
    
    
    filesManagerObj->WriteTextToFile(AllPackagesInstall_sh_path, InstallAll_text_shFile);
    EditFlag = 8;
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(AllPackagesInstall_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"All packages installer");
}
void InstallationDialog::on_pushButtonInstallAllPrerequisites_clicked()
{
    if(!QFile::exists(AllPackagesInstall_sh_path)){
        create_Xerces_and_PrerequestGeant4Install_sh_file();
        filesManagerObj->WriteTextToFile(AllPackagesInstall_sh_path,Prequest_text_shFile);
    }

    if(QFile::exists(AllPackagesInstall_sh_path)){
        //ui->XercesDirlineEdit->setText(AllPackagesInstall_sh_path);
        showResultsOutput("Installation file : " + AllPackagesInstall_sh_path + " is created." , 0);

        //QProcess* proc = new QProcess();
        //proc->execute("xterm -hold -e sh " + Prequests_sh_path);
        //proc->waitForStarted();
        //proc->waitForFinished();

        execShProcess( AllPackagesInstall_sh_path);

        //showResultsOutput("Installation terminated", 0);

    }else{
        //ui->XercesDirlineEdit->setText("");
        showResultsOutput("Installation file not found", 3);
    }

}
void InstallationDialog::on_pushButtonLoadForAll_clicked()
{
    EditFlag = 8;
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(AllPackagesInstall_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1," All packages installer");
}

// Load, edit and save configuration data to Config file
void InstallationDialog::on_ReadPathFilepushButton_clicked()
{
    getConfigurationDataForInstallation();

}
void InstallationDialog::on_EditPathsFilepushButton_clicked()
{
    EditFlag = 9;

    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(GUIPackagesAndFilesDirPath+"/"+GUIConfigFileName));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,GUIConfigFileName);

    //QString command = GUIPackagesAndFilesDirPath+"/"+GUIConfigFileName;
    //QProcess process;
    //QStringList qsl = {command};
    //process.startDetached("gedit", qsl);
}
void InstallationDialog::on_pushButton_7_clicked()
{
    

    if(ui->checkBoxRocks->isChecked()){ CENTOS_ROCKS_CLUSTER = "YES"; }else{CENTOS_ROCKS_CLUSTER = "NO";}
    if(ui->checkBox->isChecked()){ MPI_USE = "YES"; }else{MPI_USE = "NO";}
    if(ui->checkBox_2->isChecked()){ ROOT_USE = "YES"; }else{ROOT_USE = "NO";}
    if(ui->checkBox_3->isChecked()){ DCMTK_USE = "YES"; }else{DCMTK_USE = "NO";}
    if(ui->checkBox_4->isChecked()){ VERBOSE_USE = "YES"; }else{VERBOSE_USE = "NO";}
    if(ui->checkBox_5->isChecked()){ GDML_USE = "YES"; }else{GDML_USE = "NO";}

    getPathsFromEditLines();
    
    ConfigDataText =
            "\nCMAKE_INSTALL_DIR      "+ CMAKE_Lib_dir_path +"\n"+
            "GEANT4_INSTALL_DIR       "+ geant4_Lib_dir_path +"\n"+
            "MPI_INSTALL_DIR          "+ MPI_Lib_dir_path +"\n"+
            "ROOT_INSTALL_DIR         "+ Root_Lib_dir_path +"\n"+
            "DCMTK_INSTALL_DIR        "+ DCMTK_Lib_dir_path +"\n"+
            "DoseCalcs_SOURCE_DIR     "+ DoseCalcsCore_source_dir_path +"\n\n"+
            "Geant4Example_SOURCE_DIR "+ ExampleGeant4Code_source_dir_path+"\n"+
            "Geant4Example_BUILD_DIR  "+ ExampleGeant4Code_build_dir_path+"\n"+
            "Geant4Example_EXECUTION_COMMAND  "+ ExampleGeant4Code_ExecutionCommand+"\n\n"+

            "DEFAULT_DoseCalcs_INPUTS     "+ ui->lineEdit_11->text() +"\n\n"

            "CMAKE_DOWNLOAD_URL     "+ CMAKE_Url_String +"\n"+
            "GEANT4_DOWNLOAD_URL    "+ Geant4_Url_String +"\n"+
            //"XERCES_DOWNLOAD_URL    "+ Xerces_Url_String +"\n"+
            "MPI_DOWNLOAD_URL       "+ MPI_Url_String +"\n"+
            "ROOT_DOWNLOAD_URL      "+ Root_Url_String +"\n"+
            "DCMTK_DOWNLOAD_URL     "+ DCMTK_Url_String +"\n\n"+
            //"SUPPLEMENTARY_DOWNLOAD_URL    "+ Supplementary_Url_String +"\n"+


            "VERBOSE_USE            "+ VERBOSE_USE +"\n"+
            "GDML_USE               "+ GDML_USE +"\n"+
            "CENTOS_ROCKS_CLUSTER   "+ CENTOS_ROCKS_CLUSTER +"\n"+
            "MPI_USE                "+ MPI_USE +"\n"+
            "ROOT_USE               "+ ROOT_USE +"\n"+
            "DCMTK_USE              "+ DCMTK_USE +"\n\n"
            ;
    
    filesManagerObj->WriteTextToFile( GUIPackagesAndFilesDirPath+"/"+GUIConfigFileName , ConfigDataText);
    
}

void InstallationDialog::getPathsFromEditLines(){
    
    if(QFile::exists(ui->lineEdit_13->text())){CMAKE_Lib_dir_path = ui->lineEdit_13->text();}
    if(QFile::exists(CMAKE_Lib_dir_path)){ // create install dir for DCMTK
        cmakeTruePath = "\n " + CMAKE_Lib_dir_path +"/cmake ";
    }
    else{
        cmakeTruePath = "\n cmake " ;
    }
    if(QFile::exists(ui->lineEdit->text())){geant4_Lib_dir_path = ui->lineEdit->text();}
    if(QFile::exists(ui->lineEdit_2->text())){MPI_Lib_dir_path = ui->lineEdit_2->text();}
    if(QFile::exists(ui->lineEdit_3->text())){Root_Lib_dir_path = ui->lineEdit_3->text();}
    if(QFile::exists(ui->lineEdit_4->text())){DCMTK_Lib_dir_path = ui->lineEdit_4->text();}
    if(QFile::exists(ui->lineEdit_5->text())){

        DoseCalcsCore_source_dir_path = ui->lineEdit_5->text();

        /*QDir* dir = new QDir(DoseCalcsCore_source_dir_path) ; dir->cdUp();
        DoseCalcsCore_build_dir_path = dir->absolutePath()+"/"+DoseCalcs_build_dir_name;
        if(!QFile::exists(DoseCalcsCore_build_dir_path)){ // create install dir for DCMTK
            QDir* dir = new QDir(DoseCalcsCore_source_dir_path) ; dir->cdUp();
            dir->mkdir(DoseCalcs_build_dir_name);
            showResultsOutput("DoseCalcs build directory " + DoseCalcsCore_build_dir_path + " is created, but you can't run until you build DoseCalcs code", 4);
        }
        */

        QDir dir = QDir(QCoreApplication::applicationDirPath());
        DoseCalcsCore_build_dir_path = dir.absolutePath()+"/"+DoseCalcs_build_dir_name;
        if(!QFile::exists(DoseCalcsCore_build_dir_path)){ // create install dir for DCMTK
            dir.mkdir(DoseCalcs_build_dir_name);
            showResultsOutput("DoseCalcs build directory " + DoseCalcsCore_build_dir_path + " is created, but you can't run until you build DoseCalcs code", 4);
        }

    }
}

// Open install Directories
void InstallationDialog::on_pushButton_clicked()
{
    PackageName = "Geant4 Install bin ";
    QString a = GetChoosenDirFromDialog(1);
    if(a != ""){
        geant4_Lib_dir_path = a;
        ui->lineEdit->setText(geant4_Lib_dir_path);
    }

}
void InstallationDialog::on_pushButton_8_clicked()
{
    PackageName = "CMAKE Install ";
    QString a = GetChoosenDirFromDialog(1);
    if(a != ""){
        CMAKE_Lib_dir_path = a;
        ui->lineEdit_13->setText(CMAKE_Lib_dir_path);
    }

    if(QFile::exists(CMAKE_Lib_dir_path)){ // create install dir for CMAKE
        cmakeTruePath = "\n " + CMAKE_Lib_dir_path +"/cmake ";
    }
    else{
        cmakeTruePath = "\n cmake " ;
    }

}
void InstallationDialog::on_pushButton_2_clicked()
{
    PackageName = "MPI Install ";
    QString a = GetChoosenDirFromDialog(1);
    if(a != ""){
        MPI_Lib_dir_path = a;
        ui->lineEdit_2->setText(MPI_Lib_dir_path);
    }
}
void InstallationDialog::on_pushButton_3_clicked()
{
    PackageName = "ROOT Install cmake ";
    QString a = GetChoosenDirFromDialog(1);
    if(a != ""){
        Root_Lib_dir_path = a;
        ui->lineEdit_3->setText(Root_Lib_dir_path);
    }
}
void InstallationDialog::on_pushButton_4_clicked()
{
    PackageName = "DCMTK Install cmake";
    QString a = GetChoosenDirFromDialog(1);
    if(a != ""){
        DCMTK_Lib_dir_path = a;
        ui->lineEdit_4->setText(DCMTK_Lib_dir_path);
    }
}
void InstallationDialog::on_pushButton_5_clicked()
{
    PackageName = "DoseCalcs Source ";
    QString a = GetChoosenDirFromDialog(1);
    if(a != ""){
        DoseCalcsCore_source_dir_path = a;
        ui->lineEdit_5->setText(DoseCalcsCore_source_dir_path);
    }

    QDir* SourceDirectory = new QDir(DoseCalcsCore_source_dir_path);

    if(SourceDirectory->exists()){
        QDir* dir = new QDir(DoseCalcsCore_source_dir_path) ; dir->cdUp();
        DoseCalcsCore_build_dir_path = dir->absolutePath()+"/"+DoseCalcs_build_dir_name;
        if(!QFile::exists(DoseCalcsCore_build_dir_path)){ // create install dir for DCMTK
            QDir* dir = new QDir(DoseCalcsCore_source_dir_path) ; dir->cdUp();
            dir->mkdir(DoseCalcs_build_dir_name);
            showResultsOutput("DoseCalcs build directory " + DoseCalcsCore_build_dir_path + " is created", 4);
        }
        else{
            showResultsOutput("DoseCalcs build directory " + DoseCalcsCore_build_dir_path + " is already existed", 4);
        }
    }
}
void InstallationDialog::on_pushButton_6_clicked()
{
    PackageName = "DoseCalcs default input macros file ";

    //QString pp = GetChoosenDirFromDialog(0);
    //ui->lineEdit_11->setText(pp);

    chosen_Dir = QFileDialog::getOpenFileName(
                this,
                tr(PackageName.toStdString().c_str()),
                GUIPackagesAndFilesDirPath+"/"+DoseCalcsCore_build_dir_path+"/Scripts",
                "All files (*.*)"//;;Text files (*.txt)" // extentions to show
                );

    if(chosen_Dir != "" || !chosen_Dir.isEmpty()){
        ui->lineEdit_11->setText(chosen_Dir);
    }

}

void InstallationDialog::on_checkBoxUseCompiledSource_clicked()
{
    if(ui->checkBoxUseCompiledSource->isChecked()){
        ui->InstallROOTButton->setText("Choose ROOT source ");
        ui->openAndInstallROOTBtn->setEnabled(false);
    }
    else{
        ui->openAndInstallROOTBtn->setText("Generate/Install");
        ui->openAndInstallROOTBtn->setEnabled(true);
    }
}

void InstallationDialog::on_pushButtonSaveInput_clicked()
{
    if(EditFlag == 0){

    }
    else if (EditFlag == 1){ // exe.sh file text saving
        filesManagerObj->WriteTextToFile( Prequests_sh_path , ui->textEditInput->toPlainText());
    }
    else if (EditFlag == 2){ // exe.sh file text saving
        filesManagerObj->WriteTextToFile( CMAKE_sh_path , ui->textEditInput->toPlainText());
    }
    else if (EditFlag == 3){ // exe.sh file text saving
        filesManagerObj->WriteTextToFile( Geant4_sh_path , ui->textEditInput->toPlainText());
    }
    else if (EditFlag == 4){ // exe.sh file text saving
        filesManagerObj->WriteTextToFile( MPI_sh_path , ui->textEditInput->toPlainText());
    }
    else if (EditFlag == 5){ // exe.sh file text saving
        filesManagerObj->WriteTextToFile( Root_sh_path , ui->textEditInput->toPlainText());
    }
    else if (EditFlag == 6){ // exe.sh file text saving
        filesManagerObj->WriteTextToFile( DCMTK_sh_path , ui->textEditInput->toPlainText());
    }
    else if (EditFlag == 7){ // exe.sh file text saving
        filesManagerObj->WriteTextToFile( App_sh_path , ui->textEditInput->toPlainText());
    }
    else if (EditFlag == 10){ // exe.sh file text saving
        filesManagerObj->WriteTextToFile( Geant4Example_sh_path , ui->textEditInput->toPlainText());
    }
    else if (EditFlag == 11){ // exe.sh file text saving
        filesManagerObj->WriteTextToFile( Geant4Example_MacrosFile_path , ui->textEditInput->toPlainText());
    }
    else if (EditFlag == 8){ // exe.sh file text saving
        filesManagerObj->WriteTextToFile( AllPackagesInstall_sh_path , ui->textEditInput->toPlainText());
    }
    else if(EditFlag == 9){
        filesManagerObj->WriteTextToFile( GUIPackagesAndFilesDirPath+"/"+GUIConfigFileName , ui->textEditInput->toPlainText());
    }

    ui->tabWidget->setTabText(1,"Input");
    ui->textEditInput->clear();
    EditFlag = 0;
}

void InstallationDialog::on_pushButtonGenerate_clicked()
{
    EditFlag = 1;

    create_Xerces_and_PrerequestGeant4Install_sh_file();
    filesManagerObj->WriteTextToFile(Prequests_sh_path,Prequest_text_shFile);
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(Prequests_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"CMAKE and Prerequisites installer");
}
void InstallationDialog::on_pushButtonGEANTGenerate_clicked()
{
    EditFlag = 3;

    create_Geant4Install_sh_file();
    filesManagerObj->WriteTextToFile(Geant4_sh_path,Geant4_text_shFile);
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(Geant4_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"GEANT4 installer");
}
void InstallationDialog::on_pushButtonMPIGenerate_clicked()
{
    EditFlag = 4;

    create_MPIInstall_sh_file();
    filesManagerObj->WriteTextToFile(MPI_sh_path,MPI_text_shFile);
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(MPI_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"MPI installer");
}
void InstallationDialog::on_pushButtonROOTGenerate_clicked()
{
    EditFlag = 5;

    create_RootInstall_sh_file();
    filesManagerObj->WriteTextToFile(Root_sh_path,Root_text_shFile);
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(Root_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"ROOT installer");
}
void InstallationDialog::on_pushButtonDCMTKGenerate_clicked()
{
    EditFlag = 6;

    create_DCMTK_Install_sh_file();
    filesManagerObj->WriteTextToFile(DCMTK_sh_path,DCMTK_text_shFile);
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(DCMTK_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"DCMTK installer");
}
void InstallationDialog::on_pushButtonGenerateDoseCalcsBuilder_clicked()
{
    EditFlag = 7;

    create_AppInstall_sh_file();
    filesManagerObj->WriteTextToFile(App_sh_path,App_text_shFile);
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(App_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"DoseCalcs Builder");
}

void InstallationDialog::on_pushButton_9_clicked()
{
    if(ui->groupBox_2->isVisible()){
        ui->groupBox_2->setVisible(false);
        ui->pushButton_9->setText("Installations < ");
    }else{
        ui->groupBox_2->setVisible(true);
        ui->pushButton_9->setText("Installations > ");
    }
}
void InstallationDialog::on_ClearTerButton_clicked()
{
    if( ui->tabWidget->currentIndex() == 0){
        QString BashCommandsForExecuting = "#! /bin/bash \n "
                                               "cd "+DoseCalcsCore_build_dir_path+"\n"
                                               "bash \n ";

        filesManagerObj->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);
        execShProcess(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
    }else if( ui->tabWidget->currentIndex() == 1){
        ui->textEditInput->clear();
    }else if( ui->tabWidget->currentIndex() == 2){
        ui->installationDialogTextOutput->clear();
    }

}

void InstallationDialog::setCompleters(){

    completer1 = new QCompleter(this);
    completer1->setModel(new QStringListModel(CompleterWords, completer1));
    completer1->setModelSorting(QCompleter::CaseInsensitivelySortedModel);
    completer1->setCaseSensitivity(Qt::CaseInsensitive);
    completer1->setWrapAround(false);
    ui->textEditInput->setCompleter(completer1);
}


void InstallationDialog::closeEvent(QCloseEvent *event)  // show prompt when user wants to close app
{
    event->ignore();
    if (QMessageBox::Yes == QMessageBox::question(this, "Close Confirmation", "Exit?", QMessageBox::Yes | QMessageBox::No))
    {
        event->accept();
    }
}
void InstallationDialog::keyPressEvent(QKeyEvent *e) {
    if(e->key() != Qt::Key_Escape) {QDialog::keyPressEvent(e);}
}




void InstallationDialog::on_pushButtonOpenPackagesFilesDir_clicked()
{
    if(QFile::exists(GUIPackagesAndFilesDirPath)){
        QString command = GUIPackagesAndFilesDirPath;
        QProcess process;
        QStringList qsl = {command};
        process.startDetached("nautilus", qsl);
    }else{

        QString nn = "The directory "+ GUIPackagesAndFilesDirName + " will be created";
        if(QMessageBox::Yes == QMessageBox::question(this, tr("Directory not found !!"), nn)){
            QDir(QDir::currentPath()).mkdir("GUIPackagesAndFilesDirName");
        }else{

        }
    }
}
void InstallationDialog::on_pushButtonDownloadSupplement_clicked()
{
    if(QMessageBox::Yes == QMessageBox::question(this, tr("DoseCalcs Supplementaries installation"), "When you click on \"Yes\", the .tar.xz file containing the ICRPDATA and PreDefinedGeometry directories will be downloaded from Google Drive. Then, it will be unziped under the directory "+
                                                 GUIPackagesAndFilesDirPath+", the download and unzip will be done on terminal. After finishing, please check the files in "+ GUIPackagesAndFilesDirName +" directory.")){

        //DownloadFileToDirectory(GUIPackagesAndFilesDirPath, "https://drive.google.com/file/d/1arU9aJpi7M5VehO1lPKOSM1JkyRe7j2k/view?usp=sharing");

        //https://drive.google.com/file/d/1arU9aJpi7M5VehO1lPKOSM1JkyRe7j2k/view?usp=drive_link from imttarikk
        //https://drive.google.com/file/d/1SMKv5_AlP2zwMg6QOT6BldR-y02fGVXj/view?usp=drive_link from telghabzouri@uae.ac.ma
        QString text_shFile =
                "cd "+GUIPackagesAndFilesDirPath+"\n"+
                "rm -r "+ GUIPackagesAndFilesDirPath+"/DoseCalcsSupplementaries.tar.xz"+ +"\n"
                //"wget -O DoseCalcsSupplementaries.tar.xz --no-check-certificate -r 'https://drive.google.com/uc?export=download&id=1arU9aJpi7M5VehO1lPKOSM1JkyRe7j2k' \n"+
                //"wget --load-cookies /tmp/cookies.txt \"https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1arU9aJpi7M5VehO1lPKOSM1JkyRe7j2k' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\\1\\n/p')&id=1arU9aJpi7M5VehO1lPKOSM1JkyRe7j2k\" -O DoseCalcsSupplementaries.tar.xz && rm -rf /tmp/cookies.txt \n"+
                "sudo apt-get install -y curl\n"
                "curl -L \"https://drive.usercontent.google.com/download?id=1arU9aJpi7M5VehO1lPKOSM1JkyRe7j2k&export=download&confirm=y\" -o \"DoseCalcsSupplementaries.tar.xz\"\n"

                //"wget "+ Supplementary_Url_String +"\n"+
                "tar xvf "+ GUIPackagesAndFilesDirPath+"/DoseCalcsSupplementaries.tar.xz"+ +"\n"
                ;
        filesManagerObj->WriteTextToFile(Supplementary_sh_path,text_shFile);

        showResultsOutput("DoseCalcs Supplementary Download shell file : " + Supplementary_sh_path + " is created." , 0);
        execShProcess( Supplementary_sh_path);

    }else{

    }
}


void InstallationDialog::on_pushButtonBuildRunGeant4ExampleShow_clicked()
{
    if(ui->groupBox_4->isVisible()){
        ui->groupBox_4->setVisible(false);
        ui->pushButtonBuildRunGeant4ExampleShow->setText("Build, Run Geant4 application < ");
    }else{
        ui->groupBox_4->setVisible(true);
        ui->pushButtonBuildRunGeant4ExampleShow->setText("Build, Run Geant4 application > ");
    }
}
void InstallationDialog::on_pushButtonAppSourceDir_clicked()
{
    PackageName = "Geant4 Example Source ";
    QString a = GetChoosenDirFromDialog(1);
    if(a != ""){
        ExampleGeant4Code_source_dir_path = a;
        ui->pushButtonAppSourceDir->setToolTip(ExampleGeant4Code_source_dir_path);
    }

    QDir* SourceDirectory = new QDir(ExampleGeant4Code_source_dir_path);

    if(SourceDirectory->exists()){
        QDir* dir = new QDir(ExampleGeant4Code_source_dir_path) ; dir->cdUp();
        DoseCalcsCore_build_dir_path = dir->absolutePath()+"/"+DoseCalcs_build_dir_name;
        if(!QFile::exists(DoseCalcsCore_build_dir_path)){ // create install dir for DCMTK
            QDir* dir = new QDir(DoseCalcsCore_source_dir_path) ; dir->cdUp();
            dir->mkdir(DoseCalcs_build_dir_name);
            showResultsOutput("DoseCalcs build directory " + DoseCalcsCore_build_dir_path + " is created", 4);
        }
        else{
            showResultsOutput("DoseCalcs build directory " + DoseCalcsCore_build_dir_path + " is already existed", 4);
        }
    }
}
void InstallationDialog::on_pushButtonAppBuildDir_clicked()
{
    PackageName = "Geant4 Example Build ";
    QString a = GetChoosenDirFromDialog(1);
    if(a != ""){
        ExampleGeant4Code_build_dir_path = a;
        ui->pushButtonAppBuildDir->setToolTip(ExampleGeant4Code_build_dir_path);
    }
}
void InstallationDialog::on_pushButtonGenerateShellFileForExampleBuilding_clicked()
{
    EditFlag = 10;

    create_Geant4ExampleInstall_sh_file();
    filesManagerObj->WriteTextToFile(Geant4Example_sh_path,Geant4Example_text_shFile);
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(Geant4Example_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"Geant4 Example Builder");
}
void InstallationDialog::on_pushButtonLoadInstaller_clicked()
{
    EditFlag = 10;

    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(Geant4Example_sh_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"Geant4 Example Builder");
}
void InstallationDialog::on_pushButtonBuildExample_clicked()
{
    if(!QFile::exists(Geant4Example_sh_path)){
        create_Geant4ExampleInstall_sh_file();
        filesManagerObj->WriteTextToFile(Geant4Example_sh_path,Geant4Example_text_shFile);
    }

    if(QFile::exists(Geant4Example_sh_path)){
        //ui->lineEdit_5->setText(App_sh_path);
        showResultsOutput("Installation file : " + Geant4Example_sh_path + " is created." , 4);

        execShProcess( Geant4Example_sh_path);

    }

}
void InstallationDialog::on_pushButtonRunExample_clicked()
{

    QDialog * d = new QDialog(); d->setWindowTitle("Geant4Example Execution Command");

    QGridLayout* GraphLayout = new QGridLayout;

    QLineEdit * lineStyles = new QLineEdit();
    lineStyles->setText(ExampleGeant4Code_ExecutionCommand);
    lineStyles->setToolTip("Add command to Run Geant4 Example");

    int ii = 0, jj=0;
    GraphLayout->addWidget(lineStyles, jj,ii,1,1);

    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QObject::connect(buttonBox, SIGNAL(accepted()), d, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), d, SLOT(reject()));

    GraphLayout->addWidget(buttonBox);

    d->setLayout(GraphLayout);

    int result = d->exec();

    if(result == QDialog::Accepted)
    {
        ExampleGeant4Code_ExecutionCommand = lineStyles->text();
    }else{
        return;
    }


    if(ExampleGeant4Code_ExecutionCommand =="" ||ExampleGeant4Code_ExecutionCommand.isEmpty() ){
        return;
    }


    QString BashCommandsForExecuting = "#! /bin/bash \n cd "
            +geant4_Lib_dir_path +"\n"+
            + ". ./geant4.sh\n cd "
            + ExampleGeant4Code_build_dir_path + "\n"
            + ExampleGeant4Code_ExecutionCommand
            + "\n bash \n";

    showResultsOutput("Writing Run Commands : \n", 0);
    showResultsOutput(BashCommandsForExecuting , 0);
    showResultsOutput("to --> " + DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , 4);
    filesManagerObj->WriteTextToFile( DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName , BashCommandsForExecuting);

    if(QFile::exists(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName)){

        showResultsOutput("Computation Run" , 1);
        execShProcess(DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName);
    }
    else{
        showResultsOutput("Cannot find file containing execution commands "+ DoseCalcsCore_build_dir_path+"/"+DoseCalcsExecutingFileName + " , you should build example or something else.." , 3);
    }

}
void InstallationDialog::on_pushButtonOpenBuildDir_clicked()
{
    QString command = ExampleGeant4Code_build_dir_path;
    QProcess process;
    QStringList qsl = {command};
    process.startDetached("nautilus", qsl);
}

void InstallationDialog::on_pushButtonEditAFile_clicked()
{

    Geant4Example_MacrosFile_path = QFileDialog::getOpenFileName(this,tr("Geant4Example File"), ExampleGeant4Code_build_dir_path, "All files (*.*)");

    if(Geant4Example_MacrosFile_path != ""){
        ui->pushButtonEditAFile->setToolTip(Geant4Example_MacrosFile_path);
    }else{
        return;
    }

    EditFlag = 11;
    ui->textEditInput->setPlainText(filesManagerObj->ReadTextFromFileInOneString(Geant4Example_MacrosFile_path));
    ui->tabWidget->setCurrentIndex(1);
    ui->tabWidget->setTabText(1,"Geant4Example File");

}

