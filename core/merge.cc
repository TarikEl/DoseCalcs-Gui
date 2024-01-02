#include <iostream>
#include "G4TResultCalculation.hh"

#include <dirent.h>         // dirname
#include <libgen.h>         // dirname
#include <unistd.h>         // readlink
#include <linux/limits.h>   // PATH_MAX
std::string getProgpath()
{
    char result[ PATH_MAX ];
    ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
    return std::string( result, (count > 0) ? count : 0 );
}

char result[PATH_MAX];
ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
std::string getProgDirpath()
{
    std::string path;
    if (count != -1) { path = dirname(result); }
    return path;
}

G4String* CopyNumberRegionNameMap;
G4float* CopyNumberXPos;
G4float* CopyNumberYPos;
G4float* CopyNumberZPos;
G4float* CopyNumberMassSize;

G4String ResultDirectoryPath ;
G4String MacrosStartingFile ;

using namespace std;

G4String appBuildDir;

bool DirectoryExists( const char* pzPath )
{
    if ( pzPath == NULL) return false;

    DIR *pDir;
    bool bExists = false;
    pDir = opendir (pzPath);
    if (pDir != NULL)
    {
        bExists = true;
        (void) closedir (pDir);
    }
    return bExists;
}

int main(int argc,char** argv){

    appBuildDir = getProgDirpath().c_str();//(getenv("PWD"));
    ResultDirectoryPath = appBuildDir+"/Results";

    //G4cout << "\n========= DoseCalcs Build Directory "<< appBuildDir << " ========= "<< "\n" << G4endl;

    G4TResultCalculation* ResultCalculation = new G4TResultCalculation();
    ResultCalculation->setExeFromMerge(true);
    ResultCalculation->setUseAllResultsFiles(true);

    //G4cout << argv[1] << " argc " << argc << " - " << G4endl;

    if(argc > 1){
        MacrosStartingFile = argv[1];
        //if(DirectoryExists( MacrosStartingFile )){}
        //else{ G4cout << " Cannot find " << MacrosStartingFile << G4endl; return 0;}
    }else{
        G4cout << " Set the macros file as first argument " << G4endl;
        return 0;
    }

    if(argc > 2){
        G4String v = argv[2];
        if(v == "v"){ResultCalculation->setV(true);}
        else{ResultCalculation->setV(false);}
    }else{
        ResultCalculation->setV(false);
    }

    ResultCalculation->ReadSimulationData();
    ResultCalculation->Initialization();
    ResultCalculation->MergeSimulationsData();

	return 0;
}
