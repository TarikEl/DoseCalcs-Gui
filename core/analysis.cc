#include <iostream>
//#include "globals.hh"
#include "G4DoseCalcsAnalysis.hh"
#include "G4RandomTools.hh"

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

std::string MacrosStartingFile;

char result[PATH_MAX];
ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
std::string getProgDirpath()
{
    std::string path;
    if (count != -1) { path = dirname(result); }
    return path;
}
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

std::string TissueFacFilePath = "TissueRadiationFactors.mac";

using namespace std;

std::string appBuildDir;
G4String ResultDirectoryPath ;

int main(int argc,char** argv){

    G4Random::setTheEngine(new CLHEP::MixMaxRng);

    appBuildDir = getProgDirpath();
    ResultDirectoryPath = appBuildDir+"/Results";

    G4DoseCalcsAnalysis* ResultCalculation;
    ResultCalculation = new G4DoseCalcsAnalysis(ResultDirectoryPath);

    if(argc == 4){
        MacrosStartingFile = argv[3];
        ResultCalculation->GenerateGraphFromROOTGraphDATA();
        return 0;
    }

    if(argc > 1){
        MacrosStartingFile = argv[1];
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

    //std::cout << "\n========= DoseCalcs Build Directory "<< appBuildDir << " ========= "<< "\n" << G4endl;


    ResultCalculation->ReadSimulationData();
    ResultCalculation->ReadResultsAndReferenceData();

    ResultCalculation->createLatexTables();
    ResultCalculation->createDataInCSVFormat();

#ifdef ANALYSIS_USE
    ResultCalculation->GenerateROOTFileForAllResults();
    ResultCalculation->GenerateSelfCrossGraphs();
    //std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n==== GenerateCrossSectionGraphs ========= "<< "\n" << G4endl;

    ResultCalculation->GenerateCrossSectionGraphs();
    ResultCalculation->GenerateSourceComputationTimeGraph();

//#ifdef VOX_USE
    ResultCalculation->GenerateVoxel2DGraphs();
//#endif

#endif

    ResultCalculation->GenerateEventsDataHisto();

    delete ResultCalculation;

    return 0;
}
