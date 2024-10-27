
#include "G4TResultCalculation.hh"
//#include "G4SystemOfUnits.hh"
//#include "G4UnitsTable.hh"

#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <fstream>

extern G4String appBuildDir;

// this are constants between all threads the all threads have the same values
extern  G4String* CopyNumberRegionNameMap;
extern  G4double* CopyNumberXPos;
extern  G4double* CopyNumberYPos;
extern  G4double* CopyNumberZPos;
extern  G4double* CopyNumberMassSize;

extern G4String ResultDirectoryPath ;
extern G4String MacrosStartingFile ;
bool G4TResultCalculation::DirectoryExists( const char* pzPath )
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

G4TResultCalculation::G4TResultCalculation(){
    
    V = false;
    
    if(V)G4cout << "\n\n========= RESULTS MERGING =====================\n" << G4endl;
    
    SZ = 17;
    VarSZ = 8;
    DefaultSZ = std::cout.precision();

    VOX_USE = false;
    GenerateVoxelsResuls = false;
    GeometrySymbol = "phantom0";
    Physics = "EMS0";

    ExecutionMode = "MPI";
    RankID = 0;

    DoseCalcsQuantities.push_back("AE");DoseCalcsQuantities.push_back("AF");DoseCalcsQuantities.push_back("SAF"); DoseCalcsQuantities.push_back("DR");
    DoseCalcsQuantities.push_back("AD");DoseCalcsQuantities.push_back("S");DoseCalcsQuantities.push_back("H");DoseCalcsQuantities.push_back("E");
    DoseCalcsQuantities.push_back("ER");DoseCalcsQuantities.push_back("DCC");

    ErrorMessage = "\n\n========= Errors and Warnings =====================\n";

    /*
    TissueFactorMap["World"]=1;
    TissueFactorMap["Gonads"]=0.08;
    TissueFactorMap["Thyroid"]=0.04;
    TissueFactorMap["UrinaryBladder"]=0.04;
    TissueFactorMap["Oesophagus"]=0.04;
    TissueFactorMap["Liver"]=0.04;
    TissueFactorMap["Brain"]=0.01;
    TissueFactorMap["SalivaryGlands"]=0.01;
    TissueFactorMap["Skin"]=0.01;
    TissueFactorMap["ArmBone"]= 0.01;
    TissueFactorMap["LegBone"]= 0.01;
    TissueFactorMap["VOXEL"]= 0.12;
    TissueFactorMap["Others"]= 0.12;
    */

    //RadiationFactorMap["gamma"][1.]= 1.;
    //RadiationFactorMap["e-"][1.]= 1.;
    //RadiationFactorMap["muon"][1.]=  1.;
    //RadiationFactorMap["proton"][1.]= 2.;
    //RadiationFactorMap["pion"][1.]=  2.;
    //RadiationFactorMap["alpha"][1.]= 20.;
    //RadiationFactorMap["neutron"][1.]= 5.;

    MeV_to_J = 1.60218e-13;
    Gy_to_Sv = 1. ;
    Bq_to_MBq = 1e-6 ;
    
    //GenerateVoxelsResuls = true;
    VOX_USE = true;
    GETVOXELDATA = false;

    // Default code values
    // AE MeV
    // Mass kg
    // Volume cm3
    // Density g/cm3
    
    // Default quantities units
    // AE MeV
    // AF
    // SAF kg-1
    // AD Gy
    // S Gy/Bq
    // H Sv
    // E Sv
    // DCC pGy/cm2

    AEUnitFactor = 1.;
    AFUnitFactor = 1.;
    SAFUnitFactor = 1.;
    SUnitFactor = 1.;
    ADUnitFactor = 1.;
    HUnitFactor = 1.;
    EUnitFactor = 1.;
    DCCUnitFactor = MeV_to_J*1e+12 /*0.16021773*/;

    //G4cout << "DCCUnitFactor : " << DCCUnitFactor << G4endl;

    AEUnit = "MeV";
    AFUnit = "";
    SAFUnit = "kg-1";
    ADUnit = "MeV/kg";
    SUnit = "MeV/kg";
    HUnit = "MeV/kg";
    EUnit = "MeV/kg";
    DCCUnit = "pGy/cm2";

    UnitPerParticle = "Particle";
    UnitPerRadioTracerDecay = "Decay";
    UnitPerRadionuclideDecay = "Decay";

    ResidenceTimeUnitFactor = 1.;
    ResidenceTimeUnit = "s";
    
    AdministeredActivityUnitFactor = 1.;
    AdministeredActivityUnit = "Bq";
    
    intOfEneForRadFac = 1;
    //RadiationFactorMap[ParticleName][ParticleSourceEnergy] = 1.; // Name of particle, interval of energy, the radiation factor value
    PhantomEffectiveDose = 0.;
    
    RadiotracerDataFomFile = false;
    GenerateResultsForRadioTracer = false;
    GenerateResultsForRadioTracerExams = false;

}
void G4TResultCalculation::Initialization(){
    
    if(V)G4cout << "\n\n========= Initialization for new Particle-Source-energy results calculations =====================\n" << G4endl;
    
    OneEventExecutionTimeInMs = 0. ;
    ExecutionTimeInMin = 0.;
    MaxExecutionTimeInMin = 0.;
    TotalEventNumber = 0 ;
    TotalEmittedEnergy = 0. ;
    
    OrgansNameVector.clear(); // should be in Initialization() like for normal sources and not just for combinated sources
    ED_Total.clear();
    ED2_Total.clear();
    NOfValues.clear();
    Fluence.clear();
    ReadedResultFilesPaths.clear();
    
    if(SourceType == "Voxels"){
        VOX_USE = true;
    }

    if(VOX_USE && GenerateVoxelsResuls){
        
        TotVoxNum = VoxXNumber * VoxYNumber * VoxZNumber;
        
        //CNID = new unsigned int [TotVoxNum];
        
        //if (ExeFromMerge){

        if(V)G4cout << "\n-------- " << "The array CopyNumberMassSize is not defined. Execution comes from merge executable " << "\n" << G4endl ;
        if(V)G4cout << "\n-------- for " << TotVoxNum << " voxel, the array CopyNumberMassSize is not defined yet. Execution comes from ./merge executable " << "\n" << G4endl ;

        CopyNumberZPos = new G4double [TotVoxNum];
        CopyNumberYPos = new G4double [TotVoxNum];
        CopyNumberXPos = new G4double [TotVoxNum];
        CopyNumberMassSize = new G4double [TotVoxNum];
        CopyNumberRegionNameMap = new G4String [TotVoxNum];

        for( G4int C = 0 ; C < TotVoxNum ; C++ ){
            CopyNumberZPos[C] = 0;
            CopyNumberYPos[C] = 0;
            CopyNumberXPos[C] = 0 ;
            CopyNumberMassSize[C] = 0;
            CopyNumberRegionNameMap[C] = "VOXEL";
        }
        //}else{
        //    if(V)G4cout << "\n-------- " << "The array CopyNumberMassSize is already defined. Execution comes from ./simulate executable " << "\n" << G4endl ;
        //}
        
        VoxED_Total = new G4double  [TotVoxNum];
        VoxED2_Total = new G4double [TotVoxNum];
        VoxNOfValues  = new unsigned long long int[TotVoxNum];
        VoxFluence  = new unsigned long long int[TotVoxNum];

        VoxED_Mean = new G4double [TotVoxNum];
        VoxED2_Mean = new G4double [TotVoxNum];
        VoxED_Var = new G4double [TotVoxNum];
        VoxED_SDev = new G4double [TotVoxNum];
        VoxED_RelS_D = new G4double [TotVoxNum];
        
        VoxAF_Total = new G4double [TotVoxNum];
        VoxSAF_Total = new G4double [TotVoxNum];
        VoxAD_Total = new G4double [TotVoxNum];
        VoxS_Total = new G4double [TotVoxNum];
        VoxH_Total = new G4double [TotVoxNum];
        VoxE_Total = new G4double [TotVoxNum];
        
        for( G4int C = 0 ; C < TotVoxNum ; C++ ){
            VoxED_Total[C] = 0;
            VoxED2_Total[C] = 0;
            VoxNOfValues[C] = 0 ;
            VoxFluence[C] = 0 ;

            VoxED_Mean[C] = 0;
            VoxED2_Mean[C] = 0;
            VoxED_Var[C] = 0;
            VoxED_SDev[C] = 0;
            VoxED_RelS_D[C] = 0;
            
            VoxAF_Total[C] = 0;
            VoxSAF_Total[C] = 0;
            VoxAD_Total[C] = 0;
            VoxS_Total[C] = 0;
            VoxH_Total[C] = 0;
            VoxE_Total[C] = 0;
        }
    }
    //SourceParticleEnergyValues[GeometrySymbol].clear();SourceParticleEnergyValues[GeometrySymbol].empty();
}

void G4TResultCalculation::InitializeVoxelizedData(){

    TotVoxNum = VoxXNumber * VoxYNumber * VoxZNumber;

    //CNID = new unsigned int [TotVoxNum];

    if(V)G4cout << "\n-------- " << "The array CopyNumberMassSize is not defined. Execution comes from merge executable " << "\n" << G4endl ;
    if(V)G4cout << "\n-------- for " << TotVoxNum << " voxel, the array CopyNumberMassSize is not defined yet. Execution comes from ./merge executable " << "\n" << G4endl ;

    CopyNumberZPos = new G4double [TotVoxNum];
    CopyNumberYPos = new G4double [TotVoxNum];
    CopyNumberXPos = new G4double [TotVoxNum];
    CopyNumberMassSize = new G4double [TotVoxNum];
    CopyNumberRegionNameMap = new G4String [TotVoxNum];

    for( G4int C = 0 ; C < TotVoxNum ; C++ ){
        CopyNumberZPos[C] = 0;
        CopyNumberYPos[C] = 0;
        CopyNumberXPos[C] = 0 ;
        CopyNumberMassSize[C] = 0;
        CopyNumberRegionNameMap[C] = "VOXEL";
    }

    VoxED_Total = new G4double  [TotVoxNum];
    VoxED2_Total = new G4double [TotVoxNum];
    VoxNOfValues  = new unsigned long long int[TotVoxNum];
    VoxFluence  = new unsigned long long int[TotVoxNum];

    VoxED_Mean = new G4double [TotVoxNum];
    VoxED2_Mean = new G4double [TotVoxNum];
    VoxED_Var = new G4double [TotVoxNum];
    VoxED_SDev = new G4double [TotVoxNum];
    VoxED_RelS_D = new G4double [TotVoxNum];

    VoxAF_Total = new G4double [TotVoxNum];
    VoxSAF_Total = new G4double [TotVoxNum];
    VoxAD_Total = new G4double [TotVoxNum];
    VoxS_Total = new G4double [TotVoxNum];
    VoxH_Total = new G4double [TotVoxNum];
    VoxE_Total = new G4double [TotVoxNum];

    for( G4int C = 0 ; C < TotVoxNum ; C++ ){
        VoxED_Total[C] = 0;
        VoxED2_Total[C] = 0;
        VoxNOfValues[C] = 0 ;
        VoxFluence[C] = 0 ;

        VoxED_Mean[C] = 0;
        VoxED2_Mean[C] = 0;
        VoxED_Var[C] = 0;
        VoxED_SDev[C] = 0;
        VoxED_RelS_D[C] = 0;

        VoxAF_Total[C] = 0;
        VoxSAF_Total[C] = 0;
        VoxAD_Total[C] = 0;
        VoxS_Total[C] = 0;
        VoxH_Total[C] = 0;
        VoxE_Total[C] = 0;
    }

    GETVOXELDATA = true;

}

void G4TResultCalculation::G4TCoutReset(){
    std::resetiosflags( G4cout.basefield ); // clears integer manipulations
    std::resetiosflags( G4cout.floatfield  ); // clears floating-point manipulations
    std::resetiosflags( G4cout.flags() ); // clears all flags
}

void G4TResultCalculation::MergeSimulationsData(){
    
    if(UseAllResultsFiles == true){
        
        std::cout << "\n\n The resulted files are read from " << ResultDirectoryPath <<  "\n\n"<< std::endl;
        
        GeometrySourceDataRankThreadDataFile.clear();
        GeometrySourceDataRankThreadDataFile.empty();
        //SourceDataRankThreadDataFile.empty();
        //SourceDataRankThreadDataFile.clear();
        if(SourceDataRankThreadDataFile.empty()){
            //std::cout << "\n\n\n\ndddddddddddddddddddddddddddddddd" << std::endl;
        }
        
        DIR *dir; struct dirent *diread;
        std::vector<G4String> ResultsFileNames;
        
        if ((dir = opendir(ResultDirectoryPath.c_str())) != nullptr) {
            while ((diread = readdir(dir)) != nullptr) {
                std::string str1 = diread->d_name;
                std::string Idc = str1.substr(0, 3);
                //std::cout << str1 << "  " << Idc << std::endl;
                
                if(Idc == "AE@"){
                    
                    ResultsFileNames.push_back(diread->d_name);
                    
                    std::vector<std::string> result;
                    std::stringstream ss (diread->d_name);
                    std::string item;
                    while (getline (ss, item, '@')) {
                        
                        result.push_back(item);
                        //std::cout << item << " ";
                    }
                    
                    //std::cout  << result.size() << "\n";
                    
                    // forgot the files of uneded sources
                    if(IsAllSourcesToScore == false){
                        bool IsIn = false;
                        for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                            if(SourceNamesToScore[gg] == result[8]){
                                IsIn = true; break;
                            }
                        }


                        for ( auto Obeg = NeWCombinationsForSourceNamesToScore.begin(); Obeg != NeWCombinationsForSourceNamesToScore.end(); ++Obeg  )
                        {
                            for (int gg = 0 ; gg < Obeg->second.size() ; gg++) {
                                if(Obeg->second[gg] == result[8]){
                                    IsIn = true; break;
                                }
                            }
                        }

                        if(IsIn == false){
                            continue;
                        }

                    }
                    //std::cout  << result.size() << "\n";

                    if(result.size() == 10){
                        //std::cout << result[6] << " "  << result[7] << " " << result[8] << " " << result[9]<< " " << result[3] << " " << result[5] << " " << GeometrySourceDataRankThreadDataFile.size() << std::endl ;
                        GeometrySourceDataRankThreadDataFile[result[6]     ][result[7]][result[8]][std::stod(result[9])][std::atoi(result[3].c_str())].push_back(std::atoi(result[5].c_str())); //=std::atoi(result[5].c_str());
                        GeometryEnergySourceDataRankThreadDataFile[result[6]     ][result[7]][std::stod(result[9])][result[8]][std::atoi(result[3].c_str())].push_back(std::atoi(result[5].c_str())); //=std::atoi(result[5].c_str());
                    }
                    else if(result.size() == 9){
                        GeometrySourceDataRankThreadDataFile[GeometrySymbol][result[6]][result[7]][std::stod(result[8])][std::atoi(result[3].c_str())].push_back(std::atoi(result[5].c_str())); //=std::atoi(result[5].c_str());
                        GeometryEnergySourceDataRankThreadDataFile[GeometrySymbol][result[7]][std::stod(result[9])][result[8]][std::atoi(result[3].c_str())].push_back(std::atoi(result[5].c_str())); //=std::atoi(result[5].c_str());
                        //SourceDataRankThreadDataFile[result[6]][result[7]][std::stod(result[8])][std::atoi(result[3].c_str())]=std::atoi(result[5].c_str());
                        //std::cout << result[6] << " "  << result[7] << " " << result[8] << " " << result[3] << " " << result[5] << " " << GeometrySourceDataRankThreadDataFile.size() << std::endl ;
                    }
                }
            }
            closedir (dir);
        } else {
            perror ("opendir");
        }
        //for (auto file : ResultsFileNames) std::cout << file << std::endl;
    }
    
    //std::map<G4String, std::map<G4String,std::map<G4double,std::map<G4int,G4int>>>> FilesDataMap = GeometrySourceDataRankThreadDataFile[GeometrySymbol];
    //std::map<G4String, std::map<G4String,std::map<G4double,std::map<G4int,G4int>>>> FilesDataMap = SourceDataRankThreadDataFile;
    
    if(GeometrySourceDataRankThreadDataFile.size() == 0) {
        G4cout << "\n\n No data are registered from directory " << ResultDirectoryPath << G4endl ;
        G4cout << "Check the result directory, or scored source result files if exist... " << "\n\n" << G4endl ;
        return;
    }
    if(V)G4cout << "\n========= Regions Data Calculation from " << __FUNCTION__ << "\n" << G4endl ;
    
    // iterations on particle name
    for ( auto Obeg = GeometrySourceDataRankThreadDataFile.begin(); Obeg != GeometrySourceDataRankThreadDataFile.end(); ++Obeg  )
    {
        GeometrySymbol = Obeg->first;

        for ( auto Abeg = Obeg->second.begin(); Abeg != Obeg->second.end(); ++Abeg  )
        {
            G4String Particle = Abeg->first;

            // iterations on source name
            for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  )
            {

                G4String SourceReg = Bbeg->first;
                SourceRegionName = SourceReg;

                // iterations on energy
                for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg  )
                {
                    
                    G4double Energy = Cbeg->first;
                    
                    G4cout << "\n\n================================== Reading, Calculation and Generation final results from files containing data of, Particle: " << Particle << " , SourceReg: " << SourceReg << " , Energy: " << Energy << "MeV ==================================\n" << G4endl ;
                    
                    Initialization();
                    
                    G4String fm, vfm;
                    for ( auto Dbeg = Cbeg->second.begin(); Dbeg != Cbeg->second.end(); ++Dbeg  )
                    {
                        
                        G4int rankid = Dbeg->first;
                        //G4int threadid = Dbeg->second;

                        for (G4int r = 0 ; r < Dbeg->second.size() ; r++) { // calculate data just if the file is readed
                            G4int threadid = Dbeg->second[r];

                            if(V)G4cout << " Getting results from Rank " << rankid << " and Thread " << threadid << " : " << G4endl ;

                            std::ostringstream filename;
                            //filename << ResultDirectoryPath << "/AE@for@Rank@" << rankid << "@Thread@"<< threadid<<"@" << Particle << "@" << SourceReg << "@" << Energy ;
                            filename << ResultDirectoryPath << "/AE@for@Rank@" << rankid << "@Thread@"<< threadid<<"@" << GeometrySymbol <<"@" << Particle << "@" << SourceReg << "@" << Energy ;
                            fm = filename.str().c_str();

                            if(FILE* file = fopen(fm,"r")){
                                std::fclose(file);
                            }else{
                                std::ostringstream fmm; fmm << ResultDirectoryPath << "/AE@for@Rank@" << rankid << "@Thread@"<< threadid<<"@" << Particle << "@" << SourceReg << "@" << Energy ;
                                fm = fmm.str().c_str();
                            }


                            if(ReadThreadRegionResultFile(fm)){
                                G4cout << "\nThe file to read: " << fm << G4endl ;
                            }else{
                                G4cout << "\nThe file is not read: " << fm << G4endl ;
                            }


                            //if(V)
                            //G4cout << "Vox_use "<< VOX_USE<< G4endl ;
                            if(VOX_USE && GenerateVoxelsResuls){

                                std::ostringstream name2;
                                name2 << ResultDirectoryPath << "/AE@Voxel@for@Rank@" << rankid << "@Thread@"<< threadid<<"@" << GeometrySymbol <<"@" << Particle << "@" << SourceReg << "@" << Energy ;

                                vfm = name2.str().c_str();
                                if(V)G4cout << "\nThe file to read: " << vfm << G4endl ;
                                //G4cout << "Vox_use "<< VOX_USE<< G4endl ;

                                ReadThreadVoxelResultFile(vfm);
                                //G4cout << "Vox_use "<< VOX_USE<< G4endl ;

                            }
                        }
                    }
                    
                    for (G4int r = 0 ; r < ReadedResultFilesPaths[GeometrySymbol].size() ; r++) { // calculate data just if the file is readed
                        if(fm == ReadedResultFilesPaths[GeometrySymbol][r]){
                            RegionQuantitiesCalculation();
                            GenerateRegionResultFile();
                        }
                    }
                    
                    if(VOX_USE && GenerateVoxelsResuls){
                        VoxelQuantitiesCalculation();
                        GenerateVoxelsResultFiles();
                    }

                    ED_Total.clear();
                    ED_Mean.clear();
                    ED2_Total.clear();
                    ED2_Mean.clear();
                    ED_Var.clear();
                    ED_SDev.clear();
                    ED_RelS_D.clear();

                    AF_Total.clear();
                    AF2_Total.clear();
                    AF_Mean.clear();
                    AF2_Mean.clear();
                    AF_Var.clear();
                    AF_SDev.clear();
                    AF_RelS_D.clear();

                    SAF_Total.clear();
                    SAF2_Total.clear();
                    SAF_Mean.clear();
                    SAF2_Mean.clear();
                    SAF_Var.clear();
                    SAF_SDev.clear();
                    SAF_RelS_D.clear();

                    AD_Total.clear();
                    AD2_Total.clear();
                    AD_Mean.clear();
                    AD2_Mean.clear();
                    AD_Var.clear();
                    AD_SDev.clear();
                    AD_RelS_D.clear();

                    S_Total.clear();
                    S2_Total.clear();
                    S_Mean.clear();
                    S2_Mean.clear();
                    S_Var.clear();
                    S_SDev.clear();
                    S_RelS_D.clear();

                    H_Total.clear();
                    H2_Total.clear();
                    H_Mean.clear();
                    H2_Mean.clear();
                    H_Var.clear();
                    H_SDev.clear();
                    H_RelS_D.clear();

                    E_Total.clear();
                    E2_Total.clear();
                    E_Mean.clear();
                    E2_Mean.clear();
                    E_Var.clear();
                    E_SDev.clear();
                    E_RelS_D.clear();

                    DCC_Total.clear();
                    DCC2_Total.clear();
                    DCC_Mean.clear();
                    DCC2_Mean.clear();
                    DCC_Var.clear();
                    DCC_SDev.clear();
                    DCC_RelS_D.clear();



                    DR_Total.clear();
                    ER_Total.clear();

                    ChosenVariablevariance.clear();
                    ChosenVariableStandardDeviation.clear();
                    ChosenVariableMean.clear();
                    ChosenVariableTotal.clear();

                    NOfValues.clear();
                    Fluence.clear();


                }
            }
        }
        
        if(GenerateResultsForRadioTracer){
            GenerateRegionResultForRadioTracer();
        }
        
        // ////////////////////////////////////////////////////////////
        RadioTracerQuantitySourceTargetValue.clear();
        RadioTracerQuantityOrganValue.clear();
        RadioTracerQuantityOrganVAR.clear();
        TotalDoseFromRadioTracer.clear();
        TotalValueOfQuantity.clear();
        RadioTracerQuantitySourceTargetVariance.clear();
        ResultTable.clear();
        StandardDeviation.clear();
        TotalAEForRadiotracerRSD.clear();
        RadiTracerParticleEnergyDataString.clear();
        RadiTracerDataForTotalDoseString.clear();
        OrgansNameVector.clear();
        VolumeNameMassMap.clear();
        VolumeNameVolumeMap.clear();
        VolumeNameDensityMap.clear();
        // ////////////////////////////////////////////////////////////
    }

    // ////////////////////////////////////////////////////

    // for the new Combinated source region

    for ( auto Obeg = GeometryEnergySourceDataRankThreadDataFile.begin(); Obeg != GeometryEnergySourceDataRankThreadDataFile.end(); ++Obeg  )
    {
        GeometrySymbol = Obeg->first;

        for ( auto Abeg = Obeg->second.begin(); Abeg != Obeg->second.end(); ++Abeg  )
        {
            G4String Particle = Abeg->first;

            // iterations on source name
            for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  )
            {
                G4double Energy = Bbeg->first;

                for ( auto SSbeg = NeWCombinationsForSourceNamesToScore.begin(); SSbeg != NeWCombinationsForSourceNamesToScore.end(); ++SSbeg  )
                {
                    G4String NeWSourceReg = SSbeg->first;
                    bool ThisNewSourceIsAlreadyDefined = false;
                    SourceRegionName = NeWSourceReg;

                    //VolumeNameMassMap.clear();
                    //VolumeNameVolumeMap.clear();
                    //VolumeNameDensityMap.clear();

                    //OrgansNameVector.clear(); // should be in Initialization() like for normal sources and not just for combinated sources
                    Initialization();

                    G4String fm, vfm;
                    for (G4int ddd = 0 ; ddd < SSbeg->second.size() ; ddd++) { // calculate data just if the file is readed

                        G4String SourceRegC = SSbeg->second[ddd];

                        for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg  )
                        {
                            G4String SourceReg = Cbeg->first;

                            if(NeWSourceReg == SourceReg){
                                ThisNewSourceIsAlreadyDefined = true;
                                G4cout << "\n\n================================== This new combinated source region ("<<NeWSourceReg<<") is already defined in the simulation of ("<< GeometrySymbol <<", "<< Particle <<", " << NeWSourceReg <<", " << Energy <<") as a single source region, the results calculation will be ignored for this combinated source region ==================================\n" << G4endl ;

                                break;
                            }

                            if(SourceRegC != SourceReg){continue;}

                            G4cout << "\n\n================================== Reading, Calculation and Generation final results from files containing data of, Particle: " << Particle << " , SourceReg: " << SourceReg << " (in combinated source region "<< NeWSourceReg << "), Energy: " << Energy << "MeV ==================================\n" << G4endl ;

                            for ( auto Dbeg = Cbeg->second.begin(); Dbeg != Cbeg->second.end(); ++Dbeg  )
                            {

                                G4int rankid = Dbeg->first;
                                //G4int threadid = Dbeg->second;

                                for (G4int r = 0 ; r < Dbeg->second.size() ; r++) { // calculate data just if the file is readed

                                    G4int threadid = Dbeg->second[r];

                                    if(V)G4cout << " Getting results from Rank " << rankid << " and Thread " << threadid << " : " << G4endl ;

                                    std::ostringstream filename;
                                    //filename << ResultDirectoryPath << "/AE@for@Rank@" << rankid << "@Thread@"<< threadid<<"@" << Particle << "@" << SourceReg << "@" << Energy ;
                                    filename << ResultDirectoryPath << "/AE@for@Rank@" << rankid << "@Thread@"<< threadid<<"@" << GeometrySymbol <<"@" << Particle << "@" << SourceReg << "@" << Energy ;
                                    fm = filename.str().c_str();

                                    if(FILE* file = fopen(fm,"r")){
                                        std::fclose(file);
                                    }else{
                                        std::ostringstream fmm; fmm << ResultDirectoryPath << "/AE@for@Rank@" << rankid << "@Thread@"<< threadid<<"@" << Particle << "@" << SourceReg << "@" << Energy ;
                                        fm = fmm.str().c_str();
                                    }


                                    if(ReadThreadRegionResultFile(fm)){
                                        G4cout << "\nThe file to read: " << fm << G4endl ;
                                    }else{
                                        G4cout << "\nThe file is not read: " << fm << G4endl ;
                                    }


                                    //if(V)G4cout << "Vox_use "<< VOX_USE<< G4endl ;
                                    if(VOX_USE && GenerateVoxelsResuls){

                                        std::ostringstream name2;
                                        name2 << ResultDirectoryPath << "/AE@Voxel@for@Rank@" << rankid << "@Thread@"<< threadid<<"@" << GeometrySymbol <<"@" << Particle << "@" << SourceReg << "@" << Energy ;

                                        vfm = name2.str().c_str();
                                        if(V)G4cout << "\nThe file to read: " << vfm << G4endl ;

                                        ReadThreadVoxelResultFile(vfm);
                                    }
                                }
                            }
                        }
                        if(ThisNewSourceIsAlreadyDefined == true){break;}
                    }

                    if(ThisNewSourceIsAlreadyDefined == true){continue;}

                    for (G4int r = 0 ; r < ReadedResultFilesPaths[GeometrySymbol].size() ; r++) { // calculate data just if the file is readed
                        if(fm == ReadedResultFilesPaths[GeometrySymbol][r]){
                            RegionQuantitiesCalculation();
                            GenerateRegionResultFile();
                        }
                    }

                    if(VOX_USE && GenerateVoxelsResuls){
                        VoxelQuantitiesCalculation();
                        GenerateVoxelsResultFiles();
                    }

                    ED_Total.clear();
                    ED_Mean.clear();
                    ED2_Total.clear();
                    ED2_Mean.clear();
                    ED_Var.clear();
                    ED_SDev.clear();
                    ED_RelS_D.clear();

                    AF_Total.clear();
                    AF2_Total.clear();
                    AF_Mean.clear();
                    AF2_Mean.clear();
                    AF_Var.clear();
                    AF_SDev.clear();
                    AF_RelS_D.clear();

                    SAF_Total.clear();
                    SAF2_Total.clear();
                    SAF_Mean.clear();
                    SAF2_Mean.clear();
                    SAF_Var.clear();
                    SAF_SDev.clear();
                    SAF_RelS_D.clear();

                    AD_Total.clear();
                    AD2_Total.clear();
                    AD_Mean.clear();
                    AD2_Mean.clear();
                    AD_Var.clear();
                    AD_SDev.clear();
                    AD_RelS_D.clear();

                    S_Total.clear();
                    S2_Total.clear();
                    S_Mean.clear();
                    S2_Mean.clear();
                    S_Var.clear();
                    S_SDev.clear();
                    S_RelS_D.clear();

                    H_Total.clear();
                    H2_Total.clear();
                    H_Mean.clear();
                    H2_Mean.clear();
                    H_Var.clear();
                    H_SDev.clear();
                    H_RelS_D.clear();

                    E_Total.clear();
                    E2_Total.clear();
                    E_Mean.clear();
                    E2_Mean.clear();
                    E_Var.clear();
                    E_SDev.clear();
                    E_RelS_D.clear();

                    DR_Total.clear();
                    ER_Total.clear();

                    ChosenVariablevariance.clear();
                    ChosenVariableStandardDeviation.clear();
                    ChosenVariableMean.clear();
                    ChosenVariableTotal.clear();

                    NOfValues.clear();
                    Fluence.clear();

                }
            }
        }

        if(GenerateResultsForRadioTracer){
            GenerateRegionResultForRadioTracer();
        }

        // ////////////////////////////////////////////////////////////
        RadioTracerQuantitySourceTargetValue.clear();
        RadioTracerQuantityOrganValue.clear();
        RadioTracerQuantityOrganVAR.clear();
        TotalDoseFromRadioTracer.clear();
        RadioTracerQuantitySourceTargetVariance.clear();
        ResultTable.clear();
        StandardDeviation.clear();
        TotalAEForRadiotracerRSD.clear();
        RadiTracerParticleEnergyDataString.clear();
        RadiTracerDataForTotalDoseString.clear();
        OrgansNameVector.clear();
        VolumeNameMassMap.clear();
        VolumeNameVolumeMap.clear();
        VolumeNameDensityMap.clear();
        // ////////////////////////////////////////////////////////////
    }
    // ////////////////////////////////////////////////////

    G4cout << "\n-------- " << ErrorMessage << G4endl ;


}

// are called from constructor to read the data
void G4TResultCalculation::ReadSimulationData(){
    
    if(V)G4cout << "\n========= " << __FUNCTION__ << G4endl ;
    
    std::ostringstream filename1;
    std::ifstream file1(MacrosStartingFile , std::ios::binary);
    
    G4String line, word, ParameterName;

    if(file1.is_open()){
        
        double mass, volume, density;
        int first = 0;
        
        bool readgeomdata = false;
        
        if(readgeomdata){
            
            while (getline(file1, line)) {
                
                std::istringstream LineString(line);
                LineString >> ParameterName;
                
                if(ParameterName.empty() || ParameterName == ""){
                    continue;
                }
                
                if(ParameterName == "#GeometryData"){
                    first++;
                    if(first == 2){
                        break;
                    }
                    continue;
                }
                
                if(first == 1){
                    LineString >> mass >> density >> volume ;
                    
                    OrgansNameVector.push_back(ParameterName);
                    VolumeNameMassMap[ParameterName] = mass ;
                    VolumeNameVolumeMap[ParameterName] = density ;
                    VolumeNameDensityMap[ParameterName] = volume ;
                    
                    G4cout << " ParameterName " << ParameterName <<  " mass " << mass <<  " density " << density <<  " volume " << volume << G4endl ;
                    continue;
                }
            }
        }
        if(V)G4cout << "\nReading file " << filename1.str().c_str() << G4endl ;
        
        while (getline(file1, line)) {
            
            std::istringstream LineString(line);
            ParameterName = "";
            LineString >> ParameterName;
            //if(V)G4cout << "\nParameterName " << ParameterName << G4endl ;
            
            if(ParameterName =="/GeometryData/setGeometrySymbol"){// to evaluate the last 4 lines
                LineString >> GeometrySymbol ;
                if(V)G4cout << " GeometrySymbol " << GeometrySymbol << G4endl ;
            }
            else if(ParameterName =="/PhysicsData/setPhysicsData"){// to evaluate the last 4 lines
                LineString >> Physics ;
                if(V)G4cout << " Physics " << Physics << G4endl ;
            }
            else if(ParameterName =="/PhysicsData/setCutsData"){// to evaluate the last 4 lines
                LineString >> CutsEnergy  >> CutsDistance ;
                if(V)G4cout << " CutsEnergy " << CutsEnergy << " CutsDistance " << CutsDistance << G4endl ;
            }
            else if(ParameterName =="/SourceData/setEventsParticleNameData"){// to evaluate the last 4 lines
                LineString >> ParticleName ;
                if(V)G4cout << " ParticleName " << ParticleName << G4endl ;
            }
            else if(ParameterName == "/SourceData/setEventsInitialMomDirData"){// to evaluate the last 4 lines
                LineString >> ParameterName >> MomDirDistribution;
                if(V)G4cout << " MomDirDistribution " << MomDirDistribution << G4endl ;
            }
            else if(ParameterName =="/SourceData/setEventsInitialEneData"){// to evaluate the last 4 lines
                LineString >> ParameterName >> EnergyDistribution >> ParticleSourceEnergy ;
                if(V)G4cout << " EnergyDistribution " << EnergyDistribution << " ParticleSourceEnergy " << ParticleSourceEnergy << G4endl ;
            }
            else if(ParameterName =="/SourceData/setEventsInitialPosData"){// to evaluate the last 4 lines
                LineString >> ParameterName >> SourceType >>  SourceRegionName ;
                if(V)G4cout << " SourceType " << SourceType << " SourceRegionName " << SourceRegionName << G4endl ;
            }
            else if(ParameterName =="/RunAndScoreData/setSimNumOnRanks"){// to evaluate the last 4 lines
                LineString >> OneOrMultiSimulations ;
                if(V)G4cout << " OneOrMultiSimulations " << OneOrMultiSimulations << G4endl ;
            }
            else if(ParameterName =="/RunAndScoreData/setQuantitiesToScore"){
                
                G4String wordsWithSpace ;
                while(LineString >> wordsWithSpace ){
                    if(wordsWithSpace.empty() || wordsWithSpace == ""){ }
                    else if (wordsWithSpace.find("All") != std::string::npos || wordsWithSpace.find("all") != std::string::npos) {
                        //if(V)G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n       found       --- " << Organs_To_Score << G4endl;
                        QuantityNamesToScore.push_back("AE");
                        QuantityNamesToScore.push_back("SAF");
                        QuantityNamesToScore.push_back("AF");
                        QuantityNamesToScore.push_back("AD");
                        QuantityNamesToScore.push_back("S");
                        QuantityNamesToScore.push_back("H");
                        QuantityNamesToScore.push_back("E");
                        QuantityNamesToScore.push_back("DR");
                        QuantityNamesToScore.push_back("DCC");
                        break;
                    }
                    else {
                        QuantityNamesToScore.push_back(wordsWithSpace);
                    }
                }
                if(V)G4cout << "  Quantities names to score are : " << G4endl ;
                for (G4int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {
                    if(V)G4cout << QuantityNamesToScore[gg] << G4endl;
                }
            }
            else if(ParameterName =="/RunAndScoreData/setVolumesToScore"){
                G4String wordsWithSpace ;

                IsAllTargetsToScore = false;
                IsAllSourcesToScore = false;

                G4String SourceOrTarget = "";
                while(LineString >> wordsWithSpace ){
                    if(wordsWithSpace == "source"){
                        SourceOrTarget = "source";
                    }
                    else if (wordsWithSpace == "target"){
                        SourceOrTarget = "target";
                    }
                    else if(wordsWithSpace == "All" || wordsWithSpace == "all" || wordsWithSpace.empty() || wordsWithSpace == "") {
                        if(SourceOrTarget == "source"){
                            IsAllSourcesToScore = true;
                        }
                        else if(SourceOrTarget == "target"){
                            IsAllTargetsToScore = true;
                        }
                        else{
                            IsAllTargetsToScore = true;
                            IsAllSourcesToScore = true;
                        }
                    }
                    else{
                        if (wordsWithSpace.find(":") != std::string::npos) {
                            if(SourceOrTarget == "source"){
                                IsAllSourcesToScore = false;
                            }
                            else if(SourceOrTarget == "target"){
                                IsAllTargetsToScore = false;
                            }
                            std::string str1 = wordsWithSpace;
                            std::stringstream ss (str1);
                            std::string item;
                            int ii = 0;
                            G4String RegionName;
                            while (getline (ss, item, ':')) {

                                if(ii == 0){
                                    RegionName = item;

                                    //G4cout << SourceOrTarget << " RegionName " << RegionName << G4endl ;

                                    if(SourceOrTarget == "source"){
                                        SourceNamesToScore.push_back(RegionName);
                                    }
                                    else if(SourceOrTarget == "target"){
                                        TargetNamesToScore.push_back(RegionName);
                                    }
                                }else{

                                    if (item.find("@") != std::string::npos) {
                                        std::string str2 = item;
                                        std::stringstream pp (str2);
                                        std::string subitem;
                                        int jj = 0;
                                        std::string regionname;
                                        G4double Fraction = 1;
                                        while (getline (pp, subitem, '@')) {

                                            if(jj == 0){
                                                regionname = subitem;
                                                bool IsIn = false;
                                                if(SourceOrTarget == "source"){
                                                    for (int gg = 0 ; gg < RegionsVolumesNamesMap[RegionName].size() ; gg++) { if(RegionsVolumesNamesMap[RegionName][gg] == subitem){ IsIn = true; }}
                                                    if(IsIn == false){ RegionsVolumesNamesMap[RegionName].push_back(subitem);}
                                                    NeWCombinationsForSourceNamesToScore[RegionName].push_back(subitem);
                                                    OldTargetNamesToScore.push_back(subitem);
                                                }
                                                else if(SourceOrTarget == "target"){
                                                    for (int gg = 0 ; gg < RegionsVolumesNamesMap[RegionName].size() ; gg++) { if(RegionsVolumesNamesMap[RegionName][gg] == subitem){ IsIn = true; }}
                                                    if(IsIn == false){
                                                        RegionsVolumesNamesMap[RegionName].push_back(subitem);
                                                        OldTargetNamesToScore.push_back(subitem);
                                                    }
                                                }
                                                G4cout << " SourceOrTarget " << SourceOrTarget
                                                       << " jj " << jj
                                                       << " RegionName " << RegionName
                                                       << " regionname " << regionname
                                                       << G4endl ;

                                            }else{

                                                Fraction = atof(subitem.c_str());

                                                if(SourceOrTarget == "source"){
                                                    RegionRegionsVolumesNamesFractionMap[RegionName][regionname] = Fraction;
                                                }
                                                else if(SourceOrTarget == "target"){
                                                    RegionRegionsVolumesNamesFractionMap[RegionName][regionname] = Fraction;
                                                }
                                                G4cout << " SourceOrTarget " << SourceOrTarget
                                                       << " jj " << jj
                                                       << " RegionName " << RegionName
                                                       << " regionname " << regionname
                                                       << " Fraction " << Fraction
                                                       << G4endl ;

                                            }



                                            jj++;
                                        }

                                    }else{


                                        bool IsIn = false;
                                        if(SourceOrTarget == "source"){
                                            for (int gg = 0 ; gg < RegionsVolumesNamesMap[RegionName].size() ; gg++) { if(RegionsVolumesNamesMap[RegionName][gg] == item){ IsIn = true; }}
                                            if(IsIn == false){ RegionsVolumesNamesMap[RegionName].push_back(item);}
                                            NeWCombinationsForSourceNamesToScore[RegionName].push_back(item);
                                            RegionRegionsVolumesNamesFractionMap[RegionName][item] = 1;
                                            OldTargetNamesToScore.push_back(item);
                                        }
                                        else if(SourceOrTarget == "target"){
                                            for (int gg = 0 ; gg < RegionsVolumesNamesMap[RegionName].size() ; gg++) { if(RegionsVolumesNamesMap[RegionName][gg] == item){ IsIn = true; }}
                                            if(IsIn == false){
                                                RegionsVolumesNamesMap[RegionName].push_back(item);
                                                RegionRegionsVolumesNamesFractionMap[RegionName][item] = 1;
                                                OldTargetNamesToScore.push_back(item);
                                            }
                                        }
                                    }

                                }
                                ii++;
                            }
                        }
                        else {
                            if(SourceOrTarget == "source"){
                                IsAllSourcesToScore = false;
                                SourceNamesToScore.push_back(wordsWithSpace);
                            }
                            else if(SourceOrTarget == "target"){
                                IsAllTargetsToScore = false;
                                TargetNamesToScore.push_back(wordsWithSpace);
                            }

                        }
                    }
                }

                if(V)G4cout << "  Source region names to score are : " << G4endl ;
                for (G4int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                    if(V)G4cout << SourceNamesToScore[gg] << G4endl;
                }

                if(V)G4cout << "  Target region names to score are : " << G4endl ;
                for (G4int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                    if(V)G4cout << TargetNamesToScore[gg] << G4endl;
                }

                if(V)G4cout << " IsAllSourcesToScore " << IsAllSourcesToScore << G4endl;
                if(V)G4cout << " IsAllTargetsToScore " << IsAllTargetsToScore << G4endl;

            }
            else if(ParameterName =="/RunAndScoreData/setResultDirectoryPath"){
                LineString >> ResultDirectoryPath;
                if(!DirectoryExists(ResultDirectoryPath)){
                    ResultDirectoryPath = appBuildDir+"/Results";
                }
                if(V)G4cout << " ResultDirectoryPath " << ResultDirectoryPath << G4endl ;
            }
            else if(ParameterName =="/RunAndScoreData/setRadioTracerData"){// to evaluate the last 4 lines
                
                G4String filename, type;

                LineString >> RadioTracerName ;
                LineString >> type ;
                if(type == "File" || type == "ICRP107"){
                    
                    RadiotracerDataFomFile = true;
                    //LineString >> RadioTracerName ;
                    
                    LineString >> filename ;
                    
                    std::ifstream file1(filename.c_str() , std::ios::binary);

                    if(file1.is_open()){
                        
                        G4String RPar;
                        G4String S;
                        G4double total = 0;
                        G4double LastEnergy = 0.;
                        std::string line, indicator;
                        G4double EVal;
                        G4double PVal;
                        while (getline(file1, line)) {

                            //G4cout << " the line " << line << G4"\n" ;

                            std::istringstream A(line);

                            if(A.str().empty()){
                                continue;
                            }

                            //QTextStream(stdout) << word.c_str() << " "<< Mass1 <<  " " << Mass2 <<"\n";

                            A >> word ;
                            //std::cout << " Line " << A.str() << "             word = " << word << std::endl;

                            if(word == "#"){
                                continue;
                            }
                            else if(word.empty()){
                                continue;
                            }
                            else if(word == "Radionuclide"){
                                indicator = "Radionuclide";
                                A >> RadioTracerName;
                                continue;
                            }
                            else if(   word == "gamma"
                                    || word == "e-"
                                    || word == "e+"
                                    || word == "alpha"
                                    || word == "proton"
                                    || word == "neutron"){

                                RPar = word;
                                A >> S;
                                LastEnergy = 0;
                                total = 0;

                                continue;

                            }else{

                                if(type == "File"){
                                    EVal = atof(word);
                                    A >> PVal;
                                    total += PVal;
                                    if(S == "Spectrum"){
                                        RadioTracerEnergyPerCent[RadioTracerName][RPar][EVal] = (EVal-LastEnergy)*PVal*100;
                                    }else{
                                        RadioTracerEnergyPerCent[RadioTracerName][RPar][EVal] = PVal*100;
                                    }
                                }
                                else if (type == "ICRP107"){
                                    A >> PVal;
                                    A >> EVal;
                                    A >> word;
                                    //total += PVal;
                                    if(word == "G")      {RPar = "gamma";   S == "Discrete";}
                                    else if(word == "PG"){RPar = "gamma";   S == "Discrete";}
                                    else if(word == "DG"){RPar = "gamma";   S == "Discrete";}
                                    else if(word == "X") {RPar = "gamma";   S == "Discrete";}
                                    else if(word == "AQ"){RPar = "gamma";   S == "Discrete";}
                                    else if(word == "B+"){RPar = "e+";      S == "Discrete";}
                                    else if(word == "B-"){RPar = "e-";      S == "Discrete";}
                                    else if(word == "BD"){RPar = "e-";      S == "Discrete";}
                                    else if(word == "IE"){RPar = "e-";      S == "Discrete";}
                                    else if(word == "AE"){RPar = "e-";      S == "Discrete";}
                                    else if(word == "A") {RPar = "alpha";   S == "Discrete";}
                                    else if(word == "AR"){RPar = "alpha";   S == "Discrete";}
                                    else if(word == "FF"){RPar = "FF";      S == "Discrete";}
                                    else if(word == "N") {RPar = "neutron"; S == "Discrete";}

                                    if(RPar == "e+"){
                                        RadioTracerEnergyPerCent[RadioTracerName]["e-"][EVal] = PVal*100;
                                        RadioTracerEnergyPerCent[RadioTracerName]["gamma"][0.511] = 2*PVal*100;
                                    }else{
                                        RadioTracerEnergyPerCent[RadioTracerName][RPar][EVal] = PVal*100;
                                    }
                                }


                                std::cout << " RadioTracerName " << RadioTracerName
                                          << " RPar = " << RPar
                                          << " S = " << S
                                          << " LastEnergy = " << LastEnergy
                                          << " EVal = " << EVal
                                          << " PVal = " << PVal
                                          << " Yield% = " << RadioTracerEnergyPerCent[RadioTracerName][RPar][EVal] << std::endl;

                                LastEnergy = EVal;

                                continue;
                            }
                        }

                        file1.close();
                    }
                }
                else{
                    if(V)G4cout << " RadioTracerName " << RadioTracerName << G4endl ;
                    G4double PerC, Ene; G4String pn;
                    while(LineString >> pn ){
                        LineString >> Ene;
                        LineString >> PerC;
                        RadioTracerEnergyPerCent[RadioTracerName][pn][Ene] = PerC;
                        if(V)G4cout << " pn " << pn << " Ene " << Ene << " PerC " << PerC << G4endl ;
                    }
                }
                GenerateResultsForRadioTracer = true;
            }
            else if(ParameterName =="/RunAndScoreData/setRadioTracerBiokinetic"){// to evaluate the last 4 lines
                
                LineString >> RadioTracerName >> InjectedActivity;
                LineString >> AdministeredActivityUnit >> ResidenceTimeUnit;
                
                if(AdministeredActivityUnit == "kBq"){
                    AdministeredActivityUnitFactor = 1e+3;
                }else if(AdministeredActivityUnit == "MBq"){
                    AdministeredActivityUnitFactor = 1e+6;
                }else if(AdministeredActivityUnit == "GBq"){
                    AdministeredActivityUnitFactor = 1e+9;
                }else {
                    AdministeredActivityUnit = "Bq";
                    AdministeredActivityUnitFactor = 1.;
                }
                
                std::stringstream unittext ;
                unittext<<InjectedActivity<<AdministeredActivityUnit;
                UnitPerRadioTracerDecay = unittext.str().c_str();

                if(ResidenceTimeUnit == "min"){
                    ResidenceTimeUnitFactor = 60;
                }else if(ResidenceTimeUnit == "h"){
                    ResidenceTimeUnitFactor = 3600;
                }else if(ResidenceTimeUnit == "d"){
                    ResidenceTimeUnitFactor = 86400;
                }else if(ResidenceTimeUnit == "y"){
                    ResidenceTimeUnitFactor = 31536000;
                }else{
                    ResidenceTimeUnit = "s";
                    ResidenceTimeUnitFactor = 1.;
                }
                
                if(V)G4cout << " AdministeredActivityUnit " << AdministeredActivityUnit
                            << " AdministeredActivityUnitFactor " << AdministeredActivityUnitFactor
                            << " ResidenceTimeUnit " << ResidenceTimeUnit
                            << " ResidenceTimeUnitFactor " << ResidenceTimeUnitFactor
                            << G4endl ;
                
                InjectedActivity = InjectedActivity*AdministeredActivityUnitFactor;
                RadioTracerInjectedActivity[RadioTracerName] = InjectedActivity;
                if(V)G4cout << " RadioTracerName " << RadioTracerName << " InjectedActivity " << InjectedActivity << G4endl ;
                
                G4double Fs, Ti, ai, AsPerA0 ; G4String orgName;
                while(LineString >> orgName ){
                    
                    //LineString >> Fs;RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracerName][orgName].push_back(Fs);
                    //LineString >> Ti;RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracerName][orgName].push_back(Ti);
                    //LineString >> ai;RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracerName][orgName].push_back(ai);
                    LineString >> AsPerA0 ;
                    RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracerName][orgName]=AsPerA0*ResidenceTimeUnitFactor;
                    //" Fs " << Fs << " Ti " << Ti  << " ai " << ai  <<
                    if(V)G4cout << " orgName " << orgName << " AsPerA0 " << AsPerA0 << G4endl ;
                }
                GenerateResultsForRadioTracerExams = true;
                
            }
            else if(ParameterName =="/RunAndScoreData/setRadiationFactors"){// to evaluate the last 4 lines
                //G4double Ene, Fac; G4String pn;
                //while(LineString >> pn ){
                //LineString >> Ene;
                //LineString >> Fac;
                //RadiationFactorMap[pn][Ene] = Fac;
                //if(V)G4cout << " pn " << pn << " Ene " << Ene << " Fac " << Fac << G4endl ;
                //}
            }
            else if(ParameterName =="/RunAndScoreData/setTissueFactors"){// to evaluate the last 4 lines
                G4double Fac; G4String pn;
                
                while(LineString >> pn ){
                    LineString >> Fac;
                    TissueFactorMap[pn] = Fac;
                    if(V)G4cout << " pn " << pn << " Fac " << Fac << G4endl ;
                }
            }
            else if(ParameterName =="/RunAndScoreData/setQuantitiesUnits"){// to evaluate the last 4 lines
                G4String Quant; G4String uni;
                
                while(LineString >> Quant ){
                    
                    LineString >> uni;
                    
                    if(Quant == "AE"){
                        if(uni == "MeV"){
                            AEUnit = "MeV";
                            AEUnitFactor = 1;
                        }else if(uni == "J"){
                            AEUnit = "J";
                            AEUnitFactor = 1.60218e-13;
                        }
                        if(V)G4cout << " Quantity " << Quant << " unit " << uni << " AEUnitFactor " << AEUnitFactor << G4endl ;
                    }else if(Quant == "AF"){
                        if(uni == ""){
                            AFUnit = "";
                            AFUnitFactor = 1;
                        }
                        if(V)G4cout << " Quantity " << Quant << " unit " << uni << " AFUnitFactor " << AFUnitFactor << G4endl ;
                    }else if(Quant == "SAF"){
                        if(uni == "kg-1"){
                            SAFUnit = "kg-1";
                            SAFUnitFactor = 1;
                        }else if(uni == "g-1"){
                            SAFUnit = "g-1";
                            SAFUnitFactor = 1/1e+3;
                        }
                        if(V)G4cout << " Quantity " << Quant << " unit " << uni << " SAFUnitFactor " << SAFUnitFactor << G4endl ;
                    }else if(Quant == "AD"){
                        if(uni == "MeV/kg"){
                            ADUnit = "MeV/kg";
                            ADUnitFactor = 1;
                        }else if(uni == "Gy"){
                            ADUnit = "Gy";
                            ADUnitFactor = MeV_to_J;
                        }
                        else if(uni == "miGy"){
                            ADUnit = "miGy";
                            ADUnitFactor = MeV_to_J*1e+6;
                        }else if(uni == "nGy"){
                            ADUnit = "nGy";
                            ADUnitFactor = MeV_to_J*1e+9;
                        }else if(uni == "mGy"){
                            ADUnit = "mGy";
                            ADUnitFactor = MeV_to_J*1e+3;
                        }else if(uni == "kGy"){
                            ADUnit = "kGy";
                            ADUnitFactor = MeV_to_J*1e-3;
                        }else if(uni == "MGy"){
                            ADUnit = "MGy";
                            ADUnitFactor = MeV_to_J*1e-6;
                        }
                        if(V)G4cout << " Quantity " << Quant << " unit " << uni << " ADUnitFactor " << ADUnitFactor << G4endl ;
                    }else if(Quant == "S"){
                        if(uni == "MeV/kg"){
                            SUnit = "MeV/kg";
                            SUnitFactor = 1;
                        }else if(uni == "Gy"){
                            SUnit = "Gy";
                            SUnitFactor = MeV_to_J;
                        }else if(uni == "miGy"){
                            SUnit = "miGy";
                            SUnitFactor = MeV_to_J*1e+6;
                        }else if(uni == "nGy"){
                            SUnit = "nGy";
                            SUnitFactor = MeV_to_J*1e+9;
                        }else if(uni == "mGy"){
                            SUnit = "mGy";
                            SUnitFactor = MeV_to_J*1e+3;
                        }else if(uni == "MGy"){
                            SUnit = "MGy";
                            SUnitFactor = MeV_to_J*1e-6;
                        }else if(uni == "kGy"){
                            SUnit = "kGy";
                            SUnitFactor = MeV_to_J*1e-3;
                        }else if(uni == "mGy/MBq"){
                            SUnit = "mGy";
                            SUnitFactor = MeV_to_J*1e+3/Bq_to_MBq;
                            UnitPerRadionuclideDecay = "MBq";
                        }
                        if(V)G4cout << " Quantity " << Quant << " unit " << uni << " SUnitFactor " << SUnitFactor << G4endl ;
                    }else if(Quant == "H"){
                        if(uni == "MeV/kg"){
                            HUnit = "MeV/kg";
                            HUnitFactor = 1;
                        }else if(uni == "Gy"){
                            HUnit = "Gy";
                            HUnitFactor = MeV_to_J;
                        }else if(uni == "miGy"){
                            HUnit = "miGy";
                            HUnitFactor = MeV_to_J*1e+6;
                        }else if(uni == "nGy"){
                            HUnit = "nGy";
                            HUnitFactor = MeV_to_J*1e+9;
                        }else if(uni == "mGy"){
                            HUnit = "mGy";
                            HUnitFactor = MeV_to_J*1e+3;
                        }else if(uni == "MGy"){
                            HUnit = "MGy";
                            HUnitFactor = MeV_to_J*1e-6;
                        }else if(uni == "kGy"){
                            HUnit = "kGy";
                            HUnitFactor = MeV_to_J*1e-3;
                        }else if(uni == "Sv"){
                            HUnit = "Sv";
                            HUnitFactor = MeV_to_J*Gy_to_Sv;
                        }else if(uni == "mSv"){
                            HUnit = "mSv";
                            HUnitFactor = MeV_to_J*1e+3*Gy_to_Sv;
                        }else if(uni == "mGy/MBq"){
                            SUnit = "mGy";
                            SUnitFactor = MeV_to_J*1e+3/Bq_to_MBq;
                            UnitPerRadionuclideDecay = "MBq";
                        }else if(uni == "mSv/MBq"){
                            SUnit = "mSv";
                            SUnitFactor = MeV_to_J*1e+3*Gy_to_Sv/Bq_to_MBq;
                            UnitPerRadionuclideDecay = "MBq";
                        }
                        if(V)G4cout << " Quantity " << Quant << " unit " << uni << " HUnitFactor " << HUnitFactor << G4endl ;
                    }else if(Quant == "E"){
                        if(uni == "MeV/kg"){
                            EUnit = "MeV/kg";
                            EUnitFactor = 1;
                        }else if(uni == "Gy"){
                            EUnit = "Gy";
                            EUnitFactor = MeV_to_J;
                        }else if(uni == "miGy"){
                            EUnit = "miGy";
                            EUnitFactor = MeV_to_J*1e+6;
                        }else if(uni == "nGy"){
                            EUnit = "nGy";
                            EUnitFactor = MeV_to_J*1e+9;
                        }else if(uni == "mGy"){
                            EUnit = "mGy";
                            EUnitFactor = MeV_to_J*1e+3;
                        }else if(uni == "Sv"){
                            EUnit = "Sv";
                            EUnitFactor = MeV_to_J*Gy_to_Sv;
                        }else if(uni == "mSv"){
                            EUnit = "mSv";
                            EUnitFactor = MeV_to_J*1e+3*Gy_to_Sv;
                        }else if(uni == "MGy"){
                            EUnit = "MGy";
                            EUnitFactor = MeV_to_J*1e-6;
                        }else if(uni == "kGy"){
                            EUnit = "kGy";
                            EUnitFactor = MeV_to_J*1e-3;
                        }else if(uni == "mGy/MBq"){
                            SUnit = "mGy";
                            SUnitFactor = MeV_to_J*1e+3/Bq_to_MBq;
                            UnitPerRadionuclideDecay = "MBq";
                        }else if(uni == "mSv/MBq"){
                            SUnit = "mSv";
                            SUnitFactor = MeV_to_J*1e+3*Gy_to_Sv/Bq_to_MBq;
                            UnitPerRadionuclideDecay = "MBq";
                        }
                        if(V)G4cout << " Quantity " << Quant << " unit " << uni << " EUnitFactor " << EUnitFactor << G4endl ;
                    }
                    else if(Quant == "DCC"){
                        if(uni == "pGy/cm2"){
                            DCCUnit = "pGy/cm2";
                            DCCUnitFactor = MeV_to_J*1e+12 /*0.16021773*/;
                        }else {
                            DCCUnit = "pGy/cm2";
                            DCCUnitFactor = MeV_to_J*1e+12 /*0.16021773*/;
                        }
                        if(V)G4cout << " Quantity " << Quant << " unit " << uni << " SAFUnitFactor " << SAFUnitFactor << G4endl ;
                    }
                }
            }
            else if(ParameterName =="/GeometryData/setVoxelsData"){
                LineString >> VoxZNumber >> VoxYNumber >> VoxZNumber >> word >> word >> VoxXHalfSize >> VoxYHalfSize >> VoxZHalfSize ;
                if(V)G4cout << " VoxXNumber " << VoxXNumber << " VoxYNumber " << VoxYNumber << " VoxZNumber " << VoxZNumber << " VoxXHalfSize " << VoxXHalfSize << " VoxYHalfSize " << VoxYHalfSize << " VoxZHalfSize " << VoxZHalfSize << G4endl ;
            }
            else if(ParameterName =="/RunAndScoreData/generateVoxelsResults"){
                GenerateVoxelsResuls = true;
                if(V)G4cout << " GenerateVoxelsResults " << GenerateVoxelsResuls << G4endl ;
            }
            else if(ParameterName =="ExecutionMode"){
                LineString >> ExecutionMode;
                if(V)G4cout << " ExecutionMode " << ExecutionMode << G4endl ;
            }
        }
        file1.close();


        if(TissueFactorMap.size()==0){
            G4String FilePath = "../PackagesAndFiles/ICRPDATA/ICRP110RegionsData";
            std::ifstream file(FilePath.c_str() , std::ios::binary);

            TissueFactorMap.clear();

            if(file.is_open()){

                double WT, frac, frac1, Mass1 , Mass2;

                std::string line, indicator, organ, word;

                while (getline(file, line)) {

                    //G4cout << " the line " << line << G4"\n" ;

                    std::istringstream A(line);

                    if(A.str().empty()){
                        continue;
                    }

                    //QTextStream(stdout) << word.c_str() << " "<< Mass1 <<  " " << Mass2 <<"\n";

                    A >> word ;
                    if(word == "#"){
                        continue;
                    }
                    else if(word == "Source"){
                        indicator = "Source";
                        continue;
                    } else if (word == "Target"){
                        indicator = "Target";
                        continue;
                    } else if (word == "WT-factor"){
                        indicator = "WT-factor";
                        continue;
                    } else if (word == "OtherTissues"){
                        //QTextStream(stdout) << " the word " << word.c_str() << "\n" ;
                        indicator = "OtherTissues";
                        continue;
                    } else if (word == "TotalBody"){
                        //QTextStream(stdout) << " the word " << word.c_str() << "\n" ;
                        indicator = "TotalBody";
                        continue;
                    } else if (word == "NewSourceRegions"){
                        //QTextStream(stdout) << " the word " << word.c_str() << "\n" ;
                        indicator = "NewSourceRegions";
                        continue;
                    } else{

                        if (indicator == "WT-factor"){

                            organ = word;
                            A >> WT ;
                            if(organ.c_str() == "#"){
                                continue;
                            }

                            //G4cout << organ.c_str() << " " << WT << " ";
                            TissueFactorMap[organ.c_str()] = WT;
                            while (A >> word && A >> frac){

                                TissueFactorMap[word.c_str()] = WT;
                                //G4cout << word.c_str() << " "<< frac <<  " ";
                            }
                        }
                    }
                }
            }
            else{
                G4cout << "cannot open the file " << FilePath.c_str() << " to read the tissue weighting factors"<< G4endl ;
            }
        }
    }
    else{
        G4cout << "cannot open the file " << filename1.str().c_str() << G4endl ;
    }
    
}

// for step and event level for Voxelized
void G4TResultCalculation::ReadThreadVoxelResultFile(G4String fm){
    
    //G4cout << "\n------------------------------- VOXELS RESULTS MERGING -------------------------------\n" << G4endl;
    //G4cout << "\n\n\n\n\n" << __FUNCTION__<< G4endl;
    
    // this are constants between all threads the all threads have the same values
    
    //G4cout << "\n-------- " << __FUNCTION__ << "\n" << G4endl ;
    G4double nnn = 0, bbb = 0;
    std::ifstream file(fm , std::ios::binary);
    if(file.is_open()){

        //if(V)G4cout << "\nReading file  " << fm << G4endl ;
        
        unsigned int LineNum, CN, NOFVal;
        G4double val, Z, Y, X, mass;
        G4String word, line;
        
        //for( unsigned int jj = 0 ; jj < TotVoxNum ; jj++ ){
        //VoxED_Total[jj] = 0.;
        //}
        
        file >> LineNum >> VoxXNumber >> VoxYNumber >> VoxZNumber >> VoxXHalfSize >> VoxYHalfSize >> VoxZHalfSize;

        if(GETVOXELDATA==false){
            InitializeVoxelizedData();
        }
        //if(V)
        //G4cout << " LineNum " << LineNum << G4endl ;
        
        for( unsigned int jj = 0 ; jj < LineNum ; jj++ ){
            
            //if (ExeFromMerge){
            file >> CN    ;  //CNID[CN] = CN;
            file >> X     ;  CopyNumberXPos[CN]=X;
            file >> Y     ;  CopyNumberYPos[CN]=Y;
            file >> Z     ;  CopyNumberZPos[CN]=Z;
            file >> word  ;  CopyNumberRegionNameMap[CN] = word;
            file >> val   ;  VoxED_Total[CN] = VoxED_Total[CN] + val;
            file >> val   ;  VoxED2_Total[CN] = VoxED2_Total[CN] + val;
            file >> NOFVal;  VoxNOfValues[CN] = VoxNOfValues[CN] + NOFVal;
            file >> mass  ;  CopyNumberMassSize[CN] = mass;
            file >> val   ;  VoxFluence[CN] = VoxFluence[CN] + val;

            //if(V)G4cout << CN << " Steps " << NOFVal << " TotalSteps " << VoxNOfValues[CN] << " Region="<<word << " AE="<< val << " AETotal="<< VoxED_Total[CN] << G4endl ;

            //}else{
            //    file >> CN >> X >> X >> X >> word  ;
            //    file >> val   ;  VoxED_Total[CN] = VoxED_Total[CN] + val;
            //    file >> val   ;  VoxED2_Total[CN] = VoxED2_Total[CN] + val;
            //    file >> NOFVal;  VoxNOfValues[CN] = VoxNOfValues[CN] + NOFVal;
            //    file >> X ;
            //    //if(V)G4cout << CN << " Steps " << NOFVal << " TotalSteps " << VoxNOfValues[CN] << " AE Total="<< VoxED_Total[CN] << G4endl ;
            //}
            //if(V)G4cout << CN << " Steps " << NOFVal << " TotalSteps " << VoxNOfValues[CN] << " Region="<<word << " AE="<< val << " AETotal="<< VoxED_Total[CN] << G4endl ;
            //if(V)G4cout << " " << VoxNOfValues[CN] << "\n" << G4endl ;
            
            //nnn += VoxNOfValues[CN];
            //bbb += VoxED_Total[CN];
            //if(V)G4cout << "jj " << jj << "\n" << G4endl ;
            
            //if(V)G4cout << CN << " " << Z << " " << Y << " " << X << " " << CopyNumberRegionNameMap[CN] << " " << VoxED_Total[CN] << " "<< VoxED2_Total[CN] << " " << VoxNOfValues[CN] << " " << CopyNumberMassSize[CN] << "\n" << G4endl ;
        }
        
        G4cout << " LineNum " << LineNum << G4endl ;

        file.close();
    }
    else{
        G4cout << "cannot open the file " << fm << G4endl ;
    }
    //if(V)G4cout << "Steps " << nnn << " ED " << bbb << G4endl ;
}
void G4TResultCalculation::VoxelQuantitiesCalculation(){
    
    if(V)G4cout << "\n-------- " << __FUNCTION__ << "\n" << G4endl;
    
    if(V)G4cout << "TotalEventNumber " << TotalEventNumber<< G4endl;
    if(V)G4cout << "TotalEmittedEnergy " << TotalEmittedEnergy<< G4endl;
    
    TotalNumberOfSteps = 0;
    TotalAbsorbedEnergy = 0.;
    
    TotVoxNum = VoxXNumber * VoxYNumber * VoxZNumber;

    //G4double jjj = 0, TotalMass = 0.;
    for(unsigned int gg = 0 ; gg < TotVoxNum ; gg++){
        
        // || isnanl(VoxED_Total[gg]) || isinfl(VoxED_Total[gg])
        if(VoxED_Total[gg] == 0 ){
            continue;
        }
        
        VoxAFCte[gg]  = 1./(double)TotalEmittedEnergy;
        VoxSAFCte[gg] = SAFUnitFactor/((double)TotalEmittedEnergy*CopyNumberMassSize[gg]);
        VoxADCte[gg]  = (ADUnitFactor/(CopyNumberMassSize[gg])) ;
        VoxSCte[gg]   = SUnitFactor/(double)TotalEventNumber;
        VoxHCte[gg]   = HUnitFactor*GenerateRadiationFactor(ParticleName,ParticleSourceEnergy);
        VoxECte[gg]   = EUnitFactor*(double)TissueFactorMap[CopyNumberRegionNameMap[gg]] ;
        VoxDCCCte[gg] = DCCUnitFactor/(CopyNumberMassSize[gg]*VoxFluence[gg]);
        //std::fprintf(stdout,"%-15u%-15u%-15e%-15e%-15e%-15e%-15e%-15e\n", gg , VoxNOfValues[gg], VoxAFCte[gg] , VoxSAFCte[gg] , VoxADCte[gg] , VoxSCte[gg], VoxHCte[gg] , VoxECte[gg]);
        
        //if(V)G4cout << "VoxHCte[gg]: " << VoxHCte[gg]<< G4endl;
        //if(V)G4cout << "VoxECte[gg]: " << VoxECte[gg]<< G4endl;
        
        TotalNumberOfSteps = TotalNumberOfSteps + VoxNOfValues[gg] ;
        TotalAbsorbedEnergy = TotalAbsorbedEnergy + VoxED_Total[gg];
        
        //if(V)G4cout << gg << " " << CopyNumberRegionNameMap[gg] << " " << VoxNOfValues[gg] << " TotalNumberOfSteps: " << TotalNumberOfSteps << " TotalAbsorbedEnergy " << TotalAbsorbedEnergy << G4endl;
        
        //if(V)G4cout << "VoxED_Total[gg]: " << VoxED_Total[gg] << " TotalAbsorbedEnergy: " << TotalAbsorbedEnergy << G4endl;
        
        VoxED_Mean[gg]  = VoxED_Total[gg]/VoxNOfValues[gg];
        VoxED2_Mean[gg] = VoxED2_Total[gg]/VoxNOfValues[gg];
        VoxED_Var[gg]   = (VoxED2_Mean[gg]-(VoxED_Mean[gg]*VoxED_Mean[gg]))/(VoxNOfValues[gg]-1);
        VoxED_SDev[gg]  = std::sqrt(std::abs(VoxED_Var[gg]));
        VoxED_RelS_D[gg]= ((VoxED_SDev[gg]/VoxED_Mean[gg])*100.);
        
        //if(V)G4cout << "VoxED_Total[gg] " << VoxED_Total[gg] << G4endl;
        
        VoxAF_Total[gg]  = VoxAFCte[gg]  * VoxED_Total[gg] ;
        VoxSAF_Total[gg] = VoxSAFCte[gg] * VoxED_Total[gg] ;
        VoxAD_Total[gg]  = VoxADCte[gg]  * VoxED_Total[gg] ;
        VoxS_Total[gg]   = VoxSCte[gg]   * VoxAD_Total[gg] ;
        VoxH_Total[gg]   = VoxHCte[gg]   * VoxAD_Total[gg] ;
        VoxE_Total[gg]   = VoxECte[gg]   * VoxH_Total[gg]  ;
        VoxDCC_Total[gg] = VoxDCCCte[gg] * VoxED_Total[gg] ;

    }
    //std::fprintf(stdout,"----------------------------------------------------------------------------------------------------------------------\n");
    
    
    if(V)G4cout << "TotalNumberOfSteps: " << TotalNumberOfSteps<< G4endl;
    if(V)G4cout << "TotalAbsorbedEnergy: " << TotalAbsorbedEnergy<< G4endl;
    if(V)G4cout << "PhantomEffectiveDose: " << PhantomEffectiveDose<< G4endl;
    
    if(V)printf("\n");
    
}
void G4TResultCalculation::GenerateVoxelsResultFiles(){
    
    if(V)G4cout << "-------- GenerateVoxelsResultFile\n"<< G4endl;
    
    VoxMaxRel_SDv=0; VoxMinRel_SDv=100;
    
    int WW = 15;
    
    //G4cout << " TotVoxNum "<< TotVoxNum << G4endl;

    //if(V)G4cout << "Total Voxels Number " << TotVoxNum << G4endl;
    
    for(unsigned int dd = 0 ; dd < TotVoxNum; dd++){
        
        //if(V)G4cout << VoxED_RelS_D[dd]<< G4endl;
        
        if(VoxED_Total[dd] == 0. || isnanl(VoxED_Total[dd]) || isinfl(VoxED_Total[dd])){
            continue;
        }
        
        if(isinfl(VoxED_RelS_D[dd]) || isnanl(VoxED_RelS_D[dd])){
            VoxED_RelS_D[dd] = 100;
            //if(V)G4cout << VoxED_RelS_D[dd] << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;
        }
        
        if(VoxED_RelS_D[dd] < VoxMinRel_SDv){
            VoxMinRel_SDv = VoxED_RelS_D[dd];
            VoxelCNForMinRel_SDv = dd;
        }
        if(VoxED_RelS_D[dd] > VoxMaxRel_SDv){
            VoxMaxRel_SDv = VoxED_RelS_D[dd];
            VoxelCNForMaxRel_SDv = dd;
        }
    }
    
    
    //std::ostringstream fname;
    //fname << ResultDirectoryPath << "/VoxelsResults";
    //std::ofstream file1(fname.str().c_str(), std::ios_base::binary);
    //if(file1.is_open()){
    
    std::ostringstream filname ;
    filname << ResultDirectoryPath << "/VoxelsResults";
    G4String fileNameString = filname.str().c_str();
    
    if(FILE* file1 = fopen(fileNameString.c_str(),"a")){
        //if(fileR.is_open()){
        
        
        std::stringstream aastr;
        aastr << "# "
              << "Source:"<< SourceRegionName << " "
              << ParticleName << " "
              << ParticleSourceEnergy << "MeV "
              << EnergyDistribution <<  " "
              << MomDirDistribution << " "
              << SourceType << " "
              << TotalEventNumber << "Event "
              << TotalNumberOfSteps << "Step "
              << CutsDistance << "mm "
              << CutsEnergy << "MeV "
              << GenerateRadiationFactor(ParticleName,ParticleSourceEnergy) <<"Wr "
              << ExecutionMode << " "
              << RankID << " "
              << TotalEmittedEnergy << "MeV "
              << TotalAbsorbedEnergy - ED_Total["World"] << "MeV "
              << PhantomEffectiveDose << "Sv "
              << VoxelCNForMinRel_SDv << " " << VoxMinRel_SDv << "% "
              << VoxelCNForMaxRel_SDv << " " << VoxMaxRel_SDv << "% "
              << (double)OneEventExecutionTimeInMs << "ms "
              << (double)MaxExecutionTimeInMin << "min "
              << (double)ExecutionTimeInMin << "min " << "\n\n";
        
        G4cout << "\n2) Calculated Voxelized Data(just header):\n\n" << aastr.str().c_str()<< G4endl;
        
        aastr << std::setw(WW) << std::left << "CN"
              << std::setw(WW) << std::left << "# XPos(mm)"
              << std::setw(WW) << std::left << "# YPos(mm)"
              << std::setw(WW) << std::left << "# ZPos(mm)"
              << std::setw(WW) << std::left << "# Mass[kg]"
              << std::setw(WW) << std::left << "# AE(MeV)"
              << std::setw(WW) << std::left << "# AD(MeV/kg)"
              << std::setw(WW) << std::left << "# AF"
              << std::setw(WW) << std::left << "# SAF(kg-1)"
              << std::setw(WW) << std::left << "# S(Gray/ev)"
              << std::setw(WW) << std::left << "# H(Sv)"
              << std::setw(WW) << std::left << "# E(Sv)"
              << std::setw(WW) << std::left << "# RelSDv(MeV)"
              << "\n\n";

        //file1 << aastr.str().c_str();
        //if(V)G4cout << aastr.str().c_str()<< G4endl;
        
        std::fprintf(file1,"%s",aastr.str().c_str());
        
        for(G4int dd = 0 ; dd < TotVoxNum; dd++){
            
            if(VoxED_Total[dd] == 0. || isnanl(VoxED_Total[dd]) || isinfl(VoxED_Total[dd])){
                continue;
            }
            
            std::stringstream table;
            
            //if(V)G4cout << "ED " << VoxED_Total[dd] << " " << CopyNumberMassSize[dd] << G4endl;
            
            table << std::setw(WW) << std::left << dd
                     << std::setw(WW) << std::left << CopyNumberXPos[dd]
                     << std::setw(WW) << std::left << CopyNumberYPos[dd]
                        << std::setw(WW) << std::left << CopyNumberZPos[dd]
                           << std::setw(WW) << std::left << CopyNumberMassSize[dd]
                              << std::setw(WW) << std::left << VoxED_Total[dd]
                                 << std::setw(WW) << std::left << VoxAD_Total[dd]
                                    << std::setw(WW) << std::left << VoxAF_Total[dd]
                                       << std::setw(WW) << std::left << VoxSAF_Total[dd]
                                          << std::setw(WW) << std::left << VoxS_Total[dd]
                                             << std::setw(WW) << std::left << VoxH_Total[dd]
                                                << std::setw(WW) << std::left << VoxE_Total[dd]
                                                   << std::setw(WW) << std::left << VoxED_RelS_D[dd]
                                                      << "\n";
            std::fprintf(file1,"%s",table.str().c_str());
            //file1 << table.str().c_str();
            //if(V)G4cout << table.str().c_str();
            
        }
        if(V)G4cout << "\n" << G4endl;
        std::fclose(file1);
        
        //if(V)G4cout << "The CN " << VoxelCNForMinRel_SDv << ", have Min Rel SDev " << (double)VoxMinRel_SDv << G4endl;
        //if(V)G4cout << "The CN " << VoxelCNForMaxRel_SDv << ", have Max Rel SDev " << (double)VoxMaxRel_SDv << G4endl;
    }
}

// for step and event level
bool G4TResultCalculation::ReadThreadRegionResultFile(G4String fm){
    
    //if(V) G4cout << "\n-------- " << __FUNCTION__ << G4endl ;

    bool IsRegionsAreRead = false;

    std::ifstream file(fm , std::ios::binary);
    if(file.is_open()){

        // ///////////////////////////////////////////////////////////////////////////////////
        //G4cout << "\n-------- is open " << G4endl ;

        G4double val = 0; unsigned long long int ival = 0;
        G4String line , ParameterName;
        while (getline(file, line)) {
            
            std::istringstream LineString(line);
            LineString >> ParameterName;

            //G4cout << LineString.str() << G4endl ;

            if(V)G4cout << "\nParameterName " << ParameterName << G4endl ;
            if(ParameterName.empty() || ParameterName == "RegionName" || ParameterName == "RegionsData"){ IsRegionsAreRead = true ; continue; }
            
            if(ParameterName =="TotalEventNumber"){
                LineString >> ival; TotalEventNumber += ival ;
                if(V)G4cout << " TotalEventNumber " << TotalEventNumber << G4endl ;
            }
            else if(ParameterName =="EnergyEmittedPerThread"){
                LineString >> EnergyEmittedPerThread ;
                if(V)G4cout << " EnergyEmittedPerThread " << EnergyEmittedPerThread << G4endl ;
                TotalEmittedEnergy += EnergyEmittedPerThread;
            }
            else if(ParameterName =="ExecutionTimeInMin"){
                
                LineString >> val;
                if(val > MaxExecutionTimeInMin){
                    MaxExecutionTimeInMin = val;
                }
                ExecutionTimeInMin += val ;
                
                if(V)G4cout << " ExecutionTimeInMin " << ExecutionTimeInMin << G4endl ;
            }
            else if(ParameterName =="OneEventExecutionTimeInMs"){
                LineString >> val; OneEventExecutionTimeInMs+=val ;
                if(V)G4cout << " OneEventExecutionTimeInMs " << OneEventExecutionTimeInMs << G4endl ;
            }
            else if(ParameterName =="ParticleName"){
                LineString >> ParticleName ;
                if(V)G4cout << " ParticleName " << ParticleName << G4endl ;
            }
            else if(ParameterName =="SourceType"){
                LineString >> SourceType ;
                if(V)G4cout << " SourceType " << SourceType << G4endl ;
            }
            else if(ParameterName =="SourceRegionName"){
                //LineString >> SourceRegionName ;
                //if(V)G4cout << " SourceRegionName " << SourceRegionName << G4endl ;
            }
            else if(ParameterName =="EnergyDistribution"){
                LineString >> EnergyDistribution ;
                if(V)G4cout << " EnergyDistribution " << EnergyDistribution << G4endl ;
            }
            else if(ParameterName =="ParticleSourceEnergy"){
                LineString >> ParticleSourceEnergy ;
                if(V)G4cout << " ParticleSourceEnergy " << ParticleSourceEnergy << G4endl ;
            }
            else if(ParameterName =="MomDirDistribution"){
                LineString >> MomDirDistribution ;
                if(V)G4cout << " MomDirDistribution " << MomDirDistribution << G4endl ;
            }
            else if(ParameterName =="RankID"){
                LineString >> RankID ;
                if(V)G4cout << " RankID " << RankID << G4endl ;
            }
            else{
                
                G4String VolumeName = ParameterName;
                if(VolumeName == "World" || VolumeName == "VOXEL"){continue;}

                LineString >> val ; ED_Total[VolumeName] += val ;
                LineString >> val ; ED2_Total[VolumeName] += val ;
                LineString >> ival ; NOfValues[VolumeName] += ival ;

                bool IsIn = false;
                for (int gg = 0 ; gg < OrgansNameVector.size() ; gg++) { if(OrgansNameVector[gg] == VolumeName){ IsIn = true; break; } }
                if(IsIn == false){
                    OrgansNameVector.push_back(VolumeName);
                    LineString >> val ; VolumeNameMassMap[VolumeName] = val ;
                    LineString >> val ; VolumeNameVolumeMap[VolumeName] = val ;
                    LineString >> val ; VolumeNameDensityMap[VolumeName] = val ;
                    LineString >> val ; Fluence[VolumeName] += val ;

                }else{
                    LineString >> val>> val>> val>> val ; Fluence[VolumeName] += val ;

                }
                if(TissueFactorMap[VolumeName] == 0.){
                    //TissueFactorMap[VolumeName] = 0.12;
                }

                //G4cout << LineString.str() << G4endl ;

                //G4cout << VolumeName << " " << ED_Total[VolumeName] <<  " " << ED2_Total[VolumeName] <<  " " << NOfValues[VolumeName] <<  " " << VolumeNameMassMap[VolumeName] <<  " " << VolumeNameVolumeMap[VolumeName] <<  " " << VolumeNameDensityMap[VolumeName] <<  " " << TissueFactorMap[VolumeName] << G4endl ;

            }
        }
        
        //if(GenerateRadiationFactor(ParticleName,ParticleSourceEnergy) == 0.){
        //  RadiationFactorMap[ParticleName][ParticleSourceEnergy] = 1;
        //}

        bool IsIn = false;
        for (int gg = 0 ; gg < SourceParticleEnergyValues[GeometrySymbol][SourceRegionName][ParticleName].size() ; gg++) {if(SourceParticleEnergyValues[GeometrySymbol][SourceRegionName][ParticleName][gg] == 0){IsIn = true;break;}}
        if(IsIn == false){SourceParticleEnergyValues[GeometrySymbol][SourceRegionName][ParticleName].push_back(0);}
        
        IsIn = false;
        for (int gg = 0 ; gg < SourceParticleEnergyValues[GeometrySymbol][SourceRegionName][ParticleName].size() ; gg++) {if(SourceParticleEnergyValues[GeometrySymbol][SourceRegionName][ParticleName][gg] == ParticleSourceEnergy){IsIn = true;break;}}
        if(IsIn == false){
            //G4cout << " Save Energy data " << ParticleSourceEnergy << G4endl ;
            SourceParticleEnergyValues[GeometrySymbol][SourceRegionName][ParticleName].push_back(ParticleSourceEnergy);
        }

        if(IsRegionsAreRead == false){
            return false;
        }else{
            ReadedResultFilesPaths[GeometrySymbol].push_back(fm);
            return true;
        }

        if(V)G4cout << G4endl ;
        
        file.close();
    }
    else{
        return false;
        G4cout << "cannot open the file " << fm << G4endl ;
    }
}
void G4TResultCalculation::RegionQuantitiesCalculation(){
    
    if(V)G4cout << "\n-------- " << __FUNCTION__ << G4endl ;

    for (int gg = 0 ; gg < OrgansNameVector.size() ; gg++) {
        //G4cout <<  " OrgansNameVector[gg] " << OrgansNameVector[gg] <<G4endl;
    }

    // add new combinated source target region data
    for ( auto it = RegionsVolumesNamesMap.begin(); it != RegionsVolumesNamesMap.end(); ++it  ){

        //G4cout <<  " it->first " << it->first << " " <<   VolumeNameMassMap[it->first] << " " << ED_Total[it->first] <<G4endl;

        // add region name, mass, density and volume to the defined regions
        bool IsIn = false;
        for (int gg = 0 ; gg < OrgansNameVector.size() ; gg++) {
            if(OrgansNameVector[gg] == it->first){ IsIn = true; break; }
        }
        if(IsIn == false){
            OrgansNameVector.push_back(it->first);

            VolumeNameMassMap[it->first] = 0. ;
            VolumeNameVolumeMap[it->first] = 0. ;
            VolumeNameDensityMap[it->first] = 0. ;
            ED_Total[it->first] = 0. ;
            ED2_Total[it->first] = 0. ;
            NOfValues[it->first] = 0. ;
            Fluence[it->first] = 0. ;
        }else{
            continue;
        }

        if(TissueFactorMap[it->first] == 0. || __isinf(TissueFactorMap[it->first]) || __isnan(TissueFactorMap[it->first])){
            //TissueFactorMap[it->first] = 0.12;
        }

        G4double multipliedmasses = 1;
        G4double addedmasses = 0;
        G4double multipliedmasses2 = 1;

        for (int gg = 0 ; gg < it->second.size() ; gg++) {

            //G4cout << it->second[gg] <<  " " <<   VolumeNameMassMap[it->second[gg]] <<  " " << ED_Total[it->second[gg]] << " it->first " <<  it->first << " " << VolumeNameMassMap[it->first]  <<G4endl;

            multipliedmasses = 1;
            addedmasses = 0;
            multipliedmasses2 = 1;

            bool IsIn = false;
            for (int zz = 0 ; zz < OrgansNameVector.size() ; zz++) {
                if(OrgansNameVector[zz] == it->second[gg]){ IsIn = true; break; }
            }


            if(IsIn == true){ //if exist means that its data existe

                //ED_Total[it->first] += ED_Total[it->second[gg]];
                //ED2_Total[it->first] += ED2_Total[it->second[gg]];

                std::vector<G4String> subregions = it->second;
                for (int AAA = 0 ; AAA < subregions.size() ; AAA++) {
                    if( VolumeNameMassMap[subregions[AAA]] != 0 && !__isnan(VolumeNameMassMap[subregions[AAA]]) && !__isinf(VolumeNameMassMap[subregions[AAA]])){
                        multipliedmasses *= VolumeNameMassMap[subregions[AAA]];
                        addedmasses += VolumeNameMassMap[subregions[AAA]];
                    }else{
                        ErrorMessage = ErrorMessage + "\n- The mass of "+ subregions[AAA] +" in " +it->first +" region is 0 \n";
                    }
                }

                multipliedmasses2 = multipliedmasses/VolumeNameMassMap[it->second[gg]];

                if( multipliedmasses2 != 0 && !__isnan(multipliedmasses2) && !__isinf(multipliedmasses2)){
                    if( multipliedmasses != 0 && !__isnan(multipliedmasses) && !__isinf(multipliedmasses)){
                        if( addedmasses != 0 && !__isnan(addedmasses) && !__isinf(addedmasses)){
                            ED_Total[it->first] += ((ED_Total[it->second[gg]]*multipliedmasses2)/multipliedmasses)
                                    *RegionRegionsVolumesNamesFractionMap[it->first][it->second[gg]]*addedmasses;
                            ED2_Total[it->first] += ((ED2_Total[it->second[gg]]*multipliedmasses2)/multipliedmasses)
                                    *RegionRegionsVolumesNamesFractionMap[it->first][it->second[gg]]*addedmasses;
                        }
                    }
                }


                //if(it->first == "CoW"){

                //}

                if( NOfValues[it->second[gg]] != 0 && !__isnan(NOfValues[it->second[gg]]) && !__isinf(NOfValues[it->second[gg]])){
                    NOfValues[it->first] += NOfValues[it->second[gg]];
                }

                if( Fluence[it->second[gg]] != 0 && !__isnan(Fluence[it->second[gg]]) && !__isinf(Fluence[it->second[gg]])){
                    Fluence[it->first] += Fluence[it->second[gg]];
                }

                if( VolumeNameMassMap[it->second[gg]] != 0 && !__isnan(VolumeNameMassMap[it->second[gg]]) && !__isinf(VolumeNameMassMap[it->second[gg]])){
                    VolumeNameMassMap[it->first] += VolumeNameMassMap[it->second[gg]] ;
                }

                if( VolumeNameDensityMap[it->second[gg]] != 0 && !__isnan(VolumeNameDensityMap[it->second[gg]]) && !__isinf(VolumeNameDensityMap[it->second[gg]])){
                    VolumeNameDensityMap[it->first] += VolumeNameDensityMap[it->second[gg]] ;
                }

                if( VolumeNameVolumeMap[it->second[gg]] != 0 && !__isnan(VolumeNameVolumeMap[it->second[gg]]) && !__isinf(VolumeNameVolumeMap[it->second[gg]])){
                    VolumeNameVolumeMap[it->first] += VolumeNameVolumeMap[it->second[gg]] ;
                }


                // to eliminate the initial regions data which affect data in DR and ER ratios which depends on the sum of dose
                //ED_Total[it->second[gg]] = 0;
                //ED2_Total[it->second[gg]] = 0;
                //NOfValues[it->second[gg]] = 0;
                //VolumeNameMassMap[it->second[gg]] = 0;
                //VolumeNameVolumeMap[it->second[gg]] = 0;
                //VolumeNameDensityMap[it->second[gg]] = 0;

            }
        }
        if( !__isnan(VolumeNameDensityMap[it->first]) && !__isinf(VolumeNameDensityMap[it->first])){
            //VolumeNameDensityMap[it->first] = VolumeNameDensityMap[it->first]/it->second.size();
            VolumeNameDensityMap[it->first] = 1000*(VolumeNameMassMap[it->first]/VolumeNameVolumeMap[it->first]);
        }
    }

    if(V)G4cout << " \n The " << ReadedResultFilesPaths[GeometrySymbol].size() << " result files readed are : \n"<<  G4endl ;
    for(G4int fge = 0 ; fge < ReadedResultFilesPaths[GeometrySymbol].size() ; fge++ ){
        if(V)G4cout << ReadedResultFilesPaths[GeometrySymbol][fge] << G4endl ;
    }
    
    if(V)G4cout << "\n========= Values needed to calculate the dosimetric quantities =========" << G4endl;
    
    ExecutionTimeInMin = ExecutionTimeInMin/ReadedResultFilesPaths[GeometrySymbol].size();
    OneEventExecutionTimeInMs = OneEventExecutionTimeInMs/ReadedResultFilesPaths[GeometrySymbol].size();
    
    if(V)G4cout << " OneEventExecutionTimeInMs " << OneEventExecutionTimeInMs << G4endl;
    if(V)G4cout << " MaxExecutionTimeInMin " << MaxExecutionTimeInMin << G4endl;
    if(V)G4cout << " ExecutionTimeInMin " << ExecutionTimeInMin << "\n" << G4endl;
    if(V)G4cout << "\n TotalNumberOfSteps " << TotalNumberOfSteps << G4endl;
    
    if(V)G4cout << "\n Constants values used in conversion between scored quantities" << G4endl;
    if(V)G4cout << "  - Number of events : " << TotalEventNumber << G4endl;
    if(V)G4cout << "  - ParticleSourceEnergy " << ParticleSourceEnergy << G4endl;
    if(V)G4cout << "  - Total emitted energy : " << TotalEmittedEnergy << G4endl;
    
    if(V)std::fprintf(stdout,"---------------------------------------------------------------------------------------------------------------------------------------------------\n");
    if(V)std::fprintf(stdout,"%-30s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n", "| Region Volume" , "| NOfValues" , "| ED_Total", "| ED_Mean" , "| ED2_Total" ,"| ED2_Mean", "| Variance", "| SDev", "| Rel_SD" );
    if(V)std::fprintf(stdout,"---------------------------------------------------------------------------------------------------------------------------------------------------\n");
    
    TotalNumberOfSteps = 0;
    TotalAbsorbedEnergy = 0.;
    PhantomEffectiveDose = 0.;
    TotalAbsorbedDose = 0;


    // to remove all old targets used to construct new target
    for (int zz = 0 ; zz < OldTargetNamesToScore.size() ; zz++) {
        //ED_Total[OldTargetNamesToScore[zz]] = NULL;
    }



    for ( auto it = ED_Total.begin(); it != ED_Total.end(); ++it  ){
        
        if(ED_Total[it->first] == 0. ){
            continue;
        }else {
            
            if(NOfValues[it->first] == 1. ){
                NOfValues[it->first] = 2.;
            }
            
            unsigned long long int ccc = NOfValues[it->first];
            if( __isnan(ccc) && __isinf(ccc) && ccc == 0 && ccc == NULL){
                G4cout << "\n\n !!!!!!!!!!!!!!!!!!!!!!!!!! " << it->first << " number of steps is not defined (" << NOfValues[it->first] << ")\n";
                continue;
            }

            double eee = ED_Total[it->first];
            if( __isnan(eee) && __isinf(eee) && eee == 0 && eee == NULL){
                G4cout << "\n\n !!!!!!!!!!!!!!!!!!!!!!!!!! " << it->first << " energy deposited is not defined (" << ED_Total[it->first] << ")\n";
                continue;
            }

            ED_Mean[it->first] = ED_Total[it->first]/NOfValues[it->first];
            ED2_Mean[it->first] = ED2_Total[it->first]/NOfValues[it->first];
            ED_Var[it->first] = (ED2_Mean[it->first]-(ED_Mean[it->first]*ED_Mean[it->first]))/(NOfValues[it->first]-1);
            ED_SDev[it->first] = std::sqrt(std::abs(ED_Var[it->first]));
            ED_RelS_D[it->first] = ((ED_SDev[it->first]/ED_Mean[it->first])*100.);
            TotalAbsorbedEnergy += ED_Total[it->first];

            //if(it->first == "CoW"){
            //    G4cout << it->first << " " << ParticleName <<" "<< ParticleSourceEnergy << " Total for 100*" << ED_Mean[it->first] << " " << ED_Var[it->first] << " " << ED_SDev[it->first] << "*" << NOfValues[it->first] << "/" << ED_Total[it->first] << "=" << ED_RelS_D[it->first] << G4endl;
            //}
        }
        
        G4String nn = it->first;nn.resize(29);
        if(V)std::fprintf(stdout,"%-30s%-15u%-15e%-15e%-15e%-15e%-15e%-15e%-15e\n", nn.c_str() , (int)NOfValues[it->first] , (double)ED_Total[it->first] , (double)ED_Mean[it->first], (double)ED2_Total[it->first] , (double)ED2_Mean[it->first], (double)ED_Var[it->first] , (double)ED_SDev[it->first] , (double)ED_RelS_D[it->first] ,"\n");

        if( !__isnan(NOfValues[it->first]) && !__isinf(NOfValues[it->first]) && NOfValues[it->first] != 0 && NOfValues[it->first] != NULL){
            TotalNumberOfSteps += NOfValues[it->first] ;
            //G4cout << TotalNumberOfSteps << G4endl;
        }

    }
    if(V)std::fprintf(stdout,"---------------------------------------------------------------------------------------------------------------------------------------------------\n");
    
    G4double A1 = 0;
    G4double A2 = 0;

    G4double jjj = 0. ;
    for(G4int gg = 0 ; gg < (G4int)OrgansNameVector.size() ; gg++){
        
        //E_Total[OrgansNameVector[gg]] = 0. ;// to initialize the value because if we let it gives -nan val which dont give error value in phantom effective dose
        if(__isnan(ED_Total[OrgansNameVector[gg]]) || __isinf(ED_Total[OrgansNameVector[gg]]) || ED_Total[OrgansNameVector[gg]] == 0 ){
            continue;
        }
        

        //This for just the error in mass, it should be removed after the calculation for article
        //if (VolumeNameMassMap["Liver"] > 1000){
        //    for (int zz = 0 ; zz < OrgansNameVector.size() ; zz++) {
        //        VolumeNameMassMap[OrgansNameVector[zz]] = VolumeNameMassMap[OrgansNameVector[zz]]/1000;
        //    }
        //    G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n " << G4endl ;
        //}

        AFCte[OrgansNameVector[gg]]  = AFUnitFactor/TotalEmittedEnergy;
        //G4cout << TotalEmittedEnergy << "  - AFCte[OrgansNameVector[gg]] " << AFCte[OrgansNameVector[gg]] << G4endl;
        SAFCte[OrgansNameVector[gg]] = SAFUnitFactor  /(VolumeNameMassMap[OrgansNameVector[gg]]*TotalEmittedEnergy);
        ADCte[OrgansNameVector[gg]]  = ADUnitFactor   /(VolumeNameMassMap[OrgansNameVector[gg]]*(double)TotalEventNumber);
        SCte[OrgansNameVector[gg]]   = SUnitFactor    /(VolumeNameMassMap[OrgansNameVector[gg]]*(double)TotalEventNumber);
        HCte[OrgansNameVector[gg]]   = HUnitFactor*(1./(VolumeNameMassMap[OrgansNameVector[gg]]*(double)TotalEventNumber))*GenerateRadiationFactor(ParticleName,ParticleSourceEnergy);
        ECte[OrgansNameVector[gg]]   = EUnitFactor*(1./(VolumeNameMassMap[OrgansNameVector[gg]]*(double)TotalEventNumber))*GenerateRadiationFactor(ParticleName,ParticleSourceEnergy)*TissueFactorMap[OrgansNameVector[gg]];
        DCCCte[OrgansNameVector[gg]] = DCCUnitFactor  /(VolumeNameMassMap[OrgansNameVector[gg]]*Fluence[OrgansNameVector[gg]]);

        AF_Total[OrgansNameVector[gg]]  = ED_Total[OrgansNameVector[gg]]*AFCte[OrgansNameVector[gg]];
        SAF_Total[OrgansNameVector[gg]] = ED_Total[OrgansNameVector[gg]]*SAFCte[OrgansNameVector[gg]];
        AD_Total[OrgansNameVector[gg]]  = ED_Total[OrgansNameVector[gg]]*ADCte[OrgansNameVector[gg]];
        S_Total[OrgansNameVector[gg]]   = ED_Total[OrgansNameVector[gg]]*SCte[OrgansNameVector[gg]];
        H_Total[OrgansNameVector[gg]]   = ED_Total[OrgansNameVector[gg]]*HCte[OrgansNameVector[gg]];
        E_Total[OrgansNameVector[gg]]   = ED_Total[OrgansNameVector[gg]]*ECte[OrgansNameVector[gg]];
        DCC_Total[OrgansNameVector[gg]] = ED_Total[OrgansNameVector[gg]]*DCCCte[OrgansNameVector[gg]];

        //G4cout << " OrgansNameVector[gg] " << OrgansNameVector[gg]
        //       << " DCCUnitFactor " << DCCUnitFactor
        //       << " Mass " << VolumeNameMassMap[OrgansNameVector[gg]]
        //       << " Fluence " << Fluence[OrgansNameVector[gg]]

        //       << " ED " << ED_Total[OrgansNameVector[gg]]
        //       << " DCC " << DCCCte[OrgansNameVector[gg]]
        //       << " DCC_Total " << DCC_Total[OrgansNameVector[gg]]
        //       <<G4endl;

        /*
                AF_Total[OrgansNameVector[gg]]  = AFCte[OrgansNameVector[gg]] * ED_Total[OrgansNameVector[gg]] ;
                SAF_Total[OrgansNameVector[gg]] = SAFCte[OrgansNameVector[gg]] * ED_Total[OrgansNameVector[gg]] ;
                AD_Total[OrgansNameVector[gg]]  = ADCte[OrgansNameVector[gg]] * ED_Total[OrgansNameVector[gg]] ;
                S_Total[OrgansNameVector[gg]]   = SCte[OrgansNameVector[gg]] * AD_Total[OrgansNameVector[gg]] ;
                H_Total[OrgansNameVector[gg]]   = HCte[OrgansNameVector[gg]] * AD_Total[OrgansNameVector[gg]] ;
                E_Total[OrgansNameVector[gg]]   = ECte[OrgansNameVector[gg]] * H_Total[OrgansNameVector[gg]];
        */
        
        // needed in calculation of Radio-tracer S-values
        ResultTable["AE" ][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = ED_Total [OrgansNameVector[gg]];
        ResultTable["AF" ][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = AF_Total [OrgansNameVector[gg]];
        ResultTable["SAF"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = SAF_Total[OrgansNameVector[gg]];
        ResultTable["AD" ][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = AD_Total [OrgansNameVector[gg]];
        ResultTable["S"  ][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = S_Total  [OrgansNameVector[gg]];
        ResultTable["H"  ][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = H_Total  [OrgansNameVector[gg]];
        ResultTable["E"  ][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = E_Total  [OrgansNameVector[gg]];
        ResultTable["DCC"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = DCC_Total[OrgansNameVector[gg]];

        //ResultTable["EAm"][ParticleName][ParticleSourceEnergy][SourceRegionName]["All"] += E_Total[OrgansNameVector[gg]];
        //ResultTable["FluenceAm"][ParticleName][ParticleSourceEnergy][SourceRegionName]["All"] += Fluence[OrgansNameVector[gg]];
        A1 += E_Total[OrgansNameVector[gg]];
        A2 += Fluence[OrgansNameVector[gg]]*TissueFactorMap[OrgansNameVector[gg]];
        if( __isnan(A2) && __isinf(A2) && A2 == 0 && A2 == NULL){}else{
            ResultTable["DCCAm"][ParticleName][ParticleSourceEnergy][SourceRegionName]["All"] = A1/A2;
        }

        if( __isnan(NOfValues[OrgansNameVector[gg]]) && __isinf(NOfValues[OrgansNameVector[gg]]) && NOfValues[OrgansNameVector[gg]] == 0 && NOfValues[OrgansNameVector[gg]] == NULL){
            G4cout << "\n\n !!!!!!!!!!!!!!!!!!!!!!!!!! " << OrgansNameVector[gg] << " total number of steps is not defined (" << TotalNumberOfSteps << ")\n";
        }else{
            //G4cout << TotalNumberOfSteps << G4endl;
            NumberOfSteps[ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = NOfValues[OrgansNameVector[gg]];
        }

        if(SourceRegionName == "Liver" && OrgansNameVector[gg] == "Liver"){
            //            G4cout << OrgansNameVector[gg] << " " << ParticleName  << " " << ParticleSourceEnergy << " " << SourceRegionName << " " << OrgansNameVector[gg] << " SAF " <<  ResultTable["H"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] << G4endl;
        }
        
        if(!__isinf(AD_Total[OrgansNameVector[gg]]) || !__isnan(AD_Total[OrgansNameVector[gg]])){
            //G4cout << OrgansNameVector[gg] << " " << ED_Total[OrgansNameVector[gg]] << " " << AD_Total[OrgansNameVector[gg]] << " " << TotalAbsorbedDose  << " ";
            TotalAbsorbedDose = TotalAbsorbedDose + AD_Total[OrgansNameVector[gg]];
            //G4cout << TotalAbsorbedDose << G4endl;
        }

        AF_Mean[OrgansNameVector[gg]]  = AF_Total[OrgansNameVector[gg]]/NOfValues[OrgansNameVector[gg]] ;
        SAF_Mean[OrgansNameVector[gg]] = SAF_Total[OrgansNameVector[gg]]/NOfValues[OrgansNameVector[gg]] ;
        AD_Mean[OrgansNameVector[gg]]  = AD_Total[OrgansNameVector[gg]]/NOfValues[OrgansNameVector[gg]] ;
        S_Mean[OrgansNameVector[gg]]   = S_Total[OrgansNameVector[gg]]/NOfValues[OrgansNameVector[gg]] ;
        H_Mean[OrgansNameVector[gg]]   = H_Total[OrgansNameVector[gg]]/NOfValues[OrgansNameVector[gg]] ;
        E_Mean[OrgansNameVector[gg]]   = E_Total[OrgansNameVector[gg]]/NOfValues[OrgansNameVector[gg]] ;
        DCC_Mean[OrgansNameVector[gg]]   = DCC_Total[OrgansNameVector[gg]]/NOfValues[OrgansNameVector[gg]] ;

        AF_Var[OrgansNameVector[gg]]  = AFCte[OrgansNameVector[gg]] * AFCte[OrgansNameVector[gg]] * ED_Var[OrgansNameVector[gg]] ;
        SAF_Var[OrgansNameVector[gg]] = SAFCte[OrgansNameVector[gg]] * SAFCte[OrgansNameVector[gg]] * ED_Var[OrgansNameVector[gg]] ;
        AD_Var[OrgansNameVector[gg]]  = ADCte[OrgansNameVector[gg]] * ADCte[OrgansNameVector[gg]] * ED_Var[OrgansNameVector[gg]] ;
        S_Var[OrgansNameVector[gg]]   = SCte[OrgansNameVector[gg]] * SCte[OrgansNameVector[gg]] * AD_Var[OrgansNameVector[gg]] ;
        H_Var[OrgansNameVector[gg]]   = HCte[OrgansNameVector[gg]] * HCte[OrgansNameVector[gg]] * AD_Var[OrgansNameVector[gg]] ;
        E_Var[OrgansNameVector[gg]]   = ECte[OrgansNameVector[gg]] * ECte[OrgansNameVector[gg]] * H_Var[OrgansNameVector[gg]] ;
        DCC_Var[OrgansNameVector[gg]]   = DCCCte[OrgansNameVector[gg]] * DCCCte[OrgansNameVector[gg]] * H_Var[OrgansNameVector[gg]] ;

        AF_SDev[OrgansNameVector[gg]] = std::sqrt(std::abs(AF_Var[OrgansNameVector[gg]]));
        SAF_SDev[OrgansNameVector[gg]] = std::sqrt(std::abs(SAF_Var[OrgansNameVector[gg]]));
        AD_SDev[OrgansNameVector[gg]] = std::sqrt(std::abs(AD_Var[OrgansNameVector[gg]]));
        S_SDev[OrgansNameVector[gg]] = std::sqrt(std::abs(S_Var[OrgansNameVector[gg]]));
        H_SDev[OrgansNameVector[gg]] = std::sqrt(std::abs(H_Var[OrgansNameVector[gg]]));
        E_SDev[OrgansNameVector[gg]] = std::sqrt(std::abs(E_Var[OrgansNameVector[gg]]));
        DCC_SDev[OrgansNameVector[gg]] = std::sqrt(std::abs(DCC_Var[OrgansNameVector[gg]]));

        StandardDeviation["AE"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = ED_SDev[OrgansNameVector[gg]];
        StandardDeviation["AF"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = AF_SDev[OrgansNameVector[gg]];
        StandardDeviation["SAF"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = SAF_SDev[OrgansNameVector[gg]];
        StandardDeviation["AD"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = AD_SDev[OrgansNameVector[gg]];
        StandardDeviation["S"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = S_SDev[OrgansNameVector[gg]];
        StandardDeviation["H"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = H_SDev[OrgansNameVector[gg]];
        StandardDeviation["E"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = E_SDev[OrgansNameVector[gg]];
        StandardDeviation["DCC"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = DCC_SDev[OrgansNameVector[gg]];

        TotalAEForRadiotracerRSD["AE"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = ED_Total[OrgansNameVector[gg]];

        if(!__isnan(E_Total[OrgansNameVector[gg]]) || !__isinf(E_Total[OrgansNameVector[gg]])){
            PhantomEffectiveDose = PhantomEffectiveDose + E_Total[OrgansNameVector[gg]];
            //G4cout << OrgansNameVector[gg] << " " << E_Total[OrgansNameVector[gg]] << " " << PhantomEffectiveDose<<G4endl;
        }

    }
    
    for(G4int gg = 0 ; gg < (G4int)OrgansNameVector.size() ; gg++){

        if(!__isinf(AD_Total[OrgansNameVector[gg]]) || !__isnan(AD_Total[OrgansNameVector[gg]])){
            DR_Total[OrgansNameVector[gg]] = (AD_Total[OrgansNameVector[gg]]/TotalAbsorbedDose)*100;

            ER_Total[OrgansNameVector[gg]] = (E_Total[OrgansNameVector[gg]]/PhantomEffectiveDose)*100;
            StandardDeviation["DR"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = AD_SDev[OrgansNameVector[gg]];
            StandardDeviation["ER"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] = E_SDev[OrgansNameVector[gg]];

            //if(ParticleName == "gamma" && SourceRegionName == "Liver" && ParticleSourceEnergy == 0.511){
            //G4cout << "Source " << SourceRegionName << " " << OrgansNameVector[gg] << " " << TotalAEForRadiotracerRSD["AE"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] << " " << StandardDeviation["DR"][ParticleName][ParticleSourceEnergy][SourceRegionName][OrgansNameVector[gg]] << G4endl;
            //}

            //G4cout << OrgansNameVector[gg] << " " << E_Total[OrgansNameVector[gg]] << " " << TotalEffectiveDose << " " << ER_Total[OrgansNameVector[gg]] << G4endl;

        }
    }
    
    //G4cout << " ----------- " <<G4endl;
    
    printf("\n");
    
}
void G4TResultCalculation::GenerateRegionResultFile(){
    

    // i want just the results of radionuclides

    //return;

    if(QuantityNamesToScore.size() == 1){
        if(QuantityNamesToScore[0] == "S"){
            return;
        }
    }









    if(V)G4cout << "\n-------- " << __FUNCTION__ << "\n" << G4endl ;
    G4cout << "\n1) Calculated Region Data For Particles:\n"<< G4endl;
    
    std::map<G4String,std::map<G4String,std::map<G4double,G4double>>> ModifiedRadioTracerEnergyPerCent = RadioTracerEnergyPerCent ;
    
    MaxRel_SDv=0; MinRel_SDv=100;
    for ( auto it = ED_RelS_D.begin(); it != ED_RelS_D.end(); ++it  ){
        
        //if(V)G4cout << it->first << " " << ED_RelS_D[it->first] << G4endl;
        
        if(ED_RelS_D[it->first] < MinRel_SDv){
            MinRel_SDv = ED_RelS_D[it->first];
            RegionForMinRel_SDv = it->first;
        }
        if(ED_RelS_D[it->first] > MaxRel_SDv){
            MaxRel_SDv = ED_RelS_D[it->first];
            RegionForMaxRel_SDv = it->first;
        }
    }

    //G4cout << " IsAllSourcesToScore " << IsAllSourcesToScore << G4endl;
    //G4cout << " IsAllTargetsToScore " << IsAllTargetsToScore << G4endl;

    if(IsAllTargetsToScore == true){ TargetNamesToScore.clear();for (G4int gg = 0 ; gg < OrgansNameVector.size() ; gg++) {TargetNamesToScore.push_back(OrgansNameVector[gg]);
            //G4cout << " OrgansNameVector[gg] "<< OrgansNameVector[gg] << G4endl;
        }}
    
    for (G4int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {

        //G4cout << gg << " Target To Score[gg] "<< TargetNamesToScore[gg] << G4endl;
    }

    for(G4int ikl = 0 ; ikl < (G4int)QuantityNamesToScore.size() ; ikl++){
        
        //if(V)G4cout << QuantityNamesToScore[ikl] << " ikl "<< ikl << " size " << QuantityNamesToScore.size() << "              ......  " << G4endl;
        
        std::ostringstream filname ;
        filname << ResultDirectoryPath << "/ResultsData";
        G4String fileNameString = filname.str().c_str();
        std::ofstream file(fileNameString.c_str(), std::ios_base::app);

        //if(FILE* file = fopen(fileNameString.c_str(),"a")){
        if(file.is_open()){
            
            if(V)G4cout << "Creating file  " << fileNameString.c_str() << " - writing the data output : \n"<< G4endl;
            
            std::ostringstream OutS;
            OutS << "****** "
                 << QuantityNamesToScore[ikl] << " "
                 << SourceRegionName << " "
                 << ParticleName << " "
                 << ParticleSourceEnergy << " "
                 << GeometrySymbol <<  " "
                 << EnergyDistribution <<  " "
                 << MomDirDistribution << " "
                 << SourceType << " "
                 << TotalEventNumber << "Event "
                 << TotalNumberOfSteps << "Step "
                 << TotalEmittedEnergy << "MeV "
                 << TotalAbsorbedEnergy - ED_Total["World"] << "MeV "
                 << Physics << " "
                 << CutsDistance << "mm "
                 << CutsEnergy << "MeV "
                 << RegionForMinRel_SDv << " "
                 << MinRel_SDv << "% "
                 << RegionForMaxRel_SDv << " "
                 << MaxRel_SDv << "% "
                 << (double)OneEventExecutionTimeInMs << "ms "
                 << (double)MaxExecutionTimeInMin << " "
                 << (double)ExecutionTimeInMin << " min ";
            
            OutS << "\n";
            
            
            G4String headerText= OutS.str().c_str();
            
            //std::fprintf(file,"%s", headerText.c_str());
            //std::fprintf(stdout,"%s",headerText.c_str());
            
            file << headerText.c_str();
            G4cout << headerText.c_str();

            G4String Quantity_NAME_UNIT = "";
            
            if(QuantityNamesToScore[ikl] == "AE" ){
                Quantity_NAME_UNIT = "AE [" + AEUnit +"] ";
                ChosenVariableTotal = ED_Total;
                ChosenVariableMean = ED_Mean;
                ChosenVariablevariance = ED_Var;
                ChosenVariableStandardDeviation = ED_SDev;
            }
            else if(QuantityNamesToScore[ikl] == "AF"){
                Quantity_NAME_UNIT = "AF[" + AFUnit + "] ";
                ChosenVariableTotal = AF_Total;
                ChosenVariableMean = AF_Mean;
                ChosenVariablevariance = AF_Var;
                ChosenVariableStandardDeviation = AF_SDev;
            }
            else if(QuantityNamesToScore[ikl] == "SAF"){
                Quantity_NAME_UNIT = "SAF[" + SAFUnit + "] ";
                ChosenVariableTotal = SAF_Total;
                ChosenVariableMean = SAF_Mean;
                ChosenVariablevariance = SAF_Var;
                ChosenVariableStandardDeviation = SAF_SDev;
            }
            else if(QuantityNamesToScore[ikl] == "DR" ){
                Quantity_NAME_UNIT = "DR [%]";
                ChosenVariableTotal = DR_Total;
                ChosenVariableMean = AD_Mean;
                ChosenVariablevariance = AD_Var;
                ChosenVariableStandardDeviation = AD_SDev;
            }
            else if(QuantityNamesToScore[ikl] == "ER" ){
                Quantity_NAME_UNIT = "ER [%]";
                ChosenVariableTotal = ER_Total;
                ChosenVariableMean = E_Mean;
                ChosenVariablevariance = E_Var;
                ChosenVariableStandardDeviation = E_SDev;
            }
            else if(QuantityNamesToScore[ikl] == "AD" ){
                Quantity_NAME_UNIT = "AD[(" + ADUnit+ ")/"+UnitPerParticle+"]";
                ChosenVariableTotal = AD_Total;
                ChosenVariableMean = AD_Mean;
                ChosenVariablevariance = AD_Var;
                ChosenVariableStandardDeviation = AD_SDev;
            }
            else if(QuantityNamesToScore[ikl] == "S"){
                Quantity_NAME_UNIT = "S[(" + SUnit+ ")/"+UnitPerParticle+"]";
                ChosenVariableTotal = S_Total;
                ChosenVariableMean = S_Mean;
                ChosenVariablevariance = S_Var;
                ChosenVariableStandardDeviation = S_SDev;
            }
            else if(QuantityNamesToScore[ikl] == "H"){
                Quantity_NAME_UNIT = "H[(" + HUnit+ ")/"+UnitPerParticle+"]";
                ChosenVariableTotal = H_Total;
                ChosenVariableMean = H_Mean;
                ChosenVariablevariance = H_Var;
                ChosenVariableStandardDeviation = H_SDev;
            }
            else if(QuantityNamesToScore[ikl] == "E"){
                Quantity_NAME_UNIT = "E[(" + EUnit+ ")/"+UnitPerParticle+"]";
                ChosenVariableTotal = E_Total;
                ChosenVariableMean = E_Mean;
                ChosenVariablevariance = E_Var;
                ChosenVariableStandardDeviation = E_SDev;
            }
            else if(QuantityNamesToScore[ikl] == "DCC"){
                Quantity_NAME_UNIT = "DCC[(" + DCCUnit+ ")]";
                ChosenVariableTotal = DCC_Total;
                ChosenVariableMean = DCC_Mean;
                ChosenVariablevariance = DCC_Var;
                ChosenVariableStandardDeviation = DCC_SDev;
            }
            else{
                
                if(V)G4cout << " No quantity to score is chosen (default is AE)" << G4endl;
                Quantity_NAME_UNIT = "AE [" + AEUnit + "] ";
                //QuantitiesToScore = "AE";
                ChosenVariableTotal = ED_Total;
                ChosenVariableMean = ED_Mean;
                ChosenVariablevariance = ED_Var;
                ChosenVariableStandardDeviation = ED_SDev;
            }

            G4cout << std::setw(24) << std::left << "# Volume" << " "
                   << std::setw(SZ) << std::left << Quantity_NAME_UNIT.c_str() << " "
                   << std::setw(SZ) << std::left << "SDev"  << " "
                   << std::setw(SZ) << std::left << "Rel_SDev(%)" << " "
                   << std::setw(SZ) << std::left << "Values Num" << " "
                   << std::setw(SZ) << std::left << "Mass[kg]" << " "
                   << std::setw(SZ) << std::left << "Volume[cm3]" << " "
                   << std::setw(SZ) << std::left << "Density[g/cm3] " << " \n";

            file << std::setw(24) << std::left << "# Volume" << " "
                 << std::setw(SZ) << std::left << Quantity_NAME_UNIT.c_str() << " "
                 << std::setw(SZ) << std::left << "SDev"  << " "
                 << std::setw(SZ) << std::left << "Rel_SDev(%)" << " "
                 << std::setw(SZ) << std::left << "Values Num" << " "
                 << std::setw(SZ) << std::left << "Mass[kg]" << " "
                 << std::setw(SZ) << std::left << "Volume[cm3]" << " "
                 << std::setw(SZ) << std::left << "Density[g/cm3] " << " \n";

            //std::fprintf(file  ,"%-24s%-1s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n", "# Volume" ," ", Quantity_NAME_UNIT.c_str() , "SDev" , "Rel_SDev(%)" , "Values Num" , "Mass[kg]" , "Volume[cm3]", "Density[g/cm3]","\n");
            //std::fprintf(stdout  ,"%-24s%-1s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n", "# Volume" ," ", Quantity_NAME_UNIT.c_str() , "SDev" , "Rel_SDev(%)" , "Values Num" , "Mass[kg]" , "Volume[cm3]", "Density[g/cm3]","\n");
            
            for(G4int dd = 0 ; dd < (G4int)TargetNamesToScore.size(); dd++){
                if(TargetNamesToScore[dd] == "World" || TargetNamesToScore[dd] == "VOXEL"){continue;}
                
                if(!__isnan(ChosenVariableTotal[TargetNamesToScore[dd]]) || !__isinf(ChosenVariableTotal[TargetNamesToScore[dd]]) || ChosenVariableTotal[TargetNamesToScore[dd]] != NULL){

                    //G4cout.setf(std::ios::fixed,std::ios::floatfield);

                    G4cout << std::setw(24) << std::left << TargetNamesToScore[dd].c_str() << " "
                           << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << ChosenVariableTotal[TargetNamesToScore[dd]] << " "
                           << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << ChosenVariableStandardDeviation[TargetNamesToScore[dd]]  << " "
                           << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << ED_RelS_D[TargetNamesToScore[dd]] << " "
                           << std::resetiosflags(G4cout.basefield) << std::resetiosflags( G4cout.floatfield) << std::resetiosflags( G4cout.flags())
                           << std::setw(SZ) << std::left << NOfValues[TargetNamesToScore[dd]] << " "
                           << std::setw(SZ) << std::left << VolumeNameMassMap[TargetNamesToScore[dd]] << " "
                           << std::setw(SZ) << std::left << VolumeNameVolumeMap[TargetNamesToScore[dd]] << " "
                           << std::setw(SZ) << std::left << VolumeNameDensityMap[TargetNamesToScore[dd]] << " \n";

                    file << std::setw(24) << std::left << TargetNamesToScore[dd].c_str() << " "
                         << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << ChosenVariableTotal[TargetNamesToScore[dd]] << " "
                         << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << ChosenVariableStandardDeviation[TargetNamesToScore[dd]]  << " "
                         << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << ED_RelS_D[TargetNamesToScore[dd]] << " "
                         << std::resetiosflags(file.basefield) << std::resetiosflags( file.floatfield) << std::resetiosflags( file.flags())
                         << std::setw(SZ) << std::left << NOfValues[TargetNamesToScore[dd]] << " "
                         << std::setw(SZ) << std::left << VolumeNameMassMap[TargetNamesToScore[dd]] << " "
                         << std::setw(SZ) << std::left << VolumeNameVolumeMap[TargetNamesToScore[dd]] << " "
                         << std::setw(SZ) << std::left << VolumeNameDensityMap[TargetNamesToScore[dd]] << " \n";

                    //std::fprintf(file  ,"%-24s%-1s%-15e%-15e%-15.6f%-15u%-15.6f%-15.3f%-15.6f\n", TargetNamesToScore[dd].c_str() ," ", (double) ChosenVariableTotal[TargetNamesToScore[dd]] , (double)ChosenVariableStandardDeviation[TargetNamesToScore[dd]] , (double)ED_RelS_D[TargetNamesToScore[dd]] , (unsigned long long int)NOfValues[TargetNamesToScore[dd]] , (double)VolumeNameMassMap[TargetNamesToScore[dd]] , (double)VolumeNameVolumeMap[TargetNamesToScore[dd]] , (double)VolumeNameDensityMap[TargetNamesToScore[dd]],"\n");
                    //std::fprintf(stdout,"%-24s%-1s%-15e%-15e%-15.6f%-15u%-15.6f%-15.3f%-15.6f\n", TargetNamesToScore[dd].c_str() ," ", (double) ChosenVariableTotal[TargetNamesToScore[dd]] , (double)ChosenVariableStandardDeviation[TargetNamesToScore[dd]] , (double)ED_RelS_D[TargetNamesToScore[dd]] , (unsigned long long int)NOfValues[TargetNamesToScore[dd]] , (double)VolumeNameMassMap[TargetNamesToScore[dd]] , (double)VolumeNameVolumeMap[TargetNamesToScore[dd]] , (double)VolumeNameDensityMap[TargetNamesToScore[dd]],"\n");
                }
            }
            //std::fprintf(file,"* ----------------------------------------------------------------------------------------------------------------------------------------------\n");
            //std::fprintf(stdout,"* ----------------------------------------------------------------------------------------------------------------------------------------------\n");

            if(QuantityNamesToScore[ikl] == "DCC"){

                G4cout << std::setw(24) << std::left << "Ambient_DCC_ToED " << " "
                       << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << ResultTable["DCCAm"][ParticleName][ParticleSourceEnergy][SourceRegionName]["All"] << " "
                       << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << "1 "  << " "
                       << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << "1 "
                       << std::resetiosflags(G4cout.basefield) << std::resetiosflags( G4cout.floatfield) << std::resetiosflags( G4cout.flags())
                       << std::setw(SZ) << std::left << "1 "
                       << std::setw(SZ) << std::left << "1 "
                       << std::setw(SZ) << std::left << "1 "
                       << std::setw(SZ) << std::left << "1 \n";

                file << std::setw(24) << std::left << "Ambient_DCC_ToED " << " "
                     << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << ResultTable["DCCAm"][ParticleName][ParticleSourceEnergy][SourceRegionName]["All"] << " "
                     << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << "1 "  << " "
                     << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << "1 "
                     << std::resetiosflags(G4cout.basefield) << std::resetiosflags( G4cout.floatfield) << std::resetiosflags( G4cout.flags())
                     << std::setw(SZ) << std::left << "1 "
                     << std::setw(SZ) << std::left << "1 "
                     << std::setw(SZ) << std::left << "1 "
                     << std::setw(SZ) << std::left << "1 \n";
            }
            G4cout << "* ----------------------------------------------------------------------------------------------------------------------------------------------" << G4endl;
            file << "* ----------------------------------------------------------------------------------------------------------------------------------------------\n";

            //std::fclose (file);
            file.close();
        }
        else{
            
            G4cout << "Cannot open the file " << fileNameString.c_str() << G4endl;
        }
    }
}
void G4TResultCalculation::GenerateExternalDosimetryCoefficients(){

    if(V)G4cout << "\n-------- " << __FUNCTION__ << "\n" << G4endl ;

    G4String RadioTracer_NAME;
    G4String Quantity_NAME;
    G4String Particle_NAME;
    G4String Source_NAME;
    G4String Target_NAME;

    G4double Energy_Val;

    G4String Quantity_NAME_UNIT;

    // calculation of radiotracer quantity for src-trg combination and for organ
    for ( auto it = RadioTracerEnergyPerCent.begin(); it != RadioTracerEnergyPerCent.end(); ++it  ){

        RadioTracer_NAME = it->first;

        RadiTracerParticleEnergyDataString[RadioTracer_NAME].str("");
        RadiTracerParticleEnergyDataString[RadioTracer_NAME].clear();
        RadiTracerDataForTotalDoseString[RadioTracer_NAME].str("");
        RadiTracerDataForTotalDoseString[RadioTracer_NAME].clear();

        G4cout << "Generating Results for: " << RadioTracer_NAME <<G4endl;

        for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2  ){
            Particle_NAME = it2->first;
            //G4cout << "Particle_NAME " << Particle_NAME <<G4endl;

            for ( auto it3 = it2->second.begin(); it3 != it2->second.end(); ++it3  ){
                Energy_Val = it3->first;
                //G4cout << "Energy_Val " << Energy_Val <<G4endl;

                if(ResultTable["AE"][Particle_NAME][Energy_Val].size() == 0){
                    //G4cout << "\n-----> For Radiotracer " << RadioTracer_NAME
                    //       << ", the data for particle "<< Particle_NAME << " with energy "<< Energy_Val
                    //       << " are not found (this source configuration not simulated), the related quantities values will be generated by interpolating the already existed data of "
                    //       << Particle_NAME << " with energies surround the "
                    //       << Energy_Val << " from the existing source regions." << G4endl;
                    if(RadiotracerDataFomFile == false){
                        RadiTracerParticleEnergyDataString[RadioTracer_NAME] << " | (by interpolation) "<< Particle_NAME << " " << Energy_Val << "MeV " << RadioTracerEnergyPerCent[RadioTracer_NAME][Particle_NAME][Energy_Val] << "%";
                        RadiTracerDataForTotalDoseString[RadioTracer_NAME] << " | (by interpolation) "<< Particle_NAME << " " << Energy_Val << "MeV " << RadioTracerEnergyPerCent[RadioTracer_NAME][Particle_NAME][Energy_Val] << "%";
                    }
                    GenerateRadiotracerQuantitiesByInterpolation(Particle_NAME,Energy_Val);
                }else{
                    if(RadiotracerDataFomFile == false){
                        RadiTracerParticleEnergyDataString[RadioTracer_NAME] << " | (by simulation) "<< Particle_NAME << " " << Energy_Val << "MeV " << RadioTracerEnergyPerCent[RadioTracer_NAME][Particle_NAME][Energy_Val] << "%";
                        RadiTracerDataForTotalDoseString[RadioTracer_NAME] << " | (by simulation) "<< Particle_NAME << " " << Energy_Val << "MeV " << RadioTracerEnergyPerCent[RadioTracer_NAME][Particle_NAME][Energy_Val] << "%";
                    }
                }

                TotalRadioTracerEmmitedEnergy_Val[RadioTracer_NAME] += Energy_Val;

                // this is to provide results for all sources that can user simulate without taking into account the RadioTracer source Organs
                for ( auto DD = ResultTable["S"][Particle_NAME][Energy_Val].begin(); DD != ResultTable["S"][Particle_NAME][Energy_Val].end(); ++DD  ){
                    Source_NAME  = DD->first;
                    //G4cout << "Source_NAME " << Source_NAME <<G4endl;

                    for ( auto XX = RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracer_NAME].begin(); XX != RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracer_NAME].end(); ++XX  ){
                        G4String SrcOrg = XX->first;
                        if(Source_NAME == SrcOrg){
                            if(RadiotracerDataFomFile == false){
                                RadiTracerDataForTotalDoseString[RadioTracer_NAME] << " " << Source_NAME <<" As/A0(s)=" << RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracer_NAME][Source_NAME] ;
                            }
                        }
                    }

                    for ( auto CC = DD->second.begin(); CC != DD->second.end(); ++CC  ){
                        Target_NAME  = CC->first;

                        double vb= ResultTable["AE"][Particle_NAME][Energy_Val][Source_NAME][Target_NAME];
                        if(__isinf(vb) || __isnan(vb) || vb == 0.){continue;}

                        G4double RadiationPerCent = RadioTracerEnergyPerCent[RadioTracer_NAME][Particle_NAME][Energy_Val]/100;

                        G4double CumulatedActivityInSource = RadioTracerInjectedActivity[RadioTracer_NAME]
                                *RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracer_NAME][Source_NAME];

                        // values and variance calculation for chosen quantities
                        for (int rr = 0 ; rr < DoseCalcsQuantities.size() ; rr++) {

                            double ccc = ResultTable[DoseCalcsQuantities[rr]][Particle_NAME][Energy_Val][Source_NAME][Target_NAME];
                            if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                                RadioTracerQuantitySourceTargetValue[RadioTracer_NAME][DoseCalcsQuantities[rr]][Source_NAME][Target_NAME] +=
                                        RadiationPerCent
                                        *ccc;


                                // this condition just to show the regions with Wt different to 0 for AD in intake
                                if( DoseCalcsQuantities[rr] == "AD" && (__isnan(TissueFactorMap[Target_NAME])
                                                                        || __isinf(TissueFactorMap[Target_NAME])
                                                                        || TissueFactorMap[Target_NAME] == 0
                                                                        || TissueFactorMap[Target_NAME] == NULL)){

                                }else{
                                    RadioTracerQuantityOrganValue[RadioTracer_NAME][DoseCalcsQuantities[rr]][Target_NAME] +=
                                            RadiationPerCent
                                            *CumulatedActivityInSource
                                            *ccc;
                                }
                            }

                        }

                        if( !__isnan(NumberOfSteps[Particle_NAME][Energy_Val][Source_NAME][Target_NAME]) && !__isinf(NumberOfSteps[Particle_NAME][Energy_Val][Source_NAME][Target_NAME]) && NumberOfSteps[Particle_NAME][Energy_Val][Source_NAME][Target_NAME] != 0 && NumberOfSteps[Particle_NAME][Energy_Val][Source_NAME][Target_NAME] != NULL){

                            NumberOfStepsInRadiotracerSourceTarget[RadioTracer_NAME][Source_NAME][Target_NAME] +=
                                    NumberOfSteps[Particle_NAME][Energy_Val][Source_NAME][Target_NAME] ;

                            NumberOfStepsInRadiotracerOrgan[RadioTracer_NAME][Target_NAME] +=
                                    NumberOfSteps[Particle_NAME][Energy_Val][Source_NAME][Target_NAME] ;

                        }

                        double ccc = TotalAEForRadiotracerRSD["AE"][Particle_NAME][Energy_Val][Source_NAME][Target_NAME];

                        if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                            RadioTracerQuantitySourceTargetAETotalForRSD[RadioTracer_NAME]["AE"][Source_NAME][Target_NAME] +=
                                    ccc;

                            //G4cout << " Total for " << Particle_NAME << " "<< Energy_Val << " " << Source_NAME << " "<< Target_NAME << " is " << ccc << G4endl;
                            RadioTracerQuantityOrganAETotalForRSD[RadioTracer_NAME]["AE"][Target_NAME] +=
                                    ccc;
                        }

                        ccc = StandardDeviation["AE"][Particle_NAME][Energy_Val][Source_NAME][Target_NAME];
                        if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                            RadioTracerQuantitySourceTargetVariance[RadioTracer_NAME]["AE"][Source_NAME][Target_NAME] +=
                                    ((RadiationPerCent/100)*ccc)*((RadiationPerCent/100)*ccc);

                            //G4cout << " SD for AE " << Particle_NAME << " "<< Energy_Val << " " << Source_NAME << " "<< Target_NAME << " is " << ccc
                            //<< " . and Total " << RadioTracerQuantitySourceTargetVariance[RadioTracer_NAME]["AE"][Source_NAME][Target_NAME] << G4endl ;

                            RadioTracerQuantityOrganVAR[RadioTracer_NAME]["AE"][Target_NAME] +=
                                    ((RadiationPerCent/100)*ccc)*((RadiationPerCent/100)*ccc);
                        }
                    }
                }
            }
        }
        RadiTracerParticleEnergyDataString[RadioTracer_NAME] << " | ";
        RadiTracerDataForTotalDoseString[RadioTracer_NAME] << " | ";
    }


    // to calculate quantity ratio ER for each source-target
    for ( auto DoseCalcsQuantities = RadioTracerQuantitySourceTargetValue.begin(); DoseCalcsQuantities != RadioTracerQuantitySourceTargetValue.end(); ++DoseCalcsQuantities  ){
        RadioTracer_NAME = DoseCalcsQuantities->first;
        for ( auto it3 = DoseCalcsQuantities->second.begin(); it3 != DoseCalcsQuantities->second.end(); ++it3  ){
            Quantity_NAME = it3->first;
            for ( auto it4 = it3->second.begin(); it4 != it3->second.end(); ++it4  ){
                Source_NAME = it4->first;
                G4double CumulatedActivityInSource = (RadioTracerInjectedActivity[RadioTracer_NAME]
                                                      *RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracer_NAME][Source_NAME]);
                TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME] = 0;

                long double val = 0. ;
                for ( auto it5 = it4->second.begin(); it5 != it4->second.end(); ++it5  ){
                    Target_NAME = it5->first;

                    // to calculate ER for each source-target
                    if(Quantity_NAME == "E" || Quantity_NAME == "AD"){
                        double ccc = RadioTracerQuantitySourceTargetValue[RadioTracer_NAME][Quantity_NAME][Source_NAME][Target_NAME];
                        if( !__isnan(ccc) && !__isinf(ccc) && ccc != NULL){
                            TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME] = TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME] + ccc;
                        }
                    }
                    if( !__isnan(it5->second) && !__isinf(it5->second) && it5->second != 0 && it5->second != NULL){
                        val = val + it5->second;
                        //G4cout <<  " DoseCalcsQuantities->first " << DoseCalcsQuantities->first << " it3->first " << it3->first << " it4->first " << it4->first << " it5->first " << it5->first <<  " it5->second " << it5->second << " val: " << val << G4endl ;
                    }
                }
                // to calculate ER for each source-target
                if(Quantity_NAME == "E" || Quantity_NAME == "AD"){
                    if( !__isnan(TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME]) && !__isinf(TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME]) && TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME] != 0 && TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME] != NULL){
                        //double val = 0;
                        for ( auto it5 = it4->second.begin(); it5 != it4->second.end(); ++it5  ){
                            Target_NAME = it5->first;
                            if(Quantity_NAME == "E"){
                                double ccc = RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["E"][Source_NAME][Target_NAME];
                                if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                                    RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["ER"][Source_NAME][Target_NAME] =
                                            (RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["E"][Source_NAME][Target_NAME]/TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME])*100;
                                    //val = val + RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["ER"][Source_NAME][Target_NAME];
                                    //G4cout << " Source_NAME: " << Source_NAME << " Target_NAME: " << Target_NAME << " val: " << val << G4endl ;
                                }
                            }
                            if(Quantity_NAME == "AD"){
                                double ccc = RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["AD"][Source_NAME][Target_NAME];
                                if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                                    RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["DR"][Source_NAME][Target_NAME] =
                                            (RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["AD"][Source_NAME][Target_NAME]/TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME])*100;
                                    //val = val + RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["ER"][Source_NAME][Target_NAME];
                                    //G4cout << " Source_NAME: " << Source_NAME << " Target_NAME: " << Target_NAME << " val: " << val << G4endl ;
                                }
                            }
                        }
                    }
                }
                TotalDoseFromRadioTracer[RadioTracer_NAME][Quantity_NAME] += val * CumulatedActivityInSource;
                //G4cout << " CumulatedActivityInSource " << CumulatedActivityInSource << " val: " << val << " TotalDoseFromRadioTracer[RadioTracer_NAME][Quantity_NAME] " << TotalDoseFromRadioTracer[RadioTracer_NAME][Quantity_NAME] << G4endl ;
            }
        }
    }


    // to calculate quantity ratio ER for each organ
    for ( auto DoseCalcsQuantities = RadioTracerQuantityOrganValue.begin(); DoseCalcsQuantities != RadioTracerQuantityOrganValue.end(); ++DoseCalcsQuantities  ){

        RadioTracer_NAME = DoseCalcsQuantities->first;
        TotalValueOfQuantity[RadioTracer_NAME]["E"] = 0;
        TotalValueOfQuantity[RadioTracer_NAME]["AD"] = 0;

        //For ER
        for ( auto it4 = DoseCalcsQuantities->second["E"].begin(); it4 != DoseCalcsQuantities->second["E"].end(); ++it4  ){
            Target_NAME = it4->first;
            double ccc = it4->second;
            if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                TotalValueOfQuantity[RadioTracer_NAME]["E"] = TotalValueOfQuantity[RadioTracer_NAME]["E"] + ccc;
            }
        }
        //double val;
        for ( auto it4 = DoseCalcsQuantities->second["E"].begin(); it4 != DoseCalcsQuantities->second["E"].end(); ++it4  ){
            Target_NAME = it4->first;
            double ccc = it4->second;
            if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){

                RadioTracerQuantityOrganValue[RadioTracer_NAME]["ER"][Target_NAME] = (ccc/TotalValueOfQuantity[RadioTracer_NAME]["E"])*100;
                //val = val + RadioTracerQuantityOrganValue[RadioTracer_NAME]["ER"][Target_NAME];
                //G4cout << " Target_NAME: " << Target_NAME << " val: " << val << G4endl ;
            }
        }

        // For DR
        for ( auto it4 = DoseCalcsQuantities->second["AD"].begin(); it4 != DoseCalcsQuantities->second["AD"].end(); ++it4  ){
            Target_NAME = it4->first;
            double ccc = it4->second;
            if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                TotalValueOfQuantity[RadioTracer_NAME]["AD"] = TotalValueOfQuantity[RadioTracer_NAME]["AD"] + ccc;
            }
        }
        //double val;
        for ( auto it4 = DoseCalcsQuantities->second["AD"].begin(); it4 != DoseCalcsQuantities->second["AD"].end(); ++it4  ){
            Target_NAME = it4->first;
            double ccc = it4->second;
            if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){

                RadioTracerQuantityOrganValue[RadioTracer_NAME]["DR"][Target_NAME] = (ccc/TotalValueOfQuantity[RadioTracer_NAME]["AD"])*100;
                //val = val + RadioTracerQuantityOrganValue[RadioTracer_NAME]["ER"][Target_NAME];
                //G4cout << " Target_NAME: " << Target_NAME << " val: " << val << G4endl ;
            }
        }
    }


    if(IsAllTargetsToScore == true){ TargetNamesToScore.clear();for (G4int gg = 0 ; gg < OrgansNameVector.size() ; gg++) {TargetNamesToScore.push_back(OrgansNameVector[gg]);}}

    for ( auto it0 = RadioTracerQuantitySourceTargetValue.begin(); it0 != RadioTracerQuantitySourceTargetValue.end(); ++it0  ){

        RadioTracer_NAME = it0->first;

        for ( auto it1 = it0->second.begin(); it1 != it0->second.end(); ++it1  ){

            Quantity_NAME = it1->first;

            // to write just the results of chosen quantities
            bool IsIn = false;
            for (int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {if(QuantityNamesToScore[gg] == Quantity_NAME){IsIn = true; break; }}
            if(IsIn == false){continue;}


            if(Quantity_NAME == "AF"){
                Quantity_NAME_UNIT = "AF[" + AFUnit + "]";
            }
            else if(Quantity_NAME == "SAF"){
                Quantity_NAME_UNIT = "SAF[" + SAFUnit + "]";
            }
            else if(Quantity_NAME == "S"){
                Quantity_NAME_UNIT = "S[(" + SUnit+ ")/"+UnitPerRadioTracerDecay+"]";
            }
            else if(Quantity_NAME == "AD"){
                Quantity_NAME_UNIT = "AD[(" + ADUnit+ ")/"+UnitPerRadioTracerDecay+"]";
            }
            else if(Quantity_NAME == "H"){
                Quantity_NAME_UNIT = "H[(" + HUnit+ ")/"+UnitPerRadioTracerDecay+"]";
            }
            else if(Quantity_NAME == "E"){
                Quantity_NAME_UNIT = "E[(" + EUnit+ ")/"+UnitPerRadioTracerDecay+"]";
            }
            else if(Quantity_NAME == "ER"){
                Quantity_NAME_UNIT = "ER[\%/"+UnitPerRadioTracerDecay+"]";
            }

            for ( auto it2 = it1->second.begin(); it2 != it1->second.end(); ++it2  ){
                Source_NAME = it2->first;

                G4cout << "\n\n================================== Generation of " << Quantity_NAME_UNIT << " deposition data in body target organs for intake of radio-tracer: " << RadioTracer_NAME << " into simulated source organ: " << Source_NAME << " ==================================\n" << G4endl ;

                //bool IsIn = false;
                //for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                //    if(TargetNamesToScore[gg] == Source_NAME){IsIn = true; break; }
                //}
                //if(IsIn == false){continue;}

                std::ostringstream filname ;
                filname << ResultDirectoryPath << "/ResultsData";
                G4String fileNameString = filname.str().c_str();
                std::ofstream file(fileNameString.c_str(), std::ios_base::app);

                //if(FILE* file = fopen(fileNameString.c_str(),"a")){
                if(file.is_open()){

                    if(V)G4cout << "Creating file  " << fileNameString.c_str() << " - writing the data output : \n"<< G4endl;

                    std::ostringstream OutS;
                    OutS << "****** "
                         << Quantity_NAME << " "
                         << Source_NAME << " "
                         << "RadioTracer "
                         << RadioTracer_NAME << " "
                         << GeometrySymbol <<  " "
                         << Physics << " "
                         << TotalRadioTracerEmmitedEnergy_Val[RadioTracer_NAME] << " "
                         << RadiTracerParticleEnergyDataString[RadioTracer_NAME].str().c_str() << " ";
                    OutS << "\n";

                    G4String headerText= OutS.str().c_str();

                    file << headerText.c_str();
                    G4cout << headerText.c_str();

                    G4cout << std::setw(24) << std::left << "# Volume" << " "
                           << std::setw(SZ) << std::left << Quantity_NAME_UNIT.c_str() << " "
                           << std::setw(SZ) << std::left << "SDev"  << " "
                           << std::setw(SZ) << std::left << "Rel_SDev(%)" << " "
                           << std::setw(SZ) << std::left << "Values Num" << " "
                           << std::setw(SZ) << std::left << "Mass[kg]" << " "
                           << std::setw(SZ) << std::left << "Volume[cm3]" << " "
                           << std::setw(SZ) << std::left << "Density[g/cm3] " << " \n";

                    file << std::setw(24) << std::left << "# Volume" << " "
                         << std::setw(SZ) << std::left << Quantity_NAME_UNIT.c_str() << " "
                         << std::setw(SZ) << std::left << "SDev"  << " "
                         << std::setw(SZ) << std::left << "Rel_SDev(%)" << " "
                         << std::setw(SZ) << std::left << "Values Num" << " "
                         << std::setw(SZ) << std::left << "Mass[kg]" << " "
                         << std::setw(SZ) << std::left << "Volume[cm3]" << " "
                         << std::setw(SZ) << std::left << "Density[g/cm3] " << " \n";

                    for ( auto it3 = it2->second.begin(); it3 != it2->second.end(); ++it3  ){
                        Target_NAME = it3->first;
                        if(Target_NAME == "World"){continue;}

                        bool IsIn = false;
                        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                            if(TargetNamesToScore[gg] == Target_NAME){IsIn = true; break; }
                        }
                        if(IsIn == true){

                            double sdv = std::sqrt(std::abs(RadioTracerQuantitySourceTargetVariance[RadioTracer_NAME]["AE"][Source_NAME][Target_NAME]));
                            unsigned long long int sm = NumberOfStepsInRadiotracerSourceTarget[RadioTracer_NAME][Source_NAME][Target_NAME];
                            double rsdv = ((sdv * sm)/ RadioTracerQuantitySourceTargetAETotalForRSD[RadioTracer_NAME]["AE"][Source_NAME][Target_NAME])*100;

                            //G4cout << " Total for 100*" << sdv << "*"<< sm << "/" << RadioTracerQuantitySourceTargetAETotalForRSD[RadioTracer_NAME]["AE"][Source_NAME][Target_NAME] << "="<< rsdv << G4endl;

                            G4cout << std::setw(24) << std::left << Target_NAME.c_str() << " "
                                   << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << it3->second << " "
                                   << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << sdv << " "
                                   << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << rsdv << " "
                                   << std::resetiosflags(G4cout.basefield) << std::resetiosflags( G4cout.floatfield) << std::resetiosflags( G4cout.flags())
                                   << std::setw(SZ) << std::left << sm << " "
                                   << std::setw(SZ) << std::left << VolumeNameMassMap[Target_NAME] << " "
                                   << std::setw(SZ) << std::left << VolumeNameVolumeMap[Target_NAME] << " "
                                   << std::setw(SZ) << std::left << VolumeNameDensityMap[Target_NAME] << " \n";

                            file << std::setw(24) << std::left << Target_NAME.c_str() << " "
                                 << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << it3->second << " "
                                 << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << sdv << " "
                                 << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << rsdv << " "
                                 << std::resetiosflags(file.basefield) << std::resetiosflags( file.floatfield) << std::resetiosflags( file.flags())
                                 << std::setw(SZ) << std::left << sm << " "
                                 << std::setw(SZ) << std::left << VolumeNameMassMap[Target_NAME] << " "
                                 << std::setw(SZ) << std::left << VolumeNameVolumeMap[Target_NAME] << " "
                                 << std::setw(SZ) << std::left << VolumeNameDensityMap[Target_NAME] << " \n";

                        }
                        else{continue;}
                    }
                    //std::fprintf(file,"* ----------------------------------------------------------------------------------------------------------------------------------------------\n");
                    //std::fprintf(stdout,"* ----------------------------------------------------------------------------------------------------------------------------------------------\n");

                    G4cout << "* ----------------------------------------------------------------------------------------------------------------------------------------------" << G4endl;
                    file << "* ----------------------------------------------------------------------------------------------------------------------------------------------\n";

                    file.close();

                }
            }
        }
    }

    if(GenerateResultsForRadioTracerExams){

        for ( auto it0 = RadioTracerQuantityOrganValue.begin(); it0 != RadioTracerQuantityOrganValue.end(); ++it0  ){

            RadioTracer_NAME = it0->first;

            if(RadioTracerInjectedActivity[RadioTracer_NAME] == 0.){
                continue;
            }

            for ( auto it1 = it0->second.begin(); it1 != it0->second.end(); ++it1  ){

                Quantity_NAME = it1->first;

                // to write just the results of chosen quantities
                bool IsIn = false;
                for (int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {if(QuantityNamesToScore[gg] == Quantity_NAME){IsIn = true; break; }}
                if(IsIn == false){continue;}


                if(Quantity_NAME == "AF"){
                    Quantity_NAME_UNIT = "AF[" + AFUnit + "]";
                }
                else if(Quantity_NAME == "SAF"){
                    Quantity_NAME_UNIT = "SAF[" + SAFUnit + "]";
                }
                else if(Quantity_NAME == "S"){
                    Quantity_NAME_UNIT = "S[(" + SUnit+ ")/"+UnitPerRadioTracerDecay+"]";
                }
                else if(Quantity_NAME == "AD"){
                    Quantity_NAME_UNIT = "AD[(" + ADUnit+ ")/"+UnitPerRadioTracerDecay+"]";
                }
                else if(Quantity_NAME == "H"){
                    Quantity_NAME_UNIT = "H[(" + HUnit+ ")/"+UnitPerRadioTracerDecay+"]";
                }
                else if(Quantity_NAME == "E"){
                    Quantity_NAME_UNIT = "E[(" + EUnit+ ")/"+UnitPerRadioTracerDecay+"]";
                }
                else if(Quantity_NAME == "ER"){
                    Quantity_NAME_UNIT = "ER[\%/"+UnitPerRadioTracerDecay+"]";
                }
                else if(Quantity_NAME == "DR"){
                    Quantity_NAME_UNIT = "DR[\%/"+UnitPerRadioTracerDecay+"]";
                }

                G4cout << "\n\n================================== Generation of " << Quantity_NAME_UNIT << " deposition data in all simulated body organs for " << RadioTracerInjectedActivity[RadioTracer_NAME] << " Bq intake of radio-tracer: " << RadioTracer_NAME << " into " << GeometrySymbol<< " simulated body ==================================\n" << G4endl ;

                std::ostringstream filname ;
                filname << ResultDirectoryPath << "/ResultsData";
                G4String fileNameString = filname.str().c_str();
                std::ofstream file(fileNameString.c_str(), std::ios_base::app);

                //if(FILE* file = fopen(fileNameString.c_str(),"a")){
                if(file.is_open()){

                    if(V)G4cout << "Creating file  " << fileNameString.c_str() << " - writing the data output : \n"<< G4endl;

                    std::ostringstream OutS;
                    OutS << "****** "
                         << Quantity_NAME << " "
                         << "IntakeIntoBody "
                         << "RadioTracer "
                         << RadioTracer_NAME << " "
                         << GeometrySymbol <<  " "
                         << Physics << " "
                         << TotalRadioTracerEmmitedEnergy_Val[RadioTracer_NAME] << " "
                         << RadiTracerDataForTotalDoseString[RadioTracer_NAME].str().c_str() << " "
                         << "InjectedActivity(Bq)=" << RadioTracerInjectedActivity[RadioTracer_NAME] << " "
                         << Quantity_NAME_UNIT << "_Total"<< "=" << TotalDoseFromRadioTracer[RadioTracer_NAME][Quantity_NAME] << " "
                            ;

                    OutS << "\n";

                    G4String headerText= OutS.str().c_str();


                    file << headerText.c_str();
                    G4cout << headerText.c_str();

                    G4cout << std::setw(24) << std::left << "# Volume" << " "
                           << std::setw(SZ) << std::left << Quantity_NAME_UNIT.c_str() << " "
                           << std::setw(SZ) << std::left << "SDev"  << " "
                           << std::setw(SZ) << std::left << "Rel_SDev(%)" << " "
                           << std::setw(SZ) << std::left << "Values Num" << " "
                           << std::setw(SZ) << std::left << "Mass[kg]" << " "
                           << std::setw(SZ) << std::left << "Volume[cm3]" << " "
                           << std::setw(SZ) << std::left << "Density[g/cm3] " << " \n";

                    file << std::setw(24) << std::left << "# Volume" << " "
                         << std::setw(SZ) << std::left << Quantity_NAME_UNIT.c_str() << " "
                         << std::setw(SZ) << std::left << "SDev"  << " "
                         << std::setw(SZ) << std::left << "Rel_SDev(%)" << " "
                         << std::setw(SZ) << std::left << "Values Num" << " "
                         << std::setw(SZ) << std::left << "Mass[kg]" << " "
                         << std::setw(SZ) << std::left << "Volume[cm3]" << " "
                         << std::setw(SZ) << std::left << "Density[g/cm3] " << " \n";

                    for ( auto it3 = it1->second.begin(); it3 != it1->second.end(); ++it3  ){
                        Target_NAME = it3->first;
                        if(Target_NAME == "World"){continue;}

                        double sdv = std::sqrt(std::abs(RadioTracerQuantityOrganVAR[RadioTracer_NAME]["AE"][Target_NAME]));
                        unsigned long long int sm = NumberOfStepsInRadiotracerOrgan[RadioTracer_NAME][Target_NAME];
                        double rsdv = ((sdv*sm)/RadioTracerQuantityOrganAETotalForRSD[RadioTracer_NAME]["AE"][Target_NAME])*100;

                        G4cout << std::setw(24) << std::left << Target_NAME.c_str() << " "
                               << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << it3->second << " "
                               << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << sdv << " "
                               << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << rsdv << " "
                               << std::resetiosflags(G4cout.basefield) << std::resetiosflags( G4cout.floatfield) << std::resetiosflags( G4cout.flags())
                               << std::setw(SZ) << std::left << sm << " "
                               << std::setw(SZ) << std::left << VolumeNameMassMap[Target_NAME] << " "
                               << std::setw(SZ) << std::left << VolumeNameVolumeMap[Target_NAME] << " "
                               << std::setw(SZ) << std::left << VolumeNameDensityMap[Target_NAME] << " \n";

                        file << std::setw(24) << std::left << Target_NAME.c_str() << " "
                             << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << it3->second << " "
                             << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << sdv << " "
                             << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << rsdv << " "
                             << std::resetiosflags(file.basefield) << std::resetiosflags( file.floatfield) << std::resetiosflags( file.flags())
                             << std::setw(SZ) << std::left << sm << " "
                             << std::setw(SZ) << std::left << VolumeNameMassMap[Target_NAME] << " "
                             << std::setw(SZ) << std::left << VolumeNameVolumeMap[Target_NAME] << " "
                             << std::setw(SZ) << std::left << VolumeNameDensityMap[Target_NAME] << " \n";

                    }

                    G4cout << std::setw(24) << std::left << "Total" << " "
                           << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << TotalDoseFromRadioTracer[RadioTracer_NAME][Quantity_NAME] << " "
                           << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << 1 << " "
                           << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << 1 << " "
                           << std::resetiosflags(G4cout.basefield) << std::resetiosflags( G4cout.floatfield) << std::resetiosflags( G4cout.flags())
                           << std::setw(SZ) << std::left << 1 << " "
                           << std::setw(SZ) << std::left << 1 << " "
                           << std::setw(SZ) << std::left << 1 << " "
                           << std::setw(SZ) << std::left << 1 << " \n";

                    file << std::setw(24) << std::left << "Total" << " "
                         << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << TotalDoseFromRadioTracer[RadioTracer_NAME][Quantity_NAME] << " "
                         << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << 1 << " "
                         << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << 1 << " "
                         << std::resetiosflags(file.basefield) << std::resetiosflags( file.floatfield) << std::resetiosflags( file.flags())
                         << std::setw(SZ) << std::left << 1 << " "
                         << std::setw(SZ) << std::left << 1 << " "
                         << std::setw(SZ) << std::left << 1 << " "
                         << std::setw(SZ) << std::left << 1 << " \n";


                    G4cout << "* ----------------------------------------------------------------------------------------------------------------------------------------------" << G4endl;
                    file << "* ----------------------------------------------------------------------------------------------------------------------------------------------\n";

                    file.close();
                }
            }
        }
    }

}
void G4TResultCalculation::GenerateRegionResultForRadioTracer(){

    if(V)G4cout << "\n-------- " << __FUNCTION__ << "\n" << G4endl ;

    G4String RadioTracer_NAME;
    G4String Quantity_NAME;
    G4String Particle_NAME;
    G4String Source_NAME;
    G4String Target_NAME;

    G4double Energy_Val;

    G4String Quantity_NAME_UNIT;

    // calculation of radiotracer quantity for src-trg combination and for organ
    for ( auto it = RadioTracerEnergyPerCent.begin(); it != RadioTracerEnergyPerCent.end(); ++it  ){

        RadioTracer_NAME = it->first;

        RadiTracerParticleEnergyDataString[RadioTracer_NAME].str("");
        RadiTracerParticleEnergyDataString[RadioTracer_NAME].clear();
        RadiTracerDataForTotalDoseString[RadioTracer_NAME].str("");
        RadiTracerDataForTotalDoseString[RadioTracer_NAME].clear();

        G4cout << "Generating Results for: " << RadioTracer_NAME <<G4endl;

        //G4cout << "RadioTracer_NAME " << RadioTracer_NAME <<G4endl;

        for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2  ){
            Particle_NAME = it2->first;
            //G4cout << "Particle_NAME " << Particle_NAME <<G4endl;

            for ( auto it3 = it2->second.begin(); it3 != it2->second.end(); ++it3  ){
                Energy_Val = it3->first;
                //G4cout << "Energy_Val " << Energy_Val <<G4endl;

                if(ResultTable["AE"][Particle_NAME][Energy_Val].size() == 0){
                    //G4cout << "\n-----> For Radiotracer " << RadioTracer_NAME
                    //       << ", the data for particle "<< Particle_NAME << " with energy "<< Energy_Val
                    //       << " are not found (this source configuration not simulated), the related quantities values will be generated by interpolating the already existed data of "
                    //       << Particle_NAME << " with energies surround the "
                    //       << Energy_Val << " from the existing source regions." << G4endl;
                    if(RadiotracerDataFomFile == false){
                        RadiTracerParticleEnergyDataString[RadioTracer_NAME] << " | (by interpolation) "<< Particle_NAME << " " << Energy_Val << "MeV " << RadioTracerEnergyPerCent[RadioTracer_NAME][Particle_NAME][Energy_Val] << "%";
                        RadiTracerDataForTotalDoseString[RadioTracer_NAME] << " | (by interpolation) "<< Particle_NAME << " " << Energy_Val << "MeV " << RadioTracerEnergyPerCent[RadioTracer_NAME][Particle_NAME][Energy_Val] << "%";
                    }
                    GenerateRadiotracerQuantitiesByInterpolation(Particle_NAME,Energy_Val);
                }else{
                    if(RadiotracerDataFomFile == false){
                        RadiTracerParticleEnergyDataString[RadioTracer_NAME] << " | (by simulation) "<< Particle_NAME << " " << Energy_Val << "MeV " << RadioTracerEnergyPerCent[RadioTracer_NAME][Particle_NAME][Energy_Val] << "%";
                        RadiTracerDataForTotalDoseString[RadioTracer_NAME] << " | (by simulation) "<< Particle_NAME << " " << Energy_Val << "MeV " << RadioTracerEnergyPerCent[RadioTracer_NAME][Particle_NAME][Energy_Val] << "%";
                    }
                }

                TotalRadioTracerEmmitedEnergy_Val[RadioTracer_NAME] += Energy_Val;

                // this is to provide results for all sources that can user simulate without taking into account the RadioTracer source Organs
                for ( auto DD = ResultTable["S"][Particle_NAME][Energy_Val].begin(); DD != ResultTable["S"][Particle_NAME][Energy_Val].end(); ++DD  ){
                    Source_NAME  = DD->first;
                    //G4cout << "Source_NAME " << Source_NAME <<G4endl;

                    for ( auto XX = RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracer_NAME].begin(); XX != RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracer_NAME].end(); ++XX  ){
                        G4String SrcOrg = XX->first;
                        if(Source_NAME == SrcOrg){
                            if(RadiotracerDataFomFile == false){
                                RadiTracerDataForTotalDoseString[RadioTracer_NAME] << " " << Source_NAME <<" As/A0(s)=" << RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracer_NAME][Source_NAME] ;
                            }
                        }
                    }

                    for ( auto CC = DD->second.begin(); CC != DD->second.end(); ++CC  ){
                        Target_NAME  = CC->first;

                        double vb= ResultTable["AE"][Particle_NAME][Energy_Val][Source_NAME][Target_NAME];
                        if(__isinf(vb) || __isnan(vb) || vb == 0.){continue;}

                        G4double RadiationPerCent = RadioTracerEnergyPerCent[RadioTracer_NAME][Particle_NAME][Energy_Val]/100;

                        G4double CumulatedActivityInSource = RadioTracerInjectedActivity[RadioTracer_NAME]
                                *RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracer_NAME][Source_NAME];

                        // values and variance calculation for chosen quantities
                        for (int rr = 0 ; rr < DoseCalcsQuantities.size() ; rr++) {

                            double ccc = ResultTable[DoseCalcsQuantities[rr]][Particle_NAME][Energy_Val][Source_NAME][Target_NAME];
                            if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                                RadioTracerQuantitySourceTargetValue[RadioTracer_NAME][DoseCalcsQuantities[rr]][Source_NAME][Target_NAME] +=
                                        RadiationPerCent
                                        *ccc;


                                // this condition just to show the regions with Wt different to 0 for AD in intake
                                if( DoseCalcsQuantities[rr] == "AD" && (__isnan(TissueFactorMap[Target_NAME])
                                                                        || __isinf(TissueFactorMap[Target_NAME])
                                                                        || TissueFactorMap[Target_NAME] == 0
                                                                        || TissueFactorMap[Target_NAME] == NULL)){

                                }else{
                                    RadioTracerQuantityOrganValue[RadioTracer_NAME][DoseCalcsQuantities[rr]][Target_NAME] +=
                                            RadiationPerCent
                                            *CumulatedActivityInSource
                                            *ccc;
                                }
                            }

                        }

                        if( !__isnan(NumberOfSteps[Particle_NAME][Energy_Val][Source_NAME][Target_NAME]) && !__isinf(NumberOfSteps[Particle_NAME][Energy_Val][Source_NAME][Target_NAME]) && NumberOfSteps[Particle_NAME][Energy_Val][Source_NAME][Target_NAME] != 0 && NumberOfSteps[Particle_NAME][Energy_Val][Source_NAME][Target_NAME] != NULL){

                            NumberOfStepsInRadiotracerSourceTarget[RadioTracer_NAME][Source_NAME][Target_NAME] +=
                                    NumberOfSteps[Particle_NAME][Energy_Val][Source_NAME][Target_NAME] ;

                            NumberOfStepsInRadiotracerOrgan[RadioTracer_NAME][Target_NAME] +=
                                    NumberOfSteps[Particle_NAME][Energy_Val][Source_NAME][Target_NAME] ;

                        }

                        double ccc = TotalAEForRadiotracerRSD["AE"][Particle_NAME][Energy_Val][Source_NAME][Target_NAME];

                        if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                            RadioTracerQuantitySourceTargetAETotalForRSD[RadioTracer_NAME]["AE"][Source_NAME][Target_NAME] +=
                                    ccc;

                            //G4cout << " Total for " << Particle_NAME << " "<< Energy_Val << " " << Source_NAME << " "<< Target_NAME << " is " << ccc << G4endl;
                            RadioTracerQuantityOrganAETotalForRSD[RadioTracer_NAME]["AE"][Target_NAME] +=
                                    ccc;
                        }

                        ccc = StandardDeviation["AE"][Particle_NAME][Energy_Val][Source_NAME][Target_NAME];
                        //G4cout << " SD for AE " << Particle_NAME << " "<< Energy_Val << " " << Source_NAME << " "<< Target_NAME << " is " << ccc << G4endl;

                        if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                            RadioTracerQuantitySourceTargetVariance[RadioTracer_NAME]["AE"][Source_NAME][Target_NAME] +=
                                    ((RadiationPerCent/100)*ccc)*((RadiationPerCent/100)*ccc);

                            //G4cout << " SD for AE " << Particle_NAME << " "<< Energy_Val << " " << Source_NAME << " "<< Target_NAME << " is " << ccc
                            //<< " . and Total " << RadioTracerQuantitySourceTargetVariance[RadioTracer_NAME]["AE"][Source_NAME][Target_NAME] << G4endl ;

                            RadioTracerQuantityOrganVAR[RadioTracer_NAME]["AE"][Target_NAME] +=
                                    ((RadiationPerCent/100)*ccc)*((RadiationPerCent/100)*ccc);
                        }
                    }
                }
            }
        }
        RadiTracerParticleEnergyDataString[RadioTracer_NAME] << " | ";
        RadiTracerDataForTotalDoseString[RadioTracer_NAME] << " | ";
    }


    // to calculate quantity ratio ER for each source-target
    for ( auto DoseCalcsQuantities = RadioTracerQuantitySourceTargetValue.begin(); DoseCalcsQuantities != RadioTracerQuantitySourceTargetValue.end(); ++DoseCalcsQuantities  ){
        RadioTracer_NAME = DoseCalcsQuantities->first;
        for ( auto it3 = DoseCalcsQuantities->second.begin(); it3 != DoseCalcsQuantities->second.end(); ++it3  ){
            Quantity_NAME = it3->first;
            for ( auto it4 = it3->second.begin(); it4 != it3->second.end(); ++it4  ){
                Source_NAME = it4->first;
                G4double CumulatedActivityInSource = (RadioTracerInjectedActivity[RadioTracer_NAME]
                                                      *RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracer_NAME][Source_NAME]);
                TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME] = 0;

                long double val = 0. ;
                for ( auto it5 = it4->second.begin(); it5 != it4->second.end(); ++it5  ){
                    Target_NAME = it5->first;

                    // to calculate ER for each source-target
                    if(Quantity_NAME == "E" || Quantity_NAME == "AD"){
                        double ccc = RadioTracerQuantitySourceTargetValue[RadioTracer_NAME][Quantity_NAME][Source_NAME][Target_NAME];
                        if( !__isnan(ccc) && !__isinf(ccc) && ccc != NULL){
                            TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME] = TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME] + ccc;
                        }
                    }
                    if( !__isnan(it5->second) && !__isinf(it5->second) && it5->second != 0 && it5->second != NULL){
                        val = val + it5->second;
                        //G4cout <<  " DoseCalcsQuantities->first " << DoseCalcsQuantities->first << " it3->first " << it3->first << " it4->first " << it4->first << " it5->first " << it5->first <<  " it5->second " << it5->second << " val: " << val << G4endl ;
                    }
                }
                // to calculate ER for each source-target
                if(Quantity_NAME == "E" || Quantity_NAME == "AD"){
                    if( !__isnan(TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME]) && !__isinf(TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME]) && TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME] != 0 && TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME] != NULL){
                        //double val = 0;
                        for ( auto it5 = it4->second.begin(); it5 != it4->second.end(); ++it5  ){
                            Target_NAME = it5->first;
                            if(Quantity_NAME == "E"){
                                double ccc = RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["E"][Source_NAME][Target_NAME];
                                if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                                    RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["ER"][Source_NAME][Target_NAME] =
                                            (RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["E"][Source_NAME][Target_NAME]/TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME])*100;
                                    //val = val + RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["ER"][Source_NAME][Target_NAME];
                                    //G4cout << " Source_NAME: " << Source_NAME << " Target_NAME: " << Target_NAME << " val: " << val << G4endl ;
                                }
                            }
                            if(Quantity_NAME == "AD"){
                                double ccc = RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["AD"][Source_NAME][Target_NAME];
                                if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                                    RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["DR"][Source_NAME][Target_NAME] =
                                            (RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["AD"][Source_NAME][Target_NAME]/TotalValueOfQuantity[RadioTracer_NAME][Quantity_NAME])*100;
                                    //val = val + RadioTracerQuantitySourceTargetValue[RadioTracer_NAME]["ER"][Source_NAME][Target_NAME];
                                    //G4cout << " Source_NAME: " << Source_NAME << " Target_NAME: " << Target_NAME << " val: " << val << G4endl ;
                                }
                            }
                        }
                    }
                }
                TotalDoseFromRadioTracer[RadioTracer_NAME][Quantity_NAME] += val * CumulatedActivityInSource;
                //G4cout << " CumulatedActivityInSource " << CumulatedActivityInSource << " val: " << val << " TotalDoseFromRadioTracer[RadioTracer_NAME][Quantity_NAME] " << TotalDoseFromRadioTracer[RadioTracer_NAME][Quantity_NAME] << G4endl ;
            }
        }
    }


    // to calculate quantity ratio ER for each organ
    for ( auto DoseCalcsQuantities = RadioTracerQuantityOrganValue.begin(); DoseCalcsQuantities != RadioTracerQuantityOrganValue.end(); ++DoseCalcsQuantities  ){

        RadioTracer_NAME = DoseCalcsQuantities->first;
        TotalValueOfQuantity[RadioTracer_NAME]["E"] = 0;
        TotalValueOfQuantity[RadioTracer_NAME]["AD"] = 0;

        //For ER
        for ( auto it4 = DoseCalcsQuantities->second["E"].begin(); it4 != DoseCalcsQuantities->second["E"].end(); ++it4  ){
            Target_NAME = it4->first;
            double ccc = it4->second;
            if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                TotalValueOfQuantity[RadioTracer_NAME]["E"] = TotalValueOfQuantity[RadioTracer_NAME]["E"] + ccc;
            }
        }
        //double val;
        for ( auto it4 = DoseCalcsQuantities->second["E"].begin(); it4 != DoseCalcsQuantities->second["E"].end(); ++it4  ){
            Target_NAME = it4->first;
            double ccc = it4->second;
            if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){

                RadioTracerQuantityOrganValue[RadioTracer_NAME]["ER"][Target_NAME] = (ccc/TotalValueOfQuantity[RadioTracer_NAME]["E"])*100;
                //val = val + RadioTracerQuantityOrganValue[RadioTracer_NAME]["ER"][Target_NAME];
                //G4cout << " Target_NAME: " << Target_NAME << " val: " << val << G4endl ;
            }
        }

        // For DR
        for ( auto it4 = DoseCalcsQuantities->second["AD"].begin(); it4 != DoseCalcsQuantities->second["AD"].end(); ++it4  ){
            Target_NAME = it4->first;
            double ccc = it4->second;
            if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){
                TotalValueOfQuantity[RadioTracer_NAME]["AD"] = TotalValueOfQuantity[RadioTracer_NAME]["AD"] + ccc;
            }
        }
        //double val;
        for ( auto it4 = DoseCalcsQuantities->second["AD"].begin(); it4 != DoseCalcsQuantities->second["AD"].end(); ++it4  ){
            Target_NAME = it4->first;
            double ccc = it4->second;
            if( !__isnan(ccc) && !__isinf(ccc) && ccc != 0 && ccc != NULL){

                RadioTracerQuantityOrganValue[RadioTracer_NAME]["DR"][Target_NAME] = (ccc/TotalValueOfQuantity[RadioTracer_NAME]["AD"])*100;
                //val = val + RadioTracerQuantityOrganValue[RadioTracer_NAME]["ER"][Target_NAME];
                //G4cout << " Target_NAME: " << Target_NAME << " val: " << val << G4endl ;
            }
        }
    }


    if(IsAllTargetsToScore == true){ TargetNamesToScore.clear();for (G4int gg = 0 ; gg < OrgansNameVector.size() ; gg++) {TargetNamesToScore.push_back(OrgansNameVector[gg]);}}

    for ( auto it0 = RadioTracerQuantitySourceTargetValue.begin(); it0 != RadioTracerQuantitySourceTargetValue.end(); ++it0  ){

        RadioTracer_NAME = it0->first;

        for ( auto it1 = it0->second.begin(); it1 != it0->second.end(); ++it1  ){

            Quantity_NAME = it1->first;

            // to write just the results of chosen quantities
            bool IsIn = false;
            for (int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {if(QuantityNamesToScore[gg] == Quantity_NAME){IsIn = true; break; }}
            if(IsIn == false){continue;}


            if(Quantity_NAME == "AF"){
                Quantity_NAME_UNIT = "AF[" + AFUnit + "]";
            }
            else if(Quantity_NAME == "SAF"){
                Quantity_NAME_UNIT = "SAF[" + SAFUnit + "]";
            }
            else if(Quantity_NAME == "S"){
                Quantity_NAME_UNIT = "S[(" + SUnit+ ")/"+UnitPerRadionuclideDecay+"]";
            }
            else if(Quantity_NAME == "AD"){
                Quantity_NAME_UNIT = "AD[(" + ADUnit+ ")/"+UnitPerRadionuclideDecay+"]";
            }
            else if(Quantity_NAME == "H"){
                Quantity_NAME_UNIT = "H[(" + HUnit+ ")/"+UnitPerRadionuclideDecay+"]";
            }
            else if(Quantity_NAME == "E"){
                Quantity_NAME_UNIT = "E[(" + EUnit+ ")/"+UnitPerRadionuclideDecay+"]";
            }
            else if(Quantity_NAME == "ER"){
                Quantity_NAME_UNIT = "ER[\%/"+UnitPerRadionuclideDecay+"]";
            }

            for ( auto it2 = it1->second.begin(); it2 != it1->second.end(); ++it2  ){
                Source_NAME = it2->first;

                G4cout << "\n\n================================== Generation of " << Quantity_NAME_UNIT << " deposition data in body target organs for intake of radio-tracer: " << RadioTracer_NAME << " into simulated source organ: " << Source_NAME << " ==================================\n" << G4endl ;

                //bool IsIn = false;
                //for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                //    if(TargetNamesToScore[gg] == Source_NAME){IsIn = true; break; }
                //}
                //if(IsIn == false){continue;}

                std::ostringstream filname ;
                filname << ResultDirectoryPath << "/ResultsData";
                G4String fileNameString = filname.str().c_str();
                std::ofstream file(fileNameString.c_str(), std::ios_base::app);

                //if(FILE* file = fopen(fileNameString.c_str(),"a")){
                if(file.is_open()){

                    if(V)G4cout << "Creating file  " << fileNameString.c_str() << " - writing the data output : \n"<< G4endl;

                    std::ostringstream OutS;
                    OutS << "****** "
                         << Quantity_NAME << " "
                         << Source_NAME << " "
                         << "RadioTracer "
                         << RadioTracer_NAME << " "
                         << GeometrySymbol <<  " "
                         << Physics << " "
                         << TotalRadioTracerEmmitedEnergy_Val[RadioTracer_NAME] << " "
                         << RadiTracerParticleEnergyDataString[RadioTracer_NAME].str().c_str() << " ";
                    OutS << "\n";

                    G4String headerText= OutS.str().c_str();

                    file << headerText.c_str();
                    G4cout << headerText.c_str();

                    G4cout << std::setw(24) << std::left << "# Volume" << " "
                           << std::setw(SZ) << std::left << Quantity_NAME_UNIT.c_str() << " "
                           << std::setw(SZ) << std::left << "SDev"  << " "
                           << std::setw(SZ) << std::left << "Rel_SDev(%)" << " "
                           << std::setw(SZ) << std::left << "Values Num" << " "
                           << std::setw(SZ) << std::left << "Mass[kg]" << " "
                           << std::setw(SZ) << std::left << "Volume[cm3]" << " "
                           << std::setw(SZ) << std::left << "Density[g/cm3] " << " \n";

                    file << std::setw(24) << std::left << "# Volume" << " "
                         << std::setw(SZ) << std::left << Quantity_NAME_UNIT.c_str() << " "
                         << std::setw(SZ) << std::left << "SDev"  << " "
                         << std::setw(SZ) << std::left << "Rel_SDev(%)" << " "
                         << std::setw(SZ) << std::left << "Values Num" << " "
                         << std::setw(SZ) << std::left << "Mass[kg]" << " "
                         << std::setw(SZ) << std::left << "Volume[cm3]" << " "
                         << std::setw(SZ) << std::left << "Density[g/cm3] " << " \n";

                    for ( auto it3 = it2->second.begin(); it3 != it2->second.end(); ++it3  ){
                        Target_NAME = it3->first;
                        if(Target_NAME == "World"){continue;}

                        bool IsIn = false;
                        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                            if(TargetNamesToScore[gg] == Target_NAME){IsIn = true; break; }
                        }
                        if(IsIn == true){

                            double sdv = std::sqrt(std::abs(RadioTracerQuantitySourceTargetVariance[RadioTracer_NAME]["AE"][Source_NAME][Target_NAME]));
                            unsigned long long int sm = NumberOfStepsInRadiotracerSourceTarget[RadioTracer_NAME][Source_NAME][Target_NAME];
                            double rsdv = ((sdv * sm)/ RadioTracerQuantitySourceTargetAETotalForRSD[RadioTracer_NAME]["AE"][Source_NAME][Target_NAME])*100;

                            //G4cout << " Total for 100*" << sdv << "*"<< sm << "/" << RadioTracerQuantitySourceTargetAETotalForRSD[RadioTracer_NAME]["AE"][Source_NAME][Target_NAME] << "="<< rsdv << G4endl;

                            G4cout << std::setw(24) << std::left << Target_NAME.c_str() << " "
                                   << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << it3->second << " "
                                   << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << sdv << " "
                                   << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << rsdv << " "
                                   << std::resetiosflags(G4cout.basefield) << std::resetiosflags( G4cout.floatfield) << std::resetiosflags( G4cout.flags())
                                   << std::setw(SZ) << std::left << sm << " "
                                   << std::setw(SZ) << std::left << VolumeNameMassMap[Target_NAME] << " "
                                   << std::setw(SZ) << std::left << VolumeNameVolumeMap[Target_NAME] << " "
                                   << std::setw(SZ) << std::left << VolumeNameDensityMap[Target_NAME] << " \n";

                            file << std::setw(24) << std::left << Target_NAME.c_str() << " "
                                 << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << it3->second << " "
                                 << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << sdv << " "
                                 << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << rsdv << " "
                                 << std::resetiosflags(file.basefield) << std::resetiosflags( file.floatfield) << std::resetiosflags( file.flags())
                                 << std::setw(SZ) << std::left << sm << " "
                                 << std::setw(SZ) << std::left << VolumeNameMassMap[Target_NAME] << " "
                                 << std::setw(SZ) << std::left << VolumeNameVolumeMap[Target_NAME] << " "
                                 << std::setw(SZ) << std::left << VolumeNameDensityMap[Target_NAME] << " \n";

                        }
                        else{continue;}
                    }
                    //std::fprintf(file,"* ----------------------------------------------------------------------------------------------------------------------------------------------\n");
                    //std::fprintf(stdout,"* ----------------------------------------------------------------------------------------------------------------------------------------------\n");

                    G4cout << "* ----------------------------------------------------------------------------------------------------------------------------------------------" << G4endl;
                    file << "* ----------------------------------------------------------------------------------------------------------------------------------------------\n";

                    file.close();

                }
            }
        }
    }



    // i want just the results of radionuclides for articles
    return;









    if(GenerateResultsForRadioTracerExams){

        for ( auto it0 = RadioTracerQuantityOrganValue.begin(); it0 != RadioTracerQuantityOrganValue.end(); ++it0  ){

            RadioTracer_NAME = it0->first;

            if(RadioTracerInjectedActivity[RadioTracer_NAME] == 0.){
                continue;
            }

            for ( auto it1 = it0->second.begin(); it1 != it0->second.end(); ++it1  ){

                Quantity_NAME = it1->first;

                // to write just the results of chosen quantities
                bool IsIn = false;
                for (int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {if(QuantityNamesToScore[gg] == Quantity_NAME){IsIn = true; break; }}
                if(IsIn == false){continue;}


                if(Quantity_NAME == "AF"){
                    Quantity_NAME_UNIT = "AF[" + AFUnit + "]";
                }
                else if(Quantity_NAME == "SAF"){
                    Quantity_NAME_UNIT = "SAF[" + SAFUnit + "]";
                }
                else if(Quantity_NAME == "S"){
                    Quantity_NAME_UNIT = "S[(" + SUnit+ ")/"+UnitPerRadioTracerDecay+"]";
                }
                else if(Quantity_NAME == "AD"){
                    Quantity_NAME_UNIT = "AD[(" + ADUnit+ ")/"+UnitPerRadioTracerDecay+"]";
                }
                else if(Quantity_NAME == "H"){
                    Quantity_NAME_UNIT = "H[(" + HUnit+ ")/"+UnitPerRadioTracerDecay+"]";
                }
                else if(Quantity_NAME == "E"){
                    Quantity_NAME_UNIT = "E[(" + EUnit+ ")/"+UnitPerRadioTracerDecay+"]";
                }
                else if(Quantity_NAME == "ER"){
                    Quantity_NAME_UNIT = "ER[\%/"+UnitPerRadioTracerDecay+"]";
                }
                else if(Quantity_NAME == "DR"){
                    Quantity_NAME_UNIT = "DR[\%/"+UnitPerRadioTracerDecay+"]";
                }

                G4cout << "\n\n================================== Generation of " << Quantity_NAME_UNIT << " deposition data in all simulated body organs for " << RadioTracerInjectedActivity[RadioTracer_NAME] << " Bq intake of radio-tracer: " << RadioTracer_NAME << " into " << GeometrySymbol<< " simulated body ==================================\n" << G4endl ;

                std::ostringstream filname ;
                filname << ResultDirectoryPath << "/ResultsData";
                G4String fileNameString = filname.str().c_str();
                std::ofstream file(fileNameString.c_str(), std::ios_base::app);

                //if(FILE* file = fopen(fileNameString.c_str(),"a")){
                if(file.is_open()){

                    if(V)G4cout << "Creating file  " << fileNameString.c_str() << " - writing the data output : \n"<< G4endl;

                    std::ostringstream OutS;
                    OutS << "****** "
                         << Quantity_NAME << " "
                         << "IntakeIntoBody "
                         << "RadioTracer "
                         << RadioTracer_NAME << " "
                         << GeometrySymbol <<  " "
                         << Physics << " "
                         << TotalRadioTracerEmmitedEnergy_Val[RadioTracer_NAME] << " "
                         << RadiTracerDataForTotalDoseString[RadioTracer_NAME].str().c_str() << " "
                         << "InjectedActivity(Bq)=" << RadioTracerInjectedActivity[RadioTracer_NAME] << " "
                         << Quantity_NAME_UNIT << "_Total"<< "=" << TotalDoseFromRadioTracer[RadioTracer_NAME][Quantity_NAME] << " "
                            ;

                    OutS << "\n";

                    G4String headerText= OutS.str().c_str();


                    file << headerText.c_str();
                    G4cout << headerText.c_str();

                    G4cout << std::setw(24) << std::left << "# Volume" << " "
                           << std::setw(SZ) << std::left << Quantity_NAME_UNIT.c_str() << " "
                           << std::setw(SZ) << std::left << "SDev"  << " "
                           << std::setw(SZ) << std::left << "Rel_SDev(%)" << " "
                           << std::setw(SZ) << std::left << "Values Num" << " "
                           << std::setw(SZ) << std::left << "Mass[kg]" << " "
                           << std::setw(SZ) << std::left << "Volume[cm3]" << " "
                           << std::setw(SZ) << std::left << "Density[g/cm3] " << " \n";

                    file << std::setw(24) << std::left << "# Volume" << " "
                         << std::setw(SZ) << std::left << Quantity_NAME_UNIT.c_str() << " "
                         << std::setw(SZ) << std::left << "SDev"  << " "
                         << std::setw(SZ) << std::left << "Rel_SDev(%)" << " "
                         << std::setw(SZ) << std::left << "Values Num" << " "
                         << std::setw(SZ) << std::left << "Mass[kg]" << " "
                         << std::setw(SZ) << std::left << "Volume[cm3]" << " "
                         << std::setw(SZ) << std::left << "Density[g/cm3] " << " \n";

                    for ( auto it3 = it1->second.begin(); it3 != it1->second.end(); ++it3  ){
                        Target_NAME = it3->first;
                        if(Target_NAME == "World"){continue;}

                        double sdv = std::sqrt(std::abs(RadioTracerQuantityOrganVAR[RadioTracer_NAME]["AE"][Target_NAME]));
                        unsigned long long int sm = NumberOfStepsInRadiotracerOrgan[RadioTracer_NAME][Target_NAME];
                        double rsdv = ((sdv*sm)/RadioTracerQuantityOrganAETotalForRSD[RadioTracer_NAME]["AE"][Target_NAME])*100;

                        G4cout << std::setw(24) << std::left << Target_NAME.c_str() << " "
                               << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << it3->second << " "
                               << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << sdv << " "
                               << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << rsdv << " "
                               << std::resetiosflags(G4cout.basefield) << std::resetiosflags( G4cout.floatfield) << std::resetiosflags( G4cout.flags())
                               << std::setw(SZ) << std::left << sm << " "
                               << std::setw(SZ) << std::left << VolumeNameMassMap[Target_NAME] << " "
                               << std::setw(SZ) << std::left << VolumeNameVolumeMap[Target_NAME] << " "
                               << std::setw(SZ) << std::left << VolumeNameDensityMap[Target_NAME] << " \n";

                        file << std::setw(24) << std::left << Target_NAME.c_str() << " "
                             << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << it3->second << " "
                             << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << sdv << " "
                             << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << rsdv << " "
                             << std::resetiosflags(file.basefield) << std::resetiosflags( file.floatfield) << std::resetiosflags( file.flags())
                             << std::setw(SZ) << std::left << sm << " "
                             << std::setw(SZ) << std::left << VolumeNameMassMap[Target_NAME] << " "
                             << std::setw(SZ) << std::left << VolumeNameVolumeMap[Target_NAME] << " "
                             << std::setw(SZ) << std::left << VolumeNameDensityMap[Target_NAME] << " \n";

                    }

                    G4cout << std::setw(24) << std::left << "Total" << " "
                           << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << TotalDoseFromRadioTracer[RadioTracer_NAME][Quantity_NAME] << " "
                           << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << 1 << " "
                           << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << 1 << " "
                           << std::resetiosflags(G4cout.basefield) << std::resetiosflags( G4cout.floatfield) << std::resetiosflags( G4cout.flags())
                           << std::setw(SZ) << std::left << 1 << " "
                           << std::setw(SZ) << std::left << 1 << " "
                           << std::setw(SZ) << std::left << 1 << " "
                           << std::setw(SZ) << std::left << 1 << " \n";

                    file << std::setw(24) << std::left << "Total" << " "
                         << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << TotalDoseFromRadioTracer[RadioTracer_NAME][Quantity_NAME] << " "
                         << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << 1 << " "
                         << std::setw(SZ) << std::left << std::scientific << std::setprecision(VarSZ) << 1 << " "
                         << std::resetiosflags(file.basefield) << std::resetiosflags( file.floatfield) << std::resetiosflags( file.flags())
                         << std::setw(SZ) << std::left << 1 << " "
                         << std::setw(SZ) << std::left << 1 << " "
                         << std::setw(SZ) << std::left << 1 << " "
                         << std::setw(SZ) << std::left << 1 << " \n";


                    G4cout << "* ----------------------------------------------------------------------------------------------------------------------------------------------" << G4endl;
                    file << "* ----------------------------------------------------------------------------------------------------------------------------------------------\n";

                    file.close();
                }
            }
        }
    }

}
void G4TResultCalculation::GenerateRadiotracerQuantitiesByInterpolation(G4String ParticleName, G4double Energy){
    
    if(V)G4cout << "\n-------- " << __FUNCTION__ << G4endl ;
    
    for ( auto Abeg = SourceParticleEnergyValues[GeometrySymbol].begin(); Abeg != SourceParticleEnergyValues[GeometrySymbol].end(); ++Abeg  )
    {
        for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  ){
            std::sort(Bbeg->second.begin(), Bbeg->second.end());
        }
    }
    
    G4double Energy1 = 0;
    G4double Energy2 = 0;
    
    //for(G4int gg = 0 ; gg < (G4int)OrgansNameVector.size() ; gg++){
    for ( auto Abeg = SourceParticleEnergyValues[GeometrySymbol].begin(); Abeg != SourceParticleEnergyValues[GeometrySymbol].end(); ++Abeg  ){
        
        for(G4int ss = 0 ; ss < Abeg->second[ParticleName].size() ; ss++){
            
            G4int ff = ss+1;
            
            G4int da = Abeg->second[ParticleName].size()-1; if(ff == da){break;}
            
            G4double E1 = Abeg->second[ParticleName][ss];
            G4double E2 = Abeg->second[ParticleName][ff];
            //G4cout << ss << " " << E1 << " " << ff << " " << E2 << G4endl ;
            
            if(E1 < Energy && Energy < E2){
                Energy1 = E1;
                Energy2 = E2;
                //G4cout << " Energy1 " << Energy1 << " Energy " << Energy << " Energy2 " << Energy2 << G4endl ;
                break;
            }
        }
        
        for(G4int hh = 0 ; hh < (G4int)OrgansNameVector.size() ; hh++){
            
            for(G4int vc = 0 ; vc < (G4int)DoseCalcsQuantities.size() ; vc++){
                
                G4double Val1 = ResultTable[DoseCalcsQuantities[vc]][ParticleName][Energy1][Abeg->first][OrgansNameVector[hh]];
                G4double Val2 = ResultTable[DoseCalcsQuantities[vc]][ParticleName][Energy2][Abeg->first][OrgansNameVector[hh]];
                //G4double Val  = Val1 + (Energy-Energy1)*((Val2-Val1)/(Energy2-Energy1));
                G4double Val  = std::exp(std::log(Val1) + (std::log(Energy) - std::log(Energy1)) * (std::log(Val2) - std::log(Val1)) / (std::log(Energy2) - std::log(Energy1)));

                ResultTable[DoseCalcsQuantities[vc]][ParticleName][Energy][Abeg->first][OrgansNameVector[hh]] = Val;

                Val1 = StandardDeviation["AE"][ParticleName][Energy1][Abeg->first][OrgansNameVector[hh]];
                Val2 = StandardDeviation["AE"][ParticleName][Energy2][Abeg->first][OrgansNameVector[hh]];
                //Val  = Val1 + (Energy-Energy1)*((Val2-Val1)/(Energy2-Energy1));
                Val  = std::exp(std::log(Val1) + (std::log(Energy) - std::log(Energy1)) * (std::log(Val2) - std::log(Val1)) / (std::log(Energy2) - std::log(Energy1)));

                StandardDeviation["AE"][ParticleName][Energy][Abeg->first][OrgansNameVector[hh]] = Val;

                Val1 = TotalAEForRadiotracerRSD["AE"][ParticleName][Energy1][Abeg->first][OrgansNameVector[hh]];
                Val2 = TotalAEForRadiotracerRSD["AE"][ParticleName][Energy2][Abeg->first][OrgansNameVector[hh]];
                //Val  = Val1 + (Energy-Energy1)*((Val2-Val1)/(Energy2-Energy1));
                Val  = std::exp(std::log(Val1) + (std::log(Energy) - std::log(Energy1)) * (std::log(Val2) - std::log(Val1)) / (std::log(Energy2) - std::log(Energy1)));
                TotalAEForRadiotracerRSD["AE"][ParticleName][Energy][Abeg->first][OrgansNameVector[hh]] = Val;

                unsigned long long int uVal1 = NumberOfSteps[ParticleName][Energy1][Abeg->first][OrgansNameVector[hh]];
                unsigned long long int uVal2 = NumberOfSteps[ParticleName][Energy2][Abeg->first][OrgansNameVector[hh]];
                //unsigned long long int uVal  = uVal1 + (Energy-Energy1)*((uVal2-uVal1)/(Energy2-Energy1));
                unsigned long long int uVal  = std::exp(std::log(uVal1) + (std::log(Energy) - std::log(Energy1)) * (std::log(uVal2) - std::log(uVal1)) / (std::log(Energy2) - std::log(Energy1)));
                NumberOfSteps[ParticleName][Energy][Abeg->first][OrgansNameVector[hh]] = uVal;

                //G4cout << "\n\n\n\n Variable " << DoseCalcsQuantities[vc]  << " ParticleName " << ParticleName  << " Energy1 " << Energy1 << " Energy2 " << Energy2   << " Abeg->first " << Abeg->first  << " OrgansNameVector[hh] " << OrgansNameVector[hh]  << " " << Val << " ResultTable[DoseCalcsQuantities[vc]][ParticleName][Energy][Abeg->first][OrgansNameVector[hh]] " << ResultTable[DoseCalcsQuantities[vc]][ParticleName][Energy][Abeg->first][OrgansNameVector[hh]] << " " <<  OrgansNameVector.size() << G4endl ;
                
                //if(Val1 != 0. && Val2 != 0.){
                //ResultTable[DoseCalcsQuantities[vc]][ParticleName][Energy][Abeg->first][OrgansNameVector[hh]] = Val1 + (Energy-Energy1)*((Val2-Val1)/(Energy2-Energy1));
                //}
            }
        }
    }
}

double G4TResultCalculation::GenerateRadiationFactor(G4String ParticleName, double Energy){ // enerrgy in MeV

    double factor = 1.;
    //double factor = RadiationFactorMap[ParticleName][Energy];

    if(ParticleName == "gamma"){
        factor = 1; // all energies
    }
    else if (ParticleName == "e-" || ParticleName == "e+"){
        factor = 1; // all energies
    }
    else if (ParticleName == "alpha"){
        factor = 20; // all energies for alpha and heavy nuclei
    }
    else if (ParticleName == "proton"){
        if(Energy <= 2.){
            factor = 1;
        }else{
            factor = 5;
        }
    }
    else if (ParticleName == "neutron"){
        if(Energy < 0.01){
            factor = 5;
        }
        else if( 0.01 <= Energy && Energy <= 0.1){
            factor = 10;
        }
        else if( 0.1 < Energy && Energy <= 2){
            factor = 20;
        }
        else if( 2 < Energy && Energy <= 10){
            factor = 20;
        }
        else if( 20 < Energy){
            factor = 5;
        }
    }
    return factor;
}
G4String G4TResultCalculation::RemoveWordFromString(G4String str, G4String word){
    
    if(V)G4cout << str << G4endl;
    
    if (str.find(word) != std::string::npos)
    {
        size_t p = -1;
        
        // To cover the case if the word is at the beginning of the string or anywhere in the middle
        std::string tempWord = word + " ";
        while ((p = str.find(word)) != std::string::npos)
            str.replace(p, tempWord.length(), "");
        
        // To cover the edge case if the word is at the end of the string
        tempWord = " " + word;
        while ((p = str.find(word)) != std::string::npos)
            str.replace(p, tempWord.length(), "");
    }
    if(V)G4cout << str << G4endl;
    return str;
}
