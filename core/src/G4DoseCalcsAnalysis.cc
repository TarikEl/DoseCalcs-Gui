

#include "G4DoseCalcsAnalysis.hh"
#include <stdio.h>
#include <filesystem>
#include <iostream>
//#include "G4PhysicalConstants.hh"
#include "G4RandomTools.hh"

#ifdef ANALYSIS_USE
#include "Riostream.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
//#include "TFile.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TPaletteAxis.h"
#include "TText.h"
#include "TF1.h"
#include "TFitResult.h"
#endif

extern std::string TissueFacFilePath;
extern std::string appBuildDir;
extern G4String ResultDirectoryPath ;
extern std::string MacrosStartingFile;
extern bool DirectoryExists(const char* pzPath );

G4DoseCalcsAnalysis::~G4DoseCalcsAnalysis(){}
G4DoseCalcsAnalysis::G4DoseCalcsAnalysis(std::string CurrentResultsDir){
    
    DefaultErrorDistance = 0.00000;
    MinValForLog = 1e+12;
    DiffSym = "RA";
    GeometrySymbol = "phantom0";
    MeV_to_J = 1.60218e-13;
    Gy_to_Sv = 1. ;
    Bq_to_MBq = 1e-6 ;
    
    GeometryVarDigitNum = 3;
    QuantityDigitNum = 4;
    DiffDigitNum = 2;
    SDevDigitNum = 4;
    
    ResultDirectoryPath = CurrentResultsDir ;
    GraphsDirectoryPath = ResultDirectoryPath+"/PlotsAndTables/";
    
    CombinedOutFileName = GraphsDirectoryPath+"AllOuts.root";
    
    UseLogE = "yes";
    UseLogVariable = "yes";
    UseGridXY = "yes";
    PrintTitle = "yes";
    LegendPos = "RightTop";
    LegendXWidth  = 0.15;
    LegendYHeight  = 0.23;
    
    AEUnitFactor = 1.;
    AFUnitFactor = 1.;
    SAFUnitFactor = 1.;
    SUnitFactor = 1.;
    ADUnitFactor = 1.;
    HUnitFactor = 1.;
    EUnitFactor = 1.;
    
    AEUnit = "MeV";
    AFUnit = "";
    SAFUnit = "kg-1";
    ADUnit = "MeV/kg";
    SUnit = "MeV/kg-Bq";
    HUnit = "MeV/kg";
    EUnit = "MeV/kg";
    
    ResidenceTimeUnitFactor = 1.;
    ResidenceTimeUnit = "s";
    
    AdministeredActivityUnitFactor = 1.;
    AdministeredActivityUnit = "MeV";
    
    RootC = {1,2,3,4,5,6,7,8,9,11,12,28,30,32,35,38,39,41, 42, 43, 46, 47, 49, 48, 49, 22, 23 };
    RootM = {53, 54, 55, 56, 59, 65, 67, 50, 51, 57, 58, 60, 61, 62, 63, 64, 66, 29, 43, 41, 39, 52};

    GenerateResultsForRadioTracer = false;
}

// read inputs for analysis
// called from main() in graph class
void G4DoseCalcsAnalysis::ReadSimulationData(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    //std::ostringstream filename1;
    //filename1 << ResultDirectoryPath <<"/SimData" ;
    //std::ifstream file1(filename1.str().c_str() , std::ios::binary);
    
    std::ifstream file1(MacrosStartingFile.c_str() , std::ios::binary);
    
    if(file1.is_open()){
        
        std::cout << "\nReading file " << MacrosStartingFile.c_str() << std::endl ;
        
        std::string ParameterName, line, word;
        //while(file1.peek()!=EOF){
        while (getline(file1, line)) {
            
            std::istringstream LineString(line);
            ParameterName = "";
            LineString >> ParameterName;
            
            if(ParameterName =="/GeometryData/setGeometrySymbol"){// to evaluate the last 4 lines
                LineString >> GeometrySymbol ;
                G4cout << " GeometrySymbol " << GeometrySymbol << G4endl ;
            }
            else if(ParameterName =="/PhysicsData/setCutsData"){// to evaluate the last 4 lines
                //LineString >> CutsEnergy  >> CutsDistance ;
                //G4cout << " CutsEnergy " << CutsEnergy << " CutsDistance " << CutsDistance << G4endl ;
            }
            else if(ParameterName =="/SourceData/setEventsParticleNameData"){// to evaluate the last 4 lines
                LineString >> ParticleName ;
                G4cout << " ParticleName " << ParticleName << G4endl ;
            }
            else if(ParameterName == "/SourceData/setEventsInitialMomDirData"){// to evaluate the last 4 lines
                LineString >> ParameterName >> MomDirDistribution;
                G4cout << " MomDirDistribution " << MomDirDistribution << G4endl ;
            }
            else if(ParameterName =="/SourceData/setEventsInitialEneData"){// to evaluate the last 4 lines
                LineString >> ParameterName >> EnergyDistribution ;//>> ParticleSourceEnergy ;
                G4cout << " EnergyDistribution " << EnergyDistribution << G4endl ;
                if(EnergyDistribution =="Rayleigh"){// to evaluate the last 4 lines
                    LineString >> RayleighEmax ;
                    std::cout << " RayleighEmax " << RayleighEmax << std::endl ;
                }
                else if(EnergyDistribution =="Uniform"){// to evaluate the last 4 lines
                    LineString >> UniformEmin >> UniformEmax;
                    std::cout << " UniformEmin " << UniformEmin<< " UniformEmax " << UniformEmax << std::endl ;
                }
                else if(EnergyDistribution =="Gauss"){// to evaluate the last 4 lines
                    LineString >> GaussSDev >> GaussMean ;
                    std::cout << " GaussSDev " << GaussSDev << " GaussMean " << GaussMean << std::endl ;
                }
                else if(EnergyDistribution =="Mono"){// to evaluate the last 4 lines
                    LineString >> MonoEnergy ;
                    std::cout << " MonoEnergy " << MonoEnergy << std::endl ;
                }
            }
            
            else if(ParameterName =="/SourceData/setEventsInitialPosData"){// to evaluate the last 4 lines
                LineString >> ParameterName >> SourceType >>  SourceRegionName ;
                G4cout << " SourceType " << SourceType << " SourceRegionName " << SourceRegionName << G4endl ;
            }
            else if(ParameterName =="/RunAndScoreData/setSimNumOnRanks"){// to evaluate the last 4 lines
                //LineString >> OneOrMultiSimulations ;
                //G4cout << " OneOrMultiSimulations " << OneOrMultiSimulations << G4endl ;
            }
            else if(ParameterName =="/RunAndScoreData/setQuantitiesToScore"){
                
                G4String wordsWithSpace ;
                while(LineString >> wordsWithSpace ){
                    if(wordsWithSpace.empty() || wordsWithSpace == ""){ }
                    else if (wordsWithSpace.find("All") != std::string::npos || wordsWithSpace.find("all") != std::string::npos) {
                        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n       found       --- " << Organs_To_Score << G4endl;
                        QuantityNamesToScore.push_back("AE");
                        QuantityNamesToScore.push_back("SAF");
                        QuantityNamesToScore.push_back("AF");
                        QuantityNamesToScore.push_back("AD");
                        QuantityNamesToScore.push_back("S");
                        QuantityNamesToScore.push_back("H");
                        QuantityNamesToScore.push_back("E");
                        QuantityNamesToScore.push_back("ER");
                        QuantityNamesToScore.push_back("DR");
                        break;
                    }
                    else {
                        QuantityNamesToScore.push_back(wordsWithSpace);
                    }
                }
                G4cout << "  Quantities names to score are : " << G4endl ;
                for (G4int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {
                    G4cout << QuantityNamesToScore[gg] << G4endl;
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
                                    
                                    //G4cout << " RegionName " << RegionName << G4endl ;
                                    
                                    if(SourceOrTarget == "source"){
                                        SourceNamesToScore.push_back(RegionName);
                                    }
                                    else if(SourceOrTarget == "target"){
                                        TargetNamesToScore.push_back(RegionName);
                                    }
                                }else{
                                    //G4cout << " item " << item << G4endl ;
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
                    std::sort(TargetNamesToScore.begin(), TargetNamesToScore.end());
                    std::sort(SourceNamesToScore.begin(), SourceNamesToScore.end());


                }
                
                if(V)G4cout << "IsAllSourcesToScore "<< IsAllSourcesToScore << "IsAllTargetsToScore "<< IsAllTargetsToScore << G4endl;
                
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
                
                GraphsDirectoryPath = ResultDirectoryPath+"/PlotsAndTables/";
                if(!DirectoryExists(GraphsDirectoryPath.c_str())){
                    std::filesystem::create_directory(GraphsDirectoryPath);
                }
                G4cout << " ResultDirectoryPath " << ResultDirectoryPath << " GraphsDirectoryPath " << GraphsDirectoryPath << G4endl ;
                
            }
            else if(ParameterName =="/RunAndScoreData/setRadioTracerData"){// to evaluate the last 4 lines
                
                LineString >> RadioTracerName ;
                if(RadioTracerName == "File"){
                    
                    RadiotracerDataFomFile = true;
                    LineString >> RadioTracerName ;
                    
                    G4String filename;
                    LineString >> filename ;
                    
                    std::ifstream file1(filename.c_str() , std::ios::binary);
                    if(file1.is_open()){
                        
                        G4String RPar;
                        G4String Word;
                        G4String S;
                        
                        while (file1 >> Word) {
                            
                            if(Word == "gamma" || Word == "e-" || Word == "e+" || Word == "alpha" || Word == "proton" || Word == "neutron"){
                                RPar = Word;
                                file1 >> S;
                                //std::cout<< " RPar = " << RPar << " S = " << S << std::endl;
                                continue;
                            }
                            
                            G4double Pval = atof(Word);
                            G4double EVal; file1 >> EVal;
                            RadioTracerEnergyPerCent[RadioTracerName][RPar][EVal] = Pval*100;

                            GenerateResultsForRadioTracer = true;
                            //std::cout << " RadioTracerName " << RadioTracerName << " RPar = " << RPar << " Pval = " << Pval << " EVal = " << EVal << std::endl;
                        }
                        file1.close();
                    }
                }
                else{
                    G4cout << " RadioTracerName " << RadioTracerName << G4endl ;
                    G4double PerC, Ene; G4String pn;
                    while(LineString >> pn ){
                        LineString >> Ene;
                        LineString >> PerC;
                        RadioTracerEnergyPerCent[RadioTracerName][pn][Ene] = PerC;
                        G4cout << " pn " << pn << " Ene " << Ene << " PerC " << PerC << G4endl ;

                        GenerateResultsForRadioTracer = true;
                    }
                }
            }
            else if(ParameterName =="/RunAndScoreData/setQuantitiesUnits"){// to evaluate the last 4 lines
                std::string Quant; std::string uni;
                
                while(LineString >> Quant ){
                    
                    LineString >> uni;
                    
                    if(Quant == "AE"){
                        AEUnit = uni;
                    }else if(Quant == "AF"){
                        AFUnit = uni;
                    }else if(Quant == "SAF"){
                        SAFUnit = uni;
                    }else if(Quant == "AD"){
                        ADUnit = uni;
                    }else if(Quant == "S"){
                        SUnit = uni;
                    }else if(Quant == "H"){
                        HUnit = uni;
                    }else if(Quant == "E"){
                        EUnit = uni;
                    }
                    std::cout << " Quantity " << Quant << " unit " << uni << std::endl ;
                }
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
                
                G4cout << " AdministeredActivityUnit " << AdministeredActivityUnit
                       << " AdministeredActivityUnitFactor " << AdministeredActivityUnitFactor
                       << " ResidenceTimeUnit " << ResidenceTimeUnit
                       << " ResidenceTimeUnitFactor " << ResidenceTimeUnitFactor
                       << G4endl ;
                
                InjectedActivity = InjectedActivity*AdministeredActivityUnitFactor;
                RadioTracerInjectedActivity[RadioTracerName] = InjectedActivity;
                G4cout << " RadioTracerName " << RadioTracerName << " InjectedActivity " << InjectedActivity << G4endl ;
                
                G4double Fs, Ti, ai, AsPerA0 ; G4String orgName;
                while(LineString >> orgName ){
                    
                    //LineString >> Fs;RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracerName][orgName].push_back(Fs);
                    //LineString >> Ti;RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracerName][orgName].push_back(Ti);
                    //LineString >> ai;RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracerName][orgName].push_back(ai);
                    LineString >> AsPerA0 ;
                    RadioTracerSourceTi_Fs_ai_AsPerA0[RadioTracerName][orgName]=AsPerA0*ResidenceTimeUnitFactor;
                    //" Fs " << Fs << " Ti " << Ti  << " ai " << ai  <<
                    G4cout << " orgName " << orgName << " AsPerA0 " << AsPerA0 << G4endl ;
                }
                GenerateResultsForRadioTracerExams = true;
                
            }
            else if(ParameterName =="/RunAndScoreData/setRadiationFactors"){// to evaluate the last 4 lines
                G4double Ene, Fac; G4String pn;
                
                while(LineString >> pn ){
                    LineString >> Ene;
                    LineString >> Fac;
                    RadiationFactorMap[pn][Ene] = Fac;
                    G4cout << " pn " << pn << " Ene " << Ene << " Fac " << Fac << G4endl ;
                }
            }
            else if(ParameterName =="/RunAndScoreData/setTissueFactors"){// to evaluate the last 4 lines
                G4double Fac; G4String pn;
                
                while(LineString >> pn ){
                    LineString >> Fac;
                    TissueFactorMap[pn] = Fac;
                    G4cout << " pn " << pn << " Fac " << Fac << G4endl ;
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
                        G4cout << " Quantity " << Quant << " unit " << uni << " AEUnitFactor " << AEUnitFactor << G4endl ;
                    }else if(Quant == "AF"){
                        if(uni == ""){
                            AFUnit = "";
                            AFUnitFactor = 1;
                        }
                        G4cout << " Quantity " << Quant << " unit " << uni << " AFUnitFactor " << AFUnitFactor << G4endl ;
                    }else if(Quant == "SAF"){
                        if(uni == "kg-1"){
                            SAFUnit = "kg-1";
                            SAFUnitFactor = 1;
                        }else if(uni == "g-1"){
                            SAFUnit = "g-1";
                            SAFUnitFactor = 1/1e+3;
                        }
                        G4cout << " Quantity " << Quant << " unit " << uni << " SAFUnitFactor " << SAFUnitFactor << G4endl ;
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
                        G4cout << " Quantity " << Quant << " unit " << uni << " ADUnitFactor " << ADUnitFactor << G4endl ;
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
                        }
                        G4cout << " Quantity " << Quant << " unit " << uni << " SUnitFactor " << SUnitFactor << G4endl ;
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
                        }
                        G4cout << " Quantity " << Quant << " unit " << uni << " HUnitFactor " << HUnitFactor << G4endl ;
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
                        }
                        G4cout << " Quantity " << Quant << " unit " << uni << " EUnitFactor " << EUnitFactor << G4endl ;
                    }
                }
            }
            else if(ParameterName =="/GeometryData/setVoxelsData"){
                LineString >> VoxXNumber >> VoxYNumber >> VoxZNumber >> word >> word >> VoxXHalfSize >> VoxYHalfSize >> VoxZHalfSize ;
                G4cout << " VoxXNumber " << VoxXNumber << " VoxYNumber " << VoxYNumber << " VoxZNumber " << VoxZNumber << " VoxXHalfSize " << VoxXHalfSize << " VoxYHalfSize " << VoxYHalfSize << " VoxZHalfSize " << VoxZHalfSize << G4endl ;
            }
            else if(ParameterName =="/RunAndScoreData/generateVoxelsResuls"){
                GenerateVoxelsResuls = true;
                G4cout << " GenerateVoxelsResuls " << GenerateVoxelsResuls << G4endl ;
            }
            
            else if(ParameterName =="/AnalysisData/generateSelfCrossGraphs"){// to evaluate the last 4 lines
                LineString >> GraphsData >> CompareType ;
                std::cout << " GraphsData " << GraphsData << " CompareType " << CompareType << std::endl ;
                std::string aa, bb;
                while (LineString >> aa && LineString >> bb) {
                    std::cout << " CompareReferenceName " << aa << " CompareReferencePath " << bb << std::endl ;
                    CompareReferenceNames.push_back(aa);
                    RefFilePaths.push_back(bb);
                }
                
                if(CompareReferenceNames.size() > 0 && RefFilePaths.size() > 0){
                    CompareReferenceName = CompareReferenceNames[0];
                    RefFilePath = RefFilePaths[0];
                }
            }
            
            else if(ParameterName =="/AnalysisData/generateRelativeErrGraph"){// to evaluate the last 4 lines
                LineString >> DifferenceMethod ;
                if(DifferenceMethod != ""){
                    GenerateRelativeErrGraph = "yes" ;
                    std::cout << " GenerateRelativeErrGraph " << GenerateRelativeErrGraph << " DifferenceMethod " << DifferenceMethod << std::endl ;
                }
            }
            else if(ParameterName =="/AnalysisData/generateRelativeSDevGraph"){// to evaluate the last 4 lines
                GenerateRelativeSDevGraph = "yes" ;
                std::cout << " GenerateRelativeSDevGraph " << GenerateRelativeSDevGraph << std::endl ;
            }
            else if(ParameterName =="/AnalysisData/generateVariableRegionGraph"){// to evaluate the last 4 lines
                LineString >> RegionVariableName ;
                GenerateRegionsVariableGraph = "yes" ;
                std::cout << " GenerateRegionsVariableGraph " << GenerateRegionsVariableGraph << " RegionVariableName " << RegionVariableName << std::endl ;
            }
            else if(ParameterName =="/PhysicsData/generateCrossSectionFor"){// to evaluate the last 4 lines
                GenerateCrossSectionGraph = "yes";
                std::cout << " GenerateCrossSectionGraph " << GenerateCrossSectionGraph << std::endl ;
            }
            else if(ParameterName =="/AnalysisData/generateEventsDataHisto"){// to evaluate the last 4 lines
                LineString >> PositionDataFile >> EnergyDataFile >> MomDirDataFile ;
                EventsDataHistograms = "yes";
                std::cout << " EventsDataHistograms " << EventsDataHistograms << " PositionDataFile " << PositionDataFile
                          << " EnergyDataFile " << EnergyDataFile << " MomDirDataFile " << MomDirDataFile << std::endl ;
            }
            else if(ParameterName =="/AnalysisData/generateVoxelizedHistograms"){// to evaluate the last 4 lines
                LineString >> DoseProfilQuantity >> BeamAxis >> SliceFor2DGraph >> SliceID ;
                std::cout << " DoseProfilQuantity " << DoseProfilQuantity << " BeamAxis " << BeamAxis << " SliceFor2DGraph " << SliceFor2DGraph << " SliceID " << SliceID << std::endl ;
            }
            else if(ParameterName =="/AnalysisData/setGraphsParameters"){// to evaluate the last 4 lines
                LineString >> UseLogE >> UseLogVariable >> UseGridXY >> PrintTitle >> LegendPos >> LegendXWidth >> LegendYHeight >> AddErrorBarInGraphs >> GraphsExt ;
                std::cout << " UseLogE " << UseLogE <<
                             " UseLogVariable " << UseLogVariable <<
                             " UseGridXY " << UseGridXY <<
                             " PrintTitle " << PrintTitle <<
                             " LegendPos " << LegendPos <<
                             " LegendXWidth " << LegendXWidth <<
                             " LegendYHeight " << LegendYHeight <<
                             " AddErrorBarInGraphs " << AddErrorBarInGraphs <<
                             " GraphsExt " << GraphsExt <<
                             std::endl ;
            }
        }
        
        file1.close();
    }
    else{
        
        std::cout << "cannot open the file " << MacrosStartingFile.c_str() << std::endl ;
    }
}

void G4DoseCalcsAnalysis::ReadResultsAndReferenceData(){
    ReadResultFile();
    if(RefFilePaths.size() > 0){
        ReferenceTable = ReadReferenceFile(RefFilePaths[0]);
    }
    if(RefFilePaths.size() > 1){
        ReferenceTable2 = ReadReferenceFile(RefFilePaths[1]);
    }

    std::sort(OrgansNameVector.begin(), OrgansNameVector.end());

}

// called from GenerateSelfCrossGraphs() and GenerateVoxel2DGraphs()
void G4DoseCalcsAnalysis::DataInitialization(){
    
    double X1 = 0.16;
    double Y1 = 0.16;
    double Xmax1 = 0.85;
    double Ymax1 = 0.85;
    
    if(LegendPos == "RightTop"){
        X1LegPos = Xmax1-LegendXWidth;
        Y1LegPos = Ymax1-LegendYHeight;
        X2LegPos = Xmax1;
        Y2LegPos = Ymax1;
    }
    else if(LegendPos == "LeftTop"){
        X1LegPos = X1;
        Y1LegPos = Ymax1-LegendYHeight;
        X2LegPos = X1+LegendXWidth;
        Y2LegPos = Ymax1;
    }
    else if(LegendPos == "RightBottom"){
        X1LegPos = Xmax1-LegendXWidth;
        Y1LegPos = Y1;
        X2LegPos = Xmax1;
        Y2LegPos = Y1+LegendYHeight;
    }
    else if(LegendPos == "LeftBottom"){
        X1LegPos = X1;
        Y1LegPos = Y1;
        X2LegPos = X1+LegendXWidth;
        Y2LegPos = Y1+LegendYHeight;
    }
    else if(LegendPos == "MiddleTop"){
        X1LegPos = 0.4;
        Y1LegPos = Ymax1-LegendYHeight;
        X2LegPos = 0.4+LegendXWidth;
        Y2LegPos = Ymax1;
    }
    else if(LegendPos == "MiddleBottom"){
        X1LegPos = 0.4;
        Y1LegPos = Y1;
        X2LegPos = 0.4+LegendXWidth;
        Y2LegPos = Y1+LegendYHeight;
    }
    
    
    if(RegionVariableName == "Volume"){
        RegionVariableNameWithUnit="Volume(cm3)";
    }
    else if(RegionVariableName == "Mass"){
        RegionVariableNameWithUnit="Mass(Kg)";
    }
    else if(RegionVariableName == "Density"){
        RegionVariableNameWithUnit="Density(mg/cm3)";
    }
    else if(RegionVariableName == "Distance"){
        RegionVariableNameWithUnit="Distance(mm)";
    }
    else if(RegionVariableName == "Target"){
        RegionVariableNameWithUnit="Targets";
    }
    
    
    QuantityUnit["AE"] = "AE("+AEUnit+")";
    QuantityUnit["AF"] = "AF";
    QuantityUnit["SAF"] = "SAF("+SAFUnit+")";
    QuantityUnit["AD"] = "AD("+ADUnit+")";
    QuantityUnit["S"] = "S("+SUnit+")";
    QuantityUnit["H"] = "H("+HUnit+")";
    QuantityUnit["E"] = "E("+EUnit+")";
    QuantityUnit["DR"] = "DR(\\%)";
    QuantityUnit["ER"] = "ER(\\%)";
    
    if(QuantitiesToScore == "AF" ){
        //std::cout << " QuantitiesToScore "<< QuantitiesToScore << " , and QuantityUnit[QuantitiesToScore] " << QuantityUnit[QuantitiesToScore] << std::endl;
        QuantityUseLog[QuantitiesToScore] = true;
    }
    else if(QuantitiesToScore == "AE"){
        //std::cout << " QuantitiesToScore "<< QuantitiesToScore << " , and QuantityUnit[QuantitiesToScore] " << QuantityUnit[QuantitiesToScore] << std::endl;
        QuantityUseLog[QuantitiesToScore] = false;
    }
    else if(QuantitiesToScore == "SAF"){
        //std::cout << " QuantitiesToScore "<< QuantitiesToScore << " , and QuantityUnit[QuantitiesToScore] " << QuantityUnit[QuantitiesToScore] << std::endl;
        QuantityUseLog[QuantitiesToScore] = true;
        
    }
    else if(QuantitiesToScore == "H" ){
        //std::cout << " QuantitiesToScore "<< QuantitiesToScore << " , and QuantityUnit[QuantitiesToScore] " << QuantityUnit[QuantitiesToScore] << std::endl;
        QuantityUseLog[QuantitiesToScore] = true;
        
    }
    else if(QuantitiesToScore == "E" ){
        //std::cout << " QuantitiesToScore "<< QuantitiesToScore << " , and QuantityUnit[QuantitiesToScore] " << QuantityUnit[QuantitiesToScore] << std::endl;
        QuantityUseLog[QuantitiesToScore] = true;
        
    }
    else if(QuantitiesToScore == "AD" ){
        //std::cout << " QuantitiesToScore "<< QuantitiesToScore << " , and QuantityUnit[QuantitiesToScore] " << QuantityUnit[QuantitiesToScore] << std::endl;
        QuantityUseLog[QuantitiesToScore] = false;
        
    }
    else if(QuantitiesToScore == "DR" ){
        //std::cout << " QuantitiesToScore "<< QuantitiesToScore << " , and QuantityUnit[QuantitiesToScore] " << QuantityUnit[QuantitiesToScore] << std::endl;
        QuantityUseLog[QuantitiesToScore] = false;
        
    }
    else if(QuantitiesToScore == "ER" ){
        //std::cout << " QuantitiesToScore "<< QuantitiesToScore << " , and QuantityUnit[QuantitiesToScore] " << QuantityUnit[QuantitiesToScore] << std::endl;
        QuantityUseLog[QuantitiesToScore] = false;
        
    }
    else if(QuantitiesToScore == "S" ){
        //std::cout << " QuantitiesToScore "<< QuantitiesToScore << " , and QuantityUnit[QuantitiesToScore] " << QuantityUnit[QuantitiesToScore] << std::endl;
        QuantityUseLog[QuantitiesToScore] = true;
    }
    
    if(QuantityUseLog[QuantitiesToScore] == true ){
        //MinValForLog =

        // iterations on particle name
        for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  )
        {
            Abeg->first;
            // iterations on source name
            for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  )
            {
                Bbeg->first;
                // iterations on target name
                for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg  )
                {
                    Cbeg->first;
                    for(int ii = 0; ii < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size(); ii++){
                        
                        double val = ResultTable[QuantitiesToScore][GeometrySymbol][Abeg->first][Bbeg->first][Cbeg->first][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]];
                        double val2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][Abeg->first][Bbeg->first][Cbeg->first][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]];
                        
                        if(val == 0.){
                            //ResultTable[QuantitiesToScore][GeometrySymbol][Abeg->first][Bbeg->first][Cbeg->first][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]] = MinValForLog;
                        }
                        if(val2 == 0.){
                            //ReferenceTable[QuantitiesToScore][GeometrySymbol][Abeg->first][Bbeg->first][Cbeg->first][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]] = MinValForLog;
                        }
                        
                        //std::cout <<  "Abeg->first " << Abeg->first  <<  "Bbeg->first " << Abeg->first  <<  "Cbeg->first " << Cbeg->first <<  "val " << val <<std::endl ;
                    }
                }
            }
        }
    }
    else{
        
    }
    
    if(DifferenceMethod == "LRD"){
        DiffSym = "LRD (\\%)";
        DiffExp = "Logarithmic Relative Difference (\\%)";
    }else if( DifferenceMethod == "RA"){
        DiffSym = "R";
        DiffExp = "Ratio ";
    }else{
        DiffSym = "RD (\\%)";
        DiffExp = "Relative Difference (\\%)";
    }
}

// called from GenerateResultGraphs() and GenerateResultReferenceGraphs
void G4DoseCalcsAnalysis::ReadResultFile(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    ResultTable.clear();
    GeometryList.clear();
    ParticleList.clear();
    RadiotracerList.clear();
    GeometryRadiotracerSources.clear();
    
    std::string line , unit, Geometry,Quantity,SrcRegionName , organTargetname, ParticleName , RadTracerName, word;
    double particleE , mass, volume, density, Val, StandardDeviation , err, CompuTime ;
    bool isADataLines = false ;
    unsigned long long int ival;
    std::string fileNameString = ResultDirectoryPath + "/ResultsData" ;
    std::ifstream fileR(fileNameString.c_str());
    
    if(fileR.is_open()){
        
        std::cout << "Reading Result file "<< fileNameString << " ... " << std::endl ;
        
        // to fill the map array DataTables[SrcRegionName][organTargetname][ParticleName][particleE] = Val;
        while (getline(fileR, line)) {
            
            //std::cout << " the line " << line << std::endl ;
            
            std::istringstream LineString(line);
            LineString >> word;
            
            //std::cout << " the word " << word << std::endl ;
            
            if(isADataLines == true){
                
                if (word == "*") {
                    
                    isADataLines = false;
                    continue;
                }
                if (word == "#") {
                    continue;
                }
                organTargetname = word;
                
                LineString >> Val;
                LineString >> StandardDeviation;
                LineString >> err;
                LineString >> ival; // for Num Steps
                LineString >> mass;
                LineString >> volume;
                LineString >> density;

                int Score = 0;

                // to see if we will save its data or not in ResultsTable
                if(IsAllSourcesToScore){
                    if(IsAllTargetsToScore){}
                    else{
                        bool IsIn = false;
                        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) { if(TargetNamesToScore[gg] == organTargetname){ IsIn = true; break; } }
                        if(IsIn == false){Score++;}
                    }
                }else{
                    bool IsIn = false;
                    for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) { if(SourceNamesToScore[gg] == SrcRegionName){ IsIn = true; break; } }
                    if(IsIn == false){Score++;}
                    if(IsAllTargetsToScore){}
                    else{
                        bool IsIn = false;
                        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) { if(TargetNamesToScore[gg] == organTargetname){ IsIn = true; break; } }
                        if(IsIn == false){Score++;}
                    }
                }

                //std::cout << "organTargetname " << organTargetname << " SrcRegionName " << SrcRegionName  << " Score " << Score << std::endl;
                
                if(ParticleName == "RadioTracer"){
                    
                    if(Score == 0){
                        ResultQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][RadTracerName][SrcRegionName][organTargetname] = Val;
                        QuantityGeometryRadioTracerSourceTargetStandartDeviation[Quantity][Geometry][RadTracerName][SrcRegionName][organTargetname] = StandardDeviation;
                        QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[Quantity][Geometry][RadTracerName][SrcRegionName][organTargetname] = err;
                        
                        bool IsIn = false;
                        for (int gg = 0 ; gg < OrgansNameVector.size() ; gg++) { if(OrgansNameVector[gg] == organTargetname){ IsIn = true; break; } }
                        if(IsIn == false){ OrgansNameVector.push_back(organTargetname);}
                        
                        IsIn = false;
                        for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) { if(SourceNamesToScore[gg] == SrcRegionName){ IsIn = true; break; } }
                        if(IsIn == false){ SourceNamesToScore.push_back(SrcRegionName);}

                        //std::cout << " Quantity " << Quantity << " Geometry " << Geometry << " RadTracerName " << RadTracerName << " SrcRegionName " << SrcRegionName << " organTargetname " << organTargetname << " Val " << Val << std::endl ;
                    }
                    //std::cout << " Quantity " << Quantity << " Geometry " << Geometry << " RadTracerName " << RadTracerName << " SrcRegionName " << SrcRegionName << " organTargetname " << organTargetname << " Val " << Val << std::endl ;
                }
                else{
                    
                    //std::string jjjj; LineString >> jjjj; // this is the error, but because the % we have to write the double value of error in the file to save it in the map
                    
                    if(Score == 0){
                        ResultTable[Quantity][Geometry][ParticleName][SrcRegionName][organTargetname][particleE] = Val;
                        StandartDeviation[Quantity][Geometry][ParticleName][SrcRegionName][organTargetname][particleE] = StandardDeviation;
                        RelativeStandartDeviationPerCent[Quantity][Geometry][ParticleName][SrcRegionName][organTargetname][particleE] = err;
                        
                        bool IsIn = false;
                        for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) { if(SourceNamesToScore[gg] == SrcRegionName){ IsIn = true; break; } }
                        if(IsIn == false){ SourceNamesToScore.push_back(SrcRegionName);}
                        
                        IsIn = false;
                        for (int ss = 0 ; ss < ResEnergies[Quantity][Geometry][ParticleName].size(); ss++) {if(ResEnergies[Quantity][Geometry][ParticleName][ss] == particleE){IsIn = true;}}
                        if(IsIn == false){ResEnergies[Quantity][Geometry][ParticleName].push_back(particleE);}
                        
                        IsIn = false;
                        for (int gg = 0 ; gg < OrgansNameVector.size() ; gg++) { if(OrgansNameVector[gg] == organTargetname){ IsIn = true; break; } }
                        if(IsIn == false){ OrgansNameVector.push_back(organTargetname);}
                        
                        //std::cout << " Quantity " << Quantity << " Geometry " << Geometry << " ParticleName " << ParticleName << " SrcRegionName " << SrcRegionName << " organTargetname " << organTargetname << " Val " << Val << " mass " << mass << " volume " << volume << " density " << density << std::endl ;
                    }
                    
                }
                
                if(Score == 0){
                    GeometryRegionVariableValue[Geometry]["Mass"][organTargetname] = mass;
                    GeometryRegionVariableValue[Geometry]["Volume"][organTargetname] = volume;
                    GeometryRegionVariableValue[Geometry]["Density"][organTargetname] = density;
                }
                if(MinValForLog > Val && Val != 0){
                    MinValForLog = Val;
                }
                //std::cout << "-----------------------MinValForLog: " << MinValForLog << " Val: " << Val << std::endl;

                //std::cout << "----------------------- For Part:" << ParticleName << " Src:"<<  SrcRegionName << " trg:" << organTargetname << " Val:" << Val << " Relative SD(%):" << err << std::endl;
            }
            
            if (word == "******" && isADataLines == false) {
                
                LineString >> Quantity >> SrcRegionName >> ParticleName ;
                //LineString >> Quantity >> unit >> SrcRegionName >> ParticleName ;
                //QuantityUnit[Quantity] = unit;
                
                if(ParticleName == "RadioTracer"){
                    LineString >> RadTracerName >> Geometry ;
                    
                    bool isin = false;
                    for (int ss = 0 ; ss < RadiotracerList.size(); ss++) {if(RadiotracerList[ss] == RadTracerName){isin = true;}}
                    if(isin == false){RadiotracerList.push_back(RadTracerName);}
                    
                    isin = false;
                    for (int ss = 0 ; ss < GeometryList.size(); ss++) {if(GeometryList[ss] == Geometry){isin = true;}}
                    if(isin == false){GeometryList.push_back(Geometry);}
                    
                    bool ScoreSource = true;
                    if(!IsAllSourcesToScore){
                        bool IsIn = false;
                        for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) { if(SourceNamesToScore[gg] == SrcRegionName){ IsIn = true; break; } }
                        if(IsIn == false){
                            ScoreSource = false;
                        }
                    }
                    
                    if(ScoreSource){
                        isin = false;
                        for (int ss = 0 ; ss < GeometryRadiotracerSources.size(); ss++) {if(GeometryRadiotracerSources[ss] == SrcRegionName){isin = true;}}
                        if(isin == false){GeometryRadiotracerSources.push_back(SrcRegionName);}
                    }
                    
                }
                else{
                    LineString >> particleE
                            >> Geometry >> word >> word >> word >> word
                            >> word >> word >> word >> word >> word
                            >> word >> word >> word >> word >> word
                            >> word >> CompuTime ;
                    
                    bool ScoreSource = true;
                    if(!IsAllSourcesToScore){
                        bool IsIn = false;
                        for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) { if(SourceNamesToScore[gg] == SrcRegionName){ IsIn = true; break; } }
                        if(IsIn == false){
                            ScoreSource = false;
                        }
                    }
                    
                    if(ScoreSource){
                        
                        bool isin = false;
                        for (int ss = 0 ; ss < GeometryList.size(); ss++) {if(GeometryList[ss] == Geometry){isin = true;}}
                        if(isin == false){GeometryList.push_back(Geometry);}
                        
                        isin = false;
                        for (int ss = 0 ; ss < ParticleList.size(); ss++) {if(ParticleList[ss] == ParticleName){isin = true;}}
                        if(isin == false){ParticleList.push_back(ParticleName);}
                        
                        //std::cout << " Quantity "<< Quantity << " SrcRegionName " << SrcRegionName << " ParticleName " << ParticleName << " particleE " << particleE  << " TotalEventNumber " << TotalEventNumber << std::endl ;
                        ResultParticleSourceEnergyTime[Geometry][ParticleName][SrcRegionName][particleE] = CompuTime;
                        //std::cout << "----------------------- For Part:" << ParticleName << " Src:"<<  SrcRegionName << " particleE:" << particleE << " CompuTime:" << CompuTime << std::endl;
                    }
                }
                
                //std::cout << " Geometry " << Geometry << " Quantity " << Quantity << " SrcRegionName " << SrcRegionName<< " ParticleName " << ParticleName  << std::endl ;
                
                isADataLines = true;
            }
        }
        
        std::cout << "----------------------- From DoseCalcs Results, Min Value For Log is: " << MinValForLog << " . Where we will update the 0 values for logarithmic scale " << std::endl;


        for ( auto Obeg = ResEnergies.begin(); Obeg != ResEnergies.end(); ++Obeg  ){
            
            for ( auto Abeg = Obeg->second.begin(); Abeg != Obeg->second.end(); ++Abeg  )
            {
                for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  )
                {
                    std::sort(ResEnergies[Obeg->first][Abeg->first][Bbeg->first].begin(), ResEnergies[Obeg->first][Abeg->first][Bbeg->first].end());
                }
            }
        }
        fileR.close();
    }
    else{
        
        std::cout << "cannot open the file " << fileNameString.c_str() << std::endl ;
    }
    
    for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol].end(); ++Abeg  )
    {
        //std::cout << "PARTICLE_NAME " << Abeg->first  <<std::endl ;
        
        // iterations on source name
        for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  )
        {
            //std::cout << "Source_ORG " << Bbeg->first  <<std::endl ;
            
            // iterations on target name
            for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg )
            {
                //std::cout << "Target_ORG " << Cbeg->first  <<std::endl ;
            }
        }
    }
    
    if(IsAllTargetsToScore){
        TargetNamesToScore.clear();
        for (G4int gg = 0 ; gg < OrgansNameVector.size() ; gg++) {
            TargetNamesToScore.push_back(OrgansNameVector[gg]);
        }
    }
    
    // because sources canot be read if the file by merge not readed and calculated its results
    if(IsAllSourcesToScore){
        //SourceNamesToScore.clear();
        for (G4int gg = 0 ; gg < OrgansNameVector.size() ; gg++) {
            //SourceNamesToScore.push_back(OrgansNameVector[gg]);
        }
    }
    
}
std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<double,double>>>>>> G4DoseCalcsAnalysis::ReadReferenceFile(std::string flnm){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    std::cout << flnm << "\n" << std::endl;
    
    std::string line , Geometry, Quantity, SrcRegionName , organTargetname, ParticleName , RadTracerName , word, QuantityUnit, EnergyUnit;
    double Val;
    bool isADataLines = false ;
    
    std::ifstream fileR(flnm.c_str());
    std::vector< double > particleEnergies;
    
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<double,double>>>>>> RefData;
    
    if(fileR.is_open()){
        
        std::cout << "Reading Reference file... " << std::endl ;
        
        // to fill the map array DataTables[SrcRegionName][organTargetname][ParticleName][particleE] = Val;
        while (getline(fileR, line)) {
            
            //std::cout << "the line " << line << std::endl ;
            
            std::istringstream LineString(line);
            LineString >> word;
            
            //std::cout << " the word " << word << std::endl ;
            
            if (word == "******" && isADataLines == false) {
                
                particleEnergies.clear();
                
                LineString >> Quantity >> QuantityUnit >> SrcRegionName >> ParticleName;
                //std::cout << " SrcRegionName " << SrcRegionName << " ParticleName " << ParticleName  << std::endl ;
                //std::cout << LineString.str().c_str() << std::endl ;
                
                if(ParticleName == "RadioTracer"){
                    LineString >> RadTracerName;
                    LineString >> Geometry;
                }else{
                    LineString >> Geometry >> EnergyUnit ;
                    
                    double Enee;
                    int numWords = 0;
                    while (LineString >> Enee) {
                        if(EnergyUnit == "keV" || EnergyUnit == "KeV" || EnergyUnit == "KEV" || EnergyUnit == "kev") {Enee = Enee * 0.001;}
                        else if(EnergyUnit == "eV" || EnergyUnit == "EV" || EnergyUnit == "Ev" || EnergyUnit == "ev") {Enee = Enee * 0.000001;}
                        particleEnergies.push_back(Enee);
                        //std::cout << " Enee " << Enee<< std::endl ;
                        ++numWords;
                    }
                    NumOfRefEne = particleEnergies.size();
                }
                
                isADataLines = true;
            }
            
            // read line of data that contains the name of target organ and the SAF correspondants
            if(isADataLines == true && word != "******"){
                
                if (word == "*") {
                    
                    isADataLines = false;
                    continue;
                }
                
                organTargetname = word;
                //std::cout << " - SrcRegionName " << SrcRegionName << " - organTargetname " << organTargetname << " - ParticleName " << ParticleName << std::endl ;
                if(ParticleName == "RadioTracer"){
                    LineString >> Val;
                    if(Quantity == "S"){
                        if(QuantityUnit == "rad/uCi-h" || QuantityUnit == "RAD/UCI-H" || QuantityUnit == "rad/uci-h") {Val = Val/*0.532*/*0.01/(37000000000*1e-6*3600);}
                    }
                    ReferenceQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][RadTracerName][SrcRegionName][organTargetname] = Val;
                    //std::cout << " - Quantity " << Quantity << " - Geometry " << Geometry << " - RadTracerName " << RadTracerName << " - SrcRegionName " << SrcRegionName << " - organTargetname " << organTargetname << " - Val " <<  Val  <<  std::endl ;
                }
                else{
                    
                    for( int zas = 0 ; zas < NumOfRefEne ; zas++ ){
                        LineString >> Val ;
                        if(Quantity == "SAF"){ if(QuantityUnit == "g-1" || QuantityUnit == "G-1" || QuantityUnit == "gram-1" || QuantityUnit == "Gram-1" || QuantityUnit == "GRAM-1") {Val = Val * 1000;}}
                        if(Quantity == "S"){if(QuantityUnit == "rad/uCi-h" || QuantityUnit == "RAD/UCI-H" || QuantityUnit == "rad/uci-h") {Val = Val/*0.532*/*0.01/(37000000000*1e-6*3600);}}
                        RefData[Quantity][Geometry][ParticleName][SrcRegionName][organTargetname][particleEnergies[zas]] = Val;

                        if(MinValForLog > Val && Val != 0){
                            MinValForLog = Val;
                        }

                        if(SrcRegionName== "Liver" && organTargetname== "Liver"){
                        //std::cout << " - SrcRegionName " << SrcRegionName << " - organTargetname " << organTargetname << " - ParticleName " << ParticleName << " particleEnergies[zas] " << particleEnergies[zas] << " - Val " <<  Val  <<  std::endl ;
                    }
                    }
                }
            }
        }
        std::cout << " From Reference, Min Value For Log is: " << MinValForLog << " . Where we will update the 0 values for logarithmic scale " << std::endl;

    }
    else {
        std::cout << "Canno't open the file " << flnm << std::endl ;
    }
    return RefData;
}

double* G4DoseCalcsAnalysis::AccumulateThirdRefDataToGraphs(std::string Part, std::string src, std::string trg){
    
    //ReferenceTable2 = ReadReferenceFile(RefFilePaths[1]);
    std::vector<double*> data;
    
    
    double ai[ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size()], ei[ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size()], errInERef[ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size()], errRef[ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size()] ;
    
    //std::cout << CompareReferenceNames[1] << " " <<  RefFilePaths[1] << " Data for " << Part << " " << src << " " << trg << " -Number of energies " << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() << std::endl ;
    
    //for ( auto Abeg = ReferenceTable2[QuantitiesToScore][GeometrySymbol][Part][src][trg].begin(); Abeg != ReferenceTable2[QuantitiesToScore][GeometrySymbol][Part][src][trg].end(); ++Abeg  )
    //{}
    
    for(int ii = 0; ii < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size(); ii++){
        
        ai[ii] = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii];
        ei[ii] = ReferenceTable2[QuantitiesToScore][GeometrySymbol][Part][src][trg][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]];

        //if(QuantityUseLog[QuantitiesToScore] == true && ei[ii] < MinValForLog){
        //  ei[ii] = MinValForLog;
        //}

        //std::cout <<  ai[ii] << " --- " << ei[ii] <<std::endl ;
    }
    
    return ei;
}
void G4DoseCalcsAnalysis::Remove0FromList(int num,double* nk , double* ss){
    
    std::vector<double> data;
    std::vector<double> data2;
    for(int ii = 0; ii < num; ii++){
        if(nk[ii] == 0 || ss[ii] == 0){
        }else{
            data.push_back(nk[ii]);
            data2.push_back(ss[ii]);
        }
    }
    
    double* newvec = new double[data.size()];
    double* newvec2 = new double[data.size()];
    for(int ii = 0; ii < data.size(); ii++){
        newvec[ii]  =data[ii];
        newvec2[ii] =data2[ii];
    }
    
    vecdata.val=data.size();
    vecdata.vec=newvec;
    vecdata.vec2=newvec2;
}
double G4DoseCalcsAnalysis::RelativeDifferenceCalculation(double val1, double val2){
    
    double a1 = val1;
    double a2 = val2;
    double a3;
    if( DifferenceMethod == "LRD"){
        if( a1 != 0. && a2 != 0. ){
            a3 = std::log(a1/a2)*100;
        }else{
            a3= NULL;
        }
    }
    if( DifferenceMethod == "RD"){
        if( a1 != 0. && a2 != 0. ){
            a3 = (a1-a2/a2)*100;
        }else{
            a3= NULL;
        }
    }
    if( DifferenceMethod == "RA"){
        if( a1 != 0. && a2 != 0. ){
            a3 = (a1/a2);
        }
        else{
            a3 = NULL;
        }
    }else{
        if( a1 != 0. && a2 != 0. ){
            a3 = (a1-a2/a2)*100;
        }else{
            a3 = NULL;
        }
    }
    
    return a3;
}
TGraph* G4DoseCalcsAnalysis::CreateGraph(int num, double* x, double* y){
    TGraph* graph = new TGraph (num , x, y );

    for(int ml = 0; ml < num; ml++){
        //std::cout << " y[ml] " << y[ml] << " MinValForLog " << MinValForLog <<  std::endl ;

        double val = 0. ;
        val = y[ml];
        if(val == 0. || __isnan(val) || __isinf(val) || val == NULL) //val <= MinValForLog ||
        {
            //std::cout << "------------------ " <<  std::endl ;
            graph->RemovePoint(ml);
        }
    }
    return graph;
}
TGraphErrors* G4DoseCalcsAnalysis::CreateGraphErrors(int num, double* x, double* y, double* a, double* b){
    TGraphErrors* graph = new TGraphErrors (num , x, y , a, b ); //this graph is related if we want to remove the zeros of SAF from Graph, if not we have to uncomment the graph1 below

    for(int ml = 0; ml < num; ml++){if(y[ml] < MinValForLog || y[ml] == MinValForLog ||y[ml] == 0. || __isnan(y[ml]) || __isinf(y[ml]) || y[ml] == NULL){graph->RemovePoint(ml);}}
    return graph;
}
TGraph* G4DoseCalcsAnalysis::setGraphData(TGraph* graph, int MS, int ColorNum){

    if(MS >= RootM.size()){MS=0;}
    if(ColorNum >= RootC.size()){ColorNum=0;}

    int MaxMarkerSize = RootM.size();
    int MarkerSize = MaxMarkerSize - MS;if(MarkerSize <= 0){ColorNum=MaxMarkerSize;}
    graph->SetMarkerSize(MarkerSize);

    graph->SetLineWidth(1.5);
    graph->SetMarkerSize(2);
    graph->SetMarkerStyle(RootM[MS]);
    graph->SetMarkerColor(RootC[ColorNum]);
    graph->SetLineColor(RootC[ColorNum]);

    return graph;
}
TGraphErrors* G4DoseCalcsAnalysis::setGraphErrorsData(TGraphErrors* graph, int MS, int ColorNum){

    if(MS >= RootM.size()){MS=0;}
    if(ColorNum >= RootC.size()){ColorNum=0;}

    int MaxMarkerSize = 4;
    int MarkerSize = MaxMarkerSize - MS;if(MarkerSize <= 0){ColorNum=MaxMarkerSize;}
    graph->SetMarkerSize(MarkerSize);

    graph->SetLineWidth(1.5);
    graph->SetMarkerSize(2);
    graph->SetMarkerStyle(RootM[MS]);
    graph->SetMarkerColor(RootC[ColorNum]);
    graph->SetLineColor(RootC[ColorNum]);

    return graph;
}
void G4DoseCalcsAnalysis::CreateMultiGraphParametersAndCanvas(std::string mgt, std::string fn, TMultiGraph* mg, TLegend*l){

    TCanvas* ResCanvas = new TCanvas(fn.c_str(), fn.c_str());

    if(UseLogVariable == "yes"){
        gPad->SetLogy(1);
    }
    if(UseLogE == "yes"){
        gPad->SetLogx(1);
    }
    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    mg->SetTitle(mgt.c_str());

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);

    l->SetX1(X1LegPos);
    l->SetX2(X2LegPos);
    l->SetY1(Y1LegPos);
    l->SetY2(Y2LegPos);

    l->Draw();

    ResCanvas->Print(fn.c_str());
    delete ResCanvas;

}


// Tables for self and cross and regions data
void G4DoseCalcsAnalysis::createLatexTables(){

    GenerateRegionsDataLatexTable();
    GenerateLatexTableResultReference();
    GenerateLatexTableResultReferenceForOneEnergy();
    GenerateLatexTableResultReferenceOfQuantitiesGeometriesRadioTracers();
    GenerateLatexTableResultForRadioTracerGeometry();
    GenerateLatexTablesForRadioTracerResult();
}
void G4DoseCalcsAnalysis::GenerateRegionsDataLatexTable(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    DataInitialization();
    
    std::string FileName = GraphsDirectoryPath + "RegionsLatexTables";
    
    std::ostringstream LatexText;
    
    LatexText << "=============================== Geometrical Data (volume, mass, density) in all simulated geometries \n\n";
    
    // ////////////////////////////////////////////////// Geometrical data, each geometry ina distinct table

    for ( int DD = 0; DD < GeometryList.size() ; DD++  ){
        
        GeometrySymbol = GeometryList[DD];
        
        LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                  << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n"
                  << "\\caption{Regions names implemented in the in " << GeometrySymbol << " phantom geometry and the corresponding volume, mass and density.} \n"
                  << "\\begin{tabular}{llll} \n";
        LatexText << "\\hline \n";
        LatexText << "\\textbf{Region name}         & \\textbf{Volume (cm3)}      & \\textbf{Mass (kg)}   & \\textbf{Density (g/cm3)}              ";
        
        for ( int A = 0; A < OrgansNameVector.size() ; A++  ){
            if(OrgansNameVector[A] == "World"){
                continue;
            }
            LatexText << "       \\\\\\hline\n";
            LatexText << OrgansNameVector[A] <<"     & " << std::fixed << std::setprecision(GeometryVarDigitNum) << GeometryRegionVariableValue[GeometrySymbol]["Volume"][OrgansNameVector[A]] <<"      & "<< GeometryRegionVariableValue[GeometrySymbol]["Mass"][OrgansNameVector[A]] <<"     & " << GeometryRegionVariableValue[GeometrySymbol]["Density"][OrgansNameVector[A]];
        }
        
        LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
        LatexText << "\\label{"<< GeometrySymbol <<"RegionImpleData}\n";
        LatexText << "\\end{table}\n\%\\end{sidewaystable}";
        LatexText << "\n\n\n\n\n";
        
    }
    

    // ////////////////////////////////////////////////// Geometrical data for all geometries in one table

    LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
              << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n"
              << "\\caption{Regions names implemented in the in geometries and the corresponding volume, mass and density.} \n"
              << "\\begin{tabular}{lll} \n";
    LatexText << "\\hline \n";
    LatexText << "\\textbf{Region name}         & \\textbf{Abbreviation}      & \\textbf{Mass (kg) (";
                 for ( int SD = 0; SD < GeometryList.size() ; SD++  ){
                     LatexText << GeometryList[SD];
                     if(SD != GeometryList.size()-1){ LatexText <<", ";}
                 }
                 LatexText << ")}";

    for ( int A = 0; A < OrgansNameVector.size() ; A++  ){
        if(OrgansNameVector[A] == "World"){
            continue;
        }
        LatexText << "       \\\\\\hline\n";
        LatexText << OrgansNameVector[A] <<"     & "
                    << OrgansNameVector[A] << "     & ";
                    for ( int SD = 0; SD < GeometryList.size() ; SD++  ){
                        LatexText << std::fixed << std::setprecision(GeometryVarDigitNum)  << GeometryRegionVariableValue[GeometryList[SD]]["Mass"][OrgansNameVector[A]];
                        if(SD != GeometryList.size()-1){ LatexText <<", ";}
                    }
                    //LatexText <<"     & ";
                    //for ( int SD = 0; SD < GeometryList.size() ; SD++  ){
                    //    LatexText << std::fixed << std::setprecision(GeometryVarDigitNum) << GeometryRegionVariableValue[GeometryList[SD]]["Density"][OrgansNameVector[A]];
                    //    if(SD != GeometryList.size()-1){ LatexText <<", ";}
                    //}
    }

    LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
    LatexText << "\\label{GeometriesRegionImplementedData}\n";
    LatexText << "\\end{table}\n\%\\end{sidewaystable}";
    LatexText << "\n\n\n\n\n";


    // //////////////////////////////////////////////////////////////////
    
    int NumberOfGeometries = GeometryList.size();
    int NumberOfCol = NumberOfGeometries + 2;
    
    int RowNum = 2;
    std::string ValNm = "Mass";
    
    std::cout << "\n\n                                                          ========= Creation of Latex Table For IntakeIntoBody, Geometry Data and radiotracers ========= "<< "\n" << std::endl;
    
    //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << GeometryRadiotracerSources[srcInc] << "\n\n";
    
    LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
    LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
    LatexText << "\\usepackage{multirow}\n";
    LatexText << "\n\n";
    
    LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
              << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";
    
    LatexText << "\\caption{DoseCalcs implemented mass (kg) and density (g/cm$^3$) of regions in the simulated phantom geometries} \n";
    
    LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
              << "\\begin{threeparttable}\n"
              << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";
    
    for ( int A = 0; A < NumberOfCol ; A++  )
    {
        LatexText << "l";
    }
    LatexText << "} \\hline \n";
    LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Region Name}}} & \\multirow{2}{*}{\\textbf{Parameter}} & \\multicolumn{";
    
    LatexText << NumberOfGeometries << "}{c}{ ";
    LatexText << " \\textbf{ Phantoms }}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
    
    for ( int A = 0; A < NumberOfGeometries ; A++  )
    {
        LatexText << "   & \\textbf{" << GeometryList[A]<< "}";
    }
    
    for ( int A = 0; A < OrgansNameVector.size() ; A++  ){
        
        LatexText << "       \\\\\\hline \n";
        LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{" << OrgansNameVector[A] << "}}       & " <<ValNm  << "                                         ";
        for ( int B = 0; B < NumberOfGeometries ; B++  ){
            LatexText << "& " << std::fixed << std::setprecision(GeometryVarDigitNum) << GeometryRegionVariableValue[GeometryList[B]]["Mass"][OrgansNameVector[A]] << "        ";
        }
        
        LatexText << "\\\\ \n                             & "<< "Density" << "                                        ";
        
        for ( int B = 0; B < NumberOfGeometries ; B++  ){
            double a1 = GeometryRegionVariableValue[GeometryList[B]]["Density"][OrgansNameVector[A]];
            LatexText << "& " << std::fixed << std::setprecision(GeometryVarDigitNum) << a1 << "        ";
        }
    }
    
    LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
    //LatexText << "\\begin{tablenotes}\\footnotesize\n";
    //LatexText << "\\item[a] " << DiffExp << "\n";
    //LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
    //LatexText << "\\end{tablenotes}\n";
    LatexText << "\\end{threeparttable}\n";
    LatexText << "\\end{adjustbox}\n";
    LatexText << "\\label{tab:DoseCalcsRegionGeometricalData}\n";
    LatexText << "\\end{table}\n\%\\end{sidewaystable}";
    LatexText << "\n\n";
    
    std::ofstream outfile(FileName , std::ios::app);
    if(outfile.is_open()){
        
        std::cout << "\nCreating file " << FileName << std::endl ;
        outfile << LatexText.str();
        outfile.close();
    }
    
    // //////////////////////////////////////////////////////////////////
}
void G4DoseCalcsAnalysis::GenerateLatexTableResultReference(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    ReadResultFile();
    for (int SelfOfCross = 1 ; SelfOfCross < 6 ; SelfOfCross++) { // 1: one table just for Self results, 2: a table for each source just for cross results, 3: one table for self-cross data for all sources and targets  results. 4: one table for self-cross data with ratio included with ref name for all sources and targets with rsd, 5: is like 4 without rsd

        if(RefFilePaths.size() == 0 || CompareReferenceNames.size() == 0 && ReferenceTable.size() == 0 && (SelfOfCross == 4 || SelfOfCross == 5)){
            //continue;
        }

        for (int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {
            for (int vvv = 0 ; vvv < GeometryList.size() ; vvv++) {
                for (int jjj = 0 ; jjj < ParticleList.size() ; jjj++) {
                    
                    QuantitiesToScore = QuantityNamesToScore[gg];
                    GeometrySymbol = GeometryList[vvv];
                    ParticleName = ParticleList[jjj];
                    
                    if(ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size()==0){continue;}
                    
                    DataInitialization();
                    
                    std::string FileName = GraphsDirectoryPath + "ResRefLatexTables";
                    
                    std::vector<std::string> SourcesName = SourceNamesToScore;
                    std::ostringstream LatexText;
                    
                    
                    int NumberOfEne = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size();
                    int NumberOfCol = NumberOfEne + 2;
                    
                    int RowNum = 2; // 2: for Data,SD, 4: Data,Ref1,RA,SD, 3: Data,Ref1,Ref2, 5: Data,Ref1+RA,ref2+RA,SD
                    //std::string ValNm = QuantitiesToScore;
                    std::string ValNm = "DoseCalcs";
                    

                    if(GraphsData == "Reference_Result" && RefFilePaths.size() != 0 && CompareReferenceNames.size() != 0 && ReferenceTable.size() != 0 ){
                        RowNum = 4;
                        if(RefFilePaths.size() > 1 && CompareReferenceNames.size() > 1 && CompareReferenceNames.size() == RefFilePaths.size()){
                            //RowNum = RefFilePaths.size();
                            RowNum = 3;
                        }
                        //ValNm = CompareReferenceName;
                    }
                    
                    if(SelfOfCross == 4 || SelfOfCross == 5){
                        RowNum = 5;
                    }
                    
                    if(SelfOfCross == 1 || SelfOfCross == 3 || SelfOfCross == 4 || SelfOfCross == 5){
                        
                        std::cout << "\n\n                                                          ========= Creation of Latex Table For irradiation Data ========= "<< "\n" << std::endl;
                        
                        //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << SourcesName[srcInc] << "\n\n";
                        
                        LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
                        LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
                        LatexText << "\\usepackage{multirow}\n";
                        LatexText << "\n\n";
                        
                        LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                                  << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";
                        if(SelfOfCross == 1){
                            if(RowNum == 4 || RowNum == 3){
                                LatexText << "\\caption{Self-absorption " << QuantitiesToScore << " values of " << ParticleName << " calculated in " << GeometrySymbol << " by DoseCalcs and compared to the " << CompareReferenceName << " reference ";
                                if(RefFilePaths.size() > 1 && CompareReferenceNames.size() > 1){
                                    LatexText << "and the " << CompareReferenceNames[1] << " reference";
                                }
                                LatexText << "} \n";
                            }else{
                                LatexText << "\\caption{Self-absorption " << QuantitiesToScore << " values of " << ParticleName << " calculated in " << GeometrySymbol << " by DoseCalcs} \n";
                            }
                        }
                        else if(SelfOfCross == 3 || SelfOfCross == 4 || SelfOfCross == 5){
                            if(RowNum == 4 || RowNum == 3 || RowNum == 5){
                                LatexText << "\\caption{Self- and cross-irradiation " << QuantitiesToScore << " values of " << ParticleName << " for each source-target combination calculated in " << GeometrySymbol << " by DoseCalcs and compared to the " << CompareReferenceName << " reference ";
                                if(RefFilePaths.size() > 1 && CompareReferenceNames.size() > 1){
                                    LatexText << "and the " << CompareReferenceNames[1] << " reference";
                                }
                                LatexText << "} \n";
                                
                            }else{
                                LatexText << "\\caption{Self- and cross-irradiation " << QuantitiesToScore << " values of " << ParticleName << " for each source-target combination calculated in " << GeometrySymbol << " by DoseCalcs} \n";
                            }
                        }
                        LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                                  << "\\begin{threeparttable}\n"
                                  << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";
                        
                        for ( int A = 0; A < NumberOfCol ; A++  )
                        {
                            LatexText << "l";
                        }
                        LatexText << "} \\hline \n";
                        if(SelfOfCross == 1){
                            LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                        }
                        else if(SelfOfCross == 3 || SelfOfCross == 4 || SelfOfCross == 5){
                            LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                        }
                        
                        LatexText << NumberOfEne << "}{c}{ ";
                        LatexText << " \\textbf{ Energies in MeV}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                        
                        for ( int A = 0; A < NumberOfEne ; A++  )
                        {
                            LatexText << "   & \\textbf{" << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][A]<< "}";
                        }
                    }
                    
                    for ( int srcInc = 0; srcInc < SourcesName.size() ; srcInc++  ){
                        
                        if(SelfOfCross == 2){
                            
                            std::cout << "\n\n                                                          ========= Creation of Latex Table For " << SourcesName[srcInc] << " Cross-irradiation Data ========= "<< "\n" << std::endl;
                            
                            //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << SourcesName[srcInc] << "\n\n";
                            
                            LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
                            LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
                            LatexText << "\\usepackage{multirow}\n";
                            LatexText << "\n\n";
                            
                            LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                                      << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";
                            if(RowNum == 4 || RowNum == 3){
                                LatexText << "\\caption{Cross-irradiation " << QuantitiesToScore << " values of " << ParticleName << " from " << SourcesName[srcInc] << " calculated in " << GeometrySymbol << " by DoseCalcs and compared to the " << CompareReferenceName << " reference ";
                                if(RefFilePaths.size() > 1 && CompareReferenceNames.size() > 1){
                                    LatexText << "and the " << CompareReferenceNames[1] << " reference";
                                }
                                LatexText << "} \n";
                            }else{
                                LatexText << "\\caption{Cross-irradiation " << QuantitiesToScore << " values of " << ParticleName << " from " << SourcesName[srcInc] << " calculated in " << GeometrySymbol << " by DoseCalcs} \n";
                            }
                            LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                                      << "\\begin{threeparttable}\n"
                                      << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";
                            
                            for ( int A = 0; A < NumberOfCol ; A++  )
                            {
                                LatexText << "l";
                            }
                            LatexText << "} \\hline \n";
                            LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Target region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                            
                            LatexText << NumberOfEne << "}{c}{ ";
                            LatexText << " \\textbf{" << SourcesName[srcInc] << " source with "<< ParticleName << " energies in MeV}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                            
                            for ( int A = 0; A < NumberOfEne ; A++  )
                            {
                                LatexText << "   & \\textbf{" << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][A]<< "}";
                            }
                            
                            for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                                
                                if(TargetNamesToScore[A] == SourcesName[srcInc]){
                                    //continue;
                                }
                                
                                LatexText << "       \\\\\\hline \n";
                                LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<< TargetNamesToScore[A] << "}}       & " <<ValNm  << "                                         ";
                                for ( int B = 0; B < NumberOfEne ; B++  ){
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                }
                                
                                if(RowNum == 4){
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceName << "                                        ";
                                    
                                    for ( int B = 0; B < NumberOfEne ; B++  ){
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                    }
                                    
                                    if(RefFilePaths.size() > 1 && CompareReferenceNames.size() > 1){
                                        LatexText << "\\\\ \n                             & "<< CompareReferenceNames[1] << "                                        ";
                                        
                                        for ( int B = 0; B < NumberOfEne ; B++ ){
                                            LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable2[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                        }
                                        
                                        LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{a}  " << "                                        ";
                                        
                                    }else{
                                        
                                        LatexText << "\\\\ \n                             & " << DiffSym << "\\tnote{a} " << "                                        ";
                                        
                                        for ( int B = 0; B < NumberOfEne ; B++  ){
                                            double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                            double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                            double a3 = RelativeDifferenceCalculation(a1,a2);
                                            
                                            LatexText << "& " << std::fixed << std::setprecision(DiffDigitNum) << a3 << "        ";
                                        }
                                        
                                        LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                                    }
                                }
                                else if(RowNum == 3){
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceName << "                                        ";
                                    
                                    for ( int B = 0; B < NumberOfEne ; B++  ){
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                    }
                                    
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceNames[1] << "                                        ";
                                    
                                    for ( int B = 0; B < NumberOfEne ; B++ ){
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable2[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                    }
                                }
                                else{
                                    LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{a}  " << "                                        ";
                                }
                                
                                if(RowNum != 3){
                                    for ( int B = 0; B < NumberOfEne ; B++  ){
                                        double a1 = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                                    }
                                }
                            }
                            
                            LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
                            LatexText << "\\begin{tablenotes}\\footnotesize\n";
                            
                            if(RowNum == 5){
                                LatexText << "\\item[a] " << DiffExp << "\n";
                                LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                            }
                            else if(RowNum == 4){
                                if(RefFilePaths.size() > 1 && CompareReferenceNames.size() > 1){
                                    LatexText << "\\item[a] Relative Standard Deviation (\\%)\n";
                                }else{
                                    LatexText << "\\item[a] " << DiffExp << "\n";
                                    LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                                }
                            }else if(RowNum == 3){}
                            else{
                                LatexText << "\\item[a] Relative Standard Deviation (\\%)\n";
                            }
                            LatexText << "\\end{tablenotes}\n";
                            LatexText << "\\end{threeparttable}\n";
                            LatexText << "\\end{adjustbox}\n";
                            LatexText << "\\label{tab:"<<GeometrySymbol<<"CrossFrom"<< SourcesName[srcInc] <<"}\n";
                            LatexText << "\\end{table}\n\%\\end{sidewaystable}";
                            LatexText << "\n\n";
                            
                        }
                        if(SelfOfCross == 1){
                            
                            LatexText << "       \\\\\\hline \n";
                            LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<< SourcesName[srcInc] << "}}       & " <<ValNm  << "                                         ";
                            for ( int B = 0; B < NumberOfEne ; B++  ){
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                            }
                            
                            if(RowNum == 4){
                                
                                LatexText << "\\\\ \n                             & "<< CompareReferenceName << "                                        ";
                                
                                for ( int B = 0; B < NumberOfEne ; B++  ){
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                }
                                
                                if(RefFilePaths.size() > 1 && CompareReferenceNames.size() > 1){
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceNames[1] << "                                        ";
                                    
                                    for ( int B = 0; B < NumberOfEne ; B++ ){
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable2[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                    }
                                    
                                    LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{a}  " << "                                        ";
                                    
                                }else{
                                    LatexText << "\\\\ \n                             & "<< DiffSym <<"\\tnote{a} " << "                                        ";
                                    
                                    for ( int B = 0; B < NumberOfEne ; B++  ){
                                        double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                        double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                        double a3 = RelativeDifferenceCalculation(a1,a2);
                                        
                                        LatexText << "& " << std::fixed << std::setprecision(DiffDigitNum) << a3 << "        ";
                                    }
                                    
                                    LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                                }
                            }
                            else if(RowNum == 3){
                                LatexText << "\\\\ \n                             & "<< CompareReferenceName << "                                        ";
                                
                                for ( int B = 0; B < NumberOfEne ; B++  ){
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                }
                                
                                LatexText << "\\\\ \n                             & "<< CompareReferenceNames[1] << "                                        ";
                                
                                for ( int B = 0; B < NumberOfEne ; B++ ){
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable2[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                }
                            }
                            else{
                                LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{a}  " << "                                        ";
                            }
                            
                            if(RowNum != 3){
                                for ( int B = 0; B < NumberOfEne ; B++  ){
                                    double a1 = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][SourcesName[srcInc]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                                }
                            }
                        }
                        if(SelfOfCross == 3){
                            
                            for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                                
                                if(TargetNamesToScore[A] == SourcesName[srcInc]){
                                    //continue;
                                }
                                
                                LatexText << "       \\\\\\hline \n";
                                LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<SourcesName[srcInc]<<"$\\to$" << TargetNamesToScore[A] << "}}       & " <<ValNm  << "                                         ";
                                for ( int B = 0; B < NumberOfEne ; B++  ){
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                }
                                
                                if(RowNum == 4){
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceName << "                                        ";
                                    
                                    for ( int B = 0; B < NumberOfEne ; B++  ){
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                    }
                                    
                                    if(RefFilePaths.size() > 1 && CompareReferenceNames.size() > 1){
                                        LatexText << "\\\\ \n                             & "<< CompareReferenceNames[1] << "                                        ";
                                        
                                        for ( int B = 0; B < NumberOfEne ; B++  ){
                                            LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum)<< ReferenceTable2[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                        }
                                        
                                        LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{a}  " << "                                        ";
                                        
                                    }else{
                                        
                                        LatexText << "\\\\ \n                             & " << DiffSym << "\\tnote{a} " << "                                        ";
                                        
                                        for ( int B = 0; B < NumberOfEne ; B++  ){
                                            double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                            double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                            double a3 = RelativeDifferenceCalculation(a1,a2);
                                            
                                            LatexText << "& " << std::fixed << std::setprecision(DiffDigitNum) << a3 << "        ";
                                        }
                                        
                                        LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                                    }
                                    
                                }
                                else if(RowNum == 3){
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceName << "                                        ";
                                    
                                    for ( int B = 0; B < NumberOfEne ; B++  ){
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                    }
                                    
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceNames[1] << "                                        ";
                                    
                                    for ( int B = 0; B < NumberOfEne ; B++ ){
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable2[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                    }
                                }
                                else{
                                    LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{a}  " << "                                        ";
                                }
                                
                                if(RowNum != 3){
                                    for ( int B = 0; B < NumberOfEne ; B++  ){
                                        double a1 = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                                    }
                                }
                            }
                        }
                        if(SelfOfCross == 4){
                            
                            for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                                
                                if(TargetNamesToScore[A] == SourcesName[srcInc]){
                                    //continue;
                                }
                                
                                int rrrr;
                                
                                if(RefFilePaths.size() == 1 && CompareReferenceNames.size() == 1){
                                    rrrr = RowNum-2;
                                }else{
                                    rrrr = RowNum-1;
                                }
                                LatexText << "       \\\\\\hline \n";
                                LatexText << "\\multirow{"<< rrrr <<"}{*}{\\textbf{"<<SourcesName[srcInc]<<"$\\to$" << TargetNamesToScore[A] << "}}       & " <<ValNm  << "                                         ";
                                for ( int B = 0; B < NumberOfEne ; B++  ){
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                }
                                
                                LatexText << "\\\\ \n                             & "<< CompareReferenceName << " (" << DiffSym << "\\tnote{a} )                                        ";
                                
                                for ( int B = 0; B < NumberOfEne ; B++  ){
                                    
                                    double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                    double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                    double a3 = RelativeDifferenceCalculation(a1,a2);
                                    
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                }
                                
                                if(RefFilePaths.size() > 1 && CompareReferenceNames.size() > 1){
                                    
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceNames[1] << " (" << DiffSym << "\\tnote{a} )                                        ";
                                    
                                    for ( int B = 0; B < NumberOfEne ; B++  ){
                                        double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                        double a2 = ReferenceTable2[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                        double a3 = RelativeDifferenceCalculation(a1,a2);
                                        
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable2[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                    }
                                }
                                
                                LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                                for ( int B = 0; B < NumberOfEne ; B++  ){
                                    double a1 = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                    LatexText << std::fixed << std::setprecision(DiffDigitNum);
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                                }
                            }
                        }
                        if(SelfOfCross == 5){
                            
                            for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                                
                                if(TargetNamesToScore[A] == SourcesName[srcInc]){
                                    //continue;
                                }
                                
                                int rrrr;
                                
                                if(RefFilePaths.size() == 1 && CompareReferenceNames.size() == 1){
                                    rrrr = 2;
                                }else{
                                    rrrr = 1;
                                }
                                LatexText << "       \\\\\\hline \n";
                                LatexText << "\\multirow{"<< rrrr <<"}{*}{\\textbf{"<<SourcesName[srcInc]<<"$\\to$" << TargetNamesToScore[A] << "}}       & " <<ValNm  << "                                         ";
                                for ( int B = 0; B < NumberOfEne ; B++  ){
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << "        ";
                                }
                                
                                LatexText << "\\\\ \n                             & "<< CompareReferenceName << " (" << DiffSym << "\\tnote{a} )                                        ";
                                
                                for ( int B = 0; B < NumberOfEne ; B++  ){
                                    
                                    double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                    double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                    double a3 = RelativeDifferenceCalculation(a1,a2);
                                    
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                }
                                
                                if(RefFilePaths.size() > 1 && CompareReferenceNames.size() > 1){
                                    
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceNames[1] << " (" << DiffSym << "\\tnote{a} )                                        ";
                                    
                                    for ( int B = 0; B < NumberOfEne ; B++  ){
                                        double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                        double a2 = ReferenceTable2[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]];
                                        double a3 = RelativeDifferenceCalculation(a1,a2);
                                        
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable2[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][B]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                    }
                                }
                            }
                        }
                    }
                    
                    if(SelfOfCross == 1 || SelfOfCross == 3|| SelfOfCross == 4|| SelfOfCross == 5){
                        LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
                        LatexText << "\\begin{tablenotes}\\footnotesize\n";
                    }
                    
                    if(SelfOfCross == 1 || SelfOfCross == 3){
                        if(RowNum == 4){
                            if(RefFilePaths.size() > 1 && CompareReferenceNames.size() > 1){
                                LatexText << "\\item[a] Relative Standard Deviation (\\%)\n";
                            }else{
                                LatexText << "\\item[a] " << DiffExp << "\n";
                                LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                            }
                        }
                        else if(RowNum == 3){}
                        else{
                            LatexText << "\\item[a] Relative Standard Deviation (\\%)\n";
                        }
                    }
                    else if(SelfOfCross == 4){
                        LatexText << "\\item[a] " << DiffExp << "\n";
                        LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                    }
                    else if(SelfOfCross == 5){
                        LatexText << "\\item[a] " << DiffExp << "\n";
                    }
                    
                    if(SelfOfCross == 1 || SelfOfCross == 3|| SelfOfCross == 4|| SelfOfCross == 5){
                        LatexText << "\\end{tablenotes}\n";
                        LatexText << "\\end{threeparttable}\n";
                        LatexText << "\\end{adjustbox}\n";
                        if(SelfOfCross == 1){LatexText << "\\label{tab:"<<GeometrySymbol<<"SelfIrr}\n";
                        }else if(SelfOfCross == 3){LatexText << "\\label{tab:"<<GeometrySymbol<<"SelfCrossIrr}\n";
                        }else if(SelfOfCross == 4){LatexText << "\\label{tab:"<<GeometrySymbol<<"SelfCrossIrrWithCompWithRSD}\n";
                        }else if(SelfOfCross == 5){LatexText << "\\label{tab:"<<GeometrySymbol<<"SelfCrossIrrWithCompWithoutRSD}\n";}
                        
                        LatexText << "\\end{table}\n\%\\end{sidewaystable}";
                        LatexText << "\n\n";
                    }
                    
                    std::ofstream outfile(FileName , std::ios::app);
                    if(outfile.is_open()){
                        
                        std::cout << "\nCreating file " << FileName << std::endl ;
                        outfile << LatexText.str();
                        outfile.close();
                    }
                }
            }
        }
    }
    
}
void G4DoseCalcsAnalysis::GenerateLatexTableResultReferenceForOneEnergy(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    ReadResultFile();
    for (int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {
        
        QuantitiesToScore = QuantityNamesToScore[gg];
        
        for (int ss = 0 ; ss < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; ss++) {
            
            DataInitialization();

            std::string FileName = GraphsDirectoryPath + "ResRefLatexTables";
            
            std::vector<std::string> SourcesName = SourceNamesToScore;
            std::ostringstream LatexText;
            
            int NumberOfTargets = TargetNamesToScore.size();
            int NumberOfCol = NumberOfTargets + 2;
            
            int RowNum = 2;
            //std::string ValNm = QuantitiesToScore;
            std::string ValNm = "DoseCalcs";
            
            if(GraphsData == "Reference_Result"){
                RowNum = 4;
                //ValNm = CompareReferenceName;
            }
            
            std::cout << "\n\n                                                          ========= Creation of Latex Table For " << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss] << " Self-absorption and Cross-irradiation Data ========= "<< "\n" << std::endl;
            
            //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << SourcesName[A] << "\n\n";
            
            LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                      << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";
            if(RowNum == 4){
                LatexText << "\\caption{" << QuantitiesToScore << " values for " << ParticleName << " of energy " << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss] << " MeV calculated in " << GeometrySymbol << " by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
            }
            else{
                LatexText << "\\caption{" << QuantitiesToScore << " values for " << ParticleName << " of energy " << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss] << " MeV calculated in " << GeometrySymbol << " by DoseCalcs} \n";
            }
            LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                      << "\\begin{threeparttable}\n"
                      << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";

            for ( int A = 0; A < NumberOfCol ; A++  )
            {
                LatexText << "l";
            }
            LatexText << "} \\hline \n";
            LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
            
            LatexText << NumberOfTargets << "}{c}{ ";
            LatexText << " \\textbf{Target region}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
            
            for ( int A = 0; A < NumberOfTargets ; A++  )
            {

                //std::cout << A << " " << NumberOfTargets << std::endl;

                LatexText << "   & \\textbf{" << TargetNamesToScore[A]<< "}";
            }


            for ( int A = 0; A < SourcesName.size() ; A++  ){
                
                //std::cout << "0" << std::endl;

                LatexText << "       \\\\\\hline \n";
                LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<< SourcesName[A] << "}}       & " <<ValNm  << "                                         ";
                for ( int B = 0; B < NumberOfTargets ; B++  ){
                    //G4cout << " " << QuantitiesToScore << " " << ParticleName << " " << SourcesName[A] << " " << TargetNamesToScore[B] << " " << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss] << " " << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[A]][TargetNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]] << " " << std::endl ;
                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[A]][TargetNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]] << "        ";
                }

                //std::cout << "1" << std::endl;


                if(RowNum == 4){
                    LatexText << "\\\\ \n                             & "<< CompareReferenceName << "                                        ";
                    
                    for ( int B = 0; B < NumberOfTargets ; B++  ){
                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[A]][TargetNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]] << "        ";
                    }
                    
                    LatexText << "\\\\ \n                             & "<< DiffSym <<"\\tnote{a} " << "                                        ";
                    
                    for ( int B = 0; B < NumberOfTargets ; B++  ){
                        double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[A]][TargetNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]];
                        double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[A]][TargetNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]];
                        double a3 = RelativeDifferenceCalculation(a1,a2);
                        
                        LatexText << "& " << std::fixed << std::setprecision(DiffDigitNum) << a3 << "        ";
                    }
                }
                
                //std::cout << "2" << std::endl;

                LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                
                for ( int B = 0; B < NumberOfTargets ; B++  ){
                    double a1 = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[A]][TargetNamesToScore[B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]];
                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                }
            }

            std::cout << "3" << std::endl;

            LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
            LatexText << "\\begin{tablenotes}\\footnotesize\n";
            LatexText << "\\item[a] " << DiffExp << " \n";
            LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
            LatexText << "\\end{tablenotes}\n";
            LatexText << "\\end{threeparttable}\n";
            LatexText << "\\end{adjustbox}\n";
            LatexText << "\\label{"<< QuantitiesToScore << "_" <<ParticleName <<"_SelfCross_"<< ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss] <<"}\n";
            LatexText << "\\end{table}\n\%\\end{sidewaystable}";
            LatexText << "\n\n";
            
            std::ofstream outfile(FileName , std::ios::app);
            if(outfile.is_open()){
                
                std::cout << "\nCreating file " << FileName << std::endl ;
                outfile << LatexText.str();
                outfile.close();
            }
        }
    }
    
    bool IsIn = false;
    for (int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {
        if(QuantityNamesToScore[gg] == "S"){IsIn = true;break;}
    }
    if(IsIn == true && ResultQuantityGeometryRadioTracerSourceTargetValues.size() != 0){
        GenerateLatexTableResultReferenceOfQuantitiesGeometriesRadioTracers();
    }
}
void G4DoseCalcsAnalysis::GenerateLatexTableResultReferenceOfQuantitiesGeometriesRadioTracers(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    ReadResultFile();
    
    for (int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {
        
        QuantitiesToScore = QuantityNamesToScore[gg];
        if(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].size()== 0
                && ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].size()== 0){continue;}
        
        for (int fd = 0 ; fd < GeometryList.size() ; fd++) {
            
            GeometrySymbol = GeometryList[fd];
            if(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]].size()== 0
                    && ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]].size()== 0){continue;}
            
            for (int qq = 0 ; qq < RadiotracerList.size() ; qq++) {
                
                RadioTracerName = RadiotracerList[qq];
                if(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName].size()== 0
                        && ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName].size()== 0){continue;}
                
                std::string FileName = GraphsDirectoryPath + "ResRefLatexTables";
                
                std::ostringstream LatexText;
                
                int NumberOfSourceTargetColumn = GeometryRadiotracerSources.size();
                int NumberOfCol = NumberOfSourceTargetColumn + 2;
                
                int RowNum = 2;
                
                std::string ValNm = "DoseCalcs";
                
                if(GraphsData == "Reference_Result"){
                    RowNum = 4;
                    //ValNm = CompareReferenceName;
                }
                
                std::cout << "\n\n                                                          ========= Creation of Latex Table of " << QuantityUnit[QuantitiesToScore] << " of " << RadioTracerName << " calculated in " << GeometrySymbol << " ========= "<< "\n" << std::endl;
                
                //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << SourcesName[A] << "\n\n";
                
                LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                          << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";
                if(RowNum = 4){
                    LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " for " << RadioTracerName << " calculated in " << GeometrySymbol << " by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                }
                else{
                    LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " for " << RadioTracerName << " calculated in " << GeometrySymbol << " by DoseCalcs} \n";
                }
                LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                          << "\\begin{threeparttable}\n"
                          << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";
                
                for ( int A = 0; A < NumberOfCol ; A++  )
                {
                    LatexText << "l";
                }
                LatexText << "} \\hline \n";
                LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Target region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                
                LatexText << NumberOfSourceTargetColumn << "}{c}{ ";
                LatexText << " \\textbf{Source region}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                
                for ( int B = 0; B < NumberOfSourceTargetColumn ; B++  )
                {
                    LatexText << "   & \\textbf{" << GeometryRadiotracerSources[B]<< "}";
                }
                
                for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                    
                    LatexText << "       \\\\\\hline \n";
                    LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<< TargetNamesToScore[A] << "}}       & " <<ValNm  << "                                         ";
                    
                    //for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg  ){}
                    
                    for ( int B = 0; B < NumberOfSourceTargetColumn ; B++  ){
                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][GeometryRadiotracerSources[B]][TargetNamesToScore[A]] << "        ";
                        //G4cout << " " << QuantitiesToScore << " " << ParticleName << " " << TargetNamesToScore[A] << " " << TargetNamesToScore[A] << " " << ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss] << " " << ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][SourcesName[srcInc]][TargetNamesToScore[A]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ss]] << " " << std::endl ;
                    }
                    
                    if(RowNum == 4){
                        LatexText << "\\\\ \n                             & "<< CompareReferenceName << "            ";
                        
                        for ( int B = 0; B < NumberOfSourceTargetColumn ; B++  ){
                            LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][GeometryRadiotracerSources[B]][TargetNamesToScore[A]] << "        ";
                        }
                        
                        LatexText << "\\\\ \n                             & "<< DiffSym <<"\\tnote{a} " << "                                        ";
                        
                        for ( int B = 0; B < NumberOfSourceTargetColumn ; B++  ){

                            double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][GeometryRadiotracerSources[B]][TargetNamesToScore[A]];
                            double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][GeometryRadiotracerSources[B]][TargetNamesToScore[A]];
                            double a3 = RelativeDifferenceCalculation(a1,a2);
                            
                            LatexText << "& " << std::fixed << std::setprecision(DiffDigitNum) << a3 << "        ";
                        }
                    }
                    
                    LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                    
                    for ( int B = 0; B < NumberOfSourceTargetColumn ; B++  ){
                        double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometrySymbol][RadioTracerName][GeometryRadiotracerSources[B]][TargetNamesToScore[A]];
                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                    }
                    
                }
                
                LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
                LatexText << "\\begin{tablenotes}\\footnotesize\n";
                LatexText << "\\item[a] " << DiffExp << " \n";
                LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                LatexText << "\\end{tablenotes}\n";
                LatexText << "\\end{threeparttable}\n";
                LatexText << "\\end{adjustbox}\n";
                LatexText << "\\label{"<< QuantityUnit[QuantitiesToScore] << "_" << GeometrySymbol <<"_" << RadioTracerName<<"}\n";
                LatexText << "\\end{table}\n\%\\end{sidewaystable}";
                LatexText << "\n\n";
                
                std::ofstream outfile(FileName , std::ios::app);
                if(outfile.is_open()){
                    
                    std::cout << "\nCreating file " << FileName << std::endl ;
                    outfile << LatexText.str();
                    outfile.close();
                }
                
                for ( int srcInc = 0; srcInc < GeometryRadiotracerSources.size() ; srcInc++  ){
                    
                }
            }
        }
    }
    
}
void G4DoseCalcsAnalysis::GenerateLatexTableResultForRadioTracerGeometry(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    if(ResultQuantityGeometryRadioTracerSourceTargetValues.size()== 0){return;}
    
    for (int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {
        
        QuantitiesToScore = QuantityNamesToScore[gg];
        if(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].size()== 0){continue;}
        
        int RowNum = 2;
        std::string ValNm = "DoseCalcs";
        
        if(GraphsData == "Reference_Result" && ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].size() != 0){
            RowNum = 3;
            //ValNm = CompareReferenceName;
        }
        
        for (int sss = 0 ; sss < 3 ; sss++) { //0 self(sources) and 1 Cross(source-target combinations) // 2 for SourceTargets
            for (int fd = 0 ; fd < GeometryList.size() ; fd++) {
                
                if(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]].size()== 0){continue;}
                
                int NumberOfRadioTracerOrGeom = RadiotracerList.size();
                int NumberOfCol = NumberOfRadioTracerOrGeom + 2;
                
                std::cout << "\n\n                                                          ========= Creation of Latex Table For Radiotracer Data in geometries ========= "<< "\n" << std::endl;
                
                DataInitialization();
                std::string FileName = GraphsDirectoryPath + "ResRefLatexTables";
                
                if(sss==2){
                    for ( int srcInc = 0; srcInc < GeometryRadiotracerSources.size() ; srcInc++  ){
                        
                        //if(GeometryRadiotracerSources[srcInc] == "IntakeIntoBody"){ continue;}
                        
                        std::ostringstream LatexText;
                        
                        std::cout << "\n\n                                                          ========= Creation of Latex Table For Geometry Data and radiotracers ========= "<< "\n" << std::endl;
                        
                        //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << GeometryRadiotracerSources[srcInc] << "\n\n";
                        
                        LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
                        LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
                        LatexText << "\\usepackage{multirow}\n";
                        LatexText << "\n\n";
                        
                        LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                                  << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";
                        
                        LatexText << "\%for long tatbles ------------------------------ \n\%\\begin{tiny}\n\% \\begin{longtable}{@{\\extracolsep{\\fill}}*{"<<NumberOfCol<<"}{l}}\n\%------------------------------\n";
                        
                        if(RowNum == 3){
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << GeometryList[fd] << " for internal irradiation from "<<GeometryRadiotracerSources[srcInc]<<" calculated in phantoms target regions by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                        }else{
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << GeometryList[fd] << " for internal irradiation from  "<<GeometryRadiotracerSources[srcInc]<<" calculated in phantoms target regions by DoseCalcs} \n";
                        }
                        
                        LatexText << "\%for long tatbles ------------------------------ \n\%\\\\ \n\%\\hline\n\% \\endfirsthead\n\% \\caption[]{(\\textit{continued})}\\\\\n\% \\hline\n\%------------------------------\n";
                        
                        LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                                  << "\\begin{threeparttable}\n"
                                  << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";
                        
                        for ( int A = 0; A < NumberOfCol ; A++  )
                        {
                            LatexText << "l";
                        }
                        LatexText << "} \\hline \n";
                        
                        LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                        LatexText << NumberOfRadioTracerOrGeom << "}{c}{ ";
                        LatexText << " \\textbf{Phantoms}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                        for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                        {
                            LatexText << "   & \\textbf{" << RadiotracerList[A]<< "}";
                        }
                        
                        LatexText << "\%for long tatbles ------------------------------ \n\%\\\\\\hline \n\%\\endhead \n\%\\hline \n\%\\multicolumn{6}{r}{(\\textit{Continued on next page})} \\\\ \n\%\\endfoot \n\%\\hline \n\%\\endlastfoot\n";
                        
                        LatexText << "\n\%\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                        LatexText << NumberOfRadioTracerOrGeom << "}{c}{ ";
                        LatexText << "\n\% \\textbf{Phantoms}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n\%                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                        for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                        {
                            LatexText << "   & \\textbf{" << RadiotracerList[A]<< "}";
                        }
                        LatexText << "\n\%------------------------------\n";
                        
                        for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                            
                            LatexText << "       \\\\\\hline \n";
                            LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<TargetNamesToScore[A]<<"$\\leftarrow$" <<GeometryRadiotracerSources[srcInc] << "}}       & " <<ValNm  << "                                         ";
                            
                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << "        ";
                            }
                            
                            if(RowNum == 3){
                                LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";
                                
                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    double a3 = RelativeDifferenceCalculation(a1,a2);
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                }
                            }
                            
                            LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                            
                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                            }
                        }
                        
                        LatexText << "\%for long tables ------------------------------ \n\%\\end{longtable} \n\%\\footnotetext[1]{Ratio} \n\%\\footnotetext[2]{Relative standard deviation}\n\%\\end{tiny}\n\%------------------------------\n";
                        
                        LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
                        LatexText << "\\begin{tablenotes}\\footnotesize\n";
                        LatexText << "\\item[a] " << DiffExp << "\n";
                        LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                        LatexText << "\\end{tablenotes}\n";
                        LatexText << "\\end{threeparttable}\n";
                        LatexText << "\\end{adjustbox}\n";
                        LatexText << "\\label{tab:"<<GeometryRadiotracerSources[srcInc]<<"RadiotracersIn" << GeometryList[fd] << "}\n";
                        LatexText << "\\end{table}\n\%\\end{sidewaystable}";
                        LatexText << "\n\n";
                        
                        std::ofstream outfile(FileName , std::ios::app);
                        if(outfile.is_open()){
                            
                            std::cout << "\nCreating file " << FileName << std::endl ;
                            outfile << LatexText.str();
                            outfile.close();
                        }
                    }
                }
                else{
                    std::ostringstream LatexText;
                    
                    //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << GeometryRadiotracerSources[srcInc] << "\n\n";
                    
                    LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
                    LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
                    LatexText << "\\usepackage{multirow}\n";
                    LatexText << "\n\n";
                    
                    LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                              << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";
                    
                    if(sss==0){
                        if(RowNum == 3){
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of radiotracers for each source region calculated in " << GeometryList[fd] << " by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                        }else{
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of radiotracers for each source region calculated in " << GeometryList[fd] << " by DoseCalcs} \n";
                        }
                    }
                    else if(sss==1){
                        if(RowNum == 3){
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of radiotracers for each source-target combination calculated in " << GeometryList[fd] << " by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                        }else{
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of radiotracers for each source-target combination calculated in " << GeometryList[fd] << " by DoseCalcs} \n";
                        }
                    }
                    
                    
                    LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                              << "\\begin{threeparttable}\n"
                              << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";
                    
                    for ( int A = 0; A < NumberOfCol ; A++  )
                    {
                        LatexText << "l";
                    }
                    LatexText << "} \\hline \n";
                    
                    if(sss==0){
                        LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source Region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                    }
                    else if(sss==1){
                        LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                    }
                    
                    LatexText << NumberOfRadioTracerOrGeom << "}{c}{ ";
                    LatexText << " \\textbf{Radiotracers}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                    
                    for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                    {
                        LatexText << "   & \\textbf{" << RadiotracerList[A]<< "}";
                    }
                    
                    for ( int srcInc = 0; srcInc < GeometryRadiotracerSources.size() ; srcInc++  ){
                        
                        //if(GeometryRadiotracerSources[srcInc] == "IntakeIntoBody"){ continue;}
                        if(sss==0){
                            LatexText << "       \\\\\\hline \n";
                            LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<GeometryRadiotracerSources[srcInc] << "}}       & " <<ValNm  << "                                         ";
                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]] << "        ";
                            }
                            
                            if(RowNum == 3){
                                LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";
                                
                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]];
                                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]];
                                    double a3 = RelativeDifferenceCalculation(a1,a2);
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                }
                            }
                            
                            LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                            
                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]];
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                            }
                        }
                        else if(sss==1){
                            for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                                
                                LatexText << "       \\\\\\hline \n";
                                LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<TargetNamesToScore[A]<<"$\\leftarrow$" <<GeometryRadiotracerSources[srcInc] << "}}       & " <<ValNm  << "                                         ";
                                
                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << "        ";
                                }
                                
                                if(RowNum == 3){
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";
                                    
                                    for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                        double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                        double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                        double a3 = RelativeDifferenceCalculation(a1,a2);
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                    }
                                }
                                
                                LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                                
                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                                }
                            }
                        }
                    }
                    
                    LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
                    LatexText << "\\begin{tablenotes}\\footnotesize\n";
                    LatexText << "\\item[a] " << DiffExp << "\n";
                    LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                    LatexText << "\\end{tablenotes}\n";
                    LatexText << "\\end{threeparttable}\n";
                    LatexText << "\\end{adjustbox}\n";
                    if(sss==0){
                        LatexText << "\\label{tab:SourceRadiotracersIn" << GeometryList[fd] << "}\n";
                    }
                    else if(sss==1){
                        LatexText << "\\label{tab:SourceTargetRadiotracersIn" << GeometryList[fd] << "}\n";
                    }
                    LatexText << "\\end{table}\n\%\\end{sidewaystable}";
                    LatexText << "\n\n";
                    
                    std::ofstream outfile(FileName , std::ios::app);
                    if(outfile.is_open()){
                        
                        std::cout << "\nCreating file " << FileName << std::endl ;
                        outfile << LatexText.str();
                        outfile.close();
                    }
                }
            }
            for (int fd = 0 ; fd < RadiotracerList.size() ; fd++) {
                
                int NumberOfRadioTracerOrGeom = GeometryList.size();
                int NumberOfCol = NumberOfRadioTracerOrGeom + 2;
                
                std::cout << "\n\n                                                          ========= Creation of Latex Table For Radiotracer Data for radiotracers ========= "<< "\n" << std::endl;
                
                DataInitialization();
                std::string FileName = GraphsDirectoryPath + "ResRefLatexTables";
                
                if(sss==2){
                    for ( int srcInc = 0; srcInc < GeometryRadiotracerSources.size() ; srcInc++  ){
                        
                        //if(GeometryRadiotracerSources[srcInc] == "IntakeIntoBody"){ continue;}
                        
                        std::ostringstream LatexText;
                        
                        std::cout << "\n\n                                                          ========= Creation of Latex Table For Geometry Data and radiotracers ========= "<< "\n" << std::endl;
                        
                        //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << GeometryRadiotracerSources[srcInc] << "\n\n";
                        
                        LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
                        LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
                        LatexText << "\\usepackage{multirow}\n";
                        LatexText << "\n\n";
                        
                        LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                                  << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";
                        
                        LatexText << "\%for long tatbles ------------------------------ \n\%\\begin{tiny}\n\% \\begin{longtable}{@{\\extracolsep{\\fill}}*{"<<NumberOfCol<<"}{l}}\n\%------------------------------\n";
                        
                        if(RowNum == 3){
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for internal irradiation from "<<GeometryRadiotracerSources[srcInc]<<" calculated in phantoms target regions by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                        }else{
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for internal irradiation from "<<GeometryRadiotracerSources[srcInc]<<" calculated in phantoms target regions  by DoseCalcs} \n";
                        }
                        
                        LatexText << "\%for long tatbles ------------------------------ \n\%\\\\ \n\%\\hline\n\% \\endfirsthead\n\% \\caption[]{(\\textit{continued})}\\\\\n\% \\hline\n\%------------------------------\n";
                        
                        LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                                  << "\\begin{threeparttable}\n"
                                  << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";
                        
                        for ( int A = 0; A < NumberOfCol ; A++  )
                        {
                            LatexText << "l";
                        }
                        LatexText << "} \\hline \n";
                        
                        LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                        LatexText << NumberOfRadioTracerOrGeom << "}{c}{ ";
                        LatexText << " \\textbf{Phantoms}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                        for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                        {
                            LatexText << "   & \\textbf{" << GeometryList[A]<< "}";
                        }
                        
                        LatexText << "\%for long tatbles ------------------------------ \n\%\\\\\\hline \n\%\\endhead \n\%\\hline \n\%\\multicolumn{6}{r}{(\\textit{Continued on next page})} \\\\ \n\%\\endfoot \n\%\\hline \n\%\\endlastfoot\n";
                        
                        LatexText << "\n\%\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                        LatexText << NumberOfRadioTracerOrGeom << "}{c}{ ";
                        LatexText << "\n\% \\textbf{Phantoms}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n\%                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                        for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                        {
                            LatexText << "   & \\textbf{" << GeometryList[A]<< "}";
                        }
                        LatexText << "\n\%------------------------------\n";
                        
                        for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                            
                            LatexText << "       \\\\\\hline \n";
                            LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<TargetNamesToScore[A]<<"$\\leftarrow$" <<GeometryRadiotracerSources[srcInc] << "}}       & " <<ValNm  << "                                         ";
                            
                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << "        ";
                            }
                            
                            if(RowNum == 3){
                                LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";
                                
                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    double a3 = RelativeDifferenceCalculation(a1,a2);
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                }
                            }
                            
                            LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                            
                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                            }
                        }
                        
                        LatexText << "\%for long tatbles ------------------------------ \n\%\\end{longtable} \n\%\\footnotetext[1]{Ratio} \n\%\\footnotetext[2]{Relative standard deviation}\n\%\\end{tiny}\n\%------------------------------\n";
                        
                        LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
                        LatexText << "\\begin{tablenotes}\\footnotesize\n";
                        LatexText << "\\item[a] " << DiffExp << "\n";
                        LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                        LatexText << "\\end{tablenotes}\n";
                        LatexText << "\\end{threeparttable}\n";
                        LatexText << "\\end{adjustbox}\n";
                        LatexText << "\\label{tab:"<<GeometryRadiotracerSources[srcInc]<<"GeometriesFrom" << RadiotracerList[fd] << "}\n";
                        LatexText << "\\end{table}\n\%\\end{sidewaystable}";
                        LatexText << "\n\n";
                        
                        std::ofstream outfile(FileName , std::ios::app);
                        if(outfile.is_open()){
                            
                            std::cout << "\nCreating file " << FileName << std::endl ;
                            outfile << LatexText.str();
                            outfile.close();
                        }
                    }
                }
                else{
                    
                    std::ostringstream LatexText;
                    
                    std::cout << "\n\n                                                          ========= Creation of Latex Table For Geometry Data and radiotracers ========= "<< "\n" << std::endl;
                    
                    //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << GeometryRadiotracerSources[srcInc] << "\n\n";
                    
                    LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
                    LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
                    LatexText << "\\usepackage{multirow}\n";
                    LatexText << "\n\n";
                    
                    LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                              << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";
                    
                    LatexText << "\%for long tatbles ------------------------------ \n\%\\begin{tiny}\n\% \\begin{longtable}{@{\\extracolsep{\\fill}}*{"<<NumberOfCol<<"}{l}}\n\%------------------------------\n";
                    
                    if(sss==0){
                        if(RowNum == 3){
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for each source region calculated in simulated geometries by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                        }else{
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for each source region calculated in simulated geometries  by DoseCalcs} \n";
                        }
                    }
                    else if(sss==1){
                        if(RowNum == 3){
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for each source-target combination calculated in simulated geometries by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                        }else{
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for each source-target combination calculated in simulated geometries  by DoseCalcs} \n";
                        }
                    }
                    
                    LatexText << "\%for long tatbles ------------------------------ \n\%\\\\ \n\%\\hline\n\% \\endfirsthead\n\% \\caption[]{(\\textit{continued})}\\\\\n\% \\hline\n\%------------------------------\n";
                    
                    LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                              << "\\begin{threeparttable}\n"
                              << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";
                    
                    for ( int A = 0; A < NumberOfCol ; A++  )
                    {
                        LatexText << "l";
                    }
                    LatexText << "} \\hline \n";
                    
                    if(sss==0){
                        LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source Region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                    }
                    else if(sss==1){
                        LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                    }
                    LatexText << NumberOfRadioTracerOrGeom << "}{c}{ ";
                    LatexText << " \\textbf{Phantoms}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                    for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                    {
                        LatexText << "   & \\textbf{" << GeometryList[A]<< "}";
                    }
                    
                    LatexText << "\%for long tatbles ------------------------------ \n\%\\\\\\hline \n\%\\endhead \n\%\\hline \n\%\\multicolumn{6}{r}{(\\textit{Continued on next page})} \\\\ \n\%\\endfoot \n\%\\hline \n\%\\endlastfoot\n";
                    if(sss==0){
                        LatexText << "\n\%\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source Region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                    }
                    else if(sss==1){
                        LatexText << "\n\%\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                    }
                    LatexText << NumberOfRadioTracerOrGeom << "}{c}{ ";
                    LatexText << "\n\% \\textbf{Phantoms}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n\%                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                    for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                    {
                        LatexText << "   & \\textbf{" << GeometryList[A]<< "}";
                    }
                    LatexText << "\n\%------------------------------\n";
                    
                    for ( int srcInc = 0; srcInc < GeometryRadiotracerSources.size() ; srcInc++  ){
                        //if(GeometryRadiotracerSources[srcInc] == "IntakeIntoBody"){ continue;}
                        
                        if(sss==0){
                            LatexText << "       \\\\\\hline \n";
                            LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<GeometryRadiotracerSources[srcInc] << "}}       & " <<ValNm  << "                                         ";
                            
                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]] << "        ";
                            }
                            
                            if(RowNum == 3){
                                LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";
                                
                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]];
                                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]];
                                    double a3 = RelativeDifferenceCalculation(a1,a2);
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                }
                            }
                            
                            LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                            
                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]];
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                            }
                        }
                        else if(sss==2){
                            
                        }
                        else if(sss==1){
                            for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                                
                                LatexText << "       \\\\\\hline \n";
                                LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<TargetNamesToScore[A]<<"$\\leftarrow$" <<GeometryRadiotracerSources[srcInc] << "}}       & " <<ValNm  << "                                         ";
                                
                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << "        ";
                                }
                                
                                if(RowNum == 3){
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";
                                    
                                    for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                        double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                        double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                        double a3 = RelativeDifferenceCalculation(a1,a2);
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                    }
                                }
                                
                                LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                                
                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                                }
                            }
                        }
                    }
                    
                    LatexText << "\%for long tatbles ------------------------------ \n\%\\end{longtable} \n\%\\footnotetext[1]{Ratio} \n\%\\footnotetext[2]{Relative standard deviation}\n\%\\end{tiny}\n\%------------------------------\n";
                    
                    LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
                    LatexText << "\\begin{tablenotes}\\footnotesize\n";
                    LatexText << "\\item[a] " << DiffExp << "\n";
                    LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                    LatexText << "\\end{tablenotes}\n";
                    LatexText << "\\end{threeparttable}\n";
                    LatexText << "\\end{adjustbox}\n";
                    if(sss==0){
                        LatexText << "\\label{tab:SourceGeometriesFrom" << RadiotracerList[fd] << "}\n";
                    }
                    else if(sss==1){
                        LatexText << "\\label{tab:SourceTargetGeometriesFrom" << RadiotracerList[fd] << "}\n";
                    }
                    LatexText << "\\end{table}\n\%\\end{sidewaystable}";
                    LatexText << "\n\n";
                    
                    std::ofstream outfile(FileName , std::ios::app);
                    if(outfile.is_open()){
                        
                        std::cout << "\nCreating file " << FileName << std::endl ;
                        outfile << LatexText.str();
                        outfile.close();
                    }
                }
            }
        }
        
        // just for IntakeIntoDoby source
        
        for (int fd = 0 ; fd < GeometryList.size() ; fd++) {
            
            if(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]].size()== 0){continue;}
            
            DataInitialization();
            std::string FileName = GraphsDirectoryPath + "ResRefLatexTables";
            std::ostringstream LatexText;
            
            int NumberOfRadioTracerOrGeom = RadiotracerList.size();
            int NumberOfCol = NumberOfRadioTracerOrGeom + 2;
            
            int RowNum = 2;
            //std::string ValNm = QuantitiesToScore;
            std::string ValNm = "DoseCalcs";
            
            if(GraphsData == "Reference_Result" && ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].size() != 0){
                RowNum = 4;
                //ValNm = CompareReferenceName;
            }
            
            std::cout << "\n\n                                                          ========= Creation of Latex Table For Intake Data ========= "<< "\n" << std::endl;
            
            LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
            LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
            LatexText << "\\usepackage{multirow}\n";
            LatexText << "\n\n";
            
            LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                      << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";
            
            if(RowNum == 4){
                LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of radiotracers for IntakeIntoBody calculated in " << GeometryList[fd] << " by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
            }else{
                LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of radiotracers for IntakeIntoBody calculated in " << GeometryList[fd] << " by DoseCalcs} \n";
            }
            
            LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                      << "\\begin{threeparttable}\n"
                      << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";
            
            for ( int A = 0; A < NumberOfCol ; A++  )
            {
                LatexText << "l";
            }
            LatexText << "} \\hline \n";
            LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
            
            LatexText << NumberOfRadioTracerOrGeom << "}{c}{ ";
            LatexText << " \\textbf{Radiotracers}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
            
            for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
            {
                LatexText << "   & \\textbf{" << RadiotracerList[A]<< "}";
            }
            
            for ( int srcInc = 0; srcInc < 1 ; srcInc++  ){ // for just the source IntakeIntoBody
                
                for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                    
                    LatexText << "       \\\\\\hline \n";
                    LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{Body$\\to$" << TargetNamesToScore[A] << "}}       & " <<ValNm  << "                                         ";
                    for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]]["IntakeIntoBody"][TargetNamesToScore[A]] << "        ";
                    }
                    
                    if(RowNum == 4){
                        LatexText << "\\\\ \n                             & "<< CompareReferenceName << "                                        ";
                        
                        for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                            LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]]["IntakeIntoBody"][TargetNamesToScore[A]] << "        ";
                        }
                        
                        LatexText << "\\\\ \n                             & " << DiffSym << "\\tnote{a} " << "                                        ";
                        
                        for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                            double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]]["IntakeIntoBody"][TargetNamesToScore[A]];
                            double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]]["IntakeIntoBody"][TargetNamesToScore[A]];
                            double a3 = RelativeDifferenceCalculation(a1,a2);
                            
                            LatexText << "& " << std::fixed << std::setprecision(DiffDigitNum) << a3 << "        ";
                        }
                    }
                    
                    LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                    for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                        double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]]["IntakeIntoBody"][TargetNamesToScore[A]];
                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                    }
                }
            }
            
            LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
            LatexText << "\\begin{tablenotes}\\footnotesize\n";
            LatexText << "\\item[a] " << DiffExp << "\n";
            LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
            LatexText << "\\end{tablenotes}\n";
            LatexText << "\\end{threeparttable}\n";
            LatexText << "\\end{adjustbox}\n";
            LatexText << "\\label{tab:IntakeIntoBodyRadiotracersIn" << GeometryList[fd] << "}\n";
            LatexText << "\\end{table}\n\%\\end{sidewaystable}";
            LatexText << "\n\n";
            
            std::ofstream outfile(FileName , std::ios::app);
            if(outfile.is_open()){
                
                std::cout << "\nCreating file " << FileName << std::endl ;
                outfile << LatexText.str();
                outfile.close();
            }
        }
        for (int fd = 0 ; fd < RadiotracerList.size() ; fd++) {
            
            //if(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]].size()== 0){continue;}
            
            DataInitialization();
            std::string FileName = GraphsDirectoryPath + "ResRefLatexTables";
            std::ostringstream LatexText;
            
            int NumberOfRadioTracerOrGeom = GeometryList.size();
            int NumberOfCol = NumberOfRadioTracerOrGeom + 2;
            
            int RowNum = 2;
            //std::string ValNm = QuantitiesToScore;
            std::string ValNm = "DoseCalcs";
            
            if(GraphsData == "Reference_Result" && ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].size() != 0){
                RowNum = 4;
                //ValNm = CompareReferenceName;
            }
            
            std::cout << "\n\n                                                          ========= Creation of Latex Table For IntakeIntoBody, Geometry Data and radiotracers ========= "<< "\n" << std::endl;
            
            //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << GeometryRadiotracerSources[srcInc] << "\n\n";
            
            LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
            LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
            LatexText << "\\usepackage{multirow}\n";
            LatexText << "\n\n";
            
            LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                      << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";
            
            if(RowNum == 4){
                LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for IntakeIntoBody calculated in simulated geometries by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
            }else{
                LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for IntakeIntoBody calculated in simulated geometries  by DoseCalcs} \n";
            }
            
            LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                      << "\\begin{threeparttable}\n"
                      << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";
            
            for ( int A = 0; A < NumberOfCol ; A++  )
            {
                LatexText << "l";
            }
            LatexText << "} \\hline \n";
            LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
            
            LatexText << NumberOfRadioTracerOrGeom << "}{c}{ ";
            LatexText << " \\textbf{Phantoms}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
            
            for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
            {
                LatexText << "   & \\textbf{" << GeometryList[A]<< "}";
            }
            
            for ( int srcInc = 0; srcInc < 1 ; srcInc++  ){
                
                for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                    
                    LatexText << "       \\\\\\hline \n";
                    LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{Body$\\to$" << TargetNamesToScore[A] << "}}       & " <<ValNm  << "                                         ";
                    for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]]["IntakeIntoBody"][TargetNamesToScore[A]] << "        ";
                    }
                    
                    if(RowNum == 4){
                        LatexText << "\\\\ \n                             & "<< CompareReferenceName << "                                        ";
                        
                        for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                            LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]]["IntakeIntoBody"][TargetNamesToScore[A]] << "        ";
                        }
                        
                        LatexText << "\\\\ \n                             & " << DiffSym << "\\tnote{a} " << "                                        ";
                        
                        for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                            double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]]["IntakeIntoBody"][TargetNamesToScore[A]];
                            double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]]["IntakeIntoBody"][TargetNamesToScore[A]];
                            double a3 = RelativeDifferenceCalculation(a1,a2);
                            
                            LatexText << "& " << std::fixed << std::setprecision(DiffDigitNum) << a3 << "        ";
                        }
                    }
                    
                    LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";
                    
                    for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                        double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]]["IntakeIntoBody"][TargetNamesToScore[A]];
                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                    }
                }
            }
            
            LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
            LatexText << "\\begin{tablenotes}\\footnotesize\n";
            LatexText << "\\item[a] " << DiffExp << "\n";
            LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
            LatexText << "\\end{tablenotes}\n";
            LatexText << "\\end{threeparttable}\n";
            LatexText << "\\end{adjustbox}\n";
            LatexText << "\\label{tab:IntakeIntoBodyGeometriesFrom" << RadiotracerList[fd] << "}\n";
            LatexText << "\\end{table}\n\%\\end{sidewaystable}";
            LatexText << "\n\n";
            
            std::ofstream outfile(FileName , std::ios::app);
            if(outfile.is_open()){
                
                std::cout << "\nCreating file " << FileName << std::endl ;
                outfile << LatexText.str();
                outfile.close();
            }
        }
        
        // For value(RSD) in one line without method column (Tested for table for each source)

        RowNum = 1; // just one line value(RSD) and not two lines (value and RSD%)

        for (int sss = 0 ; sss < 3 ; sss++) { //0 self(sources) and 1 Cross(source-target combinations) // 2 for SourceTargets
            for (int fd = 0 ; fd < GeometryList.size() ; fd++) {

                if(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]].size()== 0){continue;}

                int NumberOfRadioTracerOrGeom = RadiotracerList.size();
                int NumberOfCol = NumberOfRadioTracerOrGeom + 1;

                std::cout << "\n\n                                                          ========= Creation of Latex Table For Radiotracer Data in geometries ========= "<< "\n" << std::endl;

                DataInitialization();
                std::string FileName = GraphsDirectoryPath + "ResRefLatexTables";

                if(sss==2){
                    for ( int srcInc = 0; srcInc < GeometryRadiotracerSources.size() ; srcInc++  ){

                        std::ostringstream LatexText;

                        std::cout << "\n\n                                                          ========= Creation of Latex Table For Geometry Data and radiotracers ========= "<< "\n" << std::endl;

                        //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << GeometryRadiotracerSources[srcInc] << "\n\n";

                        LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
                        LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
                        LatexText << "\\usepackage{multirow}\n";
                        LatexText << "\n\n";

                        LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                                  << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";

                        LatexText << "\%for long tatbles ------------------------------ \n\%\\begin{tiny}\n\% \\begin{longtable}{@{\\extracolsep{\\fill}}*{"<<NumberOfCol<<"}{l}}\n\%------------------------------\n";

                        if(RowNum == 3){
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << GeometryList[fd] << " for internal irradiation from "<<GeometryRadiotracerSources[srcInc]<<" calculated in phantoms target regions by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                        }else{
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << GeometryList[fd] << " for internal irradiation from  "<<GeometryRadiotracerSources[srcInc]<<" calculated in phantoms target regions by DoseCalcs} \n";
                        }

                        LatexText << "\%for long tatbles ------------------------------ \n\%\\\\ \n\%\\hline\n\% \\endfirsthead\n\% \\caption[]{(\\textit{continued})}\\\\\n\% \\hline\n\%------------------------------\n";

                        LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                                  << "\\begin{threeparttable}\n"
                                  << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";

                        for ( int A = 0; A < NumberOfCol ; A++  )
                        {
                            LatexText << "l";
                        }
                        LatexText << "} \\hline \n";

                        LatexText << "\\multicolumn{1}{l}{\\multirow{1}{*}{\\textbf{Target}}} & \\multicolumn{";
                        LatexText << NumberOfRadioTracerOrGeom << "}{l}{ ";
                        LatexText << " \\textbf{Radionuclides}}       \\\\ \\cline{2-"<< NumberOfCol << "}\n                 \\multicolumn{1}{l}{}       ";
                        for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                        {
                            LatexText << "   & \\textbf{" << RadiotracerList[A]<< "}";
                        }

                        LatexText << "\%for long tatbles ------------------------------ \n\%\\\\\\hline \n\%\\endhead \n\%\\hline \n\%\\multicolumn{6}{r}{(\\textit{Continued on next page})} \\\\ \n\%\\endfoot \n\%\\hline \n\%\\endlastfoot\n";

                        LatexText << "\n\%\\multicolumn{1}{l}{\\multirow{1}{*}{\\textbf{Target}}} & \\multicolumn{";
                        LatexText << NumberOfRadioTracerOrGeom << "}{l}{ ";
                        LatexText << "\n\% \\textbf{Radionuclides}}       \\\\ \\cline{2-"<< NumberOfCol << "}\n\%                 \\multicolumn{1}{l}{      ";
                        for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                        {
                            LatexText << "   & \\textbf{" << RadiotracerList[A]<< "}";
                        }
                        LatexText << "\n\%------------------------------\n";

                        for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){

                            LatexText << "       \\\\\\hline \n";
                            LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<TargetNamesToScore[A]<< "}}         ";

                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];

                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << " ("<< std::scientific << std::setprecision(QuantityDigitNum) << a1 << "\\%)";
                            }

                            if(RowNum == 3){
                                LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";

                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    double a3 = RelativeDifferenceCalculation(a1,a2);
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                }
                            }

                        }

                        LatexText << "\%for long tables ------------------------------ \n\%\\end{longtable} \n\%\\footnotetext[1]{Ratio} \n\%\\footnotetext[2]{Relative standard deviation}\n\%\\end{tiny}\n\%------------------------------\n";

                        LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
                        //LatexText << "\\begin{tablenotes}\\footnotesize\n";
                        //LatexText << "\\item[a] " << DiffExp << "\n";
                        //LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                        //LatexText << "\\end{tablenotes}\n";
                        LatexText << "\\end{threeparttable}\n";
                        LatexText << "\\end{adjustbox}\n";
                        LatexText << "\\label{tab:"<<GeometryRadiotracerSources[srcInc]<<"RadiotracersIn" << GeometryList[fd] << "}\n";
                        LatexText << "\\end{table}\n\%\\end{sidewaystable}";
                        LatexText << "\n\n";

                        std::ofstream outfile(FileName , std::ios::app);
                        if(outfile.is_open()){

                            std::cout << "\nCreating file " << FileName << std::endl ;
                            outfile << LatexText.str();
                            outfile.close();
                        }
                    }
                }
                else{
                    std::ostringstream LatexText;

                    //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << GeometryRadiotracerSources[srcInc] << "\n\n";

                    LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
                    LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
                    LatexText << "\\usepackage{multirow}\n";
                    LatexText << "\n\n";

                    LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                              << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";

                    if(sss==0){
                        if(RowNum == 3){
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of radiotracers for each source region calculated in " << GeometryList[fd] << " by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                        }else{
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of radiotracers for each source region calculated in " << GeometryList[fd] << " by DoseCalcs} \n";
                        }
                    }
                    else if(sss==1){
                        if(RowNum == 3){
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of radiotracers for each source-target combination calculated in " << GeometryList[fd] << " by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                        }else{
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of radiotracers for each source-target combination calculated in " << GeometryList[fd] << " by DoseCalcs} \n";
                        }
                    }


                    LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                              << "\\begin{threeparttable}\n"
                              << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";

                    for ( int A = 0; A < NumberOfCol ; A++  )
                    {
                        LatexText << "l";
                    }
                    LatexText << "} \\hline \n";

                    if(sss==0){
                        LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source Region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                    }
                    else if(sss==1){
                        LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                    }

                    LatexText << NumberOfRadioTracerOrGeom << "}{c}{ ";
                    LatexText << " \\textbf{Radiotracers}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";

                    for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                    {
                        LatexText << "   & \\textbf{" << RadiotracerList[A]<< "}";
                    }

                    for ( int srcInc = 0; srcInc < GeometryRadiotracerSources.size() ; srcInc++  ){

                        //if(GeometryRadiotracerSources[srcInc] == "IntakeIntoBody"){ continue;}
                        if(sss==0){
                            LatexText << "       \\\\\\hline \n";
                            LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<GeometryRadiotracerSources[srcInc] << "}}       & " <<ValNm  << "                                         ";
                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]] << "        ";
                            }

                            if(RowNum == 3){
                                LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";

                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]];
                                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]];
                                    double a3 = RelativeDifferenceCalculation(a1,a2);
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                }
                            }

                            LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";

                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]];
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                            }
                        }
                        else if(sss==1){
                            for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){

                                LatexText << "       \\\\\\hline \n";
                                LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<TargetNamesToScore[A]<<"$\\leftarrow$" <<GeometryRadiotracerSources[srcInc] << "}}       & " <<ValNm  << "                                         ";

                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << "        ";
                                }

                                if(RowNum == 3){
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";

                                    for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                        double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                        double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                        double a3 = RelativeDifferenceCalculation(a1,a2);
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                    }
                                }

                                LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";

                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                                }
                            }
                        }
                    }

                    LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
                    LatexText << "\\begin{tablenotes}\\footnotesize\n";
                    LatexText << "\\item[a] " << DiffExp << "\n";
                    LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                    LatexText << "\\end{tablenotes}\n";
                    LatexText << "\\end{threeparttable}\n";
                    LatexText << "\\end{adjustbox}\n";
                    if(sss==0){
                        LatexText << "\\label{tab:SourceRadiotracersIn" << GeometryList[fd] << "}\n";
                    }
                    else if(sss==1){
                        LatexText << "\\label{tab:SourceTargetRadiotracersIn" << GeometryList[fd] << "}\n";
                    }
                    LatexText << "\\end{table}\n\%\\end{sidewaystable}";
                    LatexText << "\n\n";

                    std::ofstream outfile(FileName , std::ios::app);
                    if(outfile.is_open()){

                        std::cout << "\nCreating file " << FileName << std::endl ;
                        outfile << LatexText.str();
                        outfile.close();
                    }
                }
            }
            for (int fd = 0 ; fd < RadiotracerList.size() ; fd++) {

                int NumberOfRadioTracerOrGeom = GeometryList.size();
                int NumberOfCol = NumberOfRadioTracerOrGeom + 1;

                std::cout << "\n\n                                                          ========= Creation of Latex Table For Radiotracer Data for radiotracers ========= "<< "\n" << std::endl;

                DataInitialization();
                std::string FileName = GraphsDirectoryPath + "ResRefLatexTables";

                if(sss==2){
                    for ( int srcInc = 0; srcInc < GeometryRadiotracerSources.size() ; srcInc++  ){

                        //if(GeometryRadiotracerSources[srcInc] == "IntakeIntoBody"){ continue;}

                        std::ostringstream LatexText;

                        std::cout << "\n\n                                                          ========= Creation of Latex Table For Geometry Data and radiotracers ========= "<< "\n" << std::endl;

                        //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << GeometryRadiotracerSources[srcInc] << "\n\n";

                        LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
                        LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
                        LatexText << "\\usepackage{multirow}\n";
                        LatexText << "\n\n";

                        LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                                  << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";

                        LatexText << "\%for long tatbles ------------------------------ \n\%\\begin{tiny}\n\% \\begin{longtable}{@{\\extracolsep{\\fill}}*{"<<NumberOfCol<<"}{l}}\n\%------------------------------\n";

                        if(RowNum == 3){
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for internal irradiation from "<<GeometryRadiotracerSources[srcInc]<<" calculated in phantoms target regions by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                        }else{
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for internal irradiation from "<<GeometryRadiotracerSources[srcInc]<<" calculated in phantoms target regions  by DoseCalcs} \n";
                        }

                        LatexText << "\%for long tatbles ------------------------------ \n\%\\\\ \n\%\\hline\n\% \\endfirsthead\n\% \\caption[]{(\\textit{continued})}\\\\\n\% \\hline\n\%------------------------------\n";

                        LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                                  << "\\begin{threeparttable}\n"
                                  << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";

                        for ( int A = 0; A < NumberOfCol ; A++  )
                        {
                            LatexText << "l";
                        }
                        LatexText << "} \\hline \n";

                        LatexText << "\\multicolumn{1}{l}{\\multirow{1}{*}{\\textbf{Target}}} & \\multicolumn{";
                        LatexText << NumberOfRadioTracerOrGeom << "}{l}{ ";
                        LatexText << " \\textbf{Phantoms}}       \\\\ \\cline{2-"<< NumberOfCol << "}\n                 \\multicolumn{1}{l}{}      ";
                        for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                        {
                            LatexText << "   & \\textbf{" << GeometryList[A]<< "}";
                        }

                        LatexText << "\%for long tatbles ------------------------------ \n\%\\\\\\hline \n\%\\endhead \n\%\\hline \n\%\\multicolumn{6}{r}{(\\textit{Continued on next page})} \\\\ \n\%\\endfoot \n\%\\hline \n\%\\endlastfoot\n";

                        LatexText << "\n\%\\multicolumn{1}{l}{\\multirow{1}{*}{\\textbf{Target}}} & \\multicolumn{";
                        LatexText << NumberOfRadioTracerOrGeom << "}{l}{ ";
                        LatexText << "\n\% \\textbf{Phantoms}}       \\\\ \\cline{2-"<< NumberOfCol << "}\n\%                 \\multicolumn{1}{l}{}    ";
                        for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                        {
                            LatexText << "   & \\textbf{" << GeometryList[A]<< "}";
                        }
                        LatexText << "\n\%------------------------------\n";

                        for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){

                            LatexText << "       \\\\\\hline \n";
                            LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<TargetNamesToScore[A]<< "}}                                            ";

                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << " ("<<std::scientific << std::setprecision(QuantityDigitNum) << a1 << "\\%)        ";
                            }

                            if(RowNum == 3){
                                LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";

                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    double a3 = RelativeDifferenceCalculation(a1,a2);
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                }
                            }

                            //LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";

                            //for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                            //    double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                            //    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                            //}
                        }

                        LatexText << "\%for long tatbles ------------------------------ \n\%\\end{longtable} \n\%\\footnotetext[1]{Ratio} \n\%\\footnotetext[2]{Relative standard deviation}\n\%\\end{tiny}\n\%------------------------------\n";

                        LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
                        //LatexText << "\\begin{tablenotes}\\footnotesize\n";
                        //LatexText << "\\item[a] " << DiffExp << "\n";
                        //LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                        //LatexText << "\\end{tablenotes}\n";
                        LatexText << "\\end{threeparttable}\n";
                        LatexText << "\\end{adjustbox}\n";
                        LatexText << "\\label{tab:"<<GeometryRadiotracerSources[srcInc]<<"GeometriesFrom" << RadiotracerList[fd] << "}\n";
                        LatexText << "\\end{table}\n\%\\end{sidewaystable}";
                        LatexText << "\n\n";

                        std::ofstream outfile(FileName , std::ios::app);
                        if(outfile.is_open()){

                            std::cout << "\nCreating file " << FileName << std::endl ;
                            outfile << LatexText.str();
                            outfile.close();
                        }
                    }
                }
                else{

                    std::ostringstream LatexText;

                    std::cout << "\n\n                                                          ========= Creation of Latex Table For Geometry Data and radiotracers ========= "<< "\n" << std::endl;

                    //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << GeometryRadiotracerSources[srcInc] << "\n\n";

                    LatexText << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
                    LatexText << "\\usepackage{threeparttable}\n\%\\usepackage{longtable}\n\%\\usepackage{rotating}\n";
                    LatexText << "\\usepackage{multirow}\n";
                    LatexText << "\n\n";

                    LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                              << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";

                    LatexText << "\%for long tatbles ------------------------------ \n\%\\begin{tiny}\n\% \\begin{longtable}{@{\\extracolsep{\\fill}}*{"<<NumberOfCol<<"}{l}}\n\%------------------------------\n";

                    if(sss==0){
                        if(RowNum == 3){
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for each source region calculated in simulated geometries by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                        }else{
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for each source region calculated in simulated geometries  by DoseCalcs} \n";
                        }
                    }
                    else if(sss==1){
                        if(RowNum == 3){
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for each source-target combination calculated in simulated geometries by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                        }else{
                            LatexText << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << RadiotracerList[fd] << " for each source-target combination calculated in simulated geometries  by DoseCalcs} \n";
                        }
                    }

                    LatexText << "\%for long tatbles ------------------------------ \n\%\\\\ \n\%\\hline\n\% \\endfirsthead\n\% \\caption[]{(\\textit{continued})}\\\\\n\% \\hline\n\%------------------------------\n";

                    LatexText << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                              << "\\begin{threeparttable}\n"
                              << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";

                    for ( int A = 0; A < NumberOfCol ; A++  )
                    {
                        LatexText << "l";
                    }
                    LatexText << "} \\hline \n";

                    if(sss==0){
                        LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source Region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                    }
                    else if(sss==1){
                        LatexText << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                    }
                    LatexText << NumberOfRadioTracerOrGeom << "}{c}{ ";
                    LatexText << " \\textbf{Phantoms}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                    for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                    {
                        LatexText << "   & \\textbf{" << GeometryList[A]<< "}";
                    }

                    LatexText << "\%for long tatbles ------------------------------ \n\%\\\\\\hline \n\%\\endhead \n\%\\hline \n\%\\multicolumn{6}{r}{(\\textit{Continued on next page})} \\\\ \n\%\\endfoot \n\%\\hline \n\%\\endlastfoot\n";
                    if(sss==0){
                        LatexText << "\n\%\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source Region}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                    }
                    else if(sss==1){
                        LatexText << "\n\%\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                    }
                    LatexText << NumberOfRadioTracerOrGeom << "}{c}{ ";
                    LatexText << "\n\% \\textbf{Phantoms}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n\%                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                    for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                    {
                        LatexText << "   & \\textbf{" << GeometryList[A]<< "}";
                    }
                    LatexText << "\n\%------------------------------\n";

                    for ( int srcInc = 0; srcInc < GeometryRadiotracerSources.size() ; srcInc++  ){
                        //if(GeometryRadiotracerSources[srcInc] == "IntakeIntoBody"){ continue;}

                        if(sss==0){
                            LatexText << "       \\\\\\hline \n";
                            LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<GeometryRadiotracerSources[srcInc] << "}}       & " <<ValNm  << "                                         ";

                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]] << "        ";
                            }

                            if(RowNum == 3){
                                LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";

                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]];
                                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]];
                                    double a3 = RelativeDifferenceCalculation(a1,a2);
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                }
                            }

                            LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";

                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][GeometryRadiotracerSources[srcInc]];
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                            }
                        }
                        else if(sss==2){

                        }
                        else if(sss==1){
                            for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){

                                LatexText << "       \\\\\\hline \n";
                                LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<TargetNamesToScore[A]<<"$\\leftarrow$" <<GeometryRadiotracerSources[srcInc] << "}}       & " <<ValNm  << "                                         ";

                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << "        ";
                                }

                                if(RowNum == 3){
                                    LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";

                                    for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                        double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                        double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                        double a3 = RelativeDifferenceCalculation(a1,a2);
                                        LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                    }
                                }

                                LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";

                                for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                    double a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                                }
                            }
                        }
                    }

                    LatexText << "\%for long tatbles ------------------------------ \n\%\\end{longtable} \n\%\\footnotetext[1]{Ratio} \n\%\\footnotetext[2]{Relative standard deviation}\n\%\\end{tiny}\n\%------------------------------\n";

                    LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
                    LatexText << "\\begin{tablenotes}\\footnotesize\n";
                    LatexText << "\\item[a] " << DiffExp << "\n";
                    LatexText << "\\item[b] Relative Standard Deviation (\\%)\n";
                    LatexText << "\\end{tablenotes}\n";
                    LatexText << "\\end{threeparttable}\n";
                    LatexText << "\\end{adjustbox}\n";
                    if(sss==0){
                        LatexText << "\\label{tab:SourceGeometriesFrom" << RadiotracerList[fd] << "}\n";
                    }
                    else if(sss==1){
                        LatexText << "\\label{tab:SourceTargetGeometriesFrom" << RadiotracerList[fd] << "}\n";
                    }
                    LatexText << "\\end{table}\n\%\\end{sidewaystable}";
                    LatexText << "\n\n";

                    std::ofstream outfile(FileName , std::ios::app);
                    if(outfile.is_open()){

                        std::cout << "\nCreating file " << FileName << std::endl ;
                        outfile << LatexText.str();
                        outfile.close();
                    }
                }
            }
        }

        // For Comparison factor simple table

        for (int fd = 0 ; fd < RadiotracerList.size() ; fd++) {

            DataInitialization();
            std::string FileName = GraphsDirectoryPath + "ResRefLatexTables";

            std::ostringstream LatexText;

            LatexText << "\\begin{table}[H] \n\%\\begin{sidewaystable}\n"
                      << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n"
                      << "\\caption{" << DifferenceMethod << " for " << RadiotracerList[fd] << " calculated in geometries by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n"
                      << "\\begin{tabular}{llllll} \n";
            LatexText << "\\hline \n";

            for ( int dss = 0; dss < 2 ; dss++  ){
                LatexText << "\\textbf{Target$\\leftarrow$Source} ";
                for ( int DD = 0; DD < GeometryList.size() ; DD++  ){
                    LatexText << " & \\textbf{" << GeometryList[DD] << "}";
                }
                if(dss != 1){
                    LatexText << " & ";
                }
            }

            int ja=1;
            for ( int A = 0; A < SourceNamesToScore.size() ; A++  ){
                for ( int B = 0; B < TargetNamesToScore.size() ; B++  ){
                    if(TargetNamesToScore[B] == "World"){
                        continue;
                    }

                    if(ja==2){
                        LatexText << " &  " ;
                        ja=1;
                    }else{
                        LatexText << "       \\\\\\hline\n";
                        ja=2;
                    }

                    LatexText << TargetNamesToScore[B] << "$\\leftarrow$" << SourceNamesToScore[A] ;//
                    //<<"     & " << std::fixed << std::setprecision(GeometryVarDigitNum) << GeometryRegionVariableValue[GeometrySymbol]["Volume"][OrgansNameVector[A]] <<"      & "<< GeometryRegionVariableValue[GeometrySymbol]["Mass"][OrgansNameVector[A]] <<"     & " << GeometryRegionVariableValue[GeometrySymbol]["Density"][OrgansNameVector[A]];

                    for ( int DD = 0; DD < GeometryList.size() ; DD++  ){
                        double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[DD]][RadiotracerList[fd]][SourceNamesToScore[A]][TargetNamesToScore[B]];
                        double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[DD]][RadiotracerList[fd]][SourceNamesToScore[A]][TargetNamesToScore[B]];
                        double a3 = RelativeDifferenceCalculation(a1,a2);

                        LatexText <<"     & " << std::fixed << std::setprecision(DiffDigitNum) << a3;

                    }
                }
            }

            LatexText << "\\\\ \\hline\n\\end{tabular} \n\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n";
            LatexText << "\\label{"<< RadiotracerList[fd] <<"GeometriesComparison}\n";
            LatexText << "\\end{table}\n\%\\end{sidewaystable}";
            LatexText << "\n\n\n\n\n";

            std::ofstream outfile(FileName , std::ios::app);
            if(outfile.is_open()){

                std::cout << "\nCreating file " << FileName << std::endl ;
                outfile << LatexText.str();
                outfile.close();
            }
        }

    }

}
void G4DoseCalcsAnalysis::GenerateLatexTablesForRadioTracerResult(){

    bool islongtable = true;
    bool islandscape = true;
    int NumberOfLines = 20 ;

    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;

    if(ResultQuantityGeometryRadioTracerSourceTargetValues.size()== 0){return;}

    for (int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {

        QuantitiesToScore = QuantityNamesToScore[gg];
        if(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].size()== 0){continue;}

        int RowNum = 2;
        std::string ValNm = "DoseCalcs";

        if(GraphsData == "Reference_Result" && ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].size() != 0){
            RowNum = 3;
            //ValNm = CompareReferenceName;
        }

        std::vector <std::string> GeomRadiotList;
        for (int CC = 0 ; CC < 2 ; CC++){ // Geometry or Radiotracer

            int NumberOfRadioTracerOrGeom;

            if(CC == 0){
                NumberOfRadioTracerOrGeom = RadiotracerList.size();
                GeomRadiotList = GeometryList;
            }else{
                NumberOfRadioTracerOrGeom = GeomRadiotList.size();
                GeomRadiotList = RadiotracerList;
            }

            for (int fd = 0 ; fd < GeomRadiotList.size() ; fd++) {

                int NumberOfCol = NumberOfRadioTracerOrGeom + 2;

                std::cout << "\n\n                                                          ========= Creation of Latex Table For Radiotracer Data in geometries ========= "<< "\n" << std::endl;

                DataInitialization();
                std::string FileName = GraphsDirectoryPath + "ResRefRadiotracerLatexTables";

                std::ostringstream PackageImp;

                PackageImp << "\n\n\n\n\\usepackage{booktabs, makecell, graphicx, caption, subcaption}\n";
                if(islongtable){
                    PackageImp << "\\usepackage{longtable}\n";
                }else{
                    PackageImp << "\\usepackage{threeparttable}\n"
                                  "\\usepackage{multirow}\n";
                }
                if(islandscape){
                    PackageImp << "\\usepackage{lscape}\n";
                }

                PackageImp << "\n\n";

                int ii = 0;
                for ( int srcInc = 0; srcInc < GeometryRadiotracerSources.size() ; srcInc++  ){

                    std::ostringstream Header1;

                    //LatexText << "=============================== Latex Table for particle " << ParticleName << " and Source " << GeometryRadiotracerSources[srcInc] << "\n\n";

                    if(islandscape){
                        Header1 << "\\begin{landscape}\n";
                    }
                    if(islongtable){
                        Header1 << "\%\\begin{tiny}\n"
                                   "\\begin{longtable}{";
                        for ( int A = 0; A < NumberOfCol ; A++  )
                        {
                            Header1 << "l";
                        }
                        Header1 << "} \n";
                    }else{
                        Header1 << "\\begin{table}[H]\n"
                                << "\\centering \n\%\\caption*{Table 3: (\\textit{continued})}\n";
                    }

                    if(RowNum == 3){
                        Header1 << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << GeomRadiotList[fd] << " for internal irradiation from "<<GeometryRadiotracerSources[srcInc]<<" calculated by DoseCalcs and compared to the " << CompareReferenceName << " reference} \n";
                    }else{
                        Header1 << "\\caption{" << QuantityUnit[QuantitiesToScore] << " values of " << GeomRadiotList[fd] << " for internal irradiation from  "<<GeometryRadiotracerSources[srcInc]<<" calculated by DoseCalcs} \n";
                    }

                    if(islongtable){
                        Header1 << "\\endfirsthead\n  "
                                   "\\hline \n"
                                   "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{"<< NumberOfRadioTracerOrGeom << "}{c}{ ";
                        Header1 << "\n \\textbf{PhantomsOrRadionuclides}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                        for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                        {
                            if(CC==0){
                                Header1 << "   & \\textbf{" << RadiotracerList[A]<< "}";
                            }
                            else{
                                Header1 << "   & \\textbf{" << GeometryList[A]<< "}";
                            }
                        }
                        Header1 << "\\\\\\hline \n \\endhead % all the lines above this will be repeated on every page \n";
                        Header1 <<"\\multicolumn{"<< NumberOfCol << "}{r}{{Continued on next page \\ldots}} \\\\ \n \\endfoot \n";
                        Header1 <<"\\hline \\\\ \n \\multicolumn{"<< NumberOfCol << "}{l}{\\textbf{RSD(\\%): Relative standard deviation (\\%)}} \n \\endlastfoot \n\n";

                        Header1 <<"\\hline \n"
                                  "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{"<< NumberOfRadioTracerOrGeom << "}{c}{ ";
                        Header1 <<"\n \\textbf{PhantomsOrRadionuclides}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                        for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                        {
                            if(CC==0){
                                Header1 << "   & \\textbf{" << RadiotracerList[A]<< "}";
                            }
                            else{
                                Header1 << "   & \\textbf{" << GeometryList[A]<< "}";
                            }
                        }
                        Header1 << "\\label{" << GeometryRadiotracerSources[srcInc] << "} ";

                    }else{
                        Header1 << "\\begin{adjustbox}{width=\\columnwidth,center}\n"
                                << "\\begin{threeparttable}\n"
                                << "\%\\tiny \%to make any table size fill one page\n\\begin{tabular}{";

                        for ( int A = 0; A < NumberOfCol ; A++  )
                        {
                            Header1 << "l";
                        }
                        Header1 << "} \\hline \n";

                        Header1 << "\\multicolumn{1}{c}{\\multirow{2}{*}{\\textbf{Source$\\to$Target}}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{";
                        Header1 << NumberOfRadioTracerOrGeom << "}{c}{ ";
                        Header1 << " \\textbf{PhantomsOrRadionuclides}}       \\\\ \\cline{3-"<< NumberOfCol << "}\n                 \\multicolumn{1}{c}{}                             & \\multicolumn{1}{c}{}                        ";
                        for ( int A = 0; A < NumberOfRadioTracerOrGeom ; A++  )
                        {
                            if(CC==0){
                                Header1 << "   & \\textbf{" << RadiotracerList[A]<< "}";
                            }
                            else{
                                Header1 << "   & \\textbf{" << GeometryList[A]<< "}";
                            }
                        }
                    }

                    //Header1.str().replace("",GeometryRadiotracerSources[srcInc]);

                    std::ostringstream Footer1;
                    if(islongtable){
                        Footer1 << "\n\\end{longtable} \n "
                                   "\%\\end{tiny}\n\n";
                    }else{
                        Footer1 << "\\\\ \\hline\n"
                                   "\\end{tabular} \n"
                                   "\%\\begin{flushright}\\textit{Continued on next page}\\end{flushright}\n"
                                   "\\begin{tablenotes}\\footnotesize\n";
                        if(RowNum == 3){
                            Footer1 << "\\item[a] " << DiffExp << "\n";
                        }
                        Footer1 << "\\item[b] Relative Standard Deviation (\\%)\n"
                                   "\\end{tablenotes}\n"
                                   "\\end{threeparttable}\n"
                                   "\\end{adjustbox}\n"
                                   "\\label{tab:"<<GeometryRadiotracerSources[srcInc]<<"RadiotracerGeometry}\n\\end{table}\n";
                    }
                    if(islandscape){
                        Footer1 << "\\end{landscape}\n";
                    }

                    Footer1 << "\n\n";

                    FileName = GraphsDirectoryPath + QuantitiesToScore+"_"+GeometryRadiotracerSources[srcInc]+"_RadiotracerGeometryLatexTable.tex";
                    std::ostringstream LatexText;
                    std::ostringstream LatexTextForAllComb;

                    //std::ostringstream nullString ; nullString  << "";

                    //std::cout << "\n\n\n\n\n\n New Source " << std::endl ;

                    for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){

                        //if(GeometryRadiotracerSources[srcInc] == "IntakeIntoBody"){ continue;}

                        LatexText << "       \\\\\\hline \n";
                        LatexText << "\\multirow{"<< RowNum <<"}{*}{\\textbf{"<<TargetNamesToScore[A]<<"$\\leftarrow$" <<GeometryRadiotracerSources[srcInc] << "}}       & " <<ValNm  << "                                         ";

                        for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                            if(CC == 0){
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeomRadiotList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << "        ";
                            }else{
                                LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][GeomRadiotList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << "        ";
                            }
                        }

                        if(RowNum == 3){
                            LatexText << "\\\\ \n                             & "<< CompareReferenceName << "(" << DiffSym << "\\tnote{a} " << ")                      ";

                            for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                                double a1;
                                double a2;
                                if(CC==0){
                                    a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeomRadiotList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeomRadiotList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                }else{
                                    a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][GeomRadiotList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                    a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][GeomRadiotList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                                }
                                double a3 = RelativeDifferenceCalculation(a1,a2);
                                if(CC==0){
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeomRadiotList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                }
                                else{
                                    LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometryList[B]][RadiotracerList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]] << " (" << std::fixed << std::setprecision(DiffDigitNum) << a3 << ")        ";
                                }
                            }
                        }

                        LatexText << "\\\\ \n                             & "<< "RSD(\\%)\\tnote{b} " << "                                        ";

                        for ( int B = 0; B < NumberOfRadioTracerOrGeom ; B++  ){
                            double a1;
                            if(CC==0){
                                a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeomRadiotList[fd]][RadiotracerList[B]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                            }else{
                                a1 = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[QuantitiesToScore][GeometryList[B]][GeomRadiotList[fd]][GeometryRadiotracerSources[srcInc]][TargetNamesToScore[A]];
                            }
                            LatexText << "& " << std::scientific << std::setprecision(QuantityDigitNum) << a1 << "        ";
                        }

                        ii++;
                        if(NumberOfLines == ii && !islongtable){

                            std::cout << " NumberOfLines " << NumberOfLines << " ii " << ii << " A " << A << std::endl ;

                            ii = 0;
                            std::ofstream outfile(FileName , std::ios::app);
                            if(outfile.is_open()){

                                std::cout << "\nCreating file " << FileName << std::endl ;
                                outfile << PackageImp.str();
                                outfile << Header1.str();
                                outfile << LatexText.str();
                                outfile << Footer1.str();
                                outfile.close();
                            }
                            LatexText.str("");
                            LatexText.clear();
                        }

                        if(NumberOfLines == -1 || islongtable){
                            LatexTextForAllComb << LatexText.str();
                            LatexText.str("");
                            LatexText.clear();
                        }
                    }
                    if( ii < NumberOfLines && ii > 0 && NumberOfLines > 0 && !islongtable){

                        std::ofstream outfile(FileName , std::ios::app);
                        if(outfile.is_open()){

                            std::cout << "\nCreating file " << FileName << std::endl ;
                            outfile << PackageImp.str();
                            outfile << Header1.str();
                            outfile << LatexText.str();
                            outfile << Footer1.str();
                            outfile.close();
                        }
                    }

                    //if(NumberOfLines == -1 || islongtable || CC==1){
                       //std::ofstream outfile(GraphsDirectoryPath +"RES", std::ios::app);

                    if(NumberOfLines == -1 || islongtable){
                        std::ofstream outfile(FileName , std::ios::app);
                        if(outfile.is_open()){

                            std::cout << "\n 11Creating file " << FileName << std::endl ;
                            outfile << PackageImp.str();
                            outfile << Header1.str();
                            outfile << LatexTextForAllComb.str();
                            outfile << Footer1.str();
                            outfile.close();
                        }
                    }

                    LatexText.str("");
                    LatexText.clear();

                    LatexTextForAllComb.str("");
                    LatexTextForAllComb.clear();
                }
            }
        }
    }
}

void G4DoseCalcsAnalysis::GenerateGraphFromROOTGraphDATA(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    std::map<std::string,std::map<double,double>> DataForROOTMap;
    std::vector<double> EnergiesVec;
    std::vector<std::string> XLabelsVec;
    std::vector<int> XLabelsIDVec;
    std::string FileToRead = MacrosStartingFile+"/ROOTGraphData";
    
    std::ifstream file1(FileToRead.c_str() , std::ios::binary);
    
    std::string title, xlbl, ylbl;
    int graphtype = 0;
    
    if(file1.is_open()){
        
        std::cout << "\nReading file " << MacrosStartingFile.c_str() << std::endl ;
        
        std::string ParameterName, line, word;
        double val;
        int dataid = 0 ;
        
        while (getline(file1, line)) {
            
            std::istringstream LineString(line);
            ParameterName = "";
            
            if(dataid == 0 ){
                
                LineString >> graphtype >> title >> xlbl >> ylbl >> UseLogE >> UseLogVariable >> UseGridXY >> LegendPos >> AddErrorBarInGraphs ;
                std::cout << " UseLogE " << UseLogE <<
                             " UseLogVariable " << UseLogVariable <<
                             " UseGridXY " << UseGridXY <<
                             " LegendPos " << LegendPos <<
                             " AddErrorBarInGraphs " << AddErrorBarInGraphs <<
                             std::endl ;
                dataid++ ;
                continue;
            }
            else if(dataid == 1){
                
                LineString >> ParameterName;
                
                if(graphtype == 0){
                    while(LineString >> val ){
                        EnergiesVec.push_back(val);
                    }
                }
                else if(graphtype == 1){
                    while(LineString >> val ){
                        EnergiesVec.push_back(val);
                        LineString >> ParameterName ;
                        XLabelsVec.push_back(ParameterName);
                        int vv = (int)val;
                        XLabelsIDVec.push_back(vv);
                    }
                }
                G4cout << "  Readed energies are : " << G4endl ;
                for (G4int gg = 0 ; gg < EnergiesVec.size() ; gg++) {
                    G4cout << EnergiesVec[gg] << G4endl;
                }
                dataid++;
                continue;
                
            }else{
                
                LineString >> ParameterName;
                int ii = 0;
                while(LineString >> val ){
                    DataForROOTMap[ParameterName][EnergiesVec[ii]] = val;
                    G4cout << " ParameterName : " << ParameterName << " EnergiesVec[ii] : " << EnergiesVec[ii] << " val : " << val << G4endl ;
                    ii++;
                }
                continue;
            }
        }
    }
    
    std::string FileName ;
    TCanvas* ResCanvas ;
    TMultiGraph *mg ;
    TLegend *leg ;
    
    FileName = MacrosStartingFile+"/Graph.root";
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    
    int jj = 0;
    // iterations on source name
    for ( auto Bbeg = DataForROOTMap.begin(); Bbeg != DataForROOTMap.end(); ++Bbeg  )
    {
        std::string GraphLabel = Bbeg->first;
        std::cout << "\n\n-------------------------------- For Graph label "<< GraphLabel << "----------------------" << std::endl ;
        
        int jh = 0; // to fill the dataArray that will be used by the graph
        
        int NumOfEne = Bbeg->second.size();
        double xi[NumOfEne], yi[NumOfEne];
        
        for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg  )
        {
            xi[jh] = Cbeg->first;
            yi[jh] = Cbeg->second;
            
            //std::cout << " PARTICLE_NAME " << PARTICLE_NAME << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;
            
            jh++;
        }
        
        TGraph* graph;
        graph = CreateGraph (NumOfEne, xi, yi);
        graph->SetName(GraphLabel.c_str());
        graph->SetTitle(graph->GetName());
        graph->SetLineWidth(1);
        graph->SetMarkerStyle(50+jj);
        graph->SetMarkerColor(jj+1);
        graph->SetLineColor(jj+1);
        
        /*
        if(graphtype == 1){
            for(int aa = 0; aa < NumOfEne; aa++){
                int i = xi[aa];
                graph->GetXaxis()->ChangeLabel(i,-1,-1,-1,-1,-1,XLabelsVec[aa].c_str());;
            }
        }
        */
        mg->Add(graph);
        leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
        
        jj++;
    }
    
    if(UseLogVariable == "yes"){
        gPad->SetLogy(1);
    }
    if(UseLogE == "yes"){
        gPad->SetLogx(1);
    }
    
    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }
    
    std::string multiGraphTitle = title + " ; " + xlbl + " ; " + ylbl;
    
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = " ; " + xlbl + " ; " + ylbl;
        mg->SetTitle(multiGraphTitle.c_str());
    }
    
    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);
    
    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();
    
    if(graphtype == 1){
        /*
        TText* t = new TText();
        t->SetTextAlign(0);
        t->SetTextAngle(90);
        t->SetTextSize(0.018);
        t->SetTextFont(72);
        for(int aa = 0; aa < EnergiesVec.size(); aa++){
            t->DrawText(EnergiesVec[aa],0,const_cast<char*>(XLabelsVec[aa].c_str()));
            //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
        }
*/
        for(int aa = 0; aa < XLabelsIDVec.size(); aa++){
            int i = XLabelsIDVec[aa];
            mg->GetXaxis()->ChangeLabel(i,-1,-1,-1,-1,-1,XLabelsVec[aa].c_str());;
        }
    }
    
    mg->Draw("ALP");
    
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);
    
    //leg->SetX1(X1LegPos);
    //leg->SetX2(X2LegPos);
    //leg->SetY1(Y1LegPos);
    //leg->SetY2(Y2LegPos);
    
    leg->Draw();
    
    ResCanvas->Print(FileName.c_str());
    delete ResCanvas;
    
}

void G4DoseCalcsAnalysis::createDataInCSVFormat(){
    
    std::string Quantity;
    std::string Geometry;
    std::string PARTICLE_NAME;
    std::string Source_ORG;
    std::string Target_ORG;
    
    // file for each quantity (a table for each geometry and particle and source and SD)
    
    for ( auto Obeg = ResultTable.begin(); Obeg != ResultTable.end(); ++Obeg  )
    {
        Quantity = Obeg->first;
        
        std::string FileName = GraphsDirectoryPath +"DoseCalcs_Particles_"+Quantity+".csv";
        std::ofstream outfile(FileName , std::ios::app);
        if(outfile.is_open()){
            std::cout << "\nCreating file " << FileName << std::endl ;
        }
        
        
        std::ostringstream LatexText;
        LatexText << QuantityUnit[Quantity] << "\n";
        
        for ( auto Mbeg = Obeg->second.begin(); Mbeg != Obeg->second.end(); ++Mbeg  )
        {
            Geometry = Mbeg->first;
            
            
            for ( auto Abeg = Mbeg->second.begin(); Abeg != Mbeg->second.end(); ++Abeg  )
            {
                PARTICLE_NAME = Abeg->first;
                
                for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  )
                {
                    Source_ORG = Bbeg->first;
                    
                    LatexText << "\n" << Geometry << "," << PARTICLE_NAME << ","<< Source_ORG << ","<<  "\n";
                    
                    LatexText << "Target/Energies-SDev" << ",";
                    for (G4int gg = 0 ; gg < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; gg++) {
                        double val = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg];
                        LatexText << val  << "," << " ,";
                    }
                    LatexText << "\n";
                    
                    for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg  )
                    {
                        
                        Target_ORG = Cbeg->first;
                        LatexText << Target_ORG  << ",";
                        for (G4int gg = 0 ; gg < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; gg++) {
                            double val = ResultTable[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg]];
                            //double error = RelativeStandartDeviationPerCent[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg]];
                            double error = StandartDeviation[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg]];
                            if(val == 0 || __isinf(val) || __isnan(val)){
                                LatexText << " ," << " ,";
                            }else{
                                LatexText << val << "," << error << ",";
                            }
                        }
                        LatexText << "\n";
                        
                    }
                    
                    LatexText << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
                    
                }
            }
            
            LatexText << "\n\n*****************************************************************************************************************************************************************************************************************************************************************\n";
            
        }
        if(outfile.is_open()){
            outfile << LatexText.str();
            outfile.close();
        }
    }
    
    // file for each quantity (a table for each geometry and particle and source and SD) with reference
    
    for ( auto Obeg = ResultTable.begin(); Obeg != ResultTable.end(); ++Obeg  )
    {
        Quantity = Obeg->first;
        
        if(ReferenceTable[Quantity].size() == 0){continue;}
        
        std::string FileName = GraphsDirectoryPath + "DoseCalcs-" + CompareReferenceName + "_Particles_"+Quantity+".csv";
        std::ofstream outfile(FileName , std::ios::app);
        if(outfile.is_open()){
            std::cout << "\nCreating file " << FileName << std::endl ;
        }
        
        
        std::ostringstream LatexText;
        LatexText << QuantityUnit[Quantity] << "\n";
        
        for ( auto Mbeg = Obeg->second.begin(); Mbeg != Obeg->second.end(); ++Mbeg  )
        {
            Geometry = Mbeg->first;
            
            if(ReferenceTable[Quantity][Geometry].size() == 0){continue;}
            
            for ( auto Abeg = Mbeg->second.begin(); Abeg != Mbeg->second.end(); ++Abeg  )
            {
                PARTICLE_NAME = Abeg->first;
                
                if(ReferenceTable[Quantity][Geometry][PARTICLE_NAME].size() == 0){continue;}
                
                for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  )
                {
                    Source_ORG = Bbeg->first;
                    
                    if(ReferenceTable[Quantity][Geometry][PARTICLE_NAME][Source_ORG].size() == 0){continue;}
                    
                    LatexText << "\n" << Geometry << "," << PARTICLE_NAME << ","<< Source_ORG << ","<<  "\n";
                    
                    LatexText << "Target/Energies(DoseCalcs-"<< CompareReferenceName <<"-Ratio)" << ",";
                    for (G4int gg = 0 ; gg < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; gg++) {
                        double val = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg];
                        LatexText << val  << "," << " ,"<< " Ratio,";
                    }
                    LatexText << "\n";
                    
                    for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg  )
                    {
                        
                        Target_ORG = Cbeg->first;
                        LatexText << Target_ORG  << ",";
                        for (G4int gg = 0 ; gg < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; gg++) {
                            double val = ResultTable[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg]];
                            double val2 = ReferenceTable[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][gg]];
                            if(val == MinValForLog || val == 0 || __isinf(val) || __isnan(val)){
                                
                                LatexText << " ,";
                                if(val2 == MinValForLog || val2 == 0 || __isinf(val2) || __isnan(val2))
                                {
                                    LatexText << " ," << " ,";
                                    
                                }else{
                                    double ratio = val/val2;
                                    LatexText << val2 << " ,"  << ratio << " ,";
                                }
                                
                            }else{
                                
                                LatexText << val << ",";
                                if(val2 == MinValForLog || val2 == 0 || __isinf(val2) || __isnan(val2))
                                {
                                    LatexText << " ," << " ,";
                                    
                                }else{
                                    double ratio = val/val2;
                                    LatexText << val2 << " ,"  << ratio << " ,";
                                }
                            }
                        }
                        LatexText << "\n";
                        
                    }
                    
                    LatexText << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
                    
                }
            }
            
            LatexText << "\n\n*****************************************************************************************************************************************************************************************************************************************************************\n";
            
        }
        if(outfile.is_open()){
            outfile << LatexText.str();
            outfile.close();
        }
    }
    
    // file for each quantity (a table for each geometry and radiotracer and source and SD)
    
    for ( auto Obeg = ResultQuantityGeometryRadioTracerSourceTargetValues.begin(); Obeg != ResultQuantityGeometryRadioTracerSourceTargetValues.end(); ++Obeg  )
    {
        Quantity = Obeg->first;
        
        std::string FileName = GraphsDirectoryPath +"DoseCalcs_RadioNuclides_"+Quantity+".csv";
        std::ofstream outfile(FileName , std::ios::app);
        if(outfile.is_open()){
            std::cout << "\nCreating file " << FileName << std::endl ;
        }
        
        
        std::ostringstream LatexText;
        LatexText << QuantityUnit[Quantity] << "\n";
        
        for ( auto Mbeg = Obeg->second.begin(); Mbeg != Obeg->second.end(); ++Mbeg  )
        {
            Geometry = Mbeg->first;
            
            for ( auto Abeg = Mbeg->second.begin(); Abeg != Mbeg->second.end(); ++Abeg  )
            {
                PARTICLE_NAME = Abeg->first;
                
                for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  )
                {
                    Source_ORG = Bbeg->first;
                    
                    LatexText << "\n" << Geometry << "," << PARTICLE_NAME << ","<< Source_ORG << ","<<  "\n";
                    
                    LatexText << "Target" << "," << QuantityUnit[Quantity] << ", Rel_SD%" << "\n";
                    
                    for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg  )
                    {
                        
                        Target_ORG = Cbeg->first;
                        LatexText << Target_ORG  << ",";
                        //LatexText << Target_ORG  << "<-" << Source_ORG << ",";

                        double val = ResultQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG];
                        double error = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG];
                        if(val == 0 || __isinf(val) || __isnan(val)){
                            LatexText << " ," << " ,";
                        }else{
                            LatexText << val << "," << error << ",";
                        }


                        //if(val == 0 || __isinf(val) || __isnan(val)){
                        //    LatexText << "," ;
                        //}else{
                        //    LatexText << val << "," ;
                        //}



                        
                        LatexText << "\n";
                        
                    }
                    
                    LatexText << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
                    
                }
            }
            
            LatexText << "\n\n*****************************************************************************************************************************************************************************************************************************************************************\n";
            
        }
        if(outfile.is_open()){
            outfile << LatexText.str();
            outfile.close();
        }
    }
    
    // file for each quantity (a table for each radiotracer (all geometries and combinations with SD)

    for ( auto Obeg = ResultQuantityGeometryRadioTracerSourceTargetValues.begin(); Obeg != ResultQuantityGeometryRadioTracerSourceTargetValues.end(); ++Obeg  )
    {
        Quantity = Obeg->first;

        std::string FileName = GraphsDirectoryPath +"DoseCalcs_RadioNuclides_"+Quantity+".csv";
        std::ofstream outfile(FileName , std::ios::app);
        if(outfile.is_open()){
            std::cout << "\nCreating file " << FileName << std::endl ;
        }


        std::ostringstream LatexText;
        LatexText << QuantityUnit[Quantity] << "\n";

        for ( auto Mbeg = Obeg->second.begin(); Mbeg != Obeg->second.end(); ++Mbeg  )
        {
            Geometry = Mbeg->first;

            for ( auto Abeg = Mbeg->second.begin(); Abeg != Mbeg->second.end(); ++Abeg  )
            {
                PARTICLE_NAME = Abeg->first;

                LatexText << "\n" << Geometry << "," << PARTICLE_NAME <<  "\n";

                LatexText << "Target<-Source" << "," << QuantityUnit[Quantity] << ", Rel_SD%" << "\n";

                for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  )
                {
                    Source_ORG = Bbeg->first;
                    if(Source_ORG == "IntakeIntoBody"){
                        continue;
                    }

                    for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg  )
                    {

                        Target_ORG = Cbeg->first;
                        //LatexText << Target_ORG  << ",";
                        LatexText << Target_ORG  << "<-" << Source_ORG << ",";

                        double val = ResultQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG];
                        double error = QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG];
                        if(val == 0 || __isinf(val) || __isnan(val)){
                            LatexText << " ," << " ,";
                        }else{
                            LatexText << val << "," << error << ",";
                        }


                        //if(val == 0 || __isinf(val) || __isnan(val)){
                        //    LatexText << "," ;
                        //}else{
                        //    LatexText << val << "," ;
                        //}




                        LatexText << "\n";

                    }

                    //LatexText << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";

                }
            }

            LatexText << "\n\n*****************************************************************************************************************************************************************************************************************************************************************\n";

        }
        if(outfile.is_open()){
            outfile << LatexText.str();
            outfile.close();
        }
    }

    // file for each quantity (a table for each geometry and radiotracer and source and SD) with reference
    
    for ( auto Obeg = ResultQuantityGeometryRadioTracerSourceTargetValues.begin(); Obeg != ResultQuantityGeometryRadioTracerSourceTargetValues.end(); ++Obeg  )
    {
        Quantity = Obeg->first;
        
        if(ReferenceQuantityGeometryRadioTracerSourceTargetValues[Quantity].size() == 0){continue;}
        
        std::string FileName = GraphsDirectoryPath + "DoseCalcs-" + CompareReferenceName + "_Radionuclide_"+Quantity+".csv";
        std::ofstream outfile(FileName , std::ios::app);
        if(outfile.is_open()){
            std::cout << "\nCreating file " << FileName << std::endl ;
        }
        
        std::ostringstream LatexText;
        LatexText << QuantityUnit[Quantity] << "\n";
        
        for ( auto Mbeg = Obeg->second.begin(); Mbeg != Obeg->second.end(); ++Mbeg  )
        {
            Geometry = Mbeg->first;
            
            if(ReferenceQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry].size() == 0){continue;}
            
            for ( auto Abeg = Mbeg->second.begin(); Abeg != Mbeg->second.end(); ++Abeg  )
            {
                PARTICLE_NAME = Abeg->first;
                
                if(ReferenceQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][PARTICLE_NAME].size() == 0){continue;}
                
                for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  )
                {
                    Source_ORG = Bbeg->first;
                    
                    if(ReferenceQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][PARTICLE_NAME][Source_ORG].size() == 0){continue;}
                    
                    LatexText << "\n" << Geometry << "," << PARTICLE_NAME << ","<< Source_ORG << ","<<  "\n";
                    
                    LatexText << "Target" << "," << QuantityUnit[Quantity] << "-DoseCalcs,"<< QuantityUnit[Quantity] <<"-"<< CompareReferenceName <<", Ratio" << "\n";
                    
                    for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg  )
                    {
                        
                        Target_ORG = Cbeg->first;
                        LatexText << Target_ORG  << ",";
                        //LatexText << Target_ORG  << "<-" << Source_ORG << ",";

                        double val = ResultQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG];
                        double val2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[Quantity][Geometry][PARTICLE_NAME][Source_ORG][Target_ORG];
                        if(val == MinValForLog || val == 0 || __isinf(val) || __isnan(val)){
                            
                            LatexText << " ,";
                            if(val2 == MinValForLog || val2 == 0 || __isinf(val2) || __isnan(val2))
                            {
                                LatexText << " ," << " ,";
                                
                            }else{
                                double ratio = val/val2;
                                LatexText << val2 << " ,"  << ratio << " ,";
                            }
                            
                        }else{
                            
                            LatexText << val << ",";
                            if(val2 == MinValForLog || val2 == 0 || __isinf(val2) || __isnan(val2))
                            {
                                LatexText << " ," << " ,";
                                
                            }else{
                                double ratio = val/val2;
                                LatexText << val2 << " ,"  << ratio << " ,";
                            }
                        }
                        LatexText << "\n";
                    }
                    
                    LatexText << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
                    
                }
            }
            
            LatexText << "\n\n*****************************************************************************************************************************************************************************************************************************************************************\n";
            
        }
        if(outfile.is_open()){
            outfile << LatexText.str();
            outfile.close();
        }
    }
    
}

#ifdef ANALYSIS_USE

// callResultParticleSourceEnergyTimeed from main() in graph class
void G4DoseCalcsAnalysis::GenerateSelfCrossGraphs(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
# ifndef __CINT__
    
    //std::cout << " No Results data for " << QuantitiesToScore << "\n" << std::endl;
    
    ReadResultFile();

    for (int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {

        QuantitiesToScore = QuantityNamesToScore[gg];

        std::cout <<"\nFor quantity: " << QuantityUnit[QuantitiesToScore] << std::endl;

        DataInitialization();

        if(GenerateRelativeSDevGraph == "yes"){
            GenerateStandardDeviationGraphs();
        }

        GenerateResultInOneGraph();
        GenerateResultForEachParticleEnergyInOneGraph();
        GenerateSourceEnegyGraph();


        if(GenerateRelativeErrGraph == "yes" && CompareReferenceNames.size() != 0 && RefFilePaths.size() != 0 && ReferenceTable.size() != 0){
            GenerateComparisonFactorGraphs();
            GenerateResultReferenceInOneGraph();
        }else{
            std::cout << "\n\n                                                          ========= Canno't generate graphs for comparison with reference ========= "<< "\n" << std::endl;
            std::cout << "\n\nCheck the name and data file path of reference, or set /AnalysisData/generateRelativeErrGraph ... command" << std::endl;

            std::cout << "\n\nCheck the name and data file path of reference, or set /AnalysisData/generateRelativeErrGraph ... command" << std::endl;
        }

        if(GenerateResultsForRadioTracer == true){
            GenerateRadioTracerResultsInOneGraphForAllComAndSrcReg();
            if(GenerateRelativeErrGraph == "yes" && CompareReferenceNames.size() != 0 && RefFilePaths.size() != 0 && ReferenceQuantityGeometryRadioTracerSourceTargetValues.size() != 0){
                GenerateRadioTracerResultsInOneGraphWithComparisonForAllComAndSrcReg();
                GenerateRadioTracerComparisonFactorGraphsForAllComAndSrcReg();
            }else{
                std::cout << "Check the name and data file path of reference, or set /AnalysisData/generateRelativeErrGraph ... command" << std::endl;
            }
        }
    }
# endif

}
void G4DoseCalcsAnalysis::GenerateVoxel2DGraphs(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    std::ostringstream filename1;
    filename1 << ResultDirectoryPath << "/VoxelsResults" ;
    std::ifstream file(filename1.str().c_str() , std::ios::binary);
    if(!file.is_open()){
        std::cout << "The voxelized results will not be analyzed, cannot open the file " << filename1.str().c_str() << std::endl ;
        return;
    }

    ReadVoxelsResult();

    if( SliceFor2DGraph != "none" ){
        for (int gg = 0 ; gg < QuantityNamesToScore.size() ; gg++) {
            
            QuantitiesToScore = QuantityNamesToScore[gg];
            DataInitialization();
            std::cout <<"\nFor the quantity: " << QuantityUnit[QuantitiesToScore] << std::endl;
            //CreateHeatMap();
        }
    }
    else{
        std::cout <<" SliceFor2DGraph = none !!!!!!!!!!!!!!"<< std::endl;
    }

    CreatePercentageDepthDoseGraph();
    CreateDoseProfile();
    
}
void G4DoseCalcsAnalysis::GenerateEventsDataHisto(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    if(EventsDataHistograms == "yes"){
        
    }
    else {
        return;
    }
    std::ifstream PosFileStream, EneFileStream, MomDirFileStream;
    PosFileStream.open(PositionDataFile , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
    EneFileStream.open(EnergyDataFile , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
    MomDirFileStream.open(MomDirDataFile , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
    
    if(PosFileStream.is_open() ){}
    else {
        std::cout <<"Canno't open data file " << appBuildDir+"/"+PositionDataFile << std::endl;
        return;
    }
    
    if(EneFileStream.is_open() ){}
    else {
        std::cout <<"Canno't open data file " << appBuildDir+"/"+EnergyDataFile << std::endl;
        return;
    }
    
    if(MomDirFileStream.is_open() ){}
    else {
        std::cout <<"Canno't open data file " << appBuildDir+"/"+MomDirDataFile << std::endl;
        return;
    }
    
    double valX, valY, valZ, valE, valMDX, valMDY, valMDZ;
    //float valE;
    std::vector<double> E;
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Z;
    std::vector<double> MDX;
    std::vector<double> MDY;
    std::vector<double> MDZ;
    
    int EventIncre = 1;
    double MinE, MaxE,MinX,MaxX,MinY,MaxY, MinZ,MaxZ, MinMDX,MaxMDX,MinMDY,MaxMDY, MinMDZ,MaxMDZ;
    
    while(PosFileStream.peek() != EOF){ // read just what we need to simulate
        
        EneFileStream >> valE ;
        //std::cout << EventIncre << " valE " << valE << std::endl;
        
        PosFileStream  >> valX >> valY >> valZ ;
        //std::cout << EventIncre << " valX " << valX << std::endl;
        
        MomDirFileStream  >> valMDX >> valMDY >> valMDZ ;
        //std::cout << EventIncre << " valMDX  " << valMDX << std::endl;
        
        if(EventIncre == 1){
            MinE = valE ;
            MaxE = valE ;
            MinX = valX ;
            MaxX = valX ;
            MinY = valY ;
            MaxY = valY ;
            MinZ = valZ ;
            MaxZ = valZ ;
            MinMDX = valMDX ;
            MaxMDX = valMDX ;
            MinMDY = valMDY ;
            MaxMDY = valMDY ;
            MinMDZ = valMDZ ;
            MaxMDZ = valMDZ ;
            //std::cout << " MinZ " << MinZ << std::endl;
        }
        
        if(valE < MinE){
            MinE = valE ;
        }
        if(valE > MaxE){
            MaxE = valE ;
        }
        E.push_back(valE);
        
        if(valX < MinX){
            MinX = valX ;
        }
        if(valX > MaxX){
            MaxX = valX ;
        }
        X.push_back(valX);
        
        if(valY < MinY){
            MinY = valY ;
        }
        if(valY > MaxY){
            MaxY = valY ;
        }
        Y.push_back(valY);
        
        if(valZ < MinZ){
            MinZ = valZ ;
        }
        if(valZ > MaxZ){
            MaxZ = valZ ;
        }
        Z.push_back(valZ);
        
        if(valMDX < MinMDX){
            MinMDX = valMDX ;
        }
        if(valMDX > MaxMDX){
            MaxMDX = valMDX ;
        }
        MDX.push_back(valMDX);
        
        if(valMDY < MinMDY){
            MinMDY = valMDY ;
        }
        if(valMDY > MaxMDY){
            MaxMDY = valMDY ;
        }
        MDY.push_back(valMDY);
        
        if(valMDZ < MinMDZ){
            MinMDZ = valMDZ ;
        }
        if(valMDZ > MaxMDZ){
            MaxMDZ = valMDZ ;
        }
        MDZ.push_back(valMDZ);
        
        EventIncre++;
    }
    
    EneFileStream.close();
    PosFileStream.close();
    MomDirFileStream.close();
    
    std::cout << " End Of Reading Data File " << appBuildDir << "/" << PositionDataFile << std::endl;
    
    // ******************************************************************************************* Energy-Spatial Histograms
    
    TCanvas * ESDHC = new TCanvas("Energy Spatial Data Histograms", "EnergySpatialDataHistograms", 600,500);
    ESDHC->Divide(2,2); // 2 columns, each is divided on 3
    
    //std::cout << "MinE: "<< MinE << " MaxE: " << MaxE << std::endl;
    
    double MaxENew, MinENew;
    if(EnergyDistribution == "Rayleigh"){
        MaxENew = RayleighEmax*2;
    }
    
    if(EnergyDistribution == "Mono"){
        MinENew = 0.5*MonoEnergy; if(MinE < 0){ MinE = 0; }
        MaxENew = 1.5*MonoEnergy;
    }
    
    std::string HistoTitle = "Events Energy Spectrum; E(MeV); N(E) ";
    
    TH1F * EnergyHisto = new TH1F(EnergyDistribution.c_str() , HistoTitle.c_str(), 500 , MinENew, MaxENew);
    //TH1F * EnergyHisto = new TH1F();
    //EnergyHisto->SetName(EnergyDistribution.c_str());
    //EnergyHisto->SetTitle(HistoTitle.c_str());
    
    //Energies = new double[E.size()];
    for (unsigned int kl = 0 ; kl <E.size() ;kl++) {
        EnergyHisto->Fill(E[kl]);
    }
    
    std::cout << " Creating Energy Histogram " << std::endl;
    ESDHC->cd(1); EnergyHisto->Draw(); //ALP
    if(UseGridXY=="yes"){
        ESDHC->SetGridx(); ESDHC->SetGridy();
    }
    MinENew = MinE - 0.5*MinE; if(MinENew < 0){ MinENew = 0; }
    MaxENew = MaxE + 0.5*MaxE;
    
    //std::cout << "MinX: "<< MinX << " MaxX: "<< MaxX << " MinY: " << MinY  << " MaxY: " << MaxY  << " MinZ: " << MinZ << " MaxZ: " << MaxZ << std::endl;
    
    unsigned int NXExBins = EventIncre/5000;
    
    MinX = MinX - (MaxX - MinX)/2; MaxX = MaxX + (MaxX - MinX)/2;
    MinY = MinY - (MaxY - MinY)/2; MaxY = MaxY + (MaxY - MinY)/2;
    MinZ = MinZ - (MaxZ - MinZ)/2; MaxZ = MaxZ + (MaxZ - MinZ)/2;
    
    std::stringstream NXEx, NYEy, NZEz ;
    NXEx <<"2D distribution for X-E"<<";X(mm)"<< ";E(Mev)"  ;
    NYEy <<"2D distribution for Y-E"<<";Y(mm)"<< ";E(Mev)"  ;
    NZEz <<"2D distribution for Z-E"<<";Z(mm)"<< ";E(Mev)"  ;
    
    TH2F *DistrHistoXEx = new TH2F("X-E",NXEx.str().c_str() ,NXExBins,MinX,MaxX,NXExBins,MinENew,MaxENew);
    TH2F *DistrHistoYEy = new TH2F("Y-E",NYEy.str().c_str() ,NXExBins,MinY,MaxY,NXExBins,MinENew,MaxENew);
    TH2F *DistrHistoZEz = new TH2F("Z-E",NZEz.str().c_str() ,NXExBins,MinZ,MaxZ,NXExBins,MinENew,MaxENew);
    
    for (unsigned int kl = 0 ; kl <X.size() ;kl++) {
        //std::cout << X[kl] << " " << Y[kl] << " "<< Z[kl] << " "<< E[kl] << std::endl;
        DistrHistoXEx->Fill(X[kl], E[kl]);
        DistrHistoYEy->Fill(Y[kl], E[kl]);
        DistrHistoZEz->Fill(Z[kl], E[kl]);
    }
    
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);
    
    ESDHC->cd(2); DistrHistoXEx->Draw(); //ALP
    ESDHC->cd(3); DistrHistoYEy->Draw(); //ALP
    ESDHC->cd(4); DistrHistoZEz->Draw(); //ALP
    
    std::stringstream file1 ; file1 << GraphsDirectoryPath << "Energy_Spatial_Histo"<<GraphsExt ;
    std::cout << " Creating Energy-Spatial Histograms and Saving to " << file1.str().c_str() << " ...." << std::endl;
    ESDHC->Print(file1.str().c_str());
    
    //delete DistrHistoXEx; delete DistrHistoYEy; delete DistrHistoZEz; delete EnergyHisto;
    delete ESDHC;
    
    // ******************************************************************************************* Spatial 2D and 3D histograms
    
    TCanvas * SDHC = new TCanvas("Spatial Data Histograms", "SpatialDataHistograms", 600,500);
    SDHC->Divide(2,2); // 2 columns, each is divided on 3
    
    //std::cout << "MinX: "<< MinX << " MaxX: "<< MaxX << " MinY: " << MinY  << " MaxY: " << MaxY  << " MinZ: " << MinZ << " MaxZ: " << MaxZ << std::endl;
    
    unsigned int NBins = EventIncre/5000;
    
    std::stringstream NXY, NXZ, NYZ , NXYZ;
    NXY <<"2D Spatial distribution for XY plane"<<";X(mm)"<< ";Y(mm)"  ;
    NXZ <<"2D Spatial distribution for XZ plane"<<";X(mm)"<< ";Z(mm)"  ;
    NYZ <<"2D Spatial distribution for YZ plane"<<";Y(mm)"<< ";Z(mm)"  ;
    NXYZ <<"3D Spatial distribution for XYZ "<<";X(mm)"<< ";Y(mm)"<< ";Z(mm)"  ;
    
    TH2F *SpatialDistrHistoXY = new TH2F("Spatial XY",NXY.str().c_str() ,NBins,MinX,MaxX,NBins,MinY,MaxY);
    TH2F *SpatialDistrHistoXZ = new TH2F("Spatial XZ",NXZ.str().c_str() ,NBins,MinX,MaxX,NBins,MinZ,MaxZ);
    TH2F *SpatialDistrHistoYZ = new TH2F("Spatial YZ",NYZ.str().c_str() ,NBins,MinY,MaxY,NBins,MinZ,MaxZ);
    
    TH3F *SpatialDistrHistoXYZ = new TH3F("Spatial XYZ",NXYZ.str().c_str() , NBins,MinX,MaxX, NBins,MinY,MaxY, NBins,MinZ,MaxZ);
    
    for (unsigned int kl = 0 ; kl <X.size() ;kl++) {
        
        //std::cout << X[kl] << " " << Y[kl] << " "<< Z[kl] << std::endl;
        
        SpatialDistrHistoXY->Fill(X[kl], Y[kl]);
        SpatialDistrHistoXZ->Fill(X[kl], Z[kl]);
        SpatialDistrHistoYZ->Fill(Y[kl], Z[kl]);
        SpatialDistrHistoXYZ->Fill(X[kl], Y[kl], Z[kl]);
        
    }
    
    SDHC->cd(1); SpatialDistrHistoXY->Draw(); //ALP
    SDHC->cd(2); SpatialDistrHistoXZ->Draw(); //ALP
    SDHC->cd(3); SpatialDistrHistoYZ->Draw(); //ALP
    SDHC->cd(4); SpatialDistrHistoXYZ->Draw(); //ALP
    
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);
    if(UseGridXY=="yes"){
        SDHC->SetGridx(); SDHC->SetGridy();
    }
    std::stringstream file2 ; file2 << GraphsDirectoryPath << "Spatial_Histo"<<GraphsExt ;
    std::cout << " Creating Spatial Histograms and Saving to " << file2.str().c_str() << " ...." << std::endl;
    SDHC->Print(file2.str().c_str());
    
    //delete SpatialDistrHistoXY; delete SpatialDistrHistoYZ; delete SpatialDistrHistoXZ ; delete SpatialDistrHistoXYZ;
    delete SDHC;
    
    // ******************************************************************************************* Angular 2D and 3D histograms
    
    TCanvas * ADHC = new TCanvas("Angular Data Histograms", "AngularDataHistograms", 600,500);     //TCanvas * c2 = new TCanvas("c2", "c2"); // should be given befor Draw()
    ADHC->Divide(2,2); // 2 columns, each is divided on 3
    
    MinMDX = MinMDX - (MaxMDX - MinMDX)/2; MaxMDX = MaxMDX + (MaxMDX - MinMDX)/2;
    MinMDY = MinMDY - (MaxMDY - MinMDY)/2; MaxMDY = MaxMDY + (MaxMDY - MinMDY)/2;
    MinMDZ = MinMDZ - (MaxMDZ - MinMDZ)/2; MaxMDZ = MaxMDZ + (MaxMDZ - MinMDZ)/2;
    
    //std::cout << "MinX: "<< MinX << " MaxX: "<< MaxX << " MinY: " << MinY  << " MaxY: " << MaxY  << " MinZ: " << MinZ << " MaxZ: " << MaxZ << std::endl;
    
    unsigned int MDNBins = EventIncre/5000;
    
    std::stringstream MDNX, MDNY, MDNZ , MDNXYZ;
    MDNX <<"Momentum Direction distribution for ux"<<";ux"<< ";N(ux)"  ;
    MDNY <<"Momentum Direction distribution for uy"<<";uy"<< ";N(uy)"  ;
    MDNZ <<"Momentum Direction distribution for uz"<<";uz"<< ";N(uz)"  ;
    MDNXYZ <<"Momentum Direction distribution for ux, uy, uz "<<";ux"<< ";uy"<< ";uz"  ;
    
    std::stringstream ux, uy, uz;
    ux << "ux in " << MomDirDistribution  ;
    uy << "uy in " << MomDirDistribution  ;
    uz << "uz in " << MomDirDistribution  ;
    
    TH1F *AngularDistrHistoX = new TH1F(ux.str().c_str() ,MDNX.str().c_str() ,MDNBins,MinMDX,MaxMDX);
    TH1F *AngularDistrHistoY = new TH1F(uy.str().c_str() ,MDNY.str().c_str() ,MDNBins,MinMDX,MaxMDX);
    TH1F *AngularDistrHistoZ = new TH1F(uz.str().c_str() ,MDNZ.str().c_str() ,MDNBins,MinMDY,MaxMDY);
    
    TH3F *AngularDistrHistoXYZ = new TH3F(MomDirDistribution.c_str() ,MDNXYZ.str().c_str() , MDNBins,MinMDX,MaxMDX, MDNBins,MinMDY,MaxMDY, MDNBins,MinMDZ,MaxMDZ);
    
    for (unsigned int kl = 0 ; kl <MDX.size() ;kl++) {
        //std::cout << MDX[kl] << " " << MDY[kl] << " "<< MDZ[kl] << std::endl;
        
        AngularDistrHistoX->Fill(MDX[kl]);
        AngularDistrHistoY->Fill(MDY[kl]);
        AngularDistrHistoZ->Fill(MDZ[kl]);
        AngularDistrHistoXYZ->Fill(MDX[kl], MDY[kl], MDZ[kl]);
    }
    
    ADHC->cd(1); AngularDistrHistoX->Draw(); //ALP
    ADHC->cd(2); AngularDistrHistoY->Draw(); //ALP
    ADHC->cd(3); AngularDistrHistoZ->Draw(); //ALP
    ADHC->cd(4); AngularDistrHistoXYZ->Draw(); //ALP
    
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);
    if(UseGridXY=="yes"){
        ADHC->SetGridx(); ADHC->SetGridy();
    }
    std::stringstream file3 ; file3 << GraphsDirectoryPath << "Momentum_Direction_Histo"<<GraphsExt ;
    std::cout << " Creating Momentum Direction Histograms and Saving to " << file3.str().c_str() << " ...." << std::endl;
    ADHC->Print(file3.str().c_str());
    
    //delete AngularDistrHistoX; delete AngularDistrHistoY; delete AngularDistrHistoZ ; delete AngularDistrHistoXYZ;
    delete ADHC;
    
}
void G4DoseCalcsAnalysis::TestGraphGeneration(){
    
    std::string PARTICLE_NAME = "gamma";
    std::string Source_ORG = "Thyroid";
    std::string Target_ORG ="Thyroid";
    QuantitiesToScore = "SAF";
    
    ReadResultFile();
    DataInitialization();
    //ReferenceTable =  ReadReferenceFile(RefFilePath);
    //if(RefFilePaths.size() > 1){
    //  ReferenceTable2 = ReadReferenceFile(RefFilePaths[1]);
    //}
    
    NumOfEne = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size();
    
    double xi[NumOfEne], yi[NumOfEne];
    double ai[NumOfEne], ei[NumOfEne], errInERef[NumOfRefEne], errRef[NumOfRefEne] ;
    
    std::cout <<  "RefFilePaths.size()" << " -> "<< RefFilePaths.size() << std::endl ;
    std::cout <<  "NumOfEne" << " -> "<< NumOfEne << std::endl ;
    std::cout <<  "QuantitiesToScore" << " -> "<< QuantitiesToScore << std::endl ;
    std::cout <<  "QuantityUseLog[QuantitiesToScore]" << " -> "<< QuantityUseLog[QuantitiesToScore] << std::endl ;
    std::cout <<  "PARTICLE_NAME" << " -> "<< PARTICLE_NAME << std::endl ;
    std::cout <<  "Source_ORG" << " -> "<< Source_ORG << std::endl ;
    std::cout <<  "Target_ORG" << " -> "<< Target_ORG << std::endl ;
    std::cout <<  "QuantityUnit[QuantitiesToScore]" << " -> "<< QuantityUnit[QuantitiesToScore] << std::endl ;
    
    std::string multigraphTitle = Target_ORG+"<-"+Source_ORG + " for " +PARTICLE_NAME+"; Energy(MeV); "+ QuantityUnit[QuantitiesToScore];
    TCanvas * Res_Ref_Canvas = new TCanvas("Test", "Test");     //TCanvas * c2 = new TCanvas("c2", "c2"); // should be given befor Draw()
    std::string FileName = GraphsDirectoryPath +"Test_"+QuantitiesToScore+"_"+GeometrySymbol+"_"+PARTICLE_NAME+"_"+Source_ORG +"_"+Target_ORG+GraphsExt;
    
    for(int ii = 0; ii < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size(); ii++){
        xi[ii] = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii];
        yi[ii] = ResultTable[QuantitiesToScore][GeometrySymbol][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]];
        ai[ii] = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii];
        ei[ii] = ReferenceTable[QuantitiesToScore][GeometrySymbol][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]];

        if(QuantityUseLog[QuantitiesToScore] == true && ei[ii] < MinValForLog){
            ei[ii] = MinValForLog;
        }
        if(QuantityUseLog[QuantitiesToScore] == true && yi[ii] < MinValForLog){
            yi[ii] = MinValForLog;
        }
        std::cout << " xi[ii] " <<  xi[ii] ;
        std::cout << " yi[ii] " <<  yi[ii] ;
        std::cout << " ai[ii] " <<  ai[ii] ;
        std::cout << " ei[ii] " <<  ei[ii] ;
        std::cout << " MinValForLog "<< MinValForLog <<std::endl ;
    }
    
    TMultiGraph *mg = new TMultiGraph();
    
    TGraph* gr1;
    if(AddErrorBarInGraphs=="yes"){
        
        double zi[ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size()], mi[ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size()];
        for(int ii = 0; ii < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size(); ii++){
            zi[ii] = 0;
            mi[ii] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]];
            std::cout << "error X " <<  zi[ii] << " error Y " << mi[ii] <<std::endl ;
        }
        
        gr1 = CreateGraphErrors (NumOfEne , xi, yi, zi, mi); //this graph is related if we want to remove the zeros of SAF from Graph, if not we have to uncomment the graph1 below
    }
    else{
        gr1 = CreateGraph (NumOfEne , xi, yi); //this graph is related if we want to remove the zeros of SAF from Graph, if not we have to uncomment the graph1 below
    }
    
    gr1->SetName("DoseCalcs");
    gr1->SetTitle(gr1->GetName());
    gr1->SetLineWidth(1);
    gr1->SetMarkerStyle(50+1);
    gr1->SetMarkerColor(1+1);
    gr1->SetLineColor(1+1);
    mg->Add(gr1);
    
    TGraph* gr2 = CreateGraph (NumOfEne , ai, ei ); //this graph is related if we want to remove the zeros of SAF from Graph, if not we have to uncomment the graph1 below
    gr2->SetName(CompareReferenceName.c_str());
    gr2->SetTitle(gr2->GetName());
    gr2->SetLineWidth(1);
    gr2->SetMarkerStyle(50+2);
    gr2->SetMarkerColor(1+2);
    gr2->SetLineColor(1+2);
    mg->Add(gr2);
    
    TGraph* gr3;
    if(RefFilePaths.size() > 1){
        
        double zi[ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size()], mi[ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size()];
        for(int ii = 0; ii < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size(); ii++){
            zi[ii] = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii];
            mi[ii] = ReferenceTable2[QuantitiesToScore][GeometrySymbol][PARTICLE_NAME][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]];
            //std::cout <<  zi[ii] << " --- " << mi[ii] <<std::endl ;
        }
        
        //double* mi = AccumulateThirdRefDataToGraphs(PARTICLE_NAME, Source_ORG, Target_ORG);
        gr3 = CreateGraph (ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() , zi, mi ); //this graph is related if we want to remove the zeros of SAF from Graph, if not we have to uncomment the graph1 below
        gr3->SetName(CompareReferenceNames[1].c_str());
        gr3->SetTitle(gr3->GetName());
        gr3->SetLineWidth(1);
        gr3->SetMarkerStyle(50+3);
        gr3->SetMarkerColor(1+3);
        gr3->SetLineColor(1+3);
        mg->Add(gr3);
    }
    
    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }
    if(UseLogE == "yes"){
        gPad->SetLogx(1);
    }
    
    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }
    
    if(PrintTitle == "yes"){
        mg->SetTitle(multigraphTitle.c_str());
    }
    else{
        multigraphTitle = "; Energy(MeV); "+ QuantityUnit[QuantitiesToScore];
        mg->SetTitle(multigraphTitle.c_str());
    }
    
    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();
    
    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);
    
    
    mg->Draw("ALP");
    //gPad->Update();
    
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);
    
    TLegend *leg = new TLegend();
    //leg->AddEntry(gr1,gr1->GetName(),"LP");  // to add the explanation of this colored line
    leg->AddEntry(gr2,gr2->GetName(),"LP");  // to add the explanation of this colored line
    if(RefFilePaths.size() > 1){
        //leg->AddEntry(gr3,gr3->GetName(),"LP");  // to add the explanation of this colored line
    }
    
    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();
    
    //PadID++;
    
    Res_Ref_Canvas->Print(FileName.c_str());
    delete Res_Ref_Canvas;
}

void G4DoseCalcsAnalysis::GenerateSourceComputationTimeGraph(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    std::string Source_ORG;
    
    //std::cout << " The Quantity " << QuantitiesToScore << "\n" << std::endl;
    
    if(ResultParticleSourceEnergyTime.size() == 0 ){
        std::cout << " No Results data for " << QuantitiesToScore << "\n" << std::endl;
        return;
    }

    std::string FileName ;
    TMultiGraph *mg ;
    TLegend *leg ;
    
    std::cout << "\n\n-------------------------------- Results Computation Time For All Energies for each source Region"<< QuantitiesToScore << "----------------------" << std::endl ;
    
    FileName = GraphsDirectoryPath+"ComputationTimeFromSourceEnergy"+GraphsExt;
    mg = new TMultiGraph();
    leg = new TLegend();
    
    int jj = 0;
    
    // iterations on particle name
    for ( auto Mbeg = ResultParticleSourceEnergyTime.begin(); Mbeg != ResultParticleSourceEnergyTime.end(); ++Mbeg  )
    {
        GeometrySymbol = Mbeg->first;

        // iterations on particle name
        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {
            ParticleName = Nbeg->first;

            int nn = 0;
            // iterations on particle name
            for ( auto Abeg = Nbeg->second.begin(); Abeg != Nbeg->second.end(); ++Abeg  )
            {
                Source_ORG = Abeg->first;

                bool IsIn = false;
                for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                    if(SourceNamesToScore[gg] == Source_ORG){
                        IsIn = true;
                        break;
                    }
                }
                if(IsIn == false){
                    continue;
                }

                int jh = 0; // to fill the dataArray that will be used by the graph

                NumOfEne = Abeg->second.size();
                double xi[NumOfEne], yi[NumOfEne];
                //, errInE[NumOfEne] , err[NumOfEne];

                for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                {
                    xi[jh] = Cbeg->first;
                    yi[jh] = Cbeg->second;

                    //std::cout << " PARTICLE_NAME " << PARTICLE_NAME << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;

                    jh++;
                }

                TGraph* graph;
                graph = CreateGraph(NumOfEne, xi, yi);
                graph->SetName((GeometrySymbol+", "+ParticleName+", "+Source_ORG).c_str());
                graph->SetTitle(graph->GetName());

                mg->Add(setGraphData(graph, nn, jj));
                leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

                nn++;
            }
            jj++;
        }
    }
    

    std::string multiGraphTitle = "Simulation time from each source region and energy ; Energy(MeV); Computation Time (min)";
    if(PrintTitle != "yes"){multiGraphTitle = "; Energy(MeV); Computation Time (min)";}
    CreateMultiGraphParametersAndCanvas(multiGraphTitle, FileName, mg, leg);
}
void G4DoseCalcsAnalysis::GenerateStandardDeviationGraphs(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    std::string Source_ORG;
    std::string Target_ORG;
    
    std::string FileName ;
    TMultiGraph *mg ;
    TLegend *leg ;

    //std::cout << " The Quantity " << QuantitiesToScore << "\n" << std::endl;
    

    if(RelativeStandartDeviationPerCent.size() == 0 ){
        std::cout << " No Results data for " << QuantitiesToScore << "\n" << std::endl;
        return;
    }
    
    std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>> GraphsMap; // Result_Reference Particle Source Target Energy Value

    // generate one graph for each quantity and geometry and particle
    for (int bb = 0 ; bb < GeometryList.size() ; bb++) {
        for (int cc = 0 ; cc < ParticleList.size() ; cc++) {

            GeometrySymbol = GeometryList[bb];
            ParticleName = ParticleList[cc];
            for ( auto Abeg = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  )
            {

                Source_ORG = Abeg->first;

                bool IsIn = false;
                for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                    if(SourceNamesToScore[gg] == Source_ORG){
                        IsIn = true;
                        break;
                    }
                }
                if(IsIn == false){
                    continue;
                }

                // iterations on target name
                for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                {

                    Target_ORG = Cbeg->first;

                    bool IsIn = false;
                    for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                        if(TargetNamesToScore[gg] == Target_ORG){
                            IsIn = true;
                            break;
                        }
                    }
                    if(IsIn == false){
                        continue;
                    }

                    int jh = 0; // to fill the dataArray that will be used by the graph

                    NumOfEne = Cbeg->second.size();
                    double xi[NumOfEne], yi[NumOfEne];
                    //, errInE[NumOfEne] , err[NumOfEne];

                    for ( auto Dbeg = Cbeg->second.begin(); Dbeg != Cbeg->second.end(); ++Dbeg  )
                    {
                        xi[jh] = Dbeg->first;
                        yi[jh] = Dbeg->second;
                        //err[jh] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][Source_ORG][Target_ORG][PARTICLE_NAME][Dbeg->first];
                        //errInE[jh] = DefaultErrorDistance;

                        if(QuantityUseLog[QuantitiesToScore] == true && yi[jh] < MinValForLog){
                            yi[jh] = MinValForLog;
                        }

                        if(Source_ORG == "Liver" && Target_ORG == "Liver"){
                            //std::cout << " PARTICLE_NAME " << PARTICLE_NAME << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;
                        }

                        jh++;
                    }

                    TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

                    GraphsMap["Result"][Source_ORG][Target_ORG] = gr1;
                }
            }

            FileName = GraphsDirectoryPath+"SDv_"+ParticleName+"_"+QuantitiesToScore+"_"+GeometrySymbol+"_Combinations"+GraphsExt;
            mg = new TMultiGraph();
            leg = new TLegend();

            int jj = 0;

            // generate just for Cross-irradiation
            for ( auto Abeg = GraphsMap["Result"].begin(); Abeg != GraphsMap["Result"].end(); ++Abeg  )
            {

                Source_ORG = Abeg->first;

                int hh = 0;

                // iterations on particle name
                for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                {
                    Target_ORG = Cbeg->first;

                    //Cross

                    TGraph* graph = GraphsMap["Result"][Source_ORG][Target_ORG];

                    std::string graph_label = Target_ORG+"<-"+Source_ORG;
                    graph->SetName(graph_label.c_str());
                    graph->SetTitle(graph->GetName());

                    mg->Add(setGraphData(graph, hh, jj));
                    leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

                    hh++;

                }

                jj++;

            }

            std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, "<< QuantitiesToScore << " " << ParticleName << " ----------------------" << std::endl ;

            std::string multiGraphTitle = "SDv for "+GeometrySymbol+", "+ParticleName+" and combinations ;Energy(MeV); "+ QuantityUnit[QuantitiesToScore];
            if(PrintTitle != "yes"){multiGraphTitle = "; Energy(MeV); "+ QuantityUnit[QuantitiesToScore];}
            CreateMultiGraphParametersAndCanvas(multiGraphTitle, FileName, mg, leg);
        }
    }


    // generate one graph for each quantity

    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>>> GraphsForSDv; // Result_Reference Particle Source Target Energy Value
    GraphsForSDv.clear();

    for ( auto Mbeg = RelativeStandartDeviationPerCent[QuantitiesToScore].begin(); Mbeg != RelativeStandartDeviationPerCent[QuantitiesToScore].end(); ++Mbeg  )
    {
        GeometrySymbol = Mbeg->first;

        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {
            ParticleName = Nbeg->first;

            for ( auto Abeg = Nbeg->second.begin(); Abeg != Nbeg->second.end(); ++Abeg  )
            {
                Source_ORG = Abeg->first;

                bool IsIn = false;
                for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                    if(SourceNamesToScore[gg] == Source_ORG){
                        IsIn = true;
                        break;
                    }
                }
                if(IsIn == false){
                    continue;
                }

                // iterations on target name
                for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                {

                    Target_ORG = Cbeg->first;

                    bool IsIn = false;
                    for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                        if(TargetNamesToScore[gg] == Target_ORG){
                            IsIn = true;
                            break;
                        }
                    }
                    if(IsIn == false){
                        continue;
                    }

                    int jh = 0; // to fill the dataArray that will be used by the graph

                    NumOfEne = Cbeg->second.size();
                    double xi[NumOfEne], yi[NumOfEne];
                    //, errInE[NumOfEne] , err[NumOfEne];

                    for ( auto Dbeg = Cbeg->second.begin(); Dbeg != Cbeg->second.end(); ++Dbeg  )
                    {
                        xi[jh] = Dbeg->first;
                        yi[jh] = Dbeg->second;

                        if(QuantityUseLog[QuantitiesToScore] == true && yi[jh] < MinValForLog){
                            yi[jh] = MinValForLog;
                        }

                        if(Source_ORG == "Liver" && Target_ORG == "Liver"){
                            //std::cout << " PARTICLE_NAME " << PARTICLE_NAME << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;
                        }

                        jh++;
                    }

                    TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
                    GraphsForSDv[GeometrySymbol][ParticleName][Source_ORG][Target_ORG] = gr1;

                }
            }
        }
    }

    FileName = GraphsDirectoryPath+"SDv_"+QuantitiesToScore+"_AllParticlesGeometriesCombinations"+GraphsExt;
    mg = new TMultiGraph();
    leg = new TLegend();

    int jj = 0;

    for ( auto Mbeg = GraphsForSDv.begin(); Mbeg != GraphsForSDv.end(); ++Mbeg  )
    {
        GeometrySymbol = Mbeg->first;

        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {
            ParticleName = Nbeg->first;

            for ( auto Abeg = Nbeg->second.begin(); Abeg != Nbeg->second.end(); ++Abeg  )
            {
                Source_ORG = Abeg->first;

                int hh = 0;

                // iterations on particle name
                for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                {
                    Target_ORG = Cbeg->first;

                    TGraph* graph = GraphsForSDv[GeometrySymbol][ParticleName][Source_ORG][Target_ORG];

                    graph->SetName((GeometrySymbol+", "+ParticleName+", "+Target_ORG+"<-"+Source_ORG).c_str());
                    graph->SetTitle(graph->GetName());

                    mg->Add(setGraphData(graph, hh, jj));
                    leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

                    hh++;

                }

                jj++;

            }
        }
    }
    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, "<< QuantitiesToScore << " " << ParticleName << " ----------------------" << std::endl ;

    std::string multiGraphTitle = "SDv for all geometries, particles and combinations for "+QuantitiesToScore+" ;Energy(MeV); "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle != "yes"){multiGraphTitle = "; Energy(MeV); "+ QuantityUnit[QuantitiesToScore];}
    CreateMultiGraphParametersAndCanvas(multiGraphTitle, FileName, mg, leg);
}
void G4DoseCalcsAnalysis::GenerateResultInOneGraph(){ // just for cross in one graph, because the self is already generated in one graph previously
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    std::string Source_ORG;
    std::string Target_ORG;
    
    std::string FileName ;
    TMultiGraph *mg ;
    TLegend *leg ;

    std::map<std::string,std::map<std::string,TGraph*>> GraphsMap; // Result_Reference Particle Source Target Energy Value

    // all combinations for each Quantity, geometry and particle
    for (int aa = 0 ; aa < QuantityNamesToScore.size() ; aa++) {
        for (int bb = 0 ; bb < GeometryList.size() ; bb++) {
            for (int cc = 0 ; cc < ParticleList.size() ; cc++) {

                GraphsMap.clear();

                QuantitiesToScore = QuantityNamesToScore[aa];
                GeometrySymbol = GeometryList[bb];
                ParticleName = ParticleList[cc];

                for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  )
                {

                    Source_ORG = Abeg->first;

                    bool IsIn = false;
                    for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                        if(SourceNamesToScore[gg] == Source_ORG){
                            IsIn = true;
                            break;
                        }
                    }
                    if(IsIn == false){
                        continue;
                    }

                    // iterations on target name
                    for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                    {

                        Target_ORG = Cbeg->first;

                        bool IsIn = false;
                        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                            if(TargetNamesToScore[gg] == Target_ORG){
                                IsIn = true;
                                break;
                            }
                        }
                        if(IsIn == false){
                            continue;
                        }

                        int jh = 0; // to fill the dataArray that will be used by the graph

                        NumOfEne = Cbeg->second.size();
                        double xi[NumOfEne], yi[NumOfEne];
                        //, errInE[NumOfEne] , err[NumOfEne];

                        for ( auto Dbeg = Cbeg->second.begin(); Dbeg != Cbeg->second.end(); ++Dbeg  )
                        {
                            xi[jh] = Dbeg->first;
                            yi[jh] = Dbeg->second;
                            //err[jh] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][Source_ORG][Target_ORG][PARTICLE_NAME][Dbeg->first];
                            //errInE[jh] = DefaultErrorDistance;

                            if(QuantityUseLog[QuantitiesToScore] == true && yi[jh] < MinValForLog){
                                yi[jh] = MinValForLog;
                            }

                            if(Source_ORG == "Liver" && Target_ORG == "Liver"){
                                //std::cout << " PARTICLE_NAME " << PARTICLE_NAME << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;
                            }

                            jh++;
                        }

                        TGraph* gr1;
                        if(AddErrorBarInGraphs == "yes"){
                            double zi[NumOfEne], mi[NumOfEne];
                            for(int ii = 0; ii < NumOfEne; ii++){
                                zi[ii] = 0;
                                mi[ii] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]];
                                //std::cout <<  zi[ii] << " --- " << mi[ii] <<std::endl ;
                            }
                            gr1 = CreateGraphErrors(NumOfEne, xi, yi, zi, mi);

                        }else{
                            gr1 = CreateGraph (NumOfEne, xi, yi);
                        }

                        //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
                        GraphsMap[Source_ORG][Target_ORG] = gr1;
                    }
                }

                FileName = GraphsDirectoryPath+ParticleName+"_"+QuantitiesToScore+"_"+GeometrySymbol+"_"+"Combinations"+GraphsExt;
                mg = new TMultiGraph();
                leg = new TLegend();

                int jj = 0;

                // generate just for Cross-irradiation
                for ( auto Abeg = GraphsMap.begin(); Abeg != GraphsMap.end(); ++Abeg  )
                {

                    Source_ORG = Abeg->first;

                    int hh = 0;

                    // iterations on particle name
                    for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                    {
                        Target_ORG = Cbeg->first;

                        TGraph* graph = GraphsMap[Source_ORG][Target_ORG];

                        std::string graph_label = Target_ORG+"<-"+Source_ORG;
                        graph->SetName(graph_label.c_str());
                        graph->SetTitle(graph->GetName());

                        mg->Add(setGraphData(graph, hh, jj));
                        leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

                        hh++;

                    }

                    jj++;

                }

                std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, "<< QuantitiesToScore << " " << ParticleName << " ----------------------" << std::endl ;

                std::string multiGraphTitle = ParticleName +", " +GeometrySymbol + " in all combinations ;Energy(MeV); "+ QuantityUnit[QuantitiesToScore];
                if(PrintTitle != "yes"){multiGraphTitle = "; Energy(MeV); "+ QuantityUnit[QuantitiesToScore];}
                CreateMultiGraphParametersAndCanvas(multiGraphTitle, FileName, mg, leg);
            }
        }
    }

    // all sources for each Quantity, geometry and particle
    for (int aa = 0 ; aa < QuantityNamesToScore.size() ; aa++) {
        for (int bb = 0 ; bb < GeometryList.size() ; bb++) {
            for (int cc = 0 ; cc < ParticleList.size() ; cc++) {

                GraphsMap.clear();

                QuantitiesToScore = QuantityNamesToScore[aa];
                GeometrySymbol = GeometryList[bb];
                ParticleName = ParticleList[cc];

                for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  )
                {

                    Source_ORG = Abeg->first;

                    bool IsIn = false;
                    for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                        if(SourceNamesToScore[gg] == Source_ORG){
                            IsIn = true;
                            break;
                        }
                    }
                    if(IsIn == false){
                        continue;
                    }

                    NumOfEne = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Source_ORG].size();
                    double xi[NumOfEne], yi[NumOfEne];
                    //, errInE[NumOfEne] , err[NumOfEne];

                    int jh = 0; // to fill the dataArray that will be used by the graph
                    for ( auto Dbeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Source_ORG].begin(); Dbeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Source_ORG].end(); ++Dbeg  )
                    {
                        xi[jh] = Dbeg->first;
                        yi[jh] = Dbeg->second;
                        //err[jh] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][Source_ORG][Source_ORG][PARTICLE_NAME][Dbeg->first];
                        //errInE[jh] = DefaultErrorDistance;

                        if(QuantityUseLog[QuantitiesToScore] == true && yi[jh] < MinValForLog){
                            yi[jh] = MinValForLog;
                        }

                        jh++;
                    }

                    TGraph* gr1;
                    if(AddErrorBarInGraphs == "yes"){
                        double zi[NumOfEne], mi[NumOfEne];
                        for(int ii = 0; ii < NumOfEne; ii++){
                            zi[ii] = 0;
                            mi[ii] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Source_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]];
                            //std::cout <<  zi[ii] << " --- " << mi[ii] <<std::endl ;
                        }
                        gr1 = CreateGraphErrors(NumOfEne, xi, yi, zi, mi);

                    }else{
                        gr1 = CreateGraph (NumOfEne, xi, yi);
                    }

                    //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
                    GraphsMap[Source_ORG][Source_ORG] = gr1;
                }

                FileName = GraphsDirectoryPath+ParticleName+"_"+QuantitiesToScore+"_"+GeometrySymbol+"_"+"Sources"+GraphsExt;
                mg = new TMultiGraph();
                leg = new TLegend();

                int jj = 0;
                int hh = 0;

                // generate just for Cross-irradiation
                for ( auto Abeg = GraphsMap.begin(); Abeg != GraphsMap.end(); ++Abeg  )
                {

                    Source_ORG = Abeg->first;

                    TGraph* graph = GraphsMap[Source_ORG][Source_ORG];

                    //QuantitiesToScore+" "+GeometrySymbol+" "+ParticleName+" " +
                    std::string graph_label = Source_ORG;
                    graph->SetName(graph_label.c_str());
                    graph->SetTitle(graph->GetName());

                    mg->Add(setGraphData(graph, hh, jj));
                    leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

                    hh++;

                    jj++;

                }

                std::cout << "\n\n-------------------------------- In One Graph Result, for all sources, "<< QuantitiesToScore << " " << ParticleName << " ----------------------" << std::endl ;

                std::string multiGraphTitle = ParticleName +", " +GeometrySymbol + " in all sources ;Energy(MeV); "+ QuantityUnit[QuantitiesToScore];
                if(PrintTitle != "yes"){multiGraphTitle = "; Energy(MeV); "+ QuantityUnit[QuantitiesToScore];}
                CreateMultiGraphParametersAndCanvas(multiGraphTitle, FileName, mg, leg);
            }
        }
    }

    // all combinations for all Quantities, geometries and particles
    int jj = 0;
    for (int aa = 0 ; aa < QuantityNamesToScore.size() ; aa++) {

        FileName = GraphsDirectoryPath+QuantitiesToScore+"_AllParticlesGeometriesCombinations"+GraphsExt;
        mg = new TMultiGraph();
        leg = new TLegend();

        for (int cc = 0 ; cc < ParticleList.size() ; cc++) {
            for (int bb = 0 ; bb < GeometryList.size() ; bb++) {

                QuantitiesToScore = QuantityNamesToScore[aa];
                GeometrySymbol = GeometryList[bb];
                ParticleName = ParticleList[cc];

                for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  )
                {

                    Source_ORG = Abeg->first;

                    bool IsIn = false;
                    for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                        if(SourceNamesToScore[gg] == Source_ORG){
                            IsIn = true;
                            break;
                        }
                    }
                    if(IsIn == false){
                        continue;
                    }

                    int hh = 0;

                    // iterations on target name
                    for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                    {

                        Target_ORG = Cbeg->first;

                        bool IsIn = false;
                        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                            if(TargetNamesToScore[gg] == Target_ORG){
                                IsIn = true;
                                break;
                            }
                        }
                        if(IsIn == false){
                            continue;
                        }

                        int jh = 0; // to fill the dataArray that will be used by the graph

                        NumOfEne = Cbeg->second.size();
                        double xi[NumOfEne], yi[NumOfEne];
                        //, errInE[NumOfEne] , err[NumOfEne];

                        for ( auto Dbeg = Cbeg->second.begin(); Dbeg != Cbeg->second.end(); ++Dbeg  )
                        {
                            xi[jh] = Dbeg->first;
                            yi[jh] = Dbeg->second;
                            //err[jh] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][Source_ORG][Target_ORG][PARTICLE_NAME][Dbeg->first];
                            //errInE[jh] = DefaultErrorDistance;

                            if(QuantityUseLog[QuantitiesToScore] == true && yi[jh] < MinValForLog){
                                yi[jh] = MinValForLog;
                            }

                            if(Source_ORG == "Liver" && Target_ORG == "Liver"){
                                //std::cout << " PARTICLE_NAME " << PARTICLE_NAME << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;
                            }

                            jh++;
                        }

                        TGraph* graph;
                        if(AddErrorBarInGraphs == "yes"){
                            double zi[NumOfEne], mi[NumOfEne];
                            for(int ii = 0; ii < NumOfEne; ii++){
                                zi[ii] = 0;
                                mi[ii] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]];
                                //std::cout <<  zi[ii] << " --- " << mi[ii] <<std::endl ;
                            }
                            graph = CreateGraphErrors(NumOfEne, xi, yi, zi, mi);

                        }else{
                            graph = CreateGraph (NumOfEne, xi, yi);
                        }

                        std::string graph_label = GeometrySymbol+", "+ParticleName+", "+Target_ORG+"<-"+Source_ORG;
                        graph->SetName(graph_label.c_str());
                        graph->SetTitle(graph->GetName());

                        mg->Add(setGraphData(graph, hh, jj));
                        leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

                        hh++;
                    }
                    jj++;
                }
            }
        }

        std::cout << "\n\n-------------------------------- In One Graph Result, for all geometries, particles, and combinations, "<< QuantitiesToScore << " " << ParticleName << " ----------------------" << std::endl ;

        std::string multiGraphTitle = " For all geometries, particles and combinations ;Energy(MeV); "+ QuantityUnit[QuantitiesToScore];
        if(PrintTitle != "yes"){multiGraphTitle = "; Energy(MeV); "+ QuantityUnit[QuantitiesToScore];}
        CreateMultiGraphParametersAndCanvas(multiGraphTitle, FileName, mg, leg);
    }

    // all sources for all Quantities, geometries and particles
    jj = 0;
    for (int aa = 0 ; aa < QuantityNamesToScore.size() ; aa++) {

        FileName = GraphsDirectoryPath+QuantitiesToScore+"_AllParticlesGeometriesSources"+GraphsExt;
        mg = new TMultiGraph();
        leg = new TLegend();

        for (int cc = 0 ; cc < ParticleList.size() ; cc++) {
            for (int bb = 0 ; bb < GeometryList.size() ; bb++) {

                QuantitiesToScore = QuantityNamesToScore[aa];
                GeometrySymbol = GeometryList[bb];
                ParticleName = ParticleList[cc];

                int hh = 0;

                for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  )
                {

                    Source_ORG = Abeg->first;

                    bool IsIn = false;
                    for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                        if(SourceNamesToScore[gg] == Source_ORG){
                            IsIn = true;
                            break;
                        }
                    }
                    if(IsIn == false){
                        continue;
                    }

                    NumOfEne = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Source_ORG].size();
                    double xi[NumOfEne], yi[NumOfEne];
                    //, errInE[NumOfEne] , err[NumOfEne];

                    int jh = 0; // to fill the dataArray that will be used by the graph
                    for ( auto Dbeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Source_ORG].begin(); Dbeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Source_ORG].end(); ++Dbeg  )
                    {
                        xi[jh] = Dbeg->first;
                        yi[jh] = Dbeg->second;
                        //err[jh] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][Source_ORG][Target_ORG][PARTICLE_NAME][Dbeg->first];
                        //errInE[jh] = DefaultErrorDistance;

                        if(QuantityUseLog[QuantitiesToScore] == true && yi[jh] < MinValForLog){
                            yi[jh] = MinValForLog;
                        }

                        jh++;
                    }

                    TGraph* graph;
                    if(AddErrorBarInGraphs == "yes"){
                        double zi[NumOfEne], mi[NumOfEne];
                        for(int ii = 0; ii < NumOfEne; ii++){
                            zi[ii] = 0;
                            mi[ii] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii]];
                            //std::cout <<  zi[ii] << " --- " << mi[ii] <<std::endl ;
                        }
                        graph = CreateGraphErrors(NumOfEne, xi, yi, zi, mi);

                    }else{
                        graph = CreateGraph (NumOfEne, xi, yi);
                    }

                    std::string graph_label = GeometrySymbol+", "+ParticleName+", "+Source_ORG;
                    graph->SetName(graph_label.c_str());
                    graph->SetTitle(graph->GetName());

                    mg->Add(setGraphData(graph, hh, jj));
                    leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

                    hh++;

                    jj++;
                }
            }
        }

        std::cout << "\n\n-------------------------------- In One Graph Result, for all geometries, particles, and sources, "<< QuantitiesToScore << " " << ParticleName << " ----------------------" << std::endl ;

        std::string multiGraphTitle = " For all geometries, particles and sources ;Energy(MeV); "+ QuantityUnit[QuantitiesToScore];
        if(PrintTitle != "yes"){multiGraphTitle = "; Energy(MeV); "+ QuantityUnit[QuantitiesToScore];}
        CreateMultiGraphParametersAndCanvas(multiGraphTitle, FileName, mg, leg);
    }

}
void G4DoseCalcsAnalysis::GenerateResultForEachParticleEnergyInOneGraph(){

    std::string Source_ORG;
    std::string Target_ORG;

    std::string FileName ;
    TCanvas* ResCanvas ;
    TMultiGraph *mg ;
    TLegend *leg ;
    TText* t = new TText();
    std::string multiGraphTitle = "";

    std::map<std::string,std::map<std::string,std::map<double,TGraph*>>> AllGeometriesRadiotracersGraphs; // Result_Reference Particle Source Target Energy Value
    std::vector<std::string> XLabelsVec;
    std::vector<int> XLabelsIDVec;

    for ( auto AA = ResEnergies[QuantitiesToScore].begin(); AA != ResEnergies[QuantitiesToScore].end(); ++AA  ){
        GeometrySymbol = AA->first;
        for ( auto BB = AA->second.begin(); BB != AA->second.end(); ++BB  ){
            ParticleName = BB->first;
            for(int cv = 0; cv < BB->second.size(); cv++){
                EnergyGraphValue = BB->second[cv];

                NumOfEne = SourceNamesToScore.size()*TargetNamesToScore.size();
                double xi[NumOfEne]; double yi[NumOfEne];

                int jj = 0;
                for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
                    for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                        Source_ORG = SourceNamesToScore[ee];
                        Target_ORG = TargetNamesToScore[gg];
                        xi[jj] = jj;
                        yi[jj] = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][EnergyGraphValue];

                        if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                            yi[jj] = MinValForLog;
                        }

                        if(Source_ORG == "Liver" && Target_ORG == "Liver"){
                            //std::cout << " PARTICLE_NAME " << PARTICLE_NAME << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;
                        }
                        jj++;
                    }
                }
                TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

                //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
                AllGeometriesRadiotracersGraphs[GeometrySymbol][ParticleName][EnergyGraphValue] = gr1;
            }
        }
    }


    FileName = GraphsDirectoryPath+QuantitiesToScore+"_Geometries_ParticleEnergies_Combinations"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    int jj = 0;
    for ( auto Mbeg = AllGeometriesRadiotracersGraphs.begin(); Mbeg != AllGeometriesRadiotracersGraphs.end(); ++Mbeg  )
    {

        GeometrySymbol = Mbeg->first;

        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {
            ParticleName = Nbeg->first;

            int hh = 0;
            for ( auto Obeg = Nbeg->second.begin(); Obeg != Nbeg->second.end(); ++Obeg  )
            {

                EnergyGraphValue = Obeg->first;

                TGraph* graph = AllGeometriesRadiotracersGraphs[GeometrySymbol][ParticleName][EnergyGraphValue];

                std::ostringstream graph_label; graph_label << GeometrySymbol<<", "<< ParticleName<< ", "<<EnergyGraphValue;
                graph->SetName(graph_label.str().c_str());
                graph->SetTitle(graph->GetName());

                mg->Add(setGraphData(graph, hh, jj));
                leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
                hh++;
                jj++;
            }
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries, particles, and energies "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " in all combinations, geometries, particles and energies ;target<-source; "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; target<-source; "+ QuantityUnit[QuantitiesToScore];
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->GetHistogram()->GetXaxis()->SetTickLength(0); // because we will print text
    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);

    double y = /*0.67; */gPad->GetUymin() ;//- 0.2*mg->GetYaxis()->GetBinWidth(1);

    int jh = 0;
    for(int aa = 0; aa < SourceNamesToScore.size(); aa++){
        for(int bb = 0; bb < TargetNamesToScore.size(); bb++){
            double z=jh;
            t->DrawText(z,y,const_cast<char*>((TargetNamesToScore[bb]+"<-"+SourceNamesToScore[aa]).c_str()));
            //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
            jh++;
        }
    }

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;

    // /////////////////////////////////////// for reference relative differences

    if(GenerateRelativeErrGraph == "yes" && CompareReferenceNames.size() != 0 && RefFilePaths.size() != 0 && ReferenceTable.size() != 0){}
    else{return;}

    for ( auto AA = ResEnergies[QuantitiesToScore].begin(); AA != ResEnergies[QuantitiesToScore].end(); ++AA  ){
        GeometrySymbol = AA->first;
        for ( auto BB = AA->second.begin(); BB != AA->second.end(); ++BB  ){
            ParticleName = BB->first;
            for(int cv = 0; cv < BB->second.size(); cv++){
                EnergyGraphValue = BB->second[cv];

                NumOfEne = SourceNamesToScore.size()*TargetNamesToScore.size();
                double xi[NumOfEne]; double yi[NumOfEne];

                int jj = 0;
                for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
                    for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                        Source_ORG = SourceNamesToScore[ee];
                        Target_ORG = TargetNamesToScore[gg];
                        xi[jj] = jj;
                        double a1 = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][EnergyGraphValue] ;
                        double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][EnergyGraphValue];
                        double a3 = RelativeDifferenceCalculation(a1,a2);
                        if(a3==NULL || a3==0.){a3 = MinValForLog;}
                        yi[jj] = a3;

                        if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                            yi[jj] = MinValForLog;
                        }

                        if(Source_ORG == "Liver" && Target_ORG == "Liver"){
                            //std::cout << " PARTICLE_NAME " << PARTICLE_NAME << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;
                        }
                        jj++;
                    }
                }
                TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

                //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
                AllGeometriesRadiotracersGraphs[GeometrySymbol][ParticleName][EnergyGraphValue] = gr1;
            }
        }
    }

    FileName = GraphsDirectoryPath+"RelativeDiff_"+ParticleName+"_"+QuantitiesToScore+"_"+CompareReferenceName+"_Geometries_ParticleEnergies_Combinations"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for ( auto Mbeg = AllGeometriesRadiotracersGraphs.begin(); Mbeg != AllGeometriesRadiotracersGraphs.end(); ++Mbeg  )
    {

        GeometrySymbol = Mbeg->first;

        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {
            ParticleName = Nbeg->first;

            int hh = 0;
            for ( auto Obeg = Nbeg->second.begin(); Obeg != Nbeg->second.end(); ++Obeg  )
            {

                EnergyGraphValue = Obeg->first;

                TGraph* graph = AllGeometriesRadiotracersGraphs[GeometrySymbol][ParticleName][EnergyGraphValue];

                std::ostringstream graph_label; graph_label << GeometrySymbol<<", "<< ParticleName<< ", "<<EnergyGraphValue;
                graph->SetName(graph_label.str().c_str());
                graph->SetTitle(graph->GetName());

                mg->Add(setGraphData(graph, hh, jj));
                leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
                hh++;
                jj++;
            }
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries, particles, and energies "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = "DoseCalcs and "+CompareReferenceName +" differences for all combinations, geometries, particles and energies ;target<-source; "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; target<-source; "+ QuantityUnit[QuantitiesToScore];
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->GetHistogram()->GetXaxis()->SetTickLength(0); // because we will print text
    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);

    y = /*0.67; */gPad->GetUymin() ;//- 0.2*mg->GetYaxis()->GetBinWidth(1);

    jh = 0;
    for(int aa = 0; aa < SourceNamesToScore.size(); aa++){
        for(int bb = 0; bb < TargetNamesToScore.size(); bb++){
            double z=jh;
            t->DrawText(z,y,const_cast<char*>((TargetNamesToScore[bb]+"<-"+SourceNamesToScore[aa]).c_str()));
            //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
            jh++;
        }
    }

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;
}
void G4DoseCalcsAnalysis::GenerateResultReferenceInOneGraph(){
    
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;
    
    std::string Source_ORG;
    std::string Target_ORG;

    TCanvas* Can = new TCanvas();
    std::string FileName ;
    std::string multiGraphTitle;
    TMultiGraph *mg ;
    TLegend *leg ;

    mg = new TMultiGraph();
    leg = new TLegend();

    std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>> GraphsMap; // Result_Reference Particle Source Target Energy Value
    
    for (int aa = 0 ; aa < QuantityNamesToScore.size() ; aa++) {

        for (int bb = 0 ; bb < GeometryList.size() ; bb++) {
            for (int cc = 0 ; cc < ParticleList.size() ; cc++) {

                GraphsMap.clear();
                QuantitiesToScore = QuantityNamesToScore[aa];
                GeometrySymbol = GeometryList[bb];
                ParticleName = ParticleList[cc];

                mg = new TMultiGraph();
                leg = new TLegend();

                // iterations on particle name
                for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  )
                {
                    Source_ORG = Abeg->first;

                    bool IsIn = false;
                    for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                        if(SourceNamesToScore[gg] == Source_ORG){
                            IsIn = true;
                            break;
                        }
                    }
                    if(IsIn == false){
                        continue;
                    }

                    // iterations on target name
                    for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                    {

                        Target_ORG = Cbeg->first;

                        bool IsIn = false;
                        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                            if(TargetNamesToScore[gg] == Target_ORG){
                                IsIn = true;
                                break;
                            }
                        }
                        if(IsIn == false){
                            continue;
                        }

                        int jh = 0; // to fill the dataArray that will be used by the graph

                        NumOfEne = Cbeg->second.size();
                        double xi[NumOfEne], yi[NumOfEne];
                        //, errInE[NumOfEne] , err[NumOfEne];

                        for ( auto Dbeg = Cbeg->second.begin(); Dbeg != Cbeg->second.end(); ++Dbeg  )
                        {
                            xi[jh] = Dbeg->first;
                            yi[jh] = Dbeg->second;
                            //err[jh] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][Source_ORG][Target_ORG][ParticleName][Dbeg->first];
                            //errInE[jh] = DefaultErrorDistance;

                            if(QuantityUseLog[QuantitiesToScore] == true && yi[jh] < MinValForLog){
                                yi[jh] = MinValForLog;
                            }

                            //if(Source_ORG == "Liver" && Target_ORG == "Liver"){
                            //std::cout << " ParticleName " << ParticleName << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;
                            //}

                            jh++;
                        }

                        TGraph* gr1;
                        if(AddErrorBarInGraphs == "yes"){
                            double zi[NumOfEne], mi[NumOfEne];
                            for(int ii = 0; ii < NumOfEne; ii++){
                                zi[ii] = 0;
                                mi[ii] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][xi[ii]];
                                //std::cout <<  zi[ii] << " --- " << mi[ii] <<std::endl ;
                            }
                            gr1 = CreateGraphErrors(NumOfEne, xi, yi, zi, mi);

                        }else{
                            gr1 = CreateGraph (NumOfEne, xi, yi);
                        }

                        GraphsMap["Result"][Source_ORG][Target_ORG] = gr1;

                        double vi[NumOfEne], wi[NumOfEne];
                        for(int ii = 0; ii < NumOfEne; ii++){
                            //vi[ii] = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii];
                            wi[ii] = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][xi[ii]];
                            if(QuantityUseLog[QuantitiesToScore] == true && wi[ii] < MinValForLog){
                                wi[ii] = MinValForLog;
                            }
                            //std::cout << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " vi[ii]= " << vi[ii] << " wi[ii] " << wi[ii] <<std::endl ;
                        }

                        TGraph* gr2 = CreateGraph (NumOfEne, xi, wi);

                        GraphsMap["Reference"][Source_ORG][Target_ORG] = gr2;

                        if(RefFilePaths.size() > 1){
                            TGraph* gr3;
                            double zi[NumOfEne], mi[NumOfEne];
                            for(int ii = 0; ii < NumOfEne; ii++){
                                //zi[ii] = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii];
                                mi[ii] = ReferenceTable2[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][xi[ii]];
                                if(QuantityUseLog[QuantitiesToScore] == true && mi[ii] < MinValForLog){
                                    mi[ii] = MinValForLog;
                                }
                                //std::cout << " zi[ii] " <<  zi[ii] << " mi[ii] "<< mi[ii] <<std::endl ;
                            }

                            //double* mi = AccumulateThirdRefDataToGraphs(ParticleName, Source_ORG, Target_ORG);
                            gr3 = CreateGraph (NumOfEne , xi, mi ); //this graph is related if we want to remove the zeros of SAF from Graph, if not we have to uncomment the graph1 below

                            GraphsMap["Reference1"][Source_ORG][Target_ORG] = gr3;
                        }

                    }

                }

                if( ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName].size() == 0 ){
                    std::cout << " No reference data to compare for " << QuantitiesToScore << "\n" << std::endl;
                    return;
                }

                int jj = 0;

                // generate just for Cross-irradiation
                for ( auto Abeg = GraphsMap["Result"].begin(); Abeg != GraphsMap["Result"].end(); ++Abeg  )
                {
                    Source_ORG = Abeg->first;

                    for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                    {
                        Target_ORG = Cbeg->first;
                        TGraph* graph1;
                        TGraph* graph2;
                        TGraph* graph3;

                        graph1 = GraphsMap["Result"][Source_ORG][Target_ORG];

                        // in case of the source is not defined in reference file
                        if(GraphsMap["Reference"][Source_ORG][Target_ORG] == NULL){
                            continue;
                        }

                        graph2 = GraphsMap["Reference"][Source_ORG][Target_ORG];
                        graph2->RemovePoint();
                        std::string graph_label1 = "DoseCalcs "+Target_ORG+"<-"+Source_ORG;
                        graph1->SetName(graph_label1.c_str());
                        graph1->SetTitle(graph1->GetName());
                        mg->Add(setGraphData(graph1, 1, jj));

                        leg->AddEntry(graph1,graph1->GetName(),"LP");  // to add the explanation of this colored line

                        std::string graph_label2 = CompareReferenceName +" "+Target_ORG+"<-"+Source_ORG;

                        //std::cout << graph_label2 << std::endl ;
                        graph2->SetName(graph_label2.c_str());
                        graph2->SetTitle(graph2->GetName());

                        mg->Add(setGraphData(graph2, 2, jj));
                        leg->AddEntry(graph2,graph2->GetName(),"LP");  // to add the explanation of this colored line

                        if(RefFilePaths.size() > 1){
                            graph3 = GraphsMap["Reference1"][Source_ORG][Target_ORG];
                            std::string graph_label3 = CompareReferenceNames[1] +" "+ Target_ORG+"<-"+Source_ORG;
                            graph3->SetName(graph_label3.c_str());
                            graph3->SetTitle(graph3->GetName());
                            mg->Add(setGraphData(graph3, 3, jj));
                            leg->AddEntry(graph3,graph3->GetName(),"LP");  // to add the explanation of this colored line
                        }

                        jj++;

                    }
                }

                FileName = GraphsDirectoryPath+ParticleName+"_"+QuantitiesToScore+"_"+GeometrySymbol+"_Combinations_ComparedTo"+CompareReferenceName+GraphsExt;
                multiGraphTitle = ParticleName +", "+GeometrySymbol+ " in all combinations ;Energy(MeV); "+ QuantityUnit[QuantitiesToScore];
                if(PrintTitle != "yes"){multiGraphTitle = "; Energy(MeV); "+ QuantityUnit[QuantitiesToScore];}
                CreateMultiGraphParametersAndCanvas(multiGraphTitle, FileName, mg, leg);
            }
        }
    }

    int jj = 0;
    for (int aa = 0 ; aa < QuantityNamesToScore.size() ; aa++) {

        mg = new TMultiGraph();
        leg = new TLegend();

        for (int bb = 0 ; bb < GeometryList.size() ; bb++) {
            for (int cc = 0 ; cc < ParticleList.size() ; cc++) {

                QuantitiesToScore = QuantityNamesToScore[aa];
                GeometrySymbol = GeometryList[bb];
                ParticleName = ParticleList[cc];

                // iterations on particle name
                for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  )
                {
                    Source_ORG = Abeg->first;

                    bool IsIn = false;
                    for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                        if(SourceNamesToScore[gg] == Source_ORG){
                            IsIn = true;
                            break;
                        }
                    }
                    if(IsIn == false){
                        continue;
                    }

                    // iterations on target name
                    for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                    {

                        Target_ORG = Cbeg->first;

                        bool IsIn = false;
                        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                            if(TargetNamesToScore[gg] == Target_ORG){
                                IsIn = true;
                                break;
                            }
                        }
                        if(IsIn == false){
                            continue;
                        }

                        int jh = 0;

                        NumOfEne = Cbeg->second.size();
                        double xi[NumOfEne], yi[NumOfEne];

                        for ( auto Dbeg = Cbeg->second.begin(); Dbeg != Cbeg->second.end(); ++Dbeg  )
                        {
                            xi[jh] = Dbeg->first;
                            yi[jh] = Dbeg->second;

                            if(QuantityUseLog[QuantitiesToScore] == true && yi[jh] < MinValForLog){
                                yi[jh] = MinValForLog;
                            }

                            jh++;
                        }

                        TGraph* graph1;
                        TGraph* graph2;
                        TGraph* graph3;

                        if(AddErrorBarInGraphs == "yes"){
                            double zi[NumOfEne], mi[NumOfEne];
                            for(int ii = 0; ii < NumOfEne; ii++){
                                zi[ii] = 0;
                                mi[ii] = RelativeStandartDeviationPerCent[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][xi[ii]];
                                //std::cout <<  zi[ii] << " --- " << mi[ii] <<std::endl ;
                            }
                            graph1 = CreateGraphErrors(NumOfEne, xi, yi, zi, mi);

                        }else{
                            graph1 = CreateGraph (NumOfEne, xi, yi);
                        }

                        double vi[NumOfEne], wi[NumOfEne];
                        for(int ii = 0; ii < NumOfEne; ii++){
                            //vi[ii] = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii];
                            wi[ii] = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][xi[ii]];
                            if(QuantityUseLog[QuantitiesToScore] == true && wi[ii] < MinValForLog){
                                wi[ii] = MinValForLog;
                            }
                            //std::cout << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " vi[ii]= " << vi[ii] << " wi[ii] " << wi[ii] <<std::endl ;
                        }

                        graph2 = CreateGraph (NumOfEne, xi, wi);

                        if(RefFilePaths.size() > 1){
                            double zi[NumOfEne], mi[NumOfEne];
                            for(int ii = 0; ii < NumOfEne; ii++){
                                //zi[ii] = ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][ii];
                                mi[ii] = ReferenceTable2[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][xi[ii]];
                                if(QuantityUseLog[QuantitiesToScore] == true && mi[ii] < MinValForLog){
                                    mi[ii] = MinValForLog;
                                }
                                //std::cout << " zi[ii] " <<  zi[ii] << " mi[ii] "<< mi[ii] <<std::endl ;
                            }

                            //double* mi = AccumulateThirdRefDataToGraphs(ParticleName, Source_ORG, Target_ORG);
                            graph3 = CreateGraph (NumOfEne , xi, mi ); //this graph is related if we want to remove the zeros of SAF from Graph, if not we have to uncomment the graph1 below
                        }

                        std::string graph_label1 = "DoseCalcs, " +GeometrySymbol +", "+ParticleName+", "+ Target_ORG+"<-"+Source_ORG;
                        graph1->SetName(graph_label1.c_str());
                        graph1->SetTitle(graph1->GetName());
                        mg->Add(setGraphData(graph1, 1, jj));
                        leg->AddEntry(graph1,graph1->GetName(),"LP");  // to add the explanation of this colored line

                        std::string graph_label2 = CompareReferenceName + ", " +GeometrySymbol +", "+ParticleName+", "+ Target_ORG+"<-"+Source_ORG;

                        //std::cout << graph_label2 << std::endl ;
                        graph2->SetName(graph_label2.c_str());
                        graph2->SetTitle(graph2->GetName());
                        mg->Add(setGraphData(graph2, 2, jj));
                        leg->AddEntry(graph2,graph2->GetName(),"LP");  // to add the explanation of this colored line

                        if(RefFilePaths.size() > 1){
                            std::string graph_label3 = CompareReferenceNames[1] + ", " +GeometrySymbol +", "+ParticleName+", "+ Target_ORG+"<-"+Source_ORG;
                            graph3->SetName(graph_label3.c_str());
                            graph3->SetTitle(graph3->GetName());
                            mg->Add(setGraphData(graph3, 3, jj));
                            leg->AddEntry(graph3,graph3->GetName(),"LP");  // to add the explanation of this colored line
                        }

                        jj++;
                    }
                }
            }
        }
        FileName = GraphsDirectoryPath+QuantitiesToScore+"_AllParticlesGeometriesCombinations_ComparedTo"+CompareReferenceName+GraphsExt;
        multiGraphTitle = " For all geometries, particles and combinations ;Energy(MeV); "+ QuantityUnit[QuantitiesToScore];
        if(PrintTitle != "yes"){multiGraphTitle = "; Energy(MeV); "+ QuantityUnit[QuantitiesToScore];}
        CreateMultiGraphParametersAndCanvas(multiGraphTitle, FileName, mg, leg);
    }

}
void G4DoseCalcsAnalysis::GenerateComparisonFactorGraphs(){

    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;

    std::string Source_ORG;
    std::string Target_ORG;

    std::string FileName ;
    TMultiGraph *mg ;
    TLegend *leg ;

    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>>> GraphsForSDv;
    std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>> GraphsMap; // Result_Reference Particle Source Target Energy Value

    for (int bb = 0 ; bb < GeometryList.size() ; bb++) {
        for (int cc = 0 ; cc < ParticleList.size() ; cc++) {

            GeometrySymbol = GeometryList[bb];
            ParticleName = ParticleList[cc];

            // iterations on particle name
            for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  )
            {
                Source_ORG = Abeg->first;

                bool IsIn = false;
                for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                    if(SourceNamesToScore[gg] == Source_ORG){
                        IsIn = true;
                        break;
                    }
                }
                if(IsIn == false){
                    continue;
                }

                // iterations on target name
                for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg )
                {
                    Target_ORG = Cbeg->first;

                    bool IsIn = false;
                    for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                        if(TargetNamesToScore[gg] == Target_ORG){
                            IsIn = true;
                            break;
                        }
                    }
                    if(IsIn == false){
                        continue;
                    }

                    // iterations on particle energies
                    for ( auto Dbeg = Cbeg->second.begin(); Dbeg != Cbeg->second.end(); ++Dbeg )
                    {
                        //ResEneValMap[Dbeg->first] = Dbeg->second;
                        //std::cout << ResEneValMap[Dbeg->first] << "" <<  std::endl ;

                        double a1 = Dbeg->second ;
                        double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][Dbeg->first];
                        double a3 = RelativeDifferenceCalculation(a1,a2);
                        if(a3==NULL || a3==0.){a3 = MinValForLog;}

                        ResRefErrCompTables[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][Dbeg->first] = a3;
                    }
                }
            }

            // iterations on source name
            for ( auto Abeg = ResRefErrCompTables[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResRefErrCompTables[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  )
            {

                Source_ORG = Abeg->first;

                // iterations on target name
                for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                {

                    Target_ORG = Cbeg->first;

                    bool IsIn = false;
                    for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                        if(TargetNamesToScore[gg] == Target_ORG){
                            IsIn = true;
                            break;
                        }
                    }
                    if(IsIn == false){
                        continue;
                    }

                    int jh = 0; // to fill the dataArray that will be used by the graph

                    NumOfEne = Cbeg->second.size();
                    double xi[NumOfEne], yi[NumOfEne];
                    //, errInE[NumOfEne] , err[NumOfEne];

                    for ( auto Dbeg = Cbeg->second.begin(); Dbeg != Cbeg->second.end(); ++Dbeg  )
                    {
                        xi[jh] = Dbeg->first;
                        yi[jh] = Dbeg->second;

                        if(Source_ORG == "Liver" && Target_ORG == "Liver"){
                            //std::cout << " PARTICLE_NAME " << PARTICLE_NAME << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;
                        }

                        jh++;
                    }

                    TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

                    //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
                    GraphsMap["Ratio"][Source_ORG][Target_ORG] = gr1;
                }
            }

            FileName = GraphsDirectoryPath+"RelativeDiff_"+ParticleName+"_"+QuantitiesToScore+"_"+CompareReferenceName+"_"+GeometrySymbol+GraphsExt;
            mg = new TMultiGraph();
            leg = new TLegend();

            int jj = 0;

            // generate just for Cross-irradiation
            for ( auto Abeg = GraphsMap["Ratio"].begin(); Abeg != GraphsMap["Ratio"].end(); ++Abeg  )
            {

                Source_ORG = Abeg->first;

                int hh = 0;

                // iterations on particle name
                for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                {
                    Target_ORG = Cbeg->first;

                    TGraph* graph = GraphsMap["Ratio"][Source_ORG][Target_ORG];

                    std::string graph_label = Target_ORG+"<-"+Source_ORG;
                    graph->SetName(graph_label.c_str());
                    graph->SetTitle(graph->GetName());

                    mg->Add(setGraphData(graph, hh, jj));
                    leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

                    hh++;

                }

                jj++;

            }

            std::cout << "\n\n-------------------------------- In One Graph ratios DoseCalcs" << CompareReferenceName <<", for all combinations, "<< QuantitiesToScore << " " << ParticleName << " ----------------------" << std::endl ;

            std::string multiGraphTitle = "DoseCalcs and "+CompareReferenceName +" differences for " + GeometrySymbol + ", " + ParticleName + " and combinations  ;Energy(MeV); "+ QuantityUnit[QuantitiesToScore];
            if(PrintTitle != "yes"){multiGraphTitle = "; Energy(MeV); "+ DiffExp;}
            CreateMultiGraphParametersAndCanvas(multiGraphTitle, FileName, mg, leg);
        }
    }

    GraphsForSDv.clear();

    for (int bb = 0 ; bb < GeometryList.size() ; bb++) {
        for (int cc = 0 ; cc < ParticleList.size() ; cc++) {

            GeometrySymbol = GeometryList[bb];
            ParticleName = ParticleList[cc];

            // iterations on particle name
            for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  )
            {
                Source_ORG = Abeg->first;

                bool IsIn = false;
                for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                    if(SourceNamesToScore[gg] == Source_ORG){
                        IsIn = true;
                        break;
                    }
                }
                if(IsIn == false){
                    continue;
                }

                // iterations on target name
                for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg )
                {
                    Target_ORG = Cbeg->first;

                    bool IsIn = false;
                    for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                        if(TargetNamesToScore[gg] == Target_ORG){
                            IsIn = true;
                            break;
                        }
                    }
                    if(IsIn == false){
                        continue;
                    }

                    int jh = 0; // to fill the dataArray that will be used by the graph

                    NumOfEne = Cbeg->second.size();
                    double xi[NumOfEne], yi[NumOfEne];

                    // iterations on particle energies
                    for ( auto Dbeg = Cbeg->second.begin(); Dbeg != Cbeg->second.end(); ++Dbeg )
                    {
                        double a1 = Dbeg->second ;
                        double a2 = ReferenceTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][Dbeg->first];
                        double a3 = RelativeDifferenceCalculation(a1,a2);
                        if(a3==NULL || a3==0.){a3 = MinValForLog;}

                        ResRefErrCompTables[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][Target_ORG][Dbeg->first] = a3;

                        xi[jh] = Dbeg->first;
                        yi[jh] = a3;

                        //std::cout << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;

                        jh++;
                    }

                    TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

                    GraphsForSDv[GeometrySymbol][ParticleName][Source_ORG][Target_ORG] = gr1;
                }
            }
        }
    }

    FileName = GraphsDirectoryPath+"RelativeDiff_"+QuantitiesToScore+"_AllParticlesGeometriesCombinations"+GraphsExt;
    mg = new TMultiGraph();
    leg = new TLegend();

    int jj = 0;

    // generate just for Cross-irradiation
    for ( auto Mbeg = GraphsForSDv.begin(); Mbeg != GraphsForSDv.end(); ++Mbeg  )
    {
        GeometrySymbol = Mbeg->first;

        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {
            ParticleName = Nbeg->first;

            for ( auto Abeg = Nbeg->second.begin(); Abeg != Nbeg->second.end(); ++Abeg  )
            {

                Source_ORG = Abeg->first;

                int hh = 0;

                // iterations on particle name
                for ( auto Cbeg = Abeg->second.begin(); Cbeg != Abeg->second.end(); ++Cbeg  )
                {
                    Target_ORG = Cbeg->first;

                    TGraph* graph = GraphsForSDv[GeometrySymbol][ParticleName][Source_ORG][Target_ORG];

                    graph->SetName((GeometrySymbol +", "+ParticleName+", "+Target_ORG+"<-"+Source_ORG).c_str());
                    graph->SetTitle(graph->GetName());

                    mg->Add(setGraphData(graph, hh, jj));
                    leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

                    hh++;
                }
                jj++;
            }
        }
    }
    std::cout << "\n\n-------------------------------- In One Graph ratios DoseCalcs" << CompareReferenceName <<", for all combinations, "<< QuantitiesToScore << " " << ParticleName << " ----------------------" << std::endl ;

    std::string multiGraphTitle = "DoseCalcs and "+CompareReferenceName +" differences for all geometries, particles and combinations for "+QuantitiesToScore+" ;Energy(MeV); "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle != "yes"){multiGraphTitle = "; Energy(MeV); "+ DiffExp;}
    CreateMultiGraphParametersAndCanvas(multiGraphTitle, FileName, mg, leg);
}
void G4DoseCalcsAnalysis::GenerateSourceEnegyGraph(){

    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;

    if(RegionVariableName == "Distance"){
        return;
    }

    TCanvas * Res_Ref_Canvas = new TCanvas("RegionEnergy_QuantityVsParameter", "RegionEnergy_QuantityVsParameter");
    std::string FileName;
    std::string multigraphTitle ;
    std::string Source_ORG;
    TMultiGraph *mg = new TMultiGraph();
    TLegend *leg = new TLegend();

    std::map<std::string,std::map<double,std::string>> OrgValNm;
    std::map<std::string,std::vector<double>> ValOrdered;
    std::map<std::string,std::vector<double>> ValToScor;
    std::map<std::string,std::vector<std::string>> TargetNamesToScoreOrdered;
    //std::vector<std::string,std::vector<std::string>> TargetNamesToScoreOrdered;

    for (int bb = 0 ; bb < GeometryList.size() ; bb++) {
        GeometrySymbol = GeometryList[bb];

        for ( auto Abeg = GeometryRegionVariableValue[GeometrySymbol][RegionVariableName].begin(); Abeg != GeometryRegionVariableValue[GeometrySymbol][RegionVariableName].end(); ++Abeg  ){
            //std::cout << Abeg->first << " " << GeometryRegionVariableValue[GeometrySymbol]["Mass"][Abeg->first] << std::endl;

            OrgValNm[GeometrySymbol][GeometryRegionVariableValue[GeometrySymbol][RegionVariableName][Abeg->first]]= Abeg->first;
            ValOrdered[GeometrySymbol].push_back(GeometryRegionVariableValue[GeometrySymbol][RegionVariableName][Abeg->first]);
        }

        sort(ValOrdered[GeometrySymbol].begin(), ValOrdered[GeometrySymbol].end());

        if(RegionVariableName == "Target" || RegionVariableName == "Targets"){
            for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                TargetNamesToScoreOrdered[GeometrySymbol].push_back(TargetNamesToScore[A]);
                ValToScor[GeometrySymbol].push_back(A+1);
                GeometryRegionVariableValue[GeometrySymbol][RegionVariableName][TargetNamesToScore[A]] = A+1;
                //std::cout << " mass " << ff << " " << OrgValNm[GeometrySymbol][ValOrdered[ff]] << " " << ValOrdered[GeometrySymbol][ff] << std::endl;
            }
        }
        else{
            for ( int ff = 0 ; ff < ValOrdered[GeometrySymbol].size() ; ff++ ){
                for ( int A = 0; A < TargetNamesToScore.size() ; A++  ){
                    if(TargetNamesToScore[A] == OrgValNm[GeometrySymbol][ValOrdered[GeometrySymbol][ff]] ){
                        TargetNamesToScoreOrdered[GeometrySymbol].push_back(OrgValNm[GeometrySymbol][ValOrdered[GeometrySymbol][ff]]);
                        ValToScor[GeometrySymbol].push_back(ValOrdered[GeometrySymbol][ff]);
                        //std::cout << " mass " << ff << " " << OrgValNm[GeometrySymbol][ValOrdered[ff]] << " " << ValOrdered[GeometrySymbol][ff] << std::endl;
                        break;
                    }
                }
            }
        }
    }

    int jj = 0;

    for (int bb = 0 ; bb < GeometryList.size() ; bb++) {
        GeometrySymbol = GeometryList[bb];
        for (int cc = 0 ; cc < ParticleList.size() ; cc++) {
            ParticleName = ParticleList[cc];

            for ( int A = 0; A < ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName].size() ; A++  ){

                int hh = 0;
                for ( auto Abeg = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].begin(); Abeg != ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName].end(); ++Abeg  )
                {

                    Source_ORG = Abeg->first;

                    bool IsIn = false;
                    for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                        if(SourceNamesToScore[gg] == Source_ORG){
                            IsIn = true;
                            break;
                        }
                    }
                    if(IsIn == false){
                        continue;
                    }

                    double xi[TargetNamesToScoreOrdered[GeometrySymbol].size()], yi[TargetNamesToScoreOrdered[GeometrySymbol].size()] ;
                    for ( int B = 0; B < TargetNamesToScoreOrdered[GeometrySymbol].size() ; B++  ){

                        if(RegionVariableName == "Distance"){

                            double X1 = GeometryRegionVariableValue[GeometrySymbol]["X"][Source_ORG];
                            double X2 = GeometryRegionVariableValue[GeometrySymbol]["X"][TargetNamesToScoreOrdered[GeometrySymbol][B]];
                            double X = std::abs(X1-X2);
                            double Y1 = GeometryRegionVariableValue[GeometrySymbol]["Y"][Source_ORG];
                            double Y2 = GeometryRegionVariableValue[GeometrySymbol]["Y"][TargetNamesToScoreOrdered[GeometrySymbol][B]];
                            double Y = std::abs(Y1-Y2);
                            double Z1 = GeometryRegionVariableValue[GeometrySymbol]["Z"][Source_ORG];
                            double Z2 = GeometryRegionVariableValue[GeometrySymbol]["Z"][TargetNamesToScoreOrdered[GeometrySymbol][B]];
                            double Z = std::abs(Z1-Z2);
                            GeometryRegionVariableValue[GeometrySymbol]["Distance"][TargetNamesToScoreOrdered[GeometrySymbol][B]] = std::sqrt((X*X)+(Y*Y)+(Z*Z));
                            //std::cout <<  " Source_ORG:" << Source_ORG  <<  " Target_ORG:" << Target_ORG << " X1:" << X1 <<  " X2:" << X2  <<  " Y1:" << Y1  <<  " Y2:" << Y2  << " Z1:" << Z1 <<  " Z2:" << Z2 <<  " X:" << X <<  " Y:" << Y <<  " Z:" << Z <<  " Distance:" << GeometryRegionVariableValue[GeometrySymbol]["Distance"][Target_ORG] << std::endl;

                        }

                        xi[B] = GeometryRegionVariableValue[GeometrySymbol][RegionVariableName][TargetNamesToScoreOrdered[GeometrySymbol][B]] ;
                        yi[B] = ResultTable[QuantitiesToScore][GeometrySymbol][ParticleName][Source_ORG][TargetNamesToScoreOrdered[GeometrySymbol][B]][ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][A]];

                        if(QuantityUseLog[QuantitiesToScore] == true && yi[B] < MinValForLog){
                            yi[B] = MinValForLog;
                        }
                    }

                    TGraph* gr1 = CreateGraph (TargetNamesToScoreOrdered[GeometrySymbol].size() , xi, yi); //this graph is related if we want to remove the zeros of SAF from Graph, if not we have to uncomment the graph1 below

                    std::ostringstream nm;
                    nm << GeometrySymbol<< ", " << ParticleName<< ", "<< Source_ORG << ", " <<  ResEnergies[QuantitiesToScore][GeometrySymbol][ParticleName][A] << " MeV" ;

                    std::string name = nm.str();
                    gr1->SetName(const_cast<char*>(name.c_str()));

                    gr1->SetTitle(gr1->GetName());

                    mg->Add(setGraphData(gr1,hh,jj));
                    leg->AddEntry(gr1,gr1->GetName(),"LP");  // to add the explanation of this colored line
                    hh++;
                }
                jj++;
            }
        }
    }

    std::cout << "\n\n-------------------------------- Parameter Graph For All geometries, particles and energies, "<< QuantitiesToScore << " " << ParticleName << " "<< RegionVariableNameWithUnit << " " << QuantityUnit[QuantitiesToScore] << "----------------------" << std::endl ;

    FileName = GraphsDirectoryPath+RegionVariableName+"_"+QuantitiesToScore+"_"+"_ForAllGeometriesParticleEnergiesSourcesTargets"+GraphsExt;
    multigraphTitle =  QuantitiesToScore + " versus " + RegionVariableName +" evolution in all simulations ; "+RegionVariableNameWithUnit+ ";" +QuantityUnit[QuantitiesToScore] ;
    if(PrintTitle == "yes"){
        mg->SetTitle(multigraphTitle.c_str());
    }
    else{
        multigraphTitle = " ; "+RegionVariableNameWithUnit+ ";" +QuantityUnit[QuantitiesToScore] ;
        mg->SetTitle(multigraphTitle.c_str());
    }

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }
    if(UseLogVariable == "yes"){
        gPad->SetLogx(1);
    }
    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    TText* t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);
    double y = /*0.67; */gPad->GetUymin() ;//- 0.2*mg->GetYaxis()->GetBinWidth(1);
    for (int bb = 0 ; bb < GeometryList.size() ; bb++) {
        GeometrySymbol = GeometryList[bb];
        for(int aa = 0; aa < ValToScor[GeometrySymbol].size(); aa++){
            t->DrawText(ValToScor[GeometrySymbol][aa],y,const_cast<char*>((TargetNamesToScoreOrdered[GeometrySymbol][aa]+"("+GeometrySymbol+")").c_str()));
        }
    }
    Res_Ref_Canvas->Print(FileName.c_str());

    delete leg ;
    delete mg ;
    delete Res_Ref_Canvas;
}

void G4DoseCalcsAnalysis::GenerateRadioTracerResultsInOneGraphForAllComAndSrcReg(){ // just for cross in one graph, because the self is already generated in one graph previously

    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;

    std::map<std::string,std::map<std::string,TGraph*>> AllGeometriesRadiotracersGraphs; // Result_Reference Particle Source Target Energy Value
    std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>> AllCombinationsGeometriesRadiotracerGraphs; // Result_Reference Particle Source Target Energy Value
    std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>> AllCombinationsRadioTracersGeometryGraphs; // Result_Reference Particle Source Target Energy Value
    std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>> SourcesGeometriesRadiotracerGraphs; // Result_Reference Particle Source Target Energy Value
    std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>> SourcesRadioTracersGeometryGraphs; // Result_Reference Particle Source Target Energy Value

    std::string Source_ORG;
    std::string Target_ORG;

    std::string FileName ;
    TCanvas* ResCanvas ;
    TMultiGraph *mg ;
    TLegend *leg ;
    TText* t = new TText();
    std::string multiGraphTitle = "";

    //std::cout << " The Quantity " << QuantitiesToScore << "\n" << std::endl;

    if(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].size() == 0 ){
        std::cout << " No Results data for " << QuantitiesToScore << "\n" << std::endl;
        return;
    }


    for(int aa = 0; aa < GeometryList.size(); aa++){

        GeometrySymbol = GeometryList[aa];
        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            NumOfEne = SourceNamesToScore.size()*TargetNamesToScore.size();
            double xi[NumOfEne]; double yi[NumOfEne];

            int jj = 0;
            for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
                for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                    Source_ORG = SourceNamesToScore[ee];
                    Target_ORG = TargetNamesToScore[gg];
                    xi[jj] = jj;
                    yi[jj] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];

                    if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                        yi[jj] = MinValForLog;
                    }

                    if(Source_ORG == "Liver" && Target_ORG == "Liver"){
                        //std::cout << " PARTICLE_NAME " << PARTICLE_NAME << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;
                    }
                    jj++;
                }
            }
            TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

            //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
            AllGeometriesRadiotracersGraphs[GeometrySymbol][RadioTracerName] = gr1;

        }
    }
    FileName = GraphsDirectoryPath+"Radionuclide_"+QuantitiesToScore+"_Geometries_Radiotracers_Combinations"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    int jj = 0;
    for ( auto Mbeg = AllGeometriesRadiotracersGraphs.begin(); Mbeg != AllGeometriesRadiotracersGraphs.end(); ++Mbeg  )
    {

        GeometrySymbol = Mbeg->first;

        int hh = 0;
        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {
            RadioTracerName = Nbeg->first;

            TGraph* graph = AllGeometriesRadiotracersGraphs[GeometrySymbol][RadioTracerName];

            std::string graph_label = GeometrySymbol+", "+ RadioTracerName;
            graph->SetName(graph_label.c_str());
            graph->SetTitle(graph->GetName());

            mg->Add(setGraphData(graph, hh, jj));
            leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
            hh++;
            jj++;
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " in all combinations, geometries and radiotracers ;target<-source; "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; target<-source; "+ QuantityUnit[QuantitiesToScore];
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->GetHistogram()->GetXaxis()->SetTickLength(0); // because we will print text
    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);

    double y = /*0.67; */gPad->GetUymin() ;//- 0.2*mg->GetYaxis()->GetBinWidth(1);

    int jh = 0;
    for(int aa = 0; aa < SourceNamesToScore.size(); aa++){
        for(int bb = 0; bb < TargetNamesToScore.size(); bb++){
            double z=jh;
            t->DrawText(z,y,const_cast<char*>((TargetNamesToScore[bb]+"<-"+SourceNamesToScore[aa]).c_str()));
            //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
            jh++;
        }
    }

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;


    // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>> AllGeometriesRadiotracersSourcesGraphs; // Result_Reference Particle Source Target Energy Value
    std::vector<std::string> XLabelsVec;
    std::vector<int> XLabelsIDVec;

    for(int aa = 0; aa < GeometryList.size(); aa++){

        GeometrySymbol = GeometryList[aa];
        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];
            // iterations on source name
            for ( auto Abeg = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName].begin(); Abeg != ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName].end(); ++Abeg  )
            {

                Source_ORG = Abeg->first;

                bool IsIn = false;
                for (int gg = 0 ; gg < SourceNamesToScore.size() ; gg++) {
                    if(SourceNamesToScore[gg] == Source_ORG){
                        IsIn = true;
                        break;
                    }
                }
                if(IsIn == false){
                    continue;
                }

                NumOfEne = TargetNamesToScore.size();
                double xi[NumOfEne]; double yi[NumOfEne];
                XLabelsVec.clear();
                XLabelsIDVec.clear();

                for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                    Target_ORG = TargetNamesToScore[gg];
                    xi[gg] = gg;
                    yi[gg] = Abeg->second[Target_ORG];

                    XLabelsVec.push_back(Target_ORG);
                    XLabelsIDVec.push_back(gg);

                    if(QuantityUseLog[QuantitiesToScore] == true && yi[gg] < MinValForLog){
                        yi[gg] = MinValForLog;
                    }

                    if(Source_ORG == "Liver" && Target_ORG == "Liver"){
                        //std::cout << " PARTICLE_NAME " << PARTICLE_NAME << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;
                    }
                }

                TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

                //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
                AllGeometriesRadiotracersSourcesGraphs[GeometrySymbol][RadioTracerName][Source_ORG] = gr1;
            }
        }
    }

    FileName = GraphsDirectoryPath+"Radionuclide_"+QuantitiesToScore+"_RadioTracer_Geometries_Source_Targets"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();

    jj = 0;

    for ( auto Mbeg = AllGeometriesRadiotracersSourcesGraphs.begin(); Mbeg != AllGeometriesRadiotracersSourcesGraphs.end(); ++Mbeg  )
    {

        GeometrySymbol = Mbeg->first;

        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {

            RadioTracerName = Nbeg->first;

            int hh = 0;

            for ( auto Abeg = Nbeg->second.begin(); Abeg != Nbeg->second.end(); ++Abeg  )
            {

                Source_ORG = Abeg->first;

                TGraph* graph = AllGeometriesRadiotracersSourcesGraphs[GeometrySymbol][RadioTracerName][Source_ORG];

                std::string graph_label = GeometrySymbol+", "+ RadioTracerName+", "+ Source_ORG;
                graph->SetName(graph_label.c_str());
                graph->SetTitle(graph->GetName());

                mg->Add(setGraphData(graph, hh, jj));
                leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

                hh++;
            }
            jj++;
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }

    gPad->SetLogx(0);


    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " For each geometries and radiotracers, sources ;Targets; "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; Targets; "+ QuantityUnit[QuantitiesToScore];
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);
    y = /*0.67; */gPad->GetUymin() ;//- 0.2*mg->GetYaxis()->GetBinWidth(1);
    for(int aa = 0; aa < TargetNamesToScore.size(); aa++){
        double z=aa;
        t->DrawText(z,y,const_cast<char*>(TargetNamesToScore[aa].c_str()));
        //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
    }

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;


    // //////////////////////////////////////////////////////////////////////////////////////////////

    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
            Source_ORG = SourceNamesToScore[ee];
            Target_ORG = TargetNamesToScore[gg];

            for(int aa = 0; aa < GeometryList.size(); aa++){

                GeometrySymbol = GeometryList[aa];

                NumOfEne = RadiotracerList.size();
                double xi[NumOfEne]; double yi[NumOfEne];

                int jj = 0;

                for(int bb = 0; bb < RadiotracerList.size(); bb++){

                    RadioTracerName = RadiotracerList[bb];

                    xi[jj] = jj;
                    yi[jj] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];

                    if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                        yi[jj] = MinValForLog;
                    }

                    jj++;
                }

                TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

                //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
                AllCombinationsGeometriesRadiotracerGraphs[Source_ORG][Target_ORG][GeometrySymbol] = gr1;
            }
        }
    }
    FileName = GraphsDirectoryPath+"Radionuclide_"+QuantitiesToScore+"_Combinations_Geometries_Radiotracer"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
            Source_ORG = SourceNamesToScore[ee];
            Target_ORG = TargetNamesToScore[gg];
            int hh = 0;

            for(int aa = 0; aa < GeometryList.size(); aa++){

                GeometrySymbol = GeometryList[aa];

                TGraph* graph = AllCombinationsGeometriesRadiotracerGraphs[Source_ORG][Target_ORG][GeometrySymbol];

                std::string graph_label = GeometrySymbol+", "+ Target_ORG + "<-" +Source_ORG;
                graph->SetName(graph_label.c_str());
                graph->SetTitle(graph->GetName());

                mg->Add(setGraphData(graph, hh, jj));
                leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
                hh++;
            }
            jj++;
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for each combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " in all combinations, geometries and radiotracers ;Radionuclides; "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; Radionuclides; "+ QuantityUnit[QuantitiesToScore];
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->GetXaxis()->SetTickLength(0); // because we will print text
    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);
    y = /*0.67; */gPad->GetUymin() ;//- 0.2*mg->GetYaxis()->GetBinWidth(1);
    jh = 0;
    for(int bb = 0; bb < RadiotracerList.size(); bb++){
        double z=jh;
        t->DrawText(z,y,const_cast<char*>((RadiotracerList[bb]).c_str()));
        //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
        jh++;
    }


    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;

    // //////////////////////////////////////////////////////////////////////////////////////////////

    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
            Source_ORG = SourceNamesToScore[ee];
            Target_ORG = TargetNamesToScore[gg];

            for(int bb = 0; bb < RadiotracerList.size(); bb++){

                RadioTracerName = RadiotracerList[bb];

                NumOfEne = GeometryList.size();
                double xi[NumOfEne]; double yi[NumOfEne];
                int jj = 0;
                for(int aa = 0; aa < GeometryList.size(); aa++){

                    GeometrySymbol = GeometryList[aa];

                    xi[jj] = jj;
                    yi[jj] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];

                    if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                        yi[jj] = MinValForLog;
                    }

                    jj++;
                }
                TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
                AllCombinationsRadioTracersGeometryGraphs[Source_ORG][Target_ORG][RadioTracerName] = gr1;
            }
        }
    }
    FileName = GraphsDirectoryPath+"Radionuclide_"+QuantitiesToScore+"_Combinations_Radiotracers_Geometry"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
            Source_ORG = SourceNamesToScore[ee];
            Target_ORG = TargetNamesToScore[gg];
            int hh = 0;

            for(int bb = 0; bb < RadiotracerList.size(); bb++){

                RadioTracerName = RadiotracerList[bb];

                TGraph* graph = AllCombinationsRadioTracersGeometryGraphs[Source_ORG][Target_ORG][RadioTracerName];

                std::string graph_label = RadioTracerName+", "+ Target_ORG + "<-" +Source_ORG;
                graph->SetName(graph_label.c_str());
                graph->SetTitle(graph->GetName());

                mg->Add(setGraphData(graph, hh, jj));
                leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
                hh++;
            }
            jj++;
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " in all combinations, geometries and radiotracers ;Geometries; "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; Geometries; "+ QuantityUnit[QuantitiesToScore];
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->GetXaxis()->SetTickLength(0); // because we will print text
    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);
    y = /*0.67; */gPad->GetUymin() ;//- 0.2*mg->GetYaxis()->GetBinWidth(1);
    jh = 0;
    for(int bb = 0; bb < GeometryList.size(); bb++){
        double z=jh;
        t->DrawText(z,y,const_cast<char*>((GeometryList[bb]).c_str()));
        //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
        jh++;
    }

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;


    // /////////////////////////////-------------For just sources------------- ///////////////////////////////////////////

    FileName = GraphsDirectoryPath+"Radionuclide_"+QuantitiesToScore+"_SourceRegions_Radiotracers_Geometry"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;

    if(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].size() == 0 ){
        std::cout << " No Results data for " << QuantitiesToScore << "\n" << std::endl;
        return;
    }


    for(int aa = 0; aa < GeometryList.size(); aa++){

        GeometrySymbol = GeometryList[aa];
        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            NumOfEne = SourceNamesToScore.size();
            double xi[NumOfEne]; double yi[NumOfEne];

            int jj = 0;
            for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
                Source_ORG = SourceNamesToScore[ee];
                xi[jj] = jj;
                yi[jj] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Source_ORG];

                if(QuantityUseLog[QuantitiesToScore] == true && (yi[jj] < MinValForLog || yi[jj] == 0)){
                    yi[jj] = MinValForLog;
                }

                jj++;
            }
            TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

            //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
            AllGeometriesRadiotracersGraphs[GeometrySymbol][RadioTracerName] = gr1;

        }
    }


    FileName = GraphsDirectoryPath+"Radionuclide_"+QuantitiesToScore+"_Geometries_Radiotracers_SourceRegions"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for ( auto Mbeg = AllGeometriesRadiotracersGraphs.begin(); Mbeg != AllGeometriesRadiotracersGraphs.end(); ++Mbeg  )
    {

        GeometrySymbol = Mbeg->first;

        int hh = 0;
        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {
            RadioTracerName = Nbeg->first;

            TGraph* graph = AllGeometriesRadiotracersGraphs[GeometrySymbol][RadioTracerName];

            std::string graph_label = GeometrySymbol+", "+ RadioTracerName;
            graph->SetName(graph_label.c_str());
            graph->SetTitle(graph->GetName());

            mg->Add(setGraphData(graph, hh, jj));
            leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
            hh++;
            jj++;
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " in all combinations, geometries and radiotracers ;Source regions; "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; target<-source; "+ QuantityUnit[QuantitiesToScore];
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->GetXaxis()->SetTickLength(0); // because we will print text
    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);
    y = /*0.67; */gPad->GetUymin() ;//- 0.2*mg->GetYaxis()->GetBinWidth(1);

    jh = 0;
    for(int aa = 0; aa < SourceNamesToScore.size(); aa++){
        double z=jh;
        t->DrawText(z,y,const_cast<char*>((SourceNamesToScore[aa]).c_str()));
        //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
        jh++;
    }

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;

    // //////////////////////////////////////////////////////////////////////////////////////////////

    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        Source_ORG = SourceNamesToScore[ee];

        for(int aa = 0; aa < GeometryList.size(); aa++){

            GeometrySymbol = GeometryList[aa];

            NumOfEne = RadiotracerList.size();
            double xi[NumOfEne]; double yi[NumOfEne];

            int jj = 0;

            for(int bb = 0; bb < RadiotracerList.size(); bb++){

                RadioTracerName = RadiotracerList[bb];

                xi[jj] = jj;
                yi[jj] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Source_ORG];

                if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                    yi[jj] = MinValForLog;
                }

                jj++;
            }

            TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

            //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
            AllCombinationsGeometriesRadiotracerGraphs[Source_ORG][Source_ORG][GeometrySymbol] = gr1;
        }
    }
    FileName = GraphsDirectoryPath+"Radionuclide_"+QuantitiesToScore+"_SourceRegions_Geometries_Radiotracer"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        Source_ORG = SourceNamesToScore[ee];
        int hh = 0;

        for(int aa = 0; aa < GeometryList.size(); aa++){

            GeometrySymbol = GeometryList[aa];

            TGraph* graph = AllCombinationsGeometriesRadiotracerGraphs[Source_ORG][Source_ORG][GeometrySymbol];

            std::string graph_label = GeometrySymbol+", "+Source_ORG;
            graph->SetName(graph_label.c_str());
            graph->SetTitle(graph->GetName());

            mg->Add(setGraphData(graph, hh, jj));
            leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
            hh++;
        }
        jj++;
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for each combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " For source regions, geometries and radiotracers ;Radionuclides; "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; Radionuclides; "+ QuantityUnit[QuantitiesToScore];
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->GetXaxis()->SetTickLength(0); // because we will print text
    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);
    y = /*0.67; */gPad->GetUymin() ;//- 0.2*mg->GetYaxis()->GetBinWidth(1);
    jh = 0;
    for(int bb = 0; bb < RadiotracerList.size(); bb++){
        double z=jh;
        t->DrawText(z,y,const_cast<char*>((RadiotracerList[bb]).c_str()));
        //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
        jh++;
    }


    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;

    // //////////////////////////////////////////////////////////////////////////////////////////////

    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        Source_ORG = SourceNamesToScore[ee];

        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            NumOfEne = GeometryList.size();
            double xi[NumOfEne]; double yi[NumOfEne];
            int jj = 0;
            for(int aa = 0; aa < GeometryList.size(); aa++){

                GeometrySymbol = GeometryList[aa];

                xi[jj] = jj;
                yi[jj] = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Source_ORG];

                if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                    yi[jj] = MinValForLog;
                }

                jj++;
            }
            TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
            AllCombinationsRadioTracersGeometryGraphs[Source_ORG][Source_ORG][RadioTracerName] = gr1;
        }
    }
    FileName = GraphsDirectoryPath+"Radionuclide_"+QuantitiesToScore+"_SourceRegions_Radiotracers_Geometry"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        Source_ORG = SourceNamesToScore[ee];
        int hh = 0;

        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            TGraph* graph = AllCombinationsRadioTracersGeometryGraphs[Source_ORG][Source_ORG][RadioTracerName];

            std::string graph_label = RadioTracerName+", "+Source_ORG;
            graph->SetName(graph_label.c_str());
            graph->SetTitle(graph->GetName());

            mg->Add(setGraphData(graph, hh, jj));
            leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
            hh++;
        }
        jj++;
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " For source regions, geometries and radiotracers ;Geometries; "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; Geometries; "+ QuantityUnit[QuantitiesToScore];
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->GetXaxis()->SetTickLength(0); // because we will print text
    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);
    y = /*0.67; */gPad->GetUymin() ;//- 0.2*mg->GetYaxis()->GetBinWidth(1);
    jh = 0;
    for(int bb = 0; bb < GeometryList.size(); bb++){
        double z=jh;
        t->DrawText(z,y,const_cast<char*>((GeometryList[bb]).c_str()));
        //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
        jh++;
    }

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;

}
void G4DoseCalcsAnalysis::GenerateRadioTracerResultsInOneGraphWithComparisonForAllComAndSrcReg(){ // just for cross in one graph, because the self is already generated in one graph previously

    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;

    std::map<std::string,std::map<std::string,TGraph*>> AllGeometriesRadiotracersGraphs; // Result_Reference Particle Source Target Energy Value
    std::map<std::string,std::map<std::string,TGraph*>> AllGeometriesRadiotracersGraphsReference; // Result_Reference Particle Source Target Energy Value

    std::string Source_ORG;
    std::string Target_ORG;

    std::vector<std::string> XLabelsVec;
    std::vector<int> XLabelsIDVec;

    //std::cout << " The Quantity " << QuantitiesToScore << "\n" << std::endl;

    if(ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].size() == 0 ){
        std::cout << " No Results data for " << QuantitiesToScore << "\n" << std::endl;
        return;
    }

    for(int aa = 0; aa < GeometryList.size(); aa++){

        GeometrySymbol = GeometryList[aa];
        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            NumOfEne = SourceNamesToScore.size()*TargetNamesToScore.size();
            double xi[NumOfEne]; double yi[NumOfEne];
            double ai[NumOfEne];
            XLabelsVec.clear();
            XLabelsIDVec.clear();

            int jj = 0;
            for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
                for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                    Source_ORG = SourceNamesToScore[ee];
                    Target_ORG = TargetNamesToScore[gg];
                    xi[jj] = jj;
                    yi[jj] = ResultQuantityGeometryRadioTracerSourceTargetValues   [QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    ai[jj] = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];

                    if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                        yi[jj] = MinValForLog;
                    }
                    if(QuantityUseLog[QuantitiesToScore] == true && ai[jj] < MinValForLog){
                        ai[jj] = MinValForLog;
                    }

                    //std::cout << " QuantitiesToScore " << QuantitiesToScore << " GeometrySymbol " << GeometrySymbol << " RadioTracerName " << RadioTracerName << " Source_ORG " << Source_ORG  << " Target_ORG " << Target_ORG << " " << xi[jj] << " " << yi[jj]  << " " << ai[jj] << std::endl ;

                    jj++;
                }
            }
            TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
            TGraph* gr2 = CreateGraph (NumOfEne, xi, ai);

            AllGeometriesRadiotracersGraphs[GeometrySymbol][RadioTracerName] = gr1;
            AllGeometriesRadiotracersGraphsReference[GeometrySymbol][RadioTracerName] = gr2;

        }
    }

    std::string FileName ;
    TCanvas* ResCanvas ;
    TMultiGraph *mg ;
    TLegend *leg ;

    FileName = GraphsDirectoryPath+"Radionuclide_"+QuantitiesToScore+"_Geometries_Radiotracers_Combinations_And_Reference"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();

    int jj = 0;
    int hh = 0;
    // generate just for Cross-irradiation
    for ( auto Mbeg = AllGeometriesRadiotracersGraphs.begin(); Mbeg != AllGeometriesRadiotracersGraphs.end(); ++Mbeg  )
    {

        GeometrySymbol = Mbeg->first;

        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {

            RadioTracerName = Nbeg->first;

            TGraph* graph = AllGeometriesRadiotracersGraphs[GeometrySymbol][RadioTracerName];
            graph->SetName((GeometrySymbol+", "+ RadioTracerName).c_str());
            graph->SetTitle(graph->GetName());
            mg->Add(setGraphData(graph, hh, jj));
            leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

            jj++;

            TGraph* graph1 = AllGeometriesRadiotracersGraphsReference[GeometrySymbol][RadioTracerName];
            graph1->SetName((CompareReferenceName+", "+ GeometrySymbol+", "+ RadioTracerName).c_str());
            graph1->SetTitle(graph1->GetName());
            mg->Add(setGraphData(graph1, hh, jj));
            leg->AddEntry(graph1,graph1->GetName(),"LP");  // to add the explanation of this colored line

            jj++;
            hh++;

        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    std::string multiGraphTitle = " In all combinations, geometries and radiotracers ;target<-source; "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; target<-source; "+ QuantityUnit[QuantitiesToScore];
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    TText* t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);
    double y = /*0.67; */gPad->GetUymin() ;//- 0.2*mg->GetYaxis()->GetBinWidth(1);
    int jh = 0;
    for(int aa = 0; aa < SourceNamesToScore.size(); aa++){
        for(int bb = 0; bb < TargetNamesToScore.size(); bb++){
            double z=jh;
            t->DrawText(z,y,const_cast<char*>((TargetNamesToScore[bb]+"<-"+SourceNamesToScore[aa]).c_str()));
            //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
            jh++;
        }
    }


    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;



    // ///////////////////////////------------ Source regions ---------/////////////////////////////////////////////////////////////


    for(int aa = 0; aa < GeometryList.size(); aa++){

        GeometrySymbol = GeometryList[aa];
        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            NumOfEne = SourceNamesToScore.size();
            double xi[NumOfEne]; double yi[NumOfEne];
            double ai[NumOfEne];
            XLabelsVec.clear();
            XLabelsIDVec.clear();

            int jj = 0;
            for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
                Source_ORG = SourceNamesToScore[ee];
                xi[jj] = jj;
                yi[jj] = ResultQuantityGeometryRadioTracerSourceTargetValues   [QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Source_ORG];
                ai[jj] = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Source_ORG];

                if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                    yi[jj] = MinValForLog;
                }
                if(QuantityUseLog[QuantitiesToScore] == true && ai[jj] < MinValForLog){
                    ai[jj] = MinValForLog;
                }

                //std::cout << " QuantitiesToScore " << QuantitiesToScore << " GeometrySymbol " << GeometrySymbol << " RadioTracerName " << RadioTracerName << " Source_ORG " << Source_ORG  << " Source_ORG " << Source_ORG << " " << xi[jj] << " " << yi[jj]  << " " << ai[jj] << std::endl ;

                jj++;
            }
            TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
            TGraph* gr2 = CreateGraph (NumOfEne, xi, ai);

            AllGeometriesRadiotracersGraphs[GeometrySymbol][RadioTracerName] = gr1;
            AllGeometriesRadiotracersGraphsReference[GeometrySymbol][RadioTracerName] = gr2;

        }
    }

    FileName = GraphsDirectoryPath+"Radionuclide_"+QuantitiesToScore+"_Geometries_Radiotracers_SourceRegions_And_Reference"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();

    jj = 0;
    hh = 0;
    // generate just for Cross-irradiation
    for ( auto Mbeg = AllGeometriesRadiotracersGraphs.begin(); Mbeg != AllGeometriesRadiotracersGraphs.end(); ++Mbeg  )
    {

        GeometrySymbol = Mbeg->first;

        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {

            RadioTracerName = Nbeg->first;

            TGraph* graph = AllGeometriesRadiotracersGraphs[GeometrySymbol][RadioTracerName];
            graph->SetName((GeometrySymbol+", "+ RadioTracerName).c_str());
            graph->SetTitle(graph->GetName());
            mg->Add(setGraphData(graph, hh, jj));
            leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

            jj++;

            TGraph* graph1 = AllGeometriesRadiotracersGraphsReference[GeometrySymbol][RadioTracerName];
            graph1->SetName((CompareReferenceName+", "+ GeometrySymbol+", "+ RadioTracerName).c_str());
            graph1->SetTitle(graph1->GetName());
            mg->Add(setGraphData(graph1, hh, jj));
            leg->AddEntry(graph1,graph1->GetName(),"LP");  // to add the explanation of this colored line

            jj++;
            hh++;

        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogy(1);
    }

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " In all source regions, geometries and radiotracers ;Source regions; "+ QuantityUnit[QuantitiesToScore];
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; target<-source; "+ QuantityUnit[QuantitiesToScore];
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);
    y = /*0.67; */gPad->GetUymin() ;//- 0.2*mg->GetYaxis()->GetBinWidth(1);
    jh = 0;
    for(int aa = 0; aa < SourceNamesToScore.size(); aa++){
        double z=jh;
        t->DrawText(z,y,const_cast<char*>((SourceNamesToScore[aa]).c_str()));
        //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
        jh++;
    }


    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;




}
void G4DoseCalcsAnalysis::GenerateRadioTracerComparisonFactorGraphsForAllComAndSrcReg(){ // just for cross in one graph, because the self is already generated in one graph previously

    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;

    std::map<std::string,std::map<std::string,TGraph*>> AllGeometriesRadiotracersGraphs; // Result_Reference Particle Source Target Energy Value
    std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>> AllCombinationsGeometriesRadiotracerGraphs; // Result_Reference Particle Source Target Energy Value
    std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>> AllCombinationsRadioTracersGeometryGraphs; // Result_Reference Particle Source Target Energy Value

    std::string Source_ORG;
    std::string Target_ORG;

    std::string FileName ;
    TCanvas* ResCanvas ;
    TMultiGraph *mg ;
    TLegend *leg ;
    TText* t = new TText();
    std::string multiGraphTitle = "";

    //std::cout << " The Quantity " << QuantitiesToScore << "\n" << std::endl;

    if(ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore].size() == 0 ){
        std::cout << " No Results data for " << QuantitiesToScore << "\n" << std::endl;
        return;
    }


    for(int aa = 0; aa < GeometryList.size(); aa++){

        GeometrySymbol = GeometryList[aa];
        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            NumOfEne = SourceNamesToScore.size()*TargetNamesToScore.size();
            double xi[NumOfEne]; double yi[NumOfEne];

            int jj = 0;
            for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
                for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                    Source_ORG = SourceNamesToScore[ee];
                    Target_ORG = TargetNamesToScore[gg];
                    xi[jj] = jj;

                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    double a3 = RelativeDifferenceCalculation(a1,a2);
                    if(a3==NULL || a3==0.){a3 = MinValForLog;}
                    yi[jj] = a3;

                    if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                        yi[jj] = MinValForLog;
                    }

                    if(Source_ORG == "Liver" && Target_ORG == "Liver"){
                        //std::cout << " PARTICLE_NAME " << PARTICLE_NAME << " Source_ORG " << Source_ORG << " Target_ORG " << Target_ORG << " Ene " << xi[jh] << " Val " << yi[jh] << std::endl ;
                    }
                    jj++;
                }
            }
            TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

            //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
            AllGeometriesRadiotracersGraphs[GeometrySymbol][RadioTracerName] = gr1;
        }
    }
    FileName = GraphsDirectoryPath+"Radionuclide_"+DifferenceMethod+"_"+QuantitiesToScore+"_Geometries_Radiotracers_Combinations"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    int jj = 0;
    for ( auto Mbeg = AllGeometriesRadiotracersGraphs.begin(); Mbeg != AllGeometriesRadiotracersGraphs.end(); ++Mbeg  )
    {

        GeometrySymbol = Mbeg->first;

        int hh = 0;
        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {

            RadioTracerName = Nbeg->first;


            TGraph* graph = AllGeometriesRadiotracersGraphs[GeometrySymbol][RadioTracerName];

            std::string graph_label = GeometrySymbol+", "+ RadioTracerName;
            graph->SetName(graph_label.c_str());
            graph->SetTitle(graph->GetName());

            mg->Add(setGraphData(graph, hh, jj));
            leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
            hh++;
            jj++;
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(DifferenceMethod == "RA"){gPad->SetLogy(0);}else{gPad->SetLogy(1);}

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " in all combinations, geometries and radiotracers ;target<-source; "+ DiffExp;
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; target<-source; "+ DiffExp;
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);

    double y ;
    //if(DifferenceMethod == "RA"){gPad->GetUymin() - 0.2*mg->GetYaxis()->GetBinWidth(1)}else{y = gPad->GetUymin();}

    y = 0;//gPad->GetUymin() - 0.2*mg->GetYaxis()->GetBinWidth(1);
    int jh = 0;
    for(int aa = 0; aa < SourceNamesToScore.size(); aa++){
        for(int bb = 0; bb < TargetNamesToScore.size(); bb++){
            double z=jh;
            //t->SetX(z);
            //t->SetY(y);
            t->DrawText(z,y,const_cast<char*>((TargetNamesToScore[bb]+"<-"+SourceNamesToScore[aa]).c_str()));
            //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
            jh++;
        }
    }

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;

    // //////////////////////////////////////////////////////////////////////////////////////////////

    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
            Source_ORG = SourceNamesToScore[ee];
            Target_ORG = TargetNamesToScore[gg];

            for(int aa = 0; aa < GeometryList.size(); aa++){

                GeometrySymbol = GeometryList[aa];

                NumOfEne = RadiotracerList.size();
                double xi[NumOfEne]; double yi[NumOfEne];

                int jj = 0;

                for(int bb = 0; bb < RadiotracerList.size(); bb++){

                    RadioTracerName = RadiotracerList[bb];

                    xi[jj] = jj;
                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    double a3 = RelativeDifferenceCalculation(a1,a2);
                    if(a3==NULL || a3==0.){a3 = MinValForLog;}
                    yi[jj] = a3;

                    if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                        yi[jj] = MinValForLog;
                    }

                    jj++;
                }

                TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

                //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
                AllCombinationsGeometriesRadiotracerGraphs[Source_ORG][Target_ORG][GeometrySymbol] = gr1;
            }
        }
    }
    FileName = GraphsDirectoryPath+"Radionuclide_"+DifferenceMethod+"_"+QuantitiesToScore+"_Combinations_Geometries_Radiotracer"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
            Source_ORG = SourceNamesToScore[ee];
            Target_ORG = TargetNamesToScore[gg];
            int hh = 0;

            for(int aa = 0; aa < GeometryList.size(); aa++){

                GeometrySymbol = GeometryList[aa];

                TGraph* graph = AllCombinationsGeometriesRadiotracerGraphs[Source_ORG][Target_ORG][GeometrySymbol];

                std::string graph_label = GeometrySymbol+", "+ Target_ORG + "<-" +Source_ORG;
                graph->SetName(graph_label.c_str());
                graph->SetTitle(graph->GetName());

                mg->Add(setGraphData(graph, hh, jj));
                leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
                hh++;
                jj++;
            }
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for each combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(DifferenceMethod == "RA"){gPad->SetLogy(0);}else{gPad->SetLogy(1);}

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " in all combinations, geometries and radiotracers ;Radionuclides; "+ DiffExp;
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; Radionuclides; "+ DiffExp;
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);

    y = 0;//gPad->GetUymin() - 0.2*mg->GetYaxis()->GetBinWidth(1);
    jh = 0;
    for(int bb = 0; bb < RadiotracerList.size(); bb++){
        double z=jh;
        t->DrawText(z,y,const_cast<char*>((RadiotracerList[bb]).c_str()));
        //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
        jh++;
    }


    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;

    // //////////////////////////////////////////////////////////////////////////////////////////////

    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
            Source_ORG = SourceNamesToScore[ee];
            Target_ORG = TargetNamesToScore[gg];

            for(int bb = 0; bb < RadiotracerList.size(); bb++){

                RadioTracerName = RadiotracerList[bb];

                NumOfEne = GeometryList.size();
                double xi[NumOfEne]; double yi[NumOfEne];
                int jj = 0;
                for(int aa = 0; aa < GeometryList.size(); aa++){

                    GeometrySymbol = GeometryList[aa];

                    xi[jj] = jj;
                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    double a3 = RelativeDifferenceCalculation(a1,a2);
                    if(a3==NULL || a3==0.){a3 = MinValForLog;}
                    yi[jj] = a3;

                    if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                        yi[jj] = MinValForLog;
                    }

                    jj++;
                }
                TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
                AllCombinationsRadioTracersGeometryGraphs[Source_ORG][Target_ORG][RadioTracerName] = gr1;
            }
        }
    }
    FileName = GraphsDirectoryPath+"Radionuclide_"+DifferenceMethod+"_"+QuantitiesToScore+"_Combinations_Radiotracers_Geometry"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
            Source_ORG = SourceNamesToScore[ee];
            Target_ORG = TargetNamesToScore[gg];
            int hh = 0;

            for(int bb = 0; bb < RadiotracerList.size(); bb++){

                RadioTracerName = RadiotracerList[bb];

                TGraph* graph = AllCombinationsRadioTracersGeometryGraphs[Source_ORG][Target_ORG][RadioTracerName];

                std::string graph_label = RadioTracerName+", "+ Target_ORG + "<-" +Source_ORG;
                graph->SetName(graph_label.c_str());
                graph->SetTitle(graph->GetName());

                mg->Add(setGraphData(graph, hh, jj));
                leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
                hh++;
                jj++;
            }
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(DifferenceMethod == "RA"){gPad->SetLogy(0);}else{gPad->SetLogy(1);}

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " in all combinations, geometries and radiotracers ;Geometries; "+ DiffExp;
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; Geometries; "+ DiffExp;
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);

    y = 0;//gPad->GetUymin() - 0.2*mg->GetYaxis()->GetBinWidth(1);
    jh = 0;
    for(int bb = 0; bb < GeometryList.size(); bb++){
        double z=jh;
        t->DrawText(z,y,const_cast<char*>((GeometryList[bb]).c_str()));
        //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
        jh++;
    }

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;


    // ////////////////////// ----------- Source regions ------------//////////////////////////////////////////


    for(int aa = 0; aa < GeometryList.size(); aa++){

        GeometrySymbol = GeometryList[aa];
        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            NumOfEne = SourceNamesToScore.size();
            double xi[NumOfEne]; double yi[NumOfEne];

            int jj = 0;
            for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
                Source_ORG = SourceNamesToScore[ee];
                xi[jj] = jj;

                double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Source_ORG];
                double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Source_ORG];
                double a3 = RelativeDifferenceCalculation(a1,a2);
                if(a3==NULL || a3==0.){a3 = MinValForLog;}
                yi[jj] = a3;

                if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                    yi[jj] = MinValForLog;
                }

                jj++;

            }
            TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

            //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
            AllGeometriesRadiotracersGraphs[GeometrySymbol][RadioTracerName] = gr1;

        }
    }
    FileName = GraphsDirectoryPath+"Radionuclide_"+DifferenceMethod+"_"+QuantitiesToScore+"_Geometries_Radiotracers_SourceRegions"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for ( auto Mbeg = AllGeometriesRadiotracersGraphs.begin(); Mbeg != AllGeometriesRadiotracersGraphs.end(); ++Mbeg  )
    {

        GeometrySymbol = Mbeg->first;

        int hh = 0;
        for ( auto Nbeg = Mbeg->second.begin(); Nbeg != Mbeg->second.end(); ++Nbeg  )
        {

            RadioTracerName = Nbeg->first;


            TGraph* graph = AllGeometriesRadiotracersGraphs[GeometrySymbol][RadioTracerName];

            std::string graph_label = GeometrySymbol+", "+ RadioTracerName;
            graph->SetName(graph_label.c_str());
            graph->SetTitle(graph->GetName());

            mg->Add(setGraphData(graph, hh, jj));
            leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
            hh++;
            jj++;
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(DifferenceMethod == "RA"){gPad->SetLogy(0);}else{gPad->SetLogy(1);}

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " in all source regions, geometries and radiotracers ;Source region; "+ DiffExp;
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; Source region; "+ DiffExp;
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);

    y = 0.80;//gPad->GetUymin() - 0.2*mg->GetYaxis()->GetBinWidth(1);
    jh = 0;
    for(int aa = 0; aa < SourceNamesToScore.size(); aa++){
        double z=jh;
        //t->SetX(z);
        //t->SetY(y);
        t->DrawText(z,y,const_cast<char*>((SourceNamesToScore[aa]).c_str()));
        //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
        jh++;
    }

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;

    // //////////////////////////////////////////////////////////////////////////////////////////////

    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        Source_ORG = SourceNamesToScore[ee];

        for(int aa = 0; aa < GeometryList.size(); aa++){

            GeometrySymbol = GeometryList[aa];

            NumOfEne = RadiotracerList.size();
            double xi[NumOfEne]; double yi[NumOfEne];

            int jj = 0;

            for(int bb = 0; bb < RadiotracerList.size(); bb++){

                RadioTracerName = RadiotracerList[bb];

                xi[jj] = jj;
                double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Source_ORG];
                double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Source_ORG];
                double a3 = RelativeDifferenceCalculation(a1,a2);
                if(a3==NULL || a3==0.){a3 = MinValForLog;}
                yi[jj] = a3;

                //if(a1 > 2 || a1 < 0.01){a1 = a2;}

                if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                    yi[jj] = MinValForLog;
                }

                jj++;
            }

            TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);

            //TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
            AllCombinationsGeometriesRadiotracerGraphs[Source_ORG][Source_ORG][GeometrySymbol] = gr1;
        }
    }
    FileName = GraphsDirectoryPath+"Radionuclide_"+DifferenceMethod+"_"+QuantitiesToScore+"_SourceRegions_Geometries_Radiotracer"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        Source_ORG = SourceNamesToScore[ee];
        int hh = 0;

        for(int aa = 0; aa < GeometryList.size(); aa++){

            GeometrySymbol = GeometryList[aa];

            TGraph* graph = AllCombinationsGeometriesRadiotracerGraphs[Source_ORG][Source_ORG][GeometrySymbol];

            std::string graph_label = GeometrySymbol+", "+ Source_ORG + "<-" +Source_ORG;
            graph->SetName(graph_label.c_str());
            graph->SetTitle(graph->GetName());

            mg->Add(setGraphData(graph, hh, jj));
            leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
            hh++;
            jj++;
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for each combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(DifferenceMethod == "RA"){gPad->SetLogy(0);}else{gPad->SetLogy(1);}

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " in all source regions, geometries and radiotracers ;Radionuclides; "+ DiffExp;
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; Radionuclides; "+ DiffExp;
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);

    y = 0.80;//gPad->GetUymin() - 0.2*mg->GetYaxis()->GetBinWidth(1);
    jh = 0;
    for(int bb = 0; bb < RadiotracerList.size(); bb++){
        double z=jh;
        t->DrawText(z,y,const_cast<char*>((RadiotracerList[bb]).c_str()));
        //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
        jh++;
    }


    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;

    // //////////////////////////////////////////////////////////////////////////////////////////////

    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        Source_ORG = SourceNamesToScore[ee];

        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            NumOfEne = GeometryList.size();
            double xi[NumOfEne]; double yi[NumOfEne];
            int jj = 0;
            for(int aa = 0; aa < GeometryList.size(); aa++){

                GeometrySymbol = GeometryList[aa];

                xi[jj] = jj;
                double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Source_ORG];
                double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Source_ORG];
                double a3 = RelativeDifferenceCalculation(a1,a2);
                if(a3==NULL || a3==0.){a3 = MinValForLog;}
                yi[jj] = a3;

                if(QuantityUseLog[QuantitiesToScore] == true && yi[jj] < MinValForLog){
                    yi[jj] = MinValForLog;
                }

                jj++;
            }
            TGraph* gr1 = CreateGraph (NumOfEne, xi, yi);
            AllCombinationsRadioTracersGeometryGraphs[Source_ORG][Source_ORG][RadioTracerName] = gr1;
        }
    }
    FileName = GraphsDirectoryPath+"Radionuclide_"+DifferenceMethod+"_"+QuantitiesToScore+"_SourceRegions_Radiotracers_Geometry"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
        Source_ORG = SourceNamesToScore[ee];
        int hh = 0;

        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            TGraph* graph = AllCombinationsRadioTracersGeometryGraphs[Source_ORG][Source_ORG][RadioTracerName];

            std::string graph_label = RadioTracerName+", "+Source_ORG;
            graph->SetName(graph_label.c_str());
            graph->SetTitle(graph->GetName());

            mg->Add(setGraphData(graph, hh, jj));
            leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line
            hh++;
            jj++;
        }
    }

    std::cout << "\n\n-------------------------------- In One Graph Result, for all combinations, geometries and radiotracers  "<< QuantitiesToScore << " ----------------------" << std::endl ;

    if(DifferenceMethod == "RA"){gPad->SetLogy(0);}else{gPad->SetLogy(1);}

    gPad->SetLogx(0);

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = " in all combinations, geometries and radiotracers ;Geometries; "+ DiffExp;
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = "; Geometries; "+ DiffExp;
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("ALP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    t = new TText();
    t->SetTextAlign(33);
    t->SetTextAngle(60);
    t->SetTextSize(0.02);
    t->SetTextFont(72);

    y = 0;//gPad->GetUymin() - 0.2*mg->GetYaxis()->GetBinWidth(1);
    jh = 0;
    for(int bb = 0; bb < GeometryList.size(); bb++){
        double z=jh;
        t->DrawText(z,y,const_cast<char*>((GeometryList[bb]).c_str()));
        //std::cout << NumOfTargets << " " << aa << " --> " << x[aa] << " *** " << Source_ORG << " --> ("<< z[aa] << " - " << c[aa] << ") -values: " << y[aa] << "-" << b[aa] << " " << w[aa] <<std::endl ;
        jh++;
    }

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;



















    std::cout << "\n\n-------------------------------- In One Graph Result-reference Scatter plot, for all combinations and radiotracers and geometries  "<< QuantitiesToScore << " ----------------------" << std::endl ;
    FileName = GraphsDirectoryPath+"Radionuclide_ScatterPlot_"+QuantitiesToScore+"_Combinations_Radionuclides_Geometries"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for(int aa = 0; aa < GeometryList.size(); aa++){

        GeometrySymbol = GeometryList[aa];
        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            std::vector<double> XSValues;
            std::vector<double> YSValues;
            for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {

                for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                    Source_ORG = SourceNamesToScore[ee];
                    Target_ORG = TargetNamesToScore[gg];

                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    if(     a1==NULL || a1==0. || a1 == MinValForLog || isnanl(a1) || isinfl(a1) ||
                            a2==NULL || a2==0. || a2 == MinValForLog || isnanl(a2) || isinfl(a2)){
                        continue;
                    }

                    //if(a1 > 2 || a1 < 0.01){a1 = a2;}

                    XSValues.push_back(a2);
                    YSValues.push_back(a1);
                }
            }

            if(XSValues.size() == 0){
                continue;
            }

            double xi[XSValues.size()]; double yi[XSValues.size()];
            for (int ee = 0 ; ee < XSValues.size() ; ee++) {
                xi[ee] = XSValues[ee];
                yi[ee] = YSValues[ee];

                //std::cout << XSValues.size() << " " << xi[ee] << " " << yi[ee] <<std::endl ;
            }

            TGraph* gr1 = CreateGraph (XSValues.size(), xi, yi);

            std::string graph_label = GeometrySymbol+", "+ RadioTracerName;
            gr1->SetName(graph_label.c_str());
            gr1->SetTitle(gr1->GetName());

            mg->Add(setGraphData(gr1, jj, jj));
            leg->AddEntry(gr1,gr1->GetName(),"LP");  // to add the explanation of this colored line
            jj++;
        }
    }

    TF1 *LineFunc = new TF1("Line", "pol1"); //1st degree polynom, pol2, pol3, pol4, ....
    LineFunc->FixParameter(0, 0.0); //0 for fisrt factor P0 in polynom, fix P0 initial point to be 0    mg->Fit(LineFunc, "S");

    //LineFunc->SetParLimits(0, 0.0, 1e-10);

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogx(1);
        gPad->SetLogy(1);
    }
    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = QuantitiesToScore +" scatter plot over line, in all combinations, Geometries, and Radionuclides;"+CompareReferenceName + ";DoseCalcs";
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = ";"+CompareReferenceName + ";DoseCalcs";
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("AP");
    //LineFunc->Draw("same");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;



    std::cout << "\n\n-------------------------------- In One Graph Result-reference Scatter plot, for all combinations and radiotracers, with one geometry  "<< QuantitiesToScore << " ----------------------" << std::endl ;
    for(int aa = 0; aa < GeometryList.size(); aa++){
        FileName = GraphsDirectoryPath+"Radionuclide_ScatterPlot_"+QuantitiesToScore+"_Combinations_Radionuclides_"+GeometrySymbol+""+GraphsExt;
        ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
        mg = new TMultiGraph();
        leg = new TLegend();
        jj = 0;

        GeometrySymbol = GeometryList[aa];
        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            std::vector<double> XSValues;
            std::vector<double> YSValues;

            for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
                for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {

                    Source_ORG = SourceNamesToScore[ee];
                    Target_ORG = TargetNamesToScore[gg];

                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    if(     a1==NULL || a1==0. || a1 == MinValForLog || isnanl(a1) || isinfl(a1) ||
                            a2==NULL || a2==0. || a2 == MinValForLog || isnanl(a2) || isinfl(a2)){
                        continue;
                    }

                    //if(a1 > 2 || a1 < 0.01){a1 = a2;}

                    XSValues.push_back(a2);
                    YSValues.push_back(a1);
                }
            }

            if(XSValues.size() == 0){
                continue;
            }

            double xi[XSValues.size()]; double yi[XSValues.size()];
            for (int ee = 0 ; ee < XSValues.size() ; ee++) {
                xi[ee] = XSValues[ee];
                yi[ee] = YSValues[ee];

                //std::cout << XSValues.size() << " " << xi[ee] << " " << yi[ee] <<std::endl ;
            }

            TGraph* gr1 = CreateGraph (XSValues.size(), xi, yi);

            std::string graph_label = GeometrySymbol+", "+ RadioTracerName;
            gr1->SetName(graph_label.c_str());
            gr1->SetTitle(gr1->GetName());

            mg->Add(setGraphData(gr1, jj, jj));
            leg->AddEntry(gr1,gr1->GetName(),"LP");  // to add the explanation of this colored line
            jj++;
        }

        TF1 *LineFunc = new TF1("Line", "pol1"); //1st degree polynom, pol2, pol3, pol4, ....
        LineFunc->FixParameter(0, 0.0); //0 for fisrt factor P0 in polynom, fix P0 initial point to be 0    mg->Fit(LineFunc, "S");
        mg->Fit(LineFunc, "S");

        if(QuantityUseLog[QuantitiesToScore]){
            gPad->SetLogx(1);
            gPad->SetLogy(1);
        }
        if(UseGridXY=="yes"){
            gPad->SetGridx();
            gPad->SetGridy();
        }

        multiGraphTitle = QuantitiesToScore +" scatter plot over line, in all combinations, " + GeometrySymbol + ", Radionuclides;"+CompareReferenceName + ";DoseCalcs";
        if(PrintTitle == "yes"){
            mg->SetTitle(multiGraphTitle.c_str());
        }
        else{
            multiGraphTitle = ";"+CompareReferenceName + ";DoseCalcs";
            mg->SetTitle(multiGraphTitle.c_str());
        }

        mg->GetXaxis()->CenterTitle(true);
        mg->GetYaxis()->CenterTitle(true);
        mg->GetXaxis()->SetTitleOffset(1.3);

        mg->GetHistogram()->SetMinimum();
        mg->GetHistogram()->SetMaximum();

        mg->Draw("AP");

        gPad->SetRightMargin(0.13);
        gPad->SetLeftMargin(0.2);

        leg->SetX1(X1LegPos);
        leg->SetX2(X2LegPos);
        leg->SetY1(Y1LegPos);
        leg->SetY2(Y2LegPos);
        leg->Draw();

        ResCanvas->Print(FileName.c_str());
        //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
        delete ResCanvas;

    }



    std::cout << "\n\n-------------------------------- In One Graph Result-reference Scatter plot, for all combinations, geometries and one radionuclide"<< QuantitiesToScore << " ----------------------" << std::endl ;
    for(int bb = 0; bb < RadiotracerList.size(); bb++){

        RadioTracerName = RadiotracerList[bb];
        FileName = GraphsDirectoryPath+"Radionuclide_ScatterPlot_"+QuantitiesToScore+"_Combinations_Geometries_"+RadioTracerName+""+GraphsExt;
        ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
        mg = new TMultiGraph();
        leg = new TLegend();
        jj = 0;

        for(int aa = 0; aa < GeometryList.size(); aa++){

            GeometrySymbol = GeometryList[aa];

            std::vector<double> XSValues;
            std::vector<double> YSValues;

            for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
                for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {

                    Source_ORG = SourceNamesToScore[ee];
                    Target_ORG = TargetNamesToScore[gg];

                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    if(     a1==NULL || a1==0. || a1 == MinValForLog || isnanl(a1) || isinfl(a1) ||
                            a2==NULL || a2==0. || a2 == MinValForLog || isnanl(a2) || isinfl(a2)){
                        continue;
                    }
                    XSValues.push_back(a2);
                    YSValues.push_back(a1);
                }
            }

            if(XSValues.size() == 0){
                continue;
            }

            double xi[XSValues.size()]; double yi[XSValues.size()];
            for (int ee = 0 ; ee < XSValues.size() ; ee++) {
                xi[ee] = XSValues[ee];
                yi[ee] = YSValues[ee];

                //std::cout << XSValues.size() << " " << xi[ee] << " " << yi[ee] <<std::endl ;
            }

            TGraph* gr1 = CreateGraph (XSValues.size(), xi, yi);

            std::string graph_label = GeometrySymbol+", "+ RadioTracerName;
            gr1->SetName(graph_label.c_str());
            gr1->SetTitle(gr1->GetName());

            mg->Add(setGraphData(gr1, jj, jj));
            leg->AddEntry(gr1,gr1->GetName(),"LP");  // to add the explanation of this colored line
            jj++;
        }

        TF1 *LineFunc = new TF1("Line", "pol1"); //1st degree polynom, pol2, pol3, pol4, ....
        LineFunc->FixParameter(0, 0.0); //0 for fisrt factor P0 in polynom, fix P0 initial point to be 0    mg->Fit(LineFunc, "S");
        mg->Fit(LineFunc, "S");

        if(QuantityUseLog[QuantitiesToScore]){
            gPad->SetLogx(1);
            gPad->SetLogy(1);
        }
        if(UseGridXY=="yes"){
            gPad->SetGridx();
            gPad->SetGridy();
        }

        multiGraphTitle = QuantitiesToScore +" scatter plot over line, in all combinations, " + RadioTracerName + ", Geometries;"+CompareReferenceName + ";DoseCalcs";
        if(PrintTitle == "yes"){
            mg->SetTitle(multiGraphTitle.c_str());
        }
        else{
            multiGraphTitle = ";"+CompareReferenceName + ";DoseCalcs";
            mg->SetTitle(multiGraphTitle.c_str());
        }

        mg->GetXaxis()->CenterTitle(true);
        mg->GetYaxis()->CenterTitle(true);
        mg->GetXaxis()->SetTitleOffset(1.3);

        mg->GetHistogram()->SetMinimum();
        mg->GetHistogram()->SetMaximum();

        mg->Draw("AP");

        gPad->SetRightMargin(0.13);
        gPad->SetLeftMargin(0.2);

        leg->SetX1(X1LegPos);
        leg->SetX2(X2LegPos);
        leg->SetY1(Y1LegPos);
        leg->SetY2(Y2LegPos);
        leg->Draw();

        ResCanvas->Print(FileName.c_str());
        //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
        delete ResCanvas;

    }



    std::cout << "\n\n-------------------------------- In One Graph Result-reference Scatter plot, for all combinations with one radiotracer and one geometry  "<< QuantitiesToScore << " ----------------------" << std::endl ;
    for(int aa = 0; aa < GeometryList.size(); aa++){

        GeometrySymbol = GeometryList[aa];
        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            FileName = GraphsDirectoryPath+"Radionuclide_ScatterPlot_"+QuantitiesToScore+"_Combinations_"+RadioTracerName+"_"+GeometrySymbol+""+GraphsExt;
            ResCanvas = new TCanvas("Scatter plots", "Scatter plots");

            std::vector<double> XSValues;
            std::vector<double> YSValues;

            for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
                for (int gg = 0 ; gg < TargetNamesToScore.size() ; gg++) {
                    Source_ORG = SourceNamesToScore[ee];
                    Target_ORG = TargetNamesToScore[gg];

                    double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                    if(     a1==NULL || a1==0. || a1 == MinValForLog || isnanl(a1) || isinfl(a1) ||
                            a2==NULL || a2==0. || a2 == MinValForLog || isnanl(a2) || isinfl(a2)){
                        continue;
                    }
                    XSValues.push_back(a2);
                    YSValues.push_back(a1);
                }
            }

            if(XSValues.size() == 0){
                continue;
            }

            double xi[XSValues.size()]; double yi[XSValues.size()];
            for (int ee = 0 ; ee < XSValues.size() ; ee++) {
                xi[ee] = XSValues[ee];
                yi[ee] = YSValues[ee];

                //std::cout << XSValues.size() << " " << xi[ee] << " " << yi[ee] <<std::endl ;
            }


            TGraph* gr1 = CreateGraph (XSValues.size(), xi, yi);

            gr1->SetMarkerSize(1.2);
            gr1->SetMarkerStyle(8);
            gr1->SetMarkerColor(4);

            TF1 *LineFunc = new TF1("Line", "pol1"); //1st degree polynom, pol2, pol3, pol4, ....
            LineFunc->FixParameter(0, 0.0); //0 for fisrt factor P0 in polynom, fix P0 initial point to be 0    mg->Fit(LineFunc, "S");
            gr1->Fit(LineFunc, "S");

            if(QuantityUseLog[QuantitiesToScore]){
                gPad->SetLogx(1);
                gPad->SetLogy(1);
            }
            if(UseGridXY=="yes"){
                gPad->SetGridx();
                gPad->SetGridy();
            }

            //std::cout << "----------------- 1 -----------------" <<std::endl ;

            multiGraphTitle = QuantitiesToScore + " scatter plot over line, in all combinations, " + GeometrySymbol + ", " + RadioTracerName +";"+CompareReferenceName + ";DoseCalcs";
            if(PrintTitle == "yes"){
                gr1->SetTitle(multiGraphTitle.c_str());
            }
            else{
                multiGraphTitle = ";"+CompareReferenceName + ";DoseCalcs";
                gr1->SetTitle(multiGraphTitle.c_str());
            }

            gr1->GetXaxis()->CenterTitle(true);
            gr1->GetYaxis()->CenterTitle(true);
            gr1->GetXaxis()->SetTitleOffset(1.3);

            gr1->GetHistogram()->SetMinimum();
            gr1->GetHistogram()->SetMaximum();

            //std::cout << "----------------- 1 -----------------" <<std::endl ;

            gr1->Draw("AP");

            gPad->SetRightMargin(0.13);
            gPad->SetLeftMargin(0.2);

            //std::cout << "----------------- 1 -----------------" <<std::endl ;

            FileName = GraphsDirectoryPath+"Radionuclide_ScatterPlot_"+QuantitiesToScore+"_Combinations_"+RadioTracerName+"_"+GeometrySymbol+""+GraphsExt;
            //FileName = GraphsDirectoryPath+"Radionuclide_"+DifferenceMethod+"_"+QuantitiesToScore+"_SourceRegions_Radiotracers_Geometry"+GraphsExt;

            ResCanvas->Print(FileName.c_str());
            //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
            delete ResCanvas;
        }
    }



    std::cout << "\n\n-------------------------------- In One Graph Result-reference Scatter plot, for source regions and radiotracers and geometries  "<< QuantitiesToScore << " ----------------------" << std::endl ;
    FileName = GraphsDirectoryPath+"Radionuclide_ScatterPlot_"+QuantitiesToScore+"_SourceRegions_Radionuclides_Geometries"+GraphsExt;
    ResCanvas = new TCanvas(FileName.c_str(), FileName.c_str());
    mg = new TMultiGraph();
    leg = new TLegend();
    jj = 0;
    for(int aa = 0; aa < GeometryList.size(); aa++){

        GeometrySymbol = GeometryList[aa];
        for(int bb = 0; bb < RadiotracerList.size(); bb++){

            RadioTracerName = RadiotracerList[bb];

            std::vector<double> XSValues;
            std::vector<double> YSValues;

            for (int ee = 0 ; ee < SourceNamesToScore.size() ; ee++) {
                Source_ORG = SourceNamesToScore[ee];
                Target_ORG = Source_ORG;

                double a1 = ResultQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                double a2 = ReferenceQuantityGeometryRadioTracerSourceTargetValues[QuantitiesToScore][GeometrySymbol][RadioTracerName][Source_ORG][Target_ORG];
                if(     a1==NULL || a1==0. || a1 == MinValForLog || isnanl(a1) || isinfl(a1) ||
                        a2==NULL || a2==0. || a2 == MinValForLog || isnanl(a2) || isinfl(a2)){
                    continue;
                }
                XSValues.push_back(a2);
                YSValues.push_back(a1);
            }

            if(XSValues.size() == 0){
                continue;
            }

            double xi[XSValues.size()]; double yi[XSValues.size()];
            for (int ee = 0 ; ee < XSValues.size() ; ee++) {
                xi[ee] = XSValues[ee];
                yi[ee] = YSValues[ee];

                //std::cout << XSValues.size() << " " << xi[ee] << " " << yi[ee] <<std::endl ;
            }

            TGraph* gr1 = CreateGraph (XSValues.size(), xi, yi);

            std::string graph_label = GeometrySymbol+", "+ RadioTracerName;
            gr1->SetName(graph_label.c_str());
            gr1->SetTitle(gr1->GetName());

            mg->Add(setGraphData(gr1, jj, jj));
            leg->AddEntry(gr1,gr1->GetName(),"LP");  // to add the explanation of this colored line
            jj++;
        }
    }

    LineFunc = new TF1("Line", "pol1"); //1st degree polynom, pol2, pol3, pol4, ....
    LineFunc->FixParameter(0, 0.0); //0 for fisrt factor P0 in polynom, fix P0 initial point to be 0    mg->Fit(LineFunc, "S");
    mg->Fit(LineFunc, "S");

    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogx(1);
        gPad->SetLogy(1);
    }
    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    multiGraphTitle = QuantitiesToScore +" scatter plot over line, in source regions, Geometries, Radionuclides;"+CompareReferenceName + ";DoseCalcs";
    if(PrintTitle == "yes"){
        mg->SetTitle(multiGraphTitle.c_str());
    }
    else{
        multiGraphTitle = ";"+CompareReferenceName + ";DoseCalcs";
        mg->SetTitle(multiGraphTitle.c_str());
    }

    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetTitleOffset(1.3);

    mg->GetHistogram()->SetMinimum();
    mg->GetHistogram()->SetMaximum();

    mg->Draw("AP");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.2);

    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);
    leg->Draw();

    ResCanvas->Print(FileName.c_str());
    //ResCanvas->Print(CombinedOutFileName.c_str(), ("Title:"+multiGraphTitle).c_str());
    delete ResCanvas;

}

// for step and event level for Voxelized
// called from GenerateVoxel2DGraphs()
void G4DoseCalcsAnalysis::ReadVoxelsResult(){

    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;

    //DoseProfilQuantity = "AE";

    std::ostringstream filename1;
    filename1 << ResultDirectoryPath << "/VoxelsResults" ;
    std::ifstream file(filename1.str().c_str() , std::ios::binary);
    if(file.is_open()){

        std::cout << "Reading file " << filename1.str().c_str() << std::endl ;

        double Z, Y, X, mass;
        std::string word, line;

        getline(file, line); // Simulation Data line
        //std::cout << " line " << line << std::endl;
        getline(file, line); // Simulation Data line
        //std::cout << " line " << line << std::endl;
        getline(file, line); // Table column comment line
        //std::cout << " line " << line << std::endl;
        getline(file, line); // Simulation Data line
        //std::cout << " line " << line << std::endl;

        while(file.peek()!=EOF){

            file >> X ; //CN
            file >> X ;
            file >> Y ;
            file >> Z ;
            file >> mass ;

            //std::cout << " Z " << Z << " Y " << Y << " X " << X << std::endl;

            file >> QuantVals["AE"][X][Y][Z];
            //std::cout << "AE " << QuantVals["AE"][X][Y][Z] << std::endl;
            file >> QuantVals["AD"][X][Y][Z];
            //std::cout << "AD " << QuantVals["AD"][X][Y][Z] << std::endl;
            file >> QuantVals["AF"][X][Y][Z];
            //std::cout << "AF " << QuantVals["AF"][X][Y][Z] << std::endl;
            file >> QuantVals["SAF"][X][Y][Z];
            //std::cout << "SAF " << QuantVals["SAF"][X][Y][Z] << std::endl;
            file >> QuantVals["S"][X][Y][Z];
            //std::cout << "S " << QuantVals["S"][X][Y][Z] << std::endl;
            file >> QuantVals["H"][X][Y][Z];
            //std::cout << "H " << QuantVals["H"][X][Y][Z] << std::endl;
            file >> QuantVals["E"][X][Y][Z];
            //std::cout << "E " << QuantVals["E"][X][Y][Z] << std::endl;

            //std::cout << X << " "  << Y << " " << Z<< " " << QuantVals[DoseProfilQuantity][X][Y][Z] << std::endl;

            file >> QuantValsErr[X][Y][Z];
            //std::cout << "RelErr " << " " << QuantValsErr[X][Y][Z] << std::endl;

            if(MaximalDoseVal < QuantVals[DoseProfilQuantity][X][Y][Z]){
                MaximalDoseVal = QuantVals[DoseProfilQuantity][X][Y][Z];
                ZMAX = Z ; YMAX = Y ; XMAX = X;
            }

        }

        file.close();
    }


    double Voxel0PosX = -((VoxXNumber*VoxXHalfSize) - VoxXHalfSize);
    double Voxel0PosY = -((VoxYNumber*VoxYHalfSize) - VoxYHalfSize);
    double Voxel0PosZ = -((VoxZNumber*VoxZHalfSize) - VoxZHalfSize);

    for (unsigned int i = 0; i < VoxXNumber; i++) {
        XPos[i] = Voxel0PosX + i * 2 * VoxXHalfSize;
        //std::cout << " i "  << i <<" XPos[i] " << XPos[i] << std::endl ;
    }
    for (unsigned int j = 0; j < VoxYNumber; j++) {
        YPos[j] = Voxel0PosY + j * 2 * VoxYHalfSize;
    }
    for (unsigned int k = 0; k < VoxZNumber; k++) {
        ZPos[k] = Voxel0PosZ + k * 2 * VoxZHalfSize;
    }


}
void G4DoseCalcsAnalysis::CreatePercentageDepthDoseGraph(){

    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;

    TCanvas* PDDCanvas1 = new TCanvas("PercentageDepthDose", "PercentageDepthDose");;
    std::string histTitle1 = "Percentage Depth("+BeamAxis+") Dose Curve; Depth(mm); Relative Absorbed Dose(%)";

    double minn=0, maxx=10; int numbint = 10;
    if(BeamAxis == "X"){ numbint = VoxXNumber; minn = -VoxXHalfSize*VoxXNumber; maxx = VoxXHalfSize*VoxXNumber;}
    if(BeamAxis == "Y"){ numbint = VoxYNumber; minn = -VoxYHalfSize*VoxYNumber; maxx = VoxYHalfSize*VoxYNumber;}
    if(BeamAxis == "Z"){ numbint = VoxZNumber; minn = -VoxZHalfSize*VoxZNumber; maxx = VoxZHalfSize*VoxZNumber;}
    TH1F* histPDD = new TH1F ("PDD", "", numbint, minn, maxx);

    //std::cout << "VoxZNumber " << VoxZNumber<<  " VoxYNumber " << VoxYNumber << " VoxXNumber " << VoxXNumber << " BeamAxis "<< BeamAxis<<  std::endl;

    double X,Y,Z, val;
    for ( auto it = QuantVals[DoseProfilQuantity].begin(); it != QuantVals[DoseProfilQuantity].end(); ++it  ){
        X = it->first;
        for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2  ){
            Y = it2->first;
            for ( auto it3 = it2->second.begin(); it3 != it2->second.end(); ++it3  ){
                Z = it3->first;
                val = it3->second;

                //std::cout << Z << " " << val << std::endl;

                if(BeamAxis == "X"){histPDD->Fill(X, val/MaximalDoseVal);}
                else if(BeamAxis == "Y"){histPDD->Fill(Y, val/MaximalDoseVal);}
                else if(BeamAxis == "Z"){histPDD->Fill(Z, val/MaximalDoseVal);}
            }
        }
    }

    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    histPDD->SetTitle(const_cast<char*>(histTitle1.c_str()));
    histPDD->GetXaxis()->CenterTitle(true);
    histPDD->GetYaxis()->CenterTitle(true);
    histPDD->GetXaxis()->SetTitleOffset(1.3);
    histPDD->SetFillStyle(3001);
    histPDD->SetFillColor(kBlue);
    histPDD->Draw("hist same");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);

    std::string FileName1 = GraphsDirectoryPath+"PercentageDepthDose"+GraphsExt;
    PDDCanvas1->Print(FileName1.c_str());
    delete PDDCanvas1;








    // ////////////////////////////////////////////// Dose Profile ///////////////////////////////////////////




    PDDCanvas1 = new TCanvas("DoseProfile", "DoseProfile");;

    minn=0, maxx=10; numbint = 10;
    if(BeamAxis == "X"){ numbint = VoxYNumber; minn = -VoxYHalfSize*VoxYNumber; maxx = VoxYHalfSize*VoxYNumber;}
    if(BeamAxis == "Y"){ numbint = VoxZNumber; minn = -VoxZHalfSize*VoxZNumber; maxx = VoxZHalfSize*VoxZNumber;}
    if(BeamAxis == "Z"){ numbint = VoxXNumber; minn = -VoxXHalfSize*VoxXNumber; maxx = VoxXHalfSize*VoxXNumber;}

    TH1F* histDP = new TH1F ("DoseProfile", "", numbint, minn, maxx);

    //std::cout << "VoxZNumber " << VoxZNumber<<  " VoxYNumber " << VoxYNumber << " VoxXNumber " << VoxXNumber << " BeamAxis "<< BeamAxis<<  std::endl;


    std::map<double,double> values;

    for ( auto it = QuantVals[DoseProfilQuantity].begin(); it != QuantVals[DoseProfilQuantity].end(); ++it  ){
        X = it->first;
        for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2  ){
            Y = it2->first;
            for ( auto it3 = it2->second.begin(); it3 != it2->second.end(); ++it3  ){
                Z = it3->first;
                values[X] += it3->second;
                values[Y] += it3->second;
                values[Z] += it3->second;
                //std::cout << Z << " " << val << std::endl;

                //if(BeamAxis == "X"){histDP->Fill(X, val/MaximalDoseVal);}
                //else if(BeamAxis == "Y"){histDP->Fill(Y, val/MaximalDoseVal);}
                //else if(BeamAxis == "Z"){histDP->Fill(Z, val/MaximalDoseVal);}
            }
        }
    }

    double voxelvol = 8*VoxXHalfSize*VoxYHalfSize*VoxZHalfSize;
    //double X,Y,Z, val;

    if(BeamAxis == "Z"){
        histTitle1 = "Dose Profile for beam axis ("+BeamAxis+"); X(mm); (Absorbed Dose)/mm";
        for (unsigned int i = 0; i < VoxXNumber; i++) {
            histDP->Fill(XPos[i], values[XPos[i]]/voxelvol);
            //std::cout << " i "  << i <<" XPos[i] " << XPos[i] << std::endl ;
        }
    }
    else if(BeamAxis == "X"){
        histTitle1 = "Dose Profile for beam axis ("+BeamAxis+"); Y(mm); (Absorbed Dose)/mm";
        for (unsigned int i = 0; i < VoxYNumber; i++) {
            histDP->Fill(YPos[i], values[YPos[i]]/voxelvol);
            //std::cout << " i "  << i <<" XPos[i] " << XPos[i] << std::endl ;
        }
    }
    else if(BeamAxis == "Y"){
        histTitle1 = "Dose Profile for beam axis ("+BeamAxis+"); Z(mm); (Absorbed Dose)/mm";
        for (unsigned int i = 0; i < VoxZNumber; i++) {
            histDP->Fill(ZPos[i], values[ZPos[i]]/voxelvol);
            //std::cout << " i "  << i <<" XPos[i] " << XPos[i] << std::endl ;
        }
    }


    if(UseGridXY=="yes"){
        gPad->SetGridx();
        gPad->SetGridy();
    }

    histDP->SetTitle(const_cast<char*>(histTitle1.c_str()));
    histDP->GetXaxis()->CenterTitle(true);
    histDP->GetYaxis()->CenterTitle(true);
    histDP->GetXaxis()->SetTitleOffset(1.3);
    histDP->SetFillStyle(3001);
    histDP->SetFillColor(kBlue);
    histDP->Draw("hist same");

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);

    FileName1 = GraphsDirectoryPath+"DoseProfile"+GraphsExt;
    PDDCanvas1->Print(FileName1.c_str());
    delete PDDCanvas1;



}
void G4DoseCalcsAnalysis::CreateDoseProfile(){

    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;

    //gStyle->SetOptStat(kFALSE); // for histograms

    std::stringstream GraphName, fileName ;

    TCanvas *DoseProMap = new TCanvas("2D Dose Heat Map","2D Dose Heat Map");
    TH2F *ProfDoseHist;

    double minnX = -VoxXHalfSize*VoxXNumber; double maxxX = VoxXHalfSize*VoxXNumber;
    double minnY = -VoxYHalfSize*VoxYNumber; double maxxY = VoxYHalfSize*VoxYNumber;
    double minnZ = -VoxZHalfSize*VoxZNumber; double maxxZ = VoxZHalfSize*VoxZNumber;

    double X,Y,Z, val;
    if(SliceFor2DGraph == "XY" || SliceFor2DGraph == "YX"){

        GraphName << "Dose Heat Map in "<< SliceFor2DGraph.c_str() << "-Slice(" << SliceID <<") Plane ; X(mm) ; Y(mm)"  ;
        ProfDoseHist = new TH2F("2D Dose", GraphName.str().c_str(),  VoxXNumber, minnX , maxxX ,  VoxYNumber, minnY, maxxY);

        for ( auto it = QuantVals[DoseProfilQuantity].begin(); it != QuantVals[DoseProfilQuantity].end(); ++it  ){
            X = it->first;
            for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2  ){
                Y = it2->first;
                for ( auto it3 = it2->second.begin(); it3 != it2->second.end(); ++it3  ){
                    Z = it3->first;
                    val = it3->second;
                    ProfDoseHist->Fill(X,Y,val);
                }
            }
        }
    }
    else if(SliceFor2DGraph == "ZX" || SliceFor2DGraph == "XZ"){
        GraphName << "Dose Heat Map in "<< SliceFor2DGraph.c_str() << "-Slice(" << SliceID <<") Plane ; Z(mm) ; X(mm)"  ;
        ProfDoseHist = new TH2F("2D Dose", GraphName.str().c_str(),  VoxZNumber, minnZ , maxxZ ,  VoxXNumber, minnX, maxxX);
        for ( auto it = QuantVals[DoseProfilQuantity].begin(); it != QuantVals[DoseProfilQuantity].end(); ++it  ){
            X = it->first;
            for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2  ){
                Y = it2->first;
                for ( auto it3 = it2->second.begin(); it3 != it2->second.end(); ++it3  ){
                    Z = it3->first;
                    val = it3->second;
                    ProfDoseHist->Fill(Z,X,val);
                }
            }
        }
    }
    else if(SliceFor2DGraph == "ZY" || SliceFor2DGraph == "YZ"){
        GraphName << "Dose Heat Map in "<< SliceFor2DGraph.c_str() << "-Slice(" << SliceID <<") Plane ; Z(mm) ; Y(mm)"  ;
        ProfDoseHist = new TH2F("2D Dose", GraphName.str().c_str(),  VoxZNumber, minnZ , maxxZ ,  VoxYNumber, minnY, maxxY);
        for ( auto it = QuantVals[DoseProfilQuantity].begin(); it != QuantVals[DoseProfilQuantity].end(); ++it  ){
            X = it->first;
            for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2  ){
                Y = it2->first;
                for ( auto it3 = it2->second.begin(); it3 != it2->second.end(); ++it3  ){
                    Z = it3->first;
                    val = it3->second;
                    ProfDoseHist->Fill(Z,Y,val);
                }
            }
        }
    }


    /*

    std::map<unsigned int,double> axe1, axe2;
    if(SliceFor2DGraph == "XY" || SliceFor2DGraph == "YX"){ axe1 = XPos; axe2 = YPos;}
    if(SliceFor2DGraph == "ZX" || SliceFor2DGraph == "XZ"){ axe1 = ZPos; axe2 = XPos;}
    if(SliceFor2DGraph == "ZY" || SliceFor2DGraph == "YZ"){ axe1 = ZPos; axe2 = YPos;}

    if(SliceFor2DGraph == "XY" || SliceFor2DGraph == "YX"){
        GraphName << "Dose Profile in "<< SliceFor2DGraph.c_str() <<" Plane ; X(mm) ; Y(mm)"  ;
        ProfDoseHist = new TH2F("2D Dose", GraphName.str().c_str(),  VoxXNumber, minnX , maxxX ,  VoxYNumber, minnY, maxxY);
        for (unsigned int k = 0; k < VoxZNumber; k++) {
            for (unsigned int j = 0; j < VoxYNumber; j++) {
                for (unsigned int i = 0; i < VoxXNumber; i++) {
                    ProfDoseHist->Fill(axe1[i],axe2[j],QuantVals[DoseProfilQuantity][k][j][i]);
                }
            }
        }
    }
    else if(SliceFor2DGraph == "ZX" || SliceFor2DGraph == "XZ"){
        GraphName << "Dose Profile in "<< SliceFor2DGraph.c_str() <<" Plane ; Z(mm) ; X(mm)"  ;
        ProfDoseHist = new TH2F("2D Dose", GraphName.str().c_str(),  VoxZNumber, minnZ , maxxZ ,  VoxXNumber, minnX, maxxX);
        for (unsigned int k = 0; k < VoxZNumber; k++) {
            for (unsigned int j = 0; j < VoxYNumber; j++) {
                for (unsigned int i = 0; i < VoxXNumber; i++) {
                    ProfDoseHist->Fill(axe1[k],axe2[i],QuantVals[DoseProfilQuantity][k][j][i]);
                }
            }
        }
    }
    else if(SliceFor2DGraph == "ZY" || SliceFor2DGraph == "YZ"){
        GraphName << "Dose Profile in "<< SliceFor2DGraph.c_str() <<" Plane ; Z(mm) ; Y(mm)"  ;
        ProfDoseHist = new TH2F("2D Dose", GraphName.str().c_str(),  VoxZNumber, minnZ , maxxZ ,  VoxYNumber, minnY, maxxY);
        for (unsigned int k = 0; k < VoxZNumber; k++) {
            for (unsigned int j = 0; j < VoxYNumber; j++) {
                for (unsigned int i = 0; i < VoxXNumber; i++) {
                    ProfDoseHist->Fill(axe1[k],axe2[j],QuantVals[DoseProfilQuantity][k][j][i]);
                }
            }
        }
    }

*/
    ProfDoseHist->SetTitle(GraphName.str().c_str());

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);
    gPad->Update();
    ProfDoseHist->GetXaxis()->CenterTitle(true);
    ProfDoseHist->GetYaxis()->CenterTitle(true);
    ProfDoseHist->GetXaxis()->SetTitleOffset(1.3);

    ProfDoseHist->Draw("colz");

    fileName << GraphsDirectoryPath.c_str() << "DoseProfileHist" << "_" << SliceFor2DGraph.c_str() << "_Slice" << SliceID << GraphsExt.c_str();
    DoseProMap->Print(fileName.str().c_str());//,GraphsExt.c_str());  // create PS file in the dir were this code executed

    //delete ProfDoseHist;
    delete DoseProMap;


    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    TCanvas * Canv3DdoseProMap = new TCanvas("3D Dose Profile","3D Dose Profile");
    TH3F *Prof3DDoseHist;
    std::stringstream GraphName3D, fileName3D ;
    GraphName3D << "3D Dose Profile ; Z(mm) ; Y(mm) ; X(mm)"  ;
    Prof3DDoseHist = new TH3F("3D Dose", GraphName3D.str().c_str(),  VoxZNumber, minnZ , maxxZ ,  VoxYNumber, minnY, maxxY ,  VoxXNumber, minnX, maxxX);

    for ( auto it = QuantVals[DoseProfilQuantity].begin(); it != QuantVals[DoseProfilQuantity].end(); ++it  ){
        X = it->first;
        for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2  ){
            Y = it2->first;
            for ( auto it3 = it2->second.begin(); it3 != it2->second.end(); ++it3  ){
                Z = it3->first;
                val = it3->second;
                Prof3DDoseHist->Fill(Z,Y,X,val);
            }
        }
    }

    Prof3DDoseHist->SetTitle(GraphName3D.str().c_str());

    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);
    gPad->Update();
    Prof3DDoseHist->GetXaxis()->CenterTitle(true);
    Prof3DDoseHist->GetYaxis()->CenterTitle(true);
    Prof3DDoseHist->GetXaxis()->SetTitleOffset(1.3);

    Prof3DDoseHist->Draw();

    fileName3D << GraphsDirectoryPath.c_str() << "3DdoseProfileHist" << GraphsExt.c_str();
    Canv3DdoseProMap->Print(fileName3D.str().c_str());//,GraphsExt.c_str());  // create PS file in the dir were this code executed

    //delete ProfDoseHist;
    delete Canv3DdoseProMap;
}
void G4DoseCalcsAnalysis::CreateHeatMap(){

    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;

    //gStyle->SetOptStat(kFALSE); // for histograms

    std::stringstream GraphName, fileName ;

    TCanvas *Can2DMap = new TCanvas(QuantitiesToScore.c_str(),QuantitiesToScore.c_str());
    TGraph2D *graph2D = new TGraph2D();

    unsigned int FirstAxisSize;
    unsigned int SecondAxisSize;
    int inc;

    SliceID = ZMAX;

    if(SliceFor2DGraph == "XY"){

        inc = 0;
        FirstAxisSize = VoxXNumber;
        SecondAxisSize = VoxYNumber;

        GraphName << QuantityUnit[QuantitiesToScore] << " in "<< SliceFor2DGraph.c_str() <<" slice for Z=" << SliceID <<";X(mm)"<< ";Y(mm)"  ;

        for (unsigned int i = 0; i < FirstAxisSize; i++) {
            for (unsigned int j = 0; j < SecondAxisSize; j++) {

                std::cout  << "............." << j << " " << i << std::endl;

                graph2D->SetPoint(inc,XPos[i],YPos[j],QuantVals[QuantitiesToScore][SliceID][j][i]);
                inc++;
            }
        }
    }
    else if(SliceFor2DGraph == "YZ"){

        GraphName << QuantityUnit[QuantitiesToScore] << " in "<< SliceFor2DGraph.c_str() <<" slice for X=" << SliceID <<";Y(mm)"<< ";Z(mm)"  ;

        inc = 0;
        FirstAxisSize = VoxYNumber;
        SecondAxisSize = VoxZNumber;

        for (unsigned int j = 0; j < FirstAxisSize; j++) {
            for (unsigned int k = 0; k < SecondAxisSize; k++) {
                //std::cout  << "............." << j << " " << k << std::endl;

                graph2D->SetPoint(inc,YPos[j],ZPos[k],QuantVals[QuantitiesToScore][k][j][SliceID]);
                inc++;
            }
        }
    }
    else if(SliceFor2DGraph == "XZ"){

        GraphName << QuantityUnit[QuantitiesToScore] << " in "<< SliceFor2DGraph.c_str() <<" slice for Y=" << SliceID <<";X(mm)"<< ";Z(mm)"  ;

        inc = 0;
        FirstAxisSize = VoxXNumber;
        SecondAxisSize = VoxZNumber;

        for (unsigned int i = 0; i < FirstAxisSize; i++) {
            for (unsigned int k = 0; k < SecondAxisSize; k++) {

                graph2D->SetPoint(inc,XPos[i],ZPos[k],QuantVals[QuantitiesToScore][k][SliceID][i]);
                inc++;
            }
        }
    }else {
        std::cout  << "You have to choose the SliceFor2DGraph by setting setSliceFor2DGraph ..." << std::endl;
        return;
    }
    if(PrintTitle == "yes"){
        graph2D->SetTitle(GraphName.str().c_str());
    }
    else{
        GraphName << " ;X(mm)"<< ";Z(mm)"  ;
        graph2D->SetTitle(GraphName.str().c_str());
    }
    gStyle->SetPalette(kRainBow);
    if(QuantityUseLog[QuantitiesToScore]){
        gPad->SetLogz(1);
    }
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);
    //gPad->Update();
    //TPaletteAxis *palette = (TPaletteAxis*)graph->GetListOfFunctions()->FindObject("palette");

    graph2D->Draw("COLZ");

    fileName << GraphsDirectoryPath << "HeatMap_"<< QuantitiesToScore << "_" << GeometrySymbol<< "_" << SliceFor2DGraph << "_Slice" << SliceID << GraphsExt;
    Can2DMap->Print(fileName.str().c_str());//,GraphsExt.c_str());  // create PS file in the dir were this code executed

    delete graph2D ;
    delete Can2DMap;
}

void G4DoseCalcsAnalysis::GenerateCrossSectionGraphs(){

    if(GenerateCrossSectionGraph != "yes"){return;}

    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;

    std::ostringstream fn;
    fn << ResultDirectoryPath << "/CrossSectionData" ;
    std::ifstream file(fn.str().c_str() , std::ios::binary);

    if(file.is_open()){

        std::cout << "\nReading file " << fn.str().c_str() << std::endl ;

        std::string line , indicator, word, MatN, PartN, ProcN;
        double Ene , CrosS;
        int procNum = 0 ;

        std::map<std::string,std::vector<double>> ParticleProcessNames;
        std::vector<std::string> ProcessNames, Words;

        while (getline(file, line)) {

            //std::cout << " the line " << line << std::endl ;

            std::istringstream A(line);

            if(A.str().empty()){
                continue;
            }

            A >> word ;
            if(word == "Table"){
                indicator = word;
                A >> PartN >> MatN;
                //std::cout << "Data for Particle:"<< PartN << " and Material:" << MatN << std::endl;
                ParticleProcessNames.clear();
                continue;
            }
            else if(word == "Energy(MeV)"){

                //std::cout << "The first is Energy(MeV)" << std::endl;

                std::istringstream C(line);
                std::string LineWords ;
                ProcessNames.clear();
                Words.clear();


                while( C >> LineWords ){
                    if(LineWords.empty() || LineWords == ""){continue;}
                    else{
                        //std::cout << LineWords << " " << std::endl;
                        Words.push_back(LineWords);
                    }
                }

                for (unsigned int k = 1; k < Words.size()-18; k++) {
                    //std::cout << LineWords << " " << std::endl;
                    ProcessNames.push_back(Words[k]);
                }

                procNum = ProcessNames.size();
                //std::cout << "Number Of Process for particle "<< PartN << " and material " << MatN << " is " << procNum << std::endl;

                continue;
            }
            else if(word == "Latex"){
                indicator = word;
                continue;
            }

            //std::cout << "\n" << std::endl;

            std::istringstream B(line);
            if(indicator == "Table"){

                //std::cout << "indicator == Table " << std::endl;

                double tt = 0. ;
                B >> Ene;
                for (unsigned int k = 0; k < ProcessNames.size(); k++) {
                    B >> ParticleMaterialProcessEneCrossSection[PartN][MatN][ProcessNames[k]][Ene];
                    tt += ParticleMaterialProcessEneCrossSection[PartN][MatN][ProcessNames[k]][Ene];
                    //std::cout << " - Energy " << Ene <<" - Process :  "<< ProcessNames[k] << " - MacroCrossSection " << ParticleMaterialProcessEneCrossSection[PartN][MatN][ProcessNames[k]][Ene] << " - Total " << tt<< std::endl;
                }
            }
        }
    }


    TCanvas * ResCanvas = new TCanvas("CrossSectionPerVolCanvas", "CrossSectionPerVolCanvas");
    TMultiGraph *mg = new TMultiGraph();
    TLegend *leg = new TLegend();
    leg->SetX1(X1LegPos);
    leg->SetX2(X2LegPos);
    leg->SetY1(Y1LegPos);
    leg->SetY2(Y2LegPos);

    for ( auto Abeg = ParticleMaterialProcessEneCrossSection.begin(); Abeg != ParticleMaterialProcessEneCrossSection.end(); ++Abeg  )
    {
        std::string PARTICLE_NAME = Abeg->first;
        if(PARTICLE_NAME != ParticleName){continue;}

        // iterations on source name
        for ( auto Bbeg = Abeg->second.begin(); Bbeg != Abeg->second.end(); ++Bbeg  )
        {

            std::string MaterialName = Bbeg->first;

            ResCanvas = new TCanvas("", "");
            mg = new TMultiGraph();
            leg = new TLegend();

            leg->SetX1(X1LegPos);
            leg->SetX2(X2LegPos);
            leg->SetY1(Y1LegPos);
            leg->SetY2(Y2LegPos);

            //std::cout << X1LegPos << " " << X2LegPos << " " << Y1LegPos << " " << Y2LegPos << " " << std::endl ;

            gPad->SetRightMargin(0.13);
            gPad->SetLeftMargin(0.13);

            if(UseGridXY=="yes"){
                gPad->SetGridx();
                gPad->SetGridy();
            }

            gPad->SetLogy(1); // for cross section always use log
            if(UseLogE == "yes"){
                gPad->SetLogx(1);
            }

            int jj = 0;

            for ( auto Cbeg = Bbeg->second.begin(); Cbeg != Bbeg->second.end(); ++Cbeg  )
            {
                std::string ProcessName = Cbeg->first;

                int num = Cbeg->second.size();
                double xi[num], yi[num], minSig = 1e-8;

                int jh = 0;
                for ( auto Dbeg = Cbeg->second.begin(); Dbeg != Cbeg->second.end(); ++Dbeg  )
                {
                    xi[jh] = Dbeg->first;
                    if(minSig > Dbeg->second ){
                        yi[jh] = minSig;
                    }else{
                        yi[jh] = Dbeg->second;
                    }

                    jh++;
                }

                TGraph* graph = CreateGraph (num, xi, yi);

                std::ostringstream nm; nm << ProcessName ;
                std::string name = nm.str();

                graph->SetName(const_cast<char*>(name.c_str()));
                graph->SetTitle(graph->GetName());

                mg->Add(setGraphData(graph, jj, jj));
                leg->AddEntry(graph,graph->GetName(),"LP");  // to add the explanation of this colored line

                jj++;
            }

            std::string multiGraphTitle = "Macroscopic Cross Section for " + PARTICLE_NAME + " in material " + MaterialName + ";Energy(MeV); Macroscopic Cross Section(cm-1)";
            std::string FileName = GraphsDirectoryPath + "Macroscopic_Cross_Section_for_" + PARTICLE_NAME + "_in_material_" + MaterialName + GraphsExt;
            if(PrintTitle != "yes"){multiGraphTitle = ";Energy(MeV); Macroscopic Cross Section(cm-1)";}
            CreateMultiGraphParametersAndCanvas(multiGraphTitle, FileName, mg, leg);

        }
    }
}

void G4DoseCalcsAnalysis::GenerateROOTFileForAllResults(){
    std::cout << "\n\n                                                          ========= "<< __FUNCTION__ << " ========= "<< "\n" << std::endl;


}

#endif
