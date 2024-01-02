//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Author: Tarik Elghalbzouri,  Abdelmalek Essaâdi University,
// faculty of sciences Tetouane, morocco. email : telghalbzouri@uae.ac.ma
//
// This application is based on code developed by :
// G. Guerrieri, University of Genova, Italy .
// S. Guatelli. University of Wollongong, Australia.
//

#ifndef G4TResultCalculation_h
#define G4TResultCalculation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"

#ifdef ANALYSIS_USE
#endif

#include <vector>

class G4TResultCalculation
{
public:

    G4TResultCalculation();
    ~G4TResultCalculation();

    void ReadSimulationData();
    void Initialization();
    void setRankID(G4int n ){ RankID = n;}
    void MergeSimulationsData();
    G4String getOneOrMultiSimulations(){return OneOrMultiSimulations;}
    //G4int getRankNum(){return NumberOfRanksThreads;}
    bool DirectoryExists( const char* pzPath );

    void setExeFromMerge(G4bool n ){ExeFromMerge = n;}
    void setUseAllResultsFiles(G4bool n ){UseAllResultsFiles = n;}
    void setV(G4bool n ){V = n;}

private:

    G4bool UseAllResultsFiles;
    G4bool ExeFromMerge;

    G4bool VOX_USE;
    G4bool GenerateVoxelsResuls;

    G4String RankParticleName;
    G4double RankParticleSourceEnergy;
    G4String RankSourceRegionName;

    G4int SZ;
    G4int VarSZ;
    G4int DefaultSZ;

    G4String GeometrySymbol;
    G4String Physics;

    // data created in DoseCalcs executing

    G4String OneOrMultiSimulations;

    std::vector<G4String> OrgansNameVector;

    std::map<G4String,G4double> VolumeNameMassMap;
    std::map<G4String,G4double> VolumeNameVolumeMap;
    std::map<G4String,G4double> VolumeNameDensityMap;

    std::vector<G4String> TargetNamesToScore;
    std::vector<G4String> OldTargetNamesToScore;
    std::vector<G4String> SourceNamesToScore;
    std::map<G4String,std::vector<G4String>> NeWCombinationsForSourceNamesToScore;
    std::map<G4String,std::vector<G4String>> RegionsVolumesNamesMap;
    std::vector<G4String> QuantityNamesToScore;

    G4double CutsDistance;
    G4double CutsEnergy;

    G4bool IsAllTargetsToScore;
    G4bool IsAllSourcesToScore;

    G4String ParticleName;
    G4double ParticleSourceEnergy;
    G4String SourceType;
    G4String SourceRegionName;
    G4String EnergyDistribution;
    G4String MomDirDistribution;

    G4String GraphsData;
    G4String CompareType;
    G4String GraphsExt;
    G4String RefFilePath;
    G4String RefName;
    G4int NumOfEne;
    G4int NumOfRefEne;

    G4String OrgansNamesToScoreString;
    G4String QuantitiesToScore;

    G4String ExecutionMode;
    G4int RankID;

    G4double ExecutionTimeInMin;
    G4double MaxExecutionTimeInMin;
    G4double OneEventExecutionTimeInMs;

    // data created By thread files data

    G4double EnergyEmittedPerThread;

    G4double TotalEmittedEnergy;
    G4double TotalAbsorbedEnergy;
    double long TotalAbsorbedDose = 0;

    unsigned int  eventNumberInEachSimulation;
    G4double EnergyEmittedPerSimulation;

    // data created By internal Merging code

    G4int TotVoxNum;
    G4int VoxXNumber;
    G4int VoxYNumber;
    G4int VoxZNumber;
    G4double VoxXHalfSize;
    G4double VoxYHalfSize;
    G4double VoxZHalfSize;

    G4String RadioTracerName;
    G4String SAFRefForEffectiveDose;
    G4double TotalEffectiveDose;
    G4double InjectedActivity;
    G4double CumulatedActivity;
    G4bool RadiotracerDataFomFile;
    G4bool GenerateResultsForRadioTracer;
    G4bool GenerateResultsForRadioTracerExams;

    std::map<G4String,std::map<G4String,std::map<G4double,G4double>>> RadioTracerEnergyPerCent ;
    //std::map<G4String,std::map<G4String,std::map<std::vector<G4double>,G4double>>> RadioTracerEnergyPerCent ;

    std::map<G4String,std::map<G4String,G4double>> RadioTracerParticleTotalProbability ;
    std::map<G4String,std::map<G4String,G4double>> RadioTracerSourceTi_Fs_ai_AsPerA0 ;
    std::map<G4String,G4double> RadioTracerInjectedActivity ;
    std::map<G4String,std::map<G4String,long double>> TotalDoseFromRadioTracer ;
    std::map<G4String,std::map<G4String,std::map<G4double,std::map<G4String,std::map<G4String,G4double>>>>> ResultTable;
    std::map<G4String,std::map<G4String,std::map<G4double,std::map<G4String,std::map<G4String,G4double>>>>> StandardDeviation;
    std::map<G4String,std::map<G4String,std::map<G4double,std::map<G4String,std::map<G4String,G4double>>>>> TotalAEForRadiotracerRSD;
    std::map<G4String,std::map<G4String,std::map<G4String,G4double>>> RadioTracerSourceTargetSValueForAbsorbedDose ;
    std::map<G4String,std::map<G4String,std::map<G4String,G4double>>> RadioTracerSourceTargetSValueForEquivalentDose ;
    std::map<G4String,std::map<G4String,std::map<G4String,G4double>>> RadioTracerSourceTargetAbsorbedDose ;
    std::map<G4String,std::map<G4String,std::map<G4String,G4double>>> RadioTracerSourceTargetEquivalentDose ;
    std::map<G4String,std::map<G4String,std::map<G4String,G4double>>> RadioTracerSourceTargetEffectiveDose ;
    std::map<G4String,std::map<G4String,std::map<G4String,G4double>>> RadioTracerQuantityOrganValue ;
    std::map<G4String,std::map<G4String,std::map<G4String,G4double>>> RadioTracerQuantityOrganVAR ;
    std::map<G4String,std::map<G4String,std::map<G4String,G4double>>> RadioTracerQuantityOrganSD ;
    std::map<G4String,std::map<G4String,std::map<G4String,G4double>>> RadioTracerQuantityOrganAETotalForRSD ;
    std::map<G4String,std::map<G4String,std::map<G4String,std::map<G4String,G4double>>>> RadioTracerQuantitySourceTargetValue ;
    std::map<G4String,std::map<G4String,std::map<G4String,std::map<G4String,G4double>>>> RadioTracerQuantitySourceTargetVariance ;
    std::map<G4String,std::map<G4String,std::map<G4String,std::map<G4String,G4double>>>> RadioTracerQuantitySourceTargetAETotalForRSD ;

    std::map<G4String,std::ostringstream> RadiTracerParticleEnergyDataString ;
    std::map<G4String,std::ostringstream> RadiTracerDataForTotalDoseString ;
    std::map<G4String,G4double> TotalRadioTracerEmmitedEnergy_Val ;

    G4double MeV_to_J;
    G4double Gy_to_Sv;
    G4double Bq_to_MBq;

    G4double AEUnitFactor;
    G4double AFUnitFactor;
    G4double SAFUnitFactor;
    G4double SUnitFactor;
    G4double ADUnitFactor;
    G4double HUnitFactor;
    G4double EUnitFactor;

    G4String AEUnit;
    G4String AFUnit;
    G4String SAFUnit;
    G4String SUnit;
    G4String ADUnit;
    G4String HUnit;
    G4String EUnit;

    G4String UnitPerParticle;
    G4String UnitPerDecay;

    G4double ResidenceTimeUnitFactor;
    G4String ResidenceTimeUnit;

    G4double AdministeredActivityUnitFactor;
    G4String AdministeredActivityUnit;

    G4bool V; //Verbose

    // for normal calculation
    bool ReadThreadRegionResultFile(G4String);
    void RegionQuantitiesCalculation();

    void G4TCoutReset();

    void GenerateRegionResultFile();
    void GenerateRegionResultForRadioTracer();
    void GenerateRadiotracerQuantitiesByInterpolation(G4String, G4double);

    G4String RemoveWordFromString(G4String,G4String);
    double GenerateRadiationFactor(G4String, double);

    G4String RegionForMaxRel_SDv, RegionForMinRel_SDv;
    G4double MaxRel_SDv, MinRel_SDv;

    G4double PhantomEffectiveDose;

    G4int VoxelCNForMaxRel_SDv, VoxelCNForMinRel_SDv;
    G4double VoxMaxRel_SDv, VoxMinRel_SDv;

    unsigned long long int * VoxNOfValues ;
    std::map<G4String,unsigned long long int> NOfValues;
    unsigned long long int TotalEventNumber;
    unsigned long long int TotalNumberOfSteps;
    std::map<G4String,std::map<G4double,std::map<G4String,std::map<G4String,unsigned long long int>>>> NumberOfSteps;
    std::map<G4String,std::map<G4String,std::map<G4String,unsigned long long int>>> NumberOfStepsInRadiotracerSourceTarget;
    std::map<G4String,std::map<G4String,unsigned long long int>> NumberOfStepsInRadiotracerOrgan;


    std::map<G4String,G4double> DR_Total ;
    std::map<G4String,G4double> ER_Total ;

    std::map<G4String,G4double> ED_Total ;
    std::map<G4String,G4double> ED_Mean ;
    std::map<G4String,G4double> ED2_Total ;
    std::map<G4String,G4double> ED2_Mean ;
    std::map<G4String,G4double> ED_Var ;
    std::map<G4String,G4double> ED_SDev ;
    std::map<G4String,G4double> ED_RelS_D ;

    std::map<G4String,G4double> AF_Total ;
    std::map<G4String,G4double> AF2_Total ;
    std::map<G4String,G4double> AF_Mean ;
    std::map<G4String,G4double> AF2_Mean ;
    std::map<G4String,G4double> AF_Var ;
    std::map<G4String,G4double> AF_SDev ;
    std::map<G4String,G4double> AF_RelS_D ;

    std::map<G4String,G4double> SAF_Total ;
    std::map<G4String,G4double> SAF2_Total ;
    std::map<G4String,G4double> SAF_Mean ;
    std::map<G4String,G4double> SAF2_Mean ;
    std::map<G4String,G4double> SAF_Var ;
    std::map<G4String,G4double> SAF_SDev ;
    std::map<G4String,G4double> SAF_RelS_D ;

    std::map<G4String,G4double> AD_Total ;
    std::map<G4String,G4double> AD2_Total ;
    std::map<G4String,G4double> AD_Mean ;
    std::map<G4String,G4double> AD2_Mean ;
    std::map<G4String,G4double> AD_Var ;
    std::map<G4String,G4double> AD_SDev ;
    std::map<G4String,G4double> AD_RelS_D ;

    std::map<G4String,G4double> S_Total ;
    std::map<G4String,G4double> S2_Total ;
    std::map<G4String,G4double> S_Mean ;
    std::map<G4String,G4double> S2_Mean ;
    std::map<G4String,G4double> S_Var ;
    std::map<G4String,G4double> S_SDev ;
    std::map<G4String,G4double> S_RelS_D ;

    std::map<G4String,G4double> H_Total ;
    std::map<G4String,G4double> H2_Total ;
    std::map<G4String,G4double> H_Mean ;
    std::map<G4String,G4double> H2_Mean ;
    std::map<G4String,G4double> H_Var ;
    std::map<G4String,G4double> H_SDev ;
    std::map<G4String,G4double> H_RelS_D ;

    std::map<G4String,G4double> E_Total ;
    std::map<G4String,G4double> E2_Total ;
    std::map<G4String,G4double> E_Mean ;
    std::map<G4String,G4double> E2_Mean ;
    std::map<G4String,G4double> E_Var ;
    std::map<G4String,G4double> E_SDev ;
    std::map<G4String,G4double> E_RelS_D ;

    std::vector<G4String> DoseCalcsQuantities;
    std::map<G4String,std::map<G4String,G4double>> TotalValueOfQuantity;

    std::map<G4String,G4double> ChosenVariablevariance;
    std::map<G4String,G4double> ChosenVariableStandardDeviation;
    std::map<G4String,G4double> ChosenVariableMean;
    std::map<G4String,G4double> ChosenVariableTotal;

    std::map<G4String, std::map<G4int,G4double>> AE_OrganBatchValueMap;

    std::map<G4String,std::map<G4String,std::map<G4String,std::map<G4double,std::map<G4int,std::vector<G4int>>>>>> GeometrySourceDataRankThreadDataFile;
    std::map<G4String,std::map<G4String,std::map<G4double,std::map<G4String,std::map<G4int,std::vector<G4int>>>>>> GeometryEnergySourceDataRankThreadDataFile;
    std::map<G4String, std::map<G4String,std::map<G4double,std::map<G4int,G4int>>>> SourceDataRankThreadDataFile;
    std::map<G4String, std::map<G4String,std::map<G4double,std::map<G4String,std::vector<G4double>>>>> SourceDataThreadDataFile;
    std::map<G4String, std::map<G4String,std::map<G4double,std::map<G4String,std::vector<G4double>>>>> SourceDataRankDataFile;

    std::map<G4String, G4double> TissueFactorMap; // organ Name, the tissue factor value
    //std::map<G4String, std::map<G4double, G4double>> RadiationFactorMap; // Name of particle, interval of energy, the radiation factor value
    G4int intOfEneForRadFac;

    std::map<G4String,std::map<G4String,std::map<G4String, std::vector<G4double>>>> SourceParticleEnergyValues;

    std::map<G4String, G4double> AFCte ;
    std::map<G4String, G4double> SAFCte ;
    std::map<G4String, G4double> ADCte ;
    std::map<G4String, G4double> SCte ;
    std::map<G4String, G4double> HCte ;
    std::map<G4String, G4double> ECte ;

    std::map<G4String,std::vector<G4String>> ReadedResultFilesPaths;

    // for Voxel Data

    //std::vector<unsigned int> CNX; std::vector<unsigned int> getCNX() const { return CNX;}
    //std::vector<unsigned int> CNY; std::vector<unsigned int> getCNY() const { return CNY;}
    //std::vector<unsigned int> CNZ; std::vector<unsigned int> getCNZ() const { return CNZ;}

    /*
    std::vector<unsigned int> CNID;
    G4double* VoxED_Total ;
    G4double* VoxED2_Total ;
    G4double* CopyNumberZPos;
    G4double* CopyNumberYPos;
    G4double* CopyNumberXPos;
    G4double* CopyNumberMassSize;
    std::map<unsigned int,G4String> CopyNumberRegionNameMap;
    std::map<unsigned int,unsigned int> VoxNOfValues ;
    */

    //G4double* CopyNumberZPos ;
    //G4double* CopyNumberYPos ;
    //G4double* CopyNumberXPos ;
    //G4double* CopyNumberMassSize ;
    //G4String* CopyNumberRegionNameMap ;

    unsigned int* CNID ;
    G4double* VoxED_Total ;
    G4double* VoxED2_Total ;

    G4double* VoxED_Mean ;
    G4double* VoxED2_Mean ;
    G4double* VoxED_Var ;
    G4double* VoxED_SDev ;
    G4double* VoxED_RelS_D ;

    G4double* VoxAF_Total ;
    G4double* VoxSAF_Total ;
    G4double* VoxAD_Total ;
    G4double* VoxS_Total ;
    G4double* VoxH_Total ;
    G4double* VoxE_Total ;

    //std::vector<unsigned int> getCNID() const { return CNID;}
    //G4double* getCopyNumberMassSize() const { return CopyNumberMassSize;}
    //std::map<unsigned int,G4String> GetCopyNumberRegionNameMap() const { return CopyNumberRegionNameMap;}

    //G4double* VoxChosenVariablevariance;
    //G4double* VoxChosenVariableStandardDeviation;
    //G4double* VoxChosenVariableMean;
    //G4double* VoxChosenVariableTotal;

/*
    G4double* VoxAF2_Total ;
    G4double* VoxAF_Mean ;
    G4double* VoxAF2_Mean ;
    G4double* VoxAF_Var ;
    G4double* VoxAF_SDev ;
    G4double* VoxAF_RelS_D ;

    G4double* VoxSAF2_Total ;
    G4double* VoxSAF_Mean ;
    G4double* VoxSAF2_Mean ;
    G4double* VoxSAF_Var ;
    G4double* VoxSAF_SDev ;
    G4double* VoxSAF_RelS_D ;

    G4double* VoxAD2_Total ;
    G4double* VoxAD_Mean ;
    G4double* VoxAD2_Mean ;
    G4double* VoxAD_Var ;
    G4double* VoxAD_SDev ;
    G4double* VoxAD_RelS_D ;

    G4double* VoxS2_Total ;
    G4double* VoxS_Mean ;
    G4double* VoxS2_Mean ;
    G4double* VoxS_Var ;
    G4double* VoxS_SDev ;
    G4double* VoxS_RelS_D ;

    G4double* VoxH2_Total ;
    G4double* VoxH_Mean ;
    G4double* VoxH2_Mean ;
    G4double* VoxH_Var ;
    G4double* VoxH_SDev ;
    G4double* VoxH_RelS_D ;

    G4double* VoxE2_Total ;
    G4double* VoxE_Mean ;
    G4double* VoxE2_Mean ;
    G4double* VoxE_Var ;
    G4double* VoxE_SDev ;
    G4double* VoxE_RelS_D ;
*/
    std::map<unsigned int, G4double> VoxAFCte ;
    std::map<unsigned int, G4double> VoxSAFCte ;
    std::map<unsigned int, G4double> VoxADCte ;
    std::map<unsigned int, G4double> VoxSCte ;
    std::map<unsigned int, G4double> VoxHCte ;
    std::map<unsigned int, G4double> VoxECte ;

    void ReadThreadVoxelResultFile(G4String);
    void VoxelQuantitiesCalculation();
    void GenerateVoxelsResultFiles();

protected:

};
#endif
