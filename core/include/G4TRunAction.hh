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

#ifndef G4TRunAction_h
#define G4TRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4TMessenger.hh"
//#include "G4TOutputText.hh"
#include "G4TResultCalculation.hh"


#ifdef G4MPI_USE
//#include "G4MPImanager.hh"
#include "mpi.h"
#endif


//class G4TOutputText;
class G4TVolumeConstruction;
class G4TPrimaryGeneratorAction;

class G4TRunAction : public G4UserRunAction
{
public:

    G4TRunAction();
    ~G4TRunAction();

private:

    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

    void CreateThreadVoxelsResultsFiles();

    void CreateSimulationDataFile();
    void TestResults();

    // for GDML, Text, Construct geometry

public:

    //void pushEventAbsEnergy(std::map<G4String,G4double>);
    void pushVolumeEnergy(std::map<G4String,G4double>);
    void pushVolumeEnergy2(std::map<G4String,G4double>);
    void pushVolumeHitsNumber(std::map<G4String, unsigned int>);

    // For Voxel Geometry

    void pushVoxelEnergy(std::map<unsigned int,G4double>);
    void pushVoxelEnergy2(std::map<unsigned int,G4double>);
    void pushVoxelHitsNumber(std::map<unsigned int,unsigned int>);

    void FillRegionStepHits(G4String, G4double);
    void FillVoxelStepHits(unsigned int, G4double);
    void FillVoxelLenghts(unsigned int, G4double);
    void FillRegionLenghts(G4String, G4double);

    //G4int getNumberOfFissionNeutrons() const {return NumberOfFissionNeutrons;}
    //G4int getNumberOfFission() const {return NumberOfFission;}
    //G4int getNumberOfCapture() const {return NumberOfCapture;}
    //G4double getKeffective() const {return Keffective;}
    //void setKeffective(G4double nn ){Keffective = nn;};

    void CountParticleProductionByNeutron(G4String);
    void SetPoxEneOfFluxParticle(G4double ,G4double ,G4double, G4double ,G4double ,G4double ,G4double );
    void CountFission();
    void CountAbsorption();
    void IntializeBatchData();
    void SetPoxEneForBatch(G4double ,G4double ,G4double,G4double ,G4double ,G4double ,G4double );


    void CountOpticalInteractions(G4String,G4String);
    void ShowOpticalPhotonInteractionData();

private:

    G4String GenerateCrossSectionGraph, ExecutionMode, OneOrMultiSimulations;

    void useTime(G4int);
    void WriteMacroscopicCrossSection();
    void CreateThreadRegionResultFile();


    std::vector<G4String>       OrgansNameVector;
    std::map<G4String,G4double> OrganNameMassMap;
    std::map<G4String,G4double> OrganNameDensityMap;

#ifdef G4MULTITHREADED
    G4ThreadLocal static std::ofstream CriticalityFile;

    G4ThreadLocal static unsigned long long int*    FluxParticlePosX   ;
    G4ThreadLocal static unsigned long long int*    FluxParticlePosY   ;
    G4ThreadLocal static unsigned long long int*    FluxParticlePosZ   ;
    G4ThreadLocal static unsigned long long *       FluxParticleEnergy ;

    G4ThreadLocal static std::map<G4int,std::map<unsigned int,G4double>>   VoxelsFluence ;
    G4ThreadLocal static std::map<G4int,std::map<unsigned int,G4double>>   VoxelsED_Total ;
    G4ThreadLocal static std::map<G4int,std::map<unsigned int,G4double>>   VoxelsED2_Total ;
    G4ThreadLocal static std::map<G4int,std::map<unsigned int,unsigned long long int>>     VoxelsNOfValues ;

    G4ThreadLocal static std::map<G4int,std::map<G4String,unsigned long long int>> NOfValues;
    G4ThreadLocal static std::map<G4int,std::map<G4String,G4double>> ED_Total ;
    G4ThreadLocal static std::map<G4int,std::map<G4String,G4double>> Fluence ;
    G4ThreadLocal static std::map<G4int,std::map<G4String,G4double>> ED2_Total ;
    G4ThreadLocal static std::map<G4int,std::map<G4String,G4int>> ParticleProduction;

    //G4ThreadLocal static G4double Keffective;
    //G4ThreadLocal static G4int NumberOfFissionNeutrons;
    //G4ThreadLocal static G4int NumberOfFission;
    //G4ThreadLocal static G4int NumberOfCapture;
    G4ThreadLocal static G4int rank;
    G4ThreadLocal static G4int thread;
    G4ThreadLocal static G4int DataID;
    G4ThreadLocal static G4int CurrentBatch;
    G4ThreadLocal static G4int NumberOfRanksThreads;
    G4ThreadLocal static G4int TotalEventNumber;
    G4ThreadLocal static G4double EnergyEmittedPerThread;
    G4ThreadLocal static G4double ParticleSourceEnergy;
    G4ThreadLocal static G4double ExecutionTimeInMin;
    G4ThreadLocal static G4double OneEventExecutionTimeInMs;
    G4ThreadLocal static std::chrono::steady_clock::time_point start;
    G4ThreadLocal static std::chrono::steady_clock::time_point end ;
    G4ThreadLocal static std::map<G4int,std::vector<G4double>> NewNeutronFissionEnergyList;
    G4ThreadLocal static std::map<G4int,std::vector<G4ThreeVector>> NewNeutronFissionPositionsList;
    G4ThreadLocal static std::map<G4int,std::vector<G4ThreeVector>> NewNeutronFissionMomDirecsList;

#else
    G4ThreadLocal std::ofstream CriticalityFile;
    unsigned long long int*    FluxParticlePosX   ;
    unsigned long long int*    FluxParticlePosY   ;
    unsigned long long int*    FluxParticlePosZ   ;
    unsigned long long *       FluxParticleEnergy ;
    std::map<G4int,std::vector<G4double>> NewNeutronFissionEnergyList;
    std::map<G4int,std::vector<G4ThreeVector>> NewNeutronFissionPositionsList;
    std::map<G4int,std::vector<G4ThreeVector>> NewNeutronFissionMomDirecsList;
    std::map<G4int,std::map<unsigned int,G4double>>   VoxelsED_Total ;
    std::map<G4int,std::map<unsigned int,G4double>>   VoxelsFluence ;
    std::map<G4int,std::map<unsigned int,G4double>>   VoxelsED2_Total ;
    std::map<G4int,std::map<unsigned int,unsigned long long int>>     VoxelsNOfValues ;

    std::map<G4int,std::map<G4String,unsigned long long int>> NOfValues;
    std::map<G4int,std::map<G4String,G4double>> ED_Total ;
    std::map<G4int,std::map<G4String,G4double>> Fluence ;
    std::map<G4int,std::map<G4String,G4double>> ED2_Total ;

    //G4double Keffective;
    //G4int NumberOfFissionNeutrons, NumberOfFission, NumberOfCapture;
    G4int rank, thread, DataID, CurrentBatch, NumberOfRanksThreads, TotalEventNumber;
    G4double EnergyEmittedPerThread,ParticleSourceEnergy, ExecutionTimeInMin, OneEventExecutionTimeInMs;
    std::chrono::steady_clock::time_point start, end ;
    std::map<G4int,std::map<G4String,G4int>> ParticleProduction;
#endif
    //G4ThreadLocal static unsigned int StepCN;
    //G4ThreadLocal static G4double StepEnergy;
    //G4ThreadLocal static G4String StepRegion;

protected:

#ifdef G4MPI_USE
    //G4MPImanager* g4MPI1 ;
#endif


};
#endif
