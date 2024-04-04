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
// TETRunAction.hh
// \file   MRCP_GEANT4/Internal/include/TETRunAction.hh
// \author Haegin Han
//

#ifndef TETRunAction_h
#define TETRunAction_h 1

#include <ostream>
#include <fstream>
#include <map>

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4UserRunAction.hh""
#include "G4SystemOfUnits.hh"

#include "TETRun.hh"

// *********************************************************************
// The main function of this UserRunAction class is to produce the result
// data and print them.
// -- GenergateRun: Generate TETRun class which will calculate the sum of
//                  energy deposition.
// -- BeginOfRunAction: Set the RunManager to print the progress at the
//                      interval of 10%.
// -- EndOfRunAction: Print the run result by G4cout and std::ofstream.
//  └-- PrintResult: Method to print the result.
// *********************************************************************

class TETRunAction : public G4UserRunAction
{
public:
    TETRunAction();
	virtual ~TETRunAction();

public:
	virtual G4Run* GenerateRun();
	virtual void BeginOfRunAction(const G4Run*);
	virtual void EndOfRunAction(const G4Run*);

    void PrintResult();
    void FillRegionLenghts(G4String, G4double);

private:

    G4String GenerateCrossSectionGraph, ExecutionMode, OneOrMultiSimulations;

    void useTime(G4int);
    void WriteMacroscopicCrossSection();
    void CreateThreadRegionResultFile();

    std::vector<G4String>       OrgansNameVector;
    std::map<G4String,G4double> OrganNameMassMap;
    std::map<G4String,G4double> OrganNameDensityMap;

#ifdef G4MULTITHREADED
    G4ThreadLocal static std::map<G4int,std::map<unsigned int,G4double>>   VoxelsED_Total ;
    G4ThreadLocal static std::map<G4int,std::map<unsigned int,G4double>>   VoxelsED2_Total ;
    G4ThreadLocal static std::map<G4int,std::map<unsigned int,unsigned long long int>>     VoxelsNOfValues ;

    G4ThreadLocal static std::map<G4int,std::map<G4String,unsigned long long int>> NOfValues;
    G4ThreadLocal static std::map<G4int,std::map<G4String,G4double>> ED_Total ;
    G4ThreadLocal static std::map<G4int,std::map<G4String,G4double>> ED2_Total ;
    G4ThreadLocal static std::map<G4int,std::map<G4String,G4double>> Fluence ;

    G4ThreadLocal static G4int rank;
    G4ThreadLocal static G4int thread;
    G4ThreadLocal static G4int DataID;
    G4ThreadLocal static G4int EventIndex;
    G4ThreadLocal static G4int NumberOfRanksThreads;
    G4ThreadLocal static G4int TotalEventNumber;
    G4ThreadLocal static G4double EnergyEmittedPerThread;
    G4ThreadLocal static G4double ParticleSourceEnergy;
    G4ThreadLocal static G4double ExecutionTimeInMin;
    G4ThreadLocal static G4double OneEventExecutionTimeInMs;
    G4ThreadLocal static std::chrono::steady_clock::time_point start;
    G4ThreadLocal static std::chrono::steady_clock::time_point end ;
    G4ThreadLocal static TETRun* fRun;

#else
    std::map<G4int,std::map<unsigned int,G4double>>   VoxelsED_Total ;
    std::map<G4int,std::map<unsigned int,G4double>>   VoxelsED2_Total ;
    std::map<G4int,std::map<unsigned int,unsigned long long int>>     VoxelsNOfValues ;

    std::map<G4int,std::map<G4String,unsigned long long int>> NOfValues;
    std::map<G4int,std::map<G4String,G4double>> ED_Total ;
    std::map<G4int,std::map<G4String,G4double>> ED2_Total ;
    std::map<G4int,std::map<G4String,G4double>> Fluence ;

    G4int rank, thread, DataID, EventIndex, NumberOfRanksThreads, TotalEventNumber;
    G4double EnergyEmittedPerThread,ParticleSourceEnergy, ExecutionTimeInMin, OneEventExecutionTimeInMs;
    std::chrono::steady_clock::time_point start, end ;
    TETRun* fRun;

#endif
    //G4ThreadLocal static unsigned int StepCN;
    //G4ThreadLocal static G4double StepEnergy;
    //G4ThreadLocal static G4String StepRegion;

    std::map<G4int, G4String> MatIDNameMap;

};

#endif





