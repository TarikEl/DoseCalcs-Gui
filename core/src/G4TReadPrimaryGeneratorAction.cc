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

#include "G4TReadPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
//#include "G4ParticleTable.hh"
//#include "G4ParticleDefinition.hh"
//#include "G4PhysicalConstants.hh"
#include "G4TVolumeConstruction.hh"
#include "G4TPointDataGeneration.hh"

#include "G4RunManager.hh"
//#include "G4ios.hh"
#include <iostream>
//#include <cmath>
#include "G4AutoLock.hh"

/*
extern  G4String* CopyNumberRegionNameMap;
extern  G4float* CopyNumberXPos;
extern  G4float* CopyNumberYPos;
extern  G4float* CopyNumberZPos;

// Source Data

extern G4String ParticleName;

extern G4String SourceType;
extern G4String SourceRegionName;
extern G4ThreeVector boxCenterPos;
extern G4int NumberOfGenPointsToSave;
extern G4ThreeVector BoxDimGene;
extern G4ThreeVector newOrgPos3Vector;
extern G4ThreeVector SourcePosition;
extern G4ThreeVector SourceRotation;
extern G4String SourceSolid;
extern G4String SourceSurface;
extern G4String SourcePlane;
extern G4String SourceAxis;
extern G4double Radius;
extern G4double RadiusIn;
extern G4double BeamSDev;
extern G4double HalfX, HalfY, HalfZ;
extern G4double ThetaMin;
extern G4double ThetaMax;
extern G4double PhiMin;
extern G4double PhiMax;
extern G4ThreeVector SourceRotVector1;
extern G4ThreeVector SourceRotVector2;

extern G4String EnergyDistribution;
extern G4double GaussSDev;
extern G4double GaussMean;
extern G4double UniformEmin;
extern G4double UniformEmax;
extern G4double RayleighEmax;
extern G4double MonoEnergy;

extern G4String MomDirDistribution;
extern G4double Theta;
extern G4double Phi;

extern G4String GeneratePosFlag ;
extern G4String GenerateEneFlag ;
extern G4String GenerateMomDirFlag ;
extern G4String UseGeneratedData ;
extern G4String ShowBox ;
extern G4String TestPointsPositionsList ;

extern G4String PositionDataFile, EnergyDataFile, MomDirDataFile;

extern G4String ScriptsDirectoryPath, DataDirectoryPath;


// Run And Score
extern G4String MPISimulationNum;
extern G4int BatchsNumber;
extern G4int ThreadsNumber;
extern G4int EventNumberInOneThreads;
extern G4String AccuracyCalculationLevel;
extern G4String organs_to_score ;
extern G4String variable_To_Score ;

// Voxelized Geometry
extern G4double VoxXHalfSize;
extern G4double VoxYHalfSize;
extern G4double VoxZHalfSize;
extern G4int VoxXNumber;
extern G4int VoxYNumber;
extern G4int VoxZNumber;

extern G4int EventsNumPerThreadRank;

extern std::vector<G4String> NewRankSourcePositionDataFiles;
extern std::vector<G4String> NewRankSourceEnergyDataFiles;
extern std::vector<G4String> NewRankSourceMomDirDataFiles;

#ifdef G4MULTITHREADED
G4ThreadLocal G4double G4TReadPrimaryGeneratorAction::TotalEmittedEnergy;
G4ThreadLocal unsigned int G4TReadPrimaryGeneratorAction::EvInc;
G4ThreadLocal unsigned int G4TReadPrimaryGeneratorAction::DataID;
G4ThreadLocal double* G4TReadPrimaryGeneratorAction::EnergyList;
G4ThreadLocal G4ThreeVector* G4TReadPrimaryGeneratorAction::PositionsList;
G4ThreadLocal G4ThreeVector* G4TReadPrimaryGeneratorAction::MomDirecsList;
#endif
*/

namespace
{
G4Mutex	mutex = G4MUTEX_INITIALIZER;
}

extern G4String appBuildDir;

G4TReadPrimaryGeneratorAction::G4TReadPrimaryGeneratorAction(){

    //particleGun = new G4ParticleGun();

    //EvInc = 0;
    //TotalEmittedEnergy = 0. ;

    particleGun = new G4ParticleGun(1);
    GunInitialize();

}



G4TReadPrimaryGeneratorAction::~G4TReadPrimaryGeneratorAction()
{
    G4MUTEXDESTROY(mutex);
    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;

    delete particleGun;

}


// it's called a number of times that you give in beamOn
void G4TReadPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    if(EvInc == 0){
        G4AutoLock l(&mutex);
        GetEventsData();
        l.unlock();
    }

    //G4cout << "\n\n\n\n"<< EvInc <<" PositionsList[EvInc] " << PositionsList[EvInc] << G4endl;

    particleGun->SetNumberOfParticles(1);
    particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(ParticleName));
    particleGun->SetParticlePosition(PositionsList[EvInc]);
    particleGun->SetParticleMomentumDirection(MomDirecsList[EvInc]);
    particleGun->SetParticleEnergy(EnergyList[EvInc]);  // Energy comes from .mac file and from GetEnergyWithDistribution() it stay in geant4 default unit

    particleGun->GeneratePrimaryVertex(anEvent);

    TotalEmittedEnergy += EnergyList[EvInc] ;
    EvInc++;

}

/*
// called from constructor(), to fill the vector of PositionsList sources
void G4TReadPrimaryGeneratorAction::GetEventsData(){

    //G4cout << __FUNCTION__ << G4endl;

    // to read the file and fill the vector of points if it's empty, else it would be filled then we dont need to read and fill ...

    //G4int TotalEventsNumPerThreadRank = G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed();

    G4int rank = 0 , thread = 0, DataID = 0;

#ifdef G4MPI_USE
    //rank = G4MPImanager::GetManager()->GetRank();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(MPISimulationNum == "m"){
        //DataID = G4MPImanager::GetManager()->GetRank();
        MPI_Comm_rank(MPI_COMM_WORLD, &DataID);

    }
#else
    if(G4Threading::IsMultithreadedApplication()){ // normal multiThreaded mode
        thread = G4Threading::G4GetThreadId();
        if(MPISimulationNum == "m"){
            DataID = thread;
        }
    }
#endif

    EnergyList  = new double[EventsNumPerThreadRank];
    PositionsList = new G4ThreeVector[EventsNumPerThreadRank];
    MomDirecsList = new G4ThreeVector[EventsNumPerThreadRank];

    //G4cout << " !!!!!!! If the number of generated points is less than the number of events to be processed, the application will reopens the file and reads the events data" << G4endl;
    std::cout << "\n---Reading " << EventsNumPerThreadRank << " Generated Points PositionsList from file " << NewRankSourcePositionDataFiles[DataID].c_str() << "..."<< std::endl;

    G4double posX,posY,posZ, val;

    //#########################################################################################################################################"

    //std::ifstream fileorganStream(NewRankSourcePositionDataFiles[DataID].c_str());
    std::ifstream PosFileStream;
    PosFileStream.open(NewRankSourcePositionDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary

    if(!PosFileStream.is_open()) {
        std::cout << " Unable to open file : " << NewRankSourcePositionDataFiles[DataID].c_str() << std::endl;
    }
    else {
        G4int Begin, First = 0;
        Begin = (rank*ThreadsNumber*EventsNumPerThreadRank) + (thread*EventsNumPerThreadRank);

        while(First < Begin){ // read just what we need to simulate

            PosFileStream >> posY >> posZ >> posY ;

            if(PosFileStream.peek()==EOF && First < Begin){
                PosFileStream.close(); // , std::ios_base::out | std::ios_base::binary
                PosFileStream.open(NewRankSourcePositionDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
            }

            First++;
        }
        std::cout << "    EventsNumPerThreadRank "<< EventsNumPerThreadRank << " Events Position Data For thread " << thread << " | Begins from line " << Begin << " to " << Begin + EventsNumPerThreadRank << std::endl;

        G4int EventIncre = 0;
        while(EventIncre < EventsNumPerThreadRank){ // read just what we need to simulate

            if(PosFileStream.peek()==EOF && EventIncre < EventsNumPerThreadRank){ // If the number of generated points is less than the number of events to be processed, the application will reopen the file and reads points PositionsList

                //std::cout << "EventIncre < EventsNumPerThreadRank && peek()==EOF " << std::endl;
                PosFileStream.close(); // , std::ios_base::out | std::ios_base::binary
                PosFileStream.open(NewRankSourcePositionDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
                //std::cout  << "##########################"<<EventsNumPerThreadRank << " " <<EventIncre << " " << EnergyList[EventIncre] << " " << PositionsList[EventIncre] << " " << MomDirecsList[EventIncre] << std::endl;
                //continue;
            }

            PosFileStream >> posX >> posY >> posZ ;
            PositionsList[EventIncre] = G4ThreeVector( posX , posY , posZ);

            //if( G4Threading::G4GetThreadId() == 3 && rank == 3){
            //  std::cout << EventsNumPerThreadRank << " " << EventIncre << " " << PositionsList[EventIncre] << std::endl;
            //}


            //std::cout <<EventsNumPerThreadRank << " " <<EventIncre  << " Pos " << PositionsList[EventIncre] << std::endl;

            EventIncre++;
        }

        PosFileStream.close();

        //std::cout << "End of PositionsList File reading"<< std::endl;
    }

    //###########################################posY##############################################################################################"

    std::ifstream EneFileStream;
    EneFileStream.open(NewRankSourceEnergyDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary

    if(!EneFileStream.is_open()) {
        std::cout << " Unable to open file : : "<< NewRankSourceEnergyDataFiles[DataID].c_str() << std::endl;
    }
    else {

        G4int Begin, First = 0;
        Begin = (rank*ThreadsNumber*EventsNumPerThreadRank) + (thread*EventsNumPerThreadRank);

        while(First < Begin){ // read just what we need to simulate

            EneFileStream >> val ;

            if(EneFileStream.peek()==EOF && First < Begin){
                EneFileStream.close();
                EneFileStream.open(NewRankSourceEnergyDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
            }

            First++;
        }
        std::cout << "    EventsNumPerThreadRank "<< EventsNumPerThreadRank << " Events Energy Data For thread " << thread << " | Begins from line " << Begin << " to " << Begin + EventsNumPerThreadRank << std::endl;

        G4int EventIncre = 0;
        while(EventIncre < EventsNumPerThreadRank){ // read just what we need to simulate

            if(EneFileStream.peek()==EOF && EventIncre < EventsNumPerThreadRank){ // If the number of generated points is less than the number of events to be processed, the application will reopen the file and reads points PositionsList

                //std::cout << "EventIncre < EventsNumPerThreadRank && peek()==EOF " << std::endl;
                EneFileStream.close();
                EneFileStream.open(NewRankSourceEnergyDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
                //std::cout  << "##########################"<<EventsNumPerThreadRank << " " <<EventIncre << " " << EnergyList[EventIncre] << " " << PositionsList[EventIncre] << " " << MomDirecsList[EventIncre] << std::endl;
                //continue;
            }

            EneFileStream >> val;
            EnergyList[EventIncre] = val ;

            //if( G4Threading::G4GetThreadId() == 3  && rank == 3){
            //  std::cout << EventsNumPerThreadRank << " " << EventIncre << " " << EnergyList[EventIncre] << std::endl;
            //}

            //std::cout <<EventsNumPerThreadRank << " " <<EventIncre << " " << EnergyList[EventIncre] << " " << PositionsList[EventIncre] << " " << MomDirecsList[EventIncre] << std::endl;

            EventIncre++;
        }

        EneFileStream.close();

        //std::cout << "End of Energy File reading"<< std::endl;
    }

    //#########################################################################################################################################"

    std::ifstream MomDirFileStream;
    MomDirFileStream.open(NewRankSourceMomDirDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary

    if(!MomDirFileStream.is_open()) {
        std::cout << " Unable to open file : "<< NewRankSourceMomDirDataFiles[DataID].c_str() << std::endl;
    }
    else {

        G4int Begin, First = 0;
        Begin = (rank*ThreadsNumber*EventsNumPerThreadRank) + (thread*EventsNumPerThreadRank);

        while(First < Begin){ // read just what we need to simulate

            MomDirFileStream >> posX >> posY >> posZ ;

            if(MomDirFileStream.peek()==EOF && First < Begin){
                MomDirFileStream.close();
                MomDirFileStream.open(NewRankSourceMomDirDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
            }

            First++;
        }

        std::cout << "    EventsNumPerThreadRank "<< EventsNumPerThreadRank << " Events Momentum Direction Data For thread " << thread << " | Begins from line " << Begin << " to " << Begin + EventsNumPerThreadRank << std::endl;

        G4int EventIncre = 0;
        while(EventIncre < EventsNumPerThreadRank){ // read just what we need to simulate

            if(MomDirFileStream.peek()==EOF && EventIncre < EventsNumPerThreadRank){ // If the number of generated points is less than the number of events to be processed, the application will reopen the file and reads points PositionsList

                //std::cout << "EventIncre < EventsNumPerThreadRank && peek()==EOF " << std::endl;
                MomDirFileStream.close();
                MomDirFileStream.open(NewRankSourceMomDirDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
                //std::cout  <<EventsNumPerThreadRank << " " <<EventIncre << " " << EnergyList[EventIncre] << " " << PositionsList[EventIncre] << " " << MomDirecsList[EventIncre] << std::endl;
                //continue;
            }

            MomDirFileStream >> posX >> posY >> posZ ;
            MomDirecsList[EventIncre] = G4ThreeVector( posX , posY , posZ);

            //if( G4Threading::G4GetThreadId() == 3  && rank == 3){
            //  std::cout << EventsNumPerThreadRank << " " << EventIncre << " " << MomDirecsList[EventIncre] << std::endl;
            //}

            EventIncre++;
        }

        MomDirFileStream.close();

        //std::cout << "End of MomDir File reading"<< std::endl;
    }

    std::cout << "---Getting Events Data has done. Now is time for Events simulation" << std::endl;

    //for ( int ff = 0 ; ff < EventsNumPerThreadRank ; ff++ ) {
    //std::cout << "EnergyList[ff] " << EnergyList[ff]<< std::endl;
    //std::cout << "PositionsList[ff] " << PositionsList[ff]<< std::endl;
    //std::cout << "MomDirecsList[ff] " << MomDirecsList[ff]<< std::endl;
    //}

}

 */
