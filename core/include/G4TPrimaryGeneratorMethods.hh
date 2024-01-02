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

#ifndef G4TPrimaryGeneratorMethods_h
#define G4TPrimaryGeneratorMethods_h 1

#include "G4Navigator.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleMomentum.hh"
#include "G4AutoLock.hh"

#ifdef G4MPI_USE
//#include "G4MPImanager.hh"
#include "mpi.h"
#endif

#include <fstream>

class G4TPrimaryGeneratorMethods {

public:
    G4TPrimaryGeneratorMethods();
    ~G4TPrimaryGeneratorMethods();

    G4double getEmittedEnergy() const { return TotalEmittedEnergy;}
    G4int getEvInc() const { return EvInc;}

protected:

#ifdef G4MULTITHREADED
    G4ThreadLocal static unsigned int EvInc;
    G4ThreadLocal static unsigned int ParticleID;
    G4ThreadLocal static unsigned int EnergyListInc;
    G4ThreadLocal static int DataID;
    G4ThreadLocal static int WriteSourceDataToFiles;

    G4ThreadLocal static G4double TotalEmittedEnergy;
    G4ThreadLocal static G4double ENERGY;
    G4ThreadLocal static G4double* EnergyList;
    G4ThreadLocal static G4ThreeVector* PositionsList;
    G4ThreadLocal static G4ThreeVector* MomDirecsList;
    G4ThreadLocal static std::map<unsigned int,G4ParticleDefinition* > particleDefinitionList;

    G4ThreadLocal static std::map<G4String,std::map<G4String,std::map<int,G4String>>> PosDataFileNames;
    G4ThreadLocal static std::map<G4String,std::map<double,std::map<int,G4String>>> EneDataFileNames;
    G4ThreadLocal static std::map<G4String,std::map<int,G4String>> MomDirDataFileNames;

    G4ThreadLocal static double* PosXList;
    G4ThreadLocal static double* PosYList;
    G4ThreadLocal static double* PosZList;
    G4ThreadLocal static double* MomDirXList;
    G4ThreadLocal static double* MomDirYList;
    G4ThreadLocal static double* MomDirZList;
    G4ThreadLocal static unsigned int* ParNameList;

    G4ThreadLocal static G4double X;
    G4ThreadLocal static G4double Y;
    G4ThreadLocal static G4double Z;
    G4ThreadLocal static G4double XMOMD;
    G4ThreadLocal static G4double YMOMD;
    G4ThreadLocal static G4double ZMOMD;
    G4ThreadLocal static G4double Voxel0PosX;
    G4ThreadLocal static G4double Voxel0PosY;
    G4ThreadLocal static G4double Voxel0PosZ;
    G4ThreadLocal static G4int VoxelsInc;
    G4ThreadLocal static G4int CummNumbInVoxelsInc;
    G4ThreadLocal static G4VPhysicalVolume* WorldPhysicalVolume;
    G4ThreadLocal static G4Navigator* aNavigator;
    G4ThreadLocal static int PositionTypeNum;
    G4ThreadLocal static int EnergyTypeNum;
    G4ThreadLocal static int MomDirTypeNum;
    G4ThreadLocal static G4ParticleDefinition* particleDefinition;
    G4ThreadLocal static std::ofstream PositionFileStream;
    G4ThreadLocal static std::ofstream MomDirFileStream;
    G4ThreadLocal static std::ofstream EneFileStream;
#else
    unsigned int EvInc;
    unsigned int EnergyListInc;
    int DataID;
    G4double TotalEmittedEnergy;
    G4double ENERGY;
    G4double* EnergyList;
    G4double* PositionsList;
    G4double* MomDirecsList;

    double* PosXList;
    double* PosYList;
    double* PosZList;
    double* MomDirXList;
    double* MomDirYList;
    double* MomDirZList;
    unsigned int* ParNameList;
    unsigned int ParticleID;
    std::map<unsigned int,G4ParticleDefinition* > particleDefinitionList;
    std::map<G4String,std::map<G4String,std::map<int,G4String>>> PosDataFileNames;
    std::map<G4String,std::map<double,std::map<int,G4String>>> EneDataFileNames;
    std::map<G4String,std::map<int,G4String>> MomDirDataFileNames;

    G4double X;
    G4double Y;
    G4double Z;
    G4double XMOMD;
    G4double YMOMD;
    G4double ZMOMD;
    G4double Voxel0PosX;
    G4double Voxel0PosY;
    G4double Voxel0PosZ;

    G4int VoxelsInc;
    G4int CummNumbInVoxelsInc;
    G4VPhysicalVolume* WorldPhysicalVolume;
    G4Navigator* aNavigator;
    int PositionTypeNum;
    int EnergyTypeNum;
    int MomDirTypeNum;
    G4ParticleDefinition* particleDefinition;
    std::ofstream PositionFileStream;
    std::ofstream MomDirFileStream;
    std::ofstream EneFileStream;
    int WriteSourceDataToFiles;
#endif

#ifdef G4MPI_USE
    //G4MPImanager* g4MPI1 ;
#endif

    void SourceInitialization();
    void GunInitialize();
    void GetEventsData();

    void GenerateEventsParticle();
    //virtual void GenerateEventsPosition();
    void GenerateEventsPosition();
    void GenerateEventsMomentumDirection();
    void GenerateEventsEnergy();
    void OpenFilesToSaveGeneratedData();
    void CloseFilesToSaveGeneratedData();
    void SaveGeneratedDataToFiles();

};
#endif


