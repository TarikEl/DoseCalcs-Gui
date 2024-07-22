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

#include "G4TNeutronCriticalityPrimaryGeneratorAction.hh"
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


namespace
{
G4Mutex	mutex = G4MUTEX_INITIALIZER;
}

extern G4String appBuildDir;

G4TNeutronCriticalityPrimaryGeneratorAction::G4TNeutronCriticalityPrimaryGeneratorAction(){

    //particleGun = new G4ParticleGun();

    //EvInc = 0;
    //TotalEmittedEnergy = 0. ;

    particleGun = new G4ParticleGun(1);

    GunInitialize();

    particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(NewRankSourceParticlesNamesValues[DataID]);
    if(particleDefinition == nullptr){ G4String msg = "Particle name (" + NewRankSourceParticlesNamesValues[DataID] + ") in rank/thread (" + std::to_string(DataID) + ") not found "; G4Exception("Source Data", "1", FatalErrorInArgument, msg.c_str());}
    particleGun->SetParticleDefinition(particleDefinition);
    particleGun->SetNumberOfParticles(1);

    if(UseGeneratedData == "read"){
        G4AutoLock l(&mutex);
        GetEventDataFromCriticalityDataFile();
        l.unlock();
    }else{
        FillDataMapsForNeutronCriticality();
    }

}



G4TNeutronCriticalityPrimaryGeneratorAction::~G4TNeutronCriticalityPrimaryGeneratorAction()
{
    G4MUTEXDESTROY(mutex);
    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;

    delete particleGun;

}


// it's called a number of times that you give in beamOn
void G4TNeutronCriticalityPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    //// For MPI mode, For MT, This is called in the Constructor
    //if(EvInc == 0){
    //    //G4cout << " First Event: Initialization"<< G4endl;
    //    SourceInitialization();
    //    particleGun->SetNumberOfParticles(1);
    //    if(UseGeneratedData == "read"){
    //        G4AutoLock l(&mutex);
    //        GetEventDataFromCriticalityDataFile();
    //        l.unlock();
    //    }else{
    //        FillDataMapsForNeutronCriticality();
    //    }
    //}





    //G4cout << NumberOfEventInBatch << " "<< EvInc << G4endl;
       //               1000-1000=0 should not simulated because the array data dont have the index 1000
    //if((NumberOfEventInBatch-EvInc) == 2){
    //    G4AutoLock l(&mutex);
    //    GetEventDataFromCriticalityDataFile();
    //    //FillDataMapsForNeutronCriticality();
    //    l.unlock();
    //    EvInc = 0;
    //}

    //G4cout << "NumberOfEventInBatch-EvInc " << NumberOfEventInBatch-EvInc << " EvInc "<< EvInc <<" PositionsList[EvInc] " << PositionsList[EvInc] <<" EnergyList[EvInc] " << EnergyList[EvInc] << G4endl;

    //particleGun->SetNumberOfParticles(1);
    particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(ParticleName));

    particleGun->SetParticleEnergy(EnergyListForCriticality[DataID][EvInc]);  // Energy comes from .mac file and from GetEnergyWithDistribution() it stay in geant4 default unit
    particleGun->SetParticlePosition(PositionsListForCriticality[DataID][EvInc]);
    particleGun->SetParticleMomentumDirection(MomDirecsListForCriticality[DataID][EvInc]);
    TotalEmittedEnergy += EnergyListForCriticality[DataID][EvInc] ;

    //particleGun->SetParticlePosition(PositionsList[EvInc]);
    //particleGun->SetParticleMomentumDirection(MomDirecsList[EvInc]);
    //particleGun->SetParticleEnergy(EnergyList[EvInc]);  // Energy comes from .mac file and from GetEnergyWithDistribution() it stay in geant4 default unit
    //TotalEmittedEnergy += EnergyList[EvInc] ;

    //G4cout << " Execute " << G4endl;

    particleGun->GeneratePrimaryVertex(anEvent);

    EvInc++;

}
