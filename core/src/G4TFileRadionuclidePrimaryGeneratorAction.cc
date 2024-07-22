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

#include "G4TFileRadionuclidePrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4TVolumeConstruction.hh"

#include "G4RunManager.hh"
//#include "G4ios.hh"
#include <iostream>
#include "G4PhysicalConstants.hh"

G4TFileRadionuclidePrimaryGeneratorAction::G4TFileRadionuclidePrimaryGeneratorAction(){


    particleGun = new G4ParticleGun(1);
    GunInitialize();

    //particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(NewRankSourceParticlesNamesValues[DataID]);
    //if(particleDefinition == nullptr){ G4String msg = "Particle name (" + NewRankSourceParticlesNamesValues[DataID] + ") in rank/thread (" + std::to_string(DataID) + ") not found "; G4Exception("Source Data", "1", FatalErrorInArgument, msg.c_str());}
    //particleGun->SetParticleDefinition(particleDefinition);
    particleGun->SetNumberOfParticles(1);

    //std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n************** The primary generator Action initialization... \n" << NewRankSourceParticlesNamesValues[DataID] << std::endl;

}

G4TFileRadionuclidePrimaryGeneratorAction::~G4TFileRadionuclidePrimaryGeneratorAction()
{
    //G4MUTEXDESTROY(mutex);
    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;

    delete particleGun;
}

void G4TFileRadionuclidePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //particleGun->SetParticlePosition(G4ThreeVector(0, 0, 0));
    //particleGun->SetParticleEnergy(1.);
    //particleGun->SetParticleMomentumDirection(G4ParticleMomentum(0, 0, 1.));

    //// For MPI mode, For MT, This is called in the Constructor
    //if(EvInc == 0){
    //    //std::cout << "begin of SourceInitialization()" << std::endl;
    //    SourceInitialization();
    //    particleGun->SetNumberOfParticles(1);
    //    if(EnergyTypeNum != 5){particleGun->SetParticleDefinition(particleDefinition);}
    //    EvInc++;
    //}

    particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(ParNameList[EnergyListInc]));

    //GenerateEventsParticle();
    GenerateEventsEnergy();
    GenerateEventsPosition();
    GenerateEventsMomentumDirection();

    particleGun->SetParticlePosition(G4ThreeVector(X, Y, Z));
    particleGun->SetParticleEnergy(ENERGY);
    particleGun->SetParticleMomentumDirection(G4ParticleMomentum(XMOMD, YMOMD, ZMOMD));

    //std::cout << EvInc << " ENERGY=" << ENERGY << " X=" << X << " Y=" << Y << " Z=" << Z << " XMOMD=" << XMOMD << " YMOMD=" << YMOMD << " ZMOMD=" << ZMOMD << " " << NewRankSourceRegionsBoxDimValues[DataID] << " " << particleDefinition->GetParticleName() << std::endl;

    //if(WriteSourceDataToFiles == 1){SaveGeneratedDataToFiles();}

    particleGun->GeneratePrimaryVertex(anEvent);
    //std::cout << " particleGun->GeneratePrimaryVertex(anEvent) " << std::endl;

}
