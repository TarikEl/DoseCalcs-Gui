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

#include "G4TDirectToFilesPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4TVolumeConstruction.hh"

#include "G4RunManager.hh"
//#include "G4ios.hh"
#include <iostream>
#include "G4PhysicalConstants.hh"

G4TDirectToFilesPrimaryGeneratorAction::G4TDirectToFilesPrimaryGeneratorAction(){
    particleGun = new G4ParticleGun(1);
    GunInitialize();
}

G4TDirectToFilesPrimaryGeneratorAction::~G4TDirectToFilesPrimaryGeneratorAction(){}

void G4TDirectToFilesPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //particleGun->SetParticlePosition(G4ThreeVector(0, 0, 0));
    //particleGun->SetParticleEnergy(1.);
    //particleGun->SetParticleMomentumDirection(G4ParticleMomentum(0, 0, 1.));

    if(EvInc == 0){
        //std::cout << "begin of SourceInitialization()" << std::endl;
        //std::cout << "\n\n\n\n\n\n\n\n 1 DataID "<< DataID << " SourceType=" << SourceType<< std::endl;

        SourceInitialization();

        particleGun->SetNumberOfParticles(1);

        if(EnergyTypeNum != 5){particleGun->SetParticleDefinition(particleDefinition);}

        OpenFilesToSaveGeneratedData();

        //std::cout << "\n\n\n\n\n\n\n\n 22222222 DataID "<< DataID << " SourceType=" << SourceType<< std::endl;
        EvInc++;
    }

    if(EnergyTypeNum == 5){particleGun->SetParticleDefinition(particleDefinitionList[ParNameList[EnergyListInc]]);}

    //GenerateEventsParticle();
    GenerateEventsEnergy();
    GenerateEventsPosition();
    GenerateEventsMomentumDirection();
    SaveGeneratedDataToFiles();

    particleGun->SetParticlePosition(G4ThreeVector(X, Y, Z));
    particleGun->SetParticleEnergy(ENERGY);
    particleGun->SetParticleMomentumDirection(G4ParticleMomentum(XMOMD, YMOMD, ZMOMD));

    //std::cout << " Particle=" << particleGun->GetParticleDefinition()->GetParticleName() << " ENERGY=" << ENERGY << " X=" << X << " Y=" << Y << " Z=" << Z << std::endl;

    particleGun->GeneratePrimaryVertex(anEvent);

}
