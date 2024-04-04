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

#ifndef G4TNeutronCriticalityPrimaryGeneratorAction_h
#define G4TNeutronCriticalityPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
//#include <vector>
//#include "G4ThreeVector.hh"
#include"G4TPrimaryGeneratorMessenger.hh"
#include"G4TPrimaryGeneratorMethods.hh"
//#include "G4TOutputText.hh"

#ifdef G4MPI_USE
#include "mpi.h"
//#include "G4MPImanager.hh"
#endif

class G4ParticleGun;
//class G4GeneralParticleSource;
class G4Event;

class G4TNeutronCriticalityPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction, public G4TPrimaryGeneratorMethods {

public:
    G4TNeutronCriticalityPrimaryGeneratorAction();
    ~G4TNeutronCriticalityPrimaryGeneratorAction();

    void GeneratePrimaries(G4Event* anEvent);
/*
#ifdef G4MULTITHREADED
    G4ThreadLocal static G4double TotalEmittedEnergy;
    G4ThreadLocal static unsigned int EvInc;
    G4ThreadLocal static unsigned int DataID;
    G4ThreadLocal static double* Energies;
    G4ThreadLocal static G4ThreeVector* Positions;
    G4ThreadLocal static G4ThreeVector* MomDirecs;
#else
    G4double TotalEmittedEnergy;
    unsigned int EvInc;
    unsigned int DataID;
    double* Energies;
    G4ThreeVector* Positions;
    G4ThreeVector* MomDirecs;
#endif
    // called from endOfRun for each thread to get the emmited energy because the energy distribution can can change its value, then we calc it
*/
    //G4double getEmittedEnergy() const { return TotalEmittedEnergy;}

private:

    G4ParticleGun* particleGun;

    //void GetEventsData();

#ifdef G4MPI_USE
    //G4MPImanager* g4MPI1 ;
#endif

};
#endif


