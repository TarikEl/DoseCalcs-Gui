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

#ifndef G4TModifiedPrimaryGeneratorAction_h
#define G4TModifiedPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4TPrimaryGeneratorMethods.hh"
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

class G4ParticleGun;
//class G4GeneralParticleSource;
class G4Event;

class G4TModifiedPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction, public G4TPrimaryGeneratorMethods {

public:
    G4TModifiedPrimaryGeneratorAction();
    ~G4TModifiedPrimaryGeneratorAction();

    void GeneratePrimaries(G4Event* anEvent);

private:

    G4ParticleGun* particleGun;
};
#endif
