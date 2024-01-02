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

#ifndef G4TPhysicsMessenger_h
#define G4TPhysicsMessenger_h 1

//class G4TUserPhysicsList;
class G4TVolumeConstruction;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3VectorAndUnit;


#include "G4UImessenger.hh"
#include "globals.hh"
#include <iostream>

class G4TPhysicsMessenger: public G4UImessenger
{
public:
    G4TPhysicsMessenger(G4TVolumeConstruction*);

    ~G4TPhysicsMessenger();

    void SetNewValue(G4UIcommand* command, G4String newValue);
    G4double UseG4Units(G4String);

    //G4String GetCurrentValue(G4UIcommand * command);

    void CommandsForPhysics();

private:

    G4TVolumeConstruction*  myUserPhysics;


    G4UIdirectory*				 physicsDataDir;
    G4UIcmdWithAString*          particle_PhysicsCMD ;
    G4UIcmdWithADoubleAndUnit*          particle_Cuts_energyCMD ;
    G4UIcmdWithADoubleAndUnit*          particle_Cuts_DistanceCMD ;

    G4UIcommand* PhysicsDataCMD;
    G4UIcommand* CutsDataCMD;
    G4UIcommand* EnergyRangeDataCMD;
    G4UIcommand* SetEnergiesForCrossSectionCMD;

    G4UIcmdWithAString* PhotoElectricEffectModelCMD;
    //G4UIcmdWithAString* PolarizedPhotoElectricEffectModelCMD;
    G4UIcmdWithAString* ComptonScatteringModelCMD;
    //G4UIcmdWithAString* PolarizedComptonModelCMD;
    G4UIcmdWithAString* GammaConversionModelCMD;
    //G4UIcmdWithAString* PolarizedGammaConversionModelCMD;
    G4UIcmdWithAString* RayleighScatteringModelCMD;
    //G4UIcmdWithAString* GammaConversionToMuonModelCMD;

    G4UIcmdWithAString* ElectronIonisationModelCMD;
    G4UIcmdWithAString* ElectronBremModelCMD;
    G4UIcmdWithAString* HadronIonisationModelCMD;
    G4UIcmdWithAString* IonIonisationModelCMD;


};

#endif
