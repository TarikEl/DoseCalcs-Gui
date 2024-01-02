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

#ifndef G4TPrimaryGeneratorMessenger_h
#define G4TPrimaryGeneratorMessenger_h 1

class G4TVolumeConstruction;
class G4TPrimaryGeneratorAction;
class G4TRunAction;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAnIntegerAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;

#include "G4UImessenger.hh"
#include "globals.hh"
#include <iostream>

class G4TPrimaryGeneratorMessenger: public G4UImessenger
{
public:

    G4TPrimaryGeneratorMessenger( G4TVolumeConstruction* t);

    ~G4TPrimaryGeneratorMessenger();

    void SetNewValue(G4UIcommand* command, G4String newValue);

    void CommandsForPointsGen();
    void CommandsForPrimaryGen();
    G4double UseG4Units(G4String);

private:

    G4TPrimaryGeneratorAction*  myPrimaryAct;
    G4TVolumeConstruction* GeoConst;

    G4UIdirectory*				 sourceDataDir;
    G4UIcmdWithAString*          particle_nameCMD ;
    G4UIcmdWithADoubleAndUnit*   MonoEnergyCMD ;

    G4UIcmdWithAString*          particle_angle_distributionCMD ;
    G4UIcmdWithAString*          particle_energy_distributionCMD ;
    G4UIcmdWithADoubleAndUnit*   GaussDistSDevCMD ;
    G4UIcmdWithADoubleAndUnit*   GaussMeanCMD ;
    G4UIcmdWithADoubleAndUnit*   UniformDistEminCMD ;
    G4UIcmdWithADoubleAndUnit*   UniformDistEmaxCMD ;
    G4UIcmdWithADoubleAndUnit*   RayleighDistEmaxCMD ;
    G4UIcmdWithoutParameter*     GeneratePositionsCMD ;
    G4UIcmdWithoutParameter*     GenerateEnergiesCMD ;
    G4UIcmdWithoutParameter*     GenerateMomDirsCMD ;

    G4UIcmdWithAString*          SourceRegionNameCMD ;
    G4UIcmdWith3VectorAndUnit*   Vector_organ_sourceToGenCMD ;
    G4UIcmdWithoutParameter*          showBoxCMD ;
    G4UIcmdWithoutParameter*          TestPointsPositionsCMD ;

    G4UIcmdWith3VectorAndUnit* 	 box_widthCMD ;
    G4UIcmdWithAString*          source_typeCMD ;
    G4UIcmdWithAString*          SourceSolidCMD;
    G4UIcmdWithAString*          SourceSurfaceCMD ;
    G4UIcmdWithAString*          SourcePlaneCMD ;
    G4UIcmdWithAString*          SourceAxisCMD ;
    G4UIcmdWith3VectorAndUnit*   SourceRotationCMD ;
    G4UIcmdWithADoubleAndUnit*   RadiusCMD ;
    G4UIcmdWithADoubleAndUnit*   RadiusInCMD ;
    G4UIcmdWithADoubleAndUnit*   BeamSDevCMD ;
    G4UIcmdWithADoubleAndUnit*   HalfXCMD ;
    G4UIcmdWithADoubleAndUnit*   HalfYCMD ;
    G4UIcmdWithADoubleAndUnit*   HalfZCMD ;

    G4UIcmdWithADoubleAndUnit*   ThetaMinCMD ;
    G4UIcmdWithADoubleAndUnit*   ThetaMaxCMD ;
    G4UIcmdWithADoubleAndUnit*   PhiMinCMD ;
    G4UIcmdWithADoubleAndUnit*   PhiMaxCMD ;

    G4UIcmdWith3VectorAndUnit*   SourcePositionCMD ;
    G4UIcmdWith3Vector*          SourceRotVector1CMD ;
    G4UIcmdWith3Vector*          SourceRotVector2CMD ;
    G4UIcmdWithAnInteger*        numberOfPointsToGenerateCMD ;
    G4UIcmdWithoutParameter*     UsePETCumulativeActCMD;
    G4UIcmdWithADoubleAndUnit*   ThetaForDirectionCMD ;
    G4UIcmdWithADoubleAndUnit*   PhiForDirectionCMD ;

    G4UIcommand*          SourceTypesCMD;
    G4UIcommand*          SourceEnergyDataCMD;
    G4UIcommand*          SourceMomDirDataCMD;
    G4UIcommand*          GenerateSourceDataCMD;

    G4UIcommand*          SimulatedParticlesNamesCMD;
    G4UIcommand*          GeneratePosForRegionAndBoxDimCMD;
    G4UIcommand*          GenerateEnergiesForCMD;
    G4UIcommand*          GenerateMomDirsForCMD;

};

#endif
