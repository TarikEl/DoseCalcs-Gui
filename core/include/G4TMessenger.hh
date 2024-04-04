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

#ifndef G4TMessenger_h
#define G4TMessenger_h 1

class G4TVolumeConstruction;
//class G4TPrimaryGeneratorAction;
//class G4TRunAction;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3VectorAndUnit;


#include "G4UImessenger.hh"
#include "globals.hh"
#include <iostream>

class G4TMessenger: public G4UImessenger
{
public:
    G4TMessenger(G4TVolumeConstruction* myUsrPhtm);

    ~G4TMessenger();

    void SetNewValue(G4UIcommand* command, G4String newValue);

    void CommandsForRunAndScore();
    void CommandsForAnalysis();
    G4double UseG4Units(G4String);


private:

    G4TVolumeConstruction*           VolumeConst;

    G4UIcmdWithAnInteger*        number_of_threadsCMD ;
    G4UIcmdWithAString*          AccuracyCalculationLevelCMD ;
    G4UIcmdWithAString*			 SimNumOnRanksCMD;
    G4UIcmdWithoutParameter*     GenerateVoxelsResulsCMD;
    G4UIcommand*     SimulationIntExtNeutDet;

    G4UIcmdWithAString*			 ResultDirectoryPathCMD;

    G4UIdirectory*				 RunAndScoreDir;
    G4UIcommand*                 RegionToScoreCMD;
    G4UIcommand*                 QuantitiesToScoreCMD;

    G4UIdirectory*				 analysisDataDir;

    /*
    G4UIcmdWithAnInteger*        number_of_batchsCMD ;
    G4UIcmdWithAnInteger*        number_events_per_threadCMD ;
    G4UIcmdWithAString*          Graphs_DataCMD ;
    G4UIcmdWithAString*          Compare_typeCMD ;
    G4UIcmdWithAString*          Graphs_ExtCMD ;
    G4UIcmdWithAString*          Ref_File_PathCMD ;
    G4UIcmdWithAString*			 Ref_NameCMD;
    G4UIcmdWithAnInteger*        Num_Of_EneCMD ;
    G4UIcmdWithAnInteger*        Num_Of_RefEneCMD ;
    G4UIcmdWithAString*          Organs_to_scoreCMD ;
    G4UIcmdWithAString*          Variable_To_ScoreCMD ;
    G4UIcmdWithAString*			 EventsPositionHistogramCMD;
    G4UIcmdWithAString*			 EventsEnergyHistogramCMD;
    G4UIcmdWithAString*			 EventsMomDirHistogramCMD;


    G4UIcmdWithAnInteger*        SliceIDCMD ;
    G4UIcmdWithAString*			 SliceFor2DGraphCMD;
    G4UIcmdWithAString*			 BeamAxisCMD;

    */


    G4UIcommand*                 GenerateRelativeErrGraphCMD;
    G4UIcmdWithoutParameter*     GenerateRelativeSDevGraphCMD;
    G4UIcommand*                 RadioTracerDataCMD;
    G4UIcommand*                 RadioTracerBiokineticCMD;
    G4UIcommand*                 RadiationFactorsCMD;
    G4UIcommand*                 QuantitiesUnitsCMD;
    G4UIcommand*                 TissueFactorsCMD;
    G4UIcmdWithoutParameter*     EventsDataHistogramsCMD;
    G4UIcommand*                 GenerateRegionsVariableGraph;
    G4UIcommand*                 GenerateSelfCrossGraphs;
    G4UIcommand*                 GenerateVoxelizedHist;
    G4UIcommand*                 SetGraphsParameters;

};

#endif
