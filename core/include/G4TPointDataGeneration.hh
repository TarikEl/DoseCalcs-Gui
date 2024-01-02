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

#ifndef G4TPointDataGeneration_h
#define G4TPointDataGeneration_h 1

#include "G4VPhysicalVolume.hh"
#include "globals.hh"
#include "G4Navigator.hh"

extern std::vector<G4String> SourceRegionsNamesToBeIgnoredValues;

extern  G4String* CopyNumberRegionNameMap;
extern  G4float* CopyNumberXPos;
extern  G4float* CopyNumberYPos;
extern  G4float* CopyNumberZPos;

// Source Data

extern G4String particleName;

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
extern G4String TestPointsPositions ;

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
extern G4bool VOX_USE;
extern G4double VoxXHalfSize;
extern G4double VoxYHalfSize;
extern G4double VoxZHalfSize;
extern G4int VoxXNumber;
extern G4int VoxYNumber;
extern G4int VoxZNumber;
extern G4bool UseDicomCumAct;

class G4VPhysicalVolume;
class G4TPointDataGeneration
{

public:

    G4TPointDataGeneration();
    ~G4TPointDataGeneration();

    void ConstructBoxVolume();
    void RemoveBoxVolume(G4LogicalVolume*);

    G4ThreeVector GetRegionPosition(G4String);

    void SetSourceInputs();

    G4bool GenerateEventsPosition();
    G4bool GenerateEventsMomentumDirection();
    G4bool GenerateEventsEnergy();

    //G4ThreeVector boxCenterPos;
    //G4ThreeVector getboxCenterPos(){return boxCenterPos; }// called from the construction::construct()

private:

    //std::map<G4String,std::map<unsigned int,G4ThreeVector>> RegionCopyNumberPositionMap;

    G4VPhysicalVolume*  WorldPhysicalVolume;
    unsigned int TotalVoxelsNumber;
    G4ThreeVector VoxelsSize;

    /*
    G4bool UseDicomCumAct;
G4int NumberOfGenPointsToSave;
    G4ThreeVector BoxDimGene;
    G4String SourceType;
    G4ThreeVector SourcePosition;
    G4ThreeVector SourceRotVector1;
    G4ThreeVector SourceRotVector2;
    G4ThreeVector SourceRotation;

    G4String SourceSolid;
    G4String SourcePlane;
    G4String SourceSurface;
    G4String SourceAxis;
    G4double Radius;
    G4double HalfX, HalfY, HalfZ;
    G4double RadiusIn ;
    G4double BeamSDev ;

    G4double ThetaMin;
    G4double ThetaMax;
    G4double PhiMin;
    G4double PhiMax;

    G4String SourceRegionName;
    G4String particleName;
    G4double MonoEnergy;
    G4String EnergyDistribution;
    G4String angleDistribution;
    G4double GaussSDev;
    G4double GaussMean;
    G4double UniformEmin;
    G4double UniformEmax;
    G4double RayleighEmax;

    G4double Theta;
    G4double Phi;

    G4String PositionDataFile;
    G4String EnergyDataFile;
    G4String MomDirDataFile;

    G4VPhysicalVolume*  MotherPhysVolume;
*/

};

#endif

