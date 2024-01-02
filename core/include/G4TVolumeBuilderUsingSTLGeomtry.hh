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

#ifndef G4TVolumeBuilderUsingSTLGeomtry_h
#define G4TVolumeBuilderUsingSTLGeomtry_h 1

#include "G4VPhysicalVolume.hh"
#include "globals.hh"
#include "G4Material.hh"
//#include "G4TOutputText.hh"

class G4VPhysicalVolume;
class G4TVolumeBuilderUsingSTLGeomtry
{

public:

    G4TVolumeBuilderUsingSTLGeomtry();
    ~G4TVolumeBuilderUsingSTLGeomtry();

    G4String OrganMassFilePath;
    G4String OrganSourceName;
    G4bool GenerateWithBox;
    G4bool generatePoints;
    G4bool UseOrganMassFilePath;
    void setGeometryData();

    void ConstructBox( G4ThreeVector , G4double , G4double, G4double);
    void RemoveBoxVolumeAndSetOrgVol();


    void ReadOrganMassFile(G4String);
    void ReadSTLFile(G4String);
    G4VPhysicalVolume* ReadGeometryFile(G4String);
    G4VPhysicalVolume* ReadSeparatedGeometryFile(G4String);



    G4ThreeVector boxCenterPos;
    G4ThreeVector getboxCenterPos(){return boxCenterPos; }// called from the construction::construct()

    // to fill the vector by organsNamesCreated and get it by the builder finally
    std::vector<G4String> createdOrgans;
    std::map<G4String,G4ThreeVector> createdPositions;
    std::map<G4String,G4ThreeVector> createdRotations;
    std::map<G4String,G4double> createdOrganDensity;
    std::map<G4String,G4double> createdOrganMass;
    // that called from G4TVolumeConstruction to get the vector filled by names of organ created and ...
    std::vector<G4String> getOrganCreated(){ return createdOrgans;}
    std::map<G4String,G4double> getOrganDensityCreated(){return createdOrganDensity; }
    std::map<G4String,G4double> getOrganMassCreated(){return createdOrganMass; }
    std::map< G4String , G4ThreeVector > getOrganPositionCreated(){return createdPositions;}
    std::map< G4String , G4ThreeVector > getOrganRotationCreated(){return createdRotations;}
    G4VPhysicalVolume* GetPhantomPhysicalVol(){ return WorldPhysicalVol;}
    G4VPhysicalVolume* GetMotherPhysicalsVol(){ return MotherPhysVolume;}

    G4VPhysicalVolume* PhysicalBoxVolume;
    G4LogicalVolume* LogicalBoxVol;
    G4VPhysicalVolume* MotherPhysVolume;

protected:

    std::map<G4String,G4VPhysicalVolume* > physicalMotherVolumes;

    G4LogicalVolume* LogicalBoxVolume;

    G4VPhysicalVolume* orgPhyVol;
    G4VPhysicalVolume* organMother;
    G4VPhysicalVolume* WorldPhysicalVol;

private:

    //G4TOutputText* out;

};

#endif

