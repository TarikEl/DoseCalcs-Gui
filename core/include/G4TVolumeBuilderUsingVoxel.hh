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

#ifndef INCLUDE_G4TVOLUMEBUILDERUSINGVOXEL_HH_
#define INCLUDE_G4TVOLUMEBUILDERUSINGVOXEL_HH_

#include "G4ThreeVector.hh"

#include <set>
#include <vector>
#include <map>

class G4Material;
class G4Box;
class G4VPhysicalVolume;
class G4LogicalVolume;

struct matInfo
{
    G4double fSumdens;
    G4int fNvoxels;
    G4int fId;
};

class G4TVolumeBuilderUsingVoxel {

public:

    G4TVolumeBuilderUsingVoxel();
    ~G4TVolumeBuilderUsingVoxel();

    void StartGettingDicomData();

    void InitialisationOfMaterials();
    // create the original materials

    void ReadVoxelHeaderAndMatIDs(std::ifstream&);
    // read the DICOM files describing the phantom
    void ReadVoxelDensities( std::ifstream& fin );
    void ReadVoxelActivities( std::ifstream& fin );

    G4Material* BuildMaterialWithChangingDensity( const G4Material* origMate, float density, G4String newMateName );
    // build a new material if the density of the voxel is different
    // to the other voxels

    std::vector<G4double> GetVoxelsDataXYZ(){

        std::vector<G4double> nn;
        nn.push_back(fNVoxelX);
        nn.push_back(fNVoxelY);
        nn.push_back(fNVoxelZ);
        nn.push_back(fVoxelHalfDimX);
        nn.push_back(fVoxelHalfDimY);
        nn.push_back(fVoxelHalfDimZ);
        return nn;
    }
    std::vector<G4double> GetPETVoxelsDataXYZ(){

        std::vector<G4double> nn;
        nn.push_back(fPETNVoxelX);
        nn.push_back(fPETNVoxelY);
        nn.push_back(fPETNVoxelZ);
        nn.push_back(fPETVoxelHalfDimX);
        nn.push_back(fPETVoxelHalfDimY);
        nn.push_back(fPETVoxelHalfDimZ);
        return nn;
    }
    size_t* GetMateIds(){
        return fMateIDs;
    }
    double* GetActivities(){
        return fActivities;
    }
    std::vector<G4Material*> GetMaterials(){
        return fMaterials;
    }

    void setDicomOutTextName(G4String n){DicomDatafileName = n;}
    void setDicomType(G4String n){DicomType = n;}

protected:
    G4Material* fAir;

    // World ...
    G4Box* fWorld_solid;
    G4LogicalVolume* fWorld_logic;
    G4VPhysicalVolume* fWorld_phys;

    G4Box* fContainer_solid;
    G4LogicalVolume* fContainer_logic;
    G4VPhysicalVolume* fContainer_phys;

    G4int fNoFiles; // number of DICOM files
    std::vector<G4Material*> fOriginalMaterials;  // list of original materials
    std::vector<G4Material*> fMaterials;
    // list of new materials created to distinguish different density
    //  voxels that have the same original materials
    size_t* fMateIDs; // index of material of each voxel
    //unsigned int* fMateIDs; // index of material of each voxel

    double* fActivities;

    std::map<G4int,G4double> fDensityDiffs;
    // Density difference to distinguish material for each original
    // material (by index)

    G4int fNVoxelX, fNVoxelY, fNVoxelZ;
    G4double fVoxelHalfDimX, fVoxelHalfDimY, fVoxelHalfDimZ;
    G4int fPETNVoxelX, fPETNVoxelY, fPETNVoxelZ;
    G4double fPETVoxelHalfDimX, fPETVoxelHalfDimY, fPETVoxelHalfDimZ;

    G4double fMinX,fMinY,fMinZ; // minimum extension of voxels (position wall)
    G4double fMaxX,fMaxY,fMaxZ; // maximum extension of voxels (position wall)

    std::map<G4int,G4Material*> thePhantomMaterialsOriginal;
    // map numberOfMaterial to G4Material. They are the list of materials as
    // built from .geom file

    G4String DicomDatafileName;
    G4String DicomType;


};

#endif /* INCLUDE_G4TVOLUMEBUILDERUSINGVOXEL_HH_ */
