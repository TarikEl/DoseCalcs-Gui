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

#ifndef INCLUDE_G4TVOLUMEBUILDERUSINGDICOM_HH_
#define INCLUDE_G4TVOLUMEBUILDERUSINGDICOM_HH_

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Colour.hh"
#include "G4Material.hh"


#include "G4TDCMZSliceHeader.hh"

class G4TVolumeBuilderUsingDICOM
{
public:

	G4TVolumeBuilderUsingDICOM();
	~G4TVolumeBuilderUsingDICOM();

	void Readg4dcmFiles();
	void ConstructPhantom();
	void MergeZSliceHeaders();
	std::vector<G4String> GetDicomFilesNames(G4String);

    G4Material* GenerateMaterialWithDensity(G4Material*, G4String , G4double);
	void InitialisationOfMaterials();
	void ReadDicomInputDataFile(G4String);


	void GenerateMaterialsFromDicomDataFile();

	std::map<G4String,std::vector<G4double>> material_composition;
	std::map<G4String,G4double> material_density;

	//std::map<G4int,std::map<G4String,std::vector<G4double>>> MediumMaterial_data_readed; // # MediumMaterial ID , MediumName, H     C     N     O    Na    Mg     P     S    Cl     K    Ca    Fe     I
	//std::map<G4int,G4Material*> MediumsMaterials; // # MediumMaterial ID , Material
	//std::map<G4String,G4Material*> OrgansMaterials; // # OrganMaterial ID , Material
	// # a MediumMaterial can be used from a number of organss
	//std::map<G4String,std::map<G4int,std::vector<G4double>>> OrganMaterial_data_readed; // # organID , OrganName , MediumMaterial indice , density , mass

	G4double VoxelHalfDimX, VoxelHalfDimY, VoxelHalfDimZ, RefVoxelHalfDimX, RefVoxelHalfDimY, RefVoxelHalfDimZ ;
	G4double NVoxelX, NVoxelY, NVoxelZ, RefNVoxelX, RefNVoxelY, RefNVoxelZ  ;
    std::vector<G4double> GetVoxelsDataXYZ(){

		std::vector<G4double> nn;
		nn.push_back(NVoxelX);
		nn.push_back(NVoxelY);
		nn.push_back(NVoxelZ);
		nn.push_back(VoxelHalfDimX);
		nn.push_back(VoxelHalfDimY);
		nn.push_back(VoxelHalfDimZ);
		return nn;
    }

    //size_t* GetMateIds(){ return MateIDs;}

    std::vector<G4Material*> GetMaterials(){
        return Materials;
    }


	G4int NumOfFiles;

	G4String DicomFilesDirPath;
	G4String DicomDataFile;
	G4String DicomInputDataFile;


private:

	std::vector<G4String> DicomFilePathsVector;

	std::vector<G4TDCMZSliceHeader*> ZSliceHeaders;
	std::map<G4int,G4double> DensityDiffs;
	std::vector<G4Material*> Materials;
	std::vector<G4Material*> OriginalMaterials;  // list of original materials
	// list of new materials created to distinguish different density voxels that have the same original materials

    //size_t* MateIDs; // index of material of each voxel unsigned int* fMateIDs; // index of material of each voxel

	G4Material* AirMaterial;

	std::map<G4int,std::map<G4int,std::vector<G4double>>> VoxelsIndexesMap;

	G4TDCMZSliceHeader* ZSliceHeaderMerged ;

};

#endif /* INCLUDE_G4TVOLUMEBUILDERUSINGDICOM_HH_ */
