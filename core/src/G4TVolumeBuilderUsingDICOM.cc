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

#include "globals.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
//#include "G4Material.hh"
#include "G4Element.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIcommand.hh"

#include "G4TPhantomParameterisation.hh"
#include "G4TVolumeBuilderUsingDICOM.hh"
#include "G4TDCMZSliceHeader.hh"

#include "G4TVolumeConstruction.hh"
#include "G4RunManager.hh"


extern  size_t* MateIDs; // index of material of each voxel unsigned int* fMateIDs; // index of material of each voxel

G4TVolumeBuilderUsingDICOM::G4TVolumeBuilderUsingDICOM() :

NumOfFiles(0),
//MateIDs(0),
ZSliceHeaderMerged(0),

NVoxelX(0),
NVoxelY(0),
NVoxelZ(0),
VoxelHalfDimX(0),
VoxelHalfDimY(0),
VoxelHalfDimZ(0)

{

    const G4TVolumeConstruction* TConstruction1 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    OriginalMaterials = TConstruction1->getDcmMaterialsVector();

	Readg4dcmFiles();

    NVoxelX = ZSliceHeaderMerged->GetNoVoxelX();
    NVoxelY = ZSliceHeaderMerged->GetNoVoxelY();
    NVoxelZ = ZSliceHeaderMerged->GetNoVoxelZ();

    VoxelHalfDimX = ZSliceHeaderMerged->GetVoxelHalfX();
    VoxelHalfDimY = ZSliceHeaderMerged->GetVoxelHalfY();
    VoxelHalfDimZ = ZSliceHeaderMerged->GetVoxelHalfZ();

    //G4cout << " NVoxelX        " << NVoxelX <<    " VoxelHalfDimX    " << VoxelHalfDimX <<G4endl;
    //G4cout << " NVoxelY        " << NVoxelY <<    " VoxelHalfDimY    " << VoxelHalfDimY <<G4endl;
    //G4cout << " NVoxelZ    	   " << NVoxelZ <<    " VoxelHalfDimZ    " << VoxelHalfDimZ <<G4endl;
    //G4cout << " TotalPixels    " << NVoxelX*NVoxelY*NVoxelZ <<  G4endl;

    //--- Place it on the world
    //G4double fOffsetX = (ZSliceHeaderMerged->GetMaxX() + ZSliceHeaderMerged->GetMinX() ) /2.;
    //G4double fOffsetY = (ZSliceHeaderMerged->GetMaxY() + ZSliceHeaderMerged->GetMinY() ) /2.;
    //G4double fOffsetZ = (ZSliceHeaderMerged->GetMaxZ() + ZSliceHeaderMerged->GetMinZ() ) /2.;
    //posCentreVoxels = G4ThreeVector(fOffsetX,fOffsetY,fOffsetZ);
}

G4TVolumeBuilderUsingDICOM::~G4TVolumeBuilderUsingDICOM(){

}



// called from Readg4dcmFiles()
#include <glob.h>
std::vector<G4String> G4TVolumeBuilderUsingDICOM::GetDicomFilesNames(G4String DirPath)
{

	G4cout << "--> The g4dcm DICOM Files : " << G4endl;

	glob_t glob_result;
	glob( DirPath , GLOB_TILDE , NULL , &glob_result );
	for ( unsigned int i=0 ; i < glob_result.gl_pathc ; ++i ){
		G4String vvv = glob_result.gl_pathv[i];
		if(vvv.contains(".g4dcm")){
			DicomFilePathsVector.push_back(glob_result.gl_pathv[i]);
			G4cout << glob_result.gl_pathv[i] << G4endl;
		}
	}

	return DicomFilePathsVector;
}

// this function open each file .g4dcm and remplis the fMateIDs[fNoFiles*nVoxels] with the indices readed from each file continusly
// it fill VoxelsIndexesMap[Z][Y][x] fMatsId[] and materials[]
// called from Readg4dcmFiles() for each file (slice) to get its data and save it in fZSliceHeaders vector
void G4TVolumeBuilderUsingDICOM::Readg4dcmFiles()
{

	//G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@ void G4TVolumeBuilderUsingDICOM::Readg4dcmFiles(const G4String& fname)  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

	G4String fname;

    const G4TVolumeConstruction* TConstruction1 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    DicomFilesDirPath = TConstruction1->getDicomDataDirPath()+"/*";

	GetDicomFilesNames(DicomFilesDirPath);
	NumOfFiles = DicomFilePathsVector.size();

	for(G4int i = 0; i < DicomFilePathsVector.size(); i++ ) {

		//--- Read one data file
		fname = DicomFilePathsVector[i];

		std::ifstream fin(fname.c_str(), std::ios_base::in);
		if( !fin.is_open() ) {
			G4Exception("G4TVolumeBuilderUsingDICOM::Readg4dcmFiles", "", FatalErrorInArgument, G4String("File not found " + fname ).c_str());
		}

		//----- Define density differences (maximum density difference to create a new material)
		char* part = getenv( "DICOM_CHANGE_MATERIAL_DENSITY" );
		G4double densityDiff = -1.;
		if( part ) densityDiff = G4UIcommand::ConvertToDouble(part);

		if( densityDiff != -1. ) {
			for( unsigned int ii = 0; ii < OriginalMaterials.size(); ii++ ){
				DensityDiffs[ii] = densityDiff; //currently all materials with
				// same difference
			}
        }
        else {
			if( Materials.size() == 0 ) { // do it only for first slice
				for( unsigned int ii = 0; ii < OriginalMaterials.size(); ii++ ){
					Materials.push_back( OriginalMaterials[ii] );
				}
			}
		}

		//----- Read data header

		// fZSliceHeaders is an object vector that will contain all sliceHeader created from all .g4dcm
		G4TDCMZSliceHeader* sliceHeader = new G4TDCMZSliceHeader( fin );
		ZSliceHeaders.push_back( sliceHeader );

		//----- Read material indices

		G4int nVoxels = sliceHeader->GetNoVoxels();

		//--- If first slice, initiliaze fMateIDs
		if( ZSliceHeaders.size() == 1 ) {
			MateIDs = new size_t[NumOfFiles*nVoxels];
		}

		// my functions
		//number_of_voxels_in_xyz = fNoFiles*nVoxels;
		//

		// read the matIDs from file and save it in fMateIDs[voxelCopyNo] array that contains the matIDs of all slices consecutively
		unsigned int mateID;

		// number of voxels from previously read slices
		// fill fMateIDs[voxelCopyNo]

		G4int voxelCopyNo = (ZSliceHeaders.size()-1)*nVoxels; // pass high of number of voxels in all previous slices to atach the matIDs to new slice voxels
		for( G4int ii = 0; ii < nVoxels; ii++, voxelCopyNo++ ){
			fin >> mateID;
			MateIDs[voxelCopyNo] = mateID;
		}

		// each slice with a new Zindex
		G4int Xindex = 0 ,Yindex = 0,  Zindex = ZSliceHeaders.size()-1 ;// pass high of number of voxels in all previous slices to atach the matIDs to new slice voxels
		while (Yindex < sliceHeader->GetNoVoxelY()){
			Xindex = 0;
			while (Xindex < sliceHeader->GetNoVoxelX()){
				VoxelsIndexesMap[Zindex][Yindex].push_back(Xindex); // we begin from z=0, y=0, x=0
				//G4cout << Zindex << " " << Yindex << " "  << Xindex << G4endl;
				Xindex++;
			}
			Yindex++;
		}

		//----- Read material densities and build new materials if two voxels have same material but its density is in a different density interval
		// (size of density intervals defined by densityDiff)

		G4double density;

		// number of voxels from previously read slices

		// to here, i have the matIDs of each pixel, but i have to build material with it's density , then first we have to get the density of each pixel
		voxelCopyNo = (ZSliceHeaders.size()-1)*nVoxels;
		for( G4int ii = 0; ii < nVoxels; ii++, voxelCopyNo++ ){
			fin >> density;

			//-- Get material from list of original materials
			mateID = MateIDs[voxelCopyNo];
			G4Material* mateOrig  = OriginalMaterials[mateID];

			//-- Get density bin: middle point of the bin in which the current density is included
			G4String newMateName = mateOrig->GetName();
			float densityBin = 0. ;
			if( densityDiff != -1.) {
				densityBin = DensityDiffs[mateID] * (G4int(density/DensityDiffs[mateID])+0.5);
				//-- Build the new material name
				newMateName += G4UIcommand::ConvertToString(densityBin);
			}

			//-- Look if a material with this name is already builded or created (because the previous voxels can be already in this density bin)
			unsigned int im;
			for( im = 0; im < Materials.size(); im++ ){
				if( Materials[im]->GetName() == newMateName ) {
					break;
				}
			}

			//-- If material is already created use index of this material and not out of index
			if( im != Materials.size() ) {
				MateIDs[voxelCopyNo] = im;
				//-- else, create the material

			} else { // if im is out of index of fMaterials,
				if( densityDiff != -1.) {
					Materials.push_back( GenerateMaterialWithDensity( mateOrig, newMateName , densityBin) );
					MateIDs[voxelCopyNo] = Materials.size()-1;
				} else {
					G4cerr << " im " << im << " < " << Materials.size() << " name " << newMateName << G4endl;
					G4Exception("G4TVolumeBuilderUsingDICOM::Readg4dcmFiles", "", FatalErrorInArgument, "Wrong index in material"); //it should never reach here
				}
			}
		}
	}

	//----- Merge data headers
	MergeZSliceHeaders();

}

// called from Readg4dcmFiles()
G4Material* G4TVolumeBuilderUsingDICOM::GenerateMaterialWithDensity(G4Material* OriginalMaterial, G4String NewMatName , G4double density){

	//G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@ G4TVolumeBuilderUsingDICOM::GenerateMaterialWithDensity()  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

	//----- Copy original material, but with new density
	G4int nelem = OriginalMaterial->GetNumberOfElements();
	G4Material* NewMaterial = new G4Material( NewMatName, density*g/cm3, nelem, kStateUndefined, STP_Temperature );

	for( G4int ii = 0; ii < nelem; ii++ ){
		G4double frac = OriginalMaterial->GetFractionVector()[ii];
		G4Element* element = const_cast<G4Element*>(OriginalMaterial->GetElement(ii));
		NewMaterial->AddElement( element, frac );
	}

	return NewMaterial;
}

// for the mergedfeaders by += ,there is the gol of operator fuction += , it add the NVoxelZ and VoxelHalfDimZ (and others !!! ) on each sliceHeader's
// called from Readg4dcmFiles() to merge the slices created by Readg4dcmFiles() and listed in fZSliceHeaders vector
void G4TVolumeBuilderUsingDICOM::MergeZSliceHeaders()
{

	//G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@ DicomDetectorConstruction::MergeZSliceHeaders()  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

	//----- Images must have the same dimension ...
	ZSliceHeaderMerged = new G4TDCMZSliceHeader( *ZSliceHeaders[0] );
	for( unsigned int ii = 1; ii < ZSliceHeaders.size(); ii++ ) {
		*ZSliceHeaderMerged += *ZSliceHeaders[ii];
	}
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
