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
// $Id: G4TDCMZSliceHeader.cc 101109 2016-11-07 08:14:53Z gcosmo $
//
/// \file G4TDCMZSliceHeader.cc
/// \brief Implementation of the G4TDCMZSliceHeader class
//

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4GeometryTolerance.hh"
#include "G4NistManager.hh"

#include "G4TDCMZSliceHeader.hh"



G4TDCMZSliceHeader::~G4TDCMZSliceHeader()
{

}


// this constuctor is used when we need we create and instance from this class and we send it the name of file when we will save the data readed from DICOM
// called from DicomHandler::readFile() to create a new ZSlice header from the .dcm file
G4TDCMZSliceHeader::G4TDCMZSliceHeader(const G4String& fname)
:   fNoVoxelX(0),fNoVoxelY(0),fNoVoxelZ(0),
	fMinX(0),fMinY(0),fMinZ(0),
	fMaxX(0),fMaxY(0),fMaxZ(0),
	fFilename(fname), // this is the name of file to write data into it its "1" or "2" and merge with suffix .g4dcm
	fSliceLocation(0)
{
	//G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  G4TDCMZSliceHeader constructor to initialise the variables and set the file name .g4dcm to fill it after  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

}


// called from DicomDetectorConstruction::MergeZSliceHeaders() , G4TDCMZSliceMerged::CheckSlices(),
// we call it when we need to see if the Z slice headers can be merged or not or configure them to merge.
G4TDCMZSliceHeader::G4TDCMZSliceHeader( const G4TDCMZSliceHeader& rhs )
{
	//G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  G4TDCMZSliceHeader  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

	fNoVoxelX = rhs.GetNoVoxelX();
	fNoVoxelY = rhs.GetNoVoxelY();
	fNoVoxelZ = rhs.GetNoVoxelZ();
	fMinX = rhs.GetMinX();
	fMaxX = rhs.GetMaxX();
	fMinY = rhs.GetMinY();
	fMaxY = rhs.GetMaxY();
	fMinZ = rhs.GetMinZ();
	fMaxZ = rhs.GetMaxZ();
	fMaterialNames = rhs.GetMaterialNames();
	fFilename = rhs.fFilename;
	fValues = rhs.fValues;
	fMateIDs = rhs.fMateIDs;
	fSliceLocation = rhs.fSliceLocation;
}


// this constructor is used when we need to read data from the text(or bin) file .g4dcm (resulting when we read the .dcm file).
// Sending the file input and when reading we construct an object of ZSliceHeader. this object by void DicomHandler::StoreData(std::ofstream& *.g4dcm);
// called from  DicomDetectorConstruction::ReadPhantomDataFile()
G4TDCMZSliceHeader::G4TDCMZSliceHeader( std::ifstream& fin )
{

    G4cout << " Reading Z Slice data from it's .g4dcm file with G4TDCMZSliceHeader::G4TDCMZSliceHeader" << G4endl;

	//----- Read material indices and names
	G4int nmate;
	G4String mateindex;
	G4String matename;
	fin >> nmate;

#ifdef G4VERBOSE
    //G4cout << " G4TDCMZSliceHeader reading number of materials " << nmate << G4endl;
#endif

	for( G4int im = 0; im < nmate; im++ ){
		fin >> mateindex >> matename;
#ifdef G4VERBOSE
		//G4cout << " G4TDCMZSliceHeader reading material " << im << " : "
		//        << mateindex << "  " << matename << G4endl;
#endif

		if( ! CheckMaterialExists( matename ) ) {
			G4Exception("G4TDCMZSliceHeader::G4TDCMZSliceHeader", "A material is found in file that is not built in the C++ code", FatalErrorInArgument, matename.c_str());
		}

		fMaterialNames.push_back(matename);
	}

	//----- Read number of voxels
	fin >> fNoVoxelX >> fNoVoxelY >> fNoVoxelZ;

#ifdef G4VERBOSE
    //G4cout << " Number of voxels x:" << fNoVoxelX << " y:" << fNoVoxelY << " z:" << fNoVoxelZ << G4endl;
#endif

	//----- Read minimal and maximal extensions (= walls of phantom)
	fin >> fMinX >> fMaxX;
	fin >> fMinY >> fMaxY;
	fin >> fMinZ >> fMaxZ;

#ifdef G4VERBOSE
    //G4cout  << " Extension in X:  " << fMinX << " " << fMaxX << G4endl
            //<< " Extension in Y:  " << fMinY << " " << fMaxY << G4endl
            //<< " Extension in Z:  " << fMinZ << " " << fMaxZ << G4endl;
#endif

	fSliceLocation = 0.5*(fMinZ + fMaxZ);
}



// called from the constructor G4TDCMZSliceHeader( std::ifstream& fin ) that mean it's will be called when an object of this class is created with this constructor in the classes
// it call G4MaterialTable::GetName() the material are builded by c++ code in DicomConstruction class, and this function is called after building materials needed in the DicomConstruction class
G4bool G4TDCMZSliceHeader::CheckMaterialExists( const G4String& mateName )
{

	//G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4bool G4TDCMZSliceHeader::CheckMaterialExists( const G4String& mateName )e  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

	const G4MaterialTable* matTab = G4Material::GetMaterialTable();
	std::vector<G4Material*>::const_iterator matite;
	for( matite = matTab->begin(); matite != matTab->end(); ++matite ) {
		if( (*matite)->GetName() == mateName ) { return true; }
	}

	G4Material* g4mate = G4NistManager::Instance()->FindOrBuildMaterial(mateName);
	if( g4mate ) {
		return false;
	} else {
		return true;
	}

}




void G4TDCMZSliceHeader::operator+=( const G4TDCMZSliceHeader& rhs ) {
	//G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function void G4TDCMZSliceHeader::operator+=( const G4TDCMZSliceHeader& rhs )  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

	*this = *this + rhs;
}



G4TDCMZSliceHeader G4TDCMZSliceHeader::operator+( const G4TDCMZSliceHeader& rhs )
{
	//G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMZSliceHeader G4TDCMZSliceHeader::operator+( const G4TDCMZSliceHeader& rhs )  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

	//----- Check that both slices has the same dimensions
	if( fNoVoxelX != rhs.GetNoVoxelX() || fNoVoxelY != rhs.GetNoVoxelY() ) {
		G4cerr << "G4TDCMZSliceHeader error adding two slice headers:\
        !!! Different number of voxels: "
				<< "  X= " << fNoVoxelX << " =? " << rhs.GetNoVoxelX()
				<< "  Y=  " << fNoVoxelY << " =? " << rhs.GetNoVoxelY()
				<< "  Z=  " << fNoVoxelZ << " =? " << rhs.GetNoVoxelZ()
				<< G4endl;
		G4Exception("G4TDCMZSliceHeader::G4TDCMZSliceHeader",
				"",FatalErrorInArgument,"");
	}
	//----- Check that both slices has the same extensions
	if( fMinX != rhs.GetMinX() || fMaxX != rhs.GetMaxX() || fMinY != rhs.GetMinY() || fMaxY != rhs.GetMaxY() ) {
		G4cerr << "G4TDCMZSliceHeader error adding two slice headers:\
        !!! Different extensions: "
				<< "  Xmin= " << fMinX << " =? " << rhs.GetMinX()
				<< "  Xmax= " << fMaxX << " =? " << rhs.GetMaxX()
				<< "  Ymin= " << fMinY << " =? " << rhs.GetMinY()
				<< "  Ymax= " << fMaxY << " =? " << rhs.GetMaxY()
				<< G4endl;
		G4Exception("G4TDCMZSliceHeader::operator+","",
				FatalErrorInArgument,"");
	}

	//----- Check that both slices have the same materials
	std::vector<G4String> fMaterialNames2 = rhs.GetMaterialNames();
	if( fMaterialNames.size() != fMaterialNames2.size() ) {
		G4cerr << "G4TDCMZSliceHeader error adding two slice headers:\
        !!! Different number of materials: " << fMaterialNames.size() << " =? "
		<< fMaterialNames2.size() << G4endl;
		G4Exception("G4TDCMZSliceHeader::operator+","",
				FatalErrorInArgument,"");
	}
	for( unsigned int ii = 0; ii < fMaterialNames.size(); ii++ ) {
		if( fMaterialNames[ii] != fMaterialNames2[ii] ) {
			G4cerr << "G4TDCMZSliceHeader error adding two slice headers:\
            !!! Different material number " << ii << " : " 
			<< fMaterialNames[ii] << " =? "
			<< fMaterialNames2[ii] << G4endl;
			G4Exception("G4TDCMZSliceHeader::operator+","",
					FatalErrorInArgument,"");
		}
	}

	//----- Check that the slices are contiguous in Z
	if( std::fabs( fMinZ - rhs.GetMaxZ() ) > G4GeometryTolerance::GetInstance()->GetRadialTolerance() && std::fabs( fMaxZ - rhs.GetMinZ() ) >
	G4GeometryTolerance::GetInstance()->GetRadialTolerance() ){
		G4cerr << "G4TDCMZSliceHeader error adding two slice headers: !!!\
        Slices are not contiguous in Z "
				<< "  Zmin= " << fMinZ << " & " << rhs.GetMinZ()
				<< "  Zmax= " << fMaxZ << " & " << rhs.GetMaxZ()
				<< G4endl;
		G4Exception("G4TDCMZSliceHeader::operator+","",
				FatalErrorInArgument,"");
	}

	//----- Build slice header copying first one
	G4TDCMZSliceHeader temp( *this );

	//----- Add data from second slice header
	temp.SetMinZ( std::min( fMinZ, rhs.GetMinZ() ) );
	temp.SetMaxZ( std::max( fMaxZ, rhs.GetMaxZ() ) );
	temp.SetNoVoxelZ( fNoVoxelZ + rhs.GetNoVoxelZ() );

	return temp;
}



// executed after function DicomHandle::readFile() and after G4TDCMZSliceMerged::CheckSlices() to write the data saved in it from DicomHandle::readFile() to .g4dcm
// it add the data matIds and after density_values to the file output .g4dcm to complete file
// called from G4TDCMZSliceMerged::CheckSlices()
// it call ReadDataFromFile() ,
void G4TDCMZSliceHeader::DumpToFile()
{

	//G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function void G4TDCMZSliceHeader::DumpToFile()  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

	// fFilename is the name of file to write data into it its "1" or "2" and merged with the suffix .g4dcm
    G4cout << " Dumping Z Slice data to " << fFilename << " with G4TDCMZSliceHeader::DumpToFile" << G4endl;
	//sleep(5);

	//  May seen counter-intuitive (dumping to file you are reading from), but
	//  the reason for this is modification slice spacing
	if(fMateIDs.size() == 0 || fValues.size() == 0) { ReadDataFromFile(); }


	std::ofstream out;
	out.open(fFilename.c_str());


	if(!out) {
		G4String descript = "G4TDCMZSliceHeader::DumpToFile: could not open " + fFilename;
		G4Exception(descript.c_str(),"", FatalException, "");
	}

	out << fMaterialNames.size() << std::endl;
	for(unsigned int i = 0; i < fMaterialNames.size(); ++i) {
		out << i << " " << fMaterialNames.at(i) << std::endl;
	}

	out << fNoVoxelX << " " << fNoVoxelY << " " << fNoVoxelZ << std::endl;
	out << fMinX << " " << fMaxX << std::endl;
	out << fMinY << " " << fMaxY << std::endl;
	out << fMinZ << " " << fMaxZ << std::endl;

	for(unsigned int i = 0; i < fMateIDs.size(); ++i) { Print(out,fMateIDs.at(i)," "); }
	for(unsigned int i = 0; i < fValues.size(); ++i) { Print(out,fValues.at(i)," ",6); }

	out.close();

}



// called from G4TDCMZSliceHeader::DumpToFile() when we need to open the file .g4dcm created before. this object sliceHeader and fill it with
// it's data (filled when we readData and storData in DicomHandler class)
// it call
void G4TDCMZSliceHeader::ReadDataFromFile()
{

	//G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function void G4TDCMZSliceHeader::ReadDataFromFile()  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

	std::ifstream in;
	in.open(fFilename.c_str());

	if(!in) {
		G4String descript = "G4TDCMZSliceHeader::DumpToFile: could not open " + fFilename;
		G4Exception(descript.c_str(),"", FatalException, "");
	}
	else {

		G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  G4TDCMZSliceHeader::ReadDataFromFile() to read data from .g4dcm file  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

	}

	G4int nMaterials;
	in >> nMaterials;

	fMaterialNames.resize(nMaterials,"");
	for(G4int i = 0; i < nMaterials; ++i) {
		G4String str1, str2;
		in >> str1 >> str2;
		if(!IsInteger(str1)) {
			G4String descript = "String : " + str1 + " supposed to be integer";
			G4Exception("G4TDCMZSliceHeader::ReadDataFromFile - error in \
       formatting: missing material index","", FatalException,descript.c_str());
		}
		G4int index = G4s2n<G4int>(str1);
		if(index > nMaterials || index < 0) {
			G4String descript = "Index : " + str1;
			G4Exception("G4TDCMZSliceHeader::ReadDataFromFile - error:\
            bad material index","", FatalException,descript.c_str());
		}
		fMaterialNames[index] = str2;
	}

	in >> fNoVoxelX >> fNoVoxelY >> fNoVoxelZ;

	G4double tmpMinX, tmpMinY, tmpMinZ;
	G4double tmpMaxX, tmpMaxY, tmpMaxZ;

	in >> tmpMinX >> tmpMaxX;
	in >> tmpMinY >> tmpMaxY;
	in >> tmpMinZ >> tmpMaxZ;

	fMinX = (CheckConsistency(tmpMinX,fMinX,"Min X value")) ?
			fMinX : ((fMinX == 0) ? tmpMinX : fMinX);
	fMaxX = (CheckConsistency(tmpMaxX,fMaxX,"Max X value")) ?
			fMaxX : ((fMaxX == 0) ? tmpMaxX : fMaxX);

	fMinY = (CheckConsistency(tmpMinY,fMinY,"Min Y value")) ?
			fMinY : ((fMinY == 0) ? tmpMinY : fMinY);
	fMaxY = (CheckConsistency(tmpMaxY,fMaxY,"Max Y value")) ?
			fMaxY : ((fMaxY == 0) ? tmpMaxY : fMaxY);

	fMinZ = (CheckConsistency(tmpMinZ,fMinZ,"Min Z value")) ?
			fMinZ : ((fMinZ == 0) ? tmpMinZ : fMinZ);
	fMaxZ = (CheckConsistency(tmpMaxZ,fMaxZ,"Max Z value")) ?
			fMaxZ : ((fMaxZ == 0) ? tmpMaxZ : fMaxZ);

	fMateIDs.clear();
	fValues.clear();
	fMateIDs.resize(fNoVoxelY*fNoVoxelZ,std::vector<G4int>(fNoVoxelX,0));
	fValues.resize(fNoVoxelY*fNoVoxelZ,std::vector<G4double>(fNoVoxelX,0.));
	for(G4int k = 0; k < fNoVoxelZ; ++k) {
		for(G4int j = 0; j < fNoVoxelY; ++j) {
			for(G4int i = 0; i < fNoVoxelX; ++i) {
				G4int tmpMateID;
				in >> tmpMateID;
				G4int row = j*(k+1);
				fMateIDs[row][i] = tmpMateID;
			}
		}
	}

	for(G4int k = 0; k < fNoVoxelZ; ++k) {
		for(G4int j = 0; j < fNoVoxelY; ++j) {
			for(G4int i = 0; i < fNoVoxelX; ++i) {
				G4double tmpValue;
				in >> tmpValue;
				G4int row = j*(k+1);
				fValues[row][i] = tmpValue;
			}
		}
	}

	in.close();
}


