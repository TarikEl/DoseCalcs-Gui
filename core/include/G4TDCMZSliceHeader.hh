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
// $Id: G4TDCMZSliceHeader.hh 92820 2015-09-17 15:22:14Z gcosmo $
//
/// \file medical/DICOM/include/G4TDCMZSliceHeader.hh
/// \brief Definition of the G4TDCMZSliceHeader class
//

#ifndef G4TDCMZSliceHeader_h
#define G4TDCMZSliceHeader_h 1

#include "globals.hh"
class G4material;
#include <fstream>
#include <vector>

//*******************************************************
/// G4TDCMZSliceHeader class
///
/// Contains the meta data information corresponding to one or several
/// Z slices (number of voxels, dimension)
///
/// History: 30.11.07  First version
/// \author  P. Arce
//*******************************************************

class G4TDCMZSliceHeader
{
public:
	G4TDCMZSliceHeader(const G4String&);

	G4TDCMZSliceHeader( const G4TDCMZSliceHeader& rhs );
	// build object copying an existing one (except Z dimensions)

	G4TDCMZSliceHeader( std::ifstream& fin );
	// build object reading data from a file

	~G4TDCMZSliceHeader();

	// Get and set methods
	G4int GetNoVoxelX() const { return fNoVoxelX; };
	G4int GetNoVoxelY() const { return fNoVoxelY; };
	G4int GetNoVoxelZ() const { return fNoVoxelZ; };
	G4int GetNoVoxels() const { return fNoVoxelX*fNoVoxelY*fNoVoxelZ; };

	G4double GetMinX() const { return fMinX; };
	G4double GetMinY() const { return fMinY; };
	G4double GetMinZ() const { return fMinZ; };
	G4double GetMaxX() const { return fMaxX; };
	G4double GetMaxY() const { return fMaxY; };
	G4double GetMaxZ() const { return fMaxZ; };

	G4double GetVoxelHalfX() const { return (fMaxX-fMinX)/fNoVoxelX/2.; };
	G4double GetVoxelHalfY() const { return (fMaxY-fMinY)/fNoVoxelY/2.; };
	G4double GetVoxelHalfZ() const { return (fMaxZ-fMinZ)/fNoVoxelZ/2.; };

	const std::vector<G4String>& GetMaterialNames() const { return fMaterialNames; };

	void SetNoVoxelX(const G4int& val) { fNoVoxelX = val; }
	void SetNoVoxelY(const G4int& val) { fNoVoxelY = val; }
	void SetNoVoxelZ(const G4int& val) { fNoVoxelZ = val; }

	void SetMinX(const G4double& val) { fMinX = val; };
	void SetMaxX(const G4double& val) { fMaxX = val; };
	void SetMinY(const G4double& val) { fMinY = val; };
	void SetMaxY(const G4double& val) { fMaxY = val; };
	void SetMinZ(const G4double& val) { fMinZ = val; };
	void SetMaxZ(const G4double& val) { fMaxZ = val; };

	void SetMaterialNames(std::vector<G4String>& mn ){ fMaterialNames = mn; }

	void operator+=( const G4TDCMZSliceHeader& rhs );
	G4TDCMZSliceHeader operator+( const G4TDCMZSliceHeader& rhs );
	// add two slices that have the same dimensions, merging them in Z

	//=================================================================
	//  NEW REVISION ( Jonathan Madsen - jonathan.madsen@cern.ch )
	//  December 2012
	//
	//  New data handling format -> movement away from file-based to class based
	//      -- If density and mate ID data is not present -> read from file
	//
	//  REASONING:
	//      --  DICOM data can contain inconsistencies, handling via class
	//    instead of via file
	//          allows safe/easy modification
	//
	//  Adding Data to densities and fMateIDs
	//
	void SetFilename(const G4String& val) { fFilename = val; }
	void SetSliceLocation(const G4double& val) { fSliceLocation = val; }
	void AddMaterial(const G4String& val) { fMaterialNames.push_back(val); }

	const G4double& GetSliceLocation() const { return fSliceLocation; }

	// called from DicomHandler::StoreData
	void AddRow() { // when it's called it create a new vector in fValues, fMateIDs (std::vector< vector(G4double)>) with 0 item in this new vector
		fValues.push_back(std::vector<G4double>(0)); // it create an empty vector in vector< vector(G4double)>
		fMateIDs.push_back(std::vector<G4int>(0)); } // it create an empty vector in vector< vector(G4int)>

	// called from DicomHandler::StoreData to store the data in Zheaderslice object
	void AddValue(G4double val) {
		(fValues.size() > 0) ?
				fValues.back().push_back(val) : // it's plein then we add to the next vector item in the last vector that is active
				fValues.push_back(std::vector<G4double>(1,val)); // it's empty then we add to 1 first vector item
	}
	void AddValue(const std::vector<G4double>& val) { fValues.push_back(val); }
	void AddValue(const std::vector<std::vector<G4double> >& val) {
		for(unsigned int i = 0; i < val.size(); ++i) { fValues.push_back(val.at(i)); }
	}

	// called from DicomHandler::StoreData to store the data in Zheaderslice object
	void AddMateID(G4int val) { (fMateIDs.size() > 0) ?
			fMateIDs.back().push_back(val) :
			fMateIDs.push_back(std::vector<G4int>(1,val));
	}
	void AddMateID(const std::vector<G4int>& val) { fMateIDs.push_back(val); }
	void AddMateID(const std::vector<std::vector<G4int> >& val) {
		for(unsigned int i = 0; i < val.size(); ++i) { fMateIDs.push_back(val.at(i)); }
	}


	const std::vector<std::vector<G4double> >& GetValues() const { return fValues; }
	const std::vector<std::vector<G4int> >& GetMateIDs() const { return fMateIDs; }

	void DumpToFile();

	void ReadDataFromFile();

	void DumpExcessMemory() {
		if(fFilename.length() != 0) { fValues.clear(); fMateIDs.clear(); } }

	void FlipData();

private:

	// all this method are used in G4TDCMZSliceHeader class pour lire bien les données en .g4dcom fichiers et verifier son dype

	inline G4bool IsInteger(const G4String&);

	template <typename T>
	inline void Print(std::ostream&, const std::vector<T>&, const G4String&,
			G4int breakLine = -1);
	template <typename T> inline T G4s2n(const G4String&);
	template <typename T> inline bool CheckConsistency(const T&, const T&, G4String);
	//
	//  END NEW REVISION
	//=======================================================================

private:
	G4bool CheckMaterialExists( const G4String& mateName );
	// check that material read exists as a G4Material

private:
	G4int fNoVoxelX, fNoVoxelY, fNoVoxelZ;  // number of voxels in each dimensions
	G4double fMinX,fMinY,fMinZ; // minimum extension of voxels (position of wall)
	G4double fMaxX,fMaxY,fMaxZ; // maximum extension of voxels (position of wall)

	std::vector<G4String> fMaterialNames; // list of material names

	G4String fFilename;
	std::vector<std::vector<G4double> > fValues;
	std::vector<std::vector<G4int> > fMateIDs;
	G4double fSliceLocation;

};

//============================================================================
// This function flips all the data, Otherwise, the image is upside-down
inline void G4TDCMZSliceHeader::FlipData()
{
	std::reverse(fValues.begin(), fValues.end());
	std::reverse(fMateIDs.begin(), fMateIDs.end());
}
//=============================================================================
inline G4bool G4TDCMZSliceHeader::IsInteger(const G4String& str)
{
	return (str.find_first_not_of("0123456789") == std::string::npos) ? true : false;
}
//============================================================================
template <typename T>
inline T G4TDCMZSliceHeader::G4s2n(const G4String& str)
{
	std::istringstream iss(str);
	T val;
	iss >> val;
	return val;
}

//============================================================================
template <typename T>
inline bool G4TDCMZSliceHeader::CheckConsistency(const T& val1, const T& val2,
		G4String category) {
	if(val1 != val2) {
		G4Exception("G4TDCMSliceZHeader::CheckConsistency",
				"Consistency Mismatch : Keeping previous value if nonzero",
				JustWarning, category.c_str());
		return false;
	}
	return true;
}
//============================================================================
template <typename T>
inline void G4TDCMZSliceHeader::Print(std::ostream& out, const std::vector<T>& val,
		const G4String& delim, G4int breakLine)
{
	for(unsigned int i = 0; i < val.size(); ++i) {
		out << val.at(i);
		if(breakLine < 0) {
			if(i+1 < val.size()) { out << delim; }
			else { out << G4endl; }
		} else {
			((i != 0 && i%breakLine == 0) ? (out << G4endl) : (out << delim)); }
	}
}
//==========================================================================

#endif
