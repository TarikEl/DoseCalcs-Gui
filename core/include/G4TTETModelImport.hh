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
// G4TTETModelImport.hh
// \file   MRCP_GEANT4/Internal/include/G4TTETModelImport.hh
// \author Haegin Han
//

#ifndef G4TTETModelImport_h
#define G4TTETModelImport_h 1

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>

#include "G4UIExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4Tet.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Colour.hh"

// *********************************************************************
// This class is to import the phantom data from *.ele, *.node, and
// *.material files.
// -- DataRead: Construct G4Tet by reading data from *.ele and *.node
//              files
// -- MaterialRead: Construct G4Material by reading material data from
//                  *.material file
// -- ColourRead: Construct std::map that contains G4Colour data
//                according to data in colour.dat file
// -- PrintMaterialInformation: Print a table that contains organ ID,
//                              number of tetrahedrons, volume, density,
//                              mass, and name for each organ
// *********************************************************************

class G4TTETModelImport
{
public:
    G4TTETModelImport();
    virtual ~G4TTETModelImport() {};

	// get methods
	G4String      GetPhantomName()           { return phantomName; };
	G4Material*   GetMaterial(G4int idx)     { return materialMap[idx];}
    G4int         GetNumTetrahedron()        { return tetVector.size();}
    G4int         GetTotalNumTetrahedronInFile()        { return numEle;}
    G4int         GetMaterialIndex(G4int idx){ return materialVector[idx]; }
	G4Tet*        GetTetrahedron(G4int idx)  { return tetVector[idx]; }
	G4double      GetVolume(G4int idx)       { return volumeMap[idx]; }
	std::map<G4int, G4double> GetMassMap()   { return massMap; }
    std::map<G4int, G4Colour> GetColourMap() { return colourMap; }

	G4ThreeVector GetPhantomSize()           { return phantomSize; }
	G4ThreeVector GetPhantomBoxMin()         { return boundingBox_Min; }
	G4ThreeVector GetPhantomBoxMax()         { return boundingBox_Max; }

    void setIdNameMap(std::map<G4int, G4String> IdName) { MatIDNameMap = IdName; }

    void setNameMatMap(std::map<G4String, G4Material*> NameMat) { NameMatMap = NameMat; }

    void GenerateDataforTET();
    void PrintMaterialInfomation();

    std::vector<G4String> GetOrganNamesVector() const { return OrganNamesVector;}
    std::map<G4String, G4double> GetOrganNameDensityMap() const {return OrganNameDensityMap;}
    std::map<G4String, G4double> GetOrganNameMassMap() const {return OrganNameMassMap;}
    std::map<G4String, G4double> GetOrganNameVolumeMap() const {return OrganNameVolumeMap;}



    std::vector<G4String> DcmRegionsNames;
    //std::vector<G4bool> UseVoxelMatForSegMap;
    //std::vector<G4int> VoxelMatForSegMap;
    std::vector<G4double> DcmRegionsMinDensityMap;
    std::vector<G4double> DcmRegionsMaxDensityMap;
    std::vector<G4bool> UseDcmRegionsMinDensityMap;
    std::vector<G4bool> UseDcmRegionsMaxDensityMap;


    void setRegionSegmentationDataVectors(
            std::vector<G4String> vec1,
            //std::vector<G4bool> vec2,
            //std::vector<G4int> vec3,
            std::vector<G4bool> vec4,
            std::vector<G4bool> vec5,
            std::vector<G4double> vec6,
            std::vector<G4double> vec7) {
        DcmRegionsNames = vec1;
        //UseVoxelMatForSegMap = vec2;
        //VoxelMatForSegMap = vec3;
        UseDcmRegionsMinDensityMap = vec4;
        UseDcmRegionsMaxDensityMap = vec5;
        DcmRegionsMinDensityMap = vec6;
        DcmRegionsMaxDensityMap = vec7;
    }

private:

	// methods
	void DataRead(G4String, G4String);
	void MaterialRead(G4String);
	void ColourRead();

	G4String phantomDataPath;
	G4String phantomName;

	G4ThreeVector boundingBox_Min;
	G4ThreeVector boundingBox_Max;
	G4ThreeVector phantomSize;

	std::vector<G4ThreeVector> vertexVector;
	std::vector<G4Tet*>        tetVector;
	std::vector<G4int*>        eleVector;
	std::vector<G4int>         materialVector;
	std::map<G4int, G4int>     numTetMap;
	std::map<G4int, G4double>  volumeMap;
	std::map<G4int, G4double>  massMap;
    std::map<G4int, G4Colour>  colourMap;

	std::map<G4int, std::vector<std::pair<G4int, G4double>>> materialIndexMap;
	std::vector<G4int>                                       materialIndex;
    std::map<G4int, G4Material*>                             materialMap;
    std::map<G4String, G4Material*>                          NameMatMap;
	std::map<G4int, G4double>                                densityMap;
    std::map<G4int, G4String>                                MatIDNameMap;

    std::vector<G4String> OrganNamesVector;
    std::map<G4String, G4double> OrganNameDensityMap;
    std::map<G4String, G4double> OrganNameMassMap;
    std::map<G4String, G4double> OrganNameVolumeMap;

    G4int numEle;

};

#endif
