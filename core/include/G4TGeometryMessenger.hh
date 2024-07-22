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

#ifndef G4TGeometryMessenger_h
#define G4TGeometryMessenger_h 1

class G4TVolumeConstruction;


class G4UIcommand;
class G4UIparameter;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3VectorAndUnit;
class G4VSolid;


#include "G4UImessenger.hh"
#include "globals.hh"
#include <iostream>

class G4TGeometryMessenger: public G4UImessenger
{
public:

    G4TGeometryMessenger(G4TVolumeConstruction* GeometryConstruction);


    ~G4TGeometryMessenger();

    void SetNewValue(G4UIcommand* command, G4String newValue);

    void CommandsForMaterials();
    void CommandsForGeometries();
    void CommandsForVoxGeometryType();
/*
    void CommandsForDcmGeometryType();
    void CommandsForSolids();
    void CommandsForVolumes();
*/
    G4double UseG4Units(G4String);

private:

    G4VSolid* sld;

    G4int MatElementNum;
    G4TVolumeConstruction*  GeometryConstruction;

    G4UIdirectory*               ConstructGeometryCMD;
    G4UIdirectory*               MaterialConstruction;

    G4UIcommand*          NewMaterialCMD;
    G4UIcommand*          NewElementCMD;
    G4UIcommand*          MatElementCMD;
    G4UIcommand*          MatMaterialsCMD;
    G4UIcommand*          NistMaterialNameCMD;

    G4UIcommand*          CreateWorldCMD;

    G4UIcommand*          VoxelsDataCMD;

    G4UIcmdWithAString*          MaterialNameAsRegionNameCMD ;
    G4UIcommand*                      TETRegionDataCMD;

    G4UIcmdWithAString*          VoxDefaultMaterialNameCMD ;
    G4UIcmdWithAString*          GeometrySymbolCMD ;
    G4UIcmdWith3VectorAndUnit*   VoxContainerPosCMD ;
    G4UIcmdWith3VectorAndUnit*   VoxContainerRotCMD ;
    G4UIcommand*                 VoxelizedPlanesCMD ;
    G4UIcommand*                 setTETPhantomLimitsCMD ;

    G4UIcommand*          VoxelizedRegionDataCMD;

    G4UIcmdWithAnInteger* DcmPixelsCompressionCMD;
    G4UIcommand*          DICOM_CtDensityCMD;
    G4UIcommand*          DICOM_MaterialNamesOrderedCMD;
    G4UIcommand*          DICOMFileType_FilesPathCMD;
    G4UIcommand*          setVolumesVisToNotForced;

    G4UIcommand*          CreateSolidCMD ;
    G4UIcommand*          CreateVolumeCMD ;
    G4UIcommand*          DumpGeomCMD ;

    G4UIcommand*          SourceVolumeData ;


/*

    G4UIcmdWith3VectorAndUnit*   WorldHalfSizeCMD ;
    G4UIcmdWithAString*          WorldMaterialNameCMD ;
    G4UIcmdWithoutParameter*     ConstructWorldCmd;
    G4UIcmdWithoutParameter*     WorldAsSDCMD;

    G4UIcmdWithAString*          GeometryFileTypeCMD;
    G4UIcmdWithAString*          GeometryFileOrDirCMD;
    G4UIcmdWithAString*          RegionsMassDataPathCMD;
    G4UIcmdWithAString*          GeometryPathCMD;

    G4UIcmdWithAString*          DcmTypeCMD;
    G4UIcmdWithAString*          DcmFilesDirCMD;
    G4UIcmdWithAString*          DcmMaterialNameCMD;
    G4UIcmdWithADouble*          DcmCtNumberCMD;
    G4UIcmdWithADoubleAndUnit*   DcmCtDensityCMD;
    G4UIcmdWithADoubleAndUnit*   DcmRegionMinDensityCMD;
    G4UIcmdWithADoubleAndUnit*   DcmRegionMaxDensityCMD;

    G4UIcmdWithAString*          SolidNameCMD ;
    G4UIcmdWithAString*          SolidTypeCMD ;
    G4UIcmdWithAString*          VolumeSolidCMD ;
    G4UIcmdWithAString*          VolumeMaterialCMD ;
    G4UIcmdWithAString*          VolumeNameCMD ;
    G4UIcmdWith3VectorAndUnit*   BoxDimCMD ;
    G4UIcmdWithAString*          MotherVolumeNameCMD ;
    G4UIcmdWith3VectorAndUnit*   PositionCMD ;
    G4UIcmdWith3VectorAndUnit*   RotationCMD ;
    G4UIcmdWithAString*          PlaceVolumeCMD ;

    G4UIcmdWithADouble*          ElementZCMD;
    G4UIcmdWithADouble*          ElementACMD;
    G4UIcmdWithAString*          ElementNameCMD;
    G4UIcmdWithAString*          MaterialNameCMD;
    G4UIcmdWithAnInteger*        MaterialIDCMD;
    G4UIcmdWithADoubleAndUnit*   MaterialDensityCMD;
    G4UIcmdWithAnInteger*        MaterialCompNumberCMD;
    G4UIcmdWithAString*          CompNameCMD;
    G4UIcmdWithAnInteger*        CompNumberCMD;
    G4UIcmdWithADouble*          CompFractionCMD;

    G4UIcommand*          CreateBoxCMD ;
    G4UIcommand*          CreateTubeCMD ;
    G4UIcommand*          CreateCutTubsCMD;
    G4UIcommand*          CreateConsCMD ;
    G4UIcommand*          CreateParaCMD ;
    G4UIcommand*          CreateTrdCMD ;
    G4UIcommand*          CreateSphereCMD ;
    G4UIcommand*          CreateOrbCMD ;
    G4UIcommand*          CreateTorusCMD ;
    G4UIcommand*          CreateEllipsoidCMD ;
    G4UIcommand*          CreateUnionCMD ;
    G4UIcommand*          CreateIntersectionCMD ;
    G4UIcommand*          CreateSubtractionCMD ;

    G4UIcmdWithADoubleAndUnit*   RmaxCMD;
    G4UIcmdWithADoubleAndUnit*   RminCMD;
    G4UIcmdWithADoubleAndUnit*   Rmin1CMD;
    G4UIcmdWithADoubleAndUnit*   Rmax1CMD;
    G4UIcmdWithADoubleAndUnit*   Rmin2CMD;
    G4UIcmdWithADoubleAndUnit*   Rmax2CMD;
    G4UIcmdWithADoubleAndUnit*   Dx1CMD;
    G4UIcmdWithADoubleAndUnit*   Dx2CMD;
    G4UIcmdWithADoubleAndUnit*   Dy1CMD;
    G4UIcmdWithADoubleAndUnit*   Dy2CMD;
    G4UIcmdWithADoubleAndUnit*   DxCMD;
    G4UIcmdWithADoubleAndUnit*   DyCMD;
    G4UIcmdWithADoubleAndUnit*   DzCMD;
    G4UIcmdWithADoubleAndUnit*   RtorCMD;
    G4UIcmdWithADoubleAndUnit*   SPhiCMD;
    G4UIcmdWithADoubleAndUnit*   DPhiCMD;
    G4UIcmdWithADoubleAndUnit*   SThetaCMD;
    G4UIcmdWithADoubleAndUnit*   DThetaCMD;
    G4UIcmdWithADoubleAndUnit*   AlphaCMD;
    G4UIcmdWithADoubleAndUnit*   Theta0CMD;
    G4UIcmdWithADoubleAndUnit*   Phi0CMD;
    G4UIcmdWithADoubleAndUnit*   xSemiAxisCMD;
    G4UIcmdWithADoubleAndUnit*   ySemiAxisCMD;
    G4UIcmdWithADoubleAndUnit*   zSemiAxisCMD;

    G4UIcmdWith3VectorAndUnit*   LowNormCMD ;
    G4UIcmdWith3VectorAndUnit*   HighNormCMD ;
    G4UIcmdWith3VectorAndUnit*   SolidTransCMD ;
    G4UIcmdWith3VectorAndUnit*   SolidRotCMD ;

    G4UIcmdWithAString*          FirstSolidCMD ;
    G4UIcmdWithAString*          SecondSolidCMD ;

    G4UIcmdWithAnInteger*        VoxXNumberCMD;
    G4UIcmdWithAnInteger*        VoxYNumberCMD;
    G4UIcmdWithAnInteger*        VoxZNumberCMD;
    G4UIcmdWithADoubleAndUnit*   VoxXHalfSizeCMD;
    G4UIcmdWithADoubleAndUnit*   VoxYHalfSizeCMD;
    G4UIcmdWithADoubleAndUnit*   VoxZHalfSizeCMD;
    G4UIcmdWithAString*          VoxRegionNameCMD ;
    G4UIcmdWithAnInteger*        VoxelRegionMaterialIDCMD ;
    G4UIcmdWithAnInteger*        VoxRegionMinXCMD;
    G4UIcmdWithAnInteger*        VoxRegionMaxXCMD;
    G4UIcmdWithAnInteger*        VoxRegionMinYCMD;
    G4UIcmdWithAnInteger*        VoxRegionMaxYCMD;
    G4UIcmdWithAnInteger*        VoxRegionMinZCMD;
    G4UIcmdWithAnInteger*        VoxRegionMaxZCMD;
    G4UIcmdWithAnInteger*        VoxelsegmentedMaterialCMD;

*/

};

#endif
