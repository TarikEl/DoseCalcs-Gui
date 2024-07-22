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
//
/// \file G4TFlashRadioTherapyGeometry.hh

#ifndef G4TFlashRadioTherapyGeometry_h
#define G4TFlashRadioTherapyGeometry_h 1

#include "globals.hh"

#include "G4Material.hh"
#include "tls.hh"
#include "G4ThreeVector.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;
class G4NistManager;
class G4Tubs;
class G4Box;
class G4Element;
class G4VisAttributes;



class G4TFlashRadioTherapyGeometry {
public:


    G4VPhysicalVolume *physicalTreatmentRoom;
    G4LogicalVolume *logicTreatmentRoom;
    G4VPhysicalVolume *ConstructPhantom(G4double CollPos);
    G4VPhysicalVolume *ConstructDetector();
    G4VPhysicalVolume *ConstructCollimator();

    G4VPhysicalVolume *ConstructGeometry();

    G4TFlashRadioTherapyGeometry();
    virtual ~G4TFlashRadioTherapyGeometry();

    G4bool  SetPhantomMaterial(G4String material);
    G4bool  SetDetectorMaterial(G4String material);
    void SetAirGap(G4double position);
    void SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ);

    void SetDetectorThickness(G4double thickness);
    void SetDetector_subThickness(G4double thickness_sub);
    void SetDetectorWidth(G4double width);
    void SetDetectorPosition(G4double position);
    void SetAirGap_water_detector(G4double spost);

    G4VisAttributes *skyBlue;
    G4VisAttributes *red;
    G4VisAttributes *blue;
    G4VisAttributes *green;


private:

    G4Material *airNist;
    G4Material *fPhantomMaterial;

    G4double fAirGap;
    G4double fPhantomSizeX, fPhantomSizeY, fPhantomSizeZ, fPhantom_coordinateX,fPosition_coefficient;
    G4ThreeVector fPhantomPosition;
    G4double fDet_thickness,fDet_width,fDet_sub_thickness,fDetectorPosition,fAirGap_phantom_det;
    G4Element *Si;
    G4Element *C;
    G4Material *SiC;
    G4Material *fDetectorMaterial;

    G4Box *fPhantom;


    G4Box *fDet_box;
    G4LogicalVolume *fDetLogicalVolume;
    G4VPhysicalVolume *fDet_phys;

    G4Box *fDet_sub;
    G4LogicalVolume *fDet_sub_LogicalVolume;
    G4VPhysicalVolume *fDet_sub_phys;


    void DefineMaterials();


    G4LogicalVolume *fPhantomLogicalVolume;
    G4VPhysicalVolume *fPhant_phys;
    G4VPhysicalVolume *fPhantom_physical;
    
    G4UserLimits *fStepLimit;
    G4bool fCheckOverlaps;

    G4NistManager *nist;




    // ///////////////////////////////////  Collimator Constructing data

    void SetOuterRadius(G4double radius);
    void SetApplicatorLength(G4double length);

private:

    void FlashBeamLineVacuumSource();

    void FlashBeamLineTitaniumWindows();
    void FlashVWAlcover();
    void FlashAlCover2();
    void FlashExitBit();
    void FlashToroid();
    void OverCover();
    void OverCover2();
    void MonitorChamber();
    void Flash_connector();
    void Bigconnector();
    void Bigconnector2();
    void Bigconnector3();
    void FlashBeamLineApplicator();


    G4double fInitial_pos;

    G4double fInnerRadiusFirstApplicatorFlash;

    G4Material* Fe;
    G4Material* PVDF;
    G4Material* FILM;
    G4Material* aluminumNist;
    G4Material* PMMA;

    G4double fOutRadiusVSFlash;
    G4double fHightVSFlash;
    G4double fXPositionVSFlash;
    G4double fHightFTFlash;

    G4double fOutRadiusFTFlash;
    G4double fOutRadius;
    G4double fToroid_outRadius;
    G4double fToroid_hight;
    G4double fToroid_XPosition;
    G4double fBigcover_hight;
    G4double fBigcover_XPosition;
    G4double fChamberpos;


    void SetDefaultDimensions();

    void ConstructApplicator();

    G4VisAttributes *gray;
    G4VisAttributes *white;
    G4VisAttributes *yellow;
    G4VisAttributes *darkGreen;
    G4VisAttributes *darkOrange3;
    G4VisAttributes *magenta;

    G4double fOuterRadiusFirstApplicatorFlash;
    G4Tubs *fSolidFirstApplicatorFlash;
    G4VPhysicalVolume *fPhysiFirstApplicatorFlash;
    G4Material *fFirstApplicatorMaterialFlash;


    //  Titanium Window
    G4Tubs *solidFTFlash;
    G4VPhysicalVolume *physiFTFlash;
    G4Material *FTFlashMaterialFlash;

    //  Vacuum Source
    G4Tubs *solidVSFlash;
    G4VPhysicalVolume *physiVSFlash;
    G4Material *VSFlashMaterialFlash;

    G4double fFinalApplicatorXPositionFlash;
    G4double fHightFinalApplicatorFlash;

};

#endif
