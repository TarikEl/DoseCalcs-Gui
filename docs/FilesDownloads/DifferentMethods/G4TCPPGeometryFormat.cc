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
//


#include "G4TCPPGeometryFormat.hh"
#include "G4VPhysicalVolume.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"

#include "G4Material.hh"
#include "G4LogicalVolume.hh"

#include "G4Ellipsoid.hh"

G4TCPPGeometryFormat::G4TCPPGeometryFormat(){}
G4TCPPGeometryFormat::~G4TCPPGeometryFormat(){}

G4VPhysicalVolume* G4TCPPGeometryFormat::ConstructPhysicalVolume(){
    return WorldPhysicalVolume;
}

void G4TCPPGeometryFormat::ConstructLogicalVolumes(){

    G4Material* matH2O;
    G4Material* soft;

    G4double A;  // atomic mass
    G4double Z;  // atomic number
    G4double d;  // density

    // General elements

    A = 1.01*g/mole;
    G4Element* elH = new G4Element ("Hydrogen","H",Z = 1.,A);

    A = 12.011*g/mole;
    G4Element* elC = new G4Element("Carbon","C",Z = 6.,A);

    A = 14.01*g/mole;
    G4Element* elN = new G4Element("Nitrogen","N",Z = 7.,A);

    A = 16.00*g/mole;
    G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A);

    A = 22.99*g/mole;
    G4Element* elNa = new G4Element("Sodium","Na",Z = 11.,A);

    A = 28.086*g/mole;
    G4Element* elSi = new G4Element("Si","Si",Z = 14.,A);

    A = 24.305*g/mole;
    G4Element* elMg = new G4Element("Magnesium","Mg",Z = 12.,A);

    A = 30.974*g/mole;
    G4Element* elP = new G4Element("Phosphorus","P",Z = 15.,A);

    A = 32.064*g/mole;
    G4Element* elS = new G4Element("Sulfur","S",Z = 16.,A);

    A = 35.453*g/mole;
    G4Element* elCl = new G4Element("Chlorine","Cl",Z = 17.,A);

    A = 39.098*g/mole;
    G4Element* elK = new G4Element("Potassium","K",Z = 19.,A);

    A = 40.08*g/mole;
    G4Element* elCa = new G4Element("Calcium","Ca",Z = 20.,A);

    A = 55.85*g/mole;
    G4Element* elFe  = new G4Element("Iron","Fe",Z = 26.,A);

    A = 65.38*g/mole;
    G4Element* elZn = new G4Element("Zinc","Zn",Z = 30.,A);

    A = 85.47 *g/mole;
    G4Element* elRb = new G4Element("Rb","Rb",Z = 37.,A);

    A = 87.62 *g/mole;
    G4Element* elSr = new G4Element("Sr","Sr",Z = 38.,A);

    A = 91.22 *g/mole;
    G4Element* elZr = new G4Element("Zr","Zr",Z = 40.,A);

    A = 207.19 *g/mole;
    G4Element* elPb = new G4Element("Lead","Pb", Z = 82.,A);

    // ORNL soft tissue
    d = 1.04 *g/cm3;
    soft = new G4Material("soft_tissue",d,18);
    soft->AddElement(elH,0.10454);
    soft->AddElement(elC,0.22663);
    soft->AddElement(elN,0.02490);
    soft->AddElement(elO,0.63525);
    soft->AddElement(elNa,0.00112);
    soft->AddElement(elMg,0.00013);
    soft->AddElement(elSi,0.00030);
    soft->AddElement(elP,0.00134);
    soft->AddElement(elS,0.00204);
    soft->AddElement(elCl,0.00133);
    soft->AddElement(elK,0.00208);
    soft->AddElement(elCa,0.00024);
    soft->AddElement(elFe,0.00005);
    soft->AddElement(elZn,0.00003);
    soft->AddElement(elRb,0.00001);
    soft->AddElement(elSr,0.0);
    soft->AddElement(elZr,0.00001);
    soft->AddElement(elPb,0.0);

    G4double ax = 6.84 * cm;
    G4double by= 8.7 * cm;
    G4double cz = 5.9 * cm;

    G4Ellipsoid* brain = new G4Ellipsoid("Brain", ax, by, cz);

    new G4LogicalVolume(brain, soft, "Brain", 0, 0, 0);
}
