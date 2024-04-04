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
#include "G4TNuclearReactorGeometry.hh"

////#include "G4Material.hh"
////#include "G4Box.hh"
////#include "G4Tubs.hh"
////#include "G4LogicalVolume.hh"
////#include "G4PVPlacement.hh"
////#include "G4VisAttributes.hh"
////#include "G4SystemOfUnits.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomTools.hh"
#include "G4Polyhedra.hh"

G4TNuclearReactorGeometry::G4TNuclearReactorGeometry()
{
}

G4TNuclearReactorGeometry::~G4TNuclearReactorGeometry(){}

G4VPhysicalVolume* G4TNuclearReactorGeometry::ConstructNuclearReactor()
{

    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* waterMaterial = nistManager->FindOrBuildMaterial("G4_WATER");
    G4Material* structuralMaterial = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");

    //G4Material* uraniumOxide = new G4Material("UraniumOxide", 92., 238.02891*g/mole, 8.);
    G4Material* water = new G4Material("Water", 1., 1.01*g/mole, 1.0*g/cm3);
    G4Material* graphite = new G4Material("Graphite", 6., 12.01*g/mole, 2.267*g/cm3);
    //G4Material* graphite = nistManager->FindOrBuildMaterial("G4_GRAPHITE");
    G4Material* concrete = new G4Material("Concrete", 14., 28.09*g/mole, 2.3*g/cm3);
    G4Material* beryllium = new G4Material("Beryllium", 4., 9.0122*g/mole, 1.85*g/cm3);
    G4Material* stainlessSteel = new G4Material("StainlessSteel", 13., 26.98*g/mole, 7.85*g/cm3);
    G4Material* zirconium = new G4Material("Zirconium", 40., 91.224*g/mole, 6.506*g/cm3);

    G4Isotope* B10 = new G4Isotope("B10", 5, 10, 10.0129370 * g/mole);
    G4Element* elBoron = new G4Element("Boron", "B", 1);
    elBoron->AddIsotope(B10, 1.0);
    G4Material* absorberMaterial = new G4Material("Bore-10", 2.34 * g/cm3, 1);
    absorberMaterial->AddElement(elBoron, 1);

    // Define oxygen
    G4Element* oxygen = nistManager->FindOrBuildElement("O");
    G4double enrichmentFraction = 0.1; // 90% enrichment
    // Define uranium-238 isotope
    G4Element* uranium = new G4Element("Uranium", "U", 2);
    G4Isotope* U235 = new G4Isotope("U235", 92, 235, 235.04393 * g/mole);
    G4Isotope* U238 = new G4Isotope("U238", 92, 238, 238.05078 * g/mole);
    uranium->AddIsotope(U235, enrichmentFraction);
    uranium->AddIsotope(U238, 1-enrichmentFraction);
    // Define UO2 material with 90% U-235 enrichment
    G4double density_UO2 = 10.97 * g/cm3; // Density of UO2
    G4Material* uraniumOxide = new G4Material("UO2", density_UO2, 2);
    uraniumOxide->AddElement(oxygen, 2);
    uraniumOxide->AddElement(uranium, 1);


    // World volume
    G4double worldSize = 2000.0 * mm; // Adjust as needed
    G4Box* worldSolid = new G4Box("World", worldSize, worldSize, worldSize);
    G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid, waterMaterial, "World", 0, 0,0);
    G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0, G4ThreeVector(),"World", worldLogical, 0 , false, 0 , false );

    //G4int numAbsorberBoreRods = 5;

    G4double moderatorRadius = 600 * mm; // Adjust as needed
    G4double moderatorLength = 300 * mm; // Adjust as needed
    // Create coolant (Water)
    G4Tubs* moderatorSolid = new G4Tubs("Moderator", 0.0, moderatorRadius, moderatorLength, 0.0 ,360.*deg);
    G4LogicalVolume* moderatorLogical = new G4LogicalVolume(moderatorSolid, water, "Moderator");
    new G4PVPlacement(0, G4ThreeVector(), moderatorLogical, "Moderator", worldLogical, false, 0);
    //G4Box* moderatorSolid = new G4Box("Moderator", 1000,1000,1000);
    //G4LogicalVolume* moderatorLogical = new G4LogicalVolume(moderatorSolid, water, "Moderator");
    //new G4PVPlacement(0, G4ThreeVector(), moderatorLogical, "Moderator", worldLogical, false, 0);

    //G4Colour colour1 = G4Colour((G4double)G4UniformRand(), (G4double)G4UniformRand(), (G4double)G4UniformRand());
    //G4Colour colour2 = G4Colour((G4double)G4UniformRand(), (G4double)G4UniformRand(), (G4double)G4UniformRand());
    //G4VisAttributes* organVisAtt1 = new G4VisAttributes(colour1);
    //G4VisAttributes* organVisAtt2 = new G4VisAttributes(colour2);

    G4double fuelRodRadius = 25 * mm; // Adjust as needed
    G4double fuelRodLength = 250 * mm; // Adjust as needed
    G4double distanceBetwwenX = 60;
    G4double distanceBetwwenY = 60;

    G4int numFuelRodsX = 9;
    G4int numFuelRodsY = 9;
    G4double FirstX = -((numFuelRodsX-1)*distanceBetwwenX)/2;
    G4double FirstY = -((numFuelRodsY-1)*distanceBetwwenY)/2;

    G4int numFuelRods = 10;
    //Create fuel rods
    G4int inc = 0;
    G4double posx = FirstX;
    for (G4int i = 1; i < numFuelRodsX+1; ++i) {
        G4double posy = FirstY;

        for (G4int j = 1; j < numFuelRodsY+1; ++j) {

            std::ostringstream name; name << "FuelRodX"<<i << "_Y"<< j;

            G4String namea = "MyPolyhedra";
            G4int numSides = 6;  // Number of sides (e.g., hexagon)
            G4double phiStart = 0.0;  // Starting angle of the first side
            G4double phiTotal = 360.0 * degree;  // Total angle covered by the polyhedra
            G4int numZPlanes = 2;  // Number of z planes (z coordinates of the corners)
            G4double zPlane[2] = {-250.0, 250.0};  // Z coordinates of the corners
            G4double rInner[2] = {0.0, 0.0};  // Inner radius at each z plane
            G4double rOuter[2] = {25.0, 25.0};  // Outer radius at each z plane

            // Create G4Polyhedra
            G4Polyhedra* fuelRodSolid = new G4Polyhedra(name.str(), phiStart, phiTotal, numSides, numZPlanes, zPlane, rInner, rOuter);
            //G4Tubs* fuelRodSolid = new G4Tubs(name.str(), 0.0, fuelRodRadius, fuelRodLength, 0.0, 360.0 * deg);
            G4LogicalVolume* fuelRodLogical = new G4LogicalVolume(fuelRodSolid, uraniumOxide, name.str());
            new G4PVPlacement(0, G4ThreeVector(posx, posy, 0), fuelRodLogical, name.str(), moderatorLogical, false, inc);
            //fuelRodLogical->SetVisAttributes(organVisAtt1);

            //G4VisAttributes* voxvisa = new G4VisAttributes();
            //voxvisa->SetVisibility(true); voxvisa->SetForceSolid(false);
            //fuelRodLogical->SetVisAttributes(voxvisa);

            posy += distanceBetwwenY;
            inc++;
        }
        posx += distanceBetwwenX;
    }

    // Create absorber bore rods
    G4double absorberBoreRadius = 2.0 * cm; // Adjust as needed
    G4int numAbsorberBoreRods = 5;
    G4double distanceBetwwenAbsX = 40;
    G4double distanceBetwwenAbsY = 40;
    G4int numAbsRodsX = 0;
    G4int numAbsRodsY = 0;
    G4double FirstAbsX = -(numFuelRodsX*distanceBetwwenAbsX)/2;
    G4double FirstAbsY = -(numFuelRodsY*distanceBetwwenAbsY)/2;

    //Create fuel rods
    posx = FirstAbsX;
    inc = 0;
    for (G4int i = 1; i < numAbsRodsX+1; ++i) {
        G4double posy = FirstAbsY;

        for (G4int j = 1; j < numAbsRodsY+1; ++j) {

            std::ostringstream name; name << "AbsRodX"<<i << "_Y"<< j;

            G4String namea = "MyPolyhedra";
            G4int numSides = 6;  // Number of sides (e.g., hexagon)
            G4double phiStart = 0.0;  // Starting angle of the first side
            G4double phiTotal = 360.0 * degree;  // Total angle covered by the polyhedra
            G4int numZPlanes = 2;  // Number of z planes (z coordinates of the corners)
            G4double zPlane[2] = {-250.0, 250.0};  // Z coordinates of the corners
            G4double rInner[2] = {0.0, 0.0};  // Inner radius at each z plane
            G4double rOuter[2] = {25.0, 25.0};  // Outer radius at each z plane

            // Create G4Polyhedra
            G4Polyhedra* absorberBoreRodSolid = new G4Polyhedra(name.str(), phiStart, phiTotal, numSides, numZPlanes, zPlane, rInner, rOuter);
            //G4Tubs* absorberBoreRodSolid = new G4Tubs(name.str(), 0.0, fuelRodRadius, fuelRodLength, 0.0, 360.0 * deg);
            G4LogicalVolume* absorberBoreRodLogical = new G4LogicalVolume(absorberBoreRodSolid, absorberMaterial, name.str());
            new G4PVPlacement(0, G4ThreeVector(posx, posy, 0), absorberBoreRodLogical, name.str(), moderatorLogical, false, inc);
            //absorberBoreRodLogical->SetVisAttributes(organVisAtt2);

            posy += distanceBetwwenY;
            inc++;
        }
        posx += distanceBetwwenX;
    }

    // Create reflector (Concrete, Beryllium, or Graphite)
    G4double reflectorRadius = 400*mm;
    G4Tubs* reflectorSolid = new G4Tubs("Reflector", moderatorRadius, moderatorRadius+reflectorRadius, moderatorLength, 0., 360.*deg);
    G4LogicalVolume* reflectorLogical = new G4LogicalVolume(reflectorSolid, concrete, "Reflector");
    G4VPhysicalVolume* reflectorPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), reflectorLogical, "Reflector", worldLogical, false, 0);

    // Create structural materials (Stainless Steel, Zirconium)
    G4double structuralRadius = 100*mm;
    G4Tubs* structuralSolid = new G4Tubs("Structural", moderatorRadius+reflectorRadius, moderatorRadius+reflectorRadius+structuralRadius, moderatorLength, 0., 360.*deg);
    G4LogicalVolume* structuralLogical = new G4LogicalVolume(structuralSolid, stainlessSteel, "Structural");
    G4VPhysicalVolume* structuralPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), structuralLogical, "Structural", worldLogical, false, 0);

    //worldLogical->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));
    //moderatorLogical->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));

    return worldPhysical;
}
