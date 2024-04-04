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
#include "G4TDetectorGeometry.hh"


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
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

G4TDetectorGeometry::G4TDetectorGeometry()
{
}
G4TDetectorGeometry::~G4TDetectorGeometry(){
}

G4VPhysicalVolume* G4TDetectorGeometry::ConstructDetector()
{

    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* waterMaterial = nistManager->FindOrBuildMaterial("G4_WATER");
    G4Material* structuralMaterial = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");

    //G4Material* uraniumOxide = new G4Material("UraniumOxide", 92., 238.02891*g/mole, 8.);
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


    // NaITl
    G4Element* elNa = new G4Element("Sodium", "Na", 11, 22.989769 * g/mole);
    G4Element* elI = new G4Element("Iodine", "I", 53, 126.90447 * g/mole);
    G4Element* elTl = new G4Element("Thallium", "Tl", 81, 204.383 * g/mole);
    G4Material* NaI = new G4Material("NaI", 3.67 * g/cm3, 2);
    NaI->AddElement(elNa, 1);
    NaI->AddElement(elI, 1);
    G4Material* NaITl = new G4Material("NaI_Tl", 3.67 * g/cm3, 2);
    NaITl->AddMaterial(NaI, 0.9999);
    NaITl->AddElement(elTl, 0.0001); // Thallium doping level (0.01%)

    // SiO2
    G4Element* elSi = new G4Element("Silicon", "Si", 14, 28.0855 * g/mole);
    G4Element* elO = new G4Element("Oxygen", "O", 8, 16.00 * g/mole);
    G4Material* SiO2 = new G4Material("SiO2", 2.65 * g/cm3, 2);
    SiO2->AddElement(elSi, 1);
    SiO2->AddElement(elO, 2);
    G4double HSSiO2 = 39.95*mm;

    // World volume
    G4Material* airMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4double worldSize = 1000.0 * mm; // Adjust as needed
    G4Box* worldSolid = new G4Box("World", worldSize, worldSize, worldSize);
    G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid, airMaterial, "World", 0, 0,0);
    G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0, G4ThreeVector(),"World", worldLogical, 0 , false, 0 , false );


    // /////////////////////////////////////////////////////// Germanium //////////////////////////////////////////////////////////////

    //G4Material* germaniumMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ge");
    //G4Box* HPGeSol = new G4Box("HPGe", 2.525*cm, 1*cm, 3*cm);
    //G4LogicalVolume* HPGeLogical = new G4LogicalVolume(HPGeSol, germaniumMaterial, "HPGe", 0, 0,0);
    //G4VPhysicalVolume* HPGePhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), HPGeLogical, "MgO", worldLogical, false, 0);

    //G4Isotope* Li6 = new G4Isotope("Li6", 3, 6, 6.015*g/mole); G4Isotope* Li7 = new G4Isotope("Li7", 3, 7, 7.016*g/mole);
    //G4Element* lithium = new G4Element("Lithium", "Li", 2); lithium->AddIsotope(Li6, 95.72*perCent); lithium->AddIsotope(Li7, 4.28*perCent);
    //G4Material* lithiumMaterial = new G4Material("LithiumMaterial", 0.534*g/cm3, 1); lithiumMaterial->AddElement(lithium, 1);
    //G4Box* Box = new G4Box("Box", 2.545*cm, 1.02*cm, 3*cm);
    //G4SubtractionSolid* LiSolid = new G4SubtractionSolid( "Li" , Box, HPGeSol, 0, G4ThreeVector(0,0,0));
    //G4LogicalVolume* LiLogical = new G4LogicalVolume(LiSolid, lithiumMaterial, "Li", 0, 0,0);
    //G4VPhysicalVolume* LiPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), LiLogical, "Li", worldLogical, false, 0);

    //G4Box* Box1 = new G4Box("Box1", 3.6*cm, 6.22*cm, 3*cm);
    //G4SubtractionSolid* AirSolid = new G4SubtractionSolid( "Li" , Box1, Box, 0, G4ThreeVector(0,0,0));
    //G4LogicalVolume* AirLogical = new G4LogicalVolume(AirSolid, airMaterial, "Air", 0, 0,0);
    //G4VPhysicalVolume* AirPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), AirLogical, "Air", worldLogical, false, 0);

    //G4Material* aluminumMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    //G4Box* Box2 = new G4Box("Box2", 3.9*cm, 6.32*cm, 3*cm);
    //G4SubtractionSolid* AlSolid = new G4SubtractionSolid( "Al" , Box2, Box1, 0, G4ThreeVector(0,0,0));
    //G4LogicalVolume* AlLogical = new G4LogicalVolume(AlSolid, aluminumMaterial, "Al", 0, 0,0);
    //G4VPhysicalVolume* AlPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), AlLogical, "Al", worldLogical, false, 0);




    //// The Water Tank
    //G4Box* waterTank_box = new G4Box("Tank", 50*cm, 50*cm, 50*cm);
    //G4LogicalVolume* waterTank_log = new G4LogicalVolume(waterTank_box, waterMaterial, "Tank", 0, 0, 0);
    //G4VPhysicalVolume* waterTank_phys = new G4PVPlacement( 0, G4ThreeVector(), waterTank_log, "Tank", worldLogical, false, 0);

    //return worldPhysical;

    // ///////////////////////////////////////////////////////////// Optical Properties /////////////////////////////////////////////////////////////


    // ------------ Generate & Add Material Properties Table ------------
    //
    std::vector<G4double> photonEnergy = {
        2.034 * eV, 2.068 * eV, 2.103 * eV, 2.139 * eV, 2.177 * eV, 2.216 * eV,
        2.256 * eV, 2.298 * eV, 2.341 * eV, 2.386 * eV, 2.433 * eV, 2.481 * eV,
        2.532 * eV, 2.585 * eV, 2.640 * eV, 2.697 * eV, 2.757 * eV, 2.820 * eV,
        2.885 * eV, 2.954 * eV, 3.026 * eV, 3.102 * eV, 3.181 * eV, 3.265 * eV,
        3.353 * eV, 3.446 * eV, 3.545 * eV, 3.649 * eV, 3.760 * eV, 3.877 * eV,
        4.002 * eV, 4.136 * eV
    };

    // Water
    std::vector<G4double> refractiveIndex1 = {
        1.3435, 1.344,  1.3445, 1.345,  1.3455, 1.346,  1.3465, 1.347,
        1.3475, 1.348,  1.3485, 1.3492, 1.35,   1.3505, 1.351,  1.3518,
        1.3522, 1.3530, 1.3535, 1.354,  1.3545, 1.355,  1.3555, 1.356,
        1.3568, 1.3572, 1.358,  1.3585, 1.359,  1.3595, 1.36,   1.3608
    };
    std::vector<G4double> absorption = {
        3.448 * m,  4.082 * m,  6.329 * m,  9.174 * m,  12.346 * m, 13.889 * m,
        15.152 * m, 17.241 * m, 18.868 * m, 20.000 * m, 26.316 * m, 35.714 * m,
        45.455 * m, 47.619 * m, 52.632 * m, 52.632 * m, 55.556 * m, 52.632 * m,
        52.632 * m, 47.619 * m, 45.455 * m, 41.667 * m, 37.037 * m, 33.333 * m,
        30.000 * m, 28.500 * m, 27.000 * m, 24.500 * m, 22.000 * m, 19.500 * m,
        17.500 * m, 14.500 * m
    };

    // Material properties can be added as arrays. However, in this case it is
    // up to the user to make sure both arrays have the same number of elements.
    G4double scintilFastArray[]{ 1.0, 1.0 };
    G4double energyArray[]{ 2.034 * eV, 4.136 * eV };
    G4int lenArray = 2;

    std::vector<G4double> scintilSlow = {
        0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
        7.00, 6.00, 4.00, 3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00, 4.00,
        5.00, 6.00, 7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 5.00, 4.00
    };

    G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

    // Values can be added to the material property table individually.
    // With this method, spline interpolation cannot be set. Arguments
    // createNewKey and spline both take their default values of false.
    // Need to specify the number of entries (1) in the arrays, as an argument
    // to AddProperty.
    G4int numEntries = 1;
    myMPT1->AddProperty("RINDEX", &photonEnergy[0], &refractiveIndex1[0], numEntries);
    for(size_t i = 1; i < photonEnergy.size(); ++i)
    {
        myMPT1->AddEntry("RINDEX", photonEnergy[i], refractiveIndex1[i]);
    }

    // Check that group velocity is calculated from RINDEX
    if(myMPT1->GetProperty("RINDEX")->GetVectorLength() != myMPT1->GetProperty("GROUPVEL")->GetVectorLength())
    {
        G4ExceptionDescription ed;
        ed << "Error calculating group velocities. Incorrect number of entries in group velocity material property vector.";
        G4Exception("OpNovice::OpNoviceDetectorConstruction", "OpNovice001", FatalException, ed);
    }

    // Adding a property from two std::vectors. Argument createNewKey is false
    // and spline is true.
    myMPT1->AddProperty("ABSLENGTH", photonEnergy, absorption, false, true);

    // Adding a property using a C-style array.
    // Spline interpolation isn't used for scintillation.
    // Arguments spline and createNewKey both take default value false.
    myMPT1->AddProperty("SCINTILLATIONCOMPONENT1", energyArray, scintilFastArray, lenArray);
    myMPT1->AddProperty("SCINTILLATIONCOMPONENT2", photonEnergy, scintilSlow, false, true);
    myMPT1->AddConstProperty("SCINTILLATIONYIELD", 50. / MeV);
    myMPT1->AddConstProperty("RESOLUTIONSCALE", 1.0);
    myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 1. * ns);
    myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 10. * ns);
    myMPT1->AddConstProperty("SCINTILLATIONYIELD1", 0.8);
    myMPT1->AddConstProperty("SCINTILLATIONYIELD2", 0.2);


    std::vector<G4double> energy_water = {
        1.56962 * eV, 1.58974 * eV, 1.61039 * eV, 1.63157 * eV, 1.65333 * eV,
        1.67567 * eV, 1.69863 * eV, 1.72222 * eV, 1.74647 * eV, 1.77142 * eV,
        1.7971 * eV,  1.82352 * eV, 1.85074 * eV, 1.87878 * eV, 1.90769 * eV,
        1.93749 * eV, 1.96825 * eV, 1.99999 * eV, 2.03278 * eV, 2.06666 * eV,
        2.10169 * eV, 2.13793 * eV, 2.17543 * eV, 2.21428 * eV, 2.25454 * eV,
        2.29629 * eV, 2.33962 * eV, 2.38461 * eV, 2.43137 * eV, 2.47999 * eV,
        2.53061 * eV, 2.58333 * eV, 2.63829 * eV, 2.69565 * eV, 2.75555 * eV,
        2.81817 * eV, 2.88371 * eV, 2.95237 * eV, 3.02438 * eV, 3.09999 * eV,
        3.17948 * eV, 3.26315 * eV, 3.35134 * eV, 3.44444 * eV, 3.54285 * eV,
        3.64705 * eV, 3.75757 * eV, 3.87499 * eV, 3.99999 * eV, 4.13332 * eV,
        4.27585 * eV, 4.42856 * eV, 4.59258 * eV, 4.76922 * eV, 4.95999 * eV,
        5.16665 * eV, 5.39129 * eV, 5.63635 * eV, 5.90475 * eV, 6.19998 * eV
    };

    // Rayleigh scattering length is calculated by G4OpRayleigh

    // Mie: assume 100 times larger than the rayleigh scattering
    std::vector<G4double> mie_water = {
        167024.4 * m, 158726.7 * m, 150742 * m,   143062.5 * m, 135680.2 * m,
        128587.4 * m, 121776.3 * m, 115239.5 * m, 108969.5 * m, 102958.8 * m,
        97200.35 * m, 91686.86 * m, 86411.33 * m, 81366.79 * m, 76546.42 * m,
        71943.46 * m, 67551.29 * m, 63363.36 * m, 59373.25 * m, 55574.61 * m,
        51961.24 * m, 48527.00 * m, 45265.87 * m, 42171.94 * m, 39239.39 * m,
        36462.50 * m, 33835.68 * m, 31353.41 * m, 29010.30 * m, 26801.03 * m,
        24720.42 * m, 22763.36 * m, 20924.88 * m, 19200.07 * m, 17584.16 * m,
        16072.45 * m, 14660.38 * m, 13343.46 * m, 12117.33 * m, 10977.70 * m,
        9920.416 * m, 8941.407 * m, 8036.711 * m, 7202.470 * m, 6434.927 * m,
        5730.429 * m, 5085.425 * m, 4496.467 * m, 3960.210 * m, 3473.413 * m,
        3032.937 * m, 2635.746 * m, 2278.907 * m, 1959.588 * m, 1675.064 * m,
        1422.710 * m, 1200.004 * m, 1004.528 * m, 833.9666 * m, 686.1063 * m
    };

    // Mie: gforward, gbackward, forward backward ratio
    G4double mie_water_const[3] = { 0.99, 0.99, 0.8 };

    myMPT1->AddProperty("MIEHG", energy_water, mie_water, false, true);
    myMPT1->AddConstProperty("MIEHG_FORWARD", mie_water_const[0]);
    myMPT1->AddConstProperty("MIEHG_BACKWARD", mie_water_const[1]);
    myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO", mie_water_const[2]);

    G4cout << "Water G4MaterialPropertiesTable:" << G4endl;
    myMPT1->DumpTable();

    NaITl->SetMaterialPropertiesTable(myMPT1);
    // Set the Birks Constant for the Water scintillator
    NaITl->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);


    //SiO2->SetMaterialPropertiesTable(myMPT1);
    //// Set the Birks Constant for the Water scintillator
    //SiO2->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

    // Air
    std::vector<G4double> refractiveIndex2 = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                               1.0, 1.0, 1.0, 1.0 };

    G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
    myMPT2->AddProperty("RINDEX", photonEnergy, refractiveIndex2);

    G4cout << "Air G4MaterialPropertiesTable:" << G4endl;
    myMPT2->DumpTable();


    // /////////////////////////////////////////////////////// NaI(Tl) Scintillator //////////////////////////////////////////////////////////////
    // NaITl
    G4double HSActiveArea = 38.1*mm;
    G4Box* ActiveAreaSolid = new G4Box("ActiveArea", 38.1*mm, 38.1*mm, 5*cm);
    //G4LogicalVolume* ActiveAreaLogical = new G4LogicalVolume(ActiveAreaSolid, NaITl, "ActiveArea", 0, 0,0);
    G4LogicalVolume* ActiveAreaLogical = new G4LogicalVolume(ActiveAreaSolid, NaITl, "ActiveArea", 0, 0,0);
    G4VPhysicalVolume* ActiveAreaPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), ActiveAreaLogical, "ActiveArea", worldLogical, false, 0);

    // SiO2
    G4Box* SiO2Solid = new G4Box("SiO2", 39.95*mm, 1.5*mm, 5*cm);
    G4LogicalVolume* SiO2Logical = new G4LogicalVolume(SiO2Solid, SiO2, "SiO2", 0, 0,0);
    G4VPhysicalVolume* SiO2Physical = new G4PVPlacement(0, G4ThreeVector(0, -39.6, 0), SiO2Logical, "SiO2", worldLogical, false, 0);

    // MgO
    G4Element* elMg = new G4Element("Magnesium", "Mg", 12, 24.305 * g/mole);
    elO = new G4Element("Oxygen", "O", 8, 16.00 * g/mole);
    G4Material* MgO = new G4Material("MgO", 3.58 * g/cm3, 2);
    MgO->AddElement(elMg, 1);
    MgO->AddElement(elO, 1);
    G4UnionSolid* ActSiO2 = new G4UnionSolid( "ActSiO2" , ActiveAreaSolid, SiO2Solid, 0, G4ThreeVector(0,-39.6*mm,0));
    G4Box* box1 = new G4Box("box1", 39.95*mm, 39.025*mm,  5*cm);
    G4RotationMatrix* rm = new G4RotationMatrix(); rm->rotateX(0); rm->rotateY(0); rm->rotateZ(0);
    G4SubtractionSolid* MgOSolid = new G4SubtractionSolid( "MgO" , box1,ActSiO2, 0, G4ThreeVector(0,0,0));
    G4LogicalVolume* MgOLogical = new G4LogicalVolume(MgOSolid, MgO, "MgO", 0, 0,0);
    G4VPhysicalVolume* MgOPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), MgOLogical, "MgO", worldLogical, false, 0);

    // Al
    G4Element* elAl = new G4Element("Aluminum", "Al", 13, 26.981539 * g/mole);
    G4Material* Aluminum = new G4Material("Aluminum", 2.70 * g/cm3, 1);
    Aluminum->AddElement(elAl, 1);
    G4double HSAl = 40.45*mm;
    G4double HSMgO = 39.45*mm;
    G4Box* box2 = new G4Box("box2", 40.45*mm, 40.775*mm, 5*cm);
    G4Box* box3 = new G4Box("box3", 39.95*mm, 39.7*mm, 5*cm);
    rm = new G4RotationMatrix(); rm->rotateX(0); rm->rotateY(0); rm->rotateZ(0);
    G4SubtractionSolid* AlSolid = new G4SubtractionSolid( "Al" , box2, box3, rm, G4ThreeVector(0,0,0));
    G4LogicalVolume* AlLogical = new G4LogicalVolume(AlSolid, Aluminum, "Al", 0, 0,0);
    G4VPhysicalVolume* AlPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), AlLogical, "Al", worldLogical, false, 0);

    //// The Water Tank
    //G4Box* waterTank_box = new G4Box("Tank", 50*cm, 50*cm, 50*cm);
    //G4LogicalVolume* waterTank_log = new G4LogicalVolume(waterTank_box, waterMaterial, "Tank", 0, 0, 0);
    //G4VPhysicalVolume* waterTank_phys = new G4PVPlacement( 0, G4ThreeVector(), waterTank_log, "Tank", worldLogical, false, 0);

    //airMaterial->SetMaterialPropertiesTable(myMPT2);



    // The Water Tank
    //G4Box* waterTank_box = new G4Box("Tank", 50*cm, 50*cm, 50*cm);
    //G4LogicalVolume* waterTank_log = new G4LogicalVolume(waterTank_box, waterMaterial, "Tank", 0, 0, 0);
    //G4VPhysicalVolume* waterTank_phys = new G4PVPlacement( 0, G4ThreeVector(), waterTank_log, "Tank", worldLogical, false, 0);

    // ///////////////////////////////////////////////////// Surfaces /////////////////////////////////////////////////////////////////

    //// Water Tank
    //G4OpticalSurface* opWaterSurface = new G4OpticalSurface("WaterSurface");
    //opWaterSurface->SetType(dielectric_LUTDAVIS);
    //opWaterSurface->SetFinish(Rough_LUT);
    //opWaterSurface->SetModel(DAVIS);

    //G4LogicalBorderSurface* waterSurface = new G4LogicalBorderSurface("WaterSurface", ActiveAreaPhysical, SiO2Physical, opWaterSurface);
    //G4OpticalSurface* opticalSurface = dynamic_cast<G4OpticalSurface*>(waterSurface->GetSurface(ActiveAreaPhysical, SiO2Physical)->GetSurfaceProperty());
    //if(opticalSurface) opticalSurface->DumpInfo();

    //// Air Bubble
    //G4OpticalSurface* opAirSurface = new G4OpticalSurface("AirSurface");
    //opAirSurface->SetType(dielectric_dielectric);
    //opAirSurface->SetFinish(polished);
    //opAirSurface->SetModel(glisur);

    //G4LogicalSkinSurface* airSurface = new G4LogicalSkinSurface("AirSurface", MgOLogical, opAirSurface);

    //opticalSurface = dynamic_cast<G4OpticalSurface*>( airSurface->GetSurface(MgOLogical)->GetSurfaceProperty());
    //if(opticalSurface) opticalSurface->DumpInfo();

    //// Generate & Add Material Properties Table attached to the optical surfaces
    ////
    //std::vector<G4double> ephoton = { 2.034 * eV, 4.136 * eV };

    //// OpticalAirSurface
    //std::vector<G4double> reflectivity = { 0.3, 0.5 };
    //std::vector<G4double> efficiency   = { 0.8, 1.0 };

    //G4MaterialPropertiesTable* myST2 = new G4MaterialPropertiesTable();
    //myST2->AddProperty("REFLECTIVITY", ephoton, reflectivity);
    //myST2->AddProperty("EFFICIENCY", ephoton, efficiency);
    //G4cout << "Air Surface G4MaterialPropertiesTable:" << G4endl;
    //myST2->DumpTable();

    //opAirSurface->SetMaterialPropertiesTable(myST2);

    return worldPhysical;

}
