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
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialTable.hh"
#include "Randomize.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Ellipsoid.hh"
#include "G4UnionSolid.hh"
#include "G4EllipticalTube.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"

#include "G4TNuclearReactorGeometry.hh"
#include "G4TDetectorGeometry.hh"
#include "G4TPassiveProtonBeamLineGeometry.hh"

G4TCPPGeometryFormat::G4TCPPGeometryFormat(){}
G4TCPPGeometryFormat::~G4TCPPGeometryFormat(){}

G4VPhysicalVolume* G4TCPPGeometryFormat::ConstructPhysicalVolume(){

    // //////////////////////////////// Required Material creation ////////////////////////



    //G4TDetectorGeometry* NewGeom = new G4TDetectorGeometry();
    //WorldPhysicalVolume = NewGeom->ConstructDetector();
    //return WorldPhysicalVolume;

    //G4TPassiveProtonBeamLineGeometry* NewGeom = new G4TPassiveProtonBeamLineGeometry();
    //WorldPhysicalVolume = NewGeom->Construct();
    //return WorldPhysicalVolume;

    G4Material* matH2O;
    G4Material* soft;
    G4Material* skeleton;
    G4Material* lung;
    G4Material* adipose;
    G4Material* glandular;
    G4Material* adipose_glandular;
    G4Material* galactic;


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


    d = 1.0e-25 * g/cm3; // Very low density
    galactic = new G4Material("G4_Galactic", d, 1, kStateGas);
    galactic->AddElement(elH,1);

    // Water
    d = 1.000*g/cm3;
    matH2O = new G4Material("Water",d,2);
    matH2O->AddElement(elH,2);
    matH2O->AddElement(elO,1);
    matH2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

    // MIRD soft tissue
    d = 0.9869 *g/cm3;
    soft = new G4Material("soft_tissue",d,16);
    soft->AddElement(elH,0.1047);
    soft->AddElement(elC,0.2302);
    soft->AddElement(elN,0.0234);
    soft->AddElement(elO,0.6321);
    soft->AddElement(elNa,0.0013);
    soft->AddElement(elMg,0.00015);
    soft->AddElement(elP,0.0024);
    soft->AddElement(elS,0.0022);
    soft->AddElement(elCl,0.0014);
    soft->AddElement(elK,0.0021);
    soft->AddElement(elFe,0.000063);
    soft->AddElement(elZn,0.000032);
    soft->AddElement(elRb,0.0000057);
    soft->AddElement(elSr,0.00000034);
    soft->AddElement(elZr,0.000008);
    soft->AddElement(elPb,0.00000016);

    // MIRD Skeleton

    d = 1.4862*g/cm3;
    skeleton = new G4Material("skeleton",d,15);
    skeleton -> AddElement(elH,0.0704);
    skeleton -> AddElement(elC,0.2279);
    skeleton -> AddElement(elN,0.0387);
    skeleton -> AddElement(elO,0.4856);
    skeleton -> AddElement(elNa,0.0032);
    skeleton -> AddElement(elMg,0.0011);
    skeleton -> AddElement(elP,0.0694);
    skeleton -> AddElement(elS,0.0017);
    skeleton -> AddElement(elCl,0.0014);
    skeleton -> AddElement(elK,0.0015);
    skeleton -> AddElement(elCa,0.0991);
    skeleton -> AddElement(elFe,0.00008);
    skeleton -> AddElement(elZn,0.000048);
    skeleton -> AddElement(elSr,0.000032);
    skeleton -> AddElement(elPb,0.000011);

    // MIRD lung material
    d = 0.2958 *g/cm3;
    lung = new G4Material("lung_material", d,16);
    lung -> AddElement(elH, 0.1021);
    lung -> AddElement(elC, 0.1001);
    lung -> AddElement(elN,0.028);
    lung -> AddElement(elO,0.7596);
    lung -> AddElement(elNa,0.0019);
    lung -> AddElement(elMg,0.000074);
    lung -> AddElement(elP,0.00081);
    lung -> AddElement(elS,0.0023);
    lung -> AddElement(elCl,0.0027);
    lung -> AddElement(elK,0.0020);
    lung -> AddElement(elCa,0.00007);
    lung -> AddElement(elFe,0.00037);
    lung -> AddElement(elZn,0.000011);
    lung -> AddElement(elRb,0.0000037);
    lung -> AddElement(elSr,0.000000059);
    lung -> AddElement(elPb,0.00000041);

    G4double density_adipose = 0.93 *g/cm3;
    adipose = new G4Material("adipose", density_adipose,8);
    adipose -> AddElement(elH, 0.112);
    adipose -> AddElement(elC, 0.619);
    adipose -> AddElement(elN, 0.017);
    adipose -> AddElement(elO, 0.251);
    adipose -> AddElement(elS, 0.00025);
    adipose -> AddElement(elP, 0.00025);
    adipose -> AddElement(elK, 0.00025);
    adipose -> AddElement(elCa,0.00025);

    G4double density_glandular = 1.04 * g/cm3;
    glandular = new G4Material("glandular", density_glandular,8);
    glandular -> AddElement(elH, 0.1);
    glandular -> AddElement(elC,0.184);
    glandular -> AddElement(elN, 0.032);
    glandular -> AddElement(elO, 0.679);
    glandular -> AddElement(elS, 0.00125);
    glandular -> AddElement(elP, 0.00125);
    glandular -> AddElement(elK, 0.00125);
    glandular -> AddElement(elCa,0.00125);


    d = (density_adipose * 0.5) + (density_glandular * 0.5);
    adipose_glandular = new G4Material("adipose_glandular", d, 2);
    adipose_glandular -> AddMaterial(adipose, 0.5);
    adipose_glandular -> AddMaterial(glandular, 0.5);

    /*

    // for C++ geometry

        // materials construction
        G4Element* elN = new G4Element("Nitrogen","N",Z = 7., 14.01*g/mole);
        G4Element* elO = new G4Element("Oxygen","O",Z = 8., A = 16.00*g/mole);
        G4Material* matAir = new G4Material("Air",1.290*mg/cm3,2);
        matAir->AddElement(elN,0.7);
        matAir->AddElement(elO,0.3);

        // world construction
        G4Box* worldSol = new G4Box("world", 0.25*m, 0.15*m, 1.*m);
        G4LogicalVolume* worldLog = new G4LogicalVolume(worldSol,matAir, "World", 0, 0,0);
        G4VPhysicalVolume* worldPhy = new G4PVPlacement(0,G4ThreeVector(), "World", worldLog,
                                                        0, false, 0);

        // brain construction
        G4Ellipsoid* brainSol = new G4Ellipsoid("Brain", 6.*cm, 9.*cm, 6.5*cm);
        G4LogicalVolume* brainLog =  new G4LogicalVolume(brainSol, soft, "Brain", 0, 0, 0);
        G4VPhysicalVolume* brainPhy = new G4PVPlacement(0,
                                                         G4ThreeVector(0.*cm, 0.*cm, 8.75 * cm),
                                                         "physicalBrain", brainLog, worldPhy,
                                                         false, 0, true);

    */


    // Air
    d = 1.290*mg/cm3;
    G4Material* matAir = new G4Material("Air",d,2);
    matAir->AddElement(elN,0.7);
    matAir->AddElement(elO,0.3);

    // //////////////////////////////// World Construction ////////////////////////

    G4double worldSize = 1.5 *m ;
    G4Box* world = new G4Box("world", 2.5*m, 2.5*m, 2.5*m);

    G4LogicalVolume* logicWorld = new G4LogicalVolume(world,
                                                      galactic,
                                                      "World", 0, 0,0);

    WorldPhysicalVolume = new G4PVPlacement(0,G4ThreeVector(),
                                            "World",
                                            logicWorld,
                                            0,
                                            false,
                                            0);

    // //////////////////////////////// Head Construction ////////////////////////

    // MIRD male model
    // Ellipsoid
    G4double ax = 7.0 * cm;
    G4double by = 10.0 * cm;
    G4double cz = 8.50 * cm;
    G4double zcut1 = 0.0 * cm;
    G4double zcut2 = 8.5 * cm;

    G4Ellipsoid* head1 = new G4Ellipsoid("Head1", ax, by, cz, zcut1, zcut2);

    G4double dx = 7.0 * cm;
    G4double dy = 10.0 * cm;
    G4double dz = 7.75 * cm;


    G4EllipticalTube* head2 = new G4EllipticalTube("Head2", dx, dy, dz);

    G4UnionSolid* head = new G4UnionSolid("Head",head2,head1,
                                          0, // Rotation
                                          G4ThreeVector(0.* cm, 0.*cm, 7.7500 * cm) );

    G4LogicalVolume* logicHead = new G4LogicalVolume(head, soft,"Head", 0, 0,0);
    // Define rotation and position here!
    G4RotationMatrix* rm = new G4RotationMatrix();
    rm->rotateX(180.*degree);
    rm->rotateY(180.*degree);

    G4VPhysicalVolume* physHead = new G4PVPlacement(rm,
                                                    //G4ThreeVector(0.* cm,0.*cm, 70.75*cm), //FA
                                                    G4ThreeVector(0.* cm,0.*cm, 77.75*cm),
                                                    "physicalHead",
                                                    logicHead,
                                                    WorldPhysicalVolume,
                                                    false,
                                                    0, true);

    //G4cout << " Head : Volume=" << logicHead->GetSolid()->GetCubicVolume()/cm3 << " cm^3"
    //     << " -Density of Material=" << logicHead->GetMaterial()->GetDensity()*cm3/g << " g/cm^3"
    //   << " -Mass=" << (logicHead->GetSolid()->GetCubicVolume()/cm3)*(logicHead->GetMaterial()->GetDensity()*cm3/g)/gram << " g"
    // << " -Material=" << logicHead->GetMaterial()->GetName() << G4endl;

    // //////////////////////////////// Brain Construction ////////////////////////

    ax = 6. * cm;
    by= 9. * cm;
    cz = 6.5 * cm;

    G4Ellipsoid* brain = new G4Ellipsoid("Brain", ax, by, cz);

    G4LogicalVolume* logicBrain =  new G4LogicalVolume(brain, soft, "Brain", 0, 0, 0);

    // Define rotation and position here!
    G4VPhysicalVolume* physBrain = new G4PVPlacement(0,
                                                     G4ThreeVector(0.*cm, 0.*cm, 8.75 * cm),
                                                     "physicalBrain",
                                                     logicBrain,
                                                     physHead,
                                                     false,
                                                     0, true);

    //G4cout << " Brain : Volume=" << logicBrain->GetSolid()->GetCubicVolume()/cm3 << " cm^3"
    //     << " -Density of Material=" << logicBrain->GetMaterial()->GetDensity()*cm3/g << " g/cm^3"
    //   << " -Mass=" << (logicBrain->GetSolid()->GetCubicVolume()/cm3)*(logicBrain->GetMaterial()->GetDensity()*cm3/g)/gram << " g"
    // << " -Material=" << logicBrain->GetMaterial()->GetName() << G4endl;

    // //////////////////////////////// Skull Construction ////////////////////////

    // Outer cranium
    ax = 6.8 * cm;//a out skull
    by = 9.8 * cm; // bout
    cz = 8.3 * cm; //cout

    G4Ellipsoid* craniumOut =  new G4Ellipsoid("CraniumOut", ax, by, cz);

    ax = 6. * cm; //a in
    by = 9. * cm; //b in
    cz= 6.5 * cm; // cin

    G4Ellipsoid* craniumIn =  new G4Ellipsoid("CraniumIn", ax, by, cz);


    G4SubtractionSolid* cranium =  new G4SubtractionSolid("Cranium",
                                                          craniumOut,
                                                          craniumIn,0,
                                                          G4ThreeVector(0.0, 0.0,1. * cm));

    G4LogicalVolume* logicSkull = new G4LogicalVolume(cranium, skeleton,
                                                      "Skull",
                                                      0, 0, 0);

    // Define rotation and position here!
    G4VPhysicalVolume* physSkull = new G4PVPlacement(0,
                                                     G4ThreeVector(0., 0.,7.75 * cm),
                                                     "physicalSkull",
                                                     logicSkull,
                                                     physHead,
                                                     false,
                                                     0, true);


    // //////////////////////////////// Thyroid Construction ////////////////////////

    G4double z= 4.20*cm; //c thickness = c,
    G4double rmin= 0. * cm;
    G4double rmax= 1.85 *cm; //Rmax
    G4double startphi = 0. * degree;
    G4double deltaphi= 180. * degree; // y< y0

    G4Tubs* LobOfThyroidOut = new G4Tubs("LobOfThyroidOut",
                                         rmin, rmax,z/2.,
                                         startphi, deltaphi);

    z= 4.50*cm; // c thickness + something
    rmax= 0.83 * cm; //r
    deltaphi= 360. * degree;
    G4Tubs* LobOfThyroidIn = new G4Tubs("LobOfThyroidIn",
                                        rmin, rmax,z/2.,
                                        startphi, deltaphi);

    G4double xx = 3.72*cm;
    G4double yy= 3.72*cm;
    G4double zz= 20.00*cm;
    G4Box* SubtrThyroid = new G4Box("SubtrThyroid",
                                    xx/2., yy/2., zz/2.);
    // subtraction of the two tubs
    G4SubtractionSolid* FirstThyroid = new G4SubtractionSolid("FirstThyroid",
                                                              LobOfThyroidOut,
                                                              LobOfThyroidIn);

    G4RotationMatrix* relative_matrix = new G4RotationMatrix();
    relative_matrix -> rotateX(-50.* degree);

    G4SubtractionSolid* SecondThyroid = new G4SubtractionSolid("SecondThyroid",
                                                               FirstThyroid,
                                                               SubtrThyroid,
                                                               relative_matrix,
                                                               G4ThreeVector(0.0 *cm,0.0 *cm, 4.20*cm));

    G4RotationMatrix* relative_matrix_2 = new G4RotationMatrix();
    relative_matrix_2 -> rotateX(50.* degree);

    //G4SubtractionSolid* thyroid = new G4SubtractionSolid("SecondThyroid", SecondThyroid, SubtrThyroid, relative_matrix_2, G4ThreeVector(0.0 *cm,0.0 *cm, -5.40*cm));



    //G4LogicalVolume* logicThyroid = new G4LogicalVolume(thyroid, soft, "Thyroid", 0, 0, 0);

    rm = new G4RotationMatrix();
    rm->rotateZ(180.*degree);
    //G4VPhysicalVolume* physThyroid = new G4PVPlacement(rm, G4ThreeVector(0.0*cm,-3.91*cm, -5.65*cm),//y0"physicalThyroid", logicThyroid, physHead, false, 0,true);


    // //////////////////////////////// upperspine Construction ////////////////////////


    dx = 2. *cm;
    dy = 2.5 *cm;
    dz = 4.25*cm;

    G4EllipticalTube* upperSpine = new G4EllipticalTube("UpperSpine",dx, dy, dz);

    xx = 20. * cm;
    yy = 10. * cm;
    zz = 5. * cm;

    G4Box* subtraction = new G4Box("box", xx/2., yy/2., zz/2.);

    G4RotationMatrix* matrix = new G4RotationMatrix();
    matrix -> rotateX(-25.* deg);

    G4SubtractionSolid* upper_spine = new G4SubtractionSolid("upperspine",upperSpine, subtraction,
                                                             matrix, G4ThreeVector(0., -2.5 * cm, 5.5* cm));

    G4LogicalVolume* logicUpperSpine = new G4LogicalVolume(upper_spine, skeleton,
                                                           "UpperSpine",
                                                           0, 0, 0);
    // Define rotation and position here!
    G4VPhysicalVolume* physUpperSpine = new G4PVPlacement(0,
                                                          G4ThreeVector(0.0, 5.5 *cm, -3.5 *cm),
                                                          "physicalUpperSpine",
                                                          logicUpperSpine,
                                                          physHead,
                                                          false,
                                                          0, true);

    // //////////////////////////////// Trunk Construction ////////////////////////

    dx = 20. * cm;
    dy = 10. * cm;
    dz = 35. * cm;

    G4EllipticalTube* trunk = new G4EllipticalTube("Trunk",dx, dy, dz);

    G4LogicalVolume* logicTrunk = new G4LogicalVolume(trunk, soft,
                                                      "Trunk",
                                                      0, 0, 0);
    // Define rotation and position here!
    rm = new G4RotationMatrix();
    rm->rotateX(180.*degree);
    rm->rotateY(180.*degree);
    G4VPhysicalVolume* physTrunk = new G4PVPlacement(rm,
                                                     //G4ThreeVector(0.* cm, 0. *cm, 28.*cm), //FA
                                                     G4ThreeVector(0.* cm, 0. *cm, 35.*cm),
                                                     "physicalTrunk",
                                                     logicTrunk,
                                                     WorldPhysicalVolume,
                                                     false,
                                                     0, true);



    // //////////////////////////////// UrinaryBladder Construction ////////////////////////

    ax = 4.958*cm;
    by= 3.458 *cm;
    cz= 3.458 *cm;

    G4Ellipsoid* bladder = new G4Ellipsoid("bladder_out",ax, by, cz);

    ax = 4.706 * cm;
    by = 3.206 * cm;
    cz = 3.206 * cm;
    G4Ellipsoid* inner = new G4Ellipsoid("innerBladder", ax, by, cz);

    G4SubtractionSolid* totalBladder = new G4SubtractionSolid("bladder", bladder, inner);

    G4LogicalVolume* logicUrinaryBladder = new G4LogicalVolume(totalBladder, soft,
                                                               "UrinaryBladder",
                                                               0, 0, 0);

    // Define rotation and position here!
    G4VPhysicalVolume* physUrinaryBladder = new G4PVPlacement(0,G4ThreeVector(0 *cm, -4.5 *cm,-27. *cm),
                                                              "physicalUrinaryBladder",
                                                              logicUrinaryBladder,
                                                              physTrunk,
                                                              false,
                                                              0, true);

    // //////////////////////////////// Thymus Construction ////////////////////////

    ax= 3. *cm;
    by= 0.5*cm;
    cz= 4.*cm;

    G4Ellipsoid* Thymus = new G4Ellipsoid("Thymus",
                                          ax, by, cz);

    G4LogicalVolume* logicThymus = new G4LogicalVolume(Thymus,
                                                       soft,
                                                       "Thymus",
                                                       0, 0, 0);

    // Define rotation and position here!
    G4VPhysicalVolume* physThymus = new G4PVPlacement(0,
                                                      G4ThreeVector(2.*cm,-6.*cm, 25.5*cm),
                                                      "physicalThymus",
                                                      logicThymus,
                                                      physTrunk,
                                                      false,
                                                      0, true);

    // //////////////////////////////// Heart Construction ////////////////////////

    ax= 4.00* cm;
    by= 4.00 *cm;
    cz= 7.00 *cm;
    zcut1= -7.00 *cm;
    zcut2= 0.0 *cm;

    G4Ellipsoid* heart1 =  new G4Ellipsoid("Heart1",ax, by, cz, zcut1, zcut2);

    rmin =0.*cm;
    rmax = 3.99*cm;
    startphi = 0. * degree;
    deltaphi = 360. * degree;
    G4double starttheta = 0. * degree;
    G4double deltatheta = 90. * degree;

    G4Sphere* heart2 = new G4Sphere("Heart2", rmin,rmax,
                                    startphi,   deltaphi,
                                    starttheta, deltatheta);

    G4UnionSolid* heart = new G4UnionSolid("Heart", heart1, heart2);

    G4LogicalVolume* logicHeart = new G4LogicalVolume(heart, soft,
                                                      "Heart",
                                                      0, 0, 0);

    // Define rotation and position here!
    rm = new G4RotationMatrix();
    rm->rotateY(25.*degree);
    G4VPhysicalVolume* physHeart = new G4PVPlacement(rm,G4ThreeVector(0.0,-3.0*cm, 15.32 *cm),
                                                     "physicalHeart",
                                                     logicHeart,
                                                     physTrunk,
                                                     false,
                                                     0,true);

    // //////////////////////////////// RightAdrenal Construction ////////////////////////

    ax= 1.5 *cm; //a
    by= 0.5 *cm; //b
    cz= 5.0 *cm; //c

    G4VSolid* rightAdrenal = new G4Ellipsoid("OneRightAdrenal",ax, by, cz, 0. *cm, cz);


    G4LogicalVolume* logicRightAdrenal = new G4LogicalVolume(rightAdrenal,
                                                             soft,
                                                             "Adrenals",
                                                             0, 0, 0);

    G4VPhysicalVolume* physRightAdrenal = new G4PVPlacement(0 ,G4ThreeVector(-4.5*cm,  // xo
                                                                             6.5 *cm, //yo
                                                                             3. *cm),//zo
                                                            "physicalRightAdrenal", logicRightAdrenal,
                                                            physTrunk,
                                                            false,
                                                            0, true);

    // //////////////////////////////// LeftAdrenal Construction ////////////////////////

    ax= 1.5 *cm; //a
    by= 0.5 *cm; //b
    cz= 5.0 *cm; //c

    G4VSolid* leftAdrenal = new G4Ellipsoid("OneLeftAdrenal",ax, by, cz, 0. *cm, cz);


    G4LogicalVolume* logicLeftAdrenal = new G4LogicalVolume(leftAdrenal,
                                                            soft,
                                                            "Adrenals",
                                                            0, 0, 0);

    G4VPhysicalVolume* physLeftAdrenal = new G4PVPlacement(0,G4ThreeVector(4.5*cm,  // xo
                                                                           6.5 *cm, //yo
                                                                           3. *cm),//zo
                                                           "physicalLeftAdrenal",
                                                           logicLeftAdrenal,
                                                           physTrunk,
                                                           false,
                                                           0, true);

    // //////////////////////////////// RightArmBone Construction ////////////////////////

    dx = 1.4 * cm;//a
    dy = 2.7 * cm;//b
    //  dz= 46. * cm;//z0

    //G4EllipticalCone* arm = new G4EllipticalCone("OneRightArmBone",dx/2.,dy/2.,dz, 34.5 *cm);
    G4EllipticalTube* rightArm = new G4EllipticalTube("OneRightArmBone",dx,dy,34.5 *cm);

    G4LogicalVolume* logicRightArmBone = new G4LogicalVolume(rightArm,
                                                             skeleton,
                                                             "RightArmBone",
                                                             0, 0,0);

    rm = new G4RotationMatrix();
    rm->rotateX(180.*degree);
    G4VPhysicalVolume* physRightArmBone = new G4PVPlacement(rm,
                                                            G4ThreeVector(-18.4 * cm, 0.0, -0.5*cm),
                                                            //-x0
                                                            "physicalRightArmBone",
                                                            logicRightArmBone,
                                                            physTrunk,
                                                            false,0,true);

    // //////////////////////////////// LeftArmBone Construction ////////////////////////

    dx = 1.4 * cm;//a
    dy = 2.7 * cm;//b
    //dz = 46. * cm;//z0

    //G4EllipticalCone* arm = new G4EllipticalCone("OneLeftArmBone",dx/2.,dy/2.,dz, 34.5 *cm);
    G4EllipticalTube* leftArm = new G4EllipticalTube("OneLeftArmBone",dx,dy,34.5 *cm);

    G4LogicalVolume* logicLeftArmBone = new G4LogicalVolume(leftArm,
                                                            skeleton,
                                                            "LeftArmBone",
                                                            0, 0,0);

    rm = new G4RotationMatrix();
    rm->rotateX(180.*degree);
    G4VPhysicalVolume* physLeftArmBone = new G4PVPlacement(rm,
                                                           G4ThreeVector(18.4 * cm, 0.0, -0.5*cm),
                                                           //-x0
                                                           "physicalLeftArmBone",
                                                           logicLeftArmBone,
                                                           physTrunk,
                                                           false,0,true);


    // //////////////////////////////// Uterus Construction ////////////////////////

     ax= 2.5*cm; //a
     by= 1.5*cm; //c
     cz= 5.*cm; //b

     zcut1= -5.* cm; //-b
     zcut2= 2.5*cm; //y1-y0

    G4Ellipsoid* uterus = new G4Ellipsoid("Uterus",
                                          ax, by, cz,
                                          zcut1, zcut2);

    G4LogicalVolume* logicUterus = new G4LogicalVolume(uterus,
                                                       soft,
                                                       "Uterus",
                                                       0, 0, 0);


    // Define rotation and position here!
    rm = new G4RotationMatrix();
    rm->rotateX(90.*degree);
    G4VPhysicalVolume* physUterus = new G4PVPlacement(rm,
                                                      G4ThreeVector(0. *cm, 2*cm,-21 *cm),
                                                      "physicalUterus", //y0
                                                      logicUterus,
                                                      physTrunk,
                                                      false,
                                                      0, true);

    // //////////////////////////////// RightOvary Construction ////////////////////////

     ax= 1. *cm;
     by= 0.5*cm;
     cz= 2.*cm;

    G4Ellipsoid* OneOvary = new G4Ellipsoid("OneOvary",
                                            ax, by, cz);

    G4LogicalVolume* logicRightOvary = new G4LogicalVolume(OneOvary,
                                                           soft,
                                                           "RightOvary",
                                                           0, 0, 0);

    // Define rotation and position here!
    G4VPhysicalVolume* physRightOvary = new G4PVPlacement(0,
                                                          G4ThreeVector(6. *cm,0.5*cm, -20*cm),
                                                          "physicalRightOvary",
                                                          logicRightOvary,
                                                          physTrunk,
                                                          false,
                                                          0, true);

    // //////////////////////////////// LeftOvary Construction ////////////////////////

     ax= 1. *cm;
     by= 0.5*cm;
     cz= 2.*cm;

    OneOvary = new G4Ellipsoid("OneOvary",
                                            ax, by, cz);


    G4LogicalVolume* logicLeftOvary = new G4LogicalVolume(OneOvary,
                                                          soft,
                                                          "LeftOvary",
                                                          0, 0, 0);

    // Define rotation and position here!
    G4VPhysicalVolume* physLeftOvary = new G4PVPlacement(0,
                                                         G4ThreeVector(-6. *cm,0.5*cm, -20*cm),
                                                         "physicalLeftOvary",
                                                         logicLeftOvary,
                                                         physTrunk,
                                                         false,
                                                         0, true);

    // //////////////////////////////// RightBreast Construction ////////////////////////

     ax= 4.95* cm;
     by= 4.35* cm;
     cz= 4.15*cm;

    G4Ellipsoid* oneRightBreast = new G4Ellipsoid("OneRightBreast",
                                                  ax, by, cz);

     dx= 20.* cm;
     dy= 10.* cm;
     dz= 35.* cm;

    G4EllipticalTube* Trunk = new G4EllipticalTube("Trunk",dx, dy, dz );

    G4RotationMatrix* rm_relative = new G4RotationMatrix();
    rm_relative -> rotateX(90. * degree);

    G4SubtractionSolid* breast = new G4SubtractionSolid("RightBreast",
                                                        oneRightBreast,
                                                        Trunk,
                                                        rm_relative,
                                                        G4ThreeVector(10.*cm,
                                                                      0.0*cm,
                                                                      -8.66*cm));


    G4LogicalVolume* logicRightBreast = new G4LogicalVolume(breast, soft,"RightBreast", 0, 0,0);


    // Define rotation and position here!
     rm = new G4RotationMatrix();
    rm->rotateX(90.*degree);
    rm->rotateY(0.*degree);
    rm->rotateZ(-16.*degree);
    G4VPhysicalVolume* physRightBreast = new G4PVPlacement(rm,
                                                           G4ThreeVector(-10.*cm,
                                                                         9.1 *cm,
                                                                         52.* cm),
                                                           "physicalRightBreast",
                                                           logicRightBreast,
                                                           WorldPhysicalVolume,
                                                           false,
                                                           0, true);


    // //////////////////////////////// LeftBreast Construction ////////////////////////

     ax= 4.95* cm;
     by= 4.35* cm;
     cz= 4.15*cm;

    G4Ellipsoid* oneLeftBreast = new G4Ellipsoid("OneLeftBreast",
                                                 ax, by, cz);

     dx= 20.* cm;
     dy= 10.* cm;
     dz= 35.*cm;

    Trunk = new G4EllipticalTube("Trunk",dx, dy, dz );



     rm_relative = new G4RotationMatrix();
    rm_relative -> rotateX(90. * degree);

    breast = new G4SubtractionSolid("LeftBreast",
                                                        oneLeftBreast,
                                                        Trunk,
                                                        rm_relative,
                                                        G4ThreeVector(-10.*cm,
                                                                      0.0*cm,
                                                                      -8.66*cm));


    G4LogicalVolume* logicLeftBreast = new G4LogicalVolume(breast, soft,
                                                           "LeftBreast", 0, 0,0);


    // Define rotation and position here!
     rm = new G4RotationMatrix();
    rm->rotateX(90.*degree);
    rm->rotateY(0.*degree);
    rm->rotateZ(16.*degree);
    G4VPhysicalVolume* physLeftBreast = new G4PVPlacement(rm,G4ThreeVector(10*cm, 9.1 *cm, 52.* cm),
                                                          "physicalLeftBreast",
                                                          logicLeftBreast,
                                                          WorldPhysicalVolume,
                                                          false,
                                                          0, true);


    // //////////////////////////////// RightClavicle Construction ////////////////////////

    G4double rMin = 0*cm;
    G4double rMax = 0.7883*cm;
    G4double rTor = 10*cm;
    G4double pSPhi = 201.75*degree;
    G4double pDPhi = 0.7*rad;

    G4Torus* clavicle = new G4Torus("Clavicle",rMin,rMax,rTor,pSPhi,pDPhi);

    G4LogicalVolume* logicRightClavicle = new G4LogicalVolume(clavicle,
                                                              skeleton,
                                                              "RightClavicle",
                                                              0, 0, 0);
    G4VPhysicalVolume* physRightClavicle = new G4PVPlacement(0,
                                                             G4ThreeVector(0.*cm,2.*cm,33.25*cm),
                                                             "physicalRightClavicle",
                                                             logicRightClavicle,
                                                             physTrunk,
                                                             false,
                                                             0, true);

    // //////////////////////////////// LeftClavicle Construction ////////////////////////

    rMin = 0*cm;
    rMax = 0.7883*cm;
    rTor = 10*cm;
    pSPhi = 298.15*degree;
    pDPhi = 0.7*rad;


    clavicle = new G4Torus("Clavicle",rMin,rMax,rTor,pSPhi,pDPhi);

    G4LogicalVolume* logicLeftClavicle = new G4LogicalVolume(clavicle,
                                                             skeleton,
                                                             "LeftClavicle",
                                                             0, 0, 0);


    G4VPhysicalVolume* physLeftClavicle = new G4PVPlacement(0,
                                                            G4ThreeVector(0.*cm,2. *cm,33.25*cm),
                                                            "physicalLeftClavicle",
                                                            logicLeftClavicle,
                                                            physTrunk,
                                                            false,
                                                            0, true);

    // //////////////////////////////// RightKidney Construction ////////////////////////

    ax= 4.5 *cm; //a
    by= 1.5 *cm; //b
    cz= 5.5 *cm; //c

    G4VSolid* oneRightKidney = new G4Ellipsoid("OneRightKidney",ax, by, cz);

    xx = 6. * cm;
    yy = 12.00*cm;
    zz = 12.00*cm;

    G4VSolid* subtrRightKidney = new G4Box("SubtrRightKidney",xx/2., yy/2., zz/2.);

    G4SubtractionSolid* kidney = new G4SubtractionSolid("RightKidney",
                                                        oneRightKidney,
                                                        subtrRightKidney,
                                                        0,
                                                        G4ThreeVector(6. *cm, // x0
                                                                      0.0 *cm,
                                                                      0.0 * cm));

    G4LogicalVolume* logicRightKidney = new G4LogicalVolume(kidney,
                                                            soft,
                                                            "Kidneys",
                                                            0, 0, 0);

    G4VPhysicalVolume* physRightKidney = new G4PVPlacement(0 ,G4ThreeVector(-6.*cm,  // xo
                                                                            6. *cm, //yo
                                                                            -2.50 *cm),//zo
                                                           "physicalRightKidney", logicRightKidney,
                                                           physTrunk,
                                                           false,
                                                           0, true);


    // //////////////////////////////// LeftKidney Construction ////////////////////////

    ax= 4.5 *cm; //a
    by= 1.5 *cm; //b
    cz= 5.5 *cm; //c

    G4VSolid* oneLeftKidney = new G4Ellipsoid("OneLeftKidney",ax, by, cz);

    xx = 6. * cm;
    yy = 12.00*cm;
    zz = 12.00*cm;
    G4VSolid* subtrLeftKidney = new G4Box("SubtrLeftKidney",xx/2., yy/2., zz/2.);

    G4SubtractionSolid* leftKidney = new G4SubtractionSolid("Kidneys",
                                                            oneLeftKidney,
                                                            subtrLeftKidney,
                                                            0,
                                                            G4ThreeVector(-6. *cm, // x0
                                                                          0.0 *cm,
                                                                          0.0 * cm));

    G4LogicalVolume* logicLeftKidney = new G4LogicalVolume(leftKidney,
                                                           soft,
                                                           "Kidneys",
                                                           0, 0, 0);

    G4VPhysicalVolume* physLeftKidney = new G4PVPlacement(0 ,G4ThreeVector(6.*cm,  // xo
                                                                           6. *cm, //yo
                                                                           -2.50 *cm),//zo
                                                          "physicalLeftKidney", logicLeftKidney,
                                                          physTrunk,
                                                          false,
                                                          0, true);

    // //////////////////////////////// RightLung Construction ////////////////////////

    ax = 5. *cm; //a
    by = 7.5 *cm; //b
    cz = 24.*cm; //c
    zcut1 = 0.0 *cm;
    zcut2=24. *cm;

    G4Ellipsoid* oneLung = new G4Ellipsoid("OneLung",ax, by, cz, zcut1,zcut2);

    ax= 5.*cm;
    by= 7.5*cm;
    cz= 24.*cm;


    G4Ellipsoid* subtrLung = new G4Ellipsoid("subtrLung",ax, by, cz);

    // y<0

    dx = 5.5* cm;
    dy = 8.5 * cm;
    dz = 24. * cm;

    G4Box* box = new G4Box("Box", dx, dy, dz);


    G4SubtractionSolid* section = new G4SubtractionSolid("BoxSub", subtrLung, box, 0, G4ThreeVector(0.*cm, 8.5* cm, 0.*cm));
    //G4SubtractionSolid* section2 = new G4SubtractionSolid("BoxSub2", subtrLung, box, 0, G4ThreeVector(0.*cm, -8.5* cm, 0.*cm));

    G4SubtractionSolid* lung1 =  new G4SubtractionSolid("Lung1", oneLung,
                                                        section,
                                                        0, G4ThreeVector(6.*cm,0*cm,0.0*cm));


    G4LogicalVolume* logicRightLung = new G4LogicalVolume(lung1,lung,
                                                          "Lungs", 0, 0, 0);

    rm = new G4RotationMatrix();
    rm->rotateZ(-360.*degree);
    G4VPhysicalVolume* physRightLung = new G4PVPlacement(rm,G4ThreeVector(-8.50 *cm, 0.0*cm, 8.5*cm),
                                                         "physicalRightLung",
                                                         logicRightLung,
                                                         physTrunk,
                                                         false,
                                                         0, true);

    // //////////////////////////////// LeftLung Construction ////////////////////////


    ax = 5. *cm; //a
    by = 7.5 *cm; //b
    cz = 24.*cm; //c
    zcut1 = 0.0 *cm;
    zcut2=24. *cm;

    oneLung = new G4Ellipsoid("OneLung",ax, by, cz, zcut1,zcut2);

    ax= 5.*cm;
    by= 7.5*cm;
    cz= 24.*cm;


    subtrLung = new G4Ellipsoid("subtrLung",ax, by, cz);

    // y<0

    dx = 5.5* cm;
    dy = 8.5 * cm;
    dz = 24. * cm;

    box = new G4Box("Box", dx, dy, dz);


    //G4SubtractionSolid* section = new G4SubtractionSolid("BoxSub", subtrLung, box, 0, G4ThreeVector(0.*cm, 8.5* cm, 0.*cm));
    G4SubtractionSolid* section2 = new G4SubtractionSolid("BoxSub2", subtrLung, box, 0, G4ThreeVector(0.*cm, 8.5* cm, 0.*cm));

    //G4SubtractionSolid* lung1 =  new G4SubtractionSolid("Lung1", oneLung,
    //				       section,
    //				       0, G4ThreeVector(6.*cm,0*cm,0.0*cm));

    G4SubtractionSolid* lung2 =  new G4SubtractionSolid("Lung2", oneLung,
                                                        section2,
                                                        0, G4ThreeVector(-6.*cm,0*cm,0.0*cm));

    // G4RotationMatrix* matrix = new G4RotationMatrix();
    // matrix->rotateX(180. * degree);
    //matrix ->rotateZ(180.*degree);
    //matrix -> rotateY(180.* degree);

    //G4UnionSolid* lungs = new G4UnionSolid("Lungs", lung1, lung2, matrix, G4ThreeVector(17*cm, 0., 0.));


    G4LogicalVolume* logicLeftLung = new G4LogicalVolume(lung2,lung,
                                                         "Lungs", 0, 0, 0);


    G4VPhysicalVolume* physLeftLung = new G4PVPlacement(0,G4ThreeVector(8.50 *cm, 0.0*cm, 8.5*cm),
                                                        "physicalLeftLung",
                                                        logicLeftLung,
                                                        physTrunk,
                                                        false, 0, true);

    // //////////////////////////////// RightScapula Construction ////////////////////////

    G4double ax_in = 17.* cm;
    G4double by_in = 9.8* cm;
    G4double ax_out = 19.*cm;
    G4double by_out = 9.8*cm;
    dz= 16.4* cm;


    G4EllipticalTube* inner_scapula = new G4EllipticalTube("ScapulaIn", ax_in, by_in, (dz+ 1.*cm)/2);
    G4EllipticalTube* outer_scapula = new G4EllipticalTube("ScapulaOut", ax_out, by_out, dz/2);


    subtraction = new G4Box("subtraction",ax_out, ax_out, ax_out);

    xx = -ax_out * 0.242 ; //(sin 14deg)
    yy  = - ax_out * 0.97; // (cos 14 deg)

    rm = new G4RotationMatrix();
    rm -> rotateZ(14.* degree);



    G4SubtractionSolid* scapula_first =  new G4SubtractionSolid("Scapula_first",
                                                                outer_scapula,
                                                                subtraction,
                                                                rm,
                                                                G4ThreeVector(xx, yy, 0. *cm));

    G4double xx2 = ax_out * 0.62470 ; //(cos 51.34deg)
    G4double yy2 = ax_out * 0.78087; // (sin 51.34 deg)

    G4RotationMatrix* rm2 = new G4RotationMatrix();
    rm2 -> rotateZ(38.6598* degree);


    G4SubtractionSolid* scapula_bone =  new G4SubtractionSolid("Scapula",
                                                               scapula_first,
                                                               subtraction,
                                                               rm2,
                                                               G4ThreeVector(xx2, yy2, 0. *cm));

    G4SubtractionSolid* scapula =  new G4SubtractionSolid("Scapula",
                                                          scapula_bone,
                                                          inner_scapula);

    G4LogicalVolume* logicRightScapula = new G4LogicalVolume(scapula,
                                                             skeleton,
                                                             "RightScapula",
                                                             0, 0, 0);

    G4VPhysicalVolume* physRightScapula = new G4PVPlacement(0,
                                                            G4ThreeVector(0. * cm, 0. * cm, 24.1 *cm),
                                                            "physicalRightScapula",
                                                            logicRightScapula,
                                                            physTrunk,
                                                            false,
                                                            0, true);

    // //////////////////////////////// LeftScapula Construction ////////////////////////

    ax_in = 17.* cm;
    by_in = 9.8* cm;
    ax_out = 19.*cm;
    by_out = 9.8*cm;
    dz= 16.4* cm;


    inner_scapula = new G4EllipticalTube("ScapulaIn", ax_in, by_in, (dz+ 1.*cm)/2);
    outer_scapula = new G4EllipticalTube("ScapulaOut", ax_out, by_out, dz/2);


    subtraction = new G4Box("subtraction",ax_out, ax_out, ax_out);

    xx = ax_out * 0.242 ; //(sin 14deg)
    yy  = - ax_out * 0.97; // (cos 14 deg)

    rm = new G4RotationMatrix();
    rm -> rotateZ(-14.* degree);

    scapula_first =  new G4SubtractionSolid("Scapula_first",
                                            outer_scapula,
                                            subtraction,
                                            rm,
                                            G4ThreeVector(xx, yy, 0. *cm));

    xx2 = -ax_out * 0.62470 ; //(cos 51.34deg)
    yy2  = ax_out * 0.78087; // (sin 51.34 deg)

    rm2 = new G4RotationMatrix();
    rm2 -> rotateZ(-38.6598* degree);


    scapula_bone =  new G4SubtractionSolid("Scapula",
                                           scapula_first,
                                           subtraction,
                                           rm2,
                                           G4ThreeVector(xx2, yy2, 0. *cm));

    scapula =  new G4SubtractionSolid("Scapula",
                                      scapula_bone,
                                      inner_scapula);

    G4LogicalVolume* logicLeftScapula = new G4LogicalVolume(scapula,
                                                            skeleton,
                                                            "LeftScapula",
                                                            0, 0, 0);

    G4VPhysicalVolume* physLeftScapula = new G4PVPlacement(0,
                                                           G4ThreeVector(0. * cm, 0. * cm, 24.1 *cm),
                                                           "physicalLeftScapula",
                                                           logicLeftScapula,
                                                           physTrunk,
                                                           false,
                                                           0, true);

/*
    // //////////////////////////////// MaleGenitalia Construction ////////////////////////

    G4double pDz=2.4*cm;
    G4double pTheta=0*degree;
    G4double pPhi=0*degree;
    G4double pDy1=4.76*cm;
    G4double pDx1=9.52*cm;
    G4double pDx2=9.52*cm;
    G4double pAlp1=0*degree;
    G4double pDy2=5*cm;
    G4double pDx3=10*cm;
    G4double pDx4=10*cm;
    G4double pAlp2=0*degree;

    G4Trap* genitaliaTrap= new G4Trap("GenitaliaTrap",
                                      pDz,pTheta,pPhi,pDy1,
                                      pDx1,pDx2,pAlp1,pDy2,
                                      pDx3,pDx4,pAlp2);


    G4double rmin1 = 0.* cm;
    G4double rmin2 = 0.* cm;
    dz= 5 * cm;
    G4double rmax1= 9.51 * cm;
    G4double rmax2= 10.01 * cm;
    startphi= 0.* degree;
    deltaphi= 360. * degree;

    G4Cons* genitaliaLegL = new G4Cons("GenitaliaLegL",
                                       rmin1, rmax1,
                                       rmin2, rmax2, dz/2.,
                                       startphi, deltaphi);

    G4Cons* genitaliaLegR = new G4Cons("GenitaliaLegR",
                                       rmin1, rmax1,
                                       rmin2, rmax2, dz/2.,
                                       startphi, deltaphi);

    G4UnionSolid* genitaliaLegs = new G4UnionSolid("GenitaliaLegs",genitaliaLegL,genitaliaLegR,
                                                   0,
                                                   G4ThreeVector(20.* cm, 0.*cm,0* cm) );



    G4SubtractionSolid* MaleGenitalia = new G4SubtractionSolid("MaleGenitalia",genitaliaTrap,genitaliaLegs,
                                                               0,//
                                                               G4ThreeVector(-10.* cm, -5.*cm,0* cm) );


    G4LogicalVolume* logicMaleGenitalia = new G4LogicalVolume(MaleGenitalia,
                                                              soft,
                                                              "MaleGenitalia",
                                                              0, 0, 0);

    // Define rotation and position here!
    G4VPhysicalVolume* physMaleGenitalia = new G4PVPlacement(0,
                                                             G4ThreeVector(0*cm,5.*cm, -2.4*cm),
                                                             "physicalMaleGenitalia",
                                                             logicMaleGenitalia,
                                                             WorldPhysicalVolume,
                                                             false,
                                                             0, true);

    // //////////////////////////////// RightTeste Construction ////////////////////////

    ax= 1.3*cm;
    by= 1.5*cm;
    cz= 2.3*cm;

    G4Ellipsoid* OneTeste = new G4Ellipsoid("OneTeste", ax, by, cz);

    G4LogicalVolume* logicRightTeste = new G4LogicalVolume(OneTeste, soft, "RightTeste", 0, 0, 0);


    // Define rotation and position here!
    G4VPhysicalVolume* physRightTeste = new G4PVPlacement(0, G4ThreeVector(-1.4*cm,3*cm, 0*cm),"physicalRightTeste", logicRightTeste, physMaleGenitalia, false, 0, true);

    // //////////////////////////////// LeftTeste Construction ////////////////////////


    ax= 1.3*cm;
    by= 1.5*cm;
    cz= 2.3*cm;

    OneTeste = new G4Ellipsoid("OneTeste", ax, by, cz);


    G4LogicalVolume* logicLeftTeste = new G4LogicalVolume(OneTeste, soft, "LeftTeste", 0, 0, 0);

    // Define rotation and position here!
    G4VPhysicalVolume* physLeftTeste = new G4PVPlacement(0, G4ThreeVector(1.4 *cm,3.*cm, 0*cm), "physicalLeftTeste", logicLeftTeste, physMaleGenitalia, false, 0, true);
*/

    // //////////////////////////////// Stomach Construction ////////////////////////


    ax = 4. * cm;
    by = 3. * cm;
    cz = 8. * cm;
    // zcut1 = -8. * cm;
    // zcut2 = 8* cm;

    G4Ellipsoid* stomach_out = new G4Ellipsoid("stomach_out", ax, by, cz);

    // zcut1, zcut2);
    /*
       ax = 3.387 * cm;
       by = 2.387 * cm;
       cz = 7.387 * cm;
       zcut1 = - 7.387 *cm;
       zcut2 = 7.387 *cm;

       G4Ellipsoid* cavity = new G4Ellipsoid ("cavity", ax, by, cz, zcut1, zcut2);

       G4SubtractionSolid* stomach = new G4SubtractionSolid("stomach",stomach_out, cavity);
     */

    G4LogicalVolume* logicStomach = new G4LogicalVolume(stomach_out, soft, "Stomach", 0, 0, 0);

    // Define rotation and position here!
    G4VPhysicalVolume* physStomach = new G4PVPlacement(0,G4ThreeVector(8. *cm,-4. * cm, 0), "physicalStomach", logicStomach, physTrunk, false, 0, true);

    // //////////////////////////////// Spleen Construction ////////////////////////

    ax= 3.5 *cm;
    by= 2. *cm;
    cz= 6. * cm;

    //G4Ellipsoid* spleen = new G4Ellipsoid("spleen", ax, by, cz);


    //G4LogicalVolume* logicSpleen = new G4LogicalVolume(spleen, soft, "Spleen", 0, 0, 0);

    // Define rotation and position here!
    //G4VPhysicalVolume* physSpleen = new G4PVPlacement(0, G4ThreeVector(11. *cm, 3. *cm, 2.*cm), "physicalSpleen", logicSpleen, physTrunk, false, 0, true);


    // //////////////////////////////// Liver Construction ////////////////////////

    dx= 14.19 *cm; //a
    dy= 7.84 *cm;  //b
    dz= 7.21* cm; //(z2-z1)/2

    G4EllipticalTube* firstLiver = new G4EllipticalTube("FirstLiver",dx, dy, dz);

    xx = 20.00 * cm;
    yy = 50.00 * cm;
    zz = 50.00 *cm;

    G4Box* subtrLiver = new G4Box("SubtrLiver", xx/2., yy/2., zz/2.);

    rm_relative = new G4RotationMatrix();
    rm_relative -> rotateY(32.* degree);
    rm_relative -> rotateZ(40.9* degree);

    //  G4double aa = (1.00/31.51);
    //  G4double bb = (1.00/44.75);
    //  G4double cc = (-1.00/38.76);
    //  G4cout<< aa << " "<< bb << " "<<cc<< G4endl;

    //  G4double dd = sqrt(aa*aa + bb*bb + cc*cc);
    //  G4cout<< "dd" << dd << G4endl;
    //  G4cout << aa/dd << "" << bb/dd << " "<< cc/dd <<G4endl;
    //  G4cout << (std::atan(1.42))/deg << G4endl;

    G4SubtractionSolid* liver = new G4SubtractionSolid("Liver",
                                                       firstLiver,subtrLiver,
                                                       rm_relative,
                                                       G4ThreeVector(10.0*cm,0.0*cm,0.0 *cm));

    G4LogicalVolume* logicLiver = new G4LogicalVolume(liver, soft, "Liver", 0, 0, 0);

    // Define rotation and position here!
    rm = new G4RotationMatrix();
    rm->rotateX(180.*degree);
    G4VPhysicalVolume* physLiver = new G4PVPlacement(rm,G4ThreeVector(0. *cm,0. *cm,0.*cm),
                                                     "physicalLiver",
                                                     logicLiver,
                                                     physTrunk,
                                                     false,
                                                     0,true);


    //G4LogicalVolume* LogicAA = new G4LogicalVolume(firstLiver, soft, "Liver1", 0, 0, 0);
    //G4LogicalVolume* LogicBB = new G4LogicalVolume(subtrLiver, soft, "Liver2", 0, 0, 0);
    //G4VPhysicalVolume* p1 = new G4PVPlacement(new G4RotationMatrix(),G4ThreeVector(0. *cm,0. *cm,0.*cm), "AAA", LogicAA, physTrunk, false, 0,true);
    //G4VPhysicalVolume* p2 = new G4PVPlacement(rm_relative,G4ThreeVector(10.0*cm,0.0*cm,0.0 *cm), "BBB", LogicBB, physTrunk, false, 0,true);

    // //////////////////////////////// UpperLargeIntestine Construction ////////////////////////

    dx = 2.5 * cm; // aU
    dy = 2.5* cm; //bU
    dz = 4.775 * cm; //dzU

    G4VSolid* AscendingColonUpperLargeIntestine = new G4EllipticalTube("AscendingColon",dx, dy, dz);

    dx = 2.5 * cm;//bt
    dy = 1.5 *cm;//ct
    dz = 10.5* cm;//x1t

    G4VSolid* TraverseColonUpperLargeIntestine = new G4EllipticalTube("TraverseColon",dx, dy, dz);

    G4RotationMatrix* relative_rm =  new G4RotationMatrix();
    relative_rm -> rotateX(90. * degree);
    relative_rm -> rotateZ(0. * degree);
    relative_rm -> rotateY(90. * degree);
    G4UnionSolid* upperLargeIntestine = new G4UnionSolid("UpperLargeIntestine",
                                                         AscendingColonUpperLargeIntestine,
                                                         TraverseColonUpperLargeIntestine,
                                                         relative_rm,
                                                         G4ThreeVector(8.0 *cm, 0.0,6.275 * cm)); //,0,dzU + ct transverse


    G4LogicalVolume* logicUpperLargeIntestine = new G4LogicalVolume(upperLargeIntestine, soft,
                                                                    "UpperLargeIntestine",
                                                                    0, 0, 0);

    G4VPhysicalVolume* physUpperLargeIntestine = new G4PVPlacement(0,
                                                                   G4ThreeVector(-8.0 * cm, -2.36 *cm,-15.775 *cm),
                                                                   "physicalUpperLargeIntestine",                 //xo, yo, zo ascending colon
                                                                   logicUpperLargeIntestine,
                                                                   physTrunk,
                                                                   false,
                                                                   0, true);

    // //////////////////////////////// SmallIntestine Construction ////////////////////////

    G4double boxX = 11.*cm;
    G4double boxY = 3.53*cm;
    G4double boxZ = 5*cm;


    G4Box* smallIntestineBox = new G4Box("smallIntestineBox",boxX,boxY,boxZ);

    G4double tubsRmin = 0*cm;
    G4double tubsRmax = 11.*cm;
    G4double tubsZ = 5*cm;
    G4double tubsSphi = 0*degree;
    G4double tubsDphi = 360*degree;


    G4Tubs* smallIntestineTubs = new G4Tubs("smallIntestineTubs",tubsRmin,tubsRmax,tubsZ,tubsSphi,tubsDphi);

    //G4IntersectionSolid* SmallIntestine = new G4IntersectionSolid("SmallIntestine",smallIntestineTubs,smallIntestineBox,
    G4IntersectionSolid* filledSmallIntestine1 = new G4IntersectionSolid("filledSmallIntestine1",smallIntestineTubs,smallIntestineBox,
                                                                         0,G4ThreeVector(0*cm,-1.33*cm, 0*cm));

    G4IntersectionSolid* filledSmallIntestine = new G4IntersectionSolid("filledSmallIntestine",filledSmallIntestine1,smallIntestineTubs,
                                                                        0,G4ThreeVector(0*cm,0.8*cm, 0*cm));

    dx = 2.50*cm; // aU
    dy = 2.50*cm; //bU
    dz = 4.775*cm; //dzU

    AscendingColonUpperLargeIntestine = new G4EllipticalTube("AscendingColon",dx, dy, dz);

    dx = 2.50 * cm;//bt
    dy = 1.50 *cm;//ct
    dz = 10.50* cm;//x1t

    TraverseColonUpperLargeIntestine = new G4EllipticalTube("TraverseColon",dx, dy, dz);

    relative_rm =  new G4RotationMatrix();
    relative_rm -> rotateX(90. * degree);
    //relative_rm -> rotateZ(180. * degree);
    relative_rm -> rotateY(90. * degree);
    upperLargeIntestine = new G4UnionSolid("UpperLargeIntestine",
                                           AscendingColonUpperLargeIntestine,
                                           TraverseColonUpperLargeIntestine,
                                           relative_rm,
                                           G4ThreeVector(-8.0 *cm, 0.0*cm,6.275* cm)); //,0,dzU + ct transverse

    dx = 1.88 * cm; //a
    dy = 2.13 *cm; //b
    dz = 7.64 *cm; //(z1-z2)/2

    G4EllipticalTube* DescendingColonLowerLargeIntestine = new G4EllipticalTube("DiscendingColon",dx, dy, dz);

    G4UnionSolid* upperlowerLargeIntestine = new G4UnionSolid("UpperLowerLargeIntestine",
                                                              upperLargeIntestine,
                                                              DescendingColonLowerLargeIntestine,
                                                              0,
                                                              G4ThreeVector(-16.72*cm, 0.0*cm,-2.865* cm)); //,0,dzU + ct t


    G4SubtractionSolid* SmallIntestine = new G4SubtractionSolid("SmallIntestine",
                                                                filledSmallIntestine,
                                                                upperlowerLargeIntestine,
                                                                0,
                                                                G4ThreeVector(8.0*cm,-0.3*cm,-2.775*cm));


    G4LogicalVolume* logicSmallIntestine = new G4LogicalVolume( SmallIntestine,
                                                                soft,
                                                                "SmallIntestine",
                                                                0, 0, 0);
    rm = new G4RotationMatrix();
    rm->rotateX(180.*degree);
    rm->rotateY(180.*degree);
    G4VPhysicalVolume* physSmallIntestine = new G4PVPlacement(rm,
                                                              G4ThreeVector(0*cm, -2.66*cm, -13*cm), // Xcheck the spina position the correct placement shuod be this one
                                                              //G4ThreeVector(0*cm, -5.13*cm, -13*cm), // Xcheck the spina position the correct placement shuod be this one
                                                              //G4ThreeVector(0*cm, -6*cm, -13*cm),
                                                              "physicalSmallIntestine",
                                                              logicSmallIntestine,
                                                              physTrunk,
                                                              false,
                                                              0, true);

    // //////////////////////////////// LowerLargeIntestine Construction ////////////////////////

    dx = 1.88 * cm; //a
    dy = 2.13 *cm; //b
    dz = 7.64 *cm; //(z1-z2)/2

    DescendingColonLowerLargeIntestine = new G4EllipticalTube("DiscendingColon",dx, dy, dz);


    rmin= 0.0 *cm;
    rmax = 1.88 * cm;//a
    G4double rtor= 5.72*cm; //R1
    startphi= 0. * degree;
    deltaphi= 90. * degree;

    G4Torus* SigmoidColonUpLowerLargeIntestine = new G4Torus("SigmoidColonUpLowerLargeIntestine",
                                                             rmin, rmax,rtor,
                                                             startphi, deltaphi);

    rtor = 3. * cm;//R2
    G4VSolid* SigmoidColonDownLowerLargeIntestine = new G4Torus("SigmoidColonDownLowerLargeIntestine",
                                                                rmin, rmax,
                                                                rtor,startphi,deltaphi);

    relative_rm =  new G4RotationMatrix();
    relative_rm -> rotateY(180. * degree);
    relative_rm -> rotateZ(90. * degree);

    G4UnionSolid*  SigmoidColonLowerLargeIntestine = new G4UnionSolid( "SigmoidColonLowerLargeIntestine",
                                                                       SigmoidColonUpLowerLargeIntestine,
                                                                       SigmoidColonDownLowerLargeIntestine,
                                                                       relative_rm,
                                                                       G4ThreeVector(0.0,8.72*cm,0.0));
    // R1 + R2

    G4RotationMatrix* relative_rm_2 =  new G4RotationMatrix();
    relative_rm_2 -> rotateX(90. * degree);

    G4UnionSolid* LowerLargeIntestine = new G4UnionSolid( "LowerLargeIntestine",
                                                          DescendingColonLowerLargeIntestine,
                                                          SigmoidColonLowerLargeIntestine,
                                                          relative_rm_2,
                                                          G4ThreeVector(-5.72*cm,0.0*cm, -7.64*cm)
                                                          ); // -rtor,0, -dz


    G4LogicalVolume* logicLowerLargeIntestine = new G4LogicalVolume( LowerLargeIntestine, soft,
                                                                     "LowerLargeIntestine",
                                                                     0, 0, 0);

    G4VPhysicalVolume* physLowerLargeIntestine = new G4PVPlacement(0,           // R1+ R2, -2.36 (y0), z0
                                                                   G4ThreeVector(8.72*cm, -2.36*cm,-18.64 *cm),
                                                                   "physicalLowerLargeIntestine",
                                                                   logicLowerLargeIntestine,
                                                                   physTrunk,
                                                                   false,
                                                                   0, true);

    // //////////////////////////////// MiddleLowerSpine Construction ////////////////////////

    dx = 2. *cm;
    dy = 2.5 *cm;
    dz = 24. *cm;

    G4VSolid* middleLowerSpine = new G4EllipticalTube("MiddleLowerSpine",dx, dy, dz);

    G4LogicalVolume* logicMiddleLowerSpine = new G4LogicalVolume( middleLowerSpine, skeleton,
                                                                  "MiddleLowerSpine",
                                                                  0, 0, 0);
    // Define rotation and position here!
    G4VPhysicalVolume* physMiddleLowerSpine = new G4PVPlacement(0,G4ThreeVector(0.0 *cm, 5.5 * cm,11. * cm),
                                                                "physicalMiddleLowerSpine",
                                                                logicMiddleLowerSpine,
                                                                physTrunk,
                                                                false,
                                                                0, true);

    // //////////////////////////////// Pancreas Construction ////////////////////////

    ax= 3.*cm; //c
    by= 1.*cm;//b
    cz= 15.*cm;//a
    zcut1= -15. *cm;// -a
    zcut2= 0.0 *cm;

    G4Ellipsoid* pancreasFirst =  new G4Ellipsoid("PancreasFirst",ax, by, cz,
                                                  zcut1, zcut2);

    xx = 6. * cm;// 2*c
    yy = 2. * cm;// 2*b
    zz = 12. * cm; // cz - x1 = 3 cm
    G4Box* subtrPancreas = new G4Box("SubtrPancreas",xx/2., yy/2., zz/2.);

    G4SubtractionSolid* pancreas = new G4SubtractionSolid("pancreas",
                                                          pancreasFirst,
                                                          subtrPancreas,
                                                          0,
                                                          G4ThreeVector(-3 * cm,0.0,-9.*cm));
    //
    G4LogicalVolume* logicPancreas = new G4LogicalVolume(pancreas, soft,
                                                         "Pancreas",
                                                         0, 0, 0);

    rm = new G4RotationMatrix();
    rm->rotateY(90.*degree);
    G4VPhysicalVolume* physPancreas = new G4PVPlacement(rm,
                                                        G4ThreeVector(-0. *cm, 0.0, 2*cm),//x0, 0, 2 cm
                                                        "physicalPancreas",
                                                        logicPancreas,
                                                        physTrunk,
                                                        false,
                                                        0, true);


    // //////////////////////////////// Pelvis Construction ////////////////////////

    dx= 12. *cm; // a2
    dy= 12. * cm; //b2
    dz= 11. * cm; // z2/2

    G4VSolid* outPelvis = new G4EllipticalTube("OutPelvis",dx, dy, dz);

    dx = 11.3 * cm; // a1
    dy = 11.3* cm; // b1
    dz = 12.0 *cm; // z2/2

    G4VSolid* inPelvis = new G4EllipticalTube("InPelvis",dx, dy, dz);

    G4double x = 28. * cm; // a2 * 2
    G4double y = 28. * cm; //b2*2
    z = 24. *cm; // z2

    G4VSolid* subPelvis = new G4Box("SubtrPelvis", x/2., y/2., z/2.);



    G4SubtractionSolid* firstPelvis = new G4SubtractionSolid("FirstPelvis",
                                                             outPelvis,
                                                             inPelvis, 0, G4ThreeVector(0.*cm, -0.8 *cm, 0. * cm));


    G4SubtractionSolid* secondPelvis = new G4SubtractionSolid("SecondPelvis",
                                                              firstPelvis,
                                                              subPelvis, 0,
                                                              G4ThreeVector(0.0,
                                                                            -14. * cm, 0.*cm));
    // half of the y size of the box


    G4SubtractionSolid* pelvis = new G4SubtractionSolid("Pelvis", secondPelvis, subPelvis,
                                                        0,
                                                        G4ThreeVector(0.0,
                                                                      22. * cm, -9. *cm));


    G4LogicalVolume* logicPelvis = new G4LogicalVolume(pelvis, skeleton,
                                                       "Pelvis", 0, 0, 0);


    G4VPhysicalVolume* physPelvis = new G4PVPlacement(0,G4ThreeVector(0.0, -3. * cm,-24. * cm),// 0, y02, z position
                                                      // with respect to the trunk
                                                      "physicalPelvis",
                                                      logicPelvis,
                                                      physTrunk,
                                                      false,
                                                      0, true);

    // //////////////////////////////// RibCage Construction ////////////////////////

    dx= 17. *cm; // a2
    dy= 9.8 * cm; //b2
    G4double thickness= 32.4 * cm; // z2/2 of cage

    G4EllipticalTube* outCage = new G4EllipticalTube("outCage",dx, dy, thickness/2.);

    dx = 16.4 * cm; // a1
    dy = 9.2 * cm; // b1
    dz = 34. *cm; // z2/2

    G4EllipticalTube* inCage = new G4EllipticalTube("inCage",dx, dy, dz/2.);

    G4SubtractionSolid* cage = new G4SubtractionSolid("Cage",
                                                      outCage,
                                                      inCage, 0, G4ThreeVector(0.*cm, 0.*cm, 0. * cm));


    G4LogicalVolume* logicRibCage = new G4LogicalVolume(cage, soft, "Cage", 0, 0, 0);

    G4VPhysicalVolume* physRibCage = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, thickness/2. + 0.1 * cm),
                                                       // with respect to the trunk
                                                       "physicalRibCage",
                                                       logicRibCage,
                                                       physTrunk,
                                                       false,
                                                       0, true);


    xx = 17.*cm;
    yy = 9.8*cm;
    G4double ribThickness = 1.4*cm;
    G4EllipticalTube* rib_out = new G4EllipticalTube("rib_out",xx, yy, ribThickness/2.);

    xx = 16.5 *cm;
    yy = 9.3 * cm;
    zz = 1.5 * cm;
    G4EllipticalTube* rib_in = new G4EllipticalTube("rib_in",xx, yy, zz/2.);
    G4SubtractionSolid* rib = new G4SubtractionSolid("rib",rib_out, rib_in);

    G4LogicalVolume* logicRib= new G4LogicalVolume(rib, skeleton, "Rib", 0, 0, 0);

    G4VPhysicalVolume* physRib1 = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, (- 32.2*cm/2. + 0.8 *cm)),
                                                    // with respect to the trunk
                                                    "physicalRib",
                                                    logicRib,
                                                    physRibCage,
                                                    false,
                                                    0, true);

    G4VPhysicalVolume* physRib2 = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, ( - 32.2*cm/2. + 0.8 *cm + 2.8 *cm)),
                                                    // with respect to the trunk
                                                    "physicalRib",
                                                    logicRib,
                                                    physRibCage,
                                                    false,
                                                    0, true);

    G4VPhysicalVolume* physRib3 = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8 * cm + 5.6 *cm)),
                                                    // with respect to the trunk
                                                    "physicalRib",
                                                    logicRib,
                                                    physRibCage,
                                                    false,
                                                    0, true);

    G4VPhysicalVolume* physRib4 = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8 * cm + 8.4 *cm)),
                                                    // with respect to the trunk
                                                    "physicalRib",
                                                    logicRib,
                                                    physRibCage,
                                                    false,
                                                    0, true);

    G4VPhysicalVolume* physRib5 = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8 * cm + 11.2 *cm)),
                                                    // with respect to the trunk
                                                    "physicalRib",
                                                    logicRib,
                                                    physRibCage,
                                                    false,
                                                    0, true);

    G4VPhysicalVolume* physRib6 = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, (-thickness/2. +  0.8 * cm + 14. *cm)),
                                                    // with respect to the trunk
                                                    "physicalRib",
                                                    logicRib,
                                                    physRibCage,
                                                    false,
                                                    0, true);

    G4VPhysicalVolume* physRib7 = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8 *cm + 16.8 *cm)),
                                                    // with respect to the trunk
                                                    "physicalRib",
                                                    logicRib,
                                                    physRibCage,
                                                    false,
                                                    0, true);

    G4VPhysicalVolume* physRib8 = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8 *cm + 19.6 *cm)),
                                                    // with respect to the trunk
                                                    "physicalRib",
                                                    logicRib,
                                                    physRibCage,
                                                    false,
                                                    0, true);

    G4VPhysicalVolume* physRib9 = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8*cm + 22.4 *cm)),
                                                    // with respect to the trunk
                                                    "physicalRib",
                                                    logicRib,
                                                    physRibCage,
                                                    false,
                                                    0, true);

    G4VPhysicalVolume* physRib10 = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8*cm + 25.2 *cm)),
                                                     // with respect to the trunk
                                                     "physicalRib",
                                                     logicRib,
                                                     physRibCage,
                                                     false,
                                                     0, true);

    G4VPhysicalVolume* physRib11 = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8*cm + 28. *cm)),
                                                     // with respect to the trunk
                                                     "physicalRib",
                                                     logicRib,
                                                     physRibCage,
                                                     false,
                                                     0, true);

    G4VPhysicalVolume* physRib12 = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, (-thickness/2. + 0.8*cm + 30.8 *cm)),
                                                     // with respect to the trunk
                                                     "physicalRib",
                                                     logicRib,
                                                     physRibCage,
                                                     false,
                                                     0, true);

    // //////////////////////////////// RighttLeg Construction ////////////////////////

    G4double rmin1 = 0.* cm;
    G4double rmin2 = 0.* cm;
    dz= 80.0 * cm;
    G4double rmax1= 2.0 * cm;
    G4double rmax2= 10. * cm;
    startphi= 0.* degree;
    deltaphi= 360. * degree;

    G4Cons* leg1 = new G4Cons("Leg1",
                              rmin1, rmax1,
                              rmin2, rmax2, dz/2.,
                              startphi, deltaphi);

    G4LogicalVolume* logicRightLeg = new G4LogicalVolume(leg1,
                                                         soft,
                                                         "RightLeg",
                                                         0, 0, 0);
    rm = new G4RotationMatrix();
    rm->rotateX(180.*degree);
    rm->rotateY(180.*degree);
    G4VPhysicalVolume* physRightLeg = new G4PVPlacement(rm,
                                                        //G4ThreeVector(-10. * cm, 0. * cm, -47. *cm), //FA
                                                        G4ThreeVector(-10. * cm, 0. * cm, -40. *cm),
                                                        "physicalRightLeg",
                                                        logicRightLeg,
                                                        WorldPhysicalVolume,
                                                        false,
                                                        0, true);

    // //////////////////////////////// LeftLeg Construction ////////////////////////

    rmin1 = 0.* cm;
    rmin2 = 0.* cm;
    dz= 80.0 * cm;
    rmax1= 2.0 * cm;
    rmax2= 10. * cm;
    startphi= 0.* degree;
    deltaphi= 360. * degree;

    leg1 = new G4Cons("Leg1",
                      rmin1, rmax1,
                      rmin2, rmax2, dz/2.,
                      startphi, deltaphi);

    G4LogicalVolume* logicLeftLeg = new G4LogicalVolume(leg1,
                                                        soft,
                                                        "LeftLeg",
                                                        0, 0, 0);

    rm = new G4RotationMatrix();
    rm->rotateX(180.*degree);
    rm->rotateY(180.*degree);
    G4VPhysicalVolume* physLeftLeg = new G4PVPlacement(rm,
                                                       //G4ThreeVector(10. * cm, 0. * cm, -47. *cm), //FA
                                                       G4ThreeVector(10. * cm, 0. * cm, -40. *cm),
                                                       "physicalLeftLeg",
                                                       logicLeftLeg,
                                                       WorldPhysicalVolume,
                                                       false,
                                                       0, true);

    // //////////////////////////////// RightLegBone Construction ////////////////////////


    dz = 79.8 * cm;
    rmin1 = 0.0 * cm;
    rmin2 = 0.0 * cm;
    rmax1 = 1. * cm;
    rmax2 = 3.5 * cm;
    startphi = 0. * degree;
    deltaphi = 360. * degree;

    G4Cons* leg_bone = new G4Cons("OneRightLegBone",
                                  rmin1, rmax1,
                                  rmin2, rmax2, dz/2.,
                                  startphi, deltaphi);

    //G4RotationMatrix* rm_relative = new G4RotationMatrix();
    //rm_relative -> rotateY(-12.5 * degree);

    // G4UnionSolid* legs_bones =  new G4UnionSolid("RightLegBone",
    //				       leg_bone, leg_bone,
    //				       0,
    //				       G4ThreeVector(20.* cm, 0.0,0. * cm));


    G4LogicalVolume* logicRightLegBone = new G4LogicalVolume(leg_bone, skeleton,"RightLegBone",
                                                             0, 0, 0);


    // Define rotation and position here!
    G4VPhysicalVolume* physRightLegBone = new G4PVPlacement(0,
                                                            G4ThreeVector(0.0 * cm, 0.0, 0.1*cm),
                                                            "physicalRightLegBone",
                                                            logicRightLegBone,
                                                            physRightLeg,
                                                            false,
                                                            0, true);



    // //////////////////////////////// LeftLegBone Construction ////////////////////////

    dz = 79.8 * cm;
    rmin1 = 0.0 * cm;
    rmin2 = 0.0 * cm;
    rmax1 = 1. * cm;
    rmax2 = 3.5 * cm;
    startphi = 0. * degree;
    deltaphi = 360. * degree;

    leg_bone = new G4Cons("OneLeftLegBone",
                          rmin1, rmax1,
                          rmin2, rmax2, dz/2.,
                          startphi, deltaphi);

    G4LogicalVolume* logicLeftLegBone = new G4LogicalVolume(leg_bone, skeleton,"LeftLegBone",
                                                            0, 0, 0);


    // Define rotation and position here!
    G4VPhysicalVolume* physLeftLegBone = new G4PVPlacement(0,
                                                           G4ThreeVector(0.0 * cm, 0.0, 0.1*cm),
                                                           "physicalLeftLegBone",
                                                           logicLeftLegBone,
                                                           physLeftLeg,
                                                           false,
                                                           0, true);

    return WorldPhysicalVolume;
}

void G4TCPPGeometryFormat::ConstructLogicalVolumes(){

    G4Material* matH2O;
    G4Material* soft;
    G4Material* skeleton;
    G4Material* lung;
    G4Material* adipose;
    G4Material* glandular;
    G4Material* adipose_glandular;

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

    // Water
    d = 1.000*g/cm3;
    matH2O = new G4Material("Water",d,2);
    matH2O->AddElement(elH,2);
    matH2O->AddElement(elO,1);
    matH2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

    /*
    // MIRD soft tissue
    d = 0.9869 *g/cm3;
    soft = new G4Material("soft_tissue",d,16);
    soft->AddElement(elH,0.1047);
    soft->AddElement(elC,0.2302);
    soft->AddElement(elN,0.0234);
    soft->AddElement(elO,0.6321);
    soft->AddElement(elNa,0.0013);
    soft->AddElement(elMg,0.00015);
    soft->AddElement(elP,0.0024);
    soft->AddElement(elS,0.0022);
    soft->AddElement(elCl,0.0014);
    soft->AddElement(elK,0.0021);
    soft->AddElement(elFe,0.000063);
    soft->AddElement(elZn,0.000032);
    soft->AddElement(elRb,0.0000057);
    soft->AddElement(elSr,0.00000034);
    soft->AddElement(elZr,0.000008);
    soft->AddElement(elPb,0.00000016);
    */

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

    // MIRD Skeleton

    d = 1.4862*g/cm3;
    skeleton = new G4Material("skeleton",d,15);
    skeleton -> AddElement(elH,0.0704);
    skeleton -> AddElement(elC,0.2279);
    skeleton -> AddElement(elN,0.0387);
    skeleton -> AddElement(elO,0.4856);
    skeleton -> AddElement(elNa,0.0032);
    skeleton -> AddElement(elMg,0.0011);
    skeleton -> AddElement(elP,0.0694);
    skeleton -> AddElement(elS,0.0017);
    skeleton -> AddElement(elCl,0.0014);
    skeleton -> AddElement(elK,0.0015);
    skeleton -> AddElement(elCa,0.0991);
    skeleton -> AddElement(elFe,0.00008);
    skeleton -> AddElement(elZn,0.000048);
    skeleton -> AddElement(elSr,0.000032);
    skeleton -> AddElement(elPb,0.000011);

    // MIRD lung material
    d = 0.2958 *g/cm3;
    lung = new G4Material("lung_material", d,16);
    lung -> AddElement(elH, 0.1021);
    lung -> AddElement(elC, 0.1001);
    lung -> AddElement(elN,0.028);
    lung -> AddElement(elO,0.7596);
    lung -> AddElement(elNa,0.0019);
    lung -> AddElement(elMg,0.000074);
    lung -> AddElement(elP,0.00081);
    lung -> AddElement(elS,0.0023);
    lung -> AddElement(elCl,0.0027);
    lung -> AddElement(elK,0.0020);
    lung -> AddElement(elCa,0.00007);
    lung -> AddElement(elFe,0.00037);
    lung -> AddElement(elZn,0.000011);
    lung -> AddElement(elRb,0.0000037);
    lung -> AddElement(elSr,0.000000059);
    lung -> AddElement(elPb,0.00000041);

    G4double density_adipose = 0.93 *g/cm3;
    adipose = new G4Material("adipose", density_adipose,8);
    adipose -> AddElement(elH, 0.112);
    adipose -> AddElement(elC, 0.619);
    adipose -> AddElement(elN, 0.017);
    adipose -> AddElement(elO, 0.251);
    adipose -> AddElement(elS, 0.00025);
    adipose -> AddElement(elP, 0.00025);
    adipose -> AddElement(elK, 0.00025);
    adipose -> AddElement(elCa,0.00025);

    G4double density_glandular = 1.04 * g/cm3;
    glandular = new G4Material("glandular", density_glandular,8);
    glandular -> AddElement(elH, 0.1);
    glandular -> AddElement(elC,0.184);
    glandular -> AddElement(elN, 0.032);
    glandular -> AddElement(elO, 0.679);
    glandular -> AddElement(elS, 0.00125);
    glandular -> AddElement(elP, 0.00125);
    glandular -> AddElement(elK, 0.00125);
    glandular -> AddElement(elCa,0.00125);


    d = (density_adipose * 0.5) + (density_glandular * 0.5);
    adipose_glandular = new G4Material("adipose_glandular", d, 2);
    adipose_glandular -> AddMaterial(adipose, 0.5);
    adipose_glandular -> AddMaterial(glandular, 0.5);

    // Air
    d = 1.290*mg/cm3;
    G4Material* matAir = new G4Material("Air",d,2);
    matAir->AddElement(elN,0.7);
    matAir->AddElement(elO,0.3);

    d = 1.04 *g/cm3;
    soft = new G4Material("Brain_soft_tissue",d,18);
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
    soft->AddElement(elZr,0.00001);

    G4double ax = 6.84 * cm;
    G4double by= 8.7 * cm;
    G4double cz = 5.9 * cm;
    G4Ellipsoid* brain = new G4Ellipsoid("Brain", ax, by, cz);
    G4LogicalVolume* logicBrain =  new G4LogicalVolume(brain, soft, "Brain", 0, 0, 0);

/*
    G4double ax = 6.84 * cm;
    G4double by= 8.7 * cm;
    G4double cz = 5.9 * cm;

    ax= 1.242*cm;
    by= 1.542*cm;
    cz= 2.341*cm;

    G4Ellipsoid* OneTeste = new G4Ellipsoid("OneTeste", ax, by, cz);

    G4UnionSolid* Testes = new G4UnionSolid("Teste",OneTeste,OneTeste, 0, G4ThreeVector(-2.8* cm, 0* cm , 0* cm));

    G4LogicalVolume* logicTestes = new G4LogicalVolume(Testes, soft, "Testes", 0, 0, 0);

*/
}
