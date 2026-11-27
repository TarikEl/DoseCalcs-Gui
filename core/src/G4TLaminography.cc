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


#include "G4TLaminography.hh"
#include "G4VPhysicalVolume.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialTable.hh"
#include "Randomize.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Trd.hh"
#include "G4AssemblyVolume.hh"

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
#include "G4NistManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4TNuclearReactorGeometry.hh"
#include "G4TDetectorGeometry.hh"
#include "G4TPassiveProtonBeamLineGeometry.hh"

extern G4ThreeVector PointSourcePositionTowardDetector;
extern G4double ThetaDirectionValueTowardDetector;
extern G4double PhiDirectionValueTowardDetector;
extern G4String UserParametersForCPPCodeImplementation;

G4TLaminography::G4TLaminography(){}
G4TLaminography::~G4TLaminography(){}

G4VPhysicalVolume* G4TLaminography::ConstructPhysicalVolume(){


    /*
     * *************************************   Multi sources and Detectors
    */

    // ==== Materials used ====

    auto nist = G4NistManager::Instance();

    auto Air            = nist->FindOrBuildMaterial("G4_AIR");
    auto Aluminum       = nist->FindOrBuildMaterial("G4_Al");
    auto Pb             = nist->FindOrBuildMaterial("G4_Pb");
    auto NaI            = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    auto Vacuum         = nist->FindOrBuildMaterial("G4_Galactic");

    G4double a = 1.008 * g/mole;
    G4double z = 1.0;
    G4Element* elH = new G4Element("Hydrogen", "H", z, a);
    G4double density = 1e-25 * g/cm3;
    G4double pressure = 3e-18 * pascal;
    G4double temperature = 2.73 * kelvin;
    G4Material* MyVacuum = new G4Material("MyVacuum", density, 1, kStateGas, temperature, pressure);
    MyVacuum->AddElement(elH, 1);

    // ==== WORLD ====
    G4double worldSize = 6.0*m;
    auto worldSolid = new G4Box("World", worldSize/2, worldSize/2, worldSize/2);
    auto worldLogic = new G4LogicalVolume(worldSolid, Air, "WorldLogic");
    WorldPhysicalVolume = new G4PVPlacement(nullptr, G4ThreeVector(), worldLogic, "WorldPhys", nullptr, false, 0, true);


    // ==== PIPE ====
    G4double pipeRadiusOuter = 16.0*cm;
    G4double pipeRadiusInner = 15.0*cm;
    G4double pipeLength      = 280 *cm;
    G4double sourceSideX   = -pipeRadiusInner;  // = –15 cm  (radiation source side)
    G4double detectorSideX = +pipeRadiusInner;  // = +15 cm  (detector side)


    G4String ControlFor = "ColumnCanals";//ColumnCanals, Simple, CokingInOilPipe, ForSpectrumGeneration, DifferentShapes, SedimentAndCorrosionFromSrcAndDetSides

        // 1 Coking in Oil Pipe	Crude Oil (G4_OCTANE)	Carbon (G4_C)	Heavy carbon deposits (coke) inside refinery pipes.
        // 2 Scale in Water PipeWater (G4_WATER)Calcium Carbonate (G4_CALCIUM_CARBONATE)"Limescale" buildup blocking industrial cooling systems.
        // 3 Sand in Gas LineNatural Gas (G4_METHANE)Silicon Dioxide (G4_SILICON_DIOXIDE)Sand accumulation in high-pressure gas pipelines.
        // 4 Corrosion ProductAir (G4_AIR)Iron Oxide (G4_FERRIC_OXIDE)Rust (Hematite) flakes blocking an empty/vented steel pipe.
        // 5 Paraffin Wax PlugCrude Oil (G4_OCTANE)Paraffin (G4_PARAFFIN)Wax precipitation in subsea pipelines due to cold temperatures.

    if(ControlFor == "Simple"){

        auto pipeSolid = new G4Tubs("Pipe", pipeRadiusInner, pipeRadiusOuter, pipeLength/2, 0, 360*deg);
        auto pipeLogic = new G4LogicalVolume(pipeSolid, Aluminum, "PipeLogic");
        new G4PVPlacement(nullptr, G4ThreeVector(), pipeLogic, "PipePhys", worldLogic, false, 0, true);


        // ==== Sphere object ====

        G4double radius2 = 10*cm;  // -70  -50

        auto bigSolid2 = new G4Sphere("BigObject2",
                                      0.0,          // inner radius (solid sphere)
                                      radius2,      // outer radius
                                      0.*deg,       // start phi
                                      360.*deg,     // delta phi (full circle)
                                      0.*deg,       // start theta
                                      180.*deg);    // delta theta (full sphere)

        auto bigLogic2 = new G4LogicalVolume(bigSolid2, Pb, "BigObjectLogic2");

        //Place it centered in pipe
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), bigLogic2,"BigObjectPhys2", worldLogic, false, 0,true);

    }
    else if(ControlFor == "ForSpectrumGeneration"){

        // ============================================================================
        //  GEOMETRY BUILDER FOR CONFIGURATION ID  (1 .. 437500)
        //  1D stack along X:  Wall | Oil | Object | Wall
        //  Object touches the detector-side internal wall.
        //  No return value. No helper functions. Just builds and places.
        //  Requires:  worldPhys  (G4VPhysicalVolume*)  and  configID  (G4int) defined.
        // ============================================================================

        // ---------------------------------------------------------------------------
        // 0. DECODE THE CONFIGURATION ID  (0-based)
        // ---------------------------------------------------------------------------
        G4int configID = 2147;

        G4int id = configID - 1;

        G4int iObjSize  = id % 28;  id /= 28;
        G4int iObjDens  = id % 5;   id /= 5;
        G4int iObjMat   = id % 5;   id /= 5;
        G4int iOilDens  = id % 5;   id /= 5;
        G4int iOilMat   = id % 5;   id /= 5;
        G4int iWallDens = id % 5;   id /= 5;
        G4int iWallMat  = id % 5;

        // ---------------------------------------------------------------------------
        // 1. NIST MANAGER
        // ---------------------------------------------------------------------------
        G4NistManager* nist = G4NistManager::Instance();

        // ---------------------------------------------------------------------------
        // 2. WALL BASE MATERIALS  (index 0 = "Acier inoxydable")
        // ---------------------------------------------------------------------------
        G4Material* wallBase[5];
        wallBase[0] = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
        wallBase[1] = nist->FindOrBuildMaterial("G4_Fe");
        wallBase[2] = nist->FindOrBuildMaterial("G4_Al");
        wallBase[3] = nist->FindOrBuildMaterial("G4_Cu");
        wallBase[4] = nist->FindOrBuildMaterial("G4_Pb");

        // ---------------------------------------------------------------------------
        // 3. OIL BASE MATERIALS
        // ---------------------------------------------------------------------------
        G4Material* oilBase[5];
        oilBase[0] = nist->FindOrBuildMaterial("G4_WATER");
        oilBase[2] = nist->FindOrBuildMaterial("G4_AIR");

        if (!G4Material::GetMaterial("StdOil", false)) {
            G4Material* m = new G4Material("StdOil", 0.85*g/cm3, 2);
            m->AddElement(nist->FindOrBuildElement("C"), 8);
            m->AddElement(nist->FindOrBuildElement("H"), 18);
            oilBase[1] = m;
        } else oilBase[1] = G4Material::GetMaterial("StdOil");

        if (!G4Material::GetMaterial("Brine", false)) {
            G4Material* m = new G4Material("Brine", 1.10*g/cm3, 4);
            m->AddElement(nist->FindOrBuildElement("H"), 2);
            m->AddElement(nist->FindOrBuildElement("O"), 1);
            m->AddElement(nist->FindOrBuildElement("Na"), 1);
            m->AddElement(nist->FindOrBuildElement("Cl"), 1);
            oilBase[3] = m;
        } else oilBase[3] = G4Material::GetMaterial("Brine");

        if (!G4Material::GetMaterial("LightHC", false)) {
            G4Material* m = new G4Material("LightHC", 0.60*g/cm3, 2);
            m->AddElement(nist->FindOrBuildElement("C"), 1);
            m->AddElement(nist->FindOrBuildElement("H"), 4);
            oilBase[4] = m;
        } else oilBase[4] = G4Material::GetMaterial("LightHC");

        // ---------------------------------------------------------------------------
        // 4. OBJECT BASE MATERIALS
        // ---------------------------------------------------------------------------
        G4Material* objBase[5];

        if (!G4Material::GetMaterial("Plaster", false)) {
            G4Material* m = new G4Material("Plaster", 2.32*g/cm3, 4);
            m->AddElement(nist->FindOrBuildElement("Ca"), 1);
            m->AddElement(nist->FindOrBuildElement("S"), 1);
            m->AddElement(nist->FindOrBuildElement("O"), 4);
            m->AddElement(nist->FindOrBuildElement("H"), 4);
            objBase[0] = m;
        } else objBase[0] = G4Material::GetMaterial("Plaster");

        if (!G4Material::GetMaterial("Sediment", false)) {
            G4Material* m = new G4Material("Sediment", 2.65*g/cm3, 2);
            m->AddElement(nist->FindOrBuildElement("Si"), 1);
            m->AddElement(nist->FindOrBuildElement("O"), 2);
            objBase[1] = m;
        } else objBase[1] = G4Material::GetMaterial("Sediment");

        if (!G4Material::GetMaterial("Corrosion", false)) {
            G4Material* m = new G4Material("Corrosion", 5.24*g/cm3, 2);
            m->AddElement(nist->FindOrBuildElement("Fe"), 2);
            m->AddElement(nist->FindOrBuildElement("O"), 3);
            objBase[2] = m;
        } else objBase[2] = G4Material::GetMaterial("Corrosion");

        if (!G4Material::GetMaterial("Wax", false)) {
            G4Material* m = new G4Material("Wax", 0.90*g/cm3, 2);
            m->AddElement(nist->FindOrBuildElement("C"), 25);
            m->AddElement(nist->FindOrBuildElement("H"), 52);
            objBase[3] = m;
        } else objBase[3] = G4Material::GetMaterial("Wax");

        if (!G4Material::GetMaterial("ScaleCaCO3", false)) {
            G4Material* m = new G4Material("ScaleCaCO3", 2.71*g/cm3, 3);
            m->AddElement(nist->FindOrBuildElement("Ca"), 1);
            m->AddElement(nist->FindOrBuildElement("C"), 1);
            m->AddElement(nist->FindOrBuildElement("O"), 3);
            objBase[4] = m;
        } else objBase[4] = G4Material::GetMaterial("ScaleCaCO3");

        // ---------------------------------------------------------------------------
        // 5. PICK BASE MATERIALS AND REFERENCE DENSITIES
        // ---------------------------------------------------------------------------
        G4Material* wallBaseMat = wallBase[iWallMat];
        G4Material* oilBaseMat  = oilBase [iOilMat ];
        G4Material* objBaseMat  = objBase [iObjMat ];

        G4double wallRhoRef = wallBaseMat->GetDensity();
        G4double oilRhoRef  = oilBaseMat ->GetDensity();
        G4double objRhoRef  = objBaseMat ->GetDensity();

        // ---------------------------------------------------------------------------
        // 6. APPLY 5-LEVEL DENSITY FACTORS  (index 2 = reference)
        // ---------------------------------------------------------------------------
        const G4double fac[5] = {0.95, 0.975, 1.00, 1.025, 1.05};
        G4double wallRho = wallRhoRef * fac[iWallDens];
        G4double oilRho  = oilRhoRef  * fac[iOilDens ];
        G4double objRho  = objRhoRef  * fac[iObjDens ];

        // ---------------------------------------------------------------------------
        // 7. CLONE MATERIALS WITH PER-CONFIG DENSITY
        //    Fix: GetElement() returns const G4Element*,
        //         AddElement() wants  G4Element*   -> use const_cast.
        // ---------------------------------------------------------------------------
        G4String sfx = "_cfg" + std::to_string(configID);

        // --- Wall clone ---
        G4Material* wallMat = new G4Material("WallMat"+sfx, wallRho,
                                             wallBaseMat->GetNumberOfElements());
        {
            const G4double* fracs = wallBaseMat->GetFractionVector();
            for (G4int ie = 0; ie < (G4int)wallBaseMat->GetNumberOfElements(); ++ie)
                wallMat->AddElement(
                    const_cast<G4Element*>(wallBaseMat->GetElement(ie)),
                    fracs[ie]);
        }

        // --- Oil clone ---
        G4Material* oilMat = new G4Material("OilMat"+sfx, oilRho,
                                            oilBaseMat->GetNumberOfElements());
        {
            const G4double* fracs = oilBaseMat->GetFractionVector();
            for (G4int ie = 0; ie < (G4int)oilBaseMat->GetNumberOfElements(); ++ie)
                oilMat->AddElement(
                    const_cast<G4Element*>(oilBaseMat->GetElement(ie)),
                    fracs[ie]);
        }

        // --- Object clone ---
        G4Material* objMat = new G4Material("ObjMat"+sfx, objRho,
                                            objBaseMat->GetNumberOfElements());
        {
            const G4double* fracs = objBaseMat->GetFractionVector();
            for (G4int ie = 0; ie < (G4int)objBaseMat->GetNumberOfElements(); ++ie)
                objMat->AddElement(
                    const_cast<G4Element*>(objBaseMat->GetElement(ie)),
                    fracs[ie]);
        }

        // ---------------------------------------------------------------------------
        // 8. FIXED DIMENSIONS
        // ---------------------------------------------------------------------------
        const G4double transverseY = 5.0*cm;
        const G4double transverseZ = 5.0*cm;
        const G4double wallThk     = 3.0*mm;    // 2–4 mm range
        const G4double internalDia = 28.0*cm;   // oil + object = 28 cm

        const G4double objSizeX = (0.5 + (G4double)iObjSize) * cm; // 0.5,1.5,...,27.5
        const G4double oilSizeX = internalDia - objSizeX;

        // ---------------------------------------------------------------------------
        // 9. LOGICAL VOLUMES
        // ---------------------------------------------------------------------------
        G4Box* wallUpSolid = new G4Box("WallUpSolid"+sfx, wallThk/2., transverseY/2., transverseZ/2.);
        G4LogicalVolume* wallUpLV = new G4LogicalVolume(wallUpSolid, wallMat, "WallUpLV"+sfx);

        G4Box* oilSolid = new G4Box("OilSolid"+sfx, oilSizeX/2., transverseY/2., transverseZ/2.);
        G4LogicalVolume* oilLV = new G4LogicalVolume(oilSolid, oilMat, "OilLV"+sfx);

        G4Box* objSolid = new G4Box("ObjSolid"+sfx, objSizeX/2., transverseY/2., transverseZ/2.);
        G4LogicalVolume* objLV = new G4LogicalVolume(objSolid, objMat, "ObjLV"+sfx);

        G4Box* wallDnSolid = new G4Box("WallDnSolid"+sfx, wallThk/2., transverseY/2., transverseZ/2.);
        G4LogicalVolume* wallDnLV = new G4LogicalVolume(wallDnSolid, wallMat, "WallDnLV"+sfx);

        // ---------------------------------------------------------------------------
        // 10. POSITIONS ALONG X
        //     Order -X -> +X :  Wall | Oil | Object | Wall
        //     Object touches the detector-side (+X) inner wall.
        // ---------------------------------------------------------------------------
        const G4double totalX = 2.0*wallThk + oilSizeX + objSizeX;

        G4double xWallUp = -totalX/2. + wallThk/2.;
        G4double xOil    =  xWallUp + wallThk/2. + oilSizeX/2.;
        G4double xObj    =  xOil    + oilSizeX/2. + objSizeX/2.;
        G4double xWallDn =  xObj    + objSizeX/2. + wallThk/2.;

        // ---------------------------------------------------------------------------
        // 11. PLACE VOLUMES
        // ---------------------------------------------------------------------------
        new G4PVPlacement(0, G4ThreeVector(xWallUp, 0, 0), wallUpLV, "WallUpPV"+sfx,
                          worldLogic, false, 0, true);
        new G4PVPlacement(0, G4ThreeVector(xOil,    0, 0), oilLV,    "OilPV"+sfx,
                          worldLogic, false, 0, true);
        new G4PVPlacement(0, G4ThreeVector(xObj,    0, 0), objLV,    "ObjPV"+sfx,
                          worldLogic, false, 0, true);
        new G4PVPlacement(0, G4ThreeVector(xWallDn, 0, 0), wallDnLV, "WallDnPV"+sfx,
                          worldLogic, false, 0, true);

        // ---------------------------------------------------------------------------
        // 12. VISUAL ATTRIBUTES
        // ---------------------------------------------------------------------------
        wallUpLV->SetVisAttributes(G4Colour(0.6, 0.6, 0.6));
        wallDnLV->SetVisAttributes(G4Colour(0.6, 0.6, 0.6));
        oilLV  ->SetVisAttributes(G4Colour(0.2, 0.4, 1.0, 0.4));
        objLV  ->SetVisAttributes(G4Colour(1.0, 0.3, 0.1));






        // ===================================================================================================================================


        // ============================================================================
        // 13. FULL DIAGNOSTIC PRINTOUT OF THE CURRENT CONFIGURATION
        // ============================================================================

        G4cout << "\n";
        G4cout << "=====================================================================\n";
        G4cout << "  CONFIGURATION REPORT  ---  configID = " << configID << "\n";
        G4cout << "=====================================================================\n";

        // ---------------------------------------------------------------------------
        // 13.1  ID DECODING
        // ---------------------------------------------------------------------------
        G4cout << "\n--- 13.1  ID DECODING (0-based indices) ----------------------------\n";
        G4cout << "  configID                = " << configID << "\n";
        G4cout << "  id (0-based)            = " << (configID - 1) << "\n";
        G4cout << "  iWallMat   (0..4)       = " << iWallMat  << "\n";
        G4cout << "  iWallDens  (0..4)       = " << iWallDens << "\n";
        G4cout << "  iOilMat    (0..4)       = " << iOilMat   << "\n";
        G4cout << "  iOilDens   (0..4)       = " << iOilDens  << "\n";
        G4cout << "  iObjMat    (0..4)       = " << iObjMat   << "\n";
        G4cout << "  iObjDens   (0..4)       = " << iObjDens  << "\n";
        G4cout << "  iObjSize   (0..27)      = " << iObjSize  << "\n";

        // ---------------------------------------------------------------------------
        // 13.2  MATERIAL TYPES CHOSEN
        // ---------------------------------------------------------------------------
        G4cout << "\n--- 13.2  MATERIAL TYPES CHOSEN ------------------------------------\n";
        G4cout << "  Wall   base material    = " << wallBaseMat->GetName()
               << "   (index " << iWallMat << ")\n";
        G4cout << "  Oil    base material    = " << oilBaseMat ->GetName()
               << "   (index " << iOilMat  << ")\n";
        G4cout << "  Object base material    = " << objBaseMat ->GetName()
               << "   (index " << iObjMat  << ")\n";

        // ---------------------------------------------------------------------------
        // 13.3  DENSITIES
        // ---------------------------------------------------------------------------
        G4cout << "\n--- 13.3  DENSITIES ------------------------------------------------\n";
        G4cout << "  Density factor table    = {0.95, 0.975, 1.00, 1.025, 1.05}\n";
        G4cout << "  Wall   ref density      = "
               << wallRhoRef/(g/cm3) << " g/cm3   (idx " << iWallDens
               << ", factor " << fac[iWallDens] << ")\n";
        G4cout << "  Wall   actual density   = "
               << wallRho/(g/cm3)    << " g/cm3\n";
        G4cout << "  Oil    ref density      = "
               << oilRhoRef/(g/cm3)  << " g/cm3   (idx " << iOilDens
               << ", factor " << fac[iOilDens]  << ")\n";
        G4cout << "  Oil    actual density   = "
               << oilRho/(g/cm3)     << " g/cm3\n";
        G4cout << "  Object ref density      = "
               << objRhoRef/(g/cm3)  << " g/cm3   (idx " << iObjDens
               << ", factor " << fac[iObjDens]  << ")\n";
        G4cout << "  Object actual density   = "
               << objRho/(g/cm3)     << " g/cm3\n";

        // ---------------------------------------------------------------------------
        // 13.4  COMPOSITIONS (elements + mass fractions)
        // ---------------------------------------------------------------------------
        auto printComposition = [](const G4String& label, G4Material* mat) {
            G4cout << "  " << label << "  (" << mat->GetName()
                   << ", rho = " << mat->GetDensity()/(g/cm3) << " g/cm3)\n";
            const G4double* fracs = mat->GetFractionVector();
            for (G4int ie = 0; ie < (G4int)mat->GetNumberOfElements(); ++ie) {
                G4cout << "      element " << ie
                       << " : " << mat->GetElement(ie)->GetName()
                       << "   mass fraction = " << fracs[ie] << "\n";
            }
        };

        G4cout << "\n--- 13.4  COMPOSITIONS ---------------------------------------------\n";
        printComposition("Wall   material", wallMat);
        printComposition("Oil    material", oilMat );
        printComposition("Object material", objMat );

        // ---------------------------------------------------------------------------
        // 13.5  GEOMETRY DIMENSIONS
        // ---------------------------------------------------------------------------
        G4cout << "\n--- 13.5  GEOMETRY DIMENSIONS --------------------------------------\n";
        G4cout << "  transverseY  (full)     = " << transverseY/cm << " cm\n";
        G4cout << "  transverseZ  (full)     = " << transverseZ/cm << " cm\n";
        G4cout << "  wallThk                 = " << wallThk/mm    << " mm\n";
        G4cout << "  internalDia (oil+obj)   = " << internalDia/cm << " cm\n";
        G4cout << "  objSizeX                = " << objSizeX/cm  << " cm\n";
        G4cout << "  oilSizeX                = " << oilSizeX/cm  << " cm\n";
        G4cout << "  totalX (stack)          = " << totalX/cm    << " cm\n";

        // ---------------------------------------------------------------------------
        // 13.6  VOLUMES AND MASSES
        // ---------------------------------------------------------------------------
        auto printBoxInfo = [](const G4String& label,
                               G4LogicalVolume* lv,
                               G4double xPos) {
            G4Box* box = dynamic_cast<G4Box*>(lv->GetSolid());
            G4double dx = 2.*box->GetXHalfLength();
            G4double dy = 2.*box->GetYHalfLength();
            G4double dz = 2.*box->GetZHalfLength();
            G4double vol = box->GetCubicVolume();
            G4double rho = lv->GetMaterial()->GetDensity();
            G4double mass = vol * rho;

            G4cout << "  " << label << "\n";
            G4cout << "      material            = " << lv->GetMaterial()->GetName()
                   << "   rho = " << rho/(g/cm3) << " g/cm3\n";
            G4cout << "      dimensions (X,Y,Z)  = "
                   << dx/cm << " x " << dy/cm << " x " << dz/cm << " cm\n";
            G4cout << "      volume              = " << vol/cm3 << " cm3\n";
            G4cout << "      mass                = " << mass/g  << " g\n";
            G4cout << "      X position (center) = " << xPos/cm << " cm\n";
        };

        G4cout << "\n--- 13.6  VOLUMES AND MASSES ---------------------------------------\n";
        printBoxInfo("Wall (upstream)  ", wallUpLV, xWallUp);
        printBoxInfo("Oil              ", oilLV,    xOil   );
        printBoxInfo("Object           ", objLV,    xObj   );
        printBoxInfo("Wall (downstream)", wallDnLV, xWallDn);

        // ---------------------------------------------------------------------------
        // 13.7  TOTAL MASS OF THE STACK
        // ---------------------------------------------------------------------------
        G4double massWallUp = wallUpLV->GetSolid()->GetCubicVolume()
                              * wallUpLV->GetMaterial()->GetDensity();
        G4double massWallDn = wallDnLV->GetSolid()->GetCubicVolume()
                              * wallDnLV->GetMaterial()->GetDensity();
        G4double massOil    = oilLV   ->GetSolid()->GetCubicVolume()
                           * oilLV   ->GetMaterial()->GetDensity();
        G4double massObj    = objLV   ->GetSolid()->GetCubicVolume()
                           * objLV   ->GetMaterial()->GetDensity();

        G4cout << "\n--- 13.7  TOTAL MASS OF THE STACK ----------------------------------\n";
        G4cout << "  Wall upstream  mass     = " << massWallUp/g << " g\n";
        G4cout << "  Wall downstream mass    = " << massWallDn/g << " g\n";
        G4cout << "  Oil mass                = " << massOil/g    << " g\n";
        G4cout << "  Object mass             = " << massObj/g    << " g\n";
        G4cout << "  TOTAL mass              = "
               << (massWallUp + massWallDn + massOil + massObj)/g << " g\n";

        // ---------------------------------------------------------------------------
        // 13.8  PHYSICAL VOLUME NAMES AND PLACEMENT SUMMARY
        // ---------------------------------------------------------------------------
        G4cout << "\n--- 13.8  PLACED VOLUMES (physical) --------------------------------\n";
        G4cout << "  Name                 Material            X center [cm]\n";
        G4cout << "  -------------------  ------------------  -------------\n";
        G4cout << "  WallUpPV" << sfx << "  " << wallMat->GetName()
               << "  " << xWallUp/cm << "\n";
        G4cout << "  OilPV"    << sfx << "  " << oilMat ->GetName()
               << "  " << xOil/cm    << "\n";
        G4cout << "  ObjPV"    << sfx << "  " << objMat ->GetName()
               << "  " << xObj/cm    << "\n";
        G4cout << "  WallDnPV" << sfx << "  " << wallMat->GetName()
               << "  " << xWallDn/cm << "\n";

        // ---------------------------------------------------------------------------
        // 13.9  SUMMARY LINE (one-liner for logs / CSV)
        // ---------------------------------------------------------------------------
        G4cout << "\n--- 13.9  SUMMARY (CSV-friendly) -----------------------------------\n";
        G4cout << "configID,"
               << "wallMat,wallRho,wallDensIdx,"
               << "oilMat,oilRho,oilDensIdx,"
               << "objMat,objRho,objDensIdx,"
               << "objSizeX,oilSizeX,totalX,"
               << "massWallUp,massWallDn,massOil,massObj,massTotal\n";

        G4cout << configID                          << ","
               << wallMat->GetName()                << ","
               << wallRho/(g/cm3)                   << ","
               << iWallDens                         << ","
               << oilMat->GetName()                 << ","
               << oilRho/(g/cm3)                    << ","
               << iOilDens                          << ","
               << objMat->GetName()                 << ","
               << objRho/(g/cm3)                    << ","
               << iObjDens                          << ","
               << objSizeX/cm                       << ","
               << oilSizeX/cm                       << ","
               << totalX/cm                         << ","
               << massWallUp/g                      << ","
               << massWallDn/g                      << ","
               << massOil/g                         << ","
               << massObj/g                         << ","
               << (massWallUp + massWallDn + massOil + massObj)/g
               << "\n";

        G4cout << "=====================================================================\n";
        G4cout << "  END OF CONFIGURATION REPORT\n";
        G4cout << "=====================================================================\n\n";
    }
    else if(ControlFor == "CokingInOilPipe" ){

        // ---- Materials ----
        auto Air     = nist->FindOrBuildMaterial("G4_AIR");
        auto Inox    = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL"); // Real refinery material
        auto Oil     = nist->FindOrBuildMaterial("G4_OCTANE");          // Crude oil surrogate
        auto Coke    = nist->FindOrBuildMaterial("G4_C");               // Carbon deposits

        auto pipeSolid = new G4Tubs("PipeShell", pipeRadiusInner, pipeRadiusOuter, pipeLength/2, 0, 360*deg);
        auto pipeLogic = new G4LogicalVolume(pipeSolid, Inox, "PipeLogic");
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), pipeLogic, "PipePhys", worldLogic, false, 0, true);

        // ==== OIL FILLING (The fluid inside the pipe) ====
        // This occupies the entire inner volume of the pipe
        auto oilSolid = new G4Tubs("OilInside", 0, pipeRadiusInner, pipeLength/2, 0, 360*deg);
        auto oilLogic = new G4LogicalVolume(oilSolid, Oil, "OilLogic");
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), oilLogic, "OilPhys", worldLogic, false, 0, true);

        G4int Complicated = 2;

        if(Complicated == 0){
            // ==== COKE OBJECT (The Blockage) ====
            // In a real case, coke builds up as an irregular crust.
            // Here we simulate a 4cm thick "Crescent" deposit at the bottom of the pipe.
            // Position: Shifted -11cm in Y so it sits against the bottom inner wall.

            G4double cokeThickness = 4.0*cm;
            G4double cokeRadius    = 10.0*cm; // A large chunk blocking part of the flow

            auto cokeSolid = new G4Sphere("CokeDeposit",
                                          0.0, cokeRadius,          // Solid sphere
                                          0.*deg, 360.*deg,         // Full circle
                                          0.*deg, 180.*deg);        // Full sphere

            auto cokeLogic = new G4LogicalVolume(cokeSolid, Coke, "CokeLogic");

            // Absolute Position: Placed inside the Oil volume.
            // To simulate a wall deposit, we offset it so it touches the bottom of the inner pipe.
            // Inner radius is 15cm, sphere radius is 10cm -> Offset of 5cm downwards.
            G4double offsetYY = -5.0*cm;

            new G4PVPlacement(nullptr,
                              G4ThreeVector(0, offsetYY, 0),
                              cokeLogic,
                              "CokePhys",
                              worldLogic,    // Mother volume is the Oil
                              false,
                              0,
                              true);
        }
        else if(Complicated == 1 ){

            // =============================================================
            // COMPLICATED COKE GEOMETRY: Irregular Bottom Crust
            // =============================================================
            // =============================================================
            // 3. THE COKE CRUST (Absolute Placement)
            // =============================================================
            // We create an irregular bottom crust (from 180 to 360 degrees)
            auto baseCrustSolid = new G4Tubs("BaseCrust", 11.0*cm, 15.0*cm, 20.0*cm, 180*deg, 180*deg);

            // Add jaggedness by subtracting 3 spheres from the base crust
            G4VSolid* irregularCrustSolid = baseCrustSolid;
            for(int j=0; j<3; j++) {
                G4double sx = (j-1)*5.0*cm; // Spread spheres along X
                auto cutter = new G4Sphere("Cutter", 0, 4.0*cm, 0, 360*deg, 0, 180*deg);
                irregularCrustSolid = new G4SubtractionSolid("IrregularCrust", irregularCrustSolid, cutter,
                                                             nullptr, G4ThreeVector(sx, -13.0*cm, 0));
            }

            auto crustLogic = new G4LogicalVolume(irregularCrustSolid, Coke, "CrustLogic");
            // Absolute Placement in World
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), crustLogic, "CrustPhys", worldLogic, false, 0, true);

            // =============================================================
            // 4. THE OIL FILLING (Boolean Subtraction for Absolute Placement)
            // =============================================================
            // To prevent overlap, the Oil must be (Full Cylinder - Coke Crust)
            auto fullOilCylinder = new G4Tubs("FullOil", 0, pipeRadiusInner, pipeLength/2, 0, 360*deg);
            auto oilWithHoleForCoke = new G4SubtractionSolid("OilWithHole", fullOilCylinder, irregularCrustSolid);

            auto oilLogic = new G4LogicalVolume(oilWithHoleForCoke, Oil, "OilLogic");
            // Absolute Placement in World
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), oilLogic, "OilPhys", worldLogic, false, 0, true);

            // =============================================================
            // 5. COMPLICATED DEBRIS (Absolute Placement)
            // =============================================================
            // Define a jagged chunk
            auto chunkBox = new G4Box("ChunkBox", 1.5*cm, 1.0*cm, 2.0*cm);
            auto chunkSphere = new G4Sphere("ChunkSphere", 0, 1.2*cm, 0, 360*deg, 0, 180*deg);
            auto jaggedChunkSolid = new G4IntersectionSolid("JaggedChunk", chunkBox, chunkSphere);

            auto chunkLogic = new G4LogicalVolume(jaggedChunkSolid, Coke, "ChunkLogic");

            // Array of specific Absolute Positions for the Debris
            std::vector debrisPositions = {
                G4ThreeVector(5*cm, 2*cm, 10*cm),
                G4ThreeVector(-4*cm, -5*cm, -15*cm),
                G4ThreeVector(0*cm, 8*cm, 5*cm),
                G4ThreeVector(-8*cm, 3*cm, -5*cm),
                G4ThreeVector(2*cm, -10*cm, 0*cm)
            };

            for(int i=0; i<5; i++) {
                G4RotationMatrix* absRot = new G4RotationMatrix();
                absRot->rotateX(i*45*deg); // Different rotation for each
                absRot->rotateY(i*30*deg);

                new G4PVPlacement(absRot,
                                  debrisPositions[i],
                                  chunkLogic,
                                  "DebrisPhys",
                                  worldLogic, // Mother is World
                                  false,
                                  i,
                                  true);
            }
        }else if (Complicated == 2){

            // =============================================================
            // 1. ORGANIC COKE CRUST (Irregular Deposition)
            // =============================================================

            G4double innerR = 11.0*cm;
            G4double outerR = 14.5*cm;
            auto baseCrust = new G4Tubs("BaseCrust", innerR, outerR, 20.0*cm, 180*deg, 180*deg);

            G4VSolid* messyCokeSolid = baseCrust;

            for(int i=0; i<12; i++) {
                G4double phi = (180 + (i * 15)) * deg;
                G4double zPos = (i - 6) * 3.0 * cm;
                G4double scale = 1.0 + (i % 3) * 0.5;

                auto bump = new G4Ellipsoid("CokeBump", 2.0*scale*cm, 1.5*scale*cm, 3.0*scale*cm);
                G4ThreeVector bumpPos(innerR * cos(phi), innerR * sin(phi), zPos);

                messyCokeSolid = new G4UnionSolid("CokeWithBump", messyCokeSolid, bump,
                                                  nullptr, bumpPos);

                if(i % 2 == 0) {
                    auto hole = new G4Sphere("CokeHole", 0, 1.5*scale*cm, 0, 360*deg, 0, 180*deg);
                    G4ThreeVector holePos((innerR + 2*cm) * cos(phi + 5*deg),
                                          (innerR + 2*cm) * sin(phi + 5*deg),
                                          zPos - 1*cm);

                    messyCokeSolid = new G4SubtractionSolid("CokeWithHole", messyCokeSolid, hole,
                                                            nullptr, holePos);
                }
            }

            // --- L-Rotation dial 90 daraja 3la Z ---
            G4RotationMatrix* rotZ90 = new G4RotationMatrix();
            rotZ90->rotateZ(90*deg);

            auto cokeLogic = new G4LogicalVolume(messyCokeSolid, Coke, "CokeLogic");
            new G4PVPlacement(rotZ90, G4ThreeVector(0,0,0), cokeLogic, "CokePhys", worldLogic, false, 0, true);

            // =============================================================
            // 2. OIL FILLING
            // =============================================================
            auto fullOilCylinder = new G4Tubs("FullOil", 0, 15.0*cm, 20.0*cm, 0, 360*deg);
            auto oilWithOrganicHole = new G4SubtractionSolid("OilWithOrganicHole", fullOilCylinder, messyCokeSolid);

            auto oilLogic = new G4LogicalVolume(oilWithOrganicHole, Oil, "OilLogic");
            // Rotation 3la Z bach t-aligna m3a l-Coke
            new G4PVPlacement(rotZ90, G4ThreeVector(0,0,0), oilLogic, "OilPhys", worldLogic, false, 0, true);

            // =============================================================
            // 3. FLOATING DEBRIS
            // =============================================================
            // auto part1 = new G4Box("part1", 1.0*cm, 1.0*cm, 1.0*cm);
            // auto part2 = new G4Ellipsoid("part2", 0.8*cm, 1.5*cm, 0.8*cm);
            // auto debrisSolid = new G4UnionSolid("JaggedDebris", part1, part2,
            //                                     new G4RotationMatrix(45*deg, 45*deg, 0),
            //                                     G4ThreeVector(0.5*cm, 0, 0));

            // auto debrisLogic = new G4LogicalVolume(debrisSolid, Coke, "DebrisLogic");

            // for(int k=0; k<8; k++) {
            //     G4RotationMatrix* randRot = new G4RotationMatrix();
            //     randRot->rotateX(k*60*deg);
            //     randRot->rotateZ(90*deg); // Rotation 90 3la Z 7ta hna

            //     G4double r_rand = 5.0 * (k % 3) * cm;
            //     G4double z_rand = (k - 4) * 4 * cm;
            //     G4ThreeVector pos(r_rand * cos(k*45*deg), r_rand * sin(k*45*deg), z_rand);

            //     new G4PVPlacement(randRot, pos, debrisLogic, "DebrisPhys", worldLogic, false, k, true);
            // }
        }
    }
    else if(ControlFor == "DifferentShapes"){

        pipeRadiusOuter = 23.1*cm;
        pipeRadiusInner = 23*cm;
        pipeLength      = 100 *cm;
        sourceSideX     = -pipeRadiusInner;  // = –15 cm  (radiation source side)
        detectorSideX   = +pipeRadiusInner;  // = +15 cm  (detector side)


        // ---- High-Strength Carbon Steel  (ASTM A106 / API 5L grade B) ----------
        // Composition by mass fraction:
        //   Fe 98.4 %  |  C 0.3 %  |  Mn 1.0 %  |  Si 0.3 %
        // Density: 7.85 g/cm³
        G4Element* elFe = nist->FindOrBuildElement("Fe");
        G4Element* elSi = nist->FindOrBuildElement("Si");
        G4Element* elO  = nist->FindOrBuildElement("O");


        // =====================  Materials for Pipe


        // 1
        // Plexiglas / PMMA — $\rho = 1.19\text{ g/cm}^3$
        G4double densityPMMA = 1.19 * g/cm3;
        G4Material* matPlexiglas = new G4Material("Plexiglas_Pipe", densityPMMA, 3);
        matPlexiglas->AddElement(nist->FindOrBuildElement("C"), 5);
        matPlexiglas->AddElement(nist->FindOrBuildElement("H"), 8);
        matPlexiglas->AddElement(nist->FindOrBuildElement("O"), 2);
        // 2
        // Carbon Fiber (Fibre de Carbone) — $\rho = 1.60\text{ g/cm}^3$
        G4double densityCarbon = 1.60 * g/cm3;
        G4Material* matCarbonFiber = new G4Material("CarbonFiber_Pipe", densityCarbon, 1);
        matCarbonFiber->AddElement(nist->FindOrBuildElement("C"), 1);
        // 3
        // Beryllium (G4_Be) — $\rho = 1.85\text{ g/cm}^3$
        G4Material* matBeryllium = nist->FindOrBuildMaterial("G4_Be");
        // 4
        // Aluminium (G4_Al) — $\rho = 2.70\text{ g/cm}^3$
        G4Material* matAluminium = nist->FindOrBuildMaterial("G4_Al");
        // 5
        // Inox 8.0
        auto Inox   = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL"); // Real refinery material



        // =====================  Materials for Objects


        // A ---- Corrosion : iron(III) oxide  Fe₂O₃  (haematite / rust) -------------  8.24 * g/cm3
        G4Material* matCorrosion = new G4Material("IronOxide_Fe2O3", 8.24 * g/cm3, 2);
        matCorrosion->AddElement(elFe, 2);   // 2 Fe atoms per formula unit
        matCorrosion->AddElement(elO,  3);   // 3 O  atoms per formula unit

        // B---- Sediment : silica-based  SiO₂  (compacted silt / mineral scale) ---  2.20 g/cm3
        G4Material* matSediment = new G4Material("SilicaSediment", 2.20 * g/cm3, 2);
        matSediment->AddElement(elSi, 1);
        matSediment->AddElement(elO,  2);

        // C---- Graded sediment   0.30 * 2.20 g/cm3   – low density (border region, 30 % of full density) Models the dilute/porous outer edge of the sediment deposit.
        //G4Material* matSedLow = new G4Material("SedimentBorder", 0.30 * 2.20 * g/cm3, 1);
        G4Material* matSedLow = new G4Material("SedimentBorder", 0.80 * g/cm3, 1);
        matSedLow->AddMaterial(matSediment, 1.0);



        // =====================  Materials for Oil

        auto nistManager = G4NistManager::Instance();

        // 1
        auto Oil     = nist->FindOrBuildMaterial("G4_OCTANE");          // Crude oil surrogate

        // 2
        // =====================================================================
        // OIL STANDARD (Density = 0.80 g/cm3)
        // =====================================================================
        // L-oil standard kiy-t-bna ghaliban men l-Carbon w l-Hydrogen (Hydrocarbons)
        G4Element* elH = nistManager->FindOrBuildElement("H");
        G4Element* elC = nistManager->FindOrBuildElement("C");
        G4double densityOil = 0.80 * g/cm3;
        G4Material* matOilStandard = new G4Material("OilStandard", densityOil, 2);
        matOilStandard->AddElement(elC, 6); // Hexane/Octane approximation
        matOilStandard->AddElement(elH, 14);

        // 3
        // =====================================================================
        // 3. CAS 1: HEAVY BARIUM SCALE (Density = 4.50 g/cm3)
        // =====================================================================
        // Barium Sulfate (BaSO4) huwa l-heavy scale l-ma3rouf f l-refineries (Baryte)
        G4Element* elBa = nistManager->FindOrBuildElement("Ba");
        G4Element* elS  = nistManager->FindOrBuildElement("S");
        //G4Element* elO  = nistManager->FindOrBuildElement("O"); // already defined

        G4double densityScale = 4.50 * g/cm3; // Fixed exact density li bghiti
        G4Material* matHeavyScale = new G4Material("BariumHeavyScale", densityScale, 3);
        matHeavyScale->AddElement(elBa, 1);
        matHeavyScale->AddElement(elS,  1);
        matHeavyScale->AddElement(elO,  4);


        G4Material* matObjects = matCorrosion;  //  matCorrosion 8.24  matHeavyScale 4.50   matSediment 2.20   matSedLow 0.30 x 2.20 = 0.66
        G4Material* matOil = Oil ; //  G4_Galactic 0   Oil 0.7026     matOilStandard 0.8      matHeavyScale 4.50
        G4Material* matPipe = Inox; //  Inox 8  matPlexiglas 1.16


        //std::cout << "Inox->GetDensity()/(g/cm3): " << Inox->GetDensity()/(g/cm3)  <<  std::endl;
        //std::cout << "Oil->GetDensity()/(g/cm3): " << Oil->GetDensity()/(g/cm3)  <<  std::endl;

        // ------------------  Pipe Shell
        auto pipeSolid = new G4Tubs("PipeShell", pipeRadiusInner, pipeRadiusOuter, pipeLength/2, 0, 360*deg);
        auto pipeLogic = new G4LogicalVolume(pipeSolid, matPipe, "PipeLogic");
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), pipeLogic, "PipePhys", worldLogic, false, 0, true);

        // ------------------  Pipe Oil
        auto pipeOilSol = new G4Tubs("PipeOil", 0, pipeRadiusInner, pipeLength/2, 0, 360*deg);
        auto pipeOilLV = new G4LogicalVolume(pipeOilSol, matOil, "PipeOil");
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), pipeOilLV, "PipeOil", worldLogic, false, 0, true);

        // ------------------  Pipe Objects

        G4RotationMatrix* rotSourceSide = new G4RotationMatrix();

        G4RotationMatrix* rotDetectorSide = new G4RotationMatrix();
        rotDetectorSide->rotateZ(180 * deg);

        G4RotationMatrix* rotCylinderX = new G4RotationMatrix();
        rotCylinderX->rotateX(90 * deg);

        G4double boxDim = 4 * cm;
        G4Box* sedimentBox = new G4Box("SedimentBox", boxDim+2*cm, boxDim, boxDim);
        G4LogicalVolume* boxLog = new G4LogicalVolume(sedimentBox, matObjects, "BoxLog");


        auto sedimentTrd = new G4Trd("SedimentTriangle", 4*cm, 0.01*mm, 3*cm, 3*cm, 3*cm);
        G4LogicalVolume* trdLog = new G4LogicalVolume(sedimentTrd, matObjects, "TrdLog");

        G4double semiRMax = 4 * cm;
        G4double semiHeight = 4 * cm;
        G4Tubs* sedimentSemi = new G4Tubs("SedimentSemiCircle", 0.0, semiRMax, semiHeight, 0.0*deg, 360.0*deg);
        G4LogicalVolume* semiLog = new G4LogicalVolume(sedimentSemi, matObjects, "SemiLog");


        G4double sphereRad = 4 * cm;
        G4Sphere* centralSphere = new G4Sphere("CentralSphere", 0.0, sphereRad, 0.0*deg, 360.0*deg, 0.0*deg, 360.0*deg);
        G4LogicalVolume* sphereLog = new G4LogicalVolume(centralSphere, matObjects, "SphereLog");


        G4double posX_sourceSide = (pipeRadiusOuter - 8*cm);
        G4double posX_detectorSide = -(pipeRadiusOuter - 8*cm);

        new G4PVPlacement(0, G4ThreeVector(posX_sourceSide, 0.0, 30.0*cm), boxLog,
                          "Sediment_Box_Source", pipeOilLV, false, 0, true);

        new G4PVPlacement(0, G4ThreeVector(posX_sourceSide, 0.0, 0.0*cm), trdLog,
                          "Sediment_Trd_Source", pipeOilLV, false, 0, true);

        //new G4PVPlacement(0, G4ThreeVector(posX_sourceSide, 0.0, -30.0*cm), semiLog,
        //                  "Sediment_Semi_Source", worldLogic, false, 0, true);
        new G4PVPlacement(rotCylinderX, G4ThreeVector(posX_sourceSide, 0.0, -30.0*cm), semiLog,
                          "Sediment_Semi_Source", pipeOilLV, false, 0, true);

        new G4PVPlacement(rotCylinderX, G4ThreeVector(posX_detectorSide, 0.0, 30.0*cm), boxLog,
                          "Sediment_Box_Detector", pipeOilLV, false, 0, true);


        new G4PVPlacement(rotDetectorSide, G4ThreeVector(posX_detectorSide, 0.0, 0*cm), trdLog,
                          "Sediment_Trd_Detector", pipeOilLV, false, 0, true);


        //new G4PVPlacement(rotDetectorSide, G4ThreeVector(posX_detectorSide, 0.0, -30.0*cm), semiLog,
        //                  "Sediment_Semi_Detector", worldLogic, false, 0, true);
        new G4PVPlacement(rotCylinderX, G4ThreeVector(posX_detectorSide, 0.0, -30.0*cm), semiLog,
                          "Sediment_Semi_Detector", pipeOilLV, false, 0, true);

        new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), sphereLog,
                          "Sediment_Sphere_Center", pipeOilLV, false, 0, true);


    }
    else if(ControlFor == "SedimentAndCorrosionFromSrcAndDetSides"){

        auto Oil     = nist->FindOrBuildMaterial("G4_OCTANE");          // Crude oil surrogate

        // ---- High-Strength Carbon Steel  (ASTM A106 / API 5L grade B) ----------
        // Composition by mass fraction:
        //   Fe 98.4 %  |  C 0.3 %  |  Mn 1.0 %  |  Si 0.3 %
        // Density: 7.85 g/cm³
        G4Element* elFe = nist->FindOrBuildElement("Fe");
        G4Element* elC  = nist->FindOrBuildElement("C");
        G4Element* elMn = nist->FindOrBuildElement("Mn");
        G4Element* elSi = nist->FindOrBuildElement("Si");
        G4Element* elO  = nist->FindOrBuildElement("O");

        G4Material* matSteel = new G4Material("HighStrengthCarbonSteel",
                                              7.85 * g/cm3, 4);
        matSteel->AddElement(elFe, 0.984);   // 98.4 % Fe
        matSteel->AddElement(elC,  0.003);   //  0.3 % C
        matSteel->AddElement(elMn, 0.010);   //  1.0 % Mn
        matSteel->AddElement(elSi, 0.003);   //  0.3 % Si

        // ---- Corrosion : iron(III) oxide  Fe₂O₃  (haematite / rust) -------------
        // Density: 5.24 g/cm³
        G4Material* matCorrosion = new G4Material("IronOxide_Fe2O3",
                                                  5.24 * g/cm3, 2);
        matCorrosion->AddElement(elFe, 2);   // 2 Fe atoms per formula unit
        matCorrosion->AddElement(elO,  3);   // 3 O  atoms per formula unit

        // ---- Sediment : silica-based  SiO₂  (compacted silt / mineral scale) ---
        // Density: 2.20 g/cm³
        G4Material* matSediment = new G4Material("SilicaSediment",
                                                 2.20 * g/cm3, 2);
        matSediment->AddElement(elSi, 1);
        matSediment->AddElement(elO,  2);

        // Graded sediment – low density (border region, 30 % of full density)
        // Models the dilute/porous outer edge of the sediment deposit.
        G4Material* matSedLow = new G4Material("SedimentBorder",
                                               0.30 * 2.20 * g/cm3, 1);
        matSedLow->AddMaterial(matSediment, 1.0);

        // Graded sediment – medium density (intermediate zone, 65 %)
        G4Material* matSedMed = new G4Material("SedimentMiddle",
                                               0.65 * 2.20 * g/cm3, 1);
        matSedMed->AddMaterial(matSediment, 1.0);

        // Full-density core uses matSediment directly (100 %)

        auto Inox    = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL"); // Real refinery material


        // ---- Object zone centres along Z -----------------------------------------
        // Derived from the zone layout described in the header above.
        G4double zObj1 = -110.0 * cm;   // Object 1 zone centre
        G4double zObj2 = -50.0 * cm;   // Object 2 zone centre
        G4double zObj3 = -10.0 * cm;   // Object 3 zone centre
        G4double zObj4 = +30.0 * cm;   // Object 4 zone centre
        G4double zObj5 = +70.0 * cm;   // Object 3 zone centre
        G4double zObj6 = +110.0 * cm;   // Object 4 zone centre











        auto pipeSolid = new G4Tubs("PipeShell", pipeRadiusInner, pipeRadiusOuter, pipeLength/2, 0, 360*deg);
        auto pipeLogic = new G4LogicalVolume(pipeSolid, Inox, "PipeLogic");
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), pipeLogic, "PipePhys", worldLogic, false, 0, true);



        // ==========================================================================
        // SECTION E – BORE FLUID  (water inside the pipe)
        // ==========================================================================
        // The bore is a solid cylinder (rMin = 0) with the same half-length as
        // the pipe. All internal objects are placed as daughters of this volume.

        G4Tubs* solidBore = new G4Tubs("Bore_Solid",
                                       0.0,           // solid cylinder – no hole
                                       pipeRadiusInner,
                                       pipeLength/2,
                                       0.0 * deg,
                                       360.0 * deg);

        G4LogicalVolume* logicBore = new G4LogicalVolume(solidBore,
                                                         Vacuum,//Air,//Oil,
                                                         "Fluid");

        //new G4PVPlacement(nullptr,
        //                  G4ThreeVector(0, 0, 0),
        //                  logicBore,
        //                  "Fluid",
        //                  worldLogic,               // ← placed in world, concentric with pipe
        //                  false, 0, true);



        // ==========================================================================
        // SECTION F – OBJECT 1
        // Corrosion spheres on the SOURCE side (–X inner wall)
        //
        // Physics:
        //   Corrosion pitting is modelled as Fe₂O₃ spheres of varying radii
        //   that grow on the pipe inner wall facing the radiation source (–X side).
        //   Each sphere is tangent to the inner wall:
        //     sphere_centre_X = sourceSideX + sphere_radius
        //                     = –15 cm + r   (pushes the sphere inward so it
        //                       just touches the –X inner wall surface)
        //   Spheres are distributed along Z within the 30-cm zone, with spacing
        //   chosen so no two spheres overlap.
        //   Radii range from 0.5 cm (small pitting) to 3.0 cm (large blister).
        // ==========================================================================

        {
            // List of {radius, Z-offset from zone centre}
            // Z-offsets satisfy: |dz_i - dz_j| > r_i + r_j  (no overlap)
            struct SphSpec { G4double r; G4double dz; };
            std::vector<SphSpec> corrSrc = {
                { 0.5*cm, -13.0*cm },   // tiny pit, zone start
                { 3.2*cm, -10.5*cm },   // small spot
                { 0.9*cm,  -8.0*cm },   // medium pit
                { 2.8*cm,  -5.0*cm },   // larger blister
                { 3.7*cm,  -1.5*cm },   // prominent node
                { 1.2*cm,   2.5*cm },   // largest blister  (biggest radius = 3 cm)
                { 3.7*cm,   6.5*cm },   // decreasing
                { 2.8*cm,   9.5*cm },   // medium
                { 0.7*cm,  12.0*cm },   // small
                { 4.7*cm,  13.5*cm },   // tiny pit, zone end
            };

            for (G4int i = 0; i < (G4int)corrSrc.size(); ++i)
            {
                G4double r  = corrSrc[i].r;
                G4double dz = corrSrc[i].dz;

                // Solid
                std::ostringstream sName; sName << "CorrSrc_S" << i;
                G4Sphere* solid = new G4Sphere(sName.str(),
                                               0.0, r,            // solid sphere
                                               0.0, 360.0*deg,
                                               0.0, 180.0*deg);

                // Logical volume
                std::ostringstream lName; lName << "CorrSrc_L" << i;
                G4LogicalVolume* lv = new G4LogicalVolume(solid, matCorrosion, lName.str());

                // Position:
                //   X = sourceSideX + r  → sphere surface flush with inner wall on –X side
                //   Y = 0
                //   Z = zone centre + offset
                G4double cx = sourceSideX + r;   // –15 + r  (inside bore, touching –X wall)
                G4double cy = 0.0;
                G4double cz = zObj1 + dz;

                std::ostringstream pName; pName << "CorrSrc_P" << i;
                new G4PVPlacement(nullptr,
                                  G4ThreeVector(cx, cy, cz),
                                  lv, pName.str(),
                                  worldLogic,     // mother = bore (water volume)
                                  false, i, true);
            }
        }   // end Object 1


        // ==========================================================================
        // SECTION G – OBJECT 2
        // Corrosion spheres on the DETECTOR side (+X inner wall)
        //
        // Mirror of Object 1 but on the +X side.
        //   sphere_centre_X = detectorSideX – sphere_radius
        //                   = +15 cm – r   (sphere surface flush with +X wall)
        // Zone centre Z = –24 cm
        // ==========================================================================

        {
            struct SphSpec { G4double r; G4double dz; };
            std::vector<SphSpec> corrDet = {
                                            { 0.8*cm, -12.5*cm },
                                            { 1.2*cm,  -9.5*cm },
                                            { 1.8*cm,  -6.5*cm },
                                            { 2.5*cm,  -3.0*cm },
                                            { 3.0*cm,   0.5*cm },   // largest, near zone centre
                                            { 2.8*cm,   4.0*cm },
                                            { 2.0*cm,   7.5*cm },
                                            { 1.5*cm,  10.5*cm },
                                            { 0.8*cm,  13.0*cm },
                                            };

            for (G4int i = 0; i < (G4int)corrDet.size(); ++i)
            {
                G4double r  = corrDet[i].r;
                G4double dz = corrDet[i].dz;

                std::ostringstream sName; sName << "CorrDet_S" << i;
                G4Sphere* solid = new G4Sphere(sName.str(),
                                               0.0, r,
                                               0.0, 360.0*deg,
                                               0.0, 180.0*deg);

                std::ostringstream lName; lName << "CorrDet_L" << i;
                G4LogicalVolume* lv = new G4LogicalVolume(solid, matCorrosion, lName.str());

                // Position:
                //   X = detectorSideX – r  → surface flush with +X inner wall
                G4double cx = detectorSideX - r;   // +15 – r
                G4double cy = 0.0;
                G4double cz = zObj2 + dz;

                std::ostringstream pName; pName << "CorrDet_P" << i;
                new G4PVPlacement(nullptr,
                                  G4ThreeVector(cx, cy, cz),
                                  lv, pName.str(),
                                  worldLogic,
                                  false, i, true);
            }
        }   // end Object 2




        // ==========================================================================
        // SECTION H – OBJECT 3
        // Graded sediment ellipsoid on the SOURCE side (–X inner wall)
        //
        // Physics:
        //   Mineral scale / sediment deposits on the inner pipe wall.
        //   The deposit is pressed against the –X wall (no void between deposit
        //   and wall).  Density is graded:
        //     Border region (outer ellipsoid) : low concentration  (30 % density)
        //     Intermediate  (middle ellipsoid): medium concentration (65 % density)
        //     Core          (inner ellipsoid) : maximum concentration (100 % density)
        //   This is implemented as three NESTED ellipsoids, each placed inside
        //   the previous one as a Geant4 daughter volume.
        //
        // Geometry:
        //   Outer ellipsoid semi-axes: aX = 14 cm, aY = 14 cm, aZ = 14 cm
        //   Centre X = –(pipeInnerR – aX) = –(15 – 14) = –1 cm
        //   → outer surface at cX – aX = –1 – 14 = –15 cm = –pipeInnerR  ✓
        //     (outer ellipsoid surface exactly coincides with inner pipe wall)
        //
        //   Middle ellipsoid semi-axes: aX = 9 cm, aY = 9 cm, aZ = 9 cm
        //   Core   ellipsoid semi-axes: aX = 4.5 cm, aY = 4.5 cm, aZ = 4.5 cm
        //   Both placed at (0,0,0) relative to outer ellipsoid centre (nested).
        //
        // Zone centre Z = +24 cm
        // ==========================================================================

        {
            // Semi-axes of the three nested ellipsoids
            G4double aOut = 14.0 * cm;   // outer shell
            G4double aMid =  9.0 * cm;   // middle shell
            G4double aIn  =  4.5 * cm;   // dense core

            // X-position of the common centre of all three ellipsoids
            // The outer ellipsoid's –X surface must coincide with sourceSideX = –15 cm.
            // Centre X = sourceSideX + aOut = –15 + 14 = –1 cm
            G4double sedCX = sourceSideX + aOut;   // = –1 cm
            G4double sedCY = 0.0;
            G4double sedCZ = zObj3;                // zone centre Z = +24 cm

            // ---- Outer ellipsoid (low-density border) ----------------------------
            // G4Ellipsoid(name, ax, ay, az, zBottomCut, zTopCut)
            // zBottomCut = –aOut, zTopCut = +aOut means no cut (full ellipsoid)
            G4Ellipsoid* solidSedOut = new G4Ellipsoid("SedSrc_Out_S",
                                                       4*cm, 4*cm, aOut,
                                                       -aOut, +aOut);
            G4LogicalVolume* lvSedOut = new G4LogicalVolume(solidSedOut,
                                                            matSedLow,
                                                            "SedSrc_Out_L");

            // Placed inside the bore, with centre at (sedCX, sedCY, sedCZ)
            new G4PVPlacement(nullptr,
                              G4ThreeVector(11*cm, 0, zObj3),
                              lvSedOut, "SedSrc_Out_P",
                              worldLogic, false, 0, true);

        }   // end Object 3


        // ==========================================================================
        // SECTION I – OBJECT 4
        // Graded sediment ellipsoid on the DETECTOR side (+X inner wall)
        //
        // Mirror of Object 3 but on the +X (detector) side.
        //   Centre X = detectorSideX – aOut = +15 – 14 = +1 cm
        //   → outer surface at cX + aOut = +1 + 14 = +15 cm = +pipeInnerR  ✓
        //
        // Zone centre Z = +72 cm
        //
        // A distinct material name (SedimentBorder2 / SedimentMiddle2) is used
        // to avoid Geant4 complaining about duplicate material names.
        // ==========================================================================

        {
            // Density-graded sediment materials for Object 4
            // (separate objects from Object 3 to avoid duplicate material names)
            G4Material* matSedLow2 = new G4Material("SedimentBorder2", 0.30 * 2.20 * g/cm3, 1);
            matSedLow2->AddMaterial(matSediment, 1.0);

            G4Material* matSedMed2 = new G4Material("SedimentMiddle2",
                                                    0.65 * 2.20 * g/cm3, 1);
            matSedMed2->AddMaterial(matSediment, 1.0);

            // Separate visualisation colours (teal palette) to distinguish from Obj 3
            G4VisAttributes* visSedOut2  = new G4VisAttributes(G4Colour(0.30, 0.80, 0.70, 0.30));
            visSedOut2->SetForceSolid(true);
            G4VisAttributes* visSedMid2  = new G4VisAttributes(G4Colour(0.10, 0.65, 0.55, 0.55));
            visSedMid2->SetForceSolid(true);
            G4VisAttributes* visSedCore2 = new G4VisAttributes(G4Colour(0.00, 0.45, 0.35, 0.88));
            visSedCore2->SetForceSolid(true);

            G4double aOut = 14.0 * cm;
            G4double aMid =  9.0 * cm;
            G4double aIn  =  4.5 * cm;

            // Centre X = detectorSideX – aOut = +15 – 14 = +1 cm
            G4double sedCX = detectorSideX - aOut;   // = +1 cm
            G4double sedCY = 0.0;
            G4double sedCZ = zObj4;                  // zone centre Z = +72 cm

            // ---- Outer ellipsoid ------------------------------------------------
            G4Ellipsoid* solidSedOut2 = new G4Ellipsoid("SedDet_Out_S",
                                                        4*cm, 4*cm, aOut,
                                                        -aOut, +aOut);
            G4LogicalVolume* lvSedOut2 = new G4LogicalVolume(solidSedOut2,
                                                             matSedLow2,
                                                             "SedDet_Out_L");
            lvSedOut2->SetVisAttributes(visSedOut2);

            new G4PVPlacement(nullptr,
                              G4ThreeVector(-11*cm, 0, zObj4),
                              lvSedOut2, "SedDet_Out_P",
                              worldLogic, false, 0, true);

        }   // end Object 4


        // ==========================================================================
        // SECTION J – OBJECT 5
        // Sediment wedge on the SOURCE side (−X inner wall)
        //
        // Simulated as a partial G4Tubs (tube sector):
        //   Rmin = 0        → solid wedge, no central hole
        //   Rmax = pipeInnerR = 15 cm  → flush with the bore wall
        //   halfZ = objHalfLen (= half of the 20 cm segment you assign to this zone)
        //   startPhi = 330°, deltaPhi = 60°  → arc from 330° to 30°
        //
        // PHI CONVENTION (Geant4 G4Tubs):
        //   Phi is measured counter-clockwise from the +X axis in the XY plane.
        //   +X direction = detector side  (phi = 0°)
        //   −X direction = source  side   (phi = 180°)
        //   The requested arc 330°→30° straddles the +X axis (phi = 0°/360°),
        //   which places the wedge on the SOURCE side of the pipe
        //   only if your source is at −X and the phi convention wraps as follows:
        //
        //     330° = −30° from the +X axis  →  upper-right quadrant boundary
        //      30°           →  lower-right quadrant boundary
        //   The 60° arc therefore fills the region centred on phi = 0° (+X).
        //
        //   *** IMPORTANT ***
        //   If your source is at −X (phi = 180°), change startPhi to 150° so
        //   the wedge is centred on phi = 180°:
        //       startPhi = 150 * deg, deltaPhi = 60 * deg
        //   The values below implement exactly what was specified (330°→30°).
        //
        // FLOW NOTE:
        //   The wedge spans only 60° of the 360° cross-section, leaving 300°
        //   fully open for fluid flow — the current of water/oil is not blocked.
        //
        // Zone centre Z = zObj5  (define this variable before this block,
        //   e.g.  G4double zObj5 = <your chosen Z centre in cm> * cm;  )
        // ==========================================================================

        {
            // ------------------------------------------------------------------
            // Material: Sediment (SiO2-based, same density as other sediment
            // objects).  A unique name "SedimentWedgeSrc" avoids Geant4 errors
            // from duplicate material names if matSediment is already defined
            // as "G4_SILICON_DIOXIDE" or "Sediment" elsewhere.
            // ------------------------------------------------------------------
            G4Material* matSedWedgeSrc = new G4Material(
                "SedimentWedgeSrc",          // unique material name
                2.20 * g/cm3,                // bulk sediment density (g/cm³)
                //   2.20 g/cm³ ≈ compacted fine sand
                //   adjust if looser (1.6) or denser (2.65 quartz)
                1                            // number of components
                );
            matSedWedgeSrc->AddMaterial(matSediment, 1.0);
            //   matSediment must already be defined in your code (e.g. G4_SILICON_DIOXIDE)

            // ------------------------------------------------------------------
            // Visualisation: orange-brown palette to represent source-side sediment
            //   G4Colour(R, G, B, Alpha)
            //   Alpha = 0 (transparent) → 1 (opaque); 0.75 gives good visibility
            // ------------------------------------------------------------------
            G4VisAttributes* visWedgeSrc = new G4VisAttributes(
                G4Colour(0.80, 0.45, 0.10, 0.75)  // warm orange-brown
                );
            visWedgeSrc->SetForceSolid(true);  // render as solid, not wire-frame

            // ------------------------------------------------------------------
            // Geometry: G4Tubs wedge sector
            //
            //   G4Tubs(name,
            //          Rmin,          // inner radius
            //          Rmax,          // outer radius
            //          halfZ,         // half-length along Z
            //          startPhi,      // start angle (CCW from +X axis)
            //          deltaPhi)      // angular span
            //
            //   halfZ = 10*cm  → 20 cm total length along Z for this segment.
            //   Adjust to match the 20 cm zone length in your layout.
            // ------------------------------------------------------------------
            G4double wedgeSrcHalfZ = 10.0 * cm;   // half of 20 cm object length
                // change to match your zoning

            G4Tubs* solidWedgeSrc = new G4Tubs(
                "SedSrc_Wedge_S",              // solid name
                10 * cm,                      // Rmin = 0  (no central hole)
                pipeRadiusInner,                    // Rmax = 15 cm (flush with bore wall)
                //   pipeInnerR must be in scope;
                //   if not, replace with 15.0*cm
                wedgeSrcHalfZ,                 // halfZ = 10 cm
                330.0 * deg,                   // startPhi = 330°
                //   arc starts 30° before the +X axis
                60.0 * deg                    // deltaPhi = 60°  → ends at 30°
                //   wedge centred on phi = 0° (+X)
                //   *** see PHI note above ***
                );

            G4LogicalVolume* lvWedgeSrc = new G4LogicalVolume(
                solidWedgeSrc,                 // solid
                matSedWedgeSrc,                // material
                "SedSrc_Wedge_L"               // logical volume name
                );
            lvWedgeSrc->SetVisAttributes(visWedgeSrc);

            // ------------------------------------------------------------------
            // Placement inside logicBore
            //   Position: X=0, Y=0, Z=zObj5
            //   The partial phi range already places the wedge material only in
            //   the 330°→30° sector; the remaining 300° remains water (logicBore).
            //   No rotation needed — the phi convention aligns automatically.
            // ------------------------------------------------------------------
            new G4PVPlacement(
                nullptr,                        // no rotation matrix
                G4ThreeVector(0.0, 0.0, zObj5), // centre at (0,0,zObj5)
                //   zObj5 must be defined before this block
                lvWedgeSrc,                     // logical volume to place
                "SedSrc_Wedge_P",               // physical volume name
                worldLogic,                      // mother logical volume (bore)
                false,                          // no boolean operation
                0,                              // copy number
                true                            // check for overlaps (prints warnings)
                );

        }   // end Object 5


        // ==========================================================================
        // SECTION K – OBJECT 6
        // Sediment wedge on the DETECTOR side (+X inner wall)
        //
        // Mirror counterpart of Object 5 but placed on the detector side.
        //   Rmin = 0, Rmax = pipeInnerR = 15 cm
        //   startPhi = 60°, deltaPhi = 60°  → arc from 60° to 120°
        //
        // PHI INTERPRETATION:
        //   60°  = upper-right quadrant (between +X and +Y)
        //   90°  = +Y axis (top of pipe)
        //   120° = upper-left  quadrant (between +Y and −X)
        //   The 60° arc is centred on phi = 90° (+Y, top of pipe).
        //
        //   If your detector is at +X (phi = 0°) and you want the wedge centred
        //   on the detector side, use startPhi = −30° (= 330°) and deltaPhi = 60°
        //   (same as Object 5).  The values below implement the specification
        //   exactly: startPhi = 60°, deltaPhi = 60° (top of the pipe).
        //
        // Zone centre Z = zObj6  (define this variable before this block)
        // ==========================================================================

        {
            // ------------------------------------------------------------------
            // Material: same composition as Object 5 but a distinct Geant4
            // material object (unique name required by Geant4 material store).
            // ------------------------------------------------------------------
            G4Material* matSedWedgeDet = new G4Material(
                "SedimentWedgeDet",          // unique name (must differ from Obj 5)
                2.20 * g/cm3,                // same bulk sediment density
                //   adjust independently if needed
                1
                );
            matSedWedgeDet->AddMaterial(matSediment, 1.0);

            // ------------------------------------------------------------------
            // Visualisation: teal-blue palette to distinguish from Object 5
            //   and to indicate the detector side (conventionally cooler colour)
            // ------------------------------------------------------------------
            G4VisAttributes* visWedgeDet = new G4VisAttributes(
                G4Colour(0.10, 0.60, 0.80, 0.75)  // teal-blue
                );
            visWedgeDet->SetForceSolid(true);

            // ------------------------------------------------------------------
            // Geometry: G4Tubs wedge sector, same radial dimensions as Object 5
            // ------------------------------------------------------------------
            G4double wedgeDetHalfZ = 10.0 * cm;   // half of 20 cm object length

            G4Tubs* solidWedgeDet = new G4Tubs(
                "SedDet_Wedge_S",              // solid name
                10 * cm,                       // Rmin = 0
                pipeRadiusInner,                    // Rmax = 15 cm
                wedgeDetHalfZ,                 // halfZ = 10 cm
                150 * deg,                   // startPhi = 60°
                60 * deg                    // deltaPhi = 60°  → ends at 120°
                //   wedge centred on phi = 90° (+Y, top)
                );

            G4LogicalVolume* lvWedgeDet = new G4LogicalVolume(
                solidWedgeDet,
                matSedWedgeDet,
                "SedDet_Wedge_L"
                );
            lvWedgeDet->SetVisAttributes(visWedgeDet);

            // ------------------------------------------------------------------
            // Placement inside logicBore at (0, 0, zObj6)
            // ------------------------------------------------------------------
            new G4PVPlacement(
                nullptr,
                G4ThreeVector(0.0, 0.0, zObj6),   // zObj6 must be defined before this block
                lvWedgeDet,
                "SedDet_Wedge_P",
                worldLogic,
                false,
                0,
                true
                );

        }   // end Object 6


    }



    else if(ControlFor == "ColumnCanals"){

        // ======================================================================
        //  LABORATORY PILOT COLUMN  –  4 Trays
        //  Tray 2: Pivot at +X (fixed), -X (weir) tilts DOWN into Tray 3.
        //  Tray 1: Broken weir fragment floating in fluid.
        // ======================================================================

        // ----------------------------------------------------------------------
        // 0. ALL PARAMETERS
        // ----------------------------------------------------------------------
        const G4double colOuterR   = 63.5  * cm;
        const G4double colWallT    = 0.5   * cm;
        const G4double colInnerR   = colOuterR - colWallT;
        const G4double colHeight   = 280.0 * cm;
        const G4double colHalfZ    = colHeight / 2.0;

        const G4double trayThk     = 1.5   * cm;
        const G4double weirHeight  = 10.0  * cm;
        const G4double weirThk     = 1.5   * cm;
        const G4double cutAngle    = 90.0  * deg;
        const G4double cutDepth    = colInnerR * (1.0 - std::cos(cutAngle/2.0));
        const G4double chordWidth  = 2.0 * colInnerR * std::sin(cutAngle/2.0);
        const G4double fluidHeight = 10.0 * cm;

        const G4double zTray1 = -colHalfZ + 34.0  * cm;
        const G4double zTray2 = -colHalfZ + 84.0  * cm;
        const G4double zTray3 = -colHalfZ + 154.5 * cm;
        const G4double zTray4 = -colHalfZ + 193.5 * cm;

        const G4double fillFrac1 = 0.0; // Tray 1 is empty in the reference, but we place a fragment anyway
        const G4double fillFrac2 = 0.8659;
        const G4double fillFrac3 = 0.8639;
        const G4double fillFrac4 = 0.8259;

        const G4int weirSide1 = +1;
        const G4int weirSide2 = -1;   // Weir on -X
        const G4int weirSide3 = +1;
        const G4int weirSide4 = -1;

        // Input parameter for Tray 2 tilt angle
        const G4double tiltAngle2 = 26.5 * deg;
        const G4double bottomPoolHeight = 10.0 * cm;
        const G4double bottomPoolZ      = -colHalfZ + bottomPoolHeight/2.0;

        // --- NEW PARAMETERS: Broken Weir Fragment for Tray 1 ---
        // Fragment dimensions (a piece broken from the weir)
        const G4double fragThkX = 1.5 * cm;   // Same as weir thickness
        const G4double fragWidthY = 10.0 * cm; // Width of the broken piece
        const G4double fragHeightZ = 5.0 * cm; // Height of the broken piece
        // Position inside Tray 1 fluid (relative to Tray 1 center)
        const G4double fragPosX = 0.0 * cm;   // Modify as needed
        const G4double fragPosY = 0.0 * cm;   // Modify as needed
        const G4double fragPosZ = 2.5 * cm;   // Modify as needed (sitting in fluid)
        // Rotation of the fragment (e.g., tilted/tumbled)
        const G4double fragRotX = 15.0 * deg;
        const G4double fragRotY = 45.0 * deg;
        const G4double fragRotZ = 0.0 * deg;

        // ----------------------------------------------------------------------
        // 1. MATERIALS
        // ----------------------------------------------------------------------
        G4Element* elFe = nist->FindOrBuildElement("Fe");
        G4Element* elC  = nist->FindOrBuildElement("C");
        G4Element* elMn = nist->FindOrBuildElement("Mn");
        G4Element* elSi = nist->FindOrBuildElement("Si");

        G4Material* matSteel = new G4Material("ColumnCarbonSteel", 7.85 * g/cm3, 4);
        matSteel->AddElement(elFe, 0.984);
        matSteel->AddElement(elC,  0.003);
        matSteel->AddElement(elMn, 0.010);
        matSteel->AddElement(elSi, 0.003);

        auto matWater = nist->FindOrBuildMaterial("G4_WATER");

        // ----------------------------------------------------------------------
        // 2. COLUMN SHELL
        // ----------------------------------------------------------------------
        G4Tubs* colShellSolid = new G4Tubs("ColumnShell", colInnerR, colOuterR, colHalfZ, 0.0*deg, 360.0*deg);
        G4LogicalVolume* colShellLV = new G4LogicalVolume(colShellSolid, matSteel, "ColumnShellLV");
        colShellLV->SetVisAttributes(G4Colour(0.6, 0.6, 0.6));
        colShellLV->SetVisAttributes(G4VisAttributes::GetInvisible());
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), colShellLV, "ColumnShellPV", worldLogic, false, 0, true);

        G4Tubs* colInnerSolid = new G4Tubs("ColumnInner", 0.0, colInnerR, colHalfZ, 0.0*deg, 360.0*deg);
        G4LogicalVolume* colInnerLV = new G4LogicalVolume(colInnerSolid, Air, "ColumnInnerLV");
        colInnerLV->SetVisAttributes(G4VisAttributes::GetInvisible());
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), colInnerLV, "ColumnInnerPV", worldLogic, false, 0, true);

        // ----------------------------------------------------------------------
        // 3. SHARED TRAY SOLID & WEIR SOLID
        // ----------------------------------------------------------------------
        G4Tubs* fullDiskSolid = new G4Tubs("FullDisk", 0.0, colInnerR, trayThk/2.0, 0.0*deg, 360.0*deg);

        G4double cutterHalfX = 10.0 * cm;
        G4double cutterHalfY = colInnerR * 2.0;
        G4double cutterHalfZ = trayThk * 2.0;
        G4Box* cutterBox = new G4Box("CutterBox", cutterHalfX, cutterHalfY, cutterHalfZ);

        G4double cutterPosX = (colInnerR - cutDepth) + cutterHalfX;
        G4ThreeVector cutPos(cutterPosX, 0, 0);

        G4SubtractionSolid* cutDiskSolid = new G4SubtractionSolid("CutDisk", fullDiskSolid, cutterBox, nullptr, cutPos);
        G4LogicalVolume* trayLV = new G4LogicalVolume(cutDiskSolid, matSteel, "TrayLV");
        trayLV->SetVisAttributes(G4Colour(0.3, 0.3, 0.3));

        G4Box* weirSolid = new G4Box("WeirBox", weirThk/2.0, chordWidth/2.0, weirHeight/2.0);
        G4LogicalVolume* weirLV = new G4LogicalVolume(weirSolid, matSteel, "WeirLV");
        weirLV->SetVisAttributes(G4Colour(0.4, 0.4, 0.4));

        // ----------------------------------------------------------------------
        // 4. SHARED FLUID CUTTER SETUP
        // ----------------------------------------------------------------------
        G4double fluidCutterHalfX = 10.0 * cm;
        G4double fluidCutterHalfY = chordWidth / 2.0;
        G4double fluidCutterPosX  = (colInnerR - cutDepth) + fluidCutterHalfX;
        G4ThreeVector fluidCutPos(fluidCutterPosX, 0, 0);

        // ----------------------------------------------------------------------
        // 5. HELPER LAMBDA : place normal trays
        // ----------------------------------------------------------------------
        auto placeTray = [&](const G4String& tag, G4double zPos, G4double fillFrac, G4int weirSide)
        {
            G4RotationMatrix* rot = nullptr;
            if (weirSide < 0) {
                rot = new G4RotationMatrix();
                rot->rotateZ(180.0 * deg);
            }

            new G4PVPlacement(rot, G4ThreeVector(0, 0, zPos), trayLV, "TrayPV_" + tag, colInnerLV, false, 0, true);

            G4double weirCenterXlocal = colInnerR - cutDepth + weirThk/2.0;
            G4double weirZ = zPos + trayThk/2.0 + weirHeight/2.0;
            G4ThreeVector weirPos = (weirSide < 0) ? G4ThreeVector(-weirCenterXlocal, 0, weirZ) : G4ThreeVector(weirCenterXlocal, 0, weirZ);
            new G4PVPlacement(rot, weirPos, weirLV, "WeirPV_" + tag, colInnerLV, false, 0, true);

            if (fillFrac > 0.0) {
                G4double localFluidHeight = fillFrac * fluidHeight;
                G4Tubs* fullFluidSolid = new G4Tubs("FullFluid_" + tag, 0.0, colInnerR, localFluidHeight/2.0, 0.0*deg, 360.0*deg);
                G4Box* localFluidCutter = new G4Box("FluidCutter_" + tag, fluidCutterHalfX, fluidCutterHalfY, localFluidHeight * 2.0);
                G4SubtractionSolid* fluidSolid = new G4SubtractionSolid("FluidSolid_" + tag, fullFluidSolid, localFluidCutter, nullptr, fluidCutPos);
                G4LogicalVolume* fluidLV = new G4LogicalVolume(fluidSolid, matWater, "FluidLV_" + tag);
                fluidLV->SetVisAttributes(G4Colour(0.2, 0.4, 1.0, 0.5));
                G4double fluidZ = zPos + trayThk/2.0 + localFluidHeight/2.0;
                new G4PVPlacement(rot, G4ThreeVector(0, 0, fluidZ), fluidLV, "FluidPV_" + tag, colInnerLV, false, 0, true);
            }
        };

        // ----------------------------------------------------------------------
        // 6. PLACE THE TRAYS
        // ----------------------------------------------------------------------
        // Place Tray 1 (but we will manually handle its fluid to add the fragment)
        {
            G4RotationMatrix* rot1 = nullptr; // weirSide1 = +1
            new G4PVPlacement(rot1, G4ThreeVector(0, 0, zTray1), trayLV, "TrayPV_Tray1", colInnerLV, false, 0, true);

            // Weir for Tray 1
            G4double weirCenterXlocal = colInnerR - cutDepth + weirThk/2.0;
            G4double weirZ = zTray1 + trayThk/2.0 + weirHeight/2.0;
            new G4PVPlacement(rot1, G4ThreeVector(weirCenterXlocal, 0, weirZ), weirLV, "WeirPV_Tray1", colInnerLV, false, 0, true);

            // Fluid for Tray 1 (we force a 50% fill to float the fragment, since fillFrac1=0)
            G4double localFluidHeight = 0.5 * fluidHeight;
            G4Tubs* fullFluidSolid1 = new G4Tubs("FullFluid_Tray1", 0.0, colInnerR, localFluidHeight/2.0, 0.0*deg, 360.0*deg);
            G4Box* localFluidCutter1 = new G4Box("FluidCutter_Tray1", fluidCutterHalfX, fluidCutterHalfY, localFluidHeight * 2.0);
            G4SubtractionSolid* fluidSolid1 = new G4SubtractionSolid("FluidSolid_Tray1", fullFluidSolid1, localFluidCutter1, nullptr, fluidCutPos);
            G4LogicalVolume* fluidLV1 = new G4LogicalVolume(fluidSolid1, matWater, "FluidLV_Tray1");
            fluidLV1->SetVisAttributes(G4Colour(0.2, 0.4, 1.0, 0.5));
            G4double fluidZ1 = zTray1 + trayThk/2.0 + localFluidHeight/2.0;
            new G4PVPlacement(rot1, G4ThreeVector(0, 0, fluidZ1), fluidLV1, "FluidPV_Tray1", colInnerLV, false, 0, true);


            // --- Add Broken Weir Fragment into Tray 1's Fluid ---
            G4Box* fragSolid = new G4Box("WeirFragment", fragThkX/2.0, fragWidthY/2.0, fragHeightZ/2.0);
            G4LogicalVolume* fragLV = new G4LogicalVolume(fragSolid, matSteel, "WeirFragmentLV");
            fragLV->SetVisAttributes(G4Colour(0.8, 0.2, 0.2)); // Red color to distinguish it

            G4RotationMatrix* fragRot = new G4RotationMatrix();
            fragRot->rotateX(fragRotX);
            fragRot->rotateY(fragRotY);
            fragRot->rotateZ(fragRotZ);

            // Position relative to the tray center (zTray1)
            G4double fragZ_global = zTray1 + trayThk/2.0 + localFluidHeight/2.0 + fragPosZ;
            G4ThreeVector fragPos_global(fragPosX, fragPosY, fragZ_global);

            new G4PVPlacement(fragRot, fragPos_global, fragLV, "WeirFragmentPV", colInnerLV, false, 0, true);
        }

        // --- TRAY 2: Tilted Assembly (Pivot on +X, Weir on -X) ---
        {
            G4AssemblyVolume* tray2Assembly = new G4AssemblyVolume();

            G4RotationMatrix* rotTray = new G4RotationMatrix();
            rotTray->rotateZ(180.0 * deg);

            G4ThreeVector posDisk(0, 0, 0);
            tray2Assembly->AddPlacedVolume(trayLV, posDisk, rotTray);

            G4double weirCenterXlocal = colInnerR - cutDepth + weirThk/2.0;
            G4double weirZ = trayThk/2.0 + weirHeight/2.0;
            G4ThreeVector posWeir(-weirCenterXlocal, 0, weirZ);
            tray2Assembly->AddPlacedVolume(weirLV, posWeir, rotTray);

            G4double localFluidHeight = fillFrac2 * fluidHeight;
            G4Tubs* fullFluidSolid = new G4Tubs("FullFluid_Tray2", 0.0, colInnerR, localFluidHeight/2.0, 0.0*deg, 360.0*deg);
            G4Box* localFluidCutter = new G4Box("FluidCutter_Tray2", fluidCutterHalfX, fluidCutterHalfY, localFluidHeight * 2.0);
            G4SubtractionSolid* fluidSolid = new G4SubtractionSolid("FluidSolid_Tray2", fullFluidSolid, localFluidCutter, nullptr, fluidCutPos);
            G4LogicalVolume* fluidLV = new G4LogicalVolume(fluidSolid, matWater, "FluidLV_Tray2");
            fluidLV->SetVisAttributes(G4Colour(0.2, 0.4, 1.0, 0.5));

            G4double fluidZ = trayThk/2.0 + localFluidHeight/2.0;
            G4ThreeVector posFluid(0, 0, fluidZ);
            tray2Assembly->AddPlacedVolume(fluidLV, posFluid, rotTray);

            G4RotationMatrix* rotTilt = new G4RotationMatrix();
            rotTilt->rotateY(-tiltAngle2);

            G4ThreeVector pivotLocal(colInnerR, 0, 0);
            G4ThreeVector pivotGlobal(colInnerR, 0, zTray2);

            G4ThreeVector rotatedPivot = (*rotTilt) * pivotLocal;
            G4ThreeVector translation = pivotGlobal - rotatedPivot;

            tray2Assembly->MakeImprint(colInnerLV, translation, rotTilt, 0, true);
        }

        placeTray("Tray3", zTray3, fillFrac3, weirSide3);
        placeTray("Tray4", zTray4, fillFrac4, weirSide4);

        // ----------------------------------------------------------------------
        // 6.b BOTTOM POOL FLUID
        // ----------------------------------------------------------------------
        {
            G4Tubs* bottomPoolSolid = new G4Tubs("BottomPool", 0.0, colInnerR, bottomPoolHeight/2.0, 0.0*deg, 360.0*deg);
            G4LogicalVolume* bottomPoolLV = new G4LogicalVolume(bottomPoolSolid, matWater, "BottomPoolLV");
            bottomPoolLV->SetVisAttributes(G4Colour(0.2, 0.4, 1.0, 0.5));
            new G4PVPlacement(nullptr, G4ThreeVector(0, 0, bottomPoolZ), bottomPoolLV, "BottomPoolPV", colInnerLV, false, 0, true);
        }

    } // end ColumnCanals


    // ==== Detector ====

    G4double detRadius = 3.0*cm;
    G4double detLength = 3.0*cm;

    auto detSolid = new G4Tubs(
        "DetectorSolid",
        0,
        detRadius,
        detLength/2,
        0,
        360*deg
        );
    auto detLogic = new G4LogicalVolume(
        detSolid,
        NaI,
        "D"
        );

    G4ThreeVector SrcPos = PointSourcePositionTowardDetector;
    G4double SourceDirectionTheta = ThetaDirectionValueTowardDetector;;
    G4double SourceDirectionPhi = PhiDirectionValueTowardDetector; //(180*ThetaDirectionValueTowardDetector)/CLHEP::pi

    int ScanningFor = 1; // 1 for laminography and gamma scanning, 2 for fan beam, and tomography RotationTranslation

    if(ScanningFor == 1){ // for laminography and gamma scanning

        G4double DETECTOR_X_PLANE_POSITION = -PointSourcePositionTowardDetector.getX();

        // ===============================
        // RAY DIRECTION
        // ===============================

        G4ThreeVector rayDir( std::sin(SourceDirectionTheta), 0, std::cos(SourceDirectionTheta) );
        rayDir = rayDir.unit();

        // =================================
        // INTERSECTION WITH DETECTOR PLANE
        // =================================

        G4double t = (DETECTOR_X_PLANE_POSITION - SrcPos.x()) / rayDir.x();

        G4ThreeVector detectorPos( DETECTOR_X_PLANE_POSITION, SrcPos.y(), SrcPos.z() + t * rayDir.z() );

        // ======================================================================
        // Automatic ROTATION toward the source (according to position and theta)
        // ======================================================================

        // detector axis must face source
        G4ThreeVector target = -rayDir;

        // build orthonormal basis
        G4ThreeVector zAxis = target.unit();
        G4ThreeVector up(0,1,0);

        G4ThreeVector xAxis = up.cross(zAxis).unit();
        G4ThreeVector yAxis = zAxis.cross(xAxis).unit();

        // rotation matrix
        G4RotationMatrix* rot = new G4RotationMatrix();

        rot->setRows(xAxis,yAxis,zAxis);

        // ===============================
        // Detector PLACEMENT
        // ===============================

        new G4PVPlacement(
            rot,
            detectorPos,
            detLogic,
            "DetectorPhys",
            worldLogic,
            false,
            0,
            true
            );

    }
    else if (ScanningFor == 2){ // Fan Beam

        // 1. Source Position
        //G4ThreeVector SrcPos = ...;

        // 2. Angle men l-input (matalan 180 degree ila bghitiha m9abla m3a source)
        //G4double detectorPhi = 260.1986;

        // 2. Vector dial l-atjah dial l-ray

        G4ThreeVector rayDir(std::cos(SourceDirectionPhi), std::sin(SourceDirectionPhi), 0);

        // 3. Calcul dial Position l-Detector (Intersection Line-Circle)
        // D-da2ira radius dyalha R. Src hwa l-point A. Dir hwa v.
        // Point P = A + d*v. Khass |P| = R.
        G4double R = std::sqrt(SrcPos.x()*SrcPos.x() + SrcPos.y()*SrcPos.y());
        G4double b = 2.0 * SrcPos.dot(rayDir);
        G4double c = SrcPos.mag2() - R*R;
        G4double delta = b*b - 4.0*c;

        // Kan-akhdou d li positive (fin ghadi l-ray)
        G4double d = (-b + std::sqrt(delta)) / 2.0;
        G4ThreeVector detectorPos = SrcPos + d * rayDir;

        // 4. Orientation: Face-to-Face m3a l-Beam
        // L-detector axe khassu i-koun dima perpendicular 3la l-ray
        G4ThreeVector target = -rayDir; // Pointing back along the ray

        G4ThreeVector zAxisDet = target.unit();
        G4ThreeVector up(0, 0, 1);
        G4ThreeVector xAxisDet = up.cross(zAxisDet).unit();
        G4ThreeVector yAxisDet = zAxisDet.cross(xAxisDet).unit();

        G4RotationMatrix* rot = new G4RotationMatrix();
        rot->setRows(xAxisDet, yAxisDet, zAxisDet);

        // 5. Placement
        new G4PVPlacement(rot, detectorPos, detLogic, "DetectorPhys", worldLogic, false, 0, true);

    }


    return WorldPhysicalVolume;

}
