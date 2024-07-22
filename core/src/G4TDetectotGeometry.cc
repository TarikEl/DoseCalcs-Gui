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



#include "globals.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"


#include "globals.hh"
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4OpticalSurface.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UImanager.hh"
#include "G4VisAttributes.hh"



G4TDetectorGeometry::G4TDetectorGeometry()
{
}
G4TDetectorGeometry::~G4TDetectorGeometry()
{
}



//******************************************** LXe Detector ****************************************

G4VPhysicalVolume* G4TDetectorGeometry::ConstructLXe(){


    fExperimentalHall_box  = nullptr;
    fExperimentalHall_log  = nullptr;
    fExperimentalHall_phys = nullptr;

    fLXe = fAl = fAir = fVacuum = fGlass = nullptr;
    fPstyrene = fPMMA = fPethylene1 = fPethylene2 = nullptr;

    fN = fO = fC = fH = nullptr;

    fSaveThreshold = 0;

    // Resets to default values
    fD_mtl = 0.0635 * cm;

    fScint_x = 17.8 * cm;
    fScint_y = 17.8 * cm;
    fScint_z = 22.6 * cm;

    fNx = 2;
    fNy = 2;
    fNz = 3;

    fOuterRadius_pmt = 2.3 * cm;

    fSphereOn = true;
    fRefl     = 1.0;

    fNfibers      = 15;
    fWLSslab      = false;
    fMainVolumeOn = true;
    //fMainVolume   = nullptr;
    fSlab_z       = 2.5 * mm;

    //G4UImanager::GetUIpointer()->ApplyCommand( "/LXe/detector/scintYieldFactor 1.");

    //if(fLXe_mt)
    //    fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", 12000. / MeV);
    //if(fMPTPStyrene)
    //    fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);


    DefineMaterials();

    //*************************** experimental hall

    // The experimental hall walls are all 1m away from housing walls
    G4double expHall_x = fScint_x + fD_mtl + 1. * m;
    G4double expHall_y = fScint_y + fD_mtl + 1. * m;
    G4double expHall_z = fScint_z + fD_mtl + 1. * m;

    // Create experimental hall
    fExperimentalHall_box = new G4Box("expHall_box", expHall_x, expHall_y, expHall_z);
    fExperimentalHall_log = new G4LogicalVolume(fExperimentalHall_box, fVacuum, "expHall_log", 0, 0, 0);
    fExperimentalHall_phys = new G4PVPlacement( 0, G4ThreeVector(), "expHall", fExperimentalHall_log, 0, false, 0);


    //*************************** housing and scintillator

    G4double housing_x = fScint_x + 2. * fD_mtl;
    G4double housing_y = fScint_y + 2. * fD_mtl;
    G4double housing_z = fScint_z + 2. * fD_mtl;

    fScint_box = new G4Box("scint_box", fScint_x / 2., fScint_y / 2., fScint_z / 2.);
    fHousing_box = new G4Box("housing_box", housing_x / 2., housing_y / 2., housing_z / 2.);

    fScint_log   = new G4LogicalVolume(fScint_box, G4Material::GetMaterial("LXe"), "scint_log", 0, 0, 0);
    fHousing_log = new G4LogicalVolume( fHousing_box, G4Material::GetMaterial("Al"), "housing_log", 0, 0, 0);

    new G4PVPlacement(0, G4ThreeVector(), "housing", fHousing_log, fExperimentalHall_phys, false, 0);
    new G4PVPlacement(0, G4ThreeVector(), fScint_log, "scintillator", fHousing_log, false, 0);

    //*************** Miscellaneous sphere to demonstrate skin surfaces

    fSphere = new G4Sphere("sphere", 0., 2. * cm, 0. * deg, 360. * deg, 0. * deg, 360. * deg);
    fSphere_log = new G4LogicalVolume(fSphere, G4Material::GetMaterial("Al"), "sphere_log");
    //if(fSphereOn)
    new G4PVPlacement(0, G4ThreeVector(5. * cm, 5. * cm, 5. * cm), fSphere_log, "sphere", fScint_log, false, 0);


    //****************** Build PMTs
    G4double innerRadius_pmt   = 0.;
    G4double height_pmt        = fD_mtl / 2.;
    G4double startAngle_pmt    = 0.;
    G4double spanningAngle_pmt = 360. * deg;

    fPmt = new G4Tubs("pmt_tube", innerRadius_pmt, fOuterRadius_pmt, height_pmt, startAngle_pmt, spanningAngle_pmt);

    // the "photocathode" is a metal slab at the back of the glass that
    // is only a very rough approximation of the real thing since it only
    // absorbs or detects the photons based on the efficiency set below
    fPhotocath = new G4Tubs("photocath_tube", innerRadius_pmt, fOuterRadius_pmt, height_pmt / 2., startAngle_pmt, spanningAngle_pmt);

    fPmt_log = new G4LogicalVolume(fPmt, G4Material::GetMaterial("Glass"), "pmt_log");
    fPhotocath_log = new G4LogicalVolume( fPhotocath, G4Material::GetMaterial("Al"), "photocath_log");

    new G4PVPlacement(0, G4ThreeVector(0., 0., -height_pmt / 2.), fPhotocath_log, "photocath", fPmt_log, false, 0);


    //*********** Arrange pmts around the outside of housing**********

    G4double dx = fScint_x / fNx;
    G4double dy = fScint_y / fNy;
    G4double dz = fScint_z / fNz;

    G4double x, y, z;
    G4double xmin = -fScint_x / 2. - dx / 2.;
    G4double ymin = -fScint_y / 2. - dy / 2.;
    G4double zmin = -fScint_z / 2. - dz / 2.;
    G4int k       = 0;

    z = -fScint_z / 2. - height_pmt;  // front
    PlacePMTs(fPmt_log, nullptr, x, y, dx, dy, xmin, ymin, fNx, fNy, x, y, z, k);

    G4RotationMatrix* rm_z = new G4RotationMatrix();
    rm_z->rotateY(180. * deg);
    z = fScint_z / 2. + height_pmt;  // back
    PlacePMTs(fPmt_log, rm_z, x, y, dx, dy, xmin, ymin, fNx, fNy, x, y, z, k);

    G4RotationMatrix* rm_y1 = new G4RotationMatrix();
    rm_y1->rotateY(-90. * deg);
    x = -fScint_x / 2. - height_pmt;  // left
    PlacePMTs(fPmt_log, rm_y1, y, z, dy, dz, ymin, zmin, fNy, fNz, x, y, z, k);

    G4RotationMatrix* rm_y2 = new G4RotationMatrix();
    rm_y2->rotateY(90. * deg);
    x = fScint_x / 2. + height_pmt;  // right
    PlacePMTs(fPmt_log, rm_y2, y, z, dy, dz, ymin, zmin, fNy, fNz, x, y, z, k);

    G4RotationMatrix* rm_x1 = new G4RotationMatrix();
    rm_x1->rotateX(90. * deg);
    y = -fScint_y / 2. - height_pmt;  // bottom
    PlacePMTs(fPmt_log, rm_x1, x, z, dx, dz, xmin, zmin, fNx, fNz, x, y, z, k);

    G4RotationMatrix* rm_x2 = new G4RotationMatrix();
    rm_x2->rotateX(-90. * deg);
    y = fScint_y / 2. + height_pmt;  // top
    PlacePMTs(fPmt_log, rm_x2, x, z, dx, dz, xmin, zmin, fNx, fNz, x, y, z, k);


    SurfaceProperties();



    //*********** LXe WLS Slab **********

    G4VPhysicalVolume* slab = LXeWLSSlab(0, G4ThreeVector(0., 0., -fScint_z / 2. - fSlab_z - 1. * cm),
                                         fExperimentalHall_log, false, 0);


    // Surface properties for the WLS slab
    G4OpticalSurface* scintWrap = new G4OpticalSurface("ScintWrap");

    new G4LogicalBorderSurface("ScintWrap", slab, fExperimentalHall_phys, scintWrap);

    scintWrap->SetType(dielectric_metal);
    scintWrap->SetFinish(polished);
    scintWrap->SetModel(glisur);

    std::vector<G4double> pp           = { 2.0 * eV, 3.5 * eV };
    std::vector<G4double> reflectivity = { 1.0, 1.0 };
    std::vector<G4double> efficiency   = { 0.0, 0.0 };

    G4MaterialPropertiesTable* scintWrapProperty =
            new G4MaterialPropertiesTable();

    scintWrapProperty->AddProperty("REFLECTIVITY", pp, reflectivity);
    scintWrapProperty->AddProperty("EFFICIENCY", pp, efficiency);
    scintWrap->SetMaterialPropertiesTable(scintWrapProperty);

    return fExperimentalHall_phys;
}
void G4TDetectorGeometry::DefineMaterials()
{
    G4double a;  // atomic mass
    G4double z;  // atomic number
    G4double density;

    G4int polyPMMA = 1;
    G4int nC_PMMA  = 3 + 2 * polyPMMA;
    G4int nH_PMMA  = 6 + 2 * polyPMMA;

    G4int polyeth = 1;
    G4int nC_eth  = 2 * polyeth;
    G4int nH_eth  = 4 * polyeth;

    //***Elements
    fH = new G4Element("H", "H", z = 1., a = 1.01 * g / mole);
    fC = new G4Element("C", "C", z = 6., a = 12.01 * g / mole);
    fN = new G4Element("N", "N", z = 7., a = 14.01 * g / mole);
    fO = new G4Element("O", "O", z = 8., a = 16.00 * g / mole);

    //***Materials
    // Liquid Xenon
    fLXe = new G4Material("LXe", z = 54., a = 131.29 * g / mole,
                          density = 3.020 * g / cm3);
    // Aluminum
    fAl = new G4Material("Al", z = 13., a = 26.98 * g / mole,
                         density = 2.7 * g / cm3);
    // Vacuum
    fVacuum = new G4Material("Vacuum", z = 1., a = 1.01 * g / mole,
                             density = universe_mean_density, kStateGas,
                             0.1 * kelvin, 1.e-19 * pascal);
    // Air
    fAir = new G4Material("Air", density = 1.29 * mg / cm3, 2);
    fAir->AddElement(fN, 70 * perCent);
    fAir->AddElement(fO, 30 * perCent);
    // Glass
    fGlass = new G4Material("Glass", density = 1.032 * g / cm3, 2);
    fGlass->AddElement(fC, 91.533 * perCent);
    fGlass->AddElement(fH, 8.467 * perCent);
    // Polystyrene
    fPstyrene = new G4Material("Polystyrene", density = 1.03 * g / cm3, 2);
    fPstyrene->AddElement(fC, 8);
    fPstyrene->AddElement(fH, 8);
    // Fiber(PMMA)
    fPMMA = new G4Material("PMMA", density = 1190. * kg / m3, 3);
    fPMMA->AddElement(fH, nH_PMMA);
    fPMMA->AddElement(fC, nC_PMMA);
    fPMMA->AddElement(fO, 2);
    // Cladding(polyethylene)
    fPethylene1 = new G4Material("Pethylene1", density = 1200. * kg / m3, 2);
    fPethylene1->AddElement(fH, nH_eth);
    fPethylene1->AddElement(fC, nC_eth);
    // Double cladding(flourinated polyethylene)
    fPethylene2 = new G4Material("Pethylene2", density = 1400. * kg / m3, 2);
    fPethylene2->AddElement(fH, nH_eth);
    fPethylene2->AddElement(fC, nC_eth);

    //***Material properties tables

    std::vector<G4double> lxe_Energy = { 7.0 * eV, 7.07 * eV, 7.14 * eV };

    std::vector<G4double> lxe_SCINT = { 0.1, 1.0, 0.1 };
    std::vector<G4double> lxe_RIND  = { 1.59, 1.57, 1.54 };
    std::vector<G4double> lxe_ABSL  = { 35. * cm, 35. * cm, 35. * cm };
    fLXe_mt = new G4MaterialPropertiesTable();
    fLXe_mt->AddProperty("SCINTILLATIONCOMPONENT1", lxe_Energy, lxe_SCINT);
    fLXe_mt->AddProperty("SCINTILLATIONCOMPONENT2", lxe_Energy, lxe_SCINT);
    fLXe_mt->AddProperty("RINDEX", lxe_Energy, lxe_RIND);
    fLXe_mt->AddProperty("ABSLENGTH", lxe_Energy, lxe_ABSL);
    fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", 12000. / MeV);
    fLXe_mt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    fLXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20. * ns);
    fLXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 45. * ns);
    fLXe_mt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    fLXe_mt->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
    fLXe->SetMaterialPropertiesTable(fLXe_mt);

    // Set the Birks Constant for the LXe scintillator
    fLXe->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

    std::vector<G4double> glass_AbsLength = { 420. * cm, 420. * cm, 420. * cm };
    G4MaterialPropertiesTable* glass_mt   = new G4MaterialPropertiesTable();
    glass_mt->AddProperty("ABSLENGTH", lxe_Energy, glass_AbsLength);
    glass_mt->AddProperty("RINDEX", "Fused Silica");
    fGlass->SetMaterialPropertiesTable(glass_mt);

    G4MaterialPropertiesTable* vacuum_mt = new G4MaterialPropertiesTable();
    vacuum_mt->AddProperty("RINDEX", "Air");
    fVacuum->SetMaterialPropertiesTable(vacuum_mt);
    fAir->SetMaterialPropertiesTable(vacuum_mt);  // Give air the same rindex

    std::vector<G4double> wls_Energy = { 2.00 * eV, 2.87 * eV, 2.90 * eV,
                                         3.47 * eV };

    std::vector<G4double> rIndexPstyrene = { 1.5, 1.5, 1.5, 1.5 };
    std::vector<G4double> absorption1    = { 2. * cm, 2. * cm, 2. * cm, 2. * cm };
    std::vector<G4double> scintilFast    = { 0.0, 0.0, 1.0, 1.0 };
    fMPTPStyrene = new G4MaterialPropertiesTable();
    fMPTPStyrene->AddProperty("RINDEX", wls_Energy, rIndexPstyrene);
    fMPTPStyrene->AddProperty("ABSLENGTH", wls_Energy, absorption1);
    fMPTPStyrene->AddProperty("SCINTILLATIONCOMPONENT1", wls_Energy, scintilFast);
    fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
    fMPTPStyrene->AddConstProperty("RESOLUTIONSCALE", 1.0);
    fMPTPStyrene->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 10. * ns);
    fPstyrene->SetMaterialPropertiesTable(fMPTPStyrene);

    // Set the Birks Constant for the Polystyrene scintillator
    fPstyrene->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

    std::vector<G4double> AbsFiber    = { 9.0 * m, 9.0 * m, 0.1 * mm, 0.1 * mm };
    std::vector<G4double> EmissionFib = { 1.0, 1.0, 0.0, 0.0 };
    G4MaterialPropertiesTable* fiberProperty = new G4MaterialPropertiesTable();
    fiberProperty->AddProperty("RINDEX", "PMMA");
    fiberProperty->AddProperty("WLSABSLENGTH", wls_Energy, AbsFiber);
    fiberProperty->AddProperty("WLSCOMPONENT", wls_Energy, EmissionFib);
    fiberProperty->AddConstProperty("WLSTIMECONSTANT", 0.5 * ns);
    fPMMA->SetMaterialPropertiesTable(fiberProperty);

    std::vector<G4double> RefractiveIndexClad1 = { 1.49, 1.49, 1.49, 1.49 };
    G4MaterialPropertiesTable* clad1Property   = new G4MaterialPropertiesTable();
    clad1Property->AddProperty("RINDEX", wls_Energy, RefractiveIndexClad1);
    clad1Property->AddProperty("ABSLENGTH", wls_Energy, AbsFiber);
    fPethylene1->SetMaterialPropertiesTable(clad1Property);

    std::vector<G4double> RefractiveIndexClad2 = { 1.42, 1.42, 1.42, 1.42 };
    G4MaterialPropertiesTable* clad2Property   = new G4MaterialPropertiesTable();
    clad2Property->AddProperty("RINDEX", wls_Energy, RefractiveIndexClad2);
    clad2Property->AddProperty("ABSLENGTH", wls_Energy, AbsFiber);
    fPethylene2->SetMaterialPropertiesTable(clad2Property);
}
void G4TDetectorGeometry::SurfaceProperties()
{
    std::vector<G4double> ephoton = { 7.0 * eV, 7.14 * eV };

    //**Scintillator housing properties
    std::vector<G4double> reflectivity     = { fRefl, fRefl };
    std::vector<G4double> efficiency       = { 0.0, 0.0 };
    G4MaterialPropertiesTable* scintHsngPT = new G4MaterialPropertiesTable();
    scintHsngPT->AddProperty("REFLECTIVITY", ephoton, reflectivity);
    scintHsngPT->AddProperty("EFFICIENCY", ephoton, efficiency);
    G4OpticalSurface* OpScintHousingSurface =
            new G4OpticalSurface("HousingSurface", unified, polished, dielectric_metal);
    OpScintHousingSurface->SetMaterialPropertiesTable(scintHsngPT);

    //**Sphere surface properties
    std::vector<G4double> sphereReflectivity = { 1.0, 1.0 };
    std::vector<G4double> sphereEfficiency   = { 0.0, 0.0 };
    G4MaterialPropertiesTable* spherePT      = new G4MaterialPropertiesTable();
    spherePT->AddProperty("REFLECTIVITY", ephoton, sphereReflectivity);
    spherePT->AddProperty("EFFICIENCY", ephoton, sphereEfficiency);
    G4OpticalSurface* OpSphereSurface =
            new G4OpticalSurface("SphereSurface", unified, polished, dielectric_metal);
    OpSphereSurface->SetMaterialPropertiesTable(spherePT);

    //**Photocathode surface properties
    std::vector<G4double> photocath_EFF     = { 1., 1. };
    std::vector<G4double> photocath_ReR     = { 1.92, 1.92 };
    std::vector<G4double> photocath_ImR     = { 1.69, 1.69 };
    G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
    photocath_mt->AddProperty("EFFICIENCY", ephoton, photocath_EFF);
    photocath_mt->AddProperty("REALRINDEX", ephoton, photocath_ReR);
    photocath_mt->AddProperty("IMAGINARYRINDEX", ephoton, photocath_ImR);
    G4OpticalSurface* photocath_opsurf = new G4OpticalSurface(
                "photocath_opsurf", glisur, polished, dielectric_metal);
    photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

    //**Create logical skin surfaces
    new G4LogicalSkinSurface("photocath_surf", fHousing_log, OpScintHousingSurface);
    new G4LogicalSkinSurface("sphere_surface", fSphere_log, OpSphereSurface);
    new G4LogicalSkinSurface("photocath_surf", fPhotocath_log, photocath_opsurf);
}
void G4TDetectorGeometry::PlacePMTs(G4LogicalVolume* pmt_log, G4RotationMatrix* rot,
                                    G4double& a, G4double& b, G4double da,
                                    G4double db, G4double amin, G4double bmin,
                                    G4int na, G4int nb, G4double& x, G4double& y,
                                    G4double& z, G4int& k)
{
    /*  PlacePMTs : a different way to parameterize placement that does not depend
   * on calculating the position from the copy number
   *
   *  pmt_log = logical volume for pmts to be placed
   *  rot = rotation matrix to apply
   *  a,b = coordinates to vary(ie. if varying in the xy plane then pass x,y)
   *  da,db = value to increment a,b by
   *  amin,bmin = start values for a,b
   *  na,nb = number of repitions in a and b
   *  x,y,z = just pass x,y, and z by reference (the same ones passed for a,b)
   *  k = copy number to start with
   *  sd = sensitive detector for pmts
   */
    a = amin;
    for(G4int j = 1; j <= na; ++j)
    {
        a += da;
        b = bmin;
        for(G4int i = 1; i <= nb; ++i)
        {
            b += db;
            new G4PVPlacement(rot, G4ThreeVector(x, y, z), pmt_log, "pmt",
                              fHousing_log, false, k);
            fPmtPositions.push_back(G4ThreeVector(x, y, z));
            ++k;
        }
    }
}
G4VPhysicalVolume* G4TDetectorGeometry::LXeWLSSlab(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                                                   G4LogicalVolume* pMotherLogical, G4bool pMany,
                                                   G4int pCopyNo)
{

    // commented by me
    //G4PVPlacement(pRot, tlate,
    //              new G4LogicalVolume(new G4Box("temp", 1., 1., 1.),
    //                                  G4Material::GetMaterial("Vacuum"), "temp",
    //                                  0, 0, 0),
    //              "Slab", pMotherLogical, pMany, pCopyNo) ;


    G4double slab_x = fScint_x / 2.;
    G4double slab_y = fScint_y / 2.;

    G4Box* ScintSlab_box = new G4Box("Slab", slab_x, slab_y, fSlab_z);

    fScintSlab_log = new G4LogicalVolume( ScintSlab_box, G4Material::GetMaterial("Polystyrene"), "Slab", 0, 0, 0);

    G4double spacing = 2. * slab_y / fNfibers;

    G4RotationMatrix* rm = new G4RotationMatrix();
    rm->rotateY(90. * deg);


    //*********** Place Fibers **********

    for(G4int i = 0; i < fNfibers; ++i)
    {
        G4double Y = -(spacing) * (fNfibers - 1) * 0.5 + i * spacing;
        LXeWLSFiber(rm, G4ThreeVector(0., Y, 0.), fScintSlab_log, false, 0);
    }

    // by mee
    G4VPhysicalVolume* slab = new G4PVPlacement(pRot, tlate,
                                                fScintSlab_log,
                                                "Slab", pMotherLogical, pMany, pCopyNo) ;
    return slab;
}
G4VPhysicalVolume* G4TDetectorGeometry::LXeWLSFiber(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                                                    G4LogicalVolume* pMotherLogical, G4bool pMany,
                                                    G4int pCopyNo)
{

    // The Fiber
    //
    fFiber_rmin = 0.0 * cm;
    fFiber_rmax = 0.1 * cm;
    fFiber_z    = fScint_x / 2.;
    fFiber_sphi = 0.0 * deg;
    fFiber_ephi = 360. * deg;

    fClad1_rmin = 0.;  // fFiber_rmax;
    fClad1_rmax = fFiber_rmax + 0.015 * fFiber_rmax;

    fClad1_z    = fFiber_z;
    fClad1_sphi = fFiber_sphi;
    fClad1_ephi = fFiber_ephi;

    fClad2_rmin = 0.;  // fClad1_rmax;
    fClad2_rmax = fClad1_rmax + 0.015 * fFiber_rmax;

    fClad2_z    = fFiber_z;
    fClad2_sphi = fFiber_sphi;
    fClad2_ephi = fFiber_ephi;

    // commented by mee
    //new G4PVPlacement(pRot,
    //              tlate,
    //              new G4LogicalVolume(new G4Box("temp", 1., 1., 1.),
    //                                  G4Material::GetMaterial("Vacuum"), "temp",
    //                                  0, 0, 0),
    //              "Cladding2", pMotherLogical, pMany, pCopyNo);


    G4Tubs* fiber_tube = new G4Tubs("Fiber", fFiber_rmin, fFiber_rmax, fFiber_z, fFiber_sphi, fFiber_ephi);

    G4LogicalVolume* fiber_log = new G4LogicalVolume(
                fiber_tube, G4Material::GetMaterial("PMMA"), "Fiber", 0, 0, 0);

    // Cladding (first layer)
    //
    G4Tubs* clad1_tube = new G4Tubs("Cladding1", fClad1_rmin, fClad1_rmax, fClad1_z, fClad1_sphi,
                                    fClad1_ephi);
    G4LogicalVolume* clad1_log = new G4LogicalVolume( clad1_tube, G4Material::GetMaterial("Pethylene1"),
                                                      "Cladding1", 0, 0, 0);

    // Cladding (second layer)
    //
    G4Tubs* clad2_tube = new G4Tubs("Cladding2", fClad2_rmin, fClad2_rmax, fClad2_z, fClad2_sphi,
                                    fClad2_ephi);
    fClad2_log = new G4LogicalVolume( clad2_tube, G4Material::GetMaterial("Pethylene2"),
                                      "Cladding2", 0, 0, 0);

    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), fiber_log, "Fiber", clad1_log,
                      false, 0);
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), clad1_log, "Cladding1",
                      fClad2_log, false, 0);



    // by me
    new G4PVPlacement(pRot,
                      tlate,
                      fClad2_log,
                      "Cladding2", pMotherLogical, pMany, pCopyNo);

}







//******************************************** Neutron Detector ****************************************
// from site https://github.com/joshhooker/icSimulation/tree/master

G4VPhysicalVolume* G4TDetectorGeometry::ConstructNeutronDetector() {

    // Overlaps flag
    G4bool checkOverlaps = false;

    G4double z, fractionmass;
    G4int nel, natoms;

    inchtocm = 2.54;
    fGasType =  "P10";
    fTemperature = 293.15;
    fPressureInTorr = 40;
    fNumGrids = 4;
    fGridSize = 20*cm;
    fGridRadius = 70*cm;
    fGridDist = 70*cm;
    fScintDist = 140*cm;
    fWireThickness = 0.02*cm;
    fUseInches = false;
    //"scintResolution"     : 0.1;
    //"extraGridResolution" : 0.1;
    //"interactive"         : false;

    ConstructMaterials();

    // Define Elements
    G4Element* H  = new G4Element("Hydrogen", "H",  z = 1.,  1.008*g/mole);
    G4Element* C  = new G4Element("Carbon",   "C",  z = 6.,  12.011*g/mole);
    G4Element* O  = new G4Element("Oxygen",   "O",  z = 8.,  15.999*g/mole);
    G4Element* F  = new G4Element("Fluorine", "F",  z = 9.,  18.998*g/mole);
    G4Element* Ar = new G4Element("Argon",    "Ar", z = 18., 39.948*g/mole);
    G4Element* Cr = new G4Element("Chrome",   "Cr", z = 25., 51.996*g/mole);
    G4Element* Fe = new G4Element("Iron",     "Fe", z = 26., 55.845*g/mole);
    G4Element* Co = new G4Element("Cobalt",   "Co", z = 27., 58.933*g/mole);
    G4Element* Ni = new G4Element("Nickel",   "Ni", z = 28., 58.693*g/mole);
    G4Element* W  = new G4Element("Tungsten", "W",  z = 74., 183.850*g/mole);
    G4Element* Au = new G4Element("Gold",     "Au", z = 79., 196.967*g/mole);

    // Define Havar
    G4Material* Havar = new G4Material("Havar", 8.3*g/cm3, nel = 5);
    Havar->AddElement(Cr, fractionmass = 0.1785);
    Havar->AddElement(Fe, fractionmass = 0.1822);
    Havar->AddElement(Co, fractionmass = 0.4452);
    Havar->AddElement(Ni, fractionmass = 0.1310);
    Havar->AddElement(W,  fractionmass = 0.0631);

    // Define Mylar
    G4Material* Mylar = new G4Material("Mylar", 1.397*g/cm3, nel = 3);
    Mylar->AddElement(H, natoms = 8);
    Mylar->AddElement(C, natoms = 10);
    Mylar->AddElement(O, natoms = 4);

    // Define BC-400 Plastic Scintillator
    G4Material* BC400 = new G4Material("BC400", 1.032*g/cm3, nel = 2);
    BC400->AddElement(H, natoms = 10);
    BC400->AddElement(C, natoms = 9);

    // Define C2D4 (Deurated Polyethylene)
    G4Isotope* iso_H2 = new G4Isotope("H2", 1, 2, 2.014*g/mole);
    G4Element* D = new G4Element("Deuterium Atom", "D", 1);
    D->AddIsotope(iso_H2, 1.);

    G4Material* C2D4 = new G4Material("C2D4", 1.06*g/cm3, 2);
    C2D4->AddElement(C, 2);
    C2D4->AddElement(D, 4);

    // Define Tungsten Wires
    G4Material* Tungsten = new G4Material("Tungsten", 19.3*g/cm3, nel = 1);
    Tungsten->AddElement(W, natoms = 1);

    // Define Gold Wires
    G4Material* Gold = new G4Material("Gold", 193*g/cm3, nel = 1);
    Gold->AddElement(Au, natoms = 1);

    G4double atmPressure = 760; // torr

    if(fGasType == "P10" || fGasType == "p10") {
        G4double p10Density = 0.00159*g/cm3;
        fGasMaterial = new G4Material("P10", p10Density*fPressureInTorr/atmPressure, nel = 3, kStateGas,
                                      fTemperature*kelvin, fPressureInTorr*1.333e-3*bar);
        fGasMaterial->AddElement(H,  fractionmass = 0.0155);
        fGasMaterial->AddElement(C,  fractionmass = 0.0623);
        fGasMaterial->AddElement(Ar, fractionmass = 0.9222);
    }
    else if(fGasType == "CO2" || fGasType == "co2") {
        fGasMaterial = G4Material::GetMaterial("CO2");
    }
    else if(fGasType == "Methane" || fGasType == "methane" || fGasType == "METHANE" || fGasType == "CH4") {
        fGasMaterial = G4Material::GetMaterial("Methane");
    }
    else if(fGasType == "CF4" || fGasType == "cf4") {
        G4double cf4Density = 0.0036586*g/cm3;
        fGasMaterial = new G4Material("CF4", cf4Density*fPressureInTorr/atmPressure, nel = 2, kStateGas,
                                      fTemperature*kelvin, fPressureInTorr*1.333e-3*bar);
        fGasMaterial->AddElement(C, natoms = 1);
        fGasMaterial->AddElement(F, natoms = 4);
    }
    else {
        G4cout << "Unknown Gas Type. Either a mistake or you need to add it yourself or ask to add" << G4endl;
        G4cout << "Choosing Methane for you" << G4endl;
        fGasMaterial = G4Material::GetMaterial("Methane");
    }
    //G4cout << GetGasMaterial() << G4endl;

    G4Material* vacuum = new G4Material("Vacuum",      //Name as String
                                        1,             //Atomic Number,  in this case we use 1 for hydrogen
                                        1.008*g/mole,  //Mass per Mole "Atomic Weight"  1.008*g/mole for Hydoren
                                        1.e-25*g/cm3,  //Density of Vaccuum  *Cant be Zero, Must be small instead
                                        kStateGas,     //kStateGas for Gas
                                        2.73*kelvin,   //Temperatuer for gas
                                        1.e-25*g/cm3); //Pressure for Vaccum

    // Create vacuum filled world
    G4VSolid* worldSolid = new G4Box("worldBox", 0.2*m, 0.15*m, 0.7*m);
    fWorldLogical = new G4LogicalVolume(worldSolid, vacuum, "worldLogical");
    G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0, G4ThreeVector(), fWorldLogical, "worldPhysical", 0,
                                                         false, 0, checkOverlaps);

    // C2D4 Target
    // 6.60 um for 0.7 mg/cm2
    // 3.77 um for 0.4 mg/cm2
    G4double targetThickness = 6.60*um;
    G4VSolid* targetSolid = new G4Tubs("targetSolid", 0., 20*mm, targetThickness/2., 0., 360.*deg);
    fTargetLogical =  new G4LogicalVolume(targetSolid, C2D4, "targetLogical");
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), fTargetLogical, "targetPhysical", fWorldLogical,
                      false, 0, checkOverlaps);

    G4double maxStep = 0.02*targetThickness;
    fStepLimit = new G4UserLimits(maxStep);
    fTargetLogical->SetUserLimits(fStepLimit);

    G4double icChamberLength = 0.3*m;
    G4double targetToWindow = 21.7*cm;

    // IC chamber
    G4VSolid* detectSolid = new G4Tubs("detectBox", 0., 0.1*m, icChamberLength/2., 0., 360.*deg);
    fDetectLogical = new G4LogicalVolume(detectSolid, fGasMaterial, "detectLogical");
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.367*m), fDetectLogical, "detectPhysical", fWorldLogical,
                      false, 0, checkOverlaps);

    G4double gridSize;
    G4double gridRadius;
    G4double gridDist;
    G4double foilRadius;
    G4double scintRadius;
    G4double scintDist;
    G4double wireThickness;

    if(fUseInches) {
        gridSize = fGridSize*inchtocm*cm;
        gridRadius = fGridRadius*inchtocm*cm;
        gridDist = fGridDist*inchtocm*cm;
        foilRadius = fGridRadius*inchtocm*cm;
        scintRadius = fGridRadius*inchtocm*cm;
        scintDist = fScintDist*inchtocm*cm;
        wireThickness = fWireThickness*inchtocm*cm;
    }
    else {
        gridSize = fGridSize*mm;
        gridRadius = fGridRadius*mm;
        gridDist = fGridDist*mm;
        foilRadius = fGridRadius*mm;
        scintRadius = fGridRadius*mm;
        scintDist = fScintDist*mm;
        wireThickness = fWireThickness*mm;
    }

    // Mylar foil
    G4double foilThickness = 13.*um;
    G4VSolid* foilSolid = foilSolid = new G4Tubs("foilSolid", 0., foilRadius, foilThickness/2., 0., 360.*deg);
    fFoilLogical = new G4LogicalVolume(foilSolid, Mylar, "foilLogical");
    new G4PVPlacement(0, G4ThreeVector(0., 0., targetToWindow - foilThickness/2.), fFoilLogical, "foilPhysical",
                      fWorldLogical, false, 0, checkOverlaps);

    G4double scintillatorDetectPos = -icChamberLength/2. + scintDist;
    G4double scintillatorThickness = 10.*mm;

    // IC grids
    char name[256];
    G4double distancePerGrid = gridSize;
    G4double distanceTotalGrid = distancePerGrid*static_cast<G4double>(fNumGrids);
    G4double midGrid = gridDist; // Picking halfway between foil and scintillator
    G4double midGridPos = -icChamberLength/2. + midGrid; // Position of middle of IC Grids
    // G4cout << distancePerGrid << '\t' << distanceTotalGrid << '\t' << foilToScintillator << '\t' << midGrid << G4endl;
    for(G4int i = 0; i < fNumGrids; i++) {
        sprintf(name, "grid%d", i + 1);
        G4VSolid* gridSolid = new G4Tubs(name, 0., gridRadius, distancePerGrid/2., 0., 360.*deg);
        sprintf(name, "gridLogical%d", i + 1);
        fGridLogical.push_back(new G4LogicalVolume(gridSolid, fGasMaterial, name));
        sprintf(name, "gridPhysical%d", i + 1);
        if(fNumGrids % 2 == 0) {
            int midGridNum = fNumGrids/2;
            if(i < midGridNum) {
                new G4PVPlacement(0, G4ThreeVector(0., 0., midGridPos + distancePerGrid/2. + (i - midGridNum)*distancePerGrid),
                                  fGridLogical[i], name, fDetectLogical, false, 0, checkOverlaps);
            }
            else {
                new G4PVPlacement(0, G4ThreeVector(0., 0., midGridPos + distancePerGrid/2. + (i - midGridNum)*distancePerGrid),
                                  fGridLogical[i], name, fDetectLogical, false, 0, checkOverlaps);
            }
        }
        else {
            int midGridNum = (fNumGrids - 1)/2 + 1;
            new G4PVPlacement(0, G4ThreeVector(0., 0., midGridPos + (i + 1 - midGridNum)*distancePerGrid), fGridLogical[i],
                              name, fDetectLogical, false, 0, checkOverlaps);
        }
    }

    // IC wires
    G4double wireRadius = wireThickness/2.;
    G4double wireSpacing = 2.794*mm; // 24 wires on each side of center for 49 total
    G4VSolid* gridSolid = new G4Tubs(name, 0., gridRadius, distancePerGrid/2., 0., 360.*deg);
    for(G4int i = 0; i < fNumGrids; i++) {
        sprintf(name, "wire_grid%d", i + 1);
        G4VSolid* wireSolid = new G4Tubs(name, 0., wireRadius, gridRadius*2., 0., 360.*deg);
        G4RotationMatrix *rm = new G4RotationMatrix();
        rm->rotateX(90.*deg);
        for(G4int j = -24; j < 25; j++) {
            sprintf(name, "wire_grid%d_%d", i + 1, j + 1);
            G4VSolid* wireIntersectSolid = new G4IntersectionSolid("Wire-Grid", gridSolid, wireSolid, rm, G4ThreeVector(wireSpacing*j, 0, 0.));
            sprintf(name, "wireGridLogical%d_%d", i + 1, j + 1);
            fWireGridLogical[j + 24].push_back(new G4LogicalVolume(wireIntersectSolid, Gold, name));
        }
        sprintf(name, "wireGridPhysical%d", i + 1);
        for(G4int j = -24; j < 25; j++) {
            G4double frontGrid = -distancePerGrid/2. + wireRadius;
            new G4PVPlacement(0, G4ThreeVector(0., 0, frontGrid), fWireGridLogical[j + 24][i],
                    name, fGridLogical[i], false, j + 24, checkOverlaps);
        }
        for(G4int j = -24; j < 25; j++) {
            new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), fWireGridLogical[j + 24][i],
                    name, fGridLogical[i], false, j + 24 + 49, checkOverlaps);
        }
        if(i == fNumGrids - 1) {
            for(G4int j = -24; j < 25; j++) {
                G4double backGrid = distancePerGrid/2. - wireRadius;
                new G4PVPlacement(0, G4ThreeVector(0., 0, backGrid), fWireGridLogical[j + 24][i],
                        name, fGridLogical[i], false, j + 24 + 49 + 49, checkOverlaps);
            }
        }
    }

    // BC-400 Plastic Scintillator
    G4VSolid* scintSolid = new G4Tubs("scintSolid", 0., scintRadius, scintillatorThickness/2., 0., 360.*deg);
    fScintLogical =  new G4LogicalVolume(scintSolid, BC400, "scintLogical");
    new G4PVPlacement(0, G4ThreeVector(0., 0., scintillatorDetectPos), fScintLogical, "scintPhysical", fDetectLogical,
                      false, 0, checkOverlaps);

    SetAttributes();

    return worldPhysical;
}
void G4TDetectorGeometry::ConstructMaterials() {
    G4NistManager* man = G4NistManager::Instance();

    man->ConstructNewGasMaterial("Methane", "G4_METHANE", fTemperature*kelvin,
                                 fPressureInTorr*1.333e-3*bar);
    man->ConstructNewGasMaterial("CO2", "G4_CARBON_DIOXIDE", fTemperature*kelvin,
                                 fPressureInTorr*1.333e-3*bar);
}
/*
void G4TDetectorGeometry::ConstructSDandField() {
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SDname;

  // Target
  char nameTarget[] = "target";
  G4VSensitiveDetector* targetDetector = new MMGenSD(SDname = nameTarget);
  SDman->AddNewDetector(targetDetector);
  fTargetLogical->SetSensitiveDetector(targetDetector);

  // Foil
  char nameFoil[] = "foil";
  G4VSensitiveDetector* foilDetector = new MMGenSD(SDname = nameFoil);
  SDman->AddNewDetector(foilDetector);
  fFoilLogical->SetSensitiveDetector(foilDetector);

  // IC Grids
  for(G4int i = 0; i < fNumGrids; i++) {
    char nameGrid[256];
    sprintf(nameGrid,"grid%d", i + 1);
    G4VSensitiveDetector* gridDetector = new MMGenSD(SDname = nameGrid);
    SDman->AddNewDetector(gridDetector);
    fGridLogical[i]->SetSensitiveDetector(gridDetector);
  }

  // BC-400 Plastic Scintillator
  char nameScint[] = "scint";
  G4VSensitiveDetector* scintDetector = new MMGenSD(SDname = nameScint);
  SDman->AddNewDetector(scintDetector);
  fScintLogical->SetSensitiveDetector(scintDetector);
}
*/
void G4TDetectorGeometry::SetAttributes() {
    G4VisAttributes* worldAttr = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    worldAttr->SetVisibility(true);
    worldAttr->SetForceWireframe(true);
    fWorldLogical->SetVisAttributes(worldAttr);

    G4VisAttributes* detectAttr = new G4VisAttributes(G4Colour(0., 0., 1.0));
    detectAttr->SetVisibility(true);
    detectAttr->SetForceWireframe(true);
    fDetectLogical->SetVisAttributes(detectAttr);

    G4VisAttributes* targetAttr = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    targetAttr->SetVisibility(true);
    targetAttr->SetForceSolid(true);
    fTargetLogical->SetVisAttributes(targetAttr);

    G4VisAttributes* foilAttr = new G4VisAttributes(G4Colour::White());
    foilAttr->SetVisibility(true);
    foilAttr->SetForceSolid(true);
    fFoilLogical->SetVisAttributes(foilAttr);

    G4Colour gridColors[10] = {G4Colour::Yellow(), G4Colour::Magenta(), G4Colour::Cyan(), G4Colour::Green(), G4Colour::White(),
                               G4Colour::Yellow(), G4Colour::Magenta(), G4Colour::Cyan(), G4Colour::Green(), G4Colour::White()};

    for(G4int i = 0; i < fNumGrids; i++) {
        G4VisAttributes* gridAttr = new G4VisAttributes(gridColors[i]);
        gridAttr->SetVisibility(true);
        gridAttr->SetForceSolid(true);
        fGridLogical[i]->SetVisAttributes(gridAttr);
    }

    for(G4int i = 0; i < fNumGrids; i++) {
        G4VisAttributes* wireAttr = new G4VisAttributes(G4Colour::Red());
        wireAttr->SetVisibility(true);
        wireAttr->SetForceSolid(true);
        for(G4int j = 0; j < 49; j++) {
            fWireGridLogical[j][i]->SetVisAttributes(wireAttr);
        }
    }

    G4VisAttributes* scintAttr = new G4VisAttributes(G4Colour::Red());
    scintAttr->SetVisibility(true);
    scintAttr->SetForceSolid(true);
    fScintLogical->SetVisAttributes(scintAttr);
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
