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

#include "G4TNuclearReactorGeometry.hh"
#include "G4TDetectorGeometry.hh"
#include "G4TPassiveProtonBeamLineGeometry.hh"

extern G4ThreeVector PointSourcePositionTowardDetector;
extern G4double ThetaDirectionValueTowardDetector;
extern G4double PhiDirectionValueTowardDetector;

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

    // ==== WORLD ====
    G4double worldSize = 6.0*m;
    auto worldSolid = new G4Box("World", worldSize/2, worldSize/2, worldSize/2);
    auto worldLogic = new G4LogicalVolume(worldSolid, Air, "WorldLogic");
    WorldPhysicalVolume = new G4PVPlacement(nullptr, G4ThreeVector(), worldLogic, "WorldPhys", nullptr, false, 0, true);

    // ==== PIPE ====
    G4double pipeRadiusOuter = 16.0*cm;
    G4double pipeRadiusInner = 15.0*cm;
    G4double pipeLength      = 50*cm;

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
        "Detector"
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



    // SOURCE POSITION
    //G4ThreeVector SrcPos(-20*cm , 0*cm , -30*cm);
    //G4double SourceDirectionTheta = 20*deg;


    // ==== Materials ====
    //auto Water          = nist->FindOrBuildMaterial("G4_WATER");
    //auto StainlessSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    //auto Glass          = nist->FindOrBuildMaterial("G4_GLASS_PLATE");
    //auto Concrete       = nist->FindOrBuildMaterial("G4_CONCRETE");
    //auto SourceMaterial = Air; // Placeholder, replace with real radioactive material

    //// ==== HIGH DENSITY MATERIAL ====

    //// Get Tungsten element
    //auto Tungsten = nist->FindOrBuildMaterial("G4_W");

    //G4Element* W = nist->FindOrBuildElement("W");

    //// Get Ultra Dense Tungsten element
    //G4Material* UltraDenseW = new G4Material("UltraDenseW", 15.0 * g/cm3, 1, kStateSolid);
    //UltraDenseW->AddElement(W, 1.0);
    //// Very High Density Tungsten
    //G4Material* SuperDense = new G4Material("SuperDense", 30.0 * g/cm3, 1, kStateSolid);
    //SuperDense->AddElement(W, 1.0);


    //// Define Osmium material manually
    //G4Element* Os = nist->FindOrBuildElement("Os");
    //G4Material* OsmiumMat = new G4Material("OsmiumMat", 22.6 * g/cm3, 1, kStateSolid);
    //OsmiumMat->AddElement(Os, 1.0);



    // ==== RANDOM VOLUMES INSIDE PIPE (HIGH ATTENUATION) ====


    // -20  20

    //G4double objSizeX = 20*cm;
    //G4double objSizeY = 20*cm;
    //G4double objSizeZ = 40*cm;   // BIG delta Z

    //auto bigSolid = new G4Box("BigObject",
    //                          objSizeX/2,
    //                          objSizeY/2,
    //                          objSizeZ/2);

    //auto bigLogic = new G4LogicalVolume(bigSolid,
    //                                    Pb,
    //                                    "BigObjectLogic");

    //// Place it centered in pipe
    //new G4PVPlacement(nullptr,
    //                  G4ThreeVector(0, 0, 0),
    //                  bigLogic,
    //                  "BigObjectPhys",
    //                  worldLogic,
    //                  false,
    //                  0,
    //                  true);



    //// 40  80

    //G4double radius = 10*cm;
    //G4double height = 40*cm;   // full height in Z

    //auto bigSolid1 = new G4Tubs("BigObject1",
    //                           0.0,          // inner radius (solid cylinder)
    //                           radius,       // outer radius
    //                           height/2.0,   // half length in Z
    //                           0.*deg,       // start angle
    //                           360.*deg);    // full circle

    //auto bigLogic1 = new G4LogicalVolume(bigSolid1, Pb, "BigObjectLogic1");

    //// Place it centered in pipe
    //new G4PVPlacement(nullptr,
    //                  G4ThreeVector(0, 0, 60*cm),//60*cm
    //                  bigLogic1,
    //                  "BigObjectPhys1",
    //                  worldLogic,
    //                  false,
    //                  0,
    //                  true);


    // ===================================================
    // =============== SOURCES & DETECTORS ===============
    // ===================================================

    // ==== SOURCE CAPSULE + DISK ====
    //auto srcCapsuleSolid = new G4Tubs("SrcCapsule", 2, 5*mm, 3*mm, 0, 360*deg);
    //auto srcDiskSolid    = new G4Tubs("SrcDisk", 0, 2*mm, 3*mm, 0, 360*deg);

    //auto srcCapsuleLogic = new G4LogicalVolume(srcCapsuleSolid, StainlessSteel, "SrcCapsuleLogic");
    //auto srcDiskLogic    = new G4LogicalVolume(srcDiskSolid, SourceMaterial, "SrcDiskLogic");

    //// Shielding layer around source capsule to block gamma outside collimators
    ////auto shieldSolid = new G4Tubs("ShieldLayer", 0, 6*mm, 3.5*mm, 0, 360*deg);
    ////auto shieldLogic = new G4LogicalVolume(shieldSolid, Pb, "ShieldLayerLogic");

    //// Place disk inside capsule slightly shifted in +Z direction (local)
    //new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 1.5*mm), srcDiskLogic, "SrcDiskPhys", srcCapsuleLogic, false, 0, true);

    //// ==== COLLIMATOR ====
    //auto colOuter = new G4Tubs("CollimatorOuter", 0, 3*mm, 5*mm, 0, 360*deg);
    //auto colHole  = new G4Tubs("CollimatorHole", 0, 1*mm, 6*mm, 0, 360*deg);
    //auto colSolid = new G4SubtractionSolid("CollimatorSolid", colOuter, colHole);

    //auto colLogic = new G4LogicalVolume(colSolid, Pb, "CollimatorLogic");

    //// ==== DETECTOR (Housing + NaI crystal + Photocathode) ====
    //auto detHousingSolid = new G4Tubs("DetHousing", 3, 3.2*cm, 5*cm, 0, 360*deg);
    //auto detNaISolid     = new G4Tubs("NaICrystal", 0, 3*cm, 4*cm, 0, 360*deg);
    //auto detPCSolid      = new G4Tubs("Photocathode", 0, 3*cm, 1*mm, 0, 360*deg);


    // Source data

    //std::cout << "\n\n\n\n\n\n\n\nn ---------------------- PointSourcePositionTowardDetector " << PointSourcePositionTowardDetector << "\n"
    //    << " ThetaDirectionValueTowardDetector "<< ThetaDirectionValueTowardDetector <<  std::endl;



    /*

    G4double R_det  = 45*cm;
    G4double Zslice = 0*cm;
    G4int    N      = 30;

    // Cylindrical detector (axis = local Z)
    auto detSolid = new G4Tubs("Detector",
                               0*cm,      // inner radius
                               2*cm,      // outer radius
                               3*cm,      // half length along axis
                               0,
                               360*deg);

    for(int i=0; i<N; i++)
    {

        auto detLogic = new G4LogicalVolume(detSolid, NaI, "Detector_"+std::to_string(i));

        G4double phi = 2*CLHEP::pi*i/N;

        // --- Position ---
        G4ThreeVector position(
            R_det * std::cos(phi),
            R_det * std::sin(phi),
            Zslice
            );

        // --- Direction toward center (disk normal) ---
        //G4ThreeVector axisDirection = -position.unit();

        // --- Individual rotation matrix ---
        G4ThreeVector zAxis(0,0,1);
        G4ThreeVector axisDirection = -position.unit();

        G4RotationMatrix rot;
        rot = G4RotationMatrix();
        rot.rotate(zAxis.angle(axisDirection),
                   zAxis.cross(axisDirection));

        // --- Individual full transform ---
        G4Transform3D transform(rot, position);

        new G4PVPlacement(transform,
                          detLogic,
                          "DetectorPhys",
                          worldLogic,
                          false,
                          i,
                          true);
    }


*/





    // Manually placed and in one geometry , 30 source 40 detector for 11 angles

    /*

    // ==== ARRANGEMENT PARAMETERS ====

    const G4int Nsrc = 30;
    const G4int NdetMain = 30;
    const G4int NdetExtra = 5;

    const G4int Ndet = NdetMain + 2*NdetExtra; // = 40


    // Pipe limits
    G4double pipeZmin = -pipeLength/2;
    G4double pipeZmax =  pipeLength/2;


    // Z range for sources (inside pipe)
    G4double zMinS = pipeZmin + 5*cm;
    G4double zMaxS = pipeZmax - 5*cm;

    G4double dzS = (zMaxS - zMinS)/(Nsrc-1);


    // Z range for main detectors (inside pipe)
    G4double zMinDmain = zMinS;
    G4double dzDmain   = dzS;


    // Z range for extra detectors
    G4double dzExtra = dzS;


    // X positions
    G4double xSource   = -(pipeRadiusOuter + 6*cm);
    G4double xDetector =  (pipeRadiusOuter + 7*cm);


     ===================================================
     ================== SOURCES ========================
     ===================================================

    for (G4int i=0;i<Nsrc;i++)
    {
        G4double zS = zMinS + i*dzS;

        // Rotation: face detectors
        auto srcRot = new G4RotationMatrix();
        srcRot->rotateY(90*deg);

        new G4PVPlacement(
            srcRot,
            G4ThreeVector(xSource,0,zS),
            srcCapsuleLogic,
            "CapsulePhys",
            worldLogic,
            false,
            i,
            true);
    }



     ===================================================
     ================== DETECTORS ======================
     ===================================================

    for (G4int j=0;j<Ndet;j++)
    {
        G4double zD;

        // ---- Bottom extra detectors
        if (j < NdetExtra)
        {
            zD = zMinS - (NdetExtra-j)*dzExtra;
        }

        // ---- Main detectors (inside pipe)
        else if (j < NdetExtra + NdetMain)
        {
            G4int k = j - NdetExtra;
            zD = zMinS + k*dzDmain;
        }

        // ---- Top extra detectors
        else
        {
            G4int k = j - (NdetExtra + NdetMain);
            zD = zMaxS + (k+1)*dzExtra;
        }


        // Rotation toward sources
        auto detRot = new G4RotationMatrix();
        detRot->rotateY(-90*deg);


        // Detector housing
        auto detHousingLogic =
            new G4LogicalVolume(detHousingSolid,
                                Aluminum,
                                "DetH_"+std::to_string(j));


        // NaI crystal
        auto detNaILogic =
            new G4LogicalVolume(detNaISolid,
                                NaI,
                                "Detector_"+std::to_string(j));


        // Place NaI
        new G4PVPlacement(
            nullptr,
            G4ThreeVector(0,0,-5*mm),
            detNaILogic,
            "NaIPhys_"+std::to_string(j),
            detHousingLogic,
            false,0,true);


        // Place detector
        new G4PVPlacement(
            detRot,
            G4ThreeVector(xDetector,0,zD),
            detHousingLogic,
            "DetectorPhys_"+std::to_string(j),
            worldLogic,
            false,
            j,
            true);
    }
*/


    return WorldPhysicalVolume;

}
