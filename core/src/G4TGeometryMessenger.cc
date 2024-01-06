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

#include "G4TGeometryMessenger.hh"
#include "G4TVolumeConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "G4Tokenizer.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G4SolidStore.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Para.hh"
#include "G4Trd.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Torus.hh"
#include "G4Ellipsoid.hh"

extern std::string getFileNameFromPath(std::string const & path, std::string const & delims = "/\\");
extern std::string getFileExt(const std::string& s);
//extern G4double UseG4Units(G4String Unit);
G4double G4TGeometryMessenger::UseG4Units(G4String Unit){

    G4double U = 1;

    if(Unit == "mm"){U = mm;}
    else if(Unit == "cm"){U = cm;}
    else if(Unit == "m"){U = m;}
    else if(Unit == "g/cm3"){U = g/cm3;}
    else if(Unit == "mg/cm3"){U = mg/cm3;}
    else if(Unit == "mg/mm3"){U = mg/mm3;}
    else if(Unit == "kg/m3"){U = kg/m3;}
    else if(Unit == "mg/mL"){U = mg/mL;}
    else if(Unit == "g/mL"){U = g/mL;}
    else if(Unit == "degree"){U = degree;}
    else if(Unit == "radian"){U = radian;}
    else if(Unit == "eV"){U = eV;}
    else if(Unit == "keV"){U = keV;}
    else if(Unit == "MeV"){U = MeV;}
    else if(Unit == "g/mole"){U = g/mole;}
    return  U;
}

G4TGeometryMessenger::G4TGeometryMessenger(G4TVolumeConstruction* GC) : GeometryConstruction(GC)
{
    //G4cout << "************** The primary generator Messenger... "<< G4endl;

    // materials
    MaterialConstruction = new G4UIdirectory("/MaterialData/");
    MaterialConstruction->SetGuidance("");

    // volumes
    ConstructGeometryCMD = new G4UIdirectory("/GeometryData/");
    ConstructGeometryCMD->SetGuidance("");

    CommandsForMaterials();
    CommandsForGeometries();
    CommandsForVoxGeometryType();
    /*
    CommandsForVolumes();
    CommandsForSolids();
    CommandsForDcmGeometryType();
    */
}


G4TGeometryMessenger::~G4TGeometryMessenger()
{


    delete ConstructGeometryCMD;
    delete MaterialConstruction;

    delete VoxDefaultMaterialNameCMD ;
    delete GeometrySymbolCMD ;
    delete VoxContainerPosCMD ;
    delete VoxContainerRotCMD ;

    delete VoxelizedRegionDataCMD;
    delete VoxelsDataCMD;
    delete NewElementCMD;
    delete NewMaterialCMD;
    delete MatElementCMD;
    delete MatMaterialsCMD;
    delete NistMaterialNameCMD;

    delete CreateWorldCMD;

    delete DcmPixelsCompressionCMD;
    delete DICOM_CtDensityCMD;
    delete DICOM_MaterialNamesOrderedCMD;
    delete DICOMFileType_FilesPathCMD;

    delete CreateSolidCMD;
    delete CreateVolumeCMD;
    delete DumpGeomCMD;
    delete VoxelizedPlanesCMD;
    delete setVolumesVisToNotForced;

    delete MaterialNameAsRegionNameCMD ;
    delete TETRegionDataCMD;
    delete setTETPhantomLimitsCMD;

}

// called automatically for each command in the begining of running to set the values entered by the user in Macro file of by typing commands from
// it call G4TVolumeConstruction::SetPhantomModel(G4String newModel), G4TVolumeConstruction::SetPhantomSex(G4String newSex)
// it call PrimaryGeneratorMessenger::AddBodyPart(G4String newBodyPartSensitivity)
void G4TGeometryMessenger::SetNewValue(G4UIcommand* command,G4String newValue){

    if( command == NewElementCMD)
    {
        G4Tokenizer next(newValue);
        G4double Z = StoD(next()), A = StoD(next()); G4String name = next(); //G4String Un1 = next();
        GeometryConstruction->createElement(name, Z, A*UseG4Units("g/mole"));
    }

    if( command == NewMaterialCMD)
    {
        G4Tokenizer next(newValue);
        G4String name = next(); G4int i = StoI(next()); MatElementNum = StoI(next()); G4double d = StoD(next()); G4String Un1 = next(), fracOr = next();
        GeometryConstruction->createMaterialFromComponents(name, i, MatElementNum, d*UseG4Units(Un1), fracOr);
    }

    if( command == MatElementCMD)
    {
        G4Tokenizer next(newValue);
        for(G4int ds = 0 ; ds < MatElementNum ; ds++){
            G4String ss = next();
            G4double vv =StoD(next()); G4int nn = vv;
            GeometryConstruction->AddMatElement(ss, vv , nn);
        }
    }

    if( command == MatMaterialsCMD)
    {
        G4Tokenizer next(newValue);
        for(G4int ds = 0 ; ds < MatElementNum ; ds++){

            G4String ss = next();
            G4double nn = StoD(next());
            GeometryConstruction->AddMatElement(ss, nn , nn);
        }
    }

    if( command == NistMaterialNameCMD)
    {
        G4Tokenizer next(newValue);

        G4String ss = next();
        G4double nn = StoI(next());
        GeometryConstruction->createNistMaterial(ss, nn);
    }


    if( command == CreateWorldCMD )
    {
        G4Tokenizer next(newValue);
        G4String mn = next() ;

        G4String ext = getFileExt(mn);

        if(ext == "gdml" || ext == "geom" || ext == "cpp" || ext == "c++"){

            //std::cout << "\n\n 1-------------- "  << std::endl;
            GeometryConstruction->placeVolume(mn, "", G4ThreeVector(), G4ThreeVector());
            //std::cout << "\n\n 2-------------- "  << std::endl;

        }else{
            //std::cout << "\n\n 1-------------- "  << std::endl;
            G4double ff = StoD(next()), ss = StoD(next()), gg = StoD(next()); G4String dc = next();
            GeometryConstruction->constructWorldVolume(mn, G4ThreeVector(ff*UseG4Units(dc),ss*UseG4Units(dc),gg*UseG4Units(dc)));
            //std::cout << "\n\n 2-------------- "  << std::endl;
        }
    }

    if( command == CreateSolidCMD){

        G4Tokenizer next(newValue);

        G4String SolidType = next();
        G4String SolidName = next();
        G4String ext = getFileExt(SolidType);
        G4String fn = getFileNameFromPath(SolidType);

        G4double A, B, C, D, E, F, G, H, I; G4String Un1, Un2, Un3;


        if(SolidType == "Box"){

            A = StoD(next()), B = StoD(next()), C = StoD(next());
            Un1 = next();
            G4VSolid* sld = new G4Box(SolidName, A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un1));
        }
        if(SolidType == "Tubs"){

            A = StoD(next()), B = StoD(next()), C = StoD(next()), D = StoD(next()), E = StoD(next());
            Un1 = next(), Un2 = next();
            sld = new G4Tubs(SolidName, A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un1), D*UseG4Units(Un2), E*UseG4Units(Un2));
            //new G4Tubs(SolidName, Rmin, Rmax, Dz, SPhi, DPhi);
        }
        if(SolidType == "CutTubs"){

            A = StoD(next()), B = StoD(next()), C = StoD(next()), D = StoD(next()), E = StoD(next());
            Un1 = next(), Un2 = next();
            //sld = new G4CutTubs( SolidName, A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un1), D*UseG4Units(Un2), E*UseG4Units(Un2), LowNorm, HighNorm );
            //new G4CutTubs( SolidName, Rmin, Rmax, Dz, SPhi, DPhi, LowNorm, HighNorm );
        }
        if(SolidType == "Cons"){

            A = StoD(next()), B = StoD(next()), C = StoD(next()), D = StoD(next()), E = StoD(next()), F = StoD(next()), G = StoD(next());
            Un1 = next(), Un2 = next();
            sld = new G4Cons( SolidName, A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un1), D*UseG4Units(Un1), E*UseG4Units(Un1), F*UseG4Units(Un2), G*UseG4Units(Un2));
            //new G4Cons( SolidName, Rmin1, Rmax1, Rmin2, Rmax2, Dz, SPhi, DPhi);
        }
        if(SolidType == "Para"){

            A = StoD(next()), B = StoD(next()), C = StoD(next()), D = StoD(next()), E = StoD(next()), F = StoD(next());
            Un1 = next(), Un2 = next();
            sld = new G4Para(SolidName , A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un1), D*UseG4Units(Un2), E*UseG4Units(Un2), F*UseG4Units(Un2));
            //new G4Para(SolidName , Dx, Dy, Dz, Alpha, Theta0, Phi0);
        }
        if(SolidType == "Trd"){

            A = StoD(next()), B = StoD(next()), C = StoD(next()), D = StoD(next()), E = StoD(next());
            Un1 = next();
            sld = new G4Trd( SolidName, A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un1), D*UseG4Units(Un1), E*UseG4Units(Un1));
            //new G4Trd( SolidName, Dx1, Dx2, Dy1, Dy2, Dz);
        }
        if(SolidType == "Sphere"){

            A = StoD(next()), B = StoD(next()), C = StoD(next()), D = StoD(next()), E = StoD(next()), F = StoD(next());
            Un1 = next(), Un2 = next();
            sld = new G4Sphere( SolidName, A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un2), D*UseG4Units(Un2), E*UseG4Units(Un2), F*UseG4Units(Un2) );
            //new G4Sphere( SolidName, Rmin, Rmax, SPhi, DPhi, STheta, DTheta );
        }
        if(SolidType == "Orb"){

            A = StoD(next());
            Un1 = next();
            sld = new G4Orb( SolidName, A*UseG4Units(Un1));
            //new G4Orb( SolidName, Rmax);
        }
        if(SolidType == "Torus"){

            A = StoD(next()), B = StoD(next()), C = StoD(next()), D = StoD(next()), E = StoD(next());
            Un1 = next(), Un2 = next();
            sld = new G4Torus( SolidName, A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un1), D*UseG4Units(Un2), E*UseG4Units(Un2));
            //new G4Torus( SolidName, Rmin, Rmax, Rtor, SPhi, DPhi);
        }
        if(SolidType == "Ellipsoid"){

            A = StoD(next()), B = StoD(next()), C = StoD(next()), D = 0, E = 0;
            Un1 = next();
            sld = new  G4Ellipsoid(SolidName , A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un1), D*UseG4Units(Un1), E*UseG4Units(Un1));
            //new  G4Ellipsoid(SolidName , xSemiAxis, ySemiAxis, zSemiAxis, zBottomCut, zBottomCut);
        }
        if(SolidType == "Union" || SolidType == "Intersection" || SolidType == "Subtraction"){

            G4String AA =next(), BB =next(); C=StoD(next()), D=StoD(next()), E=StoD(next()), F = StoD(next()), G = StoD(next()), H = StoD(next());
            Un1 = next(), Un2 = next();

            G4SolidStore* ss = G4SolidStore::GetInstance();
            G4VSolid* s1  = ss->GetSolid(AA);
            G4VSolid* s2  = ss->GetSolid(BB);
            G4RotationMatrix* rm = new G4RotationMatrix(); rm->rotateX(F*UseG4Units(Un2)); rm->rotateY(G*UseG4Units(Un2)); rm->rotateZ(H*UseG4Units(Un2));
            if(SolidType == "Union"){
                sld = new G4UnionSolid( SolidName , s1, s2, rm, G4ThreeVector(C*UseG4Units(Un1), D*UseG4Units(Un1), E*UseG4Units(Un1)));
            } else if(SolidType == "Intersection"){
                sld = new G4IntersectionSolid( SolidName , s1, s2, rm, G4ThreeVector(C*UseG4Units(Un1), D*UseG4Units(Un1), E*UseG4Units(Un1)));
            } else if(SolidType == "Subtraction"){
                sld = new G4SubtractionSolid( SolidName , s1, s2, rm, G4ThreeVector(C*UseG4Units(Un1), D*UseG4Units(Un1), E*UseG4Units(Un1)));
            }
        }
    }

    if( command == CreateVolumeCMD) {

        G4Tokenizer next(newValue);

        G4String volN = next();

        if(volN == "GDML" || volN == "TEXT" || volN == "VoxIDs" ){
            G4String volS = next();
            GeometryConstruction->setGeometryFileType(volN);
            GeometryConstruction->setGeometryPath(volS);
        }
        else if(volN == "DICOM" || volN == "VOXEL" || volN == "CPP"){
            GeometryConstruction->setGeometryFileType(volN);
        }
        else if(volN == "TET"){
            GeometryConstruction->setGeometryFileType(volN);
            G4String node = next();
            G4String ele = next();
            G4String mat = next();
            GeometryConstruction->setGeometryTETDataFiles(node, ele, mat);
        }
        else{

            if(volN == "FILEGDML"){
                G4String vp = next();
                G4String mv = next();
                GeometryConstruction->PlaceVolumeFromRegionfile(vp, mv);
                return;
            }
            else if(volN == "Generate"){
                GeometryConstruction->setGeometryFileType(volN);
                GeometryConstruction->setGeometryPath(next());
                return;
            }
            else {

                GeometryConstruction->setGeometryFileType("Construct");

                G4String ext = getFileExt(volN);
                G4String fn = getFileNameFromPath(volN);

                if(fn == "World"){
                    GeometryConstruction->placeVolume(volN, "", G4ThreeVector(), G4ThreeVector());
                    return;
                }
                if(ext == "gdml" || ext == "geom" || ext == "cpp" || ext == "c++"){

                }
                else if(ext == "stl" || ext == "ast"){
                    //GeometryConstruction->setVolumeName(volN);

                    G4String volM = next();
                    GeometryConstruction->setSTLVolumeMaterial(volM);
                }
                else{

                    G4String volS = next();
                    G4String volM = next();
                    G4MaterialTable Mates = *G4Material::GetMaterialTable();
                    for (G4int n = 0; n < Mates.size(); n++ ){
                        if(Mates[n]->GetName() == volM){
                            new G4LogicalVolume(G4SolidStore::GetInstance()->GetSolid(volS), Mates[n], volN, 0, 0,0);
                        }
                    }
                }


                G4String motherVol = next();

                G4double PosX = StoD(next());
                G4double PosY = StoD(next());
                G4double PosZ = StoD(next());

                G4double RotX = StoD(next());
                G4double RotY = StoD(next());
                G4double RotZ = StoD(next());

                G4String Un1 = next(), Un2 = next();

                //std::cout << " ---------------------- ext " << ext << " volN " << volN <<  std::endl;

                GeometryConstruction->placeVolume(volN,
                                                  motherVol,
                                                  G4ThreeVector(PosX*UseG4Units(Un1),PosY*UseG4Units(Un1),PosZ*UseG4Units(Un1)),
                                                  G4ThreeVector(RotX*UseG4Units(Un2),RotY*UseG4Units(Un2),RotZ*UseG4Units(Un2)));

            }

        }
    }

    if( command == DumpGeomCMD )
    {
        G4Tokenizer next(newValue);
        G4String fn = "";
        if(fn == next()){
            if(getFileExt(fn) == "geom"){
                fn = fn + ".geom";
            }
        }
        else{
            fn = "Geometry.geom";
        }
        GeometryConstruction->dumpGeometryToTextFile(fn);
    }

    if( command == setVolumesVisToNotForced )
    {
        G4Tokenizer next(newValue);
        G4String VolName = next();

        while (!VolName.empty()) {
            GeometryConstruction->setVolumesNotVisualized(VolName);
            VolName = next();
        }

    }

    if( command == GeometrySymbolCMD )
    {
        GeometryConstruction->setGeometrySymbol(newValue);
    }

    // voxelized Volume

    if( command == VoxDefaultMaterialNameCMD )
    {
        //GeometryConstruction->setVoxDefaultMaterialName(newValue);
    }
    if( command == VoxContainerPosCMD )
    {
        GeometryConstruction->setVoxContainerPos(VoxContainerPosCMD->GetNew3VectorValue(newValue));
    }
    if( command == VoxContainerRotCMD )
    {
        GeometryConstruction->setVoxContainerRot(VoxContainerRotCMD->GetNew3VectorValue(newValue));
    }

    if( command == VoxelsDataCMD)
    {
        G4Tokenizer next(newValue);

        GeometryConstruction->setVoxXNumber(StoI(next()));
        GeometryConstruction->setVoxYNumber(StoI(next()));
        GeometryConstruction->setVoxZNumber(StoI(next()));

        G4int Parr = StoI(next());
        GeometryConstruction->setParamType(Parr);

        //std::cout << "\n\n\n From Messenger Parr " << Parr << "\n\n\n" <<  std::endl;

        G4String DefVoxMat = next();

        G4double X = StoD(next()), Y = StoD(next()), Z = StoD(next()); G4String Un = next();
        GeometryConstruction->setVoxXHalfSize(X*UseG4Units(Un));
        GeometryConstruction->setVoxYHalfSize(Y*UseG4Units(Un));
        GeometryConstruction->setVoxZHalfSize(Z*UseG4Units(Un));

        GeometryConstruction->setVoxDefaultMaterialName(DefVoxMat);
    }

    if( command == VoxelizedRegionDataCMD)
    {
        G4Tokenizer next(newValue);
        GeometryConstruction->setVoxRegionName(next());
        G4String PosWord = next();
        if(PosWord=="pos"){
            G4String PosUnit = next();
            GeometryConstruction->setVoxRegionMinX(GeometryConstruction->ConvertVoxelPosToXYZVoxelIDs("X","Min",StoD(next())*UseG4Units(PosUnit)));
            GeometryConstruction->setVoxRegionMaxX(GeometryConstruction->ConvertVoxelPosToXYZVoxelIDs("X","Max",StoD(next())*UseG4Units(PosUnit)));
            GeometryConstruction->setVoxRegionMinY(GeometryConstruction->ConvertVoxelPosToXYZVoxelIDs("Y","Min",StoD(next())*UseG4Units(PosUnit)));
            GeometryConstruction->setVoxRegionMaxY(GeometryConstruction->ConvertVoxelPosToXYZVoxelIDs("Y","Max",StoD(next())*UseG4Units(PosUnit)));
            GeometryConstruction->setVoxRegionMinZ(GeometryConstruction->ConvertVoxelPosToXYZVoxelIDs("Z","Min",StoD(next())*UseG4Units(PosUnit)));
            GeometryConstruction->setVoxRegionMaxZ(GeometryConstruction->ConvertVoxelPosToXYZVoxelIDs("Z","Max",StoD(next())*UseG4Units(PosUnit)));
        }
        /*else if(PosWord=="add"){

            G4String VolName = next();

            while (!VolName.empty()) {
                //GeometryConstruction->setVolumesNotVisualized(VolName);
                double fraction = StoD(next());
                GeometryConstruction->setVoxRegionsFractions(VolName,fraction);
                VolName = next();
            }
        }*/

        else if(PosWord=="all"){
            GeometryConstruction->setAllGeomAsMinMaxVoxRegionLimits();
        }else {
            GeometryConstruction->setVoxRegionMinX(StoI(PosWord));
            GeometryConstruction->setVoxRegionMaxX(StoI(next()));
            GeometryConstruction->setVoxRegionMinY(StoI(next()));
            GeometryConstruction->setVoxRegionMaxY(StoI(next()));
            GeometryConstruction->setVoxRegionMinZ(StoI(next()));
            GeometryConstruction->setVoxRegionMaxZ(StoI(next()));
        }

        G4String Val = next(), Val2 = next(), Val3 = next(), ValUn = next() ;
        if(Val=="null"){
            //std::cout << "\n VoxelsegmentedMaterial Val " << Val <<  std::endl;
        }else {GeometryConstruction->setVoxelsegmentedMaterial(StoI(Val));} // used just for VoxIDs

        if(Val2=="null"){
            //std::cout << "\n DcmRegionMinDensity Val " << Val <<  std::endl;
        }else {GeometryConstruction->setDcmRegionMinDensity(StoD(Val2)*UseG4Units(ValUn));} // used just for Dcm

        if(Val3=="null"){
            //std::cout << "\n DcmRegionMaxDensity Val " << Val <<  std::endl;
        }else {GeometryConstruction->setDcmRegionMaxDensity(StoD(Val3)*UseG4Units(ValUn));} // used just for Dcm

        //std::cout << "\n VoxelsegmentedMaterial Val " << Val <<" " << Val2 << " " << Val3 << " " << ValUn << std::endl;

    }
    if( command == VoxelizedPlanesCMD )
    {
        G4Tokenizer next(newValue);
        G4String planes = next(); planes.toLower();
        GeometryConstruction->setPlanesToVisualize(planes);
        if(planes != "all"){
            if(planes == "regions"){

                G4String VolName = next();
                while (!VolName.empty()) {
                    GeometryConstruction->setRegionsToVisualize(VolName);
                    VolName = next();
                }
            }else{
                GeometryConstruction->setMinPlaneID(StoI(next()));
                GeometryConstruction->setMaxPlaneID(StoI(next()));
            }
        }

//      G4String nn = next() ; nn.toLower();
//      if(nn == "yes"){
//          GeometryConstruction->setForcedSolid(true);
//      }

        //G4cout << " - plane " << planes << G4endl;

    }

    if( command == setTETPhantomLimitsCMD )
    {
        G4Tokenizer next(newValue);
        G4String planes = next(); planes.toLower();
        GeometryConstruction->setPlanesToVisualize(planes);
        if(planes != "all"){
            if(planes == "regions"){
                G4String VolName = next();
                while (!VolName.empty()) {
                    GeometryConstruction->setRegionsToVisualize(VolName);
                    VolName = next();
                }
            }else{
                G4double a = StoD(next());
                G4double b = StoD(next());
                GeometryConstruction->setTETPhantomLimits(a,b);
            }
        }
        //      G4String nn = next() ; nn.toLower();
        //      if(nn == "yes"){
        //          GeometryConstruction->setForcedSolid(true);
        //      }

        //G4cout << " - plane " << planes << G4endl;

    }
    if(command == TETRegionDataCMD )
    {
        G4Tokenizer next(newValue);
        GeometryConstruction->setVoxRegionName(next());

        G4String Val2 = next(), Val3 = next(), ValUn = next() ;
        if(Val2=="null"){
            //std::cout << "\n DcmRegionMinDensity Val " << Val <<  std::endl;
        }else {GeometryConstruction->setDcmRegionMinDensity(StoD(Val2)*UseG4Units(ValUn));} // used just for Dcm

        if(Val3=="null"){
            //std::cout << "\n DcmRegionMaxDensity Val " << Val <<  std::endl;
        }else {GeometryConstruction->setDcmRegionMaxDensity(StoD(Val3)*UseG4Units(ValUn));} // used just for Dcm

        //std::cout << "\n VoxelsegmentedMaterial Val " << Val <<" " << Val2 << " " << Val3 << " " << ValUn << std::endl;
    }


    if( command == MaterialNameAsRegionNameCMD )
    {
        G4Tokenizer next(newValue);
        G4String val = next(); val.toLower();
        if(val == "y" || val == "yes"){
            GeometryConstruction->setMaterialNameAsRegionName(true);
        }else{
            GeometryConstruction->setMaterialNameAsRegionName(false);
        }
    }


    if( command == DICOMFileType_FilesPathCMD)
    {
        G4Tokenizer next(newValue);
        G4String type = next();
        GeometryConstruction->setDicomFileType(type);
        GeometryConstruction->setDicomDataDirPath(next());

        if(type == "PET"){
            GeometryConstruction->setDcmPETSerieResidenceTime(StoD(next()));
            GeometryConstruction->setDcmPixelsCompressionXY(StoI(next()));
            //GeometryConstruction->setDcmPixelsCompressionY(StoI(next()));
        }
        else if(type == "CT"){
            GeometryConstruction->setDcmPixelsCompressionXY(StoI(next()));
            //GeometryConstruction->setDcmPixelsCompressionY(StoI(next()));
        }
    }

    if( command == DICOM_MaterialNamesOrderedCMD)
    {
        G4Tokenizer next(newValue);
        G4int zz = StoI(next());
        for(G4int ds = 0 ; ds < zz ; ds++){
            GeometryConstruction->setDcmMaterialName(next());
        }
    }

    if( command == DICOM_CtDensityCMD)
    {
        G4Tokenizer next(newValue);
        GeometryConstruction->setDcmCtNumber(StoI(next()));
        GeometryConstruction->setDcmCtDensity(StoD(next())*UseG4Units(next()));
    }

    if( command == DcmPixelsCompressionCMD )
    {
        GeometryConstruction->setDcmPixelsCompression(DcmPixelsCompressionCMD->GetNewIntValue(newValue));
    }

}

// called from constructor to define the commands waited to fill from user
void  G4TGeometryMessenger::CommandsForMaterials(){

    NewElementCMD = new G4UIcommand("/MaterialData/createElement" ,this);
    //NewMaterialCMD->SetGuidance("");
    G4UIparameter* param;
    param = new G4UIparameter("Z",'i', false);      NewElementCMD->SetParameter(param);
    param = new G4UIparameter("A",'d', false);      NewElementCMD->SetParameter(param);
    param = new G4UIparameter("Name",'s', false);   NewElementCMD->SetParameter(param);

    NewMaterialCMD = new G4UIcommand("/MaterialData/createMaterial" ,this);
    //NewMaterialCMD->SetGuidance("");
    param = new G4UIparameter("Name",'s', false);      NewMaterialCMD->SetParameter(param);
    param = new G4UIparameter("ID",'i', false);      NewMaterialCMD->SetParameter(param);
    param = new G4UIparameter("Elements Number",'i', false);   NewMaterialCMD->SetParameter(param);
    param = new G4UIparameter("Density",'d', false);      NewMaterialCMD->SetParameter(param);
    param = new G4UIparameter("Density unit",'s', false);      NewMaterialCMD->SetParameter(param);
    param = new G4UIparameter("Fraction or Number",'s', false);   NewMaterialCMD->SetParameter(param);

    MatElementCMD = new G4UIcommand("/MaterialData/addElements" ,this);
    for(G4int ds = 0 ; ds < 100 ; ds++){

        G4String f = "Element Name " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      MatElementCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
        f = "Fraction or Number Value For Elem " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      MatElementCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    }

    MatMaterialsCMD = new G4UIcommand("/MaterialData/addMaterials" ,this);
    for(G4int ds = 0 ; ds < 100 ; ds++){

        G4String f = "MatMaterial Name " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      MatMaterialsCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
        f = "Fraction or Number Value For MatMaterial " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      MatMaterialsCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    }

    NistMaterialNameCMD = new G4UIcommand("/MaterialData/setNistMaterialNameAndID",this);
    param = new G4UIparameter("Name",'s', false);      NistMaterialNameCMD->SetParameter(param);
    param = new G4UIparameter("ID",'i', false);      NistMaterialNameCMD->SetParameter(param);

}
void  G4TGeometryMessenger::CommandsForGeometries(){

    G4UIparameter* param;

    CreateWorldCMD = new G4UIcommand("/GeometryData/createWorld" ,this);
    param = new G4UIparameter("Volume Name",'s', true);     CreateWorldCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("World Material",'s', true);  CreateWorldCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("World SD",'s', true);        CreateWorldCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Volume HX ",'s', true);      CreateWorldCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Volume HY ",'s', true);      CreateWorldCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Volume HZ ",'s', true);      CreateWorldCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("H Unit ",'s', true);         CreateWorldCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    CreateSolidCMD = new G4UIcommand("/GeometryData/createSolid" ,this);
    param = new G4UIparameter("Solid Type",'s', true);      CreateSolidCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Solid Name",'s', true);      CreateSolidCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    for(G4int ds = 0 ; ds < 20 ; ds++){

        G4String f = "Solid Parameter " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      CreateSolidCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    }

    CreateVolumeCMD = new G4UIcommand("/GeometryData/createVolume" ,this);
    param = new G4UIparameter("Logical Volume Name",'s', true);      CreateVolumeCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Logical Volume Solid",'s', true);      CreateVolumeCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Logical Volume Material",'s', true);      CreateVolumeCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Mother Volume Name",'s', true);      CreateVolumeCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    param = new G4UIparameter("Volume Pos X ",'s', true);      CreateVolumeCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Volume Pos Y ",'s', true);      CreateVolumeCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Volume Pos Z ",'s', true);      CreateVolumeCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    param = new G4UIparameter("Volume Rot X ",'s', true);      CreateVolumeCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Volume Rot Y ",'s', true);      CreateVolumeCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Volume Rot Z ",'s', true);      CreateVolumeCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    param = new G4UIparameter("Volume Pos Unit ",'s', true);      CreateVolumeCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Volume Rot Unit ",'s', true);      CreateVolumeCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    DumpGeomCMD = new G4UIcommand("/GeometryData/dumpGeometry" ,this);
    param = new G4UIparameter("Geometry Output file Name",'s', true);      DumpGeomCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    setVolumesVisToNotForced = new G4UIcommand("/GeometryData/setVolumesVisToNotForced" ,this);
    for(G4int ds = 0 ; ds < 100 ; ds++){

        G4String f = "Vol " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      setVolumesVisToNotForced->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    }

    GeometrySymbolCMD = new G4UIcmdWithAString("/GeometryData/setGeometrySymbol",this);
    GeometrySymbolCMD->SetGuidance(" ");
    GeometrySymbolCMD->SetParameterName("GeometrySymbolCMD",true);
    //GeometrySymbolCMD->SetDefaultValue("G4_AIR");
    //GeometrySymbolCMD->SetCandidates("yes no");
    GeometrySymbolCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

}
void  G4TGeometryMessenger::CommandsForVoxGeometryType(){

    VoxDefaultMaterialNameCMD = new G4UIcmdWithAString("/GeometryData/setVoxDefaultMaterialName",this);
    VoxDefaultMaterialNameCMD->SetGuidance(" ");
    VoxDefaultMaterialNameCMD->SetParameterName("VoxDefaultMaterialNameCMD",true);
    VoxDefaultMaterialNameCMD->SetDefaultValue("G4_AIR");
    //VoxDefaultMaterialNameCMD->SetCandidates("yes no");
    VoxDefaultMaterialNameCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    VoxelsDataCMD = new G4UIcommand("/GeometryData/setVoxelsData" ,this);
    //NewMaterialCMD->SetGuidance("");
    G4UIparameter* param;
    param = new G4UIparameter("NX",'i', false);   VoxelsDataCMD->SetParameter(param);
    param = new G4UIparameter("NY",'i', false);   VoxelsDataCMD->SetParameter(param);
    param = new G4UIparameter("NZ",'i', false);   VoxelsDataCMD->SetParameter(param);
    param = new G4UIparameter("Parametrization type",'i', false);   VoxelsDataCMD->SetParameter(param);
    param = new G4UIparameter("Def Vox Mat",'s', false);   VoxelsDataCMD->SetParameter(param);
    param = new G4UIparameter("HX",'d', false);   VoxelsDataCMD->SetParameter(param);
    param = new G4UIparameter("HY",'d', false);   VoxelsDataCMD->SetParameter(param);
    param = new G4UIparameter("HZ",'d', false);   VoxelsDataCMD->SetParameter(param);
    param = new G4UIparameter("Half size unit",'s', false);   VoxelsDataCMD->SetParameter(param);

    VoxelizedRegionDataCMD = new G4UIcommand("/GeometryData/setVoxelizedRegionData" ,this);
    //NewMaterialCMD->SetGuidance("");
    param = new G4UIparameter("Name",'s', false);   VoxelizedRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("Pos, all or MinX",'s', false);   VoxelizedRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("MaxX",'s', true);   VoxelizedRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("MinY",'s', true);   VoxelizedRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("MaxY",'s', true);   VoxelizedRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("MinZ",'s', true);   VoxelizedRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("MaxZ",'s', true);   VoxelizedRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("MatID",'s',true);   VoxelizedRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("MinD",'s', true);   VoxelizedRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("MaxD",'s', true);   VoxelizedRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("D Unit",'s', true);  VoxelizedRegionDataCMD->SetParameter(param);


    DICOM_CtDensityCMD = new G4UIcommand("/GeometryData/setCtDensityValues" ,this);
    param = new G4UIparameter("CT Number",'i', false);      DICOM_CtDensityCMD->SetParameter(param);
    param = new G4UIparameter("Density",'d', false);      DICOM_CtDensityCMD->SetParameter(param);
    param = new G4UIparameter("Density Unit",'s', false);      DICOM_CtDensityCMD->SetParameter(param);

    DICOM_MaterialNamesOrderedCMD = new G4UIcommand("/GeometryData/setMatNumberAndMaterials" ,this);
    param = new G4UIparameter("Mat Number",'i', false);      DICOM_MaterialNamesOrderedCMD->SetParameter(param);
    for(G4int ds = 0 ; ds < 150 ; ds++){

        G4String f = "Material Name" + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      DICOM_MaterialNamesOrderedCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

    DICOMFileType_FilesPathCMD = new G4UIcommand("/GeometryData/setDcmTypeAndDataPath" ,this);
    param = new G4UIparameter("DICOM Type",'s', false);      DICOMFileType_FilesPathCMD->SetParameter(param);
    param = new G4UIparameter("DICOM Data Files Path",'s', false);      DICOMFileType_FilesPathCMD->SetParameter(param);

    DcmPixelsCompressionCMD = new G4UIcmdWithAnInteger("/GeometryData/setDcmPixelsCompression",this);
    DcmPixelsCompressionCMD->SetGuidance("Set the number of pixels to compress in one pixel");
    DcmPixelsCompressionCMD->SetParameterName("DcmPixelsCompressionCMD",true);
    DcmPixelsCompressionCMD->SetDefaultValue(2);
    DcmPixelsCompressionCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    VoxContainerPosCMD = new G4UIcmdWith3VectorAndUnit("/GeometryData/setVoxContainerPos",this);
    VoxContainerPosCMD->SetGuidance(" ");
    VoxContainerPosCMD->SetParameterName("x","y","z",true,true);
    VoxContainerPosCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    VoxContainerPosCMD->SetDefaultValue(G4ThreeVector(0.,0.,0.));
    VoxContainerPosCMD->SetDefaultUnit("mm"); // then we dont need to add *cm in command
    VoxContainerPosCMD->SetUnitCandidates("mm cm m");

    VoxContainerRotCMD = new G4UIcmdWith3VectorAndUnit("/GeometryData/setVoxContainerRot",this);
    VoxContainerRotCMD->SetGuidance(" ");
    VoxContainerRotCMD->SetParameterName("x","y","z",true,true);
    VoxContainerRotCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    VoxContainerRotCMD->SetDefaultValue(G4ThreeVector(0.,0.,0.));
    VoxContainerRotCMD->SetDefaultUnit("degree"); // then we dont need to add *cm in command
    VoxContainerRotCMD->SetUnitCandidates("degree");

    VoxelizedPlanesCMD = new G4UIcommand("/GeometryData/visualizeVoxelizedPlanes" ,this);
    param = new G4UIparameter("planes to visualize",'s', false);      VoxelizedPlanesCMD->SetParameter(param);
    param = new G4UIparameter("Min plane",'s', true);      VoxelizedPlanesCMD->SetParameter(param);
    param = new G4UIparameter("Max plane",'s', true);      VoxelizedPlanesCMD->SetParameter(param);
    param = new G4UIparameter("ForceSolid",'s', true);      VoxelizedPlanesCMD->SetParameter(param);

    setTETPhantomLimitsCMD = new G4UIcommand("/GeometryData/setTETPhantomLimits" ,this);
    param = new G4UIparameter("planes to visualize",'s', false);      setTETPhantomLimitsCMD->SetParameter(param);
    param = new G4UIparameter("Min plane",'s', true);      setTETPhantomLimitsCMD->SetParameter(param);
    param = new G4UIparameter("Max plane",'s', true);      setTETPhantomLimitsCMD->SetParameter(param);
    param = new G4UIparameter("ForceSolid",'s', true);      setTETPhantomLimitsCMD->SetParameter(param);

    MaterialNameAsRegionNameCMD = new G4UIcmdWithAString("/GeometryData/setMaterialNameAsRegionName",this);
    MaterialNameAsRegionNameCMD->SetParameterName("MaterialNameAsRegionNameCMD",true);
    MaterialNameAsRegionNameCMD->SetDefaultValue("yes");
    MaterialNameAsRegionNameCMD->SetGuidance(" ");
    MaterialNameAsRegionNameCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    TETRegionDataCMD = new G4UIcommand("/GeometryData/setTETRegionData" ,this);
    //NewMaterialCMD->SetGuidance("");
    param = new G4UIparameter("Name",'s', false);   TETRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("MatID",'s',true);   TETRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("MinD",'s', true);   TETRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("MaxD",'s', true);   TETRegionDataCMD->SetParameter(param);
    param = new G4UIparameter("D Unit",'s', true);  TETRegionDataCMD->SetParameter(param);


}
