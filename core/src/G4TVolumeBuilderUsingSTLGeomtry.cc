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

#include "G4RunManager.hh"
#include "G4tgbVolumeMgr.hh"

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"

#include "G4TVolumeConstruction.hh"
#include "G4TVolumeBuilderUsingSTLGeomtry.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalVolumeStore.hh"


#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "G4LogicalVolumeStore.hh"





G4TVolumeBuilderUsingSTLGeomtry::G4TVolumeBuilderUsingSTLGeomtry()
{
    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  G4TVolumeBuilderUsingSTLGeomtry  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    PhysicalBoxVolume = 0;
    LogicalBoxVolume = 0;


    //TEXT_GEO_File="";
    OrganMassFilePath = "";
    OrganSourceName = "";
    generatePoints = false;
    GenerateWithBox = false;
}
G4TVolumeBuilderUsingSTLGeomtry::~G4TVolumeBuilderUsingSTLGeomtry()
{
}


// called from G4TConstruction::construct()
//void G4TVolumeBuilderUsingSTLGeomtry::SetTEXT_GEO_File(G4String t) {TEXT_GEO_File = t; }
void G4TVolumeBuilderUsingSTLGeomtry::setGeometryData() {

    const G4TVolumeConstruction* TConstruction = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    //UseOrganMassFilePath = TConstruction->getUseRegionsMassDataPath();
    //OrganMassFilePath = TConstruction->getRegionsMassDataPath();
    //OrganSourceName = TConstruction->getRegionsMassDataPath();
    //generatePoints = TConstruction->getRegionsMassDataPath();
}



void G4TVolumeBuilderUsingSTLGeomtry::ConstructBox(G4ThreeVector Center, G4double SizeX , G4double SizeY , G4double SizeZ){

    G4Box* box = new G4Box("world", SizeX, SizeY, SizeZ);
    G4double Z, A;
    A = 12.011*g/mole; G4Element* elC = new G4Element("Carbon","C",Z = 6.,A);
    A = 14.01*g/mole; G4Element* elN = new G4Element("Nitrogen","N",Z = 7.,A);
    A = 16.00*g/mole; G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A);

    G4double d = 1.290*mg/cm3;
    G4Material* AirMaterial = new G4Material("Air",d,2);
    AirMaterial->AddElement(elN,0.7); AirMaterial->AddElement(elO,0.3);

    LogicalBoxVolume = new G4LogicalVolume(box , AirMaterial , "Box" , 0 , 0 , 0);

    G4VisAttributes* boxVisAtt = new G4VisAttributes(G4Colour( 0.94 , 0.5 , 0.5 ) );
    boxVisAtt->SetForceSolid(true);
    LogicalBoxVolume->SetVisAttributes(boxVisAtt);

    G4RotationMatrix* BoxMatrixRot = new G4RotationMatrix();
    BoxMatrixRot->rotateX((G4ThreeVector(0. , 0. , 0.)).x()); BoxMatrixRot->rotateY((G4ThreeVector(0. , 0. , 0.)).y()); BoxMatrixRot->rotateZ((G4ThreeVector(0. , 0. , 0.)).z());

    PhysicalBoxVolume = new G4PVPlacement( BoxMatrixRot , Center , "physicalBox" , LogicalBoxVolume, MotherPhysVolume , false , 0 , false );
}
void G4TVolumeBuilderUsingSTLGeomtry::RemoveBoxVolumeAndSetOrgVol()
{

    // add the organVolume with data rot and pos of the box that was surrender it
    orgPhyVol = new G4PVPlacement( orgPhyVol->GetRotation() , PhysicalBoxVolume->GetTranslation() , orgPhyVol->GetLogicalVolume()->GetName() , orgPhyVol->GetLogicalVolume() , organMother , false , 0 , false );

    // now, remove the box with it's content from the phantom Volume
    G4LogicalVolume* motherLogical = PhysicalBoxVolume->GetMotherLogical();
    motherLogical->RemoveDaughter(PhysicalBoxVolume);

}



//called from ReadGeometryFile()
void G4TVolumeBuilderUsingSTLGeomtry::ReadOrganMassFile(G4String f ){

    std::map<G4String,std::vector<G4String>> DataTables;

    G4String line , organ, word;
    //G4double x,y,z;
    G4String WorldName;
    G4String namesOfMotherOrgan;
    std::vector<G4String> namesOfOrgans;

    std::ifstream fileR(f.c_str());

    if(fileR.is_open()){

        G4cout << "Reading "<< f << " ... " << G4endl ;

        while (getline(fileR, line)) {

            std::istringstream LineString(line);

            if(LineString.str().empty()){
                continue;
            }

            LineString >> organ;

            //G4cout << " the word " << word << G4endl ;
            if(organ.empty() || organ == ""){ //  || organ.isNull()
                continue;
            }

            //createdOrgans.push_back(organ);
            LineString >> createdOrganMass[organ];

            //LineString >> x >> y >> y;
            //createdPositions[organ] = G4ThreeVector(x,y,z);

            continue;
        }

        fileR.close();
    }
    else{

        G4cout << "cannot open the file " << f.c_str() << G4endl ;
    }
}



void G4TVolumeBuilderUsingSTLGeomtry::ReadSTLFile(G4String TEXT_GEO_File){



}


G4VPhysicalVolume* G4TVolumeBuilderUsingSTLGeomtry::ReadGeometryFile(G4String TEXT_GEO_File)
{

    std::ifstream ifile(TEXT_GEO_File.c_str());
    if (ifile) {

        ReadOrganMassFile(OrganMassFilePath);


        G4LogicalVolumeStore* nn = G4LogicalVolumeStore::GetInstance() ;

        G4tgbVolumeMgr* G4tgbVolumeMgrObj = G4tgbVolumeMgr::GetInstance();
        G4tgbVolumeMgrObj->AddTextFile(TEXT_GEO_File);
        WorldPhysicalVol = G4tgbVolumeMgrObj->ReadAndConstructDetector();
        MotherPhysVolume = G4tgbVolumeMgrObj->GetTopPhysVol();

        //G4cout << " \n\n\n\n\n\n\n\n\n\nNunmber of daughter volumes " << WorldPhyVolume->GetLogicalVolume()->GetNoDaughters() << G4endl;
        G4cout << " Defined phantom organs in file are : " << G4endl;

        for(G4int lk = 0 ; lk < (G4int)createdOrgans.size() ; lk++ ){ // the WorldVolume

            G4cout << createdOrgans[lk] << G4endl;
            nn->DeRegister(nn->GetVolume(createdOrgans[lk]));

            // Visualization Attributes
            G4Colour colour = G4Colour((G4double)G4UniformRand(), (G4double)G4UniformRand(), (G4double)G4UniformRand());
            G4VisAttributes* organVisAtt = new G4VisAttributes(colour);
            // Visualization Attributes
            organVisAtt->SetForceSolid(false);
            G4tgbVolumeMgrObj->FindG4LogVol(createdOrgans[lk])->SetVisAttributes(organVisAtt);


            if(G4tgbVolumeMgrObj->FindG4LogVol(createdOrgans[lk])){

                createdOrganDensity[createdOrgans[lk]] = G4tgbVolumeMgrObj->FindG4LogVol(createdOrgans[lk])->GetMaterial()->GetDensity()/(6.24151e+18); // density in g/cm3;
                createdPositions[createdOrgans[lk]] = G4tgbVolumeMgrObj->FindG4PhysVol(createdOrgans[lk])->GetTranslation();

                //G4cout << "createdPositions[createdOrgans[lk]] " << createdPositions[createdOrgans[lk]] << G4endl;

                G4RotationMatrix* rm = G4tgbVolumeMgrObj->FindG4PhysVol(createdOrgans[lk])->GetRotation();
                G4ThreeVector rot; rm->rotateX(rot.x()); rm->rotateY(rot.y()); rm->rotateZ(rot.z());
                createdRotations[createdOrgans[lk]] = rot;



                //G4cout << "generatePoints " << generatePoints << " OrganSourceName == createdOrgans[lk] " << OrganSourceName << G4endl;

                G4VPhysicalVolume* physVol = G4tgbVolumeMgrObj->FindG4PhysVol(createdOrgans[lk]) ;
                G4String worldname = G4tgbVolumeMgrObj->GetTopPhysVol()->GetLogicalVolume()->GetName();
                G4String PhyVolName = createdOrgans[lk] ;
                G4ThreeVector vectPos(0,0,0);

                //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n"<< G4endl;

                while ( PhyVolName != worldname )
                {
                    vectPos += physVol->GetObjectTranslation();
                    vectPos = vectPos.transform(physVol->GetObjectRotationValue());
                    //G4cout << "PhyVolName " << PhyVolName << " Pos " << physVol->GetObjectTranslation()<< " rot " << physVol->GetObjectRotationValue()<<  G4endl;
                    //G4cout << "vectPos " << vectPos << G4endl;

                    physVol = G4tgbVolumeMgrObj->FindG4PhysVol(physVol->GetMotherLogical()->GetName());
                    PhyVolName = physVol->GetLogicalVolume()->GetName();
                }

                if( generatePoints && OrganSourceName == createdOrgans[lk]){
                    boxCenterPos = vectPos;
                    G4cout << "boxCenterPos " << boxCenterPos<< G4endl;
                }

            }
        }

        G4String worldname = G4tgbVolumeMgrObj->GetTopPhysVol()->GetLogicalVolume()->GetName();
        nn->DeRegister(nn->GetVolume(worldname));

        createdOrgans.push_back(worldname);
        createdOrganMass[worldname] = G4tgbVolumeMgrObj->GetTopPhysVol()->GetLogicalVolume()->GetMass()/(6.24151e+24); // density in g/cm3;
        createdOrganDensity[worldname] = G4tgbVolumeMgrObj->GetTopPhysVol()->GetLogicalVolume()->GetMaterial()->GetDensity()/(6.24151e+18); // density in g/cm3;
        createdPositions[worldname] = G4tgbVolumeMgrObj->GetTopPhysVol()->GetTranslation();

        //G4RotationMatrix* rm = G4tgbVolumeMgrObj->GetTopPhysVol()->GetRotation();
        //G4ThreeVector rot; rm->rotateX(rot.x()); rm->rotateY(rot.y()); rm->rotateZ(rot.z());
        //createdRotations[worldname] = rot;

        //G4cout << " Geometry is readed " << worldname<< G4endl;
        return WorldPhysicalVol;

    }
    else{
        G4cout << " Unable to open file : " << TEXT_GEO_File.c_str() << G4endl;
        return 0;
    }

}



G4VPhysicalVolume* G4TVolumeBuilderUsingSTLGeomtry::ReadSeparatedGeometryFile(G4String TEXT_GEO_File)
{


}


