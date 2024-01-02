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
#include "G4TVolumeConstruction.hh"
#include "G4PhysicalConstants.hh"
#include "G4TPointDataGeneration.hh"

G4TPointDataGeneration::G4TPointDataGeneration(){}
G4TPointDataGeneration::~G4TPointDataGeneration(){}


// called from G4TVolumeConstruction
G4bool G4TPointDataGeneration::GenerateEventsPosition(){

    //G4cout << " SourceType : " << SourceType << G4endl;

    const G4TVolumeConstruction* VC1 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    if(FILE* file = fopen(PositionDataFile,"r")){

        if(VC1->getGeneratePostions() == "yes"){
            G4cout << " The file is already existed, it will be replaced by the new : " << PositionDataFile << G4endl;
        }
        else{
            G4cout << " The file is already existed (remove to recreate the Position data file ) : " << PositionDataFile << G4endl;
            fclose(file);
            return true;
        }
    }

    G4double posX , posY , posZ, X, Y, Z;
    G4int numbSucces = 0 ;

    std::ofstream PositionFileStream;
    PositionFileStream.open(PositionDataFile.c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary

    if(!PositionFileStream.is_open()) {
        G4cout << " Unable to open file OR is already existed : " << PositionDataFile << G4endl;
        return false;
    }
    else {

        if(SourceType == "Voxels"){

            unsigned int numbOfPointsToGenInVoxel ;

            unsigned int PointsNumInc = 0, numbSucces = 0;
            //G4double Xpos, Ypos, Zpos;

            if(UseDicomCumAct && VC1->getGeometryFileType() == "DICOM"){

                G4cout << "************** Generate " << NumberOfGenPointsToSave << " Event initial Position Data in all voxels (" << TotalVoxelsNumber <<  ") basing on the cumulated activity in each voxel - Saving the data to "<< PositionDataFile.c_str() << G4endl;

                double* CumulativeActivities ;//= VC1->getCumulativeActivities();

                for (unsigned int ii = 0 ; ii < TotalVoxelsNumber; ii++ ) {

                    numbOfPointsToGenInVoxel = CumulativeActivities[ii];
                    //G4cout << "the voxel : " << ii << " of Region " << CNRN[ii] << ", contains "<< numbOfPointsToGenInVoxel << " Particle " << G4endl;

                    numbSucces = 0;
                    while (numbSucces < numbOfPointsToGenInVoxel){

                        if(CopyNumberRegionNameMap[ii] != SourceRegionName){ continue;}

                        posX = (CopyNumberXPos[ii] - VoxelsSize.getX()) + (G4double)G4UniformRand() * (2*VoxelsSize.getX()); //generateRandom(hXmin,hXmax);
                        posY = (CopyNumberYPos[ii] - VoxelsSize.getY()) + (G4double)G4UniformRand() * (2*VoxelsSize.getY()); //generateRandom(hXmin,hXmax);
                        posZ = (CopyNumberZPos[ii] - VoxelsSize.getZ()) + (G4double)G4UniformRand() * (2*VoxelsSize.getZ()); //generateRandom(hXmin,hXmax);

                        /*
                        G4cout << " copyNumber : "  << RegionCopyNumberPositionMap[SourceRegionName][jh] << " in organ " <<  OrganIDOrganNameMapVoxel[CopyNumbersOrganIDMapVoxel[OrganNameCopyNumbersMap[SourceRegionName][jh]]] << " - in a position : " << SourceRegionNameRegionCopyNumberPositionMap[SourceRegionName][OrganNameCopyNumbersMap[SourceRegionName][jh]] << G4endl;
                        G4cout << " MinX=" << Xpos - hXmax << " < posX="<< posX << " < MaxX="<< Xpos+hXmax << G4endl;
                        G4cout << " MinY=" << Ypos - hYmax << " < posY="<< posY << " < MaxY="<< Ypos+hYmax << G4endl;
                        G4cout << " MinZ=" << Zpos - hZmax << " < posZ="<< posZ << " < MaxZ="<< Zpos+hZmax << G4endl;
                        */

                        PositionFileStream << posX << " " << posY << " " << posZ <<"\n" ;

                        numbSucces++ ;

                        //std::ostringstream text;
                        //text <<"\r"<< " " << PointsNumInc <<"/"<<TotalNumbOfPointsToGen << " ";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
                        //printf(text.str().c_str());

                        PointsNumInc++ ;
                    }
                }
            }
            else {

                G4cout << "************** Generate " << NumberOfGenPointsToSave << " Event initial Position Data in voxelized Region ("<< SourceRegionName << ") that contains " << VC1->getRegionNumberOfCNMap()[SourceRegionName] << " voxel - Saving the data to "<< PositionDataFile.c_str() << G4endl;


                if(SourceRegionName == "allregions"){

                    while (numbSucces < NumberOfGenPointsToSave){

                        for (unsigned int ii = 0 ; ii < TotalVoxelsNumber; ii++ ) {

                            bool Save = true;

                            for (int gg = 0 ; gg < SourceRegionsNamesToBeIgnoredValues.size() ; gg++) {
                                if(SourceRegionsNamesToBeIgnoredValues[gg] == CopyNumberRegionNameMap[ii]){
                                    Save=false;
                                    break;
                                }
                                //else{Save=true;}
                            }

                            if(Save==true){
                                //G4cout << " ii " << ii << " == "<< CNRN[ii] << " == "<< SourceRegionName << G4endl;

                                posX = (CopyNumberXPos[ii] - VoxelsSize.getX()) + (G4double)G4UniformRand() * (2*VoxelsSize.getX()); //generateRandom(hXmin,hXmax);
                                posY = (CopyNumberYPos[ii] - VoxelsSize.getY()) + (G4double)G4UniformRand() * (2*VoxelsSize.getY()); //generateRandom(hXmin,hXmax);
                                posZ = (CopyNumberZPos[ii] - VoxelsSize.getZ()) + (G4double)G4UniformRand() * (2*VoxelsSize.getZ()); //generateRandom(hXmin,hXmax);

                                PositionFileStream << posX << " " << posY << " " << posZ <<"\n" ;

                                numbSucces++ ;
                            }

                        }
                    }
                }else{
                    while (numbSucces < NumberOfGenPointsToSave){

                        //for(G4int jh = 0 ; jh < RegionCopyNumberPositionMap[SourceRegionName].size() ; jh++){
                        //for ( auto it = RegionCopyNumberPositionMap[SourceRegionName].begin(); it != RegionCopyNumberPositionMap[SourceRegionName].end(); ++it  ){
                        //Xpos = RegionCopyNumberPositionMap[SourceRegionName][it->first].getX();
                        //Ypos = RegionCopyNumberPositionMap[SourceRegionName][it->first].getY();
                        //Zpos = RegionCopyNumberPositionMap[SourceRegionName][it->first].getZ();
                        //G4cout << " Xpos " << Xpos << " < Ypos "<< Ypos << " < Zpos "<< Zpos << G4endl;

                        for (unsigned int ii = 0 ; ii < TotalVoxelsNumber; ii++ ) {

                            if(CopyNumberRegionNameMap[ii] != SourceRegionName){ continue;}

                            //G4cout << " ii " << ii << " == "<< CNRN[ii] << " == "<< SourceRegionName << G4endl;

                            posX = (CopyNumberXPos[ii] - VoxelsSize.getX()) + (G4double)G4UniformRand() * (2*VoxelsSize.getX()); //generateRandom(hXmin,hXmax);
                            posY = (CopyNumberYPos[ii] - VoxelsSize.getY()) + (G4double)G4UniformRand() * (2*VoxelsSize.getY()); //generateRandom(hXmin,hXmax);
                            posZ = (CopyNumberZPos[ii] - VoxelsSize.getZ()) + (G4double)G4UniformRand() * (2*VoxelsSize.getZ()); //generateRandom(hXmin,hXmax);

                            /*
                            G4cout << " copyNumber : "  << RegionCopyNumberPositionMap[SourceRegionName][jh] << " in organ " <<  OrganIDOrganNameMapVoxel[CopyNumbersOrganIDMapVoxel[OrganNameCopyNumbersMap[SourceRegionName][jh]]] << " - in a position : " << SourceRegionNameRegionCopyNumberPositionMap[SourceRegionName][OrganNameCopyNumbersMap[SourceRegionName][jh]] << G4endl;
                            G4cout << " MinX=" << Xpos - hXmax << " < posX="<< posX << " < MaxX="<< Xpos+hXmax << G4endl;
                            G4cout << " MinY=" << Ypos - hYmax << " < posY="<< posY << " < MaxY="<< Ypos+hYmax << G4endl;
                            G4cout << " MinZ=" << Zpos - hZmax << " < posZ="<< posZ << " < MaxZ="<< Zpos+hZmax << G4endl;
                            */

                            PositionFileStream << posX << " " << posY << " " << posZ <<"\n" ;

                            numbSucces++ ;

                            //std::ostringstream text;
                            //text <<"\r"<< " " << numbSucces <<"/"<<TotalNumbOfPointsToGen << " ";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
                            //printf(text.str().c_str());

                            //PointsNumInc++ ;
                        }
                    }
                }

            }
        }
        if(SourceType == "Volume"){

            G4cout << "************** Generate " << NumberOfGenPointsToSave << " Event initial Position Data in "<< SourceRegionName << " Of an absolute position "<< boxCenterPos << " mm - with Box Dim "<< BoxDimGene << " mm - Saving the data to "<< PositionDataFile.c_str() << G4endl;

            G4Navigator* aNavigator = new G4Navigator();
            aNavigator->SetWorldVolume(WorldPhysicalVolume);

            if(SourceRegionName == "allregions"){

                while (numbSucces < NumberOfGenPointsToSave){

                    posX = (boxCenterPos.getX() - BoxDimGene.getX()) + (G4double)G4UniformRand() * (2*BoxDimGene.getX()); //generateRandom(hXmin,BoxDimGene.getX());
                    posY = (boxCenterPos.getY() - BoxDimGene.getY()) + (G4double)G4UniformRand() * (2*BoxDimGene.getY()); //generateRandom(hXmin,BoxDimGene.getX());
                    posZ = (boxCenterPos.getZ() - BoxDimGene.getZ()) + (G4double)G4UniformRand() * (2*BoxDimGene.getZ()); //generateRandom(hXmin,BoxDimGene.getX());

                    bool Save = true;

                    for (int gg = 0 ; gg < SourceRegionsNamesToBeIgnoredValues.size() ; gg++) {
                        if(SourceRegionsNamesToBeIgnoredValues[gg] == aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(posX,posY,posZ))->GetLogicalVolume()->GetName()){
                            Save=false;
                            break;
                        }
                        //else{Save=true;}
                    }

                    while (Save==false){

                        Save = true;

                        posX = (boxCenterPos.getX() - BoxDimGene.getX()) + (G4double)G4UniformRand() * (2*BoxDimGene.getX()); //generateRandom(hXmin,BoxDimGene.getX());
                        posY = (boxCenterPos.getY() - BoxDimGene.getY()) + (G4double)G4UniformRand() * (2*BoxDimGene.getY()); //generateRandom(hXmin,BoxDimGene.getX());
                        posZ = (boxCenterPos.getZ() - BoxDimGene.getZ()) + (G4double)G4UniformRand() * (2*BoxDimGene.getZ()); //generateRandom(hXmin,BoxDimGene.getX());

                        //G4cout  << "Save "<< Save << " G4ThreeVector(posX,posY,posZ) " <<G4ThreeVector(posX,posY,posZ) << " SourceRegionName : " << SourceRegionName << "    Navigator find " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(posX,posY,posZ))->GetLogicalVolume()->GetName() << G4endl;

                        for (int gg = 0 ; gg < SourceRegionsNamesToBeIgnoredValues.size() ; gg++) {
                            if(SourceRegionsNamesToBeIgnoredValues[gg] == aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(posX,posY,posZ))->GetLogicalVolume()->GetName()){
                                Save=false;
                                break;
                            }
                            //else{Save=true;}
                        }
                    }

                    if (Save==true){

                        //G4cout  << "Save "<< Save << " G4ThreeVector(posX,posY,posZ) " <<G4ThreeVector(posX,posY,posZ) << " SourceRegionName : " << SourceRegionName << "    Navigator find " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(posX,posY,posZ))->GetLogicalVolume()->GetName() << G4endl;

                        PositionFileStream << posX << " " << posY << " " << posZ <<"\n" ;
                        numbSucces++ ;

                        //std::ostringstream text;
                        //text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " ";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
                        //printf(text.str().c_str());

                    }
                }

            }else{

                while (numbSucces < NumberOfGenPointsToSave){

                    posX = (boxCenterPos.getX() - BoxDimGene.getX()) + (G4double)G4UniformRand() * (2*BoxDimGene.getX()); //generateRandom(hXmin,BoxDimGene.getX());
                    posY = (boxCenterPos.getY() - BoxDimGene.getY()) + (G4double)G4UniformRand() * (2*BoxDimGene.getY()); //generateRandom(hXmin,BoxDimGene.getX());
                    posZ = (boxCenterPos.getZ() - BoxDimGene.getZ()) + (G4double)G4UniformRand() * (2*BoxDimGene.getZ()); //generateRandom(hXmin,BoxDimGene.getX());

                    //G4cout  << "position from the box  " << boxCenterPos << "  " << aNavigator->LocateGlobalPointAndSetup(boxCenterPos)->GetLogicalVolume()->GetName()<< G4endl;
                    //G4cout  << "position from the data " << organ3vector[SourceRegionName] << "  " << aNavigator->LocateGlobalPointAndSetup(organ3vector[SourceRegionName])->GetLogicalVolume()->GetName()<< G4endl;
                    //G4cout << " Zmin=" << (boxCenterPos.getZ() - 250) << " posZ=" << posZ << " Zmax=" << (boxCenterPos.getZ() + 250.) << G4endl;
                    //G4cout << " Ymin=" << (boxCenterPos.getY() - 250) << " posY=" << posY << " Ymax=" << (boxCenterPos.getY() + 0.) << G4endl;
                    //G4cout << " Xmin=" << (boxCenterPos.getX() - 250) << " posX=" << posX << " Xmax=" << (boxCenterPos.getX() + 0.) << G4endl;
                    //PointLogicalVolumeName = aNavigator->LocateGlobalPointAndSetup(pointVectorXYZ)->GetLogicalVolume()->GetName();
                    //G4VPhysicalVolume* physicalVolume = aNavigator->LocateGlobalPointAndSetup(pointVectorXYZ);
                    //G4cout  << posX << " " << posY << " " << posZ << " in " << PointLogicalVolumeName << G4endl;
                    //G4cout << " logical volume is : " << physicalVolume->GetLogicalVolume()->GetName() << "-- compared logical volume is : " << "HeartVolume" << G4endl;
                    //G4cout  << "SourceRegionName : " << SourceRegionName << "    compared with " << PointLogicalVolumeName << G4endl;

                    if (aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(posX,posY,posZ))->GetLogicalVolume()->GetName() == SourceRegionName){

                        PositionFileStream << posX << " " << posY << " " << posZ <<"\n" ;
                        numbSucces++ ;

                        //std::ostringstream text;
                        //text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " ";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
                        //printf(text.str().c_str());

                    }
                }
            }

            //G4cout << " The number of generated points in " << SourceRegionName << " is : "<< numbSucces << " , and saved in the file " << PositionDataFile.c_str() << G4endl;

        }
        else if(SourceType == "Solid"){

            G4cout << "************** Generate " << NumberOfGenPointsToSave << " Event Position Data in a " << SourceSolid << " Solid named " << SourceRegionName << " Of an absolute position "<< SourcePosition << " mm - Saving the data to "<< PositionDataFile.c_str() << G4endl;

            while (numbSucces < NumberOfGenPointsToSave){

                if(SourceSolid == "Sphere")
                {
                    X = Radius*2.;
                    Y = Radius*2.;
                    Z = Radius*2.;
                    while(((X*X)+(Y*Y)+(Z*Z)) > (Radius*Radius))
                    {
                        X = (G4UniformRand()*2.*Radius) - Radius;
                        Y = (G4UniformRand()*2.*Radius) - Radius;
                        Z = (G4UniformRand()*2.*Radius) - Radius;
                    }
                }
                else if(SourceSolid == "Ellipsoid")
                {
                    G4double temp = 100.;
                    while(temp > 1.)
                    {
                        X = (G4UniformRand()*2.*HalfX) - HalfX;
                        Y = (G4UniformRand()*2.*HalfY) - HalfY;
                        Y = (G4UniformRand()*2.*HalfZ) - HalfZ;

                        temp = ((X*X)/(HalfX*HalfX)) + ((Y*Y)/(HalfY*HalfY)) + ((Z*Z)/(HalfZ*HalfZ));
                    }
                }
                else if(SourceSolid == "Cylinder")
                {
                    X = Radius*2.;
                    Y = Radius*2.;
                    while(((X*X)+(Y*Y)) > (Radius*Radius))
                    {
                        X = (G4UniformRand()*2.*Radius) - Radius;
                        Y = (G4UniformRand()*2.*Radius) - Radius;
                        Z = (G4UniformRand()*2.*HalfZ) - HalfZ;
                    }
                }
                else if(SourceSolid == "EllipticCylinder")
                {
                    G4double expression = 20.;
                    while(expression > 1.)
                    {
                        X = (G4UniformRand()*2.*HalfX) - HalfX;
                        Y = (G4UniformRand()*2.*HalfY) - HalfY;
                        Z = (G4UniformRand()*2.*HalfZ) - HalfZ;

                        expression = ((X*X)/(HalfX*HalfX)) + ((Y*Y)/(HalfY*HalfY));
                    }
                }
                else if(SourceSolid == "Para")
                {
                    X = (G4UniformRand()*2.*HalfX) - HalfX;
                    Y = (G4UniformRand()*2.*HalfY) - HalfY;
                    Z = (G4UniformRand()*2.*HalfZ) - HalfZ;
                    //X = X + Z*std::tan(ParTheta)*std::cos(ParPhi) + y*std::tan(ParAlpha);
                    //y = y + z*std::tan(ParTheta)*std::sin(ParPhi);
                    //Z = Z;
                }

                X = X+SourcePosition.getX(), Y = Y+SourcePosition.getY() , Z = Z + SourcePosition.getZ();

                PositionFileStream << X << " " << Y << " " << Z <<"\n" ;
                numbSucces++ ;
                //std::ostringstream text;
                //text <<"\r"<< " " << numbSucces <<"/"<< NumberOfGenPointsToSave << " ";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
                //printf(text.str().c_str());
            }
        }
        else if(SourceType == "Surface"){

            G4cout << "************** Generate " << NumberOfGenPointsToSave << " Event Position Data in a " << SourceSurface << " Surface named " << SourceRegionName << " Of an absolute position "<< SourcePosition << " mm - Saving the data to "<< PositionDataFile.c_str() << G4endl;

            while (numbSucces < NumberOfGenPointsToSave){

                if(SourceSurface == "Sphere")
                {
                    G4double phi = twopi*G4UniformRand(), Theta = twopi*G4UniformRand();
                    X = Radius * std::sin(Theta) * std::cos(phi);
                    Y = Radius * std::sin(Theta) * std::sin(phi);
                    Z = Radius * std::cos(Theta);
                }
                else if(SourceSurface == "Ellipsoid"){

                    G4double minphi, maxphi, middlephi;
                    G4double answer, constant;

                    constant = pi/(HalfX*HalfX) + pi/(HalfY*HalfY) + twopi/(HalfZ*HalfZ);

                    // simplified approach
                    G4double theta = G4UniformRand();
                    G4double phi = G4UniformRand();

                    theta = std::acos(1. - 2.*theta);
                    minphi = 0.;
                    maxphi = twopi;
                    while(maxphi-minphi > 0.)
                    {
                        middlephi = (maxphi+minphi)/2.;
                        answer = (1./(HalfX*HalfX))*(middlephi/2. + std::sin(2*middlephi)/4.)
                                + (1./(HalfY*HalfY))*(middlephi/2. - std::sin(2*middlephi)/4.)
                                + middlephi/(HalfZ*HalfZ);
                        answer = answer/constant;
                        if(answer > phi) maxphi = middlephi;
                        if(answer < phi) minphi = middlephi;
                        if(std::fabs(answer-phi) <= 0.00001)
                        {
                            minphi = maxphi +1;
                            phi = middlephi;
                        }
                    }

                    X = std::sin(theta)*std::cos(phi);
                    Y = std::sin(theta)*std::sin(phi);
                    Z = std::cos(theta);
                }

                X = X+SourcePosition.getX(), Y = Y+SourcePosition.getY() , Z = Z + SourcePosition.getZ();

                PositionFileStream << X << " " << Y << " " << Z <<"\n" ;
                numbSucces++ ;
                //std::ostringstream text;
                //text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " ";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
                //printf(text.str().c_str());
            }
        }
        else if(SourceType == "Beam"){

            G4cout << "\n************** Generate " << NumberOfGenPointsToSave << " Event Position Data in a Beam named " << SourceRegionName << " Of an absolute position "<< SourcePosition << " mm and SDev " << BeamSDev <<" - Saving the data to "<< PositionDataFile.c_str() << G4endl;

            while (numbSucces < NumberOfGenPointsToSave){

                G4double A = G4RandGauss::shoot(0.0,BeamSDev);
                G4double B = G4RandGauss::shoot(0.0,BeamSDev);

                if(SourceAxis == "Z"){
                    X = A+SourcePosition.getX(), Y = B+SourcePosition.getY() , Z = SourcePosition.getZ();
                }
                else if (SourceAxis == "Y") {
                    X = A+SourcePosition.getX(), Z = B+SourcePosition.getZ() , Y = SourcePosition.getY();
                }
                else if (SourceAxis == "X") {
                    Y = A+SourcePosition.getY(), Z = B+SourcePosition.getZ() , X = SourcePosition.getX();
                }

                PositionFileStream << X << " " << Y << " " << Z <<"\n" ;
                numbSucces++ ;
                //std::ostringstream text;
                //text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " ";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
                //printf(text.str().c_str());

            }
        }
        else if(SourceType == "Plane"){

            G4cout << "************** Generate " << NumberOfGenPointsToSave << " Event Position Data in a " << SourcePlane << " plane named " << SourceRegionName << " Of an absolute position "<< SourcePosition << " mm - Saving the data to "<< PositionDataFile.c_str() << G4endl;

            G4double A, B, HalfA, HalfB;

            if(SourceAxis == "Z"){
                HalfA = HalfX, HalfB = HalfY;
            }
            else if (SourceAxis == "Y") {
                HalfA = HalfX, HalfB = HalfZ;
            }
            else if (SourceAxis == "X") {
                HalfA = HalfY, HalfB = HalfZ;
            }

            while (numbSucces < NumberOfGenPointsToSave){

                if(SourcePlane == "Circle")
                {
                    A = Radius + 100.;
                    B = Radius + 100.;
                    G4double Diameter = 2*Radius;
                    while(std::sqrt((A*A) + (B*B)) > Radius)
                    {
                        A = ((G4UniformRand()*Diameter) - Radius);
                        B = ((G4UniformRand()*Diameter) - Radius);
                    }
                }
                else if(SourcePlane == "Annulus")
                {
                    A = Radius + 100.;
                    B = Radius + 100.;
                    while(std::sqrt((A*A) + (B*B)) > Radius || std::sqrt((A*A) + (B*B)) < RadiusIn )
                    {
                        A = ((G4UniformRand()*2.*Radius) - Radius);
                        B = ((G4UniformRand()*2.*Radius) - Radius);
                    }
                }
                else if(SourcePlane == "Ellipse")
                {
                    G4double expression = 20.;
                    while(expression > 1.)
                    {
                        A = ((G4UniformRand()*2.*HalfA) - HalfA);
                        B = ((G4UniformRand()*2.*HalfB) - HalfB);
                        expression = ((A*A)/(HalfA*HalfA)) + ((B*B)/(HalfB*HalfB));
                    }
                }
                else if(SourcePlane == "Square")
                {
                    A = ((G4UniformRand()*2.*HalfA) - HalfA);
                    B = ((G4UniformRand()*2.*HalfA) - HalfA);
                }
                else if(SourcePlane == "Rectangle")
                {
                    A = ((G4UniformRand()*2.*HalfA) - HalfA);
                    B = ((G4UniformRand()*2.*HalfB) - HalfB);
                }

                if(SourceAxis == "Z"){
                    X = A+SourcePosition.getX(), Y = B+SourcePosition.getY() , Z = SourcePosition.getZ();
                }
                else if (SourceAxis == "Y") {
                    X = A+SourcePosition.getX(), Z = B+SourcePosition.getZ() , Y = SourcePosition.getY();
                }
                else if (SourceAxis == "X") {
                    Y = A+SourcePosition.getY(), Z = B+SourcePosition.getZ() , X = SourcePosition.getX();
                }


                PositionFileStream << X << " " << Y << " " << Z <<"\n" ;
                numbSucces++ ;
                //std::ostringstream text;
                //text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " ";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
                //printf(text.str().c_str());
            }
        }
        else if(SourceType == "Point"){

            while (numbSucces < NumberOfGenPointsToSave){
                PositionFileStream << SourcePosition.getX()<< " " << SourcePosition.getY() << " " << SourcePosition.getZ() <<"\n" ;
                numbSucces++ ;
                //std::ostringstream text;
                //text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " ";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
                //printf(text.str().c_str());
            }
        }
        else if(SourceType == "Rotated"){

            while (numbSucces < NumberOfGenPointsToSave){

                X = Radius + 100.;
                Y = Radius + 100.;
                Z = 0.;

                G4double Diameter = 2*Radius;
                while(std::sqrt((X*X) + (Y*Y)) > Radius)
                {
                    X = ((G4UniformRand()*Diameter) - Radius);
                    Y = ((G4UniformRand()*Diameter) - Radius);
                }


                // This takes in 2 vectors,1 in x' and 2 in the plane x'-y', and from these takes a cross product to calculate z'. Then a cross product is taken between x' and z' to give y'
                G4ThreeVector Rotz, RelPos, AbsPos;

                SourceRotVector1 = SourceRotVector1.unit(); // x'
                SourceRotVector2 = SourceRotVector2.unit(); // vector in x'y' plane
                Rotz = SourceRotVector1.cross(SourceRotVector2); // z'
                Rotz = Rotz.unit();
                SourceRotVector2 = Rotz.cross(SourceRotVector1); // y'
                SourceRotVector2 = SourceRotVector2.unit();

                RelPos.setX((X * SourceRotVector1.x()) + (Y * SourceRotVector2.x()) + (Z * Rotz.x()));
                RelPos.setY((X * SourceRotVector1.y()) + (Y * SourceRotVector2.y()) + (Z * Rotz.y()));
                RelPos.setZ((X * SourceRotVector1.z()) + (Y * SourceRotVector2.z()) + (Z * Rotz.z()));

                // RelPos = G4ThreeVector(X,Y,Z);
                // Translate
                AbsPos = SourcePosition + RelPos;

                PositionFileStream << AbsPos.getX()<< " " << AbsPos.getY() << " " << AbsPos.getZ() <<"\n" ;
                numbSucces++ ;
                //std::ostringstream text;
                //text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " ";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
                //printf(text.str().c_str());
            }
        }

        printf("\n");
        PositionFileStream.close();
        return true;
    }
}

// called from VolumeConstruction()
G4bool G4TPointDataGeneration::GenerateEventsEnergy(){

    const G4TVolumeConstruction* VC1 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    if(FILE* file = fopen(EnergyDataFile.c_str(),"r")){

        if(VC1->getGenerateEnergies() == "yes"){
            G4cout << " The file is already existed, it will be replaced by the new : " << EnergyDataFile << G4endl;
        }
        else{
            G4cout << " The file is already existed (remove to recreate the Energy data file ) : " << EnergyDataFile << G4endl;
            fclose(file);
            return true;
        }
    }

    G4cout << "************** Generate " << NumberOfGenPointsToSave << " Event initial Energy(MeV) and Saving the data to " << EnergyDataFile<< G4endl;

    G4int numbSucces = 0;
    std::ofstream EneFileStream;
    EneFileStream.open(EnergyDataFile.c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary

    if(EneFileStream.is_open()){

        //G4double pi = 4*std::arctan(1);  pi/4 = 4*arctan(1/5) - arctan(1/239) (Machin's formula) , or using towpi from G4PhysicalConstants.hh
        if(EnergyDistribution == "Mono"){
            while (numbSucces < NumberOfGenPointsToSave){
                EneFileStream << MonoEnergy <<" " ;
                numbSucces++ ;
                //std::ostringstream text; text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " "; printf(text.str().c_str());
            }
        }
        else if(EnergyDistribution == "Uniform"){
            while (numbSucces < NumberOfGenPointsToSave){
                G4double Energy = UniformEmin + (UniformEmax - UniformEmin) * (G4double)G4UniformRand() ;
                EneFileStream << Energy <<" " ;
                numbSucces++ ;
                //std::ostringstream text; text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " "; printf(text.str().c_str());
            }
        }
        else if(EnergyDistribution == "Rayleigh"){
            while (numbSucces < NumberOfGenPointsToSave){
                G4double mean = RayleighEmax/3. , sigma = mean * std::sqrt(1/twopi) ,Energy = sigma*std::sqrt(-2*std::log(1-(G4double)G4UniformRand()));
                EneFileStream << Energy <<" " ;
                numbSucces++ ;
                //std::ostringstream text; text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " "; printf(text.str().c_str());
            }
        }
        else if(EnergyDistribution == "Gauss"){
            while (numbSucces < NumberOfGenPointsToSave){
                EneFileStream << CLHEP::RandGauss::shoot(GaussMean, GaussSDev)  <<" " ;
                numbSucces++ ;
                //std::ostringstream text; text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " "; printf(text.str().c_str());
            }
        }
        printf("\n");
        EneFileStream.close();
        return true;
    }
    else {
        return true;
    }

}

// called from VolumeConstruction()
G4bool G4TPointDataGeneration::GenerateEventsMomentumDirection(){

    const G4TVolumeConstruction* VC1 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    if(FILE* file = fopen(MomDirDataFile,"r")){

        if(VC1->getGenerateMomDirs() == "yes"){
            G4cout << " The file is already existed, it will be replaced by the new : " << MomDirDataFile << G4endl;
        }
        else{
            G4cout << " The file is already existed (remove to recreate the MomDirs data file ) : " << MomDirDataFile << G4endl;
            fclose(file);
            return true;
        }
    }

    G4cout << "************** Generate " << NumberOfGenPointsToSave << " Event initial Momentum Direction, Dist " << MomDirDistribution << " and Saving the data to " << MomDirDataFile << G4endl;

    G4ThreeVector vectorDir = G4ThreeVector( 0. , 0. , 0. ) ;

    G4int numbSucces = 0;
    std::ofstream MomDirFileStream;
    MomDirFileStream.open(MomDirDataFile.c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary

    if(MomDirFileStream.is_open()){

        if(MomDirDistribution == "Isotropic"){
            //distribution uniform in solid angle
            while (numbSucces < NumberOfGenPointsToSave){
                G4double cosTheta = 2*G4UniformRand() - 1., phi = twopi*G4UniformRand();
                G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
                G4double ux = sinTheta*std::cos(phi), uy = sinTheta*std::sin(phi), uz = cosTheta;
                MomDirFileStream << ux << " " << uy << " " << uz <<"\n" ;
                numbSucces++ ;
                //std::ostringstream text; text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " "; printf(text.str().c_str());
            }
        }
        else if(MomDirDistribution == "Uniform"){
            while (numbSucces < NumberOfGenPointsToSave){
                vectorDir = G4ThreeVector( -1 + (G4double)G4UniformRand()*2 , -1 + (G4double)G4UniformRand()*2 , -1 + (G4double)G4UniformRand()*2 );
                MomDirFileStream << vectorDir.getX() << " " << vectorDir.getY() << " " << vectorDir.getZ() <<"\n" ;
                numbSucces++ ;
                //std::ostringstream text; text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " "; printf(text.str().c_str());
            }
        }
        else if(MomDirDistribution == "Directed"){
            while (numbSucces < NumberOfGenPointsToSave){
                G4double ux = std::sin(Theta)*std::cos(Phi), uy = std::sin(Theta)*std::sin(Phi), uz = std::cos(Theta);
                MomDirFileStream << ux << " " << uy << " " << uz <<"\n" ;
                numbSucces++ ;
                //std::ostringstream text; text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " "; printf(text.str().c_str());
            }
        }
        else if(MomDirDistribution == "Delimited"){

            while (numbSucces < NumberOfGenPointsToSave){

                //ThetaMin = ThetaMin/degree;
                //ThetaMax = ThetaMax/degree;
                //PhiMin = PhiMin/degree;
                //PhiMax = PhiMax/degree;

                G4double theta, phi;
                if(ThetaMin == ThetaMax){
                    theta = ThetaMax;
                }else {
                    theta = ThetaMin + G4UniformRand()*(ThetaMax-ThetaMin);
                }
                if(PhiMin == PhiMax){
                    phi = PhiMax;
                }else {
                    phi = PhiMin + G4UniformRand()*(PhiMax-PhiMin);
                }

                G4double ux = std::sin(theta)*std::cos(phi);
                G4double uy = std::sin(theta)*std::sin(phi);
                G4double uz = std::cos(theta);

                //G4cout << " ThetaMin " << ThetaMin  << " ThetaMax " << ThetaMax << " ThetaMin Degree " << ThetaMin/degree  << " ThetaMax Degree" << ThetaMax/degree  << " " << ux << " " << uy << " " << uz << G4endl;

                MomDirFileStream << ux << " " << uy << " " << uz <<"\n" ;
                numbSucces++ ;
                //std::ostringstream text; text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " "; printf(text.str().c_str());
            }
        }

        printf("\n");
        MomDirFileStream.close();
        return true;
    }
    else {
        return false;
    }
}

// called from VolumeConstruction()
void G4TPointDataGeneration::SetSourceInputs(){

    const G4TVolumeConstruction* TConstruction = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    WorldPhysicalVolume = TConstruction->getWorldPhyVolume();

    TotalVoxelsNumber = VoxXNumber * VoxYNumber * VoxZNumber;
    VoxelsSize = G4ThreeVector(VoxXHalfSize,VoxYHalfSize,VoxZHalfSize);

    /*
    G4ThreeVector VoxelsSize = TConstruction->getVoxelsSize();
    G4int TotalVoxelsNumber = TConstruction->getVoxXNumber() * TConstruction->getVoxYNumber() * TConstruction->getVoxZNumber();
    UseDicomCumAct = TConstruction->getUseDicomCumAct();
    SourceType = TConstruction->getSourceType();
    boxCenterPos = TConstruction->getboxCenterPos();
    SourcePosition = TConstruction->getSourcePosition();
    SourceRotVector1 = TConstruction->getSourceRotVector1();
    SourceRotVector2 = TConstruction->getSourceRotVector2();

    SourceSolid = TConstruction->getSourceSolid();
    SourceSurface = TConstruction->getSourceSurface();
    SourcePlane = TConstruction->getSourcePlane();
    SourceAxis = TConstruction->getSourceAxis();
    SourceRotation = TConstruction->getSourceRotation();
    Radius = TConstruction->getRadius();
    HalfX = TConstruction->getHalfX();
    HalfY = TConstruction->getHalfY();
    HalfZ = TConstruction->getHalfZ();

    RadiusIn = TConstruction->getRadiusIn();
    BeamSDev = TConstruction->getBeamSDev();

    PhiMax = TConstruction->getPhiMax();
    PhiMin = TConstruction->getPhiMin();
    ThetaMax = TConstruction->getThetaMax();
    ThetaMin = TConstruction->getThetaMin();

    BoxDimGene = TConstruction->getBoxDimGene();

    NumberOfGenPointsToSave = TConstruction->getPointNumberToSave();
    Theta = TConstruction->getTheta();
    Phi = TConstruction->getPhi();
    particleName = TConstruction->getParticleName();
    MonoEnergy = TConstruction->getMonoEnergy();
    MomDirDistribution = TConstruction->getMomDirDistribution();
    EnergyDistribution = TConstruction->getEnergyDistribution();
    SourceRegionName = TConstruction->getOrganSource();
    GaussSDev = TConstruction->getGaussSDev();;
    GaussMean = TConstruction->getGaussMean();
    UniformEmin = TConstruction->getUniformEmin();
    UniformEmax = TConstruction->getUniformEmax();
    RayleighEmax = TConstruction->getRayleighEmax();

    PositionDataFile = TConstruction->getPositionDataFile();
    EnergyDataFile = TConstruction->getEnergyDataFile();
    MomDirDataFile = TConstruction->getMomDirDataFile();
*/

}

