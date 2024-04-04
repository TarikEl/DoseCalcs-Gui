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

#include "G4TDirectPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4TVolumeConstruction.hh"

#include "G4RunManager.hh"
//#include "G4ios.hh"
#include <iostream>
#include "G4PhysicalConstants.hh"

G4TDirectPrimaryGeneratorAction::G4TDirectPrimaryGeneratorAction(){

    //std::cout << "dddddddddddddddddddddddd \n\n\n\n\n\n\n\n\n\n\n\n************** The primary generator Action initialization... "<< std::endl;

    particleGun = new G4ParticleGun(1);
    GunInitialize();
}

G4TDirectPrimaryGeneratorAction::~G4TDirectPrimaryGeneratorAction()
{
    //G4MUTEXDESTROY(mutex);
    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;

    delete particleGun;
}

void G4TDirectPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //particleGun->SetParticlePosition(G4ThreeVector(0, 0, 0));
    //particleGun->SetParticleEnergy(1.);
    //particleGun->SetParticleMomentumDirection(G4ParticleMomentum(0, 0, 1.));

    if(EvInc == 0){
        //std::cout << "begin of SourceInitialization()" << std::endl;
        SourceInitialization();
        //std::cout << "End of SourceInitialization()" << std::endl;

        particleGun->SetNumberOfParticles(1);

        if(EnergyTypeNum != 5){particleGun->SetParticleDefinition(particleDefinition);}

        EvInc++;
    }

    if(EnergyTypeNum == 5){particleGun->SetParticleDefinition(particleDefinitionList[ParNameList[EnergyListInc]]);}

    //GenerateEventsParticle();
    GenerateEventsEnergy();
    GenerateEventsPosition();
    GenerateEventsMomentumDirection();

    particleGun->SetParticlePosition(G4ThreeVector(X, Y, Z));
    particleGun->SetParticleEnergy(ENERGY);
    particleGun->SetParticleMomentumDirection(G4ParticleMomentum(XMOMD, YMOMD, ZMOMD));

    //std::cout << EvInc << " ENERGY=" << ENERGY << " X=" << X << " Y=" << Y << " Z=" << Z << " XMOMD=" << XMOMD << " YMOMD=" << YMOMD << " ZMOMD=" << ZMOMD << " " << NewRankSourceRegionsBoxDimValues[DataID] << " " << particleDefinition->GetParticleName() << std::endl;

    //if(WriteSourceDataToFiles == 1){SaveGeneratedDataToFiles();}

    particleGun->GeneratePrimaryVertex(anEvent);
    //std::cout << " particleGun->GeneratePrimaryVertex(anEvent) " << std::endl;


}

void G4TDirectPrimaryGeneratorAction::GenerateEventsPosition(){

    //std::cout << "G4TDirectPrimaryGeneratorAction::GenerateEventsPosition()" << std::endl;

    if(PositionTypeNum == 1){

        X = (NewRankSourceRegionsPosValues[DataID].getX() - NewRankSourceRegionsBoxDimValues[DataID].getX()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getX()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
        Y = (NewRankSourceRegionsPosValues[DataID].getY() - NewRankSourceRegionsBoxDimValues[DataID].getY()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getY()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
        Z = (NewRankSourceRegionsPosValues[DataID].getZ() - NewRankSourceRegionsBoxDimValues[DataID].getZ()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getZ()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());

        //std::cout << "DataID "<< DataID << " SourceType=" << SourceType << " NewRankSourceRegionsNamesValues[DataID]=" << NewRankSourceRegionsNamesValues[DataID] << " " << G4ThreeVector(X,Y,Z) << " is in " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() << std::endl;

        if(NewRankSourceRegionsNamesValues[DataID] == "allregions"){

            bool Save = true;

            for (int gg = 0 ; gg < SourceRegionsNamesToBeIgnoredValues.size() ; gg++) {
                if(SourceRegionsNamesToBeIgnoredValues[gg] == aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName()){
                    Save=false;
                    break;
                }
                //else{Save=true;}
            }

            //std::cout << "\n\n\n\n\n\nDataID "<< DataID << " SourceType=" << SourceType << " NewRankSourceRegionsNamesValues[DataID]=" << NewRankSourceRegionsNamesValues[DataID] << " " << G4ThreeVector(X,Y,Z) << " is in " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() << std::endl;

            while (Save==false){

                Save = true;

                X = (NewRankSourceRegionsPosValues[DataID].getX() - NewRankSourceRegionsBoxDimValues[DataID].getX()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getX()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
                Y = (NewRankSourceRegionsPosValues[DataID].getY() - NewRankSourceRegionsBoxDimValues[DataID].getY()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getY()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
                Z = (NewRankSourceRegionsPosValues[DataID].getZ() - NewRankSourceRegionsBoxDimValues[DataID].getZ()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getZ()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());

                for (int gg = 0 ; gg < SourceRegionsNamesToBeIgnoredValues.size() ; gg++) {
                    if(SourceRegionsNamesToBeIgnoredValues[gg] == aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName()){
                        Save=false;
                    }
                }

                //std::cout << "while:  DataID "<< DataID << " SourceType=" << SourceType << " NewRankSourceRegionsNamesValues[DataID]=" << NewRankSourceRegionsNamesValues[DataID] << " " << G4ThreeVector(X,Y,Z) << " is in " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() << std::endl;

            }

        }else{
            while (aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() != NewRankSourceRegionsNamesValues[DataID]){

                //std::cout << "DataID "<< DataID << " SourceType=" << SourceType << " NewRankSourceRegionsBoxDimValues[DataID]=" << NewRankSourceRegionsBoxDimValues[DataID] << " NewRankSourceRegionsNamesValues[DataID]=" << NewRankSourceRegionsNamesValues[DataID] << " NewRankSourceRegionsPosValues[DataID]=" << NewRankSourceRegionsPosValues[DataID] << " Generated Pos: " << G4ThreeVector(X,Y,Z) << " is in " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() << std::endl;

                X = (NewRankSourceRegionsPosValues[DataID].getX() - NewRankSourceRegionsBoxDimValues[DataID].getX()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getX()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
                Y = (NewRankSourceRegionsPosValues[DataID].getY() - NewRankSourceRegionsBoxDimValues[DataID].getY()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getY()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
                Z = (NewRankSourceRegionsPosValues[DataID].getZ() - NewRankSourceRegionsBoxDimValues[DataID].getZ()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getZ()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
            }
        }
    }
    else if(PositionTypeNum == 2){

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
        //return G4ThreeVector(X+SourcePosition.getX(), Y+SourcePosition.getY(), Z + SourcePosition.getZ());
    }
    else if(PositionTypeNum == 3){

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
        else if(SourceSurface == "Cylinder"){

            G4double phi = G4UniformRand() * twopi; // Random angle
            //G4double theta = G4UniformRand() * twopi; // Random angle
            X = HalfX * cos(phi);
            Y = HalfX * sin(phi);
            Z = ((G4UniformRand()*2.*HalfY) - HalfY); // Random z-coordinate within the height of the cylinder
        }

        X = X+SourcePosition.getX(), Y = Y+SourcePosition.getY() , Z = Z + SourcePosition.getZ();
        //return G4ThreeVector(X+SourcePosition.getX(), Y+SourcePosition.getY(), Z + SourcePosition.getZ());

    }
    else if(PositionTypeNum == 4){

        G4double A = G4RandGauss::shoot(0.0,BeamSDev);
        G4double B = G4RandGauss::shoot(0.0,BeamSDev);

        //std::cout << " BeamSDev " << BeamSDev << " A "<< A  << " B "<< B << " SourceAxis "<< SourceAxis << std::endl;

        if(SourceAxis == "Z"){
            X = A+SourcePosition.getX(), Y = B+SourcePosition.getY() , Z = SourcePosition.getZ();
        }
        else if (SourceAxis == "Y") {
            X = A+SourcePosition.getX(), Z = B+SourcePosition.getZ() , Y = SourcePosition.getY();
        }
        else if (SourceAxis == "X") {
            Y = A+SourcePosition.getY(), Z = B+SourcePosition.getZ() , X = SourcePosition.getX();
        }
        //return G4ThreeVector(X, Y, Z);

    }
    else if(PositionTypeNum == 5){

        G4double A, B, HalfA, HalfB;

        //if(SourceAxis == "Z"){
        //    HalfA = HalfX, HalfB = HalfY;
        //}
        //else if (SourceAxis == "Y") {
        //    HalfA = HalfX, HalfB = HalfZ;
        //}
        //else if (SourceAxis == "X") {
        //    HalfA = HalfY, HalfB = HalfZ;
        //}

        HalfA = HalfX, HalfB = HalfY;

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
            X = A+SourcePosition.getX(), Y = B+SourcePosition.getY(), Z = SourcePosition.getZ();
        }
        else if (SourceAxis == "Y") {
            Z = A+SourcePosition.getZ(), X = B+SourcePosition.getX(), Y = SourcePosition.getY();
        }
        else if (SourceAxis == "X") {
            Y = A+SourcePosition.getY(), Z = B+SourcePosition.getZ(), X = SourcePosition.getX();
        }
        //return G4ThreeVector(X, Y, Z);


        //// Rotation angle around Z-axis in degrees
        //double r = sqrt(X*X + Y*Y + Z*Z);
        //double theta1 = acos(Z / r); // inclination angle
        //double phi1 = atan2(Y, X); // azimuthal angle

        //double Theta2 = 15*degree;
        //double Phi2   = phi1;

        //X = r * sin(theta1) * cos(phi1+0*degree);
        //Y = r * sin(theta1) * sin(phi1+0*degree);
        //Z = r * cos(theta1+300*degree);



        //// Rotation angle around Y-axis in degrees
        //G4double angle = 135*degree;









        //double radius = 100.0;

        //// Define the range of angles for the plane
        //double minTheta = 0.0;  // Minimum inclination angle (from the z-axis)
        //double maxTheta = M_PI / 2.0;  // Maximum inclination angle (from the z-axis)
        //double minPhi = 0.0;  // Minimum azimuthal angle (around the z-axis)
        //double maxPhi = 2.0 * M_PI;  // Maximum azimuthal angle (around the z-axis)

        //// Generate a random angle within the specified range

        //double theta1 = minTheta + G4UniformRand()*(maxTheta-minTheta);
        //double phi1 = minPhi + G4UniformRand()*(maxPhi-minPhi);

        //// Convert spherical coordinates to Cartesian coordinates
        //double x, y, z;
        //X = radius * sin(theta1) * cos(phi1);
        //Y = radius * sin(theta1) * sin(phi1);
        //Z = radius * cos(theta1);






    }
    else if(PositionTypeNum == 6){
        X = SourcePosition.getX(), Y = SourcePosition.getY(), Z = SourcePosition.getZ();
        //return G4ThreeVector(SourcePosition.getX(), SourcePosition.getY(), SourcePosition.getZ());
    }
    else if(PositionTypeNum == 7){

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

        X = AbsPos.getX(), Y = AbsPos.getY(), Z = AbsPos.getZ() ;
        //return G4ThreeVector(AbsPos.getX(), AbsPos.getY(), AbsPos.getZ());


    }
    else if(PositionTypeNum == 8){

        G4bool insideChk(false);
        G4ThreeVector pos;
        do{
            pos = G4ThreeVector(NewRankTETBoxMinOfSourceRegion[DataID].getX()+NewRankTETBoxDimOfSourceRegion[DataID].getX()*G4UniformRand(),
                                NewRankTETBoxMinOfSourceRegion[DataID].getY()+NewRankTETBoxDimOfSourceRegion[DataID].getY()*G4UniformRand(),
                                NewRankTETBoxMinOfSourceRegion[DataID].getZ()+NewRankTETBoxDimOfSourceRegion[DataID].getZ()*G4UniformRand());

            for(auto tet:NewRankTETOfSourceRegion[DataID]){
                if(tet->Inside(pos) == kOutside) continue;
                insideChk = true;
                break;
            }
        }while(!insideChk);

        X = pos.getX();
        Y = pos.getY();
        Z = pos.getZ();

    }
    else if(PositionTypeNum == 9){

        G4bool insideChk(false);
        G4ThreeVector pos;
        do{
            pos = G4ThreeVector(NewRankTETBoxMinOfSourceRegion[DataID].getX()+NewRankTETBoxDimOfSourceRegion[DataID].getX()*G4UniformRand(),
                                NewRankTETBoxMinOfSourceRegion[DataID].getY()+NewRankTETBoxDimOfSourceRegion[DataID].getY()*G4UniformRand(),
                                NewRankTETBoxMinOfSourceRegion[DataID].getZ()+NewRankTETBoxDimOfSourceRegion[DataID].getZ()*G4UniformRand());

            for(auto tet:NewRankTETOfSourceRegion[DataID]){
                if(tet->Inside(pos) == kOutside) continue;
                insideChk = true;
                break;
            }
        }while(!insideChk);

        X = pos.getX();
        Y = pos.getY();
        Z = pos.getZ();
    }

    if(!RotPosAxis.empty() && RotTheta != 0.){
        RotatePosition();
    }

    particleGun->SetParticlePosition(G4ThreeVector(X, Y, Z));
    //return G4ThreeVector(X, Y, Z);
}
