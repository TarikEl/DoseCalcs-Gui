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

#include "G4TModifiedPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4TVolumeConstruction.hh"

#include "G4RunManager.hh"
//#include "G4ios.hh"
#include <iostream>
#include "G4PhysicalConstants.hh"

G4TModifiedPrimaryGeneratorAction::G4TModifiedPrimaryGeneratorAction(){

    //std::cout << " \n\n\n\n\n\n\n\n\n\n\n\n************** The primary generator Action initialization... "<< std::endl;

    particleGun = new G4ParticleGun(1);

    GunInitialize();

    particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(NewRankSourceParticlesNamesValues[DataID]);
    if(particleDefinition == nullptr){ G4String msg = "Particle name (" + NewRankSourceParticlesNamesValues[DataID] + ") in rank/thread (" + std::to_string(DataID) + ") not found "; G4Exception("Source Data", "1", FatalErrorInArgument, msg.c_str());}
    particleGun->SetParticleDefinition(particleDefinition);
    particleGun->SetNumberOfParticles(1);

    /*                          //// Save Data To file

        OpenFilesToSaveGeneratedData();
    */
}

G4TModifiedPrimaryGeneratorAction::~G4TModifiedPrimaryGeneratorAction()
{
    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;

    delete particleGun;
}

void G4TModifiedPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //particleGun->SetParticlePosition(G4ThreeVector(0, 0, 0));
    //particleGun->SetParticleEnergy(1.);
    //particleGun->SetParticleMomentumDirection(G4ParticleMomentum(0, 0, 1.));


    /*                          //// For Read Primary Generator, it is already separated from others Primary Generators
                                     Then, it is not included in the modified Primary Generator class
    */


    ////////////////////////////////////////////////// Initial Particle ///////////////////////////////////////////////////////

    /*                          //// RadioNuclide or File Particle

    particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(ParNameList[EnergyListInc]));
    // For other Sources Particle, The particles are set in the constructor

    */

    ////////////////////////////////////////////////// Initial Energy ///////////////////////////////////////////////////////

    /*                          //// mono energy distribution

        ENERGY = NewRankSourceEnergiesValues[DataID] ;
        TotalEmittedEnergy += ENERGY ;
    */

    /*                          //// uniform energy distribution

        ENERGY = UniformEmin + (NewRankSourceEnergiesValues[DataID] - UniformEmin) * (G4double)G4UniformRand() ;
        TotalEmittedEnergy += ENERGY ;
    */

    /*                          //// rayleigh energy distribution

        G4double mean = NewRankSourceEnergiesValues[DataID]/3. ;
        G4double sigma = mean * std::sqrt(1/twopi) ;
        ENERGY = sigma*std::sqrt(-2*std::log(1-(G4double)G4UniformRand()));
        TotalEmittedEnergy += ENERGY ;
    */

    /*                          //// gauss distribution

        ENERGY = CLHEP::RandGauss::shoot(NewRankSourceEnergiesValues[DataID], GaussSDev) ;
        TotalEmittedEnergy += ENERGY ;
    */

    /*                          //// Spectrum distribution
        //if(EvInc == 0){
        //    SourceInitialization();
        //    EvInc++;
        //}
        if(EnergyListInc == EventsNumPerThreadRank ){
            EnergyListInc= 0;
        }
        ENERGY = EnergyList[EnergyListInc] ;
        EnergyListInc++;
        TotalEmittedEnergy += ENERGY ;
    */

    /*                          //// distribution data from file
        //if(EvInc == 0){
        //    SourceInitialization();
        //    EvInc++;
        //}
        if(EnergyListInc == EventsNumPerThreadRank ){
            EnergyListInc= 0;
        }
        ENERGY = EnergyList[EnergyListInc] ;
        EnergyListInc++;
        TotalEmittedEnergy += ENERGY ;
    */


    ////////////////////////////////////////////////// Initial Position ///////////////////////////////////////////////////////

    /*                          //// stylized volumes
        //if(EvInc == 0){
        //    SourceInitialization();
        //    EvInc++;
        //}

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

                //std::cout << "DataID "<< DataID << " SourceType=" << SourceType << " NewRankSourceRegionsNamesValues[DataID]=" << NewRankSourceRegionsNamesValues[DataID] << " " << G4ThreeVector(X,Y,Z) << " is in " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() << std::endl;

            }

        }else{
            while (aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() != NewRankSourceRegionsNamesValues[DataID]){

                //std::cout << "DataID "<< DataID << " SourceType=" << SourceType << " NewRankSourceRegionsBoxDimValues[DataID]=" << NewRankSourceRegionsBoxDimValues[DataID] << " NewRankSourceRegionsNamesValues[DataID]=" << NewRankSourceRegionsNamesValues[DataID] << " NewRankSourceRegionsPosValues[DataID]=" << NewRankSourceRegionsPosValues[DataID] << " Generated Pos: " << G4ThreeVector(X,Y,Z) << " is in " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() << std::endl;

                X = (NewRankSourceRegionsPosValues[DataID].getX() - NewRankSourceRegionsBoxDimValues[DataID].getX()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getX()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
                Y = (NewRankSourceRegionsPosValues[DataID].getY() - NewRankSourceRegionsBoxDimValues[DataID].getY()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getY()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
                Z = (NewRankSourceRegionsPosValues[DataID].getZ() - NewRankSourceRegionsBoxDimValues[DataID].getZ()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getZ()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
            }
        }
    */

    /*                          //// solid Sphere

        X = Radius*2.;
        Y = Radius*2.;
        Z = Radius*2.;
        while(((X*X)+(Y*Y)+(Z*Z)) > (Radius*Radius))
        {
            X = (G4UniformRand()*2.*Radius) - Radius;
            Y = (G4UniformRand()*2.*Radius) - Radius;
            Z = (G4UniformRand()*2.*Radius) - Radius;
        }
        X = X+SourcePosition.getX(), Y = Y+SourcePosition.getY() , Z = Z + SourcePosition.getZ();
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// solid Ellipsoid

        G4double temp = 100.;
        while(temp > 1.)
        {
            X = (G4UniformRand()*2.*HalfX) - HalfX;
            Y = (G4UniformRand()*2.*HalfY) - HalfY;
            Y = (G4UniformRand()*2.*HalfZ) - HalfZ;

            temp = ((X*X)/(HalfX*HalfX)) + ((Y*Y)/(HalfY*HalfY)) + ((Z*Z)/(HalfZ*HalfZ));
        }
        X = X+SourcePosition.getX(), Y = Y+SourcePosition.getY() , Z = Z + SourcePosition.getZ();
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// solid Cylinder

        X = Radius*2.;
        Y = Radius*2.;
        while(((X*X)+(Y*Y)) > (Radius*Radius))
        {
            X = (G4UniformRand()*2.*Radius) - Radius;
            Y = (G4UniformRand()*2.*Radius) - Radius;
            Z = (G4UniformRand()*2.*HalfZ) - HalfZ;
        }
        X = X+SourcePosition.getX(), Y = Y+SourcePosition.getY() , Z = Z + SourcePosition.getZ();
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// solid EllipticCylinder

        G4double expression = 20.;
        while(expression > 1.)
        {
            X = (G4UniformRand()*2.*HalfX) - HalfX;
            Y = (G4UniformRand()*2.*HalfY) - HalfY;
            Z = (G4UniformRand()*2.*HalfZ) - HalfZ;

            expression = ((X*X)/(HalfX*HalfX)) + ((Y*Y)/(HalfY*HalfY));
        }
        X = X+SourcePosition.getX(), Y = Y+SourcePosition.getY() , Z = Z + SourcePosition.getZ();
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// solid Para

        X = (G4UniformRand()*2.*HalfX) - HalfX;
        Y = (G4UniformRand()*2.*HalfY) - HalfY;
        Z = (G4UniformRand()*2.*HalfZ) - HalfZ;
        //X = X + Z*std::tan(ParTheta)*std::cos(ParPhi) + y*std::tan(ParAlpha);
        //y = y + z*std::tan(ParTheta)*std::sin(ParPhi);
        //Z = Z;
    }
    X = X+SourcePosition.getX(), Y = Y+SourcePosition.getY() , Z = Z + SourcePosition.getZ();
    if(!RotPosAxis.empty() && RotTheta != 0.){
        RotatePosition();
    }
    */

    /*                          //// solid surface Sphere

        G4double phi = twopi*G4UniformRand(), Theta = twopi*G4UniformRand();
        X = Radius * std::sin(Theta) * std::cos(phi);
        Y = Radius * std::sin(Theta) * std::sin(phi);
        Z = Radius * std::cos(Theta);
        X = X+SourcePosition.getX(), Y = Y+SourcePosition.getY() , Z = Z + SourcePosition.getZ();
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// solid surface Ellipsoid

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
        X = X+SourcePosition.getX(), Y = Y+SourcePosition.getY() , Z = Z + SourcePosition.getZ();
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// solid surface Cylinder

        G4double phi = G4UniformRand() * twopi; // Random angle
        //G4double theta = G4UniformRand() * twopi; // Random angle
        X = HalfX * cos(phi);
        Y = HalfX * sin(phi);
        Z = ((G4UniformRand()*2.*HalfY) - HalfY); // Random z-coordinate within the height of the cylinder
        X = X+SourcePosition.getX(), Y = Y+SourcePosition.getY() , Z = Z + SourcePosition.getZ();
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// beam

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
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// Plane Circle

        G4double A, B, HalfA, HalfB;
        HalfA = HalfX, HalfB = HalfY;

        A = Radius + 100.;
        B = Radius + 100.;
        G4double Diameter = 2*Radius;
        while(std::sqrt((A*A) + (B*B)) > Radius)
        {
            A = ((G4UniformRand()*Diameter) - Radius);
            B = ((G4UniformRand()*Diameter) - Radius);
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
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// Plane Annulus

        G4double A, B, HalfA, HalfB;
        HalfA = HalfX, HalfB = HalfY;
        A = Radius + 100.;
        B = Radius + 100.;
        while(std::sqrt((A*A) + (B*B)) > Radius || std::sqrt((A*A) + (B*B)) < RadiusIn )
        {
            A = ((G4UniformRand()*2.*Radius) - Radius);
            B = ((G4UniformRand()*2.*Radius) - Radius);
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
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// Plane Ellipse

        G4double A, B, HalfA, HalfB;
        HalfA = HalfX, HalfB = HalfY;
        G4double expression = 20.;
        while(expression > 1.)
        {
            A = ((G4UniformRand()*2.*HalfA) - HalfA);
            B = ((G4UniformRand()*2.*HalfB) - HalfB);
            expression = ((A*A)/(HalfA*HalfA)) + ((B*B)/(HalfB*HalfB));
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
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// Plane Square

        G4double A, B, HalfA, HalfB;
        HalfA = HalfX, HalfB = HalfY;
        A = ((G4UniformRand()*2.*HalfA) - HalfA);
        B = ((G4UniformRand()*2.*HalfA) - HalfA);
        if(SourceAxis == "Z"){
            X = A+SourcePosition.getX(), Y = B+SourcePosition.getY(), Z = SourcePosition.getZ();
        }
        else if (SourceAxis == "Y") {
            Z = A+SourcePosition.getZ(), X = B+SourcePosition.getX(), Y = SourcePosition.getY();
        }
        else if (SourceAxis == "X") {
            Y = A+SourcePosition.getY(), Z = B+SourcePosition.getZ(), X = SourcePosition.getX();
        }
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// Plane Rectangle

        G4double A, B, HalfA, HalfB;
        HalfA = HalfX, HalfB = HalfY;
        A = ((G4UniformRand()*2.*HalfA) - HalfA);
        B = ((G4UniformRand()*2.*HalfB) - HalfB);
        if(SourceAxis == "Z"){
            X = A+SourcePosition.getX(), Y = B+SourcePosition.getY(), Z = SourcePosition.getZ();
        }
        else if (SourceAxis == "Y") {
            Z = A+SourcePosition.getZ(), X = B+SourcePosition.getX(), Y = SourcePosition.getY();
        }
        else if (SourceAxis == "X") {
            Y = A+SourcePosition.getY(), Z = B+SourcePosition.getZ(), X = SourcePosition.getX();
        }
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// point

        X = SourcePosition.getX(), Y = SourcePosition.getY(), Z = SourcePosition.getZ();
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// rotated form

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
        if(!RotPosAxis.empty() && RotTheta != 0.){
            RotatePosition();
        }
    */

    /*                          //// Tetrahedral geometry

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

    */

    ////////////////////////////////////////////////// Initial Momentum Direction ///////////////////////////////////////////////////////

    /*                          //// isotropic

        G4double cosTheta = 2*G4UniformRand() - 1., phi = twopi*G4UniformRand();
        G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
        XMOMD = sinTheta*std::cos(phi), YMOMD = sinTheta*std::sin(phi), ZMOMD = cosTheta;
    */

    /*                          //// uniforme

        XMOMD =  -1 + (G4double)G4UniformRand()*2; YMOMD = -1 + (G4double)G4UniformRand()*2 ; ZMOMD = -1 + (G4double)G4UniformRand()*2 ;
        //return G4ThreeVector(-1 + (G4double)G4UniformRand()*2, -1 + (G4double)G4UniformRand()*2, -1 + (G4double)G4UniformRand()*2);
    */

    /*                          //// with theta and phy limits

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
        XMOMD = std::sin(theta)*std::cos(phi);
        YMOMD = std::sin(theta)*std::sin(phi);
        ZMOMD = std::cos(theta);
    */

    /*                          //// directed according to ThetaPhi

        XMOMD = std::sin(NewRankSourceMomDirsDirectedThetaValues[DataID])*std::cos(NewRankSourceMomDirsDirectedPhiValues[DataID]);
        YMOMD = std::sin(NewRankSourceMomDirsDirectedThetaValues[DataID])*std::sin(NewRankSourceMomDirsDirectedPhiValues[DataID]);
        ZMOMD = std::cos(NewRankSourceMomDirsDirectedThetaValues[DataID]);

    */

    /*                          //// directed to point

        G4ThreeVector momentumDirection = (G4ThreeVector(DirectedToX, DirectedToX, DirectedToX) - G4ThreeVector(X, Y, Z)).unit();
        XMOMD = momentumDirection.getX();
        YMOMD = momentumDirection.getY();
        ZMOMD = momentumDirection.getZ();
    */

    /*                          //// directed parallel to Z, X, or Y plane, and with cordinate A and B for other axis

        G4ThreeVector momentumDirection;
        if(DirectedParallelAxis == "YZ"){
            momentumDirection = (G4ThreeVector(X, DirectedToX, DirectedToY) - G4ThreeVector(X, Y, Z)).unit();
        }else if(DirectedParallelAxis == "ZX"){
            momentumDirection = (G4ThreeVector(DirectedToX, Y, DirectedToY) - G4ThreeVector(X, Y, Z)).unit();
        }else{
            momentumDirection = (G4ThreeVector(DirectedToX, DirectedToY, Z) - G4ThreeVector(X, Y, Z)).unit();
        }
        XMOMD = momentumDirection.getX();
        YMOMD = momentumDirection.getY();
        ZMOMD = momentumDirection.getZ();

    */

    /*                          //// directed to a cubique colume isotropically

        G4ThreeVector momentumDirection = (G4ThreeVector(
                                               (ToVolumeX+(G4UniformRand()*2.*DirectedToX) - DirectedToX),
                                               (ToVolumeY+(G4UniformRand()*2.*DirectedToY) - DirectedToY),
                                               (ToVolumeZ+(G4UniformRand()*2.*DirectedToZ) - DirectedToZ)) -
                                           G4ThreeVector(X, Y, Z)).unit();
        XMOMD = momentumDirection.getX();
        YMOMD = momentumDirection.getY();
        ZMOMD = momentumDirection.getZ();
    */

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    particleGun->SetParticlePosition(G4ThreeVector(X, Y, Z));
    particleGun->SetParticleEnergy(ENERGY);
    particleGun->SetParticleMomentumDirection(G4ParticleMomentum(XMOMD, YMOMD, ZMOMD));

    /*                          //// Save Data To file After generating ENERGY X Y Z XMOMD YMOMD ZMOMD

        SaveGeneratedDataToFiles();
    */

    //std::cout << EvInc << " ENERGY=" << ENERGY << " X=" << X << " Y=" << Y << " Z=" << Z << " " << NewRankSourceRegionsBoxDimValues[DataID] << " " << particleDefinition->GetParticleName() << std::endl;

    particleGun->GeneratePrimaryVertex(anEvent);

}


