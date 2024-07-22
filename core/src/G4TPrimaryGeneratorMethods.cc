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

#include "G4TPrimaryGeneratorMethods.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4TVolumeConstruction.hh"

#include "G4RunManager.hh"
//#include "G4ios.hh"
#include <iostream>
#include "G4PhysicalConstants.hh"
#include <dirent.h>

extern  G4String* CopyNumberRegionNameMap;
extern  G4float* CopyNumberXPos;
extern  G4float* CopyNumberYPos;
extern  G4float* CopyNumberZPos;
extern double* CumulativeActivities;
// Source Data

extern G4int EventsNumPerThreadRank;
extern unsigned int* ParNameList ;
extern double* EnergyList        ;
extern double* MomDirXList       ;
extern double* MomDirYList       ;
extern double* MomDirZList       ;
extern double* PosXList          ;
extern double* PosYList          ;
extern double* PosZList          ;

extern G4bool UseDicomCumAct;

extern G4String particleName;
extern G4String SourceType;
extern G4String SourceRegionName;
extern G4ThreeVector boxCenterPos;
extern G4int NumberOfGenPointsToSave;
extern G4ThreeVector BoxDimGene;
extern G4ThreeVector newOrgPos3Vector;
extern G4ThreeVector SourcePosition;
extern G4ThreeVector SourceRotation;
extern G4String SourceSolid;
extern G4String SourceSurface;
extern G4String SourcePlane;
extern G4String SourceAxis;
extern G4double Radius;
extern G4double RadiusIn;
extern G4double BeamSDev;
extern G4double HalfX, HalfY, HalfZ;
extern G4double ThetaMin;
extern G4double ThetaMax;
extern G4double PhiMin;
extern G4double PhiMax;
extern G4String DirectedParallelAxis;
extern G4double ToVolumeX;
extern G4double ToVolumeY;
extern G4double ToVolumeZ;
extern G4double DirectedToX;
extern G4double DirectedToY;
extern G4double DirectedToZ;
extern G4String MomDirDirectedHow; // ThetaPhi, ToVolume, ToPoint, ParallelTo

extern G4double RotTheta;
extern G4double RotPhi;
extern G4String RotPosAxis;

extern G4ThreeVector SourceRotVector1;
extern G4ThreeVector SourceRotVector2;

extern G4String EnergyDistribution;
extern G4double GaussSDev;
extern G4double GaussMean;
extern G4double UniformEmin;
extern G4double UniformEmax;
extern G4double RayleighEmax;
extern G4double MonoEnergy;
extern std::map<G4double, G4double> EnergyValueProbability;

extern std::vector<double> RadioNuclideProbVec;
extern std::vector<unsigned int> RadioNuclidePartNameVec;
extern std::vector<double*> RadioNuclideEneVec;
extern std::vector<unsigned int> RadioNuclideSpectrumOrDiscreteVec;


extern G4String MomDirDistribution;
extern G4double Theta;
extern G4double Phi;

extern G4String DataFilesExtension ;

extern G4String GeneratePosFlag ;
extern G4String GenerateEneFlag ;
extern G4String GenerateMomDirFlag ;
extern G4String UseGeneratedData ;
extern G4String ShowBox ;
extern G4String TestPointsPositions ;

extern G4String PositionDataFile, EnergyDataFile, MomDirDataFile;

extern G4String ScriptsDirectoryPath, DataDirectoryPath;

// Run And Score
extern G4String MPISimulationNum;
extern G4int BatchsNumber;
extern G4int ThreadsNumber;
extern G4String AccuracyCalculationLevel;
extern G4String organs_to_score ;
extern G4String variable_To_Score ;

// Voxelized Geometry
extern G4double VoxXHalfSize;
extern G4double VoxYHalfSize;
extern G4double VoxZHalfSize;
extern G4int VoxXNumber;
extern G4int VoxYNumber;
extern G4int VoxZNumber;

extern std::vector<G4String> NewRankSourceParticlesNamesValues;
extern std::vector<G4String> NewRankSourceRegionsNamesValues;
extern std::vector<G4ThreeVector> NewRankSourceRegionsBoxDimValues;
extern std::vector<G4ThreeVector> NewRankSourceRegionsPosValues;
extern std::vector<std::vector<unsigned int>> NewRankVoxelsIDsOfSourceRegion;
extern std::vector<std::vector<G4Tet*>> NewRankTETOfSourceRegion;
extern std::vector<G4ThreeVector> NewRankTETBoxMinOfSourceRegion;
extern std::vector<G4ThreeVector> NewRankTETBoxDimOfSourceRegion;

extern std::vector<G4double> NewRankSourceEnergiesValues;
extern std::vector<G4String> NewRankSourceMomDirsValues;
extern std::vector<G4double> NewRankSourceMomDirsDirectedThetaValues;
extern std::vector<G4double> NewRankSourceMomDirsDirectedPhiValues;
extern std::vector<G4String> SourceRegionsNamesToBeIgnoredValues;
extern std::vector<G4Navigator*> NavigatorForVolumesInitialPosition;
//extern std::vector<G4ParticleDefinition*> ParticleDefinitions;

extern std::vector<G4String> NewRankSourceEnergyDataFiles;
extern std::vector<G4String> NewRankSourcePositionDataFiles;
extern std::vector<G4String> NewRankSourceMomDirDataFiles;

extern std::map<G4int,std::vector<G4double>> FissionCapturesOfThreadsRanks;


extern std::string getFileNameFromPath(std::string const & path, std::string const & delims = "/\\");

#ifdef G4MULTITHREADED
G4ThreadLocal G4ParticleDefinition* G4TPrimaryGeneratorMethods::particleDefinition;
G4ThreadLocal G4Navigator* G4TPrimaryGeneratorMethods::aNavigator;
G4ThreadLocal G4VPhysicalVolume* G4TPrimaryGeneratorMethods::WorldPhysicalVolume;
G4ThreadLocal unsigned int G4TPrimaryGeneratorMethods::EvInc;
G4ThreadLocal unsigned int G4TPrimaryGeneratorMethods::ParticleID;
G4ThreadLocal std::map<unsigned int,G4ParticleDefinition* > G4TPrimaryGeneratorMethods::particleDefinitionList;
G4ThreadLocal std::map<G4String,std::map<G4String,std::map<int,G4String>>> G4TPrimaryGeneratorMethods::PosDataFileNames;
G4ThreadLocal std::map<G4String,std::map<double,std::map<int,G4String>>> G4TPrimaryGeneratorMethods::EneDataFileNames;
G4ThreadLocal std::map<G4String,std::map<int,G4String>> G4TPrimaryGeneratorMethods::MomDirDataFileNames;
G4ThreadLocal std::map<G4String,std::map<int,G4String>> G4TPrimaryGeneratorMethods::CriticalityDataFileNames;

G4ThreadLocal unsigned int G4TPrimaryGeneratorMethods::EnergyListInc;
//G4ThreadLocal G4double* G4TPrimaryGeneratorMethods::EnergyList;
//G4ThreadLocal double* G4TPrimaryGeneratorMethods::PosXList;
//G4ThreadLocal double* G4TPrimaryGeneratorMethods::PosYList;
//G4ThreadLocal double* G4TPrimaryGeneratorMethods::PosZList;
//G4ThreadLocal double* G4TPrimaryGeneratorMethods::MomDirXList;
//G4ThreadLocal double* G4TPrimaryGeneratorMethods::MomDirYList;
//G4ThreadLocal double* G4TPrimaryGeneratorMethods::MomDirZList;
//G4ThreadLocal unsigned int* G4TPrimaryGeneratorMethods::ParNameList;

G4ThreadLocal int G4TPrimaryGeneratorMethods::DataID;
G4ThreadLocal int G4TPrimaryGeneratorMethods::WriteSourceDataToFiles;
G4ThreadLocal G4double G4TPrimaryGeneratorMethods::TotalEmittedEnergy;
G4ThreadLocal G4double G4TPrimaryGeneratorMethods::ENERGY;
G4ThreadLocal G4double* G4TPrimaryGeneratorMethods::ReadedEnergyList;
G4ThreadLocal G4ThreeVector* G4TPrimaryGeneratorMethods::ReadedPositionsList;
G4ThreadLocal G4ThreeVector* G4TPrimaryGeneratorMethods::ReadedMomDirecsList;

G4ThreadLocal G4double G4TPrimaryGeneratorMethods::X;
G4ThreadLocal G4double G4TPrimaryGeneratorMethods::Y;
G4ThreadLocal G4double G4TPrimaryGeneratorMethods::Z;
G4ThreadLocal G4double G4TPrimaryGeneratorMethods::XMOMD;
G4ThreadLocal G4double G4TPrimaryGeneratorMethods::YMOMD;
G4ThreadLocal G4double G4TPrimaryGeneratorMethods::ZMOMD;
G4ThreadLocal G4int G4TPrimaryGeneratorMethods::VoxelsInc;
G4ThreadLocal G4int G4TPrimaryGeneratorMethods::CummNumbInVoxelsInc;
G4ThreadLocal int G4TPrimaryGeneratorMethods::PositionTypeNum;
G4ThreadLocal int G4TPrimaryGeneratorMethods::EnergyTypeNum;
G4ThreadLocal int G4TPrimaryGeneratorMethods::MomDirTypeNum;
G4ThreadLocal std::ofstream G4TPrimaryGeneratorMethods::PositionFileStream;
G4ThreadLocal std::ofstream G4TPrimaryGeneratorMethods::MomDirFileStream;
G4ThreadLocal std::ofstream G4TPrimaryGeneratorMethods::EneFileStream;
#endif

namespace
{
//G4Mutex	mutex = G4MUTEX_INITIALIZER;
}

G4TPrimaryGeneratorMethods::G4TPrimaryGeneratorMethods(){
    
}

G4TPrimaryGeneratorMethods::~G4TPrimaryGeneratorMethods(){}

void G4TPrimaryGeneratorMethods::GenerateEventsParticle(){
    
    if(EnergyTypeNum == 5){
        particleDefinitionList[ParNameList[EnergyListInc]];
    }
}

void G4TPrimaryGeneratorMethods::GenerateEventsPosition(){
    
    //std::cout << "DataID "<< DataID << " SourceType=" << SourceType << " NewRankSourceRegionsNamesValues[DataID]=" << NewRankSourceRegionsNamesValues[DataID] << " " << G4ThreeVector(X,Y,Z) << " is in " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() << std::endl;

    if(PositionTypeNum == 1){
        
        //std::cout << " NewRankSourceRegionsNamesValues[DataID] " << std::endl;

        //// stylized volumes

        X = (NewRankSourceRegionsPosValues[DataID].getX() - NewRankSourceRegionsBoxDimValues[DataID].getX()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getX()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
        Y = (NewRankSourceRegionsPosValues[DataID].getY() - NewRankSourceRegionsBoxDimValues[DataID].getY()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getY()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
        Z = (NewRankSourceRegionsPosValues[DataID].getZ() - NewRankSourceRegionsBoxDimValues[DataID].getZ()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getZ()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
        
        //std::cout << "DataID "<< DataID << " SourceType=" << SourceType << " NewRankSourceRegionsNamesValues[DataID]=" << NewRankSourceRegionsNamesValues[DataID] << " " << G4ThreeVector(X,Y,Z) << " is in " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() << std::endl;
        
        if(NewRankSourceRegionsNamesValues[DataID] == "allregions"){
            
            bool Save = true;
            
            for (int gg = 0 ; gg < SourceRegionsNamesToBeIgnoredValues.size() ; gg++) {
                //if(SourceRegionsNamesToBeIgnoredValues[gg] == aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName()){
                if(SourceRegionsNamesToBeIgnoredValues[gg] == NavigatorForVolumesInitialPosition[DataID]->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName()){
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
                    //if(SourceRegionsNamesToBeIgnoredValues[gg] == NavigatorForVolumesInitialPosition[DataID]->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName()){
                    if(SourceRegionsNamesToBeIgnoredValues[gg] == NavigatorForVolumesInitialPosition[DataID]->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName()){
                        Save=false;
                    }
                }

                //std::cout << "DataID "<< DataID << " SourceType=" << SourceType << " NewRankSourceRegionsNamesValues[DataID]=" << NewRankSourceRegionsNamesValues[DataID] << " " << G4ThreeVector(X,Y,Z) << " is in " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() << std::endl;
                
            }
            
        }else{
            //while (NavigatorForVolumesInitialPosition[DataID]->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() != NewRankSourceRegionsNamesValues[DataID]){
            while (NavigatorForVolumesInitialPosition[DataID]->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() != NewRankSourceRegionsNamesValues[DataID]){

                //std::cout << "DataID "<< DataID << " SourceType=" << SourceType << " NewRankSourceRegionsBoxDimValues[DataID]=" << NewRankSourceRegionsBoxDimValues[DataID] << " NewRankSourceRegionsNamesValues[DataID]=" << NewRankSourceRegionsNamesValues[DataID] << " NewRankSourceRegionsPosValues[DataID]=" << NewRankSourceRegionsPosValues[DataID] << " Generated Pos: " << G4ThreeVector(X,Y,Z) << " is in " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName() << std::endl;
                
                X = (NewRankSourceRegionsPosValues[DataID].getX() - NewRankSourceRegionsBoxDimValues[DataID].getX()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getX()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
                Y = (NewRankSourceRegionsPosValues[DataID].getY() - NewRankSourceRegionsBoxDimValues[DataID].getY()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getY()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
                Z = (NewRankSourceRegionsPosValues[DataID].getZ() - NewRankSourceRegionsBoxDimValues[DataID].getZ()) + (G4double)G4UniformRand() * (2*NewRankSourceRegionsBoxDimValues[DataID].getZ()); //generateRandom(hXmin,NewRankSourceRegionsBoxDimValues[DataID].getX());
            }
        }
    }
    else if(PositionTypeNum == 2){
        
        //// position in solid
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
        
        //// position in solid surface

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
        
        //// position in beam

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
        //return G4ThreeVector(X, Y, Z);
        
    }
    else if(PositionTypeNum == 5){
        
        //// position in plane

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

    }
    else if(PositionTypeNum == 6){

        //// position in point

        X = SourcePosition.getX(), Y = SourcePosition.getY(), Z = SourcePosition.getZ();
        //return G4ThreeVector(SourcePosition.getX(), SourcePosition.getY(), SourcePosition.getZ());
    }
    else if(PositionTypeNum == 7){
        
        //// position in rotated form

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
        
        //// position in Tetrahedral geometry

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
    else if(PositionTypeNum == 0){
        if(VoxelsInc == NewRankVoxelsIDsOfSourceRegion[DataID].size()){ VoxelsInc=0; }

        //std::cout  << " VoxelsInc/TotalInSource " << VoxelsInc << "/" << NewRankVoxelsIDsOfSourceRegion[DataID].size() << " ID=" << NewRankVoxelsIDsOfSourceRegion[DataID][VoxelsInc]<< " XP=" << CopyNumberXPos[NewRankVoxelsIDsOfSourceRegion[DataID][VoxelsInc]] << " " << CopyNumberRegionNameMap[NewRankVoxelsIDsOfSourceRegion[DataID][VoxelsInc]] << " == "<< NewRankSourceRegionsNamesValues[DataID] << std::endl;
        X = (CopyNumberXPos[NewRankVoxelsIDsOfSourceRegion[DataID][VoxelsInc]] - VoxXHalfSize) + (G4double)G4UniformRand() * (2*VoxXHalfSize); //generateRandom(hXmin,hXmax);
        Y = (CopyNumberYPos[NewRankVoxelsIDsOfSourceRegion[DataID][VoxelsInc]] - VoxYHalfSize) + (G4double)G4UniformRand() * (2*VoxYHalfSize); //generateRandom(hXmin,hXmax);
        Z = (CopyNumberZPos[NewRankVoxelsIDsOfSourceRegion[DataID][VoxelsInc]] - VoxZHalfSize) + (G4double)G4UniformRand() * (2*VoxZHalfSize); //generateRandom(hXmin,hXmax);
        //std::cout  << " ------------------ " << std::endl;

        VoxelsInc++;
    }
    if(!RotPosAxis.empty() && RotTheta != 0.){
        RotatePosition();
    }

    //return G4ThreeVector(X, Y, Z);
}
void G4TPrimaryGeneratorMethods::RotatePosition(){

    G4RotationMatrix rotationMatrix;

    if(RotPosAxis == "Z"){

        // // Apply rotation around Z-axis
        //X = cosAngle * X - sinAngle * Y;
        //Y = sinAngle * X + cosAngle * Y;

        //// Rotation angle around Z-axis in degrees
        //double r = sqrt(X*X + Y*Y + Z*Z);
        //double theta1 = acos(Z / r); // inclination angle
        //double phi1 = atan2(Y, X); // azimuthal angle

        //X = r * sin(theta1) * cos(phi1+RotTheta*degree);
        //Y = r * sin(theta1) * sin(phi1+RotTheta*degree);
        // //Z = r * cos(theta1);
        rotationMatrix.rotateZ(RotTheta*degree);
        G4ThreeVector rotatedPoint = rotationMatrix * G4ThreeVector(X,Y,Z);

        //std::cout << "RotTheta "<< RotTheta << std::endl;

        X = rotatedPoint.getX();
        Y = rotatedPoint.getY();
        Z = rotatedPoint.getZ();
    }else if(RotPosAxis == "X"){

        rotationMatrix.rotateX(RotTheta*degree);
        G4ThreeVector rotatedPoint = rotationMatrix * G4ThreeVector(X,Y,Z);

        //std::cout << "RotTheta "<< RotTheta << std::endl;

        X = rotatedPoint.getX();
        Y = rotatedPoint.getY();
        Z = rotatedPoint.getZ();

    }else if(RotPosAxis == "Y"){

        rotationMatrix.rotateY(RotTheta*degree);
        G4ThreeVector rotatedPoint = rotationMatrix * G4ThreeVector(X,Y,Z);

        //std::cout << "RotTheta "<< RotTheta << std::endl;

        X = rotatedPoint.getX();
        Y = rotatedPoint.getY();
        Z = rotatedPoint.getZ();
    }

}

void G4TPrimaryGeneratorMethods::GenerateEventsEnergy(){
    
    //std::cout << __FUNCTION__ << std::endl;

    //G4cout << "************** Generate " << NumberOfGenPointsToSave << " Event Energy(MeV) and Saving the data to " << EnergyDataFile<< G4endl;
    
    //G4double pi = 4*std::arctan(1);  pi/4 = 4*arctan(1/5) - arctan(1/239) (Machin's formula) , or using towpi from G4PhysicalConstants.hh
    if(EnergyTypeNum == 0){

        //// mono energy distribution

        ENERGY = NewRankSourceEnergiesValues[DataID] ;
        //return MonoEnergy;
        //std::ostringstream text; text <<"\r"<< " " << numbSucces <<"/"<<NumberOfGenPointsToSave << " "; printf(text.str().c_str());
    }
    else if(EnergyTypeNum == 1){

        //// uniform energy distribution

        ENERGY = UniformEmin + (NewRankSourceEnergiesValues[DataID] - UniformEmin) * (G4double)G4UniformRand() ;
        //return UniformEmin + (UniformEmax - UniformEmin) * (G4double)G4UniformRand();
    }
    else if(EnergyTypeNum == 2){

        //// rayleigh energy distribution

        G4double mean = NewRankSourceEnergiesValues[DataID]/3. ;
        G4double sigma = mean * std::sqrt(1/twopi) ;
        ENERGY = sigma*std::sqrt(-2*std::log(1-(G4double)G4UniformRand()));
        //return RayleighEmax/3.* std::sqrt(1/twopi)*std::sqrt(-2*std::log(1-(G4double)G4UniformRand()));
    }
    else if(EnergyTypeNum == 3){

        //// gauss distribution

        ENERGY = CLHEP::RandGauss::shoot(NewRankSourceEnergiesValues[DataID], GaussSDev) ;
        //return CLHEP::RandGauss::shoot(GaussMean, GaussSDev);
    }
    else if(EnergyTypeNum == 4){

        //// Spectrum distribution

        if(EnergyListInc == EventsNumPerThreadRank ){
            //std::cout << "EnergyListInc="<< EnergyListInc << " EventsNumPerThreadRank=" << EventsNumPerThreadRank << G4endl;
            EnergyListInc= 0;
        }

        ENERGY = EnergyList[EnergyListInc] ;
        EnergyListInc++;
    }
    else if(EnergyTypeNum == 5){

        //// distribution data from file or radionuclide

        if(EnergyListInc == EventsNumPerThreadRank ){
            //std::cout << "EnergyListInc="<< EnergyListInc << " EventsNumPerThreadRank=" << EventsNumPerThreadRank << G4endl;
            EnergyListInc= 0;
        }

        ENERGY = EnergyList[EnergyListInc] ;
        EnergyListInc++;
    }
    //G4cout << UniformEmin << " ENERGY=" << ENERGY << " " << NewRankSourceEnergiesValues[DataID] << G4endl;
    
    TotalEmittedEnergy += ENERGY ;
    
    //G4cout << EnergyListInc << " ENERGY=" << ENERGY << " " << TotalEmittedEnergy << G4endl; EnergyListInc++;
    
    //return ENERGY;    
}

void G4TPrimaryGeneratorMethods::GenerateEventsMomentumDirection(){
    
    //std::cout << __FUNCTION__ << std::endl;

    if(MomDirTypeNum == 0){

        //// isotropic

        G4double cosTheta = 2*G4UniformRand() - 1., phi = twopi*G4UniformRand();
        G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
        XMOMD = sinTheta*std::cos(phi), YMOMD = sinTheta*std::sin(phi), ZMOMD = cosTheta;
        //return G4ThreeVector(sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta);
    }
    else if(MomDirTypeNum == 1){

        //// uniforme

        XMOMD =  -1 + (G4double)G4UniformRand()*2; YMOMD = -1 + (G4double)G4UniformRand()*2 ; ZMOMD = -1 + (G4double)G4UniformRand()*2 ;
        //return G4ThreeVector(-1 + (G4double)G4UniformRand()*2, -1 + (G4double)G4UniformRand()*2, -1 + (G4double)G4UniformRand()*2);
    }
    else if(MomDirTypeNum == 2){

        //// with theta and phy limits

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

        // G4cout << " XMOMD=" << XMOMD << " YMOMD=" << YMOMD << " ZMOMD=" << ZMOMD << G4endl;

        //return G4ThreeVector(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));

    }
    else if(MomDirTypeNum == 3){

        //// directed according to ThetaPhi

        //G4cout << " NewRankSourceMomDirsDirectedThetaValues[DataID]=" << NewRankSourceMomDirsDirectedThetaValues[DataID]
        //           << " NewRankSourceMomDirsDirectedPhiValues[DataID]=" << NewRankSourceMomDirsDirectedPhiValues[DataID] << G4endl;

        XMOMD = std::sin(NewRankSourceMomDirsDirectedThetaValues[DataID])*std::cos(NewRankSourceMomDirsDirectedPhiValues[DataID]);
        YMOMD = std::sin(NewRankSourceMomDirsDirectedThetaValues[DataID])*std::sin(NewRankSourceMomDirsDirectedPhiValues[DataID]);
        ZMOMD = std::cos(NewRankSourceMomDirsDirectedThetaValues[DataID]);

    }
    else if(MomDirTypeNum == 4){

        //// directed to point

        G4ThreeVector momentumDirection = (G4ThreeVector(DirectedToX, DirectedToX, DirectedToX) - G4ThreeVector(X, Y, Z)).unit();

        XMOMD = momentumDirection.getX();
        YMOMD = momentumDirection.getY();
        ZMOMD = momentumDirection.getZ();
    }
    else if(MomDirTypeNum == 5){

        //// directed parallel to Z, X, or Y plane, and with cordinate A and B for other axis

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

    }
    else if(MomDirTypeNum == 6){

        //// directed to a cubique colume isotropically

        G4ThreeVector momentumDirection = (G4ThreeVector(
                                               (ToVolumeX+(G4UniformRand()*2.*DirectedToX) - DirectedToX),
                                               (ToVolumeY+(G4UniformRand()*2.*DirectedToY) - DirectedToY),
                                               (ToVolumeZ+(G4UniformRand()*2.*DirectedToZ) - DirectedToZ)) -
                                           G4ThreeVector(X, Y, Z)).unit();

        XMOMD = momentumDirection.getX();
        YMOMD = momentumDirection.getY();
        ZMOMD = momentumDirection.getZ();

        //G4cout << " XMOMD=" << XMOMD << " YMOMD=" << YMOMD << " ZMOMD=" << ZMOMD << G4endl;

    }


    //return G4ParticleMomentum(XMOMD, YMOMD, ZMOMD);
}

void G4TPrimaryGeneratorMethods::OpenFilesToSaveGeneratedData(){
    
    std::ostringstream a ;
    if(SourceType == "Solid"){
        a << DataDirectoryPath  << "/Pos_" << NewRankSourceRegionsNamesValues[DataID] << "_" << SourceSolid << "_" << SourceType;
    }
    else if(SourceType == "Surface"){
        a << DataDirectoryPath  << "/Pos_" << NewRankSourceRegionsNamesValues[DataID] << "_" << SourceSurface << "_" << SourceType;
    }
    else if(SourceType == "Plane"){
        a << DataDirectoryPath  << "/Pos_" << NewRankSourceRegionsNamesValues[DataID] << "_" << SourcePlane << "_" << SourceType;
    }
    else {
        a << DataDirectoryPath  << "/Pos_" << NewRankSourceRegionsNamesValues[DataID] << "_" << SourceType;
    }
    a << "_" << EventsNumPerThreadRank << "_" <<DataID << DataFilesExtension;
    PositionFileStream.open(a.str().c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
    
    std::ostringstream b ;
    b << DataDirectoryPath  << "/Ene_" << EnergyDistribution <<"_";
    if(EnergyDistribution == "Spectrum"){
        b << SpectrumMaxEnergy;
    }
    else if(EnergyDistribution == "File"){
        b << FileEnergyCharacterizer;
    }
    else if(EnergyDistribution == "RadioNuclide"){
        b << RadioNuclideMaxEnergy;
    }
    else{
        b << NewRankSourceEnergiesValues[DataID];
    }
    b << "_" << EventsNumPerThreadRank << "_" << DataID << DataFilesExtension;
    
    EneFileStream.open(b.str().c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
    
    std::ostringstream c ;
    c << DataDirectoryPath  << "/MomDir_" << MomDirDistribution ; // or NewRankSourceMomDirsValues[DataID]
    if(MomDirDistribution == "Isotropic"){
        //a << "";
    }
    else if(MomDirDistribution == "Uniform"){
        //a << "";
    }
    else if(MomDirDistribution == "Directed"){
        c << Theta << "_" << Phi ;
    }
    c <<  "_" << EventsNumPerThreadRank << "_" <<DataID<< DataFilesExtension;
    MomDirFileStream.open(c.str().c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
    
}

void G4TPrimaryGeneratorMethods::CloseFilesToSaveGeneratedData(){
    PositionFileStream.close();
    MomDirFileStream.close();
    EneFileStream.close();
}

void G4TPrimaryGeneratorMethods::SaveGeneratedDataToFiles(){
    EneFileStream << ENERGY <<" " ;
    PositionFileStream << X << " " << Y << " " << Z <<"\n" ;
    MomDirFileStream << XMOMD << " " << YMOMD << " " << ZMOMD <<"\n" ;
}

void G4TPrimaryGeneratorMethods::SourceInitialization(){

    WriteSourceDataToFiles = 0;

    if(WriteSourceDataToFiles == 1){
        OpenFilesToSaveGeneratedData();
    }

    // for initial positions generation
    
    if(SourceType == "Voxels"){PositionTypeNum = 0;
        
        if(NewRankVoxelsIDsOfSourceRegion.size() == 0 || NewRankVoxelsIDsOfSourceRegion[DataID].size() == 0 ){
            G4String msg = "Source region (" + NewRankSourceRegionsNamesValues[DataID] + ") in rank/thread (" + std::to_string(DataID) + "), number voxels of voxels equal to 0 "; G4Exception("Source Data", "1", FatalErrorInArgument, msg.c_str());
        }
        
        //const G4TVolumeConstruction* TConstruction = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        //Voxel0PosX = -((VoxXNumber*VoxXHalfSize) - VoxXHalfSize) + TConstruction->getVoxContainerPos().getX();
        //Voxel0PosY = -((VoxYNumber*VoxYHalfSize) - VoxYHalfSize) + TConstruction->getVoxContainerPos().getY();
        //Voxel0PosZ = -((VoxZNumber*VoxZHalfSize) - VoxZHalfSize) + TConstruction->getVoxContainerPos().getZ();
    }
    else if(SourceType == "Volume"){PositionTypeNum = 1;
        
        //const G4TVolumeConstruction* TConstruction = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        //WorldPhysicalVolume = TConstruction->getWorldPhyVolume();
        //
        //aNavigator = new G4Navigator();
        //aNavigator->SetWorldVolume(WorldPhysicalVolume);




        /*
                std::cout << "DataID "<< DataID << " SourceType=" << SourceType
                          << " WorldLogName=" << WorldPhysicalVolume->GetLogicalVolume()->GetName()
                          << " WorldPyhName=" << WorldPhysicalVolume->GetName()
                          << " EvInc=" << EvInc
                          << " RegNmSize " << NewRankSourceRegionsPosValues .size()
                          << " BoxDimSize " << NewRankSourceRegionsBoxDimValues.size()
                          << " Navigator test " << aNavigator->LocateGlobalPointAndSetup(G4ThreeVector(X,Y,Z))->GetLogicalVolume()->GetName()
                          << std::endl;
                          */
    }
    else if(SourceType == "Solid"){PositionTypeNum = 2;}
    else if(SourceType == "Surface"){PositionTypeNum = 3;}
    else if(SourceType == "Beam"){PositionTypeNum = 4;}
    else if(SourceType == "Plane"){PositionTypeNum = 5;}
    else if(SourceType == "Point"){PositionTypeNum = 6;}
    else if(SourceType == "Rotated"){PositionTypeNum = 7;}
    else if(SourceType == "TET"){PositionTypeNum = 8;}
    //else if(SourceType == "FissionSites"){PositionTypeNum = 9;}

    // for initial MomDires generation
    
    if(NewRankSourceMomDirsValues[DataID] == "Isotropic"){MomDirTypeNum = 0;}
    else if(NewRankSourceMomDirsValues[DataID] == "Uniform"){MomDirTypeNum = 1;}
    else if(NewRankSourceMomDirsValues[DataID] == "Delimited"){MomDirTypeNum = 2;}
    else if(NewRankSourceMomDirsValues[DataID] == "Directed"){
        if(MomDirDirectedHow == "ThetaPhi"){
            MomDirTypeNum = 3;
        }
        else if(MomDirDirectedHow == "ToPoint"){
            MomDirTypeNum = 4;
        }
        else if(MomDirDirectedHow == "ParallelTo"){
            MomDirTypeNum = 5;
        }
        else if(MomDirDirectedHow == "ToVolume"){
            MomDirTypeNum = 6;
        }
    }
    //else if(NewRankSourceMomDirsValues[DataID] == "Delimited"){MomDirTypeNum = 3;}

    // for initial energies generation

    if(EnergyDistribution == "Mono"){EnergyTypeNum = 0;}
    else if(EnergyDistribution == "Uniform"){EnergyTypeNum = 1;}
    else if(EnergyDistribution == "Rayleigh"){EnergyTypeNum = 2;}
    else if(EnergyDistribution == "Gauss"){EnergyTypeNum = 3;}

    return;//else
        if(EnergyDistribution == "Spectrum"){EnergyTypeNum = 4;

        G4double TotalProbability = 0;
        
        /*
        G4double TotalYield = 0;
        G4double EnergyProbability;
        G4double EnergyYield;
        
        G4double TotalNumberOfRegisteredEventsByProbability = 0;
        G4double TotalNumberOfRegisteredEventsByYield = 0;
        EnergyValues = new G4double[EventsNumPerThreadRank];
        
        for ( auto Abeg = EnergyValueProbability.begin(); Abeg != EnergyValueProbability.end(); ++Abeg  ){
            EnergyProbability = Abeg->second;
            TotalProbability += EnergyProbability;
            EnergyYield = Abeg->second;;
            TotalYield += EnergyYield;
            //TotalProbability += Abeg->second;
        }
        for ( auto Abeg = EnergyValueProbability.begin(); Abeg != EnergyValueProbability.end(); ++Abeg  ){
        
            EnergyProbability = Abeg->second;
            int enevalNum = EventsNumPerThreadRank * (EnergyProbability/TotalProbability);
            int enevalNumByYield = EventsNumPerThreadRank * Abeg->second;
            
            //std::cout << " EnergyProbability " << EnergyProbability  << " enevalNum " << enevalNum << std::endl;
            
            for(int f = TotalNumberOfRegisteredEventsByProbability; f < TotalNumberOfRegisteredEventsByProbability+enevalNum ;f++ ){
                EnergyValues[f] = Abeg->first;
                //std::cout << f << "  " << EnergyValues[f] << " enevalNum " << enevalNum << " enevalNumByYield " << enevalNumByYield << std::endl;
            }
            
            std::cout << " Energy "<< Abeg->first << " EnergyYield " << Abeg->second << " enevalNumByYield " << enevalNumByYield << " EnergyProbability " << EnergyProbability << " enevalNum " << enevalNum << std::endl;
            TotalNumberOfRegisteredEventsByProbability += enevalNum;
            TotalNumberOfRegisteredEventsByYield += enevalNumByYield;
            
        }
        
        std::cout  << " TotalProbability " << TotalProbability
                   << " TotalNumberOfRegisteredEventsByProbability " << TotalNumberOfRegisteredEventsByProbability
                   << " TotalYield " << TotalYield
                      //<< " TotalNumberOfRegisteredEventsByYield " << TotalNumberOfRegisteredEventsByYield
                   << " EventsNumPerThreadRank " << EventsNumPerThreadRank << std::endl;
                   
*/
        
        std::map<G4double, G4double> EnergyProbabilityInterval;
        std::map<G4int, G4double> CumulatedProbabilityIndexEnergy;
        
        G4int index=0;
        
        for ( auto Abeg = EnergyValueProbability.begin(); Abeg != EnergyValueProbability.end(); ++Abeg  ){
            
            EnergyProbabilityInterval[Abeg->first] = Abeg->second;
            TotalProbability += EnergyProbabilityInterval[Abeg->first];
            
#if VERBOSE_USE
            //std::cout  << index << " TotalProbability "<< TotalProbability << " - this EnergyProbabilityInterval " << EnergyProbabilityInterval[Abeg->first] << " Energy " << std::endl;
#endif
            CumulatedProbabilityIndexEnergy[index] = Abeg->first;
            index++;
        }
        
        EnergyList = new G4double[EventsNumPerThreadRank];
        
        for(int f = 0; f < EventsNumPerThreadRank ;f++ ){
            
            G4double RV = TotalProbability*(G4double)G4UniformRand() ;
            G4double EnergyCumulatedProbability = 0;
            
            index=0;
            for ( auto Abeg = EnergyValueProbability.begin(); Abeg != EnergyValueProbability.end(); ++Abeg  ){
                
                EnergyCumulatedProbability += Abeg->second;
                
                if(RV < EnergyCumulatedProbability){
                    EnergyList[f] = CumulatedProbabilityIndexEnergy[index];
                    break;
                }
                index++;
            }
            
            //std::cout << f << " Energy " << EnergyValues[f] << " EnergyCumulatedProbability " << EnergyCumulatedProbability << " - this EnergyProbabilityInterval " << EnergyProbabilityInterval[EnergyValues[f]] << std::endl;
        }
    }
    else if(EnergyDistribution == "RadioNuclide" || EnergyDistribution == "File"){EnergyTypeNum = 5;
        
        ParNameList = new unsigned int[EventsNumPerThreadRank];
        EnergyList = new double[EventsNumPerThreadRank];
        //MomDirXList = new double[EventsNumPerThreadRank];
        //MomDirYList = new double[EventsNumPerThreadRank];
        //MomDirZList = new double[EventsNumPerThreadRank];
        //PosXList = new double[EventsNumPerThreadRank];
        //PosYList = new double[EventsNumPerThreadRank];
        //PosZList = new double[EventsNumPerThreadRank];
        
        G4double TotalProbability;
        for(int m = 0; m < RadioNuclideProbVec.size() ;m++ ){
            TotalProbability += RadioNuclideProbVec[m];
            //std::cout << m << " TotalProbability " << TotalProbability << " SpectrumOrDiscreteVec " << RadioNuclideSpectrumOrDiscreteVec[m] << " Proba " << RadioNuclideProbVec[m] << " EneMin " << RadioNuclideEneVec[m][0] << std::endl;
            //std::cout << m <<  " yield " << RadioNuclideProbVec[m] <<  " accumulated yield " << TotalProbability << std::endl;
        }
        
        for(int f = 0; f < EventsNumPerThreadRank ;f++ ){
            
            G4double RV = TotalProbability*(G4double)G4UniformRand() ;
            G4double EnergyCumulatedProbability = 0;
            
            for(int m = 0; m < RadioNuclideProbVec.size() ;m++ ){
                
                EnergyCumulatedProbability += RadioNuclideProbVec[m];
                
                if(RV < EnergyCumulatedProbability){
                    if(RadioNuclideSpectrumOrDiscreteVec[m] == 0){
                        
                        double RandomEne = RadioNuclideEneVec[m][0] + (RadioNuclideEneVec[m][1]-RadioNuclideEneVec[m][0])*(double)G4UniformRand();
                        EnergyList[f] = RandomEne;
                        //std::cout << f << " Particle "<< RadioNuclidePartNameVec[m] << " TotalProbability "<< TotalProbability << " RandomProb " << RV << " Prob "<< RadioNuclideProbVec[m] << " EnergyCumulatedProbability " << EnergyCumulatedProbability  << " SpectrumOrDiscreteVec " << RadioNuclideSpectrumOrDiscreteVec[m] << " Min " << RadioNuclideEneVec[m][0]   << " RandomEne " << EnergyList[f] << " Max " << RadioNuclideEneVec[m][1] << std::endl;
                        
                    }else{
                        
                        EnergyList[f] = RadioNuclideEneVec[m][0];
                        //std::cout << f << " Particle "<< RadioNuclidePartNameVec[m] << " TotalProbability "<< TotalProbability << " RandomProb " << RV << " Prob "<< RadioNuclideProbVec[m] << " EnergyCumulatedProbability " << EnergyCumulatedProbability  << " SpectrumOrDiscreteVec " << RadioNuclideSpectrumOrDiscreteVec[m] << " selected Energy " << EnergyList[f] << std::endl;
                    }
                    ParNameList[f]=RadioNuclidePartNameVec[m];
                    break;
                }
            }
        }
    }


    // for particles
    /*
    if(EnergyDistribution == "RadioNuclide" || EnergyDistribution == "File"){
        particleDefinitionList[0] = G4ParticleTable::GetParticleTable()->FindParticle("gamma");;
        particleDefinitionList[1] = G4ParticleTable::GetParticleTable()->FindParticle("e-");;
        particleDefinitionList[2] = G4ParticleTable::GetParticleTable()->FindParticle("e+");;
        particleDefinitionList[3] = G4ParticleTable::GetParticleTable()->FindParticle("alpha");;
        particleDefinitionList[4] = G4ParticleTable::GetParticleTable()->FindParticle("proton");;
        particleDefinitionList[5] = G4ParticleTable::GetParticleTable()->FindParticle("neutron");;
        particleDefinition = particleDefinitionList[0];
    }else{
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(NewRankSourceParticlesNamesValues[DataID]);

        //std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n DataID "<< DataID
        //          << " ParticleName=" << NewRankSourceParticlesNamesValues[DataID]
        //          << " ParticleName=" << particleDefinition->GetParticleName() << std::endl;
        if(particleDefinition == nullptr){ G4String msg = "Particle name (" + NewRankSourceParticlesNamesValues[DataID] + ") in rank/thread (" + std::to_string(DataID) + ") not found "; G4Exception("Source Data", "1", FatalErrorInArgument, msg.c_str());}
    }
    */


    //G4cout << " End of Initialization" << G4endl;

}

void G4TPrimaryGeneratorMethods::GunInitialize(){
    
    if(UseGeneratedData == "save"){
        std::cout << "\n For each thread, the direct data generation for each event simulation. The data will be saved to events data files " << std::endl;
    }else if(UseGeneratedData == "read"){
        std::cout << "\n For each thread, the simulated events data will be read from file, for each event " << std::endl;
    }else{
        std::cout << "\n For each thread, the direct data generation for each event simulation " << std::endl;
    }
    
    //std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n************** The primary generator Action initialization... "<< std::endl;
    
    DataID = 0;
#ifdef G4MPI_USE
    //DataID = G4MPImanager::GetManager()->GetRank() ;
    MPI_Comm_rank(MPI_COMM_WORLD, &DataID);
#else
    if(G4Threading::IsMultithreadedApplication()){
        DataID = G4Threading::G4GetThreadId();
    }
#endif
    
    EnergyListInc = 0;
    EvInc = 0;
    VoxelsInc = 0;
    TotalEmittedEnergy = 0. ;

    //G4cout << "\n Source Data before initialization "
    //          << " SourceParticleName:"<< NewRankSourceParticlesNamesValues[DataID]
    //          << " SourceEnergyDistribution:"<< EnergyDistribution << "(" << EnergyTypeNum << ")"
    //          << " SourcePositionsDistribution:"<< SourceType << "(" << PositionTypeNum << ")"
    //          << " SourceRegionName:"<< NewRankSourceRegionsNamesValues[DataID]
    //          << " SourceMomentumeDitribution:"<< NewRankSourceMomDirsValues[DataID] << "(" << MomDirTypeNum << ")"
    //          << G4endl;

    SourceInitialization();

    //G4cout << "\n Source Data after initialization "
    //          << " SourceParticleName:"<< NewRankSourceParticlesNamesValues[DataID]
    //          << " SourceEnergyDistribution:"<< EnergyDistribution << "(" << EnergyTypeNum << ")"
    //          << " SourcePositionsDistribution:"<< SourceType << "(" << PositionTypeNum << ")"
    //          << " SourceRegionName:"<< NewRankSourceRegionsNamesValues[DataID]
    //          << " SourceMomentumeDitribution:"<< NewRankSourceMomDirsValues[DataID] << "(" << MomDirTypeNum << ")"
    //          << G4endl;
}

void G4TPrimaryGeneratorMethods::FillDataMapsForNeutronCriticality(){

    EnergyListForCriticality[DataID]    = new double[EventsNumPerThreadRank];
    PositionsListForCriticality[DataID] = new G4ThreeVector[EventsNumPerThreadRank];
    MomDirecsListForCriticality[DataID] = new G4ParticleMomentum[EventsNumPerThreadRank];

    for(int f = 0; f < EventsNumPerThreadRank ;f++ ){
        GenerateEventsEnergy();
        GenerateEventsPosition();
        GenerateEventsMomentumDirection();
        EnergyListForCriticality[DataID][f]    = ENERGY;
        PositionsListForCriticality[DataID][f] = G4ThreeVector(X, Y, Z);
        MomDirecsListForCriticality[DataID][f] = G4ParticleMomentum(XMOMD, YMOMD, ZMOMD);

        //G4cout << NumberOfEventInBatch << " f  " <<f   << " " << ENERGY << " Pos " << G4ThreeVector(X, Y, Z) << " MomentumDir " << G4ParticleMomentum(XMOMD, YMOMD, ZMOMD) << G4endl;
    }
    //G4cout << "End Filling arrays" << G4endl;

}
void G4TPrimaryGeneratorMethods::GetEventDataFromCriticalityDataFile(){

    //G4cout << __FUNCTION__ << G4endl;

    // to read the file and fill the vector of points if it's empty, else it would be filled then we dont need to read and fill ...

    G4ThreadLocal G4int rank = 0 ;
    G4ThreadLocal G4int thread = 0;
    G4ThreadLocal G4int DataFileThreadRankID = 0; // used to read a new file with new ID and not 0 when we use "o" option

    DataID = 0;

#ifdef G4MPI_USE
    //rank = G4MPImanager::GetManager()->GetRank();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    DataID = rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &DataID);
#else
    if(G4Threading::IsMultithreadedApplication() || G4Threading::IsWorkerThread()){ // normal multiThreaded mode
        DataID = G4Threading::G4GetThreadId();
    }
#endif

    EnergyListForCriticality   [DataID] = new double[EventsNumPerThreadRank];
    PositionsListForCriticality[DataID] = new G4ThreeVector[EventsNumPerThreadRank];
    MomDirecsListForCriticality[DataID] = new G4ParticleMomentum[EventsNumPerThreadRank];

    // to test the data files if existes and if the file id correspond to the DataId if not use the file exists
    //#########################################################################################################################################"

    //G4ThreadLocal DIR *dir;
    //G4ThreadLocal struct dirent *diread;

    //if ((dir = opendir(DataDirectoryPath.c_str())) != nullptr) {
    //    while ((diread = readdir(dir)) != nullptr) {
    //        std::string str1 = diread->d_name;
    //        std::string Idc = str1.substr(0, 19);

    //        //std::cout << str1 << "  " << Idc << std::endl;

    //        std::vector<std::string> result;
    //        std::stringstream ss (getFileNameFromPath(diread->d_name));
    //        std::string item;
    //        while (getline (ss, item, '_')) {
    //            result.push_back(item);
    //            //std::cout << "item " << item << std::endl;
    //        }

    //        if(Idc == "CriticalityDataFile"){
    //            if(result.size() == 3){
    //                CriticalityDataFileNames[result[1]][std::atoi(result[2].c_str())] = DataDirectoryPath +"/"+diread->d_name;
    //                //std::cout << "diread->d_name " << diread->d_name << std::endl;
    //            }
    //        }
    //    }
    //    closedir (dir);
    //} else {
    //    perror ("opendir");
    //}
    G4ThreadLocal std::ostringstream c ;
    c << DataDirectoryPath  << "/CriticalityDataFile_"<< GeometrySymbol <<"_" <</*DataID*/0<< ".txt";//DataFilesExtension;

    //#########################################################################################################################################"

    //G4cout << " !!!!!!! If the number of generated points is less than the number of events to be processed, the application will reopens the file and reads the events data" << G4endl;
    std::cout << "\n---> for rank " << rank << ", thread " << thread << " (data ID:" << DataID << ") Reading " << EventsNumPerThreadRank << " Generated data from file " << c.str() << "..."<< std::endl;

    G4ThreadLocal G4double Ene, posX,posY,posZ, MX,MY,MZ;

    //std::ifstream fileorganStream(CriticalityDataFileNames[GeometrySymbol][DataID].c_str());
    G4ThreadLocal std::ifstream FileStream;
    FileStream.open(c.str(), std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
    bool isopen = true;

    if(!FileStream.is_open()) {

        G4ThreadLocal G4String msg = "Unable to open position data file : " + c.str() + " for rank/thread (" + std::to_string(DataID) + ") "; G4Exception("Source Data", "1", FatalErrorInArgument, msg.c_str());

        //if(CriticalityDataFileNames[GeometrySymbol].size() >0){
        //    FileStream.open(CriticalityDataFileNames[GeometrySymbol][0].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
        //    if(!FileStream.is_open()) {
        //        G4ThreadLocal G4String msg = "Unable to open position data file : " + CriticalityDataFileNames[GeometrySymbol][DataID] + " for rank/thread (" + std::to_string(DataID) + ") "; G4Exception("Source Data", "1", FatalErrorInArgument, msg.c_str());
        //        isopen = false;
        //    }else{
        //        G4ThreadLocal G4String msg = "Unable to open position data file for rank/thread (" + std::to_string(DataID) + ") but we have used the file of ID=0 "+CriticalityDataFileNames[GeometrySymbol][0]+" as a primary source data."; G4Exception("Source Data", "1", JustWarning, msg.c_str());
        //    }
        //}
    }

    if(isopen){

        G4ThreadLocal G4int valI; G4double valD;
        G4ThreadLocal G4String valS;

        //NumOfBatchs 10
        //NumOfEventsPerBatch 1000
        //NumOfFissionNeutrons 9780
        //k-eff 1.63815
        //NumOfTransportation 25182
        //NumOfhadElastic 624052
        //NumOfnCapture 5988
        //NumOfnFission 4002
        //NumOfneutronInelastic 2448
        //SimulatedNeutronSourceData E X Y Z MX MY MZ
        FileStream >> valS >> valI >> valS >> valI >> valS >> valI >> valS >> valI         ;
        //FissionCapturesOfThreadsRanks[0] = valI;
        FileStream >> valS >> valD ;
        while(FileStream >> valS){ // read just what we need to simulate
            if(valS == "SimulatedNeutronSourceData"){
                FileStream >> valS >> valS >> valS >> valS >> valS >> valS >> valS >> valS;
                break;
            }else{
                FileStream >> valI;
                //if(valS == "NumOfnFission"){
                //    //FissionCapturesOfThreadsRanks[2] = valI;
                //}else if(valS == "NumOfnCapture"){
                //    //FissionCapturesOfThreadsRanks[1]=valI;
                //}
            }
        }

        //G4cout << " DataID " << DataID << " FissionNuetrons " << FissionCapturesOfThreadsRanks[0]
        //          << " NumOfnCapture " << FissionCapturesOfThreadsRanks[1]
        //          << " NumOfnFission " << FissionCapturesOfThreadsRanks[2]
        //          << G4endl;

        G4ThreadLocal G4int EventIncre = 0;
        while(EventIncre != NumberOfEventInBatch){ // read just what we need to simulate in the first batch

            if(FileStream.peek()==EOF && EventIncre < NumberOfEventInBatch){ // If the number of generated points is less than the number of events to be processed, the application will reopen the file and reads points positions

                //std::cout << "EventIncre < NumberOfEventInBatch && peek()==EOF " << std::endl;
                FileStream.close(); // , std::ios_base::out | std::ios_base::binary
                //FileStream.open(CriticalityDataFileNames[GeometrySymbol][DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
                FileStream.open(c.str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
                //std::cout  << "##########################"<<NumberOfEventInBatch << " " <<EventIncre << " " << Energies[EventIncre] << " " << Positions[EventIncre] << " " << MomDirecs[EventIncre] << std::endl;
                //continue;

                FileStream >> valS >> valI >> valS >> valI >> valS >> valI >> valS >> valD ;
                while(FileStream >> valS){ // read just what we need to simulate
                    if(valS == "SimulatedNeutronSourceData"){
                        FileStream >> valS >> valS >> valS >> valS >> valS >> valS >> valS >> valS;
                        break;
                    }else{
                        FileStream >> valI;
                    }
                }
            }

            FileStream >> Ene >> posY >> posZ >> posY >> MY >> MZ >> MY ;

            EnergyListForCriticality   [DataID][EventIncre] = Ene ;
            PositionsListForCriticality[DataID][EventIncre] = G4ThreeVector( posX , posY , posZ);
            MomDirecsListForCriticality[DataID][EventIncre] = G4ParticleMomentum( MX , MY , MZ);


            //G4cout << NumberOfEventInBatch << " " <<EventIncre  << " " << Ene << " Pos " << G4ThreeVector( posX , posY , posZ) << " MomentumDir " << G4ParticleMomentum( MX , MY , MZ) << G4endl;

            EventIncre++;
        }

        FileStream.close();

        //std::cout << "End of Positions File reading"<< std::endl;
    }

    //###########################################posXYZ##############################################################################################"

    G4cout << "    ----> Getting Events Data has done. Now is time for Events simulation" << G4endl;

    //for ( int ff = 0 ; ff < NumberOfEventInBatch ; ff++ ) {
    //std::cout << "Energies[ff] " << Energies[ff]<< std::endl;
    //std::cout << "Positions[ff] " << Positions[ff]<< std::endl;
    //std::cout << "MomDirecs[ff] " << MomDirecs[ff]<< std::endl;
    //}

}

void G4TPrimaryGeneratorMethods::GetEventsData(){
    
    //G4cout << __FUNCTION__ << G4endl;
    
    // to read the file and fill the vector of points if it's empty, else it would be filled then we dont need to read and fill ...
    
    //G4int TotalEventsNumPerThreadRank = G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed();
    
    G4ThreadLocal G4int rank = 0 ;
    G4ThreadLocal G4int thread = 0;
    G4ThreadLocal G4int DataFileThreadRankID = 0; // used to read a new file with new ID and not 0 when we use "o" option
    
    DataID = 0;
    
#ifdef G4MPI_USE
    //rank = G4MPImanager::GetManager()->GetRank();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    DataFileThreadRankID = rank;
    
    if(MPISimulationNum == "m"){
        //DataID = G4MPImanager::GetManager()->GetRank();
        MPI_Comm_rank(MPI_COMM_WORLD, &DataID);
        rank = 0;
    }
#else
    if(G4Threading::IsMultithreadedApplication() || G4Threading::IsWorkerThread()){ // normal multiThreaded mode
        thread = G4Threading::G4GetThreadId();
        DataFileThreadRankID = thread;
        if(MPISimulationNum == "m"){
            DataID = thread;
            thread = 0;
        }
    }
#endif
    
    ReadedEnergyList  = new double[EventsNumPerThreadRank];
    ReadedPositionsList = new G4ThreeVector[EventsNumPerThreadRank];
    ReadedMomDirecsList = new G4ThreeVector[EventsNumPerThreadRank];
    
    // to test the data files if existes and if the file id correspond to the DataId if not use the file exists
    //#########################################################################################################################################"
    
    G4ThreadLocal DIR *dir;
    G4ThreadLocal struct dirent *diread;

    G4ThreadLocal std::map<G4String,std::map<G4String,int>> PosDataFileID;
    G4ThreadLocal std::map<G4String,std::map<double,int>> EneDataFileID;
    G4ThreadLocal std::map<G4String,int> MomDirDataFileID;
    
    if ((dir = opendir(DataDirectoryPath.c_str())) != nullptr) {
        while ((diread = readdir(dir)) != nullptr) {
            std::string str1 = diread->d_name;
            std::string Idc = str1.substr(0, 3);
            
            //std::cout << str1 << "  " << Idc << std::endl;
            
            std::vector<std::string> result;
            std::stringstream ss (getFileNameFromPath(diread->d_name));
            std::string item;
            while (getline (ss, item, '_')) {
                result.push_back(item);
                //std::cout << "item " << item << std::endl;
            }
            
            if(Idc == "Pos"){
                if(result.size() > 3){
                    PosDataFileNames[result[1]][result[2]][std::atoi(result[4].c_str())] = DataDirectoryPath +"/"+diread->d_name;
                    PosDataFileID[result[1]][result[2]] = std::atoi(result[4].c_str());
                    //std::cout << "diread->d_name " << diread->d_name << std::endl;
                }
            }
            else if(Idc == "Ene"){
                if(result.size() > 3){
                    EneDataFileNames[result[1]][std::atof(result[2].c_str())][std::atoi(result[4].c_str())] = DataDirectoryPath +"/"+diread->d_name;
                    EneDataFileID   [result[1]][std::atof(result[2].c_str())] = std::atoi(result[4].c_str()) ;
                }
            }
            else if(Idc == "Mom"){
                if(result.size() > 2){
                    MomDirDataFileNames[result[1]][std::atoi(result[3].c_str())] = DataDirectoryPath +"/"+diread->d_name;
                    MomDirDataFileID[result[1]] = std::atoi(result[3].c_str()) ;
                }
            }
        }
        closedir (dir);
    } else {
        perror ("opendir");
    }
    
    if(FILE* file = fopen(PosDataFileNames[NewRankSourceRegionsNamesValues[DataID]][SourceType][DataFileThreadRankID].c_str(),"r")){
        NewRankSourcePositionDataFiles[DataID] = PosDataFileNames[NewRankSourceRegionsNamesValues[DataID]][SourceType][DataFileThreadRankID]; std::fclose(file);
    }else{NewRankSourcePositionDataFiles[DataID] = PosDataFileNames[NewRankSourceRegionsNamesValues[DataID]][SourceType][PosDataFileID[NewRankSourceRegionsNamesValues[DataID]][SourceType]];}
    
    if(FILE* file = fopen(EneDataFileNames[EnergyDistribution][NewRankSourceEnergiesValues[DataID]][DataFileThreadRankID].c_str(),"r")){
        NewRankSourceEnergyDataFiles  [DataID] = EneDataFileNames[EnergyDistribution][NewRankSourceEnergiesValues[DataID]][DataFileThreadRankID]; std::fclose(file);
    }else{NewRankSourceEnergyDataFiles  [DataID] = EneDataFileNames[EnergyDistribution][NewRankSourceEnergiesValues[DataID]][EneDataFileID[EnergyDistribution][NewRankSourceEnergiesValues[DataID]]];}
    
    if(FILE* file = fopen(MomDirDataFileNames[MomDirDistribution][DataFileThreadRankID].c_str(),"r")){
        NewRankSourceMomDirDataFiles  [DataID] = MomDirDataFileNames[MomDirDistribution][DataFileThreadRankID];std::fclose(file);
    }else{NewRankSourceMomDirDataFiles  [DataID] = MomDirDataFileNames[MomDirDistribution][MomDirDataFileID[MomDirDistribution]];}
    
    //std::cout <<" DataFileThreadRankID " << DataFileThreadRankID << " DataID " << DataID << " NewRankSourcePositionDataFiles[DataID] " << NewRankSourcePositionDataFiles[DataID] << " " << PosDataFileNames[NewRankSourceRegionsNamesValues[DataID]][SourceType][DataID].c_str() << " " << PosDataFileID[NewRankSourceRegionsNamesValues[DataID]][SourceType] << std::endl;
    //std::cout <<" DataFileThreadRankID " << DataFileThreadRankID << " DataID " << DataID << " NewRankSourceEnergyDataFiles[DataID]   " << NewRankSourceEnergyDataFiles[DataID]   << " " << EneDataFileNames[EnergyDistribution][NewRankSourceEnergiesValues[DataID]][DataID].c_str()  << " " << EneDataFileID[EnergyDistribution][NewRankSourceEnergiesValues[DataID]] << std::endl;
    //std::cout <<" DataFileThreadRankID " << DataFileThreadRankID << " DataID " << DataID << " NewRankSourceMomDirDataFiles[DataID]   " << NewRankSourceMomDirDataFiles[DataID]   << " " << MomDirDataFileNames[MomDirDistribution][DataID].c_str() << " " << MomDirDataFileID[MomDirDistribution] << std::endl;
    
    //#########################################################################################################################################"
    
    //G4cout << " !!!!!!! If the number of generated points is less than the number of events to be processed, the application will reopens the file and reads the events data" << G4endl;
    std::cout << "\n---> for rank " << rank << "; Data ID:" << DataID << ", Reading " << EventsNumPerThreadRank << " Generated data from files " << NewRankSourcePositionDataFiles[DataID].c_str() << ", " << NewRankSourceEnergyDataFiles[DataID].c_str() << " and " << NewRankSourceMomDirDataFiles[DataID].c_str() << "..."<< std::endl;
    
    G4double posX,posY,posZ, val;
    
    //std::ifstream fileorganStream(NewRankSourcePositionDataFiles[DataID].c_str());
    std::ifstream PosFileStream;
    PosFileStream.open(NewRankSourcePositionDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
    
    if(!PosFileStream.is_open()) {
        G4String msg = "Unable to open position data file : " + NewRankSourcePositionDataFiles[DataID] + " for rank/thread (" + std::to_string(DataID) + ") "; G4Exception("Source Data", "1", FatalErrorInArgument, msg.c_str());
    }
    else {
        G4int Begin, First = 0;
        Begin = (rank*ThreadsNumber*EventsNumPerThreadRank) + (thread*EventsNumPerThreadRank);
        
        while(First < Begin){ // read just what we need to simulate
            
            PosFileStream >> posY >> posZ >> posY ;
            
            if(PosFileStream.peek()==EOF && First < Begin){
                PosFileStream.close(); // , std::ios_base::out | std::ios_base::binary
                PosFileStream.open(NewRankSourcePositionDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
            }
            
            First++;
        }
        std::cout << "    EventsNumPerThreadRank "<< EventsNumPerThreadRank << " Events Position Data For rank " << rank << ", thread " << thread << " (data ID:"<< DataID<< ") | Begins from line " << Begin << " to " << Begin + EventsNumPerThreadRank << std::endl;
        
        G4int EventIncre = 0;
        while(EventIncre < EventsNumPerThreadRank){ // read just what we need to simulate
            
            if(PosFileStream.peek()==EOF && EventIncre < EventsNumPerThreadRank){ // If the number of generated points is less than the number of events to be processed, the application will reopen the file and reads points positions
                
                //std::cout << "EventIncre < EventsNumPerThreadRank && peek()==EOF " << std::endl;
                PosFileStream.close(); // , std::ios_base::out | std::ios_base::binary
                PosFileStream.open(NewRankSourcePositionDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
                //std::cout  << "##########################"<<EventsNumPerThreadRank << " " <<EventIncre << " " << Energies[EventIncre] << " " << Positions[EventIncre] << " " << MomDirecs[EventIncre] << std::endl;
                //continue;
            }
            
            PosFileStream >> posX >> posY >> posZ ;
            ReadedPositionsList[EventIncre] = G4ThreeVector( posX , posY , posZ);
            
            //if( G4Threading::G4GetThreadId() == 3 && rank == 3){
            //  std::cout << EventsNumPerThreadRank << " " << EventIncre << " " << Positions[EventIncre] << std::endl;
            //}

            //std::cout <<EventsNumPerThreadRank << " " <<EventIncre  << " Pos " << Positions[EventIncre] << std::endl;
            
            EventIncre++;
        }
        
        PosFileStream.close();
        
        //std::cout << "End of Positions File reading"<< std::endl;
    }
    
    //###########################################posXYZ##############################################################################################"
    
    std::ifstream EneFileStream;
    EneFileStream.open(NewRankSourceEnergyDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
    
    if(!EneFileStream.is_open()) {
        G4String msg =  "Unable to open energy data file : " + NewRankSourceEnergyDataFiles[DataID] + " for rank/thread (" + std::to_string(DataID) + ") "; G4Exception("Source Data", "1", FatalErrorInArgument, msg.c_str());
    }
    else {
        
        G4int Begin, First = 0;
        Begin = (rank*ThreadsNumber*EventsNumPerThreadRank) + (thread*EventsNumPerThreadRank);
        
        while(First < Begin){ // read just what we need to simulate
            
            EneFileStream >> val ;
            
            if(EneFileStream.peek()==EOF && First < Begin){
                EneFileStream.close();
                EneFileStream.open(NewRankSourceEnergyDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
            }
            
            First++;
        }
        std::cout << "    EventsNumPerThreadRank "<< EventsNumPerThreadRank << " Events Energy Data For rank " << rank << ", thread " << thread << " (data ID:"<< DataID<< ") | Begins from line " << Begin << " to " << Begin + EventsNumPerThreadRank << std::endl;
        
        G4int EventIncre = 0;
        while(EventIncre < EventsNumPerThreadRank){ // read just what we need to simulate
            
            if(EneFileStream.peek()==EOF && EventIncre < EventsNumPerThreadRank){ // If the number of generated points is less than the number of events to be processed, the application will reopen the file and reads points positions
                
                //std::cout << "EventIncre < EventsNumPerThreadRank && peek()==EOF " << std::endl;
                EneFileStream.close();
                EneFileStream.open(NewRankSourceEnergyDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
                //std::cout  << "##########################"<<EventsNumPerThreadRank << " " <<EventIncre << " " << Energies[EventIncre] << " " << Positions[EventIncre] << " " << MomDirecs[EventIncre] << std::endl;
                //continue;
            }
            
            EneFileStream >> val;
            ReadedEnergyList[EventIncre] = val ;
            
            //if( G4Threading::G4GetThreadId() == 3  && rank == 3){
            //  std::cout << EventsNumPerThreadRank << " " << EventIncre << " " << Energies[EventIncre] << std::endl;
            //}
            
            //std::cout <<EventsNumPerThreadRank << " " <<EventIncre << " " << Energies[EventIncre] << " " << Positions[EventIncre] << " " << MomDirecs[EventIncre] << std::endl;
            
            EventIncre++;
        }
        
        EneFileStream.close();
        
        //std::cout << "End of Energy File reading"<< std::endl;
    }
    
    //#########################################################################################################################################"
    
    std::ifstream MomDirFileStream;
    MomDirFileStream.open(NewRankSourceMomDirDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
    
    if(!MomDirFileStream.is_open()) {
        G4String msg =  "Unable to open momentum direction data file : " + NewRankSourceMomDirDataFiles[DataID] + " for rank/thread (" + std::to_string(DataID) + ") "; G4Exception("Source Data", "1", FatalErrorInArgument, msg.c_str());
    }
    else {
        
        G4int Begin, First = 0;
        Begin = (rank*ThreadsNumber*EventsNumPerThreadRank) + (thread*EventsNumPerThreadRank);
        
        while(First < Begin){ // read just what we need to simulate
            
            MomDirFileStream >> posX >> posY >> posZ ;
            
            if(MomDirFileStream.peek()==EOF && First < Begin){
                MomDirFileStream.close();
                MomDirFileStream.open(NewRankSourceMomDirDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
            }
            
            First++;
        }
        
        std::cout << "    EventsNumPerThreadRank "<< EventsNumPerThreadRank << " Events Momentum Direction Data For rank " << rank << ", thread " << thread << " (data ID:"<< DataID<< ") | Begins from line " << Begin << " to " << Begin + EventsNumPerThreadRank << std::endl;
        
        G4int EventIncre = 0;
        while(EventIncre < EventsNumPerThreadRank){ // read just what we need to simulate
            
            if(MomDirFileStream.peek()==EOF && EventIncre < EventsNumPerThreadRank){ // If the number of generated points is less than the number of events to be processed, the application will reopen the file and reads points positions
                
                //std::cout << "EventIncre < EventsNumPerThreadRank && peek()==EOF " << std::endl;
                MomDirFileStream.close();
                MomDirFileStream.open(NewRankSourceMomDirDataFiles[DataID].c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary
                //std::cout  <<EventsNumPerThreadRank << " " <<EventIncre << " " << Energies[EventIncre] << " " << Positions[EventIncre] << " " << MomDirecs[EventIncre] << std::endl;
                //continue;
            }
            
            MomDirFileStream >> posX >> posY >> posZ ;
            ReadedMomDirecsList[EventIncre] = G4ThreeVector( posX , posY , posZ);
            
            //if( G4Threading::G4GetThreadId() == 3  && rank == 3){
            //  std::cout << EventsNumPerThreadRank << " " << EventIncre << " " << MomDirecs[EventIncre] << std::endl;
            //}
            
            EventIncre++;
        }
        
        MomDirFileStream.close();
        
        //std::cout << "End of MomDir File reading"<< std::endl;
    }
    
    std::cout << "    ----> Getting Events Data has done. Now is time for Events simulation" << std::endl;
    
    //for ( int ff = 0 ; ff < EventsNumPerThreadRank ; ff++ ) {
    //std::cout << "Energies[ff] " << Energies[ff]<< std::endl;
    //std::cout << "Positions[ff] " << Positions[ff]<< std::endl;
    //std::cout << "MomDirecs[ff] " << MomDirecs[ff]<< std::endl;
    //}
    
}
