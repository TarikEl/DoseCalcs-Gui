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

#include "G4TPrimaryGeneratorMessenger.hh"
#include "G4TVolumeConstruction.hh"
//#include "G4TPrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"

#include "G4SystemOfUnits.hh"

#include "G4UIcmdWithoutParameter.hh"

#include "globals.hh"
#include "G4Tokenizer.hh"

G4TPrimaryGeneratorMessenger::G4TPrimaryGeneratorMessenger(G4TVolumeConstruction* PGA):GeoConst(PGA){

    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@ " << __FUNCTION__ <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    sourceDataDir = new G4UIdirectory("/SourceData/");
    sourceDataDir->SetGuidance("Set type of particle and related data specifications");

    CommandsForPrimaryGen();

}

G4TPrimaryGeneratorMessenger::~G4TPrimaryGeneratorMessenger()
{

    delete sourceDataDir;
    delete source_typeCMD ;
    delete SourcePositionCMD ;
    delete SourceRotVector1CMD;
    delete SourceRotVector2CMD;

    delete SourceSolidCMD;
    delete SourceSurfaceCMD ;
    delete SourcePlaneCMD ;
    delete SourceAxisCMD ;
    delete SourceRotationCMD ;
    delete RadiusCMD ;
    delete HalfXCMD ;
    delete HalfYCMD ;
    delete HalfZCMD ;
    delete RadiusInCMD ;
    delete BeamSDevCMD ;

    delete SourceRegionNameCMD ;
    delete numberOfPointsToGenerateCMD ;
    delete UsePETCumulativeActCMD ;

    delete Vector_organ_sourceToGenCMD ;
    delete box_widthCMD;
    delete particle_nameCMD ;
    delete MonoEnergyCMD ;
    delete particle_angle_distributionCMD ;
    delete particle_energy_distributionCMD ;
    delete GaussDistSDevCMD ;
    delete GaussMeanCMD ;
    delete UniformDistEminCMD ;
    delete UniformDistEmaxCMD ;
    delete RayleighDistEmaxCMD ;
    delete showBoxCMD ;
    delete TestPointsPositionsCMD ;

    delete ThetaForDirectionCMD;
    delete PhiForDirectionCMD;

    delete ThetaMinCMD ;
    delete ThetaMaxCMD ;
    delete PhiMinCMD ;
    delete PhiMaxCMD ;

    delete GenerateMomDirsCMD;
    delete GenerateEnergiesCMD;
    delete GeneratePositionsCMD;

    delete SourceTypesCMD;
    delete SourceEnergyDataCMD;
    delete SourceMomDirDataCMD;
    delete GenerateSourceDataCMD;

    delete GeneratePosForRegionAndBoxDimCMD;
    delete GenerateEnergiesForCMD;
    delete GenerateMomDirsForCMD;

}

G4double G4TPrimaryGeneratorMessenger::UseG4Units(G4String Unit){

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

    return  U;
}

void G4TPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,G4String newValue){

    // Source Types commands

    if( command == SourceTypesCMD)
    {
        G4Tokenizer next(newValue);
    }

    // For Energy dist

    if( command == SourceEnergyDataCMD)
    {
        G4Tokenizer next(newValue);

        G4String SET = next();
        GeoConst->setEnergyDistribution(SET);

        G4double A, B; G4String Un3;

        if(SET == "Mono"){
            A = StoD(next()); Un3 = next();
            GeoConst->setMonoEnergy(A*UseG4Units(Un3));
        }
        else if(SET == "Gauss"){
            A = StoD(next()), B = StoD(next()); Un3 = next();
            GeoConst->setGaussMean(A*UseG4Units(Un3));
            GeoConst->setGaussSDev(B*UseG4Units(Un3));
        }
        else if(SET == "Uniform"){
            A = StoD(next()), B = StoD(next()); Un3 = next();
            GeoConst->setUniformEmin(A*UseG4Units(Un3));
            GeoConst->setUniformEmax(B*UseG4Units(Un3));
        }
        else if(SET == "Rayleigh"){
            A = StoD(next()); Un3 = next();
            GeoConst->setRayleighEmax(A*UseG4Units(Un3));
        }

    }

    // For angle dist

    if( command == SourceMomDirDataCMD)
    {
        G4Tokenizer next(newValue);

        G4String SET = next();
        GeoConst->setInitialDirectionModel(SET);

        G4double A, B, C, D; G4String Un2;

        if(SET == "Isotropic"){
        }
        else if(SET == "Uniform"){
        }
        else if(SET == "Delimited"){
            A = StoD(next()), B = StoD(next()), C = StoD(next()), D = StoD(next());
            Un2 = next();
            GeoConst->setThetaMin(A*UseG4Units(Un2));
            GeoConst->setThetaMax(B*UseG4Units(Un2));
            GeoConst->setPhiMin(C*UseG4Units(Un2));
            GeoConst->setPhiMax(D*UseG4Units(Un2));
        }
        else if(SET == "Directed"){
            A = StoD(next()), B = StoD(next());
            Un2 = next();
            GeoConst->setTheta(A*UseG4Units(Un2));
            GeoConst->setPhi(B*UseG4Units(Un2));
        }
    }

    if( command == particle_nameCMD )
    {
        GeoConst->setParticleName(newValue);
    }

    if( command == SimulatedParticlesNamesCMD )
    {
        G4Tokenizer next(newValue);

        G4String Par = next();
        //std::cout<< " Dddddd " << SN <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
        while (!Par.empty()) {

            //std::cout<< " Dddddd " << Par <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;

            GeoConst->setParticleName(Par);
            GeoConst->setSourceParticlesNamesValues(Par);

            Par = next();
        }
    }

    if( command == GeneratePosForRegionAndBoxDimCMD)
    {
        G4Tokenizer next(newValue);

        G4String Un1 = next(), ST = next(), SN = next();
        G4double XPos = 0, YPos= 0, ZPos= 0;

        GeoConst->setSourceType(ST);

        //std::cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@ " << ST << std::endl;

        if(ST == "Volume" || ST == "Voxels" || ST == "TET"){
        }else {
            XPos = StoD(next()), YPos = StoD(next()), ZPos = StoD(next());
            GeoConst->setSourcePosition(G4ThreeVector(XPos*UseG4Units(Un1), YPos*UseG4Units(Un1), ZPos*UseG4Units(Un1)));
        }

        GeoConst->setSourceType(ST);
        GeoConst->setSourceRegionName(SN);

        G4double A, B, C;

        if(ST == "Point"){
            GeoConst->setSourceRegionsNamesValues(SN);

            GeoConst->setRotPosAxis(next());
            GeoConst->setRotTheta(StoD(next()));
        }
        else if(ST == "Beam"){
            GeoConst->setSourceRegionsNamesValues(SN);

            //G4String STT = next();
            GeoConst->setSourceAxis(next());
            A = StoD(next()); //Un1 = next();
            GeoConst->setBeamSDev(A*UseG4Units(Un1));

            GeoConst->setRotPosAxis(next());
            GeoConst->setRotTheta(StoD(next()));
        }
        else if(ST == "Plane"){
            GeoConst->setSourceRegionsNamesValues(SN);

            G4String STT = next();
            GeoConst->setSourcePlane(STT);
            GeoConst->setSourceAxis(next());

            if(STT == "Square"){
                A = StoD(next()); //Un1 = next();
                GeoConst->setHalfX(A*UseG4Units(Un1));
            }
            else if(STT == "Circle"){
                A = StoD(next()); //Un1 = next();
                GeoConst->setRadius(A*UseG4Units(Un1));

            }
            else if(STT == "Rectangle"){
                A = StoD(next()), B = StoD(next()) ; //Un1 = next();

                GeoConst->setHalfX(A*UseG4Units(Un1));
                GeoConst->setHalfY(B*UseG4Units(Un1));

            }
            else if(STT == "Ellipse"){
                A = StoD(next()), B = StoD(next()) ; //Un1 = next();

                GeoConst->setHalfX(A*UseG4Units(Un1));
                GeoConst->setHalfY(B*UseG4Units(Un1));
            }
            else if(STT == "Annulus"){
                A = StoD(next()); //Un1 = next();

            }

            GeoConst->setRotPosAxis(next());
            GeoConst->setRotTheta(StoD(next()));

        }
        else if(ST == "Surface"){
            GeoConst->setSourceRegionsNamesValues(SN);

            G4String STT = next();
            GeoConst->setSourceSurface(STT);

            if(STT == "Sphere"){
                A = StoD(next()); //Un1 = next();
                GeoConst->setRadius(A*UseG4Units(Un1));

            }
            else if(STT == "Ellipsoid"){
                A = StoD(next()); B = StoD(next()); C = StoD(next()) ; //Un1 = next();

                GeoConst->setHalfX(A*UseG4Units(Un1));
                GeoConst->setHalfY(B*UseG4Units(Un1));
                GeoConst->setHalfZ(C*UseG4Units(Un1));
            }
            else if(STT == "Cylinder"){
                A = StoD(next()); B = StoD(next()); //C = StoD(next()) ; //Un1 = next();

                GeoConst->setHalfX(A*UseG4Units(Un1));
                GeoConst->setHalfY(B*UseG4Units(Un1));
                //GeoConst->setHalfZ(C*UseG4Units(Un1));
            }
            GeoConst->setRotPosAxis(next());
            GeoConst->setRotTheta(StoD(next()));
        }
        else if(ST == "Solid"){

            GeoConst->setSourceRegionsNamesValues(SN);

            G4String STT = next();
            GeoConst->setSourceSolid(STT);

            if(STT == "Para"){

                A = StoD(next()), B = StoD(next()), C = StoD(next()) ; //Un1 = next();

                GeoConst->setHalfX(A*UseG4Units(Un1));
                GeoConst->setHalfY(B*UseG4Units(Un1));
                GeoConst->setHalfZ(C*UseG4Units(Un1));

            }
            else if(STT == "EllipticCylinder"){
                A = StoD(next()), B = StoD(next()), C = StoD(next()) ; //Un1 = next();

                GeoConst->setHalfX(A*UseG4Units(Un1));
                GeoConst->setHalfY(B*UseG4Units(Un1));
                GeoConst->setHalfZ(C*UseG4Units(Un1));
            }
            else if(STT == "Cylinder"){
                A = StoD(next()), B = StoD(next()) ; //Un1 = next();

                GeoConst->setRadius(A*UseG4Units(Un1));
                GeoConst->setHalfZ(B*UseG4Units(Un1));
            }
            else if(STT == "Sphere"){
                A = StoD(next()); //Un1 = next();

                GeoConst->setRadius(A*UseG4Units(Un1));
            }
            else if(STT == "Ellipsoid"){
                A = StoD(next()), B = StoD(next()), C = StoD(next()) ; //Un1 = next();

                GeoConst->setHalfX(A*UseG4Units(Un1));
                GeoConst->setHalfY(B*UseG4Units(Un1));
                GeoConst->setHalfZ(C*UseG4Units(Un1));
            }
            GeoConst->setRotPosAxis(next());
            GeoConst->setRotTheta(StoD(next()));
        }
        else if(ST == "Volume"){

            //std::cout<< " Dddddd " << SN <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
            while (!SN.empty()) {

                G4String stLower = SN; stLower.toLower();

                if(stLower=="allregions"){

                    GeoConst->setSourceRegionsNamesValues(stLower);

                    G4int ExceptNum = StoI(next());

                    //std::cout<< " RegionName " << stLower << " ExceptNum " << ExceptNum <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;

                    for(G4int ds = 0 ; ds < ExceptNum ; ds++){

                        G4String Org = next();
                        //std::cout<< " Org " << Org <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;

                        GeoConst->setSourceRegionsNamesToBeIgnoredValues(Org);
                    }

                    A = StoD(next()); B = StoD(next()); C = StoD(next()) ; //Un1 = next();
                    //std::cout<< " A " << A <<  " B " << B << " C " << C << " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;

                    GeoConst->setSourceRegionsBoxDimValues(G4ThreeVector(A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un1)));
                    GeoConst->setXYZOfBox(                 G4ThreeVector(0., 0., 0.), true);

                }

                if(stLower!="allregions"){
                    GeoConst->setSourceRegionsNamesValues(SN);

                    //std::cout<< " RegionName " << SN << " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
                    A = StoD(next()); B = StoD(next()); C = StoD(next()) ; //Un1 = next();
                    //std::cout<< " A " << A <<  " B " << B << " C " << C << " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
                    GeoConst->setSourceRegionsBoxDimValues(G4ThreeVector(A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un1)));
                    GeoConst->setXYZOfBox(                 G4ThreeVector(A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un1)), true);
                }
                //std::cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@ " << SN << " " << G4ThreeVector(A*UseG4Units(Un1), B*UseG4Units(Un1), C*UseG4Units(Un1))  << " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;

                SN = next();
            }

        }
        else if(ST == "Voxels" || "TET"){
            while (!SN.empty()) {

                G4String stLower = SN; stLower.toLower();
                if(stLower=="allregions"){
                    G4int ExceptNum = StoI(next());

                    for(G4int ds = 0 ; ds < ExceptNum ; ds++){
                        G4String Org = next();

                        GeoConst->setSourceRegionsNamesToBeIgnoredValues(Org);
                    }

                    GeoConst->setSourceRegionsNamesValues(stLower);

                }else{

                    GeoConst->setSourceRegionsNamesValues(SN);
                    //std::cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@ " << SN << std::endl;
                }
                SN = next();
            }
        }

    }

    if( command == GenerateEnergiesForCMD)
    {
        G4Tokenizer next(newValue);

        //G4int NumOfEne = StoI(next());
        G4String Un3 = next(), SET = next();
        GeoConst->setEnergyDistribution(SET);

        if(SET == "Gauss"){
            GeoConst->setGaussSDev(StoD(next())*UseG4Units(Un3));
        }
        else if(SET == "Uniform"){
            GeoConst->setUniformEmin(StoD(next())*UseG4Units(Un3));
        }
        else if(SET == "Mono"){
            //GeoConst->setMonoEnergy(StoD(next())*UseG4Units(Un3));
        }
        else if(SET == "Rayleigh"){
            //GeoConst->setRayleighEmax(StoD(next())*UseG4Units(Un3));
        }

        if(SET == "Spectrum"){

            G4double EnergyCaracterizing = StoD(next());
            G4String E = next();

            G4double MaxEnergy = 0;
            //std::cout<< " Dddddd " << SN <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
            while (!E.empty()) {

                G4String P = next();

                G4double Eval = StoD(E)*UseG4Units(Un3);
                G4double Pval = StoD(P);

                //std::cout<< " E = " << E << " P = "<< P <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;

                if(Eval != 0.){
                    if(MaxEnergy < Eval){
                        MaxEnergy = Eval;
                    }
                    if(Pval != 0.){
                        GeoConst->setEnergyValueProbability(Eval, Pval);
                    }
                }
                E = next();
            }

            GeoConst->setSpectrumMaxEnergy(MaxEnergy);
            GeoConst->setSourceEnergiesValues(EnergyCaracterizing);

        }
        else if(SET == "RadioNuclide"){

            G4String RPar = next();
            G4double MaxEnergy = 0;

            while (RPar=next()) {
                /*
                G4String S = next();
                unsigned int SpectrumOrDiscrete;
                unsigned int ParticleIndex;
                bool isNewParticle = true;
                double Min=0;

                if(S=="Spectrum"){SpectrumOrDiscrete=0;}else {SpectrumOrDiscrete=1;}

                G4String Energy = next();

                /*while (!Energy.empty()) {

                    G4double EVal = StoD(Energy)*UseG4Units(Un3);
                    G4double Pval = StoD(next());

                    double* MinMax;

                    if(S=="Spectrum"){
                        MinMax = new double [2];
                        MinMax[0] = Min;
                        MinMax[1] = EVal;
                        std::cout<< " Pval = " << Pval << " RPar = "<< RPar <<  " SpectrumOrDiscrete = "<< SpectrumOrDiscrete <<  " MinMax[0] = "<< MinMax[0]  <<  " MinMax[1] = "<< MinMax[1] <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
                    } else {
                        MinMax = new double [1];
                        MinMax[0] = EVal;
                        std::cout<< " Pval = " << Pval << " RPar = "<< RPar <<  " SpectrumOrDiscrete = "<< SpectrumOrDiscrete <<  " EVal = "<< EVal <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
                    }

                    if(RPar=="gamma"){ParticleIndex = 0;}
                    if(RPar=="e-"){ParticleIndex = 1;}
                    if(RPar=="e+"){ParticleIndex = 2;}
                    if(RPar=="alpha"){ParticleIndex = 3;}
                    if(RPar=="proton"){ParticleIndex = 4;}
                    if(RPar=="neutron"){ParticleIndex = 5;}

                    GeoConst->setProbabilityParticleNameEnergy(Pval,ParticleIndex,SpectrumOrDiscrete,MinMax);

                    isNewParticle == false;
                    Min = EVal;

                    if(MaxEnergy < EVal){MaxEnergy = EVal;}

                    Energy = next();
                    std::cout<< " Energy = " << Energy << std::endl;

                    if(Energy == "gamma" || Energy == "e-" || Energy == "e+" || Energy == "alpha" || Energy == "proton" || Energy == "neutron"){
                        RPar = Energy;
                        break;
                    }
                }
                if(Energy.empty()){
                    break;
                }
*/

                G4String S = next();
                G4int NumberOfVal = StoI(next());

                //std::cout <<" S= " << S << " NumberOfVal = " << NumberOfVal << std::endl;

                unsigned int SpectrumOrDiscrete;
                unsigned int ParticleIndex;
                bool isNewParticle = true;
                double Min=0;

                if(S=="Spectrum"){SpectrumOrDiscrete=0;}else {SpectrumOrDiscrete=1;}

                for(G4int ds = 0 ; ds < NumberOfVal ; ds++){

                    G4double EVal = StoD(next())*UseG4Units(Un3);
                    G4double Pval = StoD(next());

                    double* MinMax;

                    if(S=="Spectrum"){
                        MinMax = new double [2];
                        MinMax[0] = Min;
                        MinMax[1] = EVal;
                        //std::cout << NumberOfVal << " Pval = " << Pval << " RPar = "<< RPar <<  " SpectrumOrDiscrete = "<< SpectrumOrDiscrete <<  " MinMax[0] = "<< MinMax[0]  <<  " MinMax[1] = "<< MinMax[1] <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
                    } else {
                        MinMax = new double [1];
                        MinMax[0] = EVal;
                        //std::cout << NumberOfVal << " Pval = " << Pval << " RPar = "<< RPar <<  " SpectrumOrDiscrete = "<< SpectrumOrDiscrete <<  " EVal = "<< EVal <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
                    }

                    if(RPar=="gamma"){ParticleIndex = 0;}
                    if(RPar=="e-"){ParticleIndex = 1;}
                    if(RPar=="e+"){ParticleIndex = 2;}
                    if(RPar=="alpha"){ParticleIndex = 3;}
                    if(RPar=="proton"){ParticleIndex = 4;}
                    if(RPar=="neutron"){ParticleIndex = 5;}

                    GeoConst->setProbabilityParticleNameEnergy(Pval,ParticleIndex,SpectrumOrDiscrete,MinMax);

                    isNewParticle == false;
                    Min = EVal;

                    if(MaxEnergy < EVal){MaxEnergy = EVal;}
                    //std::cout<< " Energy = " << Energy << std::endl;
                }

                RPar = next();
                if(RPar.empty()){ break;}

            }
            GeoConst->setRadioNuclideMaxEnergy(MaxEnergy);
            GeoConst->setSourceEnergiesValues(MaxEnergy);
        }
        else if(SET == "File"){

            G4double MaxEnergy = 0;
            G4double EnergyCaracterizing = StoD(next());

            G4String filename = next();
            std::ifstream file1(filename.c_str() , std::ios::binary);
            if(file1.is_open()){

                unsigned int SpectrumOrDiscrete;
                unsigned int ParticleIndex;
                bool isNewParticle = true;
                double Min=0;

                G4String RPar;
                G4String Word;
                G4String S;

                while (file1 >> Word) {

                    if(Word == "gamma" || Word == "e-" || Word == "e+" || Word == "alpha" || Word == "proton" || Word == "neutron"){
                        RPar = Word;
                        file1 >> S;
                        if(S=="Spectrum"){SpectrumOrDiscrete=0;}else {SpectrumOrDiscrete=1;}
                        isNewParticle = true;
                        Min=0;

                        //std::cout<< " RPar = " << RPar << " S = " << S << std::endl;

                        continue;
                    }

                    //std::cout<< " Word = " << Word << std::endl;

                    G4double EVal = StoD(Word); //Pval= StoD(Pval); EVal = EVal*UseG4Units(Un3);
                    G4double Pval ; file1 >> Pval;

                    double* MinMax;

                    if(S=="Spectrum"){
                        MinMax = new double [2];
                        MinMax[0] = Min;
                        MinMax[1] = EVal;
                        //std::cout<< " Pval = " << Pval << " RPar = "<< RPar <<  " SpectrumOrDiscrete = "<< SpectrumOrDiscrete <<  " MinMax[0] = "<< MinMax[0]  <<  " MinMax[1] = "<< MinMax[1] <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
                    }else {
                        MinMax = new double [1];
                        MinMax[0] = EVal;
                        //std::cout<< " Pval = " << Pval << " RPar = "<< RPar <<  " SpectrumOrDiscrete = "<< SpectrumOrDiscrete <<  " EVal = "<< EVal <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
                    }

                    if(RPar=="gamma"){ParticleIndex = 0;}
                    if(RPar=="e-"){ParticleIndex = 1;}
                    if(RPar=="e+"){ParticleIndex = 2;}
                    if(RPar=="alpha"){ParticleIndex = 3;}
                    if(RPar=="proton"){ParticleIndex = 4;}
                    if(RPar=="neutron"){ParticleIndex = 5;}

                    GeoConst->setProbabilityParticleNameEnergy(Pval,ParticleIndex,SpectrumOrDiscrete,MinMax);

                    isNewParticle == false;
                    Min = EVal;

                    if(MaxEnergy < EVal){MaxEnergy = EVal;}
                }
                file1.close();
            }

            GeoConst->setFileEnergyCharacterizer(EnergyCaracterizing);
            GeoConst->setSourceEnergiesValues(EnergyCaracterizing);
        }
        else{
            G4String E = next();
            //std::cout<< " Dddddd " << SN <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
            while (!E.empty()) {

                //std::cout<< " Dddddd " << E <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;

                G4double Eval = StoD(E)*UseG4Units(Un3);

                if(Eval != 0.){
                    GeoConst->setSourceEnergiesValues(Eval);
                }

                E = next();
            }
        }

        //std::cout<< " Dddddd hyyyyyyyyyyyyyyyyy @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;

    }

    if( command == GenerateMomDirsForCMD)
    {
        G4Tokenizer next(newValue);

        G4String Un2 = next(), SET = next();
        //std::cout<< " Dddddd " << SN <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;

        GeoConst->setInitialDirectionModel(SET);

        while (!SET.empty()) {

            //std::cout<< " Dddddd " << SET <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;

            //GeoConst->setSourceMomDirsValues(SET);
            //if(SET == "Directed"){
            //    G4double T = StoD(next()), P = StoD(next());
            //    GeoConst->setSourceMomDirsDirectedThetaValues(T*UseG4Units(Un2));
            //    GeoConst->setSourceMomDirsDirectedPhiValues(P*UseG4Units(Un2));
            //    //std::cout<< "MomDir Dist " << SET <<  " : "<< T << " " << P << std::endl;

            //}else {
            //    GeoConst->setSourceMomDirsDirectedThetaValues(0);
            //    GeoConst->setSourceMomDirsDirectedPhiValues(0);
            //    //std::cout<< "MomDir Dist " << SET <<  " : "<< 0 << " " << 0 << std::endl;

            //}

            GeoConst->setSourceMomDirsValues(SET);
            if(SET == "Directed"){
                G4String DirTo = next();
                if(DirTo == "ToPoint"){
                    GeoConst->setMomDirDirectedHow(DirTo);

                    if(Un2 == "rad" || Un2 == "degree"){Un2=="mm";}

                    GeoConst->setDirectedToX(StoD(next())*UseG4Units(Un2));
                    GeoConst->setDirectedToY(StoD(next())*UseG4Units(Un2));
                    GeoConst->setDirectedToZ(StoD(next())*UseG4Units(Un2));

                    GeoConst->setSourceMomDirsDirectedThetaValues(0);
                    GeoConst->setSourceMomDirsDirectedPhiValues(0);
                }
                else if(DirTo == "ParallelTo"){
                    GeoConst->setMomDirDirectedHow(DirTo);

                    if(Un2 == "rad" || Un2 == "degree"){Un2=="mm";}

                    GeoConst->setDirectedParallelAxis(next());
                    GeoConst->setDirectedToX(StoD(next())*UseG4Units(Un2));
                    GeoConst->setDirectedToY(StoD(next())*UseG4Units(Un2));

                    GeoConst->setSourceMomDirsDirectedThetaValues(0);
                    GeoConst->setSourceMomDirsDirectedPhiValues(0);
                }
                else if(DirTo == "ToVolume"){
                    GeoConst->setMomDirDirectedHow(DirTo);

                    if(Un2 == "rad" || Un2 == "degree"){Un2=="mm";}

                    GeoConst->setToVolumeX(StoD(next())*UseG4Units(Un2));
                    GeoConst->setToVolumeY(StoD(next())*UseG4Units(Un2));
                    GeoConst->setToVolumeZ(StoD(next())*UseG4Units(Un2));                    GeoConst->setDirectedToX(StoD(next())*UseG4Units(Un2));
                    GeoConst->setDirectedToX(StoD(next())*UseG4Units(Un2));
                    GeoConst->setDirectedToY(StoD(next())*UseG4Units(Un2));
                    GeoConst->setDirectedToZ(StoD(next())*UseG4Units(Un2));

                    GeoConst->setSourceMomDirsDirectedThetaValues(0);
                    GeoConst->setSourceMomDirsDirectedPhiValues(0);
                }
                else if(DirTo == "ThetaPhi"){
                    GeoConst->setMomDirDirectedHow(DirTo);

                    if(Un2 == "rad" || Un2 == "degree"){Un2=="degree";}

                    //std::cout<< "MomDir Dist " << SET <<  " : "<< T << " " << P << std::endl;

                    G4double T = StoD(next()), P = StoD(next());
                    GeoConst->setSourceMomDirsDirectedThetaValues(T*UseG4Units(Un2));
                    GeoConst->setSourceMomDirsDirectedPhiValues(P*UseG4Units(Un2));
                }
                //std::cout<< "MomDir Dist " << SET <<  " : "<< T << " " << P << std::endl;

            }else {
                GeoConst->setSourceMomDirsDirectedThetaValues(0);
                GeoConst->setSourceMomDirsDirectedPhiValues(0);
                //std::cout<< "MomDir Dist " << SET <<  " : "<< 0 << " " << 0 << std::endl;

            }




            SET = next();
        }
    }

    if( command == GenerateSourceDataCMD)
    {
        G4Tokenizer next(newValue);

        G4String first = next();

        if(first == "save"){
            GeoConst->setUseGeneratedData(first);
        }
        else if(first == "read"){
            GeoConst->setUseGeneratedData(first);
        }
        else if (first == "generate"){
            GeoConst->setUseGeneratedData(first);
            //GeoConst->setPointNumberToSave(StoI(next()));
            GeoConst->setGeneratePositions("yes");//next()
            GeoConst->setGenerateEnergies("yes");//next()
            GeoConst->setGenerateMomDirs("yes");//next()
        }
    }

    if( command == showBoxCMD )
    {
        GeoConst->setShowBox("yes");
    }
    if( command == TestPointsPositionsCMD )
    {
        GeoConst->setTestPointsPositions("yes");
    }

    if( command == source_typeCMD )
    {
        GeoConst->setSourceType(newValue);
    }
    if( command == SourcePositionCMD )
    {
        GeoConst->setSourcePosition(SourcePositionCMD->GetNew3VectorValue(newValue));
    }
    if( command == SourceRotVector1CMD )
    {
        GeoConst->setSourceRotVector1(SourceRotVector1CMD->GetNew3VectorValue(newValue));
    }
    if( command == SourceRotVector2CMD )
    {
        GeoConst->setSourceRotVector2(SourceRotVector2CMD->GetNew3VectorValue(newValue));
    }
    if( command == SourceSolidCMD )
    {
        GeoConst->setSourceSolid(newValue);
    }
    if( command == SourceSurfaceCMD )
    {
        GeoConst->setSourceSurface(newValue);
    }
    if( command == SourcePlaneCMD )
    {
        GeoConst->setSourcePlane(newValue);
    }
    if( command == SourceAxisCMD )
    {
        GeoConst->setSourceAxis(newValue);
    }
    if( command == SourceRotationCMD )
    {
        GeoConst->setSourceRotation(SourceRotationCMD->GetNew3VectorValue(newValue));
    }
    if( command == HalfZCMD )
    {
        GeoConst->setHalfZ(HalfZCMD->GetNewDoubleValue(newValue));
    }
    if( command == HalfXCMD )
    {
        GeoConst->setHalfX(HalfXCMD->GetNewDoubleValue(newValue));
    }
    if( command == HalfYCMD )
    {
        GeoConst->setHalfY(HalfYCMD->GetNewDoubleValue(newValue));
    }
    if( command == RadiusCMD )
    {
        GeoConst->setRadius(RadiusCMD->GetNewDoubleValue(newValue));
    }
    if( command == RadiusInCMD )
    {
        GeoConst->setRadiusIn(RadiusInCMD->GetNewDoubleValue(newValue));
    }
    if( command == BeamSDevCMD )
    {
        GeoConst->setBeamSDev(BeamSDevCMD->GetNewDoubleValue(newValue));
    }


    if( command == GaussMeanCMD )
    {
        GeoConst->setGaussMean(GaussMeanCMD->GetNewDoubleValue(newValue));
    }
    if( command == GaussDistSDevCMD )
    {
        GeoConst->setGaussSDev(GaussDistSDevCMD->GetNewDoubleValue(newValue));
    }
    if( command == UniformDistEminCMD )
    {
        GeoConst->setUniformEmin(UniformDistEminCMD->GetNewDoubleValue(newValue));
    }
    if( command == UniformDistEmaxCMD )
    {
        GeoConst->setUniformEmax(UniformDistEmaxCMD->GetNewDoubleValue(newValue));
    }
    if( command == RayleighDistEmaxCMD )
    {
        GeoConst->setRayleighEmax(RayleighDistEmaxCMD->GetNewDoubleValue(newValue));
    }

    if( command == GeneratePositionsCMD )
    {
        GeoConst->setGeneratePositions("yes");
    }
    if( command == GenerateEnergiesCMD )
    {
        GeoConst->setGenerateEnergies("yes");
    }
    if( command == GenerateMomDirsCMD )
    {
        GeoConst->setGenerateMomDirs("yes");
    }
    if( command == SourceRegionNameCMD )
    {
        GeoConst->setSourceRegionName(newValue);
    }
    if( command == box_widthCMD )
    {
        if(box_widthCMD->GetNew3VectorValue(newValue) == G4ThreeVector()){
            GeoConst->setXYZOfBox(box_widthCMD->GetNew3VectorValue(newValue), false);
        }else{
            GeoConst->setXYZOfBox(box_widthCMD->GetNew3VectorValue(newValue), true);
        }
    }
    if( command == Vector_organ_sourceToGenCMD )
    {
        GeoConst->setOrgPosThreeVector(Vector_organ_sourceToGenCMD->GetNew3VectorValue(newValue));
    }
    if( command == numberOfPointsToGenerateCMD )
    {
        GeoConst->setPointNumberToSave(numberOfPointsToGenerateCMD->GetNewIntValue(newValue));
    }
    if( command == UsePETCumulativeActCMD )
    {
        GeoConst->setUsePETCumulativeAct("yes");
    }

    // the values will be sent to GeoConst and will be used by primary generator
    // source commands


    if( command == MonoEnergyCMD )
    {
        GeoConst->setMonoEnergy(MonoEnergyCMD->GetNewDoubleValue(newValue));
    }
    if( command == particle_angle_distributionCMD )
    {
        GeoConst->setInitialDirectionModel(newValue);
    }
    if( command == particle_energy_distributionCMD )
    {
        GeoConst->setEnergyDistribution(newValue);
    }
    if( command == ThetaForDirectionCMD )
    {
        GeoConst->setTheta(ThetaForDirectionCMD->GetNewDoubleValue(newValue));
    }
    if( command == PhiForDirectionCMD )
    {
        GeoConst->setPhi(PhiForDirectionCMD->GetNewDoubleValue(newValue));
    }
    if( command == ThetaMaxCMD )
    {
        GeoConst->setThetaMax(ThetaMaxCMD->GetNewDoubleValue(newValue));
    }
    if( command == ThetaMinCMD )
    {
        GeoConst->setThetaMin(ThetaMinCMD->GetNewDoubleValue(newValue));
    }
    if( command == PhiMaxCMD )
    {
        GeoConst->setPhiMax(PhiMaxCMD->GetNewDoubleValue(newValue));
    }
    if( command == PhiMinCMD )
    {
        GeoConst->setPhiMin(PhiMinCMD->GetNewDoubleValue(newValue));
    }

}

void  G4TPrimaryGeneratorMessenger::CommandsForPrimaryGen(){

    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@ " << __FUNCTION__ <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    G4UIparameter* param;

    SourceTypesCMD = new G4UIcommand("/SourceData/setSourcePosData" ,this);
    param = new G4UIparameter("Source Name",'s', false);      SourceTypesCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Source Type",'s', false);      SourceTypesCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    for(G4int ds = 0 ; ds < 20 ; ds++){

        G4String f = "Source Type Parameter " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      SourceTypesCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }


    SourceEnergyDataCMD = new G4UIcommand("/SourceData/setSourceEnergyData" ,this);
    param = new G4UIparameter("Source Energy Distribution ",'s', false);      SourceEnergyDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    for(G4int ds = 0 ; ds < 1000 ; ds++){

        G4String f = "Source Energy Distribution Parameters " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      SourceEnergyDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }


    SourceMomDirDataCMD = new G4UIcommand("/SourceData/setSourceMomDirData" ,this);
    param = new G4UIparameter("Source Mom Dir Distribution ",'s', false);      SourceMomDirDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    for(G4int ds = 0 ; ds < 20 ; ds++){

        G4String f = "Source Mom Dir Distribution Parameters " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      SourceMomDirDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

    GenerateSourceDataCMD = new G4UIcommand("/SourceData/useDataGenerationFiles" ,this);
    param = new G4UIparameter("Generate Or Save ",'s', false);      GenerateSourceDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Generate Events Number ",'s', true);      GenerateSourceDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Generate Events Position ",'s', true);      GenerateSourceDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Generate Events Energies ",'s', true);      GenerateSourceDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Generate Events Mom Dirs ",'s', true);      GenerateSourceDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    //param = new G4UIparameter("Use Generated Data in Simulation ",'s', true);      GenerateSourceDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    //param = new G4UIparameter("Show Box ",'s', false);      GenerateSourceDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error


    GeneratePosForRegionAndBoxDimCMD = new G4UIcommand("/SourceData/setEventsInitialPosData" ,this);
    for(G4int ds = 0 ; ds < 100 ; ds++){
        G4String f = "Source Pos " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      GeneratePosForRegionAndBoxDimCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

    GenerateEnergiesForCMD = new G4UIcommand("/SourceData/setEventsInitialEneData" ,this);

    for(G4int ds = 0 ; ds < 100 ; ds++){
        G4String f = "Source Energy " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      GenerateEnergiesForCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

    GenerateMomDirsForCMD = new G4UIcommand("/SourceData/setEventsInitialMomDirData" ,this);

    for(G4int ds = 0 ; ds < 100 ; ds++){
        G4String f = "Source Mom Dir " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      GenerateMomDirsForCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

    SimulatedParticlesNamesCMD = new G4UIcommand("/SourceData/setEventsParticleNameData" ,this);

    for(G4int ds = 0 ; ds < 100 ; ds++){
        G4String f = "Source Particles names " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      SimulatedParticlesNamesCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

    source_typeCMD = new G4UIcmdWithAString("/SourceData/setSourceType",this);
    source_typeCMD->SetGuidance("");
    source_typeCMD->SetParameterName("source_type",true);
    source_typeCMD->SetDefaultValue("Point");
    source_typeCMD->SetCandidates("Rotated Point Beam Plane Surface Solid Volume Voxels");
    source_typeCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    SourceSurfaceCMD = new G4UIcmdWithAString("/SourceData/setSourceSurface",this);
    SourceSurfaceCMD->SetGuidance("");
    SourceSurfaceCMD->SetParameterName("SourceSurface",true);
    SourceSurfaceCMD->SetDefaultValue("Sphere");
    SourceSurfaceCMD->SetCandidates("Sphere Ellipsoid");
    SourceSurfaceCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    SourceSolidCMD = new G4UIcmdWithAString("/SourceData/setSourceSolid",this);
    SourceSolidCMD->SetGuidance("");
    SourceSolidCMD->SetParameterName("SourceSolidCMD",true);
    SourceSolidCMD->SetDefaultValue("Sphere");
    SourceSolidCMD->SetCandidates("Para EllipticCylinder Cylinder Sphere Ellipsoid");
    SourceSolidCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    SourcePlaneCMD = new G4UIcmdWithAString("/SourceData/setSourcePlane",this);
    SourcePlaneCMD->SetGuidance("");
    SourcePlaneCMD->SetParameterName("SourcePlaneCMD",true);
    SourcePlaneCMD->SetDefaultValue("Circle");
    SourcePlaneCMD->SetCandidates("Circle Square Rectangle Ellipse Annulus");
    SourcePlaneCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    SourceAxisCMD = new G4UIcmdWithAString("/SourceData/setSourceAxis",this);
    SourceAxisCMD->SetGuidance("");
    SourceAxisCMD->SetParameterName("SourceAxisCMD",true);
    SourceAxisCMD->SetDefaultValue("Z");
    SourceAxisCMD->SetCandidates("X Y Z");
    SourceAxisCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    SourceRotationCMD = new G4UIcmdWith3VectorAndUnit("/SourceData/setSourceRotation",this);
    SourceRotationCMD->SetGuidance("");
    SourceRotationCMD->SetParameterName("x","y","z",true,true);
    SourceRotationCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    SourceRotationCMD->SetDefaultValue(G4ThreeVector(20.,20.,20.));
    SourceRotationCMD->SetDefaultUnit("degree"); // then we dont need to add *cm in command
    SourceRotationCMD->SetUnitCandidates("degree");

    RadiusCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setSourceRadius",this);
    RadiusCMD->SetGuidance("");
    RadiusCMD->SetDefaultValue(1.);
    RadiusCMD->SetParameterName("RadiusCMD",true);
    RadiusCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //RadiusCMD->SetUnitCategory()
    RadiusCMD->SetDefaultUnit("mm");
    RadiusCMD->SetUnitCandidates("cm mm m");

    RadiusInCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setSourceRadiusIn",this);
    RadiusInCMD->SetGuidance("");
    RadiusInCMD->SetDefaultValue(1.);
    RadiusInCMD->SetParameterName("RadiusInCMD",true);
    RadiusInCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //RadiusInCMD->SetUnitCategory()
    RadiusInCMD->SetDefaultUnit("mm");
    RadiusInCMD->SetUnitCandidates("cm mm m");

    BeamSDevCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setSourceBeamSDev",this);
    BeamSDevCMD->SetGuidance("");
    BeamSDevCMD->SetDefaultValue(1.);
    BeamSDevCMD->SetParameterName("BeamSDevCMD",true);
    BeamSDevCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //BeamSDevCMD->SetUnitCategory()
    BeamSDevCMD->SetDefaultUnit("mm");
    BeamSDevCMD->SetUnitCandidates("cm mm m");

    HalfXCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setSourceHalfX",this);
    HalfXCMD->SetGuidance("");
    HalfXCMD->SetDefaultValue(1.);
    HalfXCMD->SetParameterName("HalfXCMD",true);
    HalfXCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //HalfXCMD->SetUnitCategory()
    HalfXCMD->SetDefaultUnit("mm");
    HalfXCMD->SetUnitCandidates("cm mm m");

    HalfYCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setSourceHalfY",this);
    HalfYCMD->SetGuidance("");
    HalfYCMD->SetDefaultValue(1.);
    HalfYCMD->SetParameterName("HalfYCMD",true);
    HalfYCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //HalfYCMD->SetUnitCategory()
    HalfYCMD->SetDefaultUnit("mm");
    HalfYCMD->SetUnitCandidates("cm mm m");

    HalfZCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setSourceHalfZ",this);
    HalfZCMD->SetGuidance("");
    HalfZCMD->SetDefaultValue(1.);
    HalfZCMD->SetParameterName("HalfZCMD",true);
    HalfZCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //HalfZCMD->SetUnitCategory()
    HalfZCMD->SetDefaultUnit("mm");
    HalfZCMD->SetUnitCandidates("cm mm m");

    showBoxCMD = new G4UIcmdWithoutParameter("/SourceData/showSourceBox",this);
    showBoxCMD->SetGuidance(" ");
    //showBoxCMD->SetParameterName("showBox",true);
    //showBoxCMD->SetDefaultValue("no");
    //showBoxCMD->SetCandidates("yes no");
    showBoxCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    TestPointsPositionsCMD = new G4UIcmdWithoutParameter("/SourceData/testEventsInitialPositions",this);
    TestPointsPositionsCMD->SetGuidance(" ");
    //TestPointsPositionsCMD->SetParameterName("TestPointsPositions",true);
    //TestPointsPositionsCMD->SetDefaultValue("no");
    //TestPointsPositionsCMD->SetCandidates("yes no");
    TestPointsPositionsCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    particle_nameCMD = new G4UIcmdWithAString("/SourceData/setParticleName",this);
    particle_nameCMD->SetGuidance("Set name of primary particles source : e- gamma proton ,... ");
    particle_nameCMD->SetParameterName("particle_name",true);
    particle_nameCMD->SetDefaultValue("gamma");
    particle_nameCMD->SetCandidates("e- gamma proton alpha neutron");
    particle_nameCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    MonoEnergyCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setMonoEnergy",this);
    MonoEnergyCMD->SetGuidance("Set energy initial in MeV of primary particles");
    MonoEnergyCMD->SetDefaultValue(1.);
    MonoEnergyCMD->SetParameterName("particle_initial_energy",true);
    MonoEnergyCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //MonoEnergyCMD->SetUnitCategory()
    MonoEnergyCMD->SetDefaultUnit("MeV");
    MonoEnergyCMD->SetUnitCandidates("eV keV MeV GeV");

    particle_angle_distributionCMD = new G4UIcmdWithAString("/SourceData/setAngleDistribution",this);
    particle_angle_distributionCMD->SetGuidance("Set direction model of primary particles : ");
    particle_angle_distributionCMD->SetParameterName("particle_energy_emission_direction",true);
    particle_angle_distributionCMD->SetDefaultValue("Uniform");
    particle_angle_distributionCMD->SetCandidates("Uniform Isotropic Directed Delimited");
    particle_angle_distributionCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    ThetaForDirectionCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setDirectionTheta",this);
    ThetaForDirectionCMD->SetGuidance(" in degree");
    ThetaForDirectionCMD->SetDefaultValue(45.);
    ThetaForDirectionCMD->SetParameterName("ThetaForDirection",true);
    ThetaForDirectionCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //ThetaForDirectionCMD->SetUnitCategory()
    ThetaForDirectionCMD->SetDefaultUnit("degree");
    ThetaForDirectionCMD->SetUnitCandidates("degree");

    ThetaMinCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setThetaMin",this);
    ThetaMinCMD->SetGuidance(" in degree");
    ThetaMinCMD->SetDefaultValue(45.);
    ThetaMinCMD->SetParameterName("ThetaMin",true);
    ThetaMinCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //ThetaMinCMD->SetUnitCategory()
    ThetaMinCMD->SetDefaultUnit("degree");
    ThetaMinCMD->SetUnitCandidates("degree");

    ThetaMaxCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setThetaMax",this);
    ThetaMaxCMD->SetGuidance(" in degree");
    ThetaMaxCMD->SetDefaultValue(45.);
    ThetaMaxCMD->SetParameterName("ThetaMax",true);
    ThetaMaxCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //ThetaMaxCMD->SetUnitCategory()
    ThetaMaxCMD->SetDefaultUnit("degree");
    ThetaMaxCMD->SetUnitCandidates("degree");

    PhiMinCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setPhiMin",this);
    PhiMinCMD->SetGuidance(" in degree");
    PhiMinCMD->SetDefaultValue(45.);
    PhiMinCMD->SetParameterName("PhiMin",true);
    PhiMinCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //PhiMinCMD->SetUnitCategory()
    PhiMinCMD->SetDefaultUnit("degree");
    PhiMinCMD->SetUnitCandidates("degree");

    PhiMaxCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setPhiMax",this);
    PhiMaxCMD->SetGuidance(" in degree");
    PhiMaxCMD->SetDefaultValue(45.);
    PhiMaxCMD->SetParameterName("PhiMax",true);
    PhiMaxCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //PhiMaxCMD->SetUnitCategory()
    PhiMaxCMD->SetDefaultUnit("degree");
    PhiMaxCMD->SetUnitCandidates("degree");


    PhiForDirectionCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setDirectionPhi",this);
    PhiForDirectionCMD->SetGuidance(" in degree");
    PhiForDirectionCMD->SetDefaultValue(45.);
    PhiForDirectionCMD->SetParameterName("PhiForDirection",true);
    PhiForDirectionCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    particle_energy_distributionCMD = new G4UIcmdWithAString("/SourceData/setEnergyDistribution",this);
    particle_energy_distributionCMD->SetGuidance("Set distribution name of primary particles : ");
    particle_energy_distributionCMD->SetParameterName("particle_energy_emission_direction",true);
    particle_energy_distributionCMD->SetDefaultValue("Mono");
    particle_energy_distributionCMD->SetCandidates("Mono Uniform Gauss Rayleigh");
    particle_energy_distributionCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    GaussMeanCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setGaussMean",this);
    GaussMeanCMD->SetGuidance("");
    GaussMeanCMD->SetDefaultValue(0.5);
    GaussMeanCMD->SetParameterName("GaussMean",true);
    GaussMeanCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //GaussMeanCMD->SetUnitCategory()
    GaussMeanCMD->SetDefaultUnit("MeV");
    GaussMeanCMD->SetUnitCandidates("eV keV MeV GeV");

    GaussDistSDevCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setGaussSDev",this);
    GaussDistSDevCMD->SetGuidance("");
    GaussDistSDevCMD->SetDefaultValue(0.01);
    GaussDistSDevCMD->SetParameterName("GaussDistSDev",true);
    GaussDistSDevCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //GaussDistSDevCMD->SetUnitCategory()
    GaussDistSDevCMD->SetDefaultUnit("MeV");
    GaussDistSDevCMD->SetUnitCandidates("eV keV MeV GeV");

    UniformDistEminCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setUniformEmin",this);
    UniformDistEminCMD->SetGuidance("");
    UniformDistEminCMD->SetDefaultValue(0.7);
    UniformDistEminCMD->SetParameterName("UniformDistEmin",true);
    UniformDistEminCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //RadifffusCMD->SetUnitCategory()
    UniformDistEminCMD->SetDefaultUnit("MeV");
    UniformDistEminCMD->SetUnitCandidates("eV keV MeV GeV");

    UniformDistEmaxCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setUniformEmax",this);
    UniformDistEmaxCMD->SetGuidance("");
    UniformDistEmaxCMD->SetDefaultValue(1.);
    UniformDistEmaxCMD->SetParameterName("UniformDistEmax",true);
    UniformDistEmaxCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //UniformDistEmaxCMD->SetUnitCategory()
    UniformDistEmaxCMD->SetDefaultUnit("MeV");
    UniformDistEmaxCMD->SetUnitCandidates("eV keV MeV GeV");

    RayleighDistEmaxCMD = new G4UIcmdWithADoubleAndUnit("/SourceData/setRayleighEmax",this);
    RayleighDistEmaxCMD->SetGuidance("");
    RayleighDistEmaxCMD->SetDefaultValue(1.);
    RayleighDistEmaxCMD->SetParameterName("RayleighDistEmax",true);
    RayleighDistEmaxCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //RayleighDistEmaxCMD->SetUnitCategory()
    RayleighDistEmaxCMD->SetDefaultUnit("MeV");
    RayleighDistEmaxCMD->SetUnitCandidates("eV keV MeV GeV");

    SourcePositionCMD = new G4UIcmdWith3VectorAndUnit("/SourceData/setSourcePosition",this);
    SourcePositionCMD->SetGuidance("Set the x, y and z positions of point source");
    SourcePositionCMD->SetParameterName("x","y","z",true,true);
    SourcePositionCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    SourcePositionCMD->SetDefaultValue(G4ThreeVector(20.,20.,20.));
    SourcePositionCMD->SetDefaultUnit("cm"); // then we dont need to add *cm in command
    SourcePositionCMD->SetUnitCandidates("cm");

    GeneratePositionsCMD = new G4UIcmdWithoutParameter("/SourceData/generatePositions",this);
    GeneratePositionsCMD->SetGuidance("");
    GeneratePositionsCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    GenerateEnergiesCMD = new G4UIcmdWithoutParameter("/SourceData/generateEnergies",this);
    GenerateEnergiesCMD->SetGuidance("");
    GenerateEnergiesCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    GenerateMomDirsCMD = new G4UIcmdWithoutParameter("/SourceData/generateMomDirs",this);
    GenerateMomDirsCMD->SetGuidance("");
    GenerateMomDirsCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    SourceRegionNameCMD = new G4UIcmdWithAString("/SourceData/setSourceName",this);
    SourceRegionNameCMD->SetGuidance("Set one organ name where we will generate the point sources : Head Skull Thyroid UpperSpine Brain LeftLeg RightLeg LeftLegBone RightLegBone Trunk UpperLargeIntestine LowerLargeIntestine	RightArmBone LeftArmBone LeftKidney RightKidney	LeftLung RightLung  LeftOvary RightOvary LeftTeste RightTeste LeftBreast RightBreast  LeftAdrenal	RightAdrenal executeCmd LeftClavicle RightClavicle	LeftScapula	RightScapula Heart  UrinaryBladder MiddleLowerSpine Pelvis	Stomach Pancreas SmallIntestine Uterus Liver Pancreas Kidney Lung Adrenal Spleen RibCage Thymus MaleGenitalia");
    SourceRegionNameCMD->SetParameterName("SourceRegionNameCMD",true);
    SourceRegionNameCMD->SetDefaultValue("Heart");
    //SourceRegionNameCMD->SetCandidates("Head Skull Thyroid UpperSpine Brain LeftLeg RightLeg LeftLegBone RightLegBone Trunk UpperLargeIntestine LowerLargeIntestine	RightArmBone LeftArmBone LeftKidney RightKidney	LeftLung RightLung LeftOvary RightOvary LeftTeste RightTeste LeftBreast RightBreast	LeftAdrenal	RightAdrenal LeftClavicle RightClavicle	LeftScapula	RightScapula Heart UrinaryBladder MiddleLowerSpine Pelvis	Stomach Pancreas SmallIntestine Uterus RibCage Thymus MaleGenitalia Legs LegBone ArmBone Breast Ovary Liver Pancreas Kidney Lung Adrenal Spleen");
    SourceRegionNameCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    Vector_organ_sourceToGenCMD = new G4UIcmdWith3VectorAndUnit("/SourceData/setOrganThreeVector",this);
    Vector_organ_sourceToGenCMD->SetGuidance("Set x, y and z vector of organ where we will generate source points, if not set the code will use the ThreeVector defined in it's database");
    Vector_organ_sourceToGenCMD->SetParameterName("x","y","z",true,true);
    Vector_organ_sourceToGenCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    Vector_organ_sourceToGenCMD->SetDefaultValue(G4ThreeVector(0.,0.,0.));
    Vector_organ_sourceToGenCMD->SetDefaultUnit("cm"); // then we dont need to add *cm in command
    Vector_organ_sourceToGenCMD->SetUnitCandidates("cm");

    SourceRotVector1CMD = new G4UIcmdWith3Vector("/SourceData/setSourceRotVector1",this);
    SourceRotVector1CMD->SetGuidance("");
    SourceRotVector1CMD->SetParameterName("x","y","z",true,true);
    SourceRotVector1CMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    SourceRotVector1CMD->SetDefaultValue(G4ThreeVector(0.,0.,0.));

    SourceRotVector2CMD = new G4UIcmdWith3Vector("/SourceData/setSourceRotVector2",this);
    SourceRotVector2CMD->SetGuidance("");
    SourceRotVector2CMD->SetParameterName("x","y","z",true,true);
    SourceRotVector2CMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    SourceRotVector2CMD->SetDefaultValue(G4ThreeVector(0.,0.,0.));

    box_widthCMD = new G4UIcmdWith3VectorAndUnit("/SourceData/setBoxXYZWidth",this);
    box_widthCMD->SetGuidance("Set the dx, dy and dz of box to generate points, if not set the code will use the dx, dy and dz defined in it's database");
    box_widthCMD->SetParameterName("x","y","z",true,true);
    box_widthCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    box_widthCMD->SetDefaultValue(G4ThreeVector(20.,20.,20.));
    box_widthCMD->SetDefaultUnit("cm"); // then we dont need to add *cm in command
    box_widthCMD->SetUnitCandidates("cm");

    numberOfPointsToGenerateCMD = new G4UIcmdWithAnInteger("/SourceData/setEventsNumForDataGen",this);
    numberOfPointsToGenerateCMD->SetGuidance("Set number of point to save in binary file([organName].bin)");
    numberOfPointsToGenerateCMD->SetParameterName("NumberOfPointsToGenerateCMD",true);
    numberOfPointsToGenerateCMD->SetDefaultValue(1000);
    numberOfPointsToGenerateCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    UsePETCumulativeActCMD = new G4UIcmdWithoutParameter("/SourceData/usePETCumulativeAct",this);
    UsePETCumulativeActCMD->SetGuidance("");
    //UsePETCumulativeActCMD->SetParameterName("UsePETCumulativeActCMD",true);
    //UsePETCumulativeActCMD->SetDefaultValue("no");
    UsePETCumulativeActCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

}
