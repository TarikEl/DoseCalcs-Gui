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

#ifndef G4TUserPhysicsList_h
#define G4TUserPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
//#include "G4TOutputText.hh"
#include "G4TPhysicsMessenger.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4TComptonProcess.hh"

// Physics
extern G4String CutInRangeData;
extern G4String EnergyThresholdsData;
extern G4bool IsEcutsSet;
extern G4bool IsDcutsSet;
extern G4double CutsEnergy;
extern G4double CutsDistance;
extern G4String ParticlePysics;
extern G4String PhotoElectricEffectModel;
extern G4String PolarizedPhotoElectricEffectModel;
extern G4String ComptonScatteringModel;
extern G4String PolarizedComptonModel;
extern G4String GammaConversionModel;
extern G4String PolarizedGammaConversionModel;
extern G4String RayleighScatteringModel;
extern G4String GammaConversionToMuonModel;
extern G4String ElectronIonisationModel;
extern G4String ElectronBremModel;
extern G4String HadronIonisationModel;
extern G4String HadronBremModel;
extern G4String IonIonisationModel;
extern G4bool GenerateCrossSectionTableFlag;
extern G4String ParticleForCrossSection;
extern std::vector<G4double> EnergiesForCrossSectionValues;


extern G4String TestPointsPositions;


class G4TUserPhysicsList: public G4VUserPhysicsList
{
public:
    G4TUserPhysicsList();
    ~G4TUserPhysicsList();

private:

    void ShowSourceParameters();


    //G4TComptonProcess* ComptonWrappedProcess;
public:

/*
    void setPhysics(G4String);
    void setCutsEnergy(G4double);
    void setCutsDistance(G4double);

    void setPhotoElectricEffectModel(G4String);
    void setPolarizedPhotoElectricEffectModel(G4String);
    void setComptonScatteringModel(G4String);
    void setPolarizedComptonModel(G4String);
    void setGammaConversionModel(G4String);
    void setPolarizedGammaConversionModel(G4String);
    void setRayleighScatteringModel(G4String);
    void setGammaConversionToMuonModel(G4String);

    void setElectronIonisationModel(G4String);
    void setElectronBremModel(G4String);
    void setHadronIonisationModel(G4String);
    void setIonIonisationModel(G4String);
    G4double getCutsEnergy() const { return CutsEnergy;}
    G4double getCutsDistance() const { return CutsDistance;}

*/

    //G4String getPhysics()const {return ParticlePysics;};
    //G4double getCutsEnergy()const {return CutsEnergy;};
    //G4double getCutsDistance()const {return CutsDistance;};

protected:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

private:

    G4TPhysicsMessenger* messengerPhyObj;

    G4VPhysicsConstructor* G4VEmPhysicsConstructorObj;
    G4VPhysicsConstructor* G4VHadromPhysicsConstructorObj;
    G4VPhysicsConstructor* decPhysicsList;
    G4VPhysicsConstructor* RadioactivedecPhysicsList;

};
#endif
