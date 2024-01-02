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

#include "G4TPhysicsMessenger.hh"
//#include "G4TUserPhysicsList.hh"
#include "G4TVolumeConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

#include "G4UIcmdWithoutParameter.hh"

#include "globals.hh"
#include "G4Tokenizer.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

G4TPhysicsMessenger::G4TPhysicsMessenger(G4TVolumeConstruction* G4TPhysicsObj):myUserPhysics(G4TPhysicsObj){
    
    
    physicsDataDir = new G4UIdirectory("/PhysicsData/");
    physicsDataDir->SetGuidance("Set physics and cuts specification values");
    
    CommandsForPhysics();
}



G4TPhysicsMessenger::~G4TPhysicsMessenger()
{
    delete physicsDataDir;
    delete particle_PhysicsCMD ;
    delete particle_Cuts_energyCMD ;
    delete particle_Cuts_DistanceCMD ;
    delete SetEnergiesForCrossSectionCMD;

    delete PhotoElectricEffectModelCMD;
    delete ComptonScatteringModelCMD;
    delete GammaConversionModelCMD;
    delete RayleighScatteringModelCMD;
    delete ElectronIonisationModelCMD;
    delete ElectronBremModelCMD;
    delete HadronIonisationModelCMD;
    delete IonIonisationModelCMD;
    delete PhysicsDataCMD;
    delete CutsDataCMD;
    delete EnergyRangeDataCMD;

    
    
    /*
delete GammaConversionToMuonModelCMD;
    delete PolarizedComptonModelCMD;
    delete PolarizedGammaConversionModelCMD;
    delete PolarizedPhotoElectricEffectModelCMD;
*/
}

G4double G4TPhysicsMessenger::UseG4Units(G4String Unit){
    
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
// called automatically for each command in the begining of running to set the values entered by the user in Macro file of by typing commands from
// it call G4TVolumeConstruction::SetPhantomModel(G4String newModel), G4TVolumeConstruction::SetPhantomSex(G4String newSex)
// it call G4TPhysicsMessenger::AddBodyPart(G4String newBodyPartSensitivity)
void G4TPhysicsMessenger::SetNewValue(G4UIcommand* command,G4String newValue){
    
    //physics commands
    
    if( command == PhysicsDataCMD )
    {
        G4Tokenizer next(newValue);
        
        G4String PN = next();
        myUserPhysics->setPhysics(PN);
        
        if(PN == "Construct"){
            
            myUserPhysics->setPhotoElectricEffectModel(next());
            myUserPhysics->setComptonScatteringModel(next());
            myUserPhysics->setGammaConversionModel(next());
            myUserPhysics->setRayleighScatteringModel(next());
            myUserPhysics->setElectronIonisationModel(next());
            myUserPhysics->setElectronBremModel(next());
            myUserPhysics->setHadronIonisationModel(next());
            myUserPhysics->setIonIonisationModel(next());
        }
    }
    
    if( command == CutsDataCMD )
    {
        /*
        G4Tokenizer next(newValue);
        G4String car1 = next(), car2 = next() ,Un1 = next() , Un3 = next() ;

        if(car1 == "null"){myUserPhysics->setCutsDistance(0);}
        else {myUserPhysics->setCutsDistance(StoD(car1)*UseG4Units(Un1));}
        
        if(car2 == "null"){myUserPhysics->setCutsEnergy(0);}
        else {myUserPhysics->setCutsEnergy(StoD(car2)*UseG4Units(Un3));}
        */

        G4Tokenizer next(newValue);

        G4String Data = "", dd = next();
        while (!dd.empty()) {

            Data += dd;Data += " ";
            dd = next();
        }

        myUserPhysics->setCutInRangeData(Data);

    }
    if( command == EnergyRangeDataCMD )
    {

        G4Tokenizer next(newValue);

        G4String Data = "", dd = next();
        while (!dd.empty()) {

            Data += dd;Data += " ";
            dd = next();
        }

        myUserPhysics->setEnergyThresholdsData(Data);

    }

    if( command == SetEnergiesForCrossSectionCMD)
    {
        G4Tokenizer next(newValue);
        myUserPhysics->setGenerateCrossSectionTableFlag(true); // that theris this commend then generate tables from runAction

        G4String PN = next();
        myUserPhysics->setParticleForCrossSection(PN);
        G4String Un1 = next();
        G4String P = next();
        //std::cout<< " Dddddd " << P <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@ " << PN << std::endl;
        while (!P.empty()) {
            G4double Energy = StoD(P)*UseG4Units(Un1);
            myUserPhysics->setEnergiesForCrossSectionValues(Energy);

            P = next();
        }
    }

    if( command == particle_Cuts_energyCMD )
    {
        myUserPhysics->setCutsEnergy(particle_Cuts_energyCMD->GetNewDoubleValue(newValue));
    }
    if( command == particle_Cuts_DistanceCMD )
    {
        myUserPhysics->setCutsDistance(particle_Cuts_DistanceCMD->GetNewDoubleValue(newValue));
    }
    if( command == particle_PhysicsCMD )
    {
        myUserPhysics->setPhysics(newValue);
    }
    if( command == PhotoElectricEffectModelCMD )
    {
        myUserPhysics->setPhotoElectricEffectModel(newValue);
    }
    if( command == ComptonScatteringModelCMD )
    {
        myUserPhysics->setComptonScatteringModel(newValue);
    }
    if( command == GammaConversionModelCMD )
    {
        myUserPhysics->setGammaConversionModel(newValue);
    }
    if( command == RayleighScatteringModelCMD )
    {
        myUserPhysics->setRayleighScatteringModel(newValue);
    }
    if( command == ElectronIonisationModelCMD )
    {
        myUserPhysics->setElectronIonisationModel(newValue);
    }
    if( command == ElectronBremModelCMD )
    {
        myUserPhysics->setElectronBremModel(newValue);
    }
    if( command == HadronIonisationModelCMD )
    {
        myUserPhysics->setHadronIonisationModel(newValue);
    }
    if( command == IonIonisationModelCMD )
    {
        myUserPhysics->setIonIonisationModel(newValue);
    }
    
    /*
    if( command == GammaConversionToMuonModelCMD )
    {
        myUserPhysics->setGammaConversionToMuonModel(newValue);
    }
    
    if( command == PolarizedPhotoElectricEffectModelCMD )
    {
        myUserPhysics->setPolarizedPhotoElectricEffectModel(newValue);
    }
    if( command == PolarizedGammaConversionModelCMD )
    {
        myUserPhysics->setPolarizedGammaConversionModel(newValue);
    }
    if( command == PolarizedComptonModelCMD )
    {
        myUserPhysics->setPolarizedComptonModel(newValue);
    }
*/
    
}

// called from constructor to define the commands waited to fill from user
void  G4TPhysicsMessenger::CommandsForPhysics(){
    
    G4UIparameter* param;
    
    PhysicsDataCMD = new G4UIcommand("/PhysicsData/setPhysicsData" ,this);
    param = new G4UIparameter("Physics Name",'s', false);     PhysicsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("PhotoElectricEffectModel",'s', true);  PhysicsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("ComptonScatteringModel",'s', true);  PhysicsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("GammaConversionModel",'s', true);  PhysicsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("RayleighScatteringModel",'s', true);  PhysicsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("ElectronIonisationModel",'s', true);  PhysicsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("ElectronBremModel",'s', true);  PhysicsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("HadronIonisationModel",'s', true);  PhysicsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("IonIonisationModel",'s', true);  PhysicsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    
    CutsDataCMD = new G4UIcommand("/PhysicsData/setCutsInRange" ,this);
    param = new G4UIparameter("particle 1",'s', true);     CutsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Cut in range 1",'s', true);  CutsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Unit 1",'s', true);  CutsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("particle 2",'s', true);     CutsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Cut in range 2",'s', true);  CutsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Unit 2",'s', true);  CutsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("particle 3",'s', true);     CutsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Cut in range 3",'s', true);  CutsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Unit 3",'s', true);  CutsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("particle 4",'s', true);     CutsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Cut in range 4",'s', true);  CutsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Unit 4",'s', true);  CutsDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    EnergyRangeDataCMD = new G4UIcommand("/PhysicsData/setEnergyRange" ,this);
    param = new G4UIparameter("Min Energy",'s', true);     EnergyRangeDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Unit",'s', true);  EnergyRangeDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Min Energy",'s', true);  EnergyRangeDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Unit",'s', true);  EnergyRangeDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error

    SetEnergiesForCrossSectionCMD = new G4UIcommand("/PhysicsData/generateCrossSectionFor" ,this);
    param = new G4UIparameter("Particle name",'s', true);      SetEnergiesForCrossSectionCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    param = new G4UIparameter("Energy Unit",'s', true);      SetEnergiesForCrossSectionCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    for(G4int ds = 0 ; ds < 200 ; ds++){
        G4String f = "Energy " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      SetEnergiesForCrossSectionCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }
    
    particle_Cuts_energyCMD = new G4UIcmdWithADoubleAndUnit("/PhysicsData/setCutsEnergy",this);
    particle_Cuts_energyCMD->SetGuidance("Set cuts energy value in MeV for gamma, e- and e+");
    particle_Cuts_energyCMD->SetParameterName("particle_Cuts_energy",true);
    particle_Cuts_energyCMD->SetDefaultValue(0.000001);
    particle_Cuts_energyCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //particle_Cuts_energyCMD->SetUnitCategory()
    particle_Cuts_energyCMD->SetDefaultUnit("MeV");
    particle_Cuts_energyCMD->SetUnitCandidates("eV keV MeV GeV");
    
    particle_Cuts_DistanceCMD = new G4UIcmdWithADoubleAndUnit("/PhysicsData/setCutsDistance",this);
    particle_Cuts_DistanceCMD->SetGuidance("Set cuts distance value in mm for gamma, e- and e+");
    particle_Cuts_DistanceCMD->SetParameterName("particle_Cuts_distance",true);
    particle_Cuts_DistanceCMD->SetDefaultValue(1.);
    particle_Cuts_DistanceCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    //RadifffusCMD->SetUnitCategory()
    particle_Cuts_DistanceCMD->SetDefaultUnit("mm");
    particle_Cuts_DistanceCMD->SetUnitCandidates("cm mm m");
    
    particle_PhysicsCMD = new G4UIcmdWithAString("/PhysicsData/setPhysicsName",this);
    particle_PhysicsCMD->SetGuidance("Set physics that will be used in interaction of particles with matter either electromagtic physics constructors by setting the Name of constructor, or build your physics by setting MyPhysics as parameters and set the appropriate models with commands  ");
    particle_PhysicsCMD->SetParameterName("particle_Physics",true);
    particle_PhysicsCMD->SetDefaultValue("EMS");
    particle_PhysicsCMD->SetCandidates("EMS EMS1 EMS2 EMS3 EMS4 Livermore Penelope Construct");
    particle_PhysicsCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    PhotoElectricEffectModelCMD = new G4UIcmdWithAString("/PhysicsData/setPhotoElectricEffectModel",this);
    PhotoElectricEffectModelCMD->SetGuidance("");
    PhotoElectricEffectModelCMD->SetParameterName("PhotoElectricEffectModelCMD",true);
    PhotoElectricEffectModelCMD->SetDefaultValue("");
    PhotoElectricEffectModelCMD->SetCandidates("");
    PhotoElectricEffectModelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    ComptonScatteringModelCMD = new G4UIcmdWithAString("/PhysicsData/setComptonScatteringModel",this);
    ComptonScatteringModelCMD->SetGuidance("");
    ComptonScatteringModelCMD->SetParameterName("ComptonScatteringModelCMD",true);
    ComptonScatteringModelCMD->SetDefaultValue("");
    ComptonScatteringModelCMD->SetCandidates("");
    ComptonScatteringModelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    GammaConversionModelCMD = new G4UIcmdWithAString("/PhysicsData/setGammaConversionModel",this);
    GammaConversionModelCMD->SetGuidance("");
    GammaConversionModelCMD->SetParameterName("GammaConversionModelCMD",true);
    GammaConversionModelCMD->SetDefaultValue("");
    GammaConversionModelCMD->SetCandidates("");
    GammaConversionModelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    RayleighScatteringModelCMD = new G4UIcmdWithAString("/PhysicsData/setRayleighScatteringModel",this);
    RayleighScatteringModelCMD->SetGuidance("");
    RayleighScatteringModelCMD->SetParameterName("RayleighScatteringModelCMD",true);
    RayleighScatteringModelCMD->SetDefaultValue("");
    RayleighScatteringModelCMD->SetCandidates("");
    RayleighScatteringModelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    ElectronIonisationModelCMD = new G4UIcmdWithAString("/PhysicsData/setElectronIonisationModel",this);
    ElectronIonisationModelCMD->SetGuidance("");
    ElectronIonisationModelCMD->SetParameterName("ElectronIonisationModelCMD",true);
    ElectronIonisationModelCMD->SetDefaultValue("");
    ElectronIonisationModelCMD->SetCandidates("");
    ElectronIonisationModelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    ElectronBremModelCMD = new G4UIcmdWithAString("/PhysicsData/setElectronBremModel",this);
    ElectronBremModelCMD->SetGuidance("");
    ElectronBremModelCMD->SetParameterName("ElectronBremModelCMD",true);
    ElectronBremModelCMD->SetDefaultValue("");
    ElectronBremModelCMD->SetCandidates("");
    ElectronBremModelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    HadronIonisationModelCMD = new G4UIcmdWithAString("/PhysicsData/setHadronIonisationModel",this);
    HadronIonisationModelCMD->SetGuidance("");
    HadronIonisationModelCMD->SetParameterName("HadronIonisationModelCMD",true);
    HadronIonisationModelCMD->SetDefaultValue("");
    HadronIonisationModelCMD->SetCandidates("");
    HadronIonisationModelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    IonIonisationModelCMD = new G4UIcmdWithAString("/PhysicsData/setIonIonisationModel",this);
    IonIonisationModelCMD->SetGuidance("");
    IonIonisationModelCMD->SetParameterName("IonIonisationModelCMD",true);
    IonIonisationModelCMD->SetDefaultValue("");
    IonIonisationModelCMD->SetCandidates("");
    IonIonisationModelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    
    /*
    GammaConversionToMuonModelCMD = new G4UIcmdWithAString("/PhysicsData/setGammaConversionToMuonModel",this);
    GammaConversionToMuonModelCMD->SetGuidance("Set physics that will be used in interaction of particles with matter either electromagtic physics constructors by setting the Name of constructor, or build your physics by setting MyPhysics as parameters and set the appropriate models with commands  ");
    GammaConversionToMuonModelCMD->SetParameterName("GammaConversionToMuonModelCMD",true);
    GammaConversionToMuonModelCMD->SetDefaultValue("");
    GammaConversionToMuonModelCMD->SetCandidates("");
    GammaConversionToMuonModelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    PolarizedPhotoElectricEffectModelCMD = new G4UIcmdWithAString("/PhysicsData/setPolarizedPhotoElectricEffectModel",this);
    PolarizedPhotoElectricEffectModelCMD->SetGuidance("Set physics that will be used in interaction of particles with matter either electromagtic physics constructors by setting the Name of constructor, or build your physics by setting MyPhysics as parameters and set the appropriate models with commands  ");
    PolarizedPhotoElectricEffectModelCMD->SetParameterName("PolarizedPhotoElectricEffectModelCMD",true);
    PolarizedPhotoElectricEffectModelCMD->SetDefaultValue("");
    PolarizedPhotoElectricEffectModelCMD->SetCandidates("");
    PolarizedPhotoElectricEffectModelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    PolarizedComptonModelCMD = new G4UIcmdWithAString("/PhysicsData/setPolarizedComptonModel",this);
    PolarizedComptonModelCMD->SetGuidance("Set physics that will be used in interaction of particles with matter either electromagtic physics constructors by setting the Name of constructor, or build your physics by setting MyPhysics as parameters and set the appropriate models with commands  ");
    PolarizedComptonModelCMD->SetParameterName("PolarizedComptonModelCMD",true);
    PolarizedComptonModelCMD->SetDefaultValue("");
    PolarizedComptonModelCMD->SetCandidates("");
    PolarizedComptonModelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    PolarizedGammaConversionModelCMD = new G4UIcmdWithAString("/PhysicsData/setPolarizedGammaConversionModel",this);
    PolarizedGammaConversionModelCMD->SetGuidance("Set physics that will be used in interaction of particles with matter either electromagtic physics constructors by setting the Name of constructor, or build your physics by setting MyPhysics as parameters and set the appropriate models with commands  ");
    PolarizedGammaConversionModelCMD->SetParameterName("PolarizedGammaConversionModelCMD",true);
    PolarizedGammaConversionModelCMD->SetDefaultValue("");
    PolarizedGammaConversionModelCMD->SetCandidates("");
    PolarizedGammaConversionModelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
*/
    
    
}

