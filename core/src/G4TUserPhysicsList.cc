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


//MultipleScattering process
#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4MscStepLimitType.hh"
#include "G4UrbanMscModel.hh"
#include "G4DummyModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4CoulombScattering.hh"

//------------------ for electron

//process
#include "G4eIonisation.hh"
//Models
#include "G4MollerBhabhaModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4PAIModel.hh"
#include "G4PAIPhotModel.hh"

//process
#include "G4eBremsstrahlung.hh"
//Models
#include "G4SeltzerBergerModel.hh"
#include "G4eBremsstrahlungRelModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4PenelopeBremsstrahlungModel.hh"
// Generator
#include "G4Generator2BS.hh"
#include "G4Generator2BN.hh"

//process
//#include "G4ePolarizedBremsstrahlung.hh"
//Models
//#include "G4ePolarizedBremsstrahlungModel.hh"

//process
#include "G4ePairProduction.hh"

//process
#include "G4eplusAnnihilation.hh"

#include "G4UAtomicDeexcitation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"

//----------------------------- for hadrons

// process
#include "G4hIonisation.hh"
// models
#include "G4BetheBlochModel.hh"
#include "G4BraggModel.hh"
#include "G4ICRU73QOModel.hh"
//#include "G4PAIModel.hh" // is defined for electron
//#include "G4PAIPhotModel.hh" // is defined for electron

// process
#include "G4hBremsstrahlung.hh"
// models
#include "G4hBremsstrahlungModel.hh"

// process
#include "G4hPairProduction.hh"
// models
#include "G4hPairProductionModel.hh"


//------------------ for ions

// process
#include "G4ionIonisation.hh"
// models
#include "G4BetheBlochModel.hh"
#include "G4BetheBlochIonGasModel.hh"
#include "G4BraggIonModel.hh"
#include "G4BraggIonGasModel.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4AtimaEnergyLossModel.hh"
#include "G4LindhardSorensenIonModel.hh"


#include "G4NuclearStopping.hh"


// neutron processes and models
#include "G4HadronElasticProcess.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPThermalScatteringData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPThermalScattering.hh"



/*#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"
#include "G4Neutron.hh"
*/


// ---------------------- Gamma electromagnetic process and models

// process
#include "G4PhotoElectricEffect.hh"
// model
#include "G4PEEffectFluoModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
//#include "G4LivermorePolarizedPhotoElectricModel.hh"
#include "G4PenelopePhotoElectricModel.hh"
//#include "G4LivermorePolarizedPhotoElectricGDModel.hh"


// process
#include "G4ComptonScattering.hh"
// model
#include "G4TKleinNishinaCompton.hh"
#include "G4KleinNishinaCompton.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LivermoreComptonModel.hh"
//#include "G4LivermoreComptonModifiedModel.hh"
#include "G4LivermorePolarizedComptonModel.hh"
#include "G4LowEPComptonModel.hh"
#include "G4PenelopeComptonModel.hh"

// process
#include "G4PolarizedCompton.hh"
// model
#include "G4PolarizedComptonModel.hh"



// process
#include "G4GammaConversion.hh"
// model
#include "G4BetheHeitlerModel.hh"
#include "G4BetheHeitler5DModel.hh"
#include "G4PairProductionRelModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4BoldyshevTripletModel.hh"
#include "G4LivermoreNuclearGammaConversionModel.hh"
#include "G4LivermorePolarizedGammaConversionModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
//#include "G4LivermoreGammaConversionModelRC.hh"

// process
#include "G4PolarizedGammaConversion.hh"
// model
#include "G4PolarizedGammaConversionModel.hh"


// process
#include "G4RayleighScattering.hh"
// model
#include "G4LivermoreRayleighModel.hh"
#include "G4LivermorePolarizedRayleighModel.hh"
#include "G4PenelopeRayleighModel.hh"


// process
//#include "G4GammaConversionToMuons.hh"


// ---------------------------- decay process
#include "G4DecayPhysics.hh"


// --------------------------- all types of particles
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"



// ------------------------- elecromagnetic constructors

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"


#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmModelActivator.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"


#include "G4TVolumeConstruction.hh"
#include "G4TUserPhysicsList.hh"
#include "G4RunManager.hh"
#include "G4BiasingHelper.hh"
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

// ------------------------- hadron physics constructors

#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT_ATL.hh"
#include "G4HadronPhysicsFTFP_BERT_TRV.hh"
#include "G4HadronPhysicsQGSP_FTFP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4HadronPhysicsShielding.hh"
#include "G4HadronPhysicsShieldingLEND.hh"


#include "G4RadioactiveDecayPhysics.hh"

G4TUserPhysicsList::G4TUserPhysicsList():  G4VUserPhysicsList()
{

    G4VEmPhysicsConstructorObj = nullptr;
    G4VHadromPhysicsConstructorObj = nullptr;
    decPhysicsList = new G4DecayPhysics(0);
    RadioactivedecPhysicsList = new G4RadioactiveDecayPhysics(0);

    const G4TVolumeConstruction* TConstruction2 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    int verbose = 0 ;

    if (ParticlePysics == "EMS") {
        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics(verbose);
        if (TConstruction2->getParticleName() == "neutron") {
            G4VHadromPhysicsConstructorObj = new G4HadronPhysicsFTFP_BERT();
        }
    }
    else if(ParticlePysics == "EMS1"){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n EMS1 Physics " << G4endl;

        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option1(verbose);
        if (TConstruction2->getParticleName() == "neutron") {
            G4VHadromPhysicsConstructorObj = new G4HadronPhysicsFTFP_BERT();
        }
    }
    else if(ParticlePysics == "EMS2"){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n EMS2 Physics " << G4endl;

        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option2(verbose);
        if (TConstruction2->getParticleName() == "neutron") {
            G4VHadromPhysicsConstructorObj = new G4HadronPhysicsFTFP_BERT();
        }
    }
    else if(ParticlePysics == "EMS3"){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n EMS3 Physics " << G4endl;

        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
        if (TConstruction2->getParticleName() == "neutron") {
            G4VHadromPhysicsConstructorObj = new G4HadronPhysicsFTFP_BERT();
        }
    }
    else if(ParticlePysics == "EMS4"){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n EMS4 Physics " << G4endl;

        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option4(verbose);
        if (TConstruction2->getParticleName() == "neutron") {
            G4VHadromPhysicsConstructorObj = new G4HadronPhysicsFTFP_BERT();
        }
    }
    else if(ParticlePysics == "Livermore"){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n Livermore Physics " << G4endl;

        G4VEmPhysicsConstructorObj = new G4EmLivermorePhysics(verbose);
        if (TConstruction2->getParticleName() == "neutron") {
            G4VHadromPhysicsConstructorObj = new G4HadronPhysicsFTFP_BERT();
        }
    }
    else if(ParticlePysics == "Penelope"){
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n penelope Physics " << G4endl;

        G4VEmPhysicsConstructorObj = new G4EmPenelopePhysics(verbose);
        if (TConstruction2->getParticleName() == "neutron") {
            G4VHadromPhysicsConstructorObj = new G4HadronPhysicsFTFP_BERT();
        }
    }
    else if(ParticlePysics == "HADRON_FTFP_BERT"){
        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
        G4VHadromPhysicsConstructorObj = new G4HadronPhysicsFTFP_BERT();
    }
    else if(ParticlePysics == "HADRON_FTFP_BERT_ATL"){
        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
        G4VHadromPhysicsConstructorObj = new G4HadronPhysicsFTFP_BERT_ATL();
    }
    else if(ParticlePysics == "HADRON_FTFP_BERT_TRV"){
        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
        G4VHadromPhysicsConstructorObj = new G4HadronPhysicsFTFP_BERT_TRV();
    }
    else if(ParticlePysics == "HADRON_QGSP_FTFP_BERT"){
        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
        G4VHadromPhysicsConstructorObj = new G4HadronPhysicsQGSP_FTFP_BERT();
    }
    else if(ParticlePysics == "HADRON_QGSP_BERT"){
        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
        G4VHadromPhysicsConstructorObj = new G4HadronPhysicsQGSP_BERT();
    }
    else if(ParticlePysics == "HADRON_QGSP_BERT_HP"){
        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
        G4VHadromPhysicsConstructorObj = new G4HadronPhysicsQGSP_BERT_HP();
    }
    else if(ParticlePysics == "HADRON_QGSP_BIC"){
        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
        G4VHadromPhysicsConstructorObj = new G4HadronPhysicsQGSP_BIC();
    }
    else if(ParticlePysics == "HADRON_QGSP_BIC_AllHP"){
        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
        G4VHadromPhysicsConstructorObj = new G4HadronPhysicsQGSP_BIC_AllHP();
    }
    else if(ParticlePysics == "HADRON_INCLXX"){
        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
        G4VHadromPhysicsConstructorObj = new G4HadronPhysicsINCLXX();
    }
    else if(ParticlePysics == "HADRON_Shielding"){
        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
        G4VHadromPhysicsConstructorObj = new G4HadronPhysicsShielding();
    }
    else if(ParticlePysics == "HADRON_ShieldingLEND"){
        G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
        G4VHadromPhysicsConstructorObj = new G4HadronPhysicsShieldingLEND();
    }
    else {
        G4Exception("MyPhysicsList::MyPhysicsList", "InvalidSetup", FatalException, "Invalid physics choice");
    }

}

G4TUserPhysicsList::~G4TUserPhysicsList()
{
    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;

    delete decPhysicsList;
    //delete messengerPhyObj;

    //if(ParticlePysics == "EMS" || ParticlePysics == "EMS1" || ParticlePysics == "EMS2" || ParticlePysics == "EMS3"|| ParticlePysics == "EMS4"|| ParticlePysics == "Livermore"|| ParticlePysics == "Penelope"){
    if(G4VEmPhysicsConstructorObj != nullptr){
        delete G4VEmPhysicsConstructorObj;
    }

    if(G4VHadromPhysicsConstructorObj != nullptr){
        delete G4VHadromPhysicsConstructorObj;
    }

}

void G4TUserPhysicsList::ConstructParticle()
{

    //G4cout << "from function : " << __FUNCTION__<< G4endl;

    G4BosonConstructor* Boson = new G4BosonConstructor();
    G4LeptonConstructor* Lepton = new G4LeptonConstructor();
    G4MesonConstructor* Meson = new G4MesonConstructor();
    G4BaryonConstructor* Baryon = new G4BaryonConstructor();
    G4IonConstructor* Ion = new G4IonConstructor();
    G4ShortLivedConstructor* ShortLived = new G4ShortLivedConstructor();

    Boson->ConstructParticle();
    Lepton->ConstructParticle();
    Meson->ConstructParticle();
    Baryon->ConstructParticle();
    Ion->ConstructParticle();
    ShortLived->ConstructParticle();

}

void G4TUserPhysicsList::ConstructProcess()
{

    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;

    if(ParticlePysics == "EMS" || ParticlePysics == "EMS1" || ParticlePysics == "EMS2" || ParticlePysics == "EMS3"|| ParticlePysics == "EMS4"|| ParticlePysics == "Livermore"|| ParticlePysics == "Penelope"){
        G4VEmPhysicsConstructorObj->ConstructProcess();

        const G4TVolumeConstruction* TConstruction2 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        if (TConstruction2->getParticleName() == "neutron") {
            G4VHadromPhysicsConstructorObj->ConstructProcess();
        }
    }
    else if( ParticlePysics == "Construct"){

        const G4TVolumeConstruction* TConstruction2 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

#if VERBOSE_USE
        G4cout << "\n\n - The process are builded by the G4PhysicsListHelper " << G4endl;
#endif

        // muon & hadron bremsstrahlung and pair production
        //G4MuBremsstrahlung* mub = new G4MuBremsstrahlung();
        //G4MuPairProduction* mup = new G4MuPairProduction();

        G4hPairProduction* hadronPairProdProcess = new G4hPairProduction();
        G4hBremsstrahlung* hadronBremProcess = new G4hBremsstrahlung();
        //G4hMultipleScattering* hmsc = new G4hMultipleScattering("ionmsc");

        G4ePairProduction* ee = new G4ePairProduction();

        // nuclear stopping
        G4NuclearStopping* NuclearStoppingProcess = new G4NuclearStopping();
        NuclearStoppingProcess->SetMaxKinEnergy(MeV);

        // Add standard EM Processes
        G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

        auto myParticleIterator=GetParticleIterator();

        /*
        // Step limitation seen as a process
        G4StepLimiter* stepLimiter = new G4StepLimiter();
        G4ProcessManager* pmanager = G4Electron::Electron()->GetProcessManager();
        pmanager->AddDiscreteProcess(stepLimiter); // The first item (max step size) is automatically taken into account by
        // the G4 kernel while the others items must be managed by the user, in this case we use G4StepLimiter
        //pmanager->AddProcess(new G4UserSpecialCuts(),-1,-1,1); // Concerning the others cuts values as ie in Volume contruction
        // if we use fStepLimit = new G4UserLimits(DBL_MAX,DBL_MAX,maxTime); and not just max_step_size, we use G4UserSpecialCuts
        pmanager = G4Positron::Positron()->GetProcessManager();
        pmanager->AddDiscreteProcess(stepLimiter);
        pmanager = G4Gamma::Gamma()->GetProcessManager();
        pmanager->AddDiscreteProcess(stepLimiter);
        pmanager = G4Proton::Proton()->GetProcessManager();
        pmanager->AddDiscreteProcess(stepLimiter);
        */


        myParticleIterator->reset();
        while( (*myParticleIterator)() )
        {

            G4ParticleDefinition* particle = myParticleIterator->value();
            G4String particleName = particle->GetParticleName();

            if (!particle) { continue; }

            if (particleName == "gamma") {

                G4PhotoElectricEffect* G4PhotoElectricEffectobj = new G4PhotoElectricEffect();
                if(PhotoElectricEffectModel == "1"){
                    //G4PEEffectFluoModel* G4PEEffectFluoModelobj = new G4PEEffectFluoModel();
                    G4PhotoElectricEffectobj->SetEmModel(new G4PEEffectFluoModel());
                }
                else if (PhotoElectricEffectModel == "2") {
                    //G4LivermorePhotoElectricModel* G4LivermorePhotoElectricModelobj = new G4LivermorePhotoElectricModel();
                    G4PhotoElectricEffectobj->SetEmModel(new G4LivermorePhotoElectricModel());
                }
                else if (PhotoElectricEffectModel == "3") {
                    //G4LivermorePolarizedPhotoElectricModel* G4LivermorePolarizedPhotoElectricModelobj = new G4LivermorePolarizedPhotoElectricModel();
                    //G4PhotoElectricEffectobj->SetEmModel(new G4LivermorePolarizedPhotoElectricModel());
                    G4PhotoElectricEffectobj->SetEmModel(new G4PenelopePhotoElectricModel());
                }
                else if (PhotoElectricEffectModel == "4") {
                    //G4PenelopePhotoElectricModel* G4PenelopePhotoElectricModelobj = new G4PenelopePhotoElectricModel();
                    G4PhotoElectricEffectobj->SetEmModel(new G4PenelopePhotoElectricModel());                }
                else{
                    //G4PEEffectFluoModel* G4PEEffectFluoModelobj = new G4PEEffectFluoModel();
                    G4PhotoElectricEffectobj->SetEmModel(new G4PEEffectFluoModel());
                }
                ph->RegisterProcess(G4PhotoElectricEffectobj,  particle);


                G4ComptonScattering* G4ComptonScatteringobj = new G4ComptonScattering();
                G4bool UseWeight = true;
                if(UseWeight == true){ //for Variance reduction using weight techniques

                    G4TComptonProcess* ComptonWrappedProcess = new G4TComptonProcess();
                    G4ProcessManager* pmanager = particle->GetProcessManager();
                    G4TKleinNishinaCompton* MyKleinNishinaCompton = new G4TKleinNishinaCompton();
                    MyKleinNishinaCompton->setSplittingNum(5);
                    G4ComptonScatteringobj->SetEmModel(MyKleinNishinaCompton);
                    ComptonWrappedProcess->RegisterProcess(G4ComptonScatteringobj);
                    pmanager->AddProcess(ComptonWrappedProcess,-1,-1, 3);
                }
                else if(UseWeight == false){

                    if(ComptonScatteringModel == "1"){
                        G4ComptonScatteringobj->SetEmModel(new G4KleinNishinaCompton());
                    }
                    else if(ComptonScatteringModel == "2"){
                        G4ComptonScatteringobj->SetEmModel(new G4KleinNishinaModel());
                    }
                    else if(ComptonScatteringModel == "3"){
                        G4ComptonScatteringobj->SetEmModel(new G4LowEPComptonModel());
                    }
                    else if(ComptonScatteringModel == "4"){
                        G4ComptonScatteringobj->SetEmModel(new G4LivermoreComptonModel());
                    }
                    else if(ComptonScatteringModel == "5"){
                        //G4ComptonScatteringobj->SetEmModel(new G4LivermoreComptonModifiedModel());
                        G4ComptonScatteringobj->SetEmModel(new G4LivermoreComptonModel());
                    }
                    else if(ComptonScatteringModel == "6"){
                        G4ComptonScatteringobj->SetEmModel(new G4LivermorePolarizedComptonModel());
                    }
                    else if(ComptonScatteringModel == "7"){
                        G4ComptonScatteringobj->SetEmModel(new G4PenelopeComptonModel());
                    }
                    else if(ComptonScatteringModel == "8"){
#if VERBOSE_USE
                        G4cout << " ComptonScatteringModel == 8" << G4endl;
#endif
                        G4TKleinNishinaCompton* MyKleinNishinaCompton = new G4TKleinNishinaCompton();
                        MyKleinNishinaCompton->setSplittingNum(4);
                        G4ComptonScatteringobj->SetEmModel(MyKleinNishinaCompton);
                    }
                    else {
                        G4ComptonScatteringobj->SetEmModel(new G4KleinNishinaCompton());
                    }
                    ph->RegisterProcess(G4ComptonScatteringobj,  particle);
                }


                //G4cout << "\n\n\n\n\n\n\n\n" << TConstruction2->BiasFlag << G4endl;
                if(TConstruction2->getBiasFlag() == "yes"){

                    G4ProcessManager* pmanager = particle->GetProcessManager();
                    G4BiasingHelper::ActivatePhysicsBiasing(pmanager, "compt");
                    G4BiasingHelper::ActivateNonPhysicsBiasing(pmanager);
                }
                //G4cout << "\n\n\n\n\n\n\n\n" << G4ComptonScatteringobj->GetProcessName() << G4endl;


                G4GammaConversion* G4GammaConversionobj = new G4GammaConversion();
                if(GammaConversionModel == "1"){
                    G4GammaConversionobj->SetEmModel(new G4BetheHeitlerModel());
                }
                else if(GammaConversionModel == "2"){
                    G4GammaConversionobj->SetEmModel(new G4BetheHeitler5DModel());
                }
                else if(GammaConversionModel == "3"){
                    G4GammaConversionobj->SetEmModel(new G4PairProductionRelModel());
                }
                else if(GammaConversionModel == "4"){
                    G4GammaConversionobj->SetEmModel(new G4LivermoreGammaConversionModel());
                }
                else if(GammaConversionModel == "5"){
                    G4GammaConversionobj->SetEmModel(new G4BoldyshevTripletModel());
                }
                else if(GammaConversionModel == "6"){
                    G4GammaConversionobj->SetEmModel(new G4LivermoreNuclearGammaConversionModel());
                }
                else if(GammaConversionModel == "7"){
                    G4GammaConversionobj->SetEmModel(new G4LivermorePolarizedGammaConversionModel());
                }
                else if(GammaConversionModel == "8"){
                    G4GammaConversionobj->SetEmModel(new G4PenelopeGammaConversionModel());
                }
                else{
                    G4GammaConversionobj->SetEmModel(new G4BetheHeitlerModel());
                }
                ph->RegisterProcess(G4GammaConversionobj,  particle);



                G4RayleighScattering* G4RayleighScatteringobj = new G4RayleighScattering();
                if(RayleighScatteringModel =="1"){
                    G4RayleighScatteringobj->SetEmModel(new G4LivermoreRayleighModel());
                }
                else if(RayleighScatteringModel == "2"){
                    G4RayleighScatteringobj->SetEmModel(new G4LivermorePolarizedRayleighModel());
                }
                else if(RayleighScatteringModel =="3"){
                    G4RayleighScatteringobj->SetEmModel(new G4PenelopeRayleighModel());
                }
                else{
                    G4RayleighScatteringobj->SetEmModel(new G4LivermoreRayleighModel());
                }
                ph->RegisterProcess(G4RayleighScatteringobj,  particle);


                //G4GammaConversionToMuons* G4GammaConversionToMuonsObj  = new G4GammaConversionToMuons();
                //ph->RegisterProcess(G4GammaConversionToMuonsObj,  particle);

            }
            else if (particleName == "e-") {

                G4eMultipleScattering* msc = new G4eMultipleScattering();

                G4eIonisation* eIoni = new G4eIonisation();

                if(ElectronIonisationModel == "1"){
                    eIoni->SetEmModel(new G4MollerBhabhaModel());
                }
                else if(ElectronIonisationModel == "2"){
                    eIoni->SetEmModel(new G4LivermoreIonisationModel());
                }
                /*
                else if(ElectronIonisationModel == "3"){
                    eIoni->SetEmModel(new G4PenelopeIonisationModel());
                }
                else if(ElectronIonisationModel == "4"){
                    eIoni->SetEmModel(new G4PAIModel());
                }
                else if(ElectronIonisationModel == "5"){
                    eIoni->SetEmModel(new G4PAIPhotModel());
                }*/
                else {
                    eIoni->SetEmModel(new G4MollerBhabhaModel());
                }

                ph->RegisterProcess(eIoni, particle);

                G4eBremsstrahlung* brem = new G4eBremsstrahlung();
                if(ElectronBremModel == "1"){
                    G4SeltzerBergerModel* bremModel = new G4SeltzerBergerModel();
                    bremModel->SetAngularDistribution(new G4Generator2BS());
                    brem->SetEmModel(bremModel);
                }
                else if(ElectronBremModel == "2"){
                    G4eBremsstrahlungRelModel* bremModel = new G4eBremsstrahlungRelModel();
                    bremModel->SetAngularDistribution(new G4Generator2BS());
                    brem->SetEmModel(bremModel);
                }
                else if(ElectronBremModel == "3"){
                    G4LivermoreBremsstrahlungModel* bremModel = new G4LivermoreBremsstrahlungModel();
                    bremModel->SetAngularDistribution(new G4Generator2BS());
                    brem->SetEmModel(bremModel);
                }
                else if(ElectronBremModel == "4"){
                    G4PenelopeBremsstrahlungModel* bremModel = new G4PenelopeBremsstrahlungModel();
                    bremModel->SetAngularDistribution(new G4Generator2BS());
                    brem->SetEmModel(bremModel);
                }
                else {
                    G4SeltzerBergerModel* bremModel = new G4SeltzerBergerModel();
                    bremModel->SetAngularDistribution(new G4Generator2BS());
                    brem->SetEmModel(bremModel);
                }

                ph->RegisterProcess(brem, particle);

                // register processes
                ph->RegisterProcess(msc, particle);
                ph->RegisterProcess(ee, particle);

            }
            else if (particleName == "e+") {

                G4eMultipleScattering* msc = new G4eMultipleScattering();

                G4eIonisation* eIoni = new G4eIonisation();
                eIoni->SetEmModel(new G4MollerBhabhaModel()); // which is applicable to the positron
                ph->RegisterProcess(eIoni, particle);

                G4eBremsstrahlung* brem = new G4eBremsstrahlung();
                if(ElectronBremModel == "1"){
                    G4SeltzerBergerModel* bremModel = new G4SeltzerBergerModel();
                    bremModel->SetAngularDistribution(new G4Generator2BS());
                    brem->SetEmModel(bremModel);
                }
                else if(ElectronBremModel == "2"){
                    G4eBremsstrahlungRelModel* bremModel = new G4eBremsstrahlungRelModel();
                    bremModel->SetAngularDistribution(new G4Generator2BS());
                    brem->SetEmModel(bremModel);
                }
                else if(ElectronBremModel == "3"){
                    G4LivermoreBremsstrahlungModel* bremModel = new G4LivermoreBremsstrahlungModel();
                    bremModel->SetAngularDistribution(new G4Generator2BS());
                    brem->SetEmModel(bremModel);
                }
                else if(ElectronBremModel == "4"){
                    G4PenelopeBremsstrahlungModel* bremModel = new G4PenelopeBremsstrahlungModel();
                    bremModel->SetAngularDistribution(new G4Generator2BS());
                    brem->SetEmModel(bremModel);
                }
                else {
                    G4SeltzerBergerModel* bremModel = new G4SeltzerBergerModel();
                    bremModel->SetAngularDistribution(new G4Generator2BS());
                    brem->SetEmModel(bremModel);
                }

                ph->RegisterProcess(brem, particle);

                //G4eBremsstrahlung* brem = new G4eBremsstrahlung();
                //G4SeltzerBergerModel* br1 = new G4SeltzerBergerModel();
                //G4eBremsstrahlungRelModel* br2 = new G4eBremsstrahlungRelModel();
                //br1->SetAngularDistribution(new G4Generator2BS());
                //br2->SetAngularDistribution(new G4Generator2BS());
                //brem->SetEmModel(br1);
                //brem->SetEmModel(br2);
                //br2->SetLowEnergyLimit(GeV);

                // register processes
                ph->RegisterProcess(msc, particle);
                ph->RegisterProcess(ee, particle);

                if (particleName == "e+"){

                    ph->RegisterProcess(new G4eplusAnnihilation(), particle);
                }
            }
            else if (particleName == "proton" || particleName == "anti_proton") {

                G4hMultipleScattering* pmsc = new G4hMultipleScattering();

                G4hIonisation* hadronIoniProcess = new G4hIonisation();
                if(HadronIonisationModel == "1"){
                    hadronIoniProcess->SetEmModel(new G4BetheBlochModel());
                }
                else if(HadronIonisationModel == "2"){
                    hadronIoniProcess->SetEmModel(new G4BraggModel());
                }
                else if(HadronIonisationModel == "3"){
                    hadronIoniProcess->SetEmModel(new G4ICRU73QOModel());
                }
                else{
                    hadronIoniProcess->SetEmModel(new G4BetheBlochModel());
                }
                ph->RegisterProcess(hadronIoniProcess, particle);

                ph->RegisterProcess(pmsc, particle);
                ph->RegisterProcess(hadronBremProcess, particle);
                ph->RegisterProcess(hadronPairProdProcess, particle);
                ph->RegisterProcess(NuclearStoppingProcess, particle);

            }
            else if (particleName == "alpha" ||  particleName == "He3") {

                G4ionIonisation* ionIoniProcess = new G4ionIonisation();
                if(IonIonisationModel == "1"){
                    ionIoniProcess->SetEmModel(new G4BetheBlochModel());
                }
                else if(IonIonisationModel == "2"){
                    ionIoniProcess->SetEmModel(new G4BetheBlochIonGasModel());
                }
                else if(IonIonisationModel == "3"){
                    ionIoniProcess->SetEmModel(new G4BraggIonModel());
                }
                else if(IonIonisationModel == "4"){
                    ionIoniProcess->SetEmModel(new G4BraggIonGasModel());
                }
                else if(IonIonisationModel == "5"){
                    ionIoniProcess->SetEmModel(new G4IonParametrisedLossModel());
                }
                else if(IonIonisationModel == "6"){
                    ionIoniProcess->SetEmModel(new G4AtimaEnergyLossModel());
                }
                else if(IonIonisationModel == "7"){
                    ionIoniProcess->SetEmModel(new G4LindhardSorensenIonModel());
                }
                else {
                    ionIoniProcess->SetEmModel(new G4BetheBlochModel());
                }
                ionIoniProcess->SetStepFunction(0.1, 10*um);

                G4hMultipleScattering* msc = new G4hMultipleScattering();

                ph->RegisterProcess(msc, particle);
                ph->RegisterProcess(NuclearStoppingProcess, particle);

            }
            /*
            else if (particleName == "mu+" || particleName == "mu-"    ) {

                G4MuMultipleScattering* mumsc = new G4MuMultipleScattering();
                G4MuIonisation* muIoni = new G4MuIonisation();

                ph->RegisterProcess(mumsc, particle);
                ph->RegisterProcess(muIoni, particle);
                ph->RegisterProcess(mub, particle);
                ph->RegisterProcess(mup, particle);

            } else if (particleName == "pi+" || particleName == "pi-" ) {

                G4hMultipleScattering* pimsc = new G4hMultipleScattering();
                G4hIonisation* hadronIoniProcess = new G4hIonisation();

                ph->RegisterProcess(pimsc, particle);
                ph->RegisterProcess(hadronIoniProcess, particle);
                ph->RegisterProcess(hadronBremProcess, particle);
                ph->RegisterProcess(hadronPairProdProcess, particle);

            } else if (particleName == "kaon+" || particleName == "kaon-" ) {

                G4hMultipleScattering* kmsc = new G4hMultipleScattering();
                G4hIonisation* hadronIoniProcess = new G4hIonisation();

                ph->RegisterProcess(kmsc, particle);
                ph->RegisterProcess(hadronIoniProcess, particle);
                ph->RegisterProcess(hadronBremProcess, particle);
                ph->RegisterProcess(hadronPairProdProcess, particle);

            } else if (particleName == "B+" ||
                       particleName == "B-" ||
                       particleName == "D+" ||
                       particleName == "D-" ||
                       particleName == "Ds+" ||
                       particleName == "Ds-" ||
                       particleName == "anti_He3" ||
                       particleName == "anti_alpha" ||
                       particleName == "anti_deuteron" ||
                       particleName == "anti_lambda_c+" ||
                       particleName == "anti_omega-" ||
                       particleName == "anti_sigma_c+" ||
                       particleName == "anti_sigma_c++" ||
                       particleName == "anti_sigma+" ||
                       particleName == "anti_sigma-" ||
                       particleName == "anti_triton" ||
                       particleName == "anti_xi_c+" ||
                       particleName == "anti_xi-" ||
                       particleName == "deuteron" ||
                       particleName == "lambda_c+" ||
                       particleName == "omega-" ||
                       particleName == "sigma_c+" ||
                       particleName == "sigma_c++" ||
                       particleName == "sigma+" ||
                       particleName == "sigma-" ||
                       particleName == "tau+" ||
                       particleName == "tau-" ||
                       particleName == "triton" ||
                       particleName == "xi_c+" ||
                       particleName == "xi-" ) {

                ph->RegisterProcess(hmsc, particle);
                ph->RegisterProcess(new G4hIonisation(), particle);
                ph->RegisterProcess(NuclearStoppingProcess, particle);
            }

            else if (particleName == "GenericIon") {

                G4ionIonisation* ionIoni = new G4ionIonisation();
                ionIoni->SetEmModel(new G4IonParametrisedLossModel());
                ionIoni->SetStepFunction(0.1, 1*um);

                ph->RegisterProcess(hmsc, particle);
                ph->RegisterProcess(ionIoni, particle);
                ph->RegisterProcess(NuclearStoppingProcess, particle);

            }
*/

            /*else if (particleName == "neutron") {

                G4HadronElasticProcess* theNeutronElasticProcess = new G4HadronElasticProcess;
                // Cross Section Data set
                G4NeutronHPElasticData* theHPElasticData = new G4NeutronHPElasticData;
                theNeutronElasticProcess->AddDataSet(theHPElasticData);
                G4NeutronHPThermalScatteringData* theHPThermalScatteringData = new G4NeutronHPThermalScatteringData;
                theNeutronElasticProcess->AddDataSet(theHPThermalScatteringData);
                // Models
                G4NeutronHPElastic* theNeutronElasticModel = new G4NeutronHPElastic;
                theNeutronElasticModel->SetMinEnergy(4.0*eV);
                theNeutronElasticProcess->RegisterMe(theNeutronElasticModel);
                G4NeutronHPThermalScattering* theNeutronThermalElasticModel = new G4NeutronHPThermalScattering;
                theNeutronThermalElasticModel->SetMaxEnergy(4.0*eV);
                theNeutronElasticProcess->RegisterMe(theNeutronThermalElasticModel);
                // Apply Processes to Process Manager of Neutron
                G4ProcessManager* pmanager = G4Neutron::Neutron()->GetProcessManager();
                pmanager->AddDiscreteProcess(theNeutronElasticProcess);

            }

            */
        }

        // Deexcitation
        G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
        G4LossTableManager::Instance()->SetAtomDeexcitation(de);

    }
    else if(ParticlePysics == "HADRON_FTFP_BERT" || ParticlePysics == "HADRON_FTFP_BERT_ATL" || ParticlePysics == "HADRON_FTFP_BERT_TRV"  || ParticlePysics == "HADRON_QGSP_FTFP_BERT" || ParticlePysics == "HADRON_QGSP_BERT" || ParticlePysics == "HADRON_QGSP_BERT_HP" || ParticlePysics == "HADRON_QGSP_BIC" || ParticlePysics == "HADRON_QGSP_BIC_AllHP" || ParticlePysics == "HADRON_INCLXX" || ParticlePysics == "HADRON_Shielding" || ParticlePysics == "HADRON_ShieldingLEND"){
        G4VEmPhysicsConstructorObj->ConstructProcess();
        G4VHadromPhysicsConstructorObj->ConstructProcess();
    }
    else {
        G4Exception("MyPhysicsList::MyPhysicsList", "InvalidSetup", FatalException, "Invalid physics choice");
    }

    if(TestPointsPositions == "yes"){

    }else{

        AddTransportation();
        decPhysicsList->ConstructProcess();
        RadioactivedecPhysicsList->ConstructProcess();
    }
}

void G4TUserPhysicsList::SetCuts()
{

    //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(CutsEnergy, 1*GeV);

    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;

    if(CutInRangeData == ""){
        SetCutsWithDefault();
    }else{
        std::istringstream LineString(CutInRangeData);
        //G4cout << CutInRangeData << G4endl;

        G4double Ene; G4String pn, Unit;
        while(LineString >> pn ){
            LineString >> Ene;
            LineString >> Unit;
            //G4cout << " pn " << pn << " Ene " << Ene << " Unit " << Unit << G4endl;
            if((pn == "e-" || pn == "e+" || pn == "gamma" || pn == "proton") && Ene*G4UnitDefinition::GetValueOf(Unit) != 0.){
                SetCutValue(Ene*G4UnitDefinition::GetValueOf(Unit), pn);
            }
        }
    }

    //G4cout << EnergyThresholdsData << G4endl;

    if(EnergyThresholdsData != ""){
        std::istringstream LineString(EnergyThresholdsData);
        //G4cout << EnergyThresholdsData << G4endl;
        std::string ww; std::vector<std::string> sl;

        while(LineString >> ww ){
            sl.push_back(ww);
        }

        G4double lowLimit = 1*keV;
        G4double highLimit = 100. * GeV;
        if(sl.size()>1){
            G4double lowLimit = std::stod(sl[0])*G4UnitDefinition::GetValueOf(sl[1]);
            if(sl.size()>3){
                highLimit = std::stod(sl[2])*G4UnitDefinition::GetValueOf(sl[3]);
            }
            //G4cout << CutInRangeData << " - lowLimit " << lowLimit << " highLimit " << highLimit << G4endl;
            if(lowLimit != 0. && highLimit > lowLimit ){
                G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);
            }
        }
    }

/*
    defaultCutValue = CutsDistance;

    // set distance cuts for secondary production

    //QString InputsVals = CutsEnergy.split(QRegExp("(\\s|\\n|\\r)+"), QString::SkipEmptyParts); // "/SourceData/setSourceGenerationData"

    if(IsDcutsSet == true && CutsDistance > 0.){
        defaultCutValue = CutsDistance;
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n IsDcutsSet == true " << CutsDistance<< G4endl;
    }

    if (IsEcutsSet == true && (IsDcutsSet == false )) {
        defaultCutValue = CutsEnergy;
        //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n IsEcutsSet == true " << CutsEnergy<< G4endl;
    }

    SetCutValue(defaultCutValue, "gamma");
    SetCutValue(defaultCutValue, "e-");
    SetCutValue(defaultCutValue, "e+");
    SetCutValue(100*keV, "proton");
    SetCutValue(defaultCutValue, "alpha");
    SetCutValue(100*keV, "neutron");

    // SetCutsWithDefault();
    // Set the secondary production cut lower than 990. eV,  Very important for high precision of lowenergy processes at low energies

    G4double lowLimit = CutsEnergy ;
    G4double highLimit = 100. * GeV;
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);

    // Set the tracking cut for neutrons
    G4double neutronCut = 0.5 * mm;
    SetCutValue(neutronCut, "neutron");
    G4Region* region = G4RegionStore::GetInstance()->GetRegion("DefaultRegionForTheWorld");
    G4ProductionCuts* cuts = new G4ProductionCuts;
    cuts->SetProductionCut(neutronCut);
    region->SetProductionCuts(cuts);
*/

    ShowSourceParameters();

#if VERBOSE_USE
    G4cout<<"\n\n========= Geant4 Cross Section Tables Generation ====================\n\n" <<G4endl;

    if (verboseLevel>0) DumpCutValuesTable();
#endif
}

void G4TUserPhysicsList::ShowSourceParameters() {

#if VERBOSE_USE
    G4cout<<"\n\n========= Physics Inputs ====================\n\n" <<G4endl;
    G4cout<<" >> Physics to be used: " << ParticlePysics <<G4endl;

    if(ParticlePysics =="Construct"){
        G4cout<<" >> New PhotoElectricEffectModel " << PhotoElectricEffectModel <<G4endl;
        G4cout<<" >> New PolarizedPhotoElectricEffectModel " << PolarizedPhotoElectricEffectModel <<G4endl;
        G4cout<<" >> New ComptonScatteringModel " << ComptonScatteringModel <<G4endl;
        G4cout<<" >> New PolarizedComptonModel " << PolarizedComptonModel <<G4endl;
        G4cout<<" >> New GammaConversionModel " << GammaConversionModel <<G4endl;
        G4cout<<" >> New PolarizedGammaConversionModel " << PolarizedGammaConversionModel <<G4endl;
        G4cout<<" >> New ParticlePysics " << RayleighScatteringModel <<G4endl;
        G4cout<<" >> New GammaConversionToMuonModel " << GammaConversionToMuonModel <<G4endl;
    }
    G4cout<<" >> Cuts In Range: " << CutInRangeData << G4endl;
    G4cout<<" >> Energy Range: " << EnergyThresholdsData << G4endl;

    if(GenerateCrossSectionTableFlag){
        G4cout<<" >> GenerateCrossSectionTableFlag " << "yes" <<G4endl;
    }
#endif

}
