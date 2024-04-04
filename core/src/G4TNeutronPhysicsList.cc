#include "G4TNeutronPhysicsList.hh"


#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
//#include "G4ProcessTable.hh"

// Processes

#include "G4HadronElasticProcess.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPThermalScatteringData.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPThermalScattering.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPInelastic.hh"

#include "G4NeutronCaptureProcess.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPCapture.hh"

#include "G4NeutronFissionProcess.hh"
#include "G4ParticleHPFissionData.hh"
#include "G4ParticleHPFission.hh"

#include "G4SystemOfUnits.hh"

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

#include "G4DecayPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"

#include "G4UserLimits.hh"
#include "G4Electron.hh"

#include "G4UserSpecialCuts.hh"
#include "G4Neutron.hh"

extern G4String ParticlePysics;

G4TNeutronPhysicsList::G4TNeutronPhysicsList(const G4String& name):G4VModularPhysicsList(), fThermal(false)
{
    //SetVerboseLevel(1);

}

G4TNeutronPhysicsList::~G4TNeutronPhysicsList()
{

}

void G4TNeutronPhysicsList::ConstructParticle()
{

    G4Neutron::NeutronDefinition();

    //G4cout << "from function : " << __FUNCTION__<< G4endl;

    //G4BosonConstructor* Boson = new G4BosonConstructor();
    //G4LeptonConstructor* Lepton = new G4LeptonConstructor();
    //G4MesonConstructor* Meson = new G4MesonConstructor();
    //G4BaryonConstructor* Baryon = new G4BaryonConstructor();
    //G4IonConstructor* Ion = new G4IonConstructor();
    //G4ShortLivedConstructor* ShortLived = new G4ShortLivedConstructor();
    //Boson->ConstructParticle();
    //Lepton->ConstructParticle();
    //Meson->ConstructParticle();
    //Baryon->ConstructParticle();
    //Ion->ConstructParticle();
    //ShortLived->ConstructParticle();

    SetCutsWithDefault();

    //SetDefaultCutValue();
    /*
    SetCutValue(500*mm,"neutron");
    SetCutValue(1*mm,"proton");
    SetCutValue(0.001*mm,"e-");
    SetCutValue(1*mm,"e+");
    SetCutValue(1*mm,"gamma");
    SetCutValue(0.05*m,"alpha");
*/

}

void G4TNeutronPhysicsList::ConstructProcess()
{

    G4int verbose = 0;

    G4ParticleDefinition* neutron = G4Neutron::Neutron();
    G4ProcessManager* pManager = neutron->GetProcessManager();

    // delete all neutron processes if already registered

    G4VProcess* process = 0;
    process = pManager->GetProcess("hadElastic");

    if (process) pManager->RemoveProcess(process);
    process = pManager->GetProcess("neutronInelastic");
    if (process) pManager->RemoveProcess(process);

    process = pManager->GetProcess("nCapture");
    if (process) pManager->RemoveProcess(process);

    process = pManager->GetProcess("nFission");
    if (process) pManager->RemoveProcess(process);

    // (re) create process: elastic
    G4HadronElasticProcess* process1 = new G4HadronElasticProcess();
    pManager->AddDiscreteProcess(process1);

    // model1a
    G4ParticleHPElastic*  model1a = new G4ParticleHPElastic();
    process1->RegisterMe(model1a);
    process1->AddDataSet(new G4ParticleHPElasticData());

    // model1b
    if (fThermal) {
        model1a->SetMinEnergy(4*eV);
        G4ParticleHPThermalScattering* model1b = new G4ParticleHPThermalScattering();
        process1->RegisterMe(model1b);
        process1->AddDataSet(new G4ParticleHPThermalScatteringData());
    }

    // (re) create process: inelastic
    G4HadronInelasticProcess* process2 = new G4HadronInelasticProcess( "neutronInelastic", G4Neutron::Definition() );
    pManager->AddDiscreteProcess(process2);

    // cross section data set
    G4ParticleHPInelasticData* dataSet2 = new G4ParticleHPInelasticData();
    process2->AddDataSet(dataSet2);

    // models
    G4ParticleHPInelastic* model2 = new G4ParticleHPInelastic();
    process2->RegisterMe(model2);


    // (re) create process: nCapture
    G4NeutronCaptureProcess* process3 = new G4NeutronCaptureProcess();
    pManager->AddDiscreteProcess(process3);

    // cross section data set
    G4ParticleHPCaptureData* dataSet3 = new G4ParticleHPCaptureData();
    process3->AddDataSet(dataSet3);

    // models
    G4ParticleHPCapture* model3 = new G4ParticleHPCapture();
    process3->RegisterMe(model3);

    // (re) create process: nFission

    G4NeutronFissionProcess* process4 = new G4NeutronFissionProcess();
    pManager->AddDiscreteProcess(process4);

    // cross section data set
    G4ParticleHPFissionData* dataSet4 = new G4ParticleHPFissionData();
    process4->AddDataSet(dataSet4);

    // models
    G4ParticleHPFission* model4 = new G4ParticleHPFission();
    process4->RegisterMe(model4);



    //G4ProcessManager* pmanager = G4Neutron::NeutronDefinition()->GetProcessManager();
    //pmanager->AddProcess(new G4UserSpecialCuts(),-1,-1,1);

    //pmanager = G4Proton::ProtonDefinition()->GetProcessManager();
    //pmanager->AddProcess(new G4UserSpecialCuts(),-1,-1,1);

    //G4VPhysicsConstructor* G4VEmPhysicsConstructorObj = new G4EmStandardPhysics(verbose);

    //if (ParticlePysics == "EMS") {
    //    G4VEmPhysicsConstructorObj = new G4EmStandardPhysics(verbose);
    //}
    //else if(ParticlePysics == "EMS1"){
    //    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n EMS1 Physics " << G4endl;

    //    G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option1(verbose);
    //}
    //else if(ParticlePysics == "EMS2"){
    //    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n EMS2 Physics " << G4endl;

    //    G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option2(verbose);
    //}
    //else if(ParticlePysics == "EMS3"){
    //    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n EMS3 Physics " << G4endl;

    //    G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option3(verbose);
    //}
    //else if(ParticlePysics == "EMS4"){
    //    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n EMS4 Physics " << G4endl;

    //    G4VEmPhysicsConstructorObj = new G4EmStandardPhysics_option4(verbose);
    //}
    //else if(ParticlePysics == "Livermore"){
    //    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n Livermore Physics " << G4endl;

    //    G4VEmPhysicsConstructorObj = new G4EmLivermorePhysics(verbose);
    //}
    //else if(ParticlePysics == "Penelope"){
    //    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n penelope Physics " << G4endl;

    //    G4VEmPhysicsConstructorObj = new G4EmPenelopePhysics(verbose);
    //}
    //G4VEmPhysicsConstructorObj->ConstructProcess();

    AddTransportation();

}

void G4TNeutronPhysicsList::SetCuts()
{

    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;
/*
    G4double electronCut = 1.0 * keV;
    G4double protonCut = 10.0 * keV;
    G4double alphaCut = 100.0 * keV;

    G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
    G4double defaultCut = 1.0 * m; // default cut value
    theCoupleTable->SetEnergyRange(defaultCut, 100. * GeV);

    // Set cuts for electrons
    G4ProductionCuts* ecut = new G4ProductionCuts();
    ecut->SetProductionCut(electronCut, G4Electron::Definition());

    // Set cuts for protons
    G4ProductionCuts* pcuts = new G4ProductionCuts();
    pcuts->SetProductionCut(protonCut, G4Proton::Definition());

    // Set cuts for alphas
    G4ProductionCuts* acuts = new G4ProductionCuts();
    acuts->SetProductionCut(alphaCut, G4Alpha::Definition());
*/

/*
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ProductionCutsTable* cutsTable = G4ProductionCutsTable::GetProductionCutsTable();
    G4ProductionCuts* cuts = new G4ProductionCuts();
    // Set the cut value for each particle type
    for (G4int j = 0; j < particleTable->size(); ++j) {
        G4ParticleDefinition* particle = particleTable->GetParticle(j);
        G4String particleName = particle->GetParticleName();
        cuts->SetProductionCut(1*m, particle);
    }
*/

/*
  SetCutValue(1*m,"neutron");  SetApplyCuts(true,"neutron");
  SetCutValue(0.5*m,"proton");
  SetCutValue(0.1*m,"e-");
  SetCutValue(0.1*m,"e+");
  SetCutValue(0.1*m,"gamma");
  SetCutValue(0.05*m,"alpha");
*/
}

