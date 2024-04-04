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
// TETRunAction.cc
// \file   MRCP_GEANT4/Internal/src/TETRunAction.cc
// \author Haegin Han
//

#include "TETRunAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Run.hh"
#include "G4TVolumeConstruction.hh"
#include "G4TUserPhysicsList.hh"
#include "G4RunManager.hh"
#include "G4TReadPrimaryGeneratorAction.hh"
#include "G4TDirectToFilesPrimaryGeneratorAction.hh"
#include "G4TDirectPrimaryGeneratorAction.hh"
#include "G4TDirectVoxelsPrimaryGeneratorAction.hh"

#include <stdio.h>
#include <string.h>

#include "G4Threading.hh"

#include "G4EmCalculator.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTable.hh"
#include "G4Timer.hh"

extern G4String ParticleName;
extern G4String GeometrySymbol;
extern G4String SourceType;
extern G4String EnergyDistribution;
extern G4String UseGeneratedData ;
extern G4String ResultDirectoryPath ;
extern G4String GenerateVoxelsResuls;
extern G4String DataDirectoryPath;

extern std::vector<G4String> NewRankSourceParticlesNamesValues;
extern std::vector<G4String> NewRankSourceRegionsNamesValues;
extern std::vector<G4double> NewRankSourceEnergiesValues;
extern std::vector<G4String> NewRankSourceMomDirsValues;

#ifdef G4MULTITHREADED
G4ThreadLocal std::map<G4int,std::map<unsigned int,G4double>> TETRunAction::VoxelsED_Total ;
G4ThreadLocal std::map<G4int,std::map<unsigned int,G4double>> TETRunAction::VoxelsED2_Total ;
G4ThreadLocal std::map<G4int,std::map<unsigned int,unsigned long long int>> TETRunAction::VoxelsNOfValues ;
G4ThreadLocal std::map<G4int,std::map<G4String,unsigned long long int>> TETRunAction::NOfValues;
G4ThreadLocal std::map<G4int,std::map<G4String,G4double>> TETRunAction::ED_Total ;
G4ThreadLocal std::map<G4int,std::map<G4String,G4double>> TETRunAction::ED2_Total ;
G4ThreadLocal std::map<G4int,std::map<G4String,G4double>> TETRunAction::Fluence ;
G4ThreadLocal G4int TETRunAction::rank, TETRunAction::thread, TETRunAction::DataID, TETRunAction::EventIndex, TETRunAction::NumberOfRanksThreads, TETRunAction::TotalEventNumber;
G4ThreadLocal G4double TETRunAction::EnergyEmittedPerThread, TETRunAction::ParticleSourceEnergy, TETRunAction::ExecutionTimeInMin, TETRunAction::OneEventExecutionTimeInMs;
G4ThreadLocal std::chrono::steady_clock::time_point TETRunAction::start, TETRunAction::end ;
G4ThreadLocal TETRun* TETRunAction::fRun;
#endif


TETRunAction::TETRunAction()
//:fRun(0)
{}

TETRunAction::~TETRunAction()
{}

G4Run* TETRunAction::GenerateRun()
{
	// generate run
	fRun = new TETRun();
	return fRun;
}

void TETRunAction::BeginOfRunAction(const G4Run* aRun)
{
    rank = 0 , thread = 0 ;

#ifdef G4MPI_USE
    thread = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &DataID);
    //DataID = G4MPImanager::GetManager()->GetRank();
    rank = DataID;
    /*
    if(G4Threading::IsMultithreadedApplication()){ // normal multiThreaded mode
        thread = G4Threading::G4GetThreadId() ;
    }
    else {
        thread = 0;
    }
    */
#else
    if(G4Threading::IsMultithreadedApplication()){ // normal multiThreaded mode
        rank = 0 ;
        thread = G4Threading::G4GetThreadId() ;
        DataID = thread;
    }
    else{ // Sequential mode
        rank = 0;
        thread = 0 ;
    }
#endif

    if(OneOrMultiSimulations == "o"){
        DataID = 0;
    }

    GenerateCrossSectionGraph = "no";

    //G4int run_number = aRun->GetRunID();
    TotalEventNumber = aRun->GetNumberOfEventToBeProcessed(); // for worker thread(number of event in worker thread), for master thread (number of event in all worker thread)
    //TotalEventNumber = aRun->GetNumberOfEvent();

    //std::ostringstream sd; sd << NewRankSourceParticlesNamesValues[DataID] << " " << SourceType << " " << NewRankSourceRegionsNamesValues[DataID] << " " << NewRankSourceEnergiesValues[DataID] << " " <<NewRankSourceMomDirsValues[DataID] ;
    //G4String srcdata = sd.str();
    //G4cout << srcdata << G4endl;

    //G4String srcdata = "";

    NumberOfRanksThreads = 1 ;
    G4Timer* t; G4String Time = t->GetClockTime();

#ifdef G4MPI_USE // if MPI are used then the multiThreaded mode automatically is activated.

    ExecutionMode = "MPI";

    MPI_Comm_size(MPI_COMM_WORLD, &NumberOfRanksThreads);
    //NumberOfRanksThreads = G4MPImanager::GetManager()->GetTotalSize();

    std::cout <<  "\n " << Time << " ========= MPI mode : "<<__FUNCTION__<<" from Rank " << DataID << "/" <<NumberOfRanksThreads << ". Start of " << TotalEventNumber << " events simulation loop using " << GeometrySymbol <<" geometry -Source Data: " << NewRankSourceParticlesNamesValues[DataID] << ", " << SourceType << ", " << NewRankSourceRegionsNamesValues[DataID] << ", "  << EnergyDistribution << " " << NewRankSourceEnergiesValues[DataID] << ", " << NewRankSourceMomDirsValues[DataID] << std::endl;

    if(G4Threading::IsMultithreadedApplication()){
        if(G4Threading::IsMasterThread()){
            WriteMacroscopicCrossSection();
        }
    }else{
        if(DataID == 0){
            WriteMacroscopicCrossSection();
        }
    }

#else

    if(G4Threading::IsMultithreadedApplication()){ // normal multiThreaded mode

        ExecutionMode = "MT";
        NumberOfRanksThreads = G4Threading::GetNumberOfRunningWorkerThreads() ;
        TotalEventNumber = TotalEventNumber/G4Threading::GetNumberOfRunningWorkerThreads();
        //TotalEventNumber = aRun->GetNumberOfEvent();

        if(G4Threading::IsWorkerThread()){
            std::cout <<  "\n " << Time << " ========= MT mode : "<<__FUNCTION__<<" from Worker Thread " << G4Threading::G4GetThreadId() << "/" <<G4Threading::GetNumberOfRunningWorkerThreads() << ". Start of " << TotalEventNumber << " events simulation loop using " << GeometrySymbol <<" geometry -Source Data: " << NewRankSourceParticlesNamesValues[DataID] << ", " << SourceType << ", " << NewRankSourceRegionsNamesValues[DataID] << ", "  << EnergyDistribution << " " << NewRankSourceEnergiesValues[DataID] << ", " << NewRankSourceMomDirsValues[DataID] << std::endl ;
        }
        else if(G4Threading::IsMasterThread()){
            std::cout <<  "\n " << Time << " ========= MT mode : "<<__FUNCTION__<<" from Master Thread " << G4Threading::G4GetThreadId() << ". Start of " << TotalEventNumber << " events simulation loop " << std::endl ;
            WriteMacroscopicCrossSection();
        }
    }
    else { // not MPI and not normal MultiThreading mode then it's sequencial mode
        ExecutionMode = "SEQ";
        NumberOfRanksThreads = 1;
        std::cout <<  "\n " << Time << " ========= SEQ mode : "<<__FUNCTION__<<" from Master" << ". Start of " << TotalEventNumber << " events simulation loop " << " -Source Data: " << NewRankSourceParticlesNamesValues[DataID] << ", " << SourceType << ", " << NewRankSourceRegionsNamesValues[DataID] << ", "  << EnergyDistribution << " " << NewRankSourceEnergiesValues[DataID] << ", " << NewRankSourceMomDirsValues[DataID] << std::endl;
    }

#endif

    // to get the values of cmd for primaryGenerator from the humainPhantomConstruction that is a global class(accept cmd to simulate all app) not like primary is Threadlocal class the it can't accept cmd
    // we need to define it in a function and not in the constructor, because always the constructors are called in the first for all classes then it's data its not the updated(values of messenger by default ) then when we use it in function , all the data are apdated ( from messenger by user cmd) and we can use the updated
    const G4TVolumeConstruction* TConstruction2 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    OneOrMultiSimulations = TConstruction2->getMPISimulationNum();


    // to initialize the arrays ED_Total ;
    OrgansNameVector = TConstruction2->GetOrganNamesVector();
    for(G4int nn = 0;nn < (G4int)OrgansNameVector.size();nn++){
        ED_Total[DataID][OrgansNameVector[nn]]=0.;
        ED2_Total[DataID][OrgansNameVector[nn]]=0.;
        Fluence[DataID][OrgansNameVector[nn]]=0.;
    }
    ED_Total[DataID]["World"]=0.;

    EventIndex = 0;

    useTime(0); // "0" for setting the begining time befor a loop, "1" it get the end time and calculate consumed time by a loop and show the difference, the loop should not have a setting or getting time inside it

    G4RunManager::GetRunManager()->SetPrintProgress(TotalEventNumber*0.1);

    //std::cout << "\n\n\n\n\n\n 11111111111111 \n\n\n\n\n\n "<< std::endl;

}



void TETRunAction::EndOfRunAction(const G4Run* aRun)
{

    const G4TVolumeConstruction* TConstruction2 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    G4Timer* t; G4String Time = t->GetClockTime();
    std::cout << __FUNCTION__ << " : " << Time << " " ;
    useTime(1); // "0" for setting the begining time befor a loop, "1" it get the end time and calculate consumed time by a loop and show the difference, the loop should not have a setting or getting time inside it

    if(UseGeneratedData == "read"){
        const G4TReadPrimaryGeneratorAction* primaryGeneratorAction = static_cast<const G4TReadPrimaryGeneratorAction*> (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
        EnergyEmittedPerThread = primaryGeneratorAction->getEmittedEnergy();
    }
    else if(UseGeneratedData == "save"){
        const G4TDirectToFilesPrimaryGeneratorAction* primaryGeneratorAction = static_cast<const G4TDirectToFilesPrimaryGeneratorAction*> (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
        EnergyEmittedPerThread = primaryGeneratorAction->getEmittedEnergy();
    }
    else{
        if(SourceType == "Voxels"){
            const G4TDirectVoxelsPrimaryGeneratorAction* primaryGeneratorAction = static_cast<const G4TDirectVoxelsPrimaryGeneratorAction*> (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
            EnergyEmittedPerThread = primaryGeneratorAction->getEmittedEnergy();
        }else {
            const G4TDirectPrimaryGeneratorAction* primaryGeneratorAction = static_cast<const G4TDirectPrimaryGeneratorAction*> (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
            EnergyEmittedPerThread = primaryGeneratorAction->getEmittedEnergy();
        }
    }

    G4String distriName = EnergyDistribution;
    if(distriName == "Mono"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[DataID] ;
    }
    else if(distriName == "Uniform"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[DataID] ;
    }
    else if(distriName == "Rayleigh"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[DataID];
    }
    else if(distriName == "Gauss"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[DataID];
    }
    else if(distriName == "Spectrum"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[DataID];
    }else if(distriName == "File"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[DataID];
    }
    else if(distriName == "RadioNuclide"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[DataID];
    }

#ifdef G4MULTITHREADED
    if(G4Threading::IsMasterThread()){
        return;
    }
#endif

    EDEPMAP edepMap = *fRun->GetEdepMap();
    NumStepMAP numstepsMap = *fRun->GetnumstepsMap();

    MatIDNameMap = TConstruction2->getMaterialIDName();

    for ( auto it3 = edepMap.begin(); it3 != edepMap.end(); ++it3  ){
        ED_Total [DataID][MatIDNameMap[it3->first]] = it3->second.first;
        ED2_Total[DataID][MatIDNameMap[it3->first]] = it3->second.second;
        NOfValues[DataID][MatIDNameMap[it3->first]] = numstepsMap[it3->first];
        //std::cout << " | " << it3->first << " | "  << MatIDNameMap[it3->first]  << " | " << ED_Total[DataID][MatIDNameMap[it3->first]]  << " | "  << ED2_Total[DataID][MatIDNameMap[it3->first]] << " | "  << NOfValues[DataID][MatIDNameMap[it3->first]] << std::endl;
    }

    CreateThreadRegionResultFile();

}
void TETRunAction::FillRegionLenghts(G4String StepRegion, G4double StepEnergy){

    //G4cout << N << " " << E << " "<< E*E<< G4endl;

    Fluence[DataID][StepRegion] += StepEnergy ;

    //G4cout << DataID << " " << StepRegion << " " << Fluence[DataID][StepRegion] << G4endl;
}
void TETRunAction::CreateThreadRegionResultFile(){

    const G4TVolumeConstruction* TConstruction2 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    std::ostringstream fname;

    fname << ResultDirectoryPath <<"/AE@for@Rank@"<<rank<<"@Thread@"<<thread << "@" << GeometrySymbol << "@" << NewRankSourceParticlesNamesValues[DataID] << "@" << NewRankSourceRegionsNamesValues[DataID] << "@" << ParticleSourceEnergy ;
    //fname << ResultDirectoryPath <<"/AE@for@Rank@"<<rank<<"@Thread@"<<thread << "@" << NewRankSourceParticlesNamesValues[DataID] << "@" << NewRankSourceRegionsNamesValues[DataID] << "@" << ParticleSourceEnergy ;
    std::ofstream file(fname.str().c_str(), std::ios_base::binary);

    std::cout << " Creating region result file of worker thread or rank " << fname.str().c_str()<< std::endl;

    G4int SZ = 30;

    //file << "RegionName" << "\n" ;
    file << std::setw(SZ) << std::left << "RegionName" << " "
    << std::setw(15) << std::left << "AE(MeV)" << " "
    << std::setw(15) << std::left << "AE^2(MeV^2)" << " "
    << std::setw(15) << std::left << "Steps Number" << " "
    << std::setw(20) << std::left << "Mass(kg)"<< " "
    << std::setw(20) << std::left << "Volume(cm3)"<< " "
    << std::setw(20) << std::left << "Density(g/cm3)" << " "
    << std::setw(20) << std::left << "Fluence(1/cm2)" << "\n";

    for(G4int jk = 0 ; jk < (G4int)OrgansNameVector.size() ; jk++){ // each line for an organ name
        file << std::setw(SZ) << std::left << OrgansNameVector[jk] << " "
             << std::setw(15) << std::left << ED_Total[DataID][OrgansNameVector[jk]] << " "
             << std::setw(15) << std::left << ED2_Total[DataID][OrgansNameVector[jk]] << " "
             << std::setw(15) << std::left << NOfValues[DataID][OrgansNameVector[jk]] << " "

             << std::setw(20) << std::left << TConstruction2->GetOrganNameMassMap()[OrgansNameVector[jk]]<< " "
             << std::setw(20) << std::left << TConstruction2->GetOrganNameVolumeMap()[OrgansNameVector[jk]]/1000<< " "
             << std::setw(20) << std::left << TConstruction2->GetOrganNameDensityMap()[OrgansNameVector[jk]]*1000<< " ";
        if(TConstruction2->GetOrganNameVolumeMap()[OrgansNameVector[jk]] == 0. || __isinf(TConstruction2->GetOrganNameVolumeMap()[OrgansNameVector[jk]]) || __isnan(TConstruction2->GetOrganNameVolumeMap()[OrgansNameVector[jk]])){
            file << std::setw(20) << std::left << 0 << " ";
        }else{
            file << std::setw(20) << std::left << Fluence[DataID][OrgansNameVector[jk]]/(TConstruction2->GetOrganNameVolumeMap()[OrgansNameVector[jk]]/1000) << " ";
        }

        file << "\n";
        //G4cout  << OrgansNameVector[jk] << "  " << ED_Total[DataID][OrgansNameVector[jk]] << "  " << G4endl;
    }

    SZ++;

    file << "\n";
    file << std::setw(SZ) << std::left << "TotalEventNumber " << TotalEventNumber << "\n" ;
    file << std::setw(SZ) << std::left << "EnergyEmittedPerThread " << EnergyEmittedPerThread << "\n" ;
    file << std::setw(SZ) << std::left << "ExecutionTimeInMin " << ExecutionTimeInMin << "\n" ;
    OneEventExecutionTimeInMs = (double)(ExecutionTimeInMin/TotalEventNumber)*(1000./60.);
    file << std::setw(SZ) << std::left << "OneEventExecutionTimeInMs " << OneEventExecutionTimeInMs << "\n" ;
    file << std::setw(SZ) << std::left << "ParticleName " << NewRankSourceParticlesNamesValues[DataID] <<  "\n" ;
    file << std::setw(SZ) << std::left << "SourceType " << SourceType <<  "\n" ;
    file << std::setw(SZ) << std::left << "SourceRegionName " << NewRankSourceRegionsNamesValues[DataID] <<  "\n" ;
    file << std::setw(SZ) << std::left << "EnergyDistribution " << EnergyDistribution <<  "\n" ;
    file << std::setw(SZ) << std::left << "ParticleSourceEnergy " << ParticleSourceEnergy <<  "\n" ;
    file << std::setw(SZ) << std::left << "MomDirDistribution " << NewRankSourceMomDirsValues[DataID] <<  "\n" ;
    file << std::setw(SZ) << std::left << "RankID " << rank <<  "\n" ;
    //file << std::setw(SZ) << std::left << "RankID " << TConstruction2->getGeometrySymbol() <<  "\n" ;

    file.close();
}
void TETRunAction::useTime(G4int set_get_for)
{

    if( set_get_for == 0){
        start = std::chrono::steady_clock::now();
    }

    else if(set_get_for == 1){

        end = std::chrono::steady_clock::now();
        auto diff = end - start;

        std::ostringstream text;

        G4double unitFactor = 1;
        G4String unit = " ms";
        G4double timeConsumption = std::chrono::duration < G4double, std::milli > (diff).count();

        if(timeConsumption <= 1000.){
            unit = " ms";
        }
        else if(1000. < timeConsumption && timeConsumption <= 60000.){ ;
            unitFactor = 0.001;
            unit = " s";
        }
        else if(60000. < timeConsumption && timeConsumption <= 3600000.){
            unitFactor = 0.001/60.;
            unit = " min";
        }
        else if(3600000. < timeConsumption){
            unitFactor = 0.001/(60.*60.);
            unit = " h";
        }

        ExecutionTimeInMin = timeConsumption*(0.001/60.);

        std::cout << " Events Simulation time was: " << timeConsumption * unitFactor << unit << ".";  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.

        //printf(" Time taken : %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
        //printf("\n");
        //printf(text.str().c_str());

    }
}
void TETRunAction::WriteMacroscopicCrossSection(){

    const G4TVolumeConstruction* PH = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    //const G4TUserPhysicsList* PH = static_cast<const G4TUserPhysicsList*> (G4RunManager::GetRunManager()->GetUserPhysicsList());
    if(PH->getGenerateCrossSectionTableFlag() == false) {
        GenerateCrossSectionGraph = "no";
        return;
    }
    GenerateCrossSectionGraph = "yes"; // used to be written in simulation file
#if VERBOSE_USE
    G4cout<<"\n\n========= Cross Section Data ====================" <<G4endl;
#endif

    G4double fRangeCut[3];
    G4double fEnergyCut[3];

    //G4cout<<"\n\n========= 1 ====================" <<G4endl;

    //instanciate EmCalculator
    G4EmCalculator emCal;
    // emCal.SetVerbose(2);
    G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();

    //G4cout<<"\n\n========= 2 ====================" <<G4endl;

    G4String FileName = ResultDirectoryPath + "/CrossSectionData";

    std::vector<G4double> Energies;
    G4String partName;

    Energies = PH->getEnergiesForCrossSectionValues();
    std::sort(Energies.begin(), Energies.end());
    partName = PH->getParticleForCrossSectionValues();

    //G4cout<<"\n\n========= 3 ====================" <<G4endl;

    G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle(partName);
    G4double charge   = particle->GetPDGCharge();
    G4double energy;
    G4int nee = 0;

    //G4cout<<"\n\n========= 4 ====================" <<G4endl;

    // get material
    G4Material* material;
    std::vector<G4Material*>* Mats = G4Material::GetMaterialTable();
    std::vector<G4String> MYMaterialsNames;
    std::vector<G4String> DoneMaterialsNames;

    //G4cout<<"\n\n========= 5 ====================" <<G4endl;

    //G4cout << "\nNumber of registered materials " << Mats->size() ;

    // to ignore all the repeated materials with the same names and compositions that can produce an error
    for (size_t iLV = 0; iLV < Mats->size(); iLV++ ) {

        material = (*Mats)[iLV];
        G4String matName     = material->GetName();

        //G4cout<<"iLV: " << iLV<<" ElemNum " << (*Mats)[iLV]->GetNumberOfElements() <<" matName: " << matName <<G4endl;

        if((*Mats)[iLV]->GetNumberOfElements() > 1 ){


            //G4cout<<"iLV: " << iLV<<" ElemNum " << (*Mats)[iLV]->GetNumberOfElements() <<" matName: " << matName <<G4endl;

            bool IsIn = false;
            for (size_t gg = 0 ; gg < MYMaterialsNames.size() ; gg++) {
                if(MYMaterialsNames[gg] == matName){
                    IsIn = true;
                    //G4cout<<"----------- IsIn: " << IsIn << " " << matName << G4endl;
                    break;
                }
            }
            if(IsIn == false){
                MYMaterialsNames.push_back(matName);
                //G4cout<<"----------- IsIn: " << IsIn << " Saved " << matName << G4endl;

            }
        }
    }

    //G4cout << "\nNumber of registered and filtered materials " << MYMaterials->size() ;
    for (size_t iLV = 0; iLV < Mats->size(); iLV++ ) {

        material = (*Mats)[iLV];
        G4String matName     = material->GetName();

        bool IsIn = false;
        for (size_t gg = 0 ; gg < MYMaterialsNames.size() ; gg++) {
            if(MYMaterialsNames[gg] == matName){
                IsIn = true;
                //G4cout<<"----------- IsIn: " << IsIn << " " << matName << G4endl;
                break;
            }
        }
        if((*Mats)[iLV]->GetNumberOfElements() > 1 && IsIn == true){}
        else{continue;}

        IsIn = false;
        for (size_t gg = 0 ; gg < DoneMaterialsNames.size() ; gg++) {
            if(DoneMaterialsNames[gg] == matName){
                IsIn = true;
                //G4cout<<"----------- IsIn: " << IsIn << " " << matName << G4endl;
                break;
            }
        }
        if(IsIn == false){
            DoneMaterialsNames.push_back(matName);
        }else{
            continue;
        }


        G4double density     = material->GetDensity();
        G4double radl        = material->GetRadlen();

        G4cout << "\nTable " << partName << " " << material->GetName() << " Density: " << G4BestUnit(density,"Volumic Mass") << " RadiationLength "<< G4BestUnit(radl, "Length") << "Sigma in (cm-1), Range Cut in (mm), Energy Cut in (MeV)" <<  "\n\n" ;

        std::ostringstream LatexText;
        std::ostringstream TableText;
#if VERBOSE_USE
        //G4cout << "\nTable " << partName << " " << material->GetName() << " Density: " << G4BestUnit(density,"Volumic Mass") << " RadiationLength "<< G4BestUnit(radl, "Length") << "Sigma in (cm-1), Range Cut in (mm), Energy Cut in (MeV)" <<  "\n\n" ;
#endif
        TableText << "\nTable " << partName << " " << material->GetName() << " Density: " << G4BestUnit(density,"Volumic Mass") << " RadiationLength "<< G4BestUnit(radl, "Length") << "Sigma in (cm-1), Range Cut in (mm), Energy Cut in (MeV)" <<  "\n\n" ;
        LatexText << "Latex Table \n\n";
        LatexText << "\\begin{table}[H] \n"
                  << "\\centering \n"
                  << "\\caption{Cross section per Volume (cm$^{-1}$) for "<< partName <<" interaction with material "<< matName << " of density " << G4BestUnit(density,"Volumic Mass") << ". } \n"
                  << "\\begin{tabular}{";


        //if( (*Mats)[iLV]->GetNumberOfMaterials() == 1){}

        nee = 0;
        for (G4int e = 0; e < Energies.size(); e++ ) {

            energy = Energies[e];

            // get cuts
            size_t numOfCouples = theCoupleTable->GetTableSize();
            const G4MaterialCutsCouple* couple = 0;
            G4int index = 0;
            for (size_t i=0; i<numOfCouples; i++) {
                couple = theCoupleTable->GetMaterialCutsCouple(i);
                if (couple->GetMaterial() == material) {index = i; break;}
            }
            fRangeCut[0] = (*(theCoupleTable->GetRangeCutsVector(idxG4GammaCut)))[index]/mm;
            fRangeCut[1] = (*(theCoupleTable->GetRangeCutsVector(idxG4ElectronCut)))[index]/mm;
            fRangeCut[2] = (*(theCoupleTable->GetRangeCutsVector(idxG4PositronCut)))[index]/mm;
            fEnergyCut[0] = (*(theCoupleTable->GetEnergyCutsVector(idxG4GammaCut)))[index]/MeV;
            fEnergyCut[1] = (*(theCoupleTable->GetEnergyCutsVector(idxG4ElectronCut)))[index]/MeV;
            fEnergyCut[2] = (*(theCoupleTable->GetEnergyCutsVector(idxG4PositronCut)))[index]/MeV;
            if (charge != 0.) {
                //G4cout << "\n Range cuts : \t gamma " << std::setw(8) << G4BestUnit(fRangeCut[0],"Length") << "\t e- " << std::setw(8) << G4BestUnit(fRangeCut[1],"Length");
                //G4cout << "\n Energy cuts : \t gamma " << std::setw(8) << G4BestUnit(fEnergyCut[0],"Energy") << "\t e- " << std::setw(8) << G4BestUnit(fEnergyCut[1],"Energy") << G4endl;
            }

            // get processList and extract EM processes (but not MultipleScattering)
            G4ProcessVector* plist = G4ParticleTable::GetParticleTable()->FindParticle(partName)->GetProcessManager()->GetProcessList();
            G4String procName;
            G4double cut;
            std::vector<G4String> emName;
            std::vector<G4double> enerCut;
            size_t length = plist->size();
            for (size_t j=0; j<length; j++) {
                procName = (*plist)[j]->GetProcessName();
                cut = fEnergyCut[1];
                if ((procName == "eBrem")||(procName == "muBrems")) ; cut = fEnergyCut[0];
                if (((*plist)[j]->GetProcessType() == fElectromagnetic) && (procName != "msc")) {
                    emName.push_back(procName);
                    enerCut.push_back(cut);
                }
            }

            // print list of processes
            if(nee == 0){

                // G4cout << nee << " then  processes names :                ";

                for ( G4int A = 0; A < emName.size()+2 ; A++  ) // 2 = one column for energies and one for total
                {
                    LatexText << "l";
                }
                LatexText << "} \\hline \n";
#if VERBOSE_USE
                G4cout << std::setw(17) << std::left << "Energy(MeV)" ;
#endif
                TableText << std::setw(17) << std::left << "Energy(MeV)" ;
                LatexText << "\\textbf{Energy(MeV)}           ";

                for (size_t j=0; j<emName.size();j++){
#if VERBOSE_USE
                    G4cout << std::setw(17) <<std::left << emName[j];
#endif
                    TableText << std::setw(17) <<std::left << emName[j];
                    LatexText << "& \\textbf{" << emName[j] << "}           ";
                }
#if VERBOSE_USE
                G4cout << std::setw(17) <<std::left <<"total"
                       << std::setw(17) <<std::left <<"Range Cut gamma"
                       << std::setw(17) <<std::left <<"Range Cut e-"
                       << std::setw(17) <<std::left <<"Range Cut e+"
                       << std::setw(17) <<std::left <<"Energy Cut gamma"
                       << std::setw(17) <<std::left <<"Energy Cut e-"
                       << std::setw(17) <<std::left <<"Energy Cut e+" << G4endl;
#endif

                TableText << std::setw(17) <<std::left <<"total"
                          << std::setw(17) <<std::left <<"Range Cut gamma"
                          << std::setw(17) <<std::left <<"Range Cut e-"
                          << std::setw(17) <<std::left <<"Range Cut e+"
                          << std::setw(17) <<std::left <<"Energy Cut gamma"
                          << std::setw(17) <<std::left <<"Energy Cut e-"
                          << std::setw(17) <<std::left <<"Energy Cut e+" << G4endl;

                //<< emName[0]<< "/"<< emName[1] <<G4endl;
                LatexText << "& \\textbf{total}           ";
                nee ++;
            }



            //get cross section per volume

            std::vector<G4double> sigma0;
            std::vector<G4double> sigma1;
            std::vector<G4double> sigma2;
            G4double Sig, SigtotComp = 0., Sigtot = 0.;

            for (size_t j=0; j<emName.size();j++) {
                Sig = emCal.ComputeCrossSectionPerVolume
                        (energy,particle,emName[j],material,enerCut[j]);
                SigtotComp += Sig;
                sigma0.push_back(Sig);
                Sig = emCal.GetCrossSectionPerVolume(energy,particle,emName[j],material);
                Sigtot += Sig;
                sigma1.push_back(Sig);
                sigma2.push_back(Sig/density);
            }
            sigma0.push_back(SigtotComp);
            sigma1.push_back(Sigtot);
            sigma2.push_back(Sigtot/density);
#if VERBOSE_USE
            G4cout << std::setw(17) <<std::left << energy;
#endif
            TableText << std::setw(17) <<std::left << energy;
            LatexText << "       \\\\\\hline\n";
            LatexText << energy << "      ";

            for (size_t j=0; j<sigma1.size();j++) {
                LatexText << "& "<< sigma1[j]*cm << "      ";
#if VERBOSE_USE
                G4cout << std::setw(17) <<std::left << sigma1[j]*cm;
#endif
                TableText << std::setw(17) <<std::left << sigma1[j]*cm;
            }
#if VERBOSE_USE
            G4cout << std::setw(17) <<std::left << fRangeCut[0]
                   << std::setw(17) <<std::left << fRangeCut[1]
                   << std::setw(17) <<std::left << fRangeCut[2]
                   << std::setw(17) <<std::left << fEnergyCut[0]
                   << std::setw(17) <<std::left << fEnergyCut[1]
                   << std::setw(17) <<std::left << fEnergyCut[2]
                      //G4cout << sigma1[0]/sigma1[1];
                   << G4endl;
#endif
            TableText << std::setw(17) <<std::left << fRangeCut[0]
                      << std::setw(17) <<std::left << fRangeCut[1]
                      << std::setw(17) <<std::left << fRangeCut[2]
                      << std::setw(17) <<std::left << fEnergyCut[0]
                      << std::setw(17) <<std::left << fEnergyCut[1]
                      << std::setw(17) <<std::left << fEnergyCut[2]
                         //G4cout << sigma1[0]/sigma1[1];
                      << G4endl;
        }

        LatexText << "\\\\ \\hline\n\\end{tabular} \n";
        LatexText << "\\label{CrossSectioPerVolumeInMaterial"<< matName<<"}\n";
        LatexText << "\\end{table}";
        LatexText << "\n\n\n\n\n";

        std::ofstream outfile(FileName , std::ios::app);
        if(outfile.is_open()){

            //G4cout << "\nCreating file " << FileName << G4endl ;
            outfile << TableText.str() << "\n" ;
            outfile << LatexText.str();
            outfile.close();
        }
        //G4cout << LatexText.str() << G4endl;
    }
}
