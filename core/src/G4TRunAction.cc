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
// *****************************************************************
// S. Guatelli. University of Wollongong, Australia.***
//
// Author: Tarik El Ghalbzouri,  Abdelmalek Essaâdi University,
// faculty of sciences Tetouane, morocco. email : telghalbzouri@uae.ac.ma
//
// This application is based on code developed by :
// G. Guerrieri, University of Genova, Italy .
// S. Guatelli. University of Wollongong, Australia.
//

#include "G4TRunAction.hh"
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

//#include "G4TOutputText.hh"

#include <stdio.h>
#include <string.h>

#include "G4Threading.hh"

//#include "G4SDManager.hh"
//#include "G4TSD.hh"
//#include "G4AutoLock.hh"
//#include "G4TRun.hh"

#include "G4EmCalculator.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTable.hh"
#include "G4Timer.hh"


extern G4String GeometryFileType;
extern G4String GenerateVoxelsResuls;

extern  G4String* CopyNumberRegionNameMap;
extern  G4float* CopyNumberXPos;
extern  G4float* CopyNumberYPos;
extern  G4float* CopyNumberZPos;
extern  G4float* CopyNumberMassSize;

// Source Data

extern G4String ParticleName;
extern G4String GeometrySymbol;
extern G4String SourceType;
extern G4String EnergyDistribution;
extern G4String UseGeneratedData ;
extern G4String ResultDirectoryPath ;
extern G4String DataDirectoryPath;

extern G4int NumberOfEventInBatch;
extern G4int NumberOfBatch;
extern G4String DataFilesExtension;

extern std::vector<G4String> NewRankSourceParticlesNamesValues;
extern std::vector<G4String> NewRankSourceRegionsNamesValues;
extern std::vector<G4double> NewRankSourceEnergiesValues;
extern std::vector<G4String> NewRankSourceMomDirsValues;


extern std::map<unsigned int,G4double* >      EnergyListForCriticality;
extern std::map<unsigned int,G4ThreeVector* > PositionsListForCriticality;
extern std::map<unsigned int,G4ParticleMomentum* > MomDirecsListForCriticality;

extern std::map<G4int,std::map<G4int,std::vector<G4double>>> FissionCapturesOfThreadsRanks;
extern std::map<G4int,std::map<G4String,std::map<G4String,G4int>>> OpticalPhotonInteractionRate;

extern std::vector<G4double> KeffectiveInEachBatch;
extern std::map<G4int,std::map<G4int,G4bool>> TerminatedThreadBatch;

#ifdef G4MULTITHREADED
G4ThreadLocal std::ofstream G4TRunAction::CriticalityFile;
G4ThreadLocal std::map<G4int,std::vector<G4double>>      G4TRunAction::NewNeutronFissionEnergyList;
G4ThreadLocal std::map<G4int,std::vector<G4ThreeVector>> G4TRunAction::NewNeutronFissionPositionsList;
G4ThreadLocal std::map<G4int,std::vector<G4ParticleMomentum>> G4TRunAction::NewNeutronFissionMomDirecsList;
G4ThreadLocal unsigned long long int*    G4TRunAction::FluxParticlePosX   ;
G4ThreadLocal unsigned long long int*    G4TRunAction::FluxParticlePosY   ;
G4ThreadLocal unsigned long long int*    G4TRunAction::FluxParticlePosZ   ;
G4ThreadLocal unsigned long long *       G4TRunAction::FluxParticleEnergy ;
G4ThreadLocal std::map<G4int,std::map<G4String,G4int>> G4TRunAction::ParticleProduction;
G4ThreadLocal std::map<G4int,std::map<unsigned int,G4double>> G4TRunAction::VoxelsFluence ;
G4ThreadLocal std::map<G4int,std::map<unsigned int,G4double>> G4TRunAction::VoxelsED_Total ;
G4ThreadLocal std::map<G4int,std::map<unsigned int,G4double>> G4TRunAction::VoxelsED2_Total ;
G4ThreadLocal std::map<G4int,std::map<unsigned int,unsigned long long int>> G4TRunAction::VoxelsNOfValues ;
G4ThreadLocal std::map<G4int,std::map<G4String,unsigned long long int>> G4TRunAction::NOfValues;
G4ThreadLocal std::map<G4int,std::map<G4String,G4double>> G4TRunAction::ED_Total ;
G4ThreadLocal std::map<G4int,std::map<G4String,G4double>> G4TRunAction::Fluence ;
G4ThreadLocal std::map<G4int,std::map<G4String,G4double>> G4TRunAction::ED2_Total ;
G4ThreadLocal G4int G4TRunAction::rank, G4TRunAction::thread, G4TRunAction::DataID, G4TRunAction::CurrentBatch, G4TRunAction::NumberOfRanksThreads, G4TRunAction::TotalEventNumber;
//G4ThreadLocal G4double G4TRunAction::Keffective;
//G4ThreadLocal G4int G4TRunAction::NumberOfFissionNeutrons, G4TRunAction::NumberOfFission, G4TRunAction::NumberOfCapture,
G4ThreadLocal G4double G4TRunAction::EnergyEmittedPerThread, G4TRunAction::ParticleSourceEnergy, G4TRunAction::ExecutionTimeInMin, G4TRunAction::OneEventExecutionTimeInMs;
G4ThreadLocal std::chrono::steady_clock::time_point G4TRunAction::start, G4TRunAction::end ;
#endif

//G4ThreadLocal unsigned int G4TRunAction::StepCN;
//G4ThreadLocal G4double G4TRunAction::StepEnergy;
//G4ThreadLocal G4String G4TRunAction::StepRegion;

//namespace
//{
//G4Mutex	mutex = G4MUTEX_INITIALIZER;
//}

G4TRunAction::G4TRunAction()
{
    // get the name of the working directory
    //G4String appBuildDir(getenv("PWD"));
    //ResultDirectoryPath = appBuildDir+"/Results";
}
G4TRunAction::~G4TRunAction(){
    G4MUTEXDESTROY(mutex);
    //G4cout << "\n\n\n\n\n\n from function : " << "G4TRunAction::~G4TRunAction()"<< G4endl;
}

//G4Run* G4TRunAction::GenerateRun(){ return new G4TRun;}

void G4TRunAction::WriteMacroscopicCrossSection(){

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

// called in the begining of each run for each thread and for the master for seq and master in the MPI mode
void G4TRunAction::BeginOfRunAction(const G4Run* aRun) {

    //std::cout << " from function : " << "\n\n\n\n\n\n G4TRunAction::BeginOfRunAction()"<< std::endl;
    //out = new G4TOutputText();

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

    //std::cout <<  " GeometryFileType " << GeometryFileType <<  " GenerateVoxelsResuls " << GenerateVoxelsResuls << std::endl ;

    if(G4Threading::IsMultithreadedApplication()){ // normal multiThreaded mode

        ExecutionMode = "MT";
        NumberOfRanksThreads = G4Threading::GetNumberOfRunningWorkerThreads() ;
        TotalEventNumber = TotalEventNumber/G4Threading::GetNumberOfRunningWorkerThreads();
        //TotalEventNumber = aRun->GetNumberOfEvent();

        if(G4Threading::IsWorkerThread()){
            std::cout <<  "\n " << Time << " ========= MT mode : "<<__FUNCTION__<<  " " << GeometryFileType  << " " << " GenerateVoxelsResuls:" << GenerateVoxelsResuls << " from Worker Thread " << G4Threading::G4GetThreadId() << "/" <<G4Threading::GetNumberOfRunningWorkerThreads() << ". Start of " << TotalEventNumber << " events simulation loop using " << GeometrySymbol <<" geometry -Source Data: " << NewRankSourceParticlesNamesValues[DataID] << ", " << SourceType << ", " << NewRankSourceRegionsNamesValues[DataID] << ", "  << EnergyDistribution << " " << NewRankSourceEnergiesValues[DataID] << ", " << NewRankSourceMomDirsValues[DataID] << std::endl ;

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
        Fluence[DataID][OrgansNameVector[nn]]=0.;
        ED_Total[DataID][OrgansNameVector[nn]]=0.;
        ED2_Total[DataID][OrgansNameVector[nn]]=0.;
    }
    ED_Total[DataID]["World"]=0.;

    useTime(0); // "0" for setting the begining time befor a loop, "1" it get the end time and calculate consumed time by a loop and show the difference, the loop should not have a setting or getting time inside it

    G4RunManager::GetRunManager()->SetPrintProgress(TotalEventNumber*0.1);

    CurrentBatch = 0;
    NumberOfEventInBatch = 300;
    NumberOfBatch = 10;

    //NumberOfFission = 0;
    //NumberOfFissionNeutrons = 0;
    //NumberOfCapture = 0;
    //Keffective = 0.;
    for (G4int iLV = 0; iLV < NumberOfBatch; iLV++ ) {
        if(FissionCapturesOfThreadsRanks[DataID][iLV].size() == 3){
        }else{
            FissionCapturesOfThreadsRanks[DataID][iLV].push_back(0.);
            FissionCapturesOfThreadsRanks[DataID][iLV].push_back(0.);
            FissionCapturesOfThreadsRanks[DataID][iLV].push_back(0.);
        }
        if(KeffectiveInEachBatch.size() == NumberOfBatch){
        }else{
            KeffectiveInEachBatch.push_back(0.);
        }
    }

    //std::cout << "\n\n\n\n\n\n 11111111111111 \n\n\n\n\n\n "<< std::endl;
}

// called in the end of a run for each thread and for the master for seq and master in the MPI mode
// it call analysisMan->save(); G4TRunAction::CreateResultFileFor_Seq_MT()
void G4TRunAction::EndOfRunAction(const G4Run* aRun)
{

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
        }else{
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

    ShowOpticalPhotonInteractionData();


    // because MyGeometry can be not voxelized or stylized and voxelized in the same time which is not one
    // of the voxelized geometries VoxIDs, VOXEL, and DICOM.
    //if(GeometryFileType == "VoxIDs" || GeometryFileType == "VOXEL" || GeometryFileType == "DICOM"){
        if(GenerateVoxelsResuls == "yes"){
            CreateThreadVoxelsResultsFiles();
        }else{
            CreateThreadRegionResultFile();
        }

    //}else{
    //    CreateThreadRegionResultFile();
    //}






    /*
    if(SourceType == "Voxels"){
        if(GenerateVoxelsResuls == "yes"){
            CreateThreadVoxelsResultsFiles();
        }
    }
*/
}

// for step and event accuracy calculation
void G4TRunAction::CreateThreadVoxelsResultsFiles(){

    //G4cout << "\n\n\n\n\n" << __FUNCTION__ << G4endl ;

    std::ostringstream fname;

    //G4cout << __FUNCTION__ << G4endl ;

    const G4TVolumeConstruction* VolumeConstruction1 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fname << ResultDirectoryPath <<"/AE@Voxel@for@Rank@"<<rank<<"@Thread@"<< thread<<"@" << GeometrySymbol <<"@" << NewRankSourceParticlesNamesValues[DataID] << "@" << NewRankSourceRegionsNamesValues[DataID] << "@" << NewRankSourceEnergiesValues[DataID];
    std::cout << " Creating voxels result file of worker thread or rank " << fname.str().c_str()<< std::endl;

    std::ofstream file(fname.str().c_str(), std::ios_base::binary);

    //G4cout << "Voxels Scored " << VoxelsED_Total[DataID].size() << G4endl ;
    file << VoxelsED_Total[DataID].size()
         << " " << VoxXNumber
         << " " << VoxYNumber
         << " " << VoxZNumber
         << " " << VoxXHalfSize
         << " " << VoxYHalfSize
         << " " << VoxZHalfSize
         << "\n";

    //G4cout << "VoxelsED_Total[DataID].size() " << VoxelsED_Total[DataID].size() << G4endl ;

    for ( auto it = VoxelsED_Total[DataID].begin(); it != VoxelsED_Total[DataID].end(); ++it  )
    {
        file << it->first << " "
             << CopyNumberXPos[it->first] << " "
             << CopyNumberYPos[it->first] << " "
             << CopyNumberZPos[it->first] << " "
             << CopyNumberRegionNameMap[it->first] << " "
             << VoxelsED_Total[DataID][it->first] << " "
             << VoxelsED2_Total[DataID][it->first] << " "
             << VoxelsNOfValues[DataID][it->first] << " "
             << CopyNumberMassSize[it->first] << " "
             << VoxelsFluence[DataID][it->first]/(8*VoxXNumber*VoxXHalfSize*VoxYNumber*VoxYHalfSize*VoxZNumber*VoxZHalfSize) << " "
             << "\n" ;
    }
    //G4cout << "7" << G4endl ;


    file << "TotalEventNumber" << " " << TotalEventNumber << "\n" ;

    //file << "EnergyEmittedPerThread" << " " << EnergyEmittedPerThread << "\n" ;
    //file << "ExecutionTimeInMin" << " " << ExecutionTimeInMin << "\n" ;
    //OneEventExecutionTimeInMs = (double)(ExecutionTimeInMin/TotalEventNumber)*(0.001*0.001/60.);
    //file << "OneEventExecutionTimeInMs" << " " << OneEventExecutionTimeInMs << "\n" ;
    file.close();

    //delete[] VCNM ;delete[] XP; delete[] YP; delete[] ZP;
    //for ( auto it = VoxelsED_Total[DataID].begin(); it != VoxelsED_Total[DataID].end(); ++it  )
    //{VoxelsED_Total[DataID].erase(it->first); VoxelsED2_Total[DataID].erase(it->first); VoxelsNOfValues[DataID].erase(it->first);}
    //VoxelsED_Total[DataID].clear(); VoxelsED2_Total[DataID].clear(); VoxelsNOfValues[DataID].clear();

}

// for step and event accuracy calculation
void G4TRunAction::CreateThreadRegionResultFile(){

    //G4cout << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ void G4TRunAction::CreateResultFileFor_Seq_MT()  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << G4endl ;
    //std::ostringstream OutStringWords; OutStringWords << "  ";
    //out->showTextMessage(OutStringWords.str() , 1);

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
             << std::setw(20) << std::left << TConstruction2->GetOrganNameVolumeMap()[OrgansNameVector[jk]]<< " "
             << std::setw(20) << std::left << TConstruction2->GetOrganNameDensityMap()[OrgansNameVector[jk]]<< " ";
        if(TConstruction2->GetOrganNameVolumeMap()[OrgansNameVector[jk]] == 0. || __isinf(TConstruction2->GetOrganNameVolumeMap()[OrgansNameVector[jk]]) || __isnan(TConstruction2->GetOrganNameVolumeMap()[OrgansNameVector[jk]])){
            file << std::setw(20) << std::left << 0 << " ";
        }else{
            file << std::setw(20) << std::left << Fluence[DataID][OrgansNameVector[jk]]/TConstruction2->GetOrganNameVolumeMap()[OrgansNameVector[jk]] << " ";
        }

        file << "\n";
        //G4cout  << OrgansNameVector[jk] << "  " << ED_Total[DataID][OrgansNameVector[jk]] << "  " << G4endl;
    }

    SZ++;

    file << "\n";
    file << std::setw(SZ) << std::left << "TotalEventNumber " << TotalEventNumber <<  "\n" ;
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

void G4TRunAction::TestResults(){

#if VERBOSE_USE
    G4cout << "\n========= Testing Results Informations ========= " << G4endl ;

    G4cout << "Deposited Energy (MeV)        = " << ED_Total[DataID]["Liver"]  << G4endl ;
    G4cout << "Number Of Events              = " << TotalEventNumber  << G4endl ;
    G4cout << "Number Of Values              = " << NOfValues[DataID]["Liver"]  << G4endl ;
    G4cout << "EnergyEmittedPerThread (MeV)  = " << EnergyEmittedPerThread  << G4endl ;
    G4cout << "SAF Value                     = " << ED_Total[DataID]["Liver"]/(NOfValues[DataID]["Liver"]*EnergyEmittedPerThread)  << G4endl ;
#endif

}

// called from the endOfEvent()
void G4TRunAction::CreateSimulationDataFile(){

    // this methods is called from master run, where the thread or rank ID is -1 , then its not in NewRankSourceEnergiesValues which befin from 0
    G4String distriName = EnergyDistribution;
    if(distriName == "Mono"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[0] ;
    }
    else if(distriName == "Uniform"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[0] ;
    }
    else if(distriName == "Rayleigh"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[0];
    }
    else if(distriName == "Gauss"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[0];
    }
    else if(distriName == "Spectrum"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[0];
    }
    else if(distriName == "File"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[0];
    }
    else if(distriName == "RadioNuclide"){
        ParticleSourceEnergy = NewRankSourceEnergiesValues[0];
    }

    std::ostringstream fname1;

    const G4TVolumeConstruction* TConstruction2 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    const G4TUserPhysicsList* Phy1 = static_cast<const G4TUserPhysicsList*> (G4RunManager::GetRunManager()->GetUserPhysicsList());

    fname1 << ResultDirectoryPath <<"/SimData" ;
    std::ofstream file1(fname1.str().c_str(), std::ios_base::binary);

    file1 << "\n>> Simulation Data" <<  "\n\n" ;

    if(file1.is_open()){

        // given from VolumeConstruction class

        file1 << "GeometrySymbol                   " << GeometrySymbol <<  "\n\n" ;

        file1 << "CutsDistance                     " << CutsDistance <<  "\n" ;
        file1 << "CutsEnergy                       " << CutsEnergy <<  "\n" ;
        file1 << "ParticleName                     " << NewRankSourceParticlesNamesValues[0] <<  "\n" ;
        file1 << "SourceType                       " << TConstruction2->getSourceType() <<  "\n" ;
        file1 << "SourceRegionName                 " << NewRankSourceRegionsNamesValues[0] <<  "\n" ;
        file1 << "EnergyDistribution               " << TConstruction2->getEnergyDistribution() <<  "\n" ;
        if(EnergyDistribution == "Mono"){
            file1 << "MonoEnergy                       " << TConstruction2->getMonoEnergy() <<  "\n" ;
        }
        else if(EnergyDistribution == "Uniform"){
            file1 << "UniformEmin                      " << TConstruction2->getUniformEmin() <<  "\n" ;
            file1 << "UniformEmax                      " << TConstruction2->getUniformEmax() <<  "\n" ;
        }
        else if(EnergyDistribution == "Rayleigh"){
            file1 << "RayleighEmax                     " << TConstruction2->getRayleighEmax() <<  "\n" ;
        }
        else if(EnergyDistribution == "Spectrum"){
            file1 << "SpectrumMaxEnergy                " << TConstruction2->getSpectrumMaxEnergy() <<  "\n" ;
        }else if(EnergyDistribution == "File"){
            file1 << "FileEnergyCharacterizer          " << TConstruction2->getFileEnergyCharacterizer() <<  "\n" ;
        }
        else if(EnergyDistribution == "RadioNuclide"){
            file1 << "RadioNuclideMaxEnergy            " << TConstruction2->getRadioNuclideMaxEnergy() <<  "\n" ;
        }
        else if(EnergyDistribution == "Gauss"){
            file1 << "GaussSDev                        " << TConstruction2->getGaussSDev() <<  "\n" ;
            file1 << "GaussMean                        " << TConstruction2->getGaussMean() <<  "\n" ;
        }
        file1 << "ParticleSourceEnergy             " << ParticleSourceEnergy <<  "\n" ;
        file1 << "MomDirDistribution               " << NewRankSourceMomDirsValues[0] <<  "\n\n" ;

        file1 << "ResultDirectoryPath              " << ResultDirectoryPath <<  "\n\n" ;

        file1 << "GraphsData                       " << TConstruction2->getGraphs_Data() <<  "\n" ;
        file1 << "CompareType                      " << TConstruction2->getcompare_type() <<  "\n" ;
        file1 << "RefName                          " << TConstruction2->getref_Name() <<  "\n" ;
        file1 << "RefFilePath                      " << TConstruction2->getref_File_Path() <<  "\n" ;

        file1 << "GenerateRegionsVariableGraph     " << TConstruction2->getGenerateRegionsVariableGraph() <<  "\n" ;
        file1 << "RegionVariableName               " << TConstruction2->getRegionVariableName() <<  "\n" ;
        file1 << "GenerateRelativeSDevGraph        " << TConstruction2->getGenerateRelativeSDevGraph() <<  "\n" ;
        file1 << "GenerateRelativeErrGraph         " << TConstruction2->getGenerateRelativeErrGraph() << " " << TConstruction2->getDifferenceMethod() << "\n" ;
        file1 << "GenerateCrossSectionGraph        " << GenerateCrossSectionGraph <<  "\n" ;
        file1 << "PositionDataFile                 " << TConstruction2->getPositionDataFile() <<  "\n" ;
        file1 << "EnergyDataFile                   " << TConstruction2->getEnergyDataFile() <<  "\n" ;
        file1 << "MomDirDataFile                   " << TConstruction2->getMomDirDataFile() <<  "\n" ;
        file1 << "EventsDataHistograms             " << TConstruction2->getEventsDataHistograms() <<  "\n\n" ;

        file1 << "QuantitiesUnits                  " << TConstruction2->getQuantitiesUnits() <<  "\n" ;
        file1 << "RadiationFactors                 " << TConstruction2->getRadiationFactors() <<  "\n" ;
        file1 << "TissueFactors                    " << TConstruction2->getTissueFactors() <<  "\n" ;
        for (int var = 0; var < TConstruction2->getRadioTracerData().size(); ++var) {
            file1 << "RadioTracerData                  " << TConstruction2->getRadioTracerData()[var] << "\n" ;
        }
        for (int var = 0; var < TConstruction2->getRadioTracerBiokinetic().size(); ++var) {
            file1 << "RadioTracerBiokinetic            " << TConstruction2->getRadioTracerBiokinetic()[var] << "\n" ;
        }
        file1 << "\n" ;

        file1 << "UseLogE                          " << TConstruction2->getUseLogE() <<  "\n" ;
        file1 << "UseLogVariable                   " << TConstruction2->getUseLogVariable() <<  "\n" ;
        file1 << "UseGridXY                        " << TConstruction2->getUseGridXY() <<  "\n" ;
        file1 << "PrintTitle                       " << TConstruction2->getPrintTitle() <<  "\n" ;
        file1 << "LegendPos                        " << TConstruction2->getLegendPos() <<  "\n" ;
        file1 << "LegendXWidth                     " << TConstruction2->getLegendXWidth() <<  "\n" ;
        file1 << "LegendYHeight                    " << TConstruction2->getLegendYHeight() <<  "\n" ;
        file1 << "AddErrorBarInGraphs              " << TConstruction2->getAddErrorBarInGraphs() <<  "\n" ;
        file1 << "GraphsExt                        " << TConstruction2->getgraphs_Ext() <<  "\n\n" ;

        if(SourceType == "Voxels"){
            file1 << "VoxXNumber                       " << TConstruction2->getVoxXNumber() <<  "\n" ;
            file1 << "VoxYNumber                       " << TConstruction2->getVoxYNumber() <<  "\n" ;
            file1 << "VoxZNumber                       " << TConstruction2->getVoxZNumber() <<  "\n" ;
            file1 << "VoxXHalfSize                     " << TConstruction2->getVoxXHalfSize() <<  "\n" ;
            file1 << "VoxYHalfSize                     " << TConstruction2->getVoxYHalfSize() <<  "\n" ;
            file1 << "VoxZHalfSize                     " << TConstruction2->getVoxZHalfSize() <<  "\n" ;
            file1 << "GenerateVoxelsResuls             " << TConstruction2->getGenerateVoxelsResuls() <<  "\n" ;
            file1 << "DoseProfilQuantity               " << TConstruction2->getDoseProfilQuantity() <<  "\n" ;
            file1 << "BeamAxis                         " << TConstruction2->getBeamAxis() <<  "\n" ;
            file1 << "SliceFor2DGraph                  " << TConstruction2->getSliceFor2DGraph() <<  "\n" ;
            file1 << "SliceID                          " << TConstruction2->getSliceID() <<  "\n" ;
        }
        file1 << "QuantitiesToScore                " << TConstruction2->getvariable_To_Score() <<  "\n" ;
        file1 << "RegionsNamesToScore              " << TConstruction2->getorgans_to_score() <<  "\n\n" ;
        //file1 << "AccuracyCalculationLevel         " << AccuracyCalculationLevel <<  "\n" ;

        file1 << "EventNumberForRankThread         " << TotalEventNumber <<  "\n" ;
        file1 << "ExecutionMode                    " << ExecutionMode <<  "\n" ;
        file1 << "OneOrMultiSimulations            " << OneOrMultiSimulations <<  "\n" ;
        file1 << TConstruction2->getRanksSourceData();
        file1 << "NumberOfRanksThreads             " << NumberOfRanksThreads <<  "\n\n" ;
        file1.close();
    }
    else{

        G4cout << "cannot open the file " << fname1.str() << G4endl ;
    }
}

void G4TRunAction::FillRegionStepHits(G4String StepRegion, G4double StepEnergy){

    //G4cout << N << " " << E << " "<< E*E<< G4endl;

    ED_Total[DataID][StepRegion] += StepEnergy ;
    ED2_Total[DataID][StepRegion] += StepEnergy*StepEnergy;
    NOfValues[DataID][StepRegion]++;

    //G4cout << DataID << " " << StepRegion << " " << ED_Total[DataID][StepRegion] << " " << ED2_Total[DataID][StepRegion] << G4endl;
}
void G4TRunAction::FillRegionLenghts(G4String StepRegion, G4double StepEnergy){

    //G4cout << N << " " << E << " "<< E*E<< G4endl;

    Fluence[DataID][StepRegion] += StepEnergy ;

    //G4cout << DataID << " " << StepRegion << " " << Fluence[DataID][StepRegion] << G4endl;
}
void G4TRunAction::FillVoxelStepHits(unsigned int StepCN, G4double StepEnergy){

    //G4cout << DataID << StepCN << " " << StepEnergy << G4endl;

    VoxelsED_Total[DataID][StepCN] += StepEnergy ;
    VoxelsED2_Total[DataID][StepCN] += StepEnergy*StepEnergy;
    VoxelsNOfValues[DataID][StepCN]++;

}
void G4TRunAction::FillVoxelLenghts(unsigned int StepCN, G4double StepEnergy){

    //G4cout << DataID << StepCN << " " << StepEnergy << G4endl;

    VoxelsFluence[DataID][StepCN] += StepEnergy ;

}
// for Voxelized geometry , called for each event
// Step Level accuracy calculation
void G4TRunAction::pushVoxelEnergy(std::map<unsigned int,G4double>EventMap){

    for ( auto it = EventMap.begin(); it != EventMap.end(); ++it  )
    {
        VoxelsED_Total[DataID][it->first] += it->second ;
        ED_Total[DataID][CopyNumberRegionNameMap[it->first]] += it->second;
    }
}
void G4TRunAction::pushVoxelEnergy2(std::map<unsigned int,G4double>EventMap){

    for ( auto it = EventMap.begin(); it != EventMap.end(); ++it  )
    {
        VoxelsED2_Total[DataID][it->first] += it->second*it->second ;
        ED2_Total[DataID][CopyNumberRegionNameMap[it->first]] += it->second;
    }
}
void G4TRunAction::pushVoxelHitsNumber(std::map<unsigned int,unsigned int>EventMap){

    for ( auto it = EventMap.begin(); it != EventMap.end(); ++it  )
    {
        VoxelsNOfValues[DataID][it->first] += it->second;
        NOfValues[DataID][CopyNumberRegionNameMap[it->first]] += it->second;
    }
}

// for GDML geometry
// Step Level accuracy calculation
void G4TRunAction::pushVolumeEnergy(std::map<G4String,G4double>StepsMap){
    for ( auto it = StepsMap.begin(); it != StepsMap.end(); ++it  ){

        //G4cout << it->first  << " " << it->second <<  G4endl ;
        ED_Total[DataID][it->first] += it->second ;
    }
}
void G4TRunAction::pushVolumeEnergy2(std::map<G4String,G4double>StepsMap){
    for ( auto it = StepsMap.begin(); it != StepsMap.end(); ++it  ){

        //G4cout << it->first  << " " << it->second <<  G4endl ;
        ED2_Total[DataID][it->first] += it->second ;
    }
}
void G4TRunAction::pushVolumeHitsNumber(std::map<G4String, unsigned int> StepsMap){
    for ( auto it = StepsMap.begin(); it != StepsMap.end(); ++it  ){

        //G4cout << it->first  << " " << it->second <<  G4endl ;
        NOfValues[DataID][it->first] += it->second ;
    }
}


void G4TRunAction::CountFission() {
    //G4cout << "Fission" << G4endl;
    //NumberOfFission++;
    FissionCapturesOfThreadsRanks[DataID][CurrentBatch][2]++ ;
}
void G4TRunAction::CountAbsorption() {
    //G4cout << "Capture" << G4endl;
    //NumberOfCapture++;
    FissionCapturesOfThreadsRanks[DataID][CurrentBatch][1]++;
}
void G4TRunAction::SetPoxEneForBatch(G4double X,G4double Y,G4double Z,G4double MX,G4double MY,G4double MZ,G4double E) {
    NewNeutronFissionEnergyList[DataID].push_back(E);
    NewNeutronFissionPositionsList[DataID].push_back(G4ThreeVector(X,Y,Z));
    NewNeutronFissionMomDirecsList[DataID].push_back(G4ParticleMomentum(MX,MY,MZ));
    //NumberOfFissionNeutrons++;
    FissionCapturesOfThreadsRanks[DataID][CurrentBatch][0]++;
    //G4cout << " DataID " << DataID  << " E " << E<< " X "<< X<< " Y "<< Y<< " Z "<< Z<< " MX "<< MX<< " MY "<< MY<< " MZ "<< MZ<< G4endl;

    //NewNeutronFissionMomDirecsList[DataID].push_back(G4ThreeVector(X,Y,Z));
}
void G4TRunAction::CountParticleProductionByNeutron(G4String n){
    ParticleProduction[DataID][n]++;
}

void G4TRunAction::IntializeBatchData() {

    TerminatedThreadBatch[DataID][CurrentBatch] = true;

    //G4AutoLock l(&mutex);

    //G4cout << " DataID " << DataID << " CurrentBatch "<< CurrentBatch << " NumberOfEventInBatch "<< NumberOfEventInBatch
    //       //<< " CurrentBatch*NumberOfEventInBatch+nn "<< CurrentBatch*NumberOfEventInBatch+nn
    //       << " NewNeutronFissionEnergyList[DataID].size() "<< NewNeutronFissionEnergyList[DataID].size()
    //       << " NewNeutronFissionPositionsList[DataID].size() "<< NewNeutronFissionPositionsList[DataID].size()
    //       << " NewNeutronFissionMomDirecsList[DataID].size() "<< NewNeutronFissionMomDirecsList[DataID].size()
    //       << G4endl;


    //G4cout << " DataID " << DataID << " CurrentBatch "<< CurrentBatch << " FissionCapturesOfThreadsRanks[DataID][CurrentBatch][0] "<< FissionCapturesOfThreadsRanks[DataID][CurrentBatch][0] << G4endl;



    //G4cout << " 0 --------- " << G4endl;


    G4double NumCap = 0., NumFiss = 0.,NumFissN = 0.;

    //for ( auto it = FissionCapturesOfThreadsRanks.begin(); it != FissionCapturesOfThreadsRanks.end(); ++it  ){
    //    if(it->first == -1){ continue;}
    //    //G4cout << "***************** DataID = " << it->first << " NumberOfFissionN " << it->second[0] << " NumCaptures " << it->second[1]  << " NumFission " << it->second[2] << G4endl;
    //    NumFissN += it->second[0];
    //    NumCap   += it->second[1];
    //    NumFiss  += it->second[2];
    //}

    CurrentBatch++;// the first batch is simulated from neutron source chosen by user

    //1- calculate batch Keff, should be the first in this function
    for (G4int  BatchInc = 0;  BatchInc < NumberOfBatch;  BatchInc++ ) {
        //G4cout << "AAA NumberOfBatch = " << NumberOfBatch << "  BatchInc " <<  BatchInc  << " KeffectiveInEachBatch[ BatchInc] " << KeffectiveInEachBatch[ BatchInc] << G4endl;

        NumFissN= 0;
        NumCap  = 0;
        NumFiss = 0;

        bool isready = true;
        if(KeffectiveInEachBatch[ BatchInc] <= 0){
            for ( auto it = FissionCapturesOfThreadsRanks.begin(); it != FissionCapturesOfThreadsRanks.end(); ++it  ){

                if(it->first == -1){ continue;}

                //G4cout << "BBB For Thread = " << it->first << " For BatchInc = " << BatchInc << " IsTerminated? "<< TerminatedThreadBatch[it->first][BatchInc] << " NumberOfFissionN " << it->second[ BatchInc][0] << " NumCaptures " << it->second[ BatchInc][1] << " NumFission " << it->second[ BatchInc][2] << G4endl;

                if(TerminatedThreadBatch[it->first][BatchInc] == false){
                    isready = false;
                    break;
                }
                NumFissN += it->second[ BatchInc][0];
                NumCap   += it->second[ BatchInc][1];
                NumFiss  += it->second[ BatchInc][2];
            }
            //G4cout << "CCC Total for :     For BatchInc = " << BatchInc << " NumberOfFissionN " << NumFissN << " NumCaptures " << NumCap  << " NumFission " << NumFiss << G4endl;

        }else{
            continue;
        }

        if(isready){

            //G4cout << "***************** Total for : NumberOfFissionN = " << NumFissN << " NumCaptures " << NumCap  << " NumFission " << NumFiss << G4endl;

            //G4double keff = static_cast<G4double>(RunAction->getNumberOfFissionNeutrons()/(RunAction->getNumberOfFission()+4.94066e-324)) * static_cast<G4double>(RunAction->getNumberOfFission()/ (RunAction->getNumberOfCapture() +4.94066e-324)); // Avoid division by zero
            G4double keff = static_cast<G4double>(NumCap/(NumFiss+4.94066e-324)) * static_cast<G4double>(NumFiss/ (NumFissN +4.94066e-324)); // Avoid division by zero

            if( BatchInc == 0){
                G4cout << "***************** Batch " << BatchInc+1<<"/"<< NumberOfBatch << " , "<< NumberOfEventInBatch << "Event/Batch" << " , " << NumCap << " Captures, " << NumFiss << " Fissions, " << NumFissN << " FissionNeutrons, "<< " k-eff = " << keff << G4endl;
            }else if( BatchInc > 0){
                G4double relative_error = abs((keff - KeffectiveInEachBatch[ BatchInc-1]) / keff);
                G4cout << "***************** Batch " << BatchInc+1<<"/"<< NumberOfBatch << " , "<< NumberOfEventInBatch << "Event/Batch" << " , " << NumCap << " Captures, " << NumFiss << " Fissions, " << NumFissN << " FissionNeutrons, "<< " k-eff = " << keff << " +/- " << relative_error << G4endl;
            }

            KeffectiveInEachBatch[ BatchInc] = keff;
        }else {
            break;
        }
    }

    bool WriteData = true;
    //2- write the criticality calculation data in the last thread and batch
    if(CurrentBatch == NumberOfBatch && CurrentBatch*NumberOfEventInBatch >= TotalEventNumber){
        if(KeffectiveInEachBatch[CurrentBatch-1] != 0){
            for (G4int  nn = 0;  nn < KeffectiveInEachBatch.size();  nn++ ) {
                if(KeffectiveInEachBatch[nn] == 0.){
                    WriteData = false;
                }
            }
        }else{
            WriteData = false;
        }

        //G4cout << " DataID " << DataID << " CurrentBatch "<< CurrentBatch << " WriteData " << WriteData << G4endl;

        //G4cout << " Last Call --------- " << G4endl;

        if(WriteData == true){

            //G4cout << " 1 --------- " << G4endl;

            std::ostringstream c ;
            c << DataDirectoryPath  << "/CriticalityDataFile_"<< GeometrySymbol <<"_" <</*DataID*/0<< ".txt";//DataFilesExtension;
            CriticalityFile.open(c.str().c_str() , std::ios_base::binary); // , std::ios_base::out | std::ios_base::binary

            CriticalityFile << "NumOfBatchs "    << NumberOfBatch << "\n";
            CriticalityFile << "NumOfEventsPerBatch "    << NumberOfEventInBatch<< "\n";
            CriticalityFile << "NumOfFissionNeutrons "    <<     FissionCapturesOfThreadsRanks[DataID][CurrentBatch-1][0] << "\n";
            //CriticalityFile << "NumOfFissions "  << NumberOfFission<< "\n";
            //CriticalityFile << "NumOfCaptures "  << NumberOfCapture<< "\n";
            CriticalityFile << "k-eff "    << KeffectiveInEachBatch[CurrentBatch-1] << "\n";
            //G4cout << " 2 --------- " << G4endl;

            //std::cout << "Secondary particle and the production number"  << std::endl ;
            for ( auto it = ParticleProduction[DataID].begin(); it != ParticleProduction[DataID].end(); ++it  ){
                if(it->first == "nCapture"){
                    CriticalityFile << "NumOf"<< it->first << " " << FissionCapturesOfThreadsRanks[DataID][CurrentBatch-1][1]  << "\n" ;
                }else if (it->first == "nFission"){
                    CriticalityFile << "NumOf"<< it->first << " " << FissionCapturesOfThreadsRanks[DataID][CurrentBatch-1][2]  << "\n" ;
                }
                else{
                    CriticalityFile << "NumOf"<< it->first << " " << it->second  << "\n" ;
                }
            }
            //G4cout << " 3 --------- " << G4endl;

            G4int max = NewNeutronFissionEnergyList[DataID].size();
            if(NewNeutronFissionEnergyList[DataID].size() < NumberOfEventInBatch ){
                max = NewNeutronFissionEnergyList[DataID].size() ;
            }else{
                max = NumberOfEventInBatch ;
            }

            CriticalityFile << "SimulatedNeutronSourceData ("<< max <<") E X Y Z MX MY MZ\n";

            //for(G4int nn = NumberOfEventInBatch ;nn < TotalEventNumber ;nn++){ // to discard the first batch which is the user source configuration and not the reel source in the reactor code
            for(G4int nn = 0 ;nn < max ;nn++){

                CriticalityFile << NewNeutronFissionEnergyList   [DataID][nn]       << " "
                                << NewNeutronFissionPositionsList[DataID][nn].getX()<< " "
                                << NewNeutronFissionPositionsList[DataID][nn].getY()<< " "
                                << NewNeutronFissionPositionsList[DataID][nn].getZ()<< " "
                                << NewNeutronFissionMomDirecsList[DataID][nn].getX()<< " "
                                << NewNeutronFissionMomDirecsList[DataID][nn].getY()<< " "
                                << NewNeutronFissionMomDirecsList[DataID][nn].getZ()<< "\n"  ;
                //MomDirecsListForCriticality[DataID][currentEvent*NumberOfEventInBatch+nn] = NewNeutronFissionMomDirecsList[DataID][nn];
            }

            CriticalityFile.close();
            return;
        }
    }

    //3- update the neutron fisson data for the next batch
    if(CurrentBatch < NumberOfBatch){// calculated previously

        //G4cout << " DataID " << DataID << " CurrentBatch "<< CurrentBatch << " CurrentBatch < NumberOfBatch , Generate Data" << G4endl;

        // update the neutron fisson data for the next batch
        for(G4int nn = 0 ;nn < NumberOfEventInBatch ;nn++){

            G4int rnd = G4UniformRand()*NewNeutronFissionEnergyList[DataID].size();
            G4int evntIc = CurrentBatch*NumberOfEventInBatch+nn;

            //G4cout << " 1 --------- " << G4endl;

            //G4cout << " 1 -------- NumberOfEventInBatch "<< NumberOfEventInBatch
            //       << " nn "<< nn
            //       << " rnd "<< rnd
            //       << " evntIc "<< evntIc
            //       << " NewNeutronFissionEnergyList[DataID].size() "<< NewNeutronFissionEnergyList[DataID].size()
            //       << " NewNeutronFissionEnergyList[DataID][rnd]    "<< NewNeutronFissionEnergyList[DataID][rnd]
            //       << " NewNeutronFissionPositionsList[DataID][rnd] "<< NewNeutronFissionPositionsList[DataID][rnd]
            //       << " NewNeutronFissionMomDirecsList[DataID][rnd] "<< NewNeutronFissionMomDirecsList[DataID][rnd]
            //       << G4endl;
            //G4cout << " 2 --------- " << evntIc << G4endl;
            EnergyListForCriticality[DataID][evntIc]    = NewNeutronFissionEnergyList[DataID][rnd];
            PositionsListForCriticality[DataID][evntIc] = NewNeutronFissionPositionsList[DataID][rnd];
            MomDirecsListForCriticality[DataID][evntIc] = NewNeutronFissionMomDirecsList[DataID][rnd];
            //G4cout << " 3 --------- " << G4endl;

        }

        ////FissionCapturesOfThreadsRanks[DataID][CurrentBatch][0] += NumberOfFissionNeutrons;
        ////FissionCapturesOfThreadsRanks[DataID][CurrentBatch][1] += NumberOfCapture;
        ////FissionCapturesOfThreadsRanks[DataID][CurrentBatch][2] += NumberOfFission;
        //NumberOfFissionNeutrons = FissionCapturesOfThreadsRanks[DataID][CurrentBatch][0];
        //NumberOfCapture         = FissionCapturesOfThreadsRanks[DataID][CurrentBatch][1];
        //NumberOfFission         = FissionCapturesOfThreadsRanks[DataID][CurrentBatch][2];
        NewNeutronFissionEnergyList[DataID].clear();
        NewNeutronFissionPositionsList[DataID].clear();
        NewNeutronFissionMomDirecsList[DataID].clear();

    }
    //G4cout << " 3 --------- " << G4endl;

    //l.unlock();
    //G4cout << "Capture" << G4endl;
}
void G4TRunAction::SetPoxEneOfFluxParticle(G4double X,G4double Y,G4double Z,G4double MX,G4double MY,G4double MZ,G4double E) {
    //G4cout << "Capture" << G4endl;
}







void G4TRunAction::CountOpticalInteractions(G4String n, G4String m ){
    OpticalPhotonInteractionRate[DataID][n][m]++;

    //for ( auto it = OpticalPhotonInteractionRate[DataID].begin(); it != OpticalPhotonInteractionRate[DataID].end(); ++it  ){
    //    G4cout << " DataID " << DataID  << " Interaction " << it->first << " Number "<<it->second << G4endl;
    //}

}

void G4TRunAction::ShowOpticalPhotonInteractionData(){

    G4cout << G4endl;

    for ( auto it = OpticalPhotonInteractionRate[DataID].begin(); it != OpticalPhotonInteractionRate[DataID].end(); ++it  ){
        for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2  ){
            G4cout << " DataID " << DataID  << " Particle " << it->first << " Interaction " << it2->first << " Number "<<it2->second << G4endl;
        }
    }

}

// called from endOfeventfrom event class , it get a vector of event absorbed energy by all organs
void G4TRunAction::useTime(G4int set_get_for)
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
