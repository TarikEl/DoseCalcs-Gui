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
#include "G4TEventActionForCriticality.hh"
#include "G4TRunActionForCriticality.hh"
#include "G4RunManager.hh"

extern G4int NumberOfEventInBatch;
extern G4int NumberOfBatch;
extern std::map<G4int,std::map<G4int,std::vector<G4double>>> FissionCapturesOfThreadsRanks;
extern std::vector<G4double> KeffectiveInEachBatch;

//namespace
//{
//G4Mutex	mutex = G4MUTEX_INITIALIZER;
//}

G4TEventActionForCriticality::G4TEventActionForCriticality(G4TRunActionForCriticality* runAction): G4UserEventAction(),RunAction(runAction)
{
}
G4TEventActionForCriticality::~G4TEventActionForCriticality()
{
    //G4MUTEXDESTROY(mutex);
}

void G4TEventActionForCriticality::BeginOfEventAction(const G4Event* event) {

    //fFissionCount = 0;
    //fAbsorptionCount = 0;
}
// 2
// 1 eventAction, Last event
// 2 IntializeBatchData RunAction
// 3 Generate Primaries
// 4 GetEventData

void G4TEventActionForCriticality::EndOfEventAction(const G4Event* event){
    //G4cout << "k-effective: " << keff << G4endl;
    //RunAction->pushNumberOfFissionAndCapture(fFissionCount,fAbsorptionCount);
    //G4AutoLock l(&mutex);
    //l.lock();
    if( (G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()+1)%NumberOfEventInBatch == 0){

        RunAction->IntializeBatchData();

        //G4cout << " End Ev --------- " << G4endl;

        //the first generation of user source not used to make the keffictive from the right source distribution of fission sites
        //if((G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()+1) != NumberOfEventInBatch)
        //{

        //    //G4double NumCap = 0., NumFiss = 0.,NumFissN = 0.;
        //    //bool isready = true;
        //    //for (G4int  BatchInc = 0;  BatchInc < NumberOfBatch;  BatchInc++ ) {
        //    //    for ( auto it = FissionCapturesOfThreadsRanks.begin(); it != FissionCapturesOfThreadsRanks.end(); ++it  ){

        //    //        if(it->first == -1){ continue;}

        //    //        G4cout << "For BatchInc = " << BatchInc << " For Thread = " << it->first << " NumberOfFissionN " << it->second[ BatchInc][0] << " NumCaptures " << it->second[ BatchInc][1] << " NumFission " << it->second[ BatchInc][2] << G4endl;

        //    //        NumFissN += it->second[ BatchInc][0];
        //    //        NumCap   += it->second[ BatchInc][1];
        //    //        NumFiss  += it->second[ BatchInc][2];
        //    //    }
        //    //    G4cout << "For BatchInc = " << BatchInc << " Total for :   NumberOfFissionN " << NumFissN << " NumCaptures " << NumCap  << " NumFission " << NumFiss << G4endl;

        //    //}
        //    //RunAction->IntializeBatchData();


        //    G4double NumCap = 0., NumFiss = 0.,NumFissN = 0.;

        //    //for ( auto it = FissionCapturesOfThreadsRanks.begin(); it != FissionCapturesOfThreadsRanks.end(); ++it  ){
        //    //    if(it->first == -1){ continue;}
        //    //    //G4cout << "***************** DataID = " << it->first << " NumberOfFissionN " << it->second[0] << " NumCaptures " << it->second[1]  << " NumFission " << it->second[2] << G4endl;
        //    //    NumFissN += it->second[0];
        //    //    NumCap   += it->second[1];
        //    //    NumFiss  += it->second[2];
        //    //}

        //    for (G4int  BatchInc = 0;  BatchInc < NumberOfBatch;  BatchInc++ ) {
        //        G4cout << "AAA NumberOfBatch = " << NumberOfBatch << "  BatchInc " <<  BatchInc  << " KeffectiveInEachBatch[ BatchInc] " << KeffectiveInEachBatch[ BatchInc] << G4endl;

        //        NumFissN= 0;
        //        NumCap  = 0;
        //        NumFiss = 0;

        //        bool isready = true;
        //        if(KeffectiveInEachBatch[ BatchInc] <= 0){
        //            for ( auto it = FissionCapturesOfThreadsRanks.begin(); it != FissionCapturesOfThreadsRanks.end(); ++it  ){

        //                if(it->first == -1){ continue;}

        //                G4cout << "BBB For Thread = " << it->first << " For BatchInc = " << BatchInc << " NumberOfFissionN " << it->second[ BatchInc][0] << " NumCaptures " << it->second[ BatchInc][1] << " NumFission " << it->second[ BatchInc][2] << G4endl;

        //                if(it->second[ BatchInc][0] == 0){
        //                    isready = false;
        //                    break;
        //                }else if(it->second[ BatchInc][1] == 0){
        //                    isready = false;
        //                    break;
        //                }else if(it->second[ BatchInc][2] == 0){
        //                    isready = false;
        //                    break;
        //                }

        //                NumFissN += it->second[ BatchInc][0];
        //                NumCap   += it->second[ BatchInc][1];
        //                NumFiss  += it->second[ BatchInc][2];
        //            }
        //            G4cout << "CCC Total for :     For BatchInc = " << BatchInc << " NumberOfFissionN " << NumFissN << " NumCaptures " << NumCap  << " NumFission " << NumFiss << G4endl;

        //        }else{
        //            continue;
        //        }

        //        if(isready){

        //            //G4cout << "***************** Total for : NumberOfFissionN = " << NumFissN << " NumCaptures " << NumCap  << " NumFission " << NumFiss << G4endl;

        //            //G4double keff = static_cast<G4double>(RunAction->getNumberOfFissionNeutrons()/(RunAction->getNumberOfFission()+4.94066e-324)) * static_cast<G4double>(RunAction->getNumberOfFission()/ (RunAction->getNumberOfCapture() +4.94066e-324)); // Avoid division by zero
        //            G4double keff = static_cast<G4double>(NumCap/(NumFiss+4.94066e-324)) * static_cast<G4double>(NumFiss/ (NumFissN +4.94066e-324)); // Avoid division by zero

        //            if( BatchInc == 0){
        //                G4cout << "***************** k-eff = " << keff << G4endl;
        //            }else if( BatchInc > 0){
        //                G4double relative_error = abs((keff - KeffectiveInEachBatch[ BatchInc-1]) / keff);
        //                G4cout << "***************** k-eff = " << keff << " +/- " << relative_error << G4endl;
        //            }

        //            RunAction->IntializeBatchData();
        //            KeffectiveInEachBatch[ BatchInc] = keff;
        //        }else {
        //            break;
        //        }
        //    }

        //    //G4double keff = static_cast<G4double>(RunAction->getNumberOfFissionNeutrons()/(RunAction->getNumberOfFission()+4.94066e-324)) * static_cast<G4double>(RunAction->getNumberOfFission()/ (RunAction->getNumberOfCapture() +4.94066e-324)); // Avoid division by zero
        //    //G4double relative_error = abs((keff - RunAction->getKeffective()) / keff);
        //    ////G4cout << " \n\n\n\n\nBatch Inialization \n TerminatedEventID " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << " NumberOfEventInBatch "<< NumberOfEventInBatch << " NumberOfFission: " << RunAction->getNumberOfFission() << " NumberOfCapture: " << RunAction->getNumberOfCapture() << " k-effective: " << keff << " +/- " << relative_error << G4endl;

        //    //if(RunAction->getKeffective() == 0){
        //    //    G4cout << "***************** k-eff = " << keff << G4endl;
        //    //}else{
        //    //    G4cout << "***************** k-eff = " << keff << " +/- " << relative_error << G4endl;
        //    //}
        //    //RunAction->IntializeBatchData();
        //    //RunAction->setKeffective(keff); // because we need keff in IntializeBatchData() method

        //}else{
        //    RunAction->IntializeBatchData();
        //}


    }
    //l.unlock();
}

