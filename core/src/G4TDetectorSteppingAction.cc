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
//
/// \file G4TDetectorSteppingAction.cc
/// \brief Implementation of the G4TDetectorSteppingAction class

#include "G4TDetectorSteppingAction.hh"
//#include "G4TDetectorRun.hh"
#include "G4Event.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4TDetectorSteppingAction::G4TDetectorSteppingAction(G4TRunAction* runaction)
    : G4UserSteppingAction()
    , RunAction(runaction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4TDetectorSteppingAction::~G4TDetectorSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4TDetectorSteppingAction::UserSteppingAction(const G4Step* step)
{

    //G4cout << " UserSteppingAction " << G4endl;

    auto edep = step->GetTotalEnergyDeposit();

    //std::cout <<" edep:" << edep << " @@@@@@@@@@@@@@@@@@@@@@@" << std::endl;

    if (edep == 0.) return;

    auto reg = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName();

    //std::cout <<" edep:" << edep <<" reg:" << reg << " @@@@@@@@@@@@@@@@@@@@@@@" << std::endl;

    if ( reg == "D")
    {

        //std::cout <<"\n\n\n\n\n\n\n\n edep:" << edep <<" reg:" << reg << " @@@@@@@@@@@@@@@@@@@@@@@" << std::endl;

        RunAction->AddEventsDataInDetectors(G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(), reg, edep);
        G4EventManager::GetEventManager()->AbortCurrentEvent();
    }
    //std::cout << 1 << std::endl;

    //std::cout << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << " edep:" << edep << " @@@@@@@@@@@@@@@@@@@@@@@" << std::endl;

    //in all cases we need just the particles that pass to detector directly
    G4EventManager::GetEventManager()->AbortCurrentEvent();















    // For Lam, eliminate all other detectors except one with a known name
    ////if(DetName == "Det1"){
    ////    EventAction->setDetectorNameMap(DetName, step->GetTotalEnergyDeposit());
    ////}

    //G4cout << " DetName " << DetName << " EDep " << step->GetTotalEnergyDeposit() << G4endl;

    //EventAction->setDetectorNameMap(DetName, step->GetTotalEnergyDeposit());

    //G4cout << " setDetectorNameMap has been called " << DetName << G4endl;


    //static G4ParticleDefinition* opticalphoton = G4OpticalPhoton::OpticalPhotonDefinition();

    //const G4ParticleDefinition* particleDef = step->GetTrack()->GetDynamicParticle()->GetParticleDefinition();

    //const G4Event* event = G4RunManager::GetRunManager()->GetCurrentEvent();

    //void setDetectorNameMap(G4String SS, double EnerDepo) { decetorEnerDepoMap[SS] += EnerDepo; }


    //G4int eventID = event->GetEventID();
    //G4RunManager::GetRunManager()->GetCurrentEvent()->setEnerDepo(step->GetTotalEnergyDeposit());
    //G4RunManager::GetRunManager()->GetCurrentEvent()->setDetectorNameMap(step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName());

    //setEnerDepo(G4double EE) { EnerDepo = EE; }
    //void setDetectorNameMap(G4String SS) { decetorEnerDepoMap[SS] += EnerDepo; }

    //G4cout << " EventID " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
    //       << " Track ID "<< step->GetTrack()->GetTrackID()
    //       << " Track Energy "<< step->GetTrack()->GetKineticEnergy()
    //       << " particleDef " << particleDef->GetParticleName()
    //       << " procname " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
    //       << " RegionPre " << step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName()
    //       //<< " RegionPost " << step->GetPostStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName()
    //       << G4endl;

    //RunAction->CountOpticalInteractions(step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName()+"_"+particleDef->GetParticleName(), step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());

    //if(particleDef == opticalphoton)
    //{
    //    //G4cout << " opticalphoton " << G4endl;

    //    G4StepPoint* endPoint = step->GetPostStepPoint();
    //    const G4VProcess* pds = endPoint->GetProcessDefinedStep();
    //    G4String procname     = pds->GetProcessName();

    //    //G4cout << " particleDef " << particleDef->GetParticleName() << " procname " << procname << G4endl;
    //    RunAction->CountOpticalInteractions(particleDef->GetParticleName(), procname);
    //    RunAction->CountOpticalInteractions(particleDef->GetParticleName(), procname+"_"+step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName());

    //    //G4cout << " EventID " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << " Track ID "<< step->GetTrack()->GetTrackID()<< " particleDef " << particleDef->GetParticleName() << " procname " << procname <<" Region Pre " << step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName()<<" Region Post " << step->GetPostStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName() << G4endl;

    //    //if( step->GetPostStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName() == "SiO2"){
    //    //}
    //    //if(procname.compare("OpRayleigh") == 0) RunAction->CountOpticalInteractions(procname);
    //    //else if(procname.compare("OpAbsorption") == 0) RunAction->CountOpticalInteractions(procname);
    //    //else if(procname.compare("OpMieHG") == 0) RunAction->CountOpticalInteractions(procname);

    //    // for boundary scattering, process name in 'transportation'.
    //    // Need to check differently:
    //    if(endPoint->GetStepStatus() == fGeomBoundary)
    //    {
    //        G4OpBoundaryProcessStatus theStatus = Undefined;
    //        G4ProcessManager* opManager         = opticalphoton->GetProcessManager();
    //        G4int n_proc = opManager->GetPostStepProcessVector(typeDoIt)->entries();
    //        G4ProcessVector* postStepDoItVector = opManager->GetPostStepProcessVector(typeDoIt);
    //        for(G4int i = 0; i < n_proc; ++i)
    //        {
    //            G4VProcess* currentProcess = (*postStepDoItVector)[i];
    //            G4OpBoundaryProcess* opProc = dynamic_cast<G4OpBoundaryProcess*>(currentProcess);
    //            if(opProc) theStatus = opProc->GetStatus();
    //        }
    //        if(theStatus != Undefined && theStatus != NotAtBoundary && theStatus != StepTooSmall)
    //        {
    //            //fEventAction->AddBoundary();
    //            if( step->GetPostStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName() == "SiO2"){
    //                RunAction->CountOpticalInteractions(particleDef->GetParticleName(), "BoundarySiO2");
    //                //step->GetTrack()->SetTrackStatus(fStopAndKill);
    //                //G4cout << " EventID "    << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
    //                //       << " Track ID "   << step->GetTrack()->GetTrackID()
    //                //       << " particleDef "<< particleDef->GetParticleName()
    //                //       << " procname "   << procname
    //                //       << " Region Pre " << step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName()
    //                //       << " Region Post "<< step->GetPostStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName()
    //                //       << " Energy "     << step->GetTrack()->GetKineticEnergy()
    //                //       << " Direction "  << step->GetTrack()->GetMomentumDirection()
    //                //       << G4endl;
    //            }

    //        }
    //    }
    //}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
