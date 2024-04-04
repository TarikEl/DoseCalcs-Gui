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
#include "G4TSteppingActionForCriticality.hh"
#include "G4TRunActionForCriticality.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4Neutron.hh"
#include "G4TEventActionForCriticality.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronicProcessType.hh"

//G4ThreadLocal G4double G4TSteppingActionForCriticality::edep;
//G4ThreadLocal G4String G4TSteppingActionForCriticality::reg;
//G4ThreadLocal unsigned int  G4TSteppingActionForCriticality::CN;

//G4TSteppingActionForCriticality::G4TSteppingActionForCriticality(G4TEventActionForCriticality* eventAction) : G4UserSteppingAction(), EventAction(eventAction) {}
G4TSteppingActionForCriticality::G4TSteppingActionForCriticality(G4TRunActionForCriticality* runAction) : G4UserSteppingAction(), RunAction(runAction) {}

G4TSteppingActionForCriticality::~G4TSteppingActionForCriticality(){}


void G4TSteppingActionForCriticality::UserSteppingAction(const G4Step* step)
{

    if (step->GetTrack()->GetDefinition() == G4Neutron::NeutronDefinition()) {

        //G4String Particlename = step->GetTrack()->GetParticleDefinition()->GetParticleName();
        //if(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "nFission"){
        //std::cout << " EventID " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
        //          << " TrackID " << step->GetTrack()->GetTrackID()
        //          << " Particlename " << Particlename
        //          << " ProcessName " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
        //          << std::endl;
        //}

        //const std::vector<const G4Track*>* secondaryTracks = step->GetSecondaryInCurrentStep();

        //for (const G4Track* secondaryTrack : *secondaryTracks) {
        //    // Access information about each secondary track
        //    G4ParticleDefinition* particleDefinition = secondaryTrack->GetDefinition();
        //    G4String particleName = particleDefinition->GetParticleName();

        //    //std::cout << " EventID " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
        //    //          << " TrackID " << step->GetTrack()->GetTrackID()
        //    //          << " PrimaryParticlename " << step->GetTrack()->GetParticleDefinition()->GetParticleName()
        //    //          << " ProcessName " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
        //    //          << " SecondaryParticleName " << particleName
        //    //          << " SecondaryParticleEnergy " << secondaryTrack->GetKineticEnergy()
        //    //          << " PrimaryParticleEnergy " << step->GetTrack()->GetKineticEnergy()
        //    //          << std::endl;

        //    RunAction->CountParticleProductionByNeutron(particleDefinition->GetParticleName());
        //}

        //RunAction->SetPoxEneOfFluxParticle(step->GetTrack()->GetPosition().getX(),
        //                                  step->GetTrack()->GetPosition().getY(),
        //                                  step->GetTrack()->GetPosition().getZ(),
        //                                  step->GetTrack()->GetKineticEnergy());

        RunAction->CountParticleProductionByNeutron(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
        if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "nFission") {

            RunAction->CountFission();

            //std::cout << " nFission " << std::endl;
            for (const G4Track* secondaryTrack : *step->GetSecondaryInCurrentStep()) {
                if(secondaryTrack->GetDefinition()->GetParticleName() == "neutron"){
                    //std::cout << " E " << secondaryTrack->GetKineticEnergy()
                    //          << " X " << secondaryTrack->GetPosition().getX()
                    //          << " Y " << secondaryTrack->GetPosition().getY()
                    //          << " Z " << secondaryTrack->GetPosition().getZ()
                    //          << " Weight " << secondaryTrack->GetWeight()

                    //          << std::endl;
                    RunAction->SetPoxEneForBatch(secondaryTrack->GetPosition().getX(),
                                                 secondaryTrack->GetPosition().getY(),
                                                 secondaryTrack->GetPosition().getZ(),
                                                 secondaryTrack->GetMomentumDirection().getX(),
                                                 secondaryTrack->GetMomentumDirection().getY(),
                                                 secondaryTrack->GetMomentumDirection().getZ(),
                                                 secondaryTrack->GetKineticEnergy());
                }
            }

            G4EventManager::GetEventManager()->AbortCurrentEvent();
        }
        else if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "nCapture") {
            RunAction->CountAbsorption();
            G4EventManager::GetEventManager()->AbortCurrentEvent();
        }
    }
    //else{// we need just neutron for criticality calculations
    //    step->GetTrack()->SetTrackStatus(fStopAndKill);
    //}
}
