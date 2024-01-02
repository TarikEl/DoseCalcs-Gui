#include "G4TComptonProcess.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include <assert.h>
#include <vector>
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

namespace {

G4Mutex bremspe_Mutex = G4MUTEX_INITIALIZER;

}  

// Initialise static data
G4int fNSplit = 1;
G4int fNSecondaries = 0;
G4bool fActive = true;

G4TComptonProcess::G4TComptonProcess() {

    //std::cout << " G4TComptonProcess initialization " << std::endl ;

    auto idthread = std::to_string(G4Threading::G4GetThreadId());


}

G4TComptonProcess::~G4TComptonProcess() {}




G4VParticleChange* G4TComptonProcess::PostStepDoIt(const G4Track& track, const G4Step& step) {

    //std::cout << __FUNCTION__ << std::endl ;

    G4AutoLock _m(&bremspe_Mutex);

    fNSplit = 4;

    G4VParticleChange* particleChange(0);

    if (!fActive) {
        particleChange = pRegProcess->PostStepDoIt(track, step);
        assert (0 != particleChange);
        //fNSecondaries += particleChange->GetNumberOfSecondaries();
        return particleChange;
    }

    //G4ThreeVector direction = track.GetMomentumDirection();
    //G4double theta  = std::acos(std::abs(direction.z()));

    //if ((track.GetParentID() == 0 and direction.z()> 0.0))
    //std::cout << " Track ID=" << track.GetTrackID() << std::endl ;

    if (track.GetWeight() == 1)
    {
        //std::cout << " ------------------------ Track weigh=" << track.GetWeight() << " with compton process --------------------------" << std::endl ;

        //    Do compton splitting
        assert (fNSplit > 0);
        G4int i(0);

        // Secondary store
        std::vector<G4Track*> secondaries;
        //secondaries.reserve(fNSplit);

        particleChange = pRegProcess->PostStepDoIt(track, step); particleChange->SetVerboseLevel(0);

        //std::cout << " Number Of Secondaries" << particleChange->GetNumberOfSecondaries() << std::endl ;
        G4double ParentWeight = 1. ;//track.GetWeight();
        G4double weight = 0.2 ; //track.GetWeight()/(fNSplit);
        //std::cout << "particle 1 " << particleChange->GetSecondary(0)->GetDefinition()->GetParticleName() << std::endl ;
        //std::cout << "particle 2 " << particleChange->GetSecondary(1)->GetDefinition()->GetParticleName() << std::endl ;

        //particleChange->SetParentWeightByProcess(true);
        //particleChange->ProposeParentWeight(weight);
        //std::cout << "2 weight " << weight << " ParentWeight " << particleChange->GetParentWeight() << " " << track.GetWeight() << std::endl ;

        secondaries.push_back(new G4Track(*(particleChange->GetSecondary(1))));

        // Configure particleChange to handle multiple secondaries. Other data is unchanged

        // Loop over PostStepDoIt method to generate multiple secondaries.
        for (i=0; i<fNSplit; i++)
        {
            //particleChange = pRegProcess->PostStepDoIt(track, step);
            //assert (0 != particleChange);
            //particleChange->SetVerboseLevel(0);
            //G4int j(0);
            //std::cout << track.GetParticleDefinition()->GetParticleName() << " Number of secondaries " << particleChange->GetNumberOfSecondaries() << std::endl ;
            //if(i==0 && particleChange->GetNumberOfSecondaries() > 0){
            //secondaries.push_back(new G4Track(*(particleChange->GetSecondary(0))));
            //}

            secondaries.push_back(new G4Track(*(particleChange->GetSecondary(0))));

            //std::cout << " secondaries " << secondaries[i]->GetDefinition()->GetParticleName() << std::endl ;

        }

        //particleChange->DumpInfo();

        particleChange->SetNumberOfSecondaries(fNSplit+1);
        particleChange->SetSecondaryWeightByProcess(true);

        // Add all secondaries
        std::vector<G4Track*>::iterator iter = secondaries.begin();

        // add electron secondary
        G4Track* myTrack = *iter;
        myTrack->SetWeight(ParentWeight);
        particleChange->AddSecondary(myTrack);
        //std::cout << " secondary " << myTrack->GetParticleDefinition()->GetParticleName();
        //std::cout << " W=" << myTrack->GetWeight();
        //std::cout << " ParentW=" << particleChange->GetParentWeight() << std::endl ;

        iter++;

        // add gamma secondary
        while (iter != secondaries.end())
        {
            G4Track* myTrack = *iter;
            myTrack->SetWeight(weight);
            particleChange->AddSecondary(myTrack);
            std::cout << " secondary " << myTrack->GetParticleDefinition()->GetParticleName()
                      << " W=" << myTrack->GetWeight()
                      << " ParticleChangeW=" << particleChange->GetWeight()
                      << " ParentW=" << particleChange->GetParentWeight()
                      << " E=" << myTrack->GetTotalEnergy()
                      << " M=" << myTrack->GetMomentumDirection()
                      << " MD=" << myTrack->GetMomentum()
                      << std::endl;

            iter++;
        }

        std::cout << " Number Of Secondaries = " << particleChange->GetNumberOfSecondaries() << std::endl ;

        return particleChange;
    }
    else
    {   particleChange = pRegProcess->PostStepDoIt(track, step);
        assert (0 != particleChange);
        //fNSecondaries += particleChange->GetNumberOfSecondaries();
        return particleChange;
    }

    _m.unlock();
}

void G4TComptonProcess::SetNSplit(G4int nSplit)
{
    fNSplit = nSplit;
}

void G4TComptonProcess::SetIsActive(G4bool active)
{
    fActive = active;
}

G4bool G4TComptonProcess::GetIsActive()
{
    return fActive;
}

G4int G4TComptonProcess::GetNSplit()
{
    return fNSplit;
}

G4int G4TComptonProcess::GetNSecondaries()
{
    return fNSecondaries;
}
