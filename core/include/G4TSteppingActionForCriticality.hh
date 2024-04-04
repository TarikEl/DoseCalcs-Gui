
#ifndef G4TSteppingActionForCriticality_h
#define G4TSteppingActionForCriticality_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4TRunActionForCriticality.hh"

//class G4TVolumeConstruction;
class G4TEventActionForCriticality;

class G4TSteppingActionForCriticality : public G4UserSteppingAction
{
public:
    G4TSteppingActionForCriticality(G4TRunActionForCriticality* RunAction);
    virtual ~G4TSteppingActionForCriticality();

    virtual void UserSteppingAction(const G4Step* step);

private:
    //const G4TVolumeConstruction* fDetConstruction;
    G4TRunActionForCriticality* RunAction;

    //G4ThreadLocal static G4double edep;
    //G4ThreadLocal static G4String reg;
    //G4ThreadLocal static unsigned int CN;
};

#endif
