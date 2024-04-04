
#ifndef G4TEventActionForCriticality_h
#define G4TEventActionForCriticality_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
//#include <map>

class G4TRunActionForCriticality;
class G4Event;
class G4TEventActionForCriticality : public G4UserEventAction
{
public:
    G4TEventActionForCriticality(G4TRunActionForCriticality* runAction);
    ~G4TEventActionForCriticality();

public:

    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);



private:

    G4TRunActionForCriticality* RunAction;
};
#endif
