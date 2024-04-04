
#ifndef G4TTackingActionForCriticality_h
#define G4TTackingActionForCriticality_h 1

#include "G4UserTrackingAction.hh"


class G4Track;
class G4TTackingActionForCriticality  : public G4UserTrackingAction {
  public:
  G4TTackingActionForCriticality();
  virtual ~G4TTackingActionForCriticality();

  void PreUserTrackingAction(G4Track*);
  void PostUserTrackingAction(G4Track*);

};


#endif
