
#ifndef G4TComptonProcess_HH
#define G4TComptonProcess_HH 1
#include "G4WrapperProcess.hh"

class G4TComptonProcess : public G4WrapperProcess {
  public:
  G4TComptonProcess();
  virtual ~G4TComptonProcess();
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);
  void SetNSplit(G4int);
  void SetIsActive(G4bool);
  G4bool GetIsActive();
  G4int GetNSplit();
  G4int GetNSecondaries();
};

#endif
