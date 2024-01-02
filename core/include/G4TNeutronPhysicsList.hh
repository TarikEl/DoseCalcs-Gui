#ifndef G4TNEUTRONPHYSICSLIST_HH
#define G4TNEUTRONPHYSICSLIST_HH

#include "G4VModularPhysicsList.hh"

class G4TNeutronPhysicsList : public G4VModularPhysicsList
{
public:
    G4TNeutronPhysicsList(const G4String& name="neutron");
    ~G4TNeutronPhysicsList();

public:
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    virtual void SetCuts();

private:
    G4bool              fThermal;
};
#endif
