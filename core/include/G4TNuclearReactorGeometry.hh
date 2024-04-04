
#ifndef G4TNuclearReactorGeometry_h
#define G4TNuclearReactorGeometry_h 1
#include "G4VPhysicalVolume.hh"

class G4TNuclearReactorGeometry {
  public:
  G4TNuclearReactorGeometry();
  virtual ~G4TNuclearReactorGeometry();

   G4VPhysicalVolume* ConstructNuclearReactor();
};

#endif
