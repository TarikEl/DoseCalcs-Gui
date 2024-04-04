
#ifndef G4TDetectorGeometry_h
#define G4TDetectorGeometry_h 1

#include "G4VPhysicalVolume.hh"

class G4TDetectorGeometry {
  public:
  G4TDetectorGeometry();
  virtual ~G4TDetectorGeometry();

  G4VPhysicalVolume* ConstructDetector();

};

#endif
