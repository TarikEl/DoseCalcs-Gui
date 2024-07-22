
#ifndef G4TDetectorGeometry_h
#define G4TDetectorGeometry_h 1

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"

class G4Box;
class G4LogicalVolume;
class G4Sphere;
class G4Tubs;

class G4Box;
class G4Element;
class G4LogicalVolume;
class G4Material;
class G4MaterialPropertiesTable;
class G4Sphere;
class G4Tubs;
class G4VPhysicalVolume;

class G4TDetectorGeometry {
public:
    G4TDetectorGeometry();
    virtual ~G4TDetectorGeometry();

    G4VPhysicalVolume* PlacePMTs();
    G4VPhysicalVolume* ConstructLXe();
    G4VPhysicalVolume* ConstructDetector();
public:

    G4VPhysicalVolume* LXeMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                       G4LogicalVolume* pMotherLogical, G4bool pMany, G4int pCopyNo);
    G4VPhysicalVolume* LXeWLSSlab(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                           G4LogicalVolume* pMotherLogical, G4bool pMany,
                           G4int pCopyNo);
    G4VPhysicalVolume* LXeWLSFiber(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                             G4LogicalVolume* pMotherLogical, G4bool pMany,
                             G4int pCopyNo);

    G4LogicalVolume* GetLogPhotoCath() { return fPhotocath_log; }
    G4LogicalVolume* GetLogScint() { return fScint_log; }

    std::vector<G4ThreeVector> GetPmtPositions() { return fPmtPositions; }

private:
    G4VPhysicalVolume* fMainVolume;

    void SurfaceProperties();

    void PlacePMTs(G4LogicalVolume* pmt_Log, G4RotationMatrix* rot, G4double& a,
                   G4double& b, G4double da, G4double db, G4double amin,
                   G4double bmin, G4int na, G4int nb, G4double& x, G4double& y,
                   G4double& z, G4int& k);

    void CopyValues();

    G4double fScint_x;
    G4double fScint_y;
    G4double fScint_z;
    G4double fD_mtl;
    G4int fNx;
    G4int fNy;
    G4int fNz;
    G4double fOuterRadius_pmt;
    //G4bool fSphereOn;
    G4double fRefl;

    // Basic Volumes
    //
    G4Box* fScint_box;
    G4Box* fHousing_box;
    G4Tubs* fPmt;
    G4Tubs* fPhotocath;
    G4Sphere* fSphere;

    // Logical volumes
    //
    G4LogicalVolume* fScint_log;
    G4LogicalVolume* fHousing_log;
    G4LogicalVolume* fPmt_log;
    G4LogicalVolume* fPhotocath_log;
    G4LogicalVolume* fSphere_log;

    // Sensitive Detectors positions
    std::vector<G4ThreeVector> fPmtPositions;







public:

    // Functions to modify the geometry
    void SetDimensions(G4ThreeVector);
    void SetHousingThickness(G4double);
    void SetNX(G4int);
    void SetNY(G4int);
    void SetNZ(G4int);
    void SetPMTRadius(G4double);
    void SetSaveThreshold(G4int);

    // Get values
    G4int GetNX() const { return fNx; };
    G4int GetNY() const { return fNy; };
    G4int GetNZ() const { return fNz; };
    G4int GetSaveThreshold() const { return fSaveThreshold; };
    G4double GetScintX() const { return fScint_x; }
    G4double GetScintY() const { return fScint_y; }
    G4double GetScintZ() const { return fScint_z; }
    G4double GetHousingThickness() const { return fD_mtl; }
    G4double GetPMTRadius() const { return fOuterRadius_pmt; }
    G4double GetSlabZ() const { return fSlab_z; }

    void SetHousingReflectivity(G4double);
    G4double GetHousingReflectivity() const { return fRefl; }

    void SetWLSSlabOn(G4bool b);
    G4bool GetWLSSlabOn() const { return fWLSslab; }

    void SetMainVolumeOn(G4bool b);
    G4bool GetMainVolumeOn() const { return fMainVolumeOn; }

    void SetNFibers(G4int n);
    G4int GetNFibers() const { return fNfibers; }

    void SetMainScintYield(G4double);
    void SetWLSScintYield(G4double);

private:
    void DefineMaterials();

    G4Box* fExperimentalHall_box;
    G4LogicalVolume* fExperimentalHall_log;
    G4VPhysicalVolume* fExperimentalHall_phys;

    // Materials & Elements
    G4Material* fLXe;
    G4Material* fAl;
    G4Element* fN;
    G4Element* fO;
    G4Material* fAir;
    G4Material* fVacuum;
    G4Element* fC;
    G4Element* fH;
    G4Material* fGlass;
    G4Material* fPstyrene;
    G4Material* fPMMA;
    G4Material* fPethylene1;
    G4Material* fPethylene2;

    // Geometry

    G4int fSaveThreshold;
    G4int fNfibers;
    G4bool fSphereOn;
    G4bool fWLSslab;
    G4bool fMainVolumeOn;
    G4double fSlab_z;


    G4MaterialPropertiesTable* fLXe_mt;
    G4MaterialPropertiesTable* fMPTPStyrene;

    G4LogicalVolume* fScintSlab_log;


    G4LogicalVolume* fClad2_log;

    G4double fFiber_rmin;
    G4double fFiber_rmax;
    G4double fFiber_z;
    G4double fFiber_sphi;
    G4double fFiber_ephi;

    G4double fClad1_rmin;
    G4double fClad1_rmax;
    G4double fClad1_z;
    G4double fClad1_sphi;
    G4double fClad1_ephi;

    G4double fClad2_rmin;
    G4double fClad2_rmax;
    G4double fClad2_z;
    G4double fClad2_sphi;
    G4double fClad2_ephi;


    // ******************************* for Neutron Detector
public:
    G4VPhysicalVolume* ConstructNeutronDetector();

private:

  G4double fPressureInTorr;
  G4double fTemperature;
  G4String fGasType;
  G4Material* fGasMaterial;

  G4double inchtocm;
  G4int fNumGrids;
  G4double fGridSize;
  G4double fGridRadius;
  G4double fGridDist;
  G4double fScintDist;
  G4double fWireThickness;
  G4bool fUseInches;

  G4LogicalVolume* fWorldLogical;
  G4LogicalVolume* fTargetLogical;
  G4LogicalVolume* fDetectLogical;
  G4LogicalVolume* fFoilLogical;
  std::vector<G4LogicalVolume*> fGridLogical;
  std::vector<G4LogicalVolume*> fWireGridLogical[49];
  G4LogicalVolume* fScintLogical;

  G4UserLimits* fStepLimit;

  void ConstructMaterials();
  void SetAttributes();
};

#endif
