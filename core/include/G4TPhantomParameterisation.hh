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
// Author: Tarik Elghalbzouri,  Abdelmalek Essaâdi University,
// faculty of sciences Tetouane, morocco. email : telghalbzouri@uae.ac.ma
//
// This application is based on code developed by :
// G. Guerrieri, University of Genova, Italy .
// S. Guatelli. University of Wollongong, Australia.
//

#ifndef G4TPhantomParameterisation_HH
#define G4TPhantomParameterisation_HH

#include <map>
#include "G4Colour.hh"
#include "G4PhantomParameterisation.hh"

class G4VisAttributes;

extern std::map<G4String, G4Colour> RegionNameColour;
extern G4String* CopyNumberRegionNameMap;
extern bool ForceSolid;

class G4TPhantomParameterisation : public G4PhantomParameterisation
{

public:  // with description

    G4TPhantomParameterisation();
    ~G4TPhantomParameterisation();

    virtual G4Material* ComputeMaterial(const G4int repNo, G4VPhysicalVolume *currentVol, const G4VTouchable *parentTouch=0);

    void setUseLogVolColour(G4bool e ){
        UseLogVolColour = e;
    }

    std::vector<unsigned int> RegionCopyNumber; void SetRegionCopyNumber(std::vector<unsigned int> hh ){ RegionCopyNumber = hh; }

    void SetVoxelXYZData(G4int NX, G4int NY, G4int NZ){
        NVoxelX = NX;
        NVoxelX = NY;
        NVoxelX = NZ;
    }


private:

    G4bool UseLogVolColour;
    //G4Colour* RegionCopyNumberColour;
    G4int NVoxelX, NVoxelY, NVoxelZ;

    void ReadColourData();

private:
    std::map<G4String,G4VisAttributes*> fColours;
};


#endif
