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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4TKleinNishinaCompton
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 15.03.2005
//
// Modifications:
// 18-04-05 Use G4ParticleChangeForGamma (V.Ivantchenko)
//
//
// Class Description:
//
// Implementation of gamma Compton scatteribg on free electron
// 

// -------------------------------------------------------------------
//

#ifndef G4TKleinNishinaCompton_h
#define G4TKleinNishinaCompton_h 1

#include "G4VEmModel.hh"

class G4ParticleChangeForGamma;

class G4TKleinNishinaCompton : public G4VEmModel
{

public:

    explicit G4TKleinNishinaCompton(const G4ParticleDefinition* p = nullptr,
                                   const G4String& nam = "Klein-Nishina");

    virtual ~G4TKleinNishinaCompton();

    virtual void Initialise(const G4ParticleDefinition*,
                            const G4DataVector&) override;

    virtual void InitialiseLocal(const G4ParticleDefinition*,
                                 G4VEmModel* masterModel) override;

    virtual G4double ComputeCrossSectionPerAtom(
            const G4ParticleDefinition*,
            G4double kinEnergy,
            G4double Z,
            G4double A,
            G4double cut,
            G4double emax) override;

    virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                   const G4MaterialCutsCouple*,
                                   const G4DynamicParticle*,
                                   G4double tmin,
                                   G4double maxEnergy) override;

    //virtual void StartTracking(G4Track*) override;
    void setSplittingNum(G4double S){ SplittingNum = S; }

protected:

    G4ParticleDefinition*     theGamma;
    G4ParticleDefinition*     theElectron;
    G4ParticleChangeForGamma* fParticleChange;
    G4double                  lowestSecondaryEnergy;

    G4int                  SplittingNum;

private:

    // hide assignment operator
    G4TKleinNishinaCompton & operator=
    (const G4TKleinNishinaCompton &right) = delete;
    G4TKleinNishinaCompton(const  G4TKleinNishinaCompton&) = delete;

    const G4Track* track;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
