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
// TETPhysicsList.hh
// \file   MRCP_GEANT4/Internal/include/TETPhysicsList.hh
// \author Haegin Han
//

#ifndef TETPhysicsList_h
#define TETPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

class G4VPhysicsConstructor;

// *********************************************************************
// Please note that only basic physics were registered in this
// ModularPhysicsList, and rather precise models were used for the
// production of dose coefficients provided in the ICRP Publication.
// -- SetCuts: cut values were set as default values. This can be
//             modified according to specific purposes or applications.
// *********************************************************************

class TETPhysicsList: public G4VModularPhysicsList
{
public:
	TETPhysicsList();
	virtual ~TETPhysicsList();

	virtual void SetCuts();
};

#endif
