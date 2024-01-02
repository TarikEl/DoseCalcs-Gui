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
// TETPSEnergyDeposit.hh
// \file   MRCP_GEANT4/Internal/include/TETPSEnergyDeposit.hh
// \author Haegin Han
//

#ifndef TETPSEnergyDeposit_h
#define TETPSEnergyDeposit_h 1

#include "G4PSEnergyDeposit.hh"
#include "G4TTETModelImport.hh"

// *********************************************************************
// This is the scorer based on G4PSEnergyDeposit class.
// -- GetIndex: Return the organ ID instead of copy number automatically
//              given by Parameterisation geometry.
// *********************************************************************

class TETPSEnergyDeposit : public G4PSEnergyDeposit
{
   public:
      TETPSEnergyDeposit(G4String name,G4TTETModelImport* _tetData);
      virtual ~TETPSEnergyDeposit();

  protected:
      virtual G4int GetIndex(G4Step*);

  private:
      G4TTETModelImport* tetData;
};

#endif

