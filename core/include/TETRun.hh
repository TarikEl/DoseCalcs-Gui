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
// TETRun.hh
// \file   MRCP_GEANT4/Internal/include/TETRun.hh
// \author Haegin Han
//

#ifndef TETRun_h
#define TETRun_h 1

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include "G4SDManager.hh"

typedef std::map<G4int, std::pair<G4double, G4double>> EDEPMAP;
typedef std::map<G4int,unsigned long long int> NumStepMAP;

// *********************************************************************
// This is G4Run class that sums up energy deposition from each event.
// The sum of the square of energy deposition was also calculated to
// produce the relative error of the dose.
// -- RecordEvent: Sum up the energy deposition and the square of it.
//                 The sums for each organ were saved as the form of
//                 std::map.
// -- Merge: Merge the data calculated in each thread.
// *********************************************************************

class TETRun : public G4Run 
{
public:
	TETRun();
	virtual ~TETRun();

	virtual void RecordEvent(const G4Event*);
	void ConstructMFD(const G4String& mfdName);
    virtual void Merge(const G4Run*);

    EDEPMAP* GetEdepMap() {return &edepMap;};
    NumStepMAP* GetnumstepsMap() {return &numstepsMap;};

private:
    EDEPMAP edepMap;
    NumStepMAP numstepsMap;
};

#endif
