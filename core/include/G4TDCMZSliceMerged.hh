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
// $Id: G4TDCMZSliceMerged.hh 101109 2016-11-07 08:14:53Z gcosmo $
//
/// \file medical/DICOM/include/G4TDCMZSliceMerged.hh
/// \brief Definition of the G4TDCMZSliceMerged class
//
//
// The code was written by :
//      * Jonathan Madsen : jonathan.madsen@cern.ch (12/18/2012)
//
//
//*******************************************************


#ifndef G4TDCMzslicemerged_hh_
#define G4TDCMzslicemerged_hh_

#include <map>
#include "globals.hh"

#include "G4TDCMZSliceHeader.hh"

class G4TDCMZSliceMerged
{
public:
    // Constructor and Destructors
    G4TDCMZSliceMerged();
    ~G4TDCMZSliceMerged();

public:
    // called from DicomHandler::ReadFile(FILE* dicom, char* filename2)
    void AddZSlice(G4TDCMZSliceHeader* val) {
        fSlices[val->GetSliceLocation()] = val;
    }

    // called from DicomHandler::CheckFileFormat()
    void CheckSlices();

    inline void DumpExcessMemory();

private:
    // Private functions

private:
    // Private variables
    std::map<G4double,G4TDCMZSliceHeader*> fSlices;


};

inline void G4TDCMZSliceMerged::DumpExcessMemory()
{
    for(std::map<G4double,G4TDCMZSliceHeader*>::iterator ite = fSlices.begin();
        ite != fSlices.end(); ++ite) {
        ite->second->DumpExcessMemory();
    }
}


#endif
