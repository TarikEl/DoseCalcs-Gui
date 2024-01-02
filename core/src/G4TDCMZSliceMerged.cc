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
// $Id: G4TDCMZSliceMerged.cc 92820 2015-09-17 15:22:14Z gcosmo $
//
/// \file medical/DICOM/src/G4TDCMZSliceMerged.cc
/// \brief Implementation of the G4TDCMZSliceMerged class
//
//
// The code was written by :
//      * Jonathan Madsen : jonathan.madsen@cern.ch (12/18/2012)
//
//  Texas A&M University
//  3133 TAMU, Zachry Building
//  College Station, TX 77843, USA
//
//*******************************************************

#include "G4TDCMZSliceMerged.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"



G4TDCMZSliceMerged::G4TDCMZSliceMerged()
{
    //G4cout<< " Z Slice Merged Object"<< G4endl;

}


G4TDCMZSliceMerged::~G4TDCMZSliceMerged()
{
	fSlices.clear();
}

// called from DicomHandler::CheckFileFormat()
// it call DumpToFile() ,
void G4TDCMZSliceMerged::CheckSlices()
{
    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  void G4TDCMZSliceMerged::CheckSlices()  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    G4cout << "\nG4TDCMZSliceMerged::CheckSlices(), " << fSlices.size() << " Slices...\n" << G4endl;

	G4bool uniformSliceThickness = true;

	if(fSlices.size() > 1) {
		if(fSlices.size() == 2) {

            G4cout << "one " << G4endl;

			G4TDCMZSliceHeader* one = fSlices.begin()->second;
            G4cout << "two " << G4endl;

			G4TDCMZSliceHeader* two = fSlices.end()->second;

            G4cout << "one->GetSliceLocation() " << one->GetSliceLocation() << G4endl;
            G4cout << "two->GetSliceLocation() " << two->GetSliceLocation() << G4endl;


			G4double real_distance = ( two->GetSliceLocation()-one->GetSliceLocation() ) / 2.;
            G4cout << "real_distance " << real_distance << G4endl;

			if(one->GetMaxZ() != two->GetMinZ()) {
				one->SetMaxZ(one->GetSliceLocation()+real_distance);
				two->SetMinZ(two->GetSliceLocation()-real_distance);
                G4cout << "one->GetSliceLocation()+real_distance " << one->GetSliceLocation()+real_distance << G4endl;
                G4cout << "two->GetSliceLocation()-real_distance " << two->GetSliceLocation()-real_distance << G4endl;

				//one->SetMinZ(one->GetSliceLocation()-real_distance);
				//two->SetMaxZ(two->GetSliceLocation()+real_distance);
				if(uniformSliceThickness) {
					one->SetMinZ(one->GetSliceLocation()-real_distance);
					two->SetMaxZ(two->GetSliceLocation()+real_distance);
				}
			}
		} else {

            //G4cout << "fSlices.size() != 2 " << G4endl;

			std::map<G4double,G4TDCMZSliceHeader*>::iterator ite0 = fSlices.begin();
			std::map<G4double,G4TDCMZSliceHeader*>::iterator ite1 = fSlices.begin();
			std::map<G4double,G4TDCMZSliceHeader*>::iterator ite2 = fSlices.begin();
			++ite1;
			++ite2; ++ite2;

			for(; ite2 != fSlices.end(); ++ite0, ++ite1, ++ite2) {
				G4TDCMZSliceHeader* prev = ite0->second;
				G4TDCMZSliceHeader* slice = ite1->second;
				G4TDCMZSliceHeader* next = ite2->second;
				G4double real_max_distance = (next->GetSliceLocation() - slice->GetSliceLocation())/2.;
				G4double real_min_distance = (slice->GetSliceLocation() - prev->GetSliceLocation())/2.;
				G4double real_distance = real_max_distance + real_min_distance;
				G4double stated_distance = slice->GetMaxZ()-slice->GetMinZ();
				if(real_distance != stated_distance) {
					unsigned int sliceNum = std::distance(fSlices.begin(),ite1);
					G4cout << "\tG4TDCMZSliceMerged::CheckSlices - \
                    Slice Distance Error in slice [" << sliceNum 
												   << "]: Real Distance = "
												   << real_distance/mm
												   << " mm, Stated Distance = " << stated_distance/mm << G4endl;
					slice->SetMinZ(slice->GetSliceLocation()-real_min_distance);
					slice->SetMaxZ(slice->GetSliceLocation()+real_max_distance);

					if(ite0 == fSlices.begin()) {
						prev->SetMaxZ(slice->GetMinZ());
						// Using below would make all slice same thickness
						//prev->SetMinZ(prev->GetSliceLocation()-real_min_distance);
						if(uniformSliceThickness) {
							prev->SetMinZ(prev->GetSliceLocation()-real_min_distance); }

					}
					if(static_cast<unsigned int>(std::distance(fSlices.begin(),ite2)+1)== fSlices.size()) {
						next->SetMinZ(slice->GetMaxZ());
						// Using below would make all slice same thickness
						//next->SetMaxZ(next->GetSliceLocation()+real_max_distance);
						if(uniformSliceThickness) {
							next->SetMaxZ(next->GetSliceLocation()+real_max_distance); }
					}
				}
			}
		}
	}

	std::map<G4double,G4TDCMZSliceHeader*>::iterator ite = fSlices.begin();
	for(ite = fSlices.begin(); ite != fSlices.end(); ++ite) {
		ite->second->DumpToFile();
	}

}
