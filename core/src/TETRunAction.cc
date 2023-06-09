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
// TETRunAction.cc
// \file   MRCP_GEANT4/Internal/src/TETRunAction.cc
// \author Haegin Han
//

#include "TETRunAction.hh"

TETRunAction::TETRunAction(G4String _output)
:fRun(0), numOfEvent(0), runID(0), outputFile(_output)
{}

TETRunAction::~TETRunAction()
{}

G4Run* TETRunAction::GenerateRun()
{
	// generate run
	fRun = new TETRun();
	return fRun;
}


void TETRunAction::BeginOfRunAction(const G4Run* aRun)
{
	// print the progress at the interval of 10%
	numOfEvent=aRun->GetNumberOfEventToBeProcessed();
	G4RunManager::GetRunManager()->SetPrintProgress(int(numOfEvent*0.1));
}

void TETRunAction::EndOfRunAction(const G4Run* aRun)
{
	// print the result only in the Master
	if(!isMaster) return;

	// get the run ID
	runID = aRun->GetRunID();

	// print the run result by G4cout and std::ofstream
	//
	// print by G4cout
	PrintResult(G4cout);

	// print by std::ofstream
	std::ofstream ofs(outputFile.c_str());
	PrintResult(ofs);
	ofs.close();
}

void TETRunAction::PrintResult(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	EDEPMAP edepMap = *fRun->GetEdepMap();

	out.precision(3);

    for ( auto it3 = edepMap.begin(); it3 != edepMap.end(); ++it3  ){
        std::cout << it3->first << " | "  << it3->first  << " | "<< it3->second.first  << " | "  << it3->second.second << std::endl;
    }

/*
    out << G4endl
        << "=====================================================================" << G4endl
        << " Run #" << runID << " / Number of event processed : "<< numOfEvent     << G4endl
        << "=====================================================================" << G4endl
        << "organ ID| "
        << setw(19) << "Organ Mass (g)"
        << setw(19) << "Dose (Gy/source)"
        << setw(19) << "Relative Error" << G4endl;

    for(auto itr : tetData->GetMassMap()){
		G4double meanDose    = edepMap[itr.first].first  / itr.second / numOfEvent;
		G4double squareDoese = edepMap[itr.first].second / (itr.second*itr.second);
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;

		out << setw(8)  << itr.first << "| "
			<< setw(19) << fixed      << itr.second/g;
		out	<< setw(19) << scientific << meanDose/(joule/kg);
		out	<< setw(19) << fixed      << relativeE << G4endl;
    }
  */
    out << "=====================================================================" << G4endl << G4endl;
}

