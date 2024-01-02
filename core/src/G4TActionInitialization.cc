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

#include "G4TActionInitialization.hh"
//#include "G4TTrackingAction.hh"
#include "G4TReadPrimaryGeneratorAction.hh"
#include "G4TRunAction.hh"
//#include "G4TEventAction.hh"
#include "G4TSteppingAction.hh"
#include "G4TDirectPrimaryGeneratorAction.hh"
#include "G4TDirectVoxelsPrimaryGeneratorAction.hh"
#include "G4TDirectToFilesPrimaryGeneratorAction.hh"
#include "G4TParamSteppingAction.hh"
#include "G4TNestedParamSteppingAction.hh"
#include "TETSteppingAction.hh"

#ifdef G4MPI_USE
//#include "G4TRunActionMaster.hh"
#endif

//#include "G4TTrackingAction.hh"
//#include "G4TEventAction.hh"
#include "G4TVolumeConstruction.hh"
#include "G4RunManager.hh"
#include "TETRunAction.hh"

extern G4bool VOXTET_USE;
extern G4int ParamType;

G4TActionInitialization::G4TActionInitialization(): G4VUserActionInitialization(){}

G4TActionInitialization::~G4TActionInitialization(){}

void G4TActionInitialization::BuildForMaster() const
{
    // G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;
    // In MT mode, to be clearer, the RunAction class for the master thread might be different than the one used for the workers.
    // This RunAction will be called before and after starting the  workers.

#ifdef G4MPI_USE
    //SetUserAction(new G4TRunActionMaster);
    //SetUserAction(new G4TRunAction);
#else
    //SetUserAction(new G4TRunAction);
#endif

    //G4TVolumeConstruction* Det;
    //G4TReadPrimaryGeneratorAction* Pri;

    if(SourceType=="Voxels" || SourceType=="TET"){VOXTET_USE = true;}

        SetUserAction(new G4TRunAction());
    }

void G4TActionInitialization::Build() const {

    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;

    G4TRunAction* runAction = new G4TRunAction;
    SetUserAction(runAction);

    const G4TVolumeConstruction* TConstruction2 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    if(SourceType=="Voxels"){VOXTET_USE = true;}

    if(SourceType=="Voxels"){
        if( ParamType == 0){
            //std::cout << "\n\n\n\n\n\n\n\n\n\n VOXTET_USE " << VOXTET_USE << " SourceType "<< SourceType << " ParamType "<< ParamType << std::endl;
            SetUserAction(new G4TParamSteppingAction(runAction));
        }else{
            SetUserAction(new G4TNestedParamSteppingAction(runAction));
        }
    }
    else{
        SetUserAction(new G4TSteppingAction(runAction));
    }    


    if( TConstruction2->getUseGeneratedData() == "read"){
        SetUserAction(new G4TReadPrimaryGeneratorAction);
    }else if(TConstruction2->getUseGeneratedData() == "save"){
        SetUserAction(new G4TDirectToFilesPrimaryGeneratorAction);
    }
    else{
        if(SourceType=="Voxels"){
            //std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n VOXTET_USE : " << VOXTET_USE << std::endl;
            SetUserAction(new G4TDirectVoxelsPrimaryGeneratorAction);
            //SetUserAction(new G4TDirectPrimaryGeneratorAction);
        }
        else if(SourceType=="TET"){
            SetUserAction(new G4TDirectPrimaryGeneratorAction);
        }
        else {
            SetUserAction(new G4TDirectPrimaryGeneratorAction);
        }
    }


    /*

    if(SourceType=="Voxels"){VOXTET_USE = true;}

    if(VOXTET_USE){
        if( ParamType == 0){
            //std::cout << "\n\n\n\n\n\n\n\n\n\n VOXTET_USE " << VOXTET_USE << " SourceType "<< SourceType << " ParamType "<< ParamType << std::endl;
            SetUserAction(new G4TParamSteppingAction(runAction));
        }else{
            SetUserAction(new G4TNestedParamSteppingAction(runAction));
        }
    }else{
        if(SourceType=="TET"){
            SetUserAction(new G4TParamSteppingAction(runAction));
        }
        else{
            SetUserAction(new G4TSteppingAction(runAction));
        }
    }

    if( TConstruction2->getUseGeneratedData() == "read"){
        SetUserAction(new G4TReadPrimaryGeneratorAction);
    }else if(TConstruction2->getUseGeneratedData() == "save"){
        SetUserAction(new G4TDirectToFilesPrimaryGeneratorAction);
    }
    else{
        if(VOXTET_USE){
            //std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n VOXTET_USE : " << VOXTET_USE << std::endl;
            SetUserAction(new G4TDirectVoxelsPrimaryGeneratorAction);
            //SetUserAction(new G4TDirectPrimaryGeneratorAction);
        }else{
            SetUserAction(new G4TDirectPrimaryGeneratorAction);
        }
    }

    */

}
