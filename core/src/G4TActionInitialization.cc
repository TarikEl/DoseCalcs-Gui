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
#include "G4TRunAction.hh"
#include "G4TRunActionForCriticality.hh"
//#include "G4TEventAction.hh"

#include "G4TDetectorEventAction.hh"
#include "G4TDetectorStackingAction.hh"

#include "G4TSteppingAction.hh"
#include "G4TVoxelsSteppingAction.hh"
#include "G4TParamSteppingAction.hh"
#include "G4TNestedParamSteppingAction.hh"
#include "TETSteppingAction.hh"
#include "G4TDetectorSteppingAction.hh"
#include "G4TSteppingActionExternal.hh"

#include "G4TDirectPrimaryGeneratorAction.hh"
#include "G4TDirectVoxelsPrimaryGeneratorAction.hh"
#include "G4TDirectToFilesPrimaryGeneratorAction.hh"
#include "G4TModifiedPrimaryGeneratorAction.hh"
#include "G4TReadPrimaryGeneratorAction.hh"

#ifdef G4MPI_USE
//#include "G4TRunActionMaster.hh"
#endif

//#include "G4TTrackingAction.hh"
//#include "G4TEventAction.hh"
#include "G4TVolumeConstruction.hh"
#include "G4RunManager.hh"
#include "TETRunAction.hh"
#include "G4TEventActionForCriticality.hh"
#include "G4TTackingActionForCriticality.hh"
#include "G4TSteppingActionForCriticality.hh"
#include "G4TNeutronCriticalityPrimaryGeneratorAction.hh"
#include "G4TSteppingActionVoxelizedExternal.hh"
#include "G4TModifiedSteppingAction.hh"

extern G4bool VOXTET_USE;
extern G4int ParamType;
extern G4String GenerateVoxelsResuls;
extern G4String SimulationIntExtNeutDet;

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
    //SetUserAction(new G4TRunActionForCriticality());

}

void G4TActionInitialization::Build() const {

    //G4cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n from function : " << __FUNCTION__<< G4endl;

    // For neutron criticality simulation

    //G4TRunActionForCriticality* runAction = new G4TRunActionForCriticality;
    //SetUserAction(runAction);
    //SetUserAction(new G4TEventActionForCriticality(runAction));
    ////SetUserAction(new G4TTackingActionForCriticality());
    //SetUserAction(new G4TSteppingActionForCriticality(runAction));
    //SetUserAction(new G4TNeutronCriticalityPrimaryGeneratorAction);

    // For detector simulation

    //G4TRunAction* runAction = new G4TRunAction;
    //SetUserAction(runAction);
    //G4TDetectorEventAction* event = new G4TDetectorEventAction();
    //SetUserAction(event);
    //SetUserAction(new G4TDetectorSteppingAction(runAction));
    //SetUserAction(new G4TDetectorStackingAction(runAction));
    //SetUserAction(new G4TDirectPrimaryGeneratorAction);

    // For Normal simulation

    G4TRunAction* runAction = new G4TRunAction;
    SetUserAction(runAction);

    //// Voxelized
    //SetUserAction(new G4TNestedParamSteppingAction(runAction));
    //SetUserAction(new G4TDirectVoxelsPrimaryGeneratorAction);



    //// Stylized
    //SetUserAction(new G4TSteppingAction(runAction));
    //SetUserAction(new G4TDirectPrimaryGeneratorAction);

    //return;

    const G4TVolumeConstruction* TConstruction2 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    if(SourceType=="Voxels"){VOXTET_USE = true;}

    //if(SourceType=="Voxels"){
    //    if( ParamType == 0){
    //        //std::cout << "\n\n\n\n\n\n\n\n\n\n VOXTET_USE " << VOXTET_USE << " SourceType "<< SourceType << " ParamType "<< ParamType << std::endl;
    //        SetUserAction(new G4TParamSteppingAction(runAction));
    //    }else{
    //        SetUserAction(new G4TNestedParamSteppingAction(runAction));
    //    }
    //}
    //else{
    //    if(TConstruction2->getGeometryFileType() == "VoxIDs" || TConstruction2->getGeometryFileType() == "VOXEL" || TConstruction2->getGeometryFileType() == "DICOM"){
    //        SetUserAction(new G4TNestedParamSteppingAction(runAction));
    //    }else{
    //        SetUserAction(new G4TSteppingAction(runAction));
    //    }
    //}

    G4cout << " \n\n\n * SimulationIntExtNeutDet " << SimulationIntExtNeutDet
           << "\n * GenerateVoxelsResuls " << GenerateVoxelsResuls
           << G4endl;

    if(TConstruction2->getGeometryFileType() == "VoxIDs" || TConstruction2->getGeometryFileType() == "VOXEL" || TConstruction2->getGeometryFileType() == "DICOM"){

        if(SimulationIntExtNeutDet == "InternalDosimetry"){

            if(GenerateVoxelsResuls=="yes"){
                G4cout << " \n  -- G4TVoxelsSteppingAction" << G4endl;

                SetUserAction(new G4TVoxelsSteppingAction(runAction));
            }
            else{
                if( ParamType == 0){
                    G4cout << " \n  -- G4TParamSteppingAction" << G4endl;

                    //std::cout << "\n\n\n\n\n\n\n\n\n\n VOXTET_USE " << VOXTET_USE << " SourceType "<< SourceType << " ParamType "<< ParamType << std::endl;
                    SetUserAction(new G4TParamSteppingAction(runAction));
                }else{

                    G4cout << " \n  -- G4TNestedParamSteppingAction" << G4endl;

                    SetUserAction(new G4TNestedParamSteppingAction(runAction));
                }
            }
        }
        else if (SimulationIntExtNeutDet == "ExternalDosimetry"){
            G4cout << " \n  -- G4TSteppingActionVoxelizedExternal" << G4endl;

            SetUserAction(new G4TSteppingActionVoxelizedExternal(runAction));
        }
        else if (SimulationIntExtNeutDet == "MyCPPSteppingAction"){
            G4cout << " \n  -- G4TModifiedSteppingAction" << G4endl;

            SetUserAction(new G4TModifiedSteppingAction(runAction));
        }
    }else{

        if(SimulationIntExtNeutDet == "InternalDosimetry"){
            G4cout << " \n  -- G4TSteppingAction" << G4endl;

            SetUserAction(new G4TSteppingAction(runAction));
        }
        else if(SimulationIntExtNeutDet == "ExternalDosimetry"){
            G4cout << " \n  -- G4TSteppingActionExternal" << G4endl;

            SetUserAction(new G4TSteppingActionExternal(runAction));
        }
        else if (SimulationIntExtNeutDet == "MyCPPSteppingAction"){
            G4cout << " \n  -- G4TModifiedSteppingAction" << G4endl;

            SetUserAction(new G4TModifiedSteppingAction(runAction));
        }

    }

    if( TConstruction2->getUseGeneratedData() == "read"){

        G4cout << " \n  -- G4TReadPrimaryGeneratorAction" << G4endl;
        SetUserAction(new G4TReadPrimaryGeneratorAction);
    }else if(TConstruction2->getUseGeneratedData() == "save"){

        G4cout << " \n  -- G4TDirectToFilesPrimaryGeneratorAction" << G4endl;
        SetUserAction(new G4TDirectToFilesPrimaryGeneratorAction);
    }else if(TConstruction2->getUseGeneratedData() == "MyCPPGenerator"){

        G4cout << " \n  -- G4TModifiedPrimaryGeneratorAction" << G4endl;
        SetUserAction(new G4TModifiedPrimaryGeneratorAction);
    }
    else{
        if(SourceType=="Voxels"){
            G4cout << " \n  -- G4TDirectVoxelsPrimaryGeneratorAction" << G4endl;

            SetUserAction(new G4TDirectVoxelsPrimaryGeneratorAction);
        }
        else if(SourceType=="TET"){
            G4cout << " \n  -- G4TDirectPrimaryGeneratorAction" << G4endl;

            SetUserAction(new G4TDirectPrimaryGeneratorAction);
        }
        else {
            G4cout << " \n  -- G4TDirectPrimaryGeneratorAction" << G4endl;

            SetUserAction(new G4TDirectPrimaryGeneratorAction);
        }
    }

}
