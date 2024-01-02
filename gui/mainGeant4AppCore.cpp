#include "mainGeant4AppCore.h"
#include <stdexcept>

#include "globals.hh"

#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4TransportationManager.hh"
#include "G4HumanPhantomConstruction.hh"
#include "G4HumanPhantomPhysicsList.hh"
#include "G4HumanPhantomActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#endif

#ifdef G4MULTITHREADED
#include "G4MPImanager.hh"
#include "G4MPIsession.hh"
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

mainGeant4AppCore::mainGeant4AppCore()
{

}


int mainGeant4AppCore::Geant4AppMain(G4String InterfaceMode , G4String MacrosStartingFile , G4String mode ){

	//G4cout << " ********************** mainGeant4AppCore::Geant4AppMain() *****************************" << G4endl;

	if( mode == "MPI"){

		G4MTRunManager* runManager = new G4MTRunManager;
		runManager->SetNumberOfThreads(2); // Is equal to 2 by default
		G4cout << " ********************** MPI mode ... *****************************" << G4endl;

		// At first, G4MPImanager/G4MPIsession should be created.
		G4MPImanager* g4MPI= new G4MPImanager(/*argc,argv*/); // if we use the parameter passed from terminal (argv) then the constructor must be empty

		// MPI session (G4MPIsession) instead of G4UIterminal
		G4MPIsession* session= g4MPI-> GetMPIsession();

		// Set mandatory initialization classes
		G4HumanPhantomConstruction* userPhantom = new G4HumanPhantomConstruction();
		runManager->SetUserInitialization(userPhantom);
		runManager->SetUserInitialization(new G4HumanPhantomPhysicsList);

		G4HumanPhantomActionInitialization* actions = new G4HumanPhantomActionInitialization();
		runManager->SetUserInitialization(actions);

		// After user application setting, just start a MPI session. MPIsession treats both interactive and batch modes.

		g4MPI->ExecuteMacroFile(MacrosStartingFile.c_str());
		G4cout << " UI session starts ..." << G4endl;

		// Finally, terminate the program
		delete g4MPI;
		delete runManager;

		return 0;
	}

	else if( mode == "Multi-Threading" ){

		G4MTRunManager* runManager = new G4MTRunManager;
		runManager->SetNumberOfThreads(2); // Is equal to 2 by default
		G4cout << " ********************** Multi-threading mode ... *****************************" << G4endl;

		// Set mandatory initialization classes
		G4HumanPhantomConstruction* userPhantom = new G4HumanPhantomConstruction();
		runManager->SetUserInitialization(userPhantom);
		runManager->SetUserInitialization(new G4HumanPhantomPhysicsList);

		G4HumanPhantomActionInitialization* actions = new G4HumanPhantomActionInitialization();
		runManager->SetUserInitialization(actions);

		G4UImanager* UImanager = G4UImanager::GetUIpointer();

		if(InterfaceMode == "B" || InterfaceMode == ""){

			G4cout << " UI session starts ..." << G4endl;
			UImanager->ExecuteMacroFile(MacrosStartingFile.c_str()); //UImanager->ApplyCommand("/control/loop 1OrgGen.mac ENERGYY 1. 10. 1.");
		}
		else if (InterfaceMode == "G" ){

#ifdef G4VIS_USE
			G4VisExecutive* visManager = new G4VisExecutive;
			visManager->Initialize();
#endif

#ifdef G4UI_USE
			char** argv ;//= {"","",""};
			G4UIExecutive* ui = new G4UIExecutive(3, argv);
			UImanager->ExecuteMacroFile(MacrosStartingFile.c_str());
			ui->SessionStart();
			delete ui;
#endif

		}

		// Finally, terminate the program
		delete runManager;

		return 0;
	}

	else if (mode == "Sequential") {

		G4RunManager* runManager = new G4RunManager;

		G4cout << " in sequential mode ..." << G4endl;
		// Set mandatory initialization classes

		G4HumanPhantomConstruction* userPhantom = new G4HumanPhantomConstruction();
		runManager->SetUserInitialization(userPhantom);

		runManager->SetUserInitialization(new G4HumanPhantomPhysicsList);

		G4HumanPhantomActionInitialization* actions = new G4HumanPhantomActionInitialization();
		runManager->SetUserInitialization(actions);

		G4UImanager* UImanager = G4UImanager::GetUIpointer();

		if(InterfaceMode == "B" || InterfaceMode == ""){

			G4cout << " UI session starts ..." << G4endl;
			UImanager->ExecuteMacroFile(MacrosStartingFile.c_str());

		}else if (InterfaceMode == "G" ){

#ifdef G4VIS_USE
			G4VisExecutive* visManager = new G4VisExecutive;
			visManager->Initialize();
#endif

#ifdef G4UI_USE
			char** argv ; //= {"","",""};
			G4UIExecutive* ui = new G4UIExecutive(3, argv);
			UImanager->ExecuteMacroFile(MacrosStartingFile.c_str());
			ui->SessionStart();
			delete ui;
#endif

			// job termination
#ifdef G4VIS_USE
			delete visManager;
#endif
		}


		delete runManager;
		return 0;
	}
	else{
		return 1;
	}
}
