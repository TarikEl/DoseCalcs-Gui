#ifndef MAINGEANT4APPCORE_H
#define MAINGEANT4APPCORE_H

#include "globals.hh"
#include "G4RunManager.hh"

class mainGeant4AppCore
{
public:
    mainGeant4AppCore();

    int Geant4AppMain( G4String nn, G4String dk, G4String oc );
    int NumbOfSim = 0;
    //G4RunManager* runManager;
};


#endif // MAINGEANT4APPCORE_H
