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
#ifndef G4TSTLTOGDML_HH
#define G4TSTLTOGDML_HH

#include "globals.hh"
//#include "G4LogicalVolume.hh"


#include <vector>

class G4TStlToGdml {
public:
    G4TStlToGdml();
    ~G4TStlToGdml();

    void stl_to_gdml(G4String, G4String);
    void GetLogicalVolFromGDMLFile();
    std::vector<double *> FacetNormal;
    std::vector<std::vector<G4String>> FacetVertexesNames;

    std::map<G4String, G4String> VertexName1Name2;
    std::map<G4String, double *> VertexName1Vertex;


    //std::map<G4String, double *> VertexNameValMap;

    double * Normal;
    double * Vertex;

    void getSolidFacetData(G4String);

    void __print_progress_bar__(G4String, float);
    double __get_orientation__(double * , double * , double * , double *);
    double * __vectr_subtr__(double *, double *);
    double * __vector_cross__(double *, double *);
    double __vector_inner__(double *, double *);
    void __get_inputname_base__();

    G4String volumename;
    G4String gdmlname;
    G4String gdmlFilePath ;

    G4String materialname;
    void setStlMaterialName(G4String n){
        materialname = n;
    }
    void setStructureVolumeName(G4String n){
        volumename = n;
    }

    std::map<G4String , G4String> MATERIALS_LIST;
    G4String ScriptsDirectoryPath, MaterialXMLFile, EMAIL, ERROR_CONTACT, HEADER, DOCTYPE, SCHEMA, FOOTER, LUNIT, AUNIT, MATERIALS;
    void DataDefinition();
    void CreateMaterialGDMLFile();

};
#endif //G4TSTLTOGDML_HH
