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

#include "G4TStlToGdml.hh"
#include "G4TVolumeConstruction.hh"
#include "G4RunManager.hh"

#if GDML_USE
#include "G4GDMLParser.hh"
#endif

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4TStlToGdml::G4TStlToGdml() {

    G4String appBuildDir(getenv("PWD"));
    //appBuildDir = "/home/tarik/Desktop/WorkSpace/Projects/QtCreator-build-DoseCalcsCore" ;
    ScriptsDirectoryPath = appBuildDir+"/";
    MaterialXMLFile = ScriptsDirectoryPath+"materials.xml";

    DataDefinition();
}
G4TStlToGdml::~G4TStlToGdml() {

}

//call getSolidFacetData(), __print_progress_bar__(), __get_orientation__(), guess_material()
void G4TStlToGdml::stl_to_gdml(G4String solidname, G4String STLSolidPath){

    gdmlname = solidname+".gdml";
    gdmlFilePath = ScriptsDirectoryPath+gdmlname;

    getSolidFacetData(STLSolidPath);

    std::ostringstream vertices;
    vertices << "\n    <define>\n";
    vertices << "      <position name=\""<< "pos0" << "\" lunit=\"" << LUNIT << "\" x=\"" << 0 << "\" y=\"" << 0 << "\" z=\"" << 0 << "\" />\n"  ;
    vertices << "      <rotation name=\""<< "rot0" << "\" aunit=\"" << AUNIT << "\" x=\"" << 0 << "\" y=\"" << 0 << "\" z=\"" << 0 << "\" />\n"  ;

    std::vector<G4String>::iterator itr;
    std::vector<G4String> VNV;
    for ( auto it = VertexName1Vertex.begin(); it != VertexName1Vertex.end(); ++it  ){
        //G4cout << VertexName1Vertex.size() << " " << it->first << " --- " << VertexName1Name2[it->first]
        //<< " " << it->second[0] << " " << it->second[1] << " " << it->second[2] << G4endl;

        itr = std::find (VNV.begin(), VNV.end(), it->first );
        if (itr != VNV.end()){ // found
            VNV.push_back(it->first);
        }else{ //not found
        }

        vertices << "      <position name=\""<< VertexName1Name2[it->first] << "\" lunit=\"" << LUNIT << "\" x=\"" << it->second[0]
                 << "\" y=\"" << it->second[1] << "\" z=\"" << it->second[2] << "\" />\n"  ;
    }

    for (G4int r = 0; r < VNV.size(); r++) {
        //G4cout << VNV[r] << " --- " << G4endl;
    }


    vertices << "\n    </define>\n";

    // write facet as pos triangle tags to gdml file
    int numSortrg = FacetVertexesNames.size();
    std::ostringstream solids;
    solids << "\n    <solids>\n";
    solids << "      <tessellated aunit=\""<< AUNIT << "\" lunit=\"" << LUNIT << "\" name=\"" << solidname << "\" >\n"  ;
    for (G4int r = 0; r < numSortrg; r++) {
        solids << "        <triangular vertex1=\"" << VertexName1Name2[FacetVertexesNames[r][0]];
        solids << "\" vertex2=\"" << VertexName1Name2[FacetVertexesNames[r][1]];
        solids << "\" vertex3=\"" <<  VertexName1Name2[FacetVertexesNames[r][2]] << "\" />\n";
    }
    solids << "\n      </tessellated>\n";
    solids << "\n    </solids>\n";

    std::ostringstream Phys;
    Phys << "\n      <physvol> \n        <volumeref ref=\""
              << volumename << "\"/>\n        <positionref ref=\"pos0\"/>\n        <rotationref ref=\"rot0\"/>\n      </physvol>" ;


    G4String Volnm, matref, solref;
    std::ostringstream structure;
    structure << "\n    <structure>\n       <volume name=\""
              << volumename << "\">\n           <materialref ref=\""
              << materialname << "\"/>\n           <solidref ref=\""
              << solidname << "\"/>\n       </volume> \n"
              //<< Phys.str() << "\n"
              << "    </structure>\n\n";



    G4String worref;
    std::ostringstream setup;
    setup << "\n    <setup name=\"Default\" version=\"1.0\">\n       <world ref=\""
          << volumename << "\"/>\n    </setup>";

    // create output gdml file

    std::ofstream outgdmlfile(gdmlFilePath , std::ios::binary);
    if(outgdmlfile.is_open()){


        G4cout << "\nCreating file " << gdmlFilePath << G4endl ;

        outgdmlfile << HEADER;
        outgdmlfile << DOCTYPE;
        outgdmlfile << SCHEMA;
        outgdmlfile << "&materials; ";
        outgdmlfile << vertices.str();
        outgdmlfile << solids.str();
        outgdmlfile << structure.str();
        outgdmlfile << setup.str();
        outgdmlfile << FOOTER;
        outgdmlfile.close();
    }
}



void G4TStlToGdml::GetLogicalVolFromGDMLFile(){

    std::ifstream ifile(gdmlFilePath.c_str());
    if (ifile) {
#if GDML_USE
        G4GDMLParser parser;

        parser.Read( gdmlFilePath , false );  //false to eliminate the xchema validation because it print a lot of lines of error validation

        //std::cout << " --------------- " << parser.GetPosition(" 17.802475") << std::endl;

        G4LogicalVolume* LV = parser.GetVolume(volumename);
#endif
        //G4cout << " Log Vol data : "<< G4endl;
        //G4cout << LV->GetSolid()->GetName() << G4endl;

        ifile.close();
        //return 0;
    }
    else{
        G4cout << " Unable to open file : " << gdmlFilePath.c_str() << G4endl;
        //return 0;
    }
}


// call __print_progress_bar__(), __get_three_values__()
void G4TStlToGdml::getSolidFacetData(G4String fname){

    std::ifstream fileR(fname.c_str());

    if(fileR.is_open()){

        G4String line , indicator, word;

        std::ostringstream OutStringWords; OutStringWords << "Reading "<< fname.c_str() << " ... "  ;
        //out->showTextMessage(OutStringWords.str() , 5);

        int NInc, VInc;
        int VertexInc = 0;
        std::vector<G4String> TVN;

        while (getline(fileR, line)) {
            std::istringstream LineString(line);
            if(LineString.str().empty()){
                continue;
            }

            LineString >> word ;
            //G4cout << " first word is : " << word << G4endl;

            if(word == "facet"){

                LineString >> word ;
                //G4cout << " is normal word : " << word << G4endl;
                if(word == "normal"){

                    Normal = new double [3];

                    LineString >> word ;Normal[0] =atof(word) ; NInc ++;
                    LineString >> word ;Normal[1] =atof(word) ; NInc ++;
                    LineString >> word ;Normal[2] =atof(word) ; NInc ++;

                    FacetNormal.push_back(Normal);
                    continue;
                }
            }
            else if(word == "vertex"){

                G4String w1, w2, w3;
                Vertex = new double [3];

                LineString >> w1 ;Vertex[0] =atof(w1) ; VInc ++;
                LineString >> w2 ;Vertex[1] =atof(w2) ; VInc ++;
                LineString >> w3 ;Vertex[2] =atof(w3) ;

                G4String nm1 = std::to_string(atof(w1))+std::to_string(atof(w2))+std::to_string(atof(w3));
                G4String nm2 = "v"+std::to_string(VertexInc);

                VertexName1Name2[nm1] = nm2;
                VertexName1Vertex[nm1] = Vertex;

                //G4cout << nm1 << " " << nm2 << " " << Vertex[0] << " " << Vertex[1] << " " << Vertex[3] << G4endl;
                VertexInc++;
                TVN.push_back(nm1); // records names1 of the three vertexes of facet and its cleared

                //for (G4int t = 0 ;t < TVN.size() ; t++) {
                //G4cout << TVN[t] << G4endl;
                //}
            }
            else if(word == "endfacet"){
                //for (G4int t = 0 ;t < TVN.size() ; t++) {
                //G4cout << TVN[t] << G4endl;
                //}
                //G4cout << FacetVertexesNames.size() << " " << TVN[0] << " " << TVN[1] << " " << TVN[2] <<G4endl;

                FacetVertexesNames.push_back(TVN);
                TVN.clear();

            }
        }

        fileR.close();
    }
}


void G4TStlToGdml::__print_progress_bar__(G4String text,float percentage){

    std::ostringstream tx;
    tx <<"\r"<< " " << text <<" "<< percentage ;  // You could try using \r and \b. The first moves the cursor back to the start of the line, the latter back one character. Neither erase the exisiting text, so you've got to supply enough chars to erase all the old ones.
    printf(tx.str().c_str());

    /*    void G4TStlToGdml::__print_progress_bar__(text,percentage){
    sys.stdout.write("\r%s %5.1f%%"%(text,percentage));
    sys.stdout.flush();
    if percentage==100.: print();
*/
}

//call __vectr_subtr__, __vector_cross__, __vector_inner__
double G4TStlToGdml::__get_orientation__(double * norm, double * v1, double * v2, double * v3){
    double * vec1 = __vectr_subtr__(v2, v1);
    double * vec2 = __vectr_subtr__(v3, v2);
    double * norm1 = __vector_cross__(vec1,vec2);
    if(__vector_inner__(norm,norm1)>0){
        return 1;
    }else {
        return -1;
    }
    /*
    void G4TStlToGdml::__get_orientation__(norm, v1, v2, v3){
    vec1 = __vectr_subtr__(v2, v1);
    vec2 = __vectr_subtr__(v3, v2);
    norm1 = __vector_cross__(vec1,vec2);
    return 1 if __vector_inner__(norm,norm1)>0 else -1 ;
    */
}

double * G4TStlToGdml::__vectr_subtr__(double * vector1, double * vector2){
    double * res = new double[3];
    res[0] = vector1[0] - vector2[0];
    res[1] = vector1[1] - vector2[1];
    res[2] = vector1[2] - vector2[2];
    return res;
    //return [x[0]-x[1] for x in zip(vector1,vector2)];
}
double * G4TStlToGdml::__vector_cross__(double * vector1, double * vector2){
    double * res = new double[3];
    res[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
    res[1] = - vector1[0]*vector2[2] + vector1[2]*vector2[0];
    res[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
    return res;

    /*return [      vector1[1]*vector2[2] - vector1[2]*vector2[1],
            - vector1[0]*vector2[2] + vector1[2]*vector2[0],
            vector1[0]*vector2[1] - vector1[1]*vector2[0],   ];
            */

}
double G4TStlToGdml::__vector_inner__(double * vector1,double * vector2){
    double * res = new double[3];
    res[0] = vector1[0] * vector2[0];
    res[1] = vector1[1] * vector2[1];
    res[2] = vector1[2] * vector2[2];
    double r = res[0]+res[1]+res[2];
    return r;
    //return sum([x[0]*x[1] for x in zip(vector1,vector2)]);
}

void G4TStlToGdml::CreateMaterialGDMLFile(){

    const G4TVolumeConstruction* TConstruction2 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    MATERIALS = TConstruction2->getAllMaterialsGDMLText();
    std::ofstream file(MaterialXMLFile , std::ios::binary);
    if(file.is_open()){

        G4cout << "\nCreating file " << MaterialXMLFile << G4endl ;

        file << HEADER;
        file << "\n\n<materials>\n\n";
        file << MATERIALS;
        file << "\n\n</materials>\n\n";
        file.close();
    }

}


// called first
void G4TStlToGdml::DataDefinition(){

    EMAIL         = "imttarikk@gmail.com";
    //ERROR_CONTACT = "\nPlease contact %s for more information"%EMAIL ;

    HEADER = "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    DOCTYPE = "<!DOCTYPE gdml [ <!ENTITY materials SYSTEM \"materials.xml\">]>\n\n";
    SCHEMA = "<gdml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" >\n\n";
    FOOTER = "\n\n</gdml>" ;

    AUNIT  = "deg";
    LUNIT  = "mm";

}
