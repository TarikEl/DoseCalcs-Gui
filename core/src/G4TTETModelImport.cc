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
// G4TTETModelImport.cc
// \file   MRCP_GEANT4/Internal/src/G4TTETModelImport.cc
// \author Haegin Han
//

#include "G4TTETModelImport.hh"
#include "G4RandomTools.hh"

extern G4String PlanesToVisualize;

extern G4String TETNodeDataFile, TETEleDataFile, TETMatDataFile;
extern G4double MinTETPhantom;
extern G4double MaxTETPhantom;
extern std::vector<G4String> RegionsToVisualize;

extern G4String* CopyNumberRegionNameMap;
extern G4float* CopyNumberMassSize;
extern std::map<G4String, G4Colour> RegionNameColour;

extern bool MaterialNameAsRegionName;

extern bool UseVoxelsColour;
extern std::map<G4int,G4String> MaterialIDName;

G4TTETModelImport::G4TTETModelImport()
{

    //G4String eleFile      =  phantomName + ".ele";
    //G4String nodeFile     =  phantomName + ".node";
    //G4String materialFile =  phantomName + ".material";

    G4String nodeFile     =  TETNodeDataFile;
    G4String eleFile      =  TETEleDataFile;

    // read phantom data files (*. ele, *.node)
    DataRead(eleFile, nodeFile);

    // read material file (*.material)
    //G4String materialFile =  TETMatDataFile;
    //MaterialRead(materialFile);

    // read colour data file (colour.dat) if this is interactive mode
    //if(ui) ColourRead();

    // print the summary of phantom information
    //PrintMaterialInfomation();
}

void G4TTETModelImport::DataRead(G4String eleFile, G4String nodeFile)
{
    G4String tempStr;
    G4int tempInt;

    // Read *.node file
    //
    std::ifstream ifpNode;

    //ifpMat.open((phantomDataPath + "/" + nodeFile).c_str());
    ifpNode.open((nodeFile).c_str());
    if(!ifpNode.is_open()) {
        // exception for the case when there is no *.node file
        G4Exception("G4TTETModelImport::DataRead","",FatalErrorInArgument, G4String("      There is no " + nodeFile ).c_str());
    }

#if VERBOSE_USE
    G4cout<<"\n\n========= Read Organs/Tissus/Regions vertex and tetrahedron data from .node and .ele files ====================\n" <<G4endl;
    G4cout << "Opening TETGEN node (vertex points: x y z) file " << nodeFile <<G4endl;
#endif

    G4int numVertex;
    G4double xPos, yPos, zPos;
    G4double xMin(DBL_MAX), yMin(DBL_MAX), zMin(DBL_MAX);
    G4double xMax(DBL_MIN), yMax(DBL_MIN), zMax(DBL_MIN);

    ifpNode >> numVertex >> tempInt >> tempInt >> tempInt;

#if VERBOSE_USE
    G4cout << "Total Number of nodes (vertex points: x y z) to Read: " << numVertex <<  " with plane to visualize " << PlanesToVisualize << " MinTETPhantom " << MinTETPhantom << " MaxTETPhantom " << MaxTETPhantom << "\n" << G4endl;
#endif

    for(G4int i=0; i<numVertex; i++)
    {
        ifpNode >> tempInt >> xPos >> yPos >> zPos;

        // set the unit
        xPos*=cm;
        yPos*=cm;
        zPos*=cm;

        vertexVector.push_back(G4ThreeVector(xPos, yPos, zPos));

        // to get the information of the bounding box of phantom
        if (xPos < xMin) xMin = xPos;
        if (xPos > xMax) xMax = xPos;
        if (yPos < yMin) yMin = yPos;
        if (yPos > yMax) yMax = yPos;
        if (zPos < zMin) zMin = zPos;
        if (zPos > zMax) zMax = zPos;
    }

    // set the variables for the bounding box and phantom size
    boundingBox_Min = G4ThreeVector(xMin,yMin,zMin);
    boundingBox_Max = G4ThreeVector(xMax,yMax,zMax);
    phantomSize = G4ThreeVector(xMax-xMin,yMax-yMin,zMax-zMin);

    ifpNode.close();

    // Read *.ele file
    //
    std::ifstream ifpEle;

    //ifpMat.open((phantomDataPath + "/" + eleFile).c_str());
    ifpEle.open((eleFile).c_str());
    if(!ifpEle.is_open()) {
        // exception for the case when there is no *.ele file
        G4Exception("G4TTETModelImport::DataRead","",FatalErrorInArgument, G4String("      There is no " + eleFile ).c_str());
    }
#if VERBOSE_USE
    G4cout << "Opening TETGEN elements (tetrahedron with node No.) file " << eleFile <<G4endl;
#endif

    ifpEle >> numEle  >> tempInt >> tempInt;

#if VERBOSE_USE
    G4cout << "Total Number of elements (tetrahedron with node No.) to Read: " << numEle << "\n" <<G4endl;
#endif

    G4int SavedEle = 0;

    for(G4int i=0; i<numEle; i++)
    {

        bool isOut = false;
        ifpEle >> tempInt;
        G4int* ele = new G4int[4];
        for(G4int j=0;j<4;j++){
            ifpEle >> tempInt;
            ele[j]=tempInt;

            if(PlanesToVisualize != "all" && PlanesToVisualize != "regions"){  // for each point a vector (xyz)
                if(PlanesToVisualize == "xy" || PlanesToVisualize == "yx"){
                    //G4cout << " MinTETPhantom " << MinTETPhantom << " zPos " << vertexVector[tempInt].getZ() << " MaxTETPhantom " << MaxTETPhantom <<G4endl;
                    if (vertexVector[tempInt].getZ() < MinTETPhantom || vertexVector[tempInt].getZ() > MaxTETPhantom){isOut = true;}
                }
                else if(PlanesToVisualize == "yz" || PlanesToVisualize == "zy"){
                    //G4cout << " MinTETPhantom " << MinTETPhantom << " xPos " << vertexVector[tempInt].getX() << " MaxTETPhantom " << MaxTETPhantom <<G4endl;
                    if (vertexVector[tempInt].getX() < MinTETPhantom || vertexVector[tempInt].getX() > MaxTETPhantom){isOut = true;}
                }
                else if(PlanesToVisualize == "zx" || PlanesToVisualize == "xz"){
                    //G4cout << " MinTETPhantom " << MinTETPhantom << " yPos " << vertexVector[tempInt].getY() << " MaxTETPhantom " << MaxTETPhantom <<G4endl;
                    if (vertexVector[tempInt].getY() < MinTETPhantom || vertexVector[tempInt].getY() > MaxTETPhantom){isOut = true;}
                }
            }
        }
        ifpEle >> tempInt;

        if(isOut == true){continue;}
        //G4cout << " Is IN " << " tempInt " << tempInt <<G4endl;

        // filter for a specific regions
        bool isInn = true;
        if(RegionsToVisualize.size() != 0){
            bool bbbb = false;
            for (int i = 0; i < RegionsToVisualize.size(); i++) {
                //G4cout << " tempInt " << tempInt << " MatIDNameMap[tempInt] " << MaterialIDName[tempInt] << " RegionsToVisualize[i] " << RegionsToVisualize[i] <<G4endl;
                if(MaterialIDName[tempInt] == RegionsToVisualize[i]){
                    bbbb = true;
                    break;
                }
            }
            if(bbbb == true){isInn = true;}else{isInn = false;}
        }
        if(isInn == false){continue;}

        tetVector.push_back(new G4Tet("Tet_Solid",vertexVector[ele[0]],vertexVector[ele[1]],vertexVector[ele[2]],vertexVector[ele[3]])); //for each line (tet) a Solid Tet
        eleVector.push_back(ele);
        materialVector.push_back(tempInt); // for each line (tet) a material id
        //G4cout << " 1 " <<G4endl;

        // calculate the total volume and the number of tetrahedrons for each organ
        std::map<G4int, G4double>::iterator FindIter = volumeMap.find(materialVector[SavedEle]);
        //G4cout << " 2 " <<G4endl;

        if(FindIter!=volumeMap.end()){
            FindIter->second += tetVector[SavedEle]->GetCubicVolume();
            numTetMap[materialVector[SavedEle]]++;
        }
        else {
            volumeMap[materialVector[SavedEle]] = tetVector[SavedEle]->GetCubicVolume();
            numTetMap[materialVector[SavedEle]] = 1;
        }

        SavedEle++;
        //G4cout << " 3 " <<G4endl;

    }

    ifpEle.close();
}
void G4TTETModelImport::GenerateDataforTET(){

    if(MaterialNameAsRegionName == true){

        // fill the maps needed from material commands instead of reding material files like we find in the example
        // where the material commands creation needs to use the materials reading function
        for ( auto it = MatIDNameMap.begin(); it != MatIDNameMap.end(); ++it  ){
            materialMap[it->first] = NameMatMap[it->second];
            densityMap[it->first] = materialMap[it->first]->GetDensity()/(g/cm3); //for each mat id the correcpondent material density

            if(UseVoxelsColour == true){
                colourMap[it->first] = G4Colour( (G4double) G4UniformRand(),(G4double)G4UniformRand(),(G4double)G4UniformRand(), 1.);
                RegionNameColour[it->second] = colourMap[it->first];
            }
            OrganNamesVector.push_back(it->second);
            OrganNameVolumeMap[it->second] = volumeMap[it->first]; // should be in cm3
            OrganNameDensityMap[it->second] = (materialMap[it->first]->GetDensity()/(g/cm3))/1000; // should be in g/cm3
            OrganNameMassMap[it->second] = (OrganNameDensityMap[it->second]*OrganNameVolumeMap[it->second])/1000; // should be in Kg
        }

        // The maps calculated and needed from materials file
        // MatIDNameMap[MatID], densityMap[MatID], colourMap[MatID], materialMap[MatID]

        CopyNumberRegionNameMap = new G4String[tetVector.size()];
        CopyNumberMassSize = new G4float [tetVector.size()];
        //RegionCopyNumberColour = new G4Colour[tetVector.size()]; not used the colors are used from parametrization from ColourMap from this class

        for(int d = 0; d < tetVector.size() ;d++ ){

            CopyNumberMassSize[d] = tetVector[d]->GetCubicVolume()*densityMap[materialVector[d]];
            CopyNumberRegionNameMap[d] = MatIDNameMap[materialVector[d]];

            //std::cout << "Vol "<< tetVector[d]->GetCubicVolume() << " density "<< densityMap[materialVector[d]] << " Region Name of TET CN " << CopyNumberRegionNameMap[d] << " Mass of TET CN " << CopyNumberMassSize[d] << std::endl;
        }
    }
    else{

        // fill the maps needed from material commands instead of reding material files like we find in the example
        // where the material commands creation needs to use the materials reading function

        CopyNumberRegionNameMap = new G4String[tetVector.size()];
        CopyNumberMassSize = new G4float [tetVector.size()];

        for ( auto it = MatIDNameMap.begin(); it != MatIDNameMap.end(); ++it  ){
            materialMap[it->first] = NameMatMap[it->second];
            densityMap[it->first] = materialMap[it->first]->GetDensity()/(g/cm3); //for each mat id the correcpondent material density
            if(UseVoxelsColour == true){
                colourMap[it->first] = G4Colour( (G4double) G4UniformRand(),(G4double)G4UniformRand(),(G4double)G4UniformRand(), 1.);
            }
        }

#if VERBOSE_USE
        //    G4cout<<"\n\n========= TET New Region Data ====================" << std::endl;
        //        std::cout << " DcmRegionsNames.size() "<< DcmRegionsNames.size() << std::endl;
        //        std::cout << " UseDcmRegionsMinDensityMap.size() "<< UseDcmRegionsMinDensityMap.size() << " DcmRegionsMinDensityMap.size() "<< DcmRegionsMinDensityMap.size()<< std::endl;
        //        std::cout << " UseDcmRegionsMaxDensityMap.size() "<< UseDcmRegionsMaxDensityMap.size() << " DcmRegionsMaxDensityMap.size() "<< DcmRegionsMaxDensityMap.size()<< std::endl;
        //        for(int a = 0; a < DcmRegionsNames.size() ;a++ ){
        //            std::cout << " DcmRegionsNames "<< DcmRegionsNames[a] << " UseMinD=" << UseDcmRegionsMinDensityMap[a]  << " MinD=" << DcmRegionsMinDensityMap[a]/(g/cm3)  << " UseMaxD=" << UseDcmRegionsMaxDensityMap[a] << " MaxD="<< DcmRegionsMaxDensityMap[a]/(g/cm3) << std::endl;
        //        }
#endif

        for(int d = 0; d < tetVector.size() ;d++ ){

            G4double density = NameMatMap[MatIDNameMap[materialVector[d]]]->GetDensity()/(g/cm3); //in g/cm3
            G4double volume = tetVector[d]->GetCubicVolume()/1000; // in cm3

            //std::cout << " materialVector[d] "<< materialVector[d] << std::endl;
            //std::cout << " MatIDNameMap[materialVector[d]] "<< MatIDNameMap[materialVector[d]] << std::endl;
            //std::cout << " NameMatMap[MatIDNameMap[materialVector[d]]]->GetDensity() "<< NameMatMap[MatIDNameMap[materialVector[d]]]->GetDensity()/(g/cm3) << std::endl;
            //std::cout << " tetVector[d]->GetCubicVolume() "<< tetVector[d]->GetCubicVolume() << std::endl;
            //std::cout << " tetVector[d]->GetCubicVolume()/1000 "<< tetVector[d]->GetCubicVolume()/1000 << std::endl;

            CopyNumberMassSize[d] = (density*volume)/1000.;

            bool isin = false;
            for(int a = 0; a < DcmRegionsNames.size() ;a++ ){

                //if(UseVoxelMatForSegMap[a] == true && VoxelMatForSegMap[a] == materialVector[d]){
                if(UseDcmRegionsMinDensityMap[a] == true ){
                    if(density >= DcmRegionsMinDensityMap[a]/(g/cm3)){
                        if(UseDcmRegionsMaxDensityMap[a] == true ){
                            if(density < DcmRegionsMaxDensityMap[a]/(g/cm3)){// // for min and max density seg
                                isin = true;
                                CopyNumberRegionNameMap[d] = DcmRegionsNames[a];
                                OrganNameVolumeMap[DcmRegionsNames[a]] += tetVector[d]->GetCubicVolume(); // should be in cm3
                                OrganNameMassMap[DcmRegionsNames[a]] += density*volume; // should be in Kg
                            }else{

                            }
                        }else{ // for just min density seg
                            isin = true;
                            CopyNumberRegionNameMap[d] = DcmRegionsNames[a];
                            OrganNameVolumeMap[DcmRegionsNames[a]] += tetVector[d]->GetCubicVolume(); // should be in cm3
                            OrganNameMassMap[DcmRegionsNames[a]] += density*volume; // should be in Kg
                        }
                    }else{

                    }
                }else{
                    if(UseDcmRegionsMaxDensityMap[a] == true ){
                        if(density < DcmRegionsMaxDensityMap[a]/(g/cm3)){ // for just max density seg
                            isin = true;
                            CopyNumberRegionNameMap[d] = DcmRegionsNames[a];
                            OrganNameVolumeMap[DcmRegionsNames[a]] += tetVector[d]->GetCubicVolume(); // should be in cm3
                            OrganNameMassMap[DcmRegionsNames[a]] += density*volume; // should be in Kg
                        }else{

                        }
                    }else{ // for just mat seg
                        isin = true;
                        CopyNumberRegionNameMap[d] = DcmRegionsNames[a];
                        OrganNameVolumeMap[DcmRegionsNames[a]] += tetVector[d]->GetCubicVolume(); // should be in cm3
                        OrganNameMassMap[DcmRegionsNames[a]] += density*volume; // should be in Kg
                    }
                }
                //}
                //std::cout << "Vol "<< tetVector[d]->GetCubicVolume() << " density "<< densityMap[materialVector[d]] << " Region Name of TET CN " << CopyNumberRegionNameMap[d] << " Mass of TET CN " << CopyNumberMassSize[d] << std::endl;
            }
            if(isin == false){
                CopyNumberRegionNameMap[d] = "TET";
                OrganNameVolumeMap["TET"] += tetVector[d]->GetCubicVolume(); // should be in cm3
                OrganNameMassMap["TET"] += density*volume; // should be in Kg
            }
        }
        for(int a = 0; a < DcmRegionsNames.size() ;a++ ){

            OrganNamesVector.push_back(DcmRegionsNames[a]);
            OrganNameDensityMap[DcmRegionsNames[a]] = (OrganNameMassMap[DcmRegionsNames[a]]*1000)/OrganNameVolumeMap[DcmRegionsNames[a]]; // should be in g/cm3
            RegionNameColour[DcmRegionsNames[a]] = G4Colour( (G4double) G4UniformRand(),(G4double)G4UniformRand(),(G4double)G4UniformRand(), 1.);
        }
        OrganNamesVector.push_back("TET");
        OrganNameDensityMap["TET"] = (OrganNameMassMap["TET"]*1000)/OrganNameVolumeMap["TET"]; // should be in g/cm3
        RegionNameColour["TET"] = G4Colour( (G4double) G4UniformRand(),(G4double)G4UniformRand(),(G4double)G4UniformRand(), 1.);
    }
}
void G4TTETModelImport::PrintMaterialInfomation()
{

    // Print the overal information for each organ
    //

#if VERBOSE_USE
    G4cout<<"\n\n========= Generated phantom data by material ====================\n" <<G4endl;
#endif

    G4cout << "---------------------------------------------------------------------------------------------"<<G4endl;

    G4cout << std::setw(20) << std::left << "Material ID "
           << std::setw(33) << std::left << "Material Name" << " "
           << std::setw(15) << std::left << "Density(g/cm3) "
           << std::setw(15) << std::left << "TET Num"
           << std::setw(15) << std::left << "Phantom %" <<G4endl;

    G4cout << "---------------------------------------------------------------------------------------------"<<G4endl;

    std::map<G4int, G4Material*>::iterator matIter;
    G4cout<<std::setiosflags(std::ios::fixed);
    G4cout.precision(3);
    for(matIter=materialMap.begin(); matIter!=materialMap.end();matIter++)
    {
        G4int idx = matIter->first;

        G4cout << std::setw(20) << std::left << idx                         // organ ID
               << std::setw(33) << std::left << MatIDNameMap[idx] << " "
               << std::setw(15) << std::left << densityMap[idx]
                  << std::setw(15) << std::left << numTetMap[idx]
                     << std::left << numTetMap[idx]<< "/"<< tetVector.size() <<G4endl;
    }

    G4cout << "---------------------------------------------------------------------------------------------"<<G4endl;

    if(MaterialNameAsRegionName == false){

#if VERBOSE_USE
        G4cout<<"\n\n========= Generated phantom data by region ====================\n" <<G4endl;
#endif

        G4cout << "--------------------------------------------------------------------------------------------"<<G4endl;

        G4cout << std::setw(20) << std::left << "Region Name " << " "
               << std::setw(33) << std::left << "Mass(kg)"
               << std::setw(15) << std::left << "Volume(cm3) "
               << std::setw(15) << std::left << "Density(g/cm3) "<<G4endl;

        G4cout << "--------------------------------------------------------------------------------------------"<<G4endl;

        G4cout<<std::setiosflags(std::ios::fixed);
        G4cout.precision(3);
        for(int a = 0; a < OrganNamesVector.size() ;a++ ){
            G4cout << std::setw(20) << std::left << OrganNamesVector[a] << " "                        // organ ID
                   << std::setw(33) << std::left << OrganNameMassMap[OrganNamesVector[a]]
                   << std::setw(15) << std::left << OrganNameVolumeMap[OrganNamesVector[a]]
                   << std::setw(15) << std::left << OrganNameDensityMap[OrganNamesVector[a]]<<G4endl;
        }

        G4cout << "--------------------------------------------------------------------------------------------"<<G4endl;
    }
}


void G4TTETModelImport::MaterialRead(G4String materialFile)
{

    // Read material file (*.material)
    //
    std::ifstream ifpMat;

    //ifpMat.open((phantomDataPath + "/" + materialFile).c_str());
    ifpMat.open((materialFile).c_str());
    if(!ifpMat.is_open()) {
        // exception for the case when there is no *.material file
        G4Exception("G4TTETModelImport::DataRead","",FatalErrorInArgument,
                    G4String("      There is no " + materialFile ).c_str());
    }

    G4cout << "  Opening material file '" << materialFile << "'" <<G4endl;

    char read_data[50];
    char* token;
    G4double zaid;
    G4double fraction;
    G4String MaterialName;
    G4double density;
    while(!ifpMat.eof())
    {
        ifpMat >> read_data;                   //ex) 'C' RBM
        ifpMat >> MaterialName;                //ex)  C 'RBM'
        ifpMat >> read_data;
        density = std::atof(read_data);        //ex) 1.30
        ifpMat >> read_data;                   //ex) g/cm3
        ifpMat >> read_data;
        token = std::strtok(read_data,"m");
        G4int matID = std::atoi(token);         //ex) m'10'
        materialIndex.push_back(matID);    //mat ids vector

        MatIDNameMap[matID] = MaterialName; //for each mat id the correcpondent material name
        densityMap[matID] = density*g/cm3; //for each mat id the correcpondent material density
        if(UseVoxelsColour == true){
            colourMap[matID] = G4Colour( (G4double) G4UniformRand(),(G4double)G4UniformRand(),(G4double)G4UniformRand(), 1.);
        }

        for(G4int i=0 ;  ; i++)
        {
            ifpMat >> read_data;
            if(std::strcmp(read_data, "C")==0 || ifpMat.eof()) break;

            zaid = (G4int)(std::atoi(read_data)/1000);
            ifpMat >> read_data;
            fraction = -1.0 * std::atof(read_data);
            materialIndexMap[matID].push_back(std::make_pair(G4int(zaid), fraction));
        }
    }
    ifpMat.close();

    // Construct materials for each organ
    //
    G4NistManager* nistManager = G4NistManager::Instance();

    for(G4int i=0;i<(G4int)materialIndex.size();i++){
        G4int idx = materialIndex[i];
        G4Material* mat = new G4Material(MatIDNameMap[idx], densityMap[idx], G4int(materialIndexMap[idx].size()), kStateSolid, NTP_Temperature, STP_Pressure);

        std::cout << "\n\n/MaterialData/createMaterial " << mat->GetName() << " " << idx << " " << materialIndexMap[idx].size() << " " << mat->GetDensity()/(g/cm3) << " g/cm3 frac" << " \n" ;
        std::cout << "/MaterialData/addElements " ;

        for(G4int j=0;j<G4int(materialIndexMap[idx].size());j++){
            mat->AddElement(nistManager->FindOrBuildElement(materialIndexMap[idx][j].first), materialIndexMap[idx][j].second);
            std::cout << nistManager->FindOrBuildElement(materialIndexMap[idx][j].first)->GetName() << " " << materialIndexMap[idx][j].second*100 << " " ;
        }
        materialMap[idx]=mat; // for each mat id we have a material
        massMap[idx]=densityMap[idx]*volumeMap[idx];
    }

}
void G4TTETModelImport::ColourRead()
{
    // Read colour data file (colour.dat)
    //
    std::ifstream ifpColour;

    ifpColour.open((phantomDataPath + "/" + "colour.dat").c_str());
    if(!ifpColour.is_open()) {
        // exception for the case when there is no colour.dat file
        G4Exception("G4TTETModelImport::DataRead","",FatalErrorInArgument, G4String("Colour data file was not found ").c_str());
    }

    G4cout << "  Opening colour data file 'colour.dat'" <<G4endl;

    G4int organID;
    G4double red, green, blue, alpha;
    while( ifpColour >> organID >> red >> green >> blue >> alpha ){
        if(UseVoxelsColour == true){
            colourMap[organID] = G4Colour(red, green, blue, alpha);
        }
    }
    ifpColour.close();

}
