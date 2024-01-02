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

#include "globals.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4UIcommand.hh"
#include "G4PhysicalConstants.hh"
#include "G4NistManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4TVolumeConstruction.hh"

#include "G4TPhantomParameterisation.hh"
#include "G4TVolumeBuilderUsingVoxel.hh"

G4TVolumeBuilderUsingVoxel::G4TVolumeBuilderUsingVoxel(){

    //G4cout  << "############################################## "  << __FUNCTION__ << G4endl;

    // this variable value changes, then VolumeConstructor is static and we cannot get the new value which throw an execption then we send it here by set method
    //const G4TVolumeConstruction* TConst = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    //DicomDatafileName = TConst->getDicomOutTextName();

}
G4TVolumeBuilderUsingVoxel::~G4TVolumeBuilderUsingVoxel(){}

void G4TVolumeBuilderUsingVoxel::StartGettingDicomData(){

    const G4TVolumeConstruction* TConst = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    std::ifstream fin(DicomDatafileName);

    if(DicomType == "CT"){

        //InitialisationOfMaterials();
        ReadVoxelHeaderAndMatIDs(fin);
    }
    else if (DicomType == "PET") {
        ReadVoxelActivities(fin);
    }

    fin.close();
}

void G4TVolumeBuilderUsingVoxel::InitialisationOfMaterials()
{
    // Creating elements :
    G4double z, a, density;
    G4String name, symbol;

    G4Element* elC = new G4Element( name = "Carbon",
                                    symbol = "C",
                                    z = 6.0, a = 12.011 * g/mole );
    G4Element* elH = new G4Element( name = "Hydrogen",
                                    symbol = "H",
                                    z = 1.0, a = 1.008  * g/mole );
    G4Element* elN = new G4Element( name = "Nitrogen",
                                    symbol = "N",
                                    z = 7.0, a = 14.007 * g/mole );
    G4Element* elO = new G4Element( name = "Oxygen",
                                    symbol = "O",
                                    z = 8.0, a = 16.00  * g/mole );
    G4Element* elNa = new G4Element( name = "Sodium",
                                     symbol = "Na",
                                     z= 11.0, a = 22.98977* g/mole );
    G4Element* elMg = new G4Element( name = "Magnesium",
                                     symbol = "Mg",
                                     z = 12.0, a = 24.3050* g/mole );
    G4Element* elP = new G4Element( name = "Phosphorus",
                                    symbol = "P",
                                    z = 15.0, a = 30.973976* g/mole );
    G4Element* elS = new G4Element( name = "Sulfur",
                                    symbol = "S",
                                    z = 16.0,a = 32.065* g/mole );
    G4Element* elCl = new G4Element( name = "Chlorine",
                                     symbol = "P",
                                     z = 17.0, a = 35.453* g/mole );
    G4Element* elK = new G4Element( name = "Potassium",
                                    symbol = "P",
                                    z = 19.0, a = 30.0983* g/mole );

    G4Element* elFe = new G4Element( name = "Iron",
                                     symbol = "Fe",
                                     z = 26, a = 56.845* g/mole );

    G4Element* elCa = new G4Element( name="Calcium",
                                     symbol = "Ca",
                                     z = 20.0, a = 40.078* g/mole );

    G4Element* elZn = new G4Element( name = "Zinc",
                                     symbol = "Zn",
                                     z = 30.0,a = 65.382* g/mole );

    // Creating Materials :
    G4int numberofElements;

    // Air
    fAir = new G4Material( "Air",
                           1.290*mg/cm3,
                           numberofElements = 2 );
    fAir->AddElement(elN, 0.7);
    fAir->AddElement(elO, 0.3);


    // Soft tissue (ICRP - NIST)
    G4Material* softTissue = new G4Material ("SoftTissue", 1.00*g/cm3,
                                             numberofElements = 13);
    softTissue->AddElement(elH, 10.4472*perCent);
    softTissue->AddElement(elC, 23.219*perCent);
    softTissue->AddElement(elN, 2.488*perCent);
    softTissue->AddElement(elO, 63.0238*perCent);
    softTissue->AddElement(elNa, 0.113*perCent);
    softTissue->AddElement(elMg, 0.0113*perCent);
    softTissue->AddElement(elP, 0.113*perCent);
    softTissue->AddElement(elS, 0.199*perCent);
    softTissue->AddElement(elCl, 0.134*perCent);
    softTissue->AddElement(elK, 0.199*perCent);
    softTissue->AddElement(elCa, 0.023*perCent);
    softTissue->AddElement(elFe, 0.005*perCent);
    softTissue->AddElement(elZn, 0.003*perCent);

    //  Lung Inhale
    G4Material* lunginhale = new G4Material( "LungInhale",
                                             density = 0.217*g/cm3,
                                             numberofElements = 9);
    lunginhale->AddElement(elH,0.103);
    lunginhale->AddElement(elC,0.105);
    lunginhale->AddElement(elN,0.031);
    lunginhale->AddElement(elO,0.749);
    lunginhale->AddElement(elNa,0.002);
    lunginhale->AddElement(elP,0.002);
    lunginhale->AddElement(elS,0.003);
    lunginhale->AddElement(elCl,0.002);
    lunginhale->AddElement(elK,0.003);

    // Lung exhale
    G4Material* lungexhale = new G4Material( "LungExhale",
                                             density = 0.508*g/cm3,
                                             numberofElements = 9 );
    lungexhale->AddElement(elH,0.103);
    lungexhale->AddElement(elC,0.105);
    lungexhale->AddElement(elN,0.031);
    lungexhale->AddElement(elO,0.749);
    lungexhale->AddElement(elNa,0.002);
    lungexhale->AddElement(elP,0.002);
    lungexhale->AddElement(elS,0.003);
    lungexhale->AddElement(elCl,0.002);
    lungexhale->AddElement(elK,0.003);

    // Adipose tissue
    G4Material* adiposeTissue = new G4Material( "AdiposeTissue",
                                                density = 0.967*g/cm3,
                                                numberofElements = 7);
    adiposeTissue->AddElement(elH,0.114);
    adiposeTissue->AddElement(elC,0.598);
    adiposeTissue->AddElement(elN,0.007);
    adiposeTissue->AddElement(elO,0.278);
    adiposeTissue->AddElement(elNa,0.001);
    adiposeTissue->AddElement(elS,0.001);
    adiposeTissue->AddElement(elCl,0.001);

    // Brain (ICRP - NIST)
    G4Material* brainTissue = new G4Material ("BrainTissue", 1.03 * g/cm3,
                                              numberofElements = 13);
    brainTissue->AddElement(elH, 11.0667*perCent);
    brainTissue->AddElement(elC, 12.542*perCent);
    brainTissue->AddElement(elN, 1.328*perCent);
    brainTissue->AddElement(elO, 73.7723*perCent);
    brainTissue->AddElement(elNa, 0.1840*perCent);
    brainTissue->AddElement(elMg, 0.015*perCent);
    brainTissue->AddElement(elP, 0.356*perCent);
    brainTissue->AddElement(elS, 0.177*perCent);
    brainTissue->AddElement(elCl, 0.236*perCent);
    brainTissue->AddElement(elK, 0.31*perCent);
    brainTissue->AddElement(elCa, 0.009*perCent);
    brainTissue->AddElement(elFe, 0.005*perCent);
    brainTissue->AddElement(elZn, 0.001*perCent);


    // Breast
    G4Material* breast = new G4Material( "Breast",
                                         density = 0.990*g/cm3,
                                         numberofElements = 8 );
    breast->AddElement(elH,0.109);
    breast->AddElement(elC,0.506);
    breast->AddElement(elN,0.023);
    breast->AddElement(elO,0.358);
    breast->AddElement(elNa,0.001);
    breast->AddElement(elP,0.001);
    breast->AddElement(elS,0.001);
    breast->AddElement(elCl,0.001);

    // Spinal Disc
    G4Material* spinalDisc = new G4Material ("SpinalDisc", 1.10 * g/cm3,
                                             numberofElements = 8);
    spinalDisc->AddElement(elH, 9.60*perCent);
    spinalDisc->AddElement(elC, 9.90*perCent);
    spinalDisc->AddElement(elN, 2.20*perCent);
    spinalDisc->AddElement(elO, 74.40*perCent);
    spinalDisc->AddElement(elNa, 0.50*perCent);
    spinalDisc->AddElement(elP, 2.20*perCent);
    spinalDisc->AddElement(elS, 0.90*perCent);
    spinalDisc->AddElement(elCl, 0.30*perCent);


    // Water
    G4Material* water = new G4Material( "Water",
                                        density = 1.0*g/cm3,
                                        numberofElements = 2 );
    water->AddElement(elH,0.112);
    water->AddElement(elO,0.888);

    // Muscle
    G4Material* muscle = new G4Material( "Muscle",
                                         density = 1.061*g/cm3,
                                         numberofElements = 9 );
    muscle->AddElement(elH,0.102);
    muscle->AddElement(elC,0.143);
    muscle->AddElement(elN,0.034);
    muscle->AddElement(elO,0.710);
    muscle->AddElement(elNa,0.001);
    muscle->AddElement(elP,0.002);
    muscle->AddElement(elS,0.003);
    muscle->AddElement(elCl,0.001);
    muscle->AddElement(elK,0.004);

    // Liver
    G4Material* liver = new G4Material( "Liver",
                                        density = 1.071*g/cm3,
                                        numberofElements = 9);
    liver->AddElement(elH,0.102);
    liver->AddElement(elC,0.139);
    liver->AddElement(elN,0.030);
    liver->AddElement(elO,0.716);
    liver->AddElement(elNa,0.002);
    liver->AddElement(elP,0.003);
    liver->AddElement(elS,0.003);
    liver->AddElement(elCl,0.002);
    liver->AddElement(elK,0.003);

    // Tooth Dentin
    G4Material* toothDentin = new G4Material ("ToothDentin", 2.14 * g/cm3,
                                              numberofElements = 10);
    toothDentin->AddElement(elH, 2.67*perCent);
    toothDentin->AddElement(elC, 12.77*perCent);
    toothDentin->AddElement(elN, 4.27*perCent);
    toothDentin->AddElement(elO, 40.40*perCent);
    toothDentin->AddElement(elNa, 0.65*perCent);
    toothDentin->AddElement(elMg, 0.59*perCent);
    toothDentin->AddElement(elP, 11.86*perCent);
    toothDentin->AddElement(elCl, 0.04*perCent);
    toothDentin->AddElement(elCa, 26.74*perCent);
    toothDentin->AddElement(elZn, 0.01*perCent);


    // Trabecular Bone
    G4Material* trabecularBone = new G4Material("TrabecularBone",
                                                density = 1.159*g/cm3,
                                                numberofElements = 12 );
    trabecularBone->AddElement(elH,0.085);
    trabecularBone->AddElement(elC,0.404);
    trabecularBone->AddElement(elN,0.058);
    trabecularBone->AddElement(elO,0.367);
    trabecularBone->AddElement(elNa,0.001);
    trabecularBone->AddElement(elMg,0.001);
    trabecularBone->AddElement(elP,0.034);
    trabecularBone->AddElement(elS,0.002);
    trabecularBone->AddElement(elCl,0.002);
    trabecularBone->AddElement(elK,0.001);
    trabecularBone->AddElement(elCa,0.044);
    trabecularBone->AddElement(elFe,0.001);

    // Trabecular bone used in the DICOM Head

    G4Material* trabecularBone_head = new G4Material ("TrabecularBone_HEAD",
                                                      1.18 * g/cm3,
                                                      numberofElements = 12);
    trabecularBone_head->AddElement(elH, 8.50*perCent);
    trabecularBone_head->AddElement(elC, 40.40*perCent);
    trabecularBone_head->AddElement(elN, 2.80*perCent);
    trabecularBone_head->AddElement(elO, 36.70*perCent);
    trabecularBone_head->AddElement(elNa, 0.10*perCent);
    trabecularBone_head->AddElement(elMg, 0.10*perCent);
    trabecularBone_head->AddElement(elP, 3.40*perCent);
    trabecularBone_head->AddElement(elS, 0.20*perCent);
    trabecularBone_head->AddElement(elCl, 0.20*perCent);
    trabecularBone_head->AddElement(elK, 0.10*perCent);
    trabecularBone_head->AddElement(elCa, 7.40*perCent);
    trabecularBone_head->AddElement(elFe, 0.10*perCent);

    // Dense Bone
    G4Material* denseBone = new G4Material( "DenseBone",
                                            density = 1.575*g/cm3,
                                            numberofElements = 11 );
    denseBone->AddElement(elH,0.056);
    denseBone->AddElement(elC,0.235);
    denseBone->AddElement(elN,0.050);
    denseBone->AddElement(elO,0.434);
    denseBone->AddElement(elNa,0.001);
    denseBone->AddElement(elMg,0.001);
    denseBone->AddElement(elP,0.072);
    denseBone->AddElement(elS,0.003);
    denseBone->AddElement(elCl,0.001);
    denseBone->AddElement(elK,0.001);
    denseBone->AddElement(elCa,0.146);

    // Cortical Bone (ICRP - NIST)
    G4Material* corticalBone = new G4Material ("CorticalBone", 1.85 * g/cm3,
                                               numberofElements = 9);
    corticalBone->AddElement(elH, 4.7234*perCent);
    corticalBone->AddElement(elC, 14.4330*perCent);
    corticalBone->AddElement(elN, 4.199*perCent);
    corticalBone->AddElement(elO, 44.6096*perCent);
    corticalBone->AddElement(elMg, 0.22*perCent);
    corticalBone->AddElement(elP, 10.497*perCent);
    corticalBone->AddElement(elS, 0.315*perCent);
    corticalBone->AddElement(elCa, 20.993*perCent);
    corticalBone->AddElement(elZn, 0.01*perCent);


    // Tooth enamel
    G4Material* toothEnamel = new G4Material ("ToothEnamel", 2.89 * g/cm3,
                                              numberofElements = 10);
    toothEnamel->AddElement(elH, 0.95*perCent);
    toothEnamel->AddElement(elC, 1.11*perCent);
    toothEnamel->AddElement(elN, 0.23*perCent);
    toothEnamel->AddElement(elO,41.66*perCent);
    toothEnamel->AddElement(elNa, 0.79*perCent);
    toothEnamel->AddElement(elMg, 0.23*perCent);
    toothEnamel->AddElement(elP, 18.71*perCent);
    toothEnamel->AddElement(elCl, 0.34*perCent);
    toothEnamel->AddElement(elCa, 35.97*perCent);
    toothEnamel->AddElement(elZn, 0.02*perCent);

//#ifdef DICOM_USE_HEAD
     //----- Put the materials in a vector HEAD PHANTOM
//    fOriginalMaterials.push_back(fAir); //0.00129 g/cm3
//    fOriginalMaterials.push_back(softTissue); // 1.055 g/cm3
//    fOriginalMaterials.push_back(brainTissue); // 1.07 g/cm3
//    fOriginalMaterials.push_back(spinalDisc); // 1.10 g/cm3
//    fOriginalMaterials.push_back(trabecularBone_head); // 1.13 g/cm3
//    fOriginalMaterials.push_back(toothDentin); // 1.66 g/cm3
//    fOriginalMaterials.push_back(corticalBone);  // 1.75 g/cm3
//    fOriginalMaterials.push_back(toothEnamel); // 2.04 g/cm3
//    G4cout << "The materials of the DICOM Head have been used" << G4endl;
//#else

    fOriginalMaterials.push_back(fAir); // rho = 0.00129
    fOriginalMaterials.push_back(lunginhale); // rho = 0.217
    fOriginalMaterials.push_back(lungexhale); // rho = 0.508
    fOriginalMaterials.push_back(adiposeTissue); // rho = 0.967
    fOriginalMaterials.push_back(breast ); // rho = 0.990
    fOriginalMaterials.push_back(water); // rho = 1.018
    fOriginalMaterials.push_back(muscle); // rho = 1.061
    fOriginalMaterials.push_back(liver); // rho = 1.071
    fOriginalMaterials.push_back(trabecularBone); // rho = 1.159 - HEAD PHANTOM
    fOriginalMaterials.push_back(denseBone); // rho = 1.575
    G4cout << "Default materials of the DICOM Extended examples have been used"           << G4endl;

    const G4TVolumeConstruction* TConstruction1 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fOriginalMaterials = TConstruction1->getDcmMaterialsVector();

}

void G4TVolumeBuilderUsingVoxel::ReadVoxelHeaderAndMatIDs(std::ifstream& fin){

    std::vector<G4String> wl;
    G4int nMaterials;
    fin >> nMaterials;
    G4String mateName;
    G4int nmate;
    for( G4int ii = 0; ii < nMaterials; ii++ ){
        fin >> nmate;
        fin >> mateName;
        if( mateName[0] == '"' && mateName[mateName.length()-1] == '"' ) {
            mateName = mateName.substr(1,mateName.length()-2);
        }
        //G4cout << "ReadVoxelsData() reading nmate " << ii << " = " << nmate << " mate " << mateName << G4endl;

        if( ii != nmate ) G4Exception("ReadVoxelsData()", "Wrong argument", FatalErrorInArgument, "Material number should be in increasing order:wrong material number");

        G4Material* mate = 0;
        const G4MaterialTable* matTab = G4Material::GetMaterialTable();
        std::vector<G4Material*>::const_iterator matite;
        for( matite = matTab->begin(); matite != matTab->end(); ++matite ) {
            if( (*matite)->GetName() == mateName ) {
                mate = *matite;
            }
        }
        if( mate == 0 ) {
            mate = G4NistManager::Instance()->FindOrBuildMaterial(mateName);
        }
        if( !mate ) G4Exception("ReadVoxelsData()", "Wrong argument", FatalErrorInArgument, ("Material not found" + mateName).c_str());
        thePhantomMaterialsOriginal[nmate] = mate;
    }

    fin >> fNVoxelX >> fNVoxelY >> fNVoxelZ >> fMinX >> fMaxX >> fMinY >> fMaxY >> fMinZ >> fMaxZ;
    fVoxelHalfDimX = (fMaxX-fMinX)/fNVoxelX/2.;
    fVoxelHalfDimY = (fMaxY-fMinY)/fNVoxelY/2.;
    fVoxelHalfDimZ = (fMaxZ-fMinZ)/fNVoxelZ/2.;

    G4cout  << "CT NX "  << fNVoxelX << G4endl;
    G4cout  << "CT NY "  << fNVoxelY << G4endl;
    G4cout  << "CT NZ "  << fNVoxelZ << G4endl;
    G4cout  << "CT MinX "  << fMinX  << " MaxX " << fMaxX << G4endl;
    G4cout  << "CT MinY "  << fMinY  << " MaxY " << fMaxY << G4endl;
    G4cout  << "CT MinZ "  << fMinZ  << " MaxZ " << fMaxZ << G4endl;

    G4cout  << "CT VoxelHalfDimX "  << fVoxelHalfDimX  <<  G4endl;
    G4cout  << "CT VoxelHalfDimY "  << fVoxelHalfDimY  <<  G4endl;
    G4cout  << "CT VoxelHalfDimZ "  << fVoxelHalfDimZ  <<  G4endl;

    fMateIDs = new size_t[fNVoxelX*fNVoxelY*fNVoxelZ];

    for( G4int iz = 0; iz < fNVoxelZ; iz++ ) {
        for( G4int iy = 0; iy < fNVoxelY; iy++ ) {
            for( G4int ix = 0; ix < fNVoxelX; ix++ ) {
                G4int mateID;
                fin >> mateID;
                G4int nnew = ix + (iy)*fNVoxelX + (iz)*fNVoxelX*fNVoxelY;
                if( mateID < 0 || mateID >= nMaterials ) {
                    G4Exception("ReadVoxelsData()", "Wrong index in phantom file", FatalException, G4String("It should be between 0 and " + G4UIcommand::ConvertToString(nMaterials-1) + ", while it is " + G4UIcommand::ConvertToString(mateID)).c_str());
                }
                fMateIDs[nnew] = mateID;
            }
        }
    }

    ReadVoxelDensities( fin );
}

// called from ReadVoxelHeaderAndMatIDs()
void G4TVolumeBuilderUsingVoxel::ReadVoxelDensities( std::ifstream& fin ) {

    G4String stemp;
    std::map<G4int, std::pair<G4double,G4double> > densiMinMax;
    std::map<G4int, std::pair<G4double,G4double> >::iterator mpite;
    for( size_t ii = 0; ii < thePhantomMaterialsOriginal.size(); ii++ ){
        densiMinMax[ii] = std::pair<G4double,G4double>(DBL_MAX,-DBL_MAX);        
    }

    char* part = getenv( "DICOM_CHANGE_MATERIAL_DENSITY" );
    G4double densityDiff = -1;
    if( part ) densityDiff = G4UIcommand::ConvertToDouble(part);

    std::map<G4int,G4double> densityDiffs;
    for( size_t ii = 0; ii < thePhantomMaterialsOriginal.size(); ii++ ){
        densityDiffs[ii] = densityDiff; //currently all materials with same step
    }

    //  densityDiffs[0] = 0.0001; //air

    //--- Calculate the average material density for each material/density bin
    std::map< std::pair<G4Material*,G4int>, matInfo* > newMateDens;

    //---- Read the material densities
    G4double dens;
    for( G4int iz = 0; iz < fNVoxelZ; iz++ ) {
        for( G4int iy = 0; iy < fNVoxelY; iy++ ) {
            for( G4int ix = 0; ix < fNVoxelX; ix++ ) {
                fin >> dens;
                G4int copyNo = ix + (iy)*fNVoxelX + (iz)*fNVoxelX*fNVoxelY;

                //G4cout  << copyNo << " - densityDiff "  << densityDiff << G4endl;

                if( densityDiff != -1. ) continue;

                //--- store the minimum and maximum density for each material
                mpite = densiMinMax.find( fMateIDs[copyNo] );
                if( dens < (*mpite).second.first ) (*mpite).second.first = dens;
                if( dens > (*mpite).second.second ) (*mpite).second.second = dens;
                //--- Get material from original list of material in file
                G4int mateID = fMateIDs[copyNo];
                std::map<G4int,G4Material*>::const_iterator imite = thePhantomMaterialsOriginal.find(mateID);

                //--- Check if density is equal to the original material density
                if(std::fabs(dens - (*imite).second->GetDensity()/CLHEP::g*CLHEP::cm3 ) < 1.e-9 ) continue;

                //G4cout  << "The density " << dens << " is not equal to any of the original material density " << densityDiff << G4endl;

                //--- Build material name with thePhantomMaterialsOriginal name+density
                G4int densityBin = (G4int(dens/densityDiffs[mateID]));

                G4String mateName = (*imite).second->GetName() + G4UIcommand::ConvertToString(densityBin);
                //--- Look if it is the first voxel with this material/densityBin
                std::pair<G4Material*,G4int> matdens((*imite).second, densityBin );

                std::map< std::pair<G4Material*,G4int>, matInfo* >::iterator mppite = newMateDens.find( matdens );
                if( mppite != newMateDens.end() ){
                    matInfo* mi = (*mppite).second;
                    mi->fSumdens += dens;
                    mi->fNvoxels++;
                    fMateIDs[copyNo] = thePhantomMaterialsOriginal.size()-1 + mi->fId;
                } else {
                    matInfo* mi = new matInfo;
                    mi->fSumdens = dens;
                    mi->fNvoxels = 1;
                    mi->fId = newMateDens.size()+1;
                    newMateDens[matdens] = mi;
                    fMateIDs[copyNo] = thePhantomMaterialsOriginal.size()-1 + mi->fId;
                }
                //G4cout << fMateIDs[copyNo] <<" " ;
            }

        }
    }

    //G4cout << G4endl ;

    if( densityDiff != -1. ) {
        for( mpite = densiMinMax.begin(); mpite != densiMinMax.end(); mpite++ ){
            G4cout << "DicomDetectorConstruction::ReadVoxelDensities" << " ORIG MATERIALS DENSITY " << (*mpite).first << " MIN " << (*mpite).second.first << " MAX " << (*mpite).second.second << G4endl;
        }
    }

    //----- Build the list of phantom materials that go to Parameterisation

    //--- Add original materials
    std::map<G4int,G4Material*>::const_iterator mimite;
    for( mimite = thePhantomMaterialsOriginal.begin(); mimite != thePhantomMaterialsOriginal.end(); mimite++ ){
        fMaterials.push_back( (*mimite).second );
    }
    //
    //---- Build and add new materials
    std::map< std::pair<G4Material*,G4int>, matInfo* >::iterator mppite;
    for( mppite= newMateDens.begin(); mppite != newMateDens.end(); mppite++ ){
        G4double averdens = (*mppite).second->fSumdens/(*mppite).second->fNvoxels;
        G4double saverdens = G4int(1000.001*averdens)/1000.;

#ifndef G4VERBOSE
        G4cout << "G4TVolumeBuilderUsingVoxel::ReadVoxelDensities AVER DENS " << averdens << " -> " << saverdens << " -> " << G4int(1000*averdens) << " "
               << G4int(1000*averdens)/1000 << " " <<  G4int(1000*averdens)/1000. << G4endl;
#endif

        G4String mateName = ((*mppite).first).first->GetName() + "_" + G4UIcommand::ConvertToString(saverdens);
        fMaterials.push_back( BuildMaterialWithChangingDensity( (*mppite).first.first, averdens, mateName ) );
    }

}

// called from ReadVoxelDensities()
G4Material* G4TVolumeBuilderUsingVoxel::BuildMaterialWithChangingDensity( const G4Material* origMate, float density, G4String newMateName ){

    //----- Copy original material, but with new density
    G4int nelem = origMate->GetNumberOfElements();
    G4Material* mate = new G4Material( newMateName, density*g/cm3, nelem, kStateUndefined, STP_Temperature );

    for( G4int ii = 0; ii < nelem; ii++ ){
        G4double frac = origMate->GetFractionVector()[ii];
        G4Element* elem = const_cast<G4Element*>(origMate->GetElement(ii));
        mate->AddElement( elem, frac );
    }

    return mate;
}

void G4TVolumeBuilderUsingVoxel::ReadVoxelActivities( std::ifstream& fin ) {

    //G4cout  << __FUNCTION__ << G4endl;

    G4double MinX, MinY, MinZ, MaxX, MaxY, MaxZ;

    fin >> fPETNVoxelX >> fPETNVoxelY >> fPETNVoxelZ >> MinX >> MaxX >> MinY >> MaxY >> MinZ >> MaxZ;

    fPETVoxelHalfDimX = (MaxX-MinX)/fPETNVoxelX/2.;
    fPETVoxelHalfDimY = (MaxY-MinY)/fPETNVoxelY/2.;
    fPETVoxelHalfDimZ = (MaxZ-MinZ)/fPETNVoxelZ/2.;

    G4cout  << "PET NX "  << fPETNVoxelX << G4endl;
    G4cout  << "PET NY "  << fPETNVoxelY << G4endl;
    G4cout  << "PET NZ "  << fPETNVoxelZ << G4endl;
    G4cout  << "PET MinX "  << MinX  << " MaxX " << MaxX << G4endl;
    G4cout  << "PET MinY "  << MinY  << " MaxY " << MaxY << G4endl;
    G4cout  << "PET MinZ "  << MinZ  << " MaxZ " << MaxZ << G4endl;

    G4cout  << "PET VoxelHalfDimX "  << fPETVoxelHalfDimX  <<  G4endl;
    G4cout  << "PET VoxelHalfDimY "  << fPETVoxelHalfDimY  <<  G4endl;
    G4cout  << "PET VoxelHalfDimZ "  << fPETVoxelHalfDimZ  <<  G4endl;

    fActivities = new double[fPETNVoxelX*fPETNVoxelY*fPETNVoxelZ];
    double Acti;

    G4int nnn = 0;
    for( G4int iz = 0; iz < fPETNVoxelZ; iz++ ) {
        for( G4int iy = 0; iy < fPETNVoxelY; iy++ ) {
            for( G4int ix = 0; ix < fPETNVoxelX; ix++ ) {

                fin >> Acti;
                //G4cout  << "Acti "  << Acti << G4endl;

                //G4int nnew = ix + (iy)*NX + (iz)*NX*NY;
                if( Acti < 0 ) {
                    Acti = 0;
                    G4Exception("ReadVoxelsData()", "Wrong activity value in phantom file", FatalException, G4String("It should be between 0 , while it is " + G4UIcommand::ConvertToString(Acti)).c_str());
                }
                fActivities[nnn ] = Acti;
                nnn++;
            }
        }
    }
    G4cout << nnn << " Readed activity value" << G4endl;
}
