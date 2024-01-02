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
// $Id: G4TDCMHandler.cc 101109 2016-11-07 08:14:53Z gcosmo $
//
/// \file medical/G4TDCM/src/G4TDCMHandler.cc
/// \brief Implementation of the G4TDCMHandler class
//
// The code was written by :
//      *Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + University Laval, Quebec (QC) Canada
//*******************************************************
//
//*******************************************************
//
/// G4TDCMHandler.cc :
///        - Handling of DICM images
///         - Reading headers and pixels
///        - Transforming pixel to density and creating *.g4dcm
///          files
//*******************************************************

#include "G4TDCMHandler.hh"
#include "globals.hh"
#include "G4ios.hh"
#include <fstream>

#include <cctype>
#include <cstring>

#include "G4TDCMZSliceHeader.hh"
#include "G4TDCMZSliceMerged.hh"

#include "G4TVolumeConstruction.hh"
#include "G4RunManager.hh"


G4TDCMHandler* G4TDCMHandler::fInstance = 0;


G4TDCMHandler* G4TDCMHandler::Instance(){return fInstance;}


G4TDCMHandler::G4TDCMHandler():DATABUFFSIZE(8192), LINEBUFFSIZE(5020), FILENAMESIZE(512),
    fCompression(0), fNFiles(0), fRows(0), fColumns(0),
    fBitAllocated(0), fMaxPixelValue(0), fMinPixelValue(0),
    fPixelSpacingX(0.), fPixelSpacingY(0.),
    fSliceThickness(0.), fSliceLocation(0.),
    fRescaleIntercept(0), fRescaleSlope(0),
    fLittleEndian(true), fImplicitEndian(false),
    fPixelRepresentation(0), fNbrequali(0),
    fValueDensity(NULL),fValueCT(NULL),fReadCalibration(false),
    fMergedSlices(NULL),
    //DicomDataFile("########################################AltData.dat"),
    //fCt2DensityFile("CT2Density.dat"),
    DicomFilesDirPath("/home/*"),
    DICOMInputFile("DICOM_DATA.dat")
{
    fMergedSlices = new G4TDCMZSliceMerged;
    //getDicomCmdArguments();
    G4cout<< "\n\nG4TDCMHandler() initialzation"<< G4endl;

}



G4TDCMHandler::~G4TDCMHandler(){}


/*
void CheckFileFormat();
    ReadFile(FILE *,char *); reading the header and reading the data pixel :
        template <class Type> void GetValue(char *, Type &);
        void read_undefined_nested(FILE *);
            void read_undefined_item(FILE *);
        G4int read_defined_nested(FILE *, G4int);
        void GetInformation(G4int &, char *);
            template <class Type> void GetValue(char *, Type &);
        G4String fnameG4DCM = G4String(filename2) + ".g4dcm";
        G4int ReadData(FILE *,char *); // note: always use readHeader before readData
            G4float Pixel2density(G4int pixel);
                void ReadCalibration();
            unsigned int GetMaterialIndex( G4float density );
        void StoreData(std::ofstream& foutG4DCM);
            G4float Pixel2density(G4int pixel);
                void ReadCalibration();
            unsigned int GetMaterialIndex( G4float density );
        void StoreData(G4TDCMZSliceHeader* dcmPZSH);
            G4float Pixel2density(G4int pixel);
                void ReadCalibration();
            unsigned int GetMaterialIndex( G4float density );

*/

// called from Constructor
void G4TDCMHandler::getDicomCmdArgumentsForHandler()
{
    const G4TVolumeConstruction* TConstruction0 = static_cast<const G4TVolumeConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    DicomFilesDirPath = TConstruction0->getDicomDataDirPath()+"/*";
    fCompression = TConstruction0->getDicomPixelsCompression();
    ValueCT = TConstruction0->getValueCTVector();
    ValueDensity = TConstruction0->getValueDensityVector();
    fMaterialIndices = TConstruction0->getMaterialIndicesMap();
}


// called from CheckFileFormat() to fill the vector DicomFilePathsVector
#include <glob.h>
std::vector<G4String> G4TDCMHandler::GetDicomFilesNames(G4String DirPath)
{
    glob_t glob_result;
    glob( DirPath , GLOB_TILDE , NULL , &glob_result );
    for ( unsigned int i=0 ; i < glob_result.gl_pathc ; ++i ){
        G4String vvv = glob_result.gl_pathv[i];
        if(vvv.contains(".g4dcm")){

        }else {
            if(vvv.contains(".dcm")){
                DicomFilePathsVector.push_back(glob_result.gl_pathv[i]);
                G4cout << glob_result.gl_pathv[i] << G4endl;
            }else {

            }
        }
    }
    return DicomFilePathsVector;
}



// called from the main function : dcmHandler->CheckFileFormat();
// it call G4TDCMHandler::ReadMaterialIndices() , ReadFile(G4TDCM,inputFile) , fMergedSlices->CheckSlices(),
void G4TDCMHandler::CheckFileFormat()
{
    G4cout<< " Checking .g4dcm files existing with G4TDCMHandler::CheckFileFormat "<< G4endl;

    getDicomCmdArgumentsForHandler();

    std::ifstream testExistence;
    G4bool existAlready = true;

    GetDicomFilesNames(DicomFilesDirPath); // it fill DicomFilePathsVector

    //fCompression = TConstruction0->getDicomPixelsCompression();
    fNFiles = DicomFilePathsVector.size();
    G4cout <<" Number of Files : " << DicomFilePathsVector.size() << G4endl;
    //G4cout << "Compression number " << fCompression << G4endl;

    for(G4int rep = 0; rep < (G4int)DicomFilePathsVector.size(); rep++) {

        //checkData.getline(oneLine,100); // get all data in the line , because is the name of files
        G4String FileName = DicomFilePathsVector[rep];
        FileName += ".g4dcm"; //  to test existing of G4TDCMFile.g4dcm

        //G4cout << "We test existing of file : " << DicomFilePathsVector[rep] << G4endl;

        testExistence.open(FileName);
        if(!testExistence.is_open()) {
            existAlready = false;
            testExistence.clear();
            testExistence.close();
        }
        testExistence.clear();
        testExistence.close();
    }

    // if there is no file .g4dcm (existAlready == false), it files exists, then it pass to int main and contains execution after calling checkFileFormat()
    if( existAlready == false  ) {

        // G4cout << "\nAll the necessary images were not found in processed form " << ", starting with .dcm images\n";

        FILE * dicomFile;
        G4int rflag;

        G4cout << "\nNo .g4dcm files are existed, then for each .dcm file : " << G4endl;
        G4cout << "- ReadFile() " << G4endl;
        G4cout << "--- GetInformation() " << G4endl;
        G4cout << "--- ReadData() " << G4endl;
        G4cout << "--- StoreData()" << G4endl;
        G4cout << "----- Pixel2density() " << G4endl;
        G4cout << "----- GetMaterialIndex() " << G4endl;


        for(G4int rep = 0; rep < (G4int)DicomFilePathsVector.size(); rep++) {

            std::string fileNm = DicomFilePathsVector[rep].c_str();
            char *fileNamePath = &fileNm[0];

            // Open input file and give it to dicom managing methods  :
            dicomFile = std::fopen(fileNamePath,"rb");

            if( dicomFile != 0 ) {

                G4cout<< "\n------------ Begin of file " << fileNamePath << " reading with ReadFile()--------------\n" << G4endl;

                //std::printf("Opening %s and reading using ReadFile(dicom_file_Stream,inputFilenameString) :\n",fileNamePath);

                // the fuction that read the data from the DICOM .dcm file and fill the variables defined of this class...
                // and create the object zHeaderSlice of each DICOM file and save the header data to the file .g4dcm and add it to the map in ZSliceMerged of slicesHeaders
                ReadFile( dicomFile , fileNamePath );

                G4cout<< "\n------------ End of file " << fileNamePath << " reading with ReadFile()----------------\n" << G4endl;

            } else {
                G4cout << "\nError opening file : " << fileNamePath << G4endl;
            }
            rflag = std::fclose(dicomFile) ;
        }

        // from the map filled by the slices created by readData(), Checks the spacing is correct for each sliceHeader object, and verifie that are continuouses, then for each object ZsliceHeader Dumps to file with DumpToFile()
        fMergedSlices->CheckSlices() ;

        if (rflag) return;

    }

    if(fValueDensity) { delete [] fValueDensity; }
    if(fValueCT) { delete [] fValueCT; }
    if(fMergedSlices) { delete fMergedSlices; }
}



// adding code to add activity reading pluging
// this function used to read the G4TDCM file and it's not create the .g4dcm for esach dcm, but send the name (ex 1.g4dcm) to zslice object to use it after in creation with dumptoFile() called from merged::checkSlice()
// called from CheckFileFormat()
// it call :
/* GetValue(buffer, elementLength4),
read_undefined_nested(dicom,elementLength4 ) ,
read_defined_nested(dicom) ,
GetInformation(tagDictionary, data) ,
zslice->AddMaterial(ite->second),
zslice->SetNoVoxelX(fColumns/fCompression);
zslice->SetNoVoxelY(fRows/fCompression);
zslice->SetNoVoxelZ(1);
zslice->SetMinX(-fPixelSpacingX*fColumns/2.);
zslice->SetMaxX(fPixelSpacingX*fColumns/2.);
zslice->SetMinY(-fPixelSpacingY*fRows/2.);
zslice->SetMaxY(fPixelSpacingY*fRows/2.);
zslice->SetMinZ(fSliceLocation-fSliceThickness/2.);
zslice->SetMaxZ(fSliceLocation+fSliceThickness/2.);
ReadData( dicom, filename2 ),
StoreData( zslice );
dcmPZSH->SetSliceLocation(fSliceLocation);
dcmPZSH->AddValue(density);
dcmPZSH->AddMateID(GetMaterialIndex(density));
fMergedSlices->AddZSlice(zslice);
 */
G4int G4TDCMHandler::ReadFile(FILE* dicom, char* filename2)
{

    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMHandler::ReadFile()  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    //G4cout << " if exist, ReadedFile .txt generated is : " << filename2 << G4endl;
    G4int returnvalue = 0; size_t rflag = 0;
    char * buffer = new char[LINEBUFFSIZE];

    fImplicitEndian = false;
    fLittleEndian = true;

    // fread function reads the data from a position, if we call it to the second time, it will begin from the last position when it terminated befor
    rflag = std::fread( buffer, 1, 128, dicom ); // The size "1 byte" 128 objects(bytes), Read block of data from stream dicom and save it to buffer

    // are not important. Reads the "DICOM" letters
    rflag = std::fread( buffer, 1, 4, dicom );
    // if there is no preamble, the FILE pointer is rewinded.
    if(std::strncmp("DICM", buffer, 4) != 0) {
        std::fseek(dicom, 0, SEEK_SET);  // dicom  FILE object that identifies the stream . 0 Number of bytes to offset from origin. SEEK_SET(mean Beginning of file) Position used as reference origin for the offset
        //after fseek , the fread will begin from SEEK_SET that mean from the begining
        fImplicitEndian = true;
    }

    short readGroupId;    // identify the kind of input data
    short readElementId;  // identify a particular type information
    short elementLength2; // deal with element length (value length of element) in 2 bytes
    unsigned long elementLength4; // deal with element length (value length of element) in 4 bytes

    char * data = new char[DATABUFFSIZE];

    int ii = 0;

    // Read information up to the pixel data
    while(true) {

        //Reading groups and elements :
        readGroupId = 0;
        readElementId = 0;

        // group ID
        rflag = std::fread(buffer, 2, 1, dicom);
        GetValue(buffer, readGroupId); // we give it the value readed (buffer) and it we give us the readGroupId

        // element ID
        rflag = std::fread(buffer, 2, 1, dicom);
        GetValue(buffer, readElementId);

        // Creating a tag to be identified afterward
        G4int tagDictionary = readGroupId*0x10000 + readElementId;


        //G4cout<< "ii : " << ii  << "  --> buffer of readElementId : "<< buffer << "  readGroupId : "  << readGroupId <<  "    readElementId : " << readElementId << " ---> tagDictionary : " <<tagDictionary << G4endl;
        ii = ii +1 ;

        // beginning of the pixels
        if(tagDictionary == 0x7FE00010) {

            // Folling 2 fread's are modifications to original DICOM example
            rflag = std::fread(buffer,2,1,dicom);   // Reserved 2 bytes
            // (not used for pixels)
            rflag = std::fread(buffer,4,1,dicom);   // Element Length
            // (not used for pixels)

            //G4cout<< " we are in Reading the VR and VL tagDictionary == 0x7FE00010    " << " -->  buffer Element Length : " << buffer << G4endl;

            break;      // Exit to ReadImageData()
        }

        // VR or element length , each VR has it own element_length
        rflag = std::fread(buffer,2,1,dicom);
        GetValue(buffer, elementLength2);

        // If value represent ation (VR) is OB, OW, SQ, UN, added OF and UT
        // the next length (of Value length) is 32 bits
        if((elementLength2 == 0x424f ||     	// "OB"
            elementLength2 == 0x574f ||     // "OW"
            elementLength2 == 0x464f ||     // "OF"
            elementLength2 == 0x5455 ||     // "UT"
            elementLength2 == 0x5153 ||     // "SQ"
            elementLength2 == 0x4e55) &&    // "UN"
                !fImplicitEndian ) {            // explicit VR

            rflag = std::fread(buffer, 2, 1, dicom); // to Skip 2 reserved bytes(because we dont use this buffer)

            // element length
            rflag = std::fread(buffer, 4, 1, dicom);
            GetValue(buffer, elementLength4);

            if(elementLength2 == 0x5153)
            {
                if(elementLength4 == 0xFFFFFFFF)
                {
                    read_undefined_nested( dicom );
                    elementLength4=0;
                }  else{
                    if(read_defined_nested( dicom, elementLength4 )==0){
                        G4Exception("DicomHandler::ReadFile", "DICOM001", FatalException, "Function read_defined_nested() failed!");
                    }
                }
            } else  { // Reading the information with data

                rflag = std::fread(data, elementLength4,1,dicom);
            }

        }  else { //  explicit VR but different than previous ones  VR 2 bytes and next(Value length) and 2 bytes or 16 bits

            if(!fImplicitEndian || readGroupId == 2) {

                //G4cout << "Reading  DICOM files with Explicit VR"<< G4endl;
                // element length (2 bytes)
                rflag = std::fread(buffer, 2, 1, dicom);
                GetValue(buffer, elementLength2);
                elementLength4 = elementLength2;

                rflag = std::fread(data, elementLength4, 1, dicom);

            } else {                                  // Implicit VR

                //G4cout << "Reading  DICOM files with Implicit VR"<< G4endl;

                // element length (4 bytes)
                if(std::fseek(dicom, -2, SEEK_CUR) != 0) {
                    G4Exception("DicomHandler::ReadFile", "DICOM001", FatalException, "fseek failed");
                }

                rflag = std::fread(buffer, 4, 1, dicom);
                GetValue(buffer, elementLength4);

                //G4cout <<  std::hex<< elementLength4 << G4endl;

                if(elementLength4 == 0xFFFFFFFF)
                {
                    read_undefined_nested(dicom);
                    elementLength4=0;
                }  else{
                    rflag = std::fread(data, elementLength4, 1, dicom);
                }

            }
        }

        // NULL termination
        data[elementLength4] = '\0';

        // analyzing information , G4int tagDictionary = readGroupId*0x10000 + readElementId; that mean the last tag in DICOM because in end of while loop
        GetInformation(tagDictionary, data);
    }

    // filename2 is "1" or "2" or ... that exist in directory and declared in cmakelist
    G4String fnameG4DCM = G4String(filename2) + ".g4dcm";

    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@ G4String fnameG4DCM = G4String(filename2) + .g4dcm; @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    // we will fill the object ZHeadSlice of this DICOM file .dcm with the header data values readed and filled above

    // Perform functions originally written straight to file, on le donne le fichier txt qui va le remplir par les données de cette slice en utilisant la fonction readData and storData
    G4TDCMZSliceHeader* zslice = new G4TDCMZSliceHeader(fnameG4DCM);

    // on va remplis le zslice objet par les premières données d'abord

    std::map<G4float,G4String>::const_iterator ite;
    for( ite = fMaterialIndices.begin(); ite != fMaterialIndices.end(); ++ite){
        zslice->AddMaterial(ite->second);
    }

    zslice->SetNoVoxelX(fColumns/fCompression);
    zslice->SetNoVoxelY(fRows/fCompression);
    zslice->SetNoVoxelZ(1);

    zslice->SetMinX(-fPixelSpacingX*fColumns/2.);
    zslice->SetMaxX(fPixelSpacingX*fColumns/2.);

    zslice->SetMinY(-fPixelSpacingY*fRows/2.);
    zslice->SetMaxY(fPixelSpacingY*fRows/2.);

    zslice->SetMinZ(fSliceLocation-fSliceThickness/2.);
    zslice->SetMaxZ(fSliceLocation+fSliceThickness/2.);

    //===

    // to read the pixels data values and save the density values in array tab[] and fill the bianry file .g4dcmb with the same data that will fill the .g4dcm after checkslices() consistency from mergedSlices and then dumptoFile()

    ReadData( dicom, filename2 );

    // DEPRECIATED
    //StoreData( foutG4DCM );
    //foutG4DCM.close();

    // on la donne l'objet zslice comme paramètre , et ConstructVoxeDcmGeometry qui remplis les variables
    // de cet objet , et enfin cet objet prend ces valeurs de ces variable et les ecris dans le fichier txt de cette slice
    // store the data of densiry and material index of each pixel (or voxel with compression) in this zslice object
    StoreData( zslice );

    // Dumped 2 file after G4TDCMZSliceMerged has checked for consistency(in z and x and y , and ... )
    // zslice->DumpToFile();

    // add the slice object from this .dcm to the map of ZSliceHeaders to use this

    fMergedSlices->AddZSlice(zslice);

    delete [] buffer;
    delete [] data;

    if (rflag) return returnvalue;
    return returnvalue;

}



// adding code to add activity reading pluging
// get DICOM data needed for this tasks, and save it's values to the class variables
// using tagDictionary sended by readFile function,is represented by 4 byte
// called from ReadFile()
// it call GetValue()
void G4TDCMHandler::GetInformation(G4int & tagDictionary, char * data)
{

    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMHandler::GetInformation  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    if(tagDictionary == 0x00280010 ) { // Number of Rows in image DICOM and each i row and i column is a pixel
        GetValue(data, fRows);
        std::printf("[0x00280010] Rows -> %i\n",fRows);

    } else if(tagDictionary == 0x00280011 ) { // Number of fColumns in image DICOM and each i row and i column is a pixel
        GetValue(data, fColumns);
        std::printf("[0x00280011] Columns -> %i\n",fColumns);

    } else if(tagDictionary == 0x00280102 ) { // High bits  ( not used )
        short highBits;
        GetValue(data, highBits);
        std::printf("[0x00280102] High bits -> %i\n",highBits);

    } else if(tagDictionary == 0x00280100 ) { // number of Bits where the pixel value is allocated
        GetValue(data, fBitAllocated);
        std::printf("[0x00280100] Bits allocated -> %i\n", fBitAllocated);

    } else if(tagDictionary == 0x00280101 ) { //  Bits stored ( not used ) for the value in the bytes of pixels(fBitAllocated)
        short bitStored;
        GetValue(data, bitStored);
        std::printf("[0x00280101] Bits stored -> %i\n",bitStored);

    } else if(tagDictionary == 0x00280106 ) { //  Min. pixel value for image
        GetValue(data, fMinPixelValue);
        std::printf("[0x00280106] Min. pixel value -> %i\n", fMinPixelValue);

    } else if(tagDictionary == 0x00280107 ) { //  Max. pixel value for image
        GetValue(data, fMaxPixelValue);
        std::printf("[0x00280107] Max. pixel value -> %i\n", fMaxPixelValue);

    } else if(tagDictionary == 0x00281053) { //  Rescale slope parameter used to convert value of pixel to density
        fRescaleSlope = atoi(data);
        std::printf("[0x00281053] Rescale Slope -> %d\n", fRescaleSlope);

    } else if(tagDictionary == 0x00281052 ) { // Rescalse intercept parameter used to convert value of pixel to density
        fRescaleIntercept = atoi(data);
        std::printf("[0x00281052] Rescale Intercept -> %d\n",
                    fRescaleIntercept );

    } else if(tagDictionary == 0x00280103 ) {
        //  Pixel representation ( functions not design to read signed bits )
        fPixelRepresentation = atoi(data); // 0: unsigned  1: signed
        std::printf("[0x00280103] Pixel Representation -> %i\n", fPixelRepresentation);
        if(fPixelRepresentation == 1 ) {
            std::printf("### PIXEL REPRESENTATION = 1, BITS ARE SIGNED, ");
            std::printf("DICOM READING SCAN FOR UNSIGNED VALUE, POSSIBLE ");
            std::printf("ERROR !!!!!! -> \n");
        }

    } else if(tagDictionary == 0x00080006 ) { //  Modality
        std::printf("[0x00080006] Modality -> %s\n", data);

    } else if(tagDictionary == 0x00080070 ) { //  Manufacturer
        std::printf("[0x00080070] Manufacturer -> %s\n", data);

    } else if(tagDictionary == 0x00080080 ) { //  Institution Name
        std::printf("[0x00080080] Institution Name -> %s\n", data);

    } else if(tagDictionary == 0x00080081 ) { //  Institution Address
        std::printf("[0x00080081] Institution Address -> %s\n", data);

    } else if(tagDictionary == 0x00081040 ) { //  Institution Department Name
        std::printf("[0x00081040] Institution Department Name -> %s\n", data);

    } else if(tagDictionary == 0x00081090 ) { //  Manufacturer's Model Name
        std::printf("[0x00081090] Manufacturer's Model Name -> %s\n", data);

    } else if(tagDictionary == 0x00181000 ) { //  Device Serial Number
        std::printf("[0x00181000] Device Serial Number -> %s\n", data);

    } else if(tagDictionary == 0x00080008 ) { //  Image type ( not used )
        std::printf("[0x00080008] Image Types -> %s\n", data);

    } else if(tagDictionary == 0x00283000 ) { //Modality LUT Sequence(not used)
        std::printf("[0x00283000] Modality LUT Sequence SQ 1 -> %s\n", data);

    } else if(tagDictionary == 0x00283002 ) { // LUT Descriptor ( not used )
        std::printf("[0x00283002] LUT Descriptor US or SS 3 -> %s\n", data);

    } else if(tagDictionary == 0x00283003 ) { // LUT Explanation ( not used )
        std::printf("[0x00283003] LUT Explanation LO 1 -> %s\n", data);

    } else if(tagDictionary == 0x00283004 ) { // Modality LUT ( not used )
        std::printf("[0x00283004] Modality LUT Type LO 1 -> %s\n", data);

    } else if(tagDictionary == 0x00283006 ) { // LUT Data ( not used )
        std::printf("[0x00283006] LUT Data US or SS -> %s\n", data);

    } else if(tagDictionary == 0x00283010 ) { // VOI LUT ( not used )
        std::printf("[0x00283010] VOI LUT Sequence SQ 1 -> %s\n", data);

    } else if(tagDictionary == 0x00280120 ) { // Pixel Padding Value (not used)
        std::printf("[0x00280120] Pixel Padding Value US or SS 1 -> %s\n", data);

    } else if(tagDictionary == 0x00280030 ) { // Pixel Spacing in X and in Y (in Z we dont need it for images but for volumes ) the longueur x is ncol*fPixelSpacingX ,the same for y
        G4String datas(data);
        int iss = datas.find('\\');
        fPixelSpacingX = atof( datas.substr(0,iss).c_str() );
        fPixelSpacingY = atof( datas.substr(iss+2,datas.length()).c_str() );

    } else if(tagDictionary == 0x00200037 ) { // Image Orientation ( not used )
        std::printf("[0x00200037] Image Orientation (Phantom) -> %s\n", data);

    } else if(tagDictionary == 0x00200032 ) { // Image Position ( not used )
        std::printf("[0x00200032] Image Position (Phantom,mm) -> %s\n", data);

    } else if(tagDictionary == 0x00180050 ) { // Slice Thickness
        fSliceThickness = atof(data);
        std::printf("[0x00180050] Slice Thickness (mm) -> %f\n", fSliceThickness);

    } else if(tagDictionary == 0x00201041 ) { // Slice Location
        fSliceLocation = atof(data);
        std::printf("[0x00201041] Slice Location -> %f\n", fSliceLocation);

    } else if(tagDictionary == 0x00280004 ) { // Photometric Interpretation
        // ( not used )
        std::printf("[0x00280004] Photometric Interpretation -> %s\n", data);

    } else if(tagDictionary == 0x00020010) { // Endian
        if(strcmp(data, "1.2.840.10008.1.2") == 0)
            fImplicitEndian = true;
        else if(strncmp(data, "1.2.840.10008.1.2.2", 19) == 0)
            fLittleEndian = false;
        //else 1.2.840..10008.1.2.1 (explicit little endian)

        std::printf("[0x00020010] Endian -> %s\n", data);
    }

    // others
    else {
        //std::printf("[0x%x] -> %s\n", tagDictionary, data);
        ;
    }

}


// cette fonction lit les valeur de pixel de chaque colone dans tout row et remplit le fTab[i][j] par les valeurs de densité après la
// par contre getInformation() read the value of othe tags but this read the pixel data with the last groupeId and elementId cause readFile break in this two values and pass to readData function
// conversion des valeurs de pixel et fRescaleSlope....
// called from ReadFile()
// it call GetValue() , Pixel2density() , GetMaterialIndex()
G4int G4TDCMHandler::ReadData(FILE *dicom,char * filename2)
{
    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function DicomHandler::ReadData  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    G4int returnvalue = 0; size_t rflag = 0;

    //  READING THE PIXELS values and and convert it to a density and save it to tab[][] then to material index and save it to the array  :

    G4int w = 0;

    // fTab contains the mean value of each point(pixel) in row-column in DICOM file image
    fTab = new G4int*[fRows];
    for ( G4int i = 0; i < fRows; i ++ ) {
        fTab[i] = new G4int[fColumns];
    }

    if(fBitAllocated == 8) { // Case 8 bits :

        std::printf("@@@ Error! Picture != 16 bits...\n");
        std::printf("@@@ Error! Picture != 16 bits...\n");
        std::printf("@@@ Error! Picture != 16 bits...\n");

        unsigned char ch = 0;

        for(G4int j = 0; j < fRows; j++) {
            for(G4int i = 0; i < fColumns; i++) {
                w++;
                rflag = std::fread( &ch, 1, 1, dicom);
                fTab[j][i] = ch*fRescaleSlope + fRescaleIntercept;
            }
        }
        returnvalue = 1;

    } else { //  from 12 to 16 bits :
        char sbuff[2];
        short pixel;
        for( G4int j = 0; j < fRows; j++) {
            for( G4int i = 0; i < fColumns; i++) {
                w++;
                rflag = std::fread(sbuff, 2, 1, dicom);
                GetValue(sbuff, pixel);
                fTab[j][i] = pixel*fRescaleSlope + fRescaleIntercept;
            }
        }
    }


    // Creation of .g4 files wich contains averaged density data

    /*char * nameProcessed = new char[FILENAMESIZE];
    FILE* fileOut;

    //--- create file filename2.g4dcmb, The same information is also used to fill a file in binary format, that contains the same information as the text format,
    // its name ends in .g4dcmb, instead of .g4dcm .
    std::sprintf(nameProcessed,"%s.g4dcmb",filename2);
    fileOut = std::fopen(nameProcessed,"w+b");
    std::printf("### Writing of %s ###\n",nameProcessed);

    //--- Write number of materials

    unsigned int nMate = fMaterialIndices.size();
    rflag = std::fwrite(&nMate, sizeof(unsigned int), 1, fileOut);

    //--- Write materials
    std::map<G4float,G4String>::const_iterator ite;
    for( ite = fMaterialIndices.begin(); ite != fMaterialIndices.end(); ite++ ){
        G4String mateName = (*ite).second;
        for( G4int ii = (*ite).second.length(); ii < 40; ii++ ) {
            mateName += " ";
        }         //mateName = const_cast<char*>(((*ite).second).c_str());

        const char* mateNameC = mateName.c_str();
        rflag = std::fwrite(mateNameC, sizeof(char),40, fileOut);
    }

    unsigned int fRowsC = fRows/fCompression;
    unsigned int fColumnsC = fColumns/fCompression;
    unsigned int planesC = 1;

    G4float pixelLocationXM = -fPixelSpacingX*fColumns/2.;
    G4float pixelLocationXP = fPixelSpacingX*fColumns/2.;
    G4float pixelLocationYM = -fPixelSpacingY*fRows/2.;
    G4float pixelLocationYP = fPixelSpacingY*fRows/2.;
    G4float fSliceLocationZM = fSliceLocation-fSliceThickness/2.;
    G4float fSliceLocationZP = fSliceLocation+fSliceThickness/2.;

    //--- Write number of voxels (assume only one voxel in Z)
    rflag = std::fwrite(&fRowsC, sizeof(unsigned int), 1, fileOut);
    rflag = std::fwrite(&fColumnsC, sizeof(unsigned int), 1, fileOut);
    rflag = std::fwrite(&planesC, sizeof(unsigned int), 1, fileOut);

    //--- Write minimum and maximum extensions
    rflag = std::fwrite(&pixelLocationXM, sizeof(G4float), 1, fileOut);
    rflag = std::fwrite(&pixelLocationXP, sizeof(G4float), 1, fileOut);
    rflag = std::fwrite(&pixelLocationYM, sizeof(G4float), 1, fileOut);
    rflag = std::fwrite(&pixelLocationYP, sizeof(G4float), 1, fileOut);
    rflag = std::fwrite(&fSliceLocationZM, sizeof(G4float), 1, fileOut);
    rflag = std::fwrite(&fSliceLocationZP, sizeof(G4float), 1, fileOut);
    // rflag = std::fwrite(&fCompression, sizeof(unsigned int), 1, fileOut);

    G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@ from G4TDCMHandler::ReadData the data are written to the .g4dcmbbb file and that is it :  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    std::printf("%8i   %8i\n",fRows,fColumns);
    std::printf("%8f   %8f\n",fPixelSpacingX,fPixelSpacingY);
    std::printf("%8f\n", fSliceThickness);
    std::printf("%8f\n", fSliceLocation);
    std::printf("%8i\n", fCompression);

    G4int compSize = fCompression;
    G4int mean;
    G4float density;
    G4bool overflow = false;


    // for writting the material index values

    //----- Write index of material for each pixel
    if(compSize == 1) { // no fCompression: each pixel has a density value)
        for( G4int ww = 0; ww < fRows; ww++) {
            for( G4int xx = 0; xx < fColumns; xx++) {
                mean = fTab[ww][xx];
                density = Pixel2density(mean);
                unsigned int mateID = GetMaterialIndex( density );
                rflag = std::fwrite(&mateID, sizeof(unsigned int), 1, fileOut);
            }
        }

    } else {
        // density value is the average of a square region of
        // fCompression*fCompression pixels
        for(G4int ww = 0; ww < fRows ;ww += compSize ) {
            for(G4int xx = 0; xx < fColumns ;xx +=compSize ) {
                overflow = false;
                mean = 0;
                for(int sumx = 0; sumx < compSize; sumx++) {
                    for(int sumy = 0; sumy < compSize; sumy++) {
                        if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
                        mean += fTab[ww+sumy][xx+sumx];
                    }
                    if(overflow) break;
                }
                mean /= compSize*compSize;

                if(!overflow) {
                    density = Pixel2density(mean);
                    unsigned int mateID = GetMaterialIndex( density );
                    rflag = std::fwrite(&mateID, sizeof(unsigned int), 1, fileOut);
                }
            }

        }
    }


    // for writting the density values

    //----- Write density for each pixel
    if(compSize == 1) { // no fCompression: each pixel has a density value)
        for( G4int ww = 0; ww < fRows; ww++) {
            for( G4int xx = 0; xx < fColumns; xx++) {
                mean = fTab[ww][xx];
                density = Pixel2density(mean);
                rflag = std::fwrite(&density, sizeof(G4float), 1, fileOut);
            }
        }

    } else {
        // density value is the average of a square region of
        // fCompression*fCompression pixels
        for(G4int ww = 0; ww < fRows ;ww += compSize ) {
            for(G4int xx = 0; xx < fColumns ;xx +=compSize ) {
                overflow = false;
                mean = 0;
                for(int sumx = 0; sumx < compSize; sumx++) {
                    for(int sumy = 0; sumy < compSize; sumy++) {
                        if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
                        mean += fTab[ww+sumy][xx+sumx];
                    }
                    if(overflow) break;
                }
                mean /= compSize*compSize;

                if(!overflow) {
                    density = Pixel2density(mean);
                    rflag = std::fwrite(&density, sizeof(G4float), 1, fileOut);
                }
            }

        }
    }

    rflag = std::fclose(fileOut);

    delete [] nameProcessed;

        for ( G4int i = 0; i < fRows; i ++ ) {
        delete [] fTab[i];
        }
     delete [] fTab;
     */

    if (rflag) return returnvalue;
    return returnvalue;
}


// cette fonction on la donne la densité , il retourne l'indice de matériau correspondant
// called from ReadData() ,
unsigned int G4TDCMHandler::GetMaterialIndex( G4float density )
{
    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMHandler::GetMaterialIndex  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    std::map<G4float,G4String>::reverse_iterator ite;
    G4int ii = fMaterialIndices.size();
    for( ite = fMaterialIndices.rbegin(); ite != fMaterialIndices.rend();
         ite++, ii-- ) {
        if( density >= (*ite).first ) {
            break;
        }
    }

    if(static_cast<unsigned int>(ii) == fMaterialIndices.size())
    { ii = fMaterialIndices.size()-1; }

    //G4cout << " GetMaterialIndex " << density << " = " << ii << G4endl;

    return  ii;

}


// this function use the compression value to calculate the density or mean density to a pixel point or n pixel point (compression)
// and send it to the G4TDCMZSliceHeader to save it in file of this slice
// utilisant G4TDCMZSliceHeader dcmPZSH comme argument envoyer par G4TDCMHandler::ReadFile() et qui contien les premièrs donnés comme nmat et le mat et ...
// pour compléter les données de l'indice de mat et les densitées pour chaque pixel , et finalement l'object dcmPZSH call the function to save all this data to file
// .g4dcm sended ti it befor from G4TDCMHandler::ReadFile();
// called from G4TDCMHandler::ReadFile(),
// it call Pixel2density(), dcmPZSH->AddValue(density) , dcmPZSH->AddMateID(GetMaterialIndex(density)) , dcmPZSH->FlipData();
void G4TDCMHandler::StoreData(G4TDCMZSliceHeader* dcmPZSH)
{
    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMHandler::StoreData  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    G4int mean; // la valeur moyenne de pixel ou un nombre de pixels(si compression est différent de 1)
    G4double density;
    G4bool overflow = false;

    if(!dcmPZSH) { return; }

    dcmPZSH->SetSliceLocation(fSliceLocation);

    //----- Print indices of material

    if(fCompression == 1) { // no fCompression: each pixel has a density value)
        for( G4int ww = 0; ww < fRows; ww++) {
            dcmPZSH->AddRow();
            for( G4int xx = 0; xx < fColumns; xx++) {
                mean = fTab[ww][xx];
                density = Pixel2density(mean);
                dcmPZSH->AddValue(density);
                dcmPZSH->AddMateID(GetMaterialIndex(density));
            }
        }

    } else {

        // density value is the average of a square region of fCompression*fCompression pixels

        for(G4int ww = 0; ww < fRows ;ww += fCompression ) {
            dcmPZSH->AddRow();
            for(G4int xx = 0; xx < fColumns ;xx +=fCompression ) {
                overflow = false;
                mean = 0;
                for(int sumx = 0; sumx < fCompression; sumx++) {
                    for(int sumy = 0; sumy < fCompression; sumy++) {
                        if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
                        mean += fTab[ww+sumy][xx+sumx];
                    }
                    if(overflow) break;
                }
                mean /= fCompression*fCompression;

                if(!overflow) {
                    density = Pixel2density(mean);
                    dcmPZSH->AddValue(density);
                    dcmPZSH->AddMateID(GetMaterialIndex(density));
                }
            }
        }
    }

    // This function flips all the data, Otherwise, the image is upside-down
    dcmPZSH->FlipData();
}


// called from ReadData() , StoreData(string), StoreData(File),
// it call ReadCalibration() ,
G4float G4TDCMHandler::Pixel2density(G4int pixel)
{
    //G4cout<< "@@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMHandler::Pixel2density  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    if(!fReadCalibration) {
        //ReadCalibration();
        fNbrequali = ValueCT.size();
        fValueDensity = new G4double[fNbrequali];
        fValueCT = new G4double[fNbrequali];

        for (G4int fg = 0 ; fg < fNbrequali ; fg++){
            fValueDensity[fg] = ValueDensity[fg];
            fValueCT[fg] = ValueCT[fg];
            //G4cout<< "ValueDensity: " <<fValueDensity[fg]<< " - ValueCT: "<< fValueCT[fg] << G4endl;

        }
        fReadCalibration = true;
    }

    G4float density = -1.;
    G4double deltaCT = 0;
    G4double deltaDensity = 0;

    for(G4int j = 1; j < fNbrequali; j++) {
        if( pixel >= fValueCT[j-1] && pixel < fValueCT[j]) {

            deltaCT = fValueCT[j] - fValueCT[j-1];
            deltaDensity = fValueDensity[j] - fValueDensity[j-1];

            // interpolating linearly
            density = fValueDensity[j] - ((fValueCT[j] - pixel)*deltaDensity/deltaCT );
            break;
        }
    }

    if(density < 0.) {
        std::printf("@@@ Error density = %f && Pixel = %i \(0x%x) && deltaDensity/deltaCT = %f\n",density,pixel,pixel, deltaDensity/deltaCT);
    }

    //G4cout<< "pixel: " <<pixel<< " - density: "<< density << G4endl;

    return density;

}



// called from various of functions in this classe
template <class Type>
void G4TDCMHandler::GetValue(char * _val, Type & _rval) {

    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMHandler::GetValue  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

#if BYTE_ORDER == BIG_ENDIAN
    if(fLittleEndian) {      // little endian
#else // BYTE_ORDER == LITTLE_ENDIAN
    if(!fLittleEndian) {     // big endian
#endif
        const int SIZE = sizeof(_rval);
        char ctemp;
        for(int i = 0; i < SIZE/2; i++) {
            ctemp = _val[i];
            _val[i] = _val[SIZE - 1 - i];
            _val[SIZE - 1 - i] = ctemp;
        }
    }
    _rval = *(Type *)_val;
}

// called from ReadFile()
// it call GetValue()
G4int G4TDCMHandler::read_defined_nested(FILE * nested,G4int SQ_Length)
{

    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMHandler::read_defined_nested  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    //      VARIABLES
    unsigned short item_GroupNumber;
    unsigned short item_ElementNumber;
    G4int item_Length;
    G4int items_array_length=0;
    char * buffer= new char[LINEBUFFSIZE];
    size_t rflag = 0;

    while(items_array_length < SQ_Length)
    {
        rflag = std::fread(buffer, 2, 1, nested);
        GetValue(buffer, item_GroupNumber);

        rflag = std::fread(buffer, 2, 1, nested);
        GetValue(buffer, item_ElementNumber);

        rflag = std::fread(buffer, 4, 1, nested);
        GetValue(buffer, item_Length);

        rflag = std::fread(buffer, item_Length, 1, nested);

        items_array_length= items_array_length+8+item_Length;
    }

    delete [] buffer;

    if( SQ_Length>items_array_length )
        return 0;
    else
        return 1;
    if (rflag) return 1;
}

// called from ReadFile()
// it call GetValue(), read_undefined_nested(),
void G4TDCMHandler::read_undefined_nested(FILE * nested)
{
    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMHandler::read_undefined_nested  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    //      VARIABLES
    unsigned short item_GroupNumber;
    unsigned short item_ElementNumber;
    unsigned int item_Length;
    char * buffer= new char[LINEBUFFSIZE];
    size_t rflag = 0;

    do
    {
        rflag = std::fread(buffer, 2, 1, nested);
        GetValue(buffer, item_GroupNumber);

        rflag = std::fread(buffer, 2, 1, nested);
        GetValue(buffer, item_ElementNumber);

        rflag = std::fread(buffer, 4, 1, nested);
        GetValue(buffer, item_Length);

        if(item_Length!=0xffffffff)
            rflag = std::fread(buffer, item_Length, 1, nested);
        else
            read_undefined_item(nested);


    } while(item_GroupNumber!=0xFFFE || item_ElementNumber!=0xE0DD
            || item_Length!=0);

    delete [] buffer;
    if (rflag) return;
}

// called from read_undefined_nested()
// it call GetValue()
void G4TDCMHandler::read_undefined_item(FILE * nested)
{
    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMHandler::read_undefined_item  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    //      VARIABLES
    unsigned short item_GroupNumber;
    unsigned short item_ElementNumber;
    G4int item_Length; size_t rflag = 0;
    char *buffer= new char[LINEBUFFSIZE];

    do
    {
        rflag = std::fread(buffer, 2, 1, nested);
        GetValue(buffer, item_GroupNumber);

        rflag = std::fread(buffer, 2, 1, nested);
        GetValue(buffer, item_ElementNumber);

        rflag = std::fread(buffer, 4, 1, nested);
        GetValue(buffer, item_Length);


        if(item_Length!=0)
            rflag = std::fread(buffer,item_Length,1,nested);

    }
    while(item_GroupNumber!=0xFFFE || item_ElementNumber!=0xE00D || item_Length!=0);

    delete [] buffer;
    if (rflag) return;
}




//################################################################################# not used for now




// This function is depreciated as it is handled by G4TDCMZSliceHeader::DumpToFile
// called from
// it call GetMaterialIndex() , Pixel2density() ,
void G4TDCMHandler::StoreData(std::ofstream& foutG4DCM)
{
    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMHandler::StoreData  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    G4int mean;
    G4double density;
    G4bool overflow = false;

    //----- Print indices of material
    if(fCompression == 1) { // no fCompression: each pixel has a density value)
        for( G4int ww = 0; ww < fRows; ww++) {
            for( G4int xx = 0; xx < fColumns; xx++) {
                mean = fTab[ww][xx];
                density = Pixel2density(mean);
                foutG4DCM << GetMaterialIndex( density ) << " ";
            }
            foutG4DCM << G4endl;
        }

    } else {
        // density value is the average of a square region of
        // fCompression*fCompression pixels
        for(G4int ww = 0; ww < fRows ;ww += fCompression ) {
            for(G4int xx = 0; xx < fColumns ;xx +=fCompression ) {
                overflow = false;
                mean = 0;
                for(int sumx = 0; sumx < fCompression; sumx++) {
                    for(int sumy = 0; sumy < fCompression; sumy++) {
                        if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
                        mean += fTab[ww+sumy][xx+sumx];
                    }
                    if(overflow) break;
                }
                mean /= fCompression*fCompression;

                if(!overflow) {
                    density = Pixel2density(mean);
                    foutG4DCM << GetMaterialIndex( density ) << " ";
                }
            }
            foutG4DCM << G4endl;
        }

    }

    //----- Print densities
    if(fCompression == 1) { // no fCompression: each pixel has a density value)
        for( G4int ww = 0; ww < fRows; ww++) {
            for( G4int xx = 0; xx < fColumns; xx++) {
                mean = fTab[ww][xx];
                density = Pixel2density(mean);
                foutG4DCM << density << " ";
                if( xx%8 == 3 ) foutG4DCM << G4endl; // just for nicer reading
            }
        }

    } else {
        // density value is the average of a square region of
        // fCompression*fCompression pixels
        for(G4int ww = 0; ww < fRows ;ww += fCompression ) {
            for(G4int xx = 0; xx < fColumns ;xx +=fCompression ) {
                overflow = false;
                mean = 0;
                for(int sumx = 0; sumx < fCompression; sumx++) {
                    for(int sumy = 0; sumy < fCompression; sumy++) {
                        if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
                        mean += fTab[ww+sumy][xx+sumx];
                    }
                    if(overflow) break;
                }
                mean /= fCompression*fCompression;

                if(!overflow) {
                    density = Pixel2density(mean);
                    foutG4DCM << density  << " ";
                    if( xx/fCompression%8 == 3 ) foutG4DCM << G4endl; // just for nicer
                    // reading
                }
            }
        }

    }

}
// uses the file data that is sended from CheckFileFormat(), if fill vector of fMaterialIndices by reading the matname and density fMaterialIndices[densityMax] = mateName;
// called from CheckFileFormat() ,
// it call
void G4TDCMHandler::ReadMaterialIndices( std::ifstream& finData)
{

    G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMHandler::ReadMaterialIndices  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    unsigned int nMate;
    G4String mateName;
    G4float densityMax;
    finData >> nMate;

    if( finData.eof() ) return;

    G4cout << " ReadMaterialIndices from data.dat with number is " << nMate << G4endl;
    for( unsigned int ii = 0; ii < nMate; ii++ ){
        finData >> mateName >> densityMax;
        fMaterialIndices[densityMax] = mateName;
        G4cout << ii << " ReadMaterialIndices " << mateName << " " << ", density Max : " << densityMax << G4endl;
    }


}
void G4TDCMHandler::ReadDicomDataFile()
{

    G4String line , indicator, word;
    std::ifstream fileR(DICOMInputFile.c_str());

    if(fileR.is_open()){

        G4cout  << " \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << G4endl;
        G4cout << "Reading "<< DICOMInputFile << " ... " << G4endl ;

        while (getline(fileR, line)) {

            //G4cout << " the line " << line << G4endl ;

            std::istringstream LineString(line);

            if(LineString.str().empty()){
                continue;
            }

            LineString >> word;  // the white spaces more than 1 are considered as 1

            //G4cout << " the word " << word << G4endl ;
            if(word.empty() || word == "" || word == "#"){ // || word.isNull()
                //G4cout <<  "word " << word << G4endl;
                continue;
            }
            if (word == ">>") {

                LineString >> indicator;
                //G4cout <<  "indicator " << indicator << G4endl;
                continue;
            }
            else{

                //LineString >> word;

                G4String organName, mediumName;
                G4int ik = 0, MediumId, organId;
                const std::string & Organs_To_Score = LineString.str();
                char delim = ' ';
                std::stringstream ss (Organs_To_Score);
                std::string item;

                if(indicator == "dicom_calibration"){

                    while (getline (ss, item, delim)) {
                        G4String ii = item;
                        if(ii.empty()){ //ii.isNull() ||
                            continue;
                        }
                        if(ik == 0 ){
                            ValueCT.push_back(atof(item.c_str()));
                        }
                        else if(ik == 1 ){
                            ValueDensity.push_back(atof(item.c_str()));
                        }

                        //G4cout  << "item : " << item << G4endl;
                        ik++;
                    }
                }

                else if(indicator == "material_density"){

                    while (getline (ss, item, delim)) {
                        G4String ii = item;
                        if(ii.empty()){ //ii.isNull() ||
                            continue;
                        }
                        if(ik==0){
                            mediumName = item;
                            //G4cout  << "mediumName : " << mediumName << G4endl;
                        }else if(ik == 1 ){
                            fMaterialIndices[atof(item.c_str())] = mediumName;
                            G4cout  << mediumName << " " << atof(item.c_str()) << G4endl;
                        }
                        //G4cout  << "item : " << item << G4endl;
                        ik++;
                    }
                }
                else{

                }
            }
        }

        G4cout  << " \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << G4endl;

        fileR.close();
    }
    else{

        G4cout << "Cannot open the file " << DICOMInputFile.c_str() << G4endl ;
    }

    if(ValueCT.size() != ValueDensity.size()){
        G4String func = "G4TDCMHandler::" ;
        func+= __FUNCTION__;
        G4Exception(func , "G4TDCM001", FatalException, "CT values != density values" );
    }

    fValueDensity = new G4double[ValueCT.size()];
    fValueCT = new G4double[ValueCT.size()];

    for (G4int fg = 0 ; fg < ValueCT.size() ; fg++){
        fValueDensity[fg] = ValueDensity[fg];
        fValueCT[fg] = ValueCT[fg];
    }

    fReadCalibration = true;


    // to make the map of material from low density to the high density
    typedef std::pair<G4int,G4double> pair;

    // input map
    std::map<G4int,G4double> map = { {5, 2.}, {9, 1.}, {6, 4.}, {2, 3.} };

    // create a empty vector of pairs
    std::vector<pair> vec;

    // copy key-value pairs from the map to the vector
    std::copy(map.begin(), map.end(), std::back_inserter<std::vector<pair>>(vec));

    // sort the vector by increasing order of its pair's second value, if second value are equal, order by the pair's first value
    std::sort(vec.begin(), vec.end(), [](const pair& l, const pair& r) {
        if (l.second != r.second) return l.second < r.second;
        return l.first < r.first;
    });

    // print the vector
    for (auto const &pair: vec) {
        std::cout << '{' << pair.first << "," << pair.second << '}' << '\n';
    }



}
// Separated out of Pixel2density
// No need to read in same calibration EVERY time
// Increases the speed of reading file by several orders of magnitude
// called from Pixel2density() ,
// it call
void G4TDCMHandler::ReadCalibration()
{
    //G4cout<< " @@@@@@@@@@@@@@@@@@@@@@@@@@  in function G4TDCMHandler::ReadCalibration  @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< G4endl;

    /*
    fNbrequali = 0;

    // CT2Density.dat contains the calibration curve to convert CT (Hounsfield)
    // number to physical density
    std::ifstream calibration(fCt2DensityFile.c_str());
    calibration >> fNbrequali;

    fValueDensity = new G4double[fNbrequali];
    fValueCT = new G4double[fNbrequali];

    if(!calibration) {
        G4Exception("G4TDCMHandler::ReadFile", "G4TDCM001", FatalException, "@@@ No value to transform pixels in density!");

    } else { // calibration was successfully opened
        for(G4int i = 0; i < fNbrequali; i++) { // Loop to store all the
            //pts in CT2Density.dat
            calibration >> fValueCT[i] >> fValueDensity[i];
        }
    }

    calibration.close();
    fReadCalibration = true;
    */

/*
    fNbrequali = ValueCT.size();
    fValueDensity = new G4double[fNbrequali];
    fValueCT = new G4double[fNbrequali];

    for (G4int fg = 0 ; fg < fNbrequali ; fg++){
        fValueDensity[fg] = ValueDensity[fg];
        fValueCT[fg] = ValueCT[fg];
    }
    fReadCalibration = true;
*/

}

