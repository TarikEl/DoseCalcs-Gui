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
#include "DicomVFileImage.hh"
#include "DicomFileStructure.hh"
#include "DicomROI.hh"

#include "G4GeometryTolerance.hh"

#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/dcmdata/dcdeftag.h"
#include "dcmtk/dcmdata/dcpixel.h"
#include "dcmtk/dcmdata/dcpxitem.h"
#include "dcmtk/dcmdata/dcpixseq.h"
#include "dcmtk/dcmrt/drtimage.h"
#include "dcmtk/dcmiod/iodutil.h"

#include <set>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomVFileImage::DicomVFileImage()
{
    theFileMgr = DicomFileMgr::GetInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomVFileImage::DicomVFileImage(DcmDataset* dset) : DicomVFile(dset)
{
    theFileMgr = DicomFileMgr::GetInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomVFileImage::ReadData()
{

    std::vector<double> dImagePositionPatient = Read1Data(theDataset, DCM_ImagePositionPatient,3);
    fLocation = dImagePositionPatient[2];
    std::vector<double> dSliceThickness = Read1Data(theDataset, DCM_SliceThickness, 1);
    std::vector<double> dPixelSpacing = Read1Data(theDataset, DCM_PixelSpacing, 2);

    std::vector<double> dRows = Read1Data(theDataset, DCM_Rows, 1);
    std::vector<double> dColumns = Read1Data(theDataset, DCM_Columns, 1);
    fNoVoxelY = dRows[0];
    fNoVoxelX = dColumns[0];
    fNoVoxelZ = 1;

    fMinX = dImagePositionPatient[0]; // center of upper corner of pixel?
    fMaxX = dImagePositionPatient[0]+dColumns[0]*dPixelSpacing[0];

    fMinY = dImagePositionPatient[1];
    fMaxY = dImagePositionPatient[1]+dRows[0]*dPixelSpacing[1];

    fMinZ = dImagePositionPatient[2]-dSliceThickness[0]/2.;
    fMaxZ = dImagePositionPatient[2]+dSliceThickness[0]/2.;
    fVoxelDimX = dPixelSpacing[0];
    fVoxelDimY = dPixelSpacing[1];
    fVoxelDimZ = dSliceThickness[0];

    if( DicomFileMgr::verbose >= debugVerb ) G4cout << " DicomVFileImage::ReadData:  fNoVoxel "
                                                    << fNoVoxelX << " " << fNoVoxelY << " " << fNoVoxelZ << G4endl;
    if( DicomFileMgr::verbose >= debugVerb ) G4cout << " DicomVFileImage::ReadData:  fMin "
                                                    << fMinX << " " << fMinY << " " << fMinZ << G4endl;
    if( DicomFileMgr::verbose >= debugVerb ) G4cout << " DicomVFileImage::ReadData:  fMax "
                                                    << fMaxX << " " << fMaxY << " " << fMaxZ << G4endl;
    if( DicomFileMgr::verbose >= debugVerb ) G4cout << " DicomVFileImage::ReadData:  fVoxelDim "
                                                    << fVoxelDimX << " " << fVoxelDimY << " " << fVoxelDimZ << G4endl;

    std::vector<double> dImageOrientationPatient =
            Read1Data(theDataset, DCM_ImageOrientationPatient,6);
    fOrientationRows = G4ThreeVector(dImageOrientationPatient[0],dImageOrientationPatient[1],
            dImageOrientationPatient[2]);
    fOrientationColumns = G4ThreeVector(dImageOrientationPatient[3],dImageOrientationPatient[4],
            dImageOrientationPatient[5]);

    if( fOrientationRows != G4ThreeVector(1,0,0)
            || fOrientationColumns != G4ThreeVector(0,1,0) ) {
        G4cerr << " OrientationRows " << fOrientationRows << " OrientationColumns "
               << fOrientationColumns << G4endl;
        G4Exception("DicomVFileImage::ReadData",
                    "DFCT0002",
                    JustWarning,
                    "OrientationRows must be (1,0,0) and OrientationColumns (0,1,0), please contact GAMOS authors");
    }
    fBitAllocated = Read1Data(theDataset, DCM_BitsAllocated, 1)[0];
    if( DicomFileMgr::verbose >= 4 ) G4cout << " BIT ALLOCATED " << fBitAllocated << G4endl;

    std::vector<double> dRescaleSlope = Read1Data(theDataset, DCM_RescaleSlope, 1);
    if( dRescaleSlope.size() == 1 ) {
        fRescaleSlope = dRescaleSlope[0];
    } else {
        fRescaleSlope = 1;
    }
    std::vector<double> dRescaleIntercept = Read1Data(theDataset, DCM_RescaleIntercept, 1);
    if( dRescaleIntercept.size() == 1 ) {
        fRescaleIntercept = dRescaleIntercept[0];
    } else {
        fRescaleIntercept = 1;
    }

    /*
    OFString dt ; double fv;

    DcmIODUtil::getStringValueFromItem(DCM_RadiopharmaceuticalStartTime, *theDataset, dt, 0);
    getDateTimeValues("DCM_RadiopharmaceuticalStartTime",dt);

    DcmIODUtil::getStringValueFromItem(DCM_RadiopharmaceuticalStopTime, *theDataset, dt, 0);
    getDateTimeValues("DCM_RadiopharmaceuticalStopTime",dt);

    DcmIODUtil::getFloat64ValueFromItem(DCM_RadionuclideTotalDose, *theDataset, fv, 0); PETGeneralData["DCM_RadionuclideTotalDose"]=fv;
    DcmIODUtil::getFloat64ValueFromItem(DCM_RadionuclidePositronFraction, *theDataset, fv, 0); PETGeneralData["DCM_RadionuclidePositronFraction"]=fv;
    DcmIODUtil::getFloat64ValueFromItem(DCM_RadiopharmaceuticalSpecificActivity, *theDataset, fv, 0); PETGeneralData["DCM_RadiopharmaceuticalSpecificActivity"]=fv;
    DcmIODUtil::getFloat64ValueFromItem(DCM_RadiopharmaceuticalVolume, *theDataset, fv, 0); PETGeneralData["DCM_RadiopharmaceuticalVolume"]=fv;

    G4cout << "\n\nGeneral PET Data from ReadData()): " << G4endl;
    for ( auto it = PETGeneralData.begin(); it != PETGeneralData.end(); ++it  ){
        G4cout << it->first << " " << it->second << G4endl;
    }
*/
    //G4cout  << " \n---- The Data readed from DICOM file are: Compression " << theFileMgr->GetCompression() << " Location " << fLocation << " NoVoxelXYZ " << " " << fNoVoxelX
      //      << " " << fNoVoxelY << " " << fNoVoxelZ << " VoxelDimXYZ " << fVoxelDimX  << " " << fVoxelDimY << " " << fVoxelDimZ
        //    << " MinXYZ " << fMinX  << " " << fMinY << " " << fMinZ << " fOrientationRows " << " " << fOrientationRows
          //  << " OrientationColumns " << fOrientationColumns << " BitAllocated " << fBitAllocated << " RescaleSlope " << fRescaleSlope
            //<< " RescaleIntercept " << fRescaleIntercept << "\n" << G4endl;

    G4cout  << " Reading Pixel Data for Slice with Location " << fLocation << " MinZ " << fMinZ << G4endl;

    ReadPixelData();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomVFileImage::ReadPixelData()
{
    //  READING THE PIXELS :
    OFCondition result = EC_Normal;
    //---- CHECK IF DATA IS COMPRESSED
    DcmElement* element = NULL;
    result = theDataset->findAndGetElement(DCM_PixelData, element);
    if (result.bad() || element == NULL) {
        G4Exception("ReadData",
                    "findAndGetElement(DCM_PixelData, ",
                    FatalException,
                    ("Element PixelData not found: " + G4String(result.text())).c_str());
    }
    DcmPixelData *dpix = NULL;
    dpix = OFstatic_cast(DcmPixelData*, element);
    // If we have compressed data, we must utilize DcmPixelSequence
    //   in order to access it in raw format, e. g. for decompressing it
    //   with an external library.
    DcmPixelSequence *dseq = NULL;
    E_TransferSyntax xferSyntax = EXS_Unknown;
    const DcmRepresentationParameter *rep = NULL;
    // Find the key that is needed to access the right representation of the data within DCMTK
    dpix->getOriginalRepresentationKey(xferSyntax, rep);
    // Access original data representation and get result within pixel sequence
    result = dpix->getEncapsulatedRepresentation(xferSyntax, rep, dseq);
    if ( result == EC_Normal ) // COMPRESSED DATA
    {
        G4Exception("DicomVFileImage::ReadData()",
                    "DFCT004",
                    FatalException,
                    "Compressed pixel data is not supported");

        if( DicomFileMgr::verbose >= debugVerb ) G4cout
                << " DicomVFileImage::ReadData:  result == EC_Normal Reading compressed data " << std::endl;
        DcmPixelItem* pixitem = NULL;
        // Access first frame (skipping offset table)
        for( int ii = 1; ii < 2; ii++ ) {
            OFCondition cond =  dseq->getItem(pixitem, ii);
            if( !cond.good()) break;
            G4cout << ii << " PIX LENGTH " << pixitem->getLength() << G4endl;
        }
        if (pixitem == NULL) {
            G4Exception("ReadData",
                        "dseq->getItem()",
                        FatalException,
                        "No DcmPixelItem in DcmPixelSequence");
        }
        Uint8* pixData = NULL;
        // Get the length of this pixel item
        // (i.e. fragment, i.e. most of the time, the lenght of the frame)
        Uint32 length = pixitem->getLength();
        if (length == 0) {
            G4Exception("ReadData",
                        "pixitem->getLength()",
                        FatalException,
                        "PixelData empty");
        }

        if( DicomFileMgr::verbose >= debugVerb ) G4cout
                << " DicomVFileImage::ReadData:  number of pixels " << length << G4endl;
        // Finally, get the compressed data for this pixel item
        result = pixitem->getUint8Array(pixData);
    } else { // UNCOMPRESSED DATA
        if(fBitAllocated == 8) { // Case 8 bits :
            Uint8* pixData = NULL;
            if(! (element->getUint8Array(pixData)).good() ) {
                G4Exception("ReadData",
                            "getUint8Array pixData, ",
                            FatalException,
                            ("PixelData not found: " + G4String(result.text())).c_str());
            }
            for( int ir = 0; ir < fNoVoxelY; ir++ ) {
                for( int ic = 0; ic < fNoVoxelX; ic++ ) {
                    fHounsfieldV.push_back(pixData[ic+ir*fNoVoxelX]*fRescaleSlope + fRescaleIntercept);
                }
            }
        } else if(fBitAllocated == 16) { // Case 16 bits :
            Uint16* pixData = NULL;
            if(! (element->getUint16Array(pixData)).good() ) {
                G4Exception("ReadData",
                            "getUint16Array pixData, ",
                            FatalException,
                            ("PixelData not found: " + G4String(result.text())).c_str());
            }
            for( int ir = 0; ir < fNoVoxelY; ir++ ) {
                for( int ic = 0; ic < fNoVoxelX; ic++ ) {
                    fHounsfieldV.push_back(pixData[ic+ir*fNoVoxelX]*fRescaleSlope + fRescaleIntercept);
                    //G4cout << " ic " << ic << " ir " << ir << " fNoVoxelX " << fNoVoxelX << " fRescaleSlope " << fRescaleSlope << " fRescaleIntercept "<< fRescaleIntercept << " fHounsfieldV " << pixData[ic+ir*fNoVoxelX]*fRescaleSlope + fRescaleIntercept << G4endl;
                }
            }
        } else if(fBitAllocated == 32) { // Case 32 bits :
            Uint32* pixData = NULL;
            if(! (element->getUint32Array(pixData)).good() ) {
                G4Exception("ReadData",
                            "getUint32Array pixData, ",
                            FatalException,
                            ("PixelData not found: " + G4String(result.text())).c_str());
            }
            for( int ir = 0; ir < fNoVoxelY; ir++ ) {
                for( int ic = 0; ic < fNoVoxelX; ic++ ) {
                    fHounsfieldV.push_back(pixData[ic+ir*fNoVoxelX]*fRescaleSlope + fRescaleIntercept);
                }
            }
        }
    }

    //G4cout << fHounsfieldV.size() << G4endl;
    //for( int jj = 0; jj < fHounsfieldV.size(); jj ++ ) {
      //  if(fHounsfieldV[jj] != 0.){
        //    G4cout << fHounsfieldV[jj] << " ";
        //}
    //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomVFileImage::operator+=( const DicomVFileImage& rhs )
{
    *this = *this + rhs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomVFileImage DicomVFileImage::operator+( const DicomVFileImage& rhs )
{
    //----- Check that both slices has the same dimensions
    if( fNoVoxelX != rhs.GetNoVoxelX()
            || fNoVoxelY != rhs.GetNoVoxelY() ) {
        G4cerr << "DicomVFileImage error adding two slice headers:\
                  !!! Different number of voxels: "
               << "  X= " << fNoVoxelX << " =? " << rhs.GetNoVoxelX()
               << "  Y=  " << fNoVoxelY << " =? " << rhs.GetNoVoxelY()
               << "  Z=  " << fNoVoxelZ << " =? " << rhs.GetNoVoxelZ()
               << G4endl;
        G4Exception("DicomVFileImage::DicomVFileImage",
                    "",FatalErrorInArgument,"");
    }
    //----- Check that both slices has the same extensions
    if( fMinX != rhs.GetMinX() || fMaxX != rhs.GetMaxX()
            || fMinY != rhs.GetMinY() || fMaxY != rhs.GetMaxY() ) {
        G4cerr << "DicomVFileImage error adding two slice headers:\
                  !!! Different extensions: "
               << "  Xmin= " << fMinX << " =? " << rhs.GetMinX()
               << "  Xmax= " << fMaxX << " =? " << rhs.GetMaxX()
               << "  Ymin= " << fMinY << " =? " << rhs.GetMinY()
               << "  Ymax= " << fMaxY << " =? " << rhs.GetMaxY()
               << G4endl;
        G4Exception("DicomVFileImage::operator+","",
                    FatalErrorInArgument,"");
    }

    //----- Check that both slices has the same orientations
    if( fOrientationRows != rhs.GetOrientationRows() ||
            fOrientationColumns != rhs.GetOrientationColumns() ) {
        G4cerr << "DicomVFileImage error adding two slice headers: !!!\
                  Slices have different orientations "
               << "  Orientation Rows = " << fOrientationRows << " & " << rhs.GetOrientationRows()
               << "  Orientation Columns " << fOrientationColumns << " & "
               << rhs.GetOrientationColumns() << G4endl;
        G4Exception("DicomVFileImage::operator+","",
                    FatalErrorInArgument,"");
    }

    //----- Check that the slices are contiguous in Z
    if( std::fabs( fMinZ - rhs.GetMaxZ() ) >
            G4GeometryTolerance::GetInstance()->GetRadialTolerance() &&
            std::fabs( fMaxZ - rhs.GetMinZ() ) >
            G4GeometryTolerance::GetInstance()->GetRadialTolerance() ){
        G4cerr << "DicomVFileImage error adding two slice headers: !!!\
                  Slices are not contiguous in Z "
               << "  Zmin= " << fMinZ << " & " << rhs.GetMinZ()
               << "  Zmax= " << fMaxZ << " & " << rhs.GetMaxZ()
               << G4endl;
        G4Exception("DicomVFileImage::operator+","",
                    JustWarning,"");
    }

    //----- Build slice header copying first one
    DicomVFileImage temp( *this );

    //----- Add data from second slice header
    temp.SetMinZ( std::min( fMinZ, rhs.GetMinZ() ) );
    temp.SetMaxZ( std::max( fMaxZ, rhs.GetMaxZ() ) );
    temp.SetNoVoxelZ( fNoVoxelZ + rhs.GetNoVoxelZ() );

    return temp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomVFileImage::DumpHeaderToTextFile(std::ofstream& fout)
{
    if( DicomFileMgr::verbose >= warningVerb ) G4cout << fLocation << " DumpHeaderToTextFile "
                                                      << fFileName << " " << fHounsfieldV.size() << G4endl;

    G4String fName = fFileName.substr(0,fFileName.length()-3) + "g4dcm";
    std::ofstream out(fName.c_str());

    if( DicomFileMgr::verbose >= warningVerb ) G4cout
            << "### DicomVFileImage::Dumping Z Slice header to Text file " << G4endl;

    G4int fCompress = theFileMgr->GetCompression();
    fout << fNoVoxelX/fCompress << " " << fNoVoxelY/fCompress << " " << fNoVoxelZ << std::endl;
    fout << fMinX << " " << fMaxX << std::endl;
    fout << fMinY << " " << fMaxY << std::endl;
    fout << fMinZ << " " << fMaxZ << std::endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomVFileImage::Print(std::ostream& out )
{
    G4int fCompress = theFileMgr->GetCompression();
    out << "@@ CT Slice " << fLocation << G4endl;

    out << "@ NoVoxels " << fNoVoxelX/fCompress << " " << fNoVoxelY/fCompress << " "
        << fNoVoxelZ << G4endl;
    out << "@ DIM X: " << fMinX << " " << fMaxX
        << " Y: " << fMinY << " " << fMaxY
        << " Z: " << fMinZ << " " << fMaxZ << G4endl;
}

void DicomVFileImage::ReadPETGeneralData(){

    PETGeneralData.clear();

    G4cout << "\n------------------------ General PET Files Data: \n" << G4endl;

    // get this data for all series
    OFString dt ; double fv;

    if( !theDataset->findAndGetOFString(DCM_AcquisitionDate,dt).good() ) { G4Exception("DicomHandler ", "",JustWarning, " Have not read DCM_RadiopharmaceuticalStartTime"); }
    getDateTimeValues("DCM_AcquisitionDate",dt);

    //if( !theDataset->findAndGetOFString(DCM_RadiopharmaceuticalStartTime,dt).good() ) { G4Exception("DicomHandler ", "",JustWarning, " Have not read DCM_RadiopharmaceuticalStartTime"); }
    //getDateTimeValues("DCM_RadiopharmaceuticalStartTime",dt);

    //if( !theDataset->findAndGetOFString(DCM_RadiopharmaceuticalStopTime,dt).good() ) { G4Exception("DicomHandler ", "",JustWarning, " Have not read DCM_RadiopharmaceuticalStopTime"); }
    //getDateTimeValues("DCM_RadiopharmaceuticalStopTime",dt);

    if( !theDataset->findAndGetFloat64(DCM_RadionuclideTotalDose,fv).good() ) { G4Exception("DicomHandler ", "",JustWarning, " Have not read DCM_RadionuclideTotalDose"); } PETGeneralData["DCM_RadionuclideTotalDose"]=fv;
    if( !theDataset->findAndGetFloat64(DCM_RadionuclideHalfLife,fv).good() ) { G4Exception("DicomHandler ", "",JustWarning, " Have not read DCM_RadionuclideHalfLife"); } PETGeneralData["DCM_RadionuclideHalfLife"]=fv;
    if( !theDataset->findAndGetFloat64(DCM_RadionuclidePositronFraction,fv).good() ) { G4Exception("DicomHandler ", "",JustWarning, " Have not read DCM_RadionuclidePositronFraction"); } PETGeneralData["DCM_RadionuclidePositronFraction"]=fv;
    if( !theDataset->findAndGetFloat64(DCM_RadiopharmaceuticalSpecificActivity,fv).good() ) { G4Exception("DicomHandler ", "",JustWarning, " Have not read DCM_RadiopharmaceuticalSpecificActivity"); } PETGeneralData["DCM_RadiopharmaceuticalSpecificActivity"]=fv;
    if( !theDataset->findAndGetFloat64(DCM_RadiopharmaceuticalVolume,fv).good() ) { G4Exception("DicomHandler ", "",JustWarning, " Have not read DCM_RadiopharmaceuticalVolume"); } PETGeneralData["DCM_RadiopharmaceuticalVolume"]=fv;

    /*
    DcmIODUtil::getFloat64ValueFromItem(DCM_RadionuclideTotalDose, *theDataset, fv, 0); PETGeneralData["DCM_RadionuclideTotalDose"]=fv;
    DcmIODUtil::getFloat64ValueFromItem(DCM_RadionuclideHalfLife, *theDataset, fv, 0); PETGeneralData["DCM_RadionuclideHalfLife"]=fv;
    DcmIODUtil::getFloat64ValueFromItem(DCM_RadionuclidePositronFraction, *theDataset, fv, 0); PETGeneralData["DCM_RadionuclidePositronFraction"]=fv;
    DcmIODUtil::getFloat64ValueFromItem(DCM_RadiopharmaceuticalSpecificActivity, *theDataset, fv, 0); PETGeneralData["DCM_RadiopharmaceuticalSpecificActivity"]=fv;
    DcmIODUtil::getFloat64ValueFromItem(DCM_RadiopharmaceuticalVolume, *theDataset, fv, 0); PETGeneralData["DCM_RadiopharmaceuticalVolume"]=fv;
    */
    for ( auto it = PETGeneralData.begin(); it != PETGeneralData.end(); ++it  ){
        G4cout << it->first << " " << it->second << G4endl;
    }

}

std::vector<double> DicomVFileImage::getSerieAcquisitionsData(){

    // get this data for each serie

    G4cout << "\n------------------------ Serie Acquisition Data \n" << G4endl;

    OFString dt ; double fv; OFDateTime dateTimeV; std::vector<double> vec;

    //DcmIODUtil::getStringValueFromItem(DCM_AcquisitionDateTime, *theDataset, dt, 0);
    //if( !theDataset->findAndGetOFString(DCM_AcquisitionDateTime,dt).good() ) { G4Exception("DicomHandler ", "",JustWarning, " Have not read DCM_AcquisitionDateTime"); }
    //DcmDateTime::getOFDateTimeFromString(dt,dateTimeV);
    //G4cout << "DCM_AcquisitionDateTime " << dateTimeV << G4endl;

    //if( !theDataset->findAndGetOFString(DCM_AcquisitionTime,dt).good() ) { G4Exception("DicomHandler ", "",JustWarning, " Have not read DCM_AcquisitionTime"); }
    //DcmDateTime::getOFDateTimeFromString(dt,dateTimeV);
    //G4cout << "DCM_AcquisitionTime " << dateTimeV << G4endl;

    if( !theDataset->findAndGetOFString(DCM_AcquisitionDate,dt).good() ) { G4Exception("DicomHandler ", "",JustWarning, " Have not read DCM_AcquisitionDate"); }
    DcmDateTime::getOFDateTimeFromString(dt,dateTimeV);
    G4cout << "DCM_AcquisitionDate: " << dateTimeV << G4endl;

    vec.push_back(dateTimeV.getDate().getYear());
    vec.push_back(dateTimeV.getDate().getMonth());
    vec.push_back(dateTimeV.getDate().getDay());
    vec.push_back(dateTimeV.getTime().getHour());
    vec.push_back(dateTimeV.getTime().getMinute());
    vec.push_back(dateTimeV.getTime().getSecond());

    //DcmIODUtil::getFloat64ValueFromItem(DCM_AcquisitionNumber, *theDataset, fv, 0);
    //if( !theDataset->findAndGetFloat64(DCM_AcquisitionNumber,fv).good() ) { G4Exception("DicomHandler ", "", JustWarning, " Have not read DCM_AcquisitionNumber"); }
    //vec.push_back(fv);

    //DcmIODUtil::getFloat64ValueFromItem(DCM_AcquisitionDuration, *theDataset, fv, 0);
    //if( !theDataset->findAndGetFloat64(DCM_AcquisitionDuration,fv).good() ) { G4Exception("DicomHandler ", "", JustWarning, " Have not read DCM_AcquisitionDuration"); }
    //vec.push_back(fv);

    for( int ii = 0; ii < vec.size(); ii++ ) {
        G4cout << " " << vec[ii];
    }
    G4cout << G4endl;

    return vec;
}

void DicomVFileImage::getDateTimeValues(G4String s , OFString dt){

    OFDateTime dateTimeV;

    DcmDateTime::getOFDateTimeFromString(dt,dateTimeV);
    G4cout << s << " :: " << dateTimeV << G4endl;

    std::vector<double> vec;
    vec.push_back(dateTimeV.getDate().getYear());
    vec.push_back(dateTimeV.getDate().getMonth());
    vec.push_back(dateTimeV.getDate().getDay());
    vec.push_back(dateTimeV.getTime().getHour());
    vec.push_back(dateTimeV.getTime().getMinute());
    vec.push_back(dateTimeV.getTime().getSecond());

    G4cout << s << " " << dateTimeV ;
    for( int ii = 0; ii < vec.size(); ii++ ) {
        G4cout << " " << vec[ii];
    }
    G4cout << G4endl;

    PETGeneralDateData[s] = vec;
}
