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

#ifndef INCLUDE_G4DoseCalcsAnalysis_HH_
#define INCLUDE_G4DoseCalcsAnalysis_HH_

//#include "globals.hh"

#include <map>
#include <vector>

#ifdef ANALYSIS_USE
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#endif
class G4DoseCalcsAnalysis {

public:

    G4DoseCalcsAnalysis(std::string dirPath);
    ~G4DoseCalcsAnalysis();

private:

    int GeometryVarDigitNum;
    int QuantityDigitNum;
    int DiffDigitNum;
    int SDevDigitNum;

    double MinValForLog ;

    double MeV_to_J;
    double Gy_to_Sv;
    double Bq_to_MBq;


    double AEUnitFactor;
    double AFUnitFactor;
    double SAFUnitFactor;
    double SUnitFactor;
    double ADUnitFactor;
    double HUnitFactor;
    double EUnitFactor;

    std::string AEUnit;
    std::string AFUnit;
    std::string SAFUnit;
    std::string SUnit;
    std::string ADUnit;
    std::string HUnit;
    std::string EUnit;

    double ResidenceTimeUnitFactor;
    std::string ResidenceTimeUnit;

    double AdministeredActivityUnitFactor;
    std::string AdministeredActivityUnit;

    std::string CombinedOutFileName;

    std::vector<int> RootC;
    std::vector<int> RootM;

    //std::vector<double> ResEnergies;

    std::map<std::string,std::map<std::string, std::map<std::string, std::vector<double>>>> ResEnergies;

    std::map<std::string,std::string> QuantityUnit;

    std::vector<std::string> GeometryRadiotracerSources;

    std::string GraphsDirectoryPath;

    std::vector<std::string> OrgansNameVector;

    std::string OrgansNamesToScoreString;
    std::string QuantitiesToScore;


    std::vector<std::string> SourceNamesToScore;
    std::vector<std::string> TargetNamesToScore;
    std::vector<std::string> QuantityNamesToScore;


    std::vector<std::string> ParticleList;
    std::vector<std::string> RadiotracerList;
    std::vector<std::string> GeometryList;

    std::map<std::string , bool> QuantityUseLog;
    std::string UseLogE ;
    std::string UseLogVariable ;
    std::string UseGridXY ;
    std::string PrintTitle ;
    std::string LegendPos ;
    double LegendXWidth ;
    double LegendYHeight ;
    double X1LegPos ;
    double X2LegPos ;
    double Y1LegPos ;
    double Y2LegPos ;
    std::string AddErrorBarInGraphs;

    std::string SourceType;
    std::string ParticleName;
    std::string SourceRegionName;
    std::string EnergyDistribution;
    std::string MomDirDistribution;

    double GaussSDev;
    double GaussMean;
    double UniformEmin;
    double UniformEmax;
    double RayleighEmax;
    double MonoEnergy;

    std::string GraphsData;
    std::string CompareType;
    std::string GraphsExt;
    std::string CompareReferenceName;
    int NumOfEne;
    int NumOfRefEne;

    std::string RefFilePath;
    std::vector<std::string> RefFilePaths;
    std::vector<std::string> CompareReferenceNames;

    std::string PositionDataFile;
    std::string EnergyDataFile;
    std::string MomDirDataFile;

    std::string SliceFor2DGraph;
    unsigned int SliceID;
    std::string BeamAxis ;
    std::string EventsDataHistograms;

    double ZMAX, YMAX, XMAX;
    double MaximalDoseVal;
    std::string DoseProfilQuantity;

    std::string GeometrySymbol;

    std::string RegionVariableName;
    std::string RegionVariableNameWithUnit;
    double EnergyGraphValue;
    std::string GenerateRegionsVariableGraph;
    std::string GenerateRelativeSDevGraph;
    std::string GenerateRelativeErrGraph;
    std::string GenerateCrossSectionGraph;

    std::string DiffExp;
    std::string DiffSym;
    std::string DifferenceMethod;

    std::map<std::string,std::map<std::string,std::map<std::string,double>>> GeometryRegionVariableValue;
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<double,double>>>>>> RelativeStandartDeviationPerCent;
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<double,double>>>>>> StandartDeviation;
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<double,double>>>>>> ResultTable;
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<double,double>>>>>> ReferenceTable;
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<double,double>>>>>> ReferenceTable2;
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<double,double>>>> ResultParticleSourceEnergyTime;
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<double,double>>>>>> ResRefErrCompTables;
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double>>>>> ReferenceQuantityGeometryRadioTracerSourceTargetValues ;
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double>>>>> ResultQuantityGeometryRadioTracerSourceTargetValues ;
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double>>>>> QuantityGeometryRadioTracerSourceTargetStandartDeviation ;
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double>>>>> QuantityGeometryRadioTracerSourceTargetRelativeStandartDeviation ;

#ifdef ANALYSIS_USE

    std::map<std::string,std::map<std::string,std::map<std::string,TGraph*>>> AllGraphs; // Result_Reference Particle Source Target Energy Value
    //std::vector<TCanvas*> GraphsVector;
#endif

    double DefaultErrorDistance ;

    unsigned int VoxZNumber;
    unsigned int VoxYNumber;
    unsigned int VoxXNumber;
    double VoxZHalfSize;
    double VoxYHalfSize;
    double VoxXHalfSize;

    //std::map<unsigned int,unsigned int> CopyNumberID;
    std::map<std::string, std::map<double ,std::map<double ,std::map<double , double>>>> QuantVals;
    std::map<unsigned int ,std::map<unsigned int ,std::map<unsigned int , double>>> QuantValsErr;

    std::map<std::string, std::map<std::string ,std::map<std::string ,std::map<double ,double>>>> ParticleMaterialProcessEneCrossSection;

    struct MatStr
    {
        double density;
        double RadiationLenght;
    };
    std::map<std::string, MatStr> MaterialsData;

    std::map<unsigned int,double> ZPos;
    std::map<unsigned int,double> YPos;
    std::map<unsigned int,double> XPos;

    std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<double,double>>>>>> ReadReferenceFile(std::string);
    void ReadResultFile();
    void GenerateResultReferenceGraphs(int);
    void GenerateResultReferenceInOneGraph();
    void GenerateResultGraphs(int);
    void GenerateResultInOneGraph();
    void GenerateResultForEachParticleEnergyInOneGraph();
    void GenerateComparisonFactorGraphs();
    void GenerateStandardDeviationGraphs();
    void DataInitialization();

    void GenerateRadioTracerResultsInOneGraphForAllComAndSrcReg();
    void GenerateRadioTracerResultsInOneGraphWithComparisonForAllComAndSrcReg();
    void GenerateRadioTracerComparisonFactorGraphsForAllComAndSrcReg();

    // for voxelized


    void CreateHeatMap();
    void CreatePercentageDepthDoseGraph();
    void CreateDoseProfile();


    // calculated
    std::map<std::string,double> AbsorbedDoseMap;
    std::map<std::string,double> EffectiveDoseMap;

    // Given from Input Factors File.
    std::map<std::string,double> TissueFactorMap;
    std::map<std::string,std::map<double,double>> RadiationFactorMap;
    std::map<std::string,std::map<std::string,std::vector<double>>> RadiotracerBiokineticDataMap;
    std::map<std::string,std::string> RadiotracerParticleMap;
    std::map<std::string,std::string> RadiotracerParticleEnergyDistMap;
    std::map<std::string,std::vector<double>> RadiotracerParticleEnergyValuesMap;
    std::string RadioTracerEmmitedParticle;
    std::string RadiationEnergyDistributionOfRadioTracer;
    double RadiationEnergyOfRadioTracer;

    // given from commands
    std::string RadioTracerName;
    std::string SAFRefForEffectiveDose;
    double TotalEffectiveDose;
    double InjectedActivity;
    double CumulatedActivity;


    bool RadiotracerDataFomFile;
    bool GenerateResultsForRadioTracer;
    bool GenerateResultsForRadioTracerExams;

    std::map<std::string,std::map<std::string,std::map<double,double>>> RadioTracerEnergyPerCent ;
    std::map<std::string,std::map<std::string,double>> RadioTracerSourceTi_Fs_ai_AsPerA0 ;
    std::map<std::string,double> RadioTracerInjectedActivity ;

    bool GenerateVoxelsResuls;
    bool IsAllTargetsToScore;
    bool IsAllSourcesToScore;

    bool V;

    struct DD{
        int val;
        double* vec;
        double* vec2;
    };

    DD vecdata;

public:

#ifdef ANALYSIS_USE
    TGraph* CreateGraph(int , double* , double*);
    TGraphErrors* CreateGraphErrors(int , double* , double* , double* , double*);
    void CreateMultiGraphParametersAndCanvas(std::string, std::string, TMultiGraph* , TLegend*);
    TGraph* setGraphData(TGraph*, int, int );
    TGraphErrors* setGraphErrorsData(TGraphErrors*, int, int );
#endif

    void GenerateGraphFromROOTGraphDATA();

    void createDataInCSVFormat();
    void createLatexTables();

    void setV(bool n ){V = n;}

    void ReadSimulationData();
    void ReadResultsAndReferenceData();
    void GenerateSelfCrossGraphs();
    double* AccumulateThirdRefDataToGraphs(std::string, std::string, std::string);
    void Remove0FromList(int,double*,double*);
    double RelativeDifferenceCalculation(double, double);
    void ReadVoxelsResult();
    void GenerateVoxel2DGraphs();
    void GenerateEventsDataHisto();
    void GenerateLatexTableResultReference();
    void GenerateLatexTableResultReferenceForOneEnergy();
    void GenerateLatexTableResultForRadioTracerGeometry();
    void GenerateLatexTablesForRadioTracerResult();
    void GenerateLatexTableResultReferenceOfQuantitiesGeometriesRadioTracers();
    void GenerateRegionsDataLatexTable();
    void GenerateSourceEnegyGraph();
    void GenerateCrossParameterGraph();
    void GenerateCrossSectionGraphs();
    void GenerateROOTFileForAllResults();

    void TestGraphGeneration();

    void GenerateSourceCrossEffectGraphs(int);

    void GenerateSourceComputationTimeGraph();

};





#endif /* INCLUDE_G4DoseCalcsAnalysis_HH_ */
