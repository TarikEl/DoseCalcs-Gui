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

#include "G4TMessenger.hh"
#include "G4TVolumeConstruction.hh"
//#include "G4TPrimaryGeneratorAction.hh"
//#include "G4TRunAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "globals.hh"

#include "G4Tokenizer.hh"
#include "G4SystemOfUnits.hh"

extern bool DirectoryExists( const char* pzPath );

G4TMessenger::G4TMessenger(G4TVolumeConstruction* myUsrPhtm):VolumeConst(myUsrPhtm){

    RunAndScoreDir = new G4UIdirectory("/RunAndScoreData/");
    RunAndScoreDir->SetGuidance("Set the quantities to score and target volumes, values related to the execution modes");

    CommandsForRunAndScore();
    CommandsForAnalysis();
}

G4double G4TMessenger::UseG4Units(G4String Unit){

    G4double U = 1;

    if(Unit == "mm"){U = mm;}
    else if(Unit == "cm"){U = cm;}
    else if(Unit == "m"){U = m;}
    else if(Unit == "g/cm3"){U = g/cm3;}
    else if(Unit == "mg/cm3"){U = mg/cm3;}
    else if(Unit == "mg/mm3"){U = mg/mm3;}
    else if(Unit == "kg/m3"){U = kg/m3;}
    else if(Unit == "mg/mL"){U = mg/mL;}
    else if(Unit == "g/mL"){U = g/mL;}
    else if(Unit == "degree"){U = degree;}
    else if(Unit == "radian"){U = radian;}
    else if(Unit == "eV"){U = eV;}
    else if(Unit == "keV"){U = keV;}
    else if(Unit == "MeV"){U = MeV;}

    return  U;
}


G4TMessenger::~G4TMessenger()
{
    delete number_of_threadsCMD ;
    delete SimNumOnRanksCMD;

    delete analysisDataDir ;
    /*
    delete Graphs_DataCMD ;
    delete Compare_typeCMD ;
    delete Graphs_ExtCMD ;
    delete Ref_File_PathCMD ;
    delete Ref_NameCMD;
    delete Organs_to_scoreCMD ;
    delete Variable_To_ScoreCMD ;
    delete Num_Of_EneCMD ;
    delete Num_Of_RefEneCMD ;
    delete EventsPositionHistogramCMD;
    delete EventsEnergyHistogramCMD;
    delete EventsMomDirHistogramCMD;
    delete number_events_per_threadCMD ;

    delete SliceIDCMD ;
    delete SliceFor2DGraphCMD;
    delete BeamAxisCMD;

    */

    delete EventsDataHistogramsCMD;

    delete GenerateSelfCrossGraphs;
    delete GenerateRegionsVariableGraph;
    delete GenerateRelativeErrGraphCMD;
    delete GenerateRelativeSDevGraphCMD;
    delete GenerateVoxelizedHist;

    delete SetGraphsParameters;

    delete RegionToScoreCMD;
    delete QuantitiesToScoreCMD;
    delete GenerateVoxelsResulsCMD;
    delete SimulationIntExtNeutDet;

    delete RadiationFactorsCMD;
    delete QuantitiesUnitsCMD;
    delete TissueFactorsCMD;
    delete RadioTracerBiokineticCMD;
    delete RadioTracerDataCMD;

    delete ResultDirectoryPathCMD;
}


// called automatically for each command in the begining of running to set the values entered by the user in Macro file of by typing commands from
// it call G4TVolumeConstruction::SetPhantomModel(G4String newModel), G4TVolumeConstruction::SetPhantomSex(G4String newSex)
// it call G4TMessenger::AddBodyPart(G4String newBodyPartSensitivity)
void G4TMessenger::SetNewValue(G4UIcommand* command,G4String newValue){


    if( command == RegionToScoreCMD)
    {
        G4Tokenizer next(newValue);

        G4String nn = "";
        for(G4int ds = 0 ; ds < 100 ; ds++){

            if( G4String ss = next()){
                if( ss == "" || ss.empty()){
                    break;
                }
                //G4cout << "\n\n\n\n ss = " << ss << G4endl ;

                nn += ss + " ";
            }
        }
        //G4cout << "\n\n\n\n nn = " << nn << G4endl ;

        VolumeConst->setOrgans_to_score(nn);
    }

    if( command == QuantitiesToScoreCMD)
    {
        G4Tokenizer next(newValue);

        G4String nn = "";
        for(G4int ds = 0 ; ds < 7 ; ds++){

            if( G4String ss = next()){
                if( ss == "" || ss.empty()){
                    break;
                }

                //G4cout << "\n\n\n\n ss = " << ss << G4endl ;

                nn += ss + " ";
            }
        }
        //G4cout << "\n\n\n\n nn = " << nn << G4endl ;
        VolumeConst->setVariable_To_Score(newValue);

    }

    if( command == SetGraphsParameters)
    {
        G4Tokenizer next(newValue);

        VolumeConst->setUseLogE(next());
        VolumeConst->setUseLogVariable(next());
        VolumeConst->setUseGridXY(next());
        VolumeConst->setPrintTitle(next());
        VolumeConst->setLegendPos(next());
        VolumeConst->setLegendXWidth(StoD(next()));
        VolumeConst->setLegendYHeight(StoD(next()));
        VolumeConst->setAddErrorBarInGraphs(next());
        VolumeConst->setGraphs_Ext(next());

        /*
        G4String nn = "";
        for(G4int ds = 0 ; ds < 6 ; ds++){

            if( G4String ss = next()){
                if( ss == "" || ss.empty()){
                    break;
                }
                if( ds == 5){
                    VolumeConst->setGraphs_Ext(next());
                }

                //G4cout << "\n\n\n\n ss = " << ss << G4endl ;

                nn += ss + " ";
            }
        }
        */
    }

    if( command == GenerateSelfCrossGraphs)
    {
        G4Tokenizer next(newValue);

        VolumeConst->setGraphs_Data(next());
        VolumeConst->setCompare_type(next());

        G4String Par = next();
        G4String RefNM = "";
        G4String RefPH = "";

        //std::cout<< " Dddddd " << SN <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
        while (!Par.empty()) {


            RefNM += Par + " ";
            RefPH += next() + " ";;

            //std::cout << RefNM << " " << RefPH <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;

            Par = next();
        }

        VolumeConst->setRef_Name(RefNM);
        VolumeConst->setRef_File_Path(RefPH);
    }


    if( command == GenerateRegionsVariableGraph)
    {
        G4Tokenizer next(newValue);
        VolumeConst->setRegionVariableName(next());
        VolumeConst->setGenerateRegionsVariableGraph("yes");
        //VolumeConst->setEnergyGraphValue(StoD(next())*UseG4Units(next()));
    }

    if( command == GenerateRelativeErrGraphCMD )
    {
        G4String DifferenceMethod = "RD";
        G4Tokenizer next(newValue);
        if(DifferenceMethod == next()){}

        VolumeConst->setGenerateRelativeErrGraph("yes");
        VolumeConst->setDifferenceMethod(DifferenceMethod);
    }
    if( command == GenerateRelativeSDevGraphCMD )
    {
        VolumeConst->setGenerateRelativeSDevGraph("yes");
    }
    if( command == EventsDataHistogramsCMD )
    {
        G4Tokenizer next(newValue);

        G4String Data = "", dd = next();
        while (!dd.empty()) {

            Data += dd;Data += " ";
            dd = next();
        }

        VolumeConst->setEventsDataHistograms(Data);
    }

    if( command == ResultDirectoryPathCMD )
    {
        //std::cout<< "\n\n\n\n newValue : " << newValue <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
        //std::cout<< "\n\n\n\n DirectoryExists(newValue) : " << DirectoryExists(newValue) <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;

        if(DirectoryExists(newValue)){
            //std::cout<< "\n\n\n\n yees new path : " << newValue <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
            VolumeConst->setResultDirectoryPath(newValue);
        }
    }

    if( command == GenerateVoxelizedHist)
    {
        G4Tokenizer next(newValue);
        VolumeConst->setDoseProfilQuantity(next());
        VolumeConst->setBeamAxis(next());
        VolumeConst->setSliceFor2DGraph(next());
        VolumeConst->setSliceID(StoI(next()));

        //VolumeConst->setEnergyGraphValue(StoD(next())*UseG4Units(next()));
    }

    if( command == RadioTracerDataCMD)
    {
        G4Tokenizer next(newValue);

        G4String Data = "", dd = next();
        while (!dd.empty()) {

            Data += dd;Data += " ";
            dd = next();
        }

        //std::cout<< "\n\n\n\n Radio Tracer Data : " << Data <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
        VolumeConst->setRadioTracerData(Data);

    }

    if( command == RadioTracerBiokineticCMD)
    {
        G4Tokenizer next(newValue);

        G4String Data = "", dd = next();
        while (!dd.empty()) {

            Data += dd;Data += " ";
            dd = next();
        }

        //std::cout<< "\n\n\n\n Radio Tracer Data : " << Data <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
        VolumeConst->setRadioTracerBiokinetic(Data);

    }

    if( command == RadiationFactorsCMD)
    {
        G4Tokenizer next(newValue);

        G4String Data = "", dd = next();
        while (!dd.empty()) {

            Data += dd;Data += " ";
            dd = next();
        }

        //std::cout<< "\n\n\n\n Radio Tracer Data : " << Data <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
        VolumeConst->setRadiationFactors(Data);

    }

    if( command == QuantitiesUnitsCMD)
    {
        G4Tokenizer next(newValue);

        G4String Data = "", dd = next();
        while (!dd.empty()) {

            Data += dd;Data += " ";
            dd = next();
        }

        //std::cout<< "\n\n\n\n Quantities Units : " << Data <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
        VolumeConst->setQuantitiesUnits(Data);

    }
    if( command == TissueFactorsCMD)
    {
        G4Tokenizer next(newValue);

        G4String Data = "", dd = next();
        while (!dd.empty()) {

            Data += dd;Data += " ";
            dd = next();
        }

        //std::cout<< "\n\n\n\n Radio Tracer Data : " << Data <<  " @@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
        VolumeConst->setTissueFactors(Data);
    }

    if( command == SimNumOnRanksCMD )
    {
        VolumeConst->setMPISimulationNum(newValue);
    }
    // Execution commands
    if( command == number_of_threadsCMD )
    {
        VolumeConst->setNumberOfThreads(number_of_threadsCMD->GetNewIntValue(newValue));
    }
    if( command == AccuracyCalculationLevelCMD )
    {
        VolumeConst->setAccuracyCalculationLevel(newValue);
    }

    if( command == SimulationIntExtNeutDet )
    {
        G4Tokenizer next(newValue);

        VolumeConst->setSimulationIntExtNeutDet(next());
    }
    if( command == GenerateVoxelsResulsCMD )
    {
        VolumeConst->setGenerateVoxelsResuls("yes");
    }

}


// called from constructor to define the commands waited to fill from user
void G4TMessenger::CommandsForRunAndScore(){

    G4UIparameter* param;

    ResultDirectoryPathCMD = new G4UIcmdWithAString("/RunAndScoreData/setResultDirectoryPath",this);
    ResultDirectoryPathCMD->SetGuidance("");
    //ResultDirectoryPathCMD->SetParameterName("ResultDirectoryName",true);
    //ResultDirectoryPathCMD->SetDefaultValue("yes");
    ResultDirectoryPathCMD->AvailableForStates(G4State_PreInit,G4State_Idle);


    RegionToScoreCMD    = new G4UIcommand("/RunAndScoreData/setVolumesToScore",this);
    for(G4int ds = 0 ; ds < 100 ; ds++){

        G4String f = "Regions where the quantities will be scored " + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      RegionToScoreCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

    QuantitiesToScoreCMD = new G4UIcommand("/RunAndScoreData/setQuantitiesToScore",this);
    for(G4int ds = 0 ; ds < 7 ; ds++){

        G4String f = "the physical quantities tha will be scored" + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      QuantitiesToScoreCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

    number_of_threadsCMD = new G4UIcmdWithAnInteger("/RunAndScoreData/setNumberOfThreads",this);
    number_of_threadsCMD->SetGuidance("Set number of workers threads)");
    number_of_threadsCMD->SetParameterName("number_of_threadsCMD",true);
    number_of_threadsCMD->SetDefaultValue(2);
    number_of_threadsCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    AccuracyCalculationLevelCMD = new G4UIcmdWithAString("/RunAndScoreData/setAccuracyCalculationLevel",this);
    AccuracyCalculationLevelCMD->SetGuidance("Set the level at which we will estimate the simulation accuracy");
    AccuracyCalculationLevelCMD->SetParameterName("AccuracyCalculationLevelCMD",true);
    AccuracyCalculationLevelCMD->SetDefaultValue("EventLevel");
    AccuracyCalculationLevelCMD->SetCandidates(" BatchLevel EventLevel StepLevel");
    AccuracyCalculationLevelCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    SimNumOnRanksCMD = new G4UIcmdWithAString("/RunAndScoreData/setSimNumOnRanks",this);
    SimNumOnRanksCMD->SetGuidance("Set the simulation number on ranks, o for one simulation per ranks, m for multi-simulations on ranks, which means, one simulation per rank ");
    SimNumOnRanksCMD->SetParameterName("SimNumOnRanksCMD",true);
    SimNumOnRanksCMD->SetDefaultValue("o");
    SimNumOnRanksCMD->SetCandidates(" o m ");
    SimNumOnRanksCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    SimulationIntExtNeutDet = new G4UIcommand("/RunAndScoreData/RunFor" ,this);
    param = new G4UIparameter("RunFor",'s', false);
    SimulationIntExtNeutDet->SetParameter(param);


    GenerateVoxelsResulsCMD = new G4UIcmdWithoutParameter("/RunAndScoreData/generateVoxelsResults",this);
    GenerateVoxelsResulsCMD->SetGuidance("set the command without parameter to generate results at the voxel level");
    GenerateVoxelsResulsCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    RadiationFactorsCMD = new G4UIcommand("/RunAndScoreData/setRadiationFactors" ,this);
    for(G4int ds = 0 ; ds < 150 ; ds++){

        G4String f = "RadiationFactors-" + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      RadiationFactorsCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

    QuantitiesUnitsCMD = new G4UIcommand("/RunAndScoreData/setQuantitiesUnits" ,this);
    for(G4int ds = 0 ; ds < 40 ; ds++){

        G4String f = "QuantitiesUnits-" + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      QuantitiesUnitsCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

    TissueFactorsCMD = new G4UIcommand("/RunAndScoreData/setTissueFactors" ,this);
    for(G4int ds = 0 ; ds < 150 ; ds++){

        G4String f = "TissueFactors-" + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      TissueFactorsCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

    RadioTracerDataCMD = new G4UIcommand("/RunAndScoreData/setRadioTracerData" ,this);
    for(G4int ds = 0 ; ds < 150 ; ds++){

        G4String f = "RadioTracer-" + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      RadioTracerDataCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

    RadioTracerBiokineticCMD = new G4UIcommand("/RunAndScoreData/setRadioTracerBiokinetic" ,this);
    for(G4int ds = 0 ; ds < 150 ; ds++){

        G4String f = "RadioTracerBiokinetics-" + std::to_string(ds);
        param = new G4UIparameter(f,'s', true);      RadioTracerBiokineticCMD->SetParameter(param); // true because if we sent just 11 fpr ex it will not throw an error
    }

}

// called from constructor to define the commands waited to fill from user
void  G4TMessenger::CommandsForAnalysis(){

    analysisDataDir = new G4UIdirectory("/AnalysisData/");
    analysisDataDir->SetGuidance("Set Your analysis data");

    G4UIparameter* param;

    GenerateSelfCrossGraphs = new G4UIcommand("/AnalysisData/generateSelfCrossGraphs" ,this);
    param = new G4UIparameter("Graph Data",'s', false);   GenerateSelfCrossGraphs->SetParameter(param);
    param = new G4UIparameter("Compare Data",'s', false);   GenerateSelfCrossGraphs->SetParameter(param);
    param = new G4UIparameter("Ref Name",'s', true);   GenerateSelfCrossGraphs->SetParameter(param);
    param = new G4UIparameter("Ref File Path",'s', true);   GenerateSelfCrossGraphs->SetParameter(param);

    GenerateRegionsVariableGraph = new G4UIcommand("/AnalysisData/generateVariableRegionGraph" ,this);
    param = new G4UIparameter("Variable Name",'s', false);   GenerateRegionsVariableGraph->SetParameter(param);


    GenerateVoxelizedHist = new G4UIcommand("/AnalysisData/generateVoxelizedHistograms" ,this);
    param = new G4UIparameter("Beam Axis",'s', true);   GenerateVoxelizedHist->SetParameter(param);
    param = new G4UIparameter("Dose Quantity",'s', true);   GenerateVoxelizedHist->SetParameter(param);
    param = new G4UIparameter("Slice Axis",'s', true);   GenerateVoxelizedHist->SetParameter(param);
    param = new G4UIparameter("Slice ID",'s', true);   GenerateVoxelizedHist->SetParameter(param);

    GenerateRelativeErrGraphCMD = new G4UIcommand("/AnalysisData/generateRelativeErrGraph",this);
    param = new G4UIparameter("Difference Type",'s', true);   GenerateRelativeErrGraphCMD->SetParameter(param);

    GenerateRelativeErrGraphCMD->SetGuidance("");
    //GenerateRelativeErrGraphCMD->SetParameterName("GenerateRelativeErrGraph",true);
    //GenerateRelativeErrGraphCMD->SetDefaultValue("yes");
    GenerateRelativeErrGraphCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    GenerateRelativeSDevGraphCMD = new G4UIcmdWithoutParameter("/AnalysisData/generateRelativeSDevGraph",this);
    GenerateRelativeSDevGraphCMD->SetGuidance("");
    //GenerateRelativeSDevGraphCMD->SetParameterName("GenerateRelativeSDevGraph",true);
    //GenerateRelativeSDevGraphCMD->SetDefaultValue("yes");
    GenerateRelativeSDevGraphCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    EventsDataHistogramsCMD = new G4UIcmdWithoutParameter("/AnalysisData/generateEventsDataHisto",this);
    EventsDataHistogramsCMD->SetGuidance("");
    //EventsDataHistogramsCMD->SetParameterName("EventsDataHistograms",true);
    //EventsDataHistogramsCMD->SetDefaultValue("yes");
    EventsDataHistogramsCMD->AvailableForStates(G4State_PreInit,G4State_Idle);


    SetGraphsParameters = new G4UIcommand("/AnalysisData/setGraphsParameters" ,this);
    param = new G4UIparameter("UseLogE",'s', true);   SetGraphsParameters->SetParameter(param);
    param = new G4UIparameter("UseLogVariable",'s', true);   SetGraphsParameters->SetParameter(param);
    param = new G4UIparameter("UseGridXY",'s', true);   SetGraphsParameters->SetParameter(param);
    param = new G4UIparameter("PrintTitle",'s', true);   SetGraphsParameters->SetParameter(param);
    param = new G4UIparameter("LegendPos",'s', true);   SetGraphsParameters->SetParameter(param);
    param = new G4UIparameter("LegendXWidth",'s', true);   SetGraphsParameters->SetParameter(param);
    param = new G4UIparameter("LegendYHeight",'s', true);   SetGraphsParameters->SetParameter(param);
    param = new G4UIparameter("GraphExt",'s', true);   SetGraphsParameters->SetParameter(param);

}

