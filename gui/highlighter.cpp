#include "gui/highlighter.h"
#include "QRegularExpression"

Highlighter::Highlighter(QTextDocument *parent) : QSyntaxHighlighter(parent) {

    HighlightingRule rule;

    singleLineCommentFormat.setForeground(Qt::gray);
    rule.pattern = QRegularExpression(QStringLiteral("#[^\n]*"));
    //rule.pattern = QRegularExpression(QStringLiteral("//[^\n]*"));
    rule.format = singleLineCommentFormat;
    highlightingRules.append(rule);

    //NumPointsFormat.setFontItalic(true);
    NumPointsFormat.setFontWeight(QFont::Bold);
    NumPointsFormat.setForeground(Qt::blue);
    //rule.pattern = QRegularExpression(QStringLiteral("\\b^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\\b"));
    rule.pattern = QRegularExpression(QStringLiteral("\\b[-+]?[0-9]*\.?[0-9]*?[eE]?[-+]?[0-9]\\b"));
    rule.format = NumPointsFormat;
    highlightingRules.append(rule);

    const QString DoseCalcsConstValuesPatterns[] = {

        QStringLiteral("\\bVOXEL\\b"), QStringLiteral("\\bVoxIDs\\b"),
        QStringLiteral("\\bGDML\\b"), QStringLiteral("\\bSTL\\b"),
        QStringLiteral("\\bC\+\+\\b"), QStringLiteral("\\bTEXT\\b"),
        QStringLiteral("\\bDICOM\\b"),

        QStringLiteral("\\b.gdml\\b"), QStringLiteral("\\b.geom\\b"),
        QStringLiteral("\\b.ast\\b"), QStringLiteral("\\b.stl\\b"),
        QStringLiteral("\\b.c\+\+\\b"),

        QStringLiteral("\\bxy\\b"), QStringLiteral("\\byz\\b"),
        QStringLiteral("\\byx\\b"), QStringLiteral("\\bzy\\b"),
        QStringLiteral("\\bxz\\b"), QStringLiteral("\\bzx\\b"),

        QStringLiteral("\\bBox\\b"), QStringLiteral("\\bTubs\\b"),
        QStringLiteral("\\bCutTubs\\b"), QStringLiteral("\\bCons\\b"),
        QStringLiteral("\\bPara\\b"), QStringLiteral("\\bTrd\\b"),
        QStringLiteral("\\bSphere\\b"), QStringLiteral("\\bOrb\\b"),
        QStringLiteral("\\bTorus\\b"), QStringLiteral("\\bEllipsoid\\b"),
        QStringLiteral("\\bUnion\\b"), QStringLiteral("\\bIntersection\\b"),
        QStringLiteral("\\bSubtraction\\b"),

        QStringLiteral("\\bLivermore\\b"), QStringLiteral("\\bPenelope\\b"),
        QStringLiteral("\\bEMS\\b"), QStringLiteral("\\bEMS1\\b"),
        QStringLiteral("\\bEMS2\\b"), QStringLiteral("\\bEMS3\\b"),
        QStringLiteral("\\bEMS4\\b"), QStringLiteral("\\bConstruct\\b"),

        QStringLiteral("\\bnumb\\b"), QStringLiteral("\\bfrac\\b"),
        QStringLiteral("\\bregions\\b"),

        QStringLiteral("\\bVoxels\\b"), QStringLiteral("\\bVolume\\b"),QStringLiteral("\\bTET\\b"),
        QStringLiteral("\\bPoint\\b"), QStringLiteral("\\bBeam\\b"),
        QStringLiteral("\\bSurface\\b"), QStringLiteral("\\bPlane\\b"), QStringLiteral("\\bSolid\\b"),
        QStringLiteral("\\bAllRegions\\b"),QStringLiteral("\\ballregions\\b"),QStringLiteral("\\bnull\\b"),

        QStringLiteral("\\bMass\\b"), QStringLiteral("\\bDensity\\b"), QStringLiteral("\\bVolume\\b"),

        QStringLiteral("\\bMono\\b"), QStringLiteral("\\bRayleigh\\b"),
        QStringLiteral("\\bFile\\b"), QStringLiteral("\\bUniform\\b"),
        QStringLiteral("\\bGauss\\b"), QStringLiteral("\\bSpectrum\\b"),

        QStringLiteral("\\bIsotropic\\b"), QStringLiteral("\\bUniform\\b"),
        QStringLiteral("\\bDirected\\b"),

        QStringLiteral("\\bgamma\\b"), QStringLiteral("\\be-\\b"),
        QStringLiteral("\\be+\\b"), QStringLiteral("\\bproton\\b"),
        QStringLiteral("\\balpha\\b"),QStringLiteral("\\bneutron\\b"),

        QStringLiteral("\\bSAF\\b"), QStringLiteral("\\bAF\\b"),
        QStringLiteral("\\bAE\\b"), QStringLiteral("\\bAD\\b"),
        QStringLiteral("\\bS\\b"), QStringLiteral("\\bAD\\b"),
        QStringLiteral("\\bDR\\b"), QStringLiteral("\\bE\\b"),
        QStringLiteral("\\bH\\b"),

        QStringLiteral("\\bsource\\b"), QStringLiteral("\\btarget\\b"),

        QStringLiteral("\\ball\\b"), QStringLiteral("\\bAll\\b"),
        QStringLiteral("\\bm\\b"), QStringLiteral("\\bo\\b"),

        QStringLiteral("\\byes\\b"), QStringLiteral("\\bno\\b"),

        QStringLiteral("\\bResult\\b"), QStringLiteral("\\bReference_Result\\b"),
        QStringLiteral("\\bSelf\\b"), QStringLiteral("\\bSelf_Cross\\b"),
        QStringLiteral("\\bCross\\b"), QStringLiteral("\\bRA\\b"),
        QStringLiteral("\\bRD\\b"), QStringLiteral("\\bLRD\\b"),

        QStringLiteral("\\bRightBottom\\b"), QStringLiteral("\\bLeftBottom\\b"),
        QStringLiteral("\\bRightTop\\b"), QStringLiteral("\\bLeftTop\\b"),
        QStringLiteral("\\bMiddleBottom\\b"), QStringLiteral("\\bMiddleTop\\b"),

        QStringLiteral("\\bXY\\b"), QStringLiteral("\\bYZ\\b"),
        QStringLiteral("\\bYX\\b"), QStringLiteral("\\bZY\\b"),
        QStringLiteral("\\bZX\\b"), QStringLiteral("\\bXZ\\b"),

        QStringLiteral("\\b.root\\b"), QStringLiteral("\\b.pdf\\b"),
        QStringLiteral("\\b.ps\\b"), QStringLiteral("\\b.png\\b"),
        QStringLiteral("\\b.jpeg\\b"),

    };

    ConstValuesFormat.setFontWeight(QFont::Bold);
    ConstValuesFormat.setForeground(Qt::darkYellow);
    for (const QString &pattern : DoseCalcsConstValuesPatterns) {
        rule.pattern = QRegularExpression(pattern);
        rule.format = ConstValuesFormat;
        highlightingRules.append(rule);
    }

    const QString UnitsPatterns[] = {
        QStringLiteral("\\bs\\b"),
        QStringLiteral("\\bmin\\b"),
        QStringLiteral("\\bh\\b"),
        QStringLiteral("\\bd\\b"),
        QStringLiteral("\\by\\b"),
        QStringLiteral("\\bmm\\b"),
        QStringLiteral("\\bcm\\b"),
        QStringLiteral("\\bm\\b"),
        QStringLiteral("\\bg/cm3\\b"),
        QStringLiteral("\\bmg/cm3\\b"),
        QStringLiteral("\\bmg/mm3\\b"),
        QStringLiteral("\\bkg/m3\\b"),
        QStringLiteral("\\bmg/mL\\b"),
        QStringLiteral("\\bg/mL\\b"),
        QStringLiteral("\\bdegree\\b"),
        QStringLiteral("\\bradian\\b"),
        QStringLiteral("\\beV\\b"),
        QStringLiteral("\\bkeV\\b"),
        QStringLiteral("\\bMeV\\b"),
        QStringLiteral("\\bGeV\\b"),
        QStringLiteral("\\bBq\\b"),
        QStringLiteral("\\bkBq\\b"),
        QStringLiteral("\\bMBq\\b"),
        QStringLiteral("\\bGBq\\b"),
        QStringLiteral("\\bJ\\b"), QStringLiteral("\\bkg-1\\b"), QStringLiteral("\\bg-1\\b"),
        QStringLiteral("\\bMeV/kg\\b"), QStringLiteral("\\bGy\\b"), QStringLiteral("\\bmGy\\b"),
        QStringLiteral("\\bMGy\\b"), QStringLiteral("\\bmiGy\\b"), QStringLiteral("\\bnGy\\b"),
        QStringLiteral("\\bSv\\b"), QStringLiteral("\\bmSv\\b"),
    };

    UnitsFormat.setFontWeight(QFont::Bold);
    UnitsFormat.setForeground(Qt::darkRed);
    for (const QString &pattern : UnitsPatterns) {
        rule.pattern = QRegularExpression(pattern);
        rule.format = UnitsFormat;
        highlightingRules.append(rule);
    }

    const QString DoseCalcsDirPatterns[] = {
        QStringLiteral("\\bGeometryData\\b"), QStringLiteral("\\bMaterialData\\b"),
        QStringLiteral("\\bSourceData\\b"), QStringLiteral("\\bRunAndScoreData\\b"),
        QStringLiteral("\\bAnalysisData\\b"), QStringLiteral("\\bPhysicsData\\b"),
    };

    DirectoryCommandFormat.setFontWeight(QFont::Bold);
    DirectoryCommandFormat.setForeground(Qt::darkGreen);
    for (const QString &pattern : DoseCalcsDirPatterns) {
        rule.pattern = QRegularExpression(pattern);
        rule.format = DirectoryCommandFormat;
        highlightingRules.append(rule);
    }

    const QString DoseCalcsCommandsPatterns[] = {
        QStringLiteral("\\bcreateElement\\b"),
        QStringLiteral("\\bcreateMaterial\\b"),
        QStringLiteral("\\baddElements\\b"),
        QStringLiteral("\\baddMaterial\\b"),
        QStringLiteral("\\bsetNistMaterialNameAndID\\b"),
        QStringLiteral("\\bcreateWorld\\b"),
        QStringLiteral("\\bcreateSolid\\b"),
        QStringLiteral("\\bcreateVolume\\b"),
        QStringLiteral("\\bsetTETPhantomLimits\\b"),

        QStringLiteral("\\bsetMaterialNameAsRegionName\\b"),
        QStringLiteral("\\bsetTETRegionData\\b"),
        QStringLiteral("\\bsetGeometrySymbol\\b"),
        QStringLiteral("\\bsetVolumesVisToNotForced\\b"),
        QStringLiteral("\\bsetVoxelsData\\b"),
        QStringLiteral("\\bsetVoxContainerPos\\b"),
        QStringLiteral("\\bsetVoxContainerRot\\b"),
        QStringLiteral("\\bsetVoxDefaultMaterialName\\b"),
        QStringLiteral("\\bsetVoxelizedRegionData\\b"),
        QStringLiteral("\\bsetCtDensityValues\\b"),
        QStringLiteral("\\bsetMatNumberAndMaterials\\b"),
        QStringLiteral("\\bsetDcmTypeAndDataPath\\b"),
        QStringLiteral("\\bsetDcmMaterialName\\b"),
        QStringLiteral("\\bvisualizeVoxelizedPlanes\\b"),
        QStringLiteral("\\bcreateSolid\\b"),
        QStringLiteral("\\bcreateVolume\\b"),
        QStringLiteral("\\bsetPhysicsData\\b"),
        QStringLiteral("\\bsetCutsInRange\\b"),
        QStringLiteral("\\bsetEnergyRange\\b"),
        QStringLiteral("\\bgenerateCrossSectionFor\\b"),
        QStringLiteral("\\bsetEventsParticleNameData\\b"),
        QStringLiteral("\\bsetEventsInitialPosData\\b"),
        QStringLiteral("\\bsetEventsInitialEneData\\b"),
        QStringLiteral("\\bsetEventsInitialMomDirData\\b"),
        QStringLiteral("\\buseDataGenerationFiles\\b"),
        QStringLiteral("\\btestEventsInitialPositions\\b"),
        QStringLiteral("\\bshowSourceBox\\b"),
        QStringLiteral("\\bsetVolumesToScore\\b"),
        QStringLiteral("\\bsetQuantitiesToScore\\b"),
        QStringLiteral("\\bsetNumberOfThreads\\b"),
        QStringLiteral("\\bsetAccuracyCalculationLevel\\b"),
        QStringLiteral("\\bsetQuantitiesUnits\\b"),
        QStringLiteral("\\bsetRadiationFactors\\b"),
        QStringLiteral("\\bsetEventNumberPerThread\\b"),
        QStringLiteral("\\bsetSimNumOnRanks\\b"),
        QStringLiteral("\\bsetRadioTracerData\\b"),
        QStringLiteral("\\bsetTissueFactors\\b"),
        QStringLiteral("\\bsetResultDirectoryPath\\b"),
        QStringLiteral("\\bsetRadioTracerBiokinetic\\b"),
        QStringLiteral("\\bgenerateSelfCrossGraphs\\b"),
        QStringLiteral("\\bgenerateRelativeErrGraph\\b"),
        QStringLiteral("\\bgenerateRelativeSDevGraph\\b"),
        QStringLiteral("\\bgenerateVariableRegionGraph\\b"),
        QStringLiteral("\\bgenerateEventsDataHisto\\b"),
        QStringLiteral("\\bsetSliceFor2DGraph\\b"),
        QStringLiteral("\\bsetBeamAxis\\b"),
        QStringLiteral("\\bsetSliceID\\b"),
        QStringLiteral("\\bsetGraphsParameters\\b"),
        QStringLiteral("\\b*****\\b"),
        QStringLiteral("\\bDWITH_GDML_USE*\\b"),
        QStringLiteral("\\bDWITH_MPI_USE*\\b"),
        QStringLiteral("\\bDWITH_ANALYSIS_USE*\\b"),
        QStringLiteral("\\bDROOT_DIR*\\b"),
        QStringLiteral("\\bDWITH_DCMTK_USE*\\b"),
        QStringLiteral("\\bDDCMTK_DIR*\\b"),
        QStringLiteral("\\bDWITH_VERBOSE_USE\\b"),
        QStringLiteral("\\bDCMAKE_CXX_COMPILER*\\b"),
        QStringLiteral("\\bDCMAKE_C_COMPILER*\\b"),
        QStringLiteral("\\b\ON\\b"),
        QStringLiteral("\\b\OFF\\b"),
        QStringLiteral("\\bCMAKE_INSTALL_DIR\\b"),
        QStringLiteral("\\bGEANT4_INSTALL_DIR\\b"),
        QStringLiteral("\\bMPI_INSTALL_DIR\\b"),
        QStringLiteral("\\bROOT_INSTALL_DIR\\b"),
        QStringLiteral("\\bDCMTK_INSTALL_DIR\\b"),
        QStringLiteral("\\bDoseCalcs_SOURCE_DIR\\b"),
        QStringLiteral("\\bDEFAULT_DoseCalcs_INPUTS\\b"),
        QStringLiteral("\\bCMAKE_DOWNLOAD_URL\\b"),
        QStringLiteral("\\bGEANT4_DOWNLOAD_URL\\b"),
        QStringLiteral("\\bXERCES_DOWNLOAD_URL\\b"),
        QStringLiteral("\\bMPI_DOWNLOAD_URL\\b"),
        QStringLiteral("\\bROOT_DOWNLOAD_URL\\b"),
        QStringLiteral("\\bDCMTK_DOWNLOAD_URL\\b"),
        QStringLiteral("\\bVERBOSE_USE\\b"),
        QStringLiteral("\\bGDML_USE\\b"),
        QStringLiteral("\\bCENTOS_ROCKS_CLUSTER\\b"),
        QStringLiteral("\\bMPI_USE\\b"),
        QStringLiteral("\\bROOT_USE\\b"),
        QStringLiteral("\\bDCMTK_USE\\b"),
    };

    CommandFormat.setFontWeight(QFont::Bold);
    CommandFormat.setForeground(Qt::darkBlue);
    for (const QString &pattern : DoseCalcsCommandsPatterns) {
        rule.pattern = QRegularExpression(pattern);
        rule.format = CommandFormat;
        highlightingRules.append(rule);
    }


    const QString ConfigPatterns[] = {

    };

    ConfigFormat.setFontWeight(QFont::Bold);
    ConfigFormat.setForeground(Qt::darkRed);
    for (const QString &pattern : ConfigPatterns) {
        rule.pattern = QRegularExpression(pattern);
        rule.format = ConfigFormat;
        highlightingRules.append(rule);
    }

    const QString GdmlPatterns[] = {

        // for gdml syntax

        QStringLiteral("\\b?\\b"), QStringLiteral("\\b!\\b"), QStringLiteral("\\b]\\b"), QStringLiteral("\\b]\\b"),
        QStringLiteral("\\b<\\b"), QStringLiteral("\\b>\\b"), QStringLiteral("\\b\"\\b"), QStringLiteral("\\b=\\b"),
        QStringLiteral("\\b/\\b"), QStringLiteral("\\b&\\b"), QStringLiteral("\\b:\\b"),QStringLiteral("\\b;\\b"),
        QStringLiteral("\\b,\\b"),

        QStringLiteral("\\bdefine\\b"), QStringLiteral("\\bposition\\b"), QStringLiteral("\\bname\\b"),
        QStringLiteral("\\bx\\b"), QStringLiteral("\\by\\b"), QStringLiteral("\\bz\\b"),
        QStringLiteral("\\bxml\\b"), QStringLiteral("\\bversion\\b"), QStringLiteral("\\bencoding\\b"),
        QStringLiteral("\\bsolid\\b"), QStringLiteral("\\bstructure\\b"), QStringLiteral("\\bDOCTYPE\\b"),
        QStringLiteral("\\bgdml\\b"), QStringLiteral("\\bENTITY\\b"), QStringLiteral("\\bSYSTEM\\b"),
        QStringLiteral("\\bvolume\\b"), QStringLiteral("\\bmaterialref\\b"), QStringLiteral("\\bsolidref\\b"),
        QStringLiteral("\\bref\\b"), QStringLiteral("\\bsetup\\b"), QStringLiteral("\\bworld\\b"),
        QStringLiteral("\\brotation\\b"),


        // for c++ syntax
        QStringLiteral("\\bchar\\b"), QStringLiteral("\\bclass\\b"), QStringLiteral("\\bconst\\b"),
        QStringLiteral("\\bdouble\\b"), QStringLiteral("\\benum\\b"), QStringLiteral("\\bexplicit\\b"),
        QStringLiteral("\\bfriend\\b"), QStringLiteral("\\binline\\b"), QStringLiteral("\\bint\\b"),
        QStringLiteral("\\blong\\b"), QStringLiteral("\\bnamespace\\b"), QStringLiteral("\\boperator\\b"),
        QStringLiteral("\\bprivate\\b"), QStringLiteral("\\bprotected\\b"), QStringLiteral("\\bpublic\\b"),
        QStringLiteral("\\bshort\\b"), QStringLiteral("\\bsignals\\b"), QStringLiteral("\\bsigned\\b"),
        QStringLiteral("\\bslots\\b"), QStringLiteral("\\bstatic\\b"), QStringLiteral("\\bstruct\\b"),
        QStringLiteral("\\btemplate\\b"), QStringLiteral("\\btypedef\\b"), QStringLiteral("\\btypename\\b"),
        QStringLiteral("\\bunion\\b"), QStringLiteral("\\bunsigned\\b"), QStringLiteral("\\bvirtual\\b"),
        QStringLiteral("\\bvoid\\b"), QStringLiteral("\\bvolatile\\b"), QStringLiteral("\\bbool\\b"),
    };


    keywordFormat.setForeground(Qt::darkBlue);
    keywordFormat.setFontWeight(QFont::Bold);
    for (const QString &pattern : GdmlPatterns) {
        rule.pattern = QRegularExpression(pattern);
        rule.format = keywordFormat;
        highlightingRules.append(rule);
    }

    /*
    classFormat.setFontWeight(QFont::Bold);
    classFormat.setForeground(Qt::darkMagenta);
    rule.pattern = QRegularExpression(QStringLiteral("\\bQ[A-Za-z]+\\b"));
    rule.format = classFormat;
    highlightingRules.append(rule);

    quotationFormat.setForeground(Qt::darkGreen);
    rule.pattern = QRegularExpression(QStringLiteral("\".*\""));
    rule.format = quotationFormat;
    highlightingRules.append(rule);

    functionFormat.setFontItalic(true);
    functionFormat.setForeground(Qt::blue);
    rule.pattern = QRegularExpression(QStringLiteral("\\b[A-Za-z0-9_]+(?=\\()"));
    rule.format = functionFormat;
    highlightingRules.append(rule);
    */

    multiLineCommentFormat.setForeground(Qt::red);

    commentStartExpression = QRegularExpression(QStringLiteral("/\\*"));
    commentEndExpression = QRegularExpression(QStringLiteral("\\*/"));
}


void Highlighter::highlightBlock(const QString &text)
{
    for (const HighlightingRule &rule : qAsConst(highlightingRules)) {
        QRegularExpressionMatchIterator matchIterator = rule.pattern.globalMatch(text);
        while (matchIterator.hasNext()) {
            QRegularExpressionMatch match = matchIterator.next();
            setFormat(match.capturedStart(), match.capturedLength(), rule.format);
        }
    }

    setCurrentBlockState(0);

    int startIndex = 0;
    if (previousBlockState() != 1)
        startIndex = text.indexOf(commentStartExpression);

    while (startIndex >= 0) {
        QRegularExpressionMatch match = commentEndExpression.match(text, startIndex);
        int endIndex = match.capturedStart();
        int commentLength = 0;
        if (endIndex == -1) {
            setCurrentBlockState(1);
            commentLength = text.length() - startIndex;
        } else {
            commentLength = endIndex - startIndex
                    + match.capturedLength();
        }
        setFormat(startIndex, commentLength, multiLineCommentFormat);
        startIndex = text.indexOf(commentStartExpression, startIndex + commentLength);
    }

}










