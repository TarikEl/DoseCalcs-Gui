DoseCalcs-Core Simulations
============================

.. GeneralConcept:

Simulation Units
----------------

Three computation modes are available for the DoseCalcs execution scenario: sequential (SEQ), multi-threading (MT), and Message Passing Interface (MPI). Only for the (SEQ) and (MT) scenarios are three run modes available: data generation mode (Gen), calculation mode (B), and graphical visualization mode (G). The simulation requires five data types after the messenger classes have initialized the commands parameters. These data types are utilized in the five simulation steps of material and geometry construction, physics setup, source defining, simulation run and scoring, and results analysis.  Each simulation step is carried out in accordance with the given commands. 

The [macros] file, which contains all input commands to construct geometry, physics, source, and define parameters of run, score, and analysis, is the first of two files that may be considered as principal code inputs of DoseCalcs [simulate] executable. Second, if the user imports geometry from a file format, the geometry file (GDML, TEXT or STL). [Simulate] as a result, produces two types of files: absorbed energy data files, which contain the calculation result produced from each thread or rank (i.e. [AE@for@Rank@0@Thread@0@alpha@Liver@0.2]); the [ResultsData] file, which contains the result tables for each defined source parameter and each quantity to score, can alternatively be generated using [merge]. 

The input files for the [merge] executable are [macros] and absorbed energy files, with [ResultsData] being the output file. The [analysis] executable uses [macros] and [ResultsData] as input files to generate graphs, histograms, and Latex tables containing results and simulation geometry data.

Note that the files paths specified in the [macros] file or DoseCalcs execution command should point to the DoseCalcs build directory or be set to the absolute file path. 

Build a Simple Example
----------------------

To use the DoseCalcs application to simulate an dosimetry problem, two user files must be used as input data: the [macros] file and the geometry description file (or files if needed).
The [macros] file includes commands for setting the required parameters for each simulation unit, including geometry (materials and volumes and/or geometry file), physics, source, run and score, and analysis.
Each of these categories has its own set of commands.

[macros] file
+++++++++++++

The [macros] file is a Geant4 macro syntax-based ASCII file.
The code loads this file to setup the simulation by selecting parameters from commands in it, avoiding the need to recompile the source code if the simulation changes. 

1. Macros file example

.. include:: /FilesExamples/macros.mac
   :literal:

.. 2. How the data are used


Geometry Input Methods
++++++++++++++++++++++

Each geometry method uses commands and\or a specific geometry file, the geometry files that can be used by DoseCalcs:

C++ method : Construct volumes by Geant4 standard solids and materials using either simple DoseCalcs commands or C++ code. 

GDML method : Construct volumes by reading .gdml file format, which contains all geometry description data, where materials, solids, logical and physical volumes parameters are given in tags format.

TEXT method : Construct volumes by reading .geom file format, which contains all geometry description data, where materials, solids, logical and physical volumes parameters are given in TEXT format.

STL method : Construct solids by reading .stl file format. Each solid is described by an STL file, adding the materials given in [macros] file to construct geometry volume.

The simulation geometry is constructed in this example using all available geometrical methods, with the most of of phantom volumes imported from GDML format file at path Phantom.gdml, brain volume created by standard Geant4 solids using C++ code added to block of G4TCPPGeometryFormat::ConstructLogicalVolumes(), thyroid volume from TEXT format file at path Thyroid.geom, spleen volume created by standard Geant4 solids by commands, and testes volume imported from STL format file.

Voxel IDs Geometry File : Contains IDs, each ID represents a specific material that fill the correspondent voxel.

DICOM Files : Each file contains data for slice (that represents either density or activity of the pixels in slice),

Tetrahedral data Files : .node file contains data for vertex, and .ele file contains data for tertrahedrons material and points,

Execution
---------

Execution Command
+++++++++++++++++

Of course, installing Geant4 is required before developing DoseCalcs code, and it is recommended that Geant4 be created for multi-threading computing mode. After that, to execute DoseCalcs, the user must specify three options when using a command line processing interface: [Run Mode], [macros file], and [Events Number Per Thread]. 

 .. code-block:: bash

    $ ./simulate [Run Mode] [macros file] [Events Number Per Thread]

[macros file] : inputs commands (macros) file

[Events Number Per Thread] : All CPU cores are considered one thread for sequential execution. In Multi-Threaded or MPI computation modes, each thread or rank simulates this number on its own thread or rank CPU core. 

Computation Modes
+++++++++++++++++

1. Sequential Execution Command

 .. code-block:: bash

    $ ./simulate B inp.mac 100000

Total number of events in simulation is 100000

2. Multi-threaded(MT) Execution Command

 .. code-block:: bash

    $ ./simulate B inp.mac 100000

100000 Events per Thread. The total number of events in the simulation will be 100000*ThreadNumber. ThreadNumber is set by command /RunAndScoreData/setNumberOfThreads

This mode is used when Geant4 is built in multi-threading.

3. MPI Execution Command

In MPI computation mode, the mpiexec or mpirun command should be set firstly as given in the following command:

 .. code-block:: bash

    $ /home/../mpirun -np [Rank Number] ./simulate B inp.mac 100000

[Rank Number] : The number of parallel simulations on the cluster, each running on a cluster unit (i.e core). For each simulation on a rank, DoseCalcs can use different source events data. This computation mode is enableb when the DoseCalcs is built with WITH_G4MPI_USE=ON

100000 Events per Rank. The total number of events in the simulation will be 100000*[Rank Number].

.. GRunModeCommand:

Run Modes
+++++++++

One of the [Run Modes] is B(batch), G(Graphical), or Gen(Generation). The Graphical mode is used to visualize(v) the geometry built by DoseCalcs and to test the events initial positions using the command /TestEventsInitialPositions. The box dimensions may also be shown by using the command /ShowSourceBox in the macros file. The graphical mode can only be utilized in sequential and multi-threading computation modes, as previously stated. The MPI computation mode, on the other hand, is only for Batch and Generation Run Modes. The Batch run mode is for simulation and scoring, whereas the Gen run mode is for generating event data including initial positions, energies, and momentum directions. 

[Run Mode] : can be : B (Batch "default if we don't set [Run Mode]), G (Graphical), T (Terminal) or Gen (Generation).

1. G Run Mode

In Sequential mode (DoseCalcs built with -DWITH_G4MPI_USE set to OFF), G mode is only utilized for geometry visualization, therefore we may use the commands /SourceData/TestEventsInitialPositions and /SourceData/ShowSourceBox. 

 .. code-block:: bash

    $ ./simulate G [macros file] [.]

[.] : can be : v or an empty value means view geometry using the Qt interface using the OGL driver; d means download the geometry image in PS format using the DAWNFILE driver. Instead of using the v or d values, the user can use his own Geant4 visualization instructions by setting the visualization macros file path. The user should provide xy, xz, or yz for d value; the default plane when d value is specified is xy plane. 

.. For voxelized geometries, the default value is token from /GeometryData/visualizeVoxelizedPlanes command plane parameter.  

The DAWNFILE executable [dawn] that Geant4 uses to generate a high-quality picture can be obtained at https://twiki.cern.ch/twiki/bin/view/CLIC/DawnVisualization; Follow the installation instructions on this page; before executing graphical run mode (:ref:'GRunModeCommand'), in the terminal, type: 

 .. code-block:: bash

    $ export G4DAWNFILE_VIEWER="dawn -d"


2. Gen Run Mode

Gen mode is only used to generate event data, which is then used in the simulation with the multi-threading and MPI modes. A set of /SourceData/ commands must be used to define the data to be generated in the [macros] file. 

In  Multi-threaded(MT) or Sequential execution modes
 .. code-block:: bash

    $ ./simulate Gen [macros file]

The data units are generated sequentially in multi-threading computation mode, starting with initial positions in source volume 1, then source volume 2,..., energy 1, energy 2,..., momentum direction. 

In  MPI execution modes
 .. code-block:: bash

    $ /home/../mpirun -np [Rank Number] ./simulate Gen [macros file]

[Rank Number] : To generate initial positions, energies, and momentum direction, the [Rank Number] must be equal to the data that will be generated, where each rank generates a data unit. For example, to generate initial positions in three source volumes, two energies, and one momentum direction, the [Rank Number] must be set to 6. 

3. B Run Mode

B mode is used to simulate the events interactions and gives the dosimetric quantities outputs scored.

Assuming A particles (Part), B source regions (Src), C energies (Ene), and D momentum directions (MomDir) as source inputs, the data generation and calculation tasks are distributed among the ranks. As a result, the number of ranks required for the generation of this event data is N=B+C+D, where each rank generates a unique data file, but the number of ranks required for calculation is (N = A * B * C * D). Another DoseCalcs option is to simulate one Part-Src-Ene-MomDir combination (A=B=C=D=1 -> N=1) on all ranks (or threads), which is one simulation utilizing the same source data on all ranks. In this scenario, the total number of simulated occurrences is for a single simulation is ([TotalEventsNumber]=N*[EventsNumber]).

 .. image:: /images/SimModeGen.png

 .. image:: /images/SimModeGen.png

The total number of simulated events in Multi-threading or MPI modes is the [EventsNumber] multiplied by the thread number or rank number, respectively. In multithreading mode, this total number must not exceed the maximum G4int value, and in MPI mode, EventsNumber must not exceed the maximum G4int value for each rank. 

Results
-------

Each thread calculation produces its own results. Note that the global results are calculated directly by [simulate] in multithreading mode (or by [merge] executable in MPI mode) and produces [ResultsData] file. The [analysis] executable produces graphs and histograms based on the data in [ResultsData], [macros], and events data files. 

Absorbed energy data files 
++++++++++++++++++++++++++

These files represent the calculation results; during simulation, each thread or rank generates its own file, which is identified by its name, which is a combination of rank ID, thread ID, particle name, source region name, and particle energy. Below is a sample of this file, which was generated by thread 2 and simulates 200000000 emitted gamma from the brain, with mono-distributed energies of 1 MeV. 

.. include:: /FilesExamples/AE@for@Rank@0@Thread@2@ORNLAdultMale@gamma@Brain@1
   :literal:
   
DoseCalcs calculates absorbed energy in this region, absorbed energy square, and number of histories recorded in each region in simulation geometry. To compute Monte Carlo statistical uncertainty, the two last values are required. 

[ResultsData]
+++++++++++++

The internal dosimetry quantities (AE, AD, AF, SAF, S, H, E, and DR) for each source region, particle name, particle energy, and target region are all included in the text file. It's organized as a series of text lines that summarize the obtained data for each target, such as quantity value, standard deviation, relative standard deviation in percent, among other things.

The generated scores are appended to the results file for each simulation (i.e. source volume, particle, and energy combination. The produced data is written in a simple manner, with the first line including simulation data such as scored quantity, source volume name, particle name, and so on, followed by data lines containing the results for each scored volume. 

.. include:: /FilesExamples/ResultsData
   :literal:

The header file reveals that for each source data, the results are produced for two scored quantities: S and SAF, followed by scored target regions data that are presented in separated lines, and each simulation is indicated by individual source data (particle name, source region name, energy). 

1. Header line

Here, the list of the most important data that are used by [analysis] executable, the first result: 

SAF 		: scored quantity
Thyroid 	: source region name
e+ 		: particle name
0.6335 		: particle initial energy in (MeV)
Rayleigh 	: particle energy distribution
Isotropic 	: particle momentum direction distribution
100000000Event 	: total number of simulated events
1455759402Step 	: total number of simulated steps
0.000104726ms 	: mean dutarion for events simulation 
628.357 min	: mean dutarion for simulation

2. Scored Data lines

The results for each scored volume are written in the format: 
 
Scored_Region Scored_Quantity Standard_Deviation Relative_Standard_Deviation Number_Of_Steps Region_Mass Region_Volume Region_Density
    
How to get this results
+++++++++++++++++++++++

1. Which Computation mode!

Sequential, multi-threading, and MPI computation modes are included in the DoseCalcs code. Only for the (SEQ) and (MT) scenarios are three run modes available for each computation mode: data generation mode (Gen), calculation mode (B), and graphical visualization mode (G). 

Sequential : simulation of total number of events on one thread (master thread)

Multi-threading : simulation of total number of events is divided on the number of threads, and each thread simulates corresponding events number with the same source data.

MPI : each rank simulates one run sequentially and the total number of runs will be the number of ranks, each rank simulates a specific source data separately.

The results are generated using calculation mode (B) even in Sequential, Multi-threading or MPI computation mode.

2. [merge] executable

The results calculation steps are based on two types of files: [macros] and thread or rank region result data files (absorbed energy files), in order to obtain all necessary data, merge threads and ranks, and calculate and write the final simulation results to the file [ResultsData]. 

 .. code-block:: bash
  
   $ ./merge [macros file] v 

The [merge] command takes two arguments, the first of which is the [macros] file which contains all needed parameters of scoring, the units to be used, and results directory path where the absorbed energy files (i.e. AE@For@Rank@i@Thread@j@SourceOrgan@Particle@Energy) exist. The second argument is "v" which allows you to output more information when reading and calculating results. 

If a simulation on a rank crashes due to a memory problem, the [merge] executable can be used to generate results from threads or ranks files that completed the simulation successfully and produced a results file for their respective thread or rank.


Analysis
--------

prequisites to analysis
+++++++++++++++++++++++

The objective of a graph or histogram is to present data that is too large or complicated to be fully expressed in words and in a limited amount of space. [analysis] program has been designed as a direct interface to ROOT Analysis System ROOT in order to produce graphs and histograms efficiently. This executable generate graphs based on the parameters of the /AnalysisData/ commands in the [macros] file. Note that the graph type, comparing type, reference name, reference file path, and other analysis command parameters must be specified in the [macros] file DoseCalcs. This file will be used as an input file for analysis tasks in addition to the final results file ([ResultsData]). However, the parameters in the [macros] file need to be adjusted for [analysis] to provide the required analysis results.

To use the [analysis] executable, the user must build the DoseCalcs code with option WITH\_ANALYSIS\_USE=ON and set the ROOT\_Dir=/../../root . This will build an executable called analysis.

 .. code-block:: bash
  
   -DWITH_ANALYSIS_USE=ON  -DROOT_DIR=/home/.../RootInstallDir

By executing the [analysis] executable after building DoseCalcs and configuring the Root install path, you can generate graphs using the /AnalysisData/generateSelfCrossGraphs command arguments. Note that the comparing type, graph type, reference name, reference file, and other analysis command arguments are included in this command. Only Latex format tables containing geometry regions data and scored results data can be generated if DoseCalcs is built with -DWITH_ANALYSIS_USE=OFF. 

When the source and target regions are the same, the graph is called a self-absorption graph. The data of all target region that receive radiation from a single source region is saved in a cross-irradiation graph, also known as a source graph. Given reference data in a file with a simple determined format, the code can automatically compare simulation results with this reference data and generate graphs and tables to make simulation accuracy measurement easier, as well as to search in the physical context for the source of the discrepancy between reference data and code results. 

The graphs and histograms are generated via a direct interface to the ROOT Analysis System, which supports many file formats for the graphs and histograms generated. In DoseCalcs, the user can select the suitable format. Self-absorption and cross-irradiation data, as well as additional graphs such as relative errors and relative standard deviation, macroscopic cross-section, and so on, can be generated. The simulation inputs, such as initial positions, energies, and momentum directions, can also be represented as 2D and 3D histograms. 

[analysis] executable
+++++++++++++++++++++

The [simulate] command takes two arguments, the first of which is the [macros] file which contains all needed parameters of analysis tasks, including the results directory path where the absorbed energy files (i.e. AE@For@Rank@i@Thread@j@SourceOrgan@Particle@Energy) exist. The second argument is "v" which allows you to output more information when reading and calculating results. 

The the simple command line to analysis is:

 .. code-block:: bash
  
   $ ./analysis [macros file] v

After the execution, the resulted files (graphs, histograms and latex tables) will be generated in the result directory path under the /Graphs_Histograms directory. The DoseCalcs and reference data can be generated in graphs and tables, which necessity providing a reference file in a specific format. 


For [merge] and [analysis] executables, if the given result directory path not found, the [Results] directory in the DoseCalcs build directory will be used.


.. Analysis output directory and files
.. +++++++++++++++++++++++++++++++++++

.. Build directory Final Structure
.. -------------------------------



