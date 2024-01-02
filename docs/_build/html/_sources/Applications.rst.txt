.. _Applications:

DoseCalcs-Core Applications
=============================

MIRD phantom in GDML format file
--------------------------------

The used input files for [simulate] executable was the macros.mac and MIRDGDMLPhantom.gdml, The files are provided to users in links below:

:download:`Phantom geometry in GDML format file </FilesDownloads/MIRD/MIRDGDMLPhantom.gdml>`.   

:download:`macros file </FilesDownloads/MIRD/macros.mac>`.   

Example inputs description:
+++++++++++++++++++++++++++

The simulation geometry will be imported from the MIRDGDMLPhantom.gdml file, in which all volumes, materials, and shapes, including world volume data, are described by gdml tags to describe MIRD phantom regions, and the physics used to track the particles is EMS3 (i.e. ElectromagneticStandard 3), with a cut in range of 1 mm and an energy threshold of 1 keV.

Initial positions in the liver, adrenals, and kidneys were used to provide data for the radiation source.
The event momentum directions were generated using an isotropic distribution with starting Mono energy of 0.01 0.015 0.02 0.03 0.05 0.1 0.2 0.5 and 1 MeV. 

In calculation mode, the scored quantity was Specific Absorbed Fraction (SAF), which is scored in three : liver, kidney, and adrenal.
The "m" value indicates that each rank will use the MPI computing method to simulate a combination of Particle-SourceRegion-Energy-MomentumDirection.
There are a total of 27 ranks that should be used.

The generation (Gen) mode was the initial run mode.
By changing the event data number and generation flags for initial positions, energies, and momentum directions in the macro file, you may speed up the calculation mode while getting event data: 

.. code-block::
    :linenos:
    :emphasize-lines: 1

    #######################################################  Source  ############################
    ....
    ..
    /SourceData/useDataGenerationFiles read yes yes yes

The calculation was carried out in MPI mode, with 13 ranks utilized to produce initial positions in three source regions, nine energies, and one momentum direction distribution.
Using the following command: 

 .. code-block:: bash

    $ /home/../mpirun -np 13 ./simulate Gen macros.mac 

The command /SourceData/useDataGenerationFiles is used to enable data generation in files, and the 40000000 initial positions for events are generated.

The second run mode was calculation mode (B), in which the computation was done in MPI mode to speed up the calculation mode while getting event data, and the used ranks were 27 ranks to conduct 27 simulations. The command used on the terminal, from the build directory, was: simulations. In terminal, from build directory, the used command was:

 .. code-block:: bash

    $ /home/../mpirun -np 27 ./simulate B macros.mac 40000000

If the user wants to produce events data without having to utilize events data files, they should comment the command in the macro file as shown below. 

.. code-block::
    :linenos:
    :emphasize-lines: 1

    #######################################################  Source  ############################
    ....
    ..
    #/SourceData/useDataGenerationFiles read yes yes yes
    
Results
+++++++

The Absorbed energy files, and [ResultsData] files produced by the [simulate] executable are available to users in the URLs below: 

:download:`The calculated results data for scored quantities and regions </FilesDownloads/MIRD/Results/ResultsData>`.  

The user can generate new results for new scored quantities or new target regions by using the [merge] program and specifying the new quantities and targets in [macros] as follows: 

 .. code-block::
    :linenos:
    
    QuantitiesToScore                AD AF S       
    RegionsNamesToScore              Liver Kidney Adrenal Thyroid Testes Spleen  

Then executing [merge]:

  .. code-block:: bash

    $ ./merge [Macros File] v

The new results will be written to [ResultsData] file.

Generate graphs, histograms and latex format Tables
+++++++++++++++++++++++++++++++++++++++++++++++++++

The various graph types were generated using the command below, which used the [macros] file as an input file. 

 .. code-block:: bash

    $ ./analysis [Macros File] v

An example of the generated graphs is shown below 

.. list-table:: 

    * - .. figure:: /images/Cross_Result_SAF_Liver_gamma.png

           Results Cross Irradiation Graph

      - .. figure:: /images/Cross_ReferenceResult_SAF_gamma_Adrenal_Kidney_DoseCalcs_vs_MIRD.png

           Results and Reference Cross Irradiation Graph
           
    * - .. figure:: /images/Self_Result_SAF_gamma.png

           Results Self Absorption Graph 

      - .. figure:: /images/Self_ReferenceResult_SAF_gamma_Liver_DoseCalcs_vs_MIRD.png

           Results and Reference Self Irradiation Graph
           
    * - .. figure:: /images/RelativeSDv_SelfSAF.png

           Self Relative SDv Graph

      - .. figure:: /images/RelativeSDv_CrossSAF_Liver.png

           Cross Relative SDv Graph
           
    * - .. figure:: /images/RelativeError_Self_SAF_DoseCalcs_vs_MIRD.png

           Self Relative error Graph

      - .. figure:: /images/RelativeError_Cross_SAF_Liver_DoseCalcs_vs_MIRD.png

           Cross Relative error Graph
           
    * - .. figure:: /images/Mass_SAFForAllEnergies_inSelfAbsorption.png

           Mass SAF in Self Absorption Graph 

      - .. figure:: /images/Macroscopic_Cross_Section_for_gamma_in_material_SoftTissue.png

           Macroscopic Cross Section Graph

The [Scripts] and [Results] can be downloaded and pasted to the build workspace to simulate with the same inputs or with user adjustment to get new desired results. 

:download:`Compressed file, contains [Scripts] and [Results] directories </FilesDownloads/MIRD/f.tar.xz>`.

Geometry from TEXT format file 
------------------------------

The simulation inputs are identical to those for GDML phantom geometry; the geometry was taken from a TEXT format file called MIRDTEXTPhantom.geom. As shown in the following command: 

.. code-block::
    :linenos:
    :emphasize-lines: 1

    #######################################################  Geometry ##########################
    #/GeometryData/createVolume Scripts/MIRDGDMLPhantom.geom
    /GeometryData/createVolume Scripts/MIRDTEXTPhantom.geom
    
The file is provided to users in link below:

:download:`Phantom geometry in TEXT format file </FilesDownloads/MIRD/MIRDTEXTPhantom.geom>`.  


Basic geometry from Geant4 Standard solids
------------------------------------------

The macros were used as inputs for the [simulate] executable. No geometry files are needed; instead, DoseCalcs' basic commands are used to build the world and geometry volumes. Users can access the file using the following link:

:download:`macros file </FilesDownloads/StandardGeant4Solids/macros.mac>`.

The [Scripts] and [Results] can be downloaded and pasted to the build workspace to simulate with the same inputs or with user adjustment to get new desired results. 

:download:`Compressed file, contains [Scripts] and [Results] directories </FilesDownloads/StandardGeant4Solids/f.tar.xz>`.

After copying [Scripts] and [Results] to the build workspace directory, the user should execute the command below in terminal and from the build directory to visualize the world geometry: 

 .. code-block:: bash

    $ ./simulate G macros.mac

To configure the desired simulation inputs, the user should also adjust the macros.mac to simulate a volume as a source and calculate various dosimetry quantities. Run the command below: 

 .. code-block:: bash

    $ ./simulate B macros.mac 1000000


Combination of different geometry methods to construct simulation phantom
-------------------------------------------------------------------------

This example's [Scripts] and [Results] may be downloaded, and the user should copy them to the build workspace in order to simulate with the same inputs or with user adjustment to obtain new desired results. 

:download:`Compressed file, contains [Scripts] and [Results] directories </FilesDownloads/DifferentMethods/f.tar.xz>`.

Testes are imported from STL file, thyroid from TEXT file, spleen from Geant4 standard solids using the DoseCalcs command, brain from adding geometrical C++ code to G4TCPPGeometryFormat::ConstructLogicalVolumes() function block (which requires re-building DoseCalcs source code), and all remaining organs, including world volume, from GDML files. The spleen volume-related DoseCalcs commands and brain-related C++ code are listed below, along with links to the geometry files: 

Spleen region related geometry commands defined in macros.mac:

 .. code-block:: bash

    /GeometryData/createSolid Ellipsoid SpleenSol 3.2 2.3 5.7 cm
    /GeometryData/createVolume Spleen SpleenSol SoftTissue Trunk 10.79 2.94 1.8 0 0 0 cm degree


:download:`Testes STL format file </FilesDownloads/DifferentMethods/Testes.ast>`.

:download:`Thyroid TEXT format file </FilesDownloads/DifferentMethods/Thyroid.geom>`.

:download:`World and remains organs GDML format file </FilesDownloads/DifferentMethods/ORNLPhantom.gdml>`.

:download:`Brain region C++ related code implemented in G4TCPPGeometryFormat C++ format file </FilesDownloads/DifferentMethods/G4TCPPGeometryFormat.cc>`.

This C++ file format cannot be utilized directly by DoseCalcs like other geometry files; instead, the C++ class file should be used to replace the existing G4TCPPGeometryFormat class file in the DoseCalcs source directory, or copy G4TCPPGeometryFormat::ConstructLogicalVolumes() block code from this file to a source file that already exists. To account for the new C++ built brain volume, the user must build the updated DoseCalcs source.


After copying [Results], geometry files and [macros] to the build workspace directory, add the C++ code for brain volume to G4TCPPGeometryFormat::ConstructLogicalVolumes() block and build DoseCalcs source to visualize the world geometry. The user should execute the command below in terminal and from the build directory: 

 .. code-block:: bash

    $ ./simulate G macros.mac

And to simulate a volume as a source and calculating some dosimetry quantities, the user should edit the macros.mac to set the desired simulation inputs and run the command below:

 .. code-block:: bash

    $ ./simulate B macros.mac 1000000
    $ ./merge macros.mac


.. ICRP Voxelized Phantom 
.. ----------------------

.. DICOM Phantom
.. -------------

.. TET Phantom 
.. -----------


