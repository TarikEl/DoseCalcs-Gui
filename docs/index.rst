
Welcome to DoseCalcs's documentation!
=====================================

DoseCalcs tool (`DoseCalcs application <https://github.com/TarikEl/DoseCalcs-Gui>`_) is built for Ubuntu and CentOS. It is one tool which comprises a package named gui (DoseCalcs-Gui) based on Qt5-C++ and a package named core (DoseCalcs-Core) based on Geant4 toolkit. This combination enables users to leverage the strength of the DoseCalcs-Core for internal dosimetry using command file interface while enjoying the convenience and ease of use provided by the graphical user interface of DoseCalcs-Gui.

DoseCalcs-Core is based on Geant4 toolkit, which is a set of physics libraries coded in C++ that serve to solve physics problems concerning the transport of radiation in matter using the Monte Carlo method. The DoseCalcs-Core implements these functionalities to simulate the internal dosimetry scenarious. It allows the study of interaction of radiation emitted from a region source with the matter filling all other body regions in order to score the related internal dosimetry quantities. The DoseCalcs-Core can be installed by following the prerequisites installation (`DoseCalcs-Core Installations <https://dosecalcs-gui.readthedocs.io/en/latest/Installation.html>`_), then you can execute it separately using text commands by following the `DoseCalcs-Core Simulations <https://dosecalcs-gui.readthedocs.io/en/latest/GeneralConcept.html>`_ and `DoseCalcs-Core Commands <https://dosecalcs-gui.readthedocs.io/en/latest/Commands.html>`_ documentation.

The DoseCalcs-Gui is a graphical user interface for DoseCalcs-Core. It is designed to simplify the installation, interaction, and execution of the DoseCalcs-Core graphically through a visual interface. The DoseCalcs-Gui offers intuitive controls, graphical representations, and interactive elements that enhance the user experience. Users can easily navigate through various DoseCalcs-Gui tabs to: install DoseCalcs-Core prerequisites, build DoseCalcs-Core; save configuration; add simulation inputs; monitor the progress of simulation; visualize results, and access additional functionalities through the graphical interface. It simplifies the usage of DoseCalcs-Core for users who are more comfortable with graphical interfaces or prefer a visual representation of their work.

DoseCalcs tool can be used for advanced research and as an educational tool in internal radiation protection field. `Download package. <https://codeload.github.com/TarikEl/DoseCalcs-Gui/zip/refs/heads/main>`_. The beginning of documentation will be with DoseCalcs-Gui, and then the DoseCalcs-Core can be manipulated graphically.

Use DoseCalcs-Core with GUI (DoseCalcs-Gui)
----------------------------------------------

When user clones or downloads DoseCalcs zip file and unzip it, he will find DoseCalcs-Core directory along with the DoseCalcs-Gui directory. To install the DoseCalcs-GUI package, which is a Qt C++ package, the user should install the prerequisites of Qt5 and xterm as mentioned in `DoseCalcs-Gui Installations <https://dosecalcs-gui.readthedocs.io/en/latest/GUI.html#installing-the-dosecalcs-gui-component>`_. Create an empty build directory and enter into, then build the DoseCalcs-GUI using qmake using DoseCalcs.pro file "qmake /home/user/../DoseCalcs.pro", which results in the generation of the [DoseCalcs] executable that user will use to open DoseCalcs-GUI application by calling it "./DoseCalcs". This is if the user will use directly the GUI version of DoseCalcs.


Use DoseCalcs-Core without GUI
-------------------------------

Also, the user can use the DoseCalcs-Core directly, by installing the prerequisites as mentioned in `DoseCalcs-Core Installations <https://dosecalcs-gui.readthedocs.io/en/latest/Installation.html>`_; creating an empty directory and entering into this directory; then building the DoseCalcs-Core using cmake command, which results in the generation of the [simulate], [merge], and [analysis] executables that user will use to launch the Monte Carlo simulation and merge and analyze the results.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   GUI
   Installation
   GeneralConcept
   Commands
   Applications

   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
