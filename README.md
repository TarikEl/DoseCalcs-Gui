DoseCalcs tool is built for Ubuntu and CentOS. It is one tool which comprises a package named DoseCalcs-Gui (DoseCalcs-Gui) based on Qt5-C++ and a package named DoseCalcs-Core (DoseCalcs-Core) based on Geant4 toolkit. This combination enables users to leverage the strength of the DoseCalcs-Core package for internal dosimetry using command file interface while enjoying the convenience and ease of use provided by the graphical user interface of DoseCalcs-Gui package.

Core package is based on Geant4 toolkit, which is a set of physics libraries coded in C++ that serve to solve physics problems concerning the transport of radiation in matter using the Monte Carlo method. The DoseCalcs-Core implements these functionalities to simulate the internal dosimetry scenarious. It allows the study of interaction of radiation emitted from a region source with the matter filling all other body regions in order to score the related internal dosimetry quantities. The DoseCalcs-Core package can be installed by following the prerequisites installation (DoseCalcs-Core Installations), then you can execute it separately using commands file by following the DoseCalcs-Core Simulations and DoseCalcs-Core Commands documentation.

The DoseCalcs-Gui package is a graphical user interface for DoseCalcs-Core package. It is designed to simplify the installation, interaction, and execution of the DoseCalcs-Core graphically through a visual interface. The DoseCalcs-Gui package offers intuitive controls, graphical representations, and interactive elements that enhance the user experience. Users can easily navigate through various DoseCalcs-Gui tabs to: install DoseCalcs-Core prerequisites, build DoseCalcs-Core package; save configuration; add simulation inputs; monitor the progress of simulation; visualize results, and access additional functionalities through the graphical interface. It simplifies the usage of DoseCalcs-Core for users who are more comfortable with graphical interfaces or prefer a visual representation of their work.