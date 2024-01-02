#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(DoseCalcs)

#to remove the warnings messages
add_definitions(-w)

#set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

if(WITH_VERBOSE_USE)
    add_definitions(-DVERBOSE_USE)
endif()

#find_package( Qt5Core )
#find_package( Qt5PrintSupport )
#find_package( Qt5Gui )
#find_package( Qt5Widgets )
#find_package( Qt5Network )
#option(WITH_GUI_USE "Build DoseCalcs with GUI" OFF)
#if(WITH_GUI_USE)
#    find_package(Qt5 COMPONENTS Core Widgets PrintSupport REQUIRED)
#    add_definitions(-DGUI_USE)
#    if(Qt5_FOUND)
#        message("\n\n================ QT_FOUND =============================================\n")
#        message("Libraries : ${Qt5_LIBRARIES}\n")
#        message("Libraries : ${Qt5Widgets_LIBRARIES}\n")
#        message("Libraries : ${Qt5PrintSupport_LIBRARIES}\n")
#        message("Include Files : ${Qt5_INCLUDE_DIRS}\n")
#        message("Include Files : ${Qt5Widgets_INCLUDE_DIRS}\n")
#        message("Include Files : ${Qt5PrintSupport_INCLUDE_DIRS}\n")
#        include_directories( ${Qt5_INCLUDE_DIRS} )
#        file(GLOB guisources ${PROJECT_SOURCE_DIR}/gui/*.cpp)
#        file(GLOB guiheaders ${PROJECT_SOURCE_DIR}/gui/*.h)
#        set(CMAKE_AUTOUIC ON)
#        add_executable(gui ${PROJECT_SOURCE_DIR}/gui/main.cpp ${guisources} ${guiheaders})
#        target_link_libraries(gui Qt5::Widgets Qt5::PrintSupport Qt5::Network ${Qt5_LIBRARIES})
#endif()
#endif()


#----------------------------------Libraries to link---------------------------------------------



#include(${QT_USE_FILE})

#-----------------------------------Geant4 with visualization--------------------------------------
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
    if(Geant4_FOUND)
        message("\n\n================ Geant4_FOUND =============================================\n")
        message("Libraries : ${Geant4_LIBRARIES}\n")
        message("Include Files : ${Geant4_INCLUDE_DIRS}\n")
    endif()
else()
    find_package(Geant4 REQUIRED)
    if(Geant4_FOUND)
        message("\n\n================ Geant4_FOUND =============================================\n")
        message("Libraries : ${Geant4_LIBRARIES}\n")
        message("Include Files : ${Geant4_INCLUDE_DIRS}\n")
    endif()
endif()
include(${Geant4_USE_FILE})

#---------------------------------- DCMTK ---------------------------------------
option(WITH_DCMTK_USE "Build with using DCMTK for DICOM CT and PET Data Reading" OFF)
if(WITH_DCMTK_USE)
    find_package(DCMTK REQUIRED)
    add_definitions(-DDCMTK_USE)
    file(GLOB sourcesDicom ${PROJECT_SOURCE_DIR}/src/Dicom*)
    file(GLOB headersDicom ${PROJECT_SOURCE_DIR}/include/Dicom*)
    if(DCMTK_FOUND)
        message("\n\n================ DCMTK_FOUND =============================================\n")
        message("Libraries : ${DCMTK_LIBRARIES}\n")
        message("Include Files : ${DCMTK_INCLUDE_DIRS}\n")
    endif()

endif()

#---------------------------------- Voxelized Geometry ---------------------------------------

#option(WITH_VOX_USE "Build with using dicom files for geomery construction" OFF)
#if(WITH_VOX_USE)
#    add_definitions(-DVOX_USE)
#endif()

#---------------------------------- G4MPI ----------------------------------------------------------
option(WITH_G4MPI_USE "Build example with Geant4 MPI library" OFF)
if(WITH_G4MPI_USE)
    find_package(MPI REQUIRED) # Now cmake will automatically find a MPI implementation in your system. If you have multiple different MPI versions, and want to specify which one to compile with, you can set the MPI_C_COMPILER and MPI_CXX_COMPILER variables to the corresponding mpicc and mpicxx compiler wrappers. The CMake module will then use those to figure out all the required compiler and linker flags itself. -DMPI_C_COMPILER=/.../  -DMPI_CXX_COMPILER=/.../  -DMPI_CXX_INCLUDE_PATH=/.../
    #find_package(G4mpi REQUIRED)
    add_definitions(-DG4MPI_USE)
    #file(GLOB sourcesG4MPI ${PROJECT_SOURCE_DIR}/src/G4MPI*
        #${PROJECT_SOURCE_DIR}/src/G4VMPI*
        #${PROJECT_SOURCE_DIR}/src/G4UImpish.cc
        #${PROJECT_SOURCE_DIR}/src/G4VUserMPIrunMerger.cc)
    #file(GLOB headersG4MPI ${PROJECT_SOURCE_DIR}/include/G4MPI*
        #${PROJECT_SOURCE_DIR}/include/G4VMPI*
        #${PROJECT_SOURCE_DIR}/include/G4UImpish.hh
        #${PROJECT_SOURCE_DIR}/include/G4VUserMPIrunMerger.hh)

    if(MPI_FOUND)
        message("\n\n================ MPI_FOUND =============================================\n")
        message("Libraries : ${MPI_LIBRARIES}\n")
        message("Include Files : ${MPI_INCLUDE_DIRS}\n")
    endif()
endif()

#----------------------------------Multi-Threading-------------------------------------------------
#option(WITH_MT_USE "Build with using multi-threading" OFF)
#if(WITH_MT_USE)
#    add_definitions(-DMT_USE)
#endif()

#----------------------------------GDML-----------------------------------------------------------
option(WITH_GDML_USE "Build example with GDML" OFF)
if(WITH_GDML_USE)
    add_definitions(-DGDML_USE)
else()
endif()

#----------------------------------ROOT-----------------------------------------------------------
option(WITH_ANALYSIS_USE "Build example with analysis objects" OFF)
if(WITH_ANALYSIS_USE)
find_package(ROOT REQUIRED)
add_definitions(-DANALYSIS_USE)
if(ROOT_FOUND)
    message("\n\n================ ROOT_FOUND =============================================\n")
    message("Libraries : ${ROOT_LIBRARIES}\n")
    message("Include Files : ${ROOT_INCLUDE_DIRS}\n")
endif()
endif()

#----------------------------------Include dirs---------------------------------------------------

include_directories(${PROJECT_SOURCE_DIR}/include  ${Geant4_INCLUDE_DIRS} ${GDML_INCLUDE_DIRS} )

if(WITH_GDML_USE)
    include_directories( ${GDML_INCLUDE_DIRS} )
endif()

if(WITH_ANALYSIS_USE)
    include_directories( ${ROOT_INCLUDE_DIRS} )
endif()

if(WITH_G4MPI_USE)
    include_directories( ${MPI_INCLUDE_DIRS} )
endif()

if(WITH_DCMTK_USE)
    include_directories( ${DCMTK_INCLUDE_DIRS} )
endif()

#----------------------------------Source and Header Files for simulate---------------------------------------------

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/G4T*)
file(GLOB sources1 ${PROJECT_SOURCE_DIR}/src/TET*)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/G4T*)
file(GLOB headers1 ${PROJECT_SOURCE_DIR}/include/TET*)

add_executable(simulate simulate.cc ${sources} ${sources1} ${sourcesDicom} ${sourcesG4MPI} ${headers} ${headers1} ${headersDicom} ${headersG4MPI})

#----------------------------------Libraries to link---------------------------------------------

target_link_libraries(simulate ${Geant4_LIBRARIES})

if(WITH_G4MPI_USE)
    #target_link_libraries(simulate ${G4mpi_LIBRARIES})
    target_link_libraries(simulate ${MPI_LIBRARIES})
endif()

if(WITH_DCMTK_USE)
    target_link_libraries(simulate ${DCMTK_LIBRARIES})
endif()

if(WITH_GDML_USE)
    target_link_libraries(simulate ${GDML_LIBRARIES})
endif()

#-----------------------------Build simulate executable----------------------------------------------

add_custom_target(IntDoseCalcs DEPENDS simulate)
install(TARGETS simulate DESTINATION bin)

#-----------------------------Build merge executable---------------------------------------------

#if(WITH_G4MPI_USE)

file(GLOB sourceFile ${PROJECT_SOURCE_DIR}/src/G4TResultCalculation.cc)
file(GLOB headerFile ${PROJECT_SOURCE_DIR}/include/G4TResultCalculation.hh)
add_executable(merge merge.cc ${sourceFile} ${headerFile} )
target_link_libraries(merge ${Geant4_LIBRARIES})

install(TARGETS merge DESTINATION bin)

#endif()

#-----------------------------Build analysis executable---------------------------------------------

option(WITH_ANALYSIS_USE "Build example with analysis objects" OFF)
if(WITH_ANALYSIS_USE)

    #find_package(ROOT REQUIRED)
    #add_definitions(-DANALYSIS_USE)

    include_directories( ${ROOT_INCLUDE_DIRS} )
    #include_directories(/home/tarik/Desktop/WorkSpace/1BuildsAndInstallations/root_v6.24.02.Linux-ubuntu20-x86_64-gcc9.3/root/include)

    file(GLOB sourceRoot ${PROJECT_SOURCE_DIR}/src/G4DoseCalcsAnalysis.cc)
    file(GLOB headerRoot ${PROJECT_SOURCE_DIR}/include/G4DoseCalcsAnalysis.hh)

    add_executable(analysis analysis.cc ${sourceRoot} ${headerRoot})

    #target_include_directories(analysis PUBLIC ${ROOT_INCLUDE_DIRS})
    target_link_libraries(analysis ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

    install(TARGETS analysis DESTINATION bin)
endif()

#-------------------------------Copy Files and scripts ------------------------------------------


if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Results)
else()
	file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Results)
endif()
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Results DESTINATION ${CMAKE_CURRENT_BINARY_DIR})



if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Scripts)
else()
	file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Scripts)
endif()
#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Scripts DESTINATION ${CMAKE_CURRENT_BINARY_DIR})



if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/EventsData)
else()
	file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/EventsData)
endif()
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/EventsData DESTINATION ${CMAKE_CURRENT_BINARY_DIR})



set(DoseCalcs_SCRIPTS ) # openGLVis.mac buildInLenovo.sh  macrosConst.mac macrosText.mac macrosGdml.mac macrosDicom.mac macrosVoxel.mac ShellFile.sh MacroCmd.mac
foreach(_script ${DoseCalcs_SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${PROJECT_BINARY_DIR}/${_script}
        COPYONLY
        )
endforeach()

