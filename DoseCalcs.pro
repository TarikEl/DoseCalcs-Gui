#-------------------------------------------------
#
# Project created by Tarik ElGhalbzouri QtCreator 2019-02-21T18:32:03
#
#-------------------------------------------------

QT += widgets network printsupport # core # gui # concurrent

#QT       += core gui network concurrent
#greaterThan(QT_MAJOf R_VERSION, 4): QT += widgets printsupport
#greaterThan(QT_MAJOR_VERSION, 5): QT += xml

QMAKE_POST_LINK += $$quote({ export QT_QPA_PLATFORM=xcb; })

TARGET = DoseCalcs
TEMPLATE = app

#LIBS += -L/usr/local/lib/Geant4-10.5.0/ -L/usr/lib/ -L/usr/local/lib/
#INCLUDEPATH += /usr/local/include/Geant4/ #/usr/include /usr/local/include

QMAKE_LFLAGS += -Wl,-rpath,"'\$$ORIGIN'"

# to disable warning messages
CONFIG += warn_off
#for xterm embeded in X11 window
#CONFIG += xcb wayland

#DESTDIR=gui #Target file directory
OBJECTS_DIR = gui #Intermediate object files directory
UI_HEADERS_DIR = gui
UI_DIR = gui
UI_SOURCES_DIR = gui
RCC_DIR = gui
QMAKE_MOC_SRC = gui
MOC_DIR = gui

SOURCES += gui/main.cpp\
    gui/customtextedit.cpp \
    gui/dicomreaderseg.cpp \
    gui/geometriesvisualization.cpp \
    gui/geometrymodellingeditor.cpp \
    gui/mainwindow.cpp\
    #GdmlFormatVolumesCreation.cpp \
    gui/highlighter.cpp \
    #organschooserdialog.cpp \
    gui/filesManager.cpp \
    #visualizationManager.cpp \
    gui/errorManager.cpp \
    gui/geant4_app_init.cpp \
    gui/installationDialog.cpp \
    gui/httpdownload.cpp \
    gui/plotDialog.cpp \
    gui/qcustomplot.cpp \
    gui/rootdialog.cpp \
    gui/runmpisystem.cpp \
    gui/terminal.cpp

HEADERS  += gui/mainwindow.h \
    #GdmlFormatVolumesCreation.h \
    #organschooserdialog.h \
    gui/customtextedit.h \
    gui/dicomreaderseg.h \
    gui/filesManager.h \
    #visualizationManager.h \
    gui/errorManager.h \
    gui/geant4_app_init.h \
    gui/geometriesvisualization.h \
    gui/geometrymodellingeditor.h \
    gui/highlighter.h \
    gui/installationDialog.h \
    gui/httpdownload.h \
    gui/plotDialog.h \
    gui/qcustomplot.h \
    gui/rootdialog.h \
    gui/runmpisystem.h \
    gui/terminal.h

FORMS    += gui/mainwindow.ui \
    #GdmlFormatVolumesCreation.ui \
    gui/dicomreaderseg.ui \
    gui/geometriesvisualization.ui \
    gui/geometrymodellingeditor.ui \
    gui/organschooserdialog.ui \
    gui/installationDialog.ui \
    gui/httpdownload.ui \
    gui/plotDialog.ui \
    gui/rootdialog.ui \
    gui/runmpisystem.ui

#QMAKE_CXX = g++
#QMAKE_CXX = mpicxx

RESOURCES += \
    gui/Resources.qrc


if(!exists( $$OUT_PWD/PackagesAndFiles/Config) ) {
    copypackages.commands = $(COPY_DIR) $$PWD/PackagesAndFiles $$OUT_PWD
    first.depends = $(first) copypackages
    export(first.depends)
    export(copypackages.commands)
    QMAKE_EXTRA_TARGETS += first copypackages
}
