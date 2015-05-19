TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

SOURCES += main.cpp \
    atom.cpp \
    readlammpsdump.cpp \
    test_atoms.cpp \
    diffusion.cpp \
    density.cpp

HEADERS += \
    atom.h \
    readlammpsdump.h \
    test_atoms.h \
    diffusion.h \
    density.h

QMAKE_CXXFLAGS += -std=c++11
