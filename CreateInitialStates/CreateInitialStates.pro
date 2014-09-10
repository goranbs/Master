include(./defaults.pri)
TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    moleculesystem.cpp \
    atom.cpp \
    bond.cpp \
    angle.cpp

HEADERS += \
    moleculesystem.h \
    atom.h \
    bond.h \
    angle.h

