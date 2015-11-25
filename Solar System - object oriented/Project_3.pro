LIBS += -L/usr/local/lib -larmadillo
INCLUDEPATH += /usr/local/include
DEPENDPATH += /usr/local/include

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    planet.cpp \
    solarsystem.cpp

HEADERS += \
    planet.h \
    solarsystem.h

