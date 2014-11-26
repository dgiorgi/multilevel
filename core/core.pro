QMAKE_CXXFLAGS += -pthread -std=c++11

QT       -= core gui

TARGET = core
TEMPLATE = lib

DEFINES += CORE_LIBRARY

SOURCES += \
    structuralparameters.cpp \
    model.cpp \
    multilevelparameters.cpp \
    functions.cpp \

HEADERS += \
    structuralparameters.hpp \
    model.hpp \
    multilevelparameters.hpp \
    montecarlo.hpp \
    functions.hpp \
    estimator.hpp \
    scheme.hpp

INCLUDEPATH += /usr/include/eigen

