QMAKE_CXXFLAGS += -fopenmp -Wall -pthread -std=c++11
QMAKE_LFLAGS *= -fopenmp

LIBS += -L$$OUT_PWD/../../core/ -lcore

SOURCES += \
    main.cpp

INCLUDEPATH += $$PWD/../../core
DEPENDPATH += $$PWD/../../core

INCLUDEPATH += /usr/include/eigen

