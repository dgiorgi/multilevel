QMAKE_CXXFLAGS += -pthread -std=c++11

LIBS += -L$$OUT_PWD/../../core/ -lcore

SOURCES += \
    main.cpp

INCLUDEPATH += $$PWD/../../core
DEPENDPATH += $$PWD/../../coreQMAKE_CXXFLAGS += -pthread -std=c++11

INCLUDEPATH += /usr/include/eigen
