QMAKE_CXXFLAGS += -std=c++11

SOURCES += \
    main.cpp \
    core/structuralparameters.cpp \
    core/model.cpp \
    core/multilevelparameters.cpp \
    core/functions.cpp \

HEADERS += \
    core/structuralparameters.hpp \
    core/model.hpp \
    core/multilevelparameters.hpp \
    core/montecarlo.hpp \
    core/functions.hpp \
    core/estimator.hpp \
    core/scheme.hpp

INCLUDEPATH += /usr/include/eigen core
