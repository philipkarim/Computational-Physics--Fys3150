TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        funksjoner.cpp \
        jacobi_ex2a.cpp \
        jacobimetode.cpp \
        main.cpp \
        test.cpp \
        unixtest.cpp

INCLUDEPATH += C:\Users\Philip\Desktop\armadillo-9.700.2\include
DEPENDPATH += C:\Users\Philip\Desktop\armadillo-9.700.2\include

LIBS += \
    -LC:\Users\Philip\Desktop\armadillo-9.700.2\examples\lib_win64 \
    -llapack_win64_MT \
    -lblas_win64_MT

HEADERS += \
    prosjekt2funksjoner.h \
    unixtest.h
