TEMPLATE = app 
GRAPHICS = x11
CONFIG += console 
CONFIG += release
QT += widgets
QT += gui
CONFIG -= debug
CONFIG -= app_bundle


contains( GRAPHICS, qt ) {
  
}	

TARGET = parallel
MAINFILE = $$join(TARGET, " ", , ".cpp" )

message( $$MAINFILE )
message( $$TARGET )
# Input
HEADERS += ca.h \
	   hull.h \
           cell.h \
           conrec.h \
           dish.h \
           graph.h \
           misc.h \
           output.h \
           parameter.h \
           parse.h \
           pde.h \
           random.h \
           sqr.h \
           sticky.h \
       	   crash.h \
	   warning.h 

        
SOURCES += ca.cpp \
	   hull.cpp \
           cell.cpp \
           conrec.cpp \
           dish.cpp \
           misc.cpp \
           output.cpp \
           parameter.cpp \
           parse.cpp \
           pde.cpp \
           random.cpp \
           crash.cpp \
           warning.cpp 

SOURCES += $$MAINFILE
       
#QMAKE_CXXFLAGS_RELEASE += -fexceptions -fopenmp
#QMAKE_CXXFLAGS_DEBUG += -fexceptions
#QMAKE_LFLAGS_RELEASE += -O4 -fopenmp
#QMAKE_CXXFLAGS_RELEASE += -O4




#The following line was inserted by qt3to4
QT +=  
