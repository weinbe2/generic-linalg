.SUFFIXES:
.SUFFIXES: .o .cpp

#============================================================
TARGET	= rzn_example

C_SOURCES = arpack_interface.cpp
C_OBJS    = arpack_interface.o
INCLUDES  = arpack_interface.h arpack_extern.h
C_LIBS    = -lm -larpack
LIB       = -L/media/data/weinbe2/bin/arpack/install/lib/
CCX = g++
CXXFLAGS =  -std=c++11 -O2 $(LIB)
#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================
all: rzn_example

%.o:%.cpp	$(INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $(*:=.cpp)

rzn_example :  $(C_OBJS)
	$(CCX) $@.cpp $(CXXFLAGS)  -o $@ $(C_OBJS) $(C_LIBS)


#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(INCLUDES)

clean:
	rm -f $(TARGET) $(C_OBJS) core 

TAGS:	$(ALL_SOURCES)
	etags $(ALL_SOURCES)

tar: $(ALL_SOURCES) 
	tar cvf $(TARGET).tar $(ALL_SOURCES) 

ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


