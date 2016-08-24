.SUFFIXES:
.SUFFIXES: .o .c

#============================================================

C_SOURCES = arpack_interface.cpp ../generic_bicgstab.cpp
C_OBJS    = arpack_interface.o generic_bicgstab.o
C_INCLUDES  = arpack_interface.h arpack_extern.h generic_vector.h
C_LIBS    = -lm -larpack
LIB       = -L/media/data/weinbe2/bin/arpack/install/lib/
INCL      = -I../
TARGET    = rzn_example

CCX = g++
CXXFLAGS =  -std=c++11 -O2 -Wall $(LIB) $(INCL)
#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================
all: rzn_example

%.o:%.cpp	$(INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $(C_SOURCES)

rzn_example : rzn_example.o $(C_OBJS)
	$(CCX) $@.cpp $(CXXFLAGS)  -o $@ $(C_OBJS) $(C_LIBS)


#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(INCLUDES) $(TARGET)

clean:
	rm -f $(TARGET) $(C_OBJS) rzn_example.o core 

TAGS:	$(ALL_SOURCES)
	etags $(ALL_SOURCES)

tar: $(ALL_SOURCES) 
	tar cvf $(TARGET).tar $(ALL_SOURCES) 

ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


