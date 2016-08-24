.SUFFIXES:
.SUFFIXES: .o .c

#============================================================

C_SOURCES =  ../generic_bicgstab.cpp ../u1_utils/u1_utils.cpp ../operator_utils/operators.cpp ../arpack_binding/arpack_interface.cpp
C_OBJS    = generic_bicgstab.o u1_utils.o operators.o arpack_interface.o
C_INCLUDES  = generic_inverters.h generic_vector.h operators.h
C_LIBS    = -lm -larpack
EXAMPLES = eigenspectrum


CCX = g++
#CXXFLAGS =   -std=c99 -g -O2  
CXXFLAGS = -O0 -g -Wall -std=c++11 -I../ -I../u1_utils/ -I../operator_utils/ -I../arpack_binding/ -L/media/data/weinbe2/bin/arpack/install/lib/


#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================

all: $(EXAMPLES)

%.o:%.cpp	$(INCLUDES)
	$(CCX) -c  $(CXXFLAGS) $(C_SOURCES)

eigenspectrum : eigenspectrum.o $(C_OBJS)
	$(CCX) $@.cpp $(CXXFLAGS)  -o $@ $(C_OBJS) $(C_LIBS)	

#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(INCLUDES) $(EXAMPLES)

clean: 
	rm -f $(TARGET) $(C_OBJS) $(EXAMPLES) $(C_EXAMPLES_OBJ) core

TAGS:	$(ALL_SOURCES)
	etags $(ALL_SOURCES)

tar: $(ALL_SOURCES) 
	tar cvf $(TARGET).tar $(ALL_SOURCES) 

ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


