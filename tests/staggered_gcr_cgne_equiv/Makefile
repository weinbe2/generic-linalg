.SUFFIXES:
.SUFFIXES: .o .c

#============================================================

C_SOURCES =  ../../generic_gelim.cpp ../../generic_gmres.cpp ../../generic_cg.cpp ../../generic_gcr.cpp ../../u1_utils/u1_utils.cpp ../../operator_utils/operators.cpp
C_OBJS    = generic_gelim.o generic_gmres.o generic_cg.o generic_gcr.o u1_utils.o operators.o
C_INCLUDES  = generic_cg.h generic_gcr.h generic_vector.h
C_LIBS    = -lm -lrt
EXAMPLES = gcr_cgne_equiv


CCX = g++
#CXXFLAGS =   -std=c99 -g -O2  
CXXFLAGS = -O0 -g -Wall -std=c++11 -I../../ -I../../u1_utils/ -I../../operator_utils/


#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================

all: $(EXAMPLES)

%.o:%.cpp	$(INCLUDES)
	$(CCX) -c  $(CXXFLAGS) $(C_SOURCES)

gcr_cgne_equiv : gcr_cgne_equiv.o $(C_OBJS)
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


