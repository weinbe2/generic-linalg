.SUFFIXES:
.SUFFIXES: .o .c

#============================================================

C_SOURCES =  ../../u1_utils/u1_utils.cpp ../../operator_utils/operators.cpp ../../operator_utils/operators_stencil.cpp ../../stencil_2d/coarse_stencil.cpp
C_OBJS    = u1_utils.o operators.o operators_stencil.o coarse_stencil.o
C_INCLUDES  = generic_cg.h generic_bicgstab.h generic_bicgstab_l.h generic_gcr.h generic_vector.h
C_LIBS    = -lm -lrt
EXAMPLES = staggered_stencil


CCX = g++
#CXXFLAGS =   -std=c99 -g -O2  
CXXFLAGS = -O0 -g -Wall -std=c++11 -I../../ -I../../u1_utils/ -I../../operator_utils/ -I../../lattice/ -I../../stencil_2d/


#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================

all: $(EXAMPLES)

%.o:%.cpp	$(INCLUDES)
	$(CCX) -c  $(CXXFLAGS) $(C_SOURCES)

staggered_stencil : staggered_stencil.o $(C_OBJS)
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


