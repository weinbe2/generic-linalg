.SUFFIXES:
.SUFFIXES: .o .c

#============================================================

C_SOURCES =  ../../generic_cg.cpp ../../generic_cg_m.cpp ../../generic_cr.cpp ../../generic_cr_m.cpp ../../generic_bicgstab.cpp ../../generic_bicgstab_l.cpp ../../generic_gcr.cpp ../../u1_utils/u1_utils.cpp ../../operator_utils/operators.cpp
C_OBJS    = generic_cg.o generic_cg_m.o generic_cr.o generic_cr_m.o generic_bicgstab.o generic_bicgstab_l.o generic_gcr.o u1_utils.o operators.o
C_INCLUDES  = generic_ch.h generic_bicgstab.h generic_bicgstab_l.h generic_gcr.h generic_vector.h
C_LIBS    = -lm -lrt
EXAMPLES = multishift


CCX = g++
#CXXFLAGS =   -std=c99 -g -O2  
CXXFLAGS = -O0 -g -Wall -std=c++11 -I../../ -I../../u1_utils/ -I../../operator_utils/


#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================

all: $(EXAMPLES)

%.o:%.cpp	$(INCLUDES)
	$(CCX) -c  $(CXXFLAGS) $(C_SOURCES)

multishift : multishift.o $(C_OBJS)
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


