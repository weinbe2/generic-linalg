.SUFFIXES:
.SUFFIXES: .o .c

#============================================================

C_SOURCES =  ../../u1_utils/u1_utils.cpp ../../operator_utils/operators.cpp ../../generic_gcr.cpp ../../generic_minres.cpp ../../generic_gcr_var_precond.cpp ../../generic_precond.cpp
C_OBJS    = u1_utils.o operators.o generic_gcr.o generic_minres.o generic_gcr_var_precond.o generic_precond.o
C_INCLUDES  = generic_cg.h generic_bicgstab_l.h generic_vector.h
C_LIBS    = -lm -lrt
EXAMPLES = staggered_w_laplace


CCX = g++
#CXXFLAGS =   -std=c99 -g -O2  
CXXFLAGS = -O0 -g -Wall -std=c++11 -I../../ -I../../u1_utils/ -I../../operator_utils/ 


#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================

all: $(EXAMPLES)

%.o:%.cpp	$(INCLUDES)
	$(CCX) -c  $(CXXFLAGS) $(C_SOURCES)

staggered_w_laplace : staggered_w_laplace.o $(C_OBJS)
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


