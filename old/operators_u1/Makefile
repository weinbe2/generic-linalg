.SUFFIXES:
.SUFFIXES: .o .c

#============================================================

C_SOURCES =  ../generic_cg.cpp ../generic_cr.cpp ../generic_bicgstab.cpp ../generic_gcr.cpp ../generic_gmres.cpp ../generic_gelim.cpp ../generic_sor.cpp ../generic_minres.cpp ../generic_precond.cpp ../generic_cg_precond.cpp ../generic_cg_flex_precond.cpp ../generic_gcr_var_precond.cpp ../u1_utils/u1_utils.cpp
C_OBJS    = generic_bicgstab.o generic_cr.o generic_cg.o generic_gcr.o generic_gmres.o generic_gelim.o generic_sor.o generic_minres.o generic_precond.o generic_cg_precond.o generic_cg_flex_precond.o generic_gcr_var_precond.o u1_utils.o
C_INCLUDES  = generic_inverters.h generic_inverters_precond.h generic_vector.h
C_LIBS    = -lm
EXAMPLES = square_laplace square_laplace_u1 square_staggered square_staggered_u1


CCX = g++
#CXXFLAGS =   -std=c99 -g -O2  
CXXFLAGS = -O2 -std=c++11 -I../ -I../u1_utils/


#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================

all: $(EXAMPLES)

%.o:%.cpp	$(INCLUDES)
	$(CCX) -c  $(CXXFLAGS) $(C_SOURCES)

square_staggered : square_staggered.o $(C_OBJS)
	$(CCX) $@.cpp $(CXXFLAGS)  -o $@ $(C_OBJS) $(C_LIBS)

square_staggered_u1 : square_staggered_u1.o $(C_OBJS)
	$(CCX) $@.cpp $(CXXFLAGS)  -o $@ $(C_OBJS) $(C_LIBS)

square_laplace : square_laplace.o $(C_OBJS)
	$(CCX) $@.cpp $(CXXFLAGS)  -o $@ $(C_OBJS) $(C_LIBS)

square_laplace_u1 : square_laplace_u1.o $(C_OBJS)
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


