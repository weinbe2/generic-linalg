.SUFFIXES:
.SUFFIXES: .o .c

#============================================================

C_SOURCES =  ../../generic_cg.cpp ../../generic_bicgstab_l.cpp ../../u1_utils/u1_utils.cpp ../../operator_utils/operators.cpp ../../operator_utils/operators_stencil.cpp ../../stencil_2d/coarse_stencil.cpp ../../multigrid/aa_mg/null_gen.cpp ../../multigrid/aa_mg/mg_complex.cpp ../../generic_cr.cpp ../../generic_bicgstab.cpp ../../generic_gmres.cpp ../../generic_gcr.cpp ../../generic_minres.cpp ../../generic_sor.cpp ../../generic_gelim.cpp ../../generic_inverter.cpp ../../generic_cg_flex_precond.cpp ../../generic_bicgstab_precond.cpp ../../generic_gcr_var_precond.cpp
C_OBJS    = generic_cg.o generic_bicgstab_l.o u1_utils.o operators.o operators_stencil.o coarse_stencil.o null_gen.o mg_complex.o generic_cr.o generic_bicgstab.o generic_gmres.o generic_gcr.o generic_minres.o generic_sor.o generic_gelim.o generic_inverter.o generic_cg_flex_precond.o generic_bicgstab_precond.o generic_gcr_var_precond.o
C_INCLUDES  = generic_cg.h generic_bicgstab_l.h generic_vector.h
C_LIBS    = -lm -lrt
EXAMPLES = staggered_pieces


CCX = g++
#CXXFLAGS =   -std=c99 -g -O2  
CXXFLAGS = -O0 -g -Wall -std=c++11 -I../../ -I../../u1_utils/ -I../../operator_utils/ -I../../lattice/ -I../../stencil_2d/ -I../../multigrid/aa_mg/


#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================

all: $(EXAMPLES)

%.o:%.cpp	$(INCLUDES)
	$(CCX) -c  $(CXXFLAGS) $(C_SOURCES)

staggered_pieces : staggered_pieces.o $(C_OBJS)
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


