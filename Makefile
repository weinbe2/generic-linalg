.SUFFIXES:
.SUFFIXES: .o .c

#============================================================

C_SOURCES =  generic_cg.cpp generic_cr.cpp generic_bicgstab.cpp generic_bicgstab_l.cpp generic_gcr.cpp generic_gmres.cpp generic_gelim.cpp generic_sor.cpp generic_minres.cpp generic_cg_precond.cpp generic_cg_flex_precond.cpp generic_gcr_var_precond.cpp generic_bicgstab_precond.cpp generic_poweriter.cpp generic_precond.cpp generic_inverter.cpp generic_inverter_precond.cpp
C_OBJS    = generic_bicgstab.o generic_bicgstab_l.o generic_cg.o generic_cr.o generic_gcr.o generic_gmres.o generic_gelim.o generic_sor.o generic_minres.o generic_cg_precond.o generic_cg_flex_precond.o generic_gcr_var_precond.o generic_bicgstab_precond.o generic_poweriter.o generic_precond.o generic_inverter.o generic_inverter_precond.o
C_INCLUDES  = generic_inverters.h generic_inverters_precond.h generic_vector.h generic_eigenvalues.h verbosity.h generic_traits.h
C_LIBS    = -lm
EXAMPLES = square_laplace imag_laplace unit_test


CCX = g++
#CXXFLAGS =   -std=c99 -g -O2  
CXXFLAGS = -O2 -std=c++11


#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================

all: $(EXAMPLES)

%.o:%.cpp	$(INCLUDES)
	$(CCX) -c  $(CXXFLAGS) $(*:=.cpp)

square_laplace :  $(C_OBJS)
	$(CCX) $@.cpp $(CXXFLAGS)  -o $@ $(C_OBJS) $(C_LIBS)

imag_laplace : $(C_OBJS)
	$(CCX) $@.cpp $(CXXFLAGS)  -o $@ $(C_OBJS) $(C_LIBS)

unit_test : $(C_OBJS)
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


