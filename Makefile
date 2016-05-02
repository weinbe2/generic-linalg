.SUFFIXES:
.SUFFIXES: .o .c

#============================================================

C_SOURCES =  generic_cg.cpp generic_bicgstab.cpp generic_gmres.cpp generic_gelim.cpp generic_sor.cpp
C_OBJS    = generic_bicgstab.o generic_cg.o generic_gmres.o generic_gelim.o generic_sor.o
C_INCLUDES  = generic_inverters.h generic_vector.h
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


