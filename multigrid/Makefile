.SUFFIXES:
.SUFFIXES: .o .c

#============================================================

C_SOURCES =  ../generic_cg.cpp ../generic_bicgstab.cpp ../generic_gmres.cpp ../generic_gelim.cpp ../generic_sor.cpp
C_OBJS    = generic_bicgstab.o generic_cg.o generic_gmres.o generic_gelim.o generic_sor.o
C_INCLUDES  = generic_inverters.h generic_vector.h
C_LIBS    = -lm
EXAMPLES = mg_square_laplace mg_square_laplace_sor


CCX = g++
#CXXFLAGS =   -std=c99 -g -O2  
CXXFLAGS = -O2 -std=c++11 -I../


#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================

all: $(EXAMPLES)

%.o:%.cpp	$(INCLUDES)
	$(CCX) -c  $(CXXFLAGS) $(C_SOURCES)

mg_square_laplace : mg_square_laplace.o $(C_OBJS)
	$(CCX) $@.cpp $(CXXFLAGS)  -o $@ $(C_OBJS) $(C_LIBS)

mg_square_laplace_sor : mg_square_laplace_sor.o $(C_OBJS)
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


