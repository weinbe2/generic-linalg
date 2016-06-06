.SUFFIXES:
.SUFFIXES: .o .cpp

#============================================================

C_SOURCES =  u1_utils.cpp
C_OBJS    =  u1_utils.o
C_INCLUDES  = u1_utils.h
C_LIBS    = -lm
EXAMPLES = u1_test


CCX = g++
#CXXFLAGS =   -std=c99 -g -O2  
CXXFLAGS = -O2 -std=c++11


#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================

all: $(EXAMPLES)

%.o:%.cpp	$(INCLUDES)
	$(CCX) -c  $(CXXFLAGS) $(*:=.cpp)

u1_test :  $(C_OBJS)
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


