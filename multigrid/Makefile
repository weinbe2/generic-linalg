.SUFFIXES:
.SUFFIXES: .o .c

#============================================================
TARGET	= mg_square_laplace

C_SOURCES =  mg_square_laplace.cpp ../generic_cg.cpp ../generic_bicgstab.cpp ../generic_gmres.cpp ../generic_gelim.cpp ../generic_sor.cpp
C_OBJS    = mg_square_laplace.o generic_bicgstab.o generic_cg.o generic_gmres.o generic_gelim.o generic_sor.o

CCX = g++
#CXXFLAGS =   -std=c99 -g -O2  
CXXFLAGS = -O2 -g -std=c++11 -I../


#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================
all: $(TARGET)

%.o:%.cpp	$(INCLUDES)
	$(CCX) -c  $(CXXFLAGS) $(C_SOURCES) 

$(TARGET) :   $(C_OBJS)
	$(CCX) $(C_OBJS) $(CXXFLAGS)  -o $(TARGET) -lm 

#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(C_INCLUDES)

clean:
	rm -f $(TARGET) $(C_OBJS) core 

TAGS:	$(ALL_SOURCES)
	etags $(ALL_SOURCES)

tar: $(ALL_SOURCES) 
	tar cvf $(TARGET).tar $(ALL_SOURCES) 

ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


