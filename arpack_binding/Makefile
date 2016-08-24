.SUFFIXES:
.SUFFIXES: .o .c

#============================================================
TARGET	= rzn_example

C_SOURCES = rzn_example.c arpack_interface.c
C_OBJS    = rzn_example.o arpack_interface.o
INCLUDES  = arpack_interface.h arpack_extern.h
LIB       = -L/media/data/weinbe2/bin/arpack/install/lib/
CCX = gcc
CXXFLAGS =  -std=c99 -g -O2 $(LIB)
#  = -DHAVE_CONFIG_H -I   -g -O2  -MT
# -g -O2  $(INC)

#============================================================
all: $(TARGET)

.o:.c	$(INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $(C_SOURCES) 

$(TARGET) :   $(C_OBJS)
	$(CCX) $(C_OBJS) $(CXXFLAGS)  -o $(TARGET) -lm -larpack

#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(INCLUDES)

clean:
	rm -f $(TARGET) $(C_OBJS) core 

TAGS:	$(ALL_SOURCES)
	etags $(ALL_SOURCES)

tar: $(ALL_SOURCES) 
	tar cvf $(TARGET).tar $(ALL_SOURCES) 

ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


