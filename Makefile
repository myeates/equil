# simple makefile for a Rappture-based program

RAPPTURE_DIR	= /home/myeates/local/rappture/20130903
INCLUDES	= -I$(RAPPTURE_DIR)/include
LIBS		= -L$(RAPPTURE_DIR)/lib -lrappture -lm

VENDORDIR = ./vendor
OBJDIR = ./obj

include Makefile.inc

all: mainc

mainc: main.cpp $(WN_OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@ $(WN_OBJ) $(LIBS) $(CLIBS)

clean_runs:
	$(RM) run*.xml driver*.xml

clean_rappture: clean_runs
	$(RM) main 

clean_all: clean_rappture
	rm -fr $(OBJDIR)
