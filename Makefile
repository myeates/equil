# simple makefile for a Rappture-based program

RAPPTURE_DIR	= /Users/bradleymeyer/Desktop/work/rappture/20120712/
INCLUDES	= -I$(RAPPTURE_DIR)/include
LIBS		= -L$(RAPPTURE_DIR)/lib -lrappture -lm

VENDORDIR = ./vendor
OBJDIR = ./obj

include Makefile.inc

all: mainc

mainc: main.cpp $(WN_OBJ)
	$(CC) $(CFLAGS) $(WN_OBJ) $(INCLUDES) $< -o $@ $(LIBS) $(CLIBS)

clean_runs:
	$(RM) run*.xml driver*.xml

clean_rappture: clean_runs
	$(RM) mainc 

clean_all: clean_rappture
	rm -fr VENDORDIR OBJDIR
