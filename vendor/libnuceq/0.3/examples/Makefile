#///////////////////////////////////////////////////////////////////////////////
# <file type="public">
#
#   <description>
#     <abstract>
#       Makefile for compiling the example codes.
#     </abstract>
#   </description>
#   <license>
#     See the src/README.txt file in this distribution for copyright and license
#     information.
#   </license>
#
# </file>
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
# The following lines must be edited according to where you have
# placed your libnucnet, wn_matrix, and libstatmech modules.  You may also
# edit to choose a different compiler (e.g., g++), to use Valgrind or
# not, or to profile:
#///////////////////////////////////////////////////////////////////////////////

GC=gcc

MODULES = ../../..

LIBNUCEQDIR = ../src/
LIBNUCNETDIR = $(MODULES)/libnucnet/0.11/src/
LIBSTATMECHDIR = $(MODULES)/libstatmech/0.5/src/
MATRIXSRCDIR = $(MODULES)/wn_matrix/0.13/src/
VALGRIND= yes
PROFILE= no

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

VPATH = $(MATRIXSRCDIR):$(LIBNUCNETDIR):$(LIBSTATMECHDIR):$(LIBNUCEQDIR)

CINCLUDE= `xml2-config --cflags` `gsl-config --cflags` -I$(LIBNUCNETDIR) -I$(LIBNUCEQDIR) -I$(MATRIXSRCDIR) -I$(LIBSTATMECHDIR)
CLIBS= `xml2-config --libs` `gsl-config --libs`

#===============================================================================
# Compiler flags.
#===============================================================================

CFLAGS1= -ansi -Werror -Wall -W -pedantic\
         -Wconversion -Wshadow \
         -Wpointer-arith -Wcast-qual -Wcast-align \
         -Wwrite-strings \
         -fshort-enums -fno-common -Dinline= -g \

ifeq ($(GC), gcc) 
	CFLAGS2= $(CFLAGS1) -Wmissing-prototypes -Wstrict-prototypes -Wnested-externs
else
	CFLAGS2= $(CFLAGS1)
endif

ifeq ($(VALGRIND), yes)
	CFLAGS3= $(CFLAGS2) -O0
else
	CFLAGS3= $(CFLAGS2) -O2
endif

ifeq ($(PROFILE), yes)
	CFLAGS= $(CFLAGS3) -pg
else
	CFLAGS= $(CFLAGS3)
endif

CC=$(GC) $(CFLAGS) $(CINCLUDE)

#===============================================================================
# TMPDIR is the temporary directory for codes compilation, this is where
# object files are created. 
#===============================================================================

TMPDIR = ./tmp/

TMP_DIR := $(shell mkdir -p $(TMPDIR))

#===============================================================================
# Compile auxiliary codes.
#===============================================================================

OBJ =   $(TMPDIR)Libnucnet__Nuc.o \
	$(TMPDIR)Libstatmech.o    \
	$(TMPDIR)WnMatrix.o       \
	$(TMPDIR)coul_corr.o      \
	$(TMPDIR)Libnuceq.o

$(OBJ): $(TMPDIR)%.o: %.c
	$(CC) -c -o $@ $<

#===============================================================================
# Compile all.
#===============================================================================

EXECS = compute_wse                    \
        compute_wse_from_lepton_number \
	compute_function_integrand_wse \
	compare_wse_screening          \
	compute_nse                    \
	compute_screened_degen_nse     \
	compute_qse                    \
	compute_many_qse               \
	compute_ng                     \
	compute_rnse

$(EXECS) : $(OBJ)
	$(CC) -c -o $(TMPDIR)$@.o $@.c
	$(CC) $(OBJ) $(TMPDIR)$@.o $(CLIBS) -o $@

.PHONY all : $(EXECS)
     
#===============================================================================
# Clean up.
#===============================================================================

.PHONY: clean cleanall

clean: 
	rm -fr $(TMPDIR)

cleanall: clean
	rm -f $(EXECS) *.exe

