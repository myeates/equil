#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2011-2012 Clemson University.
# 
#  This file was originally written by Bradley S. Meyer.
# 
#  This is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#  USA
# 
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file Makefile
#//! \brief A makefile to generate dependencies.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef BUILD_DEF

GC=g++

LIBNUCNET_VERSION = 0.16
LIBNUCEQ_VERSION = 0.3
LIBSTATMECH_VERSION = 0.8
WNMATRIX_VERSION = 0.14

SOURCEFORGE_URL = http://sourceforge.net/projects/

LIBNUCNET_URL = $(SOURCEFORGE_URL)/libnucnet/files/libnucnet
LIBNUCEQ_URL = $(SOURCEFORGE_URL)/libnuceq/files/libnuceq
LIBSTATMECH_URL = $(SOURCEFORGE_URL)/libstatmech/files/libstatmech
WNMATRIX_URL = $(SOURCEFORGE_URL)/wnmatrix/files/wn_matrix

LIBNUCNETDIR = $(VENDORDIR)/libnucnet/$(LIBNUCNET_VERSION)/src
LIBNUCEQDIR = $(VENDORDIR)/libnuceq/$(LIBNUCEQ_VERSION)/src
LIBSTATMECHDIR = $(VENDORDIR)/libstatmech/$(LIBSTATMECH_VERSION)/src
WNMATRIXDIR = $(VENDORDIR)/wn_matrix/$(WNMATRIX_VERSION)/src

BOOST_DIR = ./

VALGRIND= yes
PROFILE= no
OPENMP= yes
NODEPRECATED= yes

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

VPATH= $(WNMATRIXDIR):$(LIBNUCNETDIR):$(LIBNUCEQDIR):$(LIBSTATMECHDIR)

CINCLUDE= `xml2-config --cflags` \
          `gsl-config --cflags`  \
          -I$(LIBNUCNETDIR)      \
          -I$(WNMATRIXDIR)       \
          -I$(LIBSTATMECHDIR)    \
          -I$(LIBNUCEQDIR)       \
          -I$(NUCNET_TARGET)     \
          -I$(BOOST_DIR)         \

CLIBS= `xml2-config --libs` `gsl-config --cflags --libs`

#===============================================================================
# Compiler flags.
#===============================================================================

CFLAGS= -ansi -Werror -Wall \
         -Wpointer-arith -Wcast-qual \
         -Wwrite-strings -Wconversion \
         -fshort-enums -fno-common -g \

ifeq ($(GC), icpc)
	CFLAGS+= -wd9 -wd981 -wd1292 -wd1418 -wd1419 -wd10148 -wd10156
endif

ifeq ($(VALGRIND), yes)
	CFLAGS+= -O0
else
	CFLAGS+= -O2
endif

ifeq ($(PROFILE), yes)
	CFLAGS+= -pg
endif

ifeq ($(OPENMP), yes)
  ifeq ($(GC), g++)
	CFLAGS+= -fopenmp
  endif
  ifeq ($(GC), icpc)
	CFLAGS+= -openmp
 endif
endif

ifeq ($(NODEPRECATED), yes)
	CFLAGS+= -Wno-deprecated
endif

CC=$(GC) $(CFLAGS) $(CINCLUDE)

#===============================================================================
# OBJDIR is the temporary directory for codes compilation, this is where
# object files are created. 
#===============================================================================

ifndef OBJDIR
OBJDIR = ../../obj
endif
OBJ_DIR := $(shell mkdir -p $(OBJDIR))

#===============================================================================
# BINDIR is the temporary directory for code executables.
#===============================================================================

ifndef BINDIR
BINDIR = .
endif
BIN_DIR := $(shell mkdir -p $(BINDIR))

#===============================================================================
# VENDORDIR is the directory for storing webnucleo codes.
#===============================================================================

ifndef VENDORDIR
VENDORDIR = ../../vendor
endif
VENDOR_DIR := $(shell mkdir -p $(VENDORDIR))

#===============================================================================
# Build dependencies.
#===============================================================================

$(VENDORDIR)/wn_matrix.tar.gz:
	wget $(WNMATRIX_URL)/wn_matrix_$(WNMATRIX_VERSION).tar.gz -O $@

$(OBJDIR)/WnMatrix.o: $(VENDORDIR)/wn_matrix.tar.gz 
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(WNMATRIXDIR)/WnMatrix.c -o $@

$(VENDORDIR)/libstatmech.tar.gz:
	wget $(LIBSTATMECH_URL)/libstatmech_$(LIBSTATMECH_VERSION).tar.gz -O $@

$(OBJDIR)/Libstatmech.o: $(VENDORDIR)/libstatmech.tar.gz 
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBSTATMECHDIR)/Libstatmech.c -o $@

$(VENDORDIR)/libnuceq.tar.gz:
	wget $(LIBNUCEQ_URL)/libnuceq_$(LIBNUCEQ_VERSION).tar.gz -O $@

$(OBJDIR)/Libnuceq.o: $(VENDORDIR)/libnuceq.tar.gz $(OBJDIR)/Libnucnet.o
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBNUCEQDIR)/Libnuceq.c -o $@

$(VENDORDIR)/libnucnet.tar.gz:
	wget $(LIBNUCNET_URL)/libnucnet_$(LIBNUCNET_VERSION).tar.gz -O $@

$(OBJDIR)/Libnucnet__Nuc.o: $(VENDORDIR)/libnucnet.tar.gz $(OBJDIR)/WnMatrix.o
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBNUCNETDIR)/Libnucnet__Nuc.c -o $@

$(OBJDIR)/Libnucnet__Reac.o: $(VENDORDIR)/libnucnet.tar.gz $(OBJDIR)/WnMatrix.o
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBNUCNETDIR)/Libnucnet__Reac.c -o $@

$(OBJDIR)/Libnucnet.o: $(VENDORDIR)/libnucnet.tar.gz $(OBJDIR)/WnMatrix.o
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBNUCNETDIR)/Libnucnet.c -o $@

#--------------------------------------------------------------------------

WN_OBJ =$(OBJDIR)/WnMatrix.o		\
	$(OBJDIR)/Libnucnet__Nuc.o	\
	$(OBJDIR)/Libnucnet__Reac.o	\
	$(OBJDIR)/Libnucnet.o		\
	$(OBJDIR)/Libstatmech.o		\
	$(OBJDIR)/Libnuceq.o
#===============================================================================
# Get data.
#===============================================================================

ifndef DATA_DIR

DATA_URL=http://nucnet-tools.sourceforge.net/data_pub/2012-07-07
DATA=nnt_data.tar
DATA_DIR=./

data:
	wget $(DATA_URL)/$(DATA).gz
	tar xz -C ${DATA_DIR} -f ${DATA_DIR}/${DATA}.gz
	rm ${DATA_DIR}/${DATA}.gz

endif

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean clean_all clean_dist

clean:
	rm -fr $(OBJDIR)

clean_dist:
	rm -fr $(VENDORDIR)

clean_all: clean clean_dist

clean_data:
	rm -fr $(DATA_DIR)

#===============================================================================
# End include.
#===============================================================================

BUILD_DEF = yes

endif

