# $Id: Makefile.in 264 2010-12-08 20:41:23Z jared.crossley $
#-----------------------------------------------------------------------
#;  Copyright (C) 2004-2008
#;  Associated Universities, Inc. Washington DC, USA.
#;
#;  This program is free software; you can redistribute it and/or
#;  modify it under the terms of the GNU General Public License as
#;  published by the Free Software Foundation; either version 2 of
#;  the License, or (at your option) any later version.
#;
#;  This program is distributed in the hope that it will be useful,
#;  but WITHOUT ANY WARRANTY; without even the implied warranty of
#;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#;  GNU General Public License for more details.
#;
#
#;  You should have received a copy of the GNU General Public
#;  License along with this program; if not, write to the Free
#;  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#;  MA 02139, USA.
#;
#;  Correspondence concerning this software should be addressed as follows:
#;         Internet email: bcotton@nrao.edu
#;         Postal address: W. D. Cotton
#;                         National Radio Astronomy Observatory
#;                         520 Edgemont Road
#;                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------
#    Makefile to build library routines
#    Obit software
#
#   Expects CC, CFLAGS, to be set by either environment variable 
#   or Make argument
#     CC is the compiler to use
#     CFLAGS is compile options
#     INCDIR directories with include files
#     LIB defines libraries and search paths,
#
#------------------------------------------------------------------------

top_srcdir = ..

CC = gcc
#CPPFLAGS = 
CPPFLAGS = 
CFLAGS = -g -O3 -fPIC -msse
LDFLAGS = 

ALL_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/include -I/usr/include/cfitsio \
	 -I/usr/include    -I/usr/include \
	  -DPACKAGE_NAME=\"Obit\" -DPACKAGE_TARNAME=\"obit\" -DPACKAGE_VERSION=\"1.0\" -DPACKAGE_STRING=\"Obit\ 1.0\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obit\" -DVERSION=\"1.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_SSE=1 -DHAVE__TMP_DUMMY2_FITSIO_H=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE_FFTW3=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -DHAVE_PGPLOT=1 -DHAVE_PYTHON=1 -DHAVE_FSEEKO=1 -DOBIT_THREADS_ENABLED=1
ALL_CFLAGS = $(CFLAGS) -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include    \
	  
ALL_LDFLAGS = $(LDFLAGS)       \
	  -L/usr/lib/python2.7/site-packages 

LIBS = ../lib/libObit.a -lcfitsio    -lfftw3f -lglib-2.0  \
	-lgsl -lgslcblas   -lcpgplot -lpgplot    -lSM -lICE -lX11    -L/usr/lib/gcc/x86_64-redhat-linux/4.9.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.9.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.9.2/../../.. -lgfortran -lm -lquadmath -lgthread-2.0 -lpthread -lglib-2.0 

CLIENT_CPPFLAGS = $(ALL_CPPFLAGS) -I../dummy_xmlrpc 
CLIENT_CFLAGS = $(ALL_CFLAGS) 
CLIENT_LDFLAGS = $(ALL_LDFLAGS) -L../dummy_xmlrpc  
CLIENT_LIBS = $(LIBS) -lz  ../dummy_xmlrpc/libxmlrpc.a 

SERVER_CPPFLAGS = -I../dummy_xmlrpc
SERVER_CFLAGS = 
SERVER_LDFLAGS = $(LDFLAGS) -L../dummy_xmlrpc
SERVER_LIBS =  ../dummy_xmlrpc/libxmlrpc.a 

# Must use SWIG v1.1-883
SWIG = 

# Libraries in case they've changed
MYLIBS := $(wildcard ../lib/lib*.a)

# Do everything in one big module
TARGETS := Obit.so 

all: $(TARGETS)

# Build shared library for python interface
$(TARGETS): setupdata.py $(MYLIBS)
	rm -rf build
	python makesetup.py
	python setup.py build install --install-lib=.

# Build python/Obit interface
interface: Obit_wrap.c
	echo "rebuild Obit/python interface"

# Use Swig to build python/c interface if necessary
Obit_wrap.c: *.inc $(MYLIBS)
	rm Obit_wrap.c
	echo "%module Obit" > Obit.i
	cat ObitTypeMaps.swig >> Obit.i
	cat *.inc >> Obit.i
	$(SWIG) -python  $(INCDIR) Obit.i

# Write compile/link info to setupdata.py
setupdata.py: Makefile
	echo "CFLAGS='$(ALL_CFLAGS)'">setupdata.py
	echo "CPPFLAGS='$(CLIENT_CPPFLAGS) $(SERVER_CPPFLAGS)'">>setupdata.py
	echo "LDFLAGS='$(CLIENT_LDFLAGS) $(SERVER_LDFLAGS)'">>setupdata.py
	echo "LIBS='$(CLIENT_LIBS) $(SERVER_LIBS)'">>setupdata.py

clean:
	rm -f Obit.i *.o *.so *.pyc	
	rm -rf build
