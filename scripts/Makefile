#    A suite of routines to interpolate regional HYbrid Coordinate
#    Ocean Model (HYCOM) history (e.g., archv) direct-access binary
#    files to a user specified rectilinear grid covering the same
#    geographical region as the regional HYCOM experiment.  

#    Copyright (C) 2013 Henry R. Winterbottom

#    Email: Henry.Winterbottom@noaa.gov

#    Snail-mail:

#    Henry R. Winterbottom 
#    NOAA/OAR/PSD R/PSD1 
#    325 Broadway 
#    Boulder, CO 80303-3328

#    This program is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of
#    the License, or (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#    General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
#    02110-1301 USA.

# Define include files

include ./configure.hycomtograds

# Define modules and associated libraries

MODS = kinds.o                                       \
       constants.o                                   \
       mpi-interface.o                               \
       namelist.o                                    \
       kdtree2.o                                     \
       interpolation-interface.o                     \
       variable-interface.o                          \
       hycomtograds-interface.o

SUBS = time-utilities.o  

# Define executables

all: clean HYCOMTOGRADS

HYCOMTOGRADS:                      Driver.f90  $(MODS) $(SUBS) 
	$(LD) $(LDFLAGS)           Driver.f90  $(MODS) $(SUBS) -lm -o $@.exe

clean:
	-rm *.o *~ \#* *.mod *.bin fort.* *.exe >& /dev/null

# Define make rules

.f90.o:
	$(F90) $(CPPFLAGS) $(FCFFLAGS) $(OPTIMIZATION) $(DEBUG) -c $*.f90

.f.o:
	$(F77) $(CPPFLAGS) $(FCFFLAGS) $(OPTIMIZATION) $(DEBUG) -c $*.f

.F.o:
	$(F77) $(CPPFLAGS) $(FCFFLAGS) $(OPTIMIZATION) $(DEBUG) -c $*.F
