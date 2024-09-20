PROG = 	cmodel

#SUBDIRS = akt_mallit equilibrium
#OBJDIR1        = ./akt_mallit/
#OBJDIR2        = ./equilibrium/	

SRCS =	aerproc.f90 c2y.f90 chemset.f90 diffc.f90 diffcoef.f90 \
	driver.f90 driver_eq.f90 families.f90 coagulation.f90 \
	fileset.f90 headfile.f90 liq_ice2.f90\
	luo.f90 main.f90 output.f90 ph2oeq.f90 precisi.f90 radius.f90 \
	radiusice.f90 readfile.f90 readgas.f90 readliq.f90 readsld.f90 \
	setindex.f90 sizedist.f90 solids.f90 surface_tension.f90 \
	funs.f90 dvode_f90_m.f90 
# equilibrium/ylimaarainen.f90 \

OBJS =	aerproc.o c2y.o chemset.o diffc.o diffcoef.o driver.o \
	driver_eq.o families.o coagulation.o fileset.o headfile.o \
	liq_ice2.o luo.o main.o output.o ph2oeq.o precisi.o \
	radius.o radiusice.o readfile.o readgas.o readliq.o readsld.o \
	setindex.o sizedist.o solids.o surface_tension.o \
	funs.o dvode_f90_m.o
#ylimaarainen.o 
#SRCS =	aerproc.f90 aerproc2.f90 c2y.f90 chemset.f90 diffc.f90 diffcoef.f90 \
#	driver.f90 driver_eq.f90 equilset.f90 families.f90 fex.f90 \
#	fileset.f90 gases.f90 headfile.f90 jex.f90 liq_ice.f90 liquids.f90 \
#	luo.f90 main.f90 output.f90 ph2oeq.f90 precisi.f90 radius.f90 \
#	radiusice.f90 readfile.f90 readgas.f90 readliq.f90 readsld.f90 \
#	setindex.f90 sizedist.f90 solids.f90 solveq4.f90 surface_tension.f90 \
#	opkda1.f opkda2.f opkdmain.f thermo.f90#

#OBJS =	aerproc.o aerproc2.o c2y.o chemset.o diffc.o diffcoef.o driver.o \
#	driver_eq.o equilset.o families.o fex.o fileset.o gases.o headfile.o \
#	jex.o liq_ice.o liquids.o luo.o main.o output.o ph2oeq.o precisi.o \
#	radius.o radiusice.o readfile.o readgas.o readliq.o readsld.o \
#	setindex.o sizedist.o solids.o solveq4.o surface_tension.o opkda1.o \
#	opkda2.o opkdmain.o thermo.o


LIBS =#	/opt/intel/fc/9.0/lib/libguide.a

CC = cc
CFLAGS = -O
#FC = gfortran
#FC = g95
#FFLAGS = -O2 -Wall
#FC = f90
#FFLAGS = -O2 -cpu:host -N113 #-Rb
#F90 = f90
F90 = gfortran #ifort
F90FLAGS = -O1  -std=legacy -fbounds-check # -O2 -g -real_size 64 
F90FLAGS2 = -O1  -std=legacy -fbounds-check #-O2 -g -real_size 64  
#-parallel#-Wall
F90FLAGS_short = -O2 
#F90 = f90
#F90FLAGS = -O2 -cpu:host -N113 #-Rb
LDFLAGS = #-L/opt/intel/mkl/8.0.2/lib/32/ -lguide #-s

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod
#	rm -f /akt_mallit/ $(PROG) $(OBJS) *.mod
#	rm -f /equilibrium/ $(PROG) $(OBJS) *.mod

deepclean:
	rm -f *~
#	rm -f akt_mallit/*~
#	rm -f equilibrium/*~

.SUFFIXES: $(SUFFIXES) .f90 .f95 .o .mod 


.f90.o: .f90
	$(F90) $(F90FLAGS) -c $<


aerproc.o: aerproc.f90 headfile.o precisi.o 
	$(F90) $(F90FLAGS) -c $<


c2y.o: c2y.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

chemset.o: chemset.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

diffc.o: diffc.f90 precisi.o precisi.o
	$(F90) $(F90FLAGS) -c $<

diffcoef.o: diffcoef.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

driver.o: driver.f90 headfile.o precisi.o funs.o dvode_f90_m.o liq_ice2.o
	$(F90) $(F90FLAGS2) -c $<

driver_eq.o: driver_eq.f90 headfile.o precisi.o funs.o  dvode_f90_m.o
	$(F90) $(F90FLAGS) -c $<

dvode_f90_m.o: dvode_f90_m.f90 headfile.o
	$(F90) $(F90FLAGS) -c $<

families.o: families.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

coagulation.o: coagulation.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

#fex.o: headfile.o precisi.o
fileset.o: fileset.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

funs.o: funs.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

headfile.o: headfile.f90 precisi.o precisi.o
	$(F90) $(F90FLAGS) -c $<

liq_ice2.o: liq_ice2.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

luo.o: luo.f90 precisi.o precisi.o
	$(F90) $(F90FLAGS) -c $<

main.o: main.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

output.o: output.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

ph2oeq.o: ph2oeq.f90 precisi.o 
	$(F90) $(F90FLAGS) -c $<

precisi.o: precisi.f90
	$(F90) $(F90FLAGS) -c $<

radius.o: radius.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

radiusice.o: radiusice.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

readfile.o: readfile.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

readgas.o: readgas.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

readliq.o: readliq.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

readsld.o: readsld.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

setindex.o: setindex.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

sizedist.o: sizedist.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

solids.o: solids.f90 headfile.o precisi.o
	$(F90) $(F90FLAGS) -c $<

surface_tension.o: surface_tension.f90  precisi.o
	$(F90) $(F90FLAGS) -c $<


#ylimaarainen.o     : equilibrium/ylimaarainen.f90 headfile_equilibrium.o headfile_activity.o precisi.o
#	$(F90) $(F90FLAGS) -c $<


