FC = gfortran
#FC = ifort
POWHEGPATH = /home/felix/POWHEG-BOX-RES
#POWHEGPATH = /home/karlbergl/cernbox/MCgenerators/POWHEG-BOX-RES
ifeq ("$(shell whoami)","felixr")
POWHEGPATH = /home/felixr/POWHEG-BOX-RES
endif
ifeq ("$(shell whoami)","tu_ptiie01")
POWHEGPATH = /home/tu/tu_tu/tu_ptiie01/Software/POWHEG-BOX-RES
endif
STD=-lstdc++
ifeq ("$(shell whoami)","zanderighi")
POWHEGPATH = /Users/zanderighi/tpimac/POWHEG/POWHEG-BOX-RES
STD=-lc++
endif
ifeq ("$(shell whoami)","zanderi")
POWHEGPATH = /home/th316/zanderi/POWHEG-BOX-RES
STD=-lc++
endif
ifeq ("$(shell whoami)","ferrario")
POWHEGPATH = /home/ferrario/work/POWHEG/RES
endif
ifeq ("$(shell whoami)","sferrario")
POWHEGPATH = $(HOME)/work/POWHEG/RES
STD=-lc++
endif
ifeq ("$(shell whoami)","ferrarioravasio")
POWHEGPATH = $(HOME)/POWHEG-BOX/RES
endif
ifeq ("$(shell whoami)","karlbergl")
POWHEGPATH = /home/karlbergl/cernbox/MCgenerators/POWHEG-BOX-RES
endif
ifeq ("$(shell whoami)","karlberg")
POWHEGPATH = /home/karlberg/cernbox/MCgenerators/POWHEG-BOX-RES
endif
ifeq ("$(shell whoami)", "siferrar")
POWHEGPATH = /afs/cern.ch/work/s/siferrar/private/POWHEG-BOX/RES
endif
ifeq ("$(shell whoami)", "krack")
POWHEGPATH=/Users/krack/POWHEG-BOX-RES
endif



INCLUDE1=$(PWD)
INCLUDE2=$(PWD)/..
INCLUDE3=$(POWHEGPATH)/include

ifeq ("$(FC)","gfortran")	
FFLAGS = -g -fno-automatic -ffixed-line-length-none -O2 -I$(INCLUDE1) -I$(INCLUDE2) -I$(INCLUDE3) 	
endif
ifeq ("$(FC)","ifort")	
FFLAGS = -save -extend_source -O3 -I$(INCLUDE1) -I$(INCLUDE2)  -I$(INCLUDE3) 
endif

# ----------------------------------------------------------------------------

# object files for qq->qqHH

OBJECTS = ./qq_ee.o \
          ./qqj_ee.o \
          ./help_hh.o \
	  ./braket_vbf.o  

TARGETS = libamplitude.a


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# rules for vbfnlo_lib:
#
 libamplitudes.a:  $(OBJECTS)
	$(AR) r libamplitude.a $(OBJECTS) 

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# for g77 or gfortran: 
.f.o:
	$(FC) -c -o $*.o  $(FFLAGS) $*.f

# ----------------------------------------------------------------------------

clean:
	rm -f *.o core *~ *.x libamplitude.a



