# $^ : dependency list
# $@ : target

#----- Change options and library path according to your system ------------#
#-----------------------------------
# Compiler options and bin path    |
#---------------------------------------------------------------------------|
#OPTIONS= -fpp -DMPI -mcmodel=large # for curion2
################# Possible options ##########################################
#  MPI_USE       : if "YES" MPI paralallism activation with the k-point parallization.
#                  IF "YES", with "make vaspbaum.mpi" tbfit.mpi will be compiled, and
#  Note: possible make command
#      make vaspbaum.mpi # generate tbfit execution file
#############################################################################
 TBBIN=$(HOME)/code/bin
 TBLIB=$(HOME)/code/lib
 VERSION="0.0.1"
 
#####################
# MAC-INTEL COMPILE #
#####################
#FC     = mpiifort
 FC     = mpif90
 OPTIONS= -fPIC -fpp 
 FFLAG  = -O2 -heap-arrays -nogen-interfaces -assume byterecl
 MPI_USE= YES

BIN    = $(TBBIN)
LIB	   = $(TBLIB)
#---------------------------------------------------------------------------|
MKLPATH   = $(MKLROOT)
LAPACK    = -L$(MKLPATH)/lib/ \
            -lmkl_intel_lp64 -lmkl_sequential \
            -lmkl_core -liomp5
INCLUDE   = -I$(MKLPATH)/include
#---------------------------------------------------------------------------|

######################### Do not modify below ###############################
#-----------------------------------
# Objects                          |
#---------------------------------------------------------------------------|

MPI_MOD= blacs_basics.o mpi_basics.o mpi_setup.o 
TEST   = test.o
MODULE = mykind.o parameters.o do_math.o print_io.o directory.o utils.o $(MPI_MOD) kill.o time.o version.o \
		 wavecar.o
READER = parse.o read_kpoints.o read_poscar.o
WRITER = write_info.o write_result_cd.o write_result_sw.o write_kpoints.o 
GET    = get_spectral_weight.o get_circular_dichroism.f90

LIBTOOL= ar src

OBJECTS =  $(MODULE) vaspbaum.o $(READER) $(WRITER) $(GET) $(TEST)

ifeq ($(MPI_USE), YES)
  F90  = $(FC) $(OPTIONS) -DMPI
  F90FLAGS = $(FFLAG)
else
  F90  = $(FC) $(OPTIONS)
  F90FLAGS = $(FFLAG)
endif

#---------------------------------------------------------------------------|

#-----------------------------------
# Suffix rules                     |
#-----------------------------------
.SUFFIXES: .f .f90
%.o: %.f90
	$(F90) $(FFLAG) -c $<

#-----------------------------------
# Targets                          |
#-----------------------------------

ifeq ($(MPI_USE), YES)
vaspbaum.mpi: $(OBJECTS)
	$(F90) -o $@ $^ $(LAPACK) $(INCLUDE)
	cp $@ $(BIN)/vaspbaum.mpi
else
vaspbaum.serial: $(OBJECT)
	$(F90) -o $@ $^ $(LAPACK) $(INCLUDE)
	cp $@ $(BIN)/vaspbaum.serial
endif

clean:
	rm *.o *.mod
