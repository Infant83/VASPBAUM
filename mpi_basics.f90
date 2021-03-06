module mpi_basics
   implicit none
#ifdef MPI
   include 'mpif.h'
#endif
   integer, public :: npar
   integer, public :: nproc_per_band
   integer, public :: mpi_comm_earth
   integer, public :: myid
   integer, public :: nprocs
   integer, public :: yourid
   integer, public :: earth_group

   type mpicomm 
        integer, public :: mpi_comm
        integer, public :: myid = 0
        integer, public :: nprocs = 1
        integer, public :: npar 
        integer, public :: comm_group, comm_group_supp

        ! for MPI_COMM_SPLIT purpose
        logical, public :: flag_split = .FALSE.
        integer, public :: key = 0
        integer, public :: color = 0
        integer, public, allocatable :: group_main(:)

        ! following values are generated by mpi_divide routine
        integer, public :: dims(2) ! grid dimension, (1) : nprow, (2) : npcol ! gener
        integer, public :: mycoord(2) ! default 2-dimensional cartesian topology, only meaningful if COMM_EARTH
   endtype

endmodule mpi_basics
