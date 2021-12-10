#include "alias.inc"
!*****************************************************************************80
! Licensing:
!    This code is distributed under the GNU LGPL license.
! PROGRAM VASP-BAUM Version 1.0 (f90)
! Written by Hyun-Jung Kim (h.kim@fz-juelich.de)
!  PGI-1, Forschungszentrum Juelich
! Copyright 2021. Hyun-Jung Kim All rights reserved.
! version 1.00
! last update and bug fixes : 2021. Dec. 10. by H.-J. Kim

! VASP BAnd Unfolding Machinery (VASP-BAUM)
program vaspbaum 
    use parameters
    use mpi_setup
    use time
    use version
    use print_io
    use mykind
    use wavecar
    implicit none
    real(kind=dp)                              t_start
    type(incar  )                           :: PINPT       ! parameters for input arguments
    type(eigen  )                           :: WAVEC       ! Eigenvalues from wavecar


    call parse_very_init(PINPT)

#ifdef MPI
    call mpi_initialize(PINPT%fnamelog)
#else
    call open_log(PINPT%fnamelog,myid)
#endif

    call version_stamp(t_start)
    call parse(PINPT)
    if(.not. PINPT%flag_set_unfold) then
        write(message,'(A,A)')"# File reading... : " ; write_msg
        write(message,'(A,A)')"  --> ", trim(PINPT%filenm) ; write_msg
    endif
    call inforead(PINPT, WAVEC)
    call write_info(WAVEC, PINPT, 0)

    call get_spectral_weight(WAVEC, PINPT)

    if(PINPT%icd .eq. 1) then
        call get_circular_dichroism_1(WAVEC, PINPT)
    elseif(PINPT%icd .eq. 2) then
        call get_circular_dichroism_matrix(WAVEC, PINPT)
    endif

#ifdef MPI
    call mpi_finish()
#endif
    call close_log(myid)

    stop
end program
