#include "alias.inc"

  subroutine write_result_spectral_weight(WAVEC, PINPT, PGEOM)
    use parameters
    use utils
    use mpi_setup
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    type(poscar )    :: PGEOM
    character(len=256)  fname
    integer(kind=sp)    ie, je, ik, is
    integer(kind=sp)    pid
    character(len=20)   ie_str, je_str, ik_str, is_str
    integer(kind=sp)    max_ik(1), min_ik(1)
    real(kind=dp)       max_sw, min_sw
    real(kind=dp)       kline(PGEOM%nkpts)
    real(kind=dp)       SW(PINPT%nband, PINPT%ispin, PGEOM%nkpts)
    real(kind=dp)       ef

    if(myid .ne. 0) return

    pid     = pid_output

    ie      = PINPT%ie_init
    je      = PINPT%ie_fina
    
    ie_str  = int2str(ie)
    je_str  = int2str(je)
    kline   = 0d0
    ef      = PINPT%e_fermi

    call get_kline_dist(PGEOM%kpts_cart, PGEOM%nkpts, kline)

    do is = 1, PINPT%ispin
        is_str  = int2str(is-1)
        write(fname,*)trim(PINPT%folder_out)//'sw_spin',trim(adjustl(is_str)),'.dat'
        open(pid, file=trim(fname), status='unknown')
        
        do ie = 1, PINPT%nband
            write(pid,'(A,I0,A)')'# KPATH(A^-1)   ENERGY(eV)      SW :  ', ie, ' -th eigenvalue'
            do ik = 1, PGEOM%nkpts
                write(pid,'(F11.6,2F15.6)')kline(ik),WAVEC%E(ie,is,PGEOM%ikpt(ik))-ef,WAVEC%SW_BAND(ie,is,ik)
            enddo
            write(pid,'(A)')' '
            write(pid,'(A)')' '
        enddo

        close(pid)
    enddo

    return
  endsubroutine

  subroutine write_result_spectral_function(WAVEC, PINPT, PGEOM, erange)
    use parameters
    use utils
    use mpi_setup
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    type(poscar )    :: PGEOM
    character(len=256)  fname
    integer(kind=sp)    ie, je, ik, is
    integer(kind=sp)    pid
    character(len=20)   ie_str, je_str, ik_str, is_str
    integer(kind=sp)    max_ik(1), min_ik(1)
    real(kind=dp)       max_sw, min_sw
    real(kind=dp)       kline(PGEOM%nkpts), erange(PINPT%nediv)
    real(kind=dp)       SW(PINPT%nband, PINPT%ispin, PGEOM%nkpts)
    real(kind=dp)       ef

    if(myid .ne. 0) return

    ef = PINPT%e_fermi

    call get_kline_dist(PGEOM%kpts_cart, PGEOM%nkpts, kline)

    do is = 1, PINPT%ispin
        is_str = int2str(is-1)
        write(fname,*)trim(PINPT%folder_out)//'SW.SPIN',trim(adjustl(is_str)),'.dat'
        open(pid, file=trim(fname), status='unknown')
        
        do ik = 1, PGEOM%nline+1
            write(pid,'(A,I6,A, F10.5)')'# K SEGMENTS ',ik, ' : ', kline(PGEOM%k_name_index(ik))
        enddo

        do ik = 1, PGEOM%nkpts
            write(pid,'(A)')'# KPATH   ENERGY   SW '
            do ie = 1, PINPT%nediv
                write(pid,'(3F16.8)') kline(ik), erange(ie)-ef, WAVEC%SW(1, ie, is, ik)
            enddo
            write(pid,'(A)')' '
            write(pid,'(A)')' '
        enddo

        close(pid)
    enddo

    return
  endsubroutine
