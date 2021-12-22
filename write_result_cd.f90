#include "alias.inc"

  subroutine write_result_cd_mode1(WAVEC, PINPT)
    use parameters
    use utils
    use mykind
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    character(len=256)  fname
    integer(kind=sp)    ie, je, ik, is
    integer(kind=sp)    pid
    character(len=20)   ie_str,je_str,ik_str,is_str
    integer(kind=sp)    max_ik(1), min_ik(1)
    real(kind=dp)       max_cd   , min_cd

    pid     = pid_output

    ie      = PINPT%ie_init
    je      = PINPT%ie_fina

    ie_str  = int2str(ie)
    je_str  = int2str(je)

    do is = 1, PINPT%ispin
        is_str  = int2str(is)
        write(fname,*)trim(PINPT%folder_out)//'OPT_TRANS_RATE_IE',(trim(adjustl(ie_str))),'_JE',trim(adjustl(je_str)),'_SP',trim(adjustl(is_str)),'.dat'
        open(pid, file=trim(fname), status='unknown')
        call write_info(WAVEC,PINPT, pid)
        write(pid,'(A,I4,A,I4)')"# OPTICAL SELECTIVITY (CD,n) BETWEEN BANDS: ", ie,"    -  ",je
        write(pid,'(A)')"# n(k,w_cv)= |P(k,s,cv,+)|^2 - |P(k,s,cv,-)|^2"
        write(pid,'(A)')"#  ||        ---------------------------------"
        write(pid,'(A)')"#  ||        |P(k,s,cv,+)|^2 + |P(k,s,cv,-)|^2"
        write(pid,'(A)')"#  ||             RCD        -       LCD      "
        write(pid,'(A)')"#  ===> CD = ---------------------------------"
        write(pid,'(A)')"#                 RCD        +       LCD      "
        write(pid,'(A)')"#  The TRANSITION MATRIX ELEMENT P ="
        write(pid,'(A)')"#   P(k,s,cv,+ or -) = 1/sqrt(2)[p_x(k,cv,s) + (or -) i*p_y(k,cv,s)]"
        write(pid,'(A)')"#  THE INTERBAND TRANSITION MATRIX p_x,y ="
        write(pid,'(A)')"#   p_x,y(k,cv,s)=<psi(k,c,s)|-i*hbar*1/dx(y)|psi(k,v,s>"
        
        max_ik     = maxloc(WAVEC%CD(3,ie,je,is,:))
        min_ik     = minloc(WAVEC%CD(3,ie,je,is,:))
        max_cd     =        WAVEC%CD(3,ie,je,is,max_ik(1)) 
        min_cd     =        WAVEC%CD(3,ie,je,is,min_ik(1)) 
        write(pid,'(A,I0,A,4F11.6)')"# MAXVAL(SPIN=",is,") of SELECTIVITY at kx,ky,kz (in reci)= ",WAVEC%kpts(:,max_ik), max_cd
        write(pid,'(A,I0,A,4F11.6)')"# MINVAL(SPIN=",is,") of SELECTIVITY at kx,ky,kz (in reci)= ",WAVEC%kpts(:,min_ik), min_cd
        write(pid,'(A)')"# (cart) kx        ky        kz(A^-1)   (recip)kx        ky        kz             RCD             LCD              CD"

        do ik = 1, PINPT%nkpts
            write(pid,'(3F11.6,5X,3F11.6, 3(F20.9) )'), WAVEC%kpts(:,ik), WAVEC%kpts_cart(:,ik), WAVEC%CD(:,ie,je,is,ik)
        enddo

        close(pid)
    enddo

    return
  endsubroutine

  subroutine write_result_cd_spectral_weight(WAVEC, PINPT, erange)
    use parameters
    use utils
    use mykind
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    character(len=256)  fname
    character(len=20)   ik_str, is_str
    integer(kind=sp)    ie, je, ik, is
    integer(kind=sp)    pid
    real(kind=dp)       kline(PINPT%nkpts), erange(PINPT%nediv)

    pid     = pid_output

    call get_kline_dist(WAVEC%kpts_cart, PINPT%nkpts, kline)

    do is = 1, PINPT%ispin
        is_str  = int2str(is)
        write(fname,*)trim(PINPT%folder_out)//'OPT_TRANS_SP',trim(adjustl(is_str)),'.dat'
    
        open(pid, file=trim(fname), status = 'unknown')
        call write_info(WAVEC, PINPT, pid)

        write(pid,'(A,I0,A,2F16.5)')"# MAX and MIN CD (SPIN-",is,"): ", maxval(WAVEC%SW(3,:,is,:)), minval(WAVEC%SW(3,:,is,:))
        do ik = 1, PINPT%nkpts
            write(pid,'(A)')"# KPATH(A-1)           ENERGY(eV)         SW-RCD              SW-LCD                SW-CD"
            do ie = 1, PINPT%nediv
                write(pid,'(F11.6, F20.6, 3F20.9)')kline(ik), erange(ie), WAVEC%SW(:,ie,is,ik)
            enddo
            write(pid,'(A)')' '
            write(pid,'(A)')' '
        enddo

        close(pid)
    enddo

    return
  endsubroutine

  subroutine write_cd_spectral_function(WAVEC, PINPT, PGEOM, erange)
    use parameters
    use utils
    use mpi_setup
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    type(poscar )    :: PGEOM
    character(len=256)  fname
    integer(kind=sp)    ii
    integer(kind=sp)    ie, je, ik, is
    integer(kind=sp)    pid
    character(len=20)   ie_str, je_str, ik_str, is_str
    integer(kind=sp)    max_ik(1), min_ik(1)
    real(kind=dp)       max_sw, min_sw
    real(kind=dp)       kline(PGEOM%nkpts), erange(PINPT%nediv)

    if(myid .ne. 0) return

    call get_kline_dist(PGEOM%kpts_cart, PGEOM%nkpts, kline)

    do is = 1, PINPT%ispin
        is_str = int2str(is-1)
        write(fname,*)trim(PINPT%folder_out)//'UNFOLD_CD_SP',trim(adjustl(is_str)),'.dat'
        open(pid, file=trim(fname), status='unknown')

        write(pid,'(A,I0,A,2F16.5)')"# MAX and MIN CD (SPIN-",is,"): ", maxval(WAVEC%SW(3,:,is,:)), minval(WAVEC%SW(3,:,is,:))
        do ik = 1, PGEOM%nline+1
            write(pid,'(A,I6,A, F10.5)')'# K SEGMENTS ',ik, ' : ', kline(PGEOM%k_name_index(ik))
        enddo

        do ik = 1, PGEOM%nkpts
            write(pid,'(A)')"# KPATH(A-1)           ENERGY(eV)         SW-RCD              SW-LCD                SW-CD"
            write(pid,'(A)')'# KPATH   ENERGY   SW '
            do ie = 1, PINPT%nediv
                write(pid,'(F11.6,F20.6,3F20.9)') kline(ik), erange(ie) , WAVEC%SW(:, ie, is, ik)
            enddo
            write(pid,'(A)')' '
            write(pid,'(A)')' '
        enddo

        close(pid)
    enddo

    return
  endsubroutine
