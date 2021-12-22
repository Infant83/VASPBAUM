#include "alias.inc"

  subroutine get_band_unfold(WAVEC, PINPT)
    use mpi_setup
    use mykind
    use parameters
    use print_io
    use utils
    use do_math
    implicit none
    type(incar  )                    :: PINPT
    type(eigen  )                    :: WAVEC
    type(poscar )                    :: PGEOM_PC, PGEOM_SC

    if((.not. PINPT%flag_set_unfold .and. .not. PINPT%flag_unfold)) return
    if(PINPT%flag_set_unfold) PINPT%flag_unfold = .FALSE.

    if(PINPT%flag_set_unfold) then
        call prepare_unfold_kpoint(PINPT, PGEOM_PC, PGEOM_SC, 1, PINPT%flag_reduce)
        call write_unfold_kpoint(PINPT, PGEOM_SC)
        return
    endif

    ! calculate spectral weight P_Km(k)
    call get_spectral_weight_band(WAVEC, PINPT, PGEOM_PC)

    ! calculate spectral function A(k,E)
    call get_spectral_function(WAVEC, PINPT, PGEOM_PC)

    ! write result
    call write_spectral_function(WAVEC, PINPT, PGEOM_PC)

    return
  endsubroutine

  ! compute spectral function
  subroutine get_spectral_function(WAVEC, PINPT, PGEOM_PC)
    use mpi_setup
    use mykind
    use parameters
    use print_io
    use utils
    use do_math
    implicit none
    type(incar  )                    :: PINPT
    type(eigen  )                    :: WAVEC
    type(poscar )                    :: PGEOM_PC
    integer(kind=sp)                    ii
    integer(kind=sp)                    nmin, nmax, nediv
    real(kind=dp)                       sigma, erange(PINPT%nediv), de
    integer(kind=sp)                    ik, ie, je, is
    integer(kind=sp)                    ourjob(nprocs), ourjob_disp(0:nprocs-1)
    real(kind=dp), allocatable       :: SW(:,:,:)

    ! NOTE:
    ! spectral function SW -> A(E,k) = sum_m P_Km(k) * delta(Em - ef - E)
    ! ef : fermi level. It is already substituted when reading Em (see wavecar.f90:Enk_read)

    nmin        = PINPT%ie_init
    nmax        = PINPT%ie_fina ; if(nmax .gt. PINPT%nband) nmax = PINPT%nband
    sigma       = PINPT%sigma
    nediv       = PINPT%nediv
    erange      = PINPT%init_e + eta + &
                  dble((/(ii,ii=0,nediv-1)/))*(PINPT%fina_e-PINPT%init_e)/dble(nediv-1)

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif

    if(allocated(WAVEC%SW)) deallocate(WAVEC%SW)
       allocate(WAVEC%SW(1,nediv,PINPT%ispin,PGEOM_PC%nkpts)) ; WAVEC%SW = 0d0
       allocate(      SW(  nediv,PINPT%ispin,PGEOM_PC%nkpts)) ;       SW = 0d0

    call mpi_job_distribution_chain(PGEOM_PC%nkpts, nprocs, ourjob, ourjob_disp)
    do ik = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
        do ie = nmin, nmax
            do is = 1, PINPT%ispin
                SW(:,is,ik) = SW(:,is,ik) + &
                              florentz(sigma,WAVEC%E(ie,is,PGEOM_PC%ikpt(ik))-erange)*WAVEC%SW_BAND(ie,is,ik)
            enddo
        enddo
    enddo

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    call MPI_ALLREDUCE(SW, WAVEC%SW(1,:,:,:), size(SW), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
#else
    WAVEC%SW(1,:,:,:) = SW
#endif

    return
  endsubroutine

  ! compute spectral weight in the primitive BZ for each band.
  subroutine get_spectral_weight_band(WAVEC, PINPT, PGEOM_PC)
    use mpi_setup
    use mykind
    use parameters
    use print_io
    use wavecar, only: get_Gid, Cnk, CSnk, CSid
    use utils
    implicit none
    type(incar  )                    :: PINPT
    type(eigen  )                    :: WAVEC
    type(poscar )                    :: PGEOM_PC, PGEOM_SC
    integer(kind=sp)                    ik, ik_W, ispinor
    integer(kind=sp)                    is, ie
    integer(kind=sp), allocatable    :: ikpt(:)
    integer(kind=sp)                    Cid(WAVEC%nplw_max)
    integer(kind=sp)                    nplw_screen
    complex(kind=dp)                    CSik(WAVEC%nplw_max/PINPT%ispinor,PINPT%ispinor)
    complex(kind=dp)                    Cik(WAVEC%nplw_max/PINPT%ispinor,PINPT%ispinor)
    real(kind=dp),    allocatable    :: SW_BAND(:,:,:) ! spectral weight
    integer(kind=sp)                    nmin, nmax
    integer(kind=sp)                    ourjob(nprocs), ourjob_disp(0:nprocs-1)

    !   Get spectral weight for k0 [see eq(8) of PRB 85, 085201 (2012)]
    !
    !       P_{Km}(k0) = \sum_n |<Km | kn>|^2
    !
    !      --> adopting a plane-wave expansion, [ see eq(15) of PRB 85, 085201 (2012)]
    !
    !       P_{Km}(k0) = \sum_{g} |C_{Km}(g + k0 - K)|^2
    !
    !      --> k - K = G (shift G vector = PGEOM_SC%G)
    !      --> Summation runs over g vectors of primitive cell BZ.
    !          Therefore, one should find g vectors from the supercell G vectors.
    !          In "screen_PC_G_von_SC_G" routine, it finds g vectors from supercell G vectors
    !          and adding shifting G vectors (PGEOM_SC%G)
    !      
    !       g + k0 - K = g + PGEOM_SC%G => GS
    !       
    !       Set unique index for GS vectors using get_Gid function

    !   P_{Km}(k0) represents the probability of finding a set of PC states |k0,n>
    !   contributing to the SC state |K,m>, or, equivalently,
    !   the amount of Bloch character k0 preserved in |K,m> at the same energy En = Em.

    nmin        = PINPT%ie_init
    nmax        = PINPT%ie_fina ; if(nmax .gt. PINPT%nband) nmax = PINPT%nband

    call prepare_unfold_kpoint(PINPT, PGEOM_PC, PGEOM_SC, 0, .false.)
    write(message,'(A,I5,A)')'  --> find ',PGEOM_SC%nkpts, ' k-points along K-path' ; write_msgi
    allocate(SW_BAND(PINPT%nband, PINPT%ispin, PGEOM_PC%nkpts)) ; SW_BAND = 0d0
    if( allocated(WAVEC%SW_BAND)) deallocate(WAVEC%SW_BAND)
        allocate(WAVEC%SW_BAND(PINPT%nband, PINPT%ispin, PGEOM_PC%nkpts)) ; WAVEC%SW_BAND = 0d0
    allocate(ikpt(PGEOM_PC%nkpts))  ; ikpt = 0
    if( allocated(PGEOM_PC%ikpt)) deallocate(PGEOM_PC%ikpt) 
        allocate(PGEOM_PC%ikpt(PGEOM_PC%nkpts)) ; PGEOM_PC%ikpt = 0

    write(message,'(A)')' ' ; write_msgi
    write(message,'(A)')'#- START BAND UNFOLDING: ' ; write_msgi

    call mpi_job_distribution_chain(PGEOM_PC%nkpts, nprocs, ourjob, ourjob_disp)
    do ik = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))

        call K_index_von_WAVECAR(ik_W, WAVEC, PGEOM_SC%kpts(:,ik))
        write(message,'(A,I4,A,3F9.4,A,3F9.4,A,3I4,A,I5)')'  --> Processing ',ik,'-th k-point k0: ', PGEOM_PC%kpts(:,ik), &
              ' K0: ',WAVEC%kpts(:,ik_W),' G0: ',PGEOM_SC%G(:,ik),' ikp= ',ik_W ;call write_log(message,12,myid)

        call screen_PC_G_von_SC_G(ik_W, WAVEC, PGEOM_SC%M, PGEOM_SC%G(:,ik), PINPT%ispinor)
        call get_Gid(WAVEC%GSid(:,ik_W), WAVEC%GS(:,:,ik_W), WAVEC%nplw_screen(ik_W), WAVEC%nplw_max, WAVEC%nbmax)
    
        ikpt(ik)    = ik_W
        nplw_screen = WAVEC%nplw_screen(ik_W)
        Cid         = CSid(WAVEC, PINPT, ik_W)

        do is = 1, PINPT%ispin
            do ie = nmin, nmax
                CSik = CSnk(Cid, WAVEC, PINPT, ie, is, ik_W) ! obtain screened Coeff 
                do ispinor = 1, PINPT%ispinor 
                    SW_BAND(ie,is,ik)=SW_BAND(ie,is,ik)+dble(dot_product(CSik(1:nplw_screen,ispinor),CSik(1:nplw_screen,ispinor)))
                enddo
            enddo
        enddo
    enddo

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    call MPI_ALLREDUCE(SW_BAND, WAVEC%SW_BAND, size(SW_BAND), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
    call MPI_ALLREDUCE(ikpt, PGEOM_PC%ikpt, size(ikpt), MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
#else
    WAVEC%SW_BAND   = SW_BAND
    PGEOM_PC%ikpt   = ikpt
#endif

    call write_result_spectral_weight(WAVEC, PINPT, PGEOM_PC)

    return
  endsubroutine

! function integ_exp(init_z, fina_z, atten_z) result(f)
!   use mykind
!   implicit none
!   real(kind=dp)       init_z, fina_z, atten_z
!   complex(kind=dp)    f

!   


!   return
! endfunction
