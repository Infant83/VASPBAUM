#include "alias.inc"

  ! compute circular dichroism
  ! calculate interband optical transition rate (Eq. 27) in PRB 92, 205108 (2015)
  subroutine get_circular_dichroism_1(WAVEC, PINPT)
    use mpi_setup
    use mykind
    use parameters
    use print_io
    use wavecar
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    integer(kind=sp)    is, ik
    integer(kind=sp)    ispinor
    integer(kind=sp)    ie, je     
    integer(kind=sp)    ourjob(nprocs), ourjob_disp(0:nprocs-1)
    real(kind=dp)       CD(3,PINPT%ispin,PINPT%nband,PINPT%nband,PINPT%nkpts)

    CD          = 0d0
    ie          = PINPT%ie_init
    je          = PINPT%ie_fina
    call mpi_job_distribution_chain(PINPT%nkpts, nprocs, ourjob, ourjob_disp)
    do ik = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
        call get_cd_ij(WAVEC%CD(:,:,ie,je,ik), WAVEC, PINPT, ie, je, ik)
       !write(message,'(A, I6, 3F12.6, F20.9)')'# IK, K(reci), CD(spin,ik) : ', &
       !      ik,WAVEC%kpts(:,ik), WAVEC%CD(3,:,ie,je,ik) ; write_msg
        write(message,'(A, F9.4,A)')'# IK, K(reci), done : ',real(ik)/real(ourjob(1))*100.d0, '%' ; write_msg       
    enddo
#ifdef MPI
    call MPI_ALLREDUCE(WAVEC%CD, CD, size(CD), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
    WAVEC%CD = CD
#endif
    if_main call write_result_cd_mode1(WAVEC, PINPT)
    return
  endsubroutine

  ! compute circular dichroism
  ! calculate interband optical transition rate (Eq. 27) in PRB 92, 205108 (2015)
  subroutine get_circular_dichroism_matrix(WAVEC, PINPT)
    use mpi_setup
    use mykind
    use parameters
    use print_io
    use wavecar
    use do_math
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    integer(kind=sp)    is, ik
    integer(kind=sp)    ispinor
    integer(kind=sp)    ie, je, nmin, nmax, nelect
    integer(kind=sp)    ourjob(nprocs), ourjob_disp(0:nprocs-1)
    real(kind=dp), allocatable   :: CD(:,:,:,:,:)
    real(kind=dp), allocatable   :: SW(:,:,:,:)
    integer(kind=sp)    ii, idir, nediv
    real(kind=dp)       sigma, erange(PINPT%nediv), de
    real(kind=dp)       very_small

    very_small  = 0.000000001d0
    nmin        = PINPT%ie_init
    nmax        = PINPT%ie_fina
    nelect      = PINPT%nelect
    sigma       = PINPT%sigma 
    nediv       = PINPT%nediv
    erange      = PINPT%init_e + very_small + &
                  dble((/(ii,ii=0,nediv-1)/))*(PINPT%fina_e-PINPT%init_e)/dble(nediv-1)
    
    if(nmin .gt. nelect) then
        write(message,'(A,2I6)')' # Error! : ie_init > nelect, (ie_init, nelect) = ',nmin, nelect
        write(message,'(A,2I6)')'     ==> set -ie ie_init correctly. Exit... '
        kill_job
    endif

    if(nmax .gt. PINPT%nband) nmax = PINPT%nband
    if(nmax .le. nelect ) then
        write(message,'(A,2I6)')' # Error! : ie_fina <= nelect, (ie_fina, nelect) = ',nmax, nelect
        write(message,'(A,2I6)')'     ==> set -if ie_fina correctly. Exit... '
        kill_job
    endif
    if(nmax .le. nmin ) then
        write(message,'(A,2I6)')' # Error! : ie_fina <= ie_init, (ie_init, ie_fina) = ',nmin, nmax
        write(message,'(A,2I6)')'     ==> set -ie ie_init -if ie_fina correctly. Exit... '
        kill_job
    endif

    allocate(CD(3,PINPT%nband,PINPT%nband,PINPT%ispin,PINPT%nkpts))
    allocate(SW(3,PINPT%nediv,PINPT%ispin,PINPT%nkpts))
    CD          = 0d0
    SW          = 0d0
    WAVEC%CD    = 0d0
    WAVEC%SW    = 0d0

    write(message,'(A)')' ' ; write_msgi
    write(message,'(A)')'#- START CIRCULAR DICHROISM CALCULATION: ' ; write_msgi

    call mpi_job_distribution_chain(PINPT%nkpts, nprocs, ourjob, ourjob_disp)
    do ik = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
        do ie = nmin, nelect
            do je = nelect+1, nmax
                call get_cd_ij(WAVEC%CD(:,ie,je,:,ik), WAVEC, PINPT, ie, je, ik)
                do is = 1, PINPT%ispin
                    de = WAVEC%E(je,is,ik) - WAVEC%E(ie,is,ik)
                    do idir = 1, 3 ! 1: RCD, 2:LCD, 3: CD
                        WAVEC%SW(idir,:,is,ik)=WAVEC%SW(idir,:,is,ik)+fgauss(sigma,erange-de)*WAVEC%CD(idir,ie,je,is,ik)/(erange**2)
                       !WAVEC%SW(idir,:,is,ik)=WAVEC%SW(idir,:,is,ik)+fgauss(sigma,erange-de)*WAVEC%CD(idir,ie,je,is,ik) !/(erange**2)
                    enddo
                enddo
            enddo
        enddo
        
        write(message,'(A, F9.4,A)')'# IK, K(reci), done : ',real(ik)/real(ourjob(1))*100.d0, '%' ; write_msg
    enddo

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    call MPI_ALLREDUCE(WAVEC%CD, CD, size(CD), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
    WAVEC%CD = CD
    call MPI_ALLREDUCE(WAVEC%SW, SW, size(SW), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
    WAVEC%SW = SW
#endif

    deallocate(CD)
    deallocate(SW)
    if_main call write_result_cd_spectral_weight(WAVEC, PINPT, erange)
    return
  endsubroutine

  ! compute circular dichroism
  ! calculate interband optical transition rate (Eq. 27) in PRB 92, 205108 (2015)
  subroutine get_circular_dichroism_unfold(WAVEC, PINPT)
    use mpi_setup
    use mykind
    use parameters
    use print_io
    use wavecar
    use do_math
    implicit none
    type(incar  )                    :: PINPT
    type(eigen  )                    :: WAVEC
    type(poscar )                    :: PGEOM_PC, PGEOM_SC
    integer(kind=sp)                    is, ik
    integer(kind=sp)                    ik_W
    integer(kind=sp)                    ispinor
    integer(kind=sp)                    ie, je, nmin, nmax, nelect
    integer(kind=sp)                    ourjob(nprocs), ourjob_disp(0:nprocs-1)
    integer(kind=sp)                    ii, idir, iaxis
    integer(kind=sp)                    nediv
    integer(kind=sp)                    nplw_screen
    complex(kind=dp)                    vex(3) ! velocity expectation value
    real(kind=dp)                       sigma, erange(PINPT%nediv), de
    real(kind=dp)                       theta_r, phi_r  ! angle of incident light in rad
    complex(kind=dp)                    A_R(3) ! operator for Right Chirality ! vector potential of incident circularly polarized light
    complex(kind=dp)                    A_L(3) ! operator for Left  Chirality ! vector potential of incident circularly polarized light
    integer(kind=sp)                    Cid(WAVEC%nplw_max)
    complex(kind=dp)                    CSjk(WAVEC%nplw_max/PINPT%ispinor,PINPT%ispinor)
    complex(kind=dp)                    vCSik(WAVEC%nplw_max/PINPT%ispinor,PINPT%ispinor,3) ! -i*hbar*d/dx,y,z|Cik>
    real(kind=dp),      allocatable  :: CD(:,:,:,:,:)
    real(kind=dp),      allocatable  :: SW(:,:,:,:)
    integer(kind=sp),   allocatable  :: ikpt(:)
    real(kind=dp)                       very_small

    very_small  = 0.000000001d0
    nmin        = PINPT%ie_init
    nmax        = PINPT%ie_fina
    nelect      = PINPT%nelect
    sigma       = PINPT%sigma
    nediv       = PINPT%nediv
    erange      = PINPT%init_e + very_small + &
                  dble((/(ii,ii=0,nediv-1)/))*(PINPT%fina_e-PINPT%init_e)/dble(nediv-1)
   !write(6,*)"ZZZ ", PINPT%init_e, PINPT%fina_e
   !kill_job

    ! incident Right(Left) circularly polarized light with angle (theta,phi)
    ! for details, see eq(9) of R. Yu et al., PRB 93, 205133 (2016)
    theta_r     = PINPT%theta * pi / 180.d0
    phi_r       = PINPT%phi   * pi / 180.d0
    A_R(1)      = (cos(theta_r) * cos(phi_r) + zi * sin(phi_r))
    A_R(2)      = (cos(theta_r) * sin(phi_r) - zi * cos(phi_r))
    A_R(3)      = -sin(theta_r)
    A_L(1)      = (cos(theta_r) * cos(phi_r) - zi * sin(phi_r))
    A_L(2)      = (cos(theta_r) * sin(phi_r) + zi * cos(phi_r))
    A_L(3)      = -sin(theta_r)

    if(nmin .gt. nelect) then
        write(message,'(A,2I6)')' # Error! : ie_init > nelect, (ie_init, nelect) = ',nmin, nelect
        write(message,'(A,2I6)')'     ==> set -ie ie_init correctly. Exit... '
        kill_job
    endif

    if(nmax .gt. PINPT%nband) nmax = PINPT%nband
    if(nmax .le. nelect ) then
        write(message,'(A,2I6)')' # Error! : ie_fina <= nelect, (ie_fina, nelect) = ',nmax, nelect
        write(message,'(A,2I6)')'     ==> set -if ie_fina correctly. Exit... '
        kill_job
    endif
    if(nmax .le. nmin ) then
        write(message,'(A,2I6)')' # Error! : ie_fina <= ie_init, (ie_init, ie_fina) = ',nmin, nmax
        write(message,'(A,2I6)')'     ==> set -ie ie_init -if ie_fina correctly. Exit... '
        kill_job
    endif

    ! obtain PC k-path to be unfoled
    call prepare_unfold_kpoint(PINPT, PGEOM_PC, PGEOM_SC, 0, .false.)

    ! allocate and initialize CD&SW array and its k-point index
    allocate(CD(3,PINPT%nband,PINPT%nband,PINPT%ispin,PGEOM_PC%nkpts))
    allocate(SW(3,PINPT%nediv,PINPT%ispin,PGEOM_PC%nkpts))
    allocate(ikpt(PGEOM_PC%nkpts))
    CD          = 0d0
    SW          = 0d0
    ikpt        = 0
    if(allocated(WAVEC%CD))         deallocate(WAVEC%CD)
    if(allocated(WAVEC%SW))         deallocate(WAVEC%SW)
    if(allocated(PGEOM_PC%ikpt))    deallocate(PGEOM_PC%ikpt)
    allocate(WAVEC%CD(3,PINPT%ispin,PINPT%nband,PINPT%nband,PGEOM_PC%nkpts))
    allocate(WAVEC%SW(3,PINPT%ispin,PINPT%nediv,PGEOM_PC%nkpts))
    allocate(PGEOM_PC%ikpt(PGEOM_PC%nkpts))
    WAVEC%CD        = 0d0
    WAVEC%SW        = 0d0
    PGEOM_PC%ikpt   = 0

    ! start main
    write(message,'(A)')' ' ; write_msgi
    write(message,'(A)')'#- START CIRCULAR DICHROISM WITH UNFOLDING: ' ; write_msgi
#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif

    call mpi_job_distribution_chain(PGEOM_PC%nkpts, nprocs, ourjob, ourjob_disp)
 kp:do ik = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))

        call K_index_von_WAVECAR(ik_W, WAVEC, PGEOM_SC%kpts(:,ik))
        write(message,'(A,I4,A,3F9.4,A,3F9.4,A,3I4,A,I5)')'  --> Processing ',ik,'-th k-point k0: ', PGEOM_PC%kpts(:,ik), &
              ' K0: ',WAVEC%kpts(:,ik_W),' G0: ',PGEOM_SC%G(:,ik),' ikp= ',ik_W ;call write_log(message,12,myid)

        call screen_PC_G_von_SC_G(ik_W, WAVEC, PGEOM_SC%M, PGEOM_SC%G(:,ik), PINPT%ispinor)
        call get_Gid(WAVEC%GSid(:,ik_W), WAVEC%GS(:,:,ik_W), WAVEC%nplw_screen(ik_W), WAVEC%nplw_max, WAVEC%nbmax)

        ikpt(ik)    = ik_W
        nplw_screen = WAVEC%nplw_screen(ik_W)
        Cid         = CSid(WAVEC, PINPT, ik_W)
!if(ik .ne. 56) cycle
    isp:do is = 1, PINPT%ispin
        do ie = nmin, nelect
            vCSik = vCSnk(Cid, WAVEC, PINPT, ie, is, ik_W)
            do je = nelect+1, nmax
               CSjk = CSnk(Cid, WAVEC, PINPT, je, is, ik_W)

                ! <psi_j|A*v|psi_i>
                vex = 0d0
                do ispinor = 1, PINPT%ispinor
                    do iaxis = 1, 3
                        vex(iaxis) = vex(iaxis) + dot_product(CSjk(1:nplw_screen,ispinor),vCSik(1:nplw_screen,ispinor,iaxis))
                    enddo
                enddo
                
                CD(1,ie,je,is,ik) = abs(sum(A_R * vex))**2 ! RCD
                CD(2,ie,je,is,ik) = abs(sum(A_L * vex))**2 ! LCD
                CD(3,ie,je,is,ik) = (CD(1,ie,je,is,ik)-CD(2,ie,je,is,ik)) ! CD = RCD - LCD

                de = WAVEC%E(je,is,ik_W) - WAVEC%E(ie,is,ik_W)
                do idir = 1, 3 ! 1: RCD, 2:LCD, 3: CD
                   !SW(idir,:,is,ik)=SW(idir,:,is,ik)+fgauss(sigma,erange-de)*CD(idir,ie,je,is,ik)/erange**2
                    SW(idir,:,is,ik)=SW(idir,:,is,ik)+fgauss(sigma,erange-de)*CD(idir,ie,je,is,ik) !/erange**2
                enddo
!               if(ik .eq. 1  ) then
!                   write(6,'(A,3I5,7F9.4)')"ZZZZ",ik, ie, je, WAVEC%E(je,is,ik), WAVEC%E(ie,is,ik), de, erange(1:2), fgauss(sigma,erange(1:2)-de) !, SW(:,1:10,1,ik)
!                   kill_job1
!               endif
            enddo
        enddo
        enddo isp
    enddo kp

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    call MPI_ALLREDUCE(CD, WAVEC%CD, size(CD), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
    call MPI_ALLREDUCE(SW, WAVEC%SW, size(SW), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
    call MPI_ALLREDUCE(ikpt, PGEOM_PC%ikpt, size(ikpt), MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
#else
    WAVEC%CD = CD
    WAVEC%SW = SW
    PGEOM_PC%ikpt = ikpt
#endif

    call write_cd_spectral_function(WAVEC, PINPT, PGEOM_PC, erange)
    kill_job
    deallocate(CD)
    deallocate(SW)
    return
  endsubroutine


  subroutine get_cd_ij(CD_ij, WAVEC, PINPT, ie, je, ik)
    use wavecar
    use mykind
    use parameters
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    integer(kind=sp)    ie, je, is, ik
    real(kind=dp)       CD_ij(3,PINPT%ispin)
    complex(kind=dp)    vex(3) ! velocity expectation value
    real(kind=dp)       theta_r, phi_r  ! angle of incident light in rad
    complex(kind=dp)    A_R(3)          ! operator for Right Chirality
    complex(kind=dp)    A_L(3)          ! operator for Left  Chirality ! vector potential of incident circularly polarized light

    CD_ij = 0d0
    ! incident Right(Left) circularly polarized light with angle (theta,phi)
    ! for details, see eq(9) of R. Yu et al., PRB 93, 205133 (2016)
    theta_r     = PINPT%theta * pi / 180.d0
    phi_r       = PINPT%phi   * pi / 180.d0
    A_R(1)      = (cos(theta_r) * cos(phi_r) + zi * sin(phi_r))
    A_R(2)      = (cos(theta_r) * sin(phi_r) - zi * cos(phi_r))
    A_R(3)      = -sin(theta_r)
    A_L(1)      = (cos(theta_r) * cos(phi_r) - zi * sin(phi_r))
    A_L(2)      = (cos(theta_r) * sin(phi_r) + zi * cos(phi_r))
    A_L(3)      = -sin(theta_r)
 
    do is = 1, PINPT%ispin
        call get_vel_expectation(vex, ie, je, is, ik, WAVEC, PINPT)
        CD_ij(1,is) = abs(sum(A_R * vex))**2 ! RCD
        CD_ij(2,is) = abs(sum(A_L * vex))**2 ! LCD
        CD_ij(3,is) = (CD_ij(1,is)-CD_ij(2,is))
    enddo

    return
  endsubroutine
