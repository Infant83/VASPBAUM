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
    real(kind=dp)       CD(3,PINPT%ispin,PINPT%nband,PINPT%nband,PINPT%nkpts)
    real(kind=dp)       SW(3,PINPT%ispin,PINPT%nediv,PINPT%nkpts)
    integer(kind=sp)    ii, idir, nediv
    real(kind=dp)       sigma, erange(PINPT%nediv), de

    nmin        = PINPT%ie_init
    nmax        = PINPT%ie_fina
    nelect      = PINPT%nelect
    sigma       = PINPT%sigma 
    nediv       = PINPT%nediv
    erange      = PINPT%init_e + eta + &
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

    CD          = 0d0
    SW          = 0d0
    WAVEC%CD    = 0d0
    WAVEC%SW    = 0d0

    call mpi_job_distribution_chain(PINPT%nkpts, nprocs, ourjob, ourjob_disp)
    do ik = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
        do ie = nmin, nelect
            do je = nelect+1, nmax
                call get_cd_ij(WAVEC%CD(:,:,ie,je,ik), WAVEC, PINPT, ie, je, ik)
               !if(ie .eq. 18 .and. je .eq. 19 .and. ik .eq. 1) then
               !    write(6,*)"BBBB ", WAVEC%CD(:,:,ie,je,ik)
               !   !stop
               !endif
                do is = 1, PINPT%ispin
                    de = WAVEC%E(je,is,ik) - WAVEC%E(ie,is,ik)
                    do idir = 1, 3 ! 1: RCD, 2:LCD, 3: CD
                        WAVEC%SW(idir,is,:,ik) = WAVEC%SW(idir,is,:,ik) + &
                                                 fgauss(sigma,erange-de)*WAVEC%CD(idir,is,ie,je,ik) / erange**2
                    enddo
                enddo
            enddo
        enddo
        
        write(message,'(A, F9.4,A)')'# IK, K(reci), done : ',real(ik)/real(ourjob(1))*100.d0, '%' ; write_msg
    enddo

!de = WAVEC%E(19,1,1) - WAVEC%E(18,1,1)
!write(6,*)"EEE ", fgauss(sigma,erange(660)-de), de, WAVEC%CD(:,1,18,19,1), fgauss(sigma,erange(660)-de)*WAVEC%CD(:,1,18,19,1)
!write(6,*)"WWWW", erange(660), erange(660)**2, erange**2
!stop

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    call MPI_ALLREDUCE(WAVEC%CD, CD, size(CD), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
    WAVEC%CD = CD
    call MPI_ALLREDUCE(WAVEC%SW, SW, size(SW), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
    WAVEC%SW = SW
#endif

    if_main call write_result_cd_spectral_weight(WAVEC, PINPT, erange)
    return
  endsubroutine

! subroutine get_cd_spectrum(SW, CD_ik, WAVEC, PINPT, ik)



  subroutine get_cd_ij(CD_ij, WAVEC, PINPT, ie, je, ik)
    use wavecar
    use mykind
    use parameters
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    integer(kind=sp)    ie, je, is, ik
    real(kind=dp)       CD_ij(3,PINPT%ispin)
   !complex(kind=dp)    Cjk(WAVEC%plw_idx_max,PINPT%ispinor)
   !complex(kind=dp)    vCik(WAVEC%plw_idx_max,PINPT%ispinor,3) ! -i*hbar*d/dx,y,z|Cik>
    complex(kind=dp)    Cjk(WAVEC%nplw_max/PINPT%ispinor,PINPT%ispinor)
    complex(kind=dp)    vCik(WAVEC%nplw_max/PINPT%ispinor,PINPT%ispinor,3) ! -i*hbar*d/dx,y,z|Cik>
    complex(kind=dp)    vex(3) ! velocity expectation value
    complex(kind=dp)    RCD, LCD ! right, left cd
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
        Cjk         = Cnk( WAVEC, PINPT, je, is, ik)
        vCik        = vCnk(WAVEC, PINPT, ie, is, ik)
        call get_vel_expectation(vex, ie, je, is, ik, WAVEC, PINPT)
        CD_ij(1,is) = abs(sum(A_R * vex))**2 ! RCD
        CD_ij(2,is) = abs(sum(A_L * vex))**2 ! LCD
        CD_ij(3,is) = (CD_ij(1,is)-CD_ij(2,is)) !/(CD_ij(1,is)+CD_ij(2,is))
    enddo

    return
  endsubroutine
