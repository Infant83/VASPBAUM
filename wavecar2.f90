#include "alias.inc"

module wavecar
    use mpi_setup
    use parameters
    use mykind
    use print_io

contains

  subroutine inforead(PINPT, WAVEC)
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    integer(kind=sp)    ik, is, ie
    integer(kind=sp)    ierr
    integer(kind=sp)    nelect_prev
    real(kind=dp)       xirecl, xispin, xiprec
    real(kind=dp)       xnkpts, xnband, xnelect

    PINPT%irecl = 24
    xnelect     = 0.0d0
    nelect_prev = 0

    open(unit=pid_wavecar, file=trim(PINPT%filenm), access='direct', recl=PINPT%irecl, iostat=ierr, status='old')
    if(ierr .ne. 0) then 
        write(message,'(A,I0,A,A)') ' FILE OPEN ERROR (0) : IERR = ', ierr, ' FILE= ', trim(PINPT%filenm) ; write_msg
    endif

    read(pid_wavecar, rec=1) xirecl, xispin, xiprec
    close(pid_wavecar)

    PINPT%irecl=nint(xirecl)
    PINPT%iprec=nint(xiprec) ! set to integer
    PINPT%ispin=nint(xispin)
    if(PINPT%ispin .eq. 2) PINPT%ispinor = 1

    if(PINPT%iprec .eq. 45210) then
      write(message,'(A)') ' ' ; write_msg
      write(message,'(A)') '*** error - WAVECAR_double requires complex*16 for coeff' ; write_msg
      kill_job
    endif
    
    open(unit=pid_wavecar, file=trim(PINPT%filenm), access='direct', recl=PINPT%irecl, iostat=ierr, status='old')
    if(ierr .ne. 0) then
        write(message,'(A,I0,A,A)') ' FILE OPEN ERROR (1) : IERR = ', ierr, ' FILE= ', trim(PINPT%filenm) ; write_msg
    endif

    read(pid_wavecar, rec=2) xnkpts, xnband, PINPT%encut, PINPT%a1(1:3), PINPT%a2(1:3), PINPT%a3(1:3) 
    PINPT%A(:,1)= PINPT%a1
    PINPT%A(:,2)= PINPT%a2
    PINPT%A(:,3)= PINPT%a3

    PINPT%nkpts = nint(xnkpts)
    PINPT%nband = nint(xnband)

    if(allocated(WAVEC%kpts_cart)) deallocate(WAVEC%kpts_cart) ; allocate(WAVEC%kpts_cart(          3,            PINPT%nkpts))
    if(allocated(WAVEC%kpts))      deallocate(WAVEC%kpts)      ; allocate(WAVEC%kpts     (          3,            PINPT%nkpts))
    if(allocated(WAVEC%E   ))      deallocate(WAVEC%E   )      ; allocate(WAVEC%E        (PINPT%nband,PINPT%ispin,PINPT%nkpts))
    if(allocated(WAVEC%OCC ))      deallocate(WAVEC%OCC )      ; allocate(WAVEC%OCC      (PINPT%nband,PINPT%ispin,PINPT%nkpts))
    if(allocated(WAVEC%nplw))      deallocate(WAVEC%nplw)      ; allocate(WAVEC%nplw     (                        PINPT%nkpts))

    if(PINPT%icd .ge. 1) then
        if(allocated(WAVEC%CD))    deallocate(WAVEC%CD  )      
        allocate(WAVEC%CD(3,PINPT%ispin,PINPT%nband,PINPT%nband,PINPT%nkpts))
        WAVEC%CD = 0d0
        if(PINPT%icd .ge. 2) then
            if(allocated(WAVEC%SW))deallocate(WAVEC%SW  )
            allocate(WAVEC%SW(3,PINPT%ispin,PINPT%nediv,PINPT%nkpts))
            WAVEC%SW = 0d0
        endif
    endif

    ! read eigenvalue, kpoints, occupations, number of plane waves
    do is = 1, PINPT%ispin
        do ik = 1, PINPT%nkpts
            call Enk_read(WAVEC, PINPT, is, ik)
        enddo
    enddo

    ! estimate number of occupied electrons
    if(PINPT%nelect .eq. 0) then
        do is = 1, PINPT%ispin
            nelect_prev = 0
            do ik = 1, PINPT%nkpts
                xnelect = sum(WAVEC%OCC(1:PINPT%nband, is, ik))
                if(ik .eq. 1) then 
                    nelect_prev = nint(xnelect)
                    PINPT%nelect = nint(xnelect)
                elseif(ik .gt. 1) then
                    PINPT%nelect = nint(xnelect)
                    if(nelect_prev .ne. nint(xnelect)) then
                        write(message,'(A, I5, I5, I5)')' ERROR. !!! ne(k) /= ne(k', nint(xnelect), &
                                                        nelect_prev, ik; write_msg
                        kill_job
                    endif
                    nelect_prev = PINPT%nelect
                endif
            enddo
        enddo
    endif

    ! compute reciprocal lattice vector
    call get_reci(PINPT%b1, PINPT%b2, PINPT%b3, PINPT%a1, PINPT%a2, PINPT%a3)
    PINPT%B(:,1) = PINPT%b1
    PINPT%B(:,2) = PINPT%b2
    PINPT%B(:,3) = PINPT%b3
    ! get reciprocal properties: nbmax, nplw_max, nplwidx_max, kpts_cart
    call reci_property(WAVEC, PINPT)

    ! asign G vectors and its unique id 
    if(allocated(WAVEC%G))  deallocate(WAVEC%G)  ;  allocate(WAVEC%G(3,WAVEC%nplw_max,PINPT%nkpts))
    if(allocated(WAVEC%Gid))deallocate(WAVEC%Gid);  allocate(WAVEC%Gid(WAVEC%nplw_max,PINPT%nkpts))
    do ik = 1, PINPT%nkpts
        WAVEC%G(:,:,ik) = iG(WAVEC, PINPT, ik)
        call get_Gid(WAVEC%Gid(:,ik), WAVEC%G(:,:,ik), WAVEC%nplw(ik), WAVEC%nplw_max, WAVEC%nbmax)
    enddo

    return
  endsubroutine

  function Cnk(WAVEC, PINPT, ie, is, ik) result(coeff2)
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    integer(kind=sp)    ie, is, ik, ispinor
    integer(kind=sp)    irecl
    integer(kind=sp)    iplw, nplw, gid
    complex(kind=spr)   coeff_(WAVEC%nplw_max)
    complex(kind=dp)    coeff (WAVEC%nplw_max)
    complex(kind=dp)    coeff2(WAVEC%plw_idx_max,PINPT%ispinor)
   !complex(kind=dp)    coeff2(WAVEC%nplw(ik)/PINPT%ispinor,PINPT%ispinor)
    logical             flag_norm

    flag_norm = .FALSE.
    coeff = (0d0,0d0)
    coeff_= (0d0,0d0)
    coeff2= (0d0,0d0)
    nplw  = WAVEC%nplw(ik)/PINPT%ispinor

    irecl = 3 + (ik-1)*(PINPT%nband+1) + PINPT%nkpts*(PINPT%nband+1)*(is-1) + ie
    read(pid_wavecar, rec=irecl) coeff_(1:nplw*PINPT%ispinor)
    
    coeff = coeff_

    if(flag_norm) coeff = coeff / sqrt(dot_product(coeff, coeff))
 
    ! rearrange coefficient with G vector id
    do ispinor = 1, PINPT%ispinor
        coeff2(WAVEC%Gid(1:nplw,ik),ispinor) = coeff(1+(ispinor-1)*nplw:nplw*ispinor) 
    enddo
   !do ispinor = 1, PINPT%ispinor
   !    coeff2(:,ispinor) = coeff(1+(ispinor-1)*nplw:nplw*ispinor) 
   !enddo

    return
  endfunction

  ! k + g
  function kpg(kp, g, B) result(kg)
    implicit none
    real(kind=dp),      intent(in)   :: kp(3), B(3,3)
    integer(kind=sp),   intent(in)   :: g(3)
    real(kind=dp)                    :: kg(3)

    kg(1) = dot_product(kp(:) + real(g(:)), B(1,:))
    kg(2) = dot_product(kp(:) + real(g(:)), B(2,:))
    kg(3) = dot_product(kp(:) + real(g(:)), B(3,:))
    
    return
  endfunction

  ! velocity operator between two state <phi_j|v|phi_i>
  subroutine get_vel_expectation(vex, ie, je, is, ik, WAVEC, PINPT)
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    integer(kind=sp)    ie, je, is, ik
    complex(kind=dp)    Cjk(WAVEC%plw_idx_max,PINPT%ispinor)
    complex(kind=dp)    vCik(WAVEC%plw_idx_max,PINPT%ispinor,3) ! -i*hbar*d/dx,y,z|Cik>
    complex(kind=dp)    vex(3) ! velocity expectation value
    integer(kind=sp)    ispinor, iaxis

    vex     = 0d0
    
    Cjk     = Cnk( WAVEC, PINPT, je, is, ik)
    vCik    = vCnk(WAVEC, PINPT, ie, is, ik)

    do ispinor = 1, PINPT%ispinor
        do iaxis = 1, 3
            vex(iaxis)  = vex(iaxis) + dot_product(Cjk(:,ispinor),vCik(:,ispinor,iaxis))
        enddo
    enddo

    return
  endsubroutine
 
  ! v|Cnk> where v = (vx, vy, vz) and vx = -i*hbar * d/dx
  function vCnk(WAVEC, PINPT, ie, is, ik)
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    integer(kind=sp)    ie, is, ik
    integer(kind=sp)    iplw, ispinor
    complex(kind=dp)    vCnk(WAVEC%plw_idx_max,PINPT%ispinor,3)
    complex(kind=dp)    coeff(WAVEC%plw_idx_max,PINPT%ispinor)
    real(kind=dp)       kg(3)

    ! derivation of wavefunction
    ! psi_nk(r) = sum_G C_n,k+G exp[i(k+G)r]
    ! d/dx = i(kx+Gx) * sum_G C_n,k+G exp[i(k+G)r]
    ! => -i * hbar d/dx = hbar * (kx+Gx) * sum_G C_n,k+G exp[i(k+G)r]

    vCnk = 0d0
    coeff = Cnk(WAVEC, PINPT, ie, is, ik)
    do iplw = 1, WAVEC%nplw(ik)/PINPT%ispinor
        kg = kpg(WAVEC%kpts(:,ik), WAVEC%G(:,iplw,ik), PINPT%B) !* hbar = 1. for the simplicity
        do ispinor = 1, PINPT%ispinor
            vCnk(WAVEC%Gid(iplw,ik),ispinor,:) = kg(:) * coeff(WAVEC%Gid(iplw,ik),ispinor)
        enddo
    enddo

    return
  endfunction

  ! generate G-vectors: (G+k)^2  < ENCUT
  ! ENCUT specifies the cutoff energy for the plane-wave-basis set in eV.
  ! |G + k| < G_cut with ENCUT = G_cut^2 / c 
  ! --> VASP convention:https://www.vasp.at/wiki/index.php/ENCUT)
  function iG(WAVEC, PINPT, ik) result(G)
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    integer(kind=sp)    ik
    real(kind=dp)       kp(3) ! reciproal 
    real(kind=dp)       kG(3) ! G + kp cartesian
    real(kind=dp)       kGxp, kGyp, kGzp
   !integer(kind=sp)    G(3,WAVEC%nplw(ik))
    integer(kind=sp)    G(3,WAVEC%nplw_max)
    integer(kind=sp)    iGx, iGy, iGz
    integer(kind=sp)    iGxp, iGyp, iGzp
    integer(kind=sp)    nbmax(3)
    integer(kind=sp)    nplw, iplw

    nbmax   = WAVEC%nbmax
    kp      = WAVEC%kpts(:,ik)
    iplw    = 0
    nplw    = WAVEC%nplw(ik)
    G       = 0

    do iGz = 0, 2 * nbmax(3)
        iGzp = igp(iGz, nbmax(3))
        kGzp = real(iGzp) + kp(3)
        do iGy = 0, 2 * nbmax(2)
            iGyp = igp(iGy, nbmax(2))
            kGyp = real(iGyp) + kp(2)
            do iGx = 0, 2 * nbmax(1)
                iGxp = igp(iGx, nbmax(1))
                kGxp = real(iGxp) + kp(1)
                
                kG(1) = kGxp*PINPT%b1(1) + kGyp*PINPT%b2(1) + kGzp*PINPT%b3(1)
                kG(2) = kGxp*PINPT%b1(2) + kGyp*PINPT%b2(2) + kGzp*PINPT%b3(2)
                kG(3) = kGxp*PINPT%b1(3) + kGyp*PINPT%b2(3) + kGzp*PINPT%b3(3)

                if( dot_product(kG,kG)/c .lt. PINPT%encut ) then
                    iplw = iplw + 1
                    G(1,iplw) = iGxp
                    G(2,iplw) = iGyp
                    G(3,iplw) = iGzp
                endif

            enddo
        enddo
    enddo

    if(PINPT%ispinor * iplw .ne. nplw) then
        write(message,'(A)')   ' ERROR - computed number of plane waves are not equal to'
        write(message,'(A)')   '         the stored value in the file'
        write(message,'(A,I0)')'         computed NPLW*ISPINOR = ', iplw * PINPT%ispinor
        write(message,'(A,I0)')'         stored   NPLW*ISPINOR = ', nplw
    endif

    if(PINPT%ispinor .eq. 2) then
        G(:,iplw+1:nplw) = G(:,1:iplw)
    endif

    return
  endfunction

  function igp(ig, nbmax) 
    implicit none
    integer(kind=sp)    nbmax
    integer(kind=sp)    ig, igp

    if(ig .gt. nbmax) then
        igp = ig - 2 * nbmax - 1
    else
        igp = ig
    endif

    return
  endfunction

  subroutine get_Gid(Gid, G, nplw, nplw_max, nbmax)
    implicit none
    integer(kind=sp)    nplw, nplw_max, i
    integer(kind=sp)    nbmax(3)
    integer(kind=sp)    Gid(nplw_max) 
    integer(kind=sp)    G(3,nplw_max)

    Gid = 0
    Gid(1:nplw) = Gidx(G(1,1:nplw),G(2,1:nplw),G(3,1:nplw), &
                   nbmax(1),   nbmax(2),   nbmax(3))

    return
  endsubroutine

  ! G vector indexing
  elemental integer(kind=sp) function Gidx(iGx,iGy,iGz,ng1,ng2,ng3)
    implicit none
    integer(kind=sp), intent(in) :: ng1,ng2,ng3
    integer(kind=sp), intent(in) :: igx,igy,igz

    Gidx = (iGz+ng3)*(2*ng2+1)*(2*ng1+1) + &
           (iGx+ng2)*(2*ng1+1)           + &
           (iGy+ng1)                     + 1

    return
  endfunction

  subroutine Enk_read(WAVEC, PINPT, is, ik)
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    integer(kind=sp)    ik, is, ie
    integer(kind=sp)    irec ! record addres for "k"-point
    real(kind=dp)       xnplw
    complex(kind=dp)    ce(PINPT%nband)

    irec = 3+(ik-1)*(PINPT%nband+1)+PINPT%nkpts*(PINPT%nband+1)*(is-1)
    read(pid_wavecar,rec=irec) xnplw, WAVEC%kpts(1:3,ik), (ce(ie), WAVEC%OCC(ie,is,ik),ie=1,PINPT%nband)
    WAVEC%nplw(ik) = nint(xnplw)
    WAVEC%E(1:PINPT%nband,is,ik) = real(ce(1:PINPT%nband))

    return
  endsubroutine

  ! computing reciprocal properties
  subroutine reci_property(WAVEC, PINPT)
    implicit none
    type(incar  )    :: PINPT
    type(eigen  )    :: WAVEC
    real(kind=dp)       b1mag,b2mag,b3mag ! magnitude of b vector
    real(kind=dp)       phi12, phi13, phi23, sinphi123, vmag, vtmp(3)
    integer(kind=sp)    ik, i
    integer(kind=sp)    nb1maxA, nb2maxA, nb3maxA, npmaxA
    integer(kind=sp)    nb1maxB, nb2maxB, nb3maxB, npmaxB
    integer(kind=sp)    nb1maxC, nb2maxC, nb3maxC, npmaxC

    b1mag = dsqrt( dot_product(PINPT%b1, PINPT%b1) )
    b2mag = dsqrt( dot_product(PINPT%b2, PINPT%b2) )
    b3mag = dsqrt( dot_product(PINPT%b3, PINPT%b3) )

    phi12     = acos(  dot_product(PINPT%b1, PINPT%b2) / (b1mag*b2mag))
    call vcross(vtmp, PINPT%b1, PINPT%b2)
    vmag      = dsqrt( dot_product(vtmp, vtmp) )
    sinphi123 = dot_product(PINPT%b3, vtmp) / (vmag * b3mag)
    nb1maxA   = dsqrt(PINPT%encut * c) / (b1mag * abs(sin(phi12))) + 1
    nb2maxA   = dsqrt(PINPT%encut * c) / (b2mag * abs(sin(phi12))) + 1
    nb3maxA   = dsqrt(PINPT%encut * c) / (b3mag * abs(sinphi123))  + 1
    npmaxA    = nint(4.d0 * pi * nb1maxA * nb2maxA * nb3maxA / 3.d0)

    phi13     = acos(  dot_product(PINPT%b1, PINPT%b3) / (b1mag*b3mag))
    call vcross(vtmp, PINPT%b1, PINPT%b3)
    vmag      = dsqrt( dot_product(vtmp, vtmp) )
    sinphi123 = dot_product(PINPT%b2, vtmp) / (vmag * b2mag)
    nb1maxB   = dsqrt(PINPT%encut * c) / (b1mag * abs(sin(phi13))) + 1
    nb2maxB   = dsqrt(PINPT%encut * c) / (b2mag * abs(sinphi123))  + 1
    nb3maxB   = dsqrt(PINPT%encut * c) / (b3mag * abs(sin(phi13))) + 1
    npmaxB    = nint(4.d0 * pi * nb1maxB * nb2maxB * nb3maxB / 3.d0)

    phi23     = acos(  dot_product(PINPT%b2, PINPT%b3) / (b2mag*b3mag))
    call vcross(vtmp, PINPT%b2, PINPT%b3)
    vmag      = dsqrt( dot_product(vtmp, vtmp) )
    sinphi123 = dot_product(PINPT%b1, vtmp) / (vmag * b1mag)
    nb1maxC   = dsqrt(PINPT%encut * c) / (b1mag * abs(sinphi123))  + 1
    nb2maxC   = dsqrt(PINPT%encut * c) / (b2mag * abs(sin(phi23))) + 1
    nb3maxC   = dsqrt(PINPT%encut * c) / (b3mag * abs(sin(phi23))) + 1
    npmaxC    = nint(4.d0 * pi * nb1maxC * nb2maxC * nb3maxC / 3.d0)

    WAVEC%nbmax(1) = max(nb1maxA, nb1maxB, nb1maxC)
    WAVEC%nbmax(2) = max(nb2maxA, nb2maxB, nb2maxC)
    WAVEC%nbmax(3) = max(nb3maxA, nb3maxB, nb3maxC)

    ! multiply 'ispinor' to handle two component spinors
    WAVEC%nplw_max = PINPT%ispinor * min(npmaxA, npmaxB, npmaxC)

    do ik=1, PINPT%nkpts
        do i = 1, 3
            WAVEC%kpts_cart(i,ik) = WAVEC%kpts(1,ik) * PINPT%b1(i) + &
                                    WAVEC%kpts(2,ik) * PINPT%b2(i) + &
                                    WAVEC%kpts(3,ik) * PINPT%b3(i)
        enddo
    enddo

   !Gidx = (iGz+ng3)*(2*ng2+1)*(2*ng1+1) + &
   !       (iGx+ng2)*(2*ng1+1)           + &
   !       (iGy+ng1)                     + 1
   !WAVEC%plwidx_max ->  maximum index number (Gid) for plane wave G vector indexing
   !                     maximum is optained when
   !                  -> iGz = WAVEC%nbmax(3)
   !                  -> iGy = WAVEC%nbmax(2)
   !                  -> iGx = WAVEC%nbmax(1)
    WAVEC%plw_idx_max = 8 * WAVEC%nbmax(1) * WAVEC%nbmax(2) * WAVEC%nbmax(3) + & 
                        4 * WAVEC%nbmax(1) * WAVEC%nbmax(2)                  + &
                        4 * WAVEC%nbmax(2) * WAVEC%nbmax(3)                  + &
                        4 * WAVEC%nbmax(3) * WAVEC%nbmax(1)                  + &
                        2 * WAVEC%nbmax(1)                                   + &
                        2 * WAVEC%nbmax(2)                                   + &
                        2 * WAVEC%nbmax(3)

    
    return
  endsubroutine

endmodule
