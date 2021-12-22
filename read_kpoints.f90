#include "alias.inc"

subroutine get_kpath_from_file(fname, PGEOM)
    use parameters, only: poscar, incar, pid_kpts
    use mykind
    use mpi_setup
    use utils
    implicit none
    type(poscar )               :: PGEOM
    character(len=*)             fname
    integer(kind=sp), parameter :: max_kline=100
    integer(kind=sp)               i_continue
    integer(kind=sp)               i,linecount, i_dummy
    integer(kind=sp)               idiv_mode, ik
    integer(kind=sp)               iline
    real(kind=dp)                  kline_dummy(3,max_kline)
    character(len=132)             inputline
    character(len=40)              desc_str, k_name_dummy(max_kline)
    logical                        flag_skip
    real(kind=dp), allocatable  :: kpts_cart(:,:), kpts_reci(:,:)
    character(*), parameter     :: func = 'read_kpoint'

    open(pid_kpts, FILE=trim(fname),iostat=i_continue)
    linecount = 0

 line: do
        read(pid_kpts,'(A)',iostat=i_continue) inputline
        if(i_continue<0) exit               ! end of file reached
        if(i_continue>0) then
          write(message,'(3A)')'Unknown error reading file:',trim(fname),func  ; write_msgi
          kill_job
        endif
        linecount = linecount + 1
        call check_comment(inputline,linecount,i,flag_skip) ; if (flag_skip ) cycle
        call check_empty(inputline,linecount,i,flag_skip) ; if(flag_skip) cycle

        if(i_continue .ne. 0) cycle              ! skip empty line

        ! head
         if(linecount .eq. 1) then
           cycle

        ! dividing factor
         elseif(linecount .eq. 2) then
           read(inputline,*,iostat=i_continue) PGEOM%nkdiv
           cycle

        ! k-grid or -line mode
         elseif(linecount .eq. 3) then
           read(inputline,*,iostat=i_continue) desc_str

         ! k-vector type
         elseif(linecount .eq. 4) then
           read(inputline,*,iostat=i_continue) desc_str

         ! k-line if 'linemode'
         elseif(linecount .ge. 5 ) then
           backspace(pid_kpts)
           PGEOM%nline=0;i=0

    kline: do
             read(pid_kpts,'(A)',iostat=i_continue) inputline
             if(i_continue<0) exit               ! end of file reached
             if(i_continue>0) then
               write(message,*)'Unknown error reading file:',trim(fname),func  ; write_msgi
               kill_job
             endif
             linecount = linecount+1 ; i=i+1
             call check_comment(inputline,linecount,i,flag_skip) ; if (flag_skip ) cycle kline
             call check_empty(inputline,linecount,i,flag_skip) ; if (flag_skip) cycle kline
             read(inputline,*) kline_dummy(1:3,i),k_name_dummy(i)
             if( mod(i,2) .eq. 1 .and. i .ge. 1) then
               PGEOM%nline = PGEOM%nline + 1
             endif
           enddo kline

           if(allocated(PGEOM%kline))      deallocate(PGEOM%kline)      ; allocate(PGEOM%kline(3,PGEOM%nline*2))
           if(allocated(PGEOM%kline_cart)) deallocate(PGEOM%kline_cart) ; allocate(PGEOM%kline_cart(3,PGEOM%nline*2))
           if(allocated(PGEOM%k_name)) deallocate(PGEOM%k_name); allocate(PGEOM%k_name(PGEOM%nline*2))
           PGEOM%kline(1:3,1:PGEOM%nline*2) = kline_dummy(1:3,1:PGEOM%nline * 2)
           PGEOM%k_name(1:PGEOM%nline*2) = k_name_dummy(1:PGEOM%nline * 2)
         endif

         if(i_continue<0) exit ! I put this line for the gfortran issue, but I'm
      enddo line

    if (linecount == 0) then
      write(message,*)'Attention - empty input file: ',trim(fname) ; write_msgi
      stop
    endif
    close(pid_kpts)
  

    if(allocated(PGEOM%kpts))      deallocate(PGEOM%kpts)      ; allocate(PGEOM%kpts(3,PGEOM%nkdiv*PGEOM%nline))
    if(allocated(PGEOM%kpts_cart)) deallocate(PGEOM%kpts_cart) ; allocate(PGEOM%kpts_cart(3,PGEOM%nkdiv*PGEOM%nline))

    call get_kpath(PGEOM%kline, PGEOM%nkdiv, PGEOM%nline, PGEOM%kpts)
    PGEOM%nkpts = PGEOM%nline * PGEOM%nkdiv
    do iline = 1, PGEOM%nline * 2
        PGEOM%kline_cart(1:3, iline) = PGEOM%kline(1,iline) * PGEOM%b_latt(1:3,1) + &
                                       PGEOM%kline(2,iline) * PGEOM%b_latt(1:3,2) + &
                                       PGEOM%kline(3,iline) * PGEOM%b_latt(1:3,3)
    enddo
    call get_kpath(PGEOM%kline_cart, PGEOM%nkdiv, PGEOM%nline, PGEOM%kpts_cart)

    if(allocated(PGEOM%kdist))        deallocate(PGEOM%kdist)        ; allocate(PGEOM%kdist(PGEOM%nkpts))


    if(allocated(PGEOM%k_name_index)) deallocate(PGEOM%k_name_index) ; allocate(PGEOM%k_name_index(PGEOM%nline+1))
    if(allocated(PGEOM%k_name2))      deallocate(PGEOM%k_name2)      ; allocate(PGEOM%k_name2(PGEOM%nline+1))

    PGEOM%k_name2(1) = trim(PGEOM%k_name(1))
    do ik = 2, PGEOM%nline+1
      PGEOM%k_name2(ik) = trim( PGEOM%k_name((ik-1)*2) )
    enddo
    do ik = 1, PGEOM%nline
      PGEOM%k_name_index(ik) = PGEOM%nkdiv*(ik-1)+1
    enddo
    PGEOM%k_name_index(PGEOM%nline+1) = PGEOM%nkdiv*PGEOM%nline
  
    return
endsubroutine

subroutine get_kpath(PK, nkdiv, nline, kpath)
    use mykind
    implicit none
    integer(kind=sp)    nkdiv, nline
    real(kind=dp)       PK(3,nline*2), dk(3) ! high_symm_kpts
    real(kind=dp)       kpath(3,nkdiv*nline)
    integer(kind=sp)    iline
    integer(kind=sp)    ik, ii

    ii = 0
    do iline = 1, nline * 2, 2
        dk(1:3) = ( PK(:,iline + 1) - PK(:,iline) ) / dble(nkdiv-1) ! if ndivk=1 it crashes -> need solution
        do ik = 1, nkdiv
            ii = ii + 1
            kpath(1:3, ii) = PK(1:3, iline) + dk(1:3) * dble(ik-1)
        enddo
    enddo

   return
endsubroutine

subroutine prepare_unfold_kpoint(PINPT, PGEOM_PC, PGEOM_SC, imode, flag_reduce)
    use parameters, only: incar, poscar
    use utils
    use mykind
    use print_io
    use mpi_setup
    use utils
    implicit none
    type(incar  )        :: PINPT
    type(poscar )        :: PGEOM_PC, PGEOM_SC
    character(len=256)      kfilenm_pc, kfilenm_sc, filenm_pc, filenm_sc
    real(kind=dp)           M(3,3)
    integer(kind=sp)        imode
    logical                 flag_reduce ! if .TRUE. remove duplicated k-points

   !imode = 1 ! report k-point to be unfolded

    write(kfilenm_pc,'(A)') trim(PINPT%folder_out)//'KPOINTS_PC'
    write(kfilenm_sc,'(A)') trim(PINPT%folder_out)//'KPOINTS_SC'
    write(filenm_pc, '(A)') trim(PINPT%folder_out)//'POSCAR_PC'
    write(filenm_sc, '(A)') trim(PINPT%folder_out)//'POSCAR_SC'

    write(message,'(A)')' ' ; write_msgi
    write(message,'(A)')'#- READING POSCAR  FILE: ' ; write_msgi

    write(message,'(A)')'  --> '//trim(filenm_pc) ; write_msgi
    call read_poscar(filenm_pc, PGEOM_PC)
        write(message,'(A      )'  )"  Primitive cell lattice vector A_PC";write_msgi
        write(message,'(A,3F9.4,A)')"    A_PC = [", PGEOM_PC%a_latt(:,1),';' ; write_msgi
        write(message,'(A,3F9.4,A)')"            ", PGEOM_PC%a_latt(:,2),';' ; write_msgi
        write(message,'(A,3F9.4,A)')"            ", PGEOM_PC%a_latt(:,3),']' ; write_msgi
        write(message,'(A      )'  )"  Primitive cell reciprocal lattice vector B_PC";write_msgi
        write(message,'(A,3F9.4,A)')"    B_PC = [", PGEOM_PC%b_latt(:,1),';' ; write_msgi
        write(message,'(A,3F9.4,A)')"            ", PGEOM_PC%b_latt(:,2),';' ; write_msgi
        write(message,'(A,3F9.4,A)')"            ", PGEOM_PC%b_latt(:,3),']' ; write_msgi

    write(message,'(A)')' ' ; write_msgi
    write(message,'(A)')'  --> '//trim(filenm_sc) ; write_msgi
    call read_poscar(filenm_sc, PGEOM_SC)
        write(message,'(A      )'  )"  Primitive cell lattice vector A_PC";write_msgi
        write(message,'(A,3F9.4,A)')"    A_SC = [", PGEOM_SC%a_latt(:,1),';' ; write_msgi
        write(message,'(A,3F9.4,A)')"            ", PGEOM_SC%a_latt(:,2),';' ; write_msgi
        write(message,'(A,3F9.4,A)')"            ", PGEOM_SC%a_latt(:,3),']' ; write_msgi
        write(message,'(A      )'  )"  Primitive cell lattice vector A_PC";write_msgi
        write(message,'(A,3F9.4,A)')"    B_SC = [", PGEOM_SC%b_latt(:,1),';' ; write_msgi
        write(message,'(A,3F9.4,A)')"            ", PGEOM_SC%b_latt(:,2),';' ; write_msgi
        write(message,'(A,3F9.4,A)')"            ", PGEOM_SC%b_latt(:,3),']' ; write_msgi

    call get_transformation_matrix(PGEOM_PC%a_latt, PGEOM_SC%a_latt, PGEOM_PC%M)
    M = PGEOM_PC%M
    PGEOM_SC%M = M
        write(message,'(A        )')' ' ; write_msgi
        write(message,'(A        )')'  Transformation matrix M: M*A_PC = A_SC' ; write_msgi
        write(message,'(A,3F9.3,A)')'      M = [',M(1,:),';'               ; write_msgi
        write(message,'(A,3F9.3,A)')'           ',M(2,:),';'               ; write_msgi
        write(message,'(A,3F9.3,A)')'           ',M(3,:),']'               ; write_msgi

    write(message,'(A)')' ' ; write_msgi
    write(message,'(A)')'#- READING KPOINTS (PC) FILE: ' ; write_msgi
    write(message,'(A)')'  --> '//trim(kfilenm_pc) ; write_msgi
    call get_kpath_from_file(trim(kfilenm_pc), PGEOM_PC)

    write(message,'(A)')' ' ; write_msgi
    write(message,'(A)')'# Mapping primitive cell k-point to supercell BZ '
    call get_unfold_klist(PGEOM_PC, PGEOM_SC, M, imode, flag_reduce)

    return
endsubroutine

subroutine get_unfold_klist(KP_PC, KP_SC, M, imode, flag_reduce)
    use parameters, only: poscar
    use mykind
    use print_io
    use mpi_setup
    implicit none
    type(poscar )                :: KP_PC, KP_SC
    integer(kind=sp)                ik
    real(kind=dp)                   M(3,3), MT(3,3)
    real(kind=dp)                   B_PC(3,3), B_SC(3,3)
    real(kind=dp)                   k(3), KK(3), K_1BZ(3)
    real(kind=dp),allocatable    :: kpts(:,:)
    integer(kind=sp),allocatable :: gpts(:,:)
    integer(kind=sp)                G(3)
    integer(kind=sp)                imode ! report (1) or not (0)
    integer(kind=sp)                nkpts_unique
    logical                         flag_reduce ! remove duplicate k

    ! Note:
    ! * PC lattice vectors and its transformation into SC by transformation matrix M
    !       M * A_PC = A_SC
    !        --> M = A_SC * inv(A_PC)
    ! 
    ! * SC reciprocal lattice vectors and its transformation into PC by transformation matrix T
    !       B_PC = T * B_SC 
    !        --> T = B_PC * inv(B_SC)
    !
    ! * Relation between M and T
    !       T = transpose(M) = MT
    !        --> B_PC = MT * B_SC
    !        --> MT = B_PC * inv(B_SC)
    !
    ! * KK vector of SC and k vector of PC. 
    !       k       = (kx, ky, kz) ; 
    !       k_cart  = (kx, ky, kz) * (B1_PC, B2_PC, B3_PC) = kx*B1_PC + ky*B2_PC + kz*B3_PC
    !        --> Bi_PC  = ( Bix_PC, Biy_PC, Biz_PC ) in cart
    !        --> k_cart = ( kx*B1x_PC + ky*B2x_PC + kz*B3x_PC, 
    !                       kx*B1y_PC + ky*B2y_PC + kz*B3y_PC,
    !                       kx*B1z_PC + ky*B2z_PC + kz*B3z_PC ) = k * B_PC
    !       
    !       KK       = (Kx, Ky, Kz) ;
    !       KK_cart  = (Kx, Ky, Kz) * (B1_SC, B2_SC, B3_SC) = Kx*B1_SC + Ky*B2_SC + Kz*B3_SC
    !        --> Bi_SC   = ( Bix_SC, Biy_SC, Biz_SC ) in cart
    !        --> KK_cart = ( Kx*B1x_SC + Ky*B2x_SC + Kz*B3x_SC,
    !                        Kx*B1y_SC + Ky*B2y_SC + Kz*B3y_SC,
    !                        Kx*B1z_SC + Ky*B2z_SC + Kz*B3z_SC ) = KK * B_SC
    !           
    ! * Relation between k and K : k_cart = K_cart
    !       k_cart can be written by reciprocal lattice vector of SC (B_SC)
    !        --> k_cart = k * B_PC = KK * B_SC = K_cart
    !        --> KK * B_SC = k * MT * B_SC
    !        --> KK = k * MT
    !   
    ! * Since KK is related with the K_1BZ where K_1BZ is reciprocal coordinate of first BZ of SC
    !       KK = K_1BZ + G
    !        --> K_1BZ = KK - G
    !
    ! * Our goald is to find K_1BZ and the G vector that folds K into 1BZ of SC

    ! temporarily set nkpts with KP_PC
    allocate(kpts(3,KP_PC%nkpts))   ; kpts = 0d0
    allocate(gpts(3,KP_PC%nkpts))   ; gpts = 0

    do ik=1, KP_PC%nkpts 
        k       = KP_PC%kpts(:,ik)
        call map_k_to_K_1BZ(k, M, KK, G, K_1BZ)
        kpts(:,ik)  = K_1BZ
        gpts(:,ik)  = G   
        if(imode .eq. 1) then
            ! in reciprocal unit (fractional)
            write(message,'( 2(A,3F15.8), A, 3I5, A)')"  k-point k0(PC): ", k, &
                                                            " , K0(SC): ", K_1BZ, &
                                                            " , G0(SC): ", G ; write_msgi
        endif

    enddo

    if(flag_reduce) then
        ! removing duplictated KK-points
        call find_unique_kpoints(kpts, gpts, KP_PC%nkpts, nkpts_unique)
        KP_SC%nkpts = nkpts_unique
    elseif(.not. flag_reduce) then
        KP_SC%nkpts = KP_PC%nkpts
    endif

    if(allocated(KP_SC%kpts)) deallocate(KP_SC%kpts)    ; allocate(KP_SC%kpts(3,KP_SC%nkpts)) ; KP_SC%kpts = kpts
    if(allocated(KP_SC%G   )) deallocate(KP_SC%G   )    ; allocate(KP_SC%G   (3,KP_SC%nkpts)) ; KP_SC%G    = gpts

    return
endsubroutine

subroutine map_k_to_K_1BZ(k, M, KK, G, K_1BZ)
    use mykind
    use parameters, only: eta
    implicit none
    real(kind=dp)       k(3), KK(3), K_1BZ(3)
    real(kind=dp)       M(3,3), MT(3,3)
    integer(kind=sp)    G(3), i

    MT      = transpose(M)
    KK      = matmul(k,MT)
   !G       = nint(KK - eta)
    G       = idnint(KK )
    K_1BZ   = KK - dble(G)

    return
endsubroutine

subroutine find_unique_kpoints(kpts, gpts, nkpts, nkpts_unique)
    use mykind
    use print_io
    use mpi_setup
    implicit none
    integer(kind=sp)    ik, jk, ii
    integer(kind=sp)    nkpts, nkpts_unique
    real(kind=dp)       kpts(3,nkpts), kpts_temp(3,nkpts)
    integer(kind=sp)    gpts(3,nkpts), gpts_temp(3,nkpts)
    real(kind=dp)       dk(3), dk2(3), ki(3), kj(3), kk(3)
    real(kind=dp)       absk, absk2, very_small
    logical             flag_unique

    very_small      = 0.000000001d0
    nkpts_unique    = 0
    kpts_temp       = 0.d0

kki:do ik = 1, nkpts
        ! pick one kpoint and if it is first round, save it to kpts_temp
        ki = kpts(1:3, ik)
        if(ik .eq. 1) then
            nkpts_unique = 1
            kpts_temp(:,nkpts_unique) = ki
            gpts_temp(:,nkpts_unique) = gpts(1:3, ik)
            cycle kki
        endif

        ! if ik is .ge. 2, it should be checked with registered members
    kkj:do jk = 1, nkpts_unique
            kj      = kpts_temp(:, jk)
            dk      = kj - ki
            absk    = sqrt(dot_product(dk,dk))
            if(absk .lt. very_small) then
                ! since ki is already exist in the kpts_temp,
                ! cycle kki loop
                cycle kki
            endif
        enddo kkj

        ! if ki is not registered yet, one should arrive at this point
        ! and should be registered as a new kpoint member to kpts_temp
        nkpts_unique = nkpts_unique + 1
        kpts_temp(:, nkpts_unique) = ki
        gpts_temp(:, nkpts_unique) = gpts(1:3, ik)
    enddo kki

    kpts = 0d0
    kpts = kpts_temp
    gpts = 0
    gpts = gpts_temp

    write(message,'(A,2(I5,A))')'  --> find ',nkpts_unique, &
                                ' unique k-points out of ',nkpts, &
                                ' k-points..' ; write_msgi

    return
endsubroutine

subroutine map_KK_to_k(k, M, KK)
    use mykind
    use do_math, only: inv
    implicit none
    real(kind=dp)       k(3), KK(3)
    real(kind=dp)       M(3,3), M_(3,3)

    ! NOTE:
    !
    ! * SC reciprocal lattice vectors and its transformation into PC by transformation matrix T
    !
    !       B_PC = MT * B_SC   <-- MT = transpose(M)
    !        --> B_SC = inv(MT) * B_PC      <-- transpose(inv(M)) = inv(transpose(M))
    !        -->      = transpose(inv(M)) * B_PC
    !
    ! * relation betwen k and KK
    !       KK  = (Kx, Ky, Kz) ! reciprocal lattice vector with SC BZ
    !       k   = (kx, ky, kz) ! reciprocal lattice vector with PC BZ
    !
    !        --> B_SC = inv(MT) * B_PC
    !        --> K_cart = KK * B_SC = KK * inv(MT) * B_PC = k_cart = k * B_PC
    !        -->                      |---------|
    !                                      |--> k = KK * inv(MT)
    !
    ! * k = KK * inv(MT)

    M_   = transpose(inv(M))
    k    = matmul(KK, M_)

    return
endsubroutine

subroutine screen_PC_G_von_SC_G(ik_W, WAVEC, M, Gshift, ispinor)
    use parameters, only: eigen
    use mykind
    use mpi_setup
    use do_math, only: inv
    use wavecar, only: Gidx
    implicit none
    type(eigen  )    :: WAVEC
    integer(kind=sp)    G(3,WAVEC%nplw_max), nbmax(3)
    integer(kind=sp)    GG(3), Gshift(3), G_idx, G1(3)
    real(kind=dp)       G_(3), dG(3), abs_dG
    integer(kind=sp)    iplw, jplw, nplw, ik_W, ispinor
    integer(kind=sp)    iplw_screen, nplw_screen ! number of plane waves screened
    real(kind=dp)       M(3,3)
    real(kind=dp)       very_small

    very_small  = 0.000001d0   
    G           = 0 
    nplw        = WAVEC%nplw(ik_W) / ispinor
    G           = WAVEC%G(1:3,1:nplw,ik_W)
    iplw_screen = 0
    nbmax       = WAVEC%nbmax

    do iplw = 1, nplw
        call map_KK_to_k(G_, M, dble(G(:,iplw)))
        dG      = G_ - dble(nint(G_))
        abs_dG  = sqrt(dot_product(dG, dG))

        if(abs_dG .le. very_small) then
            G1          = G(:,iplw) + Gshift  ! shift G by PGEOM_SC%G
            G_idx       = Gidx(G1(1), G1(2), G1(3), nbmax(1), nbmax(2), nbmax(3))
       Gidj:do jplw = 1, nplw
                if(G_idx .eq. WAVEC%Gid(jplw,ik_W)) then
                    iplw_screen = iplw_screen + 1
                    WAVEC%GS(:,iplw_screen,ik_W) = G1
                    exit Gidj
                endif
            enddo Gidj
        endif

    enddo

    ! note: it is not multiplied with 2 the spinor indicator since 
    nplw_screen = iplw_screen
    WAVEC%nplw_screen(ik_W) = iplw_screen

    if(ispinor .eq. 2) then
        WAVEC%GS(:,nplw_screen+1:iplw_screen*2,ik_W) = WAVEC%GS(:,1:nplw_screen, ik_W)
    endif

    return
endsubroutine

subroutine K_index_von_WAVECAR(ik_W, WAVEC, K)
    use parameters, only: eigen
    use mykind
    use print_io
    use mpi_setup
    implicit none
    type(eigen  )    :: WAVEC
    real(kind=dp)       K(3), K_(3)
    real(kind=dp)       very_small
    integer(kind=sp)    ik, ik_W, nkpts

    very_small  = 0.000001d0 ! this is probably the limit due to the accuracy of KPOINTS
    nkpts       = size(WAVEC%kpts(1,:))
    ik_W        = 0


 kk:do ik = 1, nkpts
        K_ = WAVEC%kpts(:,ik)
        if(sqrt(dot_product(K_ - K,K_ - K)) .le. very_small) then
            ik_W = ik
            exit kk
        endif
    enddo kk

    if(ik_W .eq. 0) then
        write(message,'(A)')' ' ; call write_log(message,12,myid)
        write(message,'(A,I5,A,3F12.8,A)')' # The ',ik,'-th K-point (SC): [', K, &
                                          '] not exist in WAVECAR !' ; call write_log(message,12,myid)
        write(message,'(A         )')'   --> Exit program... ' ; call write_log(message,12,myid)
        kill_job1
    endif

    return
endsubroutine
