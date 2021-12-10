#include 'alias.inc'

subroutine read_poscar(fname, PGEOM)
    use parameters,  only : incar, poscar, pid_geom
    use print_io
    use utils
    use mpi_setup
    implicit none
    type (poscar)                    :: PGEOM
    integer(kind=sp)                    i_continue !,nitems
    integer(kind=sp)                    i,ii,linecount
    real(kind=dp)                       a_scale
    character(len=264)                  inputline
    character(len=256)                  fname
    character(len=*), parameter      :: func = 'read_poscar'
    logical                             flag_skip
 

    open(pid_geom, FILE=trim(fname), iostat=i_continue)
    linecount = 0
    ii = 0

 line: do
        read(pid_geom,'(A)',iostat=i_continue) inputline
        if(i_continue<0) exit               ! end of file reached
        if(i_continue>0) then
          write(message,*)'Unknown error reading file:',trim(fname),func  ; write_msgi
          kill_job
        endif

        if(linecount .eq. 0) then
          call check_comment(inputline,linecount,i,flag_skip)
          if(flag_skip) linecount = linecount + 1
        else
          call check_comment(inputline,linecount,i,flag_skip)
          if(flag_skip) linecount = linecount + 1
          if (flag_skip) cycle
        endif
        linecount = linecount + 1

        ! head
         if(linecount .eq. 1) then
           cycle

        ! scaling factor
         elseif(linecount .eq. 2) then
           read(inputline,*,iostat=i_continue) a_scale
           cycle

        ! lattice parameter
         elseif(linecount .eq. 3 ) then
           backspace(pid_geom)
           do i=1,3
             read(pid_geom,'(A)',iostat=i_continue) inputline
             read(inputline,*,iostat=i_continue) PGEOM%a_latt(1:3,i)
             PGEOM%a_latt(1:3,i) = PGEOM%a_latt(1:3,i) * a_scale ! rescale with SCALE factor
           enddo
           call get_reci(PGEOM%b_latt(:,1), PGEOM%b_latt(:,2), PGEOM%b_latt(:,3), &
                         PGEOM%a_latt(:,1), PGEOM%a_latt(:,2), PGEOM%a_latt(:,3))
           linecount = linecount + 2
           cycle
         endif
        enddo line

    close(pid_geom)

    return
endsubroutine
