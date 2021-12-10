#include "alias.inc"

subroutine write_info(WAVEC, PINPT, pid)
    use mpi_setup
    use mykind
    use parameters
    use print_io
    implicit none
    type(incar  )   ::  PINPT
    type(eigen  )   ::  WAVEC
    integer(kind=sp)    i
    integer(kind=sp)    pid

    if(PINPT%flag_set_unfold) return

    if(pid .eq. 0) then
        write(message,'(A,I6)')"# NELECT     : ",PINPT%nelect   ; write_msg
        write(message,'(A,I6)')"# ISPIN      : ",PINPT%ispin    ; write_msg
        if(PINPT%ispinor .eq. 2) then
            write(message,'(A,I6,A)')"# ISPINOR    : ",PINPT%ispinor," (LSORBIT =.TRUE.)" ; write_msg
        else
            write(message,'(A,I6,A)')"# ISPINOR    : ",PINPT%ispinor," (LSORBIT =.FALSE.)" ; write_msg
        endif
        
        write(message,'(A,F11.4)')  "# ENCUT (eV) : ",PINPT%encut           ; write_msg
        write(message,'(A,I6)')     "# NKPOINT    : ",PINPT%nkpts           ; write_msg
!       write(message,'(A,I6,A,I4)')"#  K-GRID    : ",nkx,"   X",nky        ; write_msg
        write(message,'(A,I6)')     "# NBANDS     : ",PINPT%nband           ; write_msg
        write(message,'(A,3F13.6)') "# LATTVEC A1 : ",(PINPT%a1(i),i=1,3)   ; write_msg
        write(message,'(A,3F13.6)') "# LATTVEC A2 : ",(PINPT%a2(i),i=1,3)   ; write_msg
        write(message,'(A,3F13.6)') "# LATTVEC A3 : ",(PINPT%a3(i),i=1,3)   ; write_msg
        
        write(message,'(A,3F13.6)') "# RECIVEC B1 : ",(PINPT%b1(i),i=1,3)   ; write_msg
        write(message,'(A,3F13.6)') "# RECIVEC B2 : ",(PINPT%b2(i),i=1,3)   ; write_msg
        write(message,'(A,3F13.6)') "# RECIVEC B3 : ",(PINPT%b3(i),i=1,3)   ; write_msg
        
        write(message,'(A,I6)')     "# NPMAX      : ",WAVEC%nplw_max        ; write_msg
        write(message,'(A,I6)')     "# NB1MAX     : ",WAVEC%nbmax(1)        ; write_msg
        write(message,'(A,I6)')     "# NB2MAX     : ",WAVEC%nbmax(2)        ; write_msg
        write(message,'(A,I6)')     "# NB3MAX     : ",WAVEC%nbmax(3)        ; write_msg

    elseif(pid .gt. 10) then
        write(pid,'(A,I6)')"# NELECT     : ",PINPT%nelect   
        write(pid,'(A,I6)')"# ISPIN      : ",PINPT%ispin    
        if(PINPT%ispinor .eq. 2) then
            write(pid,'(A,I6,A)')"# ISPINOR    : ",PINPT%ispinor," (LSORBIT =.TRUE.)" 
        else
            write(pid,'(A,I6,A)')"# ISPINOR    : ",PINPT%ispinor," (LSORBIT =.FALSE.)" 
        endif

        write(pid,'(A,F11.4)')  "# ENCUT (eV) : ",PINPT%encut           
        write(pid,'(A,I6)')     "# NKPOINT    : ",PINPT%nkpts           
!       write(pid,'(A,I6,A,I4)')"#  K-GRID    : ",nkx,"   X",nky        
        write(pid,'(A,I6)')     "# NBANDS     : ",PINPT%nband           
        write(pid,'(A,3F13.6)') "# LATTVEC A1 : ",(PINPT%a1(i),i=1,3)   
        write(pid,'(A,3F13.6)') "# LATTVEC A2 : ",(PINPT%a2(i),i=1,3)   
        write(pid,'(A,3F13.6)') "# LATTVEC A3 : ",(PINPT%a3(i),i=1,3)   

        write(pid,'(A,3F13.6)') "# RECIVEC B1 : ",(PINPT%b1(i),i=1,3)   
        write(pid,'(A,3F13.6)') "# RECIVEC B2 : ",(PINPT%b2(i),i=1,3)   
        write(pid,'(A,3F13.6)') "# RECIVEC B3 : ",(PINPT%b3(i),i=1,3)   

        write(pid,'(A,I6)')     "# NPMAX      : ",WAVEC%nplw_max        
        write(pid,'(A,I6)')     "# NB1MAX     : ",WAVEC%nbmax(1)        
        write(pid,'(A,I6)')     "# NB2MAX     : ",WAVEC%nbmax(2)        
        write(pid,'(A,I6)')     "# NB3MAX     : ",WAVEC%nbmax(3)        
    endif

    return
endsubroutine
