#include "alias.inc"

module utils
    use mykind

contains

  function str2lowcase(string)
    character(len=*), intent(in) :: string
    character(len=len(string))   :: str2lowcase
    character(len=26)               lower_case, upper_case
    integer(kind=sp)                i, n

    lower_case = 'abcdefghijklmnopqrstuvwxyz'
    upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    str2lowcase = string
    do i=1, len(string)
        n = index(upper_case, string(i:i))
        if(n .ne. 0) str2lowcase(i:i) = lower_case(n:n)
    enddo
  endfunction

  function int2str(w) result(string)
    implicit none
    character(len=20)               string
!   character(*), intent(out)   ::  string
    integer(kind=sp), intent(in)::  w
    
    write(string,*) w
    
    return
  endfunction
      
  function flag_number(string)
    implicit none
    character(*), intent(in)    ::  string
    real(kind=dp)                   x
    integer(kind=sp)                i_number
    logical                         flag_number
    
    read(string,*,iostat=i_number) x
    flag_number = i_number == 0
    
  endfunction
   
  !computing reciprocal lattice vector    
  subroutine get_reci(b1, b2, b3, a1, a2, a3)
    use parameters, only :  pi2
    implicit none
    real(kind=dp)                   a1(3), a2(3), a3(3)
    real(kind=dp)                   b1(3), b2(3), b3(3)
    real(kind=dp)                   a2xa3(3), b1xb2(3)
    real(kind=dp)                   Vcell

    call vcross(a2xa3, a2, a3)
    Vcell = dot_product(a1,a2xa3)

    call vcross(b1,a2,a3)
    call vcross(b2,a3,a1)
    call vcross(b3,a1,a2)

    b1 = pi2 * b1 / Vcell
    b2 = pi2 * b2 / Vcell
    b3 = pi2 * b3 / Vcell

    return
  endsubroutine

  ! computing vector cross-product
  subroutine vcross(a,b,c)
    implicit none
    real(kind=dp)                   a(3), b(3), c(3)
    
    a(1)=b(2)*c(3)-b(3)*c(2)
    a(2)=b(3)*c(1)-b(1)*c(3)
    a(3)=b(1)*c(2)-b(2)*c(1)

    return
  endsubroutine

  subroutine get_kline_dist(kpts, nkpts, kline)
    implicit none
    integer(kind=sp)    ik
    integer(kind=sp)    nkpts
    real(kind=dp)       kpts(3,nkpts)
    real(kind=dp)       kline(nkpts),k0(3)
    
    do ik = 1, nkpts
        if(ik .eq. 1) then
            k0 = kpts(:,1)
            kline(1) = 0d0
        else
            k0 = kpts(:,ik-1)
            kline(ik) = kline(ik-1)
        endif
        kline(ik) = kline(ik) + sqrt(dot_product(kpts(:,ik)-k0(:),kpts(:,ik)-k0(:)))
    enddo

    return
  endsubroutine

  subroutine get_transformation_matrix(A, B, M)
    use do_math
    implicit none
    real(kind=dp)       A(3,3)
    real(kind=dp)       B(3,3)
    real(kind=dp)       M(3,3)

    ! Note: Transpose is taken due to the matrix A and B is lattice vectors
    !       saved with column major order.
    !       --> A = [A1, A2, A3] = [A11 A21 A31]
    !           Ai= [Ai1]          [A12 A22 A32]
    !               [Ai2]          [A13 A23 A33]
    !               [Ai3]
    !       We should set A = [A1^T] = [A11 A12 A13]
    !                         [A2^T]   [A21 A22 A23]
    !                         [A3^T]   [A31 A32 A33]
    !       which is the row major order.
    !       M * AT = BT --> M = BT * (AT^-1)
    M  = matmul(transpose(B), inv(transpose(A)))

    return
  endsubroutine

subroutine check_comment(inputline,linecount,i,flag_skip)
  implicit none
  integer(kind=sp)  i,linecount,i_continue
  character(*)      inputline
  character(len=40) desc_str
  logical           flag_skip
  read(inputline,*,iostat=i_continue) desc_str
  if (linecount .ne. 1 .and. desc_str(1:1).eq.'#') then
    linecount=linecount - 1
    flag_skip = .true.
    i=i - 1
  else
    flag_skip = .false.
    i=i
  endif

return
endsubroutine

subroutine check_empty(inputline,linecount,i,flag_skip)
  implicit none
  integer(kind=sp)  i,linecount,i_continue
  character(*)      inputline
  character(len=40) desc_str
  logical           flag_skip

  read(inputline,*,iostat=i_continue) desc_str
  if(i_continue .ne. 0) then
   flag_skip = .true.
   i = i - 1
   linecount = linecount - 1
  else
   flag_skip = .false.
   i = i
  endif

return
endsubroutine

function nitems(string)
  implicit none
  logical                       blank
  integer(kind=sp)              nitems,l,i
  character(*),intent(in)    :: string

  nitems=0
  l=len_trim(string)
  blank = .true.
  do i=1,l
   if(string(i:i) .eq. '#' .or. string(i:i) .eq. '!') exit

   if (blank .and. string(i:i) .ne. ' ' ) then
     blank=.false.
     nitems=nitems + 1
   elseif( .not. blank .and. string(i:i) .eq. ' ') then
     blank=.true.
   endif
  enddo
  return
endfunction

endmodule
