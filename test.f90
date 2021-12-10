#include "alias.inc"
subroutine test()
  implicit none
    integer i
    real*8 r
    real*8 pi 
    pi=4.d0*atan(1.d0)

    r = 180d0 * pi / 180d0
    write(6,*)"ZXX ", pi, cos(pi), cos(180d0), cos(r)
 stop
end subroutine
