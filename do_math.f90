#include "alias.inc"
module do_math
    use parameters
    use mykind

contains

  elemental real(kind=dp) function fgauss(sigma, x)
    implicit none
    real(kind=dp),  intent(in)   :: sigma
    real(kind=dp),  intent(in)   :: x
    real(kind=dp)                   sigma2, xx

    xx      = x**2 
    sigma2  = sigma**2
    fgauss  = exp(-xx/(2d0*sigma2)) / (sigma*sqrt(pi2))

    return
  endfunction

  elemental real(kind=dp) function florentz(sigma, x)
    implicit none
    real(kind=dp),  intent(in)   :: sigma
    real(kind=dp),  intent(in)   :: x
    real(kind=dp)                   sigma2, xx
    
    ! NOTE: Lorentzian (https://en.wikipedia.org/wiki/Spectral_line_shape)
    !       A Lorentzian line shape function can be represented as
    !                 1                  p - p0
    !       L  = -----------   ,  x  = ----------
    !              1  + x^2               w/2
    !       where L signifies a Lorentzian function standardized,
    !       for spectroscopic purpose, to a maximum value of 1;
    !       x is a subsidiary variable where p0 is the position of the 
    !       maximum (corresponding to the transition energy E),
    !       p is a position, and w is the full width at half maximum (FWHM), 
    !       the width of the curve when the intensity is half the maximum
    !       intensity (--> p = p0 + w/2).
    !       The unit of p0, p and w is typically wavenumber or frequency.
    !       The variable x is dimensionless and is zero at p = p0.
    !       
    sigma2      = sigma**2
    xx          = x**2 / sigma2
    florentz    = 1d0/(1d0 + xx)/pi

    return
  endfunction

  ! Returns the inverse of a matrix calculated by finding the LU decomposition.
  ! http://fortranwiki.org/fortran/show/Matrix+inversion
  function inv(A) result(Ainv)
    implicit none
    real*8, dimension(:,:), intent(in) :: A
    real*8, dimension(size(A,1),size(A,2)) :: Ainv
    real*8, dimension(size(A,1)) :: work ! work array for LAPACK
    integer*4, dimension(size(A,1)) :: ipiv ! pipov indices
    integer*4 :: n, info
    external DGETRF, DGETRI  ! External procedures defined in LAPACK
    
    Ainv = A
    n = size(A,1)
    
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)
    
    ! DGETRI computes the inverse of a matrix using the LU factorization
    call DGETRI(n, Ainv, n, ipiv, work, n, info)
  
    return  
  endfunction
  
endmodule
