module routines
  implicit none


  integer, parameter :: rp = kind(1.0)
contains


  subroutine crist_to_cart( vec, lat, invlat )
    implicit none
    real(rp), intent(inout) :: vec
    real(rp), intent(in) :: lat(3,3)
    real(rp), intent(in) :: invlat(3,3)
  end subroutine crist_to_cart

  subroutine cart_to_crist( vec, lat, invlat )
    implicit none
    real(rp), intent(inout) :: vec
    real(rp), intent(in) :: lat(3,3)
    real(rp), intent(in) :: invlat(3,3)
  end subroutine cart_to_crist

  subroutine periodic( vec )
    implicit none
    real(rp), intent(inout) :: vec
  end subroutine periodic


  function frac_approx( x, h, k, kmax, nmin, nmax )result(ierr)
    implicit none
    real(rp), intent(in) :: x
    integer, intent(out) :: h
    integer, intent(out) :: k
    integer, intent(in), optional :: kmax
    integer, intent(in), optional :: nmin
    integer, intent(in), optional :: nmax
    integer :: ierr
    ! internal
    integer :: n, nmin_i, nmax_i
    integer, allocatable :: coeffs(:)
    integer :: kmax_i
    real(rp) :: q, r, xa
    real(rp) :: arow(3), brow(3), xrow(3)
    ierr = 0
    if( .not.present(kmax) .and. .not.present(nmax)) then
       write(*,*) "ERR: frac_approx needs one of: kmax or nmin"
       ierr = -1
       return
    end if

    kmax_i = huge(1)
    if(present(kmax))kmax_i = kmax
    nmin_i = 1
    if(present(nmin))nmin_i=nmin
    nmax_i=huge(1)
    if(present(nmax))nmax_i=nmax

    allocate(coeffs(0))
    xa = x
    arow = [1.0_rp, 0.0_rp, 1.0_rp]
    brow = [xa    , 1.0_rp, 0.0_rp]
    n = 0
    loop_: do while( n < nmax_i )
       q = floor(xa)
       r = xa - q
       if( r .le. 1e-6 ) exit loop_
       n = n + 1
       xa = 1.0/r
       xrow = [xa, arow(2)+q*brow(2), arow(3)+q*brow(3) ]
       if( n > nmin_i )then
          if( xrow(3) > real(kmax_i,kind=kind(xrow)) ) exit loop_
       end if
       arow = brow
       brow = xrow
       coeffs = [ coeffs, int(q) ]
    end do loop_
    h = int( brow(2) )
    k = int( brow(3) )

  end function frac_approx


  subroutine invmat3x3(mat,inv)
    implicit none

    real(rp), intent(in) :: mat(3,3)
    real(rp), intent(out) :: inv(3,3)
    !
    real(rp) :: det, invdet
    !
    det = 0.0_rp
    !
    ! calculate the determinant ...
    !
    det = det + mat(1,1)*mat(2,2)*mat(3,3) &
         + mat(1,2)*mat(2,3)*mat(3,1) &
         + mat(1,3)*mat(2,1)*mat(3,2) &
         - mat(1,3)*mat(2,2)*mat(3,1) &
         - mat(1,2)*mat(2,1)*mat(3,3) &
         - mat(1,1)*mat(2,3)*mat(3,2)
    ! invert the determinant
    invdet = 1/det
    ! calculate the inverse matrix
    inv(1,1) = invdet  * ( mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2) )
    inv(2,1) = -invdet * ( mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1) )
    inv(3,1) = invdet  * ( mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1) )
    inv(1,2) = -invdet * ( mat(1,2)*mat(3,3) - mat(1,3)*mat(3,2) )
    inv(2,2) = invdet  * ( mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1) )
    inv(3,2) = -invdet * ( mat(1,1)*mat(3,2) - mat(1,2)*mat(3,1) )
    inv(1,3) = invdet  * ( mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2) )
    inv(2,3) = -invdet * ( mat(1,1)*mat(2,3) - mat(1,3)*mat(2,1) )
    inv(3,3) = invdet  * ( mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1) )
    !
  end subroutine invmat3x3


end module routines
