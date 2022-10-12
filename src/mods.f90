module mods
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use stdlib_quadrature, only: gauss_legendre !g6
  implicit none
  contains

subroutine linspace(from, to, array)
  real(dp), intent(in) :: from, to
  real(dp), intent(out) :: array(:)
  real(dp) :: range
  integer :: n, i
  n = size(array)
  range = to - from
  if (n == 0) return
  if (n == 1) then
     array(1) = from
     return
  end if
  do i=1, n
     array(i) = from + range * (i - 1) / (n - 1)
  end do
end subroutine linspace

!g6 working method

!integer, parameter :: N = 6
!real(dp), dimension(N) :: x, w
!call gauss_legendre(x, w)
!integ = sum(w*f(x))

end module mods
