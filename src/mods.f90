module mods
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use stdlib_quadrature, only: gauss_legendre
  implicit none

  real(dp), parameter :: pi=3.14159265358979323846_dp
  real(dp), parameter :: e= 2.7182818284590452353_dp

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

subroutine gauss_chebyshev_1(x, w)
  real(dp), intent(out) :: w(:), x(:)
  integer :: n, i
  n = size(w)
  do i=1, n
     x(i) = cos((pi*(2*real(i)-1))/(2*real(n)))
     w(i) = pi/real(n)
  end do
end subroutine gauss_chebyshev_1

!use stdlib_quadrature, only: gauss_legendre

end module mods
