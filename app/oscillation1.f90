program occilation1
  use mods
  implicit none
  integer :: i
  integer, parameter :: N = 2
  real(dp), dimension(N) :: x, w
  call gauss_chebyshev_1(x, w)
  print *, sum(w)
end program occilation1
