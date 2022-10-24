program occilation2
  use mods
  implicit none
  integer :: i
  integer, parameter :: N = 5000
  real(dp), dimension(100) :: A
  real(dp), dimension(N) :: x, w
  call gauss_legendre(x, w)
  call linspace(0.01_dp, 1.0_dp, A)
  open(1, file='./src/oscillation/o2T.csv', status='old')
  do i=1,100
     write (1,'(1f10.7, a, 1f10.7)') A(i), " ", sum((sqrt(2._dp)*A(i)/(e**(A(i)**2)))*w*(1/sqrt(1-e**(2*(A(i)**2)*((x**2)-1)))))
  end do
  write (*, *) "$A$ $T$"
  do i=1,10
     write (*,'(1f10.7, a, 1f10.7)') A(i), " ", sum((sqrt(2._dp)*A(i)/(e**(A(i)**2)))*w*(1/sqrt(1-e**(2*(A(i)**2)*((x**2)-1)))))
  end do
end program occilation2
