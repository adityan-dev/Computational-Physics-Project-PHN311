program cobe
  use mods
  implicit none

  real(dp) :: v(200), I(200), vi(43), vbi(43), Ii(43), erri(43), chisq
  integer :: m, j

  m = 43 ! No of data points
  open(99, file="./src/cobe/cobe.dat")
  open(98, file="./src/cobe/cobefit.csv")
  open(97, file="./src/cobe/cobeplotfreq.csv")

  do j=1, m
     read(99, *) vbi(j), Ii(j), erri(j)
  end do
  vi = 30._dp*vbi
  do j=1, m
     write(97, *) vi(j), Ii(j)
  end do
  call linspace(0._dp, 1000._dp, v)
  do j=1, 200
     I(j) = blackbody(v(j))
     write(98, "(2f20.14)") v(j), I(j)
  end do

  chisq = 0._dp
  do j=1,m
     chisq = chisq + (((Ii(j) - blackbody(vi(j)))**2._dp)/(blackbody(vi(j))))
  end do
  write(*, "(1f20.14)") chisq
  
  close(99)
  close(98)
  close(97)
contains
  real(dp) function blackbody(v) result(I)
    real(dp), intent(in) :: v
    I = (0.0014745_dp*(v**3._dp))/(exp(v*0.01760962_dp)-1._dp)
  end function blackbody
end program cobe
