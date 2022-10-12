program occilation1
  use mods
  real(dp) :: A(200), T(200), xi(6), wi(6), f, x
  f(x) = x**2
  call gauss_legendre(xi, wi)
  integ = sum(wi*f(xi))
  call linspace(-10._dp, 10._dp, A)
  print *, integ
end program occilation1
