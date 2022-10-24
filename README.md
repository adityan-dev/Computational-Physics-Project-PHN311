---
author:
- Adityan S
title: Group Project
---

# Oscillation {#sec:orgcc67047}

Consider the one-dimensional motion of a particle of mass $m$ in a time
independent potential $V(x)$.Since the energy $E$ is conserved, one can
integrate the equation of motion and obtain a solution in a closed form.

$$t-C = \sqrt{\frac{m}{2}}\int_{x_i}^x \frac{dx}{\sqrt{E-V(x)}}$$

For the choice $C=0$ , the particle is at the position $x_i$ at the time
$t = 0$ and $x$ refers to its position at any arbitrary time $t$.
Consider a particular case where the particle is in bound motion between
twopoints $a$ and $b$ where $V(x) = E$ for $x = a, b$ and $V(x) < E$ for
$a < x < b$. The time period of theoscillation $T$ is given by,

$$T = 2 \sqrt{\frac{m}{2}}\int_a^b \frac{dx}{\sqrt{E-V(x)}}$$

::: center
:::

## First consider a particle with $m = 1 kg$ in the potential $V_1(x) = 0.5 \alpha x^2$ with $\alpha =4 kg sec^2$. Numerically calculate the time period of oscillation and check this against the expected value. Verify that the frequency does not depend on the amplitude of oscillation. {#sec:org85e2e80}

::: center
:::

From the above plot, we infer that the for a value of Energy $E$, the
particle is bound for specific values of $a$ and $b$ (decided by $E$).

$$E = V_1(a)=V_1(b) = 2x^2$$

$$\implies a = - \sqrt{\frac{E}{2}} = -A, \quad b = +\sqrt{\frac{E}{2}} = +A$$

Hence Amplitude of oscillation is a function of Energy.

$$A(E) = b-a = \sqrt{\frac{E}{2}} \implies E = 2A^2$$

From the above equations, the analytic solution,

$$T = \frac{1}{A}\int_{-A}^{+A} \frac{dx}{\sqrt{1-(\frac{x}{A})^2}} = sin^{-1}{\Big(\frac{x}{A}\Big)} \Big|_{-A}^{ +A} = \pi \quad (sec)$$

Using integral identities, $T$ can be written as,

$$T = \int_{-1}^{+1} \frac{dx}{\sqrt{1-x^2}}$$

Let us solve the above numerically,

``` {.f90 breaklines="true" breakanywhere="true"}
program occilation1
  use mods
  implicit none
  integer :: i
  integer, parameter :: N = 2
  real(dp), dimension(N) :: x, w
  call gauss_chebyshev_1(x, w)
  print *, sum(w)
end program occilation1
```

``` {.sh breaklines="true" breakanywhere="true"}
fpm run oscillation1
```

    3.141592653589793

The above integral $T = 3.141592653589793 \approx \pi \neq T(E,A)$ .
Hence the timeperiod is independent of Amplitude.

## Next consider a potential $V_2(x) = e^{0.5αx^2 } -1$. Numerically verify that for small amplitude oscillations you recover the same results as the simple harmonic oscillator. {#sec:orgbea0e7d}

::: center
:::

From the above plot, we can infer that for particles of very small
energies the given potential behaves as the simple harmonic occilator.

Let us verify the above statement numerically,

$$T =\frac{\sqrt2}{e^{A^2}} \int _{-A}^{A} \frac{dx}{\sqrt{1-e^{2(x^2-A^2)}}} = \frac{\sqrt2A}{e^{A^2}} \int _{-1}^{1} \frac{dx}{\sqrt{1-e^{2A^2(x^2-1)}}}$$

$$\implies A = \pm \sqrt{\frac{\ln{(E+1)}}{2}}$$

``` {.f90 breaklines="true" breakanywhere="true"}
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
```

``` {.sh breaklines="true" breakanywhere="true"}
fpm run oscillation2
```

::: center
  ----- -----------
    $A$         $T$
          3.1410088
          3.1403022
          3.1391248
          3.1374772
            3.13536
          3.1327742
          3.1297208
          3.1262011
          3.1222167
          3.1177692
        
  ----- -----------
:::

::: center
:::

Hence for Oscillation of small amplitudes, we get the same results as
simple harmonic oscillator.

# Black Body Spectrum {#sec:org05272dd}

Quantum mechanics began with Planck's fit to the spectrum of black body
radiation:

$$I(\nu , T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{\frac{h\nu}{kT}}-1}$$

Here $I(\nu,T)$ is the energy per unit time of radiation with frequency
$\nu$ emitted per unit area of emitting surface, per unit solid angle,
and per unit frequency by a black body at temperature $T$. The parameter
$h$ is Planck's constant, $c$ is the speed of light in vacuum, and $k$
is Boltzmann constant. The Cosmic Background Explorer (COBE) project
measured the cosmic background radiation and obtained the results given
in the file `cobe.dat`

``` {.sh breaklines="true" breakanywhere="true"}
cat ./src/cobe/cobe.dat
```

::: center
  ------------------- ------------- -------
    $\widetilde{\nu}$   $I(\nu ,T)$   Error
                            200.723      14
                            249.508      19
                            293.024      25
                             327.77      23
                  ...           ...     ...
                              7.087      88
                              5.801     155
                              4.523     282
                                    
  ------------------- ------------- -------
:::

## Plot the COBE data and see if it has a shape similar to the black body spectrum first explained by Planck. {#sec:orgf695b34}

::: center
:::

Yes, the given data has a shape similar to the black body spectrum first
explained by Planck.

## Use these data to deduce the temperature T of the cosmic microwave background radiation. {#sec:org9847d67}

$$I =  \frac{0.0014745 \nu^3}{{e^{b\nu}}-1} \quad in \quad \frac{MJy}{sr}$$

$$\nu = 30\widetilde{\nu}  \quad in \quad GHz, \quad
\widetilde{\nu} \quad in \quad \frac{1}{cm^{-1}}$$
$$T = 10^9 \frac{h}{bk}, \quad b \quad in \frac{1}{GHz}, \quad T \quad in \quad K$$

Using Least Square Optimization

$$b = min\{ \sum_{i=1}^{43} r_i(b) \} = min\{ \sum_{i=1}^{43} ( \frac{0.0014745 \nu_i^3}{e^{b\nu_i}-1} - I_i)^2 \}$$

Let the Initial Guess for $b$ be $b = 0.048  \implies T \approx 1K$

``` {.python breaklines="true" breakanywhere="true"}
import numpy as np
from scipy.optimize import least_squares
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

data = pd.read_csv("src/cobe/cobeplot.csv", header=None, delim_whitespace=True)
vbi = np.asarray(data[0])
Ii = np.asarray(data[1])
vi = 30 * vbi

def r(b, vi, Ii):
    return ((0.0014745*pow(vi,3))/(np.exp(b*vi)-1) - Ii)**2

b = least_squares(r, 0.048, args=(vi, Ii))
b.x
```

\[0.01760962\]

$$b = 0.01760962  \implies T \approx 2.73K$$

## Assess the accuracy of the fit by doing a Pearson's $\chi^2$ (chi-squared) analysis. {#sec:org43ae4db}

``` {.f90 breaklines="true" breakanywhere="true"}
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
```

::: center
:::

::: center
  -------------------- -----------------
  Hypothesis           Plank's Law Fit
  Expected Values      $I(\nu_i)$
  Observes Values      $I_i$
  No of Observations   43
                       
  -------------------- -----------------
:::

$$\chi^2 = \sum_{i=1}^{43} \frac{(I_i - I(\nu_i))^2}{I(\nu_i)}$$

``` {.sh breaklines="true" breakanywhere="true"}
fpm run cobe
```

    0.04425666834533

$$\chi^2 = 0.04425666834533 < 1$$

Hence it is a good fit

# Modules {#sec:orge7d1fb8}

## Module `mods` consists of all subroutines, variables and procedures that are used in this report. The codeblocks below are all part of the file `mods.f90` {#sec:orge98a5fa}

``` {.f90 breaklines="true" breakanywhere="true"}
module mods
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use stdlib_quadrature, only: gauss_legendre
  implicit none

  real(dp), parameter :: pi=3.14159265358979323846_dp
  real(dp), parameter :: e= 2.7182818284590452353_dp
```

## Interfaces {#sec:org1b6f9f0}

``` {.f90 breaklines="true" breakanywhere="true"}
contains
```

## Routines for Array Manipulation {#sec:org77fe6d6}

### Linspace Subroutine for creating an array(sequence) with the equal step value. {#sec:org9333989}

``` {.f90 breaklines="true" breakanywhere="true"}
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
```

## Routines for Numerical Integration {#sec:org097f986}

### N-point Gauss Chebyshev Quadrature of 1st kind. {#sec:org18f1106}

$$\int _{-1}^{+1}{\frac {f(x)}{\sqrt {1-x^{2}}}}\,dx\approx \sum _{i=1}^{n}w_{i}f(x_{i})$$

$$x_{i}= \cos \left({\frac {2i-1}{2n}}\pi \right) \quad w_{i}=\frac {\pi }{n}$$

``` {.f90 breaklines="true" breakanywhere="true"}
subroutine gauss_chebyshev_1(x, w)
  real(dp), intent(out) :: w(:), x(:)
  integer :: n, i
  n = size(w)
  do i=1, n
     x(i) = cos((pi*(2*real(i)-1))/(2*real(n)))
     w(i) = pi/real(n)
  end do
end subroutine gauss_chebyshev_1
```

### N-point Gauss Legendre Quadrature of 1st kind. {#sec:org62b5f27}

$$\int _{-1}^{+1}{f(x)}\,dx\approx \sum _{i=1}^{n}w_{i}f(x_{i})$$

$$w_{i}=\frac {2}{\left(1-x_{i}^{2}\right)\left[P'_{n}(x_{i})\right]^{2}} \quad   P_n(x_i)=0$$

``` {.f90 breaklines="true" breakanywhere="true"}
!use stdlib_quadrature, only: gauss_legendre
```

## Routines for Optimization of Non Linear Equation {#sec:org58523f3}

Using SCIPY

## End of module `mods` {#sec:org48b10e5}

``` {.f90 breaklines="true" breakanywhere="true"}
end module mods
```
