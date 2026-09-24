program test_color_corrections
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use comm_color_correction_mod
  implicit none
  real(real64) :: lo(3), hi(3), coeff(3,3), k(3), alpha(3), expected, pi, theta
  logical :: valid(3)

  coeff(:,1) = [0.99970_real64, -0.00052_real64, 0.00127_real64]
  coeff(:,2) = [0.99776_real64, -0.00716_real64, 0.00090_real64]
  coeff(:,3) = [1.00465_real64, 0.01650_real64, 0.00240_real64]

  ! A pure RJ power law beta=-3 has flux-density index alpha=-1.
  lo = [100.0_real64, 3.0_real64, -4.0_real64]
  hi = lo * 2.0_real64**(-3.0_real64)
  call evaluate()
  call require(all(valid), 'power-law slopes are valid')
  call near(alpha(1), -1.0_real64, 'RJ-to-flux conversion')
  call near(alpha(2), -1.0_real64, 'polarized intensity slope')
  call near(alpha(3), alpha(2), 'Q and U share alpha_P')
  call near(k(2), 1.00582_real64, 'Q polynomial at alpha=-1')
  call near(k(3), 0.99055_real64, 'U polynomial at alpha=-1')

  ! The sum of two power laws has a different slope from either component.
  lo(1) = 100.0_real64 + 80.0_real64
  hi(1) = 100.0_real64/8.0_real64 + 80.0_real64/4.0_real64
  call evaluate()
  expected = 2.0_real64 + log(32.5_real64/180.0_real64)/log(2.0_real64)
  call near(alpha(1), expected, 'total foreground mixture')
  call require(abs(alpha(1)+1.0_real64) > 0.1_real64, 'not the synchrotron index')

  ! Rotating Q,U across the band, including a Q sign flip, preserves P's slope.
  pi = acos(-1.0_real64)
  theta = 0.8_real64*pi
  lo(2:3) = [3.0_real64, 4.0_real64]
  hi(2) = (lo(2)*cos(theta)-lo(3)*sin(theta))/8.0_real64
  hi(3) = (lo(2)*sin(theta)+lo(3)*cos(theta))/8.0_real64
  call evaluate()
  call require(hi(2) < 0.0_real64, 'test exercises a Stokes sign flip')
  call near(alpha(2), -1.0_real64, 'P is invariant under rotation')

  ! Destructive interference: sum Q,U first, not individual component P.
  lo(2:3) = [10.0_real64-9.0_real64, 0.0_real64]
  hi(2:3) = [10.0_real64/8.0_real64-9.0_real64/4.0_real64, 0.0_real64]
  call evaluate()
  call near(alpha(2), 2.0_real64, 'polarization cancellation in the mixture')

  ! Missing/negative I does not invalidate a well-defined polarized spectrum.
  lo(1) = -1.0_real64
  call evaluate()
  call require(.not. valid(1) .and. all(valid(2:3)), 'independent I/P validity')
  call near(k(1), 1.0_real64, 'nonpositive I falls back to unity')
  lo(2:3) = 0.0_real64
  call evaluate()
  call require(.not. any(valid), 'zero P is undefined')
  call require(all(k == 1.0_real64), 'undefined slopes do not create NaNs')

  lo = 1.0_real64
  hi = 0.125_real64
  lo(2) = ieee_value(0.0_real64, ieee_quiet_nan)
  call evaluate()
  call require(valid(1) .and. .not. any(valid(2:3)), 'non-finite Q rejected')
  lo = 1.e200_real64
  hi = lo/8.0_real64
  call evaluate()
  call require(all(valid), 'hypot avoids overflow of Q squared')
  call near(alpha(2), -1.0_real64, 'large-amplitude polarized slope')

  call color_correction_factors(lo, hi, 1.0_real64, 1.0_real64, coeff, &
       & 1.e-12_real64, k, alpha, valid)
  call require(.not. any(valid), 'zero-width band is invalid')
  lo = 1.e-14_real64
  hi = lo/8.0_real64
  call evaluate()
  call require(.not. any(valid), 'signal floor is respected')
  lo = 1.0_real64
  hi = lo/8.0_real64
  coeff(:,2) = [-1.0_real64, 0.0_real64, 0.0_real64]
  call evaluate()
  call require(.not. valid(2) .and. valid(3), 'nonpositive polynomial rejected per Stokes')
  call near(k(2), 1.0_real64, 'invalid polynomial uses unity')
  print *, 'PASS: color correction numerical tests'

contains
  subroutine evaluate()
    call color_correction_factors(lo, hi, 4.e9_real64, 8.e9_real64, coeff, &
         & 1.e-12_real64, k, alpha, valid)
  end subroutine evaluate

  subroutine near(actual, target, label)
    real(real64), intent(in) :: actual, target
    character(len=*), intent(in) :: label
    call require(abs(actual-target) < 1.e-11_real64, label)
  end subroutine near

  subroutine require(ok, label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    if (.not. ok) then
       print *, 'FAIL: ', label
       error stop 1
    end if
  end subroutine require
end program test_color_corrections
