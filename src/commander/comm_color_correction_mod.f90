! Numerical kernel for foreground-based colour corrections, independent of MPI.
module comm_color_correction_mod
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: color_correction_factors

contains

  pure subroutine color_correction_factors(sky_low, sky_high, nu_low, nu_high, &
       & coeffs, min_signal, factors, alpha, valid)
    ! sky_* contains summed foreground I,Q,U in uK_RJ, with CMB excluded.
    ! coeffs(:,stokes) = [A,B,C] for the divisive correction K(alpha).
    ! Multiply the monochromatic forward model by K; do not divide it.
    real(real64), intent(in) :: sky_low(3), sky_high(3), nu_low, nu_high
    real(real64), intent(in) :: coeffs(3,3), min_signal
    real(real64), intent(out) :: factors(3), alpha(3)
    logical, intent(out) :: valid(3)
    real(real64) :: low(3), high(3), log_ratio, value
    integer :: s

    factors = 1.0_real64
    alpha = 0.0_real64       ! Not an estimate where valid is false.
    valid = .false.
    if (.not. ieee_is_finite(nu_low) .or. .not. ieee_is_finite(nu_high)) return
    if (nu_low <= 0.0_real64 .or. nu_high <= nu_low) return
    if (.not. ieee_is_finite(min_signal) .or. min_signal < 0.0_real64) return
    log_ratio = log(nu_high) - log(nu_low)
    if (log_ratio <= 0.0_real64) return

    low(1) = sky_low(1)
    high(1) = sky_high(1)
    low(2:3) = polarized_intensity(sky_low(2), sky_low(3))
    high(2:3) = polarized_intensity(sky_high(2), sky_high(3))
    do s = 1, 3
       if (.not. ieee_is_finite(low(s)) .or. .not. ieee_is_finite(high(s))) cycle
       if (low(s) <= min_signal .or. high(s) <= min_signal) cycle
       alpha(s) = 2.0_real64 + (log(high(s)) - log(low(s))) / log_ratio
       if (.not. ieee_is_finite(alpha(s))) cycle
       if (.not. all(ieee_is_finite(coeffs(:,s)))) cycle
       value = coeffs(1,s) + alpha(s) * (coeffs(2,s) + alpha(s) * coeffs(3,s))
       if (.not. ieee_is_finite(value) .or. value <= 0.0_real64) cycle
       factors(s) = value
       valid(s) = .true.
    end do
  end subroutine color_correction_factors

  pure function polarized_intensity(q, u) result(p)
    real(real64), intent(in) :: q, u
    real(real64) :: p, scale
    ! Scale before squaring to avoid underflow/overflow of individual squares.
    if (.not. ieee_is_finite(q)) then
       p = q
    else if (.not. ieee_is_finite(u)) then
       p = u
    else
       scale = max(abs(q), abs(u))
       p = 0.0_real64
       if (scale > 0.0_real64) p = scale * sqrt((q/scale)**2 + (u/scale)**2)
    end if
  end function polarized_intensity

end module comm_color_correction_mod
