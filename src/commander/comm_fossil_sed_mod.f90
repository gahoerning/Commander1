! Physical FOSSIL monopole spectra, evaluated directly in MJy/sr.
! This module deliberately does not use Commander's reference-frequency units.
module comm_fossil_sed_mod
  use healpix_types, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private

  integer, parameter, public :: FOSSIL_SED_OK = 0, FOSSIL_SED_INVALID_KIND = 1
  integer, parameter, public :: FOSSIL_SED_INVALID_PARAMETER = 2, FOSSIL_SED_BAD_TABLE = 3
  integer, parameter, public :: FOSSIL_SED_OUT_OF_RANGE = 4, FOSSIL_SED_IO_ERROR = 5
  ! SI constants and T0 match fossil_simulation.py (SciPy CODATA 2022).
  real(dp), parameter :: fossil_h = 6.62607015e-34_dp, fossil_k = 1.380649e-23_dp
  real(dp), parameter :: fossil_c = 299792458.0_dp, fossil_me = 9.1093837139e-31_dp
  real(dp), parameter :: fossil_ev = 1.602176634e-19_dp
  real(dp), parameter, public :: fossil_T0 = 2.725_dp
  real(dp), parameter :: fossil_I0 = 2.0_dp*(fossil_k*fossil_T0)**3/(fossil_h*fossil_c)**2*1.e20_dp

  type, public :: fossil_sed_table
     real(dp), allocatable :: nu_GHz(:), intensity_MJy_sr(:)
     logical :: logarithmic = .false.
     logical :: extrapolate = .false.
  end type fossil_sed_table

  public :: fossil_sed_value, fossil_load_sed_table, fossil_itoh_kernels

contains

  function fossil_sed_value(kind, nu_GHz, shape, table, status) result(value)
    ! Multiply this response by one native physical amplitude:
    ! delta_T = Delta T/T0; mu; tsz = Compton y; crr/co = dimensionless A;
    ! cib = coefficient A in MJy/sr, with SED x**(3+beta)/(exp(x)-1).
    ! tsz shape=[kTe_keV], cib shape=[beta,T_K]; others have shape of length 0.
    character(len=*), intent(in) :: kind
    real(dp), intent(in) :: nu_GHz, shape(:)
    type(fossil_sed_table), intent(in), optional :: table
    integer, intent(out), optional :: status
    real(dp) :: value, x, G, theta, kernel(0:4), denominator
    integer :: code, i

    code = FOSSIL_SED_OK
    value = ieee_value(0.0_dp, ieee_quiet_nan)
    if (.not. ieee_is_finite(nu_GHz) .or. nu_GHz <= 0.0_dp) then
       code = FOSSIL_SED_INVALID_PARAMETER
    else if (.not. all(ieee_is_finite(shape))) then
       code = FOSSIL_SED_INVALID_PARAMETER
    else
       select case (trim(fossil_lower(kind)))
       case ('delta_t', 'mu', 'tsz')
          if ((trim(fossil_lower(kind)) == 'tsz' .and. size(shape) /= 1) .or. &
              (trim(fossil_lower(kind)) /= 'tsz' .and. size(shape) /= 0)) then
             code = FOSSIL_SED_INVALID_PARAMETER
          else if (trim(fossil_lower(kind)) == 'tsz') then
             if (shape(1) < 0.0_dp) code = FOSSIL_SED_INVALID_PARAMETER
          end if
          if (code == FOSSIL_SED_OK) then
             x = fossil_h*(nu_GHz*1.e9_dp)/(fossil_k*fossil_T0)
             ! The FOSSIL range is 50--2000 GHz. This limit also avoids
             ! evaluating very high-order polynomials after the Wien tail vanishes.
             if (x > 700.0_dp) then
                value = 0.0_dp
             else
                denominator = -fossil_expm1(-x)
                G = x**4*exp(-x)/denominator**2
                select case (trim(fossil_lower(kind)))
                case ('delta_t')
                   value = fossil_I0*G
                case ('mu')
                   value = fossil_I0*G*(1.0_dp/2.1923_dp - 1.0_dp/x)
                case ('tsz')
                   call fossil_itoh_kernels(x, kernel)
                   theta = shape(1)*1.e3_dp*fossil_ev/(fossil_me*fossil_c*fossil_c)
                   value = kernel(4)
                   do i = 3, 0, -1
                      value = kernel(i) + theta*value
                   end do
                   value = fossil_I0*G*value
                end select
             end if
          end if
       case ('cib')
          if (size(shape) /= 2) then
             code = FOSSIL_SED_INVALID_PARAMETER
          else if (shape(2) <= 0.0_dp) then
             code = FOSSIL_SED_INVALID_PARAMETER
          else
             x = fossil_h*(nu_GHz*1.e9_dp)/(fossil_k*shape(2))
             ! Evaluate in log space to keep high-frequency tails finite.
             if (x > 50.0_dp) then
                value = exp((3.0_dp+shape(1))*log(x) - x - log(1.0_dp-exp(-x)))
             else
                value = x**(3.0_dp+shape(1))/fossil_expm1(x)
             end if
          end if
       case ('crr', 'co')
          if (size(shape) /= 0) then
             code = FOSSIL_SED_INVALID_PARAMETER
          else if (.not. present(table)) then
             code = FOSSIL_SED_BAD_TABLE
          else
             call fossil_interpolate(nu_GHz, table, value, code)
          end if
       case default
          code = FOSSIL_SED_INVALID_KIND
       end select
    end if
    if (code == FOSSIL_SED_OK) then
       if (.not. ieee_is_finite(value)) code = FOSSIL_SED_INVALID_PARAMETER
    end if
    if (code /= FOSSIL_SED_OK) value = ieee_value(0.0_dp, ieee_quiet_nan)
    if (present(status)) status = code
  end function fossil_sed_value

  subroutine fossil_itoh_kernels(x, kernel)
    ! Itoh et al. (1998), equations 2.26--2.30, reproduced from the simulation.
    ! Here Compton y includes theta_e, rather than the paper's optical-depth y.
    real(dp), intent(in) :: x
    real(dp), intent(out) :: kernel(0:4)
    real(dp) :: Xbar, S, z, denominator
    z = exp(-x)
    denominator = -fossil_expm1(-x)
    Xbar = x*(1.0_dp+z)/denominator
    S = 2.0_dp*x*exp(-x/2.0_dp)/denominator
    kernel(0) = Xbar - 4.0_dp
    kernel(1) = -10.0_dp + 47.0_dp/2*Xbar - 42.0_dp/5*Xbar**2 + 7.0_dp/10*Xbar**3 &
         + S**2*(-21.0_dp/5 + 7.0_dp/5*Xbar)
    kernel(2) = -15.0_dp/2 + 1023.0_dp/8*Xbar - 868.0_dp/5*Xbar**2 + 329.0_dp/5*Xbar**3 &
         - 44.0_dp/5*Xbar**4 + 11.0_dp/30*Xbar**5 &
         + S**2*(-434.0_dp/5 + 658.0_dp/5*Xbar - 242.0_dp/5*Xbar**2 + 143.0_dp/30*Xbar**3) &
         + S**4*(-44.0_dp/5 + 187.0_dp/60*Xbar)
    kernel(3) = 15.0_dp/2 + 2505.0_dp/8*Xbar - 7098.0_dp/5*Xbar**2 + 14253.0_dp/10*Xbar**3 &
         - 18594.0_dp/35*Xbar**4 + 12059.0_dp/140*Xbar**5 - 128.0_dp/21*Xbar**6 &
         + 16.0_dp/105*Xbar**7 &
         + S**2*(-7098.0_dp/10 + 14253.0_dp/5*Xbar - 102267.0_dp/35*Xbar**2 &
         + 156767.0_dp/140*Xbar**3 - 1216.0_dp/7*Xbar**4 + 64.0_dp/7*Xbar**5) &
         + S**4*(-18594.0_dp/35 + 205003.0_dp/280*Xbar - 1920.0_dp/7*Xbar**2 + 1024.0_dp/35*Xbar**3) &
         + S**6*(-544.0_dp/21 + 992.0_dp/105*Xbar)
    kernel(4) = -135.0_dp/32 + 30375.0_dp/128*Xbar - 62391.0_dp/10*Xbar**2 &
         + 614727.0_dp/40*Xbar**3 - 124389.0_dp/10*Xbar**4 + 355703.0_dp/80*Xbar**5 &
         - 16568.0_dp/21*Xbar**6 + 7516.0_dp/105*Xbar**7 - 22.0_dp/7*Xbar**8 + 11.0_dp/210*Xbar**9 &
         + S**2*(-62391.0_dp/20 + 614727.0_dp/20*Xbar - 1368279.0_dp/20*Xbar**2 &
         + 4624139.0_dp/80*Xbar**3 - 157396.0_dp/7*Xbar**4 + 30064.0_dp/7*Xbar**5 &
         - 2717.0_dp/7*Xbar**6 + 2761.0_dp/210*Xbar**7) &
         + S**4*(-124389.0_dp/10 + 6046951.0_dp/160*Xbar - 248520.0_dp/7*Xbar**2 &
         + 481024.0_dp/35*Xbar**3 - 15972.0_dp/7*Xbar**4 + 18689.0_dp/140*Xbar**5) &
         + S**6*(-70414.0_dp/21 + 465992.0_dp/105*Xbar - 11792.0_dp/7*Xbar**2 + 19778.0_dp/105*Xbar**3) &
         + S**8*(-682.0_dp/7 + 7601.0_dp/210*Xbar)
  end subroutine fossil_itoh_kernels

  subroutine fossil_load_sed_table(filename, table, status, logarithmic, extrapolate)
    ! Two-column text: frequency [GHz], intensity [MJy/sr]. Blank lines and
    ! comments beginning with # or ! are accepted. Input must be ordered.
    character(len=*), intent(in) :: filename
    type(fossil_sed_table), intent(out) :: table
    integer, intent(out) :: status
    logical, intent(in), optional :: logarithmic, extrapolate
    character(len=4096) :: line
    real(dp) :: frequency, intensity
    integer :: unit, ios, pass, n, i

    status = FOSSIL_SED_OK
    if (present(logarithmic)) table%logarithmic = logarithmic
    if (present(extrapolate)) table%extrapolate = extrapolate
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
       status = FOSSIL_SED_IO_ERROR
       return
    end if
    n = 0
    do pass = 1, 2
       i = 0
       do
          read(unit, '(A)', iostat=ios) line
          if (ios < 0) exit
          if (ios > 0) then
             status = FOSSIL_SED_IO_ERROR
             exit
          end if
          line = adjustl(line)
          if (len_trim(line) == 0) cycle
          if (line(1:1) == '#' .or. line(1:1) == '!') cycle
          read(line, *, iostat=ios) frequency, intensity
          if (ios /= 0) then
             status = FOSSIL_SED_BAD_TABLE
             exit
          end if
          if (.not. ieee_is_finite(frequency) .or. .not. ieee_is_finite(intensity)) then
             status = FOSSIL_SED_BAD_TABLE
             exit
          end if
          if (frequency <= 0.0_dp .or. (table%logarithmic .and. intensity <= 0.0_dp)) then
             status = FOSSIL_SED_BAD_TABLE
             exit
          end if
          i = i+1
          if (pass == 2) then
             ! Guard against a file being changed between the two passes.
             if (i > n) then
                status = FOSSIL_SED_BAD_TABLE
                exit
             end if
             table%nu_GHz(i) = frequency
             table%intensity_MJy_sr(i) = intensity
             if (i > 1) then
                if (frequency <= table%nu_GHz(i-1)) then
                   status = FOSSIL_SED_BAD_TABLE
                   exit
                end if
             end if
          end if
       end do
       if (status /= FOSSIL_SED_OK) exit
       if (pass == 1) then
          n = i
          if (n < 2) then
             status = FOSSIL_SED_BAD_TABLE
             exit
          end if
          allocate(table%nu_GHz(n), table%intensity_MJy_sr(n))
          rewind(unit, iostat=ios)
          if (ios /= 0) then
             status = FOSSIL_SED_IO_ERROR
             exit
          end if
       else if (i /= n) then
          status = FOSSIL_SED_BAD_TABLE
       end if
    end do
    close(unit)
    if (status /= FOSSIL_SED_OK) then
       if (allocated(table%nu_GHz)) deallocate(table%nu_GHz)
       if (allocated(table%intensity_MJy_sr)) deallocate(table%intensity_MJy_sr)
    end if
  end subroutine fossil_load_sed_table

  subroutine fossil_interpolate(nu_GHz, table, value, status)
    real(dp), intent(in) :: nu_GHz
    type(fossil_sed_table), intent(in) :: table
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    integer :: n, lo, hi, middle
    real(dp) :: weight

    value = ieee_value(0.0_dp, ieee_quiet_nan)
    status = FOSSIL_SED_BAD_TABLE
    if (.not. allocated(table%nu_GHz) .or. .not. allocated(table%intensity_MJy_sr)) return
    n = size(table%nu_GHz)
    if (n < 2 .or. size(table%intensity_MJy_sr) /= n) return
    ! Tables are validated once on load; evaluation uses a binary search.
    if (.not. table%extrapolate) then
       if (nu_GHz < table%nu_GHz(1) .or. nu_GHz > table%nu_GHz(n)) then
          status = FOSSIL_SED_OUT_OF_RANGE
          return
       end if
    end if
    lo = 1
    hi = n
    do while (hi-lo > 1)
       middle = (lo+hi)/2
       if (nu_GHz >= table%nu_GHz(middle)) then
          lo = middle
       else
          hi = middle
       end if
    end do
    if (table%logarithmic) then
       weight = log(nu_GHz/table%nu_GHz(lo))/log(table%nu_GHz(hi)/table%nu_GHz(lo))
       value = exp(log(table%intensity_MJy_sr(lo)) &
            + weight*log(table%intensity_MJy_sr(hi)/table%intensity_MJy_sr(lo)))
    else
       weight = (nu_GHz-table%nu_GHz(lo))/(table%nu_GHz(hi)-table%nu_GHz(lo))
       value = table%intensity_MJy_sr(lo) + weight*(table%intensity_MJy_sr(hi)-table%intensity_MJy_sr(lo))
    end if
    status = FOSSIL_SED_OK
  end subroutine fossil_interpolate

  pure function fossil_expm1(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    if (abs(x) < 1.e-5_dp) then
       value = x*(1.0_dp + x*(0.5_dp + x*(1.0_dp/6.0_dp + x*(1.0_dp/24.0_dp + x/120.0_dp))))
    else
       value = exp(x)-1.0_dp
    end if
  end function fossil_expm1

  pure function fossil_lower(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code
    lower = text
    do i = 1, len(text)
       code = iachar(text(i:i))
       if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code+32)
    end do
  end function fossil_lower

end module comm_fossil_sed_mod
