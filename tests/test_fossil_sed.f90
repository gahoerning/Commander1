program test_fossil_sed
  use healpix_types, only: dp
  use comm_fossil_sed_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(fossil_sed_table) :: crr, co, empty, malformed
  character(len=4096) :: directory, line
  real(dp) :: nu, reference(6), value, base, kernel(0:4), x
  integer :: status, unit, ios, n, i
  real(dp), parameter :: h = 6.62607015e-34_dp, k = 1.380649e-23_dp

  call get_command_argument(1, directory)
  if (len_trim(directory) == 0) error stop 'Supply the commander_spectral directory'
  call fossil_load_sed_table(trim(directory)//'/crr_MJy_sr.txt', crr, status)
  call require(status == FOSSIL_SED_OK, 'CRR load')
  call fossil_load_sed_table(trim(directory)//'/co_MJy_sr.txt', co, status, &
       logarithmic=.true., extrapolate=.true.)
  call require(status == FOSSIL_SED_OK, 'CO load')
  open(newunit=unit, file=trim(directory)//'/fiducial_spectra_MJy_sr.txt', status='old', action='read')
  n = 0
  do
     read(unit, '(A)', iostat=ios) line
     if (ios < 0) exit
     call require(ios == 0, 'reference file read')
     line = adjustl(line)
     if (len_trim(line) == 0) cycle
     if (line(1:1) == '#') cycle
     read(line, *, iostat=ios) nu, reference
     call require(ios == 0, 'reference row parse')
     value = fossil_sed_value('delta_T', nu, [real(dp)::], status=status)
     call close_to(value, reference(1), 'delta_T')
     call require(status == FOSSIL_SED_OK, 'delta_T status')
     value = fossil_sed_value('mu', nu, [real(dp)::], status=status)
     call close_to(value, reference(2), 'mu')
     value = fossil_sed_value('tsz', nu, [1.245_dp], status=status)
     call close_to(value, reference(3), 'tsz')
     value = fossil_sed_value('cib', nu, [0.86_dp,18.8_dp], status=status)
     call close_to(value, reference(4), 'cib')
     value = fossil_sed_value('crr', nu, [real(dp)::], table=crr, status=status)
     call close_to(value, reference(5), 'crr')
     value = fossil_sed_value('co', nu, [real(dp)::], table=co, status=status)
     call close_to(value, reference(6), 'co')
     base = fossil_sed_value('delta_T', nu, [real(dp)::])
     value = fossil_sed_value('tsz', nu, [0.0_dp])
     x = h*nu*1.e9_dp/(k*fossil_T0)
     call fossil_itoh_kernels(x, kernel)
     call close_to(value, base*kernel(0), 'zero-temperature nonrelativistic limit')
     n = n+1
  end do
  close(unit)
  call require(n == 131, '131 reference channels')

  value = fossil_sed_value('crr', crr%nu_GHz(1), [real(dp)::], table=crr, status=status)
  call close_to(value, crr%intensity_MJy_sr(1), 'CRR lower endpoint')
  value = fossil_sed_value('crr', crr%nu_GHz(size(crr%nu_GHz)), [real(dp)::], table=crr, status=status)
  call close_to(value, crr%intensity_MJy_sr(size(crr%nu_GHz)), 'CRR upper endpoint')
  value = fossil_sed_value('crr', 0.5_dp*crr%nu_GHz(1), [real(dp)::], table=crr, status=status)
  call require(status == FOSSIL_SED_OUT_OF_RANGE .and. .not. ieee_is_finite(value), 'CRR strict lower bound')
  value = fossil_sed_value('crr', 2.0_dp*crr%nu_GHz(size(crr%nu_GHz)), [real(dp)::], table=crr, status=status)
  call require(status == FOSSIL_SED_OUT_OF_RANGE, 'CRR strict upper bound')
  value = fossil_sed_value('co', 2000.0_dp, [real(dp)::], table=co, status=status)
  call require(status == FOSSIL_SED_OK .and. value > 0.0_dp, 'CO extrapolation')
  value = fossil_sed_value('co', 100.0_dp, [real(dp)::], table=empty, status=status)
  call require(status == FOSSIL_SED_BAD_TABLE, 'unallocated template rejected')
  value = fossil_sed_value('co', 100.0_dp, [real(dp)::], status=status)
  call require(status == FOSSIL_SED_BAD_TABLE, 'missing template rejected')
  value = fossil_sed_value('unknown', 100.0_dp, [real(dp)::], status=status)
  call require(status == FOSSIL_SED_INVALID_KIND, 'unknown kind rejected')
  value = fossil_sed_value('delta_T', -100.0_dp, [real(dp)::], status=status)
  call require(status == FOSSIL_SED_INVALID_PARAMETER, 'negative frequency rejected')
  value = fossil_sed_value('tsz', 100.0_dp, [-1.0_dp], status=status)
  call require(status == FOSSIL_SED_INVALID_PARAMETER, 'negative temperature rejected')
  value = fossil_sed_value('tsz', 100.0_dp, [real(dp)::], status=status)
  call require(status == FOSSIL_SED_INVALID_PARAMETER, 'missing shape rejected')
  value = fossil_sed_value('cib', 100.0_dp, [0.86_dp,0.0_dp], status=status)
  call require(status == FOSSIL_SED_INVALID_PARAMETER, 'zero CIB temperature rejected')

  ! Exercise the real table reader, including signed linear values and malformed inputs.
  do i = 1, 5
     open(newunit=unit, file='temporary_sed_table.txt', status='replace', action='write')
     select case(i)
     case(1)
        write(unit,'(A)') '# signed values are valid for a linear spectrum'
        write(unit,'(A)') '1 -2'
        write(unit,'(A)') '2 2'
     case(2)
        write(unit,'(A)') '1 1'
        write(unit,'(A)') '1 2'
     case(3)
        write(unit,'(A)') '1 1'
     case(4)
        write(unit,'(A)') '1 0'
        write(unit,'(A)') '2 1'
     case(5)
        write(unit,'(A)') 'this is not a data row'
     end select
     close(unit)
     call fossil_load_sed_table('temporary_sed_table.txt', malformed, status, logarithmic=(i==4))
     if (i == 1) then
        call require(status == FOSSIL_SED_OK, 'signed table load')
        value = fossil_sed_value('crr', 1.5_dp, [real(dp)::], table=malformed, status=status)
        call close_to(value, 0.0_dp, 'signed linear interpolation')
     else
        call require(status == FOSSIL_SED_BAD_TABLE, 'malformed table rejection')
     end if
  end do
  call fossil_load_sed_table('nonexistent_sed_table.txt', malformed, status)
  call require(status == FOSSIL_SED_IO_ERROR, 'missing file rejection')
  print *, 'PASS: production FOSSIL SEDs, 131 channels, physical limits and table/error policies'

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
       print *, 'FAIL: ', message
       error stop 1
    end if
  end subroutine require

  subroutine close_to(actual, expected, message)
    real(dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: message
    if (.not. ieee_is_finite(actual)) then
       print *, 'FAIL: nonfinite ', message
       error stop 1
    end if
    if (abs(actual-expected) > 1.e-10_dp*max(abs(expected),1.e-14_dp)) then
       print *, 'FAIL: ', message, ' actual/expected = ', actual, expected
       error stop 1
    end if
  end subroutine close_to

end program test_fossil_sed
