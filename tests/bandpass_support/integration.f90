program test_bandpass
  use comm_bp_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  type(planck_rng) :: handle
  character(len=1024) :: config, mode
  real(dp), allocatable :: x(:), s(:), t(:), cmb(:)
  real(dp) :: expected, original, oldref, delta(1)
  integer :: i

  call get_command_argument(1,config)
  call get_command_argument(2,mode)
  call initialize_bp_mod(0,1,0,trim(config),handle)
  if (trim(mode) /= 'valid') then
     ! Invalid configuration modes fail above; malformed responses fail here.
     select case(trim(mode))
     case('negative_gain')
        bp(1)%tau0(2)=-1.d0
     case('zero_gain')
        bp(1)%tau0=0.d0
     case('nan_gain')
        bp(1)%tau0(2)=ieee_value(0.d0,ieee_quiet_nan)
     case('nan_frequency')
        bp(1)%nu0(2)=ieee_value(0.d0,ieee_quiet_nan)
     case('descending_frequency')
        bp(1)%nu0(2)=bp(1)%nu0(1)-1.d0
     case('zero_frequency')
        bp(1)%nu0(1)=0.d0
     case('one_sample')
        bp(1)%n=1
     case('overflow_calibration')
        bp(1)%cal_index=1.d9
     end select
     delta=0.d0
     call update_tau(delta)
     error stop 'Invalid case was accepted'
  end if

  allocate(x(bp(1)%n),s(bp(1)%n),t(bp(1)%n),cmb(bp(1)%n))
  x=bp(1)%nu/bp(1)%nu_c
  call near(get_bp_avg_spectrum(1,x**(bp(1)%cal_index-2.d0)),1.d0,2.d-14,'calibrator')
  call near(get_bp_avg_spectrum(1,x**(-2.7d0)),1.00065665502205d0,2.d-13,'measured CSV vs independent result')
  call near(get_bp_avg_spectrum(1,x**(-2.7d0)),1.0006863d0,3.1d-5,'thesis polynomial approximation')
  s=2.3d0*x**(-3.1d0)
  t=0.8d0*x**(-2.1d0)+0.4d0*exp(-log(x/1.3d0)**2/0.08d0)/x**2
  call near(get_bp_avg_spectrum(1,s+t),get_bp_avg_spectrum(1,s)+get_bp_avg_spectrum(1,t),2.d-14,'curved mixture linearity')

  call compute_ant2thermo(bp(1)%nu,cmb)
  cmb=1.d0/cmb
  ! Existing CMB/SZ special-component paths must give the same response as integration.
  call near(get_bp_avg_spectrum(1,cmb),1.d0/bp(1)%a2t,2.d-14,'thermodynamic conversion consistency')
  call near(get_bp_avg_spectrum(1,cmb*sz_thermo(bp(1)%nu))*1.d6,1.d0/bp(1)%a2sz,1.d-8,'SZ conversion consistency')
  call near(ant2data(1),1.d0,0.d0,'RJ map output')
  call near(bp(1)%a2t/bp(1)%f2t,compute_bnu_prime_RJ_single(bp(1)%nu_c)*1.d14,1.d-14,'reference brightness conversion')
  call near(get_bp_line_ant(1,bp(1)%nu(1)),bp(1)%tau(1)*bp(1)%nu(1)/c*1.d9,1.d-12,'line at lower endpoint')
  call near(get_bp_line_ant(1,bp(1)%nu(bp(1)%n)), &
       bp(1)%tau(bp(1)%n)*bp(1)%nu(bp(1)%n)/c*1.d9,1.d-12,'line at upper endpoint')
  call near(get_bp_line_ant(1,100.d9),0.d0,0.d0,'out-of-band line')

  ! Scale cancels, and repeated updates cannot renormalize already normalized weights.
  original=get_bp_avg_spectrum(1,s)
  bp(1)%tau0=bp(1)%tau0*7.1d7
  delta=0.d0
  call update_tau(delta)
  call near(get_bp_avg_spectrum(1,s),original,2.d-14,'arbitrary gain scale')
  call update_tau(delta)
  call near(get_bp_avg_spectrum(1,s),original,2.d-14,'repeat update')

  ! Same physical SED, different reference frequency; no hard-coded 4.76 in production.
  oldref=bp(1)%nu_c
  bp(1)%nu_c=5.13d9
  call update_tau(delta)
  expected=original*(bp(1)%nu_c/oldref)**(bp(1)%cal_index-2.d0)
  call near(get_bp_avg_spectrum(1,s),expected,2.d-14,'arbitrary reference frequency')
  bp(1)%cal_index=0.41d0
  call update_tau(delta)
  x=bp(1)%nu/bp(1)%nu_c
  call near(get_bp_avg_spectrum(1,x**(-1.59d0)),1.d0,2.d-14,'arbitrary calibration index')

  ! Updates integrate the new bandshape and recompute conversions consistently.
  do i=1,2
     if (i==1) then
        bp_model='additive_shift'
        delta=0.025d0
     else
        bp_model='powlaw_tilt'
        delta=0.8d0
     end if
     call update_tau(delta)
     x=bp(1)%nu/bp(1)%nu_c
     call near(get_bp_avg_spectrum(1,x**(bp(1)%cal_index-2)),1.d0,3.d-14,'updated band calibration')
     call compute_ant2thermo(bp(1)%nu,cmb)
     call near(get_bp_avg_spectrum(1,1.d0/cmb),1.d0/bp(1)%a2t,2.d-14,'updated conversion')
  end do

  ! Regression checks on original conventions with their native normalizations.
  bp_model='additive_shift'
  delta=0.d0
  bp(1)%id='LFI'
  bp(1)%a2t=-1; bp(1)%f2t=-1; bp(1)%a2sz=-1
  call update_tau(delta)
  call compute_ant2thermo(bp(1)%nu,cmb)
  call near(get_bp_avg_spectrum(1,1.d0/cmb),1.d0,2.d-14,'legacy LFI')
  bp(1)%id='HFI_cmb'
  bp(1)%a2t=-1; bp(1)%f2t=-1; bp(1)%a2sz=-1
  call update_tau(delta)
  call near(get_bp_avg_spectrum(1,1.d0/cmb),1.d0,2.d-14,'legacy HFI_cmb')
  bp(1)%id='WMAP'
  bp(1)%tau0=bp(1)%tau0/sum(bp(1)%tau0)
  bp(1)%a2t=-1; bp(1)%f2t=-1; bp(1)%a2sz=-1
  call update_tau(delta)
  call near(get_bp_avg_spectrum(1,1.d0/cmb),1.d0,2.d-14,'legacy WMAP')

  ! Narrow-band limit about an arbitrary reference frequency.
  deallocate(bp(1)%nu0,bp(1)%tau0,bp(1)%nu,bp(1)%tau,x,s,t,cmb)
  bp(1)%id='CBASS'
  bp(1)%n=3
  allocate(bp(1)%nu0(3),bp(1)%tau0(3),bp(1)%nu(3),bp(1)%tau(3),x(3))
  bp(1)%nu0=bp(1)%nu_c+[-1.d0,0.d0,1.d0]
  bp(1)%tau0=[0.d0,1.d0,0.d0]
  call update_tau(delta)
  x=bp(1)%nu/bp(1)%nu_c
  call near(get_bp_avg_spectrum(1,3.2d0*x**(-3.7d0)),3.2d0,2.d-14,'narrow-band limit')
  print *, 'PASS: CBASS initialization, measured integration, mixtures, parameters, conversions, updates and legacy bands'
contains
  subroutine near(actual,wanted,tol,label)
    real(dp), intent(in) :: actual,wanted,tol
    character(len=*), intent(in) :: label
    if (.not. ieee_is_finite(actual) .or. abs(actual-wanted)>tol) then
       print *, 'FAIL: ',label,actual,wanted,tol
       error stop 1
    end if
  end subroutine
end program
