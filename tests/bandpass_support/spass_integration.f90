program test_bandpass
  use comm_bp_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  type(planck_rng) :: handle
  character(len=1024) :: config, mode
  real(dp), allocatable :: x(:), s(:), t(:), cmb(:)
  real(dp) :: expected, original, oldref, delta(1), beta
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
  call near(tsum(bp(1)%nu0,bp(1)%tau0),168.d6,1.d0,'DR1 useful bandwidth')
  call near(get_bp_avg_spectrum(1,1.d0+0.d0*x),1.d0,2.d-14,'constant RJ sky')
  call near(get_bp_avg_spectrum(1,bp(1)%nu),2.302666666666667d9,1.d0,'DR1 mean frequency')
  do i=1,4
     beta=-1.8d0-0.3d0*i
     call near(get_bp_avg_spectrum(1,x**beta),analytic(beta),8.d-8,'analytic two-window power law')
  end do
  s=2.3d0*x**(-3.1d0)
  t=0.8d0*x**(-2.1d0)
  call near(get_bp_avg_spectrum(1,s+t),2.3d0*analytic(-3.1d0)+0.8d0*analytic(-2.1d0),3.d-7, &
       'analytic foreground mixture')
  call near(get_bp_avg_spectrum(1,s+t),get_bp_avg_spectrum(1,s)+get_bp_avg_spectrum(1,t),2.d-14, &
       'foreground mixture linearity')
  call near(get_bp_line_ant(1,2.24d9),0.d0,0.d0,'excluded frequency gap')
  call near(get_bp_line_ant(1,2.1d9),0.d0,0.d0,'outside band')
  call near(get_bp_line_ant(1,bp(1)%nu(1)),0.d0,0.d0,'lower endpoint')
  call near(get_bp_line_ant(1,bp(1)%nu(bp(1)%n)),0.d0,0.d0,'upper endpoint')
  ! Four triangular 0.1 Hz edge ramps add 0.2 Hz to the tabulated area.
  call near(get_bp_line_ant(1,2.2d9),2.2d9/c*1.d9/(168.d6+0.2d0),1.d-10,'line in lower window')
  call near(get_bp_line_ant(1,2.3d9),2.3d9/c*1.d9/(168.d6+0.2d0),1.d-10,'line in upper window')

  call compute_ant2thermo(bp(1)%nu,cmb)
  call near(get_bp_avg_spectrum(1,1.d0/cmb),1.d0/bp(1)%a2t,2.d-14,'thermodynamic response consistency')
  call near(get_bp_avg_spectrum(1,sz_thermo(bp(1)%nu)/cmb)*1.d6,1.d0/bp(1)%a2sz,1.d-8,'SZ response consistency')
  call near(ant2data(1),1.d0,0.d0,'RJ map output')
  call near(bp(1)%a2t/bp(1)%f2t,compute_bnu_prime_RJ_single(bp(1)%nu_c)*1.d14,1.d-14,'reference brightness conversion')

  original=get_bp_avg_spectrum(1,s)
  bp(1)%tau0=bp(1)%tau0*7.1d7
  delta=0.d0
  call update_tau(delta)
  call near(get_bp_avg_spectrum(1,s),original,2.d-14,'response scale invariance')
  call update_tau(delta)
  call near(get_bp_avg_spectrum(1,s),original,2.d-14,'repeated update')
  bp(1)%nu_c=2.4d9
  call update_tau(delta)
  call near(get_bp_avg_spectrum(1,s),original,2.d-14,'same physical spectrum at another reference')
  call near(bp(1)%a2t/bp(1)%f2t,compute_bnu_prime_RJ_single(bp(1)%nu_c)*1.d14,1.d-14,'changed reference conversion')

  ! Both frequency and shape changes must retain RJ normalization.
  do i=1,2
     if (i==1) then
        bp_model='additive_shift'
        delta=0.025d0
     else
        bp_model='powlaw_tilt'
        delta=0.8d0
     end if
     call update_tau(delta)
     call near(get_bp_avg_spectrum(1,1.d0+0.d0*x),1.d0,2.d-14,'updated RJ normalization')
     call compute_ant2thermo(bp(1)%nu,cmb)
     call near(get_bp_avg_spectrum(1,1.d0/cmb),1.d0/bp(1)%a2t,2.d-14,'updated conversion')
  end do

  deallocate(bp(1)%nu0,bp(1)%tau0,bp(1)%nu,bp(1)%tau,x,s,t,cmb)
  bp_model='additive_shift'
  delta=0.d0
  bp(1)%n=3
  allocate(bp(1)%nu0(3),bp(1)%tau0(3),bp(1)%nu(3),bp(1)%tau(3),x(3))
  bp(1)%nu0=bp(1)%nu_c+[-1.d0,0.d0,1.d0]
  bp(1)%tau0=[0.d0,1.d0,0.d0]
  call update_tau(delta)
  x=bp(1)%nu/bp(1)%nu_c
  call near(get_bp_avg_spectrum(1,3.2d0*x**(-3.7d0)),3.2d0,2.d-14,'narrow-band limit')
  print *, 'PASS: SPASS_RJ initialization, analytic windows, gap, mixtures, units, updates and narrow-band limit'
contains
  real(dp) function analytic(beta)
    real(dp), intent(in) :: beta
    real(dp), parameter :: ref=2.303d0
    ! Analytic uniform RJ mean over the published intervals, independent of the grid.
    analytic=((2.216d0**(beta+1)-2.176d0**(beta+1)) + &
         (2.400d0**(beta+1)-2.272d0**(beta+1))) / ((beta+1)*0.168d0*ref**beta)
  end function
  subroutine near(actual,wanted,tol,label)
    real(dp), intent(in) :: actual,wanted,tol
    character(len=*), intent(in) :: label
    if (.not. ieee_is_finite(actual) .or. abs(actual-wanted)>tol) then
       print *, 'FAIL: ',label,actual,wanted,tol
       error stop 1
    end if
  end subroutine
end program
