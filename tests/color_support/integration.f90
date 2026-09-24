program integration
  use color_test_context
  implicit none
  real(dp) :: amps(0:1,3,4), index_map(0:1,3,8), sky(3), old_k, expected, f, derivative
  real(dp) :: response(0:0,3,4), data(0:1,3,2), invn(0:1,3,2), q, u, trial(0:1,3,8)
  integer :: i, b, c, pixels(0:0), n
  type(fg_params) :: pars
  type(fg_region) :: region

  allocate(bp(2), fg_components(4))
  do b=1,2
     bp(b)%use_color_corr=(b==1)
     bp(b)%gain=2.d0
     bp(b)%a2t=3.d0
     bp(b)%nu_c=6.d9
     bp(b)%delta_rms=0.d0
     bp(b)%n=2
     allocate(bp(b)%nu0(2), bp(b)%tau0(2), bp(b)%nu(2), bp(b)%tau(2))
     bp(b)%nu0=[4.d9,8.d9]
     bp(b)%nu=bp(b)%nu0
     bp(b)%tau0=1.d0
     bp(b)%tau=1.d0
  end do

  if (command_argument_count()>0) then
     call get_command_argument(1,config_mode)
     if (config_mode=='delta_without_limits') bp(1)%nu0=4.d9
     if (config_mode=='moving_bandpass') bp(1)%delta_rms=0.1d0
     call read_color_correction_parameters('fixture',1,.true.)
     error stop 'invalid configuration did not abort'
  end if
  call read_color_correction_parameters('fixture',1,.true.)
  call near(bp(1)%cc_nu_low,4.d9,'automatic lower support')
  call near(bp(1)%cc_nu_high,8.d9,'automatic upper support')
  config_mode='explicit'
  bp(1)%nu0=6.d9
  call read_color_correction_parameters('fixture',1,.true.)
  call near(bp(1)%cc_nu_low,4.d9,'explicit GHz converted to Hz for delta band')
  bp(1)%nu0=[4.d9,8.d9]
  bp(1)%cc_coeffs(:,1)=[1.d0,0.1d0,0.01d0]
  bp(1)%cc_coeffs(:,2)=[1.d0,0.2d0,0.01d0]
  bp(1)%cc_coeffs(:,3)=[1.d0,-0.1d0,0.02d0]

  do c=1,4
     fg_components(c)%type='power_law'
     fg_components(c)%nu_ref=4.d9
     fg_components(c)%npar=2
     if (c==3) then
        fg_components(c)%type='cmb'
        fg_components(c)%npar=0
     else if (c==4) then
        fg_components(c)%type='power_law_faraday_BS'
        fg_components(c)%npar=4
     end if
     n=fg_components(c)%npar
     allocate(fg_components(c)%mask(0:1,3),fg_components(c)%indmask(0:1,3))
     allocate(fg_components(c)%p_rms(n),fg_components(c)%priors(n,3),fg_components(c)%gauss_prior(n,2))
     allocate(fg_components(c)%S_tabulated(2),fg_components(c)%par(2,max(n,2)))
     allocate(fg_components(c)%S_2D(4,4,2,2,2),fg_components(c)%S_1D(2,2,2))
     fg_components(c)%mask=1.d0
     fg_components(c)%indmask=1.d0
     fg_components(c)%p_rms=0.d0
     fg_components(c)%priors(:,1)=-1.d6
     fg_components(c)%priors(:,2)=1.d6
     fg_components(c)%gauss_prior=0.d0
     fg_components(c)%gauss_prior(:,2)=-1.d0
     fg_components(c)%S_tabulated=1.d0
     fg_components(c)%par=0.d0
     fg_components(c)%S_2D=0.d0
     fg_components(c)%S_1D=0.d0
  end do
  amps=0.d0
  amps(:,:,1)=100.d0
  amps(:,:,2)=80.d0
  amps(:,:,3)=1.d20 ! A huge CMB must have no effect on the foreground slope.
  amps(:,2,1)=10.d0
  amps(:,2,2)=-9.d0
  amps(:,3,1:2)=0.d0
  index_map=0.d0
  index_map(:,:,1)=-3.d0
  index_map(:,:,3)=-2.d0
  index_map(:,:,7)=-3.d0
  call update_sky_color_corrections(amps,index_map)
  expected=2.d0+log(32.5d0/180.d0)/log(2.d0)
  call near(bp(1)%cc_alpha(1,1),expected,'mixture index excludes CMB')
  call near(bp(1)%cc_alpha(1,2),2.d0,'P of summed Q,U')
  call near(bp(1)%cc_alpha(1,3),2.d0,'same P index for U')
  call near(bp(1)%cc_factor(1,2),1.44d0,'Q polynomial')
  call near(bp(1)%cc_factor(1,3),0.88d0,'U polynomial')
  call require(.not. allocated(bp(2)%cc_factor),'disabled band has no correction map')

  old_k=bp(1)%cc_factor(1,1)
  amps(:,1,2)=0.d0
  call near(bp(1)%cc_factor(1,1),old_k,'cache is fixed within an iteration')
  call update_sky_color_corrections(amps,index_map)
  call near(bp(1)%cc_alpha(1,1),-1.d0,'new amplitudes change next iteration index')
  call require(abs(old_k-bp(1)%cc_factor(1,1))>1.d-3,'refresh changes correction')
  call reorder_fg_params(index_map(1,:,:),pars)
  f=get_effective_fg_spectrum(fg_components(1),1,pars%comp(1)%p(1,:),pixel=1,pol=1)
  expected=1.5d0**(-3.d0)*6.d0*bp(1)%cc_factor(1,1)
  call near(f,expected,'forward model multiplies K and applies units/gain once')
  derivative=get_effective_deriv_fg_spectrum(fg_components(1),1,1,pars%comp(1)%p(1,:),pixel=1,pol=1)
  call near(derivative,expected*log(1.5d0),'derivative of corrected response')
  fg_components(3)%S_tabulated=7.d0
  f=get_effective_fg_spectrum(fg_components(3),1,pars%comp(3)%p(1,:),pixel=1,pol=1)
  call near(f,14.d0,'CMB response untouched')
  f=get_effective_fg_spectrum(fg_components(1),2,pars%comp(1)%p(1,:),pixel=1,pol=1)
  call near(f,246.d0,'disabled band uses original spline response')

  ! Local row zero corresponds to sky pixel one, including its component mask.
  pixels(0)=1
  fg_components(1)%mask(0,:)=0.d0
  call update_fg_pix_response_map(1,pixels,response,index_map)
  call near(response(0,1,1),expected,'global pixel indexing in distributed response map')
  fg_components(1)%mask=1.d0

  ! Single Faraday component: endpoint P slope is independent of rotation measure.
  amps=0.d0
  amps(:,2:3,4)=1.d0
  index_map(:,:,5)=3.d0
  index_map(:,:,6)=4.d0
  index_map(:,:,8)=1000.d0
  call update_sky_color_corrections(amps,index_map)
  call near(bp(1)%cc_alpha(1,2),-1.d0,'Faraday P spectral index')
  call require(all(bp(1)%cc_factor(:,1)==1.d0),'zero intensity uses unity')

  ! Recover Q0,U0 with distinct Q/U factors in the actual block-sampler routine.
  bp(2)%use_color_corr=.true.
  allocate(bp(2)%cc_factor(0:1,3))
  bp(1)%cc_factor(:,2)=1.2d0
  bp(1)%cc_factor(:,3)=0.8d0
  bp(2)%cc_factor=1.d0
  bp(2)%cc_factor(:,2)=0.9d0
  bp(2)%cc_factor(:,3)=1.3d0
  bp(2)%nu_c=8.d9
  data=0.d0
  invn=1.d0
  do b=1,2
     call compute_faraday_rotation(bp(b)%nu_c,4.d9,3.d0,4.d0,1000.d0,-3.d0,q,u)
     data(1,2,b)=q*bp(b)%cc_factor(1,2)*6.d0
     data(1,3,b)=u*bp(b)%cc_factor(1,3)*6.d0
  end do
  region%n=1
  allocate(region%pix(1,2))
  region%pix(1,:)=[1,2]
  trial=index_map
  call sample_QU_block_region(data,invn,amps,region,5,4,trial,index_map)
  call near(trial(1,2,5),3.d0,'block sampler Q0 with distinct K_Q,K_U')
  call near(trial(1,2,6),4.d0,'block sampler U0 with distinct K_Q,K_U')
  print *, 'PASS: color correction integration fixtures'
contains
  subroutine near(actual,target,label)
    real(dp), intent(in) :: actual,target
    character(len=*), intent(in) :: label
    if (abs(actual-target)>1.d-8*max(1.d0,abs(target))) then
       print *, actual,target
       call require(.false.,label)
    end if
  end subroutine
  subroutine require(ok,label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    if (.not.ok) then
       print *, 'FAIL: ',label
       error stop 1
    end if
  end subroutine
end program integration
