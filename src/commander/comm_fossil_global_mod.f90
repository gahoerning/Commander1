module comm_fossil_global_mod
  ! Physical, spatially constant sky components. Amplitudes use native physical
  ! units; their channel responses and all FOSSIL data are in MJy/sr.
  use comm_utils
  use comm_bp_mod, only: bp, numband
  use comm_fossil_sed_mod, only: fossil_sed_table, fossil_sed_value, fossil_load_sed_table
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: num_global_components, global_signal, global_fixed_instrument
  public :: initialize_fossil_globals, broadcast_global_signal, update_fossil_globals
  public :: output_fossil_globals, save_fossil_globals, restore_fossil_globals
  public :: fossil_global_error, fossil_restart_iteration

  type global_component
     character(len=32) :: kind, label
     real(dp) :: amplitude=0.d0, scale=1.d0, prior_mean=0.d0, prior_sigma=-1.d0
     logical :: sample_amplitude=.true.
     integer :: nshape=0
     real(dp) :: shape(2)=0.d0, low(2)=0.d0, high(2)=0.d0
     real(dp) :: shape_mean(2)=0.d0, shape_sigma(2)=-1.d0
     logical :: sample_shape(2)=.false.
     type(fossil_sed_table) :: table
  end type
  integer :: num_global_components=0, slice_steps=8, fossil_restart_iteration=0
  logical :: global_fixed_instrument=.true.
  type(global_component), allocatable :: components(:), previous(:)
  real(dp), allocatable :: response(:,:), global_signal(:)
  type(planck_rng) :: global_rng, previous_rng
  character(len=16) :: global_operation
  character(len=512) :: restart_file='none'
contains
  subroutine fossil_global_error(message)
    character(len=*), intent(in) :: message
    integer :: ierr
    write(*,*) 'FOSSIL global components: ', trim(message)
    call mpi_abort(MPI_COMM_WORLD, 1, ierr)
  end subroutine

  subroutine initialize_fossil_globals(paramfile, seed)
    character(len=*), intent(in) :: paramfile
    integer, intent(in) :: seed
    character(len=2) :: tag, ptag
    character(len=128) :: key
    character(len=512) :: filename
    logical :: found, polar, fixed
    integer :: i,j,b,status,ntemp,unit,ios,version,nstored
    real(dp) :: rms
    call get_parameter(paramfile,'NUM_GLOBAL_COMPONENTS',par_int=num_global_components,par_present=found)
    if (.not. found) num_global_components=0
    allocate(global_signal(numband)); global_signal=0.d0
    if (num_global_components==0) return
    if (num_global_components<0 .or. num_global_components>99) call fossil_global_error('Invalid component count.')
    allocate(components(num_global_components),response(numband,num_global_components))
    call rand_init(global_rng,seed)
    ! rand_init leaves the unused Box--Muller cache undefined. Persist a fully
    ! defined state even when only uniform spectral draws have been requested.
    global_rng%gset=0.d0
    call get_parameter(paramfile,'OPERATION',par_string=global_operation)
    if (global_operation/='sample' .and. global_operation/='optimize') &
         call fossil_global_error('OPERATION must be sample or optimize.')
    call get_parameter(paramfile,'POLARIZATION',par_lgt=polar)
    if (polar) call fossil_global_error('The FOSSIL global block currently supports Stokes I only.')
    call get_parameter(paramfile,'NOISE_FORMAT',par_string=filename)
    if (trim(filename)/='rms') call fossil_global_error('Use NOISE_FORMAT=rms for global sufficient statistics.')
    call get_parameter(paramfile,'SAMPLE_INSIDE_MASK',par_lgt=fixed)
    if (fixed) call fossil_global_error('The FOSSIL baseline requires SAMPLE_INSIDE_MASK=false.')
    call get_parameter(paramfile,'GLOBAL_FIXED_INSTRUMENT',par_lgt=global_fixed_instrument,par_present=found)
    if (.not. global_fixed_instrument) call fossil_global_error('This release requires GLOBAL_FIXED_INSTRUMENT=true.')
    call get_parameter(paramfile,'GLOBAL_SLICE_STEPS',par_int=slice_steps,par_present=found)
    if (slice_steps<1) call fossil_global_error('GLOBAL_SLICE_STEPS must be positive.')
    call get_parameter(paramfile,'GLOBAL_RESTART_FILE',par_string=restart_file,par_present=found,path=base_path)
    do b=1,numband
       if (trim(bp(b)%unit)/='MJy/sr') call fossil_global_error('All FOSSIL FREQ_UNIT values must be MJy/sr.')
       ! Finite-band calibration conventions differ between instruments. Refuse
       ! them here until matching integration has been independently validated.
       if (trim(bp(b)%id)/='delta') call fossil_global_error('Use delta bands matching the monochromatic simulation.')
       if (bp(b)%gain/=1.d0 .or. bp(b)%delta/=0.d0 .or. bp(b)%delta_rms/=0.d0) &
            call fossil_global_error('Require fixed unit gains and zero bandpass shifts.')
       if (bp(b)%use_color_corr) call fossil_global_error('Disable empirical colour corrections for FOSSIL.')
       block
         character(len=:), allocatable :: band
         call index2string(b,band)
         call get_parameter(paramfile,'GAIN_PRIOR_RMS'//band,par_dp=rms)
         if (rms/=0.d0) call fossil_global_error('GAIN_PRIOR_RMS must be zero.')
         call get_parameter(paramfile,'NOISEAMP_PRIOR_RMS'//band,par_dp=rms)
         if (rms/=0.d0) call fossil_global_error('NOISEAMP_PRIOR_RMS must be zero.')
         call get_parameter(paramfile,'FIX_MONOPOLE'//band,par_lgt=fixed)
         if (.not. fixed) call fossil_global_error('Fix all instrumental monopoles.')
         call get_parameter(paramfile,'FIX_DIPOLE'//band,par_lgt=fixed)
         if (.not. fixed) call fossil_global_error('Fix all instrumental dipoles in the baseline.')
       end block
    end do
    call get_parameter(paramfile,'NUM_FG_TEMP',par_int=ntemp)
    if (ntemp/=0) call fossil_global_error('The FOSSIL baseline requires NUM_FG_TEMP=0.')
    do i=1,num_global_components
       call int2string(i,tag)
       call get_parameter(paramfile,'GLOBAL_TYPE'//tag,par_string=components(i)%kind)
       components(i)%label=components(i)%kind
       call get_parameter(paramfile,'GLOBAL_LABEL'//tag,par_string=components(i)%label,par_present=found)
       do j=1,i-1
          if (components(j)%label==components(i)%label) call fossil_global_error('Global labels must be unique.')
       end do
       call get_parameter(paramfile,'GLOBAL_INIT'//tag,par_dp=components(i)%amplitude)
       components(i)%prior_mean=components(i)%amplitude
       call get_parameter(paramfile,'GLOBAL_SAMPLE_AMPLITUDE'//tag, &
            par_lgt=components(i)%sample_amplitude,par_present=found)
       call get_parameter(paramfile,'GLOBAL_PRIOR_MEAN'//tag,par_dp=components(i)%prior_mean,par_present=found)
       call get_parameter(paramfile,'GLOBAL_PRIOR_SIGMA'//tag,par_dp=components(i)%prior_sigma,par_present=found)
       call get_parameter(paramfile,'GLOBAL_SCALE'//tag,par_dp=components(i)%scale,par_present=found)
       if (components(i)%scale<=0.d0) call fossil_global_error('GLOBAL_SCALE must be positive.')
       select case(trim(components(i)%kind))
       case('delta_t','delta_T','mu')
          components(i)%nshape=0
       case('tsz')
          components(i)%nshape=1
       case('cib')
          components(i)%nshape=2
       case('crr','co')
          call get_parameter(paramfile,'GLOBAL_TEMPLATE'//tag,par_string=filename,path=base_path)
          call fossil_load_sed_table(trim(filename),components(i)%table,status, &
               logarithmic=trim(components(i)%kind)=='co',extrapolate=trim(components(i)%kind)=='co')
          if (status/=0) call fossil_global_error('Cannot read GLOBAL_TEMPLATE: '//trim(filename))
       case default
          call fossil_global_error('Unknown GLOBAL_TYPE: '//trim(components(i)%kind))
       end select
       do j=1,components(i)%nshape
          call int2string(j,ptag); key=tag//'_'//ptag
          call get_parameter(paramfile,'GLOBAL_INIT_PARAM'//trim(key),par_dp=components(i)%shape(j))
          call get_parameter(paramfile,'GLOBAL_SAMPLE_PARAM'//trim(key), &
               par_lgt=components(i)%sample_shape(j),par_present=found)
          components(i)%shape_mean(j)=components(i)%shape(j)
          call get_parameter(paramfile,'GLOBAL_PARAM_PRIOR_MEAN'//trim(key), &
               par_dp=components(i)%shape_mean(j),par_present=found)
          call get_parameter(paramfile,'GLOBAL_PARAM_PRIOR_SIGMA'//trim(key), &
               par_dp=components(i)%shape_sigma(j),par_present=found)
          if (components(i)%sample_shape(j)) then
             call get_parameter(paramfile,'GLOBAL_PARAM_LOW'//trim(key),par_dp=components(i)%low(j))
             call get_parameter(paramfile,'GLOBAL_PARAM_HIGH'//trim(key),par_dp=components(i)%high(j))
             if (components(i)%high(j)<=components(i)%low(j) .or. &
                  components(i)%shape(j)<components(i)%low(j) .or. &
                  components(i)%shape(j)>components(i)%high(j)) &
                  call fossil_global_error('Invalid global spectral parameter bounds/start.')
          end if
       end do
    end do
    if (trim(restart_file)/='none') then
       unit=getlun(); open(unit,file=trim(restart_file),status='old',iostat=ios)
       if (ios/=0) call fossil_global_error('Cannot open GLOBAL_RESTART_FILE.')
       read(unit,*,iostat=ios) version,nstored,fossil_restart_iteration
       if (ios/=0) call fossil_global_error('Invalid global restart header.')
       if (version/=1 .or. nstored/=num_global_components .or. fossil_restart_iteration<0) &
            call fossil_global_error('Incompatible global restart header.')
       do i=1,num_global_components
          read(unit,*,iostat=ios) filename,components(i)%amplitude,components(i)%shape
          if (ios/=0) call fossil_global_error('Invalid global restart component.')
          if (trim(filename)/=trim(components(i)%label)) &
               call fossil_global_error('Incompatible global restart component.')
       end do
       read(unit,*,iostat=ios) global_rng%x,global_rng%y,global_rng%z,global_rng%w,global_rng%gset,global_rng%empty
       if (ios/=0) call fossil_global_error('Invalid global restart RNG state.')
       if (.not. ieee_is_finite(global_rng%gset) .or. &
            all([global_rng%x,global_rng%y,global_rng%z,global_rng%w]==0)) &
            call fossil_global_error('Invalid global restart RNG state.')
       close(unit)
    end if
    call validate_global_state()
    call refresh_global_response()
  end subroutine

  subroutine validate_global_state()
    integer :: i,p
    do i=1,num_global_components
       if (.not. all(ieee_is_finite([components(i)%amplitude,components(i)%scale, &
            components(i)%prior_mean,components(i)%prior_sigma]))) &
            call fossil_global_error('Non-finite global amplitude, scale or prior metadata.')
       if (components(i)%scale<=0.d0) call fossil_global_error('GLOBAL_SCALE must be positive.')
       if (.not. all(ieee_is_finite(components(i)%shape(:components(i)%nshape))) .or. &
            .not. all(ieee_is_finite(components(i)%shape_mean(:components(i)%nshape))) .or. &
            .not. all(ieee_is_finite(components(i)%shape_sigma(:components(i)%nshape)))) &
            call fossil_global_error('Non-finite global spectral parameter or prior metadata.')
       do p=1,components(i)%nshape
          if (.not. components(i)%sample_shape(p)) cycle
          if (.not. all(ieee_is_finite([components(i)%low(p),components(i)%high(p)]))) &
               call fossil_global_error('Non-finite global spectral parameter bounds.')
          if (components(i)%high(p)<=components(i)%low(p) .or. &
               components(i)%shape(p)<components(i)%low(p) .or. &
               components(i)%shape(p)>components(i)%high(p)) &
               call fossil_global_error('Invalid global spectral parameter bounds/start or restart.')
       end do
    end do
  end subroutine

  subroutine refresh_global_response()
    integer :: i,b,status
    do i=1,num_global_components
       do b=1,numband
          response(b,i)=fossil_sed_value(trim(components(i)%kind),bp(b)%nu_c/1.d9, &
               components(i)%shape(:components(i)%nshape),table=components(i)%table,status=status)
          if (status/=0 .or. .not. ieee_is_finite(response(b,i))) &
               call fossil_global_error('Invalid global SED or template frequency range.')
       end do
    end do
    call refresh_global_signal()
  end subroutine

  subroutine refresh_global_signal()
    integer :: i
    global_signal=0.d0
    do i=1,num_global_components
       global_signal=global_signal+response(:,i)*components(i)%amplitude
    end do
    if (.not. all(ieee_is_finite(global_signal))) call fossil_global_error('Non-finite global signal.')
  end subroutine

  subroutine broadcast_global_signal(comm)
    integer, intent(in) :: comm
    integer :: ierr
    if (num_global_components==0) return
    call mpi_bcast(global_signal,numband,MPI_DOUBLE_PRECISION,0,comm,ierr)
  end subroutine

  subroutine save_fossil_globals()
    if (num_global_components==0) return
    previous=components; previous_rng=global_rng
  end subroutine

  subroutine restore_fossil_globals()
    if (num_global_components==0) return
    components=previous; global_rng=previous_rng
    call refresh_global_response()
  end subroutine

  subroutine update_fossil_globals(weight, weighted_residual, spectral)
    ! The residual excludes the entire current local sky, but includes all
    ! physical global signals. Sums contain each observed pixel exactly once.
    real(dp), intent(in) :: weight(:), weighted_residual(:)
    logical, intent(in) :: spectral
    if (num_global_components==0) return
    if (size(weight)/=numband .or. size(weighted_residual)/=numband) &
         call fossil_global_error('Likelihood statistics must have one entry per band.')
    if (any(weight<0.d0) .or. .not. all(ieee_is_finite(weight)) .or. &
         .not. all(ieee_is_finite(weighted_residual))) call fossil_global_error('Invalid likelihood statistics.')
    if (sum(weight)<=0.d0) call fossil_global_error('No weighted pixels in global likelihood.')
    if (spectral) then
       call update_global_shapes(weight,weighted_residual)
    else
       call update_global_amplitudes(weight,weighted_residual)
    end if
    call refresh_global_signal()
  end subroutine

  subroutine update_global_amplitudes(weight, data)
    real(dp), intent(in) :: weight(:),data(:)
    integer :: n,i,j,k,ci,cj
    integer, allocatable :: active(:)
    real(dp), allocatable :: F(:,:),precision(:,:),rhs(:),L(:,:),d(:),x(:),z(:),r(:)
    real(dp) :: pivot
    n=count(components%sample_amplitude)
    if (n==0) return
    allocate(active(n),F(numband,n),precision(n,n),rhs(n),L(n,n),d(n),x(n),z(n),r(numband))
    k=0; r=data
    do i=1,num_global_components
       if (components(i)%sample_amplitude) then
          k=k+1; active(k)=i; F(:,k)=response(:,i)*components(i)%scale
       else
          r=r-weight*response(:,i)*components(i)%amplitude
       end if
    end do
    do i=1,n
       ci=active(i); rhs(i)=sum(F(:,i)*r)
       do j=1,n
          precision(i,j)=sum(weight*F(:,i)*F(:,j))
       end do
       if (components(ci)%prior_sigma>0.d0) then
          precision(i,i)=precision(i,i)+(components(ci)%scale/components(ci)%prior_sigma)**2
          rhs(i)=rhs(i)+components(ci)%scale*components(ci)%prior_mean/components(ci)%prior_sigma**2
       end if
       if (precision(i,i)<=0.d0) call fossil_global_error('Unconstrained global amplitude; supply an explicit prior.')
       d(i)=sqrt(precision(i,i))
    end do
    ! Whiten by the precision diagonal; no artificial ridge or hidden prior.
    do i=1,n
       precision(i,:)=precision(i,:)/d(i)/d(:); rhs(i)=rhs(i)/d(i)
    end do
    L=0.d0
    do i=1,n
       do j=1,i
          pivot=precision(i,j)-sum(L(i,:j-1)*L(j,:j-1))
          if (i==j) then
             if (pivot<=1.d-14 .or. .not. ieee_is_finite(pivot)) &
                  call fossil_global_error('Global amplitude precision is singular/ill-conditioned; inspect degeneracies.')
             L(i,j)=sqrt(pivot)
          else
             L(i,j)=pivot/L(j,j)
          end if
       end do
    end do
    do i=1,n
       x(i)=(rhs(i)-sum(L(i,:i-1)*x(:i-1)))/L(i,i)
    end do
    do i=n,1,-1
       x(i)=(x(i)-sum(L(i+1:n,i)*x(i+1:n)))/L(i,i)
    end do
    z=0.d0
    if (trim(global_operation)=='sample') then
       do i=1,n
          z(i)=rand_gauss(global_rng)
       end do
       do i=n,1,-1
          z(i)=(z(i)-sum(L(i+1:n,i)*z(i+1:n)))/L(i,i)
       end do
    end if
    do i=1,n
       ci=active(i); components(ci)%amplitude=components(ci)%scale*(x(i)+z(i))/d(i)
    end do
    if (.not. all(ieee_is_finite(components%amplitude))) call fossil_global_error('Non-finite amplitude draw.')
  end subroutine

  real(dp) function shape_log_probability(ci,p,value,weight,data) result(logp)
    integer, intent(in) :: ci,p
    real(dp), intent(in) :: value,weight(:),data(:)
    real(dp) :: shape(2),v,model,delta
    integer :: b,j,status
    logp=-huge(1.d0)/100.d0
    if (value<components(ci)%low(p) .or. value>components(ci)%high(p)) return
    shape=components(ci)%shape; shape(p)=value
    logp=0.d0
    do b=1,numband
       if (weight(b)==0.d0) cycle
       model=0.d0
       do j=1,num_global_components
          if (j==ci) then
             v=fossil_sed_value(trim(components(j)%kind),bp(b)%nu_c/1.d9, &
                  shape(:components(j)%nshape),table=components(j)%table,status=status)
             if (status/=0 .or. .not. ieee_is_finite(v)) then
                logp=-huge(1.d0)/100.d0; return
             end if
          else
             v=response(b,j)
          end if
          model=model+components(j)%amplitude*v
       end do
       delta=model-data(b)/weight(b)
       logp=logp-0.5d0*weight(b)*delta**2
    end do
    if (components(ci)%shape_sigma(p)>0.d0) &
         logp=logp-0.5d0*((value-components(ci)%shape_mean(p))/components(ci)%shape_sigma(p))**2
  end function

  subroutine update_global_shapes(weight,data)
    real(dp), intent(in) :: weight(:),data(:)
    real(dp) :: left,right,current,trial,height,lp,best,bestlp,x1,x2,f1,f2,step
    integer :: i,p,s,k,j
    logical :: accepted
    do i=1,num_global_components
       do p=1,components(i)%nshape
          if (.not. components(i)%sample_shape(p)) cycle
          if (trim(global_operation)=='sample') then
             do s=1,slice_steps
                current=components(i)%shape(p)
                lp=shape_log_probability(i,p,current,weight,data)
                if (.not. ieee_is_finite(lp) .or. lp<-huge(1.d0)/1000.d0) &
                     call fossil_global_error('Invalid current global spectral state.')
                height=lp+log(max(rand_uni(global_rng),tiny(1.d0)))
                left=components(i)%low(p); right=components(i)%high(p)
                accepted=.false.
                do k=1,10000
                   trial=left+(right-left)*rand_uni(global_rng)
                   if (shape_log_probability(i,p,trial,weight,data)>=height) then
                      accepted=.true.; exit
                   end if
                   if (trial<current) then
                      left=trial
                   else
                      right=trial
                   end if
                end do
                if (.not. accepted) call fossil_global_error('Global slice sampler failed; no ML substitution performed.')
                components(i)%shape(p)=trial
             end do
          else
             ! Coarse global search followed by bounded golden-section refinement.
             left=components(i)%low(p); right=components(i)%high(p); step=(right-left)/128.d0
             best=current_value(); bestlp=shape_log_probability(i,p,best,weight,data)
             do j=0,128
                trial=left+j*step; lp=shape_log_probability(i,p,trial,weight,data)
                if (lp>bestlp) then
                   best=trial; bestlp=lp
                end if
             end do
             left=max(left,best-step); right=min(right,best+step)
             do k=1,100
                x1=right-(right-left)*0.6180339887498949d0
                x2=left+(right-left)*0.6180339887498949d0
                f1=shape_log_probability(i,p,x1,weight,data); f2=shape_log_probability(i,p,x2,weight,data)
                if (f1>bestlp) then
                   best=x1; bestlp=f1
                end if
                if (f2>bestlp) then
                   best=x2; bestlp=f2
                end if
                if (f1>f2) then
                   right=x2
                else
                   left=x1
                end if
             end do
             components(i)%shape(p)=best
          end if
          call refresh_global_response()
       end do
    end do
  contains
    real(dp) function current_value()
      current_value=components(i)%shape(p)
    end function
  end subroutine

  subroutine output_fossil_globals(directory,chain,iteration)
    character(len=*), intent(in) :: directory
    integer, intent(in) :: chain,iteration
    character(len=4) :: tag
    character(len=512) :: filename
    integer :: unit,i,p
    logical :: exists
    if (num_global_components==0) return
    call int2string(chain,tag)
    filename=trim(directory)//'/global_components_c'//tag//'.txt'
    inquire(file=trim(filename),exist=exists)
    unit=getlun(); open(unit,file=trim(filename),status='unknown',position='append')
    if (.not. exists) then
       write(unit,'(a)',advance='no') '# iteration'
       do i=1,num_global_components
          write(unit,'(a)',advance='no') ' '//trim(components(i)%label)
          do p=1,components(i)%nshape
             block
               character(len=2) :: ptag
               call int2string(p,ptag)
               write(unit,'(a)',advance='no') ' '//trim(components(i)%label)//'_p'//ptag
             end block
          end do
       end do
       write(unit,*)
    end if
    write(unit,'(i10)',advance='no') iteration
    do i=1,num_global_components
       write(unit,'(es25.16)',advance='no') components(i)%amplitude
       do p=1,components(i)%nshape
          write(unit,'(es25.16)',advance='no') components(i)%shape(p)
       end do
    end do
    write(unit,*); close(unit)
    filename=trim(directory)//'/global_restart_c'//tag//'.txt'
    open(unit,file=trim(filename),status='replace')
    write(unit,*) 1,num_global_components,iteration
    do i=1,num_global_components
       write(unit,'(a,3es25.16)') trim(components(i)%label),components(i)%amplitude,components(i)%shape
    end do
    write(unit,*) global_rng%x,global_rng%y,global_rng%z,global_rng%w,global_rng%gset,global_rng%empty
    close(unit)
    ! Compact component spectra provide the full constant maps without storing
    ! another full-sky copy of every scalar channel prediction.
    filename=trim(directory)//'/global_spectra_c'//tag//'.txt'
    open(unit,file=trim(filename),status='replace')
    write(unit,'(a)') '# frequency_GHz total_MJy_sr followed by components in configuration order'
    do i=1,numband
       write(unit,*) bp(i)%nu_c/1.d9,global_signal(i),response(i,:)*components%amplitude
    end do
    close(unit)
  end subroutine
end module comm_fossil_global_mod
