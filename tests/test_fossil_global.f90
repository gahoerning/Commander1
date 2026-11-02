program test_fossil_global
  ! Observe the public channel signal, without accessing private sampler state.
  use comm_utils
  use comm_bp_mod
  use comm_fossil_sed_mod
  use comm_fossil_global_mod
  implicit none
  character(len=4096) :: config,likelihood,mode,directory,argument
  real(dp), allocatable :: weight(:),data(:),saved(:),basis(:,:)
  real(dp) :: nu,amplitude1,amplitude2,determinant,target1,target2,beta,x,amplitude3,error
  type(fossil_sed_table) :: crr
  integer :: unit,draws,i,b,status,seed
  call get_command_argument(1,config)
  call get_command_argument(2,likelihood)
  call get_command_argument(3,mode)
  call get_command_argument(4,argument); read(argument,*) draws
  call get_command_argument(5,directory)
  call get_command_argument(6,argument); read(argument,*) seed
  open(newunit=unit,file=trim(likelihood),status='old')
  read(unit,*) numband
  allocate(bp(numband),weight(numband),data(numband),saved(numband),basis(numband,3))
  do b=1,numband
     read(unit,*) nu,weight(b),data(b)
     bp(b)%nu_c=nu*1.d9
  end do
  close(unit)
  call initialize_fossil_globals(trim(config),seed)
  call broadcast_global_signal(MPI_COMM_WORLD)
  if (trim(mode)=='initialize') stop
  if (trim(mode)=='restore') then
     call save_fossil_globals()
     call update_fossil_globals(weight,data,.false.)
     call update_fossil_globals(weight,data,.true.)
     saved=global_signal
     call restore_fossil_globals()
     call update_fossil_globals(weight,data,.false.)
     call update_fossil_globals(weight,data,.true.)
     if (any(saved/=global_signal)) error stop 'Save/restore failed to preserve full RNG and component state'
     print *, 'PASS: exact in-memory state/RNG restore'
     stop
  end if
  if (trim(mode)=='amplitudes') then
     call fossil_load_sed_table(trim(directory)//'/crr_MJy_sr.txt',crr,status)
     if (status/=0) error stop 'Test CRR table load'
     do b=1,numband
        basis(b,1)=fossil_sed_value('delta_T',bp(b)%nu_c/1.d9,[real(dp)::])
        basis(b,2)=fossil_sed_value('mu',bp(b)%nu_c/1.d9,[real(dp)::])
        basis(b,3)=fossil_sed_value('crr',bp(b)%nu_c/1.d9,[real(dp)::],table=crr)
     end do
     determinant=basis(1,1)*basis(2,2)-basis(2,1)*basis(1,2)
  end if
  open(newunit=unit,file='draws.txt',status='replace')
  do i=1,draws
     select case(trim(mode))
     case('amplitudes')
        call update_fossil_globals(weight,data,.false.)
        amplitude3=1.7d0
        target1=global_signal(1)-amplitude3*basis(1,3)
        target2=global_signal(2)-amplitude3*basis(2,3)
        amplitude1=(target1*basis(2,2)-target2*basis(1,2))/determinant
        amplitude2=(basis(1,1)*target2-basis(2,1)*target1)/determinant
        error=maxval(abs(global_signal-amplitude1*basis(:,1)-amplitude2*basis(:,2)-amplitude3*basis(:,3)))
        write(unit,'(3es26.17)') amplitude1,amplitude2,error
     case('shapes')
        call update_fossil_globals(weight,data,.true.)
        x=6.62607015d-34*bp(1)%nu_c/(1.380649d-23*18.8d0)
        beta=log(global_signal(1)*(exp(x)-1.d0)/0.346d0)/log(x)-3.d0
        write(unit,'(es26.17)') beta
     case('chain')
        call update_fossil_globals(weight,data,.false.)
        call update_fossil_globals(weight,data,.true.)
        write(unit,'(*(es26.17))') global_signal
     case default
        error stop 'Unknown test mode'
     end select
  end do
  close(unit)
  call output_fossil_globals('.',1,fossil_restart_iteration+draws)
  print *, 'PASS: produced ',draws,' test draws in mode ',trim(mode)
end program test_fossil_global
