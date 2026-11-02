! Minimal single-process interfaces for the complete production global module.
! RNG implementation is the actual src/include/rngmod.f90, compiled separately.
module healpix_types
  implicit none
  integer, parameter :: dp=kind(1.d0), i4b=kind(1), lgt=kind(.true.)
end module healpix_types

module comm_utils
  use healpix_types
  use rngmod
  implicit none
  integer, parameter :: MPI_COMM_WORLD=0, MPI_DOUBLE_PRECISION=1
  character(len=1), parameter :: base_path=''
contains
  integer function getlun()
    logical :: opened
    getlun=10
    do
       inquire(unit=getlun,opened=opened)
       if (.not. opened) return
       getlun=getlun+1
    end do
  end function getlun

  subroutine int2string(value,text)
    integer, intent(in) :: value
    character(len=*), intent(out) :: text
    character(len=32) :: format
    write(format,'(a,i0,a,i0,a)') '(i',len(text),'.',len(text),')'
    write(text,format) value
  end subroutine int2string

  subroutine index2string(value,text)
    integer, intent(in) :: value
    character(len=:), allocatable, intent(out) :: text
    character(len=32) :: temporary
    if (value<100) then
       write(temporary,'(i2.2)') value
    else
       write(temporary,'(i0)') value
    end if
    text=trim(temporary)
  end subroutine index2string

  subroutine get_parameter(file,key,par_int,par_dp,par_string,par_lgt,par_present,path)
    character(len=*), intent(in) :: file,key
    integer, intent(inout), optional :: par_int
    real(dp), intent(inout), optional :: par_dp
    character(len=*), intent(inout), optional :: par_string
    logical, intent(inout), optional :: par_lgt
    logical, intent(out), optional :: par_present
    character(len=*), intent(in), optional :: path
    character(len=4096) :: line,value
    integer :: unit,ios,eq
    logical :: found
    found=.false.
    open(newunit=unit,file=file,status='old')
    do
       read(unit,'(a)',iostat=ios) line
       if (ios/=0) exit
       eq=index(line,'=')
       if (eq<=1) cycle
       if (trim(adjustl(line(:eq-1)))/=key) cycle
       value=adjustl(line(eq+1:))
       if (present(par_int)) read(value,*) par_int
       if (present(par_dp)) read(value,*) par_dp
       if (present(par_string)) read(value,*) par_string
       if (present(par_lgt)) read(value,*) par_lgt
       found=.true.
       exit
    end do
    close(unit)
    if (present(par_present)) then
       par_present=found
    else if (.not. found) then
       print *, 'MISSING_PARAMETER: ',trim(key)
       error stop 1
    end if
  end subroutine get_parameter

  subroutine mpi_abort(comm,code,ierr)
    integer, intent(in) :: comm,code
    integer, intent(out) :: ierr
    ierr=code
    error stop 1
  end subroutine mpi_abort

  subroutine mpi_comm_rank(comm,rank,ierr)
    integer, intent(in) :: comm
    integer, intent(out) :: rank,ierr
    rank=0; ierr=0
  end subroutine mpi_comm_rank

  subroutine mpi_bcast(values,n,typ,root,comm,ierr)
    real(dp), intent(inout) :: values(:)
    integer, intent(in) :: n,typ,root,comm
    integer, intent(out) :: ierr
    ierr=0
  end subroutine mpi_bcast
end module comm_utils

module comm_bp_mod
  use healpix_types
  implicit none
  type bandinfo
     character(len=14) :: unit='MJy/sr', id='delta'
     real(dp) :: nu_c, gain=1.d0, delta=0.d0, delta_rms=0.d0, f2t=1.d0
     logical :: use_color_corr=.false.
  end type bandinfo
  integer :: numband
  type(bandinfo), allocatable :: bp(:)
end module comm_bp_mod
