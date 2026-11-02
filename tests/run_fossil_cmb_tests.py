#!/usr/bin/env python3
"""Compile production CMB units, free-free SED, covariance and CG procedures.

These isolated fixtures execute actual production arithmetic. MPI, parameter
I/O and unrelated polarized matrix routines use explicit small test doubles;
this does not replace the full remote HEALPix/MPI acceptance test.
"""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

import numpy as np
from scipy.interpolate import CubicSpline

import run_fossil_local_sed_tests as local
import run_fossil_integration_tests as integration

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src/commander'


def compile_run(work, files, name):
    compiler = shlex.split(os.environ.get('FC', 'gfortran'))
    flags = shlex.split(os.environ.get('FFLAGS', ''))
    target = work / name
    subprocess.run(compiler + flags + ['-std=f2008', '-ffree-line-length-none', '-fcheck=all',
                   '-O0', '-J', str(work), '-I', str(work), *map(str, files), '-o', str(target)],
                   cwd=work, check=True)
    completed = subprocess.run([str(target)], cwd=work, text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(completed.stdout + completed.stderr)
    return completed.stdout


def unit_test(work):
    bp = (SRC / 'comm_bp_mod.f90').read_text()
    fg = (SRC / 'comm_fg_component_mod.f90').read_text()
    # Extract the production npar=0 initialization branch, rather than inventing
    # another CMB mixing-matrix implementation in the fixture.
    branch = fg.split("else if (trim(fg_components(i)%type) == 'cmb') then", 1)[1]
    branch = branch.split('else if', 1)[0]
    extra = '\n'.join(local.procedure(bp, name) for name in (
        'compute_ant2thermo', 'compute_bnu_prime', 'get_bp_avg_spectrum'))
    extra += '\n' + local.procedure(fg, 'spec2data')
    extra += '\nsubroutine initialize_cmb_response(i)\ninteger :: i,j\n' + branch + '\nend subroutine\n'
    # Keep the foreground source byte-for-byte unchanged: initialize and run
    # its standard 100-point free-free Te table in native uK_cmb band units.
    freefree = local.procedure(fg, 'update_eff_fg_spectrum')
    freefree = freefree.split('else if (fg_components(i)%npar == 1) then', 1)[1]
    freefree = freefree.split('else if (fg_components(i)%npar == 2) then', 1)[0]
    extra += ('\nsubroutine initialize_freefree_response(i)\n'
              'integer :: i,l,n,k,j\nreal(dp) :: nu\nreal(dp),allocatable :: s(:)\n'
              + freefree + '\nend subroutine\n')
    extra += ('\nreal(dp) function tsum(x,y)\nreal(dp) :: x(:),y(:)\n'
              "error stop 'Unexpected non-delta band in native unit fixture'\nend function\n")
    guard = bp[bp.index('       if (fossil_si_cmb) then'):]
    guard = guard.split('       bp(i)%delta     = 0.d0', 1)[0]
    extra += ('\nsubroutine production_si_unit_guard(i)\ninteger :: i,q\n'
              "character(len=4) :: paramfile='test'\ncharacter(len=2) :: i_text='01'\nlogical :: fixed\n"
              + guard + '\nend subroutine\n')
    extra += ('\nsubroutine get_parameter(file,key,par_lgt)\n'
              'character(len=*),intent(in)::file,key\nlogical,intent(out)::par_lgt\n'
              'par_lgt=.true.\nend subroutine\n')
    context = local.build_context().rsplit('end module', 1)[0] + extra + '\nend module\n'
    context = context.replace('integer, parameter :: num_recent_point=3',
                              'integer, parameter :: numgrid=100,num_recent_point=3')
    (work / 'units_context.f90').write_text(context)
    (work / 'units_driver.f90').write_text(r'''
program test_cmb_units
  use fossil_local_test_context
  use comm_fossil_sed_mod, only: fossil_sed_value
  implicit none
  integer :: i,mode,j
  real(dp) :: nu,legacy_rj(131),old_local(131),old_cmb(131),response(131),conversion(131)
  real(dp) :: value,derivative,scale,empty(0),te(1),delta_response,cmb_intensity
  character(len=20) :: argument
  allocate(bp(131),fg_components(2))
  call get_command_argument(1,argument)
  if(trim(argument)=='invalid_native_unit')then
     fossil_si_cmb=.true.;bp(1)%unit='uK_ant';bp(1)%id='delta'
     bp(1)%gain_rms=0.d0;bp(1)%a2t=0.d0;bp(1)%f2t=0.d0
     call production_si_unit_guard(1)
     error stop 'Invalid band unit unexpectedly passed SI guard'
  end if
  fg_components(1)%type='cmb';fg_components(1)%nu_ref=50.d9
  fg_components(1)%amp_unit='MJy/sr';fg_components(1)%ref_band=1
  allocate(fg_components(1)%S_tabulated(131))
  fg_components(2)%type='freefree';fg_components(2)%nu_ref=50.d9;fg_components(2)%ref_band=1
  fg_components(2)%amp_unit='MJy/sr';fg_components(2)%npar=1
  allocate(fg_components(2)%priors(1,3),fg_components(2)%p_rms(1))
  allocate(fg_components(2)%par(100,1),fg_components(2)%S_1D(100,2,131))
  allocate(fg_components(2)%S_tabulated(131))
  fg_components(2)%S_tabulated=1.d0
  fg_components(2)%priors(1,1:2)=[3000.d0,15000.d0];fg_components(2)%p_rms=0.d0
  do j=1,100
     fg_components(2)%par(j,1)=1800.d0+(16200.d0-1800.d0)*(j-1)/99.d0
  end do
  te=7000.d0
  do mode=0,1
     fossil_si_cmb=mode==1
     do i=1,131
        nu=(50.d0+15.d0*(i-1))*1.d9
        bp(i)%nu_c=nu;bp(i)%unit='MJy/sr';bp(i)%id='delta';bp(i)%gain=1.d0
        bp(i)%gain_rms=0.d0;bp(i)%a2t=0.d0;bp(i)%f2t=0.d0
        call production_si_unit_guard(i)
        bp(i)%unit='uK_cmb'
        call production_si_unit_guard(i)
        bp(i)%unit='MJy/sr'
        bp(i)%use_color_corr=.false.
        if(mode==0)allocate(bp(i)%tau(1),bp(i)%nu(1))
        bp(i)%n=1;bp(i)%nu(1)=nu
        bp(i)%a2t=compute_ant2thermo_single(nu)
        bp(i)%f2t=1.d-14/compute_bnu_prime_single(nu)
        bp(i)%tau(1)=bp(i)%a2t
        if(mode==0)then
           legacy_rj(i)=compute_bnu_prime_RJ_single(nu)
           old_local(i)=ant2data(i)/ant2unit('MJy/sr',nu,i)
           old_cmb(i)=compute_bnu_prime_single(nu)
        else
           if(compute_bnu_prime_RJ_single(nu)/=legacy_rj(i))error stop 'RJ definition changed'
           if(abs(ant2data(i)/ant2unit('MJy/sr',nu,i)-old_local(i))>3.d-15) &
                error stop 'Local MJy/RJ forward ratio changed'
        end if
     end do
     call initialize_cmb_response(1)
     call initialize_freefree_response(2)
     if(mode==1)then
        call compute_bnu_prime(bp%nu_c,response)
        call compute_ant2thermo(bp%nu_c,conversion)
        do i=1,131
           nu=bp(i)%nu_c
           delta_response=fossil_sed_value('delta_t',nu/1.d9,empty)
           cmb_intensity=get_effective_fg_spectrum(fg_components(1),i,empty) * &
                spec2data(1,cmb=.true.)/ant2unit('MJy/sr',50.d9,1)
           bp(i)%unit='uK_cmb'
           value=get_effective_fg_spectrum(fg_components(2),i,te)/ant2unit('MJy/sr',50.d9,1)/bp(i)%f2t
           derivative=get_effective_deriv_fg_spectrum(fg_components(2),i,1,te)/ &
                ant2unit('MJy/sr',50.d9,1)/bp(i)%f2t
           bp(i)%unit='MJy/sr'
           if(abs(response(i)/compute_bnu_prime_single(nu)-1.d0)>1.d-15)error stop 'Array derivative mismatch'
           if(abs(conversion(i)/compute_ant2thermo_single(nu)-1.d0)>1.d-15)error stop 'Array RJ conversion mismatch'
           write(*,'(*(es26.17,1x))')nu,response(i),conversion(i),legacy_rj(i),delta_response, &
                spec2data(i,cmb=.true.),cmb_intensity,value,derivative,old_cmb(i)
        end do
     end if
  end do
end program
''')
    output = compile_run(work, [work / 'units_context.f90', SRC / 'comm_fossil_sed_mod.f90',
                               work / 'units_driver.f90'], 'units_test')
    rejected = subprocess.run([str(work/'units_test'),'invalid_native_unit'], cwd=work,
                              capture_output=True,text=True)
    assert rejected.returncode != 0 and 'requires delta MJy/sr or uK_cmb bands' in rejected.stdout
    values = np.loadtxt(output.splitlines())
    assert values.shape == (131, 10), values.shape
    nu = values[:, 0]
    h, k, c, temperature = 6.62607015e-34, 1.380649e-23, 299792458., 2.725
    x = h * nu / (k * temperature)
    derivative = 2 * h**2 * nu**4 / (c**2 * k * temperature**2) * np.exp(x) / np.expm1(x)**2
    np.testing.assert_allclose(values[:, 1], derivative, rtol=8e-15)
    np.testing.assert_allclose(values[:, 2], values[:, 3] / derivative, rtol=8e-15)
    np.testing.assert_allclose(values[:, 4], temperature * derivative * 1e20, rtol=8e-15)
    np.testing.assert_allclose(values[:, 5:7], (derivative * 1e14)[:, None] * np.ones((1, 2)), rtol=8e-15)
    assert np.max(np.abs(values[:, 1] / values[:, 9] - 1)) > 1e-6
    # Independent natural cubic interpolation of the unchanged production
    # Draine table. This tests its native-unit likelihood and derivative,
    # without pretending the finite Te table is an exact analytic evaluation.
    te, reference = 7000., 50e9
    legacy_h = 1.0545726691251021e-34 * 2 * np.pi
    legacy_k = 1.3806503e-23
    coefficient = np.sqrt(3) / np.pi
    grid = np.linspace(1800.,16200.,100)
    def gaunt(freq):
        return np.log(np.exp(5.960-coefficient*np.log(freq/1e9*(grid[:,None]/1e4)**-1.5))+2.71828)
    tabulated = gaunt(nu[None,:])/gaunt(reference)*np.exp(
        -legacy_h*(nu[None,:]-reference)/(legacy_k*grid[:,None]))
    spline = CubicSpline(grid,tabulated,axis=0,bc_type='natural')
    np.testing.assert_allclose(values[:, 7], spline(te), rtol=2e-14)
    np.testing.assert_allclose(values[:, 8], spline(te,1), rtol=5e-8, atol=2e-14)
    print('PASS: production SI CMB derivative, DeltaT proportionality, fixed/pixel/harmonic units at 131 bands')
    print('PASS: unchanged legacy RJ/local units and standard native free-free Te table/derivative')


STUBS = r'''
module healpix_types
 implicit none
 integer,parameter::dp=kind(1.d0),sp=kind(1.e0),i4b=kind(1),lgt=kind(.true.)
end module
module comm_utils
 use healpix_types
 implicit none
 real(dp),parameter::pi=acos(-1.d0)
 integer,parameter::MPI_COMM_WORLD=0
contains
 subroutine get_parameter(file,key,par_int,par_present)
  character(len=*),intent(in)::file,key
  integer,intent(out),optional::par_int
  logical,intent(out),optional::par_present
  if(key/='FOSSIL_CMB_MIN_L')error stop 'Unexpected covariance parameter key'
  if(present(par_int))par_int=2
  if(present(par_present))par_present=.true.
 end subroutine
 subroutine mpi_abort(comm,code,ierr)
  integer::comm,code,ierr
  error stop 'Unexpected MPI abort'
 end subroutine
 integer function lm2ind(l,m)
  integer::l,m
  lm2ind=l*l+l+m+1
 end function
end module
module math_tools
 use healpix_types
 implicit none
contains
 subroutine cl2s(cl,s)
  real(dp)::cl(:),s(:,:)
  s=0.d0;s(1,1)=cl(1)
 end subroutine
 subroutine compute_hermitian_root(a,power)
  real(dp)::a(:,:),power
  error stop 'Unexpected polarized covariance'
 end subroutine
 subroutine invert_matrix(a)
  real(dp)::a(:,:)
  error stop 'Unexpected polarized inverse'
 end subroutine
end module
'''


def covariance_cg_test(work):
    procedures = [integration.procedure('comm_cgd_mod.f90', 'draw_constrained_realization_by_CG')]
    procedures += [integration.procedure('comm_genvec_mod.f90', name) for name in (
        'allocate_genvec', 'deallocate_genvec', 'genvec_set_equal', 'genvec_minus',
        'genvec_plus_times', 'genvec_dot_product')]
    context = (ROOT / 'tests/fossil_integration_support.f90.in').read_text().split('\nprogram test', 1)[0]
    context = context.replace('module context\n',
        'module context\n  use comm_S_mult_mod, only: multiply_by_sqrt_S,multiply_by_inv_sqrt_S\n')
    context = context.replace('numcomp=4,nspec=1,lmax=1,lmax_lowres=1', 'numcomp=9,nspec=1,lmax=2,lmax_lowres=2')
    for name in ('multiply_by_sqrt_S', 'multiply_by_inv_sqrt_S'):
        context = re.sub(rf'  subroutine {name}\b.*?  end subroutine\n', '', context, flags=re.S)
    context = context.replace('@PROCEDURES@', '\n'.join(procedures))
    driver = r'''
program covariance_cg
 use context
 use comm_S_mult_mod,only:initialize_S_mult_mod,update_S_mult_mod,clean_up_S_mult_mod,get_Cl
 implicit none
 type(genvec)::b,x,initial
 real(dp)::cls(0:2,1),z(9,1),expected
 integer::stat
 call initialize_S_mult_mod(2,1,'mode2')
 cls=100.d0
 call update_S_mult_mod(cls)
 if(get_Cl(0,1,1)/=0.d0 .or. get_Cl(1,1,1)/=0.d0)error stop 'Low-mode covariance nonzero'
 expected=100.d0/(6.d0/(2.d0*acos(-1.d0)))
 call check(get_Cl(2,1,1),expected,'D_l converts to C_l')
 z=1.d0;call multiply_by_sqrt_S(.false.,z)
 if(any(z(1:4,:)/=0.d0))error stop 'sqrtS leaves low modes'
 call check(z(5,1),sqrt(expected),'High modes preserved')
 z=1.d0;call multiply_by_inv_sqrt_S(.false.,z)
 if(any(z(1:4,:)/=0.d0))error stop 'inv sqrtS leaves low modes'
 call allocate_genvec(b);call allocate_genvec(x);call allocate_genvec(initial)
 b%cmb_amp=2.d0;b%fg_amp=0.d0;b%temp_amp=0.d0
 initial%cmb_amp=1.d5;initial%fg_amp=0.d0;initial%temp_amp=0.d0
 enforce_zero_cl=.false.;fossil_exact_local_seds=.true.;stat=0
 call draw_constrained_realization_by_CG(2,b,x,stat,initial)
 if(stat/=0)error stop 'CG failed'
 if(any(x%cmb_amp(1:4,:)/=0.d0))error stop 'CG reintroduced monopole or dipole'
 call check(x%cmb_amp(5,1),sqrt(expected),'Physical CG high mode')
 call clean_up_S_mult_mod
 call initialize_S_mult_mod(2,1)
 call update_S_mult_mod(cls)
 call check(get_Cl(0,1,1),1.d12,'Legacy low-mode covariance preserved')
 print *,'PASS: production covariance, exact l0/l1 zero and CG initialization/output; legacy default unchanged'
end program
'''
    (work / 'stubs.f90').write_text(STUBS)
    (work / 'cg_context.f90').write_text(context + driver)
    print(compile_run(work, [work / 'stubs.f90', SRC / 'comm_S_mult_mod.f90', work / 'cg_context.f90'], 'cg_test').strip())


def spectrum_initialization_test(work):
    text = (SRC / 'comm_Cl_util_mod.f90').read_text()
    actual = '\n'.join(local.procedure(text, name) for name in (
        'initialize_Cl_util_mod', 'initialize_ran_powspec', 'cleanup_Cl_util_mod'))
    source = STUBS.split('module comm_utils', 1)[0] + r'''
module rngmod
 use healpix_types
 implicit none
 type planck_rng
 integer::unused=0
 end type
contains
 real(dp) function rand_gauss(handle)
  type(planck_rng)::handle
  rand_gauss=0.d0
 end function
end module
module comm_utils
 use healpix_types
 implicit none
 integer,parameter::MPI_COMM_WORLD=0
contains
 subroutine get_parameter(file,key,par_int,par_lgt,par_dp,par_string,par_present)
  character(len=*),intent(in)::file,key
  integer,intent(out),optional::par_int
  logical,intent(out),optional::par_lgt,par_present
  real(dp),intent(out),optional::par_dp
  character(len=*),intent(out),optional::par_string
  if(present(par_present))par_present=.true.
  select case(key)
  case('FOSSIL_CMB_MIN_L')
   par_int=2
   if(file=='legacy')then
    par_int=0;par_present=.false.
   end if
  case('LMAX');par_int=2
  case('POLARIZATION','CL_BINNING','ENFORCE_ZERO_CL');par_lgt=.false.
  case('SAMPLE_TT_MODES','CONDITION_ON_CL');par_lgt=.true.
  case('DISPERSION_FACTOR');par_dp=0.d0
  case('IN_POWSPECFILE');par_string='input_dl.txt'
  case default;error stop 'Unexpected input-spectrum parameter'
  end select
 end subroutine
 integer function getlun()
  getlun=91
 end function
 logical function is_posdef(cl)
  real(dp)::cl(:)
  is_posdef=all(cl>0.d0)
 end function
 subroutine mpi_abort(comm,code,ierr)
  integer::comm,code,ierr
  error stop 'Unexpected MPI abort'
 end subroutine
 subroutine mpi_finalize(ierr)
  integer::ierr
  error stop 'Unexpected MPI finalize'
 end subroutine
end module
module spectrum_context
 use healpix_types
 use rngmod
 use comm_utils
 implicit none
 integer,parameter::MAX_NUM_BIN=10000
 integer::num_Cl_bin,nspec,lmax,fossil_cmb_min_l=0
 integer,allocatable::Cl_bins(:,:)
 character(len=1),allocatable::Cl_bin_stat(:,:)
contains
''' + actual + r'''
end module
program input_spectrum_test
 use spectrum_context
 implicit none
 type(planck_rng)::rng
 real(dp)::cl(0:2,1)
 call initialize_Cl_util_mod('mode2')
 call initialize_ran_powspec('mode2',rng,cl)
 if(any(cl(0:1,:)/=0.d0))error stop 'Input spectrum reintroduced monopole/dipole'
 if(cl(2,1)/=100.d0)error stop 'Input D_l changed'
 call cleanup_Cl_util_mod
 call initialize_Cl_util_mod('legacy')
 call initialize_ran_powspec('legacy',rng,cl)
 if(any(cl(0:1,1)/=1.d3))error stop 'Legacy input spectrum changed'
 print *,'PASS: production D_l reader and exact low-mode initialization; absent flag retains legacy behavior'
end program
'''
    (work / 'input_dl.txt').write_text('# ell D_l_TT sigma_D_l\n2 100.0 0.0\n')
    (work / 'spectrum.f90').write_text(source)
    print(compile_run(work, [work / 'spectrum.f90'], 'spectrum_test').strip())


def pixel_preconditioner_test(work):
    source = (SRC / 'comm_cgd_precond_mod.f90').read_text()
    # The complete production pixel block is retained, including its inversion
    # and fixed-column handling. Only the unused harmonic block is omitted.
    block = source.split('    ! Set up the pixel foreground preconditioner', 1)[1]
    block = block.split('  end subroutine initialize_preconditioner', 1)[0]
    fixture = r'''
module preconditioner_context
 implicit none
 integer,parameter::dp=kind(1.d0),npix=2,nmaps=1,num_fg_signal=2,map_size=2
 integer,parameter::MPI_DOUBLE_PRECISION=0,MPI_SUM=0
 integer::num_global_components=1,myid_chain=0,root=0,map_id=1,mask_state=1,comm_chain=0,ierr
 integer::pixels(0:1)=[0,1]
 logical::sample_fg_pix=.true.,enforce_zero_cl=.false.,fossil_exact_local_seds=.true.
 type component
  logical::enforce_positive_amplitude=.false.,sample_amplitudes=.true.
  character(len=16)::type='power_law'
 end type
 type(component)::fg_components(2)
 real(dp)::M_fg_pix(2,2,0:1,1),fg_pix_spec_response(0:1,1,2)
contains
 subroutine production_pixel_preconditioner
 integer::i,j,k,l
 real(dp),allocatable::my_M_fg_pix(:,:,:,:),M_fg_pix_sum(:,:,:,:),inv_N(:)
''' + block + r'''
 end subroutine
 subroutine get_inv_N_sub_matrix(band,state,pixel,inv_N)
  integer::band,state,pixel
  real(dp)::inv_N(:)
  inv_N=2.d0
 end subroutine
 subroutine mpi_reduce(a,b,n,dtype,op,root,comm,ierr)
  real(dp)::a(:,:,:,:),b(:,:,:,:)
  integer::n,dtype,op,root,comm,ierr
  b=a;ierr=0
 end subroutine
 subroutine invert_matrix(a)
  real(dp)::a(:,:),old(2,2),det
  old=a;det=old(1,1)*old(2,2)-old(1,2)*old(2,1)
  if(det<=0.d0)error stop 'Singular pixel preconditioner'
  a(1,1)=old(2,2)/det;a(2,2)=old(1,1)/det
  a(1,2)=-old(1,2)/det;a(2,1)=-old(2,1)/det
 end subroutine
end module
program preconditioner_test
 use preconditioner_context
 implicit none
 fg_components(2)%type='cmb'
 fg_pix_spec_response(:,:,1)=3.d0;fg_pix_spec_response(:,:,2)=5.d0
 call production_pixel_preconditioner
 if(any(M_fg_pix(1,2,:,:)/=0.d0).or.any(M_fg_pix(2,1,:,:)/=0.d0))error stop 'Pixel CMB cross terms survive'
 if(any(M_fg_pix(2,2,:,:)/=1.d0))error stop 'Pixel CMB null direction not isolated'
 if(any(abs(M_fg_pix(1,1,:,:)-1.d0/18.d0)>1.d-16))error stop 'Local preconditioner changed'
 enforce_zero_cl=.true.;fg_components(2)%sample_amplitudes=.false.
 call production_pixel_preconditioner
 if(any(M_fg_pix(1,2,:,:)/=0.d0))error stop 'Frozen CMB was included'
 print *,'PASS: production pixel preconditioner excludes the synchronized harmonic CMB and frozen CMB'
end program
'''
    (work / 'preconditioner.f90').write_text(fixture)
    print(compile_run(work, [work / 'preconditioner.f90'], 'preconditioner_test').strip())


def main():
    with tempfile.TemporaryDirectory(prefix='commander-fossil-cmb-') as directory:
        work = Path(directory)
        unit_test(work)
        covariance_cg_test(work)
        spectrum_initialization_test(work)
        pixel_preconditioner_test(work)


if __name__ == '__main__':
    main()
