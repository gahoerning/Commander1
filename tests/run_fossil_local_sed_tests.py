#!/usr/bin/env python3
"""Compile actual Commander local SED/units/spline procedures against FOSSIL truth.

Only unrelated component types, MPI abort and parameter I/O are replaced by
fail-fast test doubles.  The complete production direct-SED dispatcher,
derivative routine, local SEDs, ASCII reader, spline, tridiagonal solver and
unit conversion procedures are extracted verbatim.  This is an isolated
production-procedure test, not a full HEALPix/MPI Commander build.
"""
import argparse
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT/"src/commander"


def procedure(text, name):
    match = re.search(rf"^  *(?:(?:logical|pure) +)?(function|subroutine) +{name}\b.*?"
                      rf"^  *end +\1 +{name}\b[^\n]*\n", text,
                      re.MULTILINE|re.DOTALL|re.IGNORECASE)
    if match is None:
        raise RuntimeError(f"Cannot extract production procedure {name}")
    return match.group(0)


def build_context(si_cmb=False):
    fg=(SRC/"comm_fg_component_mod.f90").read_text()
    bp=(SRC/"comm_bp_mod.f90").read_text()
    spline=(ROOT/"src/include/spline_1D_mod.f90").read_text()
    math=(ROOT/"src/include/math_tools.f90").read_text()
    locate=(ROOT/"src/include/locate_mod.f90").read_text()
    types=fg[fg.index("  type fg_par_comp"):fg.index("  integer(i4b)                     :: num_fg_comp")]
    band=bp[bp.index("  type bandinfo"):bp.index("  real(dp)                                      :: ind_iras")]
    constants="\n".join(line for line in bp.splitlines() if re.match(r"\s*real\(dp\), parameter :: (k_B|h|c)\s*=",line))
    funcs=["direct_fossil_component","color_corrected_component","monochromatic_fg_spectrum",
           "foreground_color_factor","get_effective_fg_spectrum","get_effective_deriv_fg_spectrum",
           "compute_power_law_spectrum","compute_one_comp_dust_spectrum","compute_AME_freq_shift_spectrum","compute_freefree_spectrum",
           "compute_faraday_rotation","read_spectrum"]
    production="\n".join(procedure(fg,name) for name in funcs)
    ideal=procedure(fg,"get_ideal_fg_spectrum")
    # Unused spectral types are explicit fail-fast doubles; their production
    # dispatcher remains intact and tested branches retain their real routines.
    unrelated={
        "compute_AME_freq_shift_2par_spectrum":[0,0,0,0,0,2,1],
        "compute_AME_lognormal_width_spectrum":[0,0,0,0,1],
        "compute_magnetic_dust_spectrum":[0,0,0,2,1],
        "compute_physical_dust_spectrum":[0,0,1,0,2,4,1],
        "compute_lognormal_parabola_spectrum":[0,0,0,0,1],
        "compute_fixed_index_template_spectrum":[0,0,0,0,0,0,0,1],
        "compute_power_law_break_spectrum":[0,0,0,0,0,1],
        "compute_freefree_EM_spectrum":[0,0,0],
    }
    production+="\n"+ideal
    for name in ("ant2unit","ant2data","compute_bnu_prime_RJ_single","compute_bnu_prime_single",
                 "compute_ant2thermo_single","sz_thermo_single"):
        production+="\n"+procedure(bp,name)
    for name in ("spline_plain","splint_plain","splint_deriv"):
        production+="\n"+procedure(spline,name)
    production+="\n"+procedure(math,"tridag")+procedure(locate,"locate_dp")
    context=f"""module healpix_types
  implicit none
  integer, parameter :: dp=kind(1.d0),sp=kind(1.e0),i4b=kind(1),lgt=kind(.true.)
end module
module fossil_local_test_context
  use healpix_types
  implicit none
  integer, parameter :: num_recent_point=3,N_RMS_MARG=7,MPI_COMM_WORLD=0
  real(dp), parameter :: pi=acos(-1.d0)
{constants}
  real(dp) :: T_CMB=2.725d0
{types}
{band}
  type(fg_meta_data),allocatable :: fg_components(:)
  type(bandinfo),allocatable :: bp(:)
  integer :: myid_chain=0,numband=131
  logical :: fossil_exact_local_seds=.true.,fossil_exact_mjysr_units=.true.,fossil_si_cmb={'.true.' if si_cmb else '.false.'}
  interface spline
     module procedure spline_plain
  end interface
  interface splint
     module procedure splint_plain
  end interface
  interface locate
     module procedure locate_dp
  end interface
  interface sz_thermo
     module procedure sz_thermo_single
  end interface
contains
{production}
  subroutine mpi_abort(comm,code,ierr)
    integer,intent(in)::comm,code
    integer,intent(out)::ierr
    error stop 'Unexpected MPI abort in local SED fixture'
  end subroutine
  integer function getlun()
    getlun=99
  end function
  subroutine QuickSort_dp_dist(a,b)
    real(dp),intent(inout)::a(:),b(:)
    error stop 'Unexpected optional spectrum renormalization'
  end subroutine
  real(dp) function get_bp_line_ant(band,nu)
    integer,intent(in)::band
    real(dp),intent(in)::nu
    error stop 'Unexpected spectral-line conversion'
  end function
  real(dp) function splin2_full_precomp(x,y,table,a,b)
    real(dp),intent(in)::x(:),y(:),table(:,:,:,:),a,b
    error stop 'Unexpected spectral-grid interpolation in exact FOSSIL path'
  end function
  real(dp) function compute_CO_multiline_spectrum(band,bands,nu,p,deriv)
    integer,intent(in)::band,bands(:)
    real(dp),intent(in)::nu,p(:)
    logical,intent(in),optional::deriv
    error stop 'Unexpected Galactic CO component'
  end function
"""
    for name,ranks in unrelated.items():
        args=[f"a{i}" for i in range(len(ranks))]
        context+=f"\n  real(dp) function {name}({','.join(args)})\n"
        for arg,rank in zip(args,ranks):
            shape="("+",".join([":"]*rank)+")" if rank else ""
            context+=f"    real(dp),intent(in) :: {arg}{shape}\n"
        context+=f"    error stop 'Unexpected unrelated local SED: {name}'\n  end function\n"
    context+="end module fossil_local_test_context\n"
    return context


DRIVER=r"""program test_fossil_local
  use fossil_local_test_context
  implicit none
  character(len=1024) :: filename
  real(dp),pointer :: spectrum(:,:)
  real(dp) :: pars(2,5),result(5),derivative(8),nu,scale,old_scale,new_scale
  integer :: i,j,p,k,imax(1)
  call get_command_argument(1,filename)
  allocate(bp(131),fg_components(5))
  call read_spectrum(trim(filename),spectrum)
  spectrum(:,2)=spectrum(:,2)/maxval(spectrum(:,2))
  imax=maxloc(spectrum(:,2))
  pars(:,1)=[1.53d0,21.7d0]
  pars(:,2)=[-3.d0,0.d0]
  pars(:,3)=[-2.14d0,0.d0]
  pars(:,4)=[30.d0,0.d0]
  pars(:,5)=[33.35d0,0.d0]
  do j=1,5
     fg_components(j)%type='power_law'
     fg_components(j)%label='matched_local'
     fg_components(j)%amp_unit='MJy/sr'
     fg_components(j)%npar=2
     fg_components(j)%nu_ref=50.d9
     fg_components(j)%ref_band=1
     if(j==1) then
        fg_components(j)%type='one-component_dust'
        fg_components(j)%nu_ref=545.d9
        fg_components(j)%ref_band=34
     end if
     if(j>3) then
        fg_components(j)%type='AME_freq_shift'
        fg_components(j)%npar=1
        allocate(fg_components(j)%S_nu_ref(size(spectrum,1),3))
        fg_components(j)%S_nu_ref(:,1)=spectrum(:,1)
        fg_components(j)%S_nu_ref(:,2)=log(spectrum(:,2))
        call spline(fg_components(j)%S_nu_ref(:,1),fg_components(j)%S_nu_ref(:,2), &
             1.d30,1.d30,fg_components(j)%S_nu_ref(:,3))
        fg_components(j)%nu_peak=spectrum(imax(1),1)
     end if
     allocate(fg_components(j)%p_rms(2),fg_components(j)%priors(2,3))
     fg_components(j)%p_rms=0.d0
     fg_components(j)%nu_flat=100.d9
     fg_components(j)%frac_flat=0.d0
     fg_components(j)%priors(:,1)=[-6.d0,-.5d0]
     fg_components(j)%priors(:,2)=[3.d0,50.d0]
     if(j>3)then
        fg_components(j)%priors(1,1)=10.d0
        fg_components(j)%priors(1,2)=70.d0
     end if
  end do
  do i=1,131
     nu=(50.d0+15.d0*(i-1))*1.d9
     bp(i)%nu_c=nu
     bp(i)%unit='MJy/sr'
     bp(i)%id='delta'
     bp(i)%gain=1.d0
     bp(i)%use_color_corr=.false.
     bp(i)%a2t=compute_ant2thermo_single(nu)
     bp(i)%f2t=1.d0/compute_bnu_prime_single(nu)*1.d-14
  end do
  do i=1,131
     nu=bp(i)%nu_c
     k=0
     do j=1,5
        scale=ant2unit('MJy/sr',fg_components(j)%nu_ref,fg_components(j)%ref_band)
        result(j)=get_effective_fg_spectrum(fg_components(j),i,pars(:fg_components(j)%npar,j))/scale
        do p=1,fg_components(j)%npar
           k=k+1
           derivative(k)=get_effective_deriv_fg_spectrum(fg_components(j),i,p, &
                pars(:fg_components(j)%npar,j))/scale
        end do
     end do
     new_scale=ant2unit('MJy/sr',nu,i)
     fossil_exact_mjysr_units=.false.
     old_scale=ant2unit('MJy/sr',nu,i)
     fossil_exact_mjysr_units=.true.
     write(*,'(*(es26.17,1x))')nu/1.d9,result,derivative,old_scale/new_scale-1.d0
  end do
end program
"""


def run(package,require_compiler,si_cmb=False):
    reference=np.genfromtxt(package/"local_response.csv",delimiter=",",names=True)
    expected=np.column_stack([reference[name] for name in reference.dtype.names[1:]])
    assert expected.shape==(131,14), expected.shape
    compiler=shlex.split(os.environ.get("FC","gfortran"))
    if not compiler or shutil.which(compiler[0]) is None:
        if require_compiler:
            raise RuntimeError("Fortran compiler unavailable")
        print("SKIP: compiler unavailable; actual local SED procedures were not executed")
        return
    with tempfile.TemporaryDirectory(prefix="commander-fossil-local-") as tmp:
        work=Path(tmp)
        (work/"context.f90").write_text(build_context(si_cmb=si_cmb))
        (work/"driver.f90").write_text(DRIVER)
        executable=work/"local_sed_test"
        command=compiler+shlex.split(os.environ.get("FFLAGS",""))+[
            "-std=f2008","-ffree-line-length-none","-fcheck=all","-O0","-J",tmp,"-I",tmp,
            str(work/"context.f90"),str(work/"driver.f90"),"-o",str(executable)]
        subprocess.run(command,cwd=work,check=True)
        completed=subprocess.run([str(executable),str(package/"spectra/ame_emissivity_GHz.txt")],
                                 cwd=work,text=True,capture_output=True,check=True)
        actual=np.loadtxt(completed.stdout.splitlines())
    assert actual.shape==(131,15),actual.shape
    np.testing.assert_allclose(actual[:,:6],expected[:,:6],rtol=1e-12,atol=1e-30)
    # Both sides use the documented central finite-difference step.  Their
    # differencing roundoff is larger than the direct SED-value error.
    for column in range(6,14):
        norm=max(np.max(abs(expected[:,column])),1e-30)
        np.testing.assert_allclose(actual[:,column]/norm,expected[:,column]/norm,rtol=2e-8,atol=2e-9)
    legacy=float(np.float32(1e14))/1e14-1
    np.testing.assert_allclose(actual[:,14],legacy,rtol=0,atol=5e-16)
    relative=np.max(abs((actual[:,1:6]-expected[:,1:6])/np.maximum(abs(expected[:,1:6]),1e-300)))
    print(f"PASS: actual Commander direct local SEDs and MJy/sr forward units, 131 bands; max relative error={relative:.3e}")
    print("PASS: actual local shape derivatives, natural AME spline/peak convention, and preserved legacy unit path")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package",type=Path,help="Path to portable outputs/commander package")
    parser.add_argument("--require-compiler",action="store_true")
    parser.add_argument("--si-cmb",action="store_true",help="Check preserved local spectra with the SI CMB conversion enabled")
    args=parser.parse_args()
    run(args.package.resolve(),args.require_compiler,si_cmb=args.si_cmb)
