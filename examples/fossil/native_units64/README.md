# FOSSIL native-unit 64-channel configurations

These complete English `.txt` parameter files are copied from the validated
`outputs/commander_native_units` package. Transfer that entire package as
`fossil_inputs_native_units` under the cluster project. Each chain uses
64 MPI ranks with one process per band and one thread per rank.

Band data and RMS are `uK_cmb`; Galactic reference amplitudes are `uK_ant`
(the unchanged internal RJ convention); CMB reference amplitudes are
`uK_cmb`. Spectral parameters and physical global amplitudes keep their units.
The source foreground/free-free module is unchanged from commit8ffea73.
The new physical-global interface converts spectra to the common band unit.
Original master simulations remain MJy/sr and noise levels are unchanged.

Update and recompile Commander before using this suite. Start with case01 at
NSIDE64 and inspect the diagnostic notebook. Then case03 and case02 validate
coupled amplitudes and conditional sampling. Case04 tests global shapes;
case05 isolates amplitude sampling; case06 optimizes the full model;
case07 and the three case09 seeds provide full forecast chains. Case08
exercises disabled globals. Case12 requires the completed case11 checkpoint;
compare their concatenation with case02. Other files are independent runs.

Fits are full sky. GAL070 is only a postprocessing region. The TT prior is
fixed; CMB realizations use ell>=2. Instrumental monopoles, gains and bandpass
shifts are fixed. Fresh chain roots contain `native_units_fullsky` and do not
reuse previous failed-run directories.

Use your own launcher or the supplied `comm_run_native.sh`. The updated
notebook in `diagnostics/` reads verified conversion factors and displays
native component outputs and residuals in MJy/sr. Global spectral snapshots
remain MJy/sr; native trace amplitudes retain their physical meaning.
