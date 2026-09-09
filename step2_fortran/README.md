# step2_fortran: DVPT2 + GVPT2 in Fortran

Fortran translation of `step2/gvpt2.py`. It reads the same four files
produced by step 1 (`{MOL}_opt.out`, `{MOL}_nma.out`, `{MOL}_cubics.out`,
`{MOL}_quartics.out`, from either `step1_fortran` or `step1_python`) and
prints the same report: harmonic frequencies, deperturbed anharmonic
constants, VPT2 / DVPT2 / GVPT2 fundamentals, the energies of all states
in the resonant polyads, chi0 and the zero-point energies.

The optional double-harmonic IR/Raman intensities of `gvpt2.py`
(`--enable_intensity`) need the MACE-MDP PyTorch model and are not part of
the Fortran version; the flag is accepted but only prints a notice.

## Build

Requires a Fortran 2008 compiler and LAPACK (only `dsyev` is used).

```bash
make                                  # gfortran; Accelerate on macOS, -llapack -lblas elsewhere
make FC=ifx LIBS="-qmkl=sequential"   # Intel compiler + MKL
make clean
```

Do not add `-r8`: the code declares its own double-precision kind.

## Run

Copy the `*.out` files from step 1 into this directory, then

```bash
./gvpt2.x --mol_name=C2417905 --num_modes=21     # same arguments as step2/run
./gvpt2.x C2417905                               # positional; modes default to 3N-6
```

`./run` contains the first form.

## Files

| file | contents |
|------|----------|
| `gvpt2_constants.f90` | precision kind, unit conversions, resonance thresholds, ASE atomic-mass table |
| `gvpt2_io.f90` | readers for the four input files; Python/numpy-style number formatting |
| `gvpt2_polyad.f90` | states, interactions, polyads, LAPACK diagonalisation, Hungarian assignment |
| `gvpt2_core.f90` | rotational constants, Coriolis zeta, Q->q scaling, Fermi detection, chi, chi0, ZPVE, GVPT2 driver |
| `gvpt2.f90` | main program: argument parsing, workflow, report |

## Correspondence with gvpt2.py

| gvpt2.py | Fortran |
|----------|---------|
| `read_freqs_and_vecs`, `parse_cubic_tensor`, `parse_quartic_tensor`, `ase.io.read` | `gvpt2_io` |
| `rotational_constants`, `coriolis_zeta`, `convert_tensors_Q_to_dimensionless_q` | `gvpt2_core` |
| `identify_fermi_resonances`, `compute_chi_deperturbed_cm`, `fundamentals_from_chi_cm` | `gvpt2_core` |
| `compute_chi0_cm`, `compute_zpve_cm`, `build_explicit_state_energies`, `apply_gvpt2_fermi` | `gvpt2_core` |
| `State`, `Interaction`, `Polyad`, `fermi_solver` | `gvpt2_polyad` |
| `numpy.linalg.eigh` | LAPACK `dsyev` |
| `scipy.optimize.linear_sum_assignment` | `hungarian` in `gvpt2_polyad` |
| `atoms.get_masses()` (ASE) | `ase_mass` table in `gvpt2_constants` |

Design notes:

- The Fermi set is a logical array `fermi(k, i, j)`, true when
  `(k, (i, j))` is in the Python set. Fermi-2 entries are stored for both
  `(i, j)` and `(j, i)`, as `itertools.permutations` produces them.
- `fermi_solver` reproduces the Python grouping exactly, including the
  case where an interaction is added to more than one existing polyad and
  later polyads overwrite the energies of shared states.
- The unit-conversion constants and the atomic masses are those of
  `gvpt2.py` and ASE, not those of `step1_fortran/constants.f90`, so the
  rotational constants match the Python output.

## Procedure used to produce and check the translation

1. **Set up the environment.** Cloned MLP-VPT2 and created a Python 3.12
   virtual environment with `mace-torch` 0.3.16, `ase` 3.29, `scipy` and
   `numpy` (uv, Apple M2 Max, macOS).
2. **Generated the step-1 input files.** Ran `step1_python/run` unchanged
   (MACE-OFF23 large, geometry `C2417905.inp`, 9 atoms, 21 modes). Wall time
   2 min 43 s including the one-time model download. This produced
   `C2417905_opt.out`, `_nma.out`, `_cubics.out` and `_quartics.out`.
3. **Captured the reference output.** Copied the four files to `step2` and
   ran `gvpt2.py` as shipped, once with `./run` (energies only, 0.8 s) and once
   with `--enable_intensity` (4.3 s) to confirm both code paths work. The
   energies-only log is the correctness oracle for the translation.
4. **Read the Python line by line** and listed every external dependency
   that has no Fortran equivalent: `ase.io.read` (xyz parser and atomic
   masses), `numpy.linalg.eigh`, `scipy.optimize.linear_sum_assignment`,
   `itertools` set/permutation semantics, and the MACE-MDP model. Also noted
   the behaviours that must be reproduced exactly: last-block frequency
   parsing, taking the last `n_modes` entries, storing Fermi-2 pairs in both
   index orders, the `fermi_solver` grouping rule (an interaction joins every
   polyad containing one of its states) and the dict-update overwrite when
   polyads share a state.
5. **Decided the scope with the author:** target language Fortran; energies
   only, because the intensity block needs a PyTorch runtime.
6. **Wrote the Fortran in five modules** (`gvpt2_constants`, `gvpt2_io`,
   `gvpt2_polyad`, `gvpt2_core`, program `gvpt2`), keeping the Python
   function boundaries, operation order and printout format so the two
   programs can be compared with `diff`. The ASE mass table was generated
   from `ase.data.atomic_masses` with round-trip-exact literals. The
   Hungarian algorithm replaces SciPy's assignment; LAPACK `dsyev` replaces
   `numpy.linalg.eigh`. Python's `repr()` and numpy's array printing were
   reimplemented for the diagnostic lines so those match too.
7. **Built and diffed.** The first build differed from the Python only in
   two cosmetic formatting bugs (a lost separator space in the array print
   and the exponent threshold in the shortest-repr routine). After fixing
   them the output was byte-identical.
8. **Stress-tested the resonance machinery**, which the real molecule barely
   exercises (3 candidates, two 2-state polyads): scaled the cubic constants
   by 2, 3, 4, 6 and 10 to force up to 77 candidates, 17 polyads and 15-state
   Hamiltonians, and compared Python and Fortran. This exposed one real
   portability issue, fixed in step 9.
9. **Fixed field overflow.** Python's `{:12.6f}` widens the field when a
   number does not fit; Fortran `F12.6` prints asterisks. Added `pyf()`,
   which emulates the Python behaviour, and used it for all table output.
10. **Verified the assignment solver** against
    `scipy.optimize.linear_sum_assignment` on 400 random matrices.
11. **Verified the parsers** on files re-written in the exact formats
    produced by `step1_fortran/vpt2.f90`, and checked the positional and
    space-separated argument forms, the default mode count, the
    `--enable_intensity` notice and the error paths (missing file, missing
    argument, too many modes requested).
12. **Compiled with strict checks** (`-std=f2008 -pedantic -Wall -Wextra
    -fcheck=all -ffpe-trap=invalid,zero,overflow`) and re-ran the real case
    and all stress cases under that build.

## Validation results

All comparisons below are against `step2/gvpt2.py` run on the same input
files in the same directory.

### 1. C2417905 (MACE-OFF23 large, 21 modes, 3 Fermi resonances)

```
$ diff gvpt2_python.log gvpt2_fortran.log
(no output)                     exit code 0
MD5 python  = bb41e11bdb60426aa6a3d2b4a0cab0d8
MD5 fortran = bb41e11bdb60426aa6a3d2b4a0cab0d8      143 lines each
```

The complete output, every line including the diagnostic tensor statistics
and the polyad printout, is byte-identical.

### 2. Stress tests: cubic constants scaled by 2 to 10

Scaling the cubics multiplies the resonance criterion
`phi^4 / (64 d^3)` by the fourth power of the factor, so many more
Fermi candidates pass the threshold. The physics becomes meaningless
(negative fundamentals at scale 6 and above) but the code paths are the
same, and the numbers must still agree.

A raw full-file `diff` is not empty for these cases. Every differing line
is in the polyad printout and is one of two cosmetic effects: the order of
the basis states (Python enumerates a `set`, in hash order; the Fortran uses
insertion order) or the overall sign of an eigenvector. Excerpt of the
complete diff at scale 3:

```
43,44c43,44
<   0: (3, 3)            python: set order
<   1: (5,)
---
>   0: (5,)              fortran: insertion order
>   1: (3, 3)
55,56c55,56
<    (3, 3)         -0.999440     0.998880      same rows, sign flipped,
<    (5,)           -0.033462     0.001120      coef**2 unchanged
---
>    (5,)            0.033462     0.001120
>    (3, 3)          0.999440     0.998880
```

Comparing everything that is order- and sign-independent, `diff` was empty
in every case:

| scale | Fermi candidates | polyads | largest polyad | header + Fermi list | eigenvalues | (root, state, coef²) triples | final tables |
|------:|-----------------:|--------:|---------------:|---------------------|-------------|------------------------------|--------------|
| 2  | 11 | 6  | 3  | identical | 13 identical | 29 identical  | 95 lines identical  |
| 3  | 24 | 10 | 4  | identical | 26 identical | 72 identical  | 105 lines identical |
| 4  | 30 | 12 | 4  | identical | 31 identical | 85 identical  | 110 lines identical |
| 6  | 53 | 13 | 11 | identical | 53 identical | 291 identical | 122 lines identical |
| 10 | 77 | 17 | 15 | identical | 84 identical | 646 identical | 136 lines identical |

"Final tables" means everything from `Harmonic ω` to the end: harmonic
frequencies, deperturbed χ_ii, the fundamentals table (harmonic, VPT2,
DVPT2, GVPT2 and both shifts), the solved state energies of every polyad,
chi0 and the ZPVE block in three units. Scales 3 to 10 contain states that
belong to more than one polyad (3, 3, 10 and 18 such states), which
exercises the overwrite rule in `fermi_solver`. At scale 6, for example,
mode 8 moves by 243.9 cm⁻¹ between VPT2 and GVPT2 and the Fortran
reproduces every entry to all printed digits.

Before the `pyf()` fix (procedure step 9) scale 10 printed asterisks where
Python printed `-11666.534167`; the numbers themselves agreed.

### 3. Input files in the step1_fortran formats

The C2417905 data was re-written exactly as `step1_fortran/vpt2.f90` writes
it: `E24.12` force constants (`0.803149074200E-13`), headers `Mode    1`
from `(A,I5)`, vectors in `3F20.10`, only the 21 vibrational frequencies
listed, and the geometry with an `I2` atom count and an `Angstroms` comment
line.

```
diff python vs fortran on these files:  identical (143 lines)
```

Against the original 17-digit run only the two "RAW TENSOR MAGNITUDES"
lines change, because the inputs now carry 12 significant digits; every
results table is unchanged.

### 4. Command-line forms

`./gvpt2.x F` (positional, default `3N-6` modes) and
`./gvpt2.x --mol_name F --num_modes 21` (space-separated values) give output
identical to `./gvpt2.x --mol_name=F --num_modes=21`. `--enable_intensity`
prints a notice and continues. A missing file, a missing `--mol_name` and a
request for more modes than the file contains stop with exit code 1 and a
message on stderr.

### 5. Hungarian assignment versus SciPy

The Fortran `hungarian` routine, which replaces
`scipy.optimize.linear_sum_assignment`, was run on 400 random square cost
matrices of size 1 to 12: one third squared orthonormal eigenvector
overlaps (the polyad case), one third Gaussian, one third small integers
with many ties.

```
mismatches in optimal cost: 0     max |cost_fortran - cost_scipy| = 0.00e+00
```

### 6. Strict compiler checks

Built with
`-std=f2008 -pedantic -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow`
(gfortran 16.1, Accelerate LAPACK). The only warnings are seven intentional
exact comparisons of reals (zero tests and the round-trip test in the
number formatter). Under this bounds-checked, trap-enabled build the real
case and all five stress cases run to completion with exit code 0 and match
the Python tables.

### Not validated

The Intel `ifx` + MKL build path in the Makefile has not been tried; no
Intel compiler was available on the development machine.
