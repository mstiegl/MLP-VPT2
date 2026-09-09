# Validation data for the Fortran translation

Everything needed to check `gvpt2.x` against `step2/gvpt2.py` without
re-running MACE.

## inputs/

| files | what |
|-------|------|
| `C2417905_{opt,nma,cubics,quartics}.out` | step-1 output for the C2417905 example, produced by `step1_python/run` (MACE-OFF23 large, 9 atoms, 21 modes). These are the actual inputs of both codes. |
| `F_{opt,nma,cubics,quartics}.out` | the same data re-written in the file formats of `step1_fortran/vpt2.f90` (E24.12 constants, `Mode    N` headers, 21 frequencies only). |
| `S2_cubics.out` … `S10_cubics.out` | the C2417905 cubics multiplied by 2, 3, 4, 6, 10 to force many Fermi resonances (stress tests; geometry, modes and quartics are those of C2417905). |

## logs/

For every case, the complete standard output of the Python code
(`*_python.log`) and of the Fortran code (`*_fortran.log`).
`C2417905_python_with_intensity.log` and `C2417905_intensity_python.out`
are the Python run with `--enable_intensity`, which the Fortran does not do.

```
diff logs/C2417905_python.log logs/C2417905_fortran.log     -> empty
diff logs/F_python.log        logs/F_fortran.log            -> empty
```

For the `S*` cases the raw diff shows only the polyad basis order and
eigenvector signs; see `compare.sh` for the order-independent comparison,
and `../README.md` for the tabulated results.

## compare.sh

Rebuilds nothing; runs `../gvpt2.x` on all inputs and compares with the
stored Python logs (`./compare.sh`), or re-runs `gvpt2.py` too
(`./compare.sh --rerun-python`, needs a Python with ase, scipy, numpy).
