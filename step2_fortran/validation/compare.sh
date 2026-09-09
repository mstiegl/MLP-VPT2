#!/bin/bash
# Re-run the Fortran gvpt2.x (and optionally the Python gvpt2.py) on the
# validation inputs and compare the outputs.
#
#   ./compare.sh                 run gvpt2.x, compare with the stored Python logs
#   ./compare.sh --rerun-python  also re-run ../../step2/gvpt2.py (needs ase, scipy, numpy)
#
# Cases:
#   C2417905  the real example (MACE-OFF23 large, 21 modes)       -> full diff must be empty
#   F         same data in the step1_fortran file formats          -> full diff must be empty
#   S2..S10   cubics scaled by 2,3,4,6,10 (many resonances)        -> order-independent
#             comparison: Fermi list, eigenvalues, (root,state,coef^2) triples, final tables
#
# For the S cases the raw full-file diff is NOT empty: Python enumerates the
# polyad basis in set (hash) order and eigenvector signs are arbitrary.

set -u
here=$(cd "$(dirname "$0")" && pwd)
exe="$here/../gvpt2.x"
if [ ! -x "$exe" ]; then echo "build first: (cd $here/.. && make)"; exit 1; fi
work=$(mktemp -d)
cp "$here"/inputs/*.out "$work"/
cp "$here"/../../step2/gvpt2.py "$work"/
cd "$work"
status=0

tail_from_harmonic() { sed -n '/^Harmonic ω/,$p' "$1"; }
head_to_polyads()    { awk '/GVPT2 Polyad Prints/{exit} {print}' "$1" | sed '$d'; }
roots()              { grep '^Root' "$1" | sort; }
ci_triples()         { awk '/^Root/{r=$0} /^   \(/{print r, $1, $NF}' "$1" | sort; }

check() { # label file_a file_b
  if diff -q "$2" "$3" >/dev/null; then echo "    $1: identical"; else echo "    $1: DIFFERENT"; diff "$2" "$3" | head -20; status=1; fi
}

for mol in C2417905 F S2 S3 S4 S6 S10; do
  case $mol in S*) for suf in opt nma quartics; do cp C2417905_$suf.out ${mol}_$suf.out; done;; esac
  "$exe" --mol_name=$mol --num_modes=21 > ${mol}_fortran.log 2>&1 || { echo "$mol: gvpt2.x failed"; status=1; continue; }
  if [ "${1:-}" = "--rerun-python" ]; then
    python gvpt2.py --mol_name=$mol --num_modes=21 > ${mol}_python.log 2>&1 || { echo "$mol: gvpt2.py failed"; status=1; continue; }
  else
    cp "$here"/logs/${mol}_python.log .
  fi
  echo "$mol: $(grep Detected ${mol}_python.log)"
  case $mol in
    S*)
      head_to_polyads ${mol}_python.log > ph; head_to_polyads ${mol}_fortran.log > fh; check "header + Fermi list" ph fh
      roots ${mol}_python.log > pr; roots ${mol}_fortran.log > fr;                     check "polyad eigenvalues" pr fr
      ci_triples ${mol}_python.log > pc; ci_triples ${mol}_fortran.log > fc;           check "(root, state, coef^2)" pc fc
      tail_from_harmonic ${mol}_python.log > pt; tail_from_harmonic ${mol}_fortran.log > ft; check "final tables" pt ft
      ;;
    *)
      check "full output" ${mol}_python.log ${mol}_fortran.log
      ;;
  esac
done
rm -rf "$work"
[ $status -eq 0 ] && echo "ALL COMPARISONS PASSED" || echo "SOME COMPARISONS FAILED"
exit $status
