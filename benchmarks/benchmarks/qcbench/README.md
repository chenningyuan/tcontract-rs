# QC-realistic DF-RCCSD contraction benchmark (dfrccsd demo)

`qcbench` reuses `tcbench`'s execution engines, timing, TBLIS adapter, thread enforcement, and evidence format. It replaces isotropic TCCG sizing with fixed occupied/virtual/auxiliary extents at quantum-chemistry aspect ratios. Every operand is contiguous column-major, `alpha = 1`, `beta = 0`, and input generation and plan construction are outside the timed call.

## Corpus

The six cases are the two leading formal-scaling rungs of the DF-RCCSD energy/amplitude route from U. Bozkaya and C. D. Sherrill, *J. Chem. Phys.* **147**, 044104 (2017), DOI [10.1063/1.4994918](https://doi.org/10.1063/1.4994918), main Eq. (5) and supplementary Eqs. (10)–(43): the `o²v⁴` head `T2.vvvv`, and the `o³v³` ring terms `T2.ring_mi`, `T2.ring_mj`, `T2.ring_u`, `W.mbej.t2`, `W.mbje.t2`.

The three extent sets are aspect-ratio classes, not molecules. The ratios o/v = 1/3, 1/10, 1/20 are the median v/o of the cc-pVDZ, cc-pVTZ, and cc-pVQZ rungs of a 5,239-molecule NAO statistics panel (3.90, 10.17, 20.36), so each ratio stands for one basis-quality rung without pinning a molecule. Family IDs are `<ratio>/<case-id>`, giving 18 families:

| Ratio | Rung median it stands for | o | v | naux |
|---|---|---:|---:|---:|
| `ov_03` | cc-pVDZ (v/o 3.90, naux/nao 3.747) | 18 | 54 | 270 |
| `ov_10` | cc-pVTZ (v/o 10.17, naux/nao 2.568) | 10 | 100 | 282 |
| `ov_20` | cc-pVQZ (v/o 20.36, naux/nao 2.237) | 7 | 140 | 329 |

Each ratio's extents follow the rung-median rule `naux = round(k_rung · (o+v))`, and `v·o` stays nearly constant (972 / 1000 / 980): a larger basis rung affords more virtuals but fewer occupied orbitals. The extents are derived from a *whole-binary* memory budget: the largest f64 A+B+C footprint over the planned 25-case suite this demo belongs to — bound by the `o³v³` triples output of its perturbative-triples cases — stays within 8 GiB per cell, so no case is ever skipped for memory. This member carries the dfrccsd corpus only, as a reviewable sample of the pattern; the full suite (eight QC corpora, 25 cases) follows separately.

## Engines and correctness

- `plan`: tcontract's planner, with the same `TCBENCH_*` knobs as `tcbench`.
- `tblis`: actual C++ TBLIS through the direct FFI adapter (`--features tblis`).
- The forced packed route and TTGT remain available to the shared verification implementation but are not published qcbench engine rows.

`verify` first checks an all-ones known-value contraction, then evaluates all selected families. At realistic extents it compares full outputs from the planner, forced packed route, and built baselines using relative Frobenius error, requiring finite error at most `1e-10` for f64. It is a cross-provider check, not a naive-loop recomputation of the expensive QC contractions. Timed runs also compare each later engine with the first engine's full output and fail on `MISMATCH`.

Only f64 and unit batch are accepted. `--n 1` is solely the timing-v2 identity sentinel; it does not resize the tensors. Layout stress is fixed to `none`.

## Build and run

Use the same native release profile and TBLIS prefix for every arm:

```sh
export TBLIS_ROOT=/path/to/tblis-v2.0-beta2
export LD_LIBRARY_PATH="$TBLIS_ROOT/lib:$TBLIS_ROOT/lib64"
export RUSTFLAGS='-C target-cpu=native'
export CARGO_BUILD_JOBS=4
cargo build --release -p tcontract-bench --bin qcbench --features tblis
```

Correctness precedes timing. Pin every process, including 1T:

```sh
benchmarks/scripts/pinned.sh 0 -- target/release/qcbench verify \
  --threads 1 --dtype f64 --engines plan,tblis
benchmarks/scripts/pinned.sh 0-3 -- target/release/qcbench verify \
  --threads 4 --dtype f64 --engines plan,tblis
```

A small timing smoke test, not a published result:

```sh
benchmarks/scripts/pinned.sh 0 -- target/release/qcbench run \
  --threads 1 --case ov_10 --dtype f64 --engines plan,tblis \
  --prime-ms 1500 --reps 2
```

Campaign runs belong in `tcontract-benchmark`, not this source checkout. The legacy recorder accepts and labels an uncommitted harness as dirty evidence; timing-v2 requires clean committed sources and therefore cannot record a review working tree.
