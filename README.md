# END KEM and Supporting Analysis Scripts

This repository contains the reference implementation of the END KEM and the supporting scripts used to reproduce the probability and failure-estimation calculations in the accompanying paper.

The repository has two main directories:

```text
.
├── scripts/
│   ├── dfr_frontier_checked.py
│   └── estimate_dec_failures.py
└── END KEM/
    ├── END-512/
    └── END-1024/
└── FCL-ML-KEM/
```

## Repository Layout

### `scripts/`

This directory contains the analysis scripts used for the decryption-failure and frontier calculations.

#### `dfr_frontier_checked.py`

This script checks the DFR and frontier probability calculations used in the NTRU-with-Encoding and NTRU-with-Trapdoor analysis.

It includes routines for:

- coefficient-level probability distributions;
- convolution and compression-noise distributions;
- FCL exclusion probability estimates;
- rejection-probability calculations;
- DFR checks for NTRU-with-Encoding, NTRU-with-Trapdoor, FCL-Kyber, NwE, and END-style parameter sets.

The script is intended for reproducing and checking the probability estimates used in the paper. Some cases are included as commented blocks in `main()` and can be enabled as needed.

Example:

```bash
cd scripts
sage -python dfr_frontier_checked.py
```

Depending on the local environment, this script may require SageMath, Python 3, and `mpmath`.

#### `estimate_dec_failures.py`

This script estimates the verification/collision term in the decryption-failure analysis. It intentionally separates this term from the other DFR contributions, such as:

- the probability that the rounding error is outside the allowed error set;
- the probability that FCL misses the correct candidate;
- the probability that a wrong candidate passes the verification check.

It supports built-in cases for:

- the trapdoor frontier example;
- `END-512`;
- `END-1024`;
- encoding-framework examples.

Example usage:

```bash
cd scripts
python3 estimate_dec_failures.py --case all
```

Run a single case:

```bash
python3 estimate_dec_failures.py --case end-t-512
python3 estimate_dec_failures.py --case end-t-1024
```

Optional flags:

```bash
python3 estimate_dec_failures.py --case all --signed-candidates
python3 estimate_dec_failures.py --case all --full-candidate-count
```

Dependencies:

```bash
pip install mpmath
```

## `END KEM/`

This directory contains the source code for the END KEM implementation.

Two parameter sets are provided:

- `END-512`
- `END-1024`

For each parameter set, two implementations are included:

- a standard C implementation;
- an AVX2-optimized implementation.

The detailed build and benchmark instructions are given in:

```text
END KEM/README.md
```

A typical implementation directory has the form:

```text
END KEM/END-512/END-512-257-C/
END KEM/END-512/END-512-257-AVX2/
END KEM/END-1024/END-1024-257-C/
END KEM/END-1024/END-1024-257-AVX2/
```

Because the directory name contains a space, quote the path when using shell commands if needed:

```bash
cd "END KEM/END-512/END-512-257-C/"
```

## Building and Running END

To build and run a specific END implementation, enter the corresponding implementation directory and use the provided `Makefile`.

For example, for the standard C implementation of `END-512`:

```bash
cd "END KEM/END-512/END-512-257-C/"
make main -j$(nproc)
./main
```

For memory profiling:

```bash
make memory_usage -j$(nproc)
valgrind --tool=massif --heap=yes --stacks=yes \
    --massif-out-file=memory_profile.out ./memory_usage
ms_print memory_profile.out
```

For timing-leakage checks:

```bash
make timecop -j$(nproc)
valgrind --track-origins=yes ./timecop 2>&1 | grep "Conditional"
```

The same workflow applies to the `END-1024` and AVX2 implementation directories.

## `FCL-ML-KEM/`
The implementation of FCL-ML-KEM in NIST-I security level.

## Requirements

The exact dependencies depend on which part of the repository is used.

For the analysis scripts:

- Python 3;
- `mpmath`;
- SageMath for scripts that use Sage-specific arithmetic.

For the END implementation:

- a C compiler such as `gcc` or `clang`;
- `make`;
- AVX2-capable hardware for the AVX2 implementations;
- `valgrind` for memory and timing-leakage checks.

## Reproducibility Notes

The scripts in `scripts/` are intended to reproduce the probability estimates and sanity checks used in the paper. The implementation in `END KEM/` is intended to reproduce the concrete performance and correctness experiments for the END KEM parameter sets.

The frontier calculations should be interpreted according to the assumptions stated in the paper, including the coefficient-independence model, the FCL exclusion estimate, and the specified verification/collision model.
