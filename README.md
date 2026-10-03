# Explore NTRU Decryption — ASIACRYPT 2026 Artifact

This repository is the artifact for **"How Compact Can NTRU Encryption Be? Heuristic Frontiers and Practical Schemes"** (IACR ePrint 2026/1715):
https://eprint.iacr.org/2026/1715

It contains the END KEM implementations, the FCL-ML-KEM implementation, and the analysis material used for decryption-failure/frontier and concrete-security calculations.

This artifact-preparation branch primarily adds documentation and reproducibility support. After the paper-to-code review, two narrowly scoped implementation fixes were also made to align the submitted code with ePrint 2026/1715: (1) the NIST-I FCL-ML-KEM tag is 128 bits (16 bytes), yielding the Table-6 688-byte ciphertext; and (2) END-512 uses the Table-2 trapdoor tuning value gamma=4 in both C and AVX2 key generation. The scientific-analysis files under `scripts/` are unchanged.

The ASIACRYPT 2026 artifact-evaluation requirements used to prepare this branch are at:
https://asiacrypt.iacr.org/2026/artifacts.php

## 1. Artifact contents

```text
.
├── END KEM/
│   ├── END-512/
│   │   ├── END-512-257-C/
│   │   └── END-512-257-AVX2/
│   └── END-1024/
│       ├── END-1024-257-C/
│       └── END-1024-257-AVX2/
├── FCL-ML-KEM/
├── scripts/
│   ├── dfr_frontier_checked.py
│   ├── estimate_dec_failures.py
│   ├── dependencies_impact_on_wrap_errors.ipynb
│   └── concrete_security_estimates.ipynb
├── artifact/
│   ├── reproduce_dfr.py
│   ├── run_notebooks.sh
│   └── setup_estimators.sh
├── Dockerfile
├── LICENSE
├── THIRD_PARTY.md
├── AI_USE_DISCLOSURE.md
└── SUBMISSION_CHECKLIST.md
```

## 2. Paper-to-artifact mapping

The mapping below is checked against IACR ePrint 2026/1715.

| Paper result / discussion | Artifact | Reproduction route / status |
| --- | --- | --- |
| Section 2.4 coefficient-independence heuristic and the 100-trial END dependency stress tests | `scripts/dependencies_impact_on_wrap_errors.ipynb` | `bash artifact/run_notebooks.sh dependencies` |
| Section 3 FCL probability machinery | `scripts/dfr_frontier_checked.py` | used by the DFR/frontier reproduction wrapper |
| Section 5.1.2 NTRU-with-Encoding frontier | `scripts/dfr_frontier_checked.py` | `sage -python artifact/reproduce_dfr.py --case ntru-encoding-frontier` |
| Section 5.2.2 NTRU-with-Trapdoor frontier | `scripts/dfr_frontier_checked.py` | `sage -python artifact/reproduce_dfr.py --case ntru-trapdoor-frontier` |
| Section 6.1, Algorithms 4–14: END decoding/PKE/KEM | `END KEM/END-512/`, `END KEM/END-1024/` | implementation source; see `CODE_REVIEW.md` for paper/code consistency findings |
| Table 2: suggested END parameters | END `param.h`, `keygen.c`, `pke.c` | encoded PK/CT sizes match; END-512 C/AVX2 now use gamma=4 as in Table 2; the remaining sampler/candidate-family findings are documented in `CODE_REVIEW.md` |
| Table 3: END DFR components | `scripts/dfr_frontier_checked.py`, `scripts/estimate_dec_failures.py` | `sage -python artifact/reproduce_dfr.py --case all` and `sage -python scripts/estimate_dec_failures.py --case all`; candidate-family consistency with C code must be resolved before claiming full reproduction |
| Table 4: END sizes and concrete security | END `param.h` / codecs; `scripts/concrete_security_estimates.ipynb` | sizes are directly inspectable; security estimates via `bash artifact/setup_estimators.sh && bash artifact/run_notebooks.sh security` |
| Table 5: END REF/AVX2 cycle counts | END `main.c` + Makefiles | END rows can be rerun locally; external baseline rows are not self-contained in this repository |
| Section 6.3 timing-independence discussion | END `timecop.c`, `poison.h` | build the `timecop` target; see `CODE_REVIEW.md` regarding optimization-level reproducibility |
| Appendix A, Algorithms 15–16 and Table 6: FCL-ML-KEM | `FCL-ML-KEM/` | NIST-I (`KYBER_K=2`) now appends a 16-byte tag and has a 688-byte ciphertext, matching Algorithm 15 and Table 6; the inherited K=3/4 targets are not paper claims |

### 2.1 Paper/code consistency status

A paper-to-code review was performed against ePrint 2026/1715. The FCL-ML-KEM NIST-I tag-length mismatch and the END-512 gamma mismatch identified in that review have been fixed on this branch. Two findings intentionally remain unresolved: the END fixed-weight-vs-iid secret-sampling difference and the FCL shortlist/candidate-family difference between analysis and implementation. See `CODE_REVIEW.md` for details.

Because the implementation changed, Table 5 END-512 performance and Table 6 FCL-ML-KEM performance must be remeasured before archival; do not reuse pre-fix cycle counts without verification.

## 3. Reference environment and exact dependencies

The two notebooks record a **SageMath 10.1** kernel and **Python 3.10.12**. The artifact therefore uses SageMath 10.1 as the reference mathematical environment.

Pinned mathematical/software dependencies:

- SageMath: **10.1**
- Python in the committed notebooks: **3.10.12**
- mpmath: **1.3.0** (the version packaged by SageMath 10.1)
- gmpy2: **2.1.2** (the version packaged by SageMath 10.1)
- tqdm: **4.67.1** in the artifact container/setup
- standard lattice-estimator: `malb/lattice-estimator@3e48ef421ec256afddb3e7d2249a77eab6e9ba12`
- enhanced lattice-estimator: `identitymapping/enhanced_lattice-estimator@876b66173f4354a96ddafc0ce3a79767ec43c6d4`
- PrimalMeetLWE: `yonghaason/PrimalMeetLWE@61115115830c909e42758f2774606074bf98afb1`

For the C implementations the reference container installs GCC, Clang, GNU Make, Valgrind, Git, and OpenSSL development headers. Performance results should always be reported together with the CPU model, OS, compiler version and Git commit.

AVX2 implementations require an x86-64 CPU exposing the `avx2` flag. Non-AVX2 machines can still build/run the portable C implementations.

## 4. Recommended: Docker reproduction

Docker gives reviewers the closest single-command environment to the notebook metadata.

```bash
git clone https://github.com/Icarid-Liu/Explore_NTRU_Decryption.git
cd Explore_NTRU_Decryption
git checkout artifact-asiacrypt-2026

docker build -t explore-ntru-artifact .
docker run --rm -it explore-ntru-artifact bash
```

Inside the container, run the commands from Sections 6–8 directly. The artifact intentionally keeps the execution interface small instead of wrapping every experiment in additional shell scripts.

## 5. Native installation

A native Linux environment can be used instead of Docker. Install:

- SageMath 10.1;
- GCC and Clang;
- GNU Make;
- Git;
- Valgrind (only for memory/timing-leakage checks);
- OpenSSL development headers (for the FCL-ML-KEM NIST KAT targets).

Then install the one additional notebook package and fetch exact estimator revisions:

```bash
sage -pip install tqdm==4.67.1
bash artifact/setup_estimators.sh
```

The setup script clones third-party estimator repositories into `scripts/` at the exact commits recorded in the committed concrete-security notebook. Those repositories are not vendored into this artifact.

## 6. Analysis reproduction

### 6.1 Verification/collision term

This script is pure Python + mpmath and is the lightest numerical reproduction target:

```bash
sage -python scripts/estimate_dec_failures.py --case all
```

Individual cases:

```bash
sage -python scripts/estimate_dec_failures.py --case end-t-512
sage -python scripts/estimate_dec_failures.py --case end-t-1024
sage -python scripts/estimate_dec_failures.py --case trapdoor-frontier
sage -python scripts/estimate_dec_failures.py --case encoding
```

Optional sensitivity switches already present in the original script:

```bash
sage -python scripts/estimate_dec_failures.py --case all --signed-candidates
sage -python scripts/estimate_dec_failures.py --case all --full-candidate-count
```

Interpretation: the script reports the relevant terms as base-2 logarithms, e.g. `Pr[Dec. fails] = 2^(...)`. More negative exponents indicate smaller estimated failure probability. The optional switches deliberately change the candidate-count convention and should not be mixed with the paper's default configuration unless performing a sensitivity check.

### 6.2 DFR/frontier calculations

The original `scripts/dfr_frontier_checked.py` keeps the paper cases as commented blocks in `main()`. To make those calculations exercisable without editing that source file, `artifact/reproduce_dfr.py` imports the original functions and runs the committed END-512/END-1024 blocks:

```bash
sage -python artifact/reproduce_dfr.py --case all
```

The wrapper contains orchestration/parameter selection only; the probability-distribution and DFR computations remain in the original source.

### 6.3 Notebooks

Fetch exact external estimators first:

```bash
bash artifact/setup_estimators.sh
```

Then execute either notebook non-interactively:

```bash
bash artifact/run_notebooks.sh dependencies
bash artifact/run_notebooks.sh security
```

or both:

```bash
bash artifact/run_notebooks.sh all
```

Executed notebooks are written to `artifact-output/notebooks/`. The concrete-security notebook prints the three estimator commit hashes at the beginning; they should exactly match the pins listed above.

## 7. END KEM build and output interpretation

Example: portable END-512:

```bash
cd "END KEM/END-512/END-512-257-C"
make main -j"$(nproc)"
./main
```

Portable END-1024 and the AVX2 directories use the same `make main` interface.

The program first prints `END.KEM.Enc Shared Key` and `END.KEM.Dec Shared Key`. The two arrays are the correctness check and should match. It then prints average cycle counts for KeyGen, Enc and Dec over the implementation's built-in benchmark loop. Cycle counts are machine dependent; report the output of `artifact/environment_report.sh` with any performance numbers.

Memory profiling:

```bash
make memory_usage -j"$(nproc)"
valgrind --tool=massif --heap=yes --stacks=yes \
  --massif-out-file=memory_profile.out ./memory_usage
ms_print memory_profile.out
```

Timing-leakage instrumentation target:

```bash
make timecop -j"$(nproc)"
valgrind --track-origins=yes ./timecop 2>&1 | grep "Conditional"
```

These Valgrind-based checks are diagnostic; reviewer conclusions should record the Valgrind version and platform.

## 8. FCL-ML-KEM build and output interpretation

```bash
cd FCL-ML-KEM
make test -j"$(nproc)"
./test/test_kyber512

make speed -j"$(nproc)"
./test/test_speed512

make nistkat -j"$(nproc)"
./nistkat/PQCgenKAT_kem512
```

Equivalent targets exist for 768 and 1024 because the inherited build system exposes those parameterizations. The project-specific paper claim should be matched to the parameterization stated in the paper. Cycle-count output is hardware dependent.

## 9. Source-code organization and modifications

The scientific-analysis files under `scripts/` are not modified. Existing cryptographic source changes are intentionally limited to the two paper-alignment fixes authorized after review:

- `FCL-ML-KEM/params.h`, `indcpa.c`, `api.h`, and `test/test_kyber.c`: use a 16-byte FCL tag for the paper's NIST-I (`KYBER_K=2`) instantiation and assert the 688-byte Table-6 ciphertext size; inherited K=3/4 targets retain their previous 32-byte tag behavior.
- `END KEM/END-512/END-512-257-C/keygen.c` and `END KEM/END-512/END-512-257-AVX2/keygen.c`: use the paper's gamma=4 tuning coefficient in `compute_w()`.

New material under `artifact/`, the Dockerfile, and the documentation is packaging/reproducibility infrastructure. `artifact/reproduce_dfr.py` invokes existing functions with parameter blocks already present as comments in the original DFR script; it does not replace the scientific implementation.

## 10. Licensing and third-party material

See `LICENSE` for the license applied to original artifact contributions and `THIRD_PARTY.md` for third-party provenance and retained upstream terms. Source files containing their own license notices keep those notices and terms.

Before archival, the authors should perform the provenance review in `SUBMISSION_CHECKLIST.md`, especially for any unannotated implementation/assembly files whose origin is known to the authors but not stated in the source.

## 11. AI-use disclosure

Documentation, artifact-support material, and the two explicitly authorized paper-alignment implementation fixes were produced with generative-AI assistance and author direction. The disclosure is in `AI_USE_DISCLOSURE.md`. The scientific-analysis files under `scripts/` were not edited.

## 12. Archival packaging

After author review, freeze an immutable Git tag/commit. If a tarball is needed, Git itself is sufficient:

```bash
git archive --format=tar.gz --prefix=Explore_NTRU_Decryption/ -o artifact.tar.gz HEAD
sha256sum artifact.tar.gz
```

For the HotCRP artifact submission, provide the public repository URL together with the immutable commit or tag. See `SUBMISSION_CHECKLIST.md` for the remaining author-only checks.
