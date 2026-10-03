# Paper-to-Code Review — ASIACRYPT 2026 Artifact

Paper reviewed: **How Compact Can NTRU Encryption Be? Heuristic Frontiers and Practical Schemes**, IACR ePrint 2026/1715.

Paper URL: https://eprint.iacr.org/2026/1715

Review scope: compare the public paper against the current repository, with emphasis on the concrete END implementation, FCL-ML-KEM Appendix-A implementation, DFR/security scripts, and artifact reproducibility. Following the review, the authors explicitly requested implementation fixes for B1 and B2; those two fixes are now applied on the artifact branch. B3 and B4 remain intentionally unresolved.

## Status

Two originally blocking mismatches have been **resolved on this branch**: B1 (FCL-ML-KEM tag/ciphertext size) and B2 (END-512 gamma). Two paper/code mismatches remain intentionally unresolved: B3 (secret sampling distribution) and B4 (FCL candidate family). Performance numbers affected by B1/B2 must be remeasured before archival.

## Resolved findings

### B1. RESOLVED — FCL-ML-KEM NIST-I now uses the paper's 16-byte tag

Paper:
- Appendix A, Algorithm 15 appends `Hash(m)[0:lambda/8]`.
- NIST-I uses `lambda = 128`, so the tag is 16 bytes.
- Table 6 reports an FCL-ML-KEM-512 ciphertext of **688 bytes**.

Original pre-fix code:
- `FCL-ML-KEM/params.h` defined `KYBER_SYMBYTES = 32` and included a full `KYBER_SYMBYTES` in `KYBER_INDCPA_BYTES`.
- For `KYBER_K=2`, the compressed MLWE components occupied `2*288 + 96 = 672` bytes, so the implementation ciphertext was **672 + 32 = 704 bytes**, not 688.
- `FCL-ML-KEM/indcpa.c` wrote and verified the full SHA3-256 output rather than the first 16 bytes.

Impact before the fix:
- Table 6's primary compactness claim could not be reproduced by that code.
- The implementation was not the exact Algorithm-15/16 instantiation described by the paper.

Resolution applied on `artifact-asiacrypt-2026`:
- `FCL_TAGBYTES=16` for the paper's `KYBER_K=2` / NIST-I instantiation;
- encryption computes SHA3-256 and appends only its first 16 bytes;
- FCL candidate verification compares only those 16 bytes;
- `KYBER_CIPHERTEXTBYTES` is therefore 688 bytes for `KYBER_K=2`;
- `test/test_kyber.c` contains compile-time assertions for the 16-byte tag and 688-byte ciphertext;
- inherited K=3/4 build targets retain their prior 32-byte tag behavior and are not paper claims.

This resolves the source-level mismatch. Table 6 cycle counts still need to be remeasured after the implementation change.

### B2. RESOLVED — END-512 now uses the paper's gamma=4

Paper:
- Table 2 gives END-512 `gamma = 4`.
- Algorithm 7 uses `gamma` in the computation of `w`.

Analysis notebook:
- `scripts/dependencies_impact_on_wrap_errors.ipynb` sets `gamma = 4` for END-512 and `gamma = 3` for END-1024.

Original pre-fix implementation:
- both C and AVX2 END-512 `keygen.c::compute_w()` set the effective tuning constant named `gamma2` to **6**.
- END-1024 sets it to **3**, matching Table 2 numerically.

The pre-fix source comments described the expression using `gamma^2` while the paper writes `gamma`; the numerical END-512 value also differed (paper/notebook 4 vs implementation 6).

Impact before the fix:
- the END-512 code and the paper/analysis were not the same parameter point.
- key-generation/trapdoor geometry, DFR behavior, and performance measurements could correspond to different settings.

Resolution applied on `artifact-asiacrypt-2026`:
- both portable-C and AVX2 END-512 `keygen.c::compute_w()` now use the numerical tuning coefficient `gamma=4`;
- source comments now use the same `gamma` notation as Algorithm 7;
- this matches Table 2 and the dependency notebook, which already used 4.

This resolves the source-level mismatch. END-512 correctness/performance, especially the Table 5 cycle counts, must be rerun after the change.

## Remaining blocking findings

### B3. END implementation secret sampling is not the fixed-weight distribution defined in the paper

Paper:
- Section 2.1 defines `T_{n,k}` as exactly `k` coefficients +1, exactly `k` coefficients -1, and the rest zero.
- Algorithm 7 samples `g <- T_{n,kg}`, `f <- T_{n,kf}`.
- Table 2 uses `(kg,kf)=(72,72)` for END-512 and `(96,96)` for END-1024.
- Section 6.3 states that Table 4 security estimates use those sparse ternary NTRU distributions.

Current implementation:
- END-512 `poly.c::ternary_sample_poly()` samples each coefficient independently from a 6-bit value with
  `Pr[-1]=9/64`, `Pr[+1]=9/64`. It therefore has expected counts 72/72, not fixed counts.
- END-1024 `keygen.c::ternary_sample_poly()` constructs an iid-like sampler with
  `Pr[-1]=Pr[+1]=3/32`, giving expected counts 96/96, not fixed counts.
- the concrete-security notebook instead calls `SparseTernary(n=..., p=k, m=k)`, matching the paper's fixed-weight description rather than the implementation's variable-weight sampling.

Impact:
- the distribution used by the concrete implementation differs from the distribution stated in Table 2 and Appendix D's NTRU assumption.
- the concrete-security and DFR calculations need to be checked against the actual implemented distribution if the iid sampler is intentional.

Resolution choices:
1. implement exact fixed-weight sampling matching `T_{n,k}`; or
2. specify the iid sparse ternary distribution in the paper and rerun all affected security/DFR calculations.

### B4. FCL shortlist/candidate family differs between DFR analysis and implementation

Analysis source:
- the committed END blocks in `scripts/dfr_frontier_checked.py` use:
  - END-512 gf stage: `m0=3`
  - END-512 GF stage: `m0=4`
  - END-1024 gf stage: `m0=3`
  - END-1024 GF stage: `m0=1`
- `artifact/reproduce_dfr.py` deliberately follows those committed analysis blocks.

Current END implementation:
- the decoder explicitly identifies only two FCL indices, `idx0` and `idx1`, in both stages.
- its hard-coded correction sequence is therefore not the same shortlist configuration used by the DFR script.
- for END-1024, Table 2 gives `B_GF,0 = 1`, yet the implementation's GF-stage code also contains two-position corrections involving both `idx0` and `idx1`.

Paper:
- Section 3 defines `E_cand` from a public FCL shortlist of size `m0`, including all error supports up to `B0`.
- Section 6.3 states that the implementation evaluates all candidates up to the fixed public bounds.

Impact:
- the code path being benchmarked is not obviously the same candidate family whose FCL-failure probability is reported in Table 3.
- END-1024 additionally evaluates weight-two-looking GF corrections despite the paper's `B_GF,0=1`.

Resolution:
- define the exact public `m_gf` and `m_GF` used by each END parameter set in one shared specification;
- make the DFR calculation and implementation use that same family;
- add a test that enumerates the implemented candidate set and compares it with the mathematical `E_cand` definition for the public parameters.

## High-priority reproducibility gaps

### R1. Table 5 baseline performance is not self-contained in this repository

Table 5 compares END against DAWN-alpha, DAWN-beta, BAT, NEV, and ML-KEM. The current repository contains END and FCL-ML-KEM, but not the exact baseline source revisions used for those performance numbers.

The paper identifies pq-crystals/Kyber for ML-KEM and says protected/constant-time baseline paths were used, but the artifact should pin the exact revisions, compiler flags, and invocation for every Table-5 baseline if the goal is an "Artifacts Reproduced" claim for the full table.

At minimum, distinguish:
- **END rows reproducible here**, versus
- **baseline rows imported from external implementations** with exact repository+commit+build instructions.

### R2. Timing test is built at a different optimization level from the benchmark binary

END `Makefile` targets build:
- `main` at `-O3`;
- `timecop` at `-O0`.

Paper Section 6.3 uses TIMECOP as evidence for timing independence of the submitted implementation paths. Compiler optimization can affect whether apparently branchless C becomes a branch or vice versa.

Recommended artifact check:
- run the timing test at the same optimization level used for the reported implementation (and optionally also at `-O0`);
- preserve the compiler/version and complete command line with the result.

### R3. The END smoke path builds but does not automatically assert KEM correctness

The implementation `main` prints the encapsulated and decapsulated shared keys, but the artifact smoke test currently only builds END to avoid the built-in 100000-iteration benchmark.

Recommended:
- add a separate one-shot test harness that performs keygen/encaps/decaps and exits nonzero unless the two shared keys match;
- keep benchmarking separate from correctness tests.

## Code-quality / portability findings

### Q1. Several decapsulation functions are declared `int` but return no value

Examples:
- END-1024 C `END_cca_decrypt`;
- END-512 AVX2 `END_cca_decrypt`;
- END-1024 AVX2 `END_cca_decrypt`.

The END-512 portable C variant declares it `void`.

The existing callers ignore the return value, so this is primarily API/undefined-return-value debt. Normalize all variants to a single signature (preferably `void` if no status is intended, or explicitly `return 0`).

### Q2. END-512 AVX2 public-key pointer type is inconsistent

`END-512-257-AVX2/kem.h` and `kem.c` declare `END_cca_encrypt(int16_t *h, ...)`, while:
- `main.c` allocates/passes `int8_t *h`;
- the encoded public key is a byte string;
- `decode_pk` accepts a byte pointer.

This happens to use the same pointer address at the ABI level, but it is an incompatible C type and should be corrected.

### Q3. AVX2 harness uses `aligned_alloc(32, len)` with non-multiple sizes

The AVX2 `xmalloc` wrappers call `aligned_alloc(32, len)`, while allocations include lengths such as:
- END-512 public key: 514 bytes;
- END-1024 public key: 1027 bytes;
- scratch buffers with `+18` bytes.

C11 requires the size passed to `aligned_alloc` to be an integral multiple of the alignment. Round sizes up to 32 bytes or use `posix_memalign`/an equivalent helper.

### Q4. RESOLVED — `FCL-ML-KEM/api.h` ciphertext-size constants

The public API constants have been updated to the actual FCL ciphertext formats. In particular, the paper's NIST-I target now exposes 688 bytes. The inherited K=3/4 targets expose their actual sizes under the retained 32-byte-tag behavior.

## Findings that align well with the paper

The review also found several strong matches:

- END encoded key/ciphertext sizes in the parameter headers match Table 2: 514/384 for END-512 and 1027/832 for END-1024.
- END-1024's tuning constant is 3 in both C and AVX2, matching Table 2 numerically.
- the dependency notebook uses exactly 100 trials for each END parameter set, as described in Section 2.4, and uses paper values `gamma=4` / `gamma=3`.
- the concrete-security notebook pins/prints the three estimator revisions; the standard lattice-estimator revision is the `3e48ef4...` commit stated in Section 6.3.
- FCL-ML-KEM uses `du=9`, `dv=3` for `KYBER_K=2`, matching Appendix A's compression setting.
- END includes both portable-C and AVX2 variants and dedicated TIMECOP harnesses, corresponding to Section 6.3's implementation discussion.

## Recommended next steps

1. Run the clean-environment smoke test to validate the B1/B2 source changes.
2. Remeasure END-512 Table-5 performance after the gamma change.
3. Remeasure FCL-ML-KEM-512 Table-6 performance after the 16-byte-tag change.
4. Decide separately how to handle B3 (fixed-weight vs iid secret sampler).
5. Decide separately how to handle B4 (implementation vs Table-3 FCL candidate family).
6. Add deterministic END correctness regression tests and run TIMECOP at production optimization.
7. Pin external Table-5 baseline repositories/commits if reproducing the whole comparison table is an artifact goal.
8. Freeze the artifact tag/archive only after the selected remaining scope is documented.

## Review boundary

This is a source/paper consistency review, not a cryptographic proof audit. Findings marked "blocking" mean the artifact currently cannot straightforwardly demonstrate that the checked source implements the exact paper parameterization/experiment; they do not by themselves establish an exploitable vulnerability.
