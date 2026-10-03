# ASIACRYPT 2026 artifact submission checklist

This checklist separates repository work that can be automated from author-only metadata/legal checks that cannot be inferred safely from source code.

## Repository-level checks

- [x] Root `README.md` is the artifact entry point.
- [x] Root `LICENSE` is present for original contributions.
- [x] Third-party provenance is documented in `THIRD_PARTY.md`.
- [x] Exact estimator Git commits are pinned.
- [x] Notebook Sage/Python versions are recorded.
- [x] Build/run commands are documented.
- [x] Output interpretation is documented.
- [x] Artifact-only wrappers avoid edits to the main implementation/analysis source.
- [x] Docker and smoke/full reproduction helpers are present.
- [x] Generative-AI assistance is disclosed.
- [ ] Run `bash artifact/smoke_test.sh` on a clean Linux/x86-64 environment and retain the output.
- [ ] Run the full notebook reproduction with the pinned estimator commits.
- [ ] Run the performance experiments on the machine whose results will be reported.

## Paper/code blockers found in review

The detailed evidence is in `CODE_REVIEW.md`. Resolve these before freezing the artifact:

- [ ] FCL-ML-KEM: reconcile the paper's 16-byte authentication tag / 688-byte NIST-I ciphertext with the implementation's 32-byte tag.
- [ ] END-512: reconcile Table 2 / dependency-notebook `gamma=4` with C/AVX2 `compute_w()` using 6.
- [ ] END secret distribution: reconcile the paper's fixed-weight `T_{n,k}` with the implementation's variable-weight iid-style samplers.
- [ ] END DFR: make the FCL shortlist/candidate family used by the implementation identical to the family used for Table 3 calculations.
- [ ] After resolving the four items above, rerun Table 3 DFR and Table 4 concrete-security calculations and update outputs/documentation.
- [ ] Decide whether Table 5 baseline rows are in artifact scope; if yes, pin exact baseline repositories/commits/build commands.

## Author-only checks before archival

- [ ] Record the exact CPU, OS, compiler versions and relevant CPU-frequency/turbo settings for the paper's benchmark machine.
- [ ] Confirm that the root MIT grant is the license the rights holders intend for **their original contributions**.
- [ ] Review all unannotated files, especially optimized/assembly files, for third-party provenance and add missing notices if needed.
- [ ] Confirm whether any external dataset, script or private input used for a paper result is missing from this repository.

## HotCRP submission metadata required by the artifact call

Prepare the following separately from the repository:

- [ ] accepted-paper title and abstract;
- [ ] full author list and affiliations;
- [ ] contact email addresses (including at least one author);
- [ ] submitted/camera-ready paper PDF as requested by the artifact process;
- [ ] brief artifact description;
- [ ] public repository URL;
- [ ] immutable commit hash or release tag;
- [ ] archive upload if the artifact is below the conference upload limit, otherwise retrieval instructions;
- [ ] generative-AI acknowledgment (use `AI_USE_DISCLOSURE.md` as the basis).

## Freeze procedure

After all checks above:

```bash
git status --short
bash artifact/smoke_test.sh
bash artifact/package.sh
git rev-parse HEAD
```

Create an immutable tag for the reviewed commit, for example:

```bash
git tag -a asiacrypt2026-artifact-v1 -m "ASIACRYPT 2026 artifact v1"
git push origin asiacrypt2026-artifact-v1
```

Submit the tag and commit hash rather than a moving branch name.
