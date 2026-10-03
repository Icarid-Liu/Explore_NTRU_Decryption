# Generative-AI use disclosure

Generative-AI assistance was used to prepare the ASIACRYPT 2026 artifact packaging on the `artifact-asiacrypt-2026` branch.

OpenAI ChatGPT assisted with:

- restructuring and expanding the artifact README;
- preparing reproducibility/setup shell scripts;
- preparing the Docker/CI support files;
- preparing the artifact-only DFR orchestration wrapper;
- drafting the third-party provenance and submission-checklist documentation;
- reviewing ePrint 2026/1715 against the repository;
- under explicit author instruction, changing the NIST-I FCL-ML-KEM tag from 32 bytes to the paper's 16 bytes and changing the END-512 C/AVX2 trapdoor tuning coefficient from 6 to the paper's gamma=4.

The scientific-analysis files under `scripts/` were not edited. Pre-existing cryptographic implementation changes are limited to the two author-requested paper-alignment fixes listed above, plus the corresponding FCL size regression test/API constants. No change was made to the END secret sampler or FCL candidate-family logic identified separately in `CODE_REVIEW.md`.

The generated support material should be reviewed by the paper authors before archival, especially legal/provenance statements and the exact paper-section/table/figure mapping.
