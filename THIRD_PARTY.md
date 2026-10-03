# Third-party software and provenance

This file records third-party material that is visible from the committed source and the exact external dependencies used by the reproducibility workflow. It does not replace license notices embedded in individual files.

## ntrugen

Several files in the END implementations carry the header:

> Original source: ntrugen by Thomas Pornin

Upstream: https://github.com/pornin/ntrugen

The upstream `LICENSE` dedicates ntrugen to the public domain under the Unlicense terms. Files in this repository that retain the ntrugen attribution remain subject to those upstream terms for the upstream portions.

## CRYSTALS-Kyber-derived material

The `FCL-ML-KEM/` tree uses the CRYSTALS-Kyber reference implementation's file organization, API names, comments and build structure and contains adapted Kyber implementation material.

Upstream: https://github.com/pq-crystals/kyber

The upstream Kyber repository states that its code is available under CC0/Public Domain or Apache License 2.0, with Keccak/AES code under the public-domain terms identified in the corresponding source comments. This artifact retains those upstream terms for upstream-derived portions.

Apache License 2.0:
https://www.apache.org/licenses/LICENSE-2.0

CC0:
https://creativecommons.org/publicdomain/zero/1.0/

## Keccak / FIPS 202 code

The committed `fips202.c` files state that they are based on public-domain implementations from SUPERCOP and TweetFips202 and identify the original authors in their source headers. Those notices are retained.

## Valgrind client headers

Copies of `valgrind.h` and `memcheck.h` under the END implementation directories contain their own BSD-style license notices from the Valgrind project. Those file-specific notices govern those files.

Upstream project:
https://valgrind.org/

## External estimator repositories

The concrete-security notebook requires three repositories. The artifact setup script fetches them at the exact commits recorded by the notebook; they are **not vendored** into this repository:

- https://github.com/malb/lattice-estimator at `3e48ef421ec256afddb3e7d2249a77eab6e9ba12`
- https://github.com/identitymapping/enhanced_lattice-estimator at `876b66173f4354a96ddafc0ce3a79767ec43c6d4`
- https://github.com/yonghaason/PrimalMeetLWE at `61115115830c909e42758f2774606074bf98afb1`

Their own upstream licenses apply.

## OpenSSL

The FCL-ML-KEM NIST KAT targets link against the system OpenSSL crypto library (`-lcrypto`). OpenSSL is an external system dependency and is not vendored in this repository.

## Author provenance check before archival

Some implementation and assembly files do not contain an explicit provenance/license header. The repository authors should confirm before archival that those files are either original contributions covered by the root `LICENSE` or add the appropriate upstream attribution/license notice. In particular, review unannotated optimized/assembly files in the AVX2 END directories.

This explicit review is necessary because an artifact archive must not imply relicensing of third-party code merely because it is stored in the same repository.
