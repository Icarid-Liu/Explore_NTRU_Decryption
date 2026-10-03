# Explore NTRU Decryption

Artifact implementation for:

**How Compact Can NTRU Encryption Be? Heuristic Frontiers and Practical Schemes**  
IACR ePrint 2026/1715: https://eprint.iacr.org/2026/1715

This repository contains the two concrete implementations described in the paper:

```text
.
├── END KEM/
│   ├── END-512/
│   │   ├── END-512-257-C/
│   │   └── END-512-257-AVX2/
│   └── END-1024/
│       ├── END-1024-257-C/
│       └── END-1024-257-AVX2/
└── FCL-ML-KEM/
```

## Requirements

Linux/x86-64 is recommended.

Required tools:

- GCC or Clang
- GNU Make
- OpenSSL development headers for the FCL-ML-KEM NIST KAT target
- Valgrind only for the optional TIMECOP checks

The AVX2 END implementations require an AVX2-capable x86-64 CPU.

## END KEM

The repository contains portable C and AVX2 implementations for END-512 and END-1024.

Paper parameter sizes:

| Scheme | Public key | Ciphertext |
| --- | ---: | ---: |
| END-512 | 514 bytes | 384 bytes |
| END-1024 | 1027 bytes | 832 bytes |

### Portable C

END-512:

```bash
cd "END KEM/END-512/END-512-257-C"
make main
./main
```

END-1024:

```bash
cd "END KEM/END-1024/END-1024-257-C"
make main
./main
```

### AVX2

END-512:

```bash
cd "END KEM/END-512/END-512-257-AVX2"
make main
./main
```

END-1024:

```bash
cd "END KEM/END-1024/END-1024-257-AVX2"
make main
./main
```

Each `main` executable performs a KEM correctness run and reports KeyGen, Encapsulation, and Decapsulation cycle counts.

### Optional TIMECOP check

From an END implementation directory:

```bash
make timecop
valgrind --track-origins=yes ./timecop 2>&1 | grep "Conditional"
```

## FCL-ML-KEM

The paper implementation is the NIST-I / `KYBER_K=2` instantiation.

Its parameters include:

- `du = 9`
- `dv = 3`
- 128-bit verification tag
- 800-byte public key
- 688-byte ciphertext

Build and run the correctness test:

```bash
cd FCL-ML-KEM
make test/test_kyber512
./test/test_kyber512
```

Build and run the cycle benchmark:

```bash
make test/test_speed512
./test/test_speed512
```

Build the NIST KAT generator:

```bash
make nistkat/PQCgenKAT_kem512
./nistkat/PQCgenKAT_kem512
```

The inherited `KYBER_K=3` and `KYBER_K=4` build targets are not part of the paper's FCL-ML-KEM parameter set.

## Licensing

See `LICENSE` for the repository license and `THIRD_PARTY.md` for third-party source notices.
