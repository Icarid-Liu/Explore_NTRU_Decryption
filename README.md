## END Source Code

This repository contains the source code for the END KEM.

The implementation is provided for two parameter sets:
- `END-512`
- `END-1024`

For each parameter set, we provide two implementations:
- A standard C implementation, located in `END-n/END-n-q-C/`
- An AVX2-optimized implementation, located in `END-n/END-n-q-AVX2/`

## Building and Running

To build and run the code, navigate to one of the implementation directories and use the provided `Makefile`.

For example, for the C implementation of `END-512`:

test correctness and time cost
```bash
cd END-512/END-512-257-C/
make main -j[jobs]
./main
```

measure memory usage
```bash
cd END-512/END-512-257-C/
make memory_usage -j[jobs]
valgrind --tool=massif --heap=yes --stacks=yes --massif-out-file=memory_profile.out ./memory_usage
ms_print memory_profile.out
```

detect timing attack
```bash
cd END-512/END-512-257-C/
make timecop -j[jobs]
valgrind --track-origins=yes ./timecop 2>&1 | grep "Conditional"
```
