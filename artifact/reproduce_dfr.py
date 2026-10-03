#!/usr/bin/env python3
"""Run selected DFR/frontier cases without editing scripts/dfr_frontier_checked.py.

The numerical primitives and probability calculations live in the original
scripts/dfr_frontier_checked.py.  This file only selects parameter blocks that
are already present in that file's commented main() examples and prints them in
a reviewer-friendly form.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import dfr_frontier_checked as dfr  # noqa: E402


def p2(x) -> str:
    return f"2^({dfr.log2p(x):.6f})"


def report(label: str, outside, miss, total=None) -> None:
    if total is None:
        total = outside + miss
    print(label)
    print(f"  outside allowed error set : {p2(outside)}")
    print(f"  FCL exclusion / miss      : {p2(miss)}")
    print(f"  combined                  : {p2(total)}")


def ntru_encoding_frontier() -> None:
    n = 512
    q = 81
    kgf = 32
    kse = 32
    B0 = 4
    B1 = 1
    m0 = 6

    Denc = dfr.lntru_encoding_distribution(n, kgf, kse)
    p_out = dfr.prob_error_outside_set(Denc, q, n, B0=B0, B1=B1)
    p_miss = dfr.fcl_failure_probability(
        Denc, q, n, B0, B1,
        lambda _m: m0,
        lambda m: min(4 ** m, n),
    )
    report("NTRU-with-Encoding frontier (n=512, q=81)", p_out, p_miss)


def ntru_trapdoor_frontier() -> None:
    n = 512
    q = 59
    k = 26

    Dgf = dfr.trapdoor_gf_distribution(n, k)
    p_gf_out = dfr.prob_error_outside_set(Dgf, q, n, B0=3, B1=1)
    p_gf_miss = dfr.fcl_failure_probability(
        Dgf, q, n, 3, 1,
        lambda _m: 7,
        lambda _m: n,
    )

    DGF = dfr.trapdoor_GF_distribution(n, 0.256, k)
    p_GF_out = dfr.prob_error_outside_set(DGF, q, n, B0=2, B1=1)
    p_GF_miss = dfr.fcl_failure_probability(
        DGF, q, n, 2, 1,
        lambda _m: 7,
        lambda _m: n,
    )

    total = p_gf_out + p_gf_miss + p_GF_out + p_GF_miss
    print("NTRU-with-Trapdoor frontier (n=512, q=59)")
    print(f"  gf outside                : {p2(p_gf_out)}")
    print(f"  gf FCL miss               : {p2(p_gf_miss)}")
    print(f"  GF outside                : {p2(p_GF_out)}")
    print(f"  GF FCL miss               : {p2(p_GF_miss)}")
    print(f"  combined DFR              : {p2(total)}")


def end_512() -> None:
    n = 512
    q = 257 * 2
    kgf = 72

    pm = {-1: dfr.RF(1) / 2, 1: dfr.RF(1) / 2}
    mu = {
        -3: dfr.RF(1) / 4,
        -1: dfr.RF(1) / 4,
        1: dfr.RF(1) / 4,
        3: dfr.RF(1) / 4,
    }

    Dgs = dfr.calculate_convolution_pdf(
        n, dfr.ternary_coefficient_distribution(n, kgf), pm
    )
    Dfe = dfr.calculate_convolution_pdf(
        n, dfr.ternary_coefficient_distribution(n, kgf), mu
    )
    Dgf = dfr.add_distribution(Dgs, Dfe)

    p_gf_out = dfr.prob_error_outside_set(Dgf, q, n, B0=2, B1=1)
    p_gf_miss = dfr.fcl_failure_probability(
        Dgf, q, n, 2, 1,
        lambda _m: 3,
        lambda _m: n,
    )

    DGs = dfr.calculate_convolution_pdf(
        n, dfr.get_accurate_gaussian_distribution(dfr.RF("1.0641")), pm
    )
    DFe = dfr.calculate_convolution_pdf(
        n, dfr.get_accurate_gaussian_distribution(dfr.RF("0.3644")), mu
    )
    DGF = dfr.add_distribution(DGs, DFe)

    p_GF_out = dfr.prob_error_outside_set(DGF, q, n, B0=2, B1=1)
    p_GF_miss = dfr.fcl_failure_probability(
        DGF, q, n, 2, 1,
        lambda _m: 4,
        lambda _m: n,
    )

    total = p_gf_out + p_gf_miss + p_GF_out + p_GF_miss
    print("END-512 DFR block")
    print(f"  gf outside                : {p2(p_gf_out)}")
    print(f"  gf FCL miss               : {p2(p_gf_miss)}")
    print(f"  GF outside                : {p2(p_GF_out)}")
    print(f"  GF FCL miss               : {p2(p_GF_miss)}")
    print(f"  combined DFR              : {p2(total)}")


def end_1024() -> None:
    n = 1024
    q = 257 * 2
    kgf = 96

    pm = {-1: dfr.RF(1) / 2, 1: dfr.RF(1) / 2}
    mu = {
        -2: dfr.RF(77) / 257,
        0: dfr.RF(90) / 257,
        2: dfr.RF(90) / 257,
    }

    Dgs = dfr.calculate_convolution_pdf(
        n, dfr.ternary_coefficient_distribution(n, kgf), pm
    )
    Dfe = dfr.calculate_convolution_pdf(
        n, dfr.ternary_coefficient_distribution(n, kgf), mu
    )
    Dgf = dfr.add_distribution(Dgs, Dfe)

    p_gf_out = dfr.prob_error_outside_set(Dgf, q, n, B0=2, B1=1)
    p_gf_miss = dfr.fcl_failure_probability(
        Dgf, q, n, 2, 1,
        lambda _m: 3,
        lambda _m: n,
    )

    DGs = dfr.calculate_convolution_pdf(
        n, dfr.get_accurate_gaussian_distribution(dfr.RF("0.574")), pm
    )
    DFe = dfr.calculate_convolution_pdf(
        n, dfr.get_accurate_gaussian_distribution(dfr.RF("0.275")), mu
    )
    DGF = dfr.add_distribution(DGs, DFe)

    p_GF_out = dfr.prob_error_outside_set(DGF, q, n, B0=1, B1=1)
    p_GF_miss = dfr.fcl_failure_probability(
        DGF, q, n, 1, 1,
        lambda _m: 1,
        lambda _m: n,
    )

    total = p_gf_out + p_gf_miss + p_GF_out + p_GF_miss
    print("END-1024 DFR block")
    print(f"  gf outside                : {p2(p_gf_out)}")
    print(f"  gf FCL miss               : {p2(p_gf_miss)}")
    print(f"  GF outside                : {p2(p_GF_out)}")
    print(f"  GF FCL miss               : {p2(p_GF_miss)}")
    print(f"  combined DFR              : {p2(total)}")


CASES = {
    "ntru-encoding-frontier": ntru_encoding_frontier,
    "ntru-trapdoor-frontier": ntru_trapdoor_frontier,
    "end-512": end_512,
    "end-1024": end_1024,
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=["all", *CASES], default="all")
    args = parser.parse_args()

    names = list(CASES) if args.case == "all" else [args.case]
    for i, name in enumerate(names):
        CASES[name]()
        if i + 1 != len(names):
            print()


if __name__ == "__main__":
    main()
