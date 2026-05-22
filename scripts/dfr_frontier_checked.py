#!/usr/bin/env python3
"""
Checked DFR/frontier calculations for Sections 4--5.

This file is intentionally independent of Sage/lattice-estimator.  It only checks
probability calculations used in the NTRU-with-Encoding and NTRU-with-Trapdoor
frontier discussion.  It follows the coefficient-independence model stated in
Section 2 and the FCL exclusion estimate stated in Section 3.
"""
from __future__ import annotations

from collections import defaultdict
from math import comb, log2
from mpmath import mp
from sage.all import *

mp.dps = 256
ROUNDING_FACTOR = 2 ** 512
RF = RealField(prec=512)
PROB_BOUND = RF("1e-128")
Pdf = dict[int, RF]


def round_to_rational(x):
    A = ZZ(round(x * ROUNDING_FACTOR))
    return QQ(A) / QQ(ROUNDING_FACTOR)


def average_variance(D):
    mu = 0.
    s = 0.

    for (v, p) in D.items():
        mu += v * p
        s += v * v * p

    s -= mu * mu
    return round_to_rational(mu), round_to_rational(s)


def binomial(x, y):
    try:
        binom = factorial(x) // factorial(y) // factorial(x - y)
    except ValueError:
        binom = 0
    return binom


def centered_binomial_pdf(k, x):
    return RF(binomial(2 * k, x + k)) / RF(2.**(2 * k))


def build_centered_binomial_law(k):
    D = {}
    for i in range(-k, k + 1):
        D[i] = centered_binomial_pdf(k, i)
    return D


def build_gaussian_law(sigma, t=None, mean=0, prob_bound=10 ** -64):
    D = {}
    if t is None:
        t = int(sigma * 5)
    real_bound = t
    D[mean] = RF(1)
    for i in range(1, t + 1):
        tmp = RF(exp(-i ** 2 / (2 * sigma ** 2)))
        if tmp < prob_bound:
            real_bound = i - 1
            break
        D[mean + i] = tmp
        D[mean - i] = tmp
    normalization = sum([D[i] for i in D])
    for i in D:
        D[i] = RF(D[i] / normalization)
    assert abs(sum([D[i] for i in range(mean - real_bound, mean + real_bound + 1)]) - 1.) <= 2 ** -256
    return D


def normalize(D: Pdf) -> Pdf:
    s = mp.fsum(D.values())
    if s == 0:
        raise ValueError("zero distribution")
    return {k: RF(v) / RF(s) for k, v in D.items() if v != 0}


def add_distribution(D1: Pdf, D2: Pdf, prob_bound: RF = PROB_BOUND) -> Pdf:
    out: dict[int, RF] = defaultdict(RF)
    for a, pa in D1.items():
        if pa < prob_bound:
            continue
        for b, pb in D2.items():
            if pb < prob_bound:
                continue
            p = pa * pb
            if p >= prob_bound:
                out[a + b] += p
    return normalize(dict(out))


def calculate_convolution_pdf(n, D1, D2, bound0=1e-64, bound1=1e-64, bound2=1e-64, bound3=1e-10):
    new_D = {}
    for a1 in D1.keys():
        if D1[a1] < bound0:
            continue
        for b1 in D2.keys():
            if D2[b1] < bound0:
                continue
            tmp = a1*b1
            if tmp not in new_D:
                new_D[tmp] = D1[a1] * D2[b1]
            else:
                new_D[tmp] += D1[a1] * D2[b1]
    PR = PolynomialRing(RR, 'X')
    X = PR.gens()[0]
    rho_x = 0
    for j in new_D.keys():
        if new_D[j] < bound1:
            continue
        rho_x += new_D[j] * X**j
    rho_x_up = rho_x.numerator() ** n
    D_after_convolution = {}
    down_degree = rho_x.denominator().degree() * n
    for j in range(rho_x_up.degree()+1):
        tmp = rho_x_up[j]
        if tmp < bound2:
            continue
        D_after_convolution[j-down_degree] = tmp
    normalization = sum([RF(di) for di in list(D_after_convolution.values())])
    if abs(normalization-1) > bound3:
        raise ValueError('Not enough precision')
    else:
        for k1 in range(min(D_after_convolution.keys()), max(D_after_convolution.keys()) + 1):
            D_after_convolution[k1] /= normalization
            return D_after_convolution


def get_accurate_gaussian_distribution(sigma, bounded_num=None):
    if bounded_num is None:
        bounded_num = int(round(10 * sigma))
    factor_times = 0
    while True:
        tmp_D = build_gaussian_law(sigma * (1 + factor_times * 0.0001), bounded_num)
        if abs(RR(sqrt(average_variance(tmp_D)[1])) - sigma) < 0.001:
            return tmp_D
        factor_times += 1


def get_compress_distribution(q, d):
    tmp_l = [int(int(round(int(int(round(vi * d / q)) % d) * (q/d))-vi) % q) for vi in range(q)]
    tmp_l = [vi if int(vi) < int(-vi % q) else -int(-vi%q) for vi in tmp_l]
    tmp_pdf = {}
    for vi in set(tmp_l):
        tmp_pdf[vi] = RF(tmp_l.count(vi) / q)
    tmp_mean = int(round(average_variance(tmp_pdf)[0]))
    new_pdf = {}
    for vi in tmp_pdf:
        new_pdf[vi-tmp_mean] = tmp_pdf[vi]
    return new_pdf


def scale_distribution(D: Pdf, c: int) -> Pdf:
    out: dict[int, RF] = defaultdict(RF)
    for a, p in D.items():
        out[c * a] += p
    return normalize(dict(out))


def ternary_coefficient_distribution(n: int, k: int) -> Pdf:
    """One coefficient of T_{n,k}: k entries +1, k entries -1."""
    return {-1: RF(k) / n, 0: RF(n - 2 * k) / n, 1: RF(k) / n}


def fixed_weight_inner_product_distribution(n: int, ka: int, kb: int) -> Pdf:
    """
    Distribution of <a,b> for independent a in T_{n,ka}, b in T_{n,kb}.
    This is the same fixed-Hamming-weight model used by the notebooks.
    """
    if 2 * ka > n or 2 * kb > n:
        raise ValueError("invalid T_{n,k} parameters")
    z = n - 2 * ka
    total = RF(comb(n, kb)) * RF(comb(n - kb, kb))
    counts: dict[int, int] = defaultdict(int)

    for x in range(0, min(ka, kb) + 1):
        for y in range(0, min(ka, kb - x) + 1):
            w = kb - x - y
            if w < 0 or w > z:
                continue
            count_bplus = comb(ka, x) * comb(ka, y) * comb(z, w)
            rem_ap = ka - x
            rem_am = ka - y
            rem_z = z - w
            for u in range(0, min(rem_ap, kb) + 1):
                for v in range(0, min(rem_am, kb - u) + 1):
                    t = kb - u - v
                    if t < 0 or t > rem_z:
                        continue
                    count_bminus = comb(rem_ap, u) * comb(rem_am, v) * comb(rem_z, t)
                    c = x + v - y - u
                    counts[c] += count_bplus * count_bminus

    return normalize({k: RF(cnt) / total for k, cnt in counts.items()})


def nearest_integer_to_ratio(y: int, q: int) -> int:
    """Nearest integer to y/q, with symmetric half-up behavior. q is odd here."""
    if y >= 0:
        return int(mp.floor(RF(y) / q + RF("0.5")))
    return -int(mp.floor(RF(-y) / q + RF("0.5")))


def error_class_probabilities(D: Pdf, q: int, B1: int) -> tuple[RF, RF, RF]:
    """Return p0, p_good, p_bad for e_i=round(y/q)."""
    p0 = RF("0")
    p_good = RF("0")
    p_bad = RF("0")
    for y, p in D.items():
        r = abs(nearest_integer_to_ratio(y, q))
        if r == 0:
            p0 += p
        elif r <= B1:
            p_good += p
        else:
            p_bad += p
    return p0, p_good, p_bad


def prob_error_outside_set(D: Pdf, q: int, n: int, B0: int, B1: int) -> RF:
    """Exact coefficient-iid probability Pr[e notin E_{B0,B1}]."""
    p0, p_good, _ = error_class_probabilities(D, q, B1)
    inside = RF("0")
    for m in range(0, B0 + 1):
        inside += RF(comb(n, m)) * (p_good ** m) * (p0 ** (n - m))
    return max(RF("0"), 1 - inside)


def binomial_tail_at_least(N: int, p: RF, r: int) -> RF:
    if r <= 0:
        return RF("1")
    if r > N:
        return RF("0")
    if p == 0:
        return RF("0")
    if p == 1:
        return RF("1")
    # Direct summation is fine for N <= 512 in these checks.
    return mp.fsum(RF(comb(N, k)) * (p ** k) * ((1 - p) ** (N - k)) for k in range(r, N + 1))


def conditioned_wrapped_distribution(D: Pdf, q: int, B1: int) -> Pdf:
    out: dict[int, RF] = {}
    for y, p in D.items():
        r = abs(nearest_integer_to_ratio(y, q))
        if 1 <= r <= B1:
            out[y] = p
    return normalize(out)


def conditioned_unwrapped_abs_distribution(D: Pdf, q: int) -> Pdf:
    out: dict[int, RF] = defaultdict(RF)
    for y, p in D.items():
        if nearest_integer_to_ratio(y, q) == 0:
            out[abs(y)] += p
    return normalize(dict(out))


def fcl_exclusion_probability_for_weight(
    D: Pdf,
    q: int,
    N: int,
    m1: int,
    m0: int,
    B1: int,
    tie_is_failure: bool = True,
) -> RF:
    """
    Union-bound probability that at least one of m1 wrapped coordinates is not
    contained in the m0-coordinate FCL shortlist, conditioned on exactly m1 good
    wrapped coordinates inside the FCL search domain of size N.
    """
    if m0 < m1:
        return RF("1")
    Dw = conditioned_wrapped_distribution(D, q, B1)
    Du_abs = conditioned_unwrapped_abs_distribution(D, q)
    # tail_abs[t] = Pr[|Z| >= t | Z unwrapped].  If ties are not counted, use > t.
    support = sorted(Du_abs)

    def p_better(residue_abs: int) -> RF:
        if tie_is_failure:
            return mp.fsum(p for z_abs, p in Du_abs.items() if z_abs >= residue_abs)
        return mp.fsum(p for z_abs, p in Du_abs.items() if z_abs > residue_abs)

    threshold = m0 - m1 + 1
    per_coordinate = RF("0")
    for y, py in Dw.items():
        r = nearest_integer_to_ratio(y, q)
        residue_abs = abs(y - q * r)
        per_coordinate += py * binomial_tail_at_least(N - m1, p_better(residue_abs), threshold)
    return min(RF("1"), m1 * per_coordinate)


def fcl_failure_probability(
    D: Pdf,
    q: int,
    n: int,
    B0: int,
    B1: int,
    m0_by_weight,
    N_by_weight,
) -> RF:
    """Unconditional FCL-miss probability, within e in E_{B0,B1}."""
    p0, p_good, _ = error_class_probabilities(D, q, B1)
    out = RF("0")
    for m1 in range(1, B0 + 1):
        p_exact = RF(comb(n, m1)) * (p_good ** m1) * (p0 ** (n - m1))
        out += p_exact * fcl_exclusion_probability_for_weight(
            D, q, N_by_weight(m1), m1, m0_by_weight(m1), B1
        )
    return out


def raw_candidate_count(m0: int, B0: int, B1: int) -> int:
    return 1 + sum(comb(m0, i) * (2 * B1) ** i for i in range(1, B0 + 1))


def log2p(x: RF) -> float:
    return float(mp.log(x, 2)) if x > 0 else float("-inf")


def log2_bigint(x: int) -> float:
    """Approximate log2(x) for very large integers."""
    if x <= 0:
        raise ValueError("x must be positive")
    b = x.bit_length()
    if b <= 1022:
        return log2(x)
    top = x >> (b - 53)
    return (b - 53) + log2(top)


def rejection_probability(N: int, m: int, Kver: int, lam: int) -> RF:
    """
    N: half dimension, so n = 2N.
    m: total number of nonzero coefficients in mu.
       If mu <- T_{n,k_mu} as defined in the paper, then m = 2*k_mu.
    Kver: verification weight bound.
    lam: embedded message length lambda.
    """
    total = comb(2 * N, m)
    accepted = 0

    for j in range(0, min(m, N) + 1):
        if (m - j) % 2 != 0:
            continue

        t = (m - j) // 2  # number of colliding pairs

        if t < 0 or t > N - j:
            continue

        count_j = comb(N, j) * (1 << j) * comb(N - j, t)

        admissible = (
            j <= Kver and
            comb(N, j) >= (1 << lam)
        )

        if admissible:
            accepted += count_j

    return RF((total - accepted) / total)


def lntru_encoding_distribution(n, kgf, kse) -> Pdf:
    Dgs = fixed_weight_inner_product_distribution(n, kgf, kse)
    Dfe = fixed_weight_inner_product_distribution(n, kgf, kse)
    Dgs_fe = add_distribution(Dgs, Dfe)
    Dtgs_fe = add_distribution(Dgs_fe, Dgs_fe)  # t=(1,1)
    De = ternary_coefficient_distribution(n, kse)
    Dwe = add_distribution(De, De)              # w=(1,1)
    return add_distribution(Dtgs_fe, Dwe)


def rlwe_encoding_distribution(n, kse, ke1r, q, ke2) -> Pdf:
    De1s = fixed_weight_inner_product_distribution(n, ke1r, kse)
    Der = fixed_weight_inner_product_distribution(n, ke1r, kse)
    De2 = get_compress_distribution(q, ke2)
    De1s_er_e2 = add_distribution(add_distribution(De1s, Der), De2)
    Dte1s_er_e2 = add_distribution(De1s_er_e2, De1s_er_e2)
    return Dte1s_er_e2


def trapdoor_gf_distribution(n, k) -> Pdf:
    Dgs = fixed_weight_inner_product_distribution(n, k, k)
    Dfe = fixed_weight_inner_product_distribution(n, k, k)
    return add_distribution(Dgs, Dfe)


def trapdoor_GF_distribution(n, GFstar_std, k) -> Pdf:
    # Current notebook chooses k_GF=9 for G*,F* from its Gaussian-heuristic std screen.
    DGstar_s = calculate_convolution_pdf(n, get_accurate_gaussian_distribution(GFstar_std), ternary_coefficient_distribution(n, k))
    DFstar_e = calculate_convolution_pdf(n, get_accurate_gaussian_distribution(GFstar_std), ternary_coefficient_distribution(n, k))

    return add_distribution(DGstar_s, DFstar_e)


def main() -> None:
    # # NTRU with Encoding frontier
    # n = 512
    # q = 81
    # kgf = 32
    # kse = 32
    # Denc = lntru_encoding_distribution(n, kgf, kse)
    # B0 = 4
    # B1 = 1
    # m0 = 6
    # P_out = prob_error_outside_set(Denc, q, n, B0=B0, B1=B1)
    # P_fcl_m0 = fcl_failure_probability(Denc, q, n, B0, B1, lambda m: m0, lambda m: min(4 ** m, n))
    # print("NTRU-with-Encoding frontier check: n=%d, q=%d, T_{%d,%d}, t=w=1+x" % (n, q, n, kgf))
    # print(f"  Pr[e notin E_{{%d,%d}}] = 2^({log2p(P_out):.4f})" % (B0, B1))
    # print(f"  FCL miss, shortlist m0=%d, raw candidates {raw_candidate_count(m0,B0,B1)} = 2^({log2p(P_fcl_m0):.4f})" % m0)
    # print(f"  DFR = 2^({log2p(P_out + P_fcl_m0):.4f})")

    # # # RLWE with Encoding frontier
    # n = 512
    # q = 83
    # kse = 30
    # ke1r = 30
    # ke2 = 11
    # Denc = rlwe_encoding_distribution(n, kse, ke1r, q, ke2)
    # B0 = 4
    # B1 = 1
    # m0 = 6
    # P_out = prob_error_outside_set(Denc, q, n, B0=B0, B1=B1)
    # P_fcl_m0 = fcl_failure_probability(Denc, q, n, B0, B1, lambda m: m0, lambda m: min(4 ** m, n))
    # print("RLWE-with-Encoding frontier check: n=%d, q=%d, T_{%d,%d}, t=w=1+x" % (n, q, n, kgf))
    # print(f"  Pr[e notin E_{{%d,%d}}] = 2^({log2p(P_out):.4f})" % (B0, B1))
    # print(f"  FCL miss, shortlist m0=%d, raw candidates {raw_candidate_count(m0,B0,B1)} = 2^({log2p(P_fcl_m0):.4f})" % m0)
    # print(f"  DFR = 2^({log2p(P_out + P_fcl_m0):.4f})")

    # # Kyber with FCL
    # n = 512
    # q = 3329
    # eta1 = 3
    # eta2 = 2
    # du = 9
    # dv = 3
    # B0 = 2
    # B1 = 1
    # m0 = 2
    # De = build_centered_binomial_law(eta1)
    # Ds = build_centered_binomial_law(eta1)
    # Dr = build_centered_binomial_law(eta1)
    # De1 = build_centered_binomial_law(eta2)
    # De2 = build_centered_binomial_law(eta2)
    # Ddu = get_compress_distribution(q, 2**du)
    # Ddv = get_compress_distribution(q, 2**dv)
    # De1du = add_distribution(De1, Ddu)
    # De2dv = add_distribution(De2, Ddv)
    # De1du_s = calculate_convolution_pdf(n, De1du, Ds)
    # De_r = calculate_convolution_pdf(n, De, Dr)
    # Denc = add_distribution(De1du_s, add_distribution(De_r, De2dv))
    # Denc2 = {}
    # for vi in Denc:
    #     Denc2[2*vi] = Denc[vi]
    # P_out = prob_error_outside_set(Denc2, q, n//2, B0=B0, B1=B1)
    # P_fcl_m0 = fcl_failure_probability(Denc2, q, n, B0, B1, lambda m: m0, lambda m: n)
    # print("Kyber-512 with FCL: du = %d, dv = %d" % (du, dv))
    # print(f"  Pr[e notin E_{{%d,%d}}] = 2^({log2p(P_out):.4f})" % (B0, B1))
    # print(f"  FCL miss, aux N=4, shortlist m0=%d, raw candidates {raw_candidate_count(m0,B0,B1)} = 2^({log2p(P_fcl_m0):.4f})" % m0)
    # print(f"  DFR = 2^({log2p(P_out + P_fcl_m0):.4f})")

    # # NwE-512
    # n = 512
    # q = 193
    # kgf = 40
    # kse = 40
    # Denc = lntru_encoding_distribution(n, kgf, kse)
    # B0 = 1
    # B1 = 1
    # m0 = 1
    # P_out = prob_error_outside_set(Denc, q, n, B0=B0, B1=B1)
    # P_fcl_m0 = fcl_failure_probability(Denc, q, n, B0, B1, lambda m: m0, lambda m: min(3 ** m, n))
    # print("NwE check: n=%d, q=%d, T_{%d,%d}, t=w=1+x" % (n, q, n, kgf))
    # print(f"  Pr[e notin E_{{%d,%d}}] = 2^({log2p(P_out):.4f})" % (B0, B1))
    # print(f"  FCL miss, shortlist m0=%d, raw candidates {raw_candidate_count(m0,B0,B1)} = 2^({log2p(P_fcl_m0):.4f})" % m0)
    # print(f"  DFR = 2^({log2p(P_out + P_fcl_m0):.4f})")
    # print(f"  rho = 2^({log2p(rejection_probability(N=n//2, m=kse*2, Kver=171//2, lam=128)):.4f})")


    # # NwE-1024
    # n = 1024
    # q = 193
    # kgf = 80
    # kse = 64
    # Denc = lntru_encoding_distribution(n, kgf, kse)
    # B0 = 1
    # B1 = 1
    # m0 = 1
    # P_out = prob_error_outside_set(Denc, q, n, B0=B0, B1=B1)
    # P_fcl_m0 = fcl_failure_probability(Denc, q, n, B0, B1, lambda m: m0, lambda m: min(3 ** m, n))
    # print("NwE check: n=%d, q=%d, T_{%d,%d}, t=w=1+x" % (n, q, n, kgf))
    # print(f"  Pr[e notin E_{{%d,%d}}] = 2^({log2p(P_out):.4f})" % (B0, B1))
    # print(f"  FCL miss, shortlist m0=%d, raw candidates {raw_candidate_count(m0,B0,B1)} = 2^({log2p(P_fcl_m0):.4f})" % m0)
    # print(f"  DFR = 2^({log2p(P_out + P_fcl_m0):.4f})")
    # print(f"  rho = 2^({log2p(rejection_probability(N=n//2, m=kse*2, Kver=261//2, lam=256)):.4f})")

    # # NTRU with Trapdoor frontier
    # n = 512
    # q = 59
    # kgf = 26
    # kse = 26
    # B0_gf = 3
    # B1_gf = 1
    # m0_gf = 7
    # Dgf = trapdoor_gf_distribution(n, kgf)
    # P_gf_out = prob_error_outside_set(Dgf, q, n, B0=B0_gf, B1=B1_gf)
    # P_gf_fcl_m = fcl_failure_probability(Dgf, q, n, B0_gf, B1_gf, lambda m: m0_gf, lambda m: n)
    # print("\nNTRU-with-Trapdoor gf-stage check: n=%d, q=%d, T_{%d,%d}" % (n, q, n, kgf))
    # print(f"  Pr[E_gf notin E_{{%d, %d}}] = 2^({log2p(P_gf_out):.4f})" % (B0_gf, B1_gf))
    # print(f"  FCL miss with I_n and fixed shortlist m0=%d, raw candidates {raw_candidate_count(m0_gf,B0_gf,B1_gf)} = 2^({log2p(P_gf_fcl_m):.4f})" % m0_gf)

    # GF_star_std = 0.256
    # B0_GF = 2
    # B1_GF = 1
    # m0_GF = 7
    # DGF = trapdoor_GF_distribution(n, GF_star_std, kse)
    # P_GF_out = prob_error_outside_set(DGF, q, n, B0=B0_GF, B1=B1_GF)
    # P_GF_fcl_m = fcl_failure_probability(DGF, q, n, B0_GF, B1_GF, lambda m: m0_GF, lambda m: n)
    # print("\nNTRU-with-Trapdoor GF-stage check under notebook distribution:")
    # print(f"  Pr[E_GF notin E_{{%d,%d}}] = 2^({log2p(P_GF_out):.4f})" % (B0_GF, B1_GF))
    # print(f"  FCL miss with I_n, raw candidates {raw_candidate_count(m0_GF,B0_GF,B1_GF)} = 2^({log2p(P_GF_fcl_m):.4f})")
    # print(f"  DFR = 2^({log2p(P_gf_out + P_gf_fcl_m + P_GF_out + P_GF_fcl_m):.4f})")
    
    # # END-512
    # n = 512
    # q = 257 * 2
    # kgf = 72
    # B0_gf = 2
    # B1_gf = 1
    # m0_gf = 3
    # Dgs = calculate_convolution_pdf(n, ternary_coefficient_distribution(n, kgf), {-1: 1/2, 1: 1/2})
    # Dfe = calculate_convolution_pdf(n, ternary_coefficient_distribution(n, kgf), {-3: 1/4, -1: 1/4, 1: 1/4, 3: 1/4})
    # Dgf = add_distribution(Dgs, Dfe)
    # P_gf_out = prob_error_outside_set(Dgf, q, n, B0=B0_gf, B1=B1_gf)
    # P_gf_fcl_m = fcl_failure_probability(Dgf, q, n, B0_gf, B1_gf, lambda m: m0_gf, lambda m: n)
    # print("\nNTRU-with-Trapdoor gf-stage check: n=%d, q=%d, T_{%d,%d}" % (n, q, n, kgf))
    # print(f"  Pr[E_gf notin E_{{%d, %d}}] = 2^({log2p(P_gf_out):.4f})" % (B0_gf, B1_gf))
    # print(f"  FCL miss with I_n and fixed shortlist m0=%d, raw candidates {raw_candidate_count(m0_gf,B0_gf,B1_gf)} = 2^({log2p(P_gf_fcl_m):.4f})" % m0_gf)
    
    # G_star_std = 1.0641
    # F_star_std = 0.3644
    # B0_GF = 2
    # B1_GF = 1
    # m0_GF = 4
    # DGs = calculate_convolution_pdf(n, get_accurate_gaussian_distribution(G_star_std), {-1: 1/2, 1: 1/2})
    # DFe = calculate_convolution_pdf(n, get_accurate_gaussian_distribution(F_star_std), {-3: 1/4, -1: 1/4, 1: 1/4, 3: 1/4})
    # print(RR(sqrt(average_variance(DGs)[1])), RR(sqrt(average_variance(DFe)[1])))
    # DGF = add_distribution(DGs, DFe)
    # P_GF_out = prob_error_outside_set(DGF, q, n, B0=B0_GF, B1=B1_GF)
    # P_GF_fcl_m = fcl_failure_probability(DGF, q, n, B0_GF, B1_GF, lambda m: m0_GF, lambda m: n)
    # print("\nNTRU-with-Trapdoor GF-stage check under notebook distribution:")
    # print(f"  Pr[E_GF notin E_{{%d,%d}}] = 2^({log2p(P_GF_out):.4f})" % (B0_GF, B1_GF))
    # print(f"  FCL miss with I_n and fixed shortlist m0=2, raw candidates {raw_candidate_count(m0_GF,B0_GF,B1_GF)} = 2^({log2p(P_GF_fcl_m):.4f})")
    # print(f"  DFR = 2^({log2p(P_gf_out + P_gf_fcl_m + P_GF_out + P_GF_fcl_m):.4f})")
    
    # # END-1024
    # n = 1024
    # q = 257 * 2
    # kgf = 96
    # B0_gf = 2
    # B1_gf = 1
    # m0_gf = 3
    # Dgs = calculate_convolution_pdf(n, ternary_coefficient_distribution(n, kgf), {-1: 1/2, 1: 1/2})
    # Dfe = calculate_convolution_pdf(n, ternary_coefficient_distribution(n, kgf), {-2:77/257, 0:90/257, 2:90/257})
    # Dgf = add_distribution(Dgs, Dfe)
    # P_gf_out = prob_error_outside_set(Dgf, q, n, B0=B0_gf, B1=B1_gf)
    # P_gf_fcl_m = fcl_failure_probability(Dgf, q, n, B0_gf, B1_gf, lambda m: m0_gf, lambda m: n)
    # print("\nNTRU-with-Trapdoor gf-stage check: n=%d, q=%d, T_{%d,%d}" % (n, q, n, kgf))
    # print(f"  Pr[E_gf notin E_{{%d, %d}}] = 2^({log2p(P_gf_out):.4f})" % (B0_gf, B1_gf))
    # print(f"  FCL miss with I_n and fixed shortlist m0=%d, raw candidates {raw_candidate_count(m0_gf,B0_gf,B1_gf)} = 2^({log2p(P_gf_fcl_m):.4f})" % m0_gf)
    
    # G_star_std = 0.574
    # F_star_std = 0.275
    # B0_GF = 1
    # B1_GF = 1
    # m0_GF = 1
    # DGs = calculate_convolution_pdf(n, get_accurate_gaussian_distribution(G_star_std), {-1: 1/2, 1: 1/2})
    # DFe = calculate_convolution_pdf(n, get_accurate_gaussian_distribution(F_star_std), {-2:77/257, 0:90/257, 2:90/257})
    # print(RR(sqrt(average_variance(DGs)[1])), RR(sqrt(average_variance(DFe)[1])))
    # DGF = add_distribution(DGs, DFe)
    # P_GF_out = prob_error_outside_set(DGF, q, n, B0=B0_GF, B1=B1_GF)
    # P_GF_fcl_m = fcl_failure_probability(DGF, q, n, B0_GF, B1_GF, lambda m: m0_GF, lambda m: n)
    # print("\nNTRU-with-Trapdoor GF-stage check under notebook distribution:")
    # print(f"  Pr[E_GF notin E_{{%d,%d}}] = 2^({log2p(P_GF_out):.4f})" % (B0_GF, B1_GF))
    # print(f"  FCL miss with I_n and fixed shortlist m0=2, raw candidates {raw_candidate_count(m0_GF,B0_GF,B1_GF)} = 2^({log2p(P_GF_fcl_m):.4f})")
    # print(f"  DFR = 2^({log2p(P_gf_out + P_gf_fcl_m + P_GF_out + P_GF_fcl_m):.4f})")

if __name__ == "__main__":
    main()
