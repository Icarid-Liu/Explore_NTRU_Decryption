#!/usr/bin/env python3
"""
Estimate the Pr[Dec. fails] verification/collision term for the two frameworks.

This script intentionally separates this term from the two other DFR terms:
    Pr[rounding error not in E_err] and Pr[FCL misses].

Implemented terms.

1. NTRU with Encoding.
   The generic estimate is the packing heuristic from Section 4:
       omega = max(|M'| L / |R_M| - 1, 0).
   For binary fixed-weight embedding checks, the script can also estimate the
   probability that a wrong candidate shift keeps a projected message inside
   the accepted fixed-weight set.  If every nonzero candidate shift has weight
   > 2K and valid projected messages have weight <= K, this term is zero.

2. NTRU with Trapdoor.
   The script estimates epsilon_gf + epsilon_GF from the distributional
   verification model:
       epsilon_gf: norm-score collision, using a Gaussian approximation for
                   squared norm sums.
       epsilon_GF: range/fixed-distribution false accept for the final pair.

The script uses only Python + mpmath.  It does not use Sage or lattice-estimator.
It follows the notation of the paper as closely as possible.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import reduce
from math import comb
from typing import Dict, Iterable, Mapping, Sequence, Tuple

from mpmath import mp

mp.dps = 256
Pdf = Dict[int, mp.mpf]


# ---------------------------------------------------------------------------
# Generic probability utilities
# ---------------------------------------------------------------------------

def _mp(x) -> mp.mpf:
    return mp.mpf(x)


def normalize(D: Mapping[int, object]) -> Pdf:
    s = mp.fsum(_mp(v) for v in D.values())
    if s == 0:
        raise ValueError("zero distribution")
    return {int(k): _mp(v) / s for k, v in D.items() if _mp(v) != 0}


def add_distribution(D1: Mapping[int, object], D2: Mapping[int, object]) -> Pdf:
    out: Dict[int, mp.mpf] = {}
    for a, pa0 in D1.items():
        pa = _mp(pa0)
        for b, pb0 in D2.items():
            pb = _mp(pb0)
            out[int(a) + int(b)] = out.get(int(a) + int(b), mp.mpf(0)) + pa * pb
    return normalize(out)


def sum_distributions(distributions: Sequence[Mapping[int, object]]) -> Pdf:
    if not distributions:
        return {0: mp.mpf(1)}
    return reduce(add_distribution, distributions)  # type: ignore[arg-type]


def scale_distribution(D: Mapping[int, object], c: int) -> Pdf:
    out: Dict[int, mp.mpf] = {}
    for x, p0 in D.items():
        y = int(c) * int(x)
        out[y] = out.get(y, mp.mpf(0)) + _mp(p0)
    return normalize(out)


def ternary_coefficient_distribution(n: int, k_plus: int, k_minus: int | None = None) -> Pdf:
    """One coefficient of a fixed-weight ternary T_{n,k_plus,k_minus}."""
    if k_minus is None:
        k_minus = k_plus
    if k_plus < 0 or k_minus < 0 or k_plus + k_minus > n:
        raise ValueError("invalid ternary parameters")
    return normalize({1: mp.mpf(k_plus) / n, -1: mp.mpf(k_minus) / n, 0: mp.mpf(n - k_plus - k_minus) / n})


def binary_pm_distribution() -> Pdf:
    return {-1: mp.mpf("0.5"), 1: mp.mpf("0.5")}


def discrete_uniform_distribution(start: int, end: int) -> Pdf:
    if end < start:
        raise ValueError("empty uniform distribution")
    return normalize({i: 1 for i in range(start, end + 1)})


def moments_2_4_discrete(D: Mapping[int, object]) -> Tuple[mp.mpf, mp.mpf]:
    E2 = mp.mpf(0)
    E4 = mp.mpf(0)
    for x0, p0 in D.items():
        x = _mp(x0)
        p = _mp(p0)
        E2 += x * x * p
        E4 += x ** 4 * p
    return E2, E4


def moments_2_4_sum(
    discrete_terms: Sequence[Mapping[int, object]] = (),
    gaussian_sigmas: Sequence[object] = (),
) -> Tuple[mp.mpf, mp.mpf]:
    """Return E[X^2], E[X^4] for X=sum independent terms.

    This computes by convolving the discrete part exactly, then adding the
    combined zero-mean Gaussian variance analytically.
    """
    D = sum_distributions(discrete_terms)
    E2, E4 = moments_2_4_discrete(D)
    sigma2 = mp.fsum(_mp(s) ** 2 for s in gaussian_sigmas)
    if sigma2 != 0:
        E4 = E4 + 6 * sigma2 * E2 + 3 * sigma2 ** 2
        E2 = E2 + sigma2
    return E2, E4


def norm_sum_mean_std(n: int, coeff_E2: mp.mpf, coeff_E4: mp.mpf) -> Tuple[mp.mpf, mp.mpf]:
    """Mean/std of sum_i X_i^2 under the coefficient-i.i.d. model."""
    mean = n * coeff_E2
    var = n * max(mp.mpf(0), coeff_E4 - coeff_E2 ** 2)
    return mean, mp.sqrt(var)


def normal_score_collision_probability(
    n: int,
    wrong_E2: mp.mpf,
    wrong_E4: mp.mpf,
    right_E2: mp.mpf,
    right_E4: mp.mpf,
) -> mp.mpf:
    """Pr[wrong squared-norm score <= right squared-norm score]."""
    wrong_mean, wrong_std = norm_sum_mean_std(n, wrong_E2, wrong_E4)
    right_mean, right_std = norm_sum_mean_std(n, right_E2, right_E4)
    denom = mp.sqrt(2 * (wrong_std ** 2 + right_std ** 2))
    if denom == 0:
        return mp.mpf(1) if wrong_mean <= right_mean else mp.mpf(0)
    return mp.mpf("0.5") * mp.erfc((wrong_mean - right_mean) / denom)


def fixed_ternary_vector_probability(D: Mapping[int, object], n: int, k: int, scale: int = 1) -> mp.mpf:
    """iid coefficient probability that a vector is exactly scale*T_{n,k}."""
    p_minus = _mp(D.get(-scale, 0))
    p_plus = _mp(D.get(scale, 0))
    p_zero = _mp(D.get(0, 0))
    if k < 0 or 2 * k > n:
        return mp.mpf(0)
    if p_minus == 0 or p_plus == 0 or (n - 2 * k > 0 and p_zero == 0):
        return mp.mpf(0)
    return mp.mpf(comb(n, k)) * (p_minus ** k) * mp.mpf(comb(n - k, k)) * (p_plus ** k) * (p_zero ** (n - 2 * k))


def support_vector_probability(D: Mapping[int, object], n: int, valid_values: Iterable[int]) -> mp.mpf:
    """iid coefficient probability that every coefficient lies in valid_values."""
    p = mp.fsum(_mp(D.get(v, 0)) for v in valid_values)
    return p ** n


def candidate_count_exact_weight(N: int, B0: int, B1: int, signed: bool = True) -> int:
    values = (2 * B1 if signed else B1)
    return comb(N, B0) * (values ** B0)


def candidate_count_upto_weight(N: int, B0: int, B1: int, signed: bool = True) -> int:
    values = (2 * B1 if signed else B1)
    return sum(comb(N, j) * (values ** j) for j in range(1, B0 + 1))


def log2p(x: mp.mpf) -> str:
    if x <= 0:
        return "-inf"
    return f"{float(mp.log(x, 2)):.6f}"


# ---------------------------------------------------------------------------
# NTRU with Encoding: packing and fixed-weight false-accept checks
# ---------------------------------------------------------------------------

def encoding_dec_fail_packing(
    message_image_size: int | mp.mpf,
    ambient_ring_size: int | mp.mpf,
    list_size_distribution: Mapping[int, object] | None = None,
    worst_case_list_size: int | None = None,
) -> mp.mpf:
    """Heuristic packing estimate max(|M'| L / |R_M| - 1, 0).

    If list_size_distribution is supplied, average omega(L) over it.  Otherwise
    use the worst_case_list_size as a deterministic L.
    """
    M = _mp(message_image_size)
    R = _mp(ambient_ring_size)
    if R <= 0:
        raise ValueError("ambient_ring_size must be positive")

    def omega(L: int) -> mp.mpf:
        return max(M * L / R - 1, mp.mpf(0))

    if list_size_distribution is not None:
        return mp.fsum(_mp(p) * omega(int(L)) for L, p in list_size_distribution.items())
    if worst_case_list_size is None:
        raise ValueError("supply list_size_distribution or worst_case_list_size")
    return omega(int(worst_case_list_size))


def hypergeom_tail_at_least(population: int, marked: int, draws: int, threshold: int) -> mp.mpf:
    """Pr[X >= threshold] for X~Hypergeom(population, marked, draws)."""
    if threshold <= 0:
        return mp.mpf(1)
    if threshold > min(marked, draws):
        return mp.mpf(0)
    denom = mp.mpf(comb(population, draws))
    lo = max(threshold, 0)
    hi = min(marked, draws)
    return mp.fsum(mp.mpf(comb(marked, t)) * mp.mpf(comb(population - marked, draws - t)) / denom for t in range(lo, hi + 1))


def binary_fixed_weight_false_accept_for_shift(
    m: int,
    valid_weight_bound: int,
    shift_weight: int,
    message_weight: int | None = None,
) -> mp.mpf:
    """Pr[wt(mu+a)<=K] for binary mu and fixed nonzero shift a.

    The projected message mu is modeled as a uniform binary vector of fixed
    weight r.  If message_weight is None, the function returns the worst case
    over 0<=r<=K.
    """
    K = valid_weight_bound
    h = shift_weight
    if h < 0 or h > m:
        raise ValueError("invalid shift weight")

    def one_r(r: int) -> mp.mpf:
        # wt(mu+a) = r+h-2 overlap.  Need overlap >= (r+h-K)/2.
        threshold = int(mp.ceil(mp.mpf(r + h - K) / 2))
        return hypergeom_tail_at_least(m, h, r, threshold)

    if message_weight is not None:
        return one_r(message_weight)
    return max(one_r(r) for r in range(0, K + 1))


# GF(2) polynomial routines.  Polynomials are bitsets; bit i is coeff of x^i.
def gf2_deg(a: int) -> int:
    return a.bit_length() - 1


def gf2_mul(a: int, b: int) -> int:
    out = 0
    while b:
        if b & 1:
            out ^= a
        a <<= 1
        b >>= 1
    return out


def gf2_divmod(a: int, b: int) -> Tuple[int, int]:
    if b == 0:
        raise ZeroDivisionError
    q = 0
    r = a
    db = gf2_deg(b)
    while r and gf2_deg(r) >= db:
        sh = gf2_deg(r) - db
        q ^= 1 << sh
        r ^= b << sh
    return q, r


def gf2_mod(a: int, mod_poly: int) -> int:
    return gf2_divmod(a, mod_poly)[1]


def gf2_inv_mod(a: int, mod_poly: int) -> int:
    """Inverse of a modulo mod_poly over GF(2)."""
    r0, r1 = mod_poly, a
    t0, t1 = 0, 1
    while r1 != 0:
        q, r2 = gf2_divmod(r0, r1)
        r0, r1 = r1, r2
        t0, t1 = t1, t0 ^ gf2_mul(q, t1)
    if r0 != 1:
        raise ValueError("polynomial is not invertible modulo mod_poly")
    return gf2_mod(t0, mod_poly)


def poly_from_exponents_mod2(exponents: Sequence[int], modulus_degree: int) -> int:
    """Build a polynomial in GF(2)[x]/(x^m+1)."""
    mask = 0
    for e in exponents:
        mask ^= 1 << (e % modulus_degree)
    return mask


def cyclic_shift(mask: int, shift: int, m: int) -> int:
    shift %= m
    allbits = (1 << m) - 1
    return ((mask << shift) | (mask >> (m - shift))) & allbits if shift else (mask & allbits)


def min_cyclic_shift_weight(mask: int, m: int, max_terms: int = 1) -> int:
    """Minimum weight of XOR of 1..max_terms cyclic shifts of mask.

    Exact for max_terms<=2.  For larger max_terms this function raises, because
    exact enumeration becomes parameter-specific and should be done with the
    concrete auxiliary-filtered shift family rather than all shifts.
    """
    shifts = [cyclic_shift(mask, i, m) for i in range(m)]
    best = min(x.bit_count() for x in shifts)
    if max_terms == 1:
        return best
    if max_terms == 2:
        for i in range(m):
            si = shifts[i]
            for j in range(i + 1, m):
                best = min(best, (si ^ shifts[j]).bit_count())
        return best
    raise NotImplementedError("exact all-shift search is implemented only for max_terms <= 2")


def end_e_inverse_shift_report(n: int, w_exponents: Sequence[int], valid_weight_bound: int, max_terms: int = 1) -> Tuple[int, mp.mpf]:
    """Concrete binary fixed-weight check for R_{2,x^{n/2}+1}.

    Returns (minimum shift weight, worst-case false-accept probability for one
    wrong candidate shift).  w_exponents are reduced modulo n/2 over GF(2).
    """
    m = n // 2
    mod_poly = (1 << m) | 1  # x^m + 1 over GF(2)
    w = poly_from_exponents_mod2(w_exponents, m)
    winv = gf2_inv_mod(w, mod_poly)
    h_min = min_cyclic_shift_weight(winv, m, max_terms=max_terms)
    p = binary_fixed_weight_false_accept_for_shift(m, valid_weight_bound, h_min, message_weight=None)
    return h_min, p


# ---------------------------------------------------------------------------
# NTRU with Trapdoor: distributional verification collision estimates
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TrapdoorDecFailParams:
    name: str
    n: int
    q_effective: int
    Dg: Pdf
    Df: Pdf
    Ds: Pdf
    Dmu: Pdf
    Bgf0: int
    Bgf1: int
    BGF0: int
    BGF1: int
    Gstar_std: mp.mpf
    Fstar_std: mp.mpf
    adjust_fac: int = 1
    gf_use_mu_score: bool = True
    gf_use_s_score: bool = True
    GF_rule: str = "fixed_ternary"  # fixed_ternary or support
    mu_fixed_k: int | None = None
    s_fixed_k: int | None = None
    GF_use_mu_check: bool = True
    GF_use_s_check: bool = True
    GF_scale_base_by_adjust: bool = True
    candidate_domain_gf: int | None = None
    candidate_domain_GF: int | None = None
    signed_candidate_multiplier: bool = False
    exact_weight_multiplier: bool = True


def trapdoor_gf_collision(params: TrapdoorDecFailParams) -> Tuple[mp.mpf, Dict[str, mp.mpf]]:
    """Estimate epsilon_gf by the norm-score collision heuristic."""
    n = params.n
    q = params.q_effective
    a = params.adjust_fac

    Dmu_scaled = scale_distribution(params.Dmu, a)
    Ds_scaled = scale_distribution(params.Ds, a)

    # Fractional interaction from alpha*E_gf.  The paper model is uniform on
    # [-1/2,1/2).  The old notebook approximated this by i/q', but the moment
    # needed here is simply Var(U)=1/12.
    _ = q * a  # kept only to make the dependency visible in the parameter record
    uniform_E2 = mp.mpf(1) / 12
    Dg_E2, _ = moments_2_4_discrete(params.Dg)
    e_mid_std = mp.sqrt(n * uniform_E2 * Dg_E2)

    # wrong-gf-guess coefficient model
    mu_wrong_E2, mu_wrong_E4 = moments_2_4_sum(
        [Dmu_scaled], [params.Gstar_std, e_mid_std]
    )
    s_wrong_E2, s_wrong_E4 = moments_2_4_sum(
        [Ds_scaled], [params.Fstar_std, e_mid_std]
    )

    # right-gf-guess coefficient model, still with unknown E_GF contribution.
    mu_right_E2, mu_right_E4 = moments_2_4_sum(
        [Dmu_scaled, scale_distribution(params.Dg, params.BGF1)], []
    )
    s_right_E2, s_right_E4 = moments_2_4_sum(
        [Ds_scaled, scale_distribution(params.Df, params.BGF1)], []
    )

    p_mu = normal_score_collision_probability(n, mu_wrong_E2, mu_wrong_E4, mu_right_E2, mu_right_E4)
    p_s = normal_score_collision_probability(n, s_wrong_E2, s_wrong_E4, s_right_E2, s_right_E4)

    per_candidate = mp.mpf(1)
    if params.gf_use_mu_score:
        per_candidate *= p_mu
    if params.gf_use_s_score:
        per_candidate *= p_s

    N = params.candidate_domain_gf or n
    if params.exact_weight_multiplier:
        mult = candidate_count_exact_weight(N, params.Bgf0, params.Bgf1, signed=params.signed_candidate_multiplier)
    else:
        mult = candidate_count_upto_weight(N, params.Bgf0, params.Bgf1, signed=params.signed_candidate_multiplier)
    eps = min(mp.mpf(1), mp.mpf(mult) * per_candidate)
    details = {
        "p_mu_single": p_mu,
        "p_s_single": p_s,
        "candidate_multiplier": mp.mpf(mult),
        "e_mid_std": e_mid_std,
    }
    return eps, details


def trapdoor_GF_collision(params: TrapdoorDecFailParams) -> Tuple[mp.mpf, Dict[str, mp.mpf]]:
    """Estimate epsilon_GF by fixed-distribution/range false-accept heuristic."""
    n = params.n
    a = params.adjust_fac
    if params.GF_scale_base_by_adjust:
        Dmu_base = scale_distribution(params.Dmu, a)
        Ds_base = scale_distribution(params.Ds, a)
    else:
        # This reproduces the submitted frontier notebook for adjust_fac=2:
        # it checks whether mu + Delta*g lands in the scaled ternary alphabet.
        Dmu_base = params.Dmu
        Ds_base = params.Ds

    mu_final_dist = add_distribution(Dmu_base, scale_distribution(params.Dg, params.BGF1))
    s_final_dist = add_distribution(Ds_base, scale_distribution(params.Df, params.BGF1))

    if params.GF_rule == "fixed_ternary":
        if params.mu_fixed_k is None or params.s_fixed_k is None:
            raise ValueError("fixed_ternary GF_rule needs mu_fixed_k and s_fixed_k")
        p_mu = fixed_ternary_vector_probability(mu_final_dist, n, params.mu_fixed_k, scale=a)
        p_s = fixed_ternary_vector_probability(s_final_dist, n, params.s_fixed_k, scale=a)
    elif params.GF_rule == "support":
        p_mu = support_vector_probability(mu_final_dist, n, params.Dmu.keys())
        p_s = support_vector_probability(s_final_dist, n, params.Ds.keys())
    else:
        raise ValueError(f"unknown GF_rule: {params.GF_rule}")

    per_candidate = mp.mpf(1)
    if params.GF_use_mu_check:
        per_candidate *= p_mu
    if params.GF_use_s_check:
        per_candidate *= p_s

    N = params.candidate_domain_GF or n
    if params.exact_weight_multiplier:
        mult = candidate_count_exact_weight(N, params.BGF0, params.BGF1, signed=params.signed_candidate_multiplier)
    else:
        mult = candidate_count_upto_weight(N, params.BGF0, params.BGF1, signed=params.signed_candidate_multiplier)
    eps = min(mp.mpf(1), mp.mpf(mult) * per_candidate)
    details = {
        "p_mu_single": p_mu,
        "p_s_single": p_s,
        "candidate_multiplier": mp.mpf(mult),
    }
    return eps, details


def trapdoor_dec_fail(params: TrapdoorDecFailParams) -> Dict[str, object]:
    eps_gf, gf_details = trapdoor_gf_collision(params)
    eps_GF, GF_details = trapdoor_GF_collision(params)
    return {
        "name": params.name,
        "epsilon_gf": eps_gf,
        "epsilon_GF": eps_GF,
        "dec_fails": min(mp.mpf(1), eps_gf + eps_GF),
        "gf_details": gf_details,
        "GF_details": GF_details,
    }


# ---------------------------------------------------------------------------
# Built-in cases from the paper/notebooks
# ---------------------------------------------------------------------------

def trapdoor_frontier_q39_params() -> TrapdoorDecFailParams:
    n = 512
    q = 39
    k = 21
    # same sparse ternary model as the frontier notebook
    D = ternary_coefficient_distribution(n, k)
    GF_star_std = (mp.e * q / 2) / mp.sqrt(2 * k + 2 * k) / mp.sqrt(2 * n)
    # notebook chooses the first k_GF with std(T_{n,k_GF}) > GF_star_std, but
    # Dec.fail only needs the std model above.
    return TrapdoorDecFailParams(
        name="trapdoor-frontier-q39",
        n=n,
        q_effective=q,
        Dg=D,
        Df=D,
        Ds=D,
        Dmu=D,
        Bgf0=4,
        Bgf1=1,
        BGF0=2,
        BGF1=1,
        Gstar_std=GF_star_std,
        Fstar_std=GF_star_std,
        adjust_fac=2,
        gf_use_mu_score=True,
        gf_use_s_score=True,
        GF_rule="fixed_ternary",
        mu_fixed_k=k,
        s_fixed_k=k,
        GF_use_mu_check=True,
        GF_use_s_check=True,
        GF_scale_base_by_adjust=False,
        signed_candidate_multiplier=False,
        exact_weight_multiplier=True,
    )


def end_t_512_params() -> TrapdoorDecFailParams:
    n = 512
    q_effective = 257 * 2  # same q used in the submitted notebook for DFR arithmetic
    kg = kf = 72
    Ds = {-1: 1/2, 1: 1/2}
    Df = ternary_coefficient_distribution(n, kf)
    Dg = ternary_coefficient_distribution(n, kg)
    Dmu = normalize({-3: 1/4, -1: 1/4, 1: 1/4, 3: 1/4})
    return TrapdoorDecFailParams(
        name="END-T-512",
        n=n,
        q_effective=q_effective,
        Dg=Dg,
        Df=Df,
        Ds=Ds,
        Dmu=Dmu,
        Bgf0=2,
        Bgf1=1,
        BGF0=2,
        BGF1=1,
        Gstar_std=mp.mpf("1.0641"),
        Fstar_std=mp.mpf("0.3644"),
        adjust_fac=1,
        # The submitted practical section says the implementation can evaluate
        # only the s terms when this probability is tiny.  Keep that as default.
        gf_use_mu_score=False,
        gf_use_s_score=True,
        GF_rule="support",
        GF_use_mu_check=False,
        GF_use_s_check=True,
        signed_candidate_multiplier=False,
        exact_weight_multiplier=True,
    )


def end_t_1024_params() -> TrapdoorDecFailParams:
    n = 1024
    q_effective = 257 * 2
    kg = kf = 96
    Ds = binary_pm_distribution()
    Df = ternary_coefficient_distribution(n, kf)
    Dg = ternary_coefficient_distribution(n, kg)
    Dmu = normalize({-2: mp.mpf(77) / 257, 0: mp.mpf(90) / 257, 2: mp.mpf(90) / 257})
    return TrapdoorDecFailParams(
        name="END-T-1024",
        n=n,
        q_effective=q_effective,
        Dg=Dg,
        Df=Df,
        Ds=Ds,
        Dmu=Dmu,
        Bgf0=2,
        Bgf1=1,
        BGF0=1,
        BGF1=1,
        Gstar_std=mp.mpf("0.574"),
        Fstar_std=mp.mpf("0.275"),
        adjust_fac=1,
        gf_use_mu_score=False,
        gf_use_s_score=True,
        GF_rule="support",
        GF_use_mu_check=False,
        GF_use_s_check=True,
        signed_candidate_multiplier=False,
        exact_weight_multiplier=True,
    )


def print_trapdoor_case(params: TrapdoorDecFailParams) -> None:
    result = trapdoor_dec_fail(params)
    print(f"{params.name}")
    print(f"  epsilon_gf      = 2^({log2p(result['epsilon_gf'])})")
    print(f"  epsilon_GF      = 2^({log2p(result['epsilon_GF'])})")
    print(f"  Pr[Dec. fails] = 2^({log2p(result['dec_fails'])})")
    gf = result["gf_details"]  # type: ignore[index]
    GF = result["GF_details"]  # type: ignore[index]
    print(f"  gf single: p_mu=2^({log2p(gf['p_mu_single'])}), p_s=2^({log2p(gf['p_s_single'])}), mult=2^({log2p(gf['candidate_multiplier'])})")
    print(f"  GF single: p_mu=2^({log2p(GF['p_mu_single'])}), p_s=2^({log2p(GF['p_s_single'])}), mult=2^({log2p(GF['candidate_multiplier'])})")


def print_encoding_examples() -> None:
    print("encoding-packing-frontier-example")
    # The 512/q=109 frontier text claims an exact disjointness check, so the
    # verification term is zero.  The packing helper is shown for reproducibility.
    p = encoding_dec_fail_packing(message_image_size=2**128, ambient_ring_size=2**128, worst_case_list_size=1)
    print(f"  packing omega with |M'|=|R_M| and L=1: 2^({log2p(p)})")
    print("  use exact_disjointness=True in a parameter file to record Pr_encoding[Dec. fails]=0")

    print("END-E-512 binary fixed-weight check")
    # w = x^256 + 3x^159 + x^128 + 3x^97 + 1; coefficients are reduced mod 2.
    # The relevant shift in the decoder is w^{-1} times an error residue in
    # R_{2,x^{256}+1}; the inverse is dense.
    h, pfa = end_e_inverse_shift_report(n=512, w_exponents=[256, 159, 128, 97, 0], valid_weight_bound=80, max_terms=1)
    print(f"  min wt(x^i w^(-1)) in R_(2,x^256+1): {h}")
    print(f"  one-shift false accept, worst r<=80: 2^({log2p(pfa)})")

    print("END-E-1024 binary fixed-weight check")
    h, pfa_k = end_e_inverse_shift_report(n=1024, w_exponents=[512, 455, 397, 256, 0], valid_weight_bound=128, max_terms=1)
    _, pfa_2k = end_e_inverse_shift_report(n=1024, w_exponents=[512, 455, 397, 256, 0], valid_weight_bound=256, max_terms=1)
    print(f"  min wt(x^i w^(-1)) in R_(2,x^512+1): {h}")
    print(f"  one-shift false accept with residue-weight bound K=128: 2^({log2p(pfa_k)})")
    print(f"  one-shift false accept with loose bound K=256: 2^({log2p(pfa_2k)})")


CASES = {
    "trapdoor-frontier": trapdoor_frontier_q39_params,
    "end-t-512": end_t_512_params,
    "end-t-1024": end_t_1024_params,
}


def main() -> None:
    ap = argparse.ArgumentParser(description="Estimate Pr[Dec. fails] verification/collision terms.")
    ap.add_argument("--case", choices=["all", "encoding"] + sorted(CASES), default="all")
    ap.add_argument("--signed-candidates", action="store_true", help="multiply trapdoor union bounds by (2B1)^B0 instead of B1^B0")
    ap.add_argument("--full-candidate-count", action="store_true", help="use sum_{j<=B0} C(N,j)(values)^j instead of exact weight B0")
    args = ap.parse_args()

    if args.case in ("all", "encoding"):
        print_encoding_examples()
        if args.case == "all":
            print()

    selected = []
    if args.case == "all":
        selected = ["trapdoor-frontier", "end-t-512", "end-t-1024"]
    elif args.case in CASES:
        selected = [args.case]

    for idx, name in enumerate(selected):
        params = CASES[name]()
        if args.signed_candidates or args.full_candidate_count:
            params = TrapdoorDecFailParams(
                **{**params.__dict__,
                   "signed_candidate_multiplier": bool(args.signed_candidates),
                   "exact_weight_multiplier": not bool(args.full_candidate_count)}
            )
        print_trapdoor_case(params)
        if idx != len(selected) - 1:
            print()


if __name__ == "__main__":
    main()
