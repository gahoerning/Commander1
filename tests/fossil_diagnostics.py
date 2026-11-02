"""Rank-normalized split R-hat and ESS for independent FOSSIL chains.

Numerical dependency: NumPy. Diagnostic formulas follow the
rank-normalized split-chain approach; they are acceptance diagnostics, not
evidence that a model is identifiable or correctly specified.
"""
import numpy as np
from statistics import NormalDist


def split(chains):
    chains = np.asarray(chains, dtype=float)
    if chains.ndim != 2 or chains.shape[0] < 2 or chains.shape[1] < 20:
        raise ValueError("At least two independent chains with 20 post-warmup draws are required")
    if not np.isfinite(chains).all():
        raise ValueError("Nonfinite chain values")
    half = chains.shape[1] // 2
    return np.concatenate((chains[:, :half], chains[:, -half:]))


def rank_normalize(values):
    flat = values.reshape(-1)
    _, inverse, counts = np.unique(flat, return_inverse=True, return_counts=True)
    ends = np.cumsum(counts)
    ranks = (ends - (counts - 1) / 2.0)[inverse]
    probabilities = (ranks - 3.0 / 8.0) / (ranks.size + 1.0 / 4.0)
    return np.array([NormalDist().inv_cdf(p) for p in probabilities]).reshape(values.shape)


def _rhat(values):
    n = values.shape[1]
    within = np.var(values, axis=1, ddof=1).mean()
    between = n * np.var(values.mean(axis=1), ddof=1)
    if within == 0:
        return float("inf") if between else float("nan")
    return float(np.sqrt(((n - 1) / n * within + between / n) / within))


def _ess(values):
    m, n = values.shape
    centered = values - values.mean(axis=1, keepdims=True)
    transform = np.fft.rfft(centered, n=2 * n, axis=1)
    acov = np.fft.irfft(transform * np.conjugate(transform), n=2 * n, axis=1)[:, :n] / n
    within = acov[:, 0].mean() * n / (n - 1)
    var_plus = (n - 1) / n * within + np.var(values.mean(axis=1), ddof=1)
    if var_plus <= 0:
        return 0.0
    rho = 1 - (within - acov.mean(axis=0)) / var_plus
    rho[0] = 1.0
    pairs = []
    previous = float("inf")
    for lag in range(0, n - 1, 2):
        pair = rho[lag] + rho[lag + 1]
        if pair < 0:
            break
        previous = min(previous, pair)
        pairs.append(previous)
    tau = max(1.0, -1.0 + 2.0 * sum(pairs))
    return float(min(m * n, m * n / tau))


def diagnostics(chains):
    raw = split(chains)
    ranked = rank_normalize(raw)
    folded = rank_normalize(np.abs(raw - np.median(raw)))
    low, high = np.quantile(raw, [0.05, 0.95])
    return dict(rhat=max(_rhat(ranked), _rhat(folded)), bulk_ess=_ess(ranked),
                tail_ess=min(_ess((raw <= low).astype(float)),
                             _ess((raw >= high).astype(float))),
                draws=int(raw.size), mean=float(raw.mean()), std=float(raw.std(ddof=1)))
