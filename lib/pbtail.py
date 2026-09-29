#!/usr/bin/env python3
"""Canonical Poisson-binomial upper tail. Single source of truth."""
import numpy as np

def poisson_binomial_pmf(probs):
    p = np.asarray(probs, dtype=np.longdouble)
    dp = np.zeros(len(p) + 1, dtype=np.longdouble)
    dp[0] = np.longdouble(1.0)
    for i, q in enumerate(p):
        dp[1:i+2] = dp[1:i+2] * (1 - q) + dp[0:i+1] * q
        dp[0] *= (1 - q)
    return dp

def poisson_binomial_upper_tail(probs, k):
    return float(poisson_binomial_pmf(probs)[k:].sum())

def poisson_binomial_moments(probs):
    p = np.asarray(probs, dtype=np.longdouble)
    return float(p.sum()), float(np.sqrt((p * (1 - p)).sum()))
