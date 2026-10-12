# -*- coding: utf-8 -*-
"""
outils_tendance.py
Fonctions de test de tendance partagees par 16_analyse.py et 18_tendances.py.

  tendance_serie : pente de Sen et p-values (Mann-Kendall classique, PW, TFPW)
  p_3pw          : p-value 3PW (Collaud Coen et al. 2020) = la plus grande de
                   PW et TFPW : significatif seulement si les deux le sont
  fdr            : p-values ajustees de Benjamini-Hochberg

Memes calculs que le package pymannkendall (valide dans 18_tendances.py),
mais les annees manquantes sont gerees selon leur vraie date : le
pre-whitening n'utilise que des annees consecutives.
"""

import math

import numpy as np


def p_mann_kendall(y):
    """p-value bilaterale du test de Mann-Kendall (y dans l'ordre du temps).
    Variance avec correction des ex aequo, utile pour les comptages de jours."""
    n = len(y)
    i, j = np.triu_indices(n, k=1)
    s = np.sign(y[j] - y[i]).sum()
    _, t = np.unique(y, return_counts=True)
    var = (n * (n - 1) * (2 * n + 5) - np.sum(t * (t - 1) * (2 * t + 5))) / 18.0
    if s == 0 or var <= 0:
        return 1.0
    z = (s - np.sign(s)) / math.sqrt(var)
    return math.erfc(abs(z) / math.sqrt(2.0))


def rho1(r, suite):
    """Autocorrelation de lag 1, sur les seules annees consecutives."""
    rc = r - r.mean()
    if not np.any(rc):
        return 0.0
    return np.sum(rc[1:][suite] * rc[:-1][suite]) / np.sum(rc * rc)


def tendance_serie(annees, x):
    """Pente de Sen (par an), p Mann-Kendall classique, p PW, p TFPW."""
    ok = ~np.isnan(x)
    t, v = annees[ok], x[ok]
    n = len(v)
    i, j = np.triu_indices(n, k=1)
    pente = np.median((v[j] - v[i]) / (t[j] - t[i]))
    p_classique = p_mann_kendall(v)

    suite = (t[1:] - t[:-1]) == 1          # paires d'annees consecutives

    # PW : on retire l'autocorrelation de la serie brute
    rho = rho1(v, suite)
    y_pw = v[1:][suite] - rho * v[:-1][suite]

    # TFPW : meme chose sur la serie sans tendance, puis la tendance est remise
    r = v - pente * t
    rho = rho1(r, suite)
    y_tfpw = r[1:][suite] - rho * r[:-1][suite] + pente * t[1:][suite]

    if len(y_pw) < 3:
        return pente, p_classique, np.nan, np.nan
    return pente, p_classique, p_mann_kendall(y_pw), p_mann_kendall(y_tfpw)


def fdr(p):
    """p-values ajustees de Benjamini-Hochberg (NaN ignores)."""
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = len(pv)
    ordre = np.argsort(pv)
    ajuste = pv[ordre] * m / np.arange(1, m + 1)
    ajuste = np.minimum.accumulate(ajuste[::-1])[::-1].clip(max=1.0)
    res = np.empty(m)
    res[ordre] = ajuste
    out[ok] = res
    return out


def p_3pw(p_pw, p_tfpw):
    """3PW : la tendance n'est significative que si PW et TFPW le sont."""
    return np.fmax(p_pw, p_tfpw)

