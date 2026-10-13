# -*- coding: utf-8 -*-
"""
18_tendances.py
A lancer depuis la racine du projet, environnement geo :
    python 18_tendances.py

Tendances par pixel, climat (1961-2024) et NDVI estival (2001-2024) :
  - pente de Sen, exprimee par decennie ;
  - significativite par la methode 3PW (Collaud Coen et al. 2020) : le test
    de Mann-Kendall est fait deux fois, apres deux facons de retirer
    l'autocorrelation de lag 1 (PW, von Storch 1995 ; TFPW, Yue et al. 2002).
    Une tendance n'est significative que si les DEUX tests le sont.
    Raison : une autocorrelation positive fait trouver trop de tendances au
    test classique, et TFPW seul n'y remedie pas (trop de faux positifs) ;
  - correction des tests multiples par carte : Benjamini-Hochberg, FDR 10 %
    (Wilks 2016 : environ 5 % de risque global pour un champ spatialement
    correle).

Un pixel n'est retenu que s'il a au moins 75 % d'annees valides.

Fonctions de test dans outils_tendance.py (meme dossier). Une validation
les compare au package pymannkendall sur la temperature estivale (pixels
sans annee manquante).

Entrees : 03_arcgis/*.nc (scripts 09 a 11)
Sorties, dans le projet ArcGIS Belgium/05_tendances (memes noms qu'avant) :
  trend_<ind>_<saison>_pente_dec.tif         pente de Sen par decennie
  trend_<ind>_<saison>_pvalue.tif            p-value 3PW (la plus grande de PW et TFPW)
  trend_<ind>_<saison>_pvalue_fdr.tif        p-value ajustee (Benjamini-Hochberg)
  trend_<ind>_<saison>_pente_dec_signif.tif  pente, pixels avec p ajustee < 0,10
Resume : 05_resultats/tendances_resume.csv
Validation : 05_resultats/validation_tendances_pymannkendall.csv
"""

import math
import os

import numpy as np
import pandas as pd
import rioxarray  # noqa: F401  (active .rio pour ecrire les GeoTIFF)
import xarray as xr

from outils_tendance import fdr, p_3pw, tendance_serie

RACINE = os.getcwd()   # racine du projet
PORTFOLIO = os.path.dirname(RACINE)   # le projet ArcGIS "Belgium" est a cote
ENTREE = os.path.join(RACINE, "03_arcgis")
SORTIE = os.path.join(PORTFOLIO, "Belgium", "05_tendances")
RESULTATS = os.path.join(RACINE, "05_resultats")
os.makedirs(SORTIE, exist_ok=True)
os.makedirs(RESULTATS, exist_ok=True)

PART_MIN = 0.75   # part minimale d'annees valides par pixel
SEUIL_P = 0.05    # seuil d'un test isole (pour comparaison)
SEUIL_FDR = 0.10  # taux de fausses decouvertes par carte

# (fichier NetCDF dans 03_arcgis, variable, nom de sortie)
TENDANCES = [
    ("era5_BE_DJF.nc", "t2m", "trend_t2m_DJF"),
    ("era5_BE_MAM.nc", "t2m", "trend_t2m_MAM"),
    ("era5_BE_JJA.nc", "t2m", "trend_t2m_JJA"),
    ("era5_BE_SON.nc", "t2m", "trend_t2m_SON"),
    ("era5_BE_JJA.nc", "swvl2", "trend_swvl2_JJA"),
    ("eobs_fd_BE_DJF.nc", "fd", "trend_fd_DJF"),
    ("eobs_su_BE_JJA.nc", "su", "trend_su_JJA"),
    ("eobs_su30_BE_JJA.nc", "su30", "trend_su30_JJA"),
    ("eobs_txx_BE_JJA.nc", "txx", "trend_txx_JJA"),
    ("eobs_prcptot_BE_DJF.nc", "prcptot", "trend_prcptot_DJF"),
    ("eobs_prcptot_BE_MAM.nc", "prcptot", "trend_prcptot_MAM"),
    ("modis_ndvi_BE_JJA.nc", "ndvi", "trend_ndvi_JJA"),
]


def ecrire_tif(tableau, gabarit, chemin):
    da = xr.DataArray(tableau.astype("float32"), dims=("latitude", "longitude"),
                      coords={"latitude": gabarit["latitude"], "longitude": gabarit["longitude"]})
    # Nord en haut : E-OBS range ses latitudes du sud au nord ; ecrites telles
    # quelles, elles donnent un GeoTIFF "a l'envers" qu'ArcGIS lit mal
    # (decalage, trous sur la cote). Meme correction que 12_normales.py.
    da = da.sortby("latitude", ascending=False)
    da = da.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude")
    da = da.rio.write_crs("EPSG:4326").rio.write_nodata(np.nan)
    da.rio.to_raster(chemin)


def calculer(fichier, variable):
    da = xr.open_dataset(os.path.join(ENTREE, fichier))[variable]
    da = da.transpose("time", "latitude", "longitude")
    annees = da["time"].dt.year.values.astype("float64")
    v = da.values.astype("float64")
    n_min = math.ceil(PART_MIN * len(annees))
    pente = np.full(v.shape[1:], np.nan)
    p_cl = np.full(v.shape[1:], np.nan)
    p_pw = np.full(v.shape[1:], np.nan)
    p_tf = np.full(v.shape[1:], np.nan)
    valides = np.sum(~np.isnan(v), axis=0)
    for iy, ix in zip(*np.nonzero(valides >= n_min)):
        pente[iy, ix], p_cl[iy, ix], p_pw[iy, ix], p_tf[iy, ix] = tendance_serie(annees, v[:, iy, ix])
    return da, annees, pente, p_cl, p_pw, p_tf, p_3pw(p_pw, p_tf)


# ===========================================================================
# 1. Validation contre pymannkendall (temperature estivale)
# ===========================================================================
print("=== Validation : calcul maison vs pymannkendall, t2m JJA ===")
try:
    import pymannkendall as mk
    da, annees, pente, p_cl, p_pw, p_tf, _ = calculer("era5_BE_JJA.nc", "t2m")
    v = da.values.astype("float64")
    complets = np.all(~np.isnan(v), axis=0)
    e_pente, e_cl, e_pw, e_tf = [], [], [], []
    for iy, ix in zip(*np.nonzero(complets)):
        serie = v[:, iy, ix]
        ref = mk.trend_free_pre_whitening_modification_test(serie)
        e_pente.append(abs(pente[iy, ix] - ref.slope))
        e_tf.append(abs(p_tf[iy, ix] - ref.p))
        e_pw.append(abs(p_pw[iy, ix] - mk.pre_whitening_modification_test(serie).p))
        e_cl.append(abs(p_cl[iy, ix] - mk.original_test(serie).p))
    lignes = [{"grandeur": nom, "pixels_compares": len(e), "ecart_max": max(e)}
              for nom, e in [("pente", e_pente), ("p_classique", e_cl),
                             ("p_pw", e_pw), ("p_tfpw", e_tf)]]
    pd.DataFrame(lignes).to_csv(os.path.join(RESULTATS, "validation_tendances_pymannkendall.csv"),
                                sep=";", decimal=",", index=False)
    for l in lignes:
        print("  %-12s : %d pixels, ecart max %.1e" % (l["grandeur"], l["pixels_compares"], l["ecart_max"]))
except ImportError:
    print("  pymannkendall absent : validation sautee (conda install -c conda-forge pymannkendall)")

# ===========================================================================
# 2. Tendances et cartes
# ===========================================================================
print("\n=== Tendances ===")
print("  %-18s %8s %9s | %%sig MK  %%sig 3PW  %%sig 3PW+FDR" % ("", "pixels", "pente med"))
resume = []
for fichier, variable, nom in TENDANCES:
    if not os.path.exists(os.path.join(ENTREE, fichier)):
        print("  !! introuvable, ignore : %s" % fichier)
        continue
    da, annees, pente, p_cl, _, _, p_sig = calculer(fichier, variable)
    p_adj = fdr(p_sig)
    pente_dec = pente * 10
    signif = np.where(p_adj < SEUIL_FDR, pente_dec, np.nan)

    ecrire_tif(pente_dec, da, os.path.join(SORTIE, nom + "_pente_dec.tif"))
    ecrire_tif(p_sig, da, os.path.join(SORTIE, nom + "_pvalue.tif"))
    ecrire_tif(p_adj, da, os.path.join(SORTIE, nom + "_pvalue_fdr.tif"))
    ecrire_tif(signif, da, os.path.join(SORTIE, nom + "_pente_dec_signif.tif"))

    ok = ~np.isnan(pente_dec)
    n = int(ok.sum())
    ligne = {"tendance": nom, "pixels": n,
             "pente_dec_mediane": float(np.median(pente_dec[ok])),
             "pente_dec_min": float(pente_dec[ok].min()),
             "pente_dec_max": float(pente_dec[ok].max()),
             "pct_sig_mk_classique": 100.0 * np.sum(p_cl[ok] < SEUIL_P) / n,
             "pct_sig_3pw": 100.0 * np.sum(p_sig[ok] < SEUIL_P) / n,
             "pct_sig_3pw_fdr": 100.0 * np.sum(p_adj[ok] < SEUIL_FDR) / n}
    resume.append(ligne)
    print("  %-18s %8d %+9.3f | %7.1f  %8.1f  %12.1f" % (
        nom, n, ligne["pente_dec_mediane"], ligne["pct_sig_mk_classique"],
        ligne["pct_sig_3pw"], ligne["pct_sig_3pw_fdr"]))

pd.DataFrame(resume).round(4).to_csv(os.path.join(RESULTATS, "tendances_resume.csv"),
                                     sep=";", decimal=",", index=False)
print("\n  pentes par decennie (degC, jours, mm, m3/m3 ou NDVI selon l'indicateur)")
print("  -> GeoTIFF dans %s" % SORTIE)
print("  -> 05_resultats/tendances_resume.csv")
