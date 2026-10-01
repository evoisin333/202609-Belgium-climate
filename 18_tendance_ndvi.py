# -*- coding: utf-8 -*-
"""
18_tendance_ndvi.py
A lancer depuis la racine du projet, environnement geo :
    python 18_tendance_ndvi.py

Tendance du NDVI estival 2001-2024 : pente de Sen et test de Mann-Kendall,
calcules pixel par pixel en Python (le calcul equivalent sous ArcGIS
bloquait a l'ecriture des TIF).

Meme methode que Generate Trend Raster (MANN-KENDALL) sous ArcGIS :
  - pente de Sen : mediane des pentes entre toutes les paires d'annees
  - test de Mann-Kendall bilateral, avec correction de continuite
  - variance de S sans correction des ex aequo, negligeable pour le NDVI
    (valeur continue)
  - pente exprimee par decennie (x 10), comme les autres indicateurs

Ajout : un pixel n'est retenu que s'il compte au moins N_MIN annees valides.
Sur 2 ou 3 annees, la pente de Sen peut prendre n'importe quelle valeur.

Une validation compare d'abord cette implementation a la tendance ArcGIS de
la temperature, sur le pixel de Bruxelles.

Sorties (NetCDF, a convertir en TIF sous ArcGIS, voir fin de script) :
  Belgium\\05_tendances\\trend_ndvi_JJA_pente_dec.nc
  Belgium\\05_tendances\\trend_ndvi_JJA_pvalue.nc
  Belgium\\05_tendances\\trend_ndvi_JJA_pente_dec_signif.nc
"""

import os
import warnings

import numpy as np
import xarray as xr

RACINE = os.getcwd()   # racine du projet
PORTFOLIO = os.path.dirname(RACINE)   # le projet ArcGIS "Belgium" est a cote
ENTREE = os.path.join(RACINE, "03_arcgis", "modis_ndvi_BE_JJA.nc")
ERA5 = os.path.join(RACINE, "03_arcgis", "era5_BE_JJA.nc")
SORTIE = os.path.join(PORTFOLIO, "Belgium", "05_tendances")

N_MIN = 18        # annees valides minimum sur 24 (75 %)
SEUIL = 0.05      # significativite
LIGNES = 25       # lignes traitees a la fois (memoire)

CRS_ATTRS = {
    "grid_mapping_name": "latitude_longitude",
    "longitude_of_prime_meridian": 0.0,
    "semi_major_axis": 6378137.0,
    "inverse_flattening": 298.257223563,
    "epsg_code": "EPSG:4326",
}


def erfc(x):
    """Fonction d'erreur complementaire (Abramowitz et Stegun 7.1.26,
    erreur < 1.5e-7), pour x >= 0. Evite une dependance a scipy."""
    t = 1.0 / (1.0 + 0.3275911 * x)
    y = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741
            + t * (-1.453152027 + t * 1.061405429))))
    return y * np.exp(-x * x)


def mann_kendall(v, annees):
    """v : tableau (annees, lignes, colonnes). Retourne pente/an, p-value, n."""
    T = v.shape[0]
    I, J = np.triu_indices(T, k=1)
    dt = (annees[J] - annees[I]).astype("float32")

    n = np.sum(~np.isnan(v), axis=0).astype("float64")
    pente = np.full(v.shape[1:], np.nan, dtype="float32")
    pval = np.full(v.shape[1:], np.nan, dtype="float32")

    for y0 in range(0, v.shape[1], LIGNES):
        bloc = v[:, y0:y0 + LIGNES, :]
        d = bloc[J] - bloc[I]                              # paires x lignes x colonnes
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            pente[y0:y0 + LIGNES] = np.nanmedian(d / dt[:, None, None], axis=0)

        s = np.nansum(np.sign(d), axis=0)                  # paires avec NaN ignorees
        nb = n[y0:y0 + LIGNES]
        var = nb * (nb - 1) * (2 * nb + 5) / 18.0
        with np.errstate(invalid="ignore", divide="ignore"):
            z = np.where(s > 0, (s - 1) / np.sqrt(var),
                np.where(s < 0, (s + 1) / np.sqrt(var), 0.0))
            p = erfc(np.abs(z) / np.sqrt(2.0))
        p[nb < 3] = np.nan
        pval[y0:y0 + LIGNES] = p

    return pente, pval, n


def ecrire(tableau, gabarit, nom, attrs):
    da = xr.DataArray(tableau.astype("float32"),
                      coords={"latitude": gabarit["latitude"],
                              "longitude": gabarit["longitude"]},
                      dims=("latitude", "longitude"), name=nom)
    da.attrs = dict(attrs, grid_mapping="crs")
    ds = da.to_dataset()
    ds["latitude"].attrs = {"standard_name": "latitude", "units": "degrees_north", "axis": "Y"}
    ds["longitude"].attrs = {"standard_name": "longitude", "units": "degrees_east", "axis": "X"}
    ds["crs"] = xr.DataArray(np.int32(0), attrs=CRS_ATTRS)
    ds.attrs = {"Conventions": "CF-1.7",
                "source": "MODIS MOD13Q1 NDVI estival (JJA) 2001-2024",
                "history": "18_tendance_ndvi.py : Sen + Mann-Kendall, n >= %d annees" % N_MIN}
    chemin = os.path.join(SORTIE, "trend_ndvi_JJA_%s.nc" % nom)
    ds.to_netcdf(chemin, encoding={nom: {"zlib": True, "complevel": 4,
                                         "_FillValue": np.float32(np.nan)}})
    print("  -> %s" % os.path.basename(chemin))


# ===========================================================================
# 1. Validation sur la temperature (comparaison avec ArcGIS)
# ===========================================================================
print("=== Validation : t2m JJA, pixel de Bruxelles ===")
e = xr.open_dataset(ERA5)["t2m"].transpose("time", "latitude", "longitude")
pe, pp, _ = mann_kendall(e.values.astype("float32"),
                         e["time"].dt.year.values.astype("float64"))
iy = int(np.abs(e["latitude"].values - 50.85).argmin())
ix = int(np.abs(e["longitude"].values - 4.35).argmin())
print("  pixel (%.2f N, %.2f E)" % (e["latitude"].values[iy], e["longitude"].values[ix]))
print("  pente   : %.4f degC / decennie" % (pe[iy, ix] * 10))
print("  p-value : %.2e" % pp[iy, ix])
print("  -> a comparer avec trend_t2m_JJA_pente_dec.tif et _pvalue.tif au meme endroit")

# ===========================================================================
# 2. Tendance du NDVI
# ===========================================================================
print("\n=== Tendance du NDVI estival ===")
ds = xr.open_dataset(ENTREE)
nd = ds["ndvi"].transpose("time", "latitude", "longitude")
annees = nd["time"].dt.year.values.astype("float64")
print("  %d annees (%d-%d), grille %d x %d" % (
    len(annees), annees.min(), annees.max(), nd.sizes["latitude"], nd.sizes["longitude"]))
print("  calcul en cours ...")

pente, pval, n = mann_kendall(nd.values.astype("float32"), annees)
pente_dec = pente * 10

# --- Diagnostic : d'ou viennent les valeurs aberrantes ? --------------------
terre = n > 0
print("\n=== Annees valides par pixel ===")
for seuil in (24, 20, 18, 12, 6, 3):
    print("  pixels avec >= %2d annees : %5.1f %%" % (seuil, 100.0 * (n[terre] >= seuil).mean()))

peu = terre & (n < N_MIN)
assez = n >= N_MIN
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    print("\n  |pente| max, pixels < %d annees : %.3f" % (N_MIN, np.nanmax(np.abs(pente_dec[peu])) if peu.any() else 0))
    print("  |pente| max, pixels >= %d annees : %.3f" % (N_MIN, np.nanmax(np.abs(pente_dec[assez]))))

# --- Masque et significativite ----------------------------------------------
pente_dec[~assez] = np.nan
pval[~assez] = np.nan
signif = np.where(pval < SEUIL, pente_dec, np.nan)

v = pente_dec[~np.isnan(pente_dec)]
s = signif[~np.isnan(signif)]
print("\n=== Resultat (pixels >= %d annees) ===" % N_MIN)
print("  pixels retenus    : %d" % len(v))
print("  pente 2e / 50e / 98e centile : %+.4f / %+.4f / %+.4f par decennie" % tuple(
    np.percentile(v, [2, 50, 98])))
print("  significatifs     : %.0f %%  (baisse %.0f %%, hausse %.0f %%)" % (
    100.0 * len(s) / len(v), 100.0 * (s < 0).sum() / len(v), 100.0 * (s > 0).sum() / len(v)))

for nom, la, lo in [("Foret de Soignes", 50.77, 4.43), ("Cultures, Hesbaye", 50.65, 5.10),
                    ("Port d'Anvers", 51.28, 4.30)]:
    iy = int(np.abs(nd["latitude"].values - la).argmin())
    ix = int(np.abs(nd["longitude"].values - lo).argmin())
    print("  %-18s pente %+.4f   p %.3f   n %d" % (nom, pente_dec[iy, ix], pval[iy, ix], n[iy, ix]))

# ===========================================================================
# 3. Ecriture
# ===========================================================================
print("\n=== Ecriture ===")
ecrire(pente_dec, nd, "pente_dec",
       {"units": "NDVI / decennie", "long_name": "Pente de Sen du NDVI estival"})
ecrire(pval, nd, "pvalue",
       {"units": "1", "long_name": "p-value du test de Mann-Kendall"})
ecrire(signif, nd, "pente_dec_signif",
       {"units": "NDVI / decennie", "long_name": "Pente de Sen, p < %.2f" % SEUIL})

print("\nTermine.")
print("Conversion en TIF, dans la fenetre Python d'ArcGIS :")
print('  import arcpy, os')
print('  d = r"%s"' % SORTIE)
print('  for s in ("pente_dec", "pvalue", "pente_dec_signif"):')
print('      arcpy.management.CopyRaster(os.path.join(d, "trend_ndvi_JJA_%s.nc" % s),')
print('                                  os.path.join(d, "trend_ndvi_JJA_%s.tif" % s))')