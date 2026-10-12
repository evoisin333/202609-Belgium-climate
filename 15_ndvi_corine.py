# -*- coding: utf-8 -*-
"""
15_ndvi_corine.py
A lancer depuis la racine du projet, environnement geo :
    python 15_ndvi_corine.py

NDVI par province ET par classe d'occupation du sol (CORINE Land Cover 2018,
44 classes regroupees en 8), sur la grille du NDVI MODIS (~230 m), en ne
gardant que les pixels MODIS "purs" :
  1. lecture de CORINE (100 m, ETRS89-LAEA) sur l'emprise de la Belgique ;
  2. regroupement en 8 classes ; l'eau libre (codes 511-523) est une classe a
     part, exclue de l'analyse NDVI ;
  3. part de chaque classe dans chaque pixel MODIS (reprojection par moyenne
     d'un masque 0/1) ;
  4. un pixel MODIS n'est garde que si une classe couvre au moins 80 % de sa
     surface (Kulesza et Hoscilo 2024) : un pixel de 5 ha a cheval sur une
     ville et des jardins, ou sur un champ et une haie, n'est plus attribue a
     une seule classe ;
  5. provinces rasterisees depuis le GeoJSON NUTS ;
  6. NDVI moyen par province x classe, par saison et par annee.

Sorties :
  04_stats/zs_ndvi_<saison>.csv          meme format qu'avant (lu par 16)
  04_stats/corine_8classes_libelles.csv
  01_raw/corine_8classes_grille_modis.tif   classes pures (lu par 20)

Entrees :
  raster CORINE europeen U2018_CLC2018*.tif (dans le dossier parent du projet
  ou ses sous-dossiers), 01_raw/NUTS_RG_01M_2024_4326_LEVL_2.geojson,
  02_working/modis_ndvi_BE_saisonnier.nc (03_modis_ndvi.py)
"""

import glob
import json
import os

import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform_bounds
from rasterio.windows import from_bounds

RACINE = os.getcwd()   # racine du projet
NDVI_NC = os.path.join(RACINE, "02_working", "modis_ndvi_BE_saisonnier.nc")
NUTS = os.path.join(RACINE, "01_raw", "NUTS_RG_01M_2024_4326_LEVL_2.geojson")
CLASSES_TIF = os.path.join(RACINE, "01_raw", "corine_8classes_grille_modis.tif")
STATS = os.path.join(RACINE, "04_stats")
os.makedirs(STATS, exist_ok=True)

PURETE_MIN = 0.80        # part minimale de la classe dominante dans un pixel MODIS
SAISONS = {"DJF": 1, "MAM": 4, "JJA": 7, "SON": 10}   # mois des dates (09-11)

LIBELLES = {
    1: "Artificialise", 2: "Cultures", 3: "Prairies",
    4: "Foret feuillus", 5: "Foret resineux", 6: "Foret mixte",
    7: "Milieux semi-naturels", 8: "Zones humides",
}
EAU = 9   # eau libre : calculee pour la purete des pixels, puis exclue

# Codes CORINE a 3 chiffres -> groupe (l'ordre donne les codes de grille 1-44)
CODES_3 = {}
for c in [111, 112, 121, 122, 123, 124, 131, 132, 133, 141, 142]:
    CODES_3[c] = 1
for c in [211, 212, 213, 221, 222, 223, 241, 242, 243, 244]:
    CODES_3[c] = 2
CODES_3[231] = 3
CODES_3[311] = 4
CODES_3[312] = 5
CODES_3[313] = 6
for c in [321, 322, 323, 324, 331, 332, 333, 334, 335]:
    CODES_3[c] = 7
for c in [411, 412, 421, 422, 423]:
    CODES_3[c] = 8
for c in [511, 512, 521, 522, 523]:
    CODES_3[c] = EAU
CODES_GRILLE = dict((i + 1, CODES_3[c]) for i, c in enumerate(sorted(CODES_3)))

# ===========================================================================
# 1. Grille MODIS (celle du NDVI saisonnier)
# ===========================================================================
ds = xr.open_dataset(NDVI_NC)
ndvi = ds["ndvi"].transpose("saison", "annee", "y", "x")
x, y = ds["x"].values, ds["y"].values
res = float(abs(x[1] - x[0]))
transform = from_origin(x.min() - res / 2, y.max() + res / 2, res, res)
forme = (len(y), len(x))
if y[0] < y[-1]:
    raise SystemExit("Latitudes croissantes inattendues dans %s" % NDVI_NC)
print("grille MODIS : %d x %d pixels, %.6f deg" % (forme[0], forme[1], res))

# ===========================================================================
# 2. CORINE sur l'emprise de la Belgique, regroupe en classes
# ===========================================================================
motif = os.path.join(os.path.dirname(RACINE), "**", "U2018_CLC2018*.tif")
trouves = [f for f in glob.glob(motif, recursive=True) if f.lower().endswith(".tif")]
if not trouves:
    raise SystemExit("Raster CORINE europeen introuvable (%s)" % motif)
clc_tif = trouves[0]
print("CORINE : %s" % clc_tif)

with rasterio.open(clc_tif) as src:
    bornes = transform_bounds("EPSG:4326", src.crs, x.min() - res, y.min() - res,
                              x.max() + res, y.max() + res)
    fenetre = from_bounds(*bornes, transform=src.transform).round_offsets().round_lengths()
    clc = src.read(1, window=fenetre)
    clc_transform = src.window_transform(fenetre)
    clc_crs = src.crs
table = CODES_3 if np.isin(clc, list(CODES_3)).any() else CODES_GRILLE
groupe = np.zeros(clc.shape, dtype="uint8")
for code, g in table.items():
    groupe[clc == code] = g
print("  %d x %d cellules de 100 m, codes %s" % (
    clc.shape[0], clc.shape[1], "a 3 chiffres" if table is CODES_3 else "de grille 1-44"))

# ===========================================================================
# 3. Part de chaque classe dans chaque pixel MODIS
# ===========================================================================
def part_sur_grille(masque):
    """Moyenne d'un masque 0/1 de 100 m sur chaque pixel MODIS."""
    dst = np.full(forme, np.nan, dtype="float32")
    reproject(masque.astype("float32"), dst, src_transform=clc_transform, src_crs=clc_crs,
              dst_transform=transform, dst_crs="EPSG:4326",
              resampling=Resampling.average, src_nodata=np.nan, dst_nodata=np.nan)
    return dst


parts = np.stack([part_sur_grille(np.where(groupe > 0, (groupe == k) * 1.0, np.nan))
                  for k in range(1, EAU + 1)])        # (9, y, x)
parts = np.nan_to_num(parts, nan=0.0)
dominante = parts.argmax(axis=0) + 1
purete = parts.max(axis=0)
classes = np.where((purete >= PURETE_MIN) & (dominante != EAU), dominante, 0).astype("uint8")

# ===========================================================================
# 4. Provinces sur la grille MODIS ; tout ce qui est hors de Belgique est
#    retire (le rectangle lu dans CORINE deborde sur les pays voisins)
# ===========================================================================
with open(NUTS, encoding="utf-8") as f:
    nuts = [e for e in json.load(f)["features"] if e["properties"]["CNTR_CODE"] == "BE"]
codes_nuts = sorted(e["properties"]["NUTS_ID"] for e in nuts)
if len(codes_nuts) != 11:
    raise SystemExit("%d regions NUTS 2 belges (11 attendues)" % len(codes_nuts))
indice = {c: i + 1 for i, c in enumerate(codes_nuts)}
provinces = rasterize([(e["geometry"], indice[e["properties"]["NUTS_ID"]]) for e in nuts],
                      out_shape=forme, transform=transform, fill=0, dtype="uint8")
belgique = provinces > 0
classes[~belgique] = 0

print("\n=== Purete des pixels (seuil %.0f %%) ===" % (PURETE_MIN * 100))
print("  %-22s %10s %10s %7s" % ("classe", "dominante", "pure", "gardes"))
for k in range(1, EAU + 1):
    tous = int(np.sum(belgique & (dominante == k) & (parts.sum(axis=0) > 0)))
    purs = int(np.sum(belgique & (dominante == k) & (purete >= PURETE_MIN)))
    nom = LIBELLES.get(k, "Eau libre (exclue)")
    print("  %-22s %10d %10d %6.0f %%" % (nom, tous, purs, 100.0 * purs / tous if tous else 0))

with rasterio.open(CLASSES_TIF, "w", driver="GTiff", height=forme[0], width=forme[1],
                   count=1, dtype="uint8", crs="EPSG:4326", transform=transform,
                   nodata=0, compress="deflate") as dst:
    dst.write(classes, 1)
print("  -> %s" % os.path.relpath(CLASSES_TIF, RACINE))

# ===========================================================================
# 5. NDVI moyen par province x classe, saison, annee
# ===========================================================================
zones = np.where((provinces > 0) & (classes > 0), provinces * 10 + classes, 0).ravel()
print("\n=== NDVI par province et classe d'occupation du sol ===")
annees = ndvi["annee"].values
for saison, mois in SAISONS.items():
    lignes = []
    for a in annees:
        v = ndvi.sel(saison=saison, annee=a).values.ravel()
        ok = (zones > 0) & ~np.isnan(v)
        somme = np.bincount(zones[ok], weights=v[ok], minlength=120)
        nb = np.bincount(zones[ok], minlength=120)
        for z in np.nonzero(nb)[0]:
            p, c = divmod(int(z), 10)
            lignes.append({"Value": int(z), "COUNT": int(nb[z]),
                           "StdTime": "%d-%02d-15" % (a, mois),
                           "MEAN": somme[z] / nb[z], "NUTS_ID": codes_nuts[p - 1],
                           "CLC_CODE": c, "CLASSE": LIBELLES[c]})
    t = pd.DataFrame(lignes)
    t.to_csv(os.path.join(STATS, "zs_ndvi_%s.csv" % saison), index=False, sep=";", decimal=",")
    print("  ok : zs_ndvi_%s  %5d lignes  (%d zones)" % (saison, len(t), t["Value"].nunique()))

with open(os.path.join(STATS, "corine_8classes_libelles.csv"), "w") as f:
    f.write("code,libelle\n")
    for k in sorted(LIBELLES):
        f.write("%d,%s\n" % (k, LIBELLES[k]))

# ===========================================================================
# 6. Controle : NDVI ete 2003 par classe
# ===========================================================================
print("\n=== Controle : NDVI ete 2003 par classe (moyenne ponderee des provinces) ===")
t = pd.read_csv(os.path.join(STATS, "zs_ndvi_JJA.csv"), sep=";", decimal=",")
t = t[t["StdTime"].str.startswith("2003")]
for classe, g in t.groupby("CLASSE"):
    print("  %-22s %.3f  (%d provinces, %d pixels)" % (
        classe, np.average(g["MEAN"], weights=g["COUNT"]), len(g), g["COUNT"].sum()))
print("\nTermine.")
