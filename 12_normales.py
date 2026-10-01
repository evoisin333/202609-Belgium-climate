# -*- coding: utf-8 -*-
"""
12_normales.py

Normales climatiques et cartes de changement, pour toutes les saisons
et tous les indicateurs.

Entrees  : les NetCDF saisonniers de 03_arcgis (ceux qui alimentent ArcGIS),
           qui ont tous une dimension 'time' datee au milieu de saison.
Sorties  : des GeoTIFF dans 04_normales, directement lisibles par ArcGIS Pro,
           plus un CSV recapitulatif.

Trois periodes :
  N1  normale OMM 1961-1990
  N2  normale OMM 1991-2020
  D   decennie recente 2015-2024 (illustrative)

Deux cartes de changement par indicateur et par saison :
  diff_N2_N1   changement entre normales (reference du rapport)
  diff_D_N1    decennie recente face a 1961-1990 (plus parlant en StoryMap)

Pour les variables de cumul ou de stock (precipitations, humidite du sol),
un changement relatif en pourcentage est produit en plus du changement absolu.

Un pixel n'est calcule que si au moins 80 % des annees de la periode sont
disponibles. Cela laisse passer les pixels ou 4 etes ont ete masques
(26 annees sur 30) tout en ecartant les zones trop lacunaires.

Le NDVI n'est pas traite ici : sa serie commence en 2001 et releve de
l'axe impact, avec anomalies et tendances (phase 4).

A lancer depuis la racine du projet, environnement geo :
    python 12_normales.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import rioxarray  # noqa: F401  (active l'accesseur .rio)
import xarray as xr

# --------------------------------------------------------------------------
# Parametres
# --------------------------------------------------------------------------
RACINE = Path.cwd()   # racine du projet
ENTREE = RACINE / "03_arcgis"
SORTIE = RACINE / "04_normales"
SORTIE.mkdir(exist_ok=True)

SAISONS = ["DJF", "MAM", "JJA", "SON"]

PERIODES = {
    "N1": (1961, 1990),
    "N2": (1991, 2020),
    "D": (2015, 2024),
}

COUVERTURE_MIN = 0.8        # part minimale d'annees disponibles par pixel

# Prefixes de fichiers a traiter (le NDVI est exclu volontairement)
PREFIXES = [
    "era5_BE",
    "eobs_fd_BE",
    "eobs_su_BE",
    "eobs_tr_BE",
    "eobs_txx_BE",
    "eobs_prcptot_BE",
    "eobs_su30_BE",
]

# Variables pour lesquelles un changement relatif (%) a du sens
RELATIF = {"tp", "swvl1", "swvl2", "prcptot"}

CRS = "EPSG:4326"
NODATA = -9999.0


# --------------------------------------------------------------------------
# Outils
# --------------------------------------------------------------------------
def nom_court(variable):
    """suETCCDI -> su, prcptotETCCDI -> prcptot, t2m -> t2m."""
    return variable.replace("ETCCDI", "")


def moyenne_periode(da, an_min, an_max):
    """Moyenne sur une periode, avec exigence de couverture temporelle."""
    annees = da["time"].dt.year
    sel = da.sel(time=(annees >= an_min) & (annees <= an_max))
    n_attendu = int(sel.sizes["time"])
    if n_attendu == 0:
        return None, 0
    n_dispo = sel.notnull().sum("time")
    moy = sel.mean("time", skipna=True).where(n_dispo >= COUVERTURE_MIN * n_attendu)
    return moy, n_attendu


def ecrire_tif(da, chemin, nom_long, unite):
    """Ecrit un GeoTIFF georeference, latitudes du nord au sud.

    Le NoData est fixe explicitement a NODATA plutot que laisse a NaN :
    ArcGIS Pro reconnait alors la valeur sans ambiguite, y compris dans
    les statistiques zonales.
    """
    da = da.sortby("latitude", ascending=False)
    da = da.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude")
    da = da.rio.write_crs(CRS)
    da = da.astype("float32").fillna(NODATA)
    da = da.rio.write_nodata(NODATA)
    da.attrs = {"long_name": nom_long, "units": unite}
    da.rio.to_raster(chemin, compress="LZW")


# --------------------------------------------------------------------------
# Traitement
# --------------------------------------------------------------------------
recap = []
manquants = []

for prefixe in PREFIXES:
    for saison in SAISONS:
        fichier = ENTREE / f"{prefixe}_{saison}.nc"
        if not fichier.exists():
            manquants.append(fichier.name)
            continue

        print(f"\n=== {fichier.name} ===")
        with xr.open_dataset(fichier) as src:
            ds = src.load()

        for variable in ds.data_vars:
            if variable == "crs" or "time" not in ds[variable].dims:
                continue

            da = ds[variable]
            court = nom_court(variable)
            unite = da.attrs.get("units", "")
            nom_long = da.attrs.get("long_name", court)
            base = f"{court}_{saison}"

            # --- les trois normales ---
            normales = {}
            for cle, (a, b) in PERIODES.items():
                moy, n = moyenne_periode(da, a, b)
                if moy is None:
                    print(f"  {court:9s} {cle} : aucune annee dans {a}-{b}, ignore")
                    continue
                normales[cle] = moy
                ecrire_tif(moy, SORTIE / f"{base}_{cle}.tif",
                           f"{nom_long} - normale {a}-{b}", unite)
                recap.append({
                    "indicateur": court, "saison": saison, "sortie": f"{base}_{cle}",
                    "type": "normale", "periode": f"{a}-{b}", "annees": n,
                    "unite": unite, "moyenne_emprise": float(moy.mean()),
                })
                print(f"  {court:9s} {cle} : {float(moy.mean()):8.2f} {unite} "
                      f"({n} annees)")

            # --- les changements ---
            paires = [("diff_N2_N1", "N2", "N1"), ("diff_D_N1", "D", "N1")]
            for etiquette, fin, debut in paires:
                if fin not in normales or debut not in normales:
                    continue

                delta = normales[fin] - normales[debut]
                ecrire_tif(delta, SORTIE / f"{base}_{etiquette}.tif",
                           f"{nom_long} - changement {debut} vers {fin}", unite)
                recap.append({
                    "indicateur": court, "saison": saison,
                    "sortie": f"{base}_{etiquette}", "type": "changement absolu",
                    "periode": f"{debut} -> {fin}", "annees": "",
                    "unite": unite, "moyenne_emprise": float(delta.mean()),
                })
                print(f"  {court:9s} {etiquette:11s} : {float(delta.mean()):+8.2f} {unite}")

                if court in RELATIF:
                    # Seuls les pixels dont le denominateur est quasi nul sont
                    # ecartes (sinon le pourcentage explose) : 1 % de la
                    # mediane suffit, sans sacrifier de pixels valides.
                    ref = normales[debut]
                    seuil = 0.01 * float(abs(ref).median())
                    pct = (delta / ref.where(abs(ref) > seuil)) * 100.0
                    ecrire_tif(pct, SORTIE / f"{base}_{etiquette}_pct.tif",
                               f"{nom_long} - changement relatif {debut} vers {fin}", "%")
                    recap.append({
                        "indicateur": court, "saison": saison,
                        "sortie": f"{base}_{etiquette}_pct",
                        "type": "changement relatif", "periode": f"{debut} -> {fin}",
                        "annees": "", "unite": "%",
                        "moyenne_emprise": float(pct.mean()),
                    })
                    print(f"  {court:9s} {etiquette:11s} : {float(pct.mean()):+8.1f} %")

# --------------------------------------------------------------------------
# Recapitulatif
# --------------------------------------------------------------------------
if manquants:
    print("\n!! Fichiers absents, ignores :")
    for m in manquants:
        print(f"   {m}")
    print("   (normal si tu n'as exporte que les saisons pertinentes "
          "pour certains indices)")

df = pd.DataFrame(recap)
csv = SORTIE / "recapitulatif_normales.csv"
df.to_csv(csv, index=False, sep=";", decimal=",", encoding="utf-8-sig")

print(f"\n{len(df)} rasters ecrits dans {SORTIE}")
print(f"Recapitulatif : {csv}")
print("\nATTENTION : la colonne 'moyenne_emprise' porte sur le rectangle "
      "complet, France et mer comprises. Elle sert au controle de plausibilite, "
      "pas au rapport. Les chiffres du rapport viendront des statistiques "
      "zonales par province (phase 3).")

# --------------------------------------------------------------------------
# Controles de plausibilite
# --------------------------------------------------------------------------
print("\n=== Controles ===")
verif = df[(df["type"] == "changement absolu") & (df["periode"] == "N1 -> N2")]
for _, r in verif.iterrows():
    print(f"  {r['indicateur']:9s} {r['saison']} : "
          f"{r['moyenne_emprise']:+7.2f} {r['unite']}")
print("\n  Attendu : t2m entre +0,8 et +1,8 degC selon la saison ; "
      "fd en baisse ; su, tr, txx et su30 en hausse.")