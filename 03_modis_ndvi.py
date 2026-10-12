"""
03_modis_ndvi.py
NDVI MODIS MOD13Q1 v061 (250 m, composites 16 jours, extraction AppEEARS),
moyennes saisonnieres par pixel :
  - qualite : pixel_reliability 0 garde ; 1 (marginal) garde seulement si
    VI_Quality ne signale ni nuage, ni ombre, ni neige, et si l'indice
    d'utilite n'est pas "trop bas" (guide MOD13, voir outils_modis.py) ;
  - valeur de remplissage -3000 retiree, facteur d'echelle 0,0001 ;
  - chaque pixel est rattache a la saison de son jour d'acquisition reel
    (couche composite_day_of_the_year), pas a celle du milieu du composite ;
  - une moyenne saisonniere n'est calculee que si le pixel a au moins
    3 observations valides dans la saison (sur environ 6 composites) ;
  - printemps (MAM), ete (JJA) et automne (SON) seulement : l'hiver est
    exclu (neige retiree par le masque, donc moyenne biaisee vers les jours
    sans neige, peu d'observations, soleil bas).
Periode : 2001-2024.
Sortie : 02_working/modis_ndvi_BE_saisonnier.nc (dimensions saison, annee, y, x)

Les moyennes sont cumulees composite par composite : rien n'est empile en
memoire.

A lancer depuis la racine du projet, environnement geo :
    python 03_modis_ndvi.py
Le dossier 01_raw/MODIS doit contenir NDVI, pixel_reliability, VI_Quality
et composite_day_of_the_year.
"""
import os

import numpy as np
import rioxarray as rxr
import xarray as xr

from outils_modis import ECHELLE, date_des_pixels, fichiers_apparies, lire, masque_qualite

MODIS = os.path.join("01_raw", "MODIS")
OUT = "02_working"
os.makedirs(OUT, exist_ok=True)

ANNEES = list(range(2001, 2025))
SAISONS = {"MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}
OBS_MIN = 3   # observations valides minimum par pixel et par saison

# --- appariement des fichiers par date de debut ---
composites = fichiers_apparies(MODIS, "NDVI", (ANNEES[0], ANNEES[-1]))
print(f"{len(composites)} composites apparies (NDVI + 3 couches qualite)")
print(f"de {composites[0][:2]} a {composites[-1][:2]} (annee, jour)\n")

# --- controle du systeme de coordonnees (sortie AppEEARS) ---
# 10_export_ndvi_su30.py traite x/y comme longitude/latitude : la requete
# AppEEARS doit donc sortir en Geographic (EPSG:4326), pas en sinusoidal.
with rxr.open_rasterio(composites[0][2]["NDVI"]) as g:
    crs, forme, transform = g.rio.crs, g.rio.shape, g.rio.transform()
    x, y = g["x"].values, g["y"].values
if crs is None or crs.to_epsg() != 4326:
    raise SystemExit(f"NDVI MODIS en {crs} : EPSG:4326 attendu "
                     "(choisir Geographic comme projection dans AppEEARS)")
for couche, chemin in composites[0][2].items():
    with rxr.open_rasterio(chemin) as r:
        if r.rio.crs != crs or r.rio.shape != forme or r.rio.transform() != transform:
            raise SystemExit(f"{couche} n'est pas sur la meme grille que le NDVI")
print(f"CRS MODIS : EPSG:4326, 4 couches sur la meme grille ({forme[0]} x {forme[1]} px)\n")

# --- cumul des sommes et des nombres d'observations par saison et annee ---
somme = np.zeros((len(SAISONS), len(ANNEES)) + forme, dtype="float32")
nombre = np.zeros((len(SAISONS), len(ANNEES)) + forme, dtype="uint8")
details = {}            # pixels marginaux retires, par critere
deplaces = valides = 0  # pixels dont la saison reelle differe de celle du milieu du composite

for i, (an, jour, f) in enumerate(composites):
    n = lire(f["NDVI"])
    ok = masque_qualite(n, lire(f["pixel_reliability"]), lire(f["VI_Quality"]), details)
    annee, mois = date_des_pixels(an, jour, lire(f["composite_day_of_the_year"]))
    ok &= mois > 0

    # saison du milieu du composite (ancienne methode), pour le bilan
    _, mois_milieu = date_des_pixels(an, jour, np.array([jour + 8]))
    for k, (s, m) in enumerate(SAISONS.items()):
        dans_s = ok & np.isin(mois, m)
        if not dans_s.any():
            continue
        if mois_milieu[0] not in m:
            deplaces += int(dans_s.sum())
        valides += int(dans_s.sum())
        for a in np.unique(annee[dans_s]):
            if a not in ANNEES:
                continue
            sel = dans_s & (annee == a)
            somme[k, ANNEES.index(a)][sel] += n[sel] * ECHELLE
            nombre[k, ANNEES.index(a)][sel] += 1

    if (i + 1) % 50 == 0:
        print(f"  {i+1}/{len(composites)}")

# --- bilan du masque qualite ---
print("\n--- pixels marginaux (fiabilite 1) ---")
m = max(details["marginaux"], 1)
print(f"  {details['marginaux']} pixels marginaux, {100 * details['marginaux gardes'] / m:.1f} % gardes")
for nom in ("utilite trop basse", "nuage adjacent", "nuages melanges", "neige possible", "ombre possible"):
    print(f"  retires pour {nom:20s}: {100 * details[nom] / m:5.1f} %  (criteres non exclusifs)")
print(f"\n  observations valides dont la saison reelle differe de celle du milieu "
      f"du composite : {100 * deplaces / max(valides, 1):.1f} %")

# --- moyennes saisonnieres (minimum OBS_MIN observations) ---
with np.errstate(invalid="ignore", divide="ignore"):
    moyenne = np.where(nombre >= OBS_MIN, somme / nombre, np.nan).astype("float32")
del somme

terre = (nombre > 0).any(axis=(0, 1))   # pixels jamais observes : hors emprise
print(f"\n--- couverture : pixels-annees avec >= {OBS_MIN} observations ---")
for k, s in enumerate(SAISONS):
    obs = nombre[k][:, terre]
    print(f"  {s} : {100 * np.mean(obs >= OBS_MIN):5.1f} %   "
          f"(observations par pixel et par saison : mediane {np.median(obs):.0f})")

seas = xr.DataArray(moyenne, name="ndvi", dims=("saison", "annee", "y", "x"),
                    coords={"saison": list(SAISONS), "annee": ANNEES,
                            "y": y, "x": x},
                    attrs={"long_name": "NDVI MODIS MOD13Q1, moyenne saisonniere", "units": "-"})
seas = seas.rio.write_crs("EPSG:4326")
seas.to_netcdf(os.path.join(OUT, "modis_ndvi_BE_saisonnier.nc"),
               encoding={"ndvi": {"zlib": True, "complevel": 4}})
print(f"\n-> modis_ndvi_BE_saisonnier.nc  {dict(seas.sizes)}")

print("\n--- NDVI moyen par saison ---")
for s in SAISONS:
    print(f"  {s} : {float(seas.sel(saison=s).mean()):.3f}")
