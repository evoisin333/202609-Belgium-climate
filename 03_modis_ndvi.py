"""
03_modis_ndvi.py
NDVI MODIS MOD13Q1 v061 (250 m, composites 16 jours, extraction AppEEARS) :
  - masque qualite : pixel_reliability 0 (bon) et 1 (marginal) conserves
  - valeur de remplissage -3000 retiree, facteur d'echelle 0,0001
  - chaque composite est rattache a une saison par sa date mediane
  - moyenne saisonniere, meme convention DJF que 02_agregation.py
Periode conservee : 2001-2024 (hivers 2000 et 2025 incomplets).
Sortie : 02_working/modis_ndvi_BE_saisonnier.nc

A lancer depuis la racine du projet, environnement geo :
    python scripts/03_modis_ndvi.py
"""
import rioxarray as rxr
import xarray as xr
import numpy as np
import pandas as pd
import glob, os, re

RAW_NDVI = os.path.join("01_raw", "MODIS", "NDVI")
RAW_REL  = os.path.join("01_raw", "MODIS", "pixel_reliability")
OUT      = "02_working"

os.makedirs(OUT, exist_ok=True)

# --- appariement des fichiers par date ---
def date_du_fichier(chemin):
    m = re.search(r"doy(\d{4})(\d{3})", os.path.basename(chemin))
    an, jour = int(m.group(1)), int(m.group(2))
    debut = pd.Timestamp(f"{an}-01-01") + pd.Timedelta(days=jour - 1)
    return debut + pd.Timedelta(days=8)      # date mediane du composite 16 j

ndvi_f = {date_du_fichier(f): f for f in glob.glob(os.path.join(RAW_NDVI, "*.tif"))}
rel_f  = {date_du_fichier(f): f for f in glob.glob(os.path.join(RAW_REL,  "*.tif"))}

dates = sorted(set(ndvi_f) & set(rel_f))
print(f"{len(dates)} composites apparies")
print(f"de {dates[0].date()} a {dates[-1].date()}\n")

# --- lecture, masquage, mise a l'echelle ---
couches = []
for i, d in enumerate(dates):
    n = rxr.open_rasterio(ndvi_f[d], masked=False).squeeze(drop=True)
    r = rxr.open_rasterio(rel_f[d],  masked=False).squeeze(drop=True)

    valide = (r == 0) | (r == 1)
    n = n.where(valide & (n != -3000))
    n = n * 0.0001

    couches.append(n.expand_dims(time=[d]))

    if (i + 1) % 50 == 0:
        print(f"  {i+1}/{len(dates)}")

ndvi = xr.concat(couches, dim="time").rename("ndvi")
ndvi.attrs = {"long_name": "NDVI MODIS MOD13Q1", "units": "-"}
print("\nempilement termine :", dict(ndvi.sizes))

# --- couverture valide par saison ---
mois = ndvi.time.dt.month
saison = xr.where(mois.isin([12,1,2]), "DJF",
         xr.where(mois.isin([3,4,5]),  "MAM",
         xr.where(mois.isin([6,7,8]),  "JJA", "SON")))
annee = xr.where(mois == 12, ndvi.time.dt.year + 1, ndvi.time.dt.year)
ndvi = ndvi.assign_coords(saison=saison, annee=annee)

print("\n--- couverture valide (% de pixels) ---")
for s in ["DJF", "MAM", "JJA", "SON"]:
    sub = ndvi.sel(time=ndvi.saison == s)
    pct = float(sub.notnull().mean()) * 100
    print(f"  {s} : {pct:5.1f} %   ({sub.sizes['time']} composites)")

# --- agregation saisonniere ---
seas = ndvi.groupby("saison").map(lambda g: g.groupby("annee").mean("time"))
seas = seas.sel(annee=slice(2001, 2024))

seas.to_netcdf(os.path.join(OUT, "modis_ndvi_BE_saisonnier.nc"),
               encoding={"ndvi": {"zlib": True, "complevel": 4}})
print(f"\n-> modis_ndvi_BE_saisonnier.nc  {dict(seas.sizes)}")

print("\n--- NDVI moyen par saison ---")
for s in ["DJF", "MAM", "JJA", "SON"]:
    print(f"  {s} : {float(seas.sel(saison=s).mean()):.3f}")