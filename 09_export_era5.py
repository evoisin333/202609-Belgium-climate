"""
09_export_era5.py
Restructure era5_BE_saisonnier.nc (dimensions saison x annee) en 4 fichiers
NetCDF, un par saison, avec une vraie dimension temporelle lisible par ArcGIS Pro.
Dates au milieu de la saison : DJF -> 15 janvier, MAM -> 15 avril,
JJA -> 15 juillet, SON -> 15 octobre. Sorties dans 03_arcgis.

A lancer depuis la racine du projet, environnement geo :
    python scripts/09_export_era5.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# --- Chemins ---------------------------------------------------------------
RACINE = Path.cwd()   # racine du projet
ENTREE = RACINE / "02_working" / "era5_BE_saisonnier.nc"
SORTIE = RACINE / "03_arcgis"
SORTIE.mkdir(exist_ok=True)

# --- Parametres --------------------------------------------------------------
MOIS_MILIEU = {"DJF": 1, "MAM": 4, "JJA": 7, "SON": 10}

UNITES = {
    "t2m":   ("degC",    "Temperature moyenne saisonniere a 2 m"),
    "tp":    ("mm",      "Precipitations cumulees saisonnieres"),
    "swvl1": ("m3 m-3",  "Humidite volumique du sol 0-7 cm"),
    "swvl2": ("m3 m-3",  "Humidite volumique du sol 7-28 cm"),
}

CRS_ATTRS = {
    "grid_mapping_name": "latitude_longitude",
    "longitude_of_prime_meridian": 0.0,
    "semi_major_axis": 6378137.0,
    "inverse_flattening": 298.257223563,
    "epsg_code": "EPSG:4326",
}

# --- Traitement ----------------------------------------------------------------
ds = xr.open_dataset(ENTREE)
ds = ds.drop_vars("number", errors="ignore")

for saison, mois in MOIS_MILIEU.items():
    sub = ds.sel(saison=saison).drop_vars("saison")

    # annee -> time (dates au milieu de la saison)
    dates = pd.to_datetime(
        pd.DataFrame({"year": sub["annee"].values, "month": mois, "day": 15})
    )
    sub = (
        sub.assign_coords(time=("annee", dates))
        .swap_dims({"annee": "time"})
        .drop_vars("annee")
        .transpose("time", "latitude", "longitude")
    )

    # Saisons incompletes : tp est agrege avec min_count=3 dans 02_agregation.py,
    # donc un mois manquant rend tp entierement NaN (t2m, lui, reste calcule)
    vide = sub["tp"].isnull().all(dim=("latitude", "longitude"))
    if bool(vide.any()):
        annees_vides = sub["time"].dt.year.values[vide.values]
        print(f"[{saison}] annees entierement vides retirees : {list(annees_vides)}")
        sub = sub.sel(time=~vide)

    # Variables : float32, attributs propres, lien vers le CRS
    for v in UNITES:
        da = sub[v].astype("float32")
        da.attrs = {
            "units": UNITES[v][0],
            "long_name": UNITES[v][1],
            "grid_mapping": "crs",
        }
        da.encoding = {"zlib": True, "complevel": 4, "_FillValue": np.float32(np.nan)}
        sub[v] = da

    # Coordonnees
    sub["latitude"].attrs = {"standard_name": "latitude", "units": "degrees_north", "axis": "Y"}
    sub["longitude"].attrs = {"standard_name": "longitude", "units": "degrees_east", "axis": "X"}
    sub["time"].attrs = {"standard_name": "time", "axis": "T"}
    for c in ("latitude", "longitude", "time"):
        sub[c].encoding = {}
    sub["time"].encoding = {"units": "days since 1950-01-01", "calendar": "standard"}

    # Variable CRS
    sub["crs"] = xr.DataArray(np.int32(0), attrs=CRS_ATTRS)

    # Attributs globaux
    sub.attrs = {
        "Conventions": "CF-1.7",
        "title": f"ERA5-Land saisonnier Belgique - {saison}",
        "source": "ERA5-Land monthly means (Copernicus CDS)",
        "history": f"01_decoupe.py, 02_agregation.py, 09_export_era5.py (saison {saison})",
    }

    fichier = SORTIE / f"era5_BE_{saison}.nc"
    sub.to_netcdf(fichier)
    print(f"[{saison}] ecrit : {fichier.name}  ({sub.sizes['time']} annees)")

ds.close()

# --- Verification : relecture depuis le disque ---------------------------------------
print("\nControle : pixel Bruxelles (50.85 N, 4.35 E), ete 2003")
test = xr.open_dataset(SORTIE / "era5_BE_JJA.nc")
px = test.sel(latitude=50.85, longitude=4.35, method="nearest").sel(time="2003-07-15")
for v in UNITES:
    print(f"  {v:6s} = {float(px[v].values):.3f} {UNITES[v][0]}")
print(f"  dates : {str(test.time.values[0])[:10]} -> {str(test.time.values[-1])[:10]}")
test.close()