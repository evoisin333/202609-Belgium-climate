"""
10_export_ndvi_su30.py
Meme logique que 09_export_era5.py, appliquee a :
  - modis_ndvi_BE_saisonnier.nc    (03_modis_ndvi.py ; dimensions saison, y, x, annee)
  - eobs_su30_BE_saisonnier_qc.nc  (02_agregation.py, puis masques de
                                    07_controle_spatial_qc.py et
                                    08_masque_artefact_lux.py)
Sortie (03_arcgis) : un NetCDF par saison, avec une vraie dimension 'time'
pour ArcGIS Pro.

A lancer depuis la racine du projet, environnement geo :
    python scripts/10_export_ndvi_su30.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

RACINE = Path.cwd()   # racine du projet
ENTREE = RACINE / "02_working"
SORTIE = RACINE / "03_arcgis"
SORTIE.mkdir(exist_ok=True)

MOIS_MILIEU = {"DJF": 1, "MAM": 4, "JJA": 7, "SON": 10}

CRS_ATTRS = {
    "grid_mapping_name": "latitude_longitude",
    "longitude_of_prime_meridian": 0.0,
    "semi_major_axis": 6378137.0,
    "inverse_flattening": 298.257223563,
    "epsg_code": "EPSG:4326",
}


def exporter(fichier_entree, prefixe, variables, source, historique):
    ds = xr.open_dataset(fichier_entree)

    # Noms de coordonnees homogenes (MODIS utilise x / y)
    noms = {"x": "longitude", "y": "latitude", "lon": "longitude", "lat": "latitude"}
    ds = ds.rename({k: v for k, v in noms.items() if k in ds.dims})

    # Retire les coordonnees annexes (spatial_ref, number...) : un CRS propre est ajoute plus bas
    ds = ds.reset_coords(drop=True)

    for saison, mois in MOIS_MILIEU.items():
        sub = ds.sel(saison=saison, drop=True)

        dates = pd.to_datetime(pd.DataFrame({"year": sub["annee"].values, "month": mois, "day": 15}))
        sub = sub.assign_coords(time=("annee", dates))
        sub = sub.swap_dims({"annee": "time"})
        sub = sub.drop_vars("annee")
        sub = sub.transpose("time", "latitude", "longitude")

        # Annees entierement vides
        v0 = list(variables)[0]
        vide = sub[v0].isnull().all(dim=("latitude", "longitude")).values
        if vide.any():
            print(f"  [{saison}] annees vides retirees : {list(sub['time'].dt.year.values[vide])}")
            sub = sub.isel(time=~vide)

        for v, (unite, nom_long) in variables.items():
            da = sub[v].astype("float32")
            da.attrs = {"units": unite, "long_name": nom_long, "grid_mapping": "crs"}
            da.encoding = {"zlib": True, "complevel": 4, "_FillValue": np.float32(np.nan)}
            sub[v] = da

        sub["latitude"].attrs = {"standard_name": "latitude", "units": "degrees_north", "axis": "Y"}
        sub["longitude"].attrs = {"standard_name": "longitude", "units": "degrees_east", "axis": "X"}
        sub["time"].attrs = {"standard_name": "time", "axis": "T"}
        for c in ("latitude", "longitude"):
            sub[c].encoding = {}
        sub["time"].encoding = {"units": "days since 1950-01-01", "calendar": "standard"}

        sub["crs"] = xr.DataArray(np.int32(0), attrs=CRS_ATTRS)
        sub.attrs = {
            "Conventions": "CF-1.7",
            "title": f"{prefixe} saisonnier Belgique - {saison}",
            "source": source,
            "history": f"{historique}, 10_export_ndvi_su30.py (saison {saison})",
        }

        sortie = SORTIE / f"{prefixe}_BE_{saison}.nc"
        sub.to_netcdf(sortie)
        print(f"  [{saison}] ecrit : {sortie.name}  ({sub.sizes['time']} annees)")

    ds.close()


# ---------------------------------------------------------------------------
print("=== MODIS NDVI ===")
exporter(
    ENTREE / "modis_ndvi_BE_saisonnier.nc",
    "modis_ndvi",
    {"ndvi": ("1", "NDVI moyen saisonnier (MODIS)")},
    "MODIS NDVI via AppEEARS (NASA LP DAAC)",
    "03_modis_ndvi.py",
)

print("\n=== E-OBS jours >= 30 degC ===")
exporter(
    ENTREE / "eobs_su30_BE_saisonnier_qc.nc",
    "eobs_su30",
    {"su30": ("jours", "Jours avec Tmax >= 30 degC")},
    "E-OBS v31.0e, Tmax journaliere (tx_ens_mean 0.1deg)",
    "02_agregation.py, 07_controle_spatial_qc.py, 08_masque_artefact_lux.py",
)

# ---------------------------------------------------------------------------
print("\n=== Controles (ete 2003) ===")
n = xr.open_dataset(SORTIE / "modis_ndvi_BE_JJA.nc")
print(f"  NDVI JJA min / max : {float(n.ndvi.min()):.3f} / {float(n.ndvi.max()):.3f}  (attendu entre -1 et 1)")
foret = n.ndvi.sel(latitude=50.77, longitude=4.43, method="nearest").sel(time="2003-07-15")
print(f"  NDVI Foret de Soignes (50.77 N, 4.43 E) : {float(foret):.3f}  (attendu ~0.8)")
print(f"  dates : {str(n.time.values[0])[:10]} -> {str(n.time.values[-1])[:10]}")
n.close()

s = xr.open_dataset(SORTIE / "eobs_su30_BE_JJA.nc")
bxl = s.su30.sel(latitude=50.85, longitude=4.35, method="nearest").sel(time="2003-07-15")
print(f"  Jours >= 30 degC Bruxelles : {float(bxl):.0f} jours")
print(f"  dates : {str(s.time.values[0])[:10]} -> {str(s.time.values[-1])[:10]}")
s.close()
