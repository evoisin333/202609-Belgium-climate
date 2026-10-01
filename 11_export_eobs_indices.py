"""
11_export_eobs_indices.py
Separe les indices saisonniers E-OBS pre-calcules (fd, su, tr, txx, prcptot)
en un NetCDF par saison, periode 1961-2024, dates harmonisees avec ERA5
(DJF -> 15 janvier, MAM -> 15 avril, JJA -> 15 juillet, SON -> 15 octobre).
SU et TXx sont lus dans leur version _qc (07_controle_spatial_qc.py puis
08_masque_artefact_lux.py). Affiche aussi une comparaison entre les versions
d'origine et _qc. Sorties dans 03_arcgis.

A lancer depuis la racine du projet, environnement geo :
    python scripts/11_export_eobs_indices.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

RACINE = Path.cwd()   # racine du projet
ENTREE = RACINE / "02_working"
SORTIE = RACINE / "03_arcgis"
SORTIE.mkdir(exist_ok=True)

FICHIERS = {
    "fd":      "fd_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE.nc",
    "su":      "su_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE_qc.nc",   # version _qc
    "tr":      "tr_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE.nc",
    "txx":     "txx_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE_qc.nc",  # version _qc
    "prcptot": "prcptot_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE.nc",
}

VERSIONS_QC = {
    "su":  "su_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE_qc.nc",
    "txx": "txx_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE_qc.nc",
}

# Versions d'origine, pour la comparaison finale avec les versions _qc
NORMALES = {
    "su":  "su_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE.nc",
    "txx": "txx_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE.nc",
}

SAISON_DU_MOIS = {1: "DJF", 4: "MAM", 7: "JJA", 10: "SON"}
PERIODE = ("1961-01-01", "2024-12-31")  # DJF 1950 et DJF 2025 incomplets exclus

CRS_ATTRS = {
    "grid_mapping_name": "latitude_longitude",
    "longitude_of_prime_meridian": 0.0,
    "semi_major_axis": 6378137.0,
    "inverse_flattening": 298.257223563,
    "epsg_code": "EPSG:4326",
}


def ouvrir(fichier):
    """Ouvre un fichier E-OBS et renvoie sa variable principale (avec dimension time)."""
    ds = xr.open_dataset(fichier)
    noms = {"lon": "longitude", "lat": "latitude", "x": "longitude", "y": "latitude"}
    ds = ds.rename({k: v for k, v in noms.items() if k in ds.dims})
    ds = ds.reset_coords(drop=True)
    variables = [v for v in ds.data_vars if "time" in ds[v].dims]
    return ds[variables[0]]


# ---------------------------------------------------------------------------
for nom, fichier in FICHIERS.items():
    print(f"=== {nom} ===")
    amont = "07_controle_spatial_qc.py, 08_masque_artefact_lux.py, " if nom in VERSIONS_QC else ""
    da = ouvrir(ENTREE / fichier)
    unite = da.attrs.get("units", "")
    nom_long = da.attrs.get("long_name", nom)
    print(f"  variable : {da.name}   unites : {unite}")

    da = da.sel(time=slice(*PERIODE))

    mois_presents = set(np.unique(da.time.dt.month.values))
    inattendus = mois_presents - set(SAISON_DU_MOIS)
    if inattendus:
        print(f"  ATTENTION : mois inattendus dans les dates : {sorted(inattendus)}")

    for mois, saison in SAISON_DU_MOIS.items():
        sub = da.sel(time=da.time.dt.month == mois)
        annees = sub.time.dt.year.values
        dates = pd.to_datetime(pd.DataFrame({"year": annees, "month": mois, "day": 15}))
        sub = sub.assign_coords(time=dates).transpose("time", "latitude", "longitude")

        vide = sub.isnull().all(dim=("latitude", "longitude")).values
        if vide.any():
            print(f"  [{saison}] annees vides retirees : {list(annees[vide])}")
            sub = sub.isel(time=~vide)

        sub = sub.astype("float32")
        sub.attrs = {"units": unite, "long_name": nom_long, "grid_mapping": "crs"}
        sub.encoding = {"zlib": True, "complevel": 4, "_FillValue": np.float32(np.nan)}

        ds_out = sub.to_dataset(name=nom)
        ds_out["latitude"].attrs = {"standard_name": "latitude", "units": "degrees_north", "axis": "Y"}
        ds_out["longitude"].attrs = {"standard_name": "longitude", "units": "degrees_east", "axis": "X"}
        ds_out["time"].attrs = {"standard_name": "time", "axis": "T"}
        for c in ("latitude", "longitude"):
            ds_out[c].encoding = {}
        ds_out["time"].encoding = {"units": "days since 1950-01-01", "calendar": "standard"}
        ds_out["crs"] = xr.DataArray(np.int32(0), attrs=CRS_ATTRS)
        ds_out.attrs = {
            "Conventions": "CF-1.7",
            "title": f"E-OBS {nom} saisonnier Belgique - {saison}",
            "source": f"E-OBS v31.0e indices ETCCDI ({fichier})",
            "history": f"01_decoupe.py, {amont}11_export_eobs_indices.py (saison {saison})",
        }

        sortie = SORTIE / f"eobs_{nom}_BE_{saison}.nc"
        ds_out.to_netcdf(sortie)
        print(f"  [{saison}] ecrit : {sortie.name}  ({ds_out.sizes['time']} annees, "
              f"{str(ds_out.time.values[0])[:4]}-{str(ds_out.time.values[-1])[:4]})")

# ---------------------------------------------------------------------------
print("\n=== Controles (Bruxelles, 50.85 N 4.35 E) ===")
t = xr.open_dataset(SORTIE / "eobs_txx_BE_JJA.nc")
v = t.txx.sel(latitude=50.85, longitude=4.35, method="nearest").sel(time="2003-07-15")
print(f"  TXx ete 2003 : {float(v):.1f} degC  (attendu ~33-36)")
t.close()

f = xr.open_dataset(SORTIE / "eobs_fd_BE_DJF.nc")
v = f.fd.sel(latitude=50.85, longitude=4.35, method="nearest").sel(time="1963-01-15")
print(f"  Jours de gel hiver 1963 : {float(v):.0f} jours  (hiver tres froid, attendu eleve)")
f.close()

# ---------------------------------------------------------------------------
print("\n=== Comparaison versions d'origine / _qc ===")
for nom, fichier_qc in VERSIONS_QC.items():
    a = ouvrir(ENTREE / NORMALES[nom])
    b = ouvrir(ENTREE / fichier_qc)
    print(f"  {nom} :")
    print(f"    dimensions normale : {dict(a.sizes)}")
    print(f"    dimensions _qc     : {dict(b.sizes)}")
    print(f"    NaN normale / _qc  : {int(a.isnull().sum())} / {int(b.isnull().sum())}")
    print(f"    min-max normale    : {float(a.min()):.2f} / {float(a.max()):.2f}")
    print(f"    min-max _qc        : {float(b.min()):.2f} / {float(b.max()):.2f}")
    if a.sizes == b.sizes:
        diff = float(abs(a.values - b.values)[~np.isnan(a.values - b.values)].max())
        print(f"    ecart maximal      : {diff:.4f}")