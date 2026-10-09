"""
00_telechargement_era5.py
Telecharge ERA5-Land, moyennes mensuelles (t2m, tp, swvl1, swvl2), sur le
rectangle de la Belgique, 1960-2024, via l'API du Climate Data Store (CDS).
Decembre 1960 sert a l'hiver DJF 1961.

Prealable :
  - compte CDS et licence du jeu acceptee (onglet Download de la page
    https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-monthly-means)
  - cle API dans ~/.cdsapirc (voir https://cds.climate.copernicus.eu/how-to-api)

Le CDS renvoie un zip qui contient data_stream-moda.nc ; il est extrait a la
racine du projet, ou 01_decoupe.py le lit.

A lancer depuis la racine du projet, environnement geo :
    python 00_telechargement_era5.py
"""

import os
import zipfile

import cdsapi

ZIP = "era5_land_monthly.zip"
ATTENDU = "data_stream-moda.nc"

if os.path.exists(ATTENDU):
    raise SystemExit("%s existe deja : le renommer pour retelecharger" % ATTENDU)

requete = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": [
        "2m_temperature",
        "total_precipitation",
        "volumetric_soil_water_layer_1",
        "volumetric_soil_water_layer_2",
    ],
    "year": [str(a) for a in range(1960, 2025)],
    "month": ["%02d" % m for m in range(1, 13)],
    "time": ["00:00"],
    "area": [52.0, 2.0, 49.2, 6.8],   # N, O, S, E
    "data_format": "netcdf",
    "download_format": "zip",
}

cdsapi.Client().retrieve("reanalysis-era5-land-monthly-means", requete).download(ZIP)

with zipfile.ZipFile(ZIP) as z:
    print("contenu du zip :", z.namelist())
    z.extractall(".")

if not os.path.exists(ATTENDU):
    raise SystemExit("%s absent du zip : verifier le contenu affiche ci-dessus" % ATTENDU)
print("ecrit : %s (%.1f Mo)" % (ATTENDU, os.path.getsize(ATTENDU) / 1e6))
