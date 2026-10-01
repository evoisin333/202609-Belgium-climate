"""
01_decoupe.py
Decoupe tous les NetCDF europeens (ERA5-Land mensuel, indices saisonniers
E-OBS, Tmax journaliere E-OBS) sur le rectangle de la Belgique
(49,2-52,0 N ; 2,0-6,8 E). Les fichiers bruts sont lus a la racine du
projet ; les sorties *_BE.nc vont dans 02_working.

A lancer depuis la racine du projet, environnement geo :
    python 01_decoupe.py
"""
import xarray as xr
import os, glob

LAT_MIN, LAT_MAX = 49.2, 52.0
LON_MIN, LON_MAX = 2.0, 6.8

RAW = "."
OUT = "02_working"
os.makedirs(OUT, exist_ok=True)

GROS = "tx_ens_mean_0.1deg_reg_v31.0e.nc"

def decouper(chemin, chunks=None):
    nom = os.path.basename(chemin)
    ds = xr.open_dataset(chemin, chunks=chunks)

    if ds.latitude.values[0] > ds.latitude.values[-1]:
        tr_lat = slice(LAT_MAX, LAT_MIN)
    else:
        tr_lat = slice(LAT_MIN, LAT_MAX)

    be = ds.sel(latitude=tr_lat, longitude=slice(LON_MIN, LON_MAX))

    if be.latitude.size == 0 or be.longitude.size == 0:
        print(f"  ECHEC (emprise vide) : {nom}")
        return

    if "time_bnds" in be:
        be = be.drop_vars("time_bnds")

    dest = os.path.join(OUT, nom.replace(".nc", "_BE.nc"))
    enc = {v: {"zlib": True, "complevel": 4} for v in be.data_vars}
    be.to_netcdf(dest, encoding=enc)

    av = os.path.getsize(chemin) / 1e6
    ap = os.path.getsize(dest) / 1e6
    print(f"{nom[:45]:45s} {av:8.1f} Mo -> {ap:6.2f} Mo   "
          f"({be.latitude.size} x {be.longitude.size} px)")
    ds.close()

for f in sorted(glob.glob(os.path.join(RAW, "*.nc"))):
    if os.path.basename(f) == GROS:
        continue
    decouper(f)

print("\n--- fichier journalier (patience) ---")
decouper(os.path.join(RAW, GROS), chunks={"time": 500})