"""
05_diagnostic_pixel_max.py
Diagnostic, ne modifie aucun fichier : localise le maximum estival de su30,
recompte les jours >= 30 degC sur la Tmax journaliere de ce pixel et les
compare a l'indice SU officiel. C'est ce controle qui a revele l'artefact
E-OBS du sud du Grand-Duche de Luxembourg (55 jours >= 30 degC en 2003,
Tmax jusqu'a 47,2 degC).

A lancer depuis la racine du projet, environnement geo :
    python 05_diagnostic_pixel_max.py
"""
import xarray as xr

W = "02_working/"
j = xr.open_dataset(W + "eobs_su30_BE_saisonnier.nc").su30.sel(saison="JJA")
su = xr.open_dataset(W + "su_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE.nc",
                     decode_timedelta=False).suETCCDI
tx = xr.open_dataset(W + "tx_ens_mean_0.1deg_reg_v31.0e_BE.nc").tx

pos = j.argmax(dim=["annee", "latitude", "longitude"])
an = int(j.annee[pos["annee"]])
la = float(j.latitude[pos["latitude"]])
lo = float(j.longitude[pos["longitude"]])
print(f"max : {float(j.max()):.0f} jours en JJA {an}, lat {la:.2f} lon {lo:.2f}")

ete = slice(f"{an}-06-01", f"{an}-08-31")
serie = tx.sel(latitude=la, longitude=lo, method="nearest").sel(time=ete)
print(f"  recompte tx >= 30 : {int((serie >= 30).sum())} jours "
      f"sur {int(serie.notnull().sum())} valides")
print(f"  tx max de l'ete   : {float(serie.max()):.1f} degC")
su_pix = su.sel(latitude=la, longitude=lo, method="nearest").sel(time=ete)
print(f"  SU (>25) meme ete : {float(su_pix.mean()):.0f} jours")

print("\n5 etes avec le plus fort maximum :")
top = j.max(dim=["latitude", "longitude"]).to_series()
print(top.sort_values(ascending=False).head(5).round(0).astype(int).to_string())