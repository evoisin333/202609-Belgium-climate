"""
02_agregation.py
ERA5-Land : conversions d'unites (K -> degC ; tp, taux journalier moyen en
m/jour -> cumul mensuel en mm), puis agregation saisonniere (moyenne pour
t2m, swvl1 et swvl2 ; somme pour tp).
E-OBS : nombre de jours avec Tmax >= 30 degC par saison, calcule sur tx.

Convention DJF : decembre est rattache a l'hiver de l'annee suivante.
Periode conservee : 1961-2024 (hivers 1960 et 2025 incomplets).
Sorties (02_working) : era5_BE_saisonnier.nc, eobs_su30_BE_saisonnier.nc

A lancer depuis la racine du projet, environnement geo :
    python scripts/02_agregation.py
"""
import xarray as xr
import numpy as np
import os

OUT = "02_working"

# ---------------------------------------------------------------
# Attribution saison + annee-saison (DJF rattache a l'annee suivante)
# ---------------------------------------------------------------
def saison_annee(ds, axe_temps):
    mois = ds[axe_temps].dt.month
    annee = ds[axe_temps].dt.year

    saison = xr.where(mois.isin([12, 1, 2]), "DJF",
             xr.where(mois.isin([3, 4, 5]), "MAM",
             xr.where(mois.isin([6, 7, 8]), "JJA", "SON")))

    # decembre appartient a l'hiver de l'annee suivante
    annee_saison = xr.where(mois == 12, annee + 1, annee)

    return saison, annee_saison


# ===============================================================
# ETAPE 1 + 2 : ERA5 - conversions puis agregation saisonniere
# ===============================================================
print("=== ERA5 ===")
era = xr.open_dataset(os.path.join(OUT, "data_stream-moda_BE.nc"))

if "expver" in era.coords:
    era = era.drop_vars("expver")

# --- conversions ---
era["t2m"] = era.t2m - 273.15
era.t2m.attrs = {"units": "degC", "long_name": "Temperature a 2 m"}

era["swvl1"].attrs = {"units": "m3/m3", "long_name": "Humidite du sol 0-7 cm"}
era["swvl2"].attrs = {"units": "m3/m3", "long_name": "Humidite du sol 7-28 cm"}

# tp mensuel = taux journalier moyen en m/jour
# cumul mensuel en mm = tp * 1000 * nb_jours_du_mois
jours = era["valid_time"].dt.days_in_month
era["tp"] = era.tp * 1000.0 * jours
era.tp.attrs = {"units": "mm", "long_name": "Cumul de precipitations"}

print("  controle unites (juillet 2020) :")
juil = era.sel(valid_time="2020-07-01")
print(f"    t2m   {float(juil.t2m.mean()):6.1f} degC   (attendu ~18-20)")
print(f"    tp    {float(juil.tp.mean()):6.1f} mm     (attendu ~60-90)")
print(f"    swvl2 {float(juil.swvl2.mean()):6.3f} m3/m3 (attendu ~0.20-0.30)")

# --- agregation saisonniere ---
saison, annee_s = saison_annee(era, "valid_time")
era = era.assign_coords(saison=saison, annee=annee_s)

moyennes = era[["t2m", "swvl1", "swvl2"]].groupby("saison").map(
    lambda g: g.groupby("annee").mean("valid_time"))
sommes = era[["tp"]].groupby("saison").map(
    lambda g: g.groupby("annee").sum("valid_time", min_count=3))

era_seas = xr.merge([moyennes, sommes])

# hivers incomplets : 1960 (pas de dec 1959) et 2025 (inexistant)
era_seas = era_seas.sel(annee=slice(1961, 2024))

era_seas.to_netcdf(os.path.join(OUT, "era5_BE_saisonnier.nc"),
                   encoding={v: {"zlib": True, "complevel": 4}
                             for v in era_seas.data_vars})
print(f"  -> era5_BE_saisonnier.nc  {dict(era_seas.sizes)}")


# ===============================================================
# ETAPE 3 : jours de chaleur (Tmax >= 30 degC)
# ===============================================================
print("\n=== Jours de chaleur >= 30 degC ===")
tx = xr.open_dataset(os.path.join(OUT, "tx_ens_mean_0.1deg_reg_v31.0e_BE.nc"),
                     chunks={"time": 2000})

print(f"  unites declarees : {tx.tx.attrs.get('units', '?')}")

chaud = (tx.tx >= 30).astype("float32")
chaud = chaud.where(tx.tx.notnull())

saison, annee_s = saison_annee(tx, "time")
chaud = chaud.assign_coords(saison=saison, annee=annee_s)

su30 = chaud.groupby("saison").map(lambda g: g.groupby("annee").sum("time", min_count=1))
su30 = su30.sel(annee=slice(1961, 2024)).rename("su30")
su30.attrs = {"units": "jours", "long_name": "Jours avec Tmax >= 30 degC"}

su30 = su30.compute()
su30.to_netcdf(os.path.join(OUT, "eobs_su30_BE_saisonnier.nc"),
               encoding={"su30": {"zlib": True, "complevel": 4}})
print(f"  -> eobs_su30_BE_saisonnier.nc  {dict(su30.sizes)}")


# ===============================================================
# Controles finaux
# ===============================================================
print("\n=== Controles ===")
e = xr.open_dataset(os.path.join(OUT, "era5_BE_saisonnier.nc"))
for s in ["DJF", "MAM", "JJA", "SON"]:
    t = float(e.t2m.sel(saison=s).mean())
    p = float(e.tp.sel(saison=s).mean())
    print(f"  {s}  T = {t:5.1f} degC   P = {p:5.0f} mm")

print("\n  Jours >= 30 degC en JJA :")
j = su30.sel(saison="JJA").mean(dim=["latitude", "longitude"])
for per, lab in [(slice(1961, 1990), "1961-1990"),
                 (slice(1991, 2020), "1991-2020"),
                 (slice(2015, 2024), "2015-2024")]:
    print(f"    {lab} : {float(j.sel(annee=per).mean()):4.1f} jours/an")