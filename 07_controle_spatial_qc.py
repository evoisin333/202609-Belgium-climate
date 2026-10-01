"""
07_controle_spatial_qc.py
Controle qualite des indices de chaleur E-OBS par coherence spatiale.

Une cellule (pixel x saison) est rejetee si son TXx depasse de plus de 4 degC
la mediane de ses voisins dans une fenetre de +/- 0,5 degre. Une vraie
canicule rechauffe tout le voisinage et passe le test ; une station
defectueuse cree un pic isole et le declenche.
Resultat : 40 cellules rejetees (sud du Grand-Duche de Luxembourg entre 2003
et 2020, Flandre en hiver 2022, un cas limite pres de Sedan en 2006) ; les
canicules de 2019 et 2022 sont conservees. TR et FD (Tmin) ne sont pas touches.

Sorties (02_working), les fichiers d'origine ne sont pas modifies :
  txx_..._BE_qc.nc, su_..._BE_qc.nc, eobs_su30_BE_saisonnier_qc.nc

A lancer depuis la racine du projet, environnement geo :
    python 07_controle_spatial_qc.py
"""
import numpy as np
import xarray as xr

W = "02_working/"
ECART = 4.0   # degC au-dessus de la mediane des voisins
FEN = 5       # rayon de la fenetre en pixels (5 -> +/- 0,5 deg)
SAISON = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
          6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}
F_TXX = "txx_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE.nc"
F_SU = "su_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE.nc"
F_TR = "tr_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE.nc"
F_FD = "fd_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE.nc"
F_SU30 = "eobs_su30_BE_saisonnier.nc"


def ouvrir(nom):
    return xr.open_dataset(W + nom, decode_timedelta=False).load()


# ---------------------------------------------------------------
# 1. Detection : TXx anormalement eleve par rapport aux voisins
# ---------------------------------------------------------------
txx_ds = ouvrir(F_TXX)
txx = txx_ds.txxETCCDI
mediane = txx.rolling(latitude=2 * FEN + 1, longitude=2 * FEN + 1,
                      center=True, min_periods=20).median()
drapeau = (txx - mediane) > ECART

df = xr.Dataset({"txx": txx, "voisins": mediane, "flag": drapeau}) \
       .to_dataframe().reset_index()
df = df[df["flag"] == True].copy()
df["annee_s"] = df["time"].dt.year
df["saison"] = df["time"].dt.month.map(SAISON)
df["ecart"] = df["txx"] - df["voisins"]

print(f"=== Pixels-saisons signales (TXx > mediane voisins + {ECART:.0f} degC) ===")
print(f"  {len(df)} cas")
print(df[["annee_s", "saison", "latitude", "longitude", "txx", "voisins", "ecart"]]
      .sort_values(["annee_s", "latitude", "longitude"])
      .round(2).to_string(index=False))

# ---------------------------------------------------------------
# 2. Controle Tmin au pixel suspect (TR, FD) : touches ou non ?
# ---------------------------------------------------------------
tr = ouvrir(F_TR).trETCCDI
fd = ouvrir(F_FD).fdETCCDI
print("\n=== TR (ete) et FD (hiver) : pixel suspect  |  Arlon ===")
for an in range(2001, 2007):
    def val(da, mois, la, lo):
        x = da.sel(latitude=la, longitude=lo, method="nearest")
        choix = (x.time.dt.year == an) & x.time.dt.month.isin(mois)
        return float(x.where(choix, drop=True).squeeze())
    print(f"  {an} : TR {val(tr, [6, 7, 8], 49.55, 6.05):4.0f} | {val(tr, [6, 7, 8], 49.65, 5.85):4.0f}"
          f"     FD {val(fd, [12, 1, 2], 49.55, 6.05):4.0f} | {val(fd, [12, 1, 2], 49.65, 5.85):4.0f}")

# ---------------------------------------------------------------
# 3. Masquage -> fichiers _qc (les originaux ne sont pas modifies)
# ---------------------------------------------------------------
masque = ~drapeau
txx_ds["txxETCCDI"] = txx.where(masque)
txx_ds.to_netcdf(W + F_TXX.replace(".nc", "_qc.nc"))

su_ds = ouvrir(F_SU)
su_ds["suETCCDI"] = su_ds.suETCCDI.where(masque.reindex_like(su_ds.suETCCDI, method="nearest", tolerance=0.01))
su_ds.to_netcdf(W + F_SU.replace(".nc", "_qc.nc"))

su30_ds = ouvrir(F_SU30)
su30 = su30_ds.su30
avant = float(su30.sel(saison="JJA", annee=slice(1991, 2020)).mean())
n_masque = 0
for r in df.itertuples():
    if r.annee_s not in su30.annee.values:
        continue
    i = su30.indexes["latitude"].get_indexer([r.latitude], method="nearest", tolerance=0.01)[0]
    k = su30.indexes["longitude"].get_indexer([r.longitude], method="nearest", tolerance=0.01)[0]
    if i < 0 or k < 0:
        continue
    su30.loc[dict(saison=r.saison, annee=r.annee_s,
                  latitude=su30.latitude[i], longitude=su30.longitude[k])] = np.nan
    n_masque += 1
su30_ds["su30"] = su30
su30_ds.to_netcdf(W + F_SU30.replace(".nc", "_qc.nc"))
apres = float(su30.sel(saison="JJA", annee=slice(1991, 2020)).mean())

print("\n=== Fichiers corriges ===")
print(f"  {F_TXX.replace('.nc', '_qc.nc')}")
print(f"  {F_SU.replace('.nc', '_qc.nc')}")
print(f"  {F_SU30.replace('.nc', '_qc.nc')}  ({n_masque} valeurs masquees)")
print(f"  su30 JJA 1991-2020 : {avant:.2f} -> {apres:.2f} jours/an")