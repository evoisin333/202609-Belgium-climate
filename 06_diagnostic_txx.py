"""
06_diagnostic_txx.py
Diagnostic, ne modifie aucun fichier : liste les etes ou TXx depasse 40 degC,
situe le maximum su30 des etes les plus chauds et compare la serie TXx du
pixel suspect (49,55 N / 6,05 E) a celle d'Arlon (49,65 N / 5,85 E).
Resultat : 2019 et 2022 sont de vraies canicules, spatialement coherentes ;
2003 et 2004 portent la signature d'une station defectueuse.

A lancer depuis la racine du projet, environnement geo :
    python 06_diagnostic_txx.py
"""
import xarray as xr

W = "02_working/"
SEUIL = 40.0

txx = xr.open_dataset(W + "txx_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE.nc",
                      decode_timedelta=False).txxETCCDI
ete = txx.where(txx.time.dt.month.isin([6, 7, 8]), drop=True)

# 1. tous les etes ou TXx depasse le seuil
sus = ete.where(ete > SEUIL).to_dataframe(name="txx").dropna().reset_index()
sus["annee"] = sus["time"].dt.year
print(f"=== TXx estival > {SEUIL:.0f} degC ===")
print(f"  {len(sus)} cas, {sus[['latitude', 'longitude']].drop_duplicates().shape[0]} pixels")
print("\n  cas par annee :")
print(sus.groupby("annee").size().to_string())
print("\n  20 valeurs les plus fortes :")
print(sus.sort_values("txx", ascending=False)
         [["annee", "latitude", "longitude", "txx"]].head(20)
         .round(2).to_string(index=False))

# 2. ou se situe le maximum su30 des etes les plus forts
j = xr.open_dataset(W + "eobs_su30_BE_saisonnier.nc").su30.sel(saison="JJA")
print("\n=== Position du maximum su30 ===")
for an in [2003, 2022, 2004, 2006, 2015]:
    x = j.sel(annee=an)
    p = x.argmax(dim=["latitude", "longitude"])
    print(f"  {an} : {float(x.max()):3.0f} j  lat {float(x.latitude[p['latitude']]):.2f}"
          f"  lon {float(x.longitude[p['longitude']]):.2f}")

# 3. serie TXx du pixel suspect vs un pixel belge voisin (Arlon)
a = ete.sel(latitude=49.55, longitude=6.05, method="nearest")
b = ete.sel(latitude=49.65, longitude=5.85, method="nearest")
print("\n=== TXx estival : pixel suspect  |  Arlon ===")
for t, va, vb in zip(a.time.dt.year.values, a.values, b.values):
    if 1995 <= t <= 2012:
        print(f"  {t} : {float(va):5.1f}  |  {float(vb):5.1f}")