"""
04_controles.py
Relit les sorties de 02_agregation.py et 03_modis_ndvi.py et verifie leur
coherence avant l'import dans ArcGIS Pro : pixels terre et mer, faux zeros,
plages de valeurs, couverture du NDVI, coherence entre su30 et SU.
Ne modifie AUCUN fichier.

A lancer depuis la racine du projet, environnement geo :
    python scripts/04_controles.py
"""
import glob
import os
import xarray as xr

OUT = "02_working"
F_ERA5 = os.path.join(OUT, "era5_BE_saisonnier.nc")
F_SU30 = os.path.join(OUT, "eobs_su30_BE_saisonnier.nc")
F_TX = os.path.join(OUT, "tx_ens_mean_0.1deg_reg_v31.0e_BE.nc")
F_MODIS = os.path.join(OUT, "modis_ndvi_BE_saisonnier.nc")
SAISONS = ["DJF", "MAM", "JJA", "SON"]
PERIODES = [(slice(1961, 1990), "1961-1990"),
            (slice(1991, 2020), "1991-2020"),
            (slice(2015, 2024), "2015-2024")]


def titre(t):
    print(f"\n{'=' * 64}\n{t}\n{'=' * 64}")


def bloc(nom, fonction):
    """Execute un bloc de controle sans arreter le script en cas d'erreur."""
    try:
        fonction()
    except Exception as err:
        print(f"  !! bloc '{nom}' interrompu : {type(err).__name__} : {err}")


# ================================================================
# 1. ERA5-Land
# ================================================================
def controle_era5():
    titre("1. ERA5-Land saisonnier")
    e = xr.open_dataset(F_ERA5)
    print("  dimensions :", dict(e.sizes))
    print(f"  annees     : {int(e.annee.min())} -> {int(e.annee.max())}")

    # un pixel est 'terre' si la temperature existe partout
    terre = e.t2m.notnull().all(dim=["saison", "annee"])
    n_terre, n_tot = int(terre.sum()), int(terre.size)
    print(f"  pixels terre : {n_terre} / {n_tot} ({100 * n_terre / n_tot:.0f} %)")

    faux = ((e.tp == 0).all(dim=["saison", "annee"]) & ~terre)
    print(f"  pixels MER avec tp = 0 mm (faux zeros) : {int(faux.sum())}"
          "   <- doit etre 0")

    print("\n  saison |   T    | P rectangle | P terre seule | swvl1 | swvl2")
    tot_rect = tot_terre = 0.0
    for s in SAISONS:
        x = e.sel(saison=s)
        t = float(x.t2m.mean())
        pr = float(x.tp.mean())
        pt = float(x.tp.where(terre).mean())
        s1 = float(x.swvl1.where(terre).mean())
        s2 = float(x.swvl2.where(terre).mean())
        tot_rect += pr
        tot_terre += pt
        print(f"   {s}  | {t:5.1f}  |  {pr:6.0f} mm  |   {pt:6.0f} mm   | "
              f"{s1:.3f} | {s2:.3f}")
    print(f"  cumul annuel : rectangle {tot_rect:.0f} mm  |  "
          f"terre seule {tot_terre:.0f} mm")

    print("\n  plages de valeurs (pixels terre) :")
    for v, bas, haut in [("t2m", -15, 30), ("tp", 0, 1000),
                         ("swvl1", 0.0, 0.8), ("swvl2", 0.0, 0.8)]:
        da = e[v].where(terre)
        mn, mx = float(da.min()), float(da.max())
        ok = "OK" if (bas <= mn and mx <= haut) else "A VERIFIER"
        print(f"    {v:6s} min {mn:8.3f}  max {mx:8.3f}   {ok}")


# ================================================================
# 2. Jours Tmax >= 30 degC
# ================================================================
def controle_su30():
    titre("2. Jours Tmax >= 30 degC (E-OBS, calcul maison)")
    j = xr.open_dataset(F_SU30).su30
    print("  dimensions :", dict(j.sizes))
    print(f"  annees     : {int(j.annee.min())} -> {int(j.annee.max())}")

    tx = xr.open_dataset(F_TX, chunks={"time": 1000}).tx
    terre = tx.sel(time="2020").notnull().any("time").compute()
    n_terre, n_tot = int(terre.sum()), int(terre.size)
    print(f"  pixels terre : {n_terre} / {n_tot} ({100 * n_terre / n_tot:.0f} %)")

    faux = (j == 0).all(dim=["saison", "annee"]) & ~terre
    print(f"  pixels MER avec 0 jour (faux zeros) : {int(faux.sum())}"
          "   <- doit etre 0")

    print(f"  max JJA : {float(j.sel(saison='JJA').max()):.0f} jours (<= 92)")
    print(f"  max DJF : {float(j.sel(saison='DJF').max()):.0f} jours (attendu 0)")

    print("\n  JJA, jours/an |  rectangle | terre seule")
    jja = j.sel(saison="JJA")
    for per, lab in PERIODES:
        r = float(jja.sel(annee=per).mean())
        t = float(jja.where(terre).sel(annee=per).mean())
        print(f"    {lab}  |    {r:4.1f}    |    {t:4.1f}")


# ================================================================
# 3. MODIS NDVI
# ================================================================
def controle_modis():
    titre("3. MODIS NDVI saisonnier")
    m = xr.open_dataset(F_MODIS).ndvi
    print("  dimensions :", dict(m.sizes))
    annees = [int(a) for a in m.annee.values]
    print(f"  annees ({len(annees)}) : {annees}")

    mn, mx = float(m.min()), float(m.max())
    print(f"  NDVI min {mn:.3f}  max {mx:.3f}   "
          f"{'OK' if -0.2 <= mn and mx <= 1.0 else 'A VERIFIER (echelle ?)'}")

    extremes = annees[:2] + annees[-2:]
    print("\n  part de pixels valides, annees extremes vs moyenne :")
    for s in SAISONS:
        x = m.sel(saison=s)
        detail = "  ".join(f"{a}:{100 * float(x.sel(annee=a).notnull().mean()):3.0f}%"
                           for a in extremes)
        moy = 100 * float(x.notnull().mean())
        print(f"    {s} | {detail} | moyenne {moy:3.0f}%")


# ================================================================
# 4. Indices E-OBS pre-calcules
# ================================================================
def controle_indices():
    titre("4. Indices E-OBS pre-calcules")
    connus = {os.path.abspath(f) for f in (F_ERA5, F_SU30, F_TX, F_MODIS)}
    su = fd = None
    for f in sorted(glob.glob(os.path.join(OUT, "*.nc"))):
        if os.path.abspath(f) in connus:
            continue
        ds = xr.open_dataset(f, decode_timedelta=False)
        if "time" not in ds.dims or ds.sizes["time"] > 400:
            continue
        t = ds.time.values
        print(f"  {os.path.basename(f)}")
        print(f"     {ds.sizes['time']} pas : {str(t[0])[:10]} -> {str(t[-1])[:10]}"
              f"   variables {list(ds.data_vars)}")
        if "suETCCDI" in ds:
            su = ds.suETCCDI
        if "fdETCCDI" in ds:
            fd = ds.fdETCCDI

    if fd is not None:
        f9 = fd.sel(time=slice("1991", "2020"))
        jan = float(f9.where(f9.time.dt.month.isin([12, 1, 2])).mean())
        jul = float(f9.where(f9.time.dt.month.isin([6, 7, 8])).mean())
        print(f"\n  FD 1991-2020 : pas d'hiver {jan:.1f} j | pas d'ete {jul:.1f} j"
              "   (hiver eleve ; ete ~0)")

    if su is not None:
        s9 = su.sel(time=slice("1991", "2020"))
        su_jja = float(s9.where(s9.time.dt.month.isin([6, 7, 8])).mean())
        j = xr.open_dataset(F_SU30).su30.sel(saison="JJA", annee=slice(1991, 2020))
        tx = xr.open_dataset(F_TX, chunks={"time": 1000}).tx
        terre = tx.sel(time="2020").notnull().any("time").compute()
        su30 = float(j.where(terre).mean())
        print(f"  JJA 1991-2020 : SU (>25 degC) {su_jja:.1f} j | su30 terre {su30:.1f} j"
              "   (su30 doit etre nettement < SU)")


if __name__ == "__main__":
    bloc("ERA5", controle_era5)
    bloc("su30", controle_su30)
    bloc("MODIS", controle_modis)
    bloc("indices", controle_indices)
    print("\nFin des controles.")