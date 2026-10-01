# -*- coding: utf-8 -*-
"""
08_masque_artefact_lux.py

Masquage d'un artefact E-OBS v31.0e (station defectueuse, sud du Grand-Duche
de Luxembourg) sur les indices de chaleur estivaux : SU, TXx, jours >= 30 C.

Diagnostic (a reprendre dans les limites du rapport) :
- Ecart SU entre un pixel frontalier (49,65 N / 5,85 E, Arlon) et une colonne
  de reference a l'ouest (5,45 E) : +1,4 +/- 3,4 jours en temps normal,
  mais +9 a +17 jours en 2003, 2004, 2006 et 2007.
- ERA5-Land, source independante, ne montre aucun exces ces memes etes
  (< +0,1 C cote belge) -> artefact E-OBS confirme.

Methode :
- Exces moyen de SU (JJA) des etes suspects par rapport aux autres etes,
  apres retrait d'une colonne de reference de meme latitude (retire le
  signal climatique commun a toute la zone).
- Pixels dont l'exces atteint SEUIL jours (~3 ecarts-types du bruit residuel).
- Sur ces pixels, les etes suspects sont mis a NaN dans SU, TXx et su30.

A lancer depuis la racine du projet, environnement geo :
    python scripts/08_masque_artefact_lux.py

Fermer d'abord toute session Python (>>>) et retirer les couches d'ArcGIS Pro,
sinon les fichiers sont verrouilles en ecriture.
Ce masque s'applique AVANT les exports (09 a 11). S'il est applique apres :
relancer 10_export_ndvi_su30.py (su30) et 11_export_eobs_indices.py (SU, TXx),
puis refaire les CRF de SU, TXx et su30 dans ArcGIS (Copy Raster, extension
.crf, Process as Multidimensional).
"""

import glob
import os
import shutil

import xarray as xr

# --------------------------------------------------------------------------
# Parametres
# --------------------------------------------------------------------------
DOSSIER = "02_working"

FICHIER_SU = os.path.join(DOSSIER, "su_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE_qc.nc")
FICHIER_TXX = os.path.join(DOSSIER, "txx_seas_0.1deg_reg_ens_median_E-OBSv31.0e_BE_qc.nc")
MOTIFS_SU30 = [os.path.join(DOSSIER, "*su30*_qc.nc"),
               os.path.join(DOSSIER, "*su30*.nc")]

ANNEES_SUSPECTES = [2003, 2004, 2006, 2007]
SEUIL = 5.0                    # jours d'exces moyen
LON_REF = 5.45                 # colonne de reference, a l'ouest de la zone
LAT_MIN, LAT_MAX = 49.2, 50.5
LON_MIN, LON_MAX = 5.4, 6.5

FICHIER_MASQUE = os.path.join(DOSSIER, "masque_artefact_lux.nc")


# --------------------------------------------------------------------------
# Outils
# --------------------------------------------------------------------------
def etes_suspects(da):
    """Booleen vrai pour les etes (JJA) des annees suspectes.

    Gere les deux structures rencontrees dans le projet : un axe 'time'
    en dates de milieu de saison (E-OBS), ou des dimensions 'annee'/'saison'
    (sorties de 02_agregation.py).
    """
    if "time" in da.dims:
        t = da["time"]
        return (t.dt.month == 7) & t.dt.year.isin(ANNEES_SUSPECTES)
    if "annee" in da.dims:
        cond = da["annee"].isin(ANNEES_SUSPECTES)
        if "saison" in da.dims:
            cond = cond & (da["saison"] == "JJA")
        return cond
    raise ValueError("Dimension temporelle non reconnue : %s" % (da.dims,))


def copie_de_securite(chemin):
    copie = chemin.replace(".nc", "_avant_masque.nc")
    if os.path.exists(copie):
        print("  copie de securite deja presente : %s" % copie)
    else:
        shutil.copy2(chemin, copie)
        print("  copie de securite : %s" % copie)


def appliquer_masque(chemin, masque):
    print("\n=== %s ===" % os.path.basename(chemin))
    copie_de_securite(chemin)

    # Lecture complete en memoire puis fermeture : indispensable pour
    # pouvoir reecrire le fichier a la meme place.
    with xr.open_dataset(chemin) as src:
        ds = src.load()

    for nom, da in ds.data_vars.items():
        if not {"latitude", "longitude"}.issubset(set(da.dims)):
            continue

        # Masque 2D recale sur la grille de ce fichier. La tolerance absorbe
        # les petits ecarts d'arrondi entre fichiers E-OBS.
        autres = {d: 0 for d in da.dims if d not in ("latitude", "longitude")}
        gabarit = da.isel(autres, drop=True)
        m = (masque.astype("int8")
                   .reindex_like(gabarit, method="nearest", tolerance=0.01)
                   .fillna(0)
                   .astype(bool))

        a_masquer = etes_suspects(da) & m
        avant = int(da.notnull().sum())
        ds[nom] = da.where(~a_masquer)
        ds[nom].attrs = da.attrs
        apres = int(ds[nom].notnull().sum())

        print("  %-12s %5d valeurs -> NaN | maximum : %.2f avant, %.2f apres"
              % (nom, avant - apres, float(da.max()), float(ds[nom].max())))

    ds.attrs["masque_artefact_lux"] = (
        "Etes %s mis a NaN sur %d pixels (exces SU >= %.0f jours) - "
        "artefact station E-OBS, sud du Grand-Duche de Luxembourg. "
        "Voir 08_masque_artefact_lux.py"
        % (ANNEES_SUSPECTES, int(masque.sum()), SEUIL)
    )
    ds.to_netcdf(chemin)
    print("  reecrit : %s" % chemin)


# --------------------------------------------------------------------------
# 1. Calcul du masque a partir de SU
# --------------------------------------------------------------------------
print("=== Calcul du masque (source : SU) ===")

with xr.open_dataset(FICHIER_SU) as src:
    nom_su = [v for v in src.data_vars if "su" in v.lower()][0]
    su = src[nom_su].load()
print("  variable utilisee : %s" % nom_su)

jja = su.sel(time=su["time"].dt.month == 7)
ref = jja.sel(longitude=LON_REF, method="nearest").drop_vars("longitude")
zone = jja.sel(latitude=slice(LAT_MIN, LAT_MAX), longitude=slice(LON_MIN, LON_MAX))

diff = zone - ref
suspect = diff["time"].dt.year.isin(ANNEES_SUSPECTES)
exces = diff.where(suspect).mean("time") - diff.where(~suspect).mean("time")

masque_zone = (exces >= SEUIL)
# Etendu a la grille complete, faux partout ailleurs
masque = (masque_zone.astype("int8")
                     .reindex_like(su.isel(time=0, drop=True), fill_value=0)
                     .astype(bool))

print("  pixels retenus : %d (seuil %.0f jours)" % (int(masque.sum()), SEUIL))
print("\n  Exces moyen en jours (latitudes en lignes, longitudes en colonnes) :")
print(exces.round(0).to_pandas().to_string())

xr.Dataset({
    "exces_su_jours": exces,
    "masque": masque_zone.astype("int8"),
}).to_netcdf(FICHIER_MASQUE)
print("\n  masque enregistre : %s" % FICHIER_MASQUE)

# --------------------------------------------------------------------------
# 2. Application aux trois indices
# --------------------------------------------------------------------------
fichiers = [FICHIER_SU, FICHIER_TXX]

trouves = []
for motif in MOTIFS_SU30:
    trouves = sorted(f for f in glob.glob(motif) if "_avant_masque" not in f)
    if trouves:
        break

if trouves:
    fichiers += trouves
else:
    print("\n!! Aucun fichier su30 trouve dans %s. Verifie son nom et adapte "
          "MOTIFS_SU30, puis relance : les autres indices sont deja traites."
          % DOSSIER)

for f in fichiers:
    if os.path.exists(f):
        appliquer_masque(f, masque)
    else:
        print("\n!! Introuvable, ignore : %s" % f)

print("\nTermine.")
print("Etape suivante : les exports 09, 10 et 11. Si les CRF existent deja, "
      "relancer 10 et 11 puis refaire les CRF de SU, TXx et su30.")