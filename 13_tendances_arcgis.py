# -*- coding: utf-8 -*-
"""
13_tendances_arcgis.py

Tendances Mann-Kendall, de la conversion des NetCDF a la finalisation.

Etape A : chaque NetCDF de 03_arcgis (scripts 09 a 11) est copie en CRF
multidimensionnel (Copy Raster, toutes les tranches). Ces CRF servent aussi
aux scripts 14 et 15.

Etape B : Generate Trend Raster (MANN-KENDALL, dimension StdTime) sur les
combinaisons indicateur x saison de la liste TENDANCES_A_CALCULER, sorties
trend_<indicateur>_<saison>.crf dans le dossier "tendances" du projet ArcGIS.
Ignore NoData = DATA (valeur par defaut d'ArcGIS) : une annee manquante,
par exemple un ete masque par 08_masque_artefact_lux.py, est simplement
retiree de la serie du pixel au lieu de rendre le pixel NoData.

Les CRF et tendances deja presents sont conserves : supprimer un fichier
pour le recalculer.

Etape C, pour chaque trend_*.crf du dossier "tendances" :
  1. extraction de la bande Sens_Slope et de la bande P_Value ;
  2. conversion de la pente en unite par decennie (x10) ;
  3. version masquee ou les pixels non significatifs (p >= SEUIL_P) passent
     en NoData ;
  4. part de pixels significatifs PARMI LES PIXELS VALIDES.

Deux pieges evites :
  - un test du type Con((p < seuil) & (IsNull(p) == 0), 1, 0) force les
    pixels NoData a etre evalues et les compte comme non significatifs :
    le pourcentage obtenu n'est alors que la part de pixels terrestres.
    Con(p < seuil, 1, 0) laisse au contraire le NoData se propager.
  - les statistiques sont lues sur l'objet Raster, qui ignore le NoData,
    et non sur le fichier ecrit.

Sorties, dans 05_tendances :
  <nom>_pente_dec.tif         pente par decennie, tous pixels
  <nom>_pente_dec_signif.tif  idem, pixels non significatifs retires
  <nom>_pvalue.tif            p-value, pour hachurer les zones non significatives

A LANCER DANS LA FENETRE PYTHON D'ARCGIS PRO (Analysis > Python > Python Window),
projet Belgium ouvert.
"""

import glob
import os

import arcpy
from arcpy.ia import GenerateTrendRaster
from arcpy.sa import Con, ExtractBand, SetNull

arcpy.CheckOutExtension("ImageAnalyst")
arcpy.CheckOutExtension("Spatial")

# --------------------------------------------------------------------------
# Chemins
# --------------------------------------------------------------------------
# Seule ligne a adapter : dossier du depot (scripts, 03_arcgis)
DEPOT = r"C:\vers\chemin\Project Belgium"

# Le reste se deduit du projet ArcGIS Pro ouvert
RACINE = os.path.dirname(arcpy.mp.ArcGISProject("CURRENT").filePath)
ARCGIS = os.path.join(DEPOT, "03_arcgis")
TENDANCES = os.path.join(RACINE, "tendances")
SORTIE = os.path.join(RACINE, "05_tendances")
if not os.path.isdir(ARCGIS):
    raise SystemExit("03_arcgis introuvable dans %s : corriger DEPOT" % DEPOT)
os.makedirs(TENDANCES, exist_ok=True)
print("depot  : %s" % DEPOT)
print("projet : %s\n" % RACINE)

# Combinaisons analysees : (CRF source dans 03_arcgis, variable, sortie)
TENDANCES_A_CALCULER = [
    ("era5_BE_DJF.crf", "t2m", "trend_t2m_DJF"),
    ("era5_BE_MAM.crf", "t2m", "trend_t2m_MAM"),
    ("era5_BE_JJA.crf", "t2m", "trend_t2m_JJA"),
    ("era5_BE_SON.crf", "t2m", "trend_t2m_SON"),
    ("era5_BE_JJA.crf", "swvl2", "trend_swvl2_JJA"),
    ("eobs_fd_BE_DJF.crf", "fd", "trend_fd_DJF"),
    ("eobs_su_BE_JJA.crf", "su", "trend_su_JJA"),
    ("eobs_su30_BE_JJA.crf", "su30", "trend_su30_JJA"),
    ("eobs_txx_BE_JJA.crf", "txx", "trend_txx_JJA"),
    ("eobs_prcptot_BE_DJF.crf", "prcptot", "trend_prcptot_DJF"),
    ("eobs_prcptot_BE_MAM.crf", "prcptot", "trend_prcptot_MAM"),
]

# Parametres de Generate Trend Raster, ecrits en clair pour la reproductibilite
DIMENSION = "StdTime"         # nom ArcGIS de la dimension temps des NetCDF
METHODE = "MANN-KENDALL"
IGNORE_NODATA = "DATA"        # defaut ArcGIS : les annees NoData sont ignorees

# --------------------------------------------------------------------------
# Etape A : NetCDF -> CRF multidimensionnel
# --------------------------------------------------------------------------
print("=== A. NetCDF -> CRF ===")
for nc in sorted(glob.glob(os.path.join(ARCGIS, "*.nc"))):
    crf = nc[:-3] + ".crf"
    nom = os.path.basename(crf)
    if arcpy.Exists(crf):
        print("  deja present : %s" % nom)
        continue
    arcpy.management.CopyRaster(nc, crf, format="CRF",
                                process_as_multidimensional="ALL_SLICES")
    # Controle : le CRF doit garder toutes les variables du NetCDF
    variables = arcpy.Raster(crf, True).variableNames
    print("  cree : %s  (variables : %s)" % (nom, ", ".join(variables)))

# --------------------------------------------------------------------------
# Etape B : Generate Trend Raster
# --------------------------------------------------------------------------
print("\n=== B. Generate Trend Raster ===")
for source, variable, nom in TENDANCES_A_CALCULER:
    sortie_crf = os.path.join(TENDANCES, nom + ".crf")
    if arcpy.Exists(sortie_crf):
        print("  deja present : %s" % nom)
        continue
    entree = os.path.join(ARCGIS, source)
    if not arcpy.Exists(entree):
        print("  !! source introuvable, ignore : %s" % source)
        continue
    tendance = GenerateTrendRaster(entree, DIMENSION, variable, METHODE,
                                   ignore_nodata=IGNORE_NODATA)
    tendance.save(sortie_crf)
    print("  cree : %s" % nom)

# --------------------------------------------------------------------------
# Etape C : finalisation
# --------------------------------------------------------------------------
print("\n=== C. Finalisation ===")
SEUIL_P = 0.05
FACTEUR = 10.0        # par an -> par decennie

# Ordre des bandes de Generate Trend Raster (Mann-Kendall) :
# 1 Sens_Slope, 2 P_Value, 3 Score, 4 Score_Variance, 5 Z_Score
BANDE_PENTE = 1
BANDE_P = 2

os.makedirs(SORTIE, exist_ok=True)
arcpy.env.workspace = TENDANCES
arcpy.env.overwriteOutput = True
arcpy.env.outputCoordinateSystem = None
arcpy.env.pyramid = "NONE"
arcpy.env.rasterStatistics = "STATISTICS"


# --------------------------------------------------------------------------
def stats(raster):
    """min, moyenne, max d'un objet Raster, NoData exclu."""
    if raster.minimum is None:
        arcpy.management.CalculateStatistics(raster)
    return raster.minimum, raster.mean, raster.maximum


def enregistrer(raster, chemin):
    """Ecrit un raster, avec repli sur CopyRaster si .save echoue."""
    try:
        raster.save(chemin)
    except Exception as e:
        print(f"    .save a echoue ({e}), tentative via CopyRaster")
        arcpy.management.CopyRaster(raster, chemin, pixel_type="32_BIT_FLOAT")
    return chemin


# --------------------------------------------------------------------------
fichiers = sorted(arcpy.ListRasters("trend_*") or [])
if not fichiers:
    raise RuntimeError(f"Aucun raster trend_* trouve dans {TENDANCES}")

print(f"{len(fichiers)} tendances trouvees")
print(f"sorties : {SORTIE}\n")
resultats = []

for f in fichiers:
    nom = os.path.splitext(f)[0]           # trend_t2m_JJA
    chemin = os.path.join(TENDANCES, f)
    print(f"=== {nom} ===")

    try:
        # is_multidimensional=False : le CRF est lu comme un raster multibande,
        # ce qui rend les 5 statistiques accessibles bande par bande.
        r = arcpy.Raster(chemin, False)
        if r.bandCount < BANDE_P:
            print(f"  !! {r.bandCount} bande(s) seulement, ignore "
                  f"(attendu 5 : Sens_Slope, P_Value, Score, "
                  f"Score_Variance, Z_Score)")
            continue

        pente = ExtractBand(r, band_ids=[BANDE_PENTE]) * FACTEUR
        pval = ExtractBand(r, band_ids=[BANDE_P])

        # --- statistiques sur les pixels valides ---
        mini, moy, maxi = stats(pente)

        # Con sans IsNull : le NoData se propage, donc la moyenne de ce
        # raster binaire est bien la part de pixels VALIDES significatifs.
        signif = Con(pval < SEUIL_P, 1, 0)
        _, part, _ = stats(signif)
        pct = part * 100.0

        # --- ecriture ---
        enregistrer(pente, os.path.join(SORTIE, f"{nom}_pente_dec.tif"))
        enregistrer(pval, os.path.join(SORTIE, f"{nom}_pvalue.tif"))
        enregistrer(SetNull(pval >= SEUIL_P, pente),
                    os.path.join(SORTIE, f"{nom}_pente_dec_signif.tif"))

        print(f"  pente/decennie : min {mini:+.3f} | moyenne {moy:+.3f} "
              f"| max {maxi:+.3f}")
        print(f"  pixels significatifs (p < {SEUIL_P}) : {pct:.1f} %\n")
        resultats.append((nom, mini, moy, maxi, pct))

    except Exception as e:
        print(f"  !! echec sur {nom} : {e}\n")

# --------------------------------------------------------------------------
print("\n=== Recapitulatif ===")
print(f"{'tendance':28s} {'min':>9s} {'moyenne':>9s} {'max':>9s} {'signif.':>9s}")
for nom, mini, moy, maxi, pct in resultats:
    print(f"{nom:28s} {mini:+9.3f} {moy:+9.3f} {maxi:+9.3f} {pct:8.1f} %")

print("\nUnites : par decennie (degC, jours ou mm selon l'indicateur).")
print("Controle : des pourcentages identiques d'un indicateur a l'autre")
print("trahiraient un comptage des pixels NoData et non un vrai resultat.")
print("\nRappel pour la section Limites : le test de Mann-Kendall applique ici")
print("ne corrige pas l'autocorrelation temporelle (pas de pre-whitening), et")
print("il est repete sur plusieurs milliers de pixels sans correction pour")
print("tests multiples. Les deux points sont a documenter.")