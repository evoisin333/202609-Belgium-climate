# -*- coding: utf-8 -*-
"""
13_tendances_arcgis.py

Finalisation des tendances Mann-Kendall produites par Generate Trend Raster.

Prealable, a la main dans ArcGIS Pro : Generate Trend Raster (methode
MANN-KENDALL) sur chaque CRF de 03_arcgis a analyser, sorties nommees
trend_<indicateur>_<saison>.crf dans le dossier "tendances" du projet ArcGIS.

Pour chaque trend_*.crf du dossier "tendances" :
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

import os

import arcpy
from arcpy.sa import Con, ExtractBand, SetNull

arcpy.CheckOutExtension("ImageAnalyst")
arcpy.CheckOutExtension("Spatial")

# --------------------------------------------------------------------------
# Chemins : deduits du projet ArcGIS Pro ouvert, rien a adapter
# --------------------------------------------------------------------------
RACINE = os.path.dirname(arcpy.mp.ArcGISProject("CURRENT").filePath)
TENDANCES = os.path.join(RACINE, "tendances")
SORTIE = os.path.join(RACINE, "05_tendances")
if not os.path.isdir(TENDANCES):
    raise SystemExit("Dossier 'tendances' introuvable dans %s" % RACINE)

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