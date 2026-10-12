# -*- coding: utf-8 -*-
"""
13_conversion_crf.py

Copie chaque NetCDF de 03_arcgis (scripts 09 a 11) en CRF multidimensionnel
(Copy Raster, toutes les tranches). Ces CRF servent aux statistiques zonales
(14_stats_zonales.py et 15_ndvi_corine.py).

Les tendances ne sont plus calculees sous ArcGIS : voir 18_tendances.py
(Mann-Kendall corrige de l'autocorrelation, que Generate Trend Raster ne
propose pas).

Les CRF deja presents sont conserves : supprimer un fichier pour le recreer.

A LANCER DANS LA FENETRE PYTHON D'ARCGIS PRO (Analysis > Python > Python Window),
projet Belgium ouvert.
"""

import glob
import os

import arcpy

# Seule ligne a adapter : dossier du depot (scripts, 03_arcgis)
DEPOT = r"C:\chemin\vers\Project Belgium"

ARCGIS = os.path.join(DEPOT, "03_arcgis")
if not os.path.isdir(ARCGIS):
    raise SystemExit("03_arcgis introuvable dans %s : corriger DEPOT" % DEPOT)
print("depot : %s\n" % DEPOT)

print("=== NetCDF -> CRF ===")
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

print("\nTermine.")
