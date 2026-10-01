# -*- coding: utf-8 -*-
"""
14_stats_zonales.py
A LANCER DANS LA FENETRE PYTHON D'ARCGIS PRO (Analysis > Python > Python Window).
Necessite l'extension Spatial Analyst.

Statistiques zonales par region NUTS 2 (Bruxelles-Capitale et les 10
provinces), a partir des CRF de 03_arcgis :
  - ERA5-Land : 4 saisons, 4 variables par fichier (t2m, tp, swvl1, swvl2)
  - E-OBS     : 10 combinaisons indice x saison (liste EOBS ci-dessous)
Les tables sont exportees en CSV dans 04_stats.
Le croisement NDVI x occupation du sol est traite par 15_ndvi_corine.py.

Prealable : CRF crees dans 03_arcgis a partir des NetCDF des scripts 09 a 11
(Copy Raster, extension .crf, Process as Multidimensional).
"""

import os
import arcpy
import arcpy.sa

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

# --- Chemins : a adapter ------------------------------------------------------
PROJET = r"C:\chemin\vers\Project Belgium"           # donnees et scripts Python
GDB_DATA = r"C:\chemin\vers\Belgium\Belgium.gdb"     # geodatabase du projet ArcGIS

ARCGIS = os.path.join(PROJET, "03_arcgis")
STATS = os.path.join(PROJET, "04_stats")
GDB = arcpy.mp.ArcGISProject("CURRENT").defaultGeodatabase   # sorties

# Provinces en WGS84 : meme systeme que les rasters, pas de reprojection inutile.
# Elles ne sont pas dans la geodatabase par defaut : on cherche dans les deux.
CHAMP_ZONE = "NUTS_ID"

PROVINCES = os.path.join(GDB_DATA, "provinces_BE_wgs84")
if not arcpy.Exists(PROVINCES):
    PROVINCES = os.path.join(GDB, "provinces_BE_wgs84")
if not arcpy.Exists(PROVINCES):
    raise SystemExit("provinces_BE_wgs84 introuvable : corriger GDB_DATA")
print("provinces : %s" % PROVINCES)
print("sorties   : %s" % GDB)

if not os.path.isdir(STATS):
    os.makedirs(STATS)

# Taille de cellule fine : sans cela, Bruxelles (BE10) est trop petite pour etre
# correctement echantillonnee sur une grille a 0.1 degre
CELLULE_FINE = 0.01

SAISONS = ["DJF", "MAM", "JJA", "SON"]


def exporter_csv(table, nom):
    """Exporte une table de la geodatabase en CSV dans 04_stats."""
    csv = os.path.join(STATS, nom + ".csv")
    if os.path.exists(csv):
        os.remove(csv)
    try:
        arcpy.conversion.ExportTable(table, csv)
    except AttributeError:
        arcpy.conversion.TableToTable(table, STATS, nom + ".csv")
    return csv


def zonal_provinces(raster, nom_sortie):
    """Statistiques zonales multidimensionnelles par province.

    arcpy.Raster(chemin, True) ouvre le CRF en mode multidimensionnel :
    plus fiable que MakeMultidimensionalRasterLayer, qui refuse variables="".
    """
    table = os.path.join(GDB, nom_sortie)
    if arcpy.Exists(table):
        print("  deja present : %s" % nom_sortie)
        return table

    md = arcpy.Raster(raster, True)
    with arcpy.EnvManager(cellSize=CELLULE_FINE):
        arcpy.sa.ZonalStatisticsAsTable(
            PROVINCES, CHAMP_ZONE, md, table,
            ignore_nodata="DATA", statistics_type="MEAN",
            process_as_multidimensional="ALL_SLICES",
        )
    print("  ok : %-22s %6s lignes" % (nom_sortie, arcpy.management.GetCount(table)[0]))
    return table


# ===========================================================================
# ERA5-Land et indices E-OBS par province
# ===========================================================================
print("=== ERA5-Land (4 variables ensemble) ===")
for s in SAISONS:
    nom = "zs_era5_%s" % s
    t = zonal_provinces(os.path.join(ARCGIS, "era5_BE_%s.crf" % s), nom)
    exporter_csv(t, nom)

print("\n=== Indices E-OBS ===")
EOBS = [
    ("txx", "JJA"), ("su", "JJA"), ("tr", "JJA"), ("su30", "JJA"),
    ("fd", "DJF"), ("fd", "MAM"),
    ("prcptot", "DJF"), ("prcptot", "MAM"), ("prcptot", "JJA"), ("prcptot", "SON"),
]
for indice, s in EOBS:
    nom = "zs_%s_%s" % (indice, s)
    t = zonal_provinces(os.path.join(ARCGIS, "eobs_%s_BE_%s.crf" % (indice, s)), nom)
    exporter_csv(t, nom)

print("\nTermine. CSV dans : %s" % STATS)
