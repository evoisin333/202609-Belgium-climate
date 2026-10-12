# -*- coding: utf-8 -*-
"""
19_couche_signif_ndvi.py
A LANCER DANS LA FENETRE PYTHON D'ARCGIS PRO, carte NDVI active.

Cree la couche des zones ou la tendance du NDVI estival est significative,
d'apres 18_tendances.py : p-value 3PW ajustee par Benjamini-Hochberg
(trend_ndvi_JJA_pvalue_fdr.tif) inferieure a 0,10. La correction des tests
multiples remplace l'ancien filtre de surface (20 km2), qui n'avait pas de
fondement statistique.
"""

import os

import arcpy
from arcpy.sa import Con, Raster

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

# Chemins deduits du projet ArcGIS Pro ouvert, rien a adapter
RACINE = os.path.dirname(arcpy.mp.ArcGISProject("CURRENT").filePath)
TENDANCES = os.path.join(RACINE, "05_tendances")
GDB = os.path.join(RACINE, "Belgium.gdb")
if not os.path.isdir(TENDANCES):
    raise SystemExit("Dossier '05_tendances' introuvable dans %s" % RACINE)

# Adapte ce nom si besoin : le script liste le dossier s'il ne trouve pas.
PVALUE = os.path.join(TENDANCES, "trend_ndvi_JJA_pvalue_fdr.tif")

SORTIE = os.path.join(GDB, "signif_ndvi_JJA")
NOM_COUCHE = "Significant trends (FDR 10 %)"

SEUIL_FDR = 0.10          # p ajustee (Benjamini-Hochberg)
KM2_PAR_DEGRE2 = 7900.0   # approximation a 50 deg N, pour les surfaces affichees

if not arcpy.Exists(PVALUE):
    print("Introuvable : %s" % PVALUE)
    print("Rasters presents dans %s :" % TENDANCES)
    for f in sorted(os.listdir(TENDANCES)):
        if f.lower().endswith(".tif"):
            print("   ", f)
    raise SystemExit("Corrige la variable PVALUE puis relance.")

# ===========================================================================
# 1. Polygones des cellules significatives
# ===========================================================================
print("=== Conversion ===")
arcpy.conversion.RasterToPolygon(Con(Raster(PVALUE) < SEUIL_FDR, 1), SORTIE, "NO_SIMPLIFY")
total = int(arcpy.management.GetCount(SORTIE)[0])
print("  %d polygones" % total)

# ===========================================================================
# 2. Distribution des surfaces
# ===========================================================================
champ = "Shape_Area"
surfaces = sorted(r[0] * KM2_PAR_DEGRE2
                  for r in arcpy.da.SearchCursor(SORTIE, [champ]))
if not surfaces:
    raise SystemExit("Aucun pixel significatif apres correction FDR : pas de couche a tracer.")


def centile(p):
    return surfaces[min(len(surfaces) - 1, int(len(surfaces) * p / 100.0))]


print("\n=== Surfaces des taches (km2) ===")
print("  mediane : %.3f" % centile(50))
print("  90e     : %.3f" % centile(90))
print("  99e     : %.3f" % centile(99))
print("  maximum : %.1f" % surfaces[-1])
print("  total significatif : %.0f km2" % sum(surfaces))

# ===========================================================================
# 3. Ajout a la carte
# ===========================================================================
aprx = arcpy.mp.ArcGISProject("CURRENT")
m = aprx.activeMap
if m is None:
    raise SystemExit("Ouvre l'onglet de la carte NDVI, puis relance.")

for c in list(m.listLayers()):
    if c.name == NOM_COUCHE or "signif_ndvi" in c.name.lower():
        m.removeLayer(c)
lyr = m.addDataFromPath(SORTIE)
lyr.name = NOM_COUCHE
m.moveLayer(m.listLayers()[0], lyr, "BEFORE")
print("\n  couche ajoutee : %s" % NOM_COUCHE)

print("\n=== A finir a la main ===")
print("  1. symbologie : Fill -> No color, Outline -> rouge, 1 pt")
print("  2. supprimer l'ancienne couche rouge (copie des hachures)")
print("  3. legende : 'Significant trends (Mann-Kendall 3PW, FDR 10 %)'")
print("\nTermine.")