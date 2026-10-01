# -*- coding: utf-8 -*-
"""
19_couche_signif_ndvi.py
A LANCER DANS LA FENETRE PYTHON D'ARCGIS PRO, carte NDVI active.

Cree la VRAIE couche des zones significatives du NDVI (p < 0.05), puis n'en
garde que les taches assez grandes pour etre lues comme des regions.

Pourquoi le filtre de surface :
  a 230 m, des milliers de pixels isoles passent le seuil de 5 % par hasard
  (environ 48 000 sur un million). Les entourer en rouge donnerait du bruit et
  contredirait la mise en garde de la StoryMap. Seules les taches etendues ne
  peuvent pas s'expliquer par le hasard.

Le script imprime la distribution des surfaces : ajuste SEUIL_KM2 ensuite si
le resultat est trop charge ou trop vide. Seuil retenu pour la StoryMap :
20 km2.
"""

import os

import arcpy
from arcpy.sa import Con, Raster

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

RACINE = r"C:\chemin\vers\Belgium"   # projet ArcGIS Pro, a adapter
TENDANCES = os.path.join(RACINE, "05_tendances")
GDB = os.path.join(RACINE, "Belgium.gdb")

# Adapte ce nom si besoin : le script liste le dossier s'il ne trouve pas.
PVALUE = os.path.join(TENDANCES, "trend_ndvi_JJA_pvalue.tif")

BRUT = os.path.join(GDB, "signif_ndvi_JJA_brut")
SORTIE = os.path.join(GDB, "signif_ndvi_JJA")
NOM_COUCHE = "Significant trends (p < 0.05)"

SEUIL_KM2 = 20.0          # surface minimale d'une tache conservee
KM2_PAR_DEGRE2 = 7900.0  # approximation a 50 deg N, suffisante pour un seuil

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
arcpy.conversion.RasterToPolygon(Con(Raster(PVALUE) < 0.05, 1), BRUT, "NO_SIMPLIFY")
total = int(arcpy.management.GetCount(BRUT)[0])
print("  %d polygones bruts" % total)

# ===========================================================================
# 2. Distribution des surfaces
# ===========================================================================
champ = "Shape_Area"
surfaces = sorted(r[0] * KM2_PAR_DEGRE2
                  for r in arcpy.da.SearchCursor(BRUT, [champ]))
if not surfaces:
    raise SystemExit("Aucun polygone : verifie que le raster de p-value est le bon.")


def centile(p):
    return surfaces[min(len(surfaces) - 1, int(len(surfaces) * p / 100.0))]


print("\n=== Surfaces des taches (km2) ===")
print("  mediane : %.3f" % centile(50))
print("  90e     : %.3f" % centile(90))
print("  99e     : %.3f" % centile(99))
print("  maximum : %.1f" % surfaces[-1])
print("  total significatif : %.0f km2" % sum(surfaces))

# ===========================================================================
# 3. Filtre
# ===========================================================================
seuil_deg2 = SEUIL_KM2 / KM2_PAR_DEGRE2
arcpy.analysis.Select(BRUT, SORTIE, "%s > %.10f" % (champ, seuil_deg2))
garde = int(arcpy.management.GetCount(SORTIE)[0])
retenu = sum(s for s in surfaces if s > SEUIL_KM2)
print("\n=== Filtre a %.1f km2 ===" % SEUIL_KM2)
print("  %d taches conservees sur %d (%.1f %%)" % (garde, total, 100.0 * garde / total))
print("  soit %.0f km2, %.0f %% de la surface significative totale" % (
    retenu, 100.0 * retenu / sum(surfaces)))
if garde > 60:
    print("  -> encore charge : augmente SEUIL_KM2 et relance.")
elif garde < 5:
    print("  -> presque vide : diminue SEUIL_KM2 et relance.")

# ===========================================================================
# 4. Ajout a la carte
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
print("  3. legende : 'Significant trends (p < 0.05, areas above %.0f km2)'" % SEUIL_KM2)
print("  4. note sous la carte : les taches isolees sous ce seuil ne sont pas tracees")
print("\nTermine.")