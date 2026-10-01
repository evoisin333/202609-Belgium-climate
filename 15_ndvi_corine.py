# -*- coding: utf-8 -*-
"""
15_ndvi_corine.py
A LANCER DANS LA FENETRE PYTHON D'ARCGIS PRO. Necessite Spatial Analyst.

NDVI par province ET par classe d'occupation du sol (CORINE Land Cover 2018,
44 classes regroupees en 8), sur la grille du NDVI MODIS :
  1. decoupe du raster CORINE europeen complet sur les provinces, projetees
     d'abord dans le systeme de CORINE (ETRS89-LAEA)
  2. reclassification en 8 classes
  3. zones province x classe, identifiant = province x 10 + classe
     (ex. 72 = province 7, classe 2) ; resultat obtenu : 86 zones, 11 provinces
  4. statistiques zonales du NDVI pour les 4 saisons ; le code NUTS et le
     libelle de classe sont ecrits dans les tables de sortie (la table
     attributaire du raster de zones refuse d'etre modifiee)
Sorties (04_stats) : zs_ndvi_<saison>.csv, corine_8classes_libelles.csv

Deux pieges rencontres :
  - un premier decoupage de CORINE ne couvrait qu'environ 80 x 95 km autour
    de Namur : la decoupe repart donc du raster europeen complet ;
  - reprojection de CORINE, rasterisation des provinces et croisement doivent
    partager le meme cadre de reference (snapRaster, taille de cellule,
    emprise), sinon le croisement ne couvre qu'une partie du pays.

Prealable : CRF modis_ndvi_BE_<saison>.crf dans 03_arcgis.
"""

import os
import glob
import arcpy
from arcpy.sa import Reclassify, RemapValue, Int

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

# --- Chemins : a adapter ------------------------------------------------------
PROJET = r"C:\chemin\vers\Project Belgium"           # donnees et scripts Python
GDB_DATA = r"C:\chemin\vers\Belgium\Belgium.gdb"     # geodatabase du projet ArcGIS
# Raster CORINE europeen complet (U2018_CLC2018*.tif). Laisser None pour le
# chercher dans les cartes du projet ouvert, puis dans les dossiers voisins.
CLC_FULL = None

ARCGIS = os.path.join(PROJET, "03_arcgis")
STATS = os.path.join(PROJET, "04_stats")
GDB = arcpy.mp.ArcGISProject("CURRENT").defaultGeodatabase

PROVINCES = os.path.join(GDB_DATA, "provinces_BE_wgs84")
CHAMP_ZONE = "NUTS_ID"
NDVI_REF = os.path.join(ARCGIS, "modis_ndvi_BE_JJA.crf")

CLC_DECOUPE = os.path.join(PROJET, "01_raw", "CLC2018_BE2.tif")
CORINE_RECLASS = os.path.join(GDB, "corine_8classes")
CORINE_GRILLE = os.path.join(GDB, "corine_grille_ndvi")
PROV_GRILLE = os.path.join(GDB, "provinces_grille_ndvi")
ZONES = os.path.join(GDB, "zones_prov_corine")

SAISONS = ["DJF", "MAM", "JJA", "SON"]

LIBELLES = {
    1: "Artificialise", 2: "Cultures", 3: "Prairies",
    4: "Foret feuillus", 5: "Foret resineux", 6: "Foret mixte",
    7: "Milieux semi-naturels", 8: "Zones humides et eau",
}

# Codes CORINE a 3 chiffres -> groupe (l'ordre donne les codes de grille 1-44)
CODES_3 = {}
for c in [111, 112, 121, 122, 123, 124, 131, 132, 133, 141, 142]:
    CODES_3[c] = 1
for c in [211, 212, 213, 221, 222, 223, 241, 242, 243, 244]:
    CODES_3[c] = 2
CODES_3[231] = 3
CODES_3[311] = 4
CODES_3[312] = 5
CODES_3[313] = 6
for c in [321, 322, 323, 324, 331, 332, 333, 334, 335]:
    CODES_3[c] = 7
for c in [411, 412, 421, 422, 423, 511, 512, 521, 522, 523]:
    CODES_3[c] = 8
CODES_GRILLE = dict((i + 1, CODES_3[c]) for i, c in enumerate(sorted(CODES_3)))

os.makedirs(STATS, exist_ok=True)


def exporter_csv(table, nom):
    csv = os.path.join(STATS, nom + ".csv")
    if os.path.exists(csv):
        os.remove(csv)
    try:
        arcpy.conversion.ExportTable(table, csv)
    except AttributeError:
        arcpy.conversion.TableToTable(table, STATS, nom + ".csv")


# --- 0. Retrouver le raster CORINE europeen complet --------------------------
if CLC_FULL is None:
    projet = arcpy.mp.ArcGISProject("CURRENT")
    for m in projet.listMaps():
        for l in m.listLayers():
            if l.isRasterLayer and l.name.upper().startswith("U2018_CLC2018"):
                CLC_FULL = l.dataSource
                break
if CLC_FULL is None:
    motif = os.path.join(os.path.dirname(PROJET), "**", "U2018_CLC2018*.tif")
    trouves = glob.glob(motif, recursive=True)
    CLC_FULL = trouves[0] if trouves else None
if CLC_FULL is None:
    raise SystemExit("Raster CORINE europeen introuvable : renseigner CLC_FULL a la main")
print("CORINE complet : %s" % CLC_FULL)

# --- 1. Nettoyage des sorties precedentes ------------------------------------
print("\n=== Suppression des sorties precedentes ===")
for c in [ZONES, CORINE_GRILLE, PROV_GRILLE, CORINE_RECLASS] + \
         [os.path.join(GDB, "zs_ndvi_%s" % s) for s in SAISONS]:
    if arcpy.Exists(c):
        arcpy.management.Delete(c)
        print("  supprime : %s" % os.path.basename(c))

# --- 2. Decoupage sur la Belgique --------------------------------------------
print("\n=== Decoupage de CORINE sur la Belgique ===")
sr_clc = arcpy.Describe(CLC_FULL).spatialReference
prov_laea = os.path.join(GDB, "prov_laea_tmp")
if arcpy.Exists(prov_laea):
    arcpy.management.Delete(prov_laea)
arcpy.management.Project(PROVINCES, prov_laea, sr_clc)
e = arcpy.Describe(prov_laea).extent
print("  emprise des provinces en %s :" % sr_clc.name)
print("    X %.0f -> %.0f  (%.0f km)" % (e.XMin, e.XMax, (e.XMax - e.XMin) / 1000))
print("    Y %.0f -> %.0f  (%.0f km)" % (e.YMin, e.YMax, (e.YMax - e.YMin) / 1000))

arcpy.management.Clip(CLC_FULL, "#", CLC_DECOUPE, prov_laea,
                      "255", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
print("  ok : %s" % os.path.basename(CLC_DECOUPE))

# --- 3. Reclassification ------------------------------------------------------
print("\n=== Reclassification en 8 classes ===")
r = arcpy.Raster(CLC_DECOUPE)
if r.maximum is None:
    arcpy.management.CalculateStatistics(CLC_DECOUPE)
    r = arcpy.Raster(CLC_DECOUPE)
print("  valeur maximale : %s" % r.maximum)
# Le raster CORINE 100 m stocke souvent un code de grille 1-44 au lieu des
# codes a 3 chiffres ; tout code absent de la table (48 = hors zone,
# 255 = NoData) devient NoData
table_codes = CODES_3 if r.maximum > 100 else CODES_GRILLE
print("  -> %s" % ("codes a 3 chiffres" if r.maximum > 100 else "codes de grille 1-44"))

remap = RemapValue([[v, table_codes[v]] for v in sorted(table_codes)])
Reclassify(CLC_DECOUPE, "Value", remap, "NODATA").save(CORINE_RECLASS)
print("  ok : corine_8classes")

# --- 4. Zones province x classe sur la grille NDVI ---------------------------
desc = arcpy.Describe(NDVI_REF)
sr = desc.spatialReference
taille = desc.meanCellWidth
etendue = arcpy.Describe(PROVINCES).extent

# Les trois operations partagent le meme cadre de reference
with arcpy.EnvManager(snapRaster=NDVI_REF, cellSize=taille,
                      extent=etendue, outputCoordinateSystem=sr):
    print("\n=== Reprojection de CORINE sur la grille NDVI ===")
    # NEAREST obligatoire : donnee categorielle
    arcpy.management.ProjectRaster(CORINE_RECLASS, CORINE_GRILLE, sr,
                                   "NEAREST", str(taille))
    e2 = arcpy.Raster(CORINE_GRILLE).extent
    print("  emprise : X %.3f-%.3f  Y %.3f-%.3f" % (e2.XMin, e2.XMax, e2.YMin, e2.YMax))

    print("\n=== Rasterisation des provinces ===")
    arcpy.conversion.PolygonToRaster(PROVINCES, CHAMP_ZONE, PROV_GRILLE,
                                     "CELL_CENTER", cellsize=str(taille))

    print("\n=== Croisement ===")
    # Combine echoue a l'ecriture en geodatabase : identifiant de zone calcule
    # directement, province x 10 + classe (provinces 1-11, classes 1-8)
    zones_r = Int(arcpy.Raster(PROV_GRILLE) * 10 + arcpy.Raster(CORINE_GRILLE))
    zones_r.save(ZONES)

print("  zones creees : %s  (attendu 86)" % arcpy.management.GetCount(ZONES)[0])

# --- 5. Correspondance code numerique -> code NUTS ---------------------------
correspondance = dict((int(x[0]), x[1]) for x in
                      arcpy.da.SearchCursor(PROV_GRILLE, ["Value", CHAMP_ZONE]))
print("\n=== Provinces ===")
for k in sorted(correspondance):
    print("  %2d -> %s" % (k, correspondance[k]))

zones_vues = sorted(int(x[0]) for x in arcpy.da.SearchCursor(ZONES, ["Value"]))
provinces = sorted(set(correspondance.get(z // 10, "?") for z in zones_vues))
print("\n=== Zones ===")
print("  %d zones, %d provinces : %s" % (len(zones_vues), len(provinces), provinces))

# --- 6. Statistiques zonales du NDVI -----------------------------------------
print("\n=== NDVI par province et classe d'occupation du sol ===")
for s in SAISONS:
    nom = "zs_ndvi_%s" % s
    table = os.path.join(GDB, nom)
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

    md = arcpy.Raster(os.path.join(ARCGIS, "modis_ndvi_BE_%s.crf" % s), True)
    arcpy.sa.ZonalStatisticsAsTable(
        ZONES, "Value", md, table,
        ignore_nodata="DATA", statistics_type="MEAN",
        process_as_multidimensional="ALL_SLICES",
    )

    # Libelles ajoutes a la table de sortie, pas au raster
    arcpy.management.AddField(table, CHAMP_ZONE, "TEXT", field_length=5)
    arcpy.management.AddField(table, "CLC_CODE", "SHORT")
    arcpy.management.AddField(table, "CLASSE", "TEXT", field_length=30)
    with arcpy.da.UpdateCursor(table, ["Value", CHAMP_ZONE, "CLC_CODE", "CLASSE"]) as cur:
        for row in cur:
            p, c = divmod(int(row[0]), 10)
            row[1] = correspondance.get(p, "?")
            row[2] = c
            row[3] = LIBELLES.get(c, "?")
            cur.updateRow(row)

    n = int(arcpy.management.GetCount(table)[0])
    print("  ok : %-14s %5d lignes  (%d zones x 24 ans)" % (nom, n, n / 24))
    exporter_csv(table, nom)

# --- Table de correspondance des classes, pour le depot ----------------------
with open(os.path.join(STATS, "corine_8classes_libelles.csv"), "w") as f:
    f.write("code,libelle\n")
    for k in sorted(LIBELLES):
        f.write("%d,%s\n" % (k, LIBELLES[k]))

# --- Controle de plausibilite ------------------------------------------------
print("\n=== Controle : NDVI ete 2003 par classe (moyenne des provinces) ===")
table = os.path.join(GDB, "zs_ndvi_JJA")
champs = [f.name for f in arcpy.ListFields(table)]
champ_date = [c for c in champs if c.upper() in ("STDTIME", "DIMENSIONVALUE")]
champ_date = champ_date[0] if champ_date else None

cumul = {}
with arcpy.da.SearchCursor(table, ["CLASSE", "MEAN"] + ([champ_date] if champ_date else [])) as cur:
    for row in cur:
        if champ_date and str(row[2])[:4] != "2003":
            continue
        cumul.setdefault(row[0], []).append(row[1])
for k in sorted(cumul):
    v = cumul[k]
    print("  %-22s %.3f  (%d provinces)" % (k, sum(v) / len(v), len(v)))

print("\nTermine. CSV dans : %s" % STATS)
