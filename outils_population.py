# outils_population.py — répartition de la population sur les unités d'adresse (module sans numéro)
# Règle fixée à l'étape 2.4 (26/11/2026), écrite UNE fois ici pour tous les scripts :
#   1. unité = adresse sans boîte OU boîte (registre UrbIS – Adresses) ;
#   2. poids d'une unité = population de sa maille (grille Statbel 2025) / nombre d'unités de la maille ;
#      unité hors de toute maille habitée -> poids 0 ;
#   3. maille habitée sans unité -> un pseudo-point à son centre, de poids = population de la maille ;
#   4. dans chaque secteur 2026, les poids sont remis à l'échelle pour que leur somme = population Statbel 2026 ;
#      secteur habité dont tous les poids sont nuls -> parts égales entre ses unités, ou point intérieur
#      du secteur s'il n'en a aucune ;
#   5. secteurs « inconnus » (ZZZZ) : non placés.

import numpy as np
import pandas as pd
import geopandas as gpd

CRS_GRILLE = "EPSG:3035"


def charger_secteurs(config):
    """Secteurs 2026 de la Région (EPSG:31370), code en texte."""
    fichier = next(p for p in (config.DATA_BRUT / "statbel" / "secteurs_2026").rglob("*.geojson") if p.is_file())
    s = gpd.read_file(fichier, columns=["cd_sector", "cd_rgn_refnis"])
    s["cd_rgn_refnis"] = s["cd_rgn_refnis"].astype(str).str.zfill(5)
    s["cd_sector"] = s["cd_sector"].astype(str)
    return s[s["cd_rgn_refnis"] == "04000"].to_crs(config.CRS_CALCUL)[["cd_sector", "geometry"]].reset_index(drop=True)


def charger_population(config):
    """Population Statbel 2026 par secteur bruxellois (séparateur |), codes en texte."""
    fichier = config.DATA_BRUT / "statbel" / "population_2026" / "OPENDATA_SECTOREN_2026.txt"
    pop = pd.read_csv(fichier, sep="|", encoding="utf-8-sig", dtype={"CD_REFNIS": str, "CD_SECTOR": str})
    pop["CD_REFNIS"] = pop["CD_REFNIS"].str.zfill(5)
    pop = pop[pop["CD_REFNIS"].str.startswith("21")]
    return pop.set_index("CD_SECTOR")["TOTAL"].astype(float)


def charger_unites(config):
    """Unités d'adresse (dernier fichier téléchargé), en points EPSG:31370."""
    fichier = sorted((config.DATA_BRUT / "urbis_adresses").glob("urbis_adresses_*.csv"))[-1]
    a = pd.read_csv(fichier, dtype=str)
    a_des_boites = set(a.loc[a["PARENTID"].notna(), "PARENTID"])
    unites = a[a["PARENTID"].notna() | ~a["INSPIRE_ID"].isin(a_des_boites)]
    geom = gpd.points_from_xy(unites["XL72"].astype(float), unites["YL72"].astype(float))
    return gpd.GeoDataFrame(unites[["INSPIRE_ID"]].reset_index(drop=True), geometry=geom, crs=config.CRS_CALCUL)


def rattacher_secteur(points, secteurs):
    """Code de secteur de chaque point : secteur qui le contient, sinon le plus proche (points sur une limite)."""
    lien = gpd.sjoin(points[["geometry"]], secteurs, predicate="within", how="left")
    lien = lien[~lien.index.duplicated()]["cd_sector"]
    manque = lien.isna()
    if manque.any():
        proche = gpd.sjoin_nearest(points.loc[manque, ["geometry"]], secteurs, how="left")
        lien.loc[manque] = proche[~proche.index.duplicated()]["cd_sector"]
    return lien, int(manque.sum())


def repartir_population(config):
    """Renvoie (points, diagnostic). points : GeoDataFrame EPSG:31370 avec cd_sector, type, pop.
    type ∈ {unite, pseudo_maille, part_egale, point_secteur}."""
    secteurs = charger_secteurs(config)
    pop_secteur = charger_population(config)
    region = secteurs.geometry.union_all()

    # Grille : toutes les mailles qui touchent la Région (une maille à cheval sur la limite garde sa population ;
    # la remise à l'échelle par secteur corrige ensuite)
    grille = gpd.read_file(config.DATA_BRUT / "statbel" / "grille_2025" / "POP_GRID_2025_3035.geojson",
                           bbox=secteurs, columns=["ms_len", "ms_pop"])
    grille = grille[grille.intersects(gpd.GeoSeries([region], crs=config.CRS_CALCUL).to_crs(grille.crs).iloc[0])]
    grille = grille.reset_index(drop=True)
    grille["id_maille"] = grille.index
    grille["ms_pop"] = grille["ms_pop"].astype(float)

    # 1-2. Unités et leur poids
    unites = charger_unites(config)
    unites = unites[unites.within(region)].reset_index(drop=True)
    lien = gpd.sjoin(unites.to_crs(grille.crs)[["geometry"]], grille[["id_maille", "geometry"]], predicate="within")
    unites["id_maille"] = lien[~lien.index.duplicated()]["id_maille"]
    n_par_maille = unites.groupby("id_maille").size()
    pop_maille = grille.set_index("id_maille")["ms_pop"]
    unites["poids"] = (unites["id_maille"].map(pop_maille) / unites["id_maille"].map(n_par_maille)).fillna(0.0)
    unites["type"] = "unite"

    # 3. Mailles habitées sans unité -> pseudo-point au centre (dans la Région seulement)
    vides = grille[(grille["ms_pop"] > 0) & ~grille["id_maille"].isin(n_par_maille.index)].copy()
    centres = gpd.GeoDataFrame(vides[["id_maille", "ms_pop"]], geometry=vides.geometry.centroid,
                               crs=grille.crs).to_crs(config.CRS_CALCUL)
    centres = centres[centres.within(region)]
    pseudo = gpd.GeoDataFrame({"INSPIRE_ID": None, "id_maille": centres["id_maille"].values,
                               "poids": centres["ms_pop"].values, "type": "pseudo_maille"},
                              geometry=centres.geometry.values, crs=config.CRS_CALCUL)

    points = pd.concat([unites, pseudo], ignore_index=True)
    points["cd_sector"], n_limite = rattacher_secteur(points, secteurs)

    # 4. Remise à l'échelle par secteur
    somme = points.groupby("cd_sector")["poids"].sum()
    facteur = points["cd_sector"].map(pop_secteur / somme.replace(0, np.nan))   # population 2026 / somme des poids
    points["pop"] = (points["poids"] * facteur).fillna(0.0)

    # Secteurs habités dont tous les poids sont nuls
    habites = pop_secteur[pop_secteur > 0]
    nuls = [c for c in habites.index if c in set(secteurs["cd_sector"]) and somme.get(c, 0) == 0]
    ajouts = []
    for code in nuls:
        dans = points["cd_sector"] == code
        if dans.any():
            points.loc[dans, "pop"] = habites[code] / dans.sum()
            points.loc[dans, "type"] = "part_egale"
        else:
            geom = secteurs.loc[secteurs["cd_sector"] == code, "geometry"].iloc[0].representative_point()
            ajouts.append({"INSPIRE_ID": None, "id_maille": np.nan, "poids": 0.0, "type": "point_secteur",
                           "cd_sector": code, "pop": habites[code], "geometry": geom})
    if ajouts:
        points = pd.concat([points, gpd.GeoDataFrame(ajouts, crs=config.CRS_CALCUL)], ignore_index=True)

    inconnus = pop_secteur[pop_secteur.index.str.contains("ZZZZ")]
    diagnostic = {
        "unites": len(unites), "unites_poids_nul": int((unites["poids"] == 0).sum()),
        "pseudo_maille": len(pseudo), "pop_pseudo_maille": float(pseudo["poids"].sum()),
        "points_rattaches_par_proximite": n_limite,
        "secteurs_poids_nuls": len(nuls), "zzzz_habitants": float(inconnus.sum()),
        "mailles": len(grille),
    }
    return gpd.GeoDataFrame(points, geometry="geometry", crs=config.CRS_CALCUL), diagnostic, pop_secteur, grille
